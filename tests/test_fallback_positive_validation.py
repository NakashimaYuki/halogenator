# -*- coding: ascii -*-
"""Test positive validation for fallback paths (AtomMap and heuristic success)."""

import unittest
from unittest.mock import patch, MagicMock
from rdkit import Chem
from src.halogenator.enumerate_k import _detect_site_from_product_with_mapping_fallback, _detect_site_from_product


class TestFallbackPositiveValidation(unittest.TestCase):
    """Test positive coverage for fallback paths."""
    
    def setUp(self):
        """Set up test molecules."""
        # Simple benzene parent
        self.parent_mol = Chem.MolFromSmiles('c1ccccc1')
        if self.parent_mol:
            Chem.SanitizeMol(self.parent_mol)
    
    def test_atommap_success_positive(self):
        """Test AtomMap fallback success path."""
        # Create a product molecule with atom map number = 1 on the substitution site
        product_smiles = 'Fc1ccccc1'
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol:
            # Set atom map number on the carbon that's connected to F
            for atom in product_mol.GetAtoms():
                if atom.GetSymbol() == 'C':
                    # Find the carbon connected to fluorine
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'F':
                            atom.SetAtomMapNum(1)  # Mark this as the substitution site
                            break
                    if atom.GetAtomMapNum() == 1:
                        break
            
            Chem.SanitizeMol(product_mol)
        
        fallback_stats = {}
        
        # Call the fallback detection function
        site_atom_idx, ring_tag = _detect_site_from_product_with_mapping_fallback(
            self.parent_mol, product_mol, 'R1', 'F', fallback_stats
        )
        
        # Should have succeeded via AtomMap
        self.assertIsNotNone(site_atom_idx)
        self.assertEqual(fallback_stats.get('atommap_used', 0), 1)
        self.assertEqual(fallback_stats.get('heuristic_used', 0), 0)
    
    def test_heuristic_success_positive(self):
        """Test heuristic fallback success path."""
        # Create a product molecule with no atom map numbers (all 0)
        product_smiles = 'Fc1ccccc1'
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol:
            # Ensure all atoms have atom map number = 0 (this is usually the default)
            for atom in product_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            Chem.SanitizeMol(product_mol)
        
        fallback_stats = {}
        
        # Mock the heuristic detection to succeed
        with patch('src.halogenator.enumerate_k._detect_site_from_product') as mock_heuristic:
            mock_heuristic.return_value = (0, 'benzene')  # Return valid site and tag
            
            # Call the fallback detection function
            site_atom_idx, ring_tag = _detect_site_from_product_with_mapping_fallback(
                self.parent_mol, product_mol, 'R1', 'F', fallback_stats
            )
            
            # Should have succeeded via heuristic
            self.assertIsNotNone(site_atom_idx)
            self.assertEqual(site_atom_idx, 0)
            self.assertEqual(ring_tag, 'benzene')
            self.assertEqual(fallback_stats.get('atommap_used', 0), 0)
            self.assertEqual(fallback_stats.get('heuristic_used', 0), 1)
            
            # Verify heuristic was called
            mock_heuristic.assert_called_once_with(self.parent_mol, product_mol, 'R1', 'F')
    
    def test_atommap_vs_heuristic_mutual_exclusion(self):
        """Test that atommap_used and heuristic_used are mutually exclusive."""
        # Test case 1: AtomMap succeeds (should not call heuristic)
        product_mol_with_map = Chem.MolFromSmiles('Fc1ccccc1')
        if product_mol_with_map:
            # Set atom map on the substitution site
            for atom in product_mol_with_map.GetAtoms():
                if atom.GetSymbol() == 'C':
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'F':
                            atom.SetAtomMapNum(1)
                            break
                    if atom.GetAtomMapNum() == 1:
                        break
            Chem.SanitizeMol(product_mol_with_map)
        
        fallback_stats1 = {}
        
        with patch('src.halogenator.enumerate_k._detect_site_from_product') as mock_heuristic:
            _detect_site_from_product_with_mapping_fallback(
                self.parent_mol, product_mol_with_map, 'R1', 'F', fallback_stats1
            )
            
            # AtomMap should succeed, heuristic should not be called
            self.assertEqual(fallback_stats1.get('atommap_used', 0), 1)
            self.assertEqual(fallback_stats1.get('heuristic_used', 0), 0)
            mock_heuristic.assert_not_called()
        
        # Test case 2: AtomMap fails, heuristic succeeds
        product_mol_no_map = Chem.MolFromSmiles('Fc1ccccc1')
        if product_mol_no_map:
            for atom in product_mol_no_map.GetAtoms():
                atom.SetAtomMapNum(0)  # No atom maps
            Chem.SanitizeMol(product_mol_no_map)
        
        fallback_stats2 = {}
        
        with patch('src.halogenator.enumerate_k._detect_site_from_product') as mock_heuristic:
            mock_heuristic.return_value = (1, 'benzene')
            
            _detect_site_from_product_with_mapping_fallback(
                self.parent_mol, product_mol_no_map, 'R1', 'F', fallback_stats2
            )
            
            # Heuristic should succeed, AtomMap should not
            self.assertEqual(fallback_stats2.get('atommap_used', 0), 0)
            self.assertEqual(fallback_stats2.get('heuristic_used', 0), 1)
            mock_heuristic.assert_called_once()
    
    def test_accumulative_semantics_at_product_level(self):
        """Test that fallback flags accumulate correctly at the product detection level."""
        
        # Test AtomMap accumulation
        product_mol_with_map = Chem.MolFromSmiles('Fc1ccccc1')
        if product_mol_with_map:
            for atom in product_mol_with_map.GetAtoms():
                if atom.GetSymbol() == 'C':
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'F':
                            atom.SetAtomMapNum(1)
                            break
                    if atom.GetAtomMapNum() == 1:
                        break
            Chem.SanitizeMol(product_mol_with_map)
        
        atommap_stats = {}
        
        # Call the function multiple times (simulating multiple products in an attempt)
        # Note: We need fresh molecules each time because the function clears atom maps
        for i in range(3):
            # Create a fresh molecule with atom map each time
            fresh_product = Chem.MolFromSmiles('Fc1ccccc1')
            if fresh_product:
                for atom in fresh_product.GetAtoms():
                    if atom.GetSymbol() == 'C':
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetSymbol() == 'F':
                                atom.SetAtomMapNum(1)
                                break
                        if atom.GetAtomMapNum() == 1:
                            break
                Chem.SanitizeMol(fresh_product)
            
            _detect_site_from_product_with_mapping_fallback(
                self.parent_mol, fresh_product, 'R1', 'F', atommap_stats
            )
        
        # AtomMap should accumulate
        self.assertGreaterEqual(atommap_stats.get('atommap_used', 0), 1)
        self.assertEqual(atommap_stats.get('heuristic_used', 0), 0)
        
        # Test heuristic accumulation separately
        product_mol_no_map = Chem.MolFromSmiles('Fc1ccccc1')
        if product_mol_no_map:
            for atom in product_mol_no_map.GetAtoms():
                atom.SetAtomMapNum(0)
            Chem.SanitizeMol(product_mol_no_map)
        
        heuristic_stats = {}
        
        with patch('src.halogenator.enumerate_k._detect_site_from_product') as mock_heuristic:
            mock_heuristic.return_value = (1, 'benzene')
            
            # Call multiple times
            for i in range(2):
                _detect_site_from_product_with_mapping_fallback(
                    self.parent_mol, product_mol_no_map, 'R1', 'F', heuristic_stats
                )
            
            # Heuristic should accumulate
            self.assertGreaterEqual(heuristic_stats.get('heuristic_used', 0), 1)
            self.assertEqual(heuristic_stats.get('atommap_used', 0), 0)
    
    def test_end_to_end_atommap_integration(self):
        """Test AtomMap success in the context of enumerate_with_stats."""
        from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig
        
        # This is a minimal end-to-end test using the public API
        # We'll use a real parent molecule and verify the fallback paths work correctly
        
        cfg = EnumConfig(k_max=1, halogens=('F',))
        
        # Use a simple aromatic parent that should work with the real enumeration
        parent_smiles = 'c1ccccc1'  # benzene
        
        # Call the real enumeration function
        products, qa_stats = enumerate_with_stats(parent_smiles, cfg)
        
        # The real enumeration should produce some results
        # We can't guarantee specific fallback paths without mocking, but we can verify
        # that the QA stats structure is correct and consistent
        
        self.assertIn('qa_paths', qa_stats)
        self.assertIn('atommap_used', qa_stats['qa_paths'])
        self.assertIn('heuristic_used', qa_stats['qa_paths'])
        
        # Verify that atommap_used and heuristic_used are integers (not accumulated accidentally)
        self.assertIsInstance(qa_stats['qa_paths']['atommap_used'], int)
        self.assertIsInstance(qa_stats['qa_paths']['heuristic_used'], int)
        
        # Verify they're non-negative
        self.assertGreaterEqual(qa_stats['qa_paths']['atommap_used'], 0)
        self.assertGreaterEqual(qa_stats['qa_paths']['heuristic_used'], 0)
        
        # If we have pivots, verify consistency
        if 'pivots' in qa_stats and qa_stats['pivots'].get('by_rule_halogen_k'):
            pivots = qa_stats['pivots']['by_rule_halogen_k']
            
            # Sum up fallback usage across all attempts
            total_atommap = sum(events.get('atommap_used', 0) for events in pivots.values())
            total_heuristic = sum(events.get('heuristic_used', 0) for events in pivots.values())
            
            # Should match totals
            self.assertEqual(qa_stats['qa_paths']['atommap_used'], total_atommap)
            self.assertEqual(qa_stats['qa_paths']['heuristic_used'], total_heuristic)


if __name__ == '__main__':
    unittest.main()