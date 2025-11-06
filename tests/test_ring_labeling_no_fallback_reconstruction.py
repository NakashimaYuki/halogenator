# -*- coding: ascii -*-
"""Tests to ensure no molecule reconstruction fallback in ring labeling."""

import unittest
from unittest.mock import patch, MagicMock
from rdkit import Chem
from src.halogenator.sites import flavonoid_ring_label, clear_flavonoid_ring_label_cache, _compute_molecule_signature


class TestRingLabelingNoFallbackReconstruction(unittest.TestCase):
    """Test that ring labeling never reconstructs molecules from InChI/SMILES."""
    
    def setUp(self):
        """Clear cache before each test."""
        clear_flavonoid_ring_label_cache()
    
    def test_no_reconstruction_functions_called(self):
        """Test that no molecule reconstruction happens during ring labeling."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction function called - this violates no-reconstruction guarantee!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            # Get ring labels multiple times - should not call reconstruction
            for atom_idx in range(min(5, mol.GetNumAtoms())):
                label = flavonoid_ring_label(mol, atom_idx)
                self.assertIsInstance(label, str)
                self.assertIn(label, ['', 'A', 'B', 'C'])
            
            # Cache hit path should also not call reconstruction
            for atom_idx in range(min(5, mol.GetNumAtoms())):
                label = flavonoid_ring_label(mol, atom_idx)
                self.assertIsInstance(label, str)
    
    def test_weakref_failure_uses_live_molecule(self):
        """Test that WeakKeyDictionary failure leads to direct computation on live molecule."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction called when it should compute on live molecule!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            # Force weak reference failure by mocking dictionary operations
            with patch('src.halogenator.sites._mol_ring_labels_cache.__contains__', side_effect=TypeError("Not weak-referenceable")), \
                 patch('src.halogenator.sites._mol_ring_labels_cache.__setitem__', side_effect=TypeError("Not weak-referenceable")):
                
                # This should still work by computing directly on live molecule
                label = flavonoid_ring_label(mol, 0)
                self.assertIsInstance(label, str)
                self.assertIn(label, ['', 'A', 'B', 'C'])
    
    def test_signature_computation_uses_live_molecule(self):
        """Test that signature computation operates on live molecule without reconstruction."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction called during signature computation!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            # Compute signature - should use original molecule only
            signature = _compute_molecule_signature(mol)
            self.assertIsInstance(signature, tuple)
            self.assertEqual(len(signature), 2)
            self.assertIsInstance(signature[0], int)  # num_atoms
            self.assertGreaterEqual(signature[0], 1)  # Should have atoms
    
    def test_different_molecules_no_cross_contamination(self):
        """Test that different molecules are processed independently without reconstruction."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        kaempferol_smi = "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"
        
        mol1 = Chem.MolFromSmiles(quercetin_smi)
        mol2 = Chem.MolFromSmiles(kaempferol_smi)
        
        self.assertIsNotNone(mol1)
        self.assertIsNotNone(mol2)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction called - should use live molecules!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            clear_flavonoid_ring_label_cache()
            
            # Process both molecules - should work without reconstruction
            labels1 = []
            for i in range(min(5, mol1.GetNumAtoms())):
                label = flavonoid_ring_label(mol1, i)
                labels1.append(label)
                self.assertIsInstance(label, str)
            
            labels2 = []
            for i in range(min(5, mol2.GetNumAtoms())):
                label = flavonoid_ring_label(mol2, i)
                labels2.append(label)
                self.assertIsInstance(label, str)
    
    def test_error_paths_no_reconstruction(self):
        """Test that error handling paths don't trigger molecule reconstruction."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # Simple benzene
        self.assertIsNotNone(mol)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction called during error handling!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            # Mock implementation to cause error, should fallback without reconstruction
            with patch('src.halogenator.sites._flavonoid_ring_label_impl', side_effect=Exception("Test error")):
                label = flavonoid_ring_label(mol, 0)
                self.assertEqual(label, '')  # Should return empty on error
    
    def test_cache_operations_no_reconstruction(self):
        """Test that cache operations (hit/miss/invalidation) never trigger reconstruction."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # Mock reconstruction functions to raise if called
        def raise_if_called(*args, **kwargs):
            raise AssertionError("Molecule reconstruction called during cache operations!")
        
        with patch('rdkit.Chem.MolFromInchi', side_effect=raise_if_called), \
             patch('rdkit.Chem.MolFromSmiles', side_effect=raise_if_called):
            
            clear_flavonoid_ring_label_cache()
            
            # Cache miss
            label1 = flavonoid_ring_label(mol, 0)
            self.assertIsInstance(label1, str)
            
            # Cache hit
            label2 = flavonoid_ring_label(mol, 0)
            self.assertEqual(label1, label2)
            
            # Force invalidation by mocking signature change
            with patch('src.halogenator.sites._compute_molecule_signature', side_effect=[(999, 999), (999, 999)]):
                label3 = flavonoid_ring_label(mol, 0)
                self.assertIsInstance(label3, str)


if __name__ == '__main__':
    unittest.main()