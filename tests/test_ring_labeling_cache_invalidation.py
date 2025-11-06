# -*- coding: ascii -*-
"""Tests for ring labeling cache invalidation on topology changes."""

import unittest
from unittest.mock import patch, MagicMock
from rdkit import Chem
from src.halogenator.sites import flavonoid_ring_label, clear_flavonoid_ring_label_cache, get_flavonoid_ring_label_cache_info, _compute_molecule_signature


class TestRingLabelingCacheInvalidation(unittest.TestCase):
    """Test cache invalidation when molecular topology changes."""
    
    def setUp(self):
        """Clear cache before each test."""
        clear_flavonoid_ring_label_cache()
    
    def test_cache_hit_same_topology(self):
        """Test that cache hits work when topology stays the same."""
        # Use quercetin as test molecule
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # First call - cache miss
        info_before = get_flavonoid_ring_label_cache_info()
        label1 = flavonoid_ring_label(mol, 0)
        info_after_first = get_flavonoid_ring_label_cache_info()
        
        # Verify cache miss occurred
        self.assertEqual(info_after_first.misses, info_before.misses + 1)
        self.assertEqual(info_after_first.hits, info_before.hits)
        self.assertEqual(info_after_first.invalidations, info_before.invalidations)
        
        # Second call with same molecule - should be cache hit
        label2 = flavonoid_ring_label(mol, 0)
        info_after_second = get_flavonoid_ring_label_cache_info()
        
        # Verify cache hit occurred
        self.assertEqual(info_after_second.hits, info_after_first.hits + 1)
        self.assertEqual(info_after_second.misses, info_after_first.misses)
        self.assertEqual(info_after_second.invalidations, info_after_first.invalidations)
        
        # Labels should be identical
        self.assertEqual(label1, label2)
        self.assertIsInstance(label1, str)
    
    def test_topology_change_triggers_invalidation(self):
        """Test topology change via signature mocking triggers invalidation."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # First call - populate cache
        clear_flavonoid_ring_label_cache()
        label_before = flavonoid_ring_label(mol, 0)
        info_after_first = get_flavonoid_ring_label_cache_info()
        
        self.assertEqual(info_after_first.misses, 1)
        self.assertEqual(info_after_first.invalidations, 0)
        
        # Mock _compute_molecule_signature to return different signature on next call
        original_signature = _compute_molecule_signature(mol)
        different_signature = (original_signature[0], original_signature[1] + 1)
        
        with patch('src.halogenator.sites._compute_molecule_signature', side_effect=[different_signature, different_signature]):
            # This should trigger invalidation due to signature change
            label_after = flavonoid_ring_label(mol, 0)
            info_after_change = get_flavonoid_ring_label_cache_info()
            
            # Verify invalidation and recomputation occurred
            self.assertEqual(info_after_change.invalidations, info_after_first.invalidations + 1)
            self.assertEqual(info_after_change.misses, info_after_first.misses + 1)
    
    def test_isotope_tagging_does_not_trigger_invalidation(self):
        """Test that isotope changes do not trigger cache invalidation (by design)."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # First call - populate cache
        clear_flavonoid_ring_label_cache()
        label_before = flavonoid_ring_label(mol, 0)
        info_after_first = get_flavonoid_ring_label_cache_info()
        
        # Set isotope on atom (this should NOT change topology signature)
        original_isotope = mol.GetAtomWithIdx(0).GetIsotope()
        mol.GetAtomWithIdx(0).SetIsotope(2)
        
        try:
            # Verify signature is unchanged after isotope tagging
            sig1 = _compute_molecule_signature(mol)
            mol.GetAtomWithIdx(0).SetIsotope(0)  # Reset
            sig2 = _compute_molecule_signature(mol)
            mol.GetAtomWithIdx(0).SetIsotope(2)  # Set again
            sig3 = _compute_molecule_signature(mol)
            
            # Isotope changes should not affect topology signature
            self.assertEqual(sig1, sig2)
            self.assertEqual(sig2, sig3)
            
            # Second call with isotope - should be cache hit, no invalidation
            label_after = flavonoid_ring_label(mol, 0)
            info_after_isotope = get_flavonoid_ring_label_cache_info()
            
            # Should be cache hit, not invalidation
            self.assertEqual(info_after_isotope.hits, info_after_first.hits + 1)
            self.assertEqual(info_after_isotope.invalidations, info_after_first.invalidations)
            self.assertEqual(label_before, label_after)
            
        finally:
            # Clean up isotope
            mol.GetAtomWithIdx(0).SetIsotope(original_isotope)
    
    def test_different_molecules_separate_cache_entries(self):
        """Test that different molecules get separate cache entries."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        kaempferol_smi = "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"
        
        mol1 = Chem.MolFromSmiles(quercetin_smi)
        mol2 = Chem.MolFromSmiles(kaempferol_smi)
        
        self.assertIsNotNone(mol1)
        self.assertIsNotNone(mol2)
        
        clear_flavonoid_ring_label_cache()
        
        # Process molecule 1
        label1 = flavonoid_ring_label(mol1, 0)
        info_after_mol1 = get_flavonoid_ring_label_cache_info()
        
        # Process molecule 2 - should not interfere with mol1 cache
        label2 = flavonoid_ring_label(mol2, 0)
        info_after_mol2 = get_flavonoid_ring_label_cache_info()
        
        # Both should have valid labels
        self.assertIsInstance(label1, str)
        self.assertIsInstance(label2, str)
        
        # Should have 2 cache misses (one for each molecule)
        self.assertEqual(info_after_mol2.misses, 2)
        self.assertEqual(info_after_mol2.invalidations, 0)
        
        # Accessing mol1 again should be cache hit
        label1_again = flavonoid_ring_label(mol1, 0)
        info_final = get_flavonoid_ring_label_cache_info()
        
        self.assertEqual(info_final.hits, 1)
        self.assertEqual(label1, label1_again)
    
    def test_keyerror_safety_in_cache_deletion(self):
        """Test that KeyError in cache deletion is handled gracefully."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        # Populate cache first
        clear_flavonoid_ring_label_cache()
        label1 = flavonoid_ring_label(mol, 0)
        info_before = get_flavonoid_ring_label_cache_info()
        
        # Import weakref to patch the correct method
        import weakref
        
        # Mock the deletion to raise KeyError - patch the type method, not instance
        with patch.object(weakref.WeakKeyDictionary, '__delitem__', side_effect=KeyError("test")):
            # Mock signature to trigger invalidation path
            with patch('src.halogenator.sites._compute_molecule_signature', return_value=(999, 999)):
                # This should not crash despite KeyError in deletion
                label2 = flavonoid_ring_label(mol, 0)
                info_after = get_flavonoid_ring_label_cache_info()
                
                # Verify the function completed without error and invalidation count increased
                self.assertIsInstance(label2, str)
                self.assertEqual(info_after.invalidations, info_before.invalidations + 1)
    
    def test_topology_or_aromatic_change_triggers_invalidation(self):
        """Test that modifying bond order or toggling aromaticity triggers invalidation."""
        # Use a simple aromatic molecule 
        benzene_smi = "c1ccccc1"
        mol = Chem.MolFromSmiles(benzene_smi)
        self.assertIsNotNone(mol)
        
        # Populate cache first
        clear_flavonoid_ring_label_cache()
        label_before = flavonoid_ring_label(mol, 0)
        info_after_first = get_flavonoid_ring_label_cache_info()
        self.assertEqual(info_after_first.misses, 1)
        self.assertEqual(info_after_first.invalidations, 0)
        
        # Modify aromaticity - create non-aromatic version
        mol_copy = Chem.RWMol(mol)
        
        # Get original signature
        original_signature = _compute_molecule_signature(mol)
        
        # Kekulize to remove aromaticity (benzene -> cyclohexatriene)
        try:
            Chem.Kekulize(mol_copy, clearAromaticFlags=True)
            modified_signature = _compute_molecule_signature(mol_copy)
            
            # Signatures should be different due to aromaticity change
            self.assertNotEqual(original_signature, modified_signature, 
                              "Signature should change when aromaticity changes")
            
            # Mock the signature to simulate this topology change during lookup
            with patch('src.halogenator.sites._compute_molecule_signature', return_value=modified_signature):
                # This should trigger invalidation due to signature change
                label_after = flavonoid_ring_label(mol, 0)
                info_after_change = get_flavonoid_ring_label_cache_info()
                
                # Verify invalidation occurred
                self.assertEqual(info_after_change.invalidations, info_after_first.invalidations + 1)
                self.assertEqual(info_after_change.misses, info_after_first.misses + 1)
                
        except Exception:
            # If kekulization fails, just verify that different SMILES would trigger invalidation
            different_signature = (original_signature[0], "C1=CC=CC=C1")  # Different SMILES
            
            with patch('src.halogenator.sites._compute_molecule_signature', return_value=different_signature):
                label_after = flavonoid_ring_label(mol, 0)
                info_after_change = get_flavonoid_ring_label_cache_info()
                
                # Verify invalidation occurred
                self.assertEqual(info_after_change.invalidations, info_after_first.invalidations + 1)
                self.assertEqual(info_after_change.misses, info_after_first.misses + 1)
    
    def test_signature_fallback_safety(self):
        """Test that signature fallback handles MolToSmiles exceptions gracefully."""
        quercetin_smi = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smi)
        self.assertIsNotNone(mol)
        
        clear_flavonoid_ring_label_cache()
        
        # Mock MolToSmiles to raise exception to trigger fallback path
        with patch('src.halogenator.sites.Chem.MolToSmiles', side_effect=Exception("mock failure")):
            # This should not crash and should use fallback signature
            label = flavonoid_ring_label(mol, 0)
            info_after = get_flavonoid_ring_label_cache_info()
            
            # Should get some valid label (might be empty string for non-flavonoids in fallback path)
            self.assertIsInstance(label, str)
            
            # Cache should have recorded a miss (successful fallback computation)
            self.assertEqual(info_after.misses, 1)
            self.assertEqual(info_after.invalidations, 0)
            
            # Second call should hit cache even with fallback signature
            label2 = flavonoid_ring_label(mol, 0)
            info_after_second = get_flavonoid_ring_label_cache_info()
            
            # Should be cache hit this time
            self.assertEqual(info_after_second.hits, 1)
            self.assertEqual(info_after_second.misses, 1)
            self.assertEqual(label, label2)


if __name__ == '__main__':
    unittest.main()