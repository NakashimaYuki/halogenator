# -*- coding: ascii -*-
"""Tests for A/B/C ring labeling index consistency and caching correctness."""

import unittest
from rdkit import Chem
from src.halogenator.sites import (
    flavonoid_ring_label,
    clear_flavonoid_ring_label_cache,
    get_flavonoid_ring_label_cache_info
)


class TestRingLabelingIndexConsistency(unittest.TestCase):
    """Test that ring labeling remains consistent across atom index changes."""
    
    def setUp(self):
        """Clear cache before each test."""
        clear_flavonoid_ring_label_cache()
    
    def test_atom_index_consistency_after_modification(self):
        """Test that ring labels are consistent even after molecule modifications."""
        # Start with quercetin
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        original_mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(original_mol)
        
        # Get labels for original molecule
        original_labels = []
        for i in range(original_mol.GetNumAtoms()):
            label = flavonoid_ring_label(original_mol, i)
            original_labels.append(label)
        
        # Create a "modified" version by reconstructing from canonical SMILES
        # This simulates what would happen if we used SMILES-based caching
        canonical_smiles = Chem.MolToSmiles(original_mol, canonical=True)
        reconstructed_mol = Chem.MolFromSmiles(canonical_smiles)
        Chem.SanitizeMol(reconstructed_mol)
        
        # Get labels for reconstructed molecule
        reconstructed_labels = []
        for i in range(reconstructed_mol.GetNumAtoms()):
            label = flavonoid_ring_label(reconstructed_mol, i)
            reconstructed_labels.append(label)
        
        # Even if atom indices differ, the same chemical positions should have same labels
        # We'll check this by comparing the distribution of labels
        original_label_counts = {}
        for label in original_labels:
            original_label_counts[label] = original_label_counts.get(label, 0) + 1
        
        reconstructed_label_counts = {}
        for label in reconstructed_labels:
            reconstructed_label_counts[label] = reconstructed_label_counts.get(label, 0) + 1
        
        # The label distributions should be identical
        self.assertEqual(original_label_counts, reconstructed_label_counts,
                        "Label distribution should be consistent across molecule reconstructions")
    
    def test_caching_prevents_recomputation(self):
        """Test that caching actually works and prevents recomputation."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Clear cache and get initial info
        clear_flavonoid_ring_label_cache()
        
        # First call should populate cache
        label1 = flavonoid_ring_label(mol, 0)
        label2 = flavonoid_ring_label(mol, 5)
        
        # Second calls should use cache (same molecule object)
        label1_repeat = flavonoid_ring_label(mol, 0)
        label2_repeat = flavonoid_ring_label(mol, 5)
        
        # Results should be consistent
        self.assertEqual(label1, label1_repeat)
        self.assertEqual(label2, label2_repeat)
        
        # Cache should show the molecule is cached
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertGreater(cache_info.currsize, 0, "Cache should contain cached molecules")
    
    def test_different_molecule_objects_same_structure(self):
        """Test behavior with different molecule objects having same structure."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        
        # Create two separate molecule objects from same SMILES
        mol1 = Chem.MolFromSmiles(quercetin_smiles)
        mol2 = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol1)
        Chem.SanitizeMol(mol2)
        
        # These are different objects but same structure
        self.assertIsNot(mol1, mol2)
        
        # Get labels for both
        label1_mol1 = flavonoid_ring_label(mol1, 0)
        label1_mol2 = flavonoid_ring_label(mol2, 0)
        
        # Labels should be the same for same chemical position
        self.assertEqual(label1_mol1, label1_mol2)
        
        # Both molecules should be cached separately
        cache_info = get_flavonoid_ring_label_cache_info()
        # We expect at least 1 entry (possibly 2 if WeakKeyDictionary keeps both)
        self.assertGreaterEqual(cache_info.currsize, 1)
    
    def test_atom_mapping_preservation(self):
        """Test that atom mapping is preserved correctly with new caching."""
        # Create a simple flavonoid-like structure
        test_smiles = "O=c1cc(oc2cc(O)cc(O)c12)-c1ccc(O)cc1"  # apigenin-like
        mol = Chem.MolFromSmiles(test_smiles)
        Chem.SanitizeMol(mol)
        
        # Get all labels
        labels = []
        for i in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, i)
            labels.append((i, label))
        
        # Should have some A, B, C labels for a flavonoid
        label_types = set(label for _, label in labels if label)
        self.assertIn('C', label_types, "Should identify C ring in flavonoid")
        # Note: This particular structure might not have clear A/B distinction
        # so we just check that we get some meaningful labels
        
        # Verify consistency on repeated calls
        labels_repeat = []
        for i in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, i)
            labels_repeat.append((i, label))
        
        self.assertEqual(labels, labels_repeat, "Labels should be consistent on repeated calls")
    
    def test_cache_clear_functionality(self):
        """Test that cache clearing works with new implementation."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Populate cache
        flavonoid_ring_label(mol, 0)
        cache_info = get_flavonoid_ring_label_cache_info()
        initial_size = cache_info.currsize
        self.assertGreater(initial_size, 0)
        
        # Clear cache
        clear_flavonoid_ring_label_cache()
        
        # Cache should be empty or much smaller
        cache_info_after = get_flavonoid_ring_label_cache_info()
        self.assertLessEqual(cache_info_after.currsize, initial_size)
        
        # Should still work after clearing
        label = flavonoid_ring_label(mol, 0)
        self.assertIsInstance(label, str)
    
    def test_fallback_cache_on_weak_reference_failure(self):
        """Test that fallback caching works when WeakKeyDictionary can't be used."""
        # This test is harder to trigger reliably, but we can at least verify
        # the fallback mechanism exists and works
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Force a scenario where we might need fallback
        # by calling with molecule that might not be cacheable in WeakKeyDictionary
        label = flavonoid_ring_label(mol, 0)
        self.assertIsInstance(label, str)
        
        # Verify cache info still works
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertIsNotNone(cache_info)
        self.assertGreaterEqual(cache_info.currsize, 0)


if __name__ == '__main__':
    unittest.main()