# -*- coding: ascii -*-
"""Tests for A/B/C ring labeling caching performance optimization."""

import unittest
from rdkit import Chem
from src.halogenator.sites import (
    flavonoid_ring_label, 
    get_flavonoid_ring_label_cache_info,
    clear_flavonoid_ring_label_cache
)


class TestFlavonoidRingLabelingCache(unittest.TestCase):
    """Test LRU caching for flavonoid ring labeling performance."""
    
    def setUp(self):
        """Clear cache before each test."""
        clear_flavonoid_ring_label_cache()
    
    def test_caching_reduces_computation(self):
        """Test that repeated calls use cache instead of recomputing."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Initial cache should be empty
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.hits, 0)
        self.assertEqual(cache_info.misses, 0)
        
        # First call should be a cache miss (computes all atoms for the molecule)
        label1 = flavonoid_ring_label(mol, 0)
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 1)
        self.assertEqual(cache_info.hits, 0)
        
        # Second call to different atom of same molecule should be cache hit
        label2 = flavonoid_ring_label(mol, 5)
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 1)  # Still 1 miss
        self.assertEqual(cache_info.hits, 1)   # Now 1 hit
        
        # Repeated calls should all be cache hits
        label1_repeat = flavonoid_ring_label(mol, 0)
        label2_repeat = flavonoid_ring_label(mol, 5)
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 1)  # Still 1 miss
        self.assertEqual(cache_info.hits, 3)    # Now 3 hits total
        
        # Results should be consistent
        self.assertEqual(label1, label1_repeat)
        self.assertEqual(label2, label2_repeat)
    
    def test_different_molecules_separate_cache_entries(self):
        """Test that different molecules have separate cache entries."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        apigenin_smiles = "O=c1cc(oc2cc(O)cc(O)c12)-c1ccc(O)cc1"
        
        mol1 = Chem.MolFromSmiles(quercetin_smiles)
        mol2 = Chem.MolFromSmiles(apigenin_smiles)
        
        Chem.SanitizeMol(mol1)
        Chem.SanitizeMol(mol2)
        
        # Each molecule should create separate cache entries
        flavonoid_ring_label(mol1, 0)  # Miss for quercetin (computes all atoms)
        flavonoid_ring_label(mol2, 0)  # Miss for apigenin (computes all atoms)
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 2)
        
        # Repeated calls should hit cache
        flavonoid_ring_label(mol1, 0)  # Hit for quercetin
        flavonoid_ring_label(mol2, 0)  # Hit for apigenin
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.hits, 2)
        self.assertEqual(cache_info.misses, 2)
    
    def test_cache_handles_non_flavonoids(self):
        """Test that cache works correctly with non-flavonoid molecules."""
        benzene_smiles = "c1ccccc1"
        mol = Chem.MolFromSmiles(benzene_smiles)
        Chem.SanitizeMol(mol)
        
        # Should return empty labels and cache the result
        label1 = flavonoid_ring_label(mol, 0)
        label2 = flavonoid_ring_label(mol, 0)  # Should hit cache
        
        self.assertEqual(label1, '')
        self.assertEqual(label2, '')
        
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.hits, 1)
        self.assertEqual(cache_info.misses, 1)
    
    def test_cache_consistency_with_multiple_atoms(self):
        """Test cache consistency when calling for multiple atoms in same molecule."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Call for all atoms in the molecule
        labels_first = []
        for i in range(mol.GetNumAtoms()):
            labels_first.append(flavonoid_ring_label(mol, i))
        
        # Call again for all atoms - should hit cache
        labels_second = []
        for i in range(mol.GetNumAtoms()):
            labels_second.append(flavonoid_ring_label(mol, i))
        
        # Results should be identical
        self.assertEqual(labels_first, labels_second)
        
        # Should have 1 miss (first call computes all labels) and (2*num_atoms - 1) hits total
        # First pass: 1 miss + (num_atoms-1) hits, Second pass: num_atoms hits
        cache_info = get_flavonoid_ring_label_cache_info()
        expected_hits = (mol.GetNumAtoms() - 1) + mol.GetNumAtoms()  # First pass hits + second pass hits
        self.assertEqual(cache_info.hits, expected_hits)
        self.assertEqual(cache_info.misses, 1)  # Only first call is a miss
    
    def test_cache_clear_functionality(self):
        """Test that cache clearing works correctly."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        mol = Chem.MolFromSmiles(quercetin_smiles)
        Chem.SanitizeMol(mol)
        
        # Populate cache
        flavonoid_ring_label(mol, 0)
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 1)
        
        # Clear cache
        clear_flavonoid_ring_label_cache()
        
        # Cache should be empty
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.hits, 0)
        self.assertEqual(cache_info.misses, 0)
        
        # Next call should be a miss again
        flavonoid_ring_label(mol, 0)
        cache_info = get_flavonoid_ring_label_cache_info()
        self.assertEqual(cache_info.misses, 1)
        self.assertEqual(cache_info.hits, 0)


if __name__ == '__main__':
    unittest.main()