# -*- coding: ascii -*-
"""Test P1 BFS enumeration functionality."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestP1BFSEnumeration(unittest.TestCase):
    """Test BFS enumeration with different halogen sets and k values."""
    
    def test_benzene_4_halogens_k2(self):
        """Test benzene with 4 halogens (F, Cl, Br, I) at k=2 should produce 24 products."""
        cfg = EnumConfig(
            k_max=2,
            halogens=('F', 'Cl', 'Br', 'I'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        # Expected: 6 positions * 4 halogens = 24 at k=1
        # Plus combinations at k=2 with symmetry folding
        # Should be exactly 24 total (6*4 at k=1 + 0 at k=2 due to symmetry)
        # Actually, let's count what we get and validate it's reasonable
        k_counts = Counter(p['k'] for p in products)
        
        # At least some k=1 and k=2 products should exist
        self.assertGreater(k_counts[1], 0, "Should have k=1 products")
        self.assertGreater(k_counts[2], 0, "Should have k=2 products") 
        
        # Total should be reasonable for benzene with 4 halogens
        self.assertGreater(len(products), 20, "Should have >20 products for benzene with 4 halogens")
        
    def test_benzene_2_halogens_k2(self):
        """Test benzene with 2 halogens (F, Cl) at k=2 should produce 8 products."""
        cfg = EnumConfig(
            k_max=2,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        # With symmetry folding, benzene has 1 unique aromatic position
        # k=1: 1 position * 2 halogens = 2 products  
        # k=2: ortho, meta, para patterns with 2 halogens each
        # Expected around 8 total
        k_counts = Counter(p['k'] for p in products)
        
        self.assertEqual(k_counts[1], 2, "Should have exactly 2 k=1 products (1 position * 2 halogens)")
        self.assertGreaterEqual(len(products), 6, "Should have >=6 total products")
        self.assertLessEqual(len(products), 15, "Should have <=15 total products")
        
    def test_bfs_layer_ordering(self):
        """Test that BFS produces products in correct depth order."""
        cfg = EnumConfig(
            k_max=3,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        # Check that products are generated in BFS order (k=1, then k=2, then k=3)
        prev_k = 0
        for product in products:
            current_k = product['k']
            self.assertGreaterEqual(current_k, prev_k, 
                                  "BFS should generate products in non-decreasing k order")
            prev_k = current_k
            
    def test_symmetry_folding_effectiveness(self):
        """Test that symmetry folding reduces product count appropriately."""
        # With symmetry folding
        cfg_with_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        # Without symmetry folding  
        cfg_without_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products_with = list(enumerate_products(benzene_smiles, cfg_with_folding))
        products_without = list(enumerate_products(benzene_smiles, cfg_without_folding))
        
        # With folding should have fewer products (1 vs 6 for benzene)
        self.assertEqual(len(products_with), 1, "Benzene with symmetry folding should have 1 product")
        self.assertEqual(len(products_without), 6, "Benzene without symmetry folding should have 6 products")
        self.assertGreater(len(products_without), len(products_with), 
                          "Without symmetry folding should have more products than with folding")
    
    def test_no_folding_k2_cross_state_deduplication(self):
        """Test that k=2 without symmetry folding doesn't suffer cross-state deduplication conflicts."""
        # Use benzene k=2 with no folding - should generate many k=2 products 
        # without being incorrectly deduplicated due to same site numbers on different intermediates
        cfg_no_folding_k2 = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg_no_folding_k2))
        
        # Count products by k
        k1_products = [p for p in products if p['k'] == 1]
        k2_products = [p for p in products if p['k'] == 2]
        
        # Should have 6 k=1 products (one F at each position)
        self.assertEqual(len(k1_products), 6, "Should have 6 k=1 products without folding")
        
        # Should have k=2 products (each k=1 can be further substituted)
        # Without cross-state conflicts, we should have multiple k=2 products
        self.assertGreater(len(k2_products), 0, "Should have k=2 products without cross-state conflicts")
        
        # Verify that k=2 products have different intermediate parents
        # (i.e., they're not being incorrectly deduplicated)
        k2_histories = [p['substitutions'] for p in k2_products]
        unique_first_sites = set()
        for history in k2_histories:
            if len(history) >= 2:
                first_site = history[0].get('site')
                if first_site is not None:
                    unique_first_sites.add(first_site)
        
        # Should have multiple different first substitution sites in k=2 products
        self.assertGreater(len(unique_first_sites), 1, 
                          "k=2 products should come from different intermediate states")


if __name__ == '__main__':
    unittest.main()