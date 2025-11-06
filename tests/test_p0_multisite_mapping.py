# -*- coding: ascii -*-
"""Test P0-A/P0-B: Multi-site R3/R4/R5 site-to-product mapping and deduplication."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestMultiSiteMapping(unittest.TestCase):
    """Test correct site mapping for multi-site substrates."""
    
    def test_r3_dihydroxybenzene_site_mapping(self):
        """Test R3 on 1,4-dihydroxybenzene with correct site mapping."""
        # 1,4-dihydroxybenzene (hydroquinone) - has two symmetric OH groups
        dihydroxy_smiles = "c1cc(O)ccc1O"
        
        # Test with symmetry folding enabled (should get 1 product)
        cfg_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        products_folding = list(enumerate_products(dihydroxy_smiles, cfg_folding))
        r3_products_folding = [p for p in products_folding if p['rule'] == 'R3']
        
        # With folding, should have only 1 R3 product (OH groups are symmetric)
        self.assertEqual(len(r3_products_folding), 1, 
                        "With folding, should have 1 R3 product for symmetric diOH")
        
        # Test with symmetry folding disabled (should get 2 products)
        cfg_no_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products_no_folding = list(enumerate_products(dihydroxy_smiles, cfg_no_folding))
        r3_products_no_folding = [p for p in products_no_folding if p['rule'] == 'R3']
        
        # Without folding, should have 2 R3 products (one for each OH site)
        self.assertEqual(len(r3_products_no_folding), 2, 
                        "Without folding, should have 2 R3 products for two OH sites")
        
        # Verify that the two products have different site indices
        sites = [p['substitutions'][0]['site'] for p in r3_products_no_folding]
        unique_sites = set(sites)
        self.assertEqual(len(unique_sites), 2, 
                        "The two R3 products should have different site indices")
        
    def test_r4_diaminobenzene_site_mapping(self):
        """Test R4 on 1,3-diaminobenzene with correct site mapping."""
        # 1,3-diaminobenzene - has two NH2 groups
        diamino_smiles = "c1cc(N)cc(N)c1"
        
        cfg_no_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(diamino_smiles, cfg_no_folding))
        r4_products = [p for p in products if p['rule'] == 'R4']
        
        # Should have 2 R4 products (one for each NH2 site)
        self.assertEqual(len(r4_products), 2, 
                        "Should have 2 R4 products for two NH2 sites")
        
        # Verify different site indices
        sites = [p['substitutions'][0]['site'] for p in r4_products]
        unique_sites = set(sites)
        self.assertEqual(len(unique_sites), 2, 
                        "The two R4 products should have different site indices")
        
    def test_deduplication_strategy_alignment(self):
        """Test that R3 deduplication strategy aligns with R1 (site rules)."""
        # Use phenol which has both aromatic CH (R1) and OH (R3) sites
        phenol_smiles = "c1ccccc1O"
        
        cfg_no_folding = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(phenol_smiles, cfg_no_folding))
        
        # Check that we get both k=1 and k=2 products
        k_counts = Counter(p['k'] for p in products)
        self.assertGreater(k_counts[1], 0, "Should have k=1 products")
        self.assertGreater(k_counts[2], 0, "Should have k=2 products")
        
        # For k=2, should have products from different k=1 intermediates
        k2_products = [p for p in products if p['k'] == 2]
        
        # Check that k=2 products have different first substitution sites
        first_sites = []
        for product in k2_products:
            if len(product['substitutions']) >= 2:
                first_site = product['substitutions'][0].get('site')
                if first_site is not None:
                    first_sites.append(first_site)
        
        unique_first_sites = set(first_sites)
        self.assertGreater(len(unique_first_sites), 1, 
                          "k=2 products should come from different k=1 intermediates")
    
    def test_r5_dicarboxylic_acid_site_mapping(self):
        """Test R5 on terephthalic acid (1,4-dicarboxylic acid) with correct site mapping."""
        # Terephthalic acid: 1,4-benzenedicarboxylic acid - has two COOH groups
        dicarboxylic_smiles = "O=C(O)c1ccc(C(=O)O)cc1"
        
        cfg_no_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(dicarboxylic_smiles, cfg_no_folding))
        r5_products = [p for p in products if p['rule'] == 'R5']
        
        # Should have 2 R5 products (one for each COOH group)
        self.assertEqual(len(r5_products), 2, 
                        "Should have 2 R5 products for two COOH sites")
        
        # Verify different site indices
        sites = [p['substitutions'][0]['site'] for p in r5_products]
        unique_sites = set(sites)
        self.assertEqual(len(unique_sites), 2, 
                        "The two R5 products should have different site indices")
        
        # Test with folding - should get only 1 product due to symmetry
        cfg_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        products_folding = list(enumerate_products(dicarboxylic_smiles, cfg_folding))
        r5_products_folding = [p for p in products_folding if p['rule'] == 'R5']
        
        self.assertEqual(len(r5_products_folding), 1, 
                        "With folding, should have 1 R5 product for symmetric diCOOH")

    def test_mixed_multisite_substrate(self):
        """Test substrate with multiple R3, R4, R5 sites simultaneously."""
        # 4-aminosalicylic acid: has both NH2 (R4) and COOH (R5) and OH (R3)
        mixed_smiles = "Nc1ccc(C(=O)O)c(O)c1"  # 4-amino-2-hydroxybenzoic acid
        
        cfg_no_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(mixed_smiles, cfg_no_folding))
        
        # Count products by rule
        rule_counts = {}
        for product in products:
            rule = product['rule']
            rule_counts[rule] = rule_counts.get(rule, 0) + 1
        
        # Should have R1 (aromatic CH), R3 (OH), R4 (NH2), R5 (COOH) products
        self.assertGreater(rule_counts.get('R1', 0), 0, "Should have R1 products")
        self.assertGreater(rule_counts.get('R3', 0), 0, "Should have R3 products") 
        self.assertGreater(rule_counts.get('R4', 0), 0, "Should have R4 products")
        self.assertGreater(rule_counts.get('R5', 0), 0, "Should have R5 products")
        
        # Verify each rule type has correct site assignments
        r3_products = [p for p in products if p['rule'] == 'R3']
        r4_products = [p for p in products if p['rule'] == 'R4'] 
        r5_products = [p for p in products if p['rule'] == 'R5']
        
        # Each should have exactly 1 product (one site each)
        self.assertEqual(len(r3_products), 1, "Should have 1 R3 product (1 OH site)")
        self.assertEqual(len(r4_products), 1, "Should have 1 R4 product (1 NH2 site)")
        self.assertEqual(len(r5_products), 1, "Should have 1 R5 product (1 COOH site)")
        
        # Verify sites are different
        all_sites = []
        for products_list in [r3_products, r4_products, r5_products]:
            for p in products_list:
                all_sites.append(p['substitutions'][0]['site'])
        
        unique_sites = set(all_sites)
        self.assertEqual(len(unique_sites), len(all_sites), 
                        "All functional group sites should be different")

    def test_per_ring_quota_with_correct_sites(self):
        """Test that per_ring_quota counting uses correct ring tags from proper sites."""
        # Use a molecule where incorrect site mapping would affect quota counting
        dihydroxy_smiles = "c1cc(O)ccc1O"  # 1,4-dihydroxybenzene
        
        # Set per_ring_quota to 1 to test quota enforcement
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 1, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(dihydroxy_smiles, cfg))
        
        # With per_ring_quota=1, we should still get products from different rules
        # but not exceed quota per ring
        r1_products = [p for p in products if p['rule'] == 'R1']
        r3_products = [p for p in products if p['rule'] == 'R3']
        
        # Should have some R1 products and some R3 products
        self.assertGreater(len(r1_products), 0, "Should have R1 products")
        self.assertGreater(len(r3_products), 0, "Should have R3 products")
        
        # Total per ring should not exceed quota
        # (This test mainly ensures the system doesn't crash with correct site mapping)
        self.assertLessEqual(len(products), 10, "Should respect quota constraints")


if __name__ == '__main__':
    unittest.main()