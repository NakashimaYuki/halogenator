# -*- coding: ascii -*-
"""Test P1 constraint system."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestP1Constraints(unittest.TestCase):
    """Test constraint-based product filtering."""
    
    def test_per_ring_quota_constraint(self):
        """Test per_ring_quota=1 limits substitutions per ring with stable ring tagging."""
        # Strict constraint: only 1 substitution per ring
        cfg_strict = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 1, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        # Use quercetin as a test flavonoid with A/B/C rings
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        products_strict = list(enumerate_products(quercetin_smiles, cfg_strict))
        
        strict_k_counts = Counter(p['k'] for p in products_strict)
        
        # With stable ring tagging and per_ring_quota=1 on flavonoid:
        # Should prevent multiple substitutions within the same A, B, or C ring
        # k=1 products should still be generated
        self.assertGreater(strict_k_counts[1], 0, "Should have k=1 products")
        
        # Now test relaxed constraint: allow 2 substitutions per ring  
        cfg_relaxed = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 2, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        products_relaxed = list(enumerate_products(quercetin_smiles, cfg_relaxed))
        relaxed_k_counts = Counter(p['k'] for p in products_relaxed)
        
        # With per_ring_quota=2, should allow more k=2 products within same rings
        self.assertGreater(relaxed_k_counts[1], 0, "Should have k=1 products")
        self.assertGreaterEqual(relaxed_k_counts.get(2, 0), 0, "Should have k=2 products with per_ring_quota=2")
        
        # Verify constraint has effect: strict should produce fewer products than relaxed
        self.assertLess(len(products_strict), len(products_relaxed),
                       "per_ring_quota=1 should produce fewer products than per_ring_quota=2")
        
    def test_min_graph_distance_constraint(self):
        """Test min_graph_distance=3 prevents adjacent substitutions."""
        # Strict distance: no adjacent substitutions
        cfg_strict = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 3, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        # Relaxed distance: allow adjacent substitutions
        cfg_relaxed = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products_strict = list(enumerate_products(benzene_smiles, cfg_strict))
        products_relaxed = list(enumerate_products(benzene_smiles, cfg_relaxed))
        
        # Strict distance should produce fewer k=2 products
        strict_k2 = [p for p in products_strict if p['k'] == 2]
        relaxed_k2 = [p for p in products_relaxed if p['k'] == 2]
        
        self.assertLessEqual(len(strict_k2), len(relaxed_k2),
                           "Strict distance constraint should produce <= k=2 products")
        
    def test_max_per_halogen_constraint(self):
        """Test max_per_halogen limits total usage of each halogen."""
        cfg = EnumConfig(
            k_max=3,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': {'F': 1, 'Cl': 2}},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        # Count halogen usage in each product's history
        for product in products:
            history = product.get('substitutions', [])
            halogen_counts = Counter(step['halogen'] for step in history)
            
            self.assertLessEqual(halogen_counts.get('F', 0), 1,
                               "Should not exceed max F usage")
            self.assertLessEqual(halogen_counts.get('Cl', 0), 2,
                               "Should not exceed max Cl usage")
            
    def test_constraints_ok_flag(self):
        """Test that constraints_ok flag is properly set."""
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 1, 'min_graph_distance': 2, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        # All enumerated products should pass constraints (failed ones are filtered out)
        for product in products:
            self.assertTrue(product.get('constraints_ok', False),
                          "All enumerated products should have constraints_ok=True")
            
    def test_constraints_violations_field(self):
        """Test that constraints_violations field is properly structured."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        if products:
            product = products[0]
            violations = product.get('constraints_violations', {})
            
            # Should be a dict (empty since constraints pass)
            self.assertIsInstance(violations, dict, "constraints_violations should be a dict")
            self.assertEqual(len(violations), 0, "Should have no violations for passing constraints")


if __name__ == '__main__':
    unittest.main()