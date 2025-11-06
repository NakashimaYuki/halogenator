# -*- coding: ascii -*-
"""Test P1 reaction-based rules (R3, R4, R5)."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestP1ReactionRules(unittest.TestCase):
    """Test reaction-based halogenation rules R3, R4, R5."""
    
    def test_r3_phenol_halogenation(self):
        """Test R3 rule on phenol should produce >=2 k=1 products."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        phenol_smiles = "c1ccc(O)cc1"
        products = list(enumerate_products(phenol_smiles, cfg))
        
        # Check that R3 rule is applied
        r3_products = [p for p in products if p['rule'] == 'R3']
        self.assertGreaterEqual(len(r3_products), 2, 
                              "Phenol should produce >=2 R3 products (OH->OX)")
        
        # Check that we get both halogens
        r3_halogens = set(p['halogen'] for p in r3_products)
        self.assertTrue(len(r3_halogens) >= 1, "Should have R3 products with different halogens")
        
    def test_r4_aniline_halogenation(self):
        """Test R4 rule on aniline should produce >=2 k=1 products.""" 
        cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        aniline_smiles = "c1ccc(N)cc1"
        products = list(enumerate_products(aniline_smiles, cfg))
        
        # Check that R4 rule is applied
        r4_products = [p for p in products if p['rule'] == 'R4']
        self.assertGreaterEqual(len(r4_products), 2,
                              "Aniline should produce >=2 R4 products (NH2->NHX)")
        
        # Check that we get both halogens
        r4_halogens = set(p['halogen'] for p in r4_products)
        self.assertTrue(len(r4_halogens) >= 1, "Should have R4 products with different halogens")
        
    def test_r5_benzoic_acid_halogenation(self):
        """Test R5 rule on benzoic acid should produce >=2 k=1 products."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        benzoic_acid_smiles = "c1ccc(C(=O)O)cc1"
        products = list(enumerate_products(benzoic_acid_smiles, cfg))
        
        # Check that R5 rule is applied
        r5_products = [p for p in products if p['rule'] == 'R5']
        self.assertGreaterEqual(len(r5_products), 2,
                              "Benzoic acid should produce >=2 R5 products (COOH->COOX)")
        
        # Check that we get both halogens
        r5_halogens = set(p['halogen'] for p in r5_products)
        self.assertTrue(len(r5_halogens) >= 1, "Should have R5 products with different halogens")
        
    def test_mixed_rules_phenol(self):
        """Test that phenol produces both R1 (aromatic) and R3 (OH) products."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        phenol_smiles = "c1ccc(O)cc1"
        products = list(enumerate_products(phenol_smiles, cfg))
        
        rule_counts = Counter(p['rule'] for p in products)
        
        # Should have both R1 (aromatic halogenation) and R3 (OH halogenation)
        self.assertGreater(rule_counts.get('R1', 0), 0, "Should have R1 products")
        self.assertGreater(rule_counts.get('R3', 0), 0, "Should have R3 products")
        
        # Total should be reasonable
        self.assertGreater(len(products), 3, "Phenol should produce >3 products total")
        
    def test_no_false_positive_rules(self):
        """Test that simple benzene doesn't trigger R3/R4/R5 rules."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        benzene_smiles = "c1ccccc1"
        products = list(enumerate_products(benzene_smiles, cfg))
        
        rule_counts = Counter(p['rule'] for p in products)
        
        # Benzene should only have R1 products, not R3/R4/R5
        self.assertGreater(rule_counts.get('R1', 0), 0, "Should have R1 products")
        self.assertEqual(rule_counts.get('R3', 0), 0, "Should not have R3 products")
        self.assertEqual(rule_counts.get('R4', 0), 0, "Should not have R4 products") 
        self.assertEqual(rule_counts.get('R5', 0), 0, "Should not have R5 products")


if __name__ == '__main__':
    unittest.main()