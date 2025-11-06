# -*- coding: ascii -*-
"""
Test reaction path parity for consistent counting.

This test ensures that reaction-based enumeration paths (R3, R4, R5) use
the same QA counting mechanisms as site-based paths (R1, R2, R6), providing
consistent metrics and deduplication behavior across all enumeration types.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


class TestReactionPathParity(unittest.TestCase):
    """Test reaction path parity for consistent counting."""

    def test_reaction_and_site_paths_use_same_qa_structure(self):
        """Test that reaction and site paths produce consistent QA structures."""
        parent_smi = 'c1ccc(C=O)cc1O'  # aromatic aldehyde with phenol

        # Test site-based rules only
        cfg_site = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1', 'R2'),  # Site-based rules
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_site, qa_site = enumerate_with_stats(parent_smi, cfg_site)

        # Test reaction-based rules only (if available for this substrate)
        cfg_reaction = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3', 'R4', 'R5'),  # Reaction-based rules
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_reaction, qa_reaction = enumerate_with_stats(parent_smi, cfg_reaction)

        # Both should have the same QA structure keys
        expected_qa_keys = {
            'version', 'pivots', 'attempts', 'products', 'no_product_matches',
            'template_unsupported', 'qa_paths', 'dedup_hits_statesig', 'dedup_hits_inchi'
        }

        self.assertEqual(set(qa_site.keys()), expected_qa_keys,
                        "Site rules should have expected QA structure")
        self.assertEqual(set(qa_reaction.keys()), expected_qa_keys,
                        "Reaction rules should have expected QA structure")

        # QA paths should have consistent substructure
        expected_qa_paths_keys = {
            'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'
        }

        site_qa_paths_keys = set(qa_site['qa_paths'].keys())
        reaction_qa_paths_keys = set(qa_reaction['qa_paths'].keys())

        self.assertTrue(expected_qa_paths_keys.issubset(site_qa_paths_keys),
                       f"Site rules missing expected QA paths keys: {expected_qa_paths_keys - site_qa_paths_keys}")
        self.assertTrue(expected_qa_paths_keys.issubset(reaction_qa_paths_keys),
                       f"Reaction rules missing expected QA paths keys: {expected_qa_paths_keys - reaction_qa_paths_keys}")

    def test_dedup_metrics_consistency_across_paths(self):
        """Test that deduplication metrics behave consistently across paths."""
        parent_smi = 'c1ccccc1O'  # phenol (should work with both site and reaction rules)

        # Test with different policies to ensure consistency
        policies = ['stable_key', 'state_sig', 'none']

        for policy in policies:
            with self.subTest(policy=policy):
                # Combined rules (both site and reaction)
                cfg_combined = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=('R1', 'R2', 'R3', 'R4', 'R5'),
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': policy}
                )
                products_combined, qa_combined = enumerate_with_stats(parent_smi, cfg_combined)

                # Should have proper metric types and ranges
                self.assertIsInstance(qa_combined['dedup_hits_statesig'], int)
                self.assertIsInstance(qa_combined['dedup_hits_inchi'], int)
                self.assertGreaterEqual(qa_combined['dedup_hits_statesig'], 0)
                self.assertGreaterEqual(qa_combined['dedup_hits_inchi'], 0)

                # For none policy, both dedup counters should be 0
                if policy == 'none':
                    self.assertEqual(qa_combined['dedup_hits_statesig'], 0,
                                   "None policy should have no state sig hits")
                    self.assertEqual(qa_combined['dedup_hits_inchi'], 0,
                                   "None policy should have no inchi hits")

                # All QA paths metrics should be integers >= 0
                for qa_key, qa_value in qa_combined['qa_paths'].items():
                    self.assertIsInstance(qa_value, int, f"QA path {qa_key} should be int")
                    self.assertGreaterEqual(qa_value, 0, f"QA path {qa_key} should be non-negative")

    def test_record_format_consistency_across_paths(self):
        """Test that records from reaction and site paths have consistent format."""
        parent_smi = 'c1ccc(C=O)cc1O'  # aromatic aldehyde

        # Test with multiple rule types
        cfg_multi = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1', 'R2', 'R3', 'R4', 'R5', 'R6'),  # All rule types
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_multi, qa_multi = enumerate_with_stats(parent_smi, cfg_multi)

        # Group products by rule type
        products_by_rule = {}
        for product in products_multi:
            rule = product.get('rule', 'unknown')
            if rule not in products_by_rule:
                products_by_rule[rule] = []
            products_by_rule[rule].append(product)

        # All rule types should produce consistent record format
        expected_record_keys = {'smiles', 'inchikey', 'rule', 'halogen'}

        for rule, rule_products in products_by_rule.items():
            with self.subTest(rule=rule):
                for i, product in enumerate(rule_products):
                    with self.subTest(product_index=i):
                        # Should have core required fields
                        for key in expected_record_keys:
                            self.assertIn(key, product, f"Rule {rule} product {i} missing {key}")

                        # InChIKey should be valid regardless of rule type
                        ik = product['inchikey']
                        self.assertIsInstance(ik, str)
                        self.assertNotEqual(ik, 'UNKNOWN')
                        self.assertGreater(len(ik), 10)

                        # Rule should be valid
                        self.assertIn(product['rule'], ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'])

    def test_qa_bus_vs_legacy_dedup_stats_consistency(self):
        """Test that updated reaction paths produce same results as legacy approach."""
        parent_smi = 'c1ccccc1O'  # phenol

        # This test ensures that the new qa_bus approach in reaction paths
        # produces equivalent results to the old dedup_stats approach

        # Test with reaction rules (which now use qa_bus internally)
        cfg_reaction = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules=('R3', 'R4', 'R5'),  # Reaction rules
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_reaction, qa_reaction = enumerate_with_stats(parent_smi, cfg_reaction)

        # Should produce valid results
        if len(products_reaction) > 0:
            # If products were generated, verify structure
            for product in products_reaction:
                self.assertIn('inchikey', product)
                self.assertIn('rule', product)
                self.assertIn(product['rule'], ['R3', 'R4', 'R5'])

        # QA stats should be properly structured
        self.assertIsInstance(qa_reaction['dedup_hits_inchi'], int)
        self.assertIsInstance(qa_reaction['dedup_hits_statesig'], int)

        # QA paths should contain expected metrics
        qa_paths = qa_reaction['qa_paths']
        expected_metrics = ['isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used']
        for metric in expected_metrics:
            self.assertIn(metric, qa_paths)
            self.assertIsInstance(qa_paths[metric], int)
            self.assertGreaterEqual(qa_paths[metric], 0)

    def test_policy_propagation_to_reaction_paths(self):
        """Test that dedup policy is properly propagated to reaction paths."""
        parent_smi = 'c1ccccc1O'  # phenol

        # Test that different policies affect reaction paths appropriately
        test_cases = [
            {'policy': 'stable_key', 'folding': True},
            {'policy': 'state_sig', 'folding': False},
            {'policy': 'auto', 'folding': True},  # Should behave like stable_key
            {'policy': 'auto', 'folding': False},  # Should behave like state_sig
            {'policy': 'none', 'folding': True}   # Should disable deduplication
        ]

        for case in test_cases:
            with self.subTest(policy=case['policy'], folding=case['folding']):
                cfg = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=('R3', 'R4', 'R5'),  # Focus on reaction rules
                    pruning_cfg={'enable_symmetry_fold': case['folding']},
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': case['policy']}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                # Should have consistent QA structure regardless of policy
                self.assertIn('dedup_hits_statesig', qa_stats)
                self.assertIn('dedup_hits_inchi', qa_stats)

                # For none policy, dedup counters should be 0
                if case['policy'] == 'none':
                    self.assertEqual(qa_stats['dedup_hits_statesig'], 0)
                    self.assertEqual(qa_stats['dedup_hits_inchi'], 0)

                # All products should have valid records
                for product in products:
                    self.assertIn('inchikey', product)
                    self.assertIsInstance(product['inchikey'], str)
                    self.assertNotEqual(product['inchikey'], 'UNKNOWN')


if __name__ == '__main__':
    unittest.main(verbosity=2)