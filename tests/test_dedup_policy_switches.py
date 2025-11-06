# -*- coding: ascii -*-
"""
Test deduplication policy switches.

This test ensures that different dedup policies (auto/stable_key/state_sig/none)
produce the expected product sets and QA metrics, and that the policies work
consistently across all enumeration paths.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


class TestDedupPolicySwitches(unittest.TestCase):
    """Test deduplication policy switches and their effects."""

    def test_policy_ordering_stable_key_most_restrictive(self):
        """Test that stable_key <= state_sig <= none in terms of product count."""
        # Use a molecule that should generate some duplicates for meaningful testing
        parent_smi = 'Cc1ccc(C)cc1'  # p-xylene - symmetric

        # Test with all policies
        policies = ['stable_key', 'state_sig', 'none']
        results = {}

        for policy in policies:
            cfg = EnumConfig(
                k_max=1,
                halogens=('F', 'Cl'),
                rules=('R6',),
                rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F', 'Cl'], 'allow_on_methoxy': False}},
                engine_cfg={'budget_mode': 'ops', 'dedup_policy': policy}
            )
            products, qa_stats = enumerate_with_stats(parent_smi, cfg)
            results[policy] = {
                'product_count': len(products),
                'qa_stats': qa_stats
            }

        # Verify ordering: stable_key should produce <= state_sig <= none products
        stable_count = results['stable_key']['product_count']
        state_sig_count = results['state_sig']['product_count']
        none_count = results['none']['product_count']

        self.assertLessEqual(stable_count, state_sig_count,
                           "stable_key should produce <= products than state_sig")
        self.assertLessEqual(state_sig_count, none_count,
                           "state_sig should produce <= products than none")

        # All should generate at least some products
        for policy in policies:
            self.assertGreater(results[policy]['product_count'], 0,
                             f"Policy {policy} should generate at least one product")

    def test_dedup_metrics_behavior_by_policy(self):
        """Test that dedup metrics are correctly counted based on policy."""
        parent_smi = 'c1ccccc1O'  # phenol

        # Test each policy's metrics behavior
        test_cases = [
            {
                'policy': 'stable_key',
                'expected_metric': 'dedup_hits_inchi',
                'unexpected_metric': 'dedup_hits_statesig'
            },
            {
                'policy': 'state_sig',
                'expected_metric': 'dedup_hits_statesig',
                'unexpected_metric': None  # Both counters might be incremented
            },
            {
                'policy': 'none',
                'expected_metric': None,  # No dedup hits expected
                'unexpected_metric': None
            }
        ]

        for case in test_cases:
            with self.subTest(policy=case['policy']):
                cfg = EnumConfig(
                    k_max=2,  # Higher k to increase chance of duplicates
                    halogens=('F', 'Cl'),
                    rules=('R1', 'R2'),
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': case['policy']}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                # Should generate products
                self.assertGreater(len(products), 0, f"Policy {case['policy']} should generate products")

                # Check metrics structure
                self.assertIn('dedup_hits_statesig', qa_stats)
                self.assertIn('dedup_hits_inchi', qa_stats)
                self.assertIsInstance(qa_stats['dedup_hits_statesig'], int)
                self.assertIsInstance(qa_stats['dedup_hits_inchi'], int)

                # For 'none' policy, both counters should be 0
                if case['policy'] == 'none':
                    self.assertEqual(qa_stats['dedup_hits_statesig'], 0,
                                   "None policy should have no state sig dedup hits")
                    self.assertEqual(qa_stats['dedup_hits_inchi'], 0,
                                   "None policy should have no inchi dedup hits")

    def test_auto_policy_consistency_with_folding_setting(self):
        """Test that auto policy correctly adapts to folding settings."""
        parent_smi = 'Cc1ccccc1'  # toluene

        # Test auto with folding on (should behave like stable_key)
        cfg_auto_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': True},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'auto'}
        )
        products_auto_folding, qa_auto_folding = enumerate_with_stats(parent_smi, cfg_auto_folding)

        # Test explicit stable_key with folding on
        cfg_stable_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': True},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_stable_folding, qa_stable_folding = enumerate_with_stats(parent_smi, cfg_stable_folding)

        # Auto with folding should behave like explicit stable_key
        self.assertEqual(len(products_auto_folding), len(products_stable_folding),
                        "Auto with folding should produce same count as explicit stable_key")

        # Test auto with folding off (should behave like state_sig)
        cfg_auto_nonfolding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': False},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'auto'}
        )
        products_auto_nonfolding, qa_auto_nonfolding = enumerate_with_stats(parent_smi, cfg_auto_nonfolding)

        # Test explicit state_sig with folding off
        cfg_statesig_nonfolding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': False},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'state_sig'}
        )
        products_statesig_nonfolding, qa_statesig_nonfolding = enumerate_with_stats(parent_smi, cfg_statesig_nonfolding)

        # Auto with non-folding should behave like explicit state_sig
        self.assertEqual(len(products_auto_nonfolding), len(products_statesig_nonfolding),
                        "Auto with non-folding should produce same count as explicit state_sig")

    def test_all_policies_produce_valid_records(self):
        """Test that all policies produce valid product records with proper InChIKeys."""
        parent_smi = 'c1ccccc1O'  # phenol

        policies = ['auto', 'stable_key', 'state_sig', 'none']

        for policy in policies:
            with self.subTest(policy=policy):
                cfg = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=('R1', 'R2'),
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': policy}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                # Should generate products
                self.assertGreater(len(products), 0, f"Policy {policy} should generate products")

                # All products should have valid records
                for i, product in enumerate(products):
                    with self.subTest(product_index=i):
                        # Should have required fields
                        self.assertIn('inchikey', product, f"Product {i} should have inchikey")
                        self.assertIn('smiles', product, f"Product {i} should have smiles")

                        # InChIKey should be valid (not state signature)
                        ik = product['inchikey']
                        self.assertIsNotNone(ik, f"Product {i} inchikey should not be None")
                        self.assertNotEqual(ik, 'UNKNOWN', f"Product {i} inchikey should not be UNKNOWN")
                        self.assertIsInstance(ik, str, f"Product {i} inchikey should be string")
                        self.assertGreater(len(ik), 10, f"Product {i} inchikey should be proper length")

    def test_policy_across_multiple_enumeration_paths(self):
        """Test that policy is consistently applied across R1/R2, R6, and reaction paths."""
        parent_smi = 'c1ccc(C)cc1O'  # p-cresol (has multiple reactive sites)

        for policy in ['stable_key', 'state_sig', 'none']:
            with self.subTest(policy=policy):
                # Test with multiple rule types to exercise different paths
                cfg = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=('R1', 'R2', 'R6'),  # Mix of site and methyl rules
                    rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': policy}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                # Should generate products from multiple paths
                self.assertGreater(len(products), 0, f"Policy {policy} should generate products")

                # All products should have consistent record format
                for product in products:
                    self.assertIn('inchikey', product)
                    self.assertIn('rule', product)
                    self.assertIsInstance(product['inchikey'], str)
                    self.assertNotEqual(product['inchikey'], 'UNKNOWN')

                # Metrics should be properly structured
                self.assertIn('dedup_hits_statesig', qa_stats)
                self.assertIn('dedup_hits_inchi', qa_stats)
                self.assertIsInstance(qa_stats['dedup_hits_statesig'], int)
                self.assertIsInstance(qa_stats['dedup_hits_inchi'], int)


if __name__ == '__main__':
    unittest.main(verbosity=2)