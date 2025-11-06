# -*- coding: ascii -*-
"""
Test record InChIKey decoupled from deduplication key.

This test ensures that record['inchikey'] always contains a valid sanitized InChIKey
string, regardless of the deduplication policy being used (stable_key, state_sig, none).
This prevents pollution of record fields when deduplication strategy switches.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


class TestRecordInchikeyDecoupledFromDedupKey(unittest.TestCase):
    """Test that record InChIKey is decoupled from deduplication key."""

    def test_record_inchikey_always_stable_key_regardless_of_policy(self):
        """Test that record['inchikey'] is always sanitized InChIKey, not dedup key."""
        parent_smi = 'c1ccccc1O'  # phenol

        # Test all deduplication policies
        policies = ['stable_key', 'state_sig', 'none', 'auto']

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

                # Every product should have valid InChIKey record
                for i, product in enumerate(products):
                    with self.subTest(product_index=i):
                        ik = product.get('inchikey')

                        # Should be present and valid
                        self.assertIsNotNone(ik, f"Product {i} with policy {policy} should have inchikey")
                        self.assertIsInstance(ik, str, f"Product {i} inchikey should be string")
                        self.assertNotEqual(ik, 'UNKNOWN', f"Product {i} inchikey should not be UNKNOWN")

                        # Should be a proper InChIKey format (not state signature)
                        self.assertGreater(len(ik), 10, f"Product {i} inchikey should be proper length")

                        # InChIKey should start with standard prefix pattern
                        # (This is a basic check - real InChIKeys have specific format)
                        self.assertTrue(ik.replace('-', '').replace('_', '').isalnum(),
                                      f"Product {i} inchikey should be alphanumeric with separators")

    def test_state_sig_policy_produces_stable_inchikey_records(self):
        """Test that even with state_sig policy, records contain stable InChIKeys."""
        parent_smi = 'Cc1ccccc1'  # toluene

        # Use state_sig policy explicitly
        cfg_state_sig = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F', 'Cl'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'state_sig'}
        )
        products_state_sig, qa_state_sig = enumerate_with_stats(parent_smi, cfg_state_sig)

        # Compare with stable_key policy to ensure records are equivalent
        cfg_stable_key = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F', 'Cl'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'stable_key'}
        )
        products_stable_key, qa_stable_key = enumerate_with_stats(parent_smi, cfg_stable_key)

        # Both should generate products
        self.assertGreater(len(products_state_sig), 0, "State sig should generate products")
        self.assertGreater(len(products_stable_key), 0, "Stable key should generate products")

        # Record InChIKeys should have same format characteristics
        for products, policy_name in [(products_state_sig, "state_sig"), (products_stable_key, "stable_key")]:
            for i, product in enumerate(products):
                with self.subTest(product_index=i, policy=policy_name):
                    ik = product['inchikey']

                    # Should be proper InChIKey string format
                    self.assertIsInstance(ik, str)
                    self.assertNotEqual(ik, 'UNKNOWN')
                    self.assertGreater(len(ik), 10)

                    # Should not contain state signature patterns (which are typically shorter)
                    # This is a heuristic check since we don't know exact state sig format
                    # But InChIKeys are typically much longer than state signatures
                    self.assertGreater(len(ik), 20,
                                      f"InChIKey should be longer than typical state signature")

    def test_nonfolding_auto_policy_produces_stable_inchikey_records(self):
        """Test that auto policy with folding=False still produces stable InChIKey records."""
        parent_smi = 'Cc1ccc(C)cc1'  # p-xylene

        # Non-folding mode with auto policy (should use state_sig for dedup)
        cfg_nonfolding_auto = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': False},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'auto'}
        )
        products_nonfolding, qa_nonfolding = enumerate_with_stats(parent_smi, cfg_nonfolding_auto)

        # Should generate products
        self.assertGreater(len(products_nonfolding), 0, "Non-folding auto should generate products")

        # Even though deduplication uses state_sig, records should have stable InChIKeys
        for i, product in enumerate(products_nonfolding):
            with self.subTest(product_index=i):
                ik = product['inchikey']

                # Should be valid InChIKey (not state signature)
                self.assertIsInstance(ik, str)
                self.assertNotEqual(ik, 'UNKNOWN')
                self.assertGreater(len(ik), 10)

                # Should contain SMILES-derived product info
                self.assertIn('smiles', product)
                self.assertIsInstance(product['smiles'], str)
                self.assertGreater(len(product['smiles']), 0)

    def test_across_all_enumeration_paths(self):
        """Test InChIKey consistency across R1/R2, R6, and reaction paths."""
        parent_smi = 'c1ccc(C)cc1O'  # p-cresol (reactive on both ring and methyl)

        # Use state_sig policy to test decoupling across all paths
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1', 'R2', 'R3', 'R4', 'R5', 'R6'),  # All rule types
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'state_sig'}
        )
        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Should generate products from multiple paths
        self.assertGreater(len(products), 0, "Should generate products from multiple paths")

        # Group products by rule type
        products_by_rule = {}
        for product in products:
            rule = product.get('rule', 'unknown')
            if rule not in products_by_rule:
                products_by_rule[rule] = []
            products_by_rule[rule].append(product)

        # Verify InChIKey consistency across all rule types
        for rule, rule_products in products_by_rule.items():
            with self.subTest(rule=rule):
                for i, product in enumerate(rule_products):
                    with self.subTest(product_index=i):
                        ik = product['inchikey']

                        # Should always be stable InChIKey format
                        self.assertIsInstance(ik, str, f"Rule {rule} product {i} inchikey should be string")
                        self.assertNotEqual(ik, 'UNKNOWN', f"Rule {rule} product {i} should have valid inchikey")
                        self.assertGreater(len(ik), 10, f"Rule {rule} product {i} inchikey should be proper length")

    def test_none_policy_still_produces_valid_inchikeys(self):
        """Test that even with no deduplication, records have valid InChIKeys."""
        parent_smi = 'c1ccccc1O'  # phenol

        cfg_none = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1', 'R2'),
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'none'}
        )
        products_none, qa_none = enumerate_with_stats(parent_smi, cfg_none)

        # Should generate products (potentially with duplicates)
        self.assertGreater(len(products_none), 0, "None policy should generate products")

        # No dedup hits expected
        self.assertEqual(qa_none.get('dedup_hits_statesig', 0), 0)
        self.assertEqual(qa_none.get('dedup_hits_inchi', 0), 0)

        # But all records should still have valid InChIKeys
        for i, product in enumerate(products_none):
            with self.subTest(product_index=i):
                ik = product['inchikey']

                self.assertIsInstance(ik, str)
                self.assertNotEqual(ik, 'UNKNOWN')
                self.assertGreater(len(ik), 10)


if __name__ == '__main__':
    unittest.main(verbosity=2)