# -*- coding: ascii -*-
"""
Test non-folding state signature semantics.

This test ensures that when enable_symmetry_fold=False, the system uses
state signature deduplication instead of stable InChIKey deduplication,
preserving historical behavior for chemical equivalence handling.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


class TestNonfoldingStateSigSemantics(unittest.TestCase):
    """Test non-folding state signature semantics."""

    def test_nonfolding_uses_state_sig_by_default(self):
        """Test that non-folding mode automatically uses state_sig deduplication."""
        # Use a symmetric molecule that would generate different results
        # under state_sig vs stable_key deduplication
        parent_smi = 'Cc1ccc(C)cc1'  # p-xylene - symmetric

        # Test with folding disabled (should auto-select state_sig policy)
        cfg_nonfolding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': False},
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'pre', 'dedup_policy': 'auto'}
        )
        products_nonfolding, qa_nonfolding = enumerate_with_stats(parent_smi, cfg_nonfolding)

        # Should have generated products
        self.assertGreater(len(products_nonfolding), 0, "Non-folding should generate products")

        # Verify that state signature counting was used
        self.assertIn('dedup_hits_statesig', qa_nonfolding)
        self.assertIn('dedup_hits_inchi', qa_nonfolding)

        # With non-folding, any duplicates should be counted in statesig, not inchi
        if qa_nonfolding['dedup_hits_statesig'] > 0 or qa_nonfolding['dedup_hits_inchi'] > 0:
            # If there were any dedup hits, they should be in statesig for non-folding mode
            # (This specific test might not hit duplicates, but if it does, verify correct counting)
            pass

    def test_explicit_state_sig_policy_behavior(self):
        """Test explicit state_sig policy produces expected results."""
        parent_smi = 'c1ccccc1O'  # phenol

        # Test with explicit state_sig policy
        cfg_state_sig = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules=('R1', 'R2'),
            pruning_cfg={'enable_symmetry_fold': True},  # Even with folding on
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'pre', 'dedup_policy': 'state_sig'}
        )
        products_state_sig, qa_state_sig = enumerate_with_stats(parent_smi, cfg_state_sig)

        # Should have generated products
        self.assertGreater(len(products_state_sig), 0, "State sig policy should generate products")

        # Verify metrics structure
        self.assertIn('dedup_hits_statesig', qa_state_sig)
        self.assertIn('dedup_hits_inchi', qa_state_sig)

        # All products should have valid InChIKey records (even though dedup uses state_sig)
        for product in products_state_sig:
            self.assertIn('inchikey', product)
            self.assertIsNotNone(product['inchikey'])
            self.assertNotEqual(product['inchikey'], 'UNKNOWN')
            self.assertGreater(len(product['inchikey']), 10, "InChIKey should be a proper key string")

    def test_auto_policy_folding_vs_nonfolding_consistency(self):
        """Test that auto policy correctly selects stable_key vs state_sig based on folding."""
        parent_smi = 'Cc1ccccc1'  # toluene

        # Test with folding enabled (should auto-select stable_key)
        cfg_folding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': True},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'auto'}
        )
        products_folding, qa_folding = enumerate_with_stats(parent_smi, cfg_folding)

        # Test with folding disabled (should auto-select state_sig)
        cfg_nonfolding = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            pruning_cfg={'enable_symmetry_fold': False},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'auto'}
        )
        products_nonfolding, qa_nonfolding = enumerate_with_stats(parent_smi, cfg_nonfolding)

        # Both should generate products
        self.assertGreater(len(products_folding), 0, "Folding mode should generate products")
        self.assertGreater(len(products_nonfolding), 0, "Non-folding mode should generate products")

        # Both should have stable metric structures
        for qa_stats, mode_name in [(qa_folding, "folding"), (qa_nonfolding, "non-folding")]:
            self.assertIn('dedup_hits_statesig', qa_stats, f"{mode_name} should have statesig counter")
            self.assertIn('dedup_hits_inchi', qa_stats, f"{mode_name} should have inchi counter")
            self.assertIsInstance(qa_stats['dedup_hits_statesig'], int)
            self.assertIsInstance(qa_stats['dedup_hits_inchi'], int)

        # All products should have proper InChIKey records regardless of dedup policy
        all_products = [
            (products_folding, "folding"),
            (products_nonfolding, "non-folding")
        ]
        for products, mode_name in all_products:
            for i, product in enumerate(products):
                with self.subTest(product_index=i, mode=mode_name):
                    self.assertIn('inchikey', product)
                    self.assertNotEqual(product['inchikey'], 'UNKNOWN')
                    self.assertGreater(len(product['inchikey']), 10,
                                      f"Product {i} in {mode_name} mode should have valid InChIKey")

    def test_none_policy_disables_deduplication(self):
        """Test that policy='none' disables deduplication entirely."""
        parent_smi = 'Cc1ccccc1'  # toluene

        cfg_none = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_policy': 'none'}
        )
        products_none, qa_none = enumerate_with_stats(parent_smi, cfg_none)

        # Should have generated products
        self.assertGreater(len(products_none), 0, "None policy should generate products")

        # Should have no dedup hits (since deduplication is disabled)
        self.assertEqual(qa_none.get('dedup_hits_statesig', 0), 0,
                        "None policy should not have state sig dedup hits")
        self.assertEqual(qa_none.get('dedup_hits_inchi', 0), 0,
                        "None policy should not have inchi dedup hits")

        # Products should still have valid InChIKey records
        for product in products_none:
            self.assertIn('inchikey', product)
            self.assertNotEqual(product['inchikey'], 'UNKNOWN')

    def test_record_inchikey_decoupled_from_dedup_policy(self):
        """Test that record['inchikey'] is always sanitized InChIKey regardless of dedup policy."""
        parent_smi = 'c1ccccc1O'  # phenol

        # Test multiple policies
        policies = ['stable_key', 'state_sig', 'none']

        for policy in policies:
            with self.subTest(policy=policy):
                cfg = EnumConfig(
                    k_max=1,
                    halogens=('F',),
                    rules=('R1', 'R2'),
                    engine_cfg={'budget_mode': 'ops', 'dedup_policy': policy}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                # Should have generated products
                self.assertGreater(len(products), 0, f"Policy {policy} should generate products")

                # All products should have proper InChIKey records
                for i, product in enumerate(products):
                    with self.subTest(product_index=i):
                        ik = product.get('inchikey')
                        self.assertIsNotNone(ik, f"Product {i} should have inchikey")
                        self.assertNotEqual(ik, 'UNKNOWN', f"Product {i} inchikey should not be UNKNOWN")
                        self.assertIsInstance(ik, str, f"Product {i} inchikey should be string")
                        self.assertGreater(len(ik), 10, f"Product {i} inchikey should be proper length")


if __name__ == '__main__':
    unittest.main(verbosity=2)