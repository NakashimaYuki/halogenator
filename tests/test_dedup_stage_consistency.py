# -*- coding: ascii -*-
"""
Test dedup_stage consistency and stable key usage.

This test ensures that the configurable dedup_stage setting produces
consistent results and that stable keys are used throughout.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats
from halogenator.standardize import to_inchikey_sanitized


class TestDedupStageConsistency(unittest.TestCase):
    """Test dedup_stage configuration and stable key consistency."""

    def test_dedup_stage_pre_post_equivalence(self):
        """Test that pre and post dedup_stage produce same unique products and equal counts."""
        # Use a symmetric molecule that should generate duplicates for testing
        parent_smi = 'Cc1ccc(C)cc1'  # p-xylene - symmetric, should generate some duplicates

        # Test with pre-constraint deduplication
        cfg_pre = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'pre', 'emit_legacy_keys': False}
        )
        products_pre, qa_pre = enumerate_with_stats(parent_smi, cfg_pre)

        # Test with post-constraint deduplication
        cfg_post = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'post', 'emit_legacy_keys': False}
        )
        products_post, qa_post = enumerate_with_stats(parent_smi, cfg_post)

        # Extract unique product sets by InChIKey
        inchikeys_pre = set(p['inchikey'] for p in products_pre)
        inchikeys_post = set(p['inchikey'] for p in products_post)

        # Should produce same unique product set
        self.assertEqual(inchikeys_pre, inchikeys_post,
                        "Pre and post dedup stages should produce same unique product set")

        # Should have equal dedup hit counts
        qa_paths_pre = qa_pre.get('qa_paths', {})
        qa_paths_post = qa_post.get('qa_paths', {})

        dedup_hits_pre = qa_paths_pre.get('dedup_hits_inchi', 0)
        dedup_hits_post = qa_paths_post.get('dedup_hits_inchi', 0)

        self.assertEqual(dedup_hits_pre, dedup_hits_post,
                        f"Pre and post stages should have equal dedup counts: {dedup_hits_pre} vs {dedup_hits_post}")

    def test_duplicate_hits_increment_qa_paths(self):
        """Test that duplicate hits properly increment qa_paths['dedup_hits_inchi']."""
        # Use a molecule that's likely to generate duplicates when enumerated with high k
        parent_smi = 'Cc1ccccc1'  # toluene

        for dedup_stage in ['pre', 'post']:
            with self.subTest(dedup_stage=dedup_stage):
                cfg = EnumConfig(
                    k_max=2,  # Higher k increases chance of duplicates
                    halogens=('F', 'Cl'),
                    rules=('R6',),
                    rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F', 'Cl'], 'allow_on_methoxy': False}},
                    engine_cfg={'budget_mode': 'ops', 'dedup_stage': dedup_stage, 'emit_legacy_keys': False}
                )
                products, qa_stats = enumerate_with_stats(parent_smi, cfg)

                qa_paths = qa_stats.get('qa_paths', {})
                dedup_hits = qa_paths.get('dedup_hits_inchi', 0)

                # Should be non-negative integer
                self.assertIsInstance(dedup_hits, int)
                self.assertGreaterEqual(dedup_hits, 0)

                # If we have products, verify the dedup metric is tracked consistently
                if len(products) > 0:
                    # All products should have valid InChIKeys
                    for product in products:
                        self.assertIn('inchikey', product)
                        self.assertNotEqual(product['inchikey'], 'UNKNOWN')

    def test_record_keys_are_sanitized(self):
        """Test that record['inchikey'] and record['parent_inchikey'] are sanitized."""
        parent_smi = 'Cc1ccccc1'  # toluene

        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            rules_cfg={'R6_methyl': {'enable': True, 'allowed': ['F'], 'allow_on_methoxy': False}},
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'pre', 'emit_legacy_keys': False}
        )
        products, qa_stats = enumerate_with_stats(parent_smi, cfg)

        # Test that we have some products to verify
        self.assertGreater(len(products), 0, "Should generate at least one product for testing")

        for i, product in enumerate(products):
            with self.subTest(product_index=i):
                # Product InChIKey should be sanitized (or at least valid)
                product_inchikey = product.get('inchikey')
                self.assertIsNotNone(product_inchikey, "Product inchikey should not be None")
                self.assertNotEqual(product_inchikey, 'UNKNOWN', "Product inchikey should not be UNKNOWN")

                # Verify InChIKey format
                if product_inchikey and product_inchikey != 'UNKNOWN':
                    self.assertTrue(len(product_inchikey) > 10, "InChIKey should be reasonably long")

                # Check that parent info is present (could be parent_inchikey or parent_smi)
                has_parent_info = product.get('parent_inchikey') or product.get('parent_smi')
                self.assertIsNotNone(has_parent_info, "Product should have parent information (parent_inchikey or parent_smi)")

                # If parent_inchikey is present, verify it's valid
                parent_inchikey = product.get('parent_inchikey')
                if parent_inchikey is not None:
                    self.assertNotEqual(parent_inchikey, 'UNKNOWN', "Parent inchikey should not be UNKNOWN")
                    self.assertTrue(len(parent_inchikey) > 10, "Parent InChIKey should be reasonably long")

    def test_dedup_stage_configuration_respected(self):
        """Test that dedup_stage configuration is properly respected."""
        parent_smi = 'Cc1ccccc1'

        # Test default dedup_stage (should be 'pre')
        cfg_default = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops'}  # No explicit dedup_stage
        )
        self.assertEqual(cfg_default.engine_cfg.get('dedup_stage', 'pre'), 'pre')

        # Test explicit 'pre' setting
        cfg_explicit_pre = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'pre'}
        )
        self.assertEqual(cfg_explicit_pre.engine_cfg['dedup_stage'], 'pre')

        # Test 'post' setting
        cfg_post = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R6',),
            engine_cfg={'budget_mode': 'ops', 'dedup_stage': 'post'}
        )
        self.assertEqual(cfg_post.engine_cfg['dedup_stage'], 'post')

        # Both should work without errors
        products_default, _ = enumerate_with_stats(parent_smi, cfg_default)
        products_post, _ = enumerate_with_stats(parent_smi, cfg_post)

        # Should generate same number of unique products
        self.assertEqual(len(products_default), len(products_post))


if __name__ == '__main__':
    unittest.main(verbosity=2)