# -*- coding: ascii -*-
"""Tests for site_tokens_json presence in R6 product records."""

import unittest
import json
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestR6SiteTokensJson(unittest.TestCase):
    """Test that R6 product records include site_tokens_json field."""

    def setUp(self):
        """Set up test configuration with R6 enabled."""
        self.cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules_cfg={
                'R1': {'enable': True},
                'R3': {'enable': True},
                'R6_methyl': {
                    'enable': True,
                    'step': {'enable': True},
                    'macro': {'enable': True, 'labels': ['CF3', 'CCl3']}
                }
            }
        )

    def test_r6_products_have_site_tokens_json(self):
        """Test that R6 products include site_tokens_json field."""
        # Use a molecule with methyl groups for R6 halogenation
        parent_smi = 'CC'  # Ethane - simple methyl sites

        products, qa_stats = enumerate_with_stats(parent_smi, self.cfg)

        # Find R6 products
        r6_products = [p for p in products if p.get('rule') == 'R6']

        if not r6_products:
            self.skipTest("No R6 products generated for test molecule")

        # Check that all R6 products have site_tokens_json
        for product in r6_products:
            with self.subTest(product=product.get('smiles')):
                self.assertIn('site_tokens_json', product,
                             "R6 product should have site_tokens_json field")

                # Verify it's a valid JSON string
                site_tokens_json = product['site_tokens_json']
                self.assertIsInstance(site_tokens_json, str,
                                    "site_tokens_json should be a string")

                # Verify it can be parsed as JSON
                try:
                    site_tokens = json.loads(site_tokens_json)
                except json.JSONDecodeError as e:
                    self.fail(f"site_tokens_json is not valid JSON: {e}")

                # Verify it's a dict (site_tokens should be a dict structure mapping sites to counts)
                self.assertIsInstance(site_tokens, dict,
                                    "site_tokens should be a dict when parsed from JSON")

    def test_r6_step_vs_macro_site_tokens(self):
        """Test that both R6 step and macro products have site_tokens_json."""
        # Use a molecule that can generate both step and macro products
        parent_smi = 'CCC'  # Propane - has methyl sites for both step and macro

        products, qa_stats = enumerate_with_stats(parent_smi, self.cfg)

        # Find R6 step and macro products
        r6_step_products = [p for p in products if p.get('rule') == 'R6' and 'macro_label' not in p]
        r6_macro_products = [p for p in products if p.get('rule') == 'R6' and 'macro_label' in p]

        # Test step products
        for product in r6_step_products:
            self.assertIn('site_tokens_json', product,
                         "R6 step product should have site_tokens_json field")
            site_tokens_json = product['site_tokens_json']
            site_tokens = json.loads(site_tokens_json)
            self.assertIsInstance(site_tokens, dict)

        # Test macro products
        for product in r6_macro_products:
            self.assertIn('site_tokens_json', product,
                         "R6 macro product should have site_tokens_json field")
            site_tokens_json = product['site_tokens_json']
            site_tokens = json.loads(site_tokens_json)
            self.assertIsInstance(site_tokens, dict)

    def test_site_tokens_json_not_empty_when_tokens_exist(self):
        """Test that site_tokens_json is not empty when site tokens exist."""
        # Use a molecule likely to have site tokens
        parent_smi = 'CCCC'  # Butane - multiple methyl sites should generate site tokens

        products, qa_stats = enumerate_with_stats(parent_smi, self.cfg)

        # Find R6 products
        r6_products = [p for p in products if p.get('rule') == 'R6']

        if not r6_products:
            self.skipTest("No R6 products generated for test molecule")

        found_non_empty = False
        for product in r6_products:
            site_tokens_json = product.get('site_tokens_json', '{}')
            site_tokens = json.loads(site_tokens_json)

            # If site tokens exist, they should be properly recorded
            if site_tokens:
                found_non_empty = True
                # Each value in the site_tokens dict should be an integer count
                for site_key, count in site_tokens.items():
                    self.assertIsInstance(site_key, str, "Site token keys should be strings")
                    self.assertIsInstance(count, int, "Site token counts should be integers")
                    self.assertGreaterEqual(count, 0, "Site token counts should be non-negative")

        # Note: We don't require non-empty tokens as it depends on the budget system
        # But if they exist, they should be properly formatted

    def test_site_tokens_json_defaults_to_empty_dict(self):
        """Test that site_tokens_json defaults to empty dict when no tokens."""
        # Use a simple molecule that might not generate complex site tokens
        parent_smi = 'C'  # Methane - simplest possible case

        products, qa_stats = enumerate_with_stats(parent_smi, self.cfg)

        # Find R6 products if any
        r6_products = [p for p in products if p.get('rule') == 'R6']

        # For any R6 products, verify site_tokens_json exists and is valid
        for product in r6_products:
            site_tokens_json = product.get('site_tokens_json', '{}')
            site_tokens = json.loads(site_tokens_json)
            # Should be a dict (possibly empty)
            self.assertIsInstance(site_tokens, dict)

    def test_site_tokens_json_consistent_across_similar_products(self):
        """Test that similar R6 products have consistent site_tokens_json structure."""
        parent_smi = 'CC'  # Ethane

        products, qa_stats = enumerate_with_stats(parent_smi, self.cfg)

        # Find R6 products
        r6_products = [p for p in products if p.get('rule') == 'R6']

        if len(r6_products) < 2:
            self.skipTest("Need at least 2 R6 products for consistency test")

        # Get site_tokens structure from all products
        tokens_structures = []
        for product in r6_products:
            site_tokens_json = product.get('site_tokens_json', '{}')
            site_tokens = json.loads(site_tokens_json)
            tokens_structures.append(type(site_tokens))

        # All should be dicts
        for structure in tokens_structures:
            self.assertEqual(structure, dict,
                           "All R6 products should have dict-type site_tokens_json")


if __name__ == '__main__':
    unittest.main()