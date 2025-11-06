# -*- coding: ascii -*-
"""Test R2a/R2b halogen counts QA consistency."""

import unittest
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestR2HalogenCountsQA(unittest.TestCase):
    """Test that R2a/R2b QA counters respect allowed halogens configuration."""

    def test_r2_halogen_counts_respect_config(self):
        """Test that by_rule_halogen counters only include configured halogens."""
        # Use a simple flavonoid that should trigger R2 rules
        quercetin_smiles = "c1cc(c2c(c1)oc(c(c2=O)O)c3ccc(c(c3)O)O)O"

        # Configure with limited halogen set
        cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),  # Only F and Cl allowed
            rules=('R2',)  # Only R2 rules
        )

        products, qa_stats = enumerate_with_stats(quercetin_smiles, cfg)

        # Check if QA stats contain by_rule_halogen_counts (depends on implementation)
        if 'by_rule_halogen_counts' in qa_stats:
            by_rule_halogen_counts = qa_stats['by_rule_halogen_counts']

            # Both R2a and R2b should be present
            self.assertIn('R2a', by_rule_halogen_counts, "R2a should be present in halogen counts")
            self.assertIn('R2b', by_rule_halogen_counts, "R2b should be present in halogen counts")

            # Each subrule should only have keys for configured halogens
            for subrule in ('R2a', 'R2b'):
                halogen_keys = set(by_rule_halogen_counts[subrule].keys())
                expected_halogens = set(cfg.halogens)
                self.assertEqual(halogen_keys, expected_halogens,
                               f"{subrule} halogen keys {halogen_keys} should match config {expected_halogens}")

                # All counts should be non-negative integers
                for halogen, count in by_rule_halogen_counts[subrule].items():
                    self.assertIsInstance(count, int, f"{subrule}[{halogen}] should be int")
                    self.assertGreaterEqual(count, 0, f"{subrule}[{halogen}] should be non-negative")

    def test_r2_halogen_counts_no_extra_halogens(self):
        """Test that halogen counts don't include extra halogens not in config."""
        # Use a molecule that might naturally generate products with various halogens
        benzene_smiles = "c1ccccc1"

        # Configure with very limited halogen set - only Fluorine
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),  # Only F allowed
            rules=('R1', 'R2')  # Include R1 and R2 to test intersection
        )

        products, qa_stats = enumerate_with_stats(benzene_smiles, cfg)

        # Verify no products contain disallowed halogens
        for product in products:
            halogen = product.get('halogen')
            if halogen:
                self.assertIn(halogen, cfg.halogens,
                            f"Product halogen '{halogen}' not in config halogens {cfg.halogens}")

        # Check QA stats halogen consistency if available
        if 'by_rule_halogen_counts' in qa_stats:
            by_rule_halogen_counts = qa_stats['by_rule_halogen_counts']

            for subrule, halogen_counts in by_rule_halogen_counts.items():
                if subrule in ('R2a', 'R2b'):
                    halogen_keys = set(halogen_counts.keys())
                    self.assertEqual(halogen_keys, {'F'},
                                   f"{subrule} should only have F, got {halogen_keys}")

    def test_r2_halogen_counts_structure_consistency(self):
        """Test that R2a/R2b structure is consistent regardless of products."""
        # Test with a molecule that may not generate R2 products
        simple_alcohol = "CCO"  # Simple alcohol, unlikely to trigger R2

        cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl', 'Br'),
            rules=('R2', 'R3')  # Include R3 to generate some products
        )

        products, qa_stats = enumerate_with_stats(simple_alcohol, cfg)

        # Even if no R2 products are generated, the QA structure should be consistent
        if 'by_rule_halogen_counts' in qa_stats:
            by_rule_halogen_counts = qa_stats['by_rule_halogen_counts']

            # Both R2a and R2b should be present with zero counts
            for subrule in ('R2a', 'R2b'):
                self.assertIn(subrule, by_rule_halogen_counts,
                             f"{subrule} should be present even with zero products")

                # Should have entries for all configured halogens
                halogen_keys = set(by_rule_halogen_counts[subrule].keys())
                expected_halogens = set(cfg.halogens)
                self.assertEqual(halogen_keys, expected_halogens,
                               f"{subrule} should have all config halogens even with zero counts")


if __name__ == '__main__':
    unittest.main()