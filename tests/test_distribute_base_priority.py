"""
Tests for base priority enforcement in 2D distribution.

This module tests that when both by_rule and by_halogen marginals exist,
by_rule takes priority as the base for 2D distribution.
"""
import unittest
from src.halogenator.report import distribute_marginals_to_2d_with_warnings


class TestDistributeBasePriority(unittest.TestCase):
    def setUp(self):
        self.rules = ['R1', 'R2']
        self.halogens = ['F', 'Cl']
        self.metrics = ['products', 'no_product_matches']

    def test_by_rule_priority_over_by_halogen(self):
        """Test that by_rule takes priority when both marginals exist with different totals."""
        # Set up conflicting marginals
        by_rule = {
            'R1': {'products': 10, 'no_product_matches': 2},
            'R2': {'products': 8, 'no_product_matches': 1}
        }
        by_halogen = {
            'F': {'products': 12, 'no_product_matches': 2},
            'Cl': {'products': 4, 'no_product_matches': 3}  # Different totals than by_rule
        }

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, self.rules, self.halogens, self.metrics
        )

        # Calculate expected totals from by_rule (should match 2D aggregation)
        expected_by_rule_total_products = 10 + 8  # 18
        expected_by_rule_total_no_product = 2 + 1  # 3

        # Calculate actual totals from 2D result
        actual_total_products = sum(
            result[r][h]['products'] for r in self.rules for h in self.halogens
        )
        actual_total_no_product = sum(
            result[r][h]['no_product_matches'] for r in self.rules for h in self.halogens
        )

        # Assert that 2D totals match by_rule totals (not by_halogen)
        self.assertEqual(actual_total_products, expected_by_rule_total_products)
        self.assertEqual(actual_total_no_product, expected_by_rule_total_no_product)

        # Verify that by_halogen totals are different (confirming conflict existed)
        by_halogen_total_products = 12 + 4  # 16 != 18
        by_halogen_total_no_product = 2 + 3  # 5 != 3
        self.assertNotEqual(by_halogen_total_products, expected_by_rule_total_products)
        self.assertNotEqual(by_halogen_total_no_product, expected_by_rule_total_no_product)

    def test_by_rule_only_uses_by_rule(self):
        """Test that when only by_rule exists, it's used correctly."""
        by_rule = {
            'R1': {'products': 6, 'no_product_matches': 1},
            'R2': {'products': 4, 'no_product_matches': 2}
        }

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, None, self.rules, self.halogens, self.metrics
        )

        # Check totals match
        actual_total_products = sum(
            result[r][h]['products'] for r in self.rules for h in self.halogens
        )
        actual_total_no_product = sum(
            result[r][h]['no_product_matches'] for r in self.rules for h in self.halogens
        )

        self.assertEqual(actual_total_products, 6 + 4)  # 10
        self.assertEqual(actual_total_no_product, 1 + 2)  # 3

    def test_by_halogen_only_uses_by_halogen(self):
        """Test that when only by_halogen exists, it's used correctly."""
        by_halogen = {
            'F': {'products': 7, 'no_product_matches': 1},
            'Cl': {'products': 3, 'no_product_matches': 2}
        }

        result, warnings = distribute_marginals_to_2d_with_warnings(
            None, by_halogen, self.rules, self.halogens, self.metrics
        )

        # Check totals match
        actual_total_products = sum(
            result[r][h]['products'] for r in self.rules for h in self.halogens
        )
        actual_total_no_product = sum(
            result[r][h]['no_product_matches'] for r in self.rules for h in self.halogens
        )

        self.assertEqual(actual_total_products, 7 + 3)  # 10
        self.assertEqual(actual_total_no_product, 1 + 2)  # 3

    def test_neither_marginal_returns_zeros(self):
        """Test that when neither marginal exists, all-zero 2D is returned."""
        result, warnings = distribute_marginals_to_2d_with_warnings(
            None, None, self.rules, self.halogens, self.metrics
        )

        # Check all values are zero
        for rule in self.rules:
            for halogen in self.halogens:
                for metric in self.metrics:
                    self.assertEqual(result[rule][halogen][metric], 0)

    def test_missing_metric_no_cross_side_fallback(self):
        """
        Test that when base=by_rule, metrics missing from by_rule are not filled
        from by_halogen values. Missing metrics are not distributed at all.
        """
        # Use base selector to ensure only by_rule is passed to distribution
        from src.halogenator.report import select_distribution_base

        by_rule = {
            'R1': {'isotope_unavailable': 10}  # Only has 'isotope_unavailable', missing 'no_product_matches'
        }
        by_halogen = {
            'F': {'isotope_unavailable': 8, 'no_product_matches': 5},  # Has both metrics
            'Cl': {'isotope_unavailable': 6, 'no_product_matches': 3}
        }

        # Use the unified base selector (should pick by_rule and ignore by_halogen)
        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Verify base selection worked as expected
        self.assertEqual(base_str, 'by_rule')
        self.assertEqual(source_by_rule, by_rule)
        self.assertIsNone(source_by_halogen)

        result, warnings = distribute_marginals_to_2d_with_warnings(
            source_by_rule, source_by_halogen, self.rules, self.halogens, self.metrics
        )

        # isotope_unavailable should be distributed from by_rule
        if 'isotope_unavailable' in result[self.rules[0]][self.halogens[0]]:
            actual_isotope_total = sum(
                result[r][h]['isotope_unavailable'] for r in self.rules for h in self.halogens
            )
            self.assertEqual(actual_isotope_total, 10)  # from by_rule

        # no_product_matches should either be 0 or missing (not taken from by_halogen)
        if 'no_product_matches' in result[self.rules[0]][self.halogens[0]]:
            actual_no_product_total = sum(
                result[r][h]['no_product_matches'] for r in self.rules for h in self.halogens
            )
            self.assertEqual(actual_no_product_total, 0)  # NOT 8 from by_halogen

        # Contract: when using by_rule as base, by_halogen values are completely ignored

    def test_empty_rules_or_halogens_returns_warning(self):
        """Test that empty rules or halogens returns appropriate warning."""
        by_rule = {'R1': {'products': 5}}

        # Test empty rules
        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, None, [], self.halogens, self.metrics
        )
        self.assertEqual(len(warnings), 1)
        self.assertEqual(warnings[0]['type'], 'empty_rules_or_halogens')

        # Test empty halogens
        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, None, self.rules, [], self.metrics
        )
        self.assertEqual(len(warnings), 1)
        self.assertEqual(warnings[0]['type'], 'empty_rules_or_halogens')

    def test_2d_complete_metrics_keys_consistency(self):
        """Test that 2D structure includes complete metric keys, with missing keys set to 0."""
        # Base with incomplete metrics
        by_rule = {
            'R1': {'isotope_unavailable': 10},  # Missing 'no_product_matches'
            'R2': {'isotope_unavailable': 8, 'no_product_matches': 2}  # Complete metrics
        }

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, None, self.rules, self.halogens, self.metrics
        )

        # All 2D cells should have complete metric key sets
        expected_metrics = set(self.metrics)
        for rule in self.rules:
            for halogen in self.halogens:
                cell_metrics = set(result[rule][halogen].keys())
                self.assertEqual(cell_metrics, expected_metrics,
                               f"Cell [{rule}][{halogen}] should have complete metrics")

                # Missing metrics should be 0
                if rule == 'R1':  # R1 was missing 'no_product_matches' in base
                    self.assertEqual(result[rule][halogen]['no_product_matches'], 0,
                                   f"Missing metric should be 0 for cell [{rule}][{halogen}]")

    def test_2d_metrics_completeness_all_cells_identical_keys(self):
        """Test that all 2D cells have identical metric key sets."""
        # Mix of complete and incomplete base metrics
        by_rule = {
            'R1': {'isotope_unavailable': 5, 'no_product_matches': 1},  # Complete
            'R2': {'isotope_unavailable': 3}  # Missing 'no_product_matches'
        }

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, None, self.rules, self.halogens, self.metrics
        )

        # Collect metric key sets from all cells
        all_metric_sets = []
        for rule in self.rules:
            for halogen in self.halogens:
                cell_metrics = set(result[rule][halogen].keys())
                all_metric_sets.append(cell_metrics)

        # All cells should have identical metric key sets
        first_set = all_metric_sets[0]
        for i, metric_set in enumerate(all_metric_sets):
            self.assertEqual(metric_set, first_set,
                           f"Cell {i} has different metric keys: {metric_set} vs {first_set}")


if __name__ == '__main__':
    unittest.main()