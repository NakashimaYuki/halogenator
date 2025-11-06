#!/usr/bin/env python3
"""
Test mixed warning type sorting with stable, deterministic results.

Verifies that the _warn_sort_key function correctly handles all warning types
regardless of their internal structure (fields at top level vs in 'where' field).
"""
import sys
import os
import unittest

# Ensure testing bootstrap is available for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import testing_bootstrap

from halogenator.report import _warn_sort_key


class TestWarningSortKeyMixed(unittest.TestCase):
    """Test warning sorting with mixed warning types."""

    def test_warn_sort_key_extraction(self):
        """Test _warn_sort_key extracts fields consistently across warning types."""

        # Test non_integer_value_detected (fields in 'where')
        warning_1 = {
            "type": "non_integer_value_detected",
            "where": {
                "dimension": "by_rule",
                "key": "R1",
                "metric": "products"
            },
            "value": 2.7,
            "coerced_to_int": True
        }

        # Test unknown_metric_dropped (fields at top level)
        warning_2 = {
            "type": "unknown_metric_dropped",
            "dimension": "by_rule",
            "key": "R1",
            "metrics": ["unknown_metric", "another_unknown"]
        }

        # Test marginal_conflict (different structure)
        warning_3 = {
            "type": "marginal_conflict",
            "base": "by_rule",
            "tolerance": 2,
            "delta": {"R1": 1}
        }

        # Extract sort keys
        key_1 = _warn_sort_key(warning_1)
        key_2 = _warn_sort_key(warning_2)
        key_3 = _warn_sort_key(warning_3)

        # Verify structure and contents
        self.assertEqual(key_1, ("non_integer_value_detected", "by_rule", "R1", "products"))
        self.assertEqual(key_2, ("unknown_metric_dropped", "by_rule", "R1", "another_unknown"))  # sorted metrics
        self.assertEqual(key_3, ("marginal_conflict", "by_rule", "", ""))  # uses base field for dimension

        # Verify sorting order is stable
        all_keys = [key_1, key_2, key_3]
        sorted_keys = sorted(all_keys)
        expected_order = [
            ("marginal_conflict", "by_rule", "", ""),  # marginal_conflict with base
            ("non_integer_value_detected", "by_rule", "R1", "products"),
            ("unknown_metric_dropped", "by_rule", "R1", "another_unknown")
        ]
        self.assertEqual(sorted_keys, expected_order)

    def test_mixed_warnings_stable_sort(self):
        """Test that mixed warning types sort stably and deterministically."""

        # Create mixed warning list with predictable sort order
        warnings = [
            # This should be 3rd (unknown_metric_dropped, by_rule, R2, metric_a)
            {
                "type": "unknown_metric_dropped",
                "dimension": "by_rule",
                "key": "R2",
                "metrics": ["metric_a", "metric_b"]
            },
            # This should be 1st (marginal_conflict has empty dimension/key/metric)
            {
                "type": "marginal_conflict",
                "base": "by_rule",
                "tolerance": 2,
                "delta": {"R1": 1}
            },
            # This should be 2nd (non_integer, by_rule, R1, products)
            {
                "type": "non_integer_value_detected",
                "where": {
                    "dimension": "by_rule",
                    "key": "R1",
                    "metric": "products"
                },
                "value": 3.14,
                "coerced_to_int": True
            },
            # This should be 4th (unknown_metric_dropped, by_rule, R2, metric_z)
            {
                "type": "unknown_metric_dropped",
                "dimension": "by_rule",
                "key": "R2",
                "metrics": ["metric_z"]
            }
        ]

        # Sort using _warn_sort_key
        sorted_warnings = sorted(warnings, key=_warn_sort_key)

        # Verify exact order
        expected_types = [
            "marginal_conflict",
            "non_integer_value_detected",
            "unknown_metric_dropped",  # metric_a comes before metric_z
            "unknown_metric_dropped"
        ]

        actual_types = [w["type"] for w in sorted_warnings]
        self.assertEqual(actual_types, expected_types)

        # Verify detailed sort keys
        actual_keys = [_warn_sort_key(w) for w in sorted_warnings]
        expected_keys = [
            ("marginal_conflict", "by_rule", "", ""),  # uses base field for dimension
            ("non_integer_value_detected", "by_rule", "R1", "products"),
            ("unknown_metric_dropped", "by_rule", "R2", "metric_a"),  # sorted metrics: metric_a < metric_b
            ("unknown_metric_dropped", "by_rule", "R2", "metric_z")
        ]
        self.assertEqual(actual_keys, expected_keys)

    def test_missing_fields_handle_gracefully(self):
        """Test that missing fields are handled gracefully with empty string defaults."""

        # Warning with minimal fields
        minimal_warning = {
            "type": "custom_warning"
        }

        # Warning with partial fields
        partial_warning = {
            "type": "partial_warning",
            "dimension": "by_rule"
        }

        # Test extraction
        key_1 = _warn_sort_key(minimal_warning)
        key_2 = _warn_sort_key(partial_warning)

        self.assertEqual(key_1, ("custom_warning", "", "", ""))
        self.assertEqual(key_2, ("partial_warning", "by_rule", "", ""))

        # Verify they sort predictably
        self.assertTrue(key_1 < key_2)  # "custom" < "partial"

    def test_metrics_array_sorted_first_element(self):
        """Test that sorted first element of metrics array is used for sort key."""

        warning_with_metrics = {
            "type": "unknown_metric_dropped",
            "dimension": "by_rule",
            "key": "R1",
            "metrics": ["zebra_metric", "alpha_metric", "beta_metric"]
        }

        key = _warn_sort_key(warning_with_metrics)
        # Should use sorted first element: alpha_metric (not zebra_metric)
        self.assertEqual(key, ("unknown_metric_dropped", "by_rule", "R1", "alpha_metric"))

        # Test empty metrics array
        warning_empty_metrics = {
            "type": "unknown_metric_dropped",
            "dimension": "by_rule",
            "key": "R1",
            "metrics": []
        }

        key_empty = _warn_sort_key(warning_empty_metrics)
        self.assertEqual(key_empty, ("unknown_metric_dropped", "by_rule", "R1", ""))

    def test_deterministic_sorting_under_random_input_order(self):
        """Test that warnings sort deterministically regardless of input order."""
        import random

        # Create comprehensive test warnings covering all identified edge cases
        warnings = [
            # unknown_metric_dropped with metrics in random order
            {
                "type": "unknown_metric_dropped",
                "dimension": "by_rule",
                "key": "R1",
                "metrics": ["charlie_metric", "alpha_metric", "bravo_metric"]
            },
            # non_integer_value_detected with fields in 'where'
            {
                "type": "non_integer_value_detected",
                "where": {
                    "dimension": "by_halogen",
                    "key": "Cl",
                    "metric": "products"
                },
                "value": 2.7,
                "coerced_to_int": True
            },
            # marginal_conflict with base=by_rule
            {
                "type": "marginal_conflict",
                "base": "by_rule",
                "tolerance": 2,
                "delta": {"R1": 1}
            },
            # marginal_conflict with base=by_halogen (different base)
            {
                "type": "marginal_conflict",
                "base": "by_halogen",
                "tolerance": 1,
                "delta": {"Cl": 2}
            },
            # unknown_metric_dropped with different metrics order
            {
                "type": "unknown_metric_dropped",
                "dimension": "by_halogen",
                "key": "Br",
                "metrics": ["zulu_metric", "delta_metric"]
            }
        ]

        # Compute expected sorted keys
        expected_sorted_keys = [_warn_sort_key(w) for w in sorted(warnings, key=_warn_sort_key)]

        # Test deterministic sorting under random shuffling (20 iterations as specified)
        for i in range(20):
            # Randomly shuffle the input warnings
            shuffled_warnings = warnings.copy()
            random.shuffle(shuffled_warnings)

            # Sort and extract keys
            actual_sorted_keys = [_warn_sort_key(w) for w in sorted(shuffled_warnings, key=_warn_sort_key)]

            # Should always produce the same result regardless of input order
            self.assertEqual(actual_sorted_keys, expected_sorted_keys,
                           f"Sorting failed to be deterministic on iteration {i+1}")

        # Verify expected key structure matches specification
        expected_detailed_keys = [
            ("marginal_conflict", "by_halogen", "", ""),      # base=by_halogen sorts before by_rule
            ("marginal_conflict", "by_rule", "", ""),         # base=by_rule
            ("non_integer_value_detected", "by_halogen", "Cl", "products"),
            ("unknown_metric_dropped", "by_halogen", "Br", "delta_metric"),  # delta < zulu after sorting
            ("unknown_metric_dropped", "by_rule", "R1", "alpha_metric")      # alpha < bravo < charlie after sorting
        ]

        self.assertEqual(expected_sorted_keys, expected_detailed_keys)


if __name__ == '__main__':
    unittest.main()