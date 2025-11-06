"""
Tests for value coercion semantics and metadata consistency.

This module tests that the value coercion method annotation matches
the actual implementation behavior, particularly for edge cases like
negative numbers.
"""
import unittest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from halogenator.report import sanitize_granular_fields_with_warnings


class TestValueCoercionSemantics(unittest.TestCase):
    """Test value coercion behavior and metadata consistency."""

    def setUp(self):
        """Set up test fixtures."""
        self.distributable_set = {'isotope_unavailable', 'no_product_matches'}
        self.rules = ['R1', 'R2']
        self.halogens = ['F', 'Cl']

    def test_positive_number_coercion(self):
        """Test coercion of positive float to int (truncation toward zero)."""
        by_rule = {'R1': {'isotope_unavailable': 5.7}}

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # Should generate warning for non-integer value
        non_int_warnings = [w for w in warnings if w.get('type') == 'non_integer_value_detected']
        self.assertEqual(len(non_int_warnings), 1)

        warning = non_int_warnings[0]
        self.assertEqual(warning['value'], 5.7)
        self.assertEqual(warning['coerced_to_int'], True)
        self.assertEqual(warning['coercion_method'], 'trunc')

        # Check that value was actually truncated to 5 (not floored to 5)
        # For positive numbers, trunc and floor are the same
        self.assertEqual(by_rule['R1']['isotope_unavailable'], 5)

    def test_negative_number_coercion(self):
        """Test coercion of negative float shows truncation, not floor behavior."""
        by_halogen = {'F': {'no_product_matches': -2.7}}

        warnings = sanitize_granular_fields_with_warnings(
            None, by_halogen, None, self.distributable_set, self.rules, self.halogens
        )

        # Should generate warning for non-integer value
        non_int_warnings = [w for w in warnings if w.get('type') == 'non_integer_value_detected']
        self.assertEqual(len(non_int_warnings), 1)

        warning = non_int_warnings[0]
        self.assertEqual(warning['value'], -2.7)
        self.assertEqual(warning['coerced_to_int'], True)
        self.assertEqual(warning['coercion_method'], 'trunc')

        # Critical test: int(-2.7) = -2 (truncate toward zero)
        # floor(-2.7) would be -3 (floor toward negative infinity)
        # Our implementation uses int() which truncates toward zero
        self.assertEqual(by_halogen['F']['no_product_matches'], -2)

    def test_large_number_coercion(self):
        """Test coercion of large float values."""
        by_rule_halogen = {
            'R1': {'F': {'isotope_unavailable': 999.9}}
        }

        warnings = sanitize_granular_fields_with_warnings(
            None, None, by_rule_halogen, self.distributable_set, self.rules, self.halogens
        )

        non_int_warnings = [w for w in warnings if w.get('type') == 'non_integer_value_detected']
        self.assertEqual(len(non_int_warnings), 1)

        warning = non_int_warnings[0]
        self.assertEqual(warning['coercion_method'], 'trunc')
        self.assertEqual(by_rule_halogen['R1']['F']['isotope_unavailable'], 999)

    def test_edge_case_coercion_values(self):
        """Test edge cases for coercion behavior."""
        test_cases = [
            (0.9, 0),    # positive fractional
            (-0.9, 0),   # negative fractional
            (1.0, 1),    # exact integer as float
            (-1.0, -1),  # exact negative integer as float
        ]

        for input_val, expected_output in test_cases:
            with self.subTest(input_value=input_val):
                by_rule = {'R1': {'isotope_unavailable': input_val}}

                warnings = sanitize_granular_fields_with_warnings(
                    by_rule, None, None, self.distributable_set, self.rules, self.halogens
                )

                # All non-integer values should generate warnings
                if input_val != int(input_val):
                    non_int_warnings = [w for w in warnings if w.get('type') == 'non_integer_value_detected']
                    self.assertEqual(len(non_int_warnings), 1)
                    self.assertEqual(non_int_warnings[0]['coercion_method'], 'trunc')

                # Check actual coerced value matches int() behavior (truncation)
                self.assertEqual(by_rule['R1']['isotope_unavailable'], expected_output)

    def test_integer_values_no_coercion_warning(self):
        """Test that integer values don't generate coercion warnings."""
        by_rule = {'R1': {'isotope_unavailable': 42}}

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # No coercion warnings for integer values
        non_int_warnings = [w for w in warnings if w.get('type') == 'non_integer_value_detected']
        self.assertEqual(len(non_int_warnings), 0)

        # Value should remain unchanged
        self.assertEqual(by_rule['R1']['isotope_unavailable'], 42)


if __name__ == '__main__':
    unittest.main()