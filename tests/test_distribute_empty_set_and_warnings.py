# -*- coding: ascii -*-
"""Test distribute_marginals_to_2d empty set protection and warnings."""

import unittest
from src.halogenator.report import distribute_marginals_to_2d, distribute_marginals_to_2d_with_warnings


class TestDistributeEmptySetAndWarnings(unittest.TestCase):
    """Test empty set protection and warning functionality."""

    def test_distribute_signature_backward_compatible_returns_dict(self):
        """Test that original distribute function returns dict only."""
        by_rule = {'R1': {'metric1': 10}}
        by_halogen = None
        rules = ['R1']
        halogens = ['F', 'Cl']
        metrics = ['metric1']

        result = distribute_marginals_to_2d(by_rule, by_halogen, rules, halogens, metrics)

        # Should return dict, not tuple
        self.assertIsInstance(result, dict)
        self.assertIn('R1', result)
        self.assertIn('F', result['R1'])
        self.assertIn('metric1', result['R1']['F'])

    def test_distribute_with_warnings_returns_tuple_and_emits_empty_set_warning(self):
        """Test that new wrapper function returns tuple and handles empty sets."""
        by_rule = {'R1': {'metric1': 10}}
        by_halogen = None
        rules = []  # Empty rules set
        halogens = ['F', 'Cl']
        metrics = ['metric1']

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, rules, halogens, metrics
        )

        # Should return tuple
        self.assertIsInstance(result, dict)
        self.assertIsInstance(warnings, list)

        # Should have empty set warning
        self.assertEqual(len(warnings), 1)
        self.assertEqual(warnings[0]['type'], 'empty_rules_or_halogens')

        # Result should be empty (correct shape but all zeros)
        self.assertEqual(result, {})

    def test_distribute_handles_empty_rules_or_halogens_without_crash(self):
        """Test that distribution handles empty rules or halogens gracefully."""
        by_rule = {'R1': {'metric1': 10}}
        by_halogen = None
        metrics = ['metric1']

        # Test empty rules
        result = distribute_marginals_to_2d(by_rule, by_halogen, [], ['F', 'Cl'], metrics)
        self.assertEqual(result, {})

        # Test empty halogens
        result = distribute_marginals_to_2d(by_rule, by_halogen, ['R1'], [], metrics)
        expected = {'R1': {}}  # Correct shape for empty halogens
        self.assertEqual(result, expected)

        # Test both empty
        result = distribute_marginals_to_2d(by_rule, by_halogen, [], [], metrics)
        self.assertEqual(result, {})

    def test_distribute_with_warnings_object_structure(self):
        """Test that warnings have consistent object structure."""
        by_rule = None
        by_halogen = None
        rules = []  # Empty to trigger warning
        halogens = ['F']
        metrics = ['metric1']

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, rules, halogens, metrics
        )

        # Check warning structure
        self.assertEqual(len(warnings), 1)
        warning = warnings[0]
        self.assertIsInstance(warning, dict)
        self.assertIn('type', warning)
        self.assertEqual(warning['type'], 'empty_rules_or_halogens')

    def test_distribute_normal_operation_no_warnings(self):
        """Test that normal operation produces no warnings."""
        by_rule = {'R1': {'metric1': 10}}
        by_halogen = None
        rules = ['R1']
        halogens = ['F', 'Cl']
        metrics = ['metric1']

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, rules, halogens, metrics
        )

        # Should have no warnings for normal operation
        self.assertEqual(len(warnings), 0)

        # Should produce expected distribution
        self.assertEqual(result['R1']['F']['metric1'], 5)  # 10 // 2
        self.assertEqual(result['R1']['Cl']['metric1'], 5)  # 10 // 2

    def test_distribute_preserves_distribution_logic(self):
        """Test that wrapper preserves original distribution logic."""
        by_rule = {'R1': {'metric1': 7}}  # Odd number to test remainder
        by_halogen = None
        rules = ['R1']
        halogens = ['F', 'Cl', 'Br']  # 3 halogens
        metrics = ['metric1']

        # Test original function
        result_orig = distribute_marginals_to_2d(by_rule, by_halogen, rules, halogens, metrics)

        # Test wrapper function
        result_wrapper, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, rules, halogens, metrics
        )

        # Results should be identical
        self.assertEqual(result_orig, result_wrapper)

        # Should have no warnings
        self.assertEqual(len(warnings), 0)

        # Check distribution: 7 // 3 = 2 remainder 1
        # First halogen in alphabetical order (Br) gets the remainder
        self.assertEqual(result_wrapper['R1']['Br']['metric1'], 3)  # 2 + 1
        self.assertEqual(result_wrapper['R1']['Cl']['metric1'], 2)  # 2
        self.assertEqual(result_wrapper['R1']['F']['metric1'], 2)   # 2


if __name__ == '__main__':
    unittest.main()