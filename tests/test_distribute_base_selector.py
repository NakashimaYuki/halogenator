"""
Tests for unified base selection helper function.

This module tests the select_distribution_base function that provides
consistent base selection logic across all distribution code paths.
"""
import unittest
from src.halogenator.report import select_distribution_base


class TestDistributeBaseSelector(unittest.TestCase):
    """Test unified base selection logic."""

    def test_both_marginals_selects_by_rule_base(self):
        """Test that when both marginals exist, by_rule is selected as base."""
        by_rule = {'R1': {'isotope_unavailable': 10}}
        by_halogen = {'F': {'isotope_unavailable': 8}}

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should select by_rule as base
        self.assertEqual(source_by_rule, by_rule)
        self.assertIsNone(source_by_halogen)
        self.assertEqual(base_str, 'by_rule')

    def test_only_by_rule_selects_by_rule_base(self):
        """Test that when only by_rule exists, it's selected as base."""
        by_rule = {'R1': {'isotope_unavailable': 10}}
        by_halogen = None

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should select by_rule as base
        self.assertEqual(source_by_rule, by_rule)
        self.assertIsNone(source_by_halogen)
        self.assertEqual(base_str, 'by_rule')

    def test_only_by_halogen_selects_by_halogen_base(self):
        """Test that when only by_halogen exists, it's selected as base."""
        by_rule = None
        by_halogen = {'F': {'isotope_unavailable': 8}}

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should select by_halogen as base
        self.assertIsNone(source_by_rule)
        self.assertEqual(source_by_halogen, by_halogen)
        self.assertEqual(base_str, 'by_halogen')

    def test_neither_marginal_returns_none_base(self):
        """Test that when neither marginal exists, 'none' base is returned."""
        by_rule = None
        by_halogen = None

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should return none for both sources
        self.assertIsNone(source_by_rule)
        self.assertIsNone(source_by_halogen)
        self.assertEqual(base_str, 'none')

    def test_empty_dict_marginals_treated_as_none(self):
        """Test that empty dictionaries are treated as None."""
        by_rule = {}
        by_halogen = {}

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Empty dicts should be treated as falsy
        self.assertIsNone(source_by_rule)
        self.assertIsNone(source_by_halogen)
        self.assertEqual(base_str, 'none')

    def test_by_rule_priority_invariant_with_different_sizes(self):
        """Test that by_rule priority is maintained regardless of data size."""
        # by_halogen has more data, but by_rule should still take priority
        by_rule = {'R1': {'isotope_unavailable': 1}}
        by_halogen = {
            'F': {'isotope_unavailable': 5},
            'Cl': {'isotope_unavailable': 10},
            'Br': {'isotope_unavailable': 15}
        }

        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should still select by_rule despite by_halogen having more data
        self.assertEqual(source_by_rule, by_rule)
        self.assertIsNone(source_by_halogen)
        self.assertEqual(base_str, 'by_rule')

    def test_consistency_with_distribution_call(self):
        """Test that the helper produces inputs suitable for distribution function."""
        from src.halogenator.report import distribute_marginals_to_2d_with_warnings

        by_rule = {'R1': {'isotope_unavailable': 6}}
        by_halogen = {'F': {'isotope_unavailable': 10}}
        rules = ['R1']
        halogens = ['F', 'Cl']
        metrics = ['isotope_unavailable']

        # Get base selection
        source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

        # Should be able to call distribution function without error
        result, warnings = distribute_marginals_to_2d_with_warnings(
            source_by_rule, source_by_halogen, rules, halogens, metrics
        )

        # Result should be valid 2D structure
        self.assertIn('R1', result)
        self.assertIn('F', result['R1'])
        self.assertIn('Cl', result['R1'])
        self.assertIn('isotope_unavailable', result['R1']['F'])

        # Total should match base (by_rule in this case)
        actual_total = sum(
            result[r][h]['isotope_unavailable'] for r in rules for h in halogens
        )
        self.assertEqual(actual_total, 6)  # from by_rule, not by_halogen (10)

    def test_function_return_type_consistency(self):
        """Test that function always returns the correct types."""
        test_cases = [
            (None, None),
            ({'R1': {'m': 1}}, None),
            (None, {'F': {'m': 1}}),
            ({'R1': {'m': 1}}, {'F': {'m': 1}})
        ]

        for by_rule, by_halogen in test_cases:
            source_by_rule, source_by_halogen, base_str = select_distribution_base(by_rule, by_halogen)

            # Check types
            self.assertTrue(source_by_rule is None or isinstance(source_by_rule, dict))
            self.assertTrue(source_by_halogen is None or isinstance(source_by_halogen, dict))
            self.assertIsInstance(base_str, str)
            self.assertIn(base_str, ['by_rule', 'by_halogen', 'none'])

            # Check exclusivity: exactly one of source_by_rule/source_by_halogen should be non-None (or both None)
            if base_str == 'none':
                self.assertIsNone(source_by_rule)
                self.assertIsNone(source_by_halogen)
            else:
                self.assertTrue(
                    (source_by_rule is not None) != (source_by_halogen is not None),
                    f"Exactly one source should be non-None for base '{base_str}'"
                )


if __name__ == '__main__':
    unittest.main()