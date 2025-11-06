"""
Tests for data quality warnings in sanitize system.

This module tests the enhanced sanitize_granular_fields_with_warnings function
that provides structured warnings for data quality issues.
"""
import unittest
import tempfile
import os
import json
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from halogenator.report import write_qa_summary_json, sanitize_granular_fields_with_warnings


class TestWarningsDataQuality(unittest.TestCase):
    """Test data quality validation and warning generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.distributable_set = {'isotope_unavailable', 'no_product_matches', 'template_unsupported'}
        self.rules = ['R1', 'R2']
        self.halogens = ['F', 'Cl']

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_negative_value_detection(self):
        """Test detection of negative values with structured warning."""
        by_rule = {'R1': {'isotope_unavailable': -5}}  # Negative value

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have negative value warning
        negative_warnings = [w for w in warnings if w['type'] == 'negative_value_detected']
        self.assertEqual(len(negative_warnings), 1)

        warning = negative_warnings[0]
        self.assertEqual(warning['type'], 'negative_value_detected')
        self.assertEqual(warning['value'], -5)
        self.assertEqual(warning['where']['dimension'], 'by_rule')
        self.assertEqual(warning['where']['key'], 'R1')
        self.assertEqual(warning['where']['metric'], 'isotope_unavailable')

    def test_non_integer_value_coercion(self):
        """Test detection and coercion of non-integer values."""
        by_halogen = {'F': {'isotope_unavailable': 5.7}}  # Float value

        warnings = sanitize_granular_fields_with_warnings(
            None, by_halogen, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have non-integer warning
        non_int_warnings = [w for w in warnings if w['type'] == 'non_integer_value_detected']
        self.assertEqual(len(non_int_warnings), 1)

        warning = non_int_warnings[0]
        self.assertEqual(warning['type'], 'non_integer_value_detected')
        self.assertEqual(warning['value'], 5.7)
        self.assertEqual(warning['coerced_to_int'], True)
        self.assertEqual(warning['coercion_method'], 'trunc')
        self.assertEqual(warning['where']['dimension'], 'by_halogen')
        self.assertEqual(warning['where']['key'], 'F')
        self.assertEqual(warning['where']['metric'], 'isotope_unavailable')

        # Value should be coerced to integer
        self.assertEqual(by_halogen['F']['isotope_unavailable'], 5)

    def test_unknown_metric_dropped(self):
        """Test detection and removal of unknown metrics."""
        by_rule = {'R1': {'isotope_unavailable': 5, 'unknown_metric1': 10, 'unknown_metric2': 15}}

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have unknown metrics warning
        unknown_warnings = [w for w in warnings if w['type'] == 'unknown_metric_dropped']
        self.assertEqual(len(unknown_warnings), 1)

        warning = unknown_warnings[0]
        self.assertEqual(warning['type'], 'unknown_metric_dropped')
        self.assertEqual(warning['dimension'], 'by_rule')
        self.assertEqual(warning['key'], 'R1')
        self.assertEqual(set(warning['metrics']), {'unknown_metric1', 'unknown_metric2'})

        # Unknown metrics should be removed
        self.assertNotIn('unknown_metric1', by_rule['R1'])
        self.assertNotIn('unknown_metric2', by_rule['R1'])
        self.assertIn('isotope_unavailable', by_rule['R1'])  # Known metric preserved

    def test_unknown_dimension_key_dropped(self):
        """Test detection and removal of unknown dimension keys."""
        by_rule = {
            'R1': {'isotope_unavailable': 5},   # Valid rule
            'R99': {'isotope_unavailable': 10}  # Invalid rule
        }

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have unknown dimension key warning
        unknown_key_warnings = [w for w in warnings if w['type'] == 'unknown_dimension_key_dropped']
        self.assertEqual(len(unknown_key_warnings), 1)

        warning = unknown_key_warnings[0]
        self.assertEqual(warning['type'], 'unknown_dimension_key_dropped')
        self.assertEqual(warning['dimension'], 'by_rule')
        self.assertEqual(warning['key'], 'R99')

        # Unknown key should be removed
        self.assertNotIn('R99', by_rule)
        self.assertIn('R1', by_rule)  # Valid key preserved

    def test_2d_structure_validation_comprehensive(self):
        """Test validation of 2D structure with various issues."""
        by_rule_halogen = {
            'R1': {
                'F': {'isotope_unavailable': -5, 'bad_metric': 20},  # Negative value + unknown metric
                'Br': {'isotope_unavailable': 3}                    # Unknown halogen (valid value)
            },
            'R99': {  # Unknown rule
                'F': {'isotope_unavailable': 10}
            }
        }

        warnings = sanitize_granular_fields_with_warnings(
            None, None, by_rule_halogen, self.distributable_set, self.rules, self.halogens
        )

        # Should have multiple types of warnings
        warning_types = [w['type'] for w in warnings]
        self.assertIn('unknown_metric_dropped', warning_types)
        self.assertIn('unknown_dimension_key_dropped', warning_types)
        self.assertIn('negative_value_detected', warning_types)

        # Verify structure cleanup
        self.assertNotIn('R99', by_rule_halogen)  # Unknown rule removed
        self.assertNotIn('Br', by_rule_halogen['R1'])  # Unknown halogen removed
        self.assertNotIn('bad_metric', by_rule_halogen['R1']['F'])  # Unknown metric removed

    def test_multiple_data_quality_issues_single_cell(self):
        """Test multiple issues in a single cell."""
        by_rule = {'R1': {'isotope_unavailable': -5.7}}  # Negative AND float

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, None, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have both types of warnings
        warning_types = [w['type'] for w in warnings]
        self.assertIn('negative_value_detected', warning_types)
        self.assertIn('non_integer_value_detected', warning_types)

        # Both warnings should reference the same location
        for warning in warnings:
            self.assertEqual(warning['where']['dimension'], 'by_rule')
            self.assertEqual(warning['where']['key'], 'R1')
            self.assertEqual(warning['where']['metric'], 'isotope_unavailable')

    def test_warnings_structure_consistency(self):
        """Test that all warnings have consistent object structure."""
        by_rule = {
            'R1': {'isotope_unavailable': -5.7, 'bad_metric': 10},
            'R99': {'isotope_unavailable': 5}
        }
        by_halogen = {'Zz': {'isotope_unavailable': 8.5}}  # Unknown halogen

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, by_halogen, None, self.distributable_set, self.rules, self.halogens
        )

        # All warnings should be objects with 'type' field
        for warning in warnings:
            self.assertIsInstance(warning, dict)
            self.assertIn('type', warning)
            self.assertIsInstance(warning['type'], str)

            # Type-specific structure validation
            if warning['type'] == 'negative_value_detected':
                self.assertIn('where', warning)
                self.assertIn('value', warning)
            elif warning['type'] == 'non_integer_value_detected':
                self.assertIn('where', warning)
                self.assertIn('value', warning)
                self.assertIn('coerced_to_int', warning)
            elif warning['type'] == 'unknown_metric_dropped':
                self.assertIn('dimension', warning)
                self.assertIn('key', warning)
                self.assertIn('metrics', warning)
            elif warning['type'] == 'unknown_dimension_key_dropped':
                self.assertIn('dimension', warning)
                self.assertIn('key', warning)

    def test_end_to_end_data_quality_warnings_in_qa_summary(self):
        """Test that data quality warnings appear in final QA summary."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': -5.5, 'bad_metric': 10},
                'R99': {'isotope_unavailable': 8}  # Unknown rule
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have data quality warnings in metadata
        warnings = result['metadata'].get('warnings', [])
        warning_types = [w['type'] for w in warnings]

        self.assertIn('negative_value_detected', warning_types)
        self.assertIn('non_integer_value_detected', warning_types)
        self.assertIn('unknown_metric_dropped', warning_types)
        self.assertIn('unknown_dimension_key_dropped', warning_types)

        # All warnings should be objects
        for warning in warnings:
            self.assertIsInstance(warning, dict)
            self.assertIn('type', warning)

        # Metadata should contain value coercion strategy
        self.assertIn('value_coercion_method', result['metadata'])
        self.assertEqual(result['metadata']['value_coercion_method'], 'trunc')

    def test_no_warnings_for_clean_data(self):
        """Test that clean data produces no warnings."""
        by_rule = {'R1': {'isotope_unavailable': 10}}
        by_halogen = {'F': {'isotope_unavailable': 8}}

        warnings = sanitize_granular_fields_with_warnings(
            by_rule, by_halogen, None, self.distributable_set, self.rules, self.halogens
        )

        # Should have no warnings for clean data
        self.assertEqual(len(warnings), 0)


if __name__ == '__main__':
    unittest.main()