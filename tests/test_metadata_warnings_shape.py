# -*- coding: ascii -*-
"""Test metadata warnings structure consistency."""

import unittest
import tempfile
import os
import json
from src.halogenator.report import write_qa_summary_json, distribute_marginals_to_2d_with_warnings


class TestMetadataWarningsShape(unittest.TestCase):
    """Test that metadata warnings have consistent object structure."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_metadata_warnings_are_objects_with_type_field(self):
        """Test that all warnings are objects with type field."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10}  # Will conflict with by_halogen
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 15}   # Different value
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have warnings
        self.assertIn('warnings', result['metadata'])
        warnings = result['metadata']['warnings']
        self.assertIsInstance(warnings, list)

        # All warnings should be objects with type field
        for warning in warnings:
            self.assertIsInstance(warning, dict)
            self.assertIn('type', warning)
            self.assertIsInstance(warning['type'], str)

            # Check specific warning types have expected structure
            if warning['type'] == 'marginal_conflict':
                self.assertIn('base', warning)
                self.assertIn('tolerance', warning)
                self.assertIn('delta', warning)
            elif warning['type'] == 'empty_rules_or_halogens':
                # Should be minimal structure
                self.assertTrue(True)  # Just checking it has 'type'

    def test_distribute_with_warnings_returns_object_warnings(self):
        """Test that distribute_marginals_to_2d_with_warnings returns object warnings."""
        by_rule = {'R1': {'metric1': 10}}
        by_halogen = None
        rules = []  # Empty to trigger warning
        halogens = ['F']
        metrics = ['metric1']

        result, warnings = distribute_marginals_to_2d_with_warnings(
            by_rule, by_halogen, rules, halogens, metrics
        )

        # Should have warnings as objects
        self.assertIsInstance(warnings, list)
        self.assertEqual(len(warnings), 1)

        warning = warnings[0]
        self.assertIsInstance(warning, dict)
        self.assertEqual(warning['type'], 'empty_rules_or_halogens')

    def test_warnings_structure_unified_across_scenarios(self):
        """Test that warnings structure is consistent across different scenarios."""
        # Test both marginal conflict and empty set scenarios
        scenarios = [
            # Empty set scenario
            {
                'name': 'empty_set',
                'qa_stats': {
                    'version': '1',
                    'by_rule': {'R1': {'metric1': 5}},
                    'by_halogen': {'F': {'metric1': 5}},
                    'metadata': {'rules': [], 'halogens': ['F']}  # Empty rules
                },
                'mode': 'distribute',
                'expected_types': ['empty_rules_or_halogens']
            },
            # Marginal conflict scenario
            {
                'name': 'conflict',
                'qa_stats': {
                    'version': '1',
                    'by_rule': {'R1': {'isotope_unavailable': 10}},
                    'by_halogen': {'F': {'isotope_unavailable': 15}},
                    'metadata': {'rules': ['R1'], 'halogens': ['F']}
                },
                'mode': 'zero_fill',
                'expected_types': ['marginal_conflict']
            }
        ]

        for scenario in scenarios:
            with self.subTest(scenario=scenario['name']):
                qa_summary_path = write_qa_summary_json(
                    scenario['qa_stats'],
                    self.temp_dir,
                    completion_mode=scenario['mode']
                )

                # Read the output
                with open(qa_summary_path, 'r', encoding='utf-8') as f:
                    result = json.load(f)

                warnings = result['metadata'].get('warnings', [])

                # Check that expected warning types are present
                warning_types = [w['type'] for w in warnings]
                for expected_type in scenario['expected_types']:
                    self.assertIn(expected_type, warning_types)

                # All warnings should be objects
                for warning in warnings:
                    self.assertIsInstance(warning, dict)
                    self.assertIn('type', warning)

    def test_warnings_list_is_top_level_metadata_field(self):
        """Test that warnings is a top-level metadata field, not nested."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 13}},  # Conflict
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # warnings should be directly under metadata
        self.assertIn('metadata', result)
        self.assertIn('warnings', result['metadata'])
        self.assertIsInstance(result['metadata']['warnings'], list)

        # warnings should coexist with other metadata fields
        self.assertIn('description', result['metadata'])
        self.assertIn('generated_at', result['metadata'])

    def test_no_warnings_produces_empty_list_or_no_field(self):
        """Test behavior when there are no warnings to report."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 10}},  # Same value, no conflict
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should either have no warnings field or empty warnings list
        warnings = result['metadata'].get('warnings', [])
        self.assertIsInstance(warnings, list)
        self.assertEqual(len(warnings), 0)


if __name__ == '__main__':
    unittest.main()