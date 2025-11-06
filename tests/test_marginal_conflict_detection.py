# -*- coding: ascii -*-
"""Test marginal conflict detection in passthrough pipeline."""

import unittest
import tempfile
import os
import json
from src.halogenator.report import write_qa_summary_json


class TestMarginalConflictDetection(unittest.TestCase):
    """Test marginal conflict detection and warning generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_marginal_conflict_detects_metrics_present_in_only_one_dimension(self):
        """Test that conflict detection catches metrics in only one dimension."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10, 'no_product_matches': 5}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 12, 'template_unsupported': 3}
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='zero_fill')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have warnings
        self.assertIn('warnings', result['metadata'])
        warnings = result['metadata']['warnings']
        self.assertTrue(len(warnings) > 0)

        # Find marginal conflict warning
        conflict_warning = None
        for warning in warnings:
            if warning.get('type') == 'marginal_conflict':
                conflict_warning = warning
                break

        self.assertIsNotNone(conflict_warning, "Should have marginal conflict warning")

        # Check that conflicts include metrics from both dimensions
        delta = conflict_warning['delta']

        # isotope_unavailable has different values (10 vs 12, diff=2 > tolerance=1)
        self.assertIn('isotope_unavailable', delta)

        # no_product_matches has values (5 vs 0, diff=5 > tolerance=1)
        self.assertIn('no_product_matches', delta)

        # template_unsupported has values (0 vs 3, diff=3 > tolerance=1)
        self.assertIn('template_unsupported', delta)

    def test_marginal_conflict_payload_contains_lhs_rhs_diff_and_tolerance(self):
        """Test that conflict warning contains detailed payload."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 13}  # diff=3 > tolerance=1
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='zero_fill')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Find marginal conflict warning
        conflict_warning = None
        for warning in result['metadata']['warnings']:
            if warning.get('type') == 'marginal_conflict':
                conflict_warning = warning
                break

        self.assertIsNotNone(conflict_warning)

        # Check warning structure
        self.assertEqual(conflict_warning['type'], 'marginal_conflict')
        self.assertEqual(conflict_warning['base'], 'by_rule')
        self.assertEqual(conflict_warning['tolerance'], 1)
        self.assertIn('delta', conflict_warning)

        # Check delta structure for isotope_unavailable
        delta = conflict_warning['delta']
        self.assertIn('isotope_unavailable', delta)

        metric_conflict = delta['isotope_unavailable']
        self.assertEqual(metric_conflict['lhs'], 10)  # from by_rule
        self.assertEqual(metric_conflict['rhs'], 13)  # from by_halogen
        self.assertEqual(metric_conflict['diff'], 3)  # abs(10-13)

    def test_marginal_conflict_only_checks_distributable_metrics(self):
        """Test that conflict detection only checks distributable metrics."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {
                    'isotope_unavailable': 10,  # distributable
                    'attempts': 100  # overview counter - should be ignored
                }
            },
            'by_halogen': {
                'F': {
                    'isotope_unavailable': 10,  # same value, no conflict
                    'attempts': 200  # overview counter with different value - should be ignored
                }
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='zero_fill')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should not have marginal conflict warnings
        # because isotope_unavailable values are the same (10 vs 10)
        # and attempts is not a distributable metric
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 0)

    def test_marginal_conflict_respects_tolerance_threshold(self):
        """Test that conflict detection respects tolerance threshold."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 11}  # diff=1, equals tolerance
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='zero_fill')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should not have conflict warning because diff=1 is not > tolerance=1
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 0)

    def test_marginal_conflict_triggers_on_tolerance_plus_one(self):
        """Test that conflict detection triggers when diff > tolerance."""
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 12}  # diff=2 > tolerance=1
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='zero_fill')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have conflict warning because diff=2 > tolerance=1
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 1)


if __name__ == '__main__':
    unittest.main()