"""
Tests for warnings truncation and throttling in QA summary generation.

This module tests that warnings are properly truncated when they exceed
the max_warnings threshold and that appropriate metadata is recorded.
"""
import unittest
import tempfile
import os
import json
import shutil
from halogenator.report import write_qa_summary_json


class TestWarningsTruncation(unittest.TestCase):
    """Test warnings truncation and metadata recording."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _create_deterministic_warnings_qa_stats(self, num_warnings):
        """Create QA stats with exact number of deterministic warnings."""
        # Each rule/dimension combo with unknown metrics creates exactly one warning
        # So to get num_warnings, we need num_warnings different rule keys
        by_rule = {}
        rules = []

        for i in range(num_warnings):
            rule_key = f'R{i:02d}'  # R00, R01, R02, etc.
            rules.append(rule_key)
            by_rule[rule_key] = {
                'isotope_unavailable': 5,  # Valid metric
                f'bad_metric_{i:03d}': 10  # Unknown metric (creates 1 warning per rule)
            }

        return {
            'version': '1',
            'by_rule': by_rule,
            'metadata': {'rules': rules, 'halogens': ['F']}
        }

    def test_no_truncation_below_threshold(self):
        """Test that warnings below threshold are not truncated."""
        qa_stats = self._create_deterministic_warnings_qa_stats(5)

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=10)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        metadata = result['metadata']

        # Should not be truncated
        self.assertEqual(metadata['warnings_truncated'], False)
        self.assertEqual(metadata['warnings_count'], 5)
        self.assertEqual(metadata['warnings_returned'], 5)
        self.assertEqual(len(metadata['warnings']), 5)

    def test_truncation_above_threshold(self):
        """Test that warnings are truncated when exceeding threshold."""
        qa_stats = self._create_deterministic_warnings_qa_stats(15)

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=10)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        metadata = result['metadata']

        # Should be truncated
        self.assertEqual(metadata['warnings_truncated'], True)
        self.assertEqual(metadata['warnings_count'], 15)  # Original count
        self.assertEqual(metadata['warnings_returned'], 10)  # Truncated count
        self.assertEqual(len(metadata['warnings']), 10)  # Actual warnings in output

    def test_truncation_preserves_first_warnings_deterministically(self):
        """Test that truncation preserves the first N warnings in stable sorted order."""
        qa_stats = self._create_deterministic_warnings_qa_stats(8)

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=3)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        warnings = result['metadata']['warnings']
        self.assertEqual(len(warnings), 3)

        # Warnings should be in stable sorted order (type, dimension, key, metric)
        # All warnings are "unknown_metric_dropped" for "by_rule" dimension
        # Keys should be sorted: R00, R01, R02
        expected_keys = ['R00', 'R01', 'R02']

        for i, warning in enumerate(warnings):
            self.assertEqual(warning['type'], 'unknown_metric_dropped')
            self.assertEqual(warning['dimension'], 'by_rule')
            self.assertEqual(warning['key'], expected_keys[i])
            # Each warning should contain the corresponding bad metric
            expected_metric = f'bad_metric_{i:03d}'
            self.assertIn(expected_metric, warning['metrics'])

    def test_custom_max_warnings_parameter(self):
        """Test different max_warnings values with deterministic data."""
        test_cases = [
            {'max_warnings': 5, 'expected_returned': 5, 'expected_truncated': True},
            {'max_warnings': 15, 'expected_returned': 15, 'expected_truncated': True},
            {'max_warnings': 20, 'expected_returned': 20, 'expected_truncated': False},  # Exact match
            {'max_warnings': 25, 'expected_returned': 20, 'expected_truncated': False},  # More than available
        ]

        for case in test_cases:
            with self.subTest(max_warnings=case['max_warnings']):
                # Create fresh data for each subtest to avoid mutation issues
                qa_stats = self._create_deterministic_warnings_qa_stats(20)

                qa_summary_path = write_qa_summary_json(
                    qa_stats, self.temp_dir, max_warnings=case['max_warnings']
                )

                with open(qa_summary_path, 'r', encoding='utf-8') as f:
                    result = json.load(f)

                metadata = result['metadata']

                # Precise deterministic assertions
                self.assertEqual(metadata['warnings_count'], 20)  # Always 20 original warnings
                self.assertEqual(metadata['warnings_returned'], case['expected_returned'])
                self.assertEqual(metadata['warnings_truncated'], case['expected_truncated'])
                self.assertEqual(len(metadata['warnings']), case['expected_returned'])

    def test_zero_max_warnings(self):
        """Test that max_warnings=0 removes all warnings."""
        qa_stats = self._create_deterministic_warnings_qa_stats(5)

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=0)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        metadata = result['metadata']

        # Should be truncated to 0
        self.assertEqual(metadata['warnings_truncated'], True)
        self.assertEqual(metadata['warnings_count'], 5)  # Original count
        self.assertEqual(metadata['warnings_returned'], 0)  # Truncated to 0
        self.assertEqual(len(metadata['warnings']), 0)

    def test_no_warnings_case(self):
        """Test behavior when there are no warnings to begin with."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 5}},  # Clean data, no warnings
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=10)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        metadata = result['metadata']

        # With centralized finalize_metadata, these fields are always set
        self.assertEqual(metadata['warnings_truncated'], False)
        self.assertEqual(metadata['warnings_count'], 0)
        self.assertEqual(metadata['warnings_returned'], 0)
        self.assertEqual(metadata['warnings'], [])

    def test_mixed_warning_types_truncation_deterministic(self):
        """Test truncation with mixed warning types using deterministic data."""
        # Create exactly 4 warnings: 2 unknown metrics + 1 marginal conflict
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {
                    'isotope_unavailable': 10,
                    'bad_metric_a': 1,      # Unknown metric -> warning
                    'bad_metric_b': 2       # Unknown metric -> warning
                }
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 15}  # Conflict: 10 != 15 (diff = 5 > tolerance 1)
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(
            qa_stats, self.temp_dir,
            max_warnings=1,  # Truncate to test truncation behavior
            conflict_tolerance=1  # Ensure conflict is detected
        )

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        metadata = result['metadata']

        # Should have exactly 2 warnings total, truncated to 1
        self.assertEqual(metadata['warnings_count'], 2)
        self.assertEqual(metadata['warnings_returned'], 1)
        self.assertEqual(metadata['warnings_truncated'], True)
        self.assertEqual(len(metadata['warnings']), 1)

        # Verify warnings are sorted deterministically and truncated correctly
        warnings = metadata['warnings']
        # After sorting by (type, dimension, key, metric), the order should be predictable:
        # marginal_conflict has ("marginal_conflict", "", "", "") - empty strings sort first
        # unknown_metric_dropped has ("unknown_metric_dropped", "by_rule", "R1", "bad_metric_a")
        # So marginal_conflict should always be first
        self.assertEqual(warnings[0]['type'], 'marginal_conflict')

    def test_warnings_stable_sorting(self):
        """Test that warnings are sorted in a stable, predictable order."""
        qa_stats = self._create_deterministic_warnings_qa_stats(6)

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, max_warnings=100)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        warnings = result['metadata']['warnings']
        self.assertEqual(len(warnings), 6)

        # All warnings should be "unknown_metric_dropped" type, same dimension
        # They should be sorted by rule key: R00, R01, R02, R03, R04, R05
        expected_keys = ['R00', 'R01', 'R02', 'R03', 'R04', 'R05']

        for i, warning in enumerate(warnings):
            self.assertEqual(warning['type'], 'unknown_metric_dropped')
            self.assertEqual(warning['dimension'], 'by_rule')
            self.assertEqual(warning['key'], expected_keys[i])
            # Each warning should contain the corresponding bad metric
            expected_metric = f'bad_metric_{i:03d}'
            self.assertIn(expected_metric, warning['metrics'])


if __name__ == '__main__':
    unittest.main()