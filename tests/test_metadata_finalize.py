"""
Tests for metadata finalization centralization.

This module tests the finalize_metadata function and verifies that
all code paths in write_qa_summary_json produce consistent metadata.
"""
import unittest
import tempfile
import os
import json
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from halogenator.report import write_qa_summary_json, finalize_metadata


class TestMetadataFinalize(unittest.TestCase):
    """Test centralized metadata finalization."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_finalize_metadata_basic(self):
        """Test basic finalize_metadata functionality."""
        metadata = {}
        warnings_list = [
            {"type": "test_warning", "message": "test"},
            {"type": "another_warning", "value": 42}
        ]

        result_metadata, result_warnings = finalize_metadata(
            metadata,
            conflict_tolerance=2,
            coercion_method="trunc",
            warnings_list=warnings_list,
            max_warnings=10
        )

        # Check that standard fields were set
        self.assertEqual(result_metadata['marginal_conflict_tolerance'], 2)
        self.assertEqual(result_metadata['value_coercion_method'], 'trunc')

        # Check warnings processing
        self.assertEqual(result_metadata['warnings_count'], 2)
        self.assertEqual(result_metadata['warnings_returned'], 2)
        self.assertEqual(result_metadata['warnings_truncated'], False)
        self.assertEqual(len(result_metadata['warnings']), 2)

        # Check that warnings are sorted
        self.assertEqual(result_metadata['warnings'][0]['type'], 'another_warning')
        self.assertEqual(result_metadata['warnings'][1]['type'], 'test_warning')

    def test_finalize_metadata_truncation(self):
        """Test warnings truncation behavior."""
        metadata = {}
        warnings_list = [
            {"type": f"warning_{i}", "index": i} for i in range(5)
        ]

        result_metadata, result_warnings = finalize_metadata(
            metadata,
            conflict_tolerance=1,
            coercion_method="trunc",
            warnings_list=warnings_list,
            max_warnings=3
        )

        # Check truncation metadata
        self.assertEqual(result_metadata['warnings_count'], 5)
        self.assertEqual(result_metadata['warnings_returned'], 3)
        self.assertEqual(result_metadata['warnings_truncated'], True)
        self.assertEqual(len(result_metadata['warnings']), 3)

        # Check deterministic sorting (first 3 after sorting)
        warning_types = [w['type'] for w in result_metadata['warnings']]
        self.assertEqual(warning_types, ['warning_0', 'warning_1', 'warning_2'])

    def test_finalize_metadata_empty_warnings(self):
        """Test behavior with no warnings."""
        metadata = {'existing_field': 'value'}
        warnings_list = []

        result_metadata, result_warnings = finalize_metadata(
            metadata,
            conflict_tolerance=3,
            coercion_method="trunc",
            warnings_list=warnings_list,
            max_warnings=100
        )

        # Check that existing fields are preserved
        self.assertEqual(result_metadata['existing_field'], 'value')

        # Check that standard fields were set
        self.assertEqual(result_metadata['marginal_conflict_tolerance'], 3)
        self.assertEqual(result_metadata['value_coercion_method'], 'trunc')

        # Check empty warnings handling
        self.assertEqual(result_metadata['warnings_count'], 0)
        self.assertEqual(result_metadata['warnings_returned'], 0)
        self.assertEqual(result_metadata['warnings_truncated'], False)
        self.assertEqual(result_metadata['warnings'], [])

    def test_write_qa_summary_v2_metadata_consistency(self):
        """Test that v2 format produces consistent metadata."""
        qa_stats = {
            'version': '2',
            'total': {'attempts': 100, 'isotope_unavailable': 10},
            'pivots': {
                'by_rule': {'R1': {'isotope_unavailable': 10}},
                'by_halogen': {'F': {'isotope_unavailable': 10}}
            }
        }

        qa_path = write_qa_summary_json(
            qa_stats, self.temp_dir,
            completion_mode='zero_fill',
            conflict_tolerance=2,
            max_warnings=500
        )

        # Verify file was created and read it
        self.assertTrue(os.path.exists(qa_path))
        with open(qa_path, 'r') as f:
            result = json.load(f)

        # Check that centralized metadata fields are present and correct
        metadata = result['metadata']
        self.assertEqual(metadata['marginal_conflict_tolerance'], 2)
        self.assertEqual(metadata['value_coercion_method'], 'trunc')
        self.assertIn('warnings_count', metadata)
        self.assertIn('warnings_returned', metadata)
        self.assertIn('warnings_truncated', metadata)

    def test_write_qa_summary_v1_metadata_consistency(self):
        """Test that v1 format produces consistent metadata."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 5.7}},  # Non-integer to trigger warning
            'by_halogen': {'F': {'isotope_unavailable': 5}},
            'metadata': {
                'rules': ['R1'],
                'halogens': ['F']
            }
        }

        qa_path = write_qa_summary_json(
            qa_stats, self.temp_dir,
            completion_mode='distribute',
            conflict_tolerance=1,
            max_warnings=100
        )

        with open(qa_path, 'r') as f:
            result = json.load(f)

        # Check centralized metadata
        metadata = result['metadata']
        self.assertEqual(metadata['marginal_conflict_tolerance'], 1)
        self.assertEqual(metadata['value_coercion_method'], 'trunc')

        # Should have warnings from non-integer value
        self.assertGreater(metadata['warnings_count'], 0)
        self.assertEqual(metadata['warnings_returned'], metadata['warnings_count'])
        self.assertEqual(metadata['warnings_truncated'], False)

    def test_write_qa_summary_legacy_metadata_consistency(self):
        """Test that legacy format produces consistent metadata."""
        qa_stats = {
            'attempts': 50,
            'isotope_unavailable': 5,
            'metadata': {
                'halogens': ['F', 'Cl'],
                'rules': ['R1', 'R2']
            }
        }

        qa_path = write_qa_summary_json(
            qa_stats, self.temp_dir,
            completion_mode='zero_fill',
            conflict_tolerance=3,
            max_warnings=200
        )

        with open(qa_path, 'r') as f:
            result = json.load(f)

        # Check centralized metadata for legacy path
        metadata = result['metadata']
        self.assertEqual(metadata['marginal_conflict_tolerance'], 3)
        self.assertEqual(metadata['value_coercion_method'], 'trunc')
        self.assertEqual(metadata['warnings_count'], 0)  # No warnings expected
        self.assertEqual(metadata['warnings_returned'], 0)
        self.assertEqual(metadata['warnings_truncated'], False)

    def test_metadata_fields_identical_across_paths(self):
        """Test that all code paths produce identical standard metadata field sets."""
        test_cases = [
            # v2 format
            {
                'version': '2',
                'total': {'attempts': 10}
            },
            # v1 with granular
            {
                'version': '1',
                'by_rule': {'R1': {'isotope_unavailable': 1}},
                'metadata': {'rules': ['R1'], 'halogens': ['F']}
            },
            # Legacy with halogens
            {
                'attempts': 10,
                'metadata': {'halogens': ['F'], 'rules': ['R1']}
            },
            # Legacy totals-only
            {
                'attempts': 10
            }
        ]

        required_fields = {
            'marginal_conflict_tolerance',
            'value_coercion_method',
            'warnings_count',
            'warnings_returned',
            'warnings_truncated'
        }

        for i, qa_stats in enumerate(test_cases):
            with self.subTest(case=i):
                qa_path = write_qa_summary_json(
                    qa_stats, self.temp_dir,
                    completion_mode='zero_fill',
                    conflict_tolerance=1,
                    max_warnings=50
                )

                with open(qa_path, 'r') as f:
                    result = json.load(f)

                metadata = result['metadata']

                # Check that all required fields are present
                for field in required_fields:
                    self.assertIn(field, metadata,
                                 f"Field {field} missing in case {i}")

                # Check field values are correct
                self.assertEqual(metadata['marginal_conflict_tolerance'], 1)
                self.assertEqual(metadata['value_coercion_method'], 'trunc')


if __name__ == '__main__':
    unittest.main()