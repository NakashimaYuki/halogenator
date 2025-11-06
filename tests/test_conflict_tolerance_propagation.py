"""
Tests for conflict tolerance propagation through all write_qa_summary_json call paths.

This module ensures that the conflict_tolerance parameter is properly propagated
and recorded in metadata across all invocation scenarios.
"""
import unittest
import tempfile
import os
import json
import shutil
from src.halogenator.report import write_qa_summary_json


class TestConflictTolerancePropagation(unittest.TestCase):
    """Test conflict tolerance parameter propagation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_default_tolerance_recorded_in_metadata(self):
        """Test that default conflict tolerance is recorded in metadata."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 5}},
            'by_halogen': {'F': {'isotope_unavailable': 5}},
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        # Call without explicit conflict_tolerance (should use default)
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have default tolerance value in metadata
        self.assertIn('marginal_conflict_tolerance', result['metadata'])
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 1)

    def test_explicit_tolerance_recorded_in_metadata(self):
        """Test that explicit conflict tolerance is recorded in metadata."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 12}},  # Conflict > tolerance
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        # Call with explicit conflict_tolerance
        custom_tolerance = 5
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir,
                                                conflict_tolerance=custom_tolerance)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should have custom tolerance value in metadata
        self.assertIn('marginal_conflict_tolerance', result['metadata'])
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], custom_tolerance)

    def test_tolerance_propagation_with_completion_modes(self):
        """Test tolerance propagation works with all completion modes."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 8}},
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        # Test both completion modes
        for mode in ['zero_fill', 'distribute']:
            with self.subTest(completion_mode=mode):
                custom_tolerance = 3
                qa_summary_path = write_qa_summary_json(
                    qa_stats, self.temp_dir,
                    completion_mode=mode,
                    conflict_tolerance=custom_tolerance
                )

                with open(qa_summary_path, 'r', encoding='utf-8') as f:
                    result = json.load(f)

                self.assertIn('marginal_conflict_tolerance', result['metadata'])
                self.assertEqual(result['metadata']['marginal_conflict_tolerance'], custom_tolerance)

    def test_tolerance_propagation_with_v2_format(self):
        """Test tolerance propagation with v2 format input."""
        qa_stats = {
            'version': '2',
            'total': {'isotope_unavailable': 15},
            'pivots': {
                'by_rule': {'R1': {'isotope_unavailable': 15}},
                'by_halogen': {'F': {'isotope_unavailable': 15}}
            }
        }

        custom_tolerance = 2
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir,
                                                conflict_tolerance=custom_tolerance)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # V2 format should also record tolerance in metadata
        self.assertIn('marginal_conflict_tolerance', result['metadata'])
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], custom_tolerance)

    def test_tolerance_propagation_with_legacy_format(self):
        """Test tolerance propagation with legacy totals-only format."""
        qa_stats = {
            'isotope_unavailable': 10,
            'no_product_matches': 2
        }

        custom_tolerance = 4
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir,
                                                conflict_tolerance=custom_tolerance)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Legacy format should also record tolerance in metadata
        self.assertIn('marginal_conflict_tolerance', result['metadata'])
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], custom_tolerance)

    def test_non_cli_call_path_coverage(self):
        """Test that non-CLI call paths properly handle tolerance parameter."""
        # Simulate a direct library usage (non-CLI call)
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 20}, 'R2': {'isotope_unavailable': 15}},
            'by_halogen': {'F': {'isotope_unavailable': 18}, 'Cl': {'isotope_unavailable': 17}},
            'metadata': {'rules': ['R1', 'R2'], 'halogens': ['F', 'Cl']}
        }

        # This represents how the function might be called from other modules
        tolerance = 7
        qa_summary_path = write_qa_summary_json(
            qa_stats_dict=qa_stats,
            output_dir=self.temp_dir,
            completion_mode='distribute',
            conflict_tolerance=tolerance
        )

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should record the tolerance in metadata
        self.assertIn('marginal_conflict_tolerance', result['metadata'])
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], tolerance)

        # Should also detect conflicts since marginal totals don't match exactly
        # R1+R2 = 35, F+Cl = 35 (no conflict), but individual cases may trigger conflicts
        # based on distribution logic
        if 'warnings' in result['metadata']:
            conflict_warnings = [w for w in result['metadata']['warnings']
                               if w.get('type') == 'marginal_conflict']
            for warning in conflict_warnings:
                # Conflict warnings should reference the tolerance used
                self.assertEqual(warning.get('tolerance'), tolerance)


if __name__ == '__main__':
    unittest.main()