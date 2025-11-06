"""
Tests for conflict tolerance CLI parameter and metadata recording.

This module tests the --qa-conflict-tolerance CLI parameter and ensures
the tolerance value is properly recorded in metadata.
"""
import unittest
import tempfile
import os
import json
from src.halogenator.report import write_qa_summary_json


class TestMarginalConflictToleranceCLI(unittest.TestCase):
    """Test conflict tolerance parameter handling and metadata recording."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_default_conflict_tolerance_is_one(self):
        """Test that default conflict tolerance is 1."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 12}},  # diff=2 > 1, should trigger
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should record default tolerance in metadata
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 1)

        # Should have conflict warning because diff=2 > tolerance=1
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 1)
        self.assertEqual(conflict_warnings[0]['tolerance'], 1)

    def test_custom_conflict_tolerance_zero(self):
        """Test that tolerance=0 triggers warnings for any difference."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 11}},  # diff=1 > 0, should trigger
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, conflict_tolerance=0)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should record custom tolerance in metadata
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 0)

        # Should have conflict warning because diff=1 > tolerance=0
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 1)
        self.assertEqual(conflict_warnings[0]['tolerance'], 0)

    def test_custom_conflict_tolerance_two(self):
        """Test that tolerance=2 prevents warnings for smaller differences."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 12}},  # diff=2 == 2, should NOT trigger
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, conflict_tolerance=2)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should record custom tolerance in metadata
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 2)

        # Should NOT have conflict warning because diff=2 is not > tolerance=2
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 0)

    def test_conflict_tolerance_boundary_plus_one(self):
        """Test that tolerance+1 triggers warning at exact boundary."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'by_halogen': {'F': {'isotope_unavailable': 13}},  # diff=3 > tolerance=2
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, conflict_tolerance=2)

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should record tolerance in metadata
        self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 2)

        # Should have conflict warning because diff=3 > tolerance=2
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 1)
        self.assertEqual(conflict_warnings[0]['tolerance'], 2)

        # Check detailed warning payload
        delta = conflict_warnings[0]['delta']
        self.assertIn('isotope_unavailable', delta)
        self.assertEqual(delta['isotope_unavailable']['lhs'], 10)
        self.assertEqual(delta['isotope_unavailable']['rhs'], 13)
        self.assertEqual(delta['isotope_unavailable']['diff'], 3)

    def test_tolerance_metadata_consistency_across_modes(self):
        """Test that tolerance is recorded in metadata for both completion modes."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        for mode in ['zero_fill', 'distribute']:
            with self.subTest(mode=mode):
                qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode=mode, conflict_tolerance=5)

                with open(qa_summary_path, 'r', encoding='utf-8') as f:
                    result = json.load(f)

                # Should record tolerance regardless of completion mode
                self.assertEqual(result['metadata']['marginal_conflict_tolerance'], 5)
                self.assertEqual(result['metadata']['completion']['strategy'], mode)

    def test_function_signature_compatibility(self):
        """Test that function can be called with positional and keyword arguments."""
        qa_stats = {
            'version': '1',
            'by_rule': {'R1': {'isotope_unavailable': 10}},
            'metadata': {'rules': ['R1'], 'halogens': ['F']}
        }

        # Test positional arguments
        try:
            qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, 'zero_fill', 3)
            self.assertTrue(os.path.exists(qa_summary_path))
        except Exception as e:
            self.fail(f"Positional arguments should work: {e}")

        # Test keyword arguments
        try:
            qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute', conflict_tolerance=4)
            self.assertTrue(os.path.exists(qa_summary_path))
        except Exception as e:
            self.fail(f"Keyword arguments should work: {e}")

        # Test mixed arguments
        try:
            qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, 'zero_fill', conflict_tolerance=5)
            self.assertTrue(os.path.exists(qa_summary_path))
        except Exception as e:
            self.fail(f"Mixed arguments should work: {e}")


if __name__ == '__main__':
    unittest.main()