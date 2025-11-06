"""
Integration test for base priority and conflict warning interaction.

This module tests the integration between base selection logic and
conflict detection when both by_rule and by_halogen marginals exist.
"""
import unittest
import tempfile
import os
import json
from src.halogenator.report import write_qa_summary_json


class TestIntegrationBaseAndConflict(unittest.TestCase):
    """Test integration of base priority and conflict detection."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_base_priority_with_conflict_detection(self):
        """
        Test that by_rule takes priority for distribution while conflict is detected.

        This test verifies:
        1. 2D distribution uses by_rule as base (ignores by_halogen totals)
        2. Conflict detection still reports differences between marginals
        3. metadata.completion.base = 'by_rule'
        4. Conflict warning contains detailed delta information
        """
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10, 'no_product_matches': 2},
                'R2': {'isotope_unavailable': 8, 'no_product_matches': 1}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 12, 'no_product_matches': 2},  # Different totals
                'Cl': {'isotope_unavailable': 4, 'no_product_matches': 3}   # Different totals
            },
            'metadata': {'rules': ['R1', 'R2'], 'halogens': ['F', 'Cl']}
        }

        # Use distribute mode to trigger base priority logic
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute')

        # Read the output
        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # 1. Verify base priority: 2D totals match by_rule
        by_rule_halogen = result['by_rule_halogen']
        actual_isotope_total = sum(
            by_rule_halogen[r][h]['isotope_unavailable'] for r in ['R1', 'R2'] for h in ['F', 'Cl']
        )
        actual_no_product_total = sum(
            by_rule_halogen[r][h]['no_product_matches'] for r in ['R1', 'R2'] for h in ['F', 'Cl']
        )

        expected_by_rule_isotope = 10 + 8  # 18
        expected_by_rule_no_product = 2 + 1   # 3

        # Should match by_rule totals (not by_halogen)
        self.assertEqual(actual_isotope_total, expected_by_rule_isotope)
        self.assertEqual(actual_no_product_total, expected_by_rule_no_product)

        # Verify by_halogen totals are different (confirming conflict existed)
        by_halogen_isotope = 12 + 4  # 16 != 18
        by_halogen_no_product = 2 + 3    # 5 != 3
        self.assertNotEqual(by_halogen_isotope, expected_by_rule_isotope)
        self.assertNotEqual(by_halogen_no_product, expected_by_rule_no_product)

        # 2. Verify metadata.completion.base is by_rule
        self.assertIn('completion', result['metadata'])
        self.assertEqual(result['metadata']['completion']['base'], 'by_rule')
        self.assertEqual(result['metadata']['completion']['strategy'], 'distribute')

        # 3. Verify conflict warning is present with correct details
        self.assertIn('warnings', result['metadata'])
        warnings = result['metadata']['warnings']

        # Find marginal conflict warning
        conflict_warning = None
        for warning in warnings:
            if warning.get('type') == 'marginal_conflict':
                conflict_warning = warning
                break

        self.assertIsNotNone(conflict_warning, "Should have marginal conflict warning")

        # Check warning structure
        self.assertEqual(conflict_warning['type'], 'marginal_conflict')
        self.assertEqual(conflict_warning['base'], 'by_rule')
        self.assertEqual(conflict_warning['tolerance'], 1)
        self.assertIn('delta', conflict_warning)

        # Check delta contains both conflicting metrics
        delta = conflict_warning['delta']

        # isotope_unavailable conflict: by_rule=18, by_halogen=16, diff=2 > tolerance=1
        self.assertIn('isotope_unavailable', delta)
        isotope_conflict = delta['isotope_unavailable']
        self.assertEqual(isotope_conflict['lhs'], 18)  # by_rule total
        self.assertEqual(isotope_conflict['rhs'], 16)  # by_halogen total
        self.assertEqual(isotope_conflict['diff'], 2)

        # no_product_matches conflict: by_rule=3, by_halogen=5, diff=2 > tolerance=1
        self.assertIn('no_product_matches', delta)
        no_product_conflict = delta['no_product_matches']
        self.assertEqual(no_product_conflict['lhs'], 3)   # by_rule total
        self.assertEqual(no_product_conflict['rhs'], 5)   # by_halogen total
        self.assertEqual(no_product_conflict['diff'], 2)

    def test_by_rule_priority_invariant_under_marginal_swap(self):
        """
        Test that by_rule priority is maintained even when by_halogen has 'better' totals.
        """
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 5}  # Lower total
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 10}, # Higher total
                'Cl': {'isotope_unavailable': 8}  # Higher total
            },
            'metadata': {'rules': ['R1'], 'halogens': ['F', 'Cl']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute')

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # 2D should still match by_rule (5), not by_halogen (18)
        by_rule_halogen = result['by_rule_halogen']
        actual_total = sum(
            by_rule_halogen['R1'][h]['isotope_unavailable'] for h in ['F', 'Cl']
        )

        self.assertEqual(actual_total, 5)  # by_rule total
        self.assertNotEqual(actual_total, 18) # by_halogen total

        # Base should still be by_rule
        self.assertEqual(result['metadata']['completion']['base'], 'by_rule')

    def test_single_marginal_no_conflict_warning(self):
        """
        Test that when only one marginal exists, no conflict warning is generated.
        """
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 10}
            },
            # No by_halogen
            'metadata': {'rules': ['R1'], 'halogens': ['F', 'Cl']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute')

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Should not have marginal conflict warning
        warnings = result['metadata'].get('warnings', [])
        conflict_warnings = [w for w in warnings if w.get('type') == 'marginal_conflict']
        self.assertEqual(len(conflict_warnings), 0)

        # Base should be by_rule
        self.assertEqual(result['metadata']['completion']['base'], 'by_rule')

        # 2D should be distributed from by_rule
        by_rule_halogen = result['by_rule_halogen']
        actual_total = sum(
            by_rule_halogen['R1'][h]['isotope_unavailable'] for h in ['F', 'Cl']
        )
        self.assertEqual(actual_total, 10)

    def test_marginal_consistency_after_base_priority_distribution(self):
        """
        Test that marginals derived from 2D are consistent after base priority distribution.
        """
        qa_stats = {
            'version': '1',
            'by_rule': {
                'R1': {'isotope_unavailable': 7},  # Odd number for remainder test
                'R2': {'isotope_unavailable': 5}
            },
            'by_halogen': {
                'F': {'isotope_unavailable': 20},  # Different total (ignored)
                'Cl': {'isotope_unavailable': 30}  # Different total (ignored)
            },
            'metadata': {'rules': ['R1', 'R2'], 'halogens': ['F', 'Cl']}
        }

        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='distribute')

        with open(qa_summary_path, 'r', encoding='utf-8') as f:
            result = json.load(f)

        # Check marginal consistency
        by_rule_derived = result['by_rule']
        by_halogen_derived = result['by_halogen']
        by_rule_halogen = result['by_rule_halogen']

        # by_rule marginals should sum to original by_rule totals
        self.assertEqual(by_rule_derived['R1']['isotope_unavailable'], 7)
        self.assertEqual(by_rule_derived['R2']['isotope_unavailable'], 5)

        # by_halogen marginals should be consistent with 2D
        expected_f_total = by_rule_halogen['R1']['F']['isotope_unavailable'] + by_rule_halogen['R2']['F']['isotope_unavailable']
        expected_cl_total = by_rule_halogen['R1']['Cl']['isotope_unavailable'] + by_rule_halogen['R2']['Cl']['isotope_unavailable']

        self.assertEqual(by_halogen_derived['F']['isotope_unavailable'], expected_f_total)
        self.assertEqual(by_halogen_derived['Cl']['isotope_unavailable'], expected_cl_total)

        # Total should equal original by_rule total (12), not by_halogen total (50)
        total_from_2d = sum(
            by_rule_halogen[r][h]['isotope_unavailable']
            for r in ['R1', 'R2'] for h in ['F', 'Cl']
        )
        self.assertEqual(total_from_2d, 12)  # 7 + 5 from by_rule
        self.assertNotEqual(total_from_2d, 50)  # 20 + 30 from by_halogen


if __name__ == '__main__':
    unittest.main()