# -*- coding: ascii -*-
"""Tests for report generation with ops/atoms pivot CSVs."""

import unittest
import os
import tempfile
import shutil
import sys
from pathlib import Path
import pandas as pd

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig, QAAggregator
from halogenator.report import write_qa_summary_json


class TestReportOpsAtomsPivots(unittest.TestCase):
    """Test that report generation includes ops/atoms pivot CSV files."""

    def setUp(self):
        """Set up test environment with temporary directory."""
        self.temp_dir = tempfile.mkdtemp()
        self.cfg = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),
            rules_cfg={
                'R1': {'enable': True},
                'R3': {'enable': True},
                'R6_methyl': {
                    'enable': True,
                    'step': {'enable': True},
                    'macro': {'enable': True}
                }
            }
        )

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_pivot_csvs_generated_from_v2_data(self):
        """Test that pivot CSVs are generated when QA data has version 2 format with pivots."""
        # Generate some products with QA stats
        products, qa_stats = enumerate_with_stats('CC', self.cfg)

        # Ensure we have version 2 data with pivots
        self.assertEqual(qa_stats.get('version'), '2')
        self.assertIn('pivots', qa_stats)

        # Write QA summary to temp directory
        qa_summary_path = write_qa_summary_json(qa_stats, self.temp_dir)

        # Check that the main pivot CSV exists
        main_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_halogen_k.csv')
        self.assertTrue(os.path.exists(main_pivot_path),
                       "pivot_rule_halogen_k.csv should be generated")

        # Check that ops/atoms pivot CSVs exist
        ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_halogen_ops_atoms.csv')
        rule_ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_ops_atoms.csv')

        self.assertTrue(os.path.exists(ops_atoms_pivot_path),
                       "pivot_rule_halogen_ops_atoms.csv should be generated")
        self.assertTrue(os.path.exists(rule_ops_atoms_pivot_path),
                       "pivot_rule_ops_atoms.csv should be generated")

    def test_ops_atoms_pivot_csv_structure(self):
        """Test that ops/atoms pivot CSV has correct structure."""
        # Generate products
        products, qa_stats = enumerate_with_stats('CC', self.cfg)

        # Write QA summary
        write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the ops/atoms pivot CSV
        ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_halogen_ops_atoms.csv')

        if os.path.exists(ops_atoms_pivot_path):
            try:
                df = pd.read_csv(ops_atoms_pivot_path)

                # Check expected columns
                expected_columns = ['rule', 'halogen', 'k_ops', 'k_atoms', 'attempts', 'products',
                                  'no_product_matches', 'isotope_unavailable', 'isotope_miss',
                                  'atommap_used', 'heuristic_used', 'template_unsupported']

                for col in expected_columns:
                    self.assertIn(col, df.columns,
                                f"Column '{col}' should be present in pivot_rule_halogen_ops_atoms.csv")

                # Check data types
                if not df.empty:
                    self.assertTrue(df['rule'].dtype == 'object', "Rule column should be string")
                    self.assertTrue(df['halogen'].dtype == 'object', "Halogen column should be string")
                    # k_ops and k_atoms might be strings if they come from parsing keys
                    for metric_col in ['attempts', 'products', 'no_product_matches']:
                        if metric_col in df.columns:
                            self.assertTrue(pd.api.types.is_numeric_dtype(df[metric_col]),
                                          f"{metric_col} should be numeric")

            except Exception as e:
                self.fail(f"Failed to read/validate pivot_rule_halogen_ops_atoms.csv: {e}")

    def test_rule_ops_atoms_pivot_csv_structure(self):
        """Test that rule ops/atoms pivot CSV has correct structure."""
        # Generate products
        products, qa_stats = enumerate_with_stats('CC', self.cfg)

        # Write QA summary
        write_qa_summary_json(qa_stats, self.temp_dir)

        # Read the rule ops/atoms pivot CSV
        rule_ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_ops_atoms.csv')

        if os.path.exists(rule_ops_atoms_pivot_path):
            try:
                df = pd.read_csv(rule_ops_atoms_pivot_path)

                # Check expected columns (no halogen column in this one)
                expected_columns = ['rule', 'k_ops', 'k_atoms', 'attempts', 'products',
                                  'no_product_matches', 'isotope_unavailable', 'isotope_miss',
                                  'atommap_used', 'heuristic_used', 'template_unsupported']

                for col in expected_columns:
                    self.assertIn(col, df.columns,
                                f"Column '{col}' should be present in pivot_rule_ops_atoms.csv")

                # Check that halogen column is NOT present
                self.assertNotIn('halogen', df.columns,
                                "Halogen column should NOT be present in pivot_rule_ops_atoms.csv")

                # Check data types
                if not df.empty:
                    self.assertTrue(df['rule'].dtype == 'object', "Rule column should be string")
                    for metric_col in ['attempts', 'products', 'no_product_matches']:
                        if metric_col in df.columns:
                            self.assertTrue(pd.api.types.is_numeric_dtype(df[metric_col]),
                                          f"{metric_col} should be numeric")

            except Exception as e:
                self.fail(f"Failed to read/validate pivot_rule_ops_atoms.csv: {e}")

    def test_pivot_csvs_empty_when_no_ops_atoms_data(self):
        """Test that pivot CSVs are created (possibly empty) even when no ops/atoms data exists."""
        # Create minimal QA stats with version 2 but no ops/atoms pivot data
        minimal_qa_stats = {
            'version': '2',
            'pivots': {
                'by_rule': {},
                'by_halogen': {},
                'by_k': {},
                'by_rule_halogen': {},
                'by_rule_halogen_k': {},
                'by_rule_halogen_ops_atoms': {},  # Empty
                'by_rule_ops_atoms': {}  # Empty
            },
            'attempts': 0,
            'products': 0,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {},
            'dedup_hits_statesig': 0,
            'dedup_hits_inchi': 0
        }

        # Write QA summary
        write_qa_summary_json(minimal_qa_stats, self.temp_dir)

        # Check that CSV files are created
        ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_halogen_ops_atoms.csv')
        rule_ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_ops_atoms.csv')

        self.assertTrue(os.path.exists(ops_atoms_pivot_path),
                       "pivot_rule_halogen_ops_atoms.csv should be created even when empty")
        self.assertTrue(os.path.exists(rule_ops_atoms_pivot_path),
                       "pivot_rule_ops_atoms.csv should be created even when empty")

        # Check that they have headers even if no data
        try:
            df1 = pd.read_csv(ops_atoms_pivot_path)
            df2 = pd.read_csv(rule_ops_atoms_pivot_path)

            # Should have column headers even if no rows
            self.assertGreater(len(df1.columns), 0, "Should have column headers")
            self.assertGreater(len(df2.columns), 0, "Should have column headers")

        except Exception as e:
            self.fail(f"Failed to read empty pivot CSV files: {e}")

    def test_no_pivot_csvs_for_non_v2_data(self):
        """Test that ops/atoms pivot CSVs are NOT generated for non-version-2 data."""
        # Create legacy QA stats without version 2 format
        legacy_qa_stats = {
            'attempts': 5,
            'products': 3,
            'no_product_matches': 2,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 1},
            'dedup_hits_statesig': 0,
            'dedup_hits_inchi': 0
            # No 'version' or 'pivots' keys
        }

        # Write QA summary
        write_qa_summary_json(legacy_qa_stats, self.temp_dir)

        # Check that ops/atoms pivot CSV files are NOT created
        ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_halogen_ops_atoms.csv')
        rule_ops_atoms_pivot_path = os.path.join(self.temp_dir, 'pivot_rule_ops_atoms.csv')

        self.assertFalse(os.path.exists(ops_atoms_pivot_path),
                        "pivot_rule_halogen_ops_atoms.csv should NOT be created for non-v2 data")
        self.assertFalse(os.path.exists(rule_ops_atoms_pivot_path),
                        "pivot_rule_ops_atoms.csv should NOT be created for non-v2 data")


if __name__ == '__main__':
    unittest.main()