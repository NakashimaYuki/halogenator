# -*- coding: ascii -*-
"""Tests for business CSV comparison functionality."""

import unittest
import tempfile
import os
import pandas as pd
from src.halogenator.comparison import CSVComparison


class TestCompareBusinessCSVs(unittest.TestCase):
    """Test business CSV comparison functionality."""
    
    def setUp(self):
        """Set up test data."""
        self.test_dir = tempfile.mkdtemp()
        
        # Mock halogenation reference data
        self.ref_data = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},
            {'parent_inchikey': 'IK123', 'rule': 'R3', 'halogen': 'Cl', 'k': 1, 'product_smiles': 'ClC1=CC=CC=C1', 'site_id': 2},
            {'parent_inchikey': 'IK456', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=C(O)C=C1', 'site_id': 1}
        ])
        
        # Mock comparison data with partial overlap
        self.comp_data = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},  # Match
            {'parent_inchikey': 'IK123', 'rule': 'R4', 'halogen': 'Br', 'k': 1, 'product_smiles': 'BrC1=CC=CC=C1', 'site_id': 3},  # Extra
            {'parent_inchikey': 'IK456', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=C(O)C=C1', 'site_id': 1}   # Match
        ])
        
        self.ref_path = os.path.join(self.test_dir, 'ref.csv')
        self.comp_path = os.path.join(self.test_dir, 'comp.csv')
        
        self.ref_data.to_csv(self.ref_path, index=False)
        self.comp_data.to_csv(self.comp_path, index=False)
    
    def tearDown(self):
        """Clean up test directory."""
        import shutil
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_pair_matches_csv_structure(self):
        """Test pair_matches.csv has correct column structure."""
        comparison = CSVComparison(self.ref_path, self.comp_path)
        results = comparison.save_business_csv_reports(self.test_dir, 'test')
        
        self.assertIsNotNone(results)
        
        # Check pair_matches.csv
        pair_matches_df = pd.read_csv(results['pair_matches_file'])
        expected_cols = ['parent_id', 'product_id', 'site_id', 'rule', 'halogen', 'k', 'matched', 'ref_only', 'comp_only']
        
        for col in expected_cols:
            self.assertIn(col, pair_matches_df.columns, f"Missing column: {col}")
        
        # Should have 4 records (3 from ref + 1 extra from comp)
        self.assertEqual(len(pair_matches_df), 4)
        
        # Check match counts
        matches = pair_matches_df['matched'].sum()
        self.assertEqual(matches, 2)  # 2 matches
    
    def test_rule_confusion_csv_structure(self):
        """Test rule_confusion.csv has correct column structure."""
        comparison = CSVComparison(self.ref_path, self.comp_path)
        results = comparison.save_business_csv_reports(self.test_dir, 'test')
        
        # Check rule_confusion.csv
        rule_confusion_df = pd.read_csv(results['rule_confusion_file'])
        expected_cols = ['rule', 'halogen', 'k', 'tp', 'fp', 'fn']
        
        for col in expected_cols:
            self.assertIn(col, rule_confusion_df.columns, f"Missing column: {col}")
        
        # Should have records for each rule/halogen/k combination that appears
        self.assertGreater(len(rule_confusion_df), 0)
    
    def test_parent_summary_csv_structure(self):
        """Test parent_summary.csv has correct column structure."""
        comparison = CSVComparison(self.ref_path, self.comp_path)
        results = comparison.save_business_csv_reports(self.test_dir, 'test')
        
        # Check parent_summary.csv
        parent_summary_df = pd.read_csv(results['parent_summary_file'])
        expected_cols = ['parent_id', 'tp', 'fp', 'fn', 'coverage_pct']
        
        for col in expected_cols:
            self.assertIn(col, parent_summary_df.columns, f"Missing column: {col}")
        
        # Should have 2 parents (IK123, IK456)
        self.assertEqual(len(parent_summary_df), 2)
        
        # Check that coverage percentages are computed
        for idx, row in parent_summary_df.iterrows():
            self.assertGreaterEqual(row['coverage_pct'], 0)
            self.assertLessEqual(row['coverage_pct'], 100)
    
    def test_missing_columns_handled_gracefully(self):
        """Test that missing required columns are handled gracefully."""
        # Create data without required columns
        bad_ref_data = pd.DataFrame([
            {'smiles': 'FC1=CC=CC=C1', 'name': 'fluorobenzene'}
        ])
        
        bad_ref_path = os.path.join(self.test_dir, 'bad_ref.csv')
        bad_ref_data.to_csv(bad_ref_path, index=False)
        
        comparison = CSVComparison(bad_ref_path, self.comp_path)
        results = comparison.save_business_csv_reports(self.test_dir, 'test_bad')
        
        # Should return None when required columns are missing
        self.assertIsNone(results)
    
    def test_deduplication_and_key_design(self):
        """Test that deduplication works correctly and key design handles multiple products."""
        # Create data with duplicate records (same key) 
        ref_data_with_dupes = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},  # Duplicate
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC2=CC=CC=C2', 'site_id': 2}   # Different site_id
        ])
        
        comp_data_with_dupes = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},  # Match
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1}   # Duplicate
        ])
        
        ref_path_dupes = os.path.join(self.test_dir, 'ref_dupes.csv')
        comp_path_dupes = os.path.join(self.test_dir, 'comp_dupes.csv')
        
        ref_data_with_dupes.to_csv(ref_path_dupes, index=False)
        comp_data_with_dupes.to_csv(comp_path_dupes, index=False)
        
        comparison = CSVComparison(ref_path_dupes, comp_path_dupes)
        results = comparison.save_business_csv_reports(self.test_dir, 'test_dedup')
        
        # Read pair_matches results
        pair_matches_df = pd.read_csv(results['pair_matches_file'])
        
        # Should have 2 unique records (after deduplication):
        # 1. IK123|R1|F|1|1 - matched (in both ref and comp)
        # 2. IK123|R1|F|1|2 - ref_only (in ref but not comp)
        self.assertEqual(len(pair_matches_df), 2)
        
        # Verify one match and one ref_only
        matches = pair_matches_df['matched'].sum()
        ref_only = pair_matches_df['ref_only'].sum()
        self.assertEqual(matches, 1)
        self.assertEqual(ref_only, 1)
    
    def test_coverage_calculation_denominator(self):
        """Test that coverage calculation uses correct denominator: tp/(tp+fn) * 100."""
        # Create test data where coverage is clear
        ref_data = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'Cl', 'k': 1, 'product_smiles': 'ClC1=CC=CC=C1', 'site_id': 2},
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'Br', 'k': 1, 'product_smiles': 'BrC1=CC=CC=C1', 'site_id': 3}
        ])
        
        comp_data = pd.DataFrame([
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'F', 'k': 1, 'product_smiles': 'FC1=CC=CC=C1', 'site_id': 1},  # Match (tp=1)
            {'parent_inchikey': 'IK123', 'rule': 'R1', 'halogen': 'I', 'k': 1, 'product_smiles': 'IC1=CC=CC=C1', 'site_id': 4}    # Extra (fp=1)
            # Missing Cl and Br from ref (fn=2)
        ])
        
        ref_path_coverage = os.path.join(self.test_dir, 'ref_coverage.csv')
        comp_path_coverage = os.path.join(self.test_dir, 'comp_coverage.csv')
        
        ref_data.to_csv(ref_path_coverage, index=False)
        comp_data.to_csv(comp_path_coverage, index=False)
        
        comparison = CSVComparison(ref_path_coverage, comp_path_coverage)
        results = comparison.save_business_csv_reports(self.test_dir, 'test_coverage')
        
        # Read parent_summary results
        parent_summary_df = pd.read_csv(results['parent_summary_file'])
        
        # Should have one parent
        self.assertEqual(len(parent_summary_df), 1)
        
        row = parent_summary_df.iloc[0]
        self.assertEqual(row['tp'], 1)      # 1 match
        self.assertEqual(row['fn'], 2)      # 2 missing from comp (Cl, Br)
        self.assertEqual(row['fp'], 1)      # 1 extra in comp (I)
        
        # Coverage should be tp/(tp+fn) * 100 = 1/(1+2) * 100 = 33.33%
        # Note: fp is NOT included in denominator
        expected_coverage = (1 / (1 + 2)) * 100
        self.assertAlmostEqual(row['coverage_pct'], expected_coverage, places=2)


if __name__ == '__main__':
    unittest.main()