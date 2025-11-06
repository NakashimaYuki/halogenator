# -*- coding: ascii -*-
"""Tests for QA statistics integration with reports."""

import json
import os
import tempfile
import unittest
import pandas as pd
from src.halogenator.report import write_qa_summary_json, _write_summary_csv


class TestQAReportIntegration(unittest.TestCase):
    """Test QA statistics are properly written to reports."""
    
    def test_write_qa_summary_json(self):
        """Test QA summary JSON writing functionality."""
        qa_stats = {
            'no_product_matches': 5,
            'dedup_hits_statesig': 12,
            'dedup_hits_inchi': 3,
            'template_unsupported': 2,
            'qa_paths': {
                'isotope_unavailable': 8,
                'isotope_miss': 1,
                'atommap_used': 15,
                'heuristic_used': 4
            }
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write QA summary JSON
            json_path = write_qa_summary_json(qa_stats, tmpdir)
            
            # Verify file exists
            self.assertTrue(os.path.exists(json_path))
            self.assertEqual(os.path.basename(json_path), 'qa_summary.json')
            
            # Verify JSON content
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Check structure
            self.assertIn('total', data)
            self.assertIn('metadata', data)
            
            # Check total stats match input
            self.assertEqual(data['total'], qa_stats)
            
            # Check metadata structure
            metadata = data['metadata']
            self.assertIn('description', metadata)
            self.assertIn('semantics', metadata)
            self.assertIn('generated_at', metadata)
            
            # Verify specific QA fields
            total_stats = data['total']
            self.assertEqual(total_stats['no_product_matches'], 5)
            self.assertEqual(total_stats['template_unsupported'], 2)
            
            qa_paths = total_stats['qa_paths']
            self.assertEqual(qa_paths['isotope_unavailable'], 8)
            self.assertEqual(qa_paths['atommap_used'], 15)
    
    def test_qa_json_schema(self):
        """Test QA summary JSON has expected schema."""
        qa_stats = {
            'no_product_matches': 0,
            'dedup_hits_statesig': 0,
            'dedup_hits_inchi': 0,
            'template_unsupported': 0,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = write_qa_summary_json(qa_stats, tmpdir)
            
            # Load and validate schema
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Check required top-level keys
            required_keys = {'total', 'metadata'}
            self.assertEqual(set(data.keys()), required_keys)
            
            # Check total structure
            total = data['total']
            expected_total_keys = {
                'no_product_matches', 'dedup_hits_statesig', 'dedup_hits_inchi',
                'template_unsupported', 'qa_paths'
            }
            self.assertEqual(set(total.keys()), expected_total_keys)
            
            # Check qa_paths structure
            qa_paths = total['qa_paths']
            expected_qa_keys = {
                'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'
            }
            self.assertEqual(set(qa_paths.keys()), expected_qa_keys)
            
            # Check all values are integers
            for key, value in total.items():
                if key == 'qa_paths':
                    for qa_key, qa_value in value.items():
                        self.assertIsInstance(qa_value, int, f"qa_paths.{qa_key} should be int")
                else:
                    self.assertIsInstance(value, int, f"{key} should be int")
    
    def test_report_contains_qa_diagnostics_rows(self):
        """Test that CSV report contains QA diagnostics when JSON is present."""
        qa_stats = {
            'no_product_matches': 3,
            'dedup_hits_statesig': 8,
            'dedup_hits_inchi': 2,
            'template_unsupported': 1,
            'qa_paths': {
                'isotope_unavailable': 5,
                'isotope_miss': 0,
                'atommap_used': 7,
                'heuristic_used': 2
            }
        }
        
        # Mock stats for CSV generation
        mock_stats = {
            'parent_count': 10,
            'product_count': 50,
            'sanitize_ok_count': 45,
            'unique_parents_with_products': 8,
            'parents_without_products': ['parent1', 'parent2'],
            'avg_products_per_parent': 5.0,
            'max_products_per_parent': 15,
            'min_products_per_parent': 1,
            'k_counts': {1: 50},
            'rule_counts': {'R1': 20, 'R3': 30},
            'halogen_counts': {'F': 25, 'Cl': 25},
            'rule_halogen_counts': {'R1': {'F': 10, 'Cl': 10}, 'R3': {'F': 15, 'Cl': 15}},
            'rule_halogen_k_counts': {
                'R1': {'F': {1: 10}, 'Cl': {1: 10}},
                'R3': {'F': {1: 15}, 'Cl': {1: 15}}
            },
            'has_qc_errors': False,
            'diagnostics': {
                'diff_smiles_vs_inchikey': 2,
                'smiles_unique_parents': 8,
                'inchikey_unique_parents': 6,
                'missing_in_smiles_examples': ['example1'],
                'missing_in_ikeys_examples': ['example2'],
                'missing_in_smiles_count': 1,
                'missing_in_ikeys_count': 1
            },
            'config': {'subset': 'test'}
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write QA summary JSON first
            qa_json_path = write_qa_summary_json(qa_stats, tmpdir)
            self.assertTrue(os.path.exists(qa_json_path))
            
            # Generate CSV report
            csv_path = os.path.join(tmpdir, 'summary.csv')
            _write_summary_csv(mock_stats, csv_path)
            
            # Read and verify CSV content
            df = pd.read_csv(csv_path)
            
            # Check that QA diagnostics rows are present
            diagnostics_df = df[df['Category'] == 'Diagnostics']
            diagnostic_metrics = diagnostics_df['Metric'].tolist()
            
            expected_qa_metrics = [
                'Isotope unavailable (attempts)',
                'Isotope miss (matches)', 
                'AtomMap used (attempts)',
                'Heuristic used (attempts)',
                'No-product matches',
                'Template unsupported (attempts)'
            ]
            
            for metric in expected_qa_metrics:
                self.assertIn(metric, diagnostic_metrics, f"Missing QA metric: {metric}")
            
            # Check specific values (CSV values are strings)
            for _, row in diagnostics_df.iterrows():
                metric = row['Metric']
                value = int(row['Value'])  # Convert from CSV string to int
                
                if metric == 'Isotope unavailable (attempts)':
                    self.assertEqual(value, 5)
                elif metric == 'AtomMap used (attempts)':
                    self.assertEqual(value, 7)
                elif metric == 'Template unsupported (attempts)':
                    self.assertEqual(value, 1)
                elif metric == 'No-product matches':
                    self.assertEqual(value, 3)


if __name__ == '__main__':
    unittest.main()