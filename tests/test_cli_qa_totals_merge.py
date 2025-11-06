# -*- coding: ascii -*-
"""Tests for CLI QA totals merging and printing consistency."""

import unittest
import json
import tempfile
import os
import sys
from io import StringIO
from unittest.mock import patch, MagicMock
from src.halogenator.cli import cmd_enum


class TestCLIQATotalsMerge(unittest.TestCase):
    """Test that QA totals are properly merged and printed."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.addCleanup(lambda: self._cleanup_temp_dir())
        
    def _cleanup_temp_dir(self):
        """Clean up temporary directory."""
        import shutil
        try:
            shutil.rmtree(self.temp_dir)
        except Exception:
            pass
    
    def create_test_config_and_files(self):
        """Create test configuration and SMILES file."""
        # Create test SMILES file
        smiles_file = os.path.join(self.temp_dir, 'test_parents.smi')
        with open(smiles_file, 'w') as f:
            f.write("c1ccc(O)cc1\tphenol\n")
            f.write("CCN\tethylamine\n")
        
        # Create test config
        test_config = {
            'halogens': ['F', 'Cl'],
            'constraints': {'per_ring_quota': 2, 'min_graph_distance': 2},
            'standardize': {'do_tautomer': False},
            'qc': {'sanitize_strict': True},
            'pruning': {'enable_symmetry_fold': True},
            'io': {
                'smiles_file': smiles_file,
                'products_table': os.path.join(self.temp_dir, 'test_products.parquet')
            }
        }
        
        return test_config, smiles_file
    
    @patch('src.halogenator.enumerate_k.enumerate_with_stats')
    @patch('src.halogenator.cli.write_table')
    @patch('src.halogenator.enumerate_k.print_reaction_warning_summary')  # Silence this for cleaner test output
    def test_qa_merge_and_print_consistency(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test that QA statistics are properly merged and printed values match JSON."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return known QA stats for each molecule
        mock_enumerate.side_effect = [
            # First molecule (phenol)
            (
                [{'smiles': 'c1ccc(F)cc1', 'rule': 'R1', 'halogen': 'F'}],
                {
                    'no_product_matches': 1,
                    'template_unsupported': 2,
                    'dedup_hits_statesig': 3,
                    'dedup_hits_inchi': 4,
                    'qa_paths': {
                        'isotope_unavailable': 5,
                        'isotope_miss': 6,
                        'atommap_used': 7,
                        'heuristic_used': 8
                    }
                }
            ),
            # Second molecule (ethylamine) 
            (
                [{'smiles': 'CCFl', 'rule': 'R4', 'halogen': 'F'}],
                {
                    'no_product_matches': 10,
                    'template_unsupported': 20,
                    'dedup_hits_statesig': 30,
                    'dedup_hits_inchi': 40,
                    'qa_paths': {
                        'isotope_unavailable': 50,
                        'isotope_miss': 60,
                        'atommap_used': 70,
                        'heuristic_used': 80
                    }
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        # Create mock args with --quiet to reduce noise
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = 'flavonoids'
        mock_args.outdir = None
        
        # Capture stdout to verify printed values
        captured_output = StringIO()
        
        with patch('sys.stdout', captured_output), \
             patch('src.halogenator.report.write_qa_summary_json') as mock_write_qa:
            
            # Mock QA summary JSON write to capture the data
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            
            # Verify merged totals are correct
            expected_totals = {
                'no_product_matches': 11,        # 1 + 10
                'template_unsupported': 22,      # 2 + 20
                'dedup_hits_statesig': 33,       # 3 + 30
                'dedup_hits_inchi': 44,          # 4 + 40
                'qa_paths': {
                    'isotope_unavailable': 55,   # 5 + 50
                    'isotope_miss': 66,          # 6 + 60
                    'atommap_used': 77,          # 7 + 70
                    'heuristic_used': 88         # 8 + 80
                }
            }
            
            # Check that all expected values are in the written QA stats
            for key, expected_value in expected_totals.items():
                if key == 'qa_paths':
                    for path_key, path_value in expected_value.items():
                        self.assertEqual(written_qa_stats['qa_paths'][path_key], path_value, 
                                       f"Mismatch in qa_paths[{path_key}]: expected {path_value}, got {written_qa_stats['qa_paths'].get(path_key)}")
                else:
                    self.assertEqual(written_qa_stats[key], expected_value, 
                                   f"Mismatch in {key}: expected {expected_value}, got {written_qa_stats.get(key)}")
        
        # Verify printed output matches the merged values
        output_lines = captured_output.getvalue().split('\n')
        
        # Check specific printed values match our expected merged totals
        isotope_unavailable_line = [line for line in output_lines if 'Isotope unavailable:' in line]
        self.assertEqual(len(isotope_unavailable_line), 1)
        self.assertIn('55', isotope_unavailable_line[0])
        
        isotope_miss_line = [line for line in output_lines if 'Isotope misses:' in line]
        self.assertEqual(len(isotope_miss_line), 1)
        self.assertIn('66', isotope_miss_line[0])
        
        atommap_line = [line for line in output_lines if 'AtomMap fallback:' in line]
        self.assertEqual(len(atommap_line), 1)
        self.assertIn('77', atommap_line[0])
        
        heuristic_line = [line for line in output_lines if 'Heuristic fallback:' in line]
        self.assertEqual(len(heuristic_line), 1)
        self.assertIn('88', heuristic_line[0])
        
        template_unsupported_line = [line for line in output_lines if 'Template unsupported:' in line]
        self.assertEqual(len(template_unsupported_line), 1)
        self.assertIn('22', template_unsupported_line[0])
        
        no_product_line = [line for line in output_lines if 'No product matches:' in line]
        self.assertEqual(len(no_product_line), 1)
        self.assertIn('11', no_product_line[0])
        
        dedup_statesig_line = [line for line in output_lines if 'Dedup hits (state sig):' in line]
        self.assertEqual(len(dedup_statesig_line), 1)
        self.assertIn('33', dedup_statesig_line[0])
        
        dedup_inchi_line = [line for line in output_lines if 'Dedup hits (InChI):' in line]
        self.assertEqual(len(dedup_inchi_line), 1)
        self.assertIn('44', dedup_inchi_line[0])
        
        # Check that the correct label is used
        qa_path_events_line = [line for line in output_lines if 'Total QA path events:' in line]
        self.assertEqual(len(qa_path_events_line), 1)
        expected_total_qa_path_events = 55 + 66 + 77 + 88  # sum of qa_paths values
        self.assertIn(str(expected_total_qa_path_events), qa_path_events_line[0])
        
        # Verify old incorrect label is not used
        total_attempts_line = [line for line in output_lines if 'Total attempts:' in line]
        self.assertEqual(len(total_attempts_line), 0, "Old 'Total attempts' label should not be used")


if __name__ == '__main__':
    unittest.main()