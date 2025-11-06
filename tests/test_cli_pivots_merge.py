# -*- coding: ascii -*-
"""Tests for CLI pivots merging functionality."""

import unittest
import json
import tempfile
import os
from unittest.mock import patch, MagicMock
from src.halogenator.cli import cmd_enum


class TestCLIPivotsMerge(unittest.TestCase):
    """Test that pivots data is properly merged in CLI."""
    
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
        """Create test configuration and SMILES files."""
        # Create temporary SMILES file with 2 molecules
        smiles_file = os.path.join(self.temp_dir, 'parents.smi')
        with open(smiles_file, 'w') as f:
            f.write('c1ccccc1O\tphenol\n')
            f.write('CCN\tethylamine\n')
        
        test_config = {
            'io': {
                'smiles_file': smiles_file,
                'products_table': os.path.join(self.temp_dir, 'products_k2.parquet')
            },
            'halogens': ['F', 'Cl'],
            'k_max': 2,
            'constraints': {'per_ring_quota': 2, 'min_graph_distance': 2},
            'standardize': {'do_tautomer': False},
            'qc': {'sanitize_strict': True},
            'pruning': {'enable_symmetry_fold': True}
        }
        
        return test_config, smiles_file
    
    @patch('src.halogenator.enumerate_k.enumerate_with_stats')
    @patch('src.halogenator.cli.write_table')
    @patch('src.halogenator.enumerate_k.print_reaction_warning_summary')
    def test_pivots_merge_correctly(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test that pivots from multiple parents are properly merged."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return known pivots for each molecule
        mock_enumerate.side_effect = [
            # First molecule (phenol) - version 2 with small pivots
            (
                [{'smiles': 'c1ccc(F)cc1O', 'rule': 'R1', 'halogen': 'F', 'depth': 1}],
                {
                    'no_product_matches': 0,
                    'template_unsupported': 1,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 0,
                    'qa_paths': {'isotope_unavailable': 2, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
                    'version': '2',
                    'pivots': {
                        'by_rule': {'R1': {'products': 1, 'attempts': 2}},
                        'by_halogen': {'F': {'products': 1, 'attempts': 2}},
                        'by_k': {1: {'products': 1, 'attempts': 2}},
                        'by_rule_halogen': {'R1_F': {'products': 1, 'attempts': 2}},
                        'by_rule_halogen_k': {'R1_F_1': {'products': 1, 'attempts': 2}}
                    }
                }
            ),
            # Second molecule (ethylamine) - version 2 with different pivots
            (
                [{'smiles': 'CCF', 'rule': 'R4', 'halogen': 'F', 'depth': 1}],
                {
                    'no_product_matches': 1,
                    'template_unsupported': 0,
                    'dedup_hits_statesig': 0,
                    'dedup_hits_inchi': 1,
                    'qa_paths': {'isotope_unavailable': 1, 'isotope_miss': 1, 'atommap_used': 0, 'heuristic_used': 0},
                    'version': '2',
                    'pivots': {
                        'by_rule': {'R4': {'products': 1, 'attempts': 3}},
                        'by_halogen': {'F': {'products': 1, 'attempts': 3}},
                        'by_k': {1: {'products': 1, 'attempts': 3}},
                        'by_rule_halogen': {'R4_F': {'products': 1, 'attempts': 3}},
                        'by_rule_halogen_k': {'R4_F_1': {'products': 1, 'attempts': 3}}
                    }
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        # Create mock args
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = 'flavonoids'
        mock_args.outdir = None
        
        with patch('src.halogenator.cli.write_qa_summary_json') as mock_write_qa:
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            # Verify version is set to '2'
            self.assertEqual(written_qa_stats.get('version'), '2')
            
            # Verify pivots are properly merged
            pivots = written_qa_stats.get('pivots', {})
            self.assertIn('by_rule', pivots)
            self.assertIn('by_halogen', pivots)
            
            # Check by_rule merging: R1 and R4 should both be present
            by_rule = pivots['by_rule']
            self.assertIn('R1', by_rule)
            self.assertIn('R4', by_rule)
            self.assertEqual(by_rule['R1']['products'], 1)  # From first molecule
            self.assertEqual(by_rule['R1']['attempts'], 2)
            self.assertEqual(by_rule['R4']['products'], 1)  # From second molecule
            self.assertEqual(by_rule['R4']['attempts'], 3)
            
            # Check by_halogen merging: F appears in both, should be summed
            by_halogen = pivots['by_halogen']
            self.assertIn('F', by_halogen)
            self.assertEqual(by_halogen['F']['products'], 2)  # 1 + 1
            self.assertEqual(by_halogen['F']['attempts'], 5)  # 2 + 3
            
            # Check by_k merging: k=1 appears in both, should be summed
            by_k = pivots['by_k']
            self.assertIn(1, by_k)
            self.assertEqual(by_k[1]['products'], 2)  # 1 + 1
            self.assertEqual(by_k[1]['attempts'], 5)  # 2 + 3
            
            # Check by_rule_halogen merging: different combinations
            by_rule_halogen = pivots['by_rule_halogen']
            self.assertIn('R1_F', by_rule_halogen)
            self.assertIn('R4_F', by_rule_halogen)
            self.assertEqual(by_rule_halogen['R1_F']['products'], 1)
            self.assertEqual(by_rule_halogen['R4_F']['products'], 1)
            
            # Check by_rule_halogen_k merging: different combinations
            by_rule_halogen_k = pivots['by_rule_halogen_k']
            self.assertIn('R1_F_1', by_rule_halogen_k)
            self.assertIn('R4_F_1', by_rule_halogen_k)
            
            # Verify base totals are also correct
            self.assertEqual(written_qa_stats['no_product_matches'], 1)  # 0 + 1
            self.assertEqual(written_qa_stats['template_unsupported'], 1)  # 1 + 0
            self.assertEqual(written_qa_stats['qa_paths']['isotope_unavailable'], 3)  # 2 + 1
            self.assertEqual(written_qa_stats['qa_paths']['isotope_miss'], 1)  # 0 + 1
    
    @patch('src.halogenator.enumerate_k.enumerate_with_stats')
    @patch('src.halogenator.cli.write_table')
    @patch('src.halogenator.enumerate_k.print_reaction_warning_summary')
    def test_pivots_disjoint_and_overlapping_keys(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test pivot merging with both disjoint and overlapping keys across all dimensions."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return complex pivot patterns
        mock_enumerate.side_effect = [
            # First molecule: multiple rules, halogens, events
            (
                [{'smiles': 'c1ccc(F)cc1O', 'rule': 'R1', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '2',
                    'pivots': {
                        'by_rule': {
                            'R1': {'products': 2, 'attempts': 4, 'template_unsupported': 1},
                            'R3': {'products': 1, 'attempts': 3, 'no_product_matches': 2}
                        },
                        'by_halogen': {
                            'F': {'products': 2, 'attempts': 5, 'template_unsupported': 1},
                            'Cl': {'products': 1, 'attempts': 2}
                        },
                        'by_k': {
                            1: {'products': 3, 'attempts': 7, 'template_unsupported': 1}
                        },
                        'by_rule_halogen': {
                            'R1_F': {'products': 2, 'attempts': 4, 'template_unsupported': 1},
                            'R3_Cl': {'products': 1, 'attempts': 3}
                        },
                        'by_rule_halogen_k': {
                            'R1_F_1': {'products': 2, 'attempts': 4, 'template_unsupported': 1},
                            'R3_Cl_1': {'products': 1, 'attempts': 3}
                        }
                    }
                }
            ),
            # Second molecule: overlapping and disjoint keys
            (
                [{'smiles': 'CCF', 'rule': 'R4', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '2',
                    'pivots': {
                        'by_rule': {
                            'R1': {'products': 1, 'attempts': 2},  # overlaps with first
                            'R4': {'products': 3, 'attempts': 5, 'isotope_unavailable': 2}  # disjoint
                        },
                        'by_halogen': {
                            'F': {'products': 3, 'attempts': 4, 'isotope_unavailable': 1},  # overlaps
                            'Br': {'products': 1, 'attempts': 3}  # disjoint
                        },
                        'by_k': {
                            1: {'products': 4, 'attempts': 7, 'isotope_unavailable': 2}  # overlaps
                        },
                        'by_rule_halogen': {
                            'R1_F': {'products': 1, 'attempts': 2},  # overlaps with first
                            'R4_F': {'products': 2, 'attempts': 3, 'isotope_unavailable': 1},  # disjoint
                            'R4_Br': {'products': 1, 'attempts': 2}  # disjoint
                        },
                        'by_rule_halogen_k': {
                            'R1_F_1': {'products': 1, 'attempts': 2},  # overlaps with first
                            'R4_F_1': {'products': 2, 'attempts': 3, 'isotope_unavailable': 1},  # disjoint
                            'R4_Br_1': {'products': 1, 'attempts': 2}  # disjoint
                        }
                    }
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = 'flavonoids'
        mock_args.outdir = None
        
        with patch('src.halogenator.cli.write_qa_summary_json') as mock_write_qa:
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            # Verify version is set to '2'
            self.assertEqual(written_qa_stats.get('version'), '2')
            
            # Verify pivots are properly merged
            pivots = written_qa_stats.get('pivots', {})
            
            # Test by_rule dimension: overlapping R1 should be summed, disjoint R3/R4 should be copied
            by_rule = pivots['by_rule']
            # R1 appears in both - should be summed
            self.assertEqual(by_rule['R1']['products'], 3)  # 2 + 1
            self.assertEqual(by_rule['R1']['attempts'], 6)   # 4 + 2
            self.assertEqual(by_rule['R1']['template_unsupported'], 1)  # 1 + 0
            # R3 only in first - should be copied exactly
            self.assertEqual(by_rule['R3']['products'], 1)
            self.assertEqual(by_rule['R3']['attempts'], 3)
            self.assertEqual(by_rule['R3']['no_product_matches'], 2)
            # R4 only in second - should be copied exactly
            self.assertEqual(by_rule['R4']['products'], 3)
            self.assertEqual(by_rule['R4']['attempts'], 5)
            self.assertEqual(by_rule['R4']['isotope_unavailable'], 2)
            
            # Test by_halogen dimension: F overlaps, Cl/Br are disjoint
            by_halogen = pivots['by_halogen']
            # F appears in both - should be summed
            self.assertEqual(by_halogen['F']['products'], 5)  # 2 + 3
            self.assertEqual(by_halogen['F']['attempts'], 9)  # 5 + 4
            self.assertEqual(by_halogen['F']['template_unsupported'], 1)  # 1 + 0
            self.assertEqual(by_halogen['F']['isotope_unavailable'], 1)   # 0 + 1
            # Cl only in first - should be copied exactly
            self.assertEqual(by_halogen['Cl']['products'], 1)
            self.assertEqual(by_halogen['Cl']['attempts'], 2)
            # Br only in second - should be copied exactly
            self.assertEqual(by_halogen['Br']['products'], 1)
            self.assertEqual(by_halogen['Br']['attempts'], 3)
            
            # Test by_k dimension: k=1 overlaps
            by_k = pivots['by_k']
            self.assertEqual(by_k[1]['products'], 7)  # 3 + 4
            self.assertEqual(by_k[1]['attempts'], 14) # 7 + 7
            self.assertEqual(by_k[1]['template_unsupported'], 1)  # 1 + 0
            self.assertEqual(by_k[1]['isotope_unavailable'], 2)   # 0 + 2
            
            # Test by_rule_halogen dimension: R1_F overlaps, others disjoint
            by_rule_halogen = pivots['by_rule_halogen']
            # R1_F appears in both - should be summed
            self.assertEqual(by_rule_halogen['R1_F']['products'], 3)  # 2 + 1
            self.assertEqual(by_rule_halogen['R1_F']['attempts'], 6)  # 4 + 2
            self.assertEqual(by_rule_halogen['R1_F']['template_unsupported'], 1)  # 1 + 0
            # Disjoint keys should be copied exactly
            self.assertEqual(by_rule_halogen['R3_Cl']['products'], 1)
            self.assertEqual(by_rule_halogen['R3_Cl']['attempts'], 3)
            self.assertEqual(by_rule_halogen['R4_F']['products'], 2)
            self.assertEqual(by_rule_halogen['R4_F']['attempts'], 3)
            self.assertEqual(by_rule_halogen['R4_F']['isotope_unavailable'], 1)
            self.assertEqual(by_rule_halogen['R4_Br']['products'], 1)
            self.assertEqual(by_rule_halogen['R4_Br']['attempts'], 2)
            
            # Test by_rule_halogen_k dimension: R1_F_1 overlaps, others disjoint
            by_rule_halogen_k = pivots['by_rule_halogen_k']
            # R1_F_1 appears in both - should be summed
            self.assertEqual(by_rule_halogen_k['R1_F_1']['products'], 3)  # 2 + 1
            self.assertEqual(by_rule_halogen_k['R1_F_1']['attempts'], 6)  # 4 + 2
            self.assertEqual(by_rule_halogen_k['R1_F_1']['template_unsupported'], 1)  # 1 + 0
            # Disjoint keys should be copied exactly
            self.assertEqual(by_rule_halogen_k['R3_Cl_1']['products'], 1)
            self.assertEqual(by_rule_halogen_k['R3_Cl_1']['attempts'], 3)
            self.assertEqual(by_rule_halogen_k['R4_F_1']['products'], 2)
            self.assertEqual(by_rule_halogen_k['R4_F_1']['attempts'], 3)
            self.assertEqual(by_rule_halogen_k['R4_F_1']['isotope_unavailable'], 1)
            self.assertEqual(by_rule_halogen_k['R4_Br_1']['products'], 1)
            self.assertEqual(by_rule_halogen_k['R4_Br_1']['attempts'], 2)


if __name__ == '__main__':
    unittest.main()