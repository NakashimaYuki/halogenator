# -*- coding: ascii -*-
"""Tests for version handling logic in CLI."""

import unittest
import tempfile
import os
from unittest.mock import patch, MagicMock
from src.halogenator.cli import cmd_enum


class TestVersionHandling(unittest.TestCase):
    """Test that version handling follows correct precedence rules."""
    
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
            f.write('c1ccccc1O\\tphenol\\n')
            f.write('CCN\\tethylamine\\n')
        
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
    def test_version_2_takes_precedence(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test that version '2' takes precedence over any other version."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return different versions
        mock_enumerate.side_effect = [
            # First molecule returns version '1'
            (
                [{'smiles': 'c1ccc(F)cc1O', 'rule': 'R1', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '1',
                    'no_product_matches': 1,
                    'template_unsupported': 0
                }
            ),
            # Second molecule returns version '2'
            (
                [{'smiles': 'CCF', 'rule': 'R4', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '2',
                    'pivots': {'by_rule': {'R4': {'products': 1}}},
                    'no_product_matches': 0,
                    'template_unsupported': 1
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = None
        mock_args.outdir = None
        
        with patch('src.halogenator.cli.write_qa_summary_json') as mock_write_qa:
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            # Verify final version is '2' (takes precedence over '1')
            self.assertEqual(written_qa_stats.get('version'), '2')
    
    @patch('src.halogenator.enumerate_k.enumerate_with_stats')
    @patch('src.halogenator.cli.write_table')
    @patch('src.halogenator.enumerate_k.print_reaction_warning_summary')
    def test_version_1_preserved_when_no_version_2(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test that version '1' is preserved when no molecule has version '2'."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return version '1' for both
        mock_enumerate.side_effect = [
            # First molecule returns version '1'
            (
                [{'smiles': 'c1ccc(F)cc1O', 'rule': 'R1', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '1',
                    'no_product_matches': 1,
                    'template_unsupported': 0
                }
            ),
            # Second molecule also returns version '1'
            (
                [{'smiles': 'CCF', 'rule': 'R4', 'halogen': 'F', 'depth': 1}],
                {
                    'version': '1',
                    'no_product_matches': 0,
                    'template_unsupported': 1
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = None
        mock_args.outdir = None
        
        with patch('src.halogenator.cli.write_qa_summary_json') as mock_write_qa:
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            # Verify final version is '1' (no version '2' to override)
            self.assertEqual(written_qa_stats.get('version'), '1')
    
    @patch('src.halogenator.enumerate_k.enumerate_with_stats')
    @patch('src.halogenator.cli.write_table')
    @patch('src.halogenator.enumerate_k.print_reaction_warning_summary')
    def test_absent_version_preserved(self, mock_warning_summary, mock_write_table, mock_enumerate):
        """Test that absent version is preserved when no molecules specify version."""
        test_config, smiles_file = self.create_test_config_and_files()
        
        # Mock enumerate_with_stats to return no version for both
        mock_enumerate.side_effect = [
            # First molecule returns no version
            (
                [{'smiles': 'c1ccc(F)cc1O', 'rule': 'R1', 'halogen': 'F', 'depth': 1}],
                {
                    'no_product_matches': 1,
                    'template_unsupported': 0
                }
            ),
            # Second molecule also returns no version
            (
                [{'smiles': 'CCF', 'rule': 'R4', 'halogen': 'F', 'depth': 1}],
                {
                    'no_product_matches': 0,
                    'template_unsupported': 1
                }
            )
        ]
        
        mock_write_table.return_value = None
        mock_warning_summary.return_value = None
        
        mock_args = MagicMock()
        mock_args.k = 2
        mock_args.subset = None
        mock_args.outdir = None
        
        with patch('src.halogenator.cli.write_qa_summary_json') as mock_write_qa:
            qa_json_path = os.path.join(self.temp_dir, 'qa_summary.json')
            mock_write_qa.return_value = qa_json_path
            
            # Run the command
            cmd_enum(test_config, mock_args)
            
            # Verify write_qa_summary_json was called once
            self.assertEqual(mock_write_qa.call_count, 1)
            written_qa_stats = mock_write_qa.call_args[0][0]
            
            # Verify final stats has no version key (absent version preserved)
            self.assertNotIn('version', written_qa_stats)


if __name__ == '__main__':
    unittest.main()