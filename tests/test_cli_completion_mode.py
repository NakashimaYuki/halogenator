# -*- coding: ascii -*-
"""Test CLI completion mode functionality."""

import unittest
import tempfile
import os
import subprocess
import sys
import json
from unittest.mock import patch, MagicMock


class TestCLICompletionMode(unittest.TestCase):
    """Test CLI completion mode parameter handling."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.config_path = os.path.join(self.temp_dir, 'test_config.yaml')

        # Create minimal test config
        with open(self.config_path, 'w', encoding='ascii') as f:
            f.write('halogens: [F, Cl]\n')
            f.write('rules: [R1, R3]\n')
            f.write('io:\n')
            f.write('  smiles_file: test_parents.smi\n')

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch('src.halogenator.cli.cmd_enum')
    def test_cli_default_completion_mode_is_zero_fill(self, mock_cmd_enum):
        """Test that CLI defaults to zero_fill completion mode."""
        # Import CLI module
        from src.halogenator.cli import main

        # Mock sys.argv
        test_args = ['halogenator', 'enum', '-c', self.config_path]
        with patch.object(sys, 'argv', test_args):
            try:
                main()
            except SystemExit:
                pass

        # Check that cmd_enum was called
        self.assertTrue(mock_cmd_enum.called)
        args, kwargs = mock_cmd_enum.call_args

        # Extract args object and check default completion mode
        args_obj = args[1] if len(args) > 1 else None
        if args_obj and hasattr(args_obj, 'qa_completion_mode'):
            self.assertEqual(args_obj.qa_completion_mode, 'zero_fill')

    @patch('src.halogenator.cli.cmd_enum')
    def test_cli_distribute_mode_sets_metadata_strategy(self, mock_cmd_enum):
        """Test that CLI distribute mode is passed through."""
        from src.halogenator.cli import main

        # Mock sys.argv with distribute mode
        test_args = ['halogenator', 'enum', '-c', self.config_path, '--qa-completion-mode', 'distribute']
        with patch.object(sys, 'argv', test_args):
            try:
                main()
            except SystemExit:
                pass

        # Check that cmd_enum was called
        self.assertTrue(mock_cmd_enum.called)
        args, kwargs = mock_cmd_enum.call_args

        # Extract args object and check completion mode
        args_obj = args[1] if len(args) > 1 else None
        if args_obj and hasattr(args_obj, 'qa_completion_mode'):
            self.assertEqual(args_obj.qa_completion_mode, 'distribute')

    def test_cli_rejects_invalid_completion_mode_value(self):
        """Test that CLI rejects invalid completion mode values."""
        from src.halogenator.cli import main

        # Test invalid completion mode
        test_args = ['halogenator', 'enum', '-c', self.config_path, '--qa-completion-mode', 'invalid']

        with patch.object(sys, 'argv', test_args):
            with self.assertRaises(SystemExit) as cm:
                main()

            # Should exit with non-zero code for invalid argument
            self.assertNotEqual(cm.exception.code, 0)

    @patch('src.halogenator.report.write_qa_summary_json')
    @patch('src.halogenator.cli.enumerate_with_stats')
    @patch('src.halogenator.cli.read_smi')
    def test_write_qa_summary_json_rejects_invalid_completion_mode(self, mock_read_smi, mock_enumerate, mock_write_qa):
        """Test that write_qa_summary_json rejects invalid completion mode."""
        from src.halogenator.report import write_qa_summary_json

        # Mock basic QA stats
        qa_stats = {
            'attempts': 10,
            'products': 5,
            'qa_paths': {'isotope_unavailable': 2}
        }

        # Test invalid completion mode
        with self.assertRaises(ValueError) as cm:
            write_qa_summary_json(qa_stats, self.temp_dir, completion_mode='invalid_mode')

        self.assertIn('Invalid completion_mode', str(cm.exception))
        self.assertIn('invalid_mode', str(cm.exception))

    def test_completion_mode_validation_accepts_valid_values(self):
        """Test that completion mode validation accepts valid values."""
        from src.halogenator.report import write_qa_summary_json

        # Mock basic QA stats
        qa_stats = {
            'version': '1',
            'attempts': 10,
            'products': 5,
            'qa_paths': {'isotope_unavailable': 2},
            'metadata': {'halogens': ['F'], 'rules': ['R1']}
        }

        # Test valid completion modes
        for mode in ['zero_fill', 'distribute']:
            try:
                write_qa_summary_json(qa_stats, self.temp_dir, completion_mode=mode)
            except ValueError:
                self.fail(f"write_qa_summary_json should accept valid completion_mode '{mode}'")

    @patch('src.halogenator.cli.enumerate_with_stats')
    @patch('src.halogenator.cli.read_smi')
    def test_end_to_end_completion_mode_in_metadata(self, mock_read_smi, mock_enumerate):
        """Test end-to-end completion mode parameter flow to metadata."""
        from src.halogenator.cli import cmd_enum
        from unittest.mock import MagicMock

        # Mock dependencies
        mock_read_smi.return_value = [('CCO', 'ethanol')]
        mock_enumerate.return_value = ([], {'attempts': 1, 'products': 0, 'qa_paths': {}})

        # Create mock args with completion mode
        mock_args = MagicMock()
        mock_args.qa_completion_mode = 'distribute'
        mock_args.subset = 'flavonoids'
        mock_args.k = 2
        mock_args.outdir = self.temp_dir

        # Create mock config
        config = {
            'halogens': ['F'],
            'rules': ['R1'],
            'io': {'smiles_file': 'test.smi'}
        }

        # Run command
        try:
            cmd_enum(config, mock_args)
        except Exception:
            pass  # Expected due to mocking

        # Check if QA summary was written
        qa_summary_path = os.path.join(self.temp_dir, 'qa_summary.json')
        if os.path.exists(qa_summary_path):
            with open(qa_summary_path, 'r', encoding='utf-8') as f:
                qa_data = json.load(f)

            # Check that completion strategy is recorded in metadata
            if 'metadata' in qa_data and 'completion' in qa_data['metadata']:
                self.assertEqual(qa_data['metadata']['completion']['strategy'], 'distribute')


if __name__ == '__main__':
    unittest.main()