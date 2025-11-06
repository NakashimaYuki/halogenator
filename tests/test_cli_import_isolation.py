# -*- coding: ascii -*-
"""Test CLI import isolation without RDKit dependency."""

import unittest
import sys
import importlib
from unittest.mock import patch


class TestCLIImportIsolation(unittest.TestCase):
    """Test that CLI can be imported without RDKit and provides good error messages."""
    
    def test_cli_import_without_rdkit(self):
        """Test that CLI module can be imported even when RDKit is not available."""
        # The CLI should already be importable since we're running this test
        # But let's test the specific isolation behavior
        
        # Import CLI module - this should work even without RDKit
        import src.halogenator.cli as cli
        
        # Check that availability flags exist
        self.assertTrue(hasattr(cli, 'RULES_AVAILABLE'))
        self.assertTrue(hasattr(cli, 'IO_UTILS_AVAILABLE'))
        self.assertTrue(hasattr(cli, 'STANDARDIZE_AVAILABLE'))
        self.assertTrue(hasattr(cli, 'RING_TAG_AVAILABLE'))
        
        # Check that error variables exist when modules are not available
        if not cli.RULES_AVAILABLE:
            self.assertTrue(hasattr(cli, '_rules_error'))
            self.assertIsNotNone(cli._rules_error)
    
    def test_cli_help_works_without_rdkit_deps(self):
        """Test that CLI help can be displayed without RDKit functionality."""
        import src.halogenator.cli as cli
        
        # The main() function should exist and be callable for help
        self.assertTrue(hasattr(cli, 'main'))
        self.assertTrue(callable(cli.main))
        
        # Test that argument parser can be created (this is what --help uses)
        # We can't easily test --help output without subprocess, but we can
        # test that the basic CLI structure is intact
        
        # Check that command functions exist
        self.assertTrue(hasattr(cli, 'cmd_k1'))
        self.assertTrue(hasattr(cli, 'cmd_enum'))
        self.assertTrue(hasattr(cli, 'cmd_report'))
        self.assertTrue(hasattr(cli, 'cmd_benchmark'))
    
    def test_cmd_functions_check_rdkit_availability(self):
        """Test that command functions check RDKit availability and give helpful errors."""
        import src.halogenator.cli as cli
        
        # When RDKit-dependent modules are not available, commands should check and error gracefully
        with patch('src.halogenator.cli.RULES_AVAILABLE', False):
            with patch('src.halogenator.cli._rules_error', 'Mock RDKit unavailable'):
                # cmd_k1 should check RULES_AVAILABLE and exit gracefully
                with self.assertRaises(SystemExit):
                    cli.cmd_k1({})
                
                # cmd_enum should check RULES_AVAILABLE and exit gracefully  
                with self.assertRaises(SystemExit):
                    cli.cmd_enum({})
                
                # cmd_benchmark should check RULES_AVAILABLE and exit gracefully
                with self.assertRaises(SystemExit):
                    cli.cmd_benchmark(type('args', (), {})())
    
    def test_local_imports_in_functions(self):
        """Test that functions use local imports for RDKit-dependent modules."""
        import src.halogenator.cli as cli
        
        # Check that the top-level imports comment exists (indicating imports were moved)
        import inspect
        cli_source = inspect.getsource(cli)
        self.assertIn("Imports moved to function-local scope", cli_source)
        
        # Check that cmd_k1 has local import
        cmd_k1_source = inspect.getsource(cli.cmd_k1)
        self.assertIn("from .enumerate_k1 import", cmd_k1_source)
        
        # Check that cmd_report has local import  
        cmd_report_source = inspect.getsource(cli.cmd_report)
        self.assertIn("from .report import", cmd_report_source)
        
        # Check that cmd_enum has local import
        cmd_enum_source = inspect.getsource(cli.cmd_enum)
        self.assertIn("from .report import write_qa_summary_json", cmd_enum_source)


if __name__ == '__main__':
    unittest.main()