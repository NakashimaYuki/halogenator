# -*- coding: ascii -*-
"""Test RDKit module isolation at CLI routing layer."""

import unittest
from unittest.mock import patch, MagicMock
import sys
import io


class TestRDKitModuleIsolation(unittest.TestCase):
    """Test RDKit module isolation at CLI routing layer."""
    
    def setUp(self):
        """Set up test with mocked availability flags."""
        # Clear any cached modules
        modules_to_clear = [
            'src.halogenator.cli',
            'src.halogenator.rules', 
            'src.halogenator.ring_tag'
        ]
        for module in modules_to_clear:
            if module in sys.modules:
                del sys.modules[module]
    
    def test_cmd_enum_fails_gracefully_without_rules(self):
        """cmd_enum should exit gracefully when rules module is unavailable."""
        
        # Mock rules module import failure
        with patch.dict('sys.modules', {'src.halogenator.rules': None}):
            with patch('src.halogenator.cli.RULES_AVAILABLE', False):
                with patch('src.halogenator.cli._rules_error', 'Mock RDKit unavailable'):
                    
                    # Capture stdout and test sys.exit
                    with self.assertRaises(SystemExit) as context:
                        with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                            from src.halogenator.cli import cmd_enum
                            
                            config = {}
                            args = MagicMock()
                            cmd_enum(config, args)
                    
                    # Should exit with code 1
                    self.assertEqual(context.exception.code, 1)
    
    def test_cmd_k1_fails_gracefully_without_rules(self):
        """cmd_k1 should exit gracefully when rules module is unavailable."""
        
        with patch('src.halogenator.cli.RULES_AVAILABLE', False):
            with patch('src.halogenator.cli._rules_error', 'Mock RDKit unavailable'):
                
                with self.assertRaises(SystemExit) as context:
                    with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                        from src.halogenator.cli import cmd_k1
                        
                        config = {}
                        args = MagicMock()
                        cmd_k1(config, args)
                        
                        # Check that error message was printed
                        output = mock_stdout.getvalue()
                        self.assertIn("k=1 enumeration requires RDKit", output)
                        self.assertIn("Please install RDKit", output)
                
                # Should exit with code 1  
                self.assertEqual(context.exception.code, 1)
    
    def test_cmd_benchmark_fails_gracefully_without_rules(self):
        """cmd_benchmark should exit gracefully when rules module is unavailable."""
        
        with patch('src.halogenator.cli.RULES_AVAILABLE', False):
            with patch('src.halogenator.cli._rules_error', 'Mock RDKit unavailable'):
                
                with self.assertRaises(SystemExit) as context:
                    with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                        from src.halogenator.cli import cmd_benchmark
                        
                        args = MagicMock()
                        args.standard = True
                        cmd_benchmark(args)
                        
                        # Check that error message was printed
                        output = mock_stdout.getvalue()
                        self.assertIn("Benchmarks require RDKit", output)
                        self.assertIn("Please install RDKit", output)
                
                # Should exit with code 1
                self.assertEqual(context.exception.code, 1)
    
    def test_ring_tag_unavailable_shows_warning_not_error(self):
        """cmd_enum should show warning but continue when ring_tag is unavailable."""
        
        with patch('src.halogenator.cli.RULES_AVAILABLE', True):
            with patch('src.halogenator.cli.RING_TAG_AVAILABLE', False):
                with patch('src.halogenator.cli._ring_tag_error', 'Mock ring tag unavailable'):
                    
                    # Mock other dependencies to avoid deeper execution
                    with patch('src.halogenator.cli.read_smi', return_value=[]):
                        with patch('src.halogenator.cli.reset_reaction_warning_counts'):
                            with patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
                                try:
                                    from src.halogenator.cli import cmd_enum
                                    
                                    config = {'halogens': ['F'], 'io': {}}
                                    args = MagicMock()
                                    args.subset = 'flavonoids'
                                    args.k = 1
                                    cmd_enum(config, args)
                                    
                                except Exception:
                                    # Expected to fail deeper in execution, but should pass ring_tag check
                                    pass
                                
                                # Check that warning (not error) was printed
                                output = mock_stdout.getvalue()
                                self.assertIn("WARNING: Ring tagging functionality unavailable", output)
                                self.assertNotIn("ERROR:", output)
    
    def test_modules_load_successfully_when_rdkit_available(self):
        """Modules should load successfully when RDKit is available."""
        
        # This test verifies the normal case where RDKit is available
        with patch('src.halogenator.cli.RULES_AVAILABLE', True):
            with patch('src.halogenator.cli.RING_TAG_AVAILABLE', True):
                
                # Should be able to import the CLI module without errors
                try:
                    from src.halogenator.cli import cmd_enum, cmd_k1, cmd_benchmark
                    # If we get here, imports succeeded
                    self.assertTrue(True)
                except ImportError:
                    self.fail("CLI commands should be importable when RDKit is available")


if __name__ == '__main__':
    unittest.main()