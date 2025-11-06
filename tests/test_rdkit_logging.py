# -*- coding: ascii -*-
"""Tests for RDKit logging functionality."""

import unittest
import logging
import builtins
import sys
import os
import tempfile
from unittest.mock import patch, MagicMock
from io import StringIO


class TestRDKitLogging(unittest.TestCase):
    """Test RDKit logging threshold and CLI switch functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Clear any existing handlers
        self.rdkit_logger = logging.getLogger('rdkit')
        self.rdkit_logger.handlers = []
        
        # Create a string buffer to capture log messages
        self.log_capture = StringIO()
        self.handler = logging.StreamHandler(self.log_capture)
        self.handler.setLevel(logging.DEBUG)
        
        # Set up logger to capture all RDKit messages
        self.rdkit_logger.addHandler(self.handler)
        self.rdkit_logger.setLevel(logging.DEBUG)
        
    def tearDown(self):
        """Clean up test fixtures."""
        self.rdkit_logger.removeHandler(self.handler)
        self.handler.close()
    
    def test_rdkit_logging_levels(self):
        """Test that different log levels control RDKit message visibility."""
        from src.halogenator.cli import setup_rdkit_logging
        
        # Mock RDKit classes
        mock_rdlogger = MagicMock()
        
        with patch.dict('sys.modules', {'rdkit': MagicMock()}):
            with patch('src.halogenator.cli.RDLogger', mock_rdlogger):
                # Test WARNING level (should disable debug and info)
                setup_rdkit_logging('WARNING')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.debug')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.info')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.warning')
                
                # Reset mock
                mock_rdlogger.reset_mock()
                
                # Test ERROR level (should disable debug, info, warning)
                setup_rdkit_logging('ERROR')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.debug')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.info')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.warning')
                # Verify rdApp.error is NEVER disabled
                disable_calls = [call[0][0] for call in mock_rdlogger.DisableLog.call_args_list]
                self.assertNotIn('rdApp.error', disable_calls)
                
                # Reset mock
                mock_rdlogger.reset_mock()
                
                # Test INFO level (should disable only debug)
                setup_rdkit_logging('INFO')
                mock_rdlogger.DisableLog.assert_any_call('rdApp.debug')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.info')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.warning')
                
                # Reset mock
                mock_rdlogger.reset_mock()
                
                # Test DEBUG level (should enable all)
                setup_rdkit_logging('DEBUG')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.debug')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.info')
                mock_rdlogger.EnableLog.assert_any_call('rdApp.warning')
    
    def test_env_variable_priority(self):
        """Test that HALO_RDKIT_LOG_LEVEL environment variable is respected."""
        from src.halogenator.cli import configure_logging
        
        # Mock arguments without rdkit_log_level
        mock_args = MagicMock()
        mock_args.quiet = False
        mock_args.log_level = 'INFO'
        mock_args.rdkit_log_level = None
        delattr(mock_args, 'rdkit_log_level')  # Remove the attribute entirely
        
        # Test environment variable takes effect
        with patch.dict(os.environ, {'HALO_RDKIT_LOG_LEVEL': 'ERROR'}):
            with patch('src.halogenator.cli.setup_rdkit_logging') as mock_setup:
                configure_logging(mock_args)
                mock_setup.assert_called_with('ERROR')
        
        # Test that main log level is used as fallback
        with patch.dict(os.environ, {}, clear=True):
            with patch('src.halogenator.cli.setup_rdkit_logging') as mock_setup:
                configure_logging(mock_args)
                mock_setup.assert_called_with('INFO')  # Should fallback to main log level
    
    def test_cli_precedence_over_env(self):
        """Test that CLI --rdkit-log-level takes precedence over environment variable."""
        from src.halogenator.cli import configure_logging
        
        # Mock arguments with rdkit_log_level
        mock_args = MagicMock()
        mock_args.quiet = False
        mock_args.log_level = 'INFO'
        mock_args.rdkit_log_level = 'DEBUG'
        
        # Set environment variable to different value
        with patch.dict(os.environ, {'HALO_RDKIT_LOG_LEVEL': 'ERROR'}):
            with patch('src.halogenator.cli.setup_rdkit_logging') as mock_setup:
                configure_logging(mock_args)
                # Should use CLI value, not environment variable
                mock_setup.assert_called_with('DEBUG')
    
    def test_graceful_degradation_no_rdkit(self):
        """Test that RDKit logging setup degrades gracefully when RDKit is not available."""
        from src.halogenator.cli import setup_rdkit_logging
        
        # Mock ImportError when importing RDKit
        original_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name.startswith('rdkit'):
                raise ImportError("No module named 'rdkit'")
            return original_import(name, *args, **kwargs)

        with patch('builtins.__import__', side_effect=fake_import):
            sys.modules.pop('src.halogenator.chem_compat', None)
            sys.modules.pop('halogenator.chem_compat', None)
            with patch('logging.getLogger') as mock_get_logger:
                mock_logger = MagicMock()
                mock_get_logger.return_value = mock_logger

                setup_rdkit_logging('WARNING')

                mock_logger.debug.assert_called_once()
                self.assertIn("RDKit not available", mock_logger.debug.call_args[0][0])
    
    def test_graceful_degradation_configuration_failure(self):
        """Test graceful degradation when RDKit configuration fails."""
        from src.halogenator.cli import setup_rdkit_logging
        
        # Mock RDKit classes that raise an exception
        mock_rdlogger = MagicMock()
        mock_rdlogger.DisableLog.side_effect = Exception("Configuration failed")
        
        with patch.dict('sys.modules', {'rdkit': MagicMock()}):
            with patch('src.halogenator.cli.RDLogger', mock_rdlogger):
                with patch('logging.getLogger') as mock_get_logger:
                    mock_logger = MagicMock()
                    mock_get_logger.return_value = mock_logger
                    
                    # Should not raise an exception
                    setup_rdkit_logging('WARNING')
                    
                    # Should log a debug message about configuration failure
                    mock_logger.debug.assert_called_once()
                    self.assertIn("Failed to configure RDKit logging", mock_logger.debug.call_args[0][0])
    
    def test_invalid_env_variable_ignored(self):
        """Test that invalid HALO_RDKIT_LOG_LEVEL values are ignored."""
        from src.halogenator.cli import configure_logging
        
        # Mock arguments without rdkit_log_level
        mock_args = MagicMock()
        mock_args.quiet = False
        mock_args.log_level = 'WARNING'
        mock_args.rdkit_log_level = None
        delattr(mock_args, 'rdkit_log_level')
        
        # Test invalid environment variable is ignored
        with patch.dict(os.environ, {'HALO_RDKIT_LOG_LEVEL': 'INVALID'}):
            with patch('src.halogenator.cli.setup_rdkit_logging') as mock_setup:
                configure_logging(mock_args)
                # Should fallback to main log level since env var is invalid
                mock_setup.assert_called_with('WARNING')
        
        # Test case-insensitive valid values work
        with patch.dict(os.environ, {'HALO_RDKIT_LOG_LEVEL': 'warning'}):
            with patch('src.halogenator.cli.setup_rdkit_logging') as mock_setup:
                configure_logging(mock_args)
                mock_setup.assert_called_with('WARNING')  # Should be normalized to uppercase
    
    def test_error_logging_never_disabled(self):
        """Test that rdApp.error logging is never disabled."""
        from src.halogenator.cli import setup_rdkit_logging
        
        mock_rdlogger = MagicMock()
        
        with patch.dict('sys.modules', {'rdkit': MagicMock()}):
            with patch('src.halogenator.cli.RDLogger', mock_rdlogger):
                # Test all levels to ensure rdApp.error is never disabled
                for level in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
                    mock_rdlogger.reset_mock()
                    setup_rdkit_logging(level)
                    
                    # Verify rdApp.error is never disabled
                    disable_calls = [call[0][0] for call in mock_rdlogger.DisableLog.call_args_list]
                    self.assertNotIn('rdApp.error', disable_calls, 
                                   f"rdApp.error should never be disabled at level {level}")


if __name__ == '__main__':
    unittest.main()