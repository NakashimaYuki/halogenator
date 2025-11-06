# -*- coding: ascii -*-
"""Simple test for RDKit module isolation at CLI routing layer."""

import unittest
from unittest.mock import patch


class TestRDKitCLIIsolation(unittest.TestCase):
    """Test RDKit module isolation at CLI routing layer with simple approach."""
    
    def test_rules_available_flag_exists(self):
        """Test that RULES_AVAILABLE flag is defined in CLI module."""
        from src.halogenator.cli import RULES_AVAILABLE
        self.assertIsInstance(RULES_AVAILABLE, bool)
    
    def test_ring_tag_available_flag_exists(self):
        """Test that RING_TAG_AVAILABLE flag is defined in CLI module."""
        from src.halogenator.cli import RING_TAG_AVAILABLE  
        self.assertIsInstance(RING_TAG_AVAILABLE, bool)
    
    def test_io_utils_available_flag_exists(self):
        """Test that IO_UTILS_AVAILABLE flag is defined in CLI module."""
        from src.halogenator.cli import IO_UTILS_AVAILABLE
        self.assertIsInstance(IO_UTILS_AVAILABLE, bool)
    
    def test_standardize_available_flag_exists(self):
        """Test that STANDARDIZE_AVAILABLE flag is defined in CLI module."""
        from src.halogenator.cli import STANDARDIZE_AVAILABLE
        self.assertIsInstance(STANDARDIZE_AVAILABLE, bool)
    
    def test_cmd_enum_has_rdkit_availability_check(self):
        """Test that cmd_enum function contains RDKit availability checks."""
        import inspect
        from src.halogenator.cli import cmd_enum
        
        # Get source code of the function
        source = inspect.getsource(cmd_enum)
        
        # Should contain checks for RULES_AVAILABLE
        self.assertIn('RULES_AVAILABLE', source)
        self.assertIn('Please install RDKit', source)
    
    def test_cmd_k1_has_rdkit_availability_check(self):
        """Test that cmd_k1 function contains RDKit availability checks."""
        import inspect
        from src.halogenator.cli import cmd_k1
        
        # Get source code of the function
        source = inspect.getsource(cmd_k1)
        
        # Should contain checks for RULES_AVAILABLE
        self.assertIn('RULES_AVAILABLE', source)
        self.assertIn('Please install RDKit', source)
    
    def test_cmd_benchmark_has_rdkit_availability_check(self):
        """Test that cmd_benchmark function contains RDKit availability checks."""
        import inspect
        from src.halogenator.cli import cmd_benchmark
        
        # Get source code of the function
        source = inspect.getsource(cmd_benchmark)
        
        # Should contain checks for RULES_AVAILABLE
        self.assertIn('RULES_AVAILABLE', source)
        self.assertIn('Please install RDKit', source)
    
    def test_availability_flags_are_true_when_rdkit_available(self):
        """Test that availability flags are True when RDKit is available (integration test)."""
        from src.halogenator.cli import (
            RULES_AVAILABLE, 
            RING_TAG_AVAILABLE, 
            IO_UTILS_AVAILABLE, 
            STANDARDIZE_AVAILABLE
        )
        
        # Since we're running in an environment with RDKit installed, these should be True
        # This is an integration test that verifies our conditional imports work correctly
        self.assertTrue(RULES_AVAILABLE, "rules module should be available with RDKit installed")
        self.assertTrue(RING_TAG_AVAILABLE, "ring_tag module should be available with RDKit installed") 
        self.assertTrue(IO_UTILS_AVAILABLE, "io_utils module should be available with RDKit installed")
        self.assertTrue(STANDARDIZE_AVAILABLE, "standardize module should be available with RDKit installed")


if __name__ == '__main__':
    unittest.main()