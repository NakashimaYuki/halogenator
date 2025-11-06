# -*- coding: ascii -*-
"""Test RDKit guard functionality."""

import os
import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.guard import rdkit_guard, rdkit_guard_enabled
from src.halogenator.enumerate_k1 import enumerate_k1_with_stats


class TestRDKitGuard(unittest.TestCase):
    """Test RDKit guard context manager and wiring."""
    
    def setUp(self):
        # Clear environment before each test
        os.environ.pop('HALO_RDKIT_GUARD', None)
    
    def test_rdkit_guard_disabled_by_default(self):
        """Guard should be disabled by default."""
        self.assertFalse(rdkit_guard_enabled())
    
    def test_rdkit_guard_enabled_with_env_var(self):
        """Guard should be enabled when environment variable is set."""
        os.environ['HALO_RDKIT_GUARD'] = '1'
        self.assertTrue(rdkit_guard_enabled())
        
        os.environ['HALO_RDKIT_GUARD'] = '0'
        self.assertFalse(rdkit_guard_enabled())
    
    def test_guard_catches_exception_when_enabled(self):
        """Guard should catch exceptions and record QA stats when enabled."""
        os.environ['HALO_RDKIT_GUARD'] = '1'
        
        qa_stats = {'qa_paths': {}, 'template_unsupported': 0}
        
        with rdkit_guard(qa_stats):
            raise ValueError("Simulated RDKit error")
        
        # Should have recorded the error
        self.assertEqual(qa_stats['qa_paths']['rdkit_error'], 1)
        self.assertEqual(qa_stats['template_unsupported'], 1)
    
    def test_guard_reraises_exception_when_disabled(self):
        """Guard should reraise exceptions when disabled."""
        os.environ['HALO_RDKIT_GUARD'] = '0'
        
        qa_stats = {'qa_paths': {}, 'template_unsupported': 0}
        
        with self.assertRaises(ValueError):
            with rdkit_guard(qa_stats):
                raise ValueError("Simulated RDKit error")
        
        # Should not have recorded anything
        self.assertEqual(len(qa_stats['qa_paths']), 0)
        self.assertEqual(qa_stats['template_unsupported'], 0)
    
    def test_guard_integration_with_enumerate_k1(self):
        """Test guard integration with k=1 enumeration."""
        # Mock config object
        mock_config = MagicMock()
        mock_config.halogens = ['Cl']
        mock_config.constraints = {}
        mock_config.qc_cfg = {}
        mock_config.std_cfg = {}
        
        # Test with guard enabled and invalid SMILES
        os.environ['HALO_RDKIT_GUARD'] = '1'
        
        products, qa_stats = enumerate_k1_with_stats('invalid_smiles', mock_config)
        
        # Should not crash and should have QA stats
        self.assertEqual(len(products), 0)
        self.assertIn('qa_paths', qa_stats)
        # May have rdkit_error recorded depending on where the error occurs
        
        # Reset and test with guard disabled 
        os.environ['HALO_RDKIT_GUARD'] = '0'
        
        # Should still work (graceful handling of invalid SMILES)
        products, qa_stats = enumerate_k1_with_stats('invalid_smiles', mock_config)
        self.assertEqual(len(products), 0)
    
    def test_guard_does_not_interfere_with_valid_operations(self):
        """Guard should not interfere with valid RDKit operations."""
        os.environ['HALO_RDKIT_GUARD'] = '1'
        
        qa_stats = {'qa_paths': {}, 'template_unsupported': 0}
        
        # This should work normally
        result = None
        with rdkit_guard(qa_stats):
            result = "normal_operation"
        
        self.assertEqual(result, "normal_operation")
        # Should not have recorded any errors
        self.assertEqual(len(qa_stats['qa_paths']), 0)
        self.assertEqual(qa_stats['template_unsupported'], 0)


if __name__ == '__main__':
    unittest.main()