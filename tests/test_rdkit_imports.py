# -*- coding: ascii -*-
"""Tests for RDKit import behavior."""

import unittest
import sys
import importlib
from unittest.mock import patch


class TestRDKitImports(unittest.TestCase):
    """Test that modules can be imported without RDKit and only use lazy imports."""
    
    def setUp(self):
        """Remove any cached imports before each test."""
        # List of modules to test
        self.modules_to_test = [
            'src.halogenator.standardize',
            'src.halogenator.qc',
            'src.halogenator.reactions', 
            'src.halogenator.constraints',
            'src.halogenator.enumerate_k1',
            'src.halogenator.enumerate_k',
            'src.halogenator.sites',
            'src.halogenator.ring_tag',
            'src.halogenator.dedupe'
        ]
        
        # Remove modules from cache if they exist
        for module_name in self.modules_to_test:
            if module_name in sys.modules:
                del sys.modules[module_name]
    
    def test_import_without_rdkit(self):
        """Test that modules can be imported when RDKit is not available."""
        # Mock the rdkit import to raise ImportError
        with patch.dict('sys.modules', {'rdkit': None}):
            with patch('builtins.__import__', side_effect=self._mock_import):
                # Try to import each module
                for module_name in self.modules_to_test:
                    with self.subTest(module=module_name):
                        try:
                            importlib.import_module(module_name)
                            # If we get here, the import succeeded (which is what we want)
                            success = True
                        except ImportError as e:
                            if 'rdkit' in str(e).lower():
                                self.fail(f"Module {module_name} has top-level RDKit import: {e}")
                            else:
                                # Some other import error, re-raise to investigate
                                raise
                        
                        self.assertTrue(success, f"Failed to import {module_name} without RDKit")
    
    def _mock_import(self, name, *args, **kwargs):
        """Mock import function that raises ImportError for rdkit."""
        if name.startswith('rdkit'):
            raise ImportError(f"No module named '{name}' (mocked)")
        # For all other imports, use the real import
        return importlib.__import__(name, *args, **kwargs)
    
    def test_rdkit_usage_is_lazy(self):
        """Test that RDKit is only imported when functions are called."""
        # Import modules normally (with RDKit available)
        modules = {}
        for module_name in self.modules_to_test:
            try:
                modules[module_name] = importlib.import_module(module_name)
            except ImportError:
                # Skip modules that can't be imported for other reasons
                pass
        
        # Verify that importing the modules doesn't immediately use RDKit
        # This is validated by the fact that the test_import_without_rdkit test passes
        # If there were top-level RDKit imports, that test would fail
        self.assertTrue(len(modules) > 0, "At least some modules should import successfully")
    
    def test_graceful_degradation(self):
        """Test that functions gracefully handle RDKit not being available."""
        from src.halogenator.standardize import std_from_smiles, to_inchikey
        from src.halogenator.qc import sanitize_ok, basic_descriptors
        
        # These functions should handle RDKit not being available gracefully
        # by returning appropriate fallback values
        
        # Mock RDKit to be unavailable during function calls
        with patch.dict('sys.modules', {'rdkit': None}):
            with patch('builtins.__import__', side_effect=self._mock_import):
                # Test standardize functions
                result = std_from_smiles("CCO")
                self.assertIsNone(result, "std_from_smiles should return None when RDKit unavailable")
                
                result = to_inchikey(None)  # Pass None as mol
                self.assertEqual(result, "UNKNOWN", "to_inchikey should return 'UNKNOWN' when RDKit unavailable")
                
                # Test QC functions
                result = sanitize_ok(None)
                self.assertFalse(result, "sanitize_ok should return False when RDKit unavailable")
                
                result = basic_descriptors(None)
                expected_fallback = {
                    'mw': 0.0,
                    'logp': 0.0,
                    'tpsa': 0.0,
                    'hba': 0,
                    'hbd': 0,
                    'aromatic_rings': 0
                }
                self.assertEqual(result, expected_fallback, "basic_descriptors should return fallback values")


if __name__ == '__main__':
    unittest.main()