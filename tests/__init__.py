# -*- coding: ascii -*-
"""Test package for halogenator."""

import unittest
import warnings
import os


class CleanLogsTestCase(unittest.TestCase):
    """Base test case that suppresses RDKit logging noise."""
    
    @classmethod
    def setUpClass(cls):
        """Set up class-level test fixtures - suppress RDKit logging noise."""
        super().setUpClass()
        
        # Suppress RDKit warnings for cleaner test output
        warnings.filterwarnings("ignore", category=DeprecationWarning, module="rdkit")
        warnings.filterwarnings("ignore", message=".*Kekulization.*", module="rdkit")
        
        # Set RDKit logger to WARNING level to reduce noise
        try:
            from src.halogenator.chem_compat import RDLogger
            RDLogger.DisableLog('rdApp.warning')
            RDLogger.DisableLog('rdApp.info')  
            RDLogger.DisableLog('rdApp.debug')
            # Keep error and critical levels enabled
        except ImportError:
            # RDKit not available, skip logging setup
            pass