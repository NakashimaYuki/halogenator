"""
Test configuration for halogenator test suite.

This module sets up common test fixtures and ensures consistent import paths
across all test modules by adding the src directory to sys.path.
"""
import os
import sys

# Add src directory to Python path for consistent imports
test_dir = os.path.dirname(__file__)
src_dir = os.path.join(test_dir, '..', 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)