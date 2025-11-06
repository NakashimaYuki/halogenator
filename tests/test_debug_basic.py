# -*- coding: ascii -*-
"""Debug basic enumeration functionality."""

import unittest
import os

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import enumerate_with_stats, EnumConfig


class TestDebugBasic(unittest.TestCase):
    """Debug basic enumeration to understand current behavior."""
    
    def test_basic_k1_legacy(self):
        """Test basic k=1 enumeration without guard."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
        )
        
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
        
        print(f"K=1 Legacy - QA stats keys: {list(qa_stats.keys())}")
        print(f"K=1 Legacy - Version: {qa_stats.get('version')}")
        print(f"K=1 Legacy - QA stats: {qa_stats}")
        
        # Basic assertions
        self.assertIsInstance(products, list)
        self.assertIsInstance(qa_stats, dict)
    
    def test_basic_k1_v2(self):
        """Test basic k=1 enumeration v2 without guard."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
        )
        
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='v2')
        
        print(f"K=1 V2 - QA stats keys: {list(qa_stats.keys())}")
        print(f"K=1 V2 - Version: {qa_stats.get('version')}")
        print(f"K=1 V2 - Has pivots: {'pivots' in qa_stats}")
        
        # Basic assertions
        self.assertIsInstance(products, list)
        self.assertIsInstance(qa_stats, dict)
    
    def test_basic_k_gt_1(self):
        """Test basic k>1 enumeration without guard."""
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            rules=('R1',),
        )
        
        products, qa_stats = enumerate_with_stats("c1ccccc1", cfg)
        
        print(f"K>1 - QA stats keys: {list(qa_stats.keys())}")
        print(f"K>1 - Version: {qa_stats.get('version')}")
        print(f"K>1 - Has pivots: {'pivots' in qa_stats}")
        
        # Basic assertions
        self.assertIsInstance(products, list)
        self.assertIsInstance(qa_stats, dict)


if __name__ == '__main__':
    unittest.main()