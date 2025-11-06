# -*- coding: ascii -*-
"""Test that legacy stream_shape outputs version:'1'."""

import unittest
import os
import tempfile
import shutil

from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import EnumConfig


class TestLegacyStreamShapeVersion(unittest.TestCase):
    """Test that legacy stream_shape emits version:'1'."""
    
    def test_legacy_stream_shape_emits_version_1(self):
        """Test that --stream-shape legacy outputs version:'1'."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Test legacy format
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
        
        # Should have version:'1' for legacy
        self.assertEqual(qa_stats.get('version'), '1')
        
        # Should have basic QA structure
        self.assertIn('no_product_matches', qa_stats)
        self.assertIn('template_unsupported', qa_stats)
        self.assertIn('qa_paths', qa_stats)
        
        # Should not have pivots (legacy format doesn't include them)
        self.assertNotIn('pivots', qa_stats)
    
    def test_v2_stream_shape_emits_version_2(self):
        """Test that --stream-shape v2 outputs version:'2'."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Test v2 format
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='v2')
        
        # Should have version:'2' for v2
        self.assertEqual(qa_stats.get('version'), '2')
        
        # Should have basic QA structure
        self.assertIn('no_product_matches', qa_stats)
        self.assertIn('template_unsupported', qa_stats)
        self.assertIn('qa_paths', qa_stats)
        
        # Should have pivots (v2 format includes them)
        self.assertIn('pivots', qa_stats)
    
    def test_legacy_format_has_basic_consistency_check(self):
        """Test that legacy format performs basic consistency validation."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            constraints={'per_ring_quota': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': False}
        )
        
        # Test with a molecule that should produce products
        products, qa_stats = enumerate_k1_with_stats("c1ccccc1", cfg, stream_shape='legacy')
        
        # Basic consistency: should have non-negative counts
        self.assertGreaterEqual(qa_stats.get('no_product_matches', 0), 0)
        self.assertGreaterEqual(qa_stats.get('template_unsupported', 0), 0)
        
        # QA paths should be present and non-negative
        qa_paths = qa_stats.get('qa_paths', {})
        for metric in ['isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used']:
            self.assertGreaterEqual(qa_paths.get(metric, 0), 0, f"{metric} should be >= 0")


if __name__ == '__main__':
    unittest.main()