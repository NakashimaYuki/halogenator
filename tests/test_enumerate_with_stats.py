# -*- coding: ascii -*-
"""Tests for stable QA stats interface."""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_products, enumerate_with_stats, EnumConfig


class TestEnumerateWithStats(unittest.TestCase):
    """Test stable QA statistics interface."""
    
    def test_enumerate_with_stats_returns_final_totals(self):
        """Test that enumerate_with_stats returns consistent final totals."""
        cfg = EnumConfig(k_max=1, halogens=('F', 'Cl'))
        
        # Test with simple substrate
        products, final_stats = enumerate_with_stats('c1ccccc1O', cfg)
        
        # Should have products
        self.assertTrue(len(products) > 0, "Should generate products")
        
        # Check QA stats structure (v2 stable shape)
        expected_keys = {'version', 'pivots', 'attempts', 'products', 'no_product_matches',
                         'template_unsupported', 'qa_paths', 'dedup_hits_statesig', 'dedup_hits_inchi'}
        self.assertEqual(set(final_stats.keys()), expected_keys)
        
        # Check qa_paths substructure - should contain at least the core QA keys
        qa_paths = final_stats['qa_paths']
        required_qa_keys = {'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'}
        self.assertTrue(required_qa_keys.issubset(set(qa_paths.keys())),
                       f"Missing required QA keys. Expected subset: {required_qa_keys}, got: {set(qa_paths.keys())}")
        
        # All stats should be integers >= 0
        for key, value in final_stats.items():
            if key == 'qa_paths':
                for qa_key, qa_value in value.items():
                    self.assertIsInstance(qa_value, int, f"qa_paths.{qa_key} should be int")
                    self.assertTrue(qa_value >= 0, f"qa_paths.{qa_key} should be non-negative")
            elif key == 'version':
                self.assertIsInstance(value, str)
            elif key == 'pivots':
                self.assertIsInstance(value, dict)
            else:
                self.assertIsInstance(value, int, f"{key} should be int")
                self.assertTrue(value >= 0, f"{key} should be non-negative")
    
    def test_enumerate_with_stats_vs_streaming_consistency(self):
        """Test that enumerate_with_stats final result matches streaming final snapshot."""
        from halogenator.enumerate_k import QAAggregator

        cfg = EnumConfig(k_max=1, halogens=('F',))

        # Get results from both interfaces
        products_stable, final_stats_stable = enumerate_with_stats('c1ccccc1O', cfg)

        # Get final snapshot from streaming interface with explicit aggregator
        aggregator = QAAggregator()
        streaming_results = list(enumerate_products('c1ccccc1O', cfg, return_qa_stats=True, stream_shape='v2', aggregator=aggregator))
        # Filter out QA summary markers from streaming products
        products_streaming = [item[0] for item in streaming_results
                             if not item[0].get('is_qa_summary_marker', False)]
        final_stats_streaming = streaming_results[-1][1] if streaming_results else {}

        # Product counts should match
        self.assertEqual(len(products_stable), len(products_streaming))

        # Final QA stats should match
        self.assertEqual(final_stats_stable, final_stats_streaming)

        # Product content should match (compare SMILES)
        stable_smiles = [p.get('product_smiles', '') for p in products_stable]
        streaming_smiles = [p.get('product_smiles', '') for p in products_streaming]
        self.assertEqual(sorted(stable_smiles), sorted(streaming_smiles))
    
    def test_enumerate_with_stats_no_products(self):
        """Test enumerate_with_stats behavior when no products are generated."""
        cfg = EnumConfig(k_max=1, halogens=('F',))
        
        # Use a substrate that won't generate products (no reactive sites)
        products, final_stats = enumerate_with_stats('C', cfg)  # Methane - no reactive sites
        
        # Should have no products but valid stats structure
        self.assertEqual(len(products), 0)
        
        # Stats structure should still be valid (v2 stable shape)
        expected_keys = {'version', 'pivots', 'attempts', 'products', 'no_product_matches',
                         'template_unsupported', 'qa_paths', 'dedup_hits_statesig', 'dedup_hits_inchi'}
        self.assertEqual(set(final_stats.keys()), expected_keys)
        
        # Some attempts should have been made (isotope_unavailable > 0)
        self.assertTrue(final_stats['qa_paths']['isotope_unavailable'] > 0,
                       "Should have attempts even with no products")


if __name__ == '__main__':
    unittest.main()
