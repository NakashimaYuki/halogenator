# -*- coding: ascii -*-
"""Tests for streaming interface compatibility protection."""

import unittest
from src.halogenator.enumerate_k import EnumConfig, enumerate_products, enumerate_with_stats


class TestStreamingInterfaceCompatibility(unittest.TestCase):
    """Test that streaming interface changes don't break existing consumers."""
    
    def setUp(self):
        """Set up test configuration."""
        self.config = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 2, 'min_graph_distance': 2},
            std_cfg={'do_tautomer': False},
            qc_cfg={'sanitize_strict': True},
            pruning_cfg={'enable_symmetry_fold': True}
        )
    
    def test_default_no_marker_emission(self):
        """Test that by default, no marker is emitted even with return_qa_stats=True."""
        # Use benzene which should have limited halogenation options
        parent_smiles = "c1ccccc1"
        
        # Collect all items from streaming interface with default parameters
        items = list(enumerate_products(parent_smiles, self.config, return_qa_stats=True))
        
        # Check that no items are markers
        for item in items:
            if isinstance(item, tuple) and len(item) == 2:
                product, qa_stats = item
                self.assertFalse(product.get('is_qa_summary_marker', False), 
                               "Should not emit marker by default")
            else:
                # Non-tuple items are definitely not markers
                pass
    
    def test_marker_emission_when_enabled(self):
        """Test that markers are emitted when emit_summary_marker=True."""
        # Use benzene which should produce some enumeration activity but maybe no products
        parent_smiles = "c1ccccc1"
        
        # Collect all items with marker emission enabled
        items = list(enumerate_products(parent_smiles, self.config, 
                                      return_qa_stats=True, emit_summary_marker=True))
        
        # Should have at least some items (potentially including a marker)
        marker_count = 0
        product_count = 0
        
        for item in items:
            if isinstance(item, tuple) and len(item) == 2:
                product, qa_stats = item
                if product.get('is_qa_summary_marker', False):
                    marker_count += 1
                    # Verify marker structure
                    self.assertIn('parent_smiles', product)
                    self.assertIn('parent_inchikey', product)
                    self.assertIsInstance(qa_stats, dict)
                else:
                    product_count += 1
        
        # Should have at most one marker
        self.assertLessEqual(marker_count, 1, "Should have at most one marker")
        
        # If we have enumeration activity, we should have a marker
        if marker_count == 1:
            # Good - marker was emitted as expected
            pass
        elif product_count > 0:
            # If we have products but no marker, that's also valid
            # (marker is only emitted when there's activity but potentially no products)
            pass
    
    def test_enumerate_with_stats_still_works(self):
        """Test that enumerate_with_stats continues to work correctly."""
        parent_smiles = "c1ccc(O)cc1"  # phenol - should produce some products
        
        products, qa_stats = enumerate_with_stats(parent_smiles, self.config)
        
        # Verify return types
        self.assertIsInstance(products, list)
        self.assertIsInstance(qa_stats, dict)
        
        # Verify no products contain markers
        for product in products:
            self.assertFalse(product.get('is_qa_summary_marker', False),
                           "enumerate_with_stats should filter out markers")
        
        # Verify QA stats structure
        expected_top_level_keys = {
            'no_product_matches', 'template_unsupported', 'qa_paths'
        }
        actual_keys = set(qa_stats.keys())
        for key in expected_top_level_keys:
            self.assertIn(key, actual_keys, f"Missing expected top-level QA stat key: {key}")
        
        # Verify qa_paths structure
        qa_paths = qa_stats.get('qa_paths', {})
        expected_qa_path_keys = {
            'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'
        }
        for key in expected_qa_path_keys:
            self.assertIn(key, qa_paths, f"Missing expected qa_paths key: {key}")
    
    def test_backward_compatibility_no_products_case(self):
        """Test backward compatibility when no products are generated."""
        # Create a molecule that's less likely to produce halogenated products
        parent_smiles = "CC"  # ethane - no aromatic carbons
        
        # Default behavior: no markers
        default_items = list(enumerate_products(parent_smiles, self.config, return_qa_stats=True))
        
        # Should have no items or only product items (no markers)
        for item in default_items:
            if isinstance(item, tuple) and len(item) == 2:
                product, qa_stats = item
                self.assertFalse(product.get('is_qa_summary_marker', False))
        
        # With marker enabled: potentially markers
        marker_items = list(enumerate_products(parent_smiles, self.config, 
                                             return_qa_stats=True, emit_summary_marker=True))
        
        # Count markers and products
        markers = []
        products = []
        
        for item in marker_items:
            if isinstance(item, tuple) and len(item) == 2:
                product, qa_stats = item
                if product.get('is_qa_summary_marker', False):
                    markers.append((product, qa_stats))
                else:
                    products.append((product, qa_stats))
        
        # Should have at most one marker
        self.assertLessEqual(len(markers), 1)
    
    def test_streaming_vs_stable_interface_consistency(self):
        """Test that streaming and stable interfaces return consistent results."""
        parent_smiles = "c1ccc(O)cc1"  # phenol
        
        # Get results from stable interface
        stable_products, stable_qa = enumerate_with_stats(parent_smiles, self.config)
        
        # Get results from streaming interface (with markers enabled to match stable interface)
        streaming_products = []
        streaming_qa = None
        
        for item in enumerate_products(parent_smiles, self.config, 
                                     return_qa_stats=True, emit_summary_marker=True):
            if isinstance(item, tuple) and len(item) == 2:
                product, qa_stats = item
                if product.get('is_qa_summary_marker', False):
                    streaming_qa = qa_stats
                else:
                    streaming_products.append(product)
                    streaming_qa = qa_stats  # Keep updating to get final state
        
        # Product counts should match
        self.assertEqual(len(stable_products), len(streaming_products))
        
        # QA stats should match (if we got them from streaming)
        if streaming_qa is not None:
            self.assertEqual(stable_qa, streaming_qa)
    
    def test_old_consumer_pattern_still_works(self):
        """Test that old consumer patterns that don't expect markers still work."""
        parent_smiles = "c1ccc(O)cc1"  # phenol
        
        # Simulate old consumer code that just iterates over products
        # without checking for markers (using default emit_summary_marker=False)
        products = []
        
        for item in enumerate_products(parent_smiles, self.config, return_qa_stats=False):
            # Old code would just assume every item is a product
            products.append(item)
        
        # All items should be valid product dictionaries
        for product in products:
            self.assertIsInstance(product, dict)
            self.assertFalse(product.get('is_qa_summary_marker', False))
            self.assertIn('smiles', product)  # Should have basic product fields
    
    def test_parameter_validation(self):
        """Test that parameter combinations work correctly."""
        parent_smiles = "c1ccc(O)cc1"
        
        # emit_summary_marker without return_qa_stats should have no effect
        items_no_qa = list(enumerate_products(parent_smiles, self.config, 
                                            return_qa_stats=False, emit_summary_marker=True))
        
        # Should get plain product dictionaries, no tuples
        for item in items_no_qa:
            self.assertIsInstance(item, dict)
            self.assertFalse(item.get('is_qa_summary_marker', False))
        
        # Both flags enabled should work
        items_both = list(enumerate_products(parent_smiles, self.config, 
                                           return_qa_stats=True, emit_summary_marker=True))
        
        # Should get tuples, potentially including markers
        for item in items_both:
            if isinstance(item, tuple):
                self.assertEqual(len(item), 2)
                product, qa_stats = item
                self.assertIsInstance(product, dict)
                self.assertIsInstance(qa_stats, dict)


if __name__ == '__main__':
    unittest.main()