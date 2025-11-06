# -*- coding: ascii -*-
"""Tests for totals-pivots equality validation."""

import unittest
import os
from src.halogenator.enumerate_k import _validate_totals_pivots_consistency, QAAggregator


class TestTotalsPivotsEquality(unittest.TestCase):
    """Test that totals exactly equal sum of pivots."""
    
    def setUp(self):
        """Set up test environment."""
        # Enable debug assertions for consistency checks
        os.environ['HALO_ASSERT_PIVOT_CONSISTENCY'] = '1'
    
    def tearDown(self):
        """Clean up test environment."""
        # Reset debug env
        if 'HALO_ASSERT_PIVOT_CONSISTENCY' in os.environ:
            del os.environ['HALO_ASSERT_PIVOT_CONSISTENCY']
    
    def test_equality_constraint_pass(self):
        """Test that equal totals and pivots pass validation."""
        # Create aggregator and add some test data
        aggregator = QAAggregator()
        
        # Add some test events
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=2)
        aggregator.record('products', rule='R1', halogen='F', k=1, amount=1)
        aggregator.record('no_product_matches', rule='R1', halogen='F', k=1, amount=1)
        aggregator.record('atommap_used', rule='R1', halogen='F', k=1, amount=1)
        
        # Build consistent QA stats dict
        qa_stats = {
            'version': '2',
            'attempts': 2,
            'products': 1,
            'no_product_matches': 1,
            'template_unsupported': 0,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 1,
                'heuristic_used': 0
            },
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should not raise exception
        try:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        except AssertionError:
            self.fail("Validation should pass for consistent data")
    
    def test_equality_constraint_fail_attempts(self):
        """Test that unequal attempts totals fail validation."""
        aggregator = QAAggregator()
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=3)
        
        # Inconsistent qa_stats (attempts mismatch)
        qa_stats = {
            'version': '2',
            'attempts': 2,  # Different from pivot sum (3)
            'products': 0,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should raise AssertionError
        with self.assertRaises(AssertionError) as cm:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        
        self.assertIn("attempts", str(cm.exception))
    
    def test_equality_constraint_fail_products(self):
        """Test that unequal products totals fail validation."""
        aggregator = QAAggregator()
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=2)
        aggregator.record('products', rule='R1', halogen='F', k=1, amount=2)
        
        # Inconsistent qa_stats (products mismatch)
        qa_stats = {
            'version': '2',
            'attempts': 2,
            'products': 1,  # Different from pivot sum (2)
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should raise AssertionError
        with self.assertRaises(AssertionError) as cm:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        
        self.assertIn("products", str(cm.exception))
    
    def test_equality_constraint_fail_qa_paths(self):
        """Test that unequal qa_paths fail validation."""
        aggregator = QAAggregator()
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=1)
        aggregator.record('atommap_used', rule='R1', halogen='F', k=1, amount=1)
        
        # Inconsistent qa_stats (atommap_used mismatch)
        qa_stats = {
            'version': '2',
            'attempts': 1,
            'products': 0,
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 2,  # Different from pivot sum (1)
                'heuristic_used': 0
            },
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should raise AssertionError
        with self.assertRaises(AssertionError) as cm:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        
        self.assertIn("atommap_used", str(cm.exception))
    
    def test_attempt_invariant_still_works(self):
        """Test that attempt invariants still work alongside equality checks."""
        aggregator = QAAggregator()
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=2)
        aggregator.record('products', rule='R1', halogen='F', k=1, amount=3)  # More products than attempts!
        
        qa_stats = {
            'version': '2',
            'attempts': 2,
            'products': 3,  # This violates attempts >= products
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should raise AssertionError due to invariant violation
        with self.assertRaises(AssertionError) as cm:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        
        self.assertIn("attempts=2 < products=3", str(cm.exception))
    
    def test_multi_dimensional_consistency(self):
        """Test consistency with multi-rule, multi-halogen data."""
        aggregator = QAAggregator()
        
        # R1_F_1: 2 attempts, 1 product, 1 no_match
        aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=2)
        aggregator.record('products', rule='R1', halogen='F', k=1, amount=1)
        aggregator.record('no_product_matches', rule='R1', halogen='F', k=1, amount=1)
        
        # R1_Cl_1: 3 attempts, 0 products, 3 no_matches  
        aggregator.record('attempts', rule='R1', halogen='Cl', k=1, amount=3)
        aggregator.record('no_product_matches', rule='R1', halogen='Cl', k=1, amount=3)
        
        # R2_F_1: 1 attempt, 1 product, 0 no_match
        aggregator.record('attempts', rule='R2', halogen='F', k=1, amount=1)
        aggregator.record('products', rule='R2', halogen='F', k=1, amount=1)
        
        # Build consistent totals
        qa_stats = {
            'version': '2',
            'attempts': 6,  # 2+3+1
            'products': 2,  # 1+0+1
            'no_product_matches': 4,  # 1+3+0
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0},
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Should pass validation
        try:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
        except AssertionError:
            self.fail("Multi-dimensional consistency check should pass")


if __name__ == '__main__':
    unittest.main()