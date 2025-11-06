# -*- coding: ascii -*-
"""Tests for totals-pivots consistency validation."""

import unittest
import logging
from unittest.mock import patch
from src.halogenator.enumerate_k import QAAggregator, _validate_totals_pivots_consistency


class TestTotalsPivotsConsistency(unittest.TestCase):
    """Test that totals-pivots consistency validation works correctly."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
        # Capture log messages for testing
        self.log_messages = []
        self.log_handler = logging.Handler()
        self.log_handler.emit = lambda record: self.log_messages.append(record.getMessage())
        
        # Add handler to the enumerate_k logger
        enumerate_k_logger = logging.getLogger('src.halogenator.enumerate_k')
        enumerate_k_logger.addHandler(self.log_handler)
        enumerate_k_logger.setLevel(logging.DEBUG)
    
    def tearDown(self):
        """Clean up logging handler."""
        enumerate_k_logger = logging.getLogger('src.halogenator.enumerate_k')
        enumerate_k_logger.removeHandler(self.log_handler)
    
    def test_consistent_qa_paths_validation_passes(self):
        """Test that consistent qa_paths totals and pivots pass validation."""
        # Record events that should be consistent
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 2, 'atommap_used': 1})
        self.aggregator.record_attempt_result('R3', 'Cl', 1, 0, {'isotope_miss': 1, 'heuristic_used': 3})
        
        # Create qa_stats dict with matching totals
        qa_stats = {
            'no_product_matches': 1,  # One failed attempt (R3_Cl_1)
            'qa_paths': {
                'isotope_unavailable': 2,
                'isotope_miss': 1,
                'atommap_used': 1,
                'heuristic_used': 3
            }
        }
        
        # Should not raise any exceptions or log errors
        _validate_totals_pivots_consistency(qa_stats, self.aggregator)
        
        # Check that debug message was logged
        debug_messages = [msg for msg in self.log_messages if 'consistency validation passed' in msg]
        self.assertEqual(len(debug_messages), 1)
    
    def test_inconsistent_qa_paths_validation_fails(self):
        """Test that inconsistent qa_paths totals trigger validation errors."""
        # Record events
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 2})
        self.aggregator.record_attempt_result('R3', 'Cl', 1, 0, {'isotope_miss': 1})
        
        # Create qa_stats dict with mismatched totals
        qa_stats = {
            'no_product_matches': 1,
            'qa_paths': {
                'isotope_unavailable': 5,  # Should be 2 (mismatch!)
                'isotope_miss': 1,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        # Should log error about mismatch
        _validate_totals_pivots_consistency(qa_stats, self.aggregator)
        
        # Check that error message was logged
        error_messages = [msg for msg in self.log_messages if 'consistency violations' in msg]
        self.assertEqual(len(error_messages), 1)
        self.assertIn('isotope_unavailable: total=5 != pivot_sum=2', error_messages[0])
    
    def test_no_product_matches_consistency(self):
        """Test that no_product_matches totals are validated against pivots."""
        # Record one successful and one failed attempt
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {})  # Success
        self.aggregator.record_attempt_result('R3', 'Cl', 1, 0, {})  # Failure
        
        # Create qa_stats with mismatched no_product_matches
        qa_stats = {
            'no_product_matches': 3,  # Should be 1 (mismatch!)
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        _validate_totals_pivots_consistency(qa_stats, self.aggregator)
        
        # Check that error about no_product_matches mismatch was logged
        error_messages = [msg for msg in self.log_messages if 'consistency violations' in msg]
        self.assertEqual(len(error_messages), 1)
        self.assertIn('no_product_matches: total=3 != pivot_sum=1', error_messages[0])
    
    def test_attempt_invariant_validation(self):
        """Test that attempt invariants are validated."""
        # Record events that violate invariants
        self.aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=1)
        self.aggregator.record('products', rule='R1', halogen='F', k=1, amount=2)  # More products than attempts!
        
        qa_stats = {
            'no_product_matches': 0,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        _validate_totals_pivots_consistency(qa_stats, self.aggregator)
        
        # Check that invariant violation was logged
        error_messages = [msg for msg in self.log_messages if 'consistency violations' in msg]
        self.assertEqual(len(error_messages), 1)
        self.assertIn('attempts=1 < products=2 (invariant violation)', error_messages[0])
    
    def test_debug_assertion_enabled_new_var(self):
        """Test that assertion is raised when HALO_ASSERT_PIVOT_CONSISTENCY is set."""
        with patch.dict('os.environ', {'HALO_ASSERT_PIVOT_CONSISTENCY': '1'}):
            # Create mismatch scenario
            self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 1})
            
            qa_stats = {
                'no_product_matches': 0,
                'qa_paths': {
                    'isotope_unavailable': 5,  # Mismatch
                    'isotope_miss': 0,
                    'atommap_used': 0,
                    'heuristic_used': 0
                }
            }
            
            # Should raise AssertionError in debug mode
            with self.assertRaises(AssertionError) as cm:
                _validate_totals_pivots_consistency(qa_stats, self.aggregator)
            
            self.assertIn('isotope_unavailable: total=5 != pivot_sum=1', str(cm.exception))
    
    def test_debug_assertion_enabled_backward_compat(self):
        """Test that assertion is raised when HALOGENATOR_DEBUG is set (backward compatibility)."""
        with patch.dict('os.environ', {'HALOGENATOR_DEBUG': '1'}):
            # Create mismatch scenario
            self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 1})
            
            qa_stats = {
                'no_product_matches': 0,
                'qa_paths': {
                    'isotope_unavailable': 5,  # Mismatch
                    'isotope_miss': 0,
                    'atommap_used': 0,
                    'heuristic_used': 0
                }
            }
            
            # Should raise AssertionError in debug mode
            with self.assertRaises(AssertionError) as cm:
                _validate_totals_pivots_consistency(qa_stats, self.aggregator)
            
            self.assertIn('isotope_unavailable: total=5 != pivot_sum=1', str(cm.exception))
    
    def test_no_assertion_when_neither_env_var_set(self):
        """Test that no assertion is raised when neither env var is set."""
        with patch.dict('os.environ', {}, clear=True):  # Clear all env vars
            # Create mismatch scenario
            self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 1})
            
            qa_stats = {
                'no_product_matches': 0,
                'qa_paths': {
                    'isotope_unavailable': 5,  # Mismatch
                    'isotope_miss': 0,
                    'atommap_used': 0,
                    'heuristic_used': 0
                }
            }
            
            # Should NOT raise AssertionError when no debug env vars set
            _validate_totals_pivots_consistency(qa_stats, self.aggregator)  # Should not raise
    
    def test_validation_with_empty_data(self):
        """Test validation with empty qa_paths or pivots."""
        # Empty qa_paths
        qa_stats_empty = {}
        _validate_totals_pivots_consistency(qa_stats_empty, self.aggregator)  # Should not crash
        
        # Empty pivots
        qa_stats_with_paths = {
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        empty_aggregator = QAAggregator()
        _validate_totals_pivots_consistency(qa_stats_with_paths, empty_aggregator)  # Should not crash
    
    def test_multiple_rule_halogen_k_combinations(self):
        """Test validation across multiple rule-halogen-k combinations."""
        # Record events across different combinations
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'atommap_used': 1})
        self.aggregator.record_attempt_result('R1', 'Cl', 1, 0, {'isotope_miss': 1})
        self.aggregator.record_attempt_result('R3', 'F', 2, 1, {'heuristic_used': 2})
        self.aggregator.record_attempt_result('R3', 'Cl', 2, 0, {'isotope_unavailable': 3})
        
        # Create consistent qa_stats
        qa_stats = {
            'no_product_matches': 2,  # Two failed attempts
            'qa_paths': {
                'isotope_unavailable': 3,
                'isotope_miss': 1,
                'atommap_used': 1,
                'heuristic_used': 2
            }
        }
        
        _validate_totals_pivots_consistency(qa_stats, self.aggregator)
        
        # Should pass validation
        debug_messages = [msg for msg in self.log_messages if 'consistency validation passed' in msg]
        self.assertEqual(len(debug_messages), 1)


if __name__ == '__main__':
    unittest.main()