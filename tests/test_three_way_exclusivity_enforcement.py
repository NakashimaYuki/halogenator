# -*- coding: ascii -*-
"""Tests to ensure three-way exclusivity is enforced and prevent regression."""

import unittest
from src.halogenator.enumerate_k import _validate_totals_pivots_consistency, QAAggregator


class TestThreeWayExclusivityEnforcement(unittest.TestCase):
    """Test that three-way exclusivity violations are properly detected."""
    
    def test_mutually_exclusive_violation_detected(self):
        """Test that overlapping no_product_matches and template_unsupported triggers validation error."""
        # Create QAAggregator with violating data (same attempt counted in both categories)
        aggregator = QAAggregator()
        
        # Simulate a single attempt that incorrectly gets counted in both categories
        # This should NEVER happen in real code, but we test that validation catches it
        # NOTE: record_attempt_result automatically adds no_product_matches=1 when produced_count=0
        # So adding template_unsupported=1 in qa_events creates the violation
        aggregator.record_attempt_result(
            'R1', 'F', 1,
            produced_count=0,  # This causes automatic no_product_matches=1
            qa_events_dict={
                'template_unsupported': 1  # VIOLATION: same attempt in both categories
            }
        )
        
        # Build QA stats dict
        qa_stats = {
            'attempts': 1,
            'products': 0,
            'no_product_matches': 1,
            'template_unsupported': 1,  # This creates: 1 != 0 + 1 + 1 (attempts != products + no_match + template_unsupported)
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 0,
                'heuristic_used': 0
            }
        }
        
        # Validation should detect the violation and log errors
        import os
        os.environ['HALO_ASSERT_PIVOT_CONSISTENCY'] = '1'
        
        try:
            with self.assertRaises(AssertionError) as context:
                _validate_totals_pivots_consistency(qa_stats, aggregator)
            
            # Check that the error message mentions the invariant violation
            error_msg = str(context.exception)
            self.assertIn('invariant violation', error_msg.lower())
            # Don't check exact numbers since they depend on implementation details
            # Just verify that it's detecting the three-way invariant violation
            self.assertIn('attempts=', error_msg)
            self.assertIn('products=', error_msg)
            self.assertIn('template_unsupported=', error_msg)
            
        finally:
            # Clean up environment
            if 'HALO_ASSERT_PIVOT_CONSISTENCY' in os.environ:
                del os.environ['HALO_ASSERT_PIVOT_CONSISTENCY']
    
    def test_correct_three_way_exclusivity_passes(self):
        """Test that correct three-way exclusive categorization passes validation."""
        aggregator = QAAggregator()
        
        # Three separate attempts with mutually exclusive outcomes
        # Attempt 1: Success
        aggregator.record_attempt_result('R1', 'F', 1, produced_count=1, qa_events_dict={'atommap_used': 1})
        
        # Attempt 2: Failed but template supported  
        aggregator.record_attempt_result('R1', 'Cl', 1, produced_count=0, qa_events_dict={})
        
        # Attempt 3: Template unsupported
        aggregator.record_attempt_result('R3', 'Br', 1, produced_count=0, qa_events_dict={'template_unsupported': 1})
        
        qa_stats = {
            'attempts': 3,
            'products': 1,
            'no_product_matches': 1, 
            'template_unsupported': 1,
            'qa_paths': {
                'isotope_unavailable': 0,
                'isotope_miss': 0,
                'atommap_used': 1,
                'heuristic_used': 0
            }
        }
        
        # This should pass validation (3 = 1 + 1 + 1)
        try:
            _validate_totals_pivots_consistency(qa_stats, aggregator)
            # If we reach here, validation passed (no exception raised)
        except Exception as e:
            self.fail(f"Valid three-way exclusive data should not trigger validation error: {e}")
    
    def test_attempts_less_than_sum_detected(self):
        """Test that attempts < (products + no_matches + template_unsupported) is detected."""
        aggregator = QAAggregator()
        
        # Record 2 attempts in pivots but claim only 1 attempt in totals
        aggregator.record_attempt_result('R1', 'F', 1, produced_count=1, qa_events_dict={})
        aggregator.record_attempt_result('R1', 'Cl', 1, produced_count=0, qa_events_dict={})
        
        qa_stats = {
            'attempts': 1,  # VIOLATION: Should be 2 based on pivot data
            'products': 1,
            'no_product_matches': 1,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0}
        }
        
        import os
        os.environ['HALO_ASSERT_PIVOT_CONSISTENCY'] = '1'
        
        try:
            with self.assertRaises(AssertionError):
                _validate_totals_pivots_consistency(qa_stats, aggregator)
        finally:
            if 'HALO_ASSERT_PIVOT_CONSISTENCY' in os.environ:
                del os.environ['HALO_ASSERT_PIVOT_CONSISTENCY']
    
    def test_missing_template_unsupported_in_totals(self):
        """Test that missing template_unsupported from totals but present in pivots is detected."""
        aggregator = QAAggregator()
        
        aggregator.record_attempt_result('R3', 'F', 1, produced_count=0, qa_events_dict={'template_unsupported': 1})
        
        # qa_stats missing template_unsupported field (old format)
        qa_stats = {
            'attempts': 1,
            'products': 0, 
            'no_product_matches': 0,
            # 'template_unsupported': 1,  # Missing this field
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0}
        }
        
        import os
        os.environ['HALO_ASSERT_PIVOT_CONSISTENCY'] = '1'
        
        try:
            with self.assertRaises(AssertionError) as context:
                _validate_totals_pivots_consistency(qa_stats, aggregator)
            
            # Should detect template_unsupported: total=0 != pivot_sum=1 mismatch
            error_msg = str(context.exception)
            self.assertIn('template_unsupported', error_msg)
            
        finally:
            if 'HALO_ASSERT_PIVOT_CONSISTENCY' in os.environ:
                del os.environ['HALO_ASSERT_PIVOT_CONSISTENCY']


if __name__ == '__main__':
    unittest.main()