# -*- coding: ascii -*-
"""Tests for QA paths dimensionality in pivot reporting."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator


class TestQAPathsDimensionality(unittest.TestCase):
    """Test that qa_paths events are properly recorded with halogen dimension."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
    
    def test_qa_paths_events_have_halogen_dimension(self):
        """Test that qa_paths events are recorded with rule+halogen+k dimensions."""
        # Simulate different qa_paths events for different halogens under same rule/k
        qa_events_f = {
            'isotope_unavailable': 2,
            'atommap_used': 1
        }
        qa_events_cl = {
            'isotope_miss': 1,
            'heuristic_used': 3
        }
        
        # Record attempts with different halogen-specific QA events
        self.aggregator.record_attempt_result('R3', 'F', 1, produced_count=1, qa_events_dict=qa_events_f)
        self.aggregator.record_attempt_result('R3', 'Cl', 1, produced_count=0, qa_events_dict=qa_events_cl)
        
        # Get pivots
        pivots = self.aggregator.to_pivots_dict()
        
        # Verify events are separated by halogen in by_rule_halogen_k dimension
        f_events = pivots['by_rule_halogen_k']['R3_F_1']
        cl_events = pivots['by_rule_halogen_k']['R3_Cl_1']
        
        # F events should only contain F-specific qa events
        self.assertEqual(f_events['isotope_unavailable'], 2)
        self.assertEqual(f_events['atommap_used'], 1)
        self.assertEqual(f_events.get('isotope_miss', 0), 0)  # Should not have Cl events
        self.assertEqual(f_events.get('heuristic_used', 0), 0)  # Should not have Cl events
        
        # Cl events should only contain Cl-specific qa events
        self.assertEqual(cl_events['isotope_miss'], 1)
        self.assertEqual(cl_events['heuristic_used'], 3)
        self.assertEqual(cl_events.get('isotope_unavailable', 0), 0)  # Should not have F events
        self.assertEqual(cl_events.get('atommap_used', 0), 0)  # Should not have F events
    
    def test_pivots_sum_equals_totals_for_qa_paths(self):
        """Test that sum over pivots equals totals for each qa_path key."""
        # Record various events across different rule/halogen/k combinations
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 2, 'atommap_used': 1})
        self.aggregator.record_attempt_result('R1', 'Cl', 1, 0, {'isotope_miss': 1})
        self.aggregator.record_attempt_result('R3', 'F', 2, 1, {'heuristic_used': 3, 'isotope_unavailable': 1})
        self.aggregator.record_attempt_result('R3', 'Cl', 2, 0, {'isotope_miss': 2, 'atommap_used': 1})
        
        # Get pivots
        pivots = self.aggregator.to_pivots_dict()
        
        # Calculate sums across all rule_halogen_k combinations
        isotope_unavailable_sum = 0
        isotope_miss_sum = 0
        atommap_used_sum = 0
        heuristic_used_sum = 0
        
        for key, events in pivots['by_rule_halogen_k'].items():
            isotope_unavailable_sum += events.get('isotope_unavailable', 0)
            isotope_miss_sum += events.get('isotope_miss', 0)
            atommap_used_sum += events.get('atommap_used', 0)
            heuristic_used_sum += events.get('heuristic_used', 0)
        
        # Expected totals based on what we recorded
        expected_isotope_unavailable = 2 + 1  # R1_F_1 + R3_F_2
        expected_isotope_miss = 1 + 2  # R1_Cl_1 + R3_Cl_2
        expected_atommap_used = 1 + 1  # R1_F_1 + R3_Cl_2
        expected_heuristic_used = 3  # R3_F_2
        
        # Verify sums match expected totals
        self.assertEqual(isotope_unavailable_sum, expected_isotope_unavailable)
        self.assertEqual(isotope_miss_sum, expected_isotope_miss)
        self.assertEqual(atommap_used_sum, expected_atommap_used)
        self.assertEqual(heuristic_used_sum, expected_heuristic_used)
    
    def test_by_halogen_aggregation_correct(self):
        """Test that by_halogen dimension correctly aggregates across rules and k."""
        # Record events for same halogen across different rules and k values
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 2})
        self.aggregator.record_attempt_result('R3', 'F', 1, 1, {'isotope_unavailable': 1, 'atommap_used': 3})
        self.aggregator.record_attempt_result('R3', 'F', 2, 0, {'isotope_miss': 4})
        
        # Record events for different halogen
        self.aggregator.record_attempt_result('R1', 'Cl', 1, 1, {'heuristic_used': 5})
        self.aggregator.record_attempt_result('R3', 'Cl', 2, 0, {'isotope_miss': 2})
        
        pivots = self.aggregator.to_pivots_dict()
        
        # Check F aggregation (should sum across all R1,R3 and k=1,2)
        f_events = pivots['by_halogen']['F']
        self.assertEqual(f_events['isotope_unavailable'], 3)  # 2 + 1
        self.assertEqual(f_events['atommap_used'], 3)  # 0 + 3 + 0
        self.assertEqual(f_events['isotope_miss'], 4)  # 0 + 0 + 4
        self.assertEqual(f_events.get('heuristic_used', 0), 0)  # Should not include Cl events
        
        # Check Cl aggregation
        cl_events = pivots['by_halogen']['Cl']
        self.assertEqual(cl_events['heuristic_used'], 5)  # 5 + 0
        self.assertEqual(cl_events['isotope_miss'], 2)  # 0 + 2
        self.assertEqual(cl_events.get('isotope_unavailable', 0), 0)  # Should not include F events
        self.assertEqual(cl_events.get('atommap_used', 0), 0)  # Should not include F events
    
    def test_mixed_qa_events_in_single_attempt(self):
        """Test that multiple qa_paths events in a single attempt are recorded correctly."""
        # Single attempt with multiple qa events
        qa_events = {
            'isotope_unavailable': 1,
            'isotope_miss': 2,
            'atommap_used': 1,
            'heuristic_used': 3
        }
        
        self.aggregator.record_attempt_result('R1', 'F', 1, 1, qa_events)
        
        pivots = self.aggregator.to_pivots_dict()
        events = pivots['by_rule_halogen_k']['R1_F_1']
        
        # All events should be recorded for this specific rule-halogen-k combination
        self.assertEqual(events['isotope_unavailable'], 1)
        self.assertEqual(events['isotope_miss'], 2)
        self.assertEqual(events['atommap_used'], 1)
        self.assertEqual(events['heuristic_used'], 3)
        
        # Also verify they appear in aggregated dimensions
        self.assertEqual(pivots['by_rule']['R1']['isotope_unavailable'], 1)
        self.assertEqual(pivots['by_halogen']['F']['isotope_miss'], 2)
        self.assertEqual(pivots['by_k'][1]['atommap_used'], 1)
        self.assertEqual(pivots['by_rule_halogen']['R1_F']['heuristic_used'], 3)
    
    def test_atommap_heuristic_binary_semantics(self):
        """Test that atommap_used/heuristic_used use binary per-attempt semantics."""
        # Simulate attempt where multiple products use atommap
        # In product-level semantics: might count 2-3 times
        # In attempt-level semantics: should count exactly 1
        
        # This simulates the normalization that happens during product QA event merging
        # Multiple products using atommap in one attempt -> attempt records 1 not sum
        qa_events = {
            'atommap_used': 1,  # Binary: was atommap used in this attempt? Yes=1
            'heuristic_used': 1  # Binary: was heuristic used in this attempt? Yes=1
        }
        
        self.aggregator.record_attempt_result('R3', 'F', 1, 1, qa_events)
        
        pivots = self.aggregator.to_pivots_dict()
        events = pivots['by_rule_halogen_k']['R3_F_1']
        
        # Should be exactly 1 for binary semantics
        self.assertEqual(events['atommap_used'], 1)
        self.assertEqual(events['heuristic_used'], 1)
        
        # Even if we record another attempt with same rule-halogen-k, each attempt is independent
        self.aggregator.record_attempt_result('R3', 'F', 1, 0, {'atommap_used': 1})
        
        updated_pivots = self.aggregator.to_pivots_dict()
        updated_events = updated_pivots['by_rule_halogen_k']['R3_F_1']
        
        # Now should have 2 (one per attempt)
        self.assertEqual(updated_events['atommap_used'], 2)
        self.assertEqual(updated_events['heuristic_used'], 1)  # Only first attempt used heuristic


if __name__ == '__main__':
    unittest.main()