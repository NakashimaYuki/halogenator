# -*- coding: ascii -*-
"""Tests for QAAggregator amount parameter functionality."""

import unittest
from src.halogenator.enumerate_k import QAAggregator


class TestAggregatorAmount(unittest.TestCase):
    """Test that QAAggregator handles amount parameter correctly."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
    
    def test_record_with_amount_parameter(self):
        """Test that record(amount=7) increments exactly 7 without loops."""
        # Record an event with amount=7
        self.aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=7)
        
        # Verify all dimensions got exactly 7
        self.assertEqual(self.aggregator.by_rule['R1']['attempts'], 7)
        self.assertEqual(self.aggregator.by_halogen['F']['attempts'], 7)
        self.assertEqual(self.aggregator.by_k[1]['attempts'], 7)
        self.assertEqual(self.aggregator.by_rule_halogen['R1_F']['attempts'], 7)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['attempts'], 7)
    
    def test_record_with_amount_zero(self):
        """Test that record(amount=0) does not increment."""
        # Record an event with amount=0
        self.aggregator.record('products', rule='R2', halogen='Cl', k=2, amount=0)
        
        # Verify no dimensions were incremented
        self.assertEqual(self.aggregator.by_rule['R2'].get('products', 0), 0)
        self.assertEqual(self.aggregator.by_halogen['Cl'].get('products', 0), 0)
        self.assertEqual(self.aggregator.by_k[2].get('products', 0), 0)
        self.assertEqual(self.aggregator.by_rule_halogen['R2_Cl'].get('products', 0), 0)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R2_Cl_2'].get('products', 0), 0)
    
    def test_record_default_amount_is_one(self):
        """Test that record() without amount parameter defaults to 1."""
        # Record an event without amount parameter
        self.aggregator.record('template_unsupported', rule='R3', halogen='Br', k=1)
        
        # Verify all dimensions got exactly 1
        self.assertEqual(self.aggregator.by_rule['R3']['template_unsupported'], 1)
        self.assertEqual(self.aggregator.by_halogen['Br']['template_unsupported'], 1)
        self.assertEqual(self.aggregator.by_k[1]['template_unsupported'], 1)
        self.assertEqual(self.aggregator.by_rule_halogen['R3_Br']['template_unsupported'], 1)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R3_Br_1']['template_unsupported'], 1)
    
    def test_record_attempt_result_helper(self):
        """Test the record_attempt_result helper method."""
        # Successful attempt with 2 products and some QA events
        qa_events = {
            'isotope_unavailable': 1,
            'atommap_used': 2
        }
        self.aggregator.record_attempt_result('R1', 'F', 1, produced_count=2, qa_events_dict=qa_events)
        
        # Verify attempt semantics
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['attempts'], 1)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['products'], 1)  # 1 because produced_count > 0
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1'].get('no_product_matches', 0), 0)  # Should not be set
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['isotope_unavailable'], 1)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['atommap_used'], 2)
        
        # Failed attempt with 0 products
        qa_events_fail = {
            'template_unsupported': 1
        }
        self.aggregator.record_attempt_result('R2', 'Cl', 1, produced_count=0, qa_events_dict=qa_events_fail)
        
        # Verify failed attempt semantics
        self.assertEqual(self.aggregator.by_rule_halogen_k['R2_Cl_1']['attempts'], 1)
        self.assertEqual(self.aggregator.by_rule_halogen_k['R2_Cl_1'].get('products', 0), 0)  # Should not be set
        self.assertEqual(self.aggregator.by_rule_halogen_k['R2_Cl_1']['no_product_matches'], 1)  # 1 because produced_count == 0
        self.assertEqual(self.aggregator.by_rule_halogen_k['R2_Cl_1']['template_unsupported'], 1)
    
    def test_multiple_amounts_accumulate(self):
        """Test that multiple record calls with amounts accumulate correctly."""
        # Record same event multiple times with different amounts
        self.aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=3)
        self.aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=5)
        self.aggregator.record('attempts', rule='R1', halogen='F', k=1, amount=2)
        
        # Verify total is sum of all amounts
        self.assertEqual(self.aggregator.by_rule_halogen_k['R1_F_1']['attempts'], 10)  # 3 + 5 + 2


if __name__ == '__main__':
    unittest.main()