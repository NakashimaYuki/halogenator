# -*- coding: ascii -*-
"""Tests for reaction guarding and logging functionality."""

import unittest
import logging
from unittest.mock import patch, MagicMock
from rdkit import Chem
from src.halogenator.enumerate_k import _run_reaction_safely, reset_reaction_warning_counts, _max_warnings_per_rule_halogen


class TestReactionGuardAndLogging(unittest.TestCase):
    """Test that reactions are properly guarded and logged."""
    
    def setUp(self):
        """Set up test environment with controlled logging."""
        # Clear warning counts before each test using the provided reset function
        reset_reaction_warning_counts()
        
        # Set up logging to capture messages
        self.log_messages = []
        self.handler = logging.Handler()
        self.handler.emit = lambda record: self.log_messages.append(record.getMessage())
        
        logger = logging.getLogger('src.halogenator.enumerate_k')
        logger.addHandler(self.handler)
        logger.setLevel(logging.WARNING)
    
    def tearDown(self):
        """Clean up logging handler."""
        logger = logging.getLogger('src.halogenator.enumerate_k')
        logger.removeHandler(self.handler)
    
    def test_none_reaction_handling(self):
        """Test that None reactions are handled gracefully."""
        stats_dict = {'template_unsupported': 0}
        
        # Call with None reaction
        result = _run_reaction_safely(None, (MagicMock(),), 'R3', 'F', stats_dict)
        
        # Should return empty list and increment counter
        self.assertEqual(result, [])
        self.assertEqual(stats_dict['template_unsupported'], 1)
        
        # Should log a warning
        self.assertEqual(len(self.log_messages), 1)
        self.assertIn('Rule R3 with F: reaction is None', self.log_messages[0])
    
    def test_reaction_exception_handling(self):
        """Test that reaction exceptions are handled gracefully."""
        stats_dict = {'template_unsupported': 0}
        
        # Create a mock reaction that raises an exception
        mock_rxn = MagicMock()
        mock_rxn.RunReactants.side_effect = Exception("Test reaction failure")
        
        result = _run_reaction_safely(mock_rxn, (MagicMock(),), 'R4', 'Cl', stats_dict)
        
        # Should return empty list and increment counter
        self.assertEqual(result, [])
        self.assertEqual(stats_dict['template_unsupported'], 1)
        
        # Should log a warning
        self.assertEqual(len(self.log_messages), 1)
        self.assertIn('Rule R4 with Cl: reaction failed - Test reaction failure', self.log_messages[0])
    
    def test_successful_reaction_execution(self):
        """Test that successful reactions work normally."""
        stats_dict = {'template_unsupported': 0}
        
        # Create a mock reaction that succeeds
        mock_rxn = MagicMock()
        mock_products = [['product1'], ['product2']]
        mock_rxn.RunReactants.return_value = mock_products
        
        result = _run_reaction_safely(mock_rxn, (MagicMock(),), 'R1', 'Br', stats_dict)
        
        # Should return the products and not increment counter
        self.assertEqual(result, mock_products)
        self.assertEqual(stats_dict['template_unsupported'], 0)
        
        # Should not log any warnings
        self.assertEqual(len(self.log_messages), 0)
    
    def test_warning_deduplication_per_rule_halogen(self):
        """Test that repetitive warnings are deduplicated per rule/halogen combination."""
        stats_dict = {'template_unsupported': 0}
        
        # Call the same rule/halogen combination multiple times with failures
        for i in range(5):  # More than _max_warnings_per_rule_halogen
            _run_reaction_safely(None, (MagicMock(),), 'R5', 'I', stats_dict)
        
        # Should only log one warning despite multiple failures
        self.assertEqual(len(self.log_messages), _max_warnings_per_rule_halogen)
        self.assertIn('Rule R5 with I: reaction is None', self.log_messages[0])
        
        # But should still increment counters for all failures
        self.assertEqual(stats_dict['template_unsupported'], 5)
    
    def test_different_rule_halogen_combinations_separate_quotas(self):
        """Test that different rule/halogen combinations get separate warning quotas."""
        stats_dict = {'template_unsupported': 0}
        
        # Test different combinations - each should get its own warning
        combinations = [('R1', 'F'), ('R1', 'Cl'), ('R3', 'F')]
        
        for rule, halogen in combinations:
            _run_reaction_safely(None, (MagicMock(),), rule, halogen, stats_dict)
        
        # Should log one warning for each unique combination
        self.assertEqual(len(self.log_messages), len(combinations))
        self.assertEqual(stats_dict['template_unsupported'], len(combinations))
        
        # Verify warning content includes correct rule/halogen info
        for i, (rule, halogen) in enumerate(combinations):
            self.assertIn(f'Rule {rule} with {halogen}', self.log_messages[i])
    
    def test_no_stats_dict_provided(self):
        """Test that function works when no stats dict is provided."""
        # Should not crash when stats_dict is None
        result = _run_reaction_safely(None, (MagicMock(),), 'R2', 'Br', None)
        
        self.assertEqual(result, [])
        # Should still log warning
        self.assertEqual(len(self.log_messages), 1)
        self.assertIn('Rule R2 with Br: reaction is None', self.log_messages[0])
    
    def test_reset_function_clears_warning_counts(self):
        """Test that reset_reaction_warning_counts() properly clears the global state."""
        stats_dict = {'template_unsupported': 0}
        
        # Generate some warnings
        _run_reaction_safely(None, (MagicMock(),), 'R1', 'F', stats_dict)
        _run_reaction_safely(None, (MagicMock(),), 'R1', 'F', stats_dict)
        
        # Should only get one warning due to dedup
        self.assertEqual(len(self.log_messages), 1)
        
        # Reset and try again
        reset_reaction_warning_counts()
        self.log_messages.clear()  # Clear captured logs
        
        # Should get another warning since reset cleared the counts
        _run_reaction_safely(None, (MagicMock(),), 'R1', 'F', stats_dict)
        self.assertEqual(len(self.log_messages), 1)
        self.assertIn('Rule R1 with F: reaction is None', self.log_messages[0])
    
    def test_mixed_success_failure_scenarios(self):
        """Test mixed scenarios with both successful and failing reactions."""
        stats_dict = {'template_unsupported': 0}
        
        # Successful reaction
        mock_rxn_success = MagicMock()
        mock_rxn_success.RunReactants.return_value = [['product']]
        result1 = _run_reaction_safely(mock_rxn_success, (MagicMock(),), 'R1', 'F', stats_dict)
        
        # Failed reaction (None)
        result2 = _run_reaction_safely(None, (MagicMock(),), 'R2', 'Cl', stats_dict)
        
        # Failed reaction (exception)
        mock_rxn_fail = MagicMock()
        mock_rxn_fail.RunReactants.side_effect = RuntimeError("Test failure")
        result3 = _run_reaction_safely(mock_rxn_fail, (MagicMock(),), 'R3', 'Br', stats_dict)
        
        # Verify results
        self.assertEqual(result1, [['product']])
        self.assertEqual(result2, [])
        self.assertEqual(result3, [])
        
        # Should have 2 failures in stats
        self.assertEqual(stats_dict['template_unsupported'], 2)
        
        # Should have 2 warning messages
        self.assertEqual(len(self.log_messages), 2)
        self.assertIn('Rule R2 with Cl: reaction is None', self.log_messages[0])
        self.assertIn('Rule R3 with Br: reaction failed - Test failure', self.log_messages[1])
    
    def test_guard_prevents_runreactants_exceptions_from_propagating(self):
        """Test that RunReactants exceptions don't crash the calling code."""
        stats_dict = {'template_unsupported': 0}
        
        # Create a reaction that raises various exception types
        exception_types = [ValueError("bad value"), RuntimeError("runtime error"), TypeError("type error")]
        
        for i, exc in enumerate(exception_types):
            mock_rxn = MagicMock()
            mock_rxn.RunReactants.side_effect = exc
            
            # This should not raise - should be caught and handled
            result = _run_reaction_safely(mock_rxn, (MagicMock(),), f'R{i+1}', 'F', stats_dict)
            
            self.assertEqual(result, [])
        
        # Should have incremented for each failure
        self.assertEqual(stats_dict['template_unsupported'], len(exception_types))
        
        # Should have logged each failure
        self.assertEqual(len(self.log_messages), len(exception_types))


if __name__ == '__main__':
    unittest.main()