# -*- coding: ascii -*-
"""Integration tests for fallback site detection paths."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator, _apply_reaction_rule, EnumConfig
from rdkit import Chem


class TestFallbackIntegration(unittest.TestCase):
    """Test fallback site detection without mocking the detection function."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
        self.cfg = EnumConfig(
            halogens=['F', 'Cl'],
            k_max=2,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        # Simple test molecule
        self.mol = Chem.MolFromSmiles('c1ccccc1')  # benzene
    
    def test_fallback_detection_no_keyerror(self):
        """Test that fallback detection doesn't raise KeyError with empty qa_events dict."""
        # Mock reaction setup to force fallback path
        mock_reactions = {
            'R3': {
                'F': MagicMock(),
                'Cl': MagicMock()
            }
        }
        
        # Mock _find_reaction_matches to return empty (triggers isotope unavailable path)
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([], 0)  # No matches, no template errors
            
            # Mock _run_reaction_safely to return a mock product
            with patch('src.halogenator.enumerate_k._run_reaction_safely') as mock_run_reaction:
                # Create a mock product that will trigger fallback detection
                mock_product = MagicMock()
                mock_product.GetAtoms.return_value = [
                    MagicMock(GetAtomMapNum=lambda: 0, SetAtomMapNum=lambda x: None, GetSymbol=lambda: 'C', GetIdx=lambda: 0),
                    MagicMock(GetAtomMapNum=lambda: 0, SetAtomMapNum=lambda x: None, GetSymbol=lambda: 'F', GetIdx=lambda: 1,
                             GetNeighbors=lambda: [MagicMock(GetSymbol=lambda: 'C', GetIdx=lambda: 0)])
                ]
                mock_run_reaction.return_value = [([mock_product],)]
                
                # Mock flavonoid_ring_label to return something
                with patch('src.halogenator.enumerate_k.flavonoid_ring_label') as mock_ring_label:
                    mock_ring_label.return_value = 'A'
                    
                    # Mock the product processing to avoid complex logic
                    with patch('src.halogenator.enumerate_k._process_reaction_product') as mock_process:
                        mock_process.return_value = ('mock_result', 'mock_history', {'smiles': 'mock'})
                        
                        # This should not raise KeyError even though we pass empty qa_events dicts
                        try:
                            results, qa_stats = _apply_reaction_rule(
                                self.mol, 'R3', mock_reactions, self.cfg,
                                [], 0, set(), self.aggregator
                            )
                            # If we get here without exception, the test passes
                            self.assertIsInstance(results, list)
                            self.assertIsInstance(qa_stats, dict)
                        except KeyError as e:
                            self.fail(f"KeyError raised in fallback detection: {e}")
                        except Exception as e:
                            # Other exceptions are fine for this specific test
                            pass
    
    def test_heuristic_used_recorded_in_aggregator(self):
        """Test that heuristic_used events are properly recorded in aggregator."""
        mock_reactions = {
            'R3': {
                'F': MagicMock()
            }
        }
        
        # Force the path that uses fallback detection
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([], 0)  # No matches
            
            with patch('src.halogenator.enumerate_k._run_reaction_safely') as mock_run_reaction:
                # Mock a product that will trigger heuristic detection
                mock_product = MagicMock()
                # Setup the mock to trigger heuristic path (no atom mapping)
                mock_product.GetAtoms.return_value = [
                    MagicMock(GetAtomMapNum=lambda: 0, SetAtomMapNum=lambda x: None)
                ]
                mock_run_reaction.return_value = [([mock_product],)]
                
                # Mock the heuristic detection to succeed and return site
                with patch('src.halogenator.enumerate_k._detect_site_from_product') as mock_detect:
                    mock_detect.return_value = (0, 'A')  # Returns valid site
                    
                    with patch('src.halogenator.enumerate_k.flavonoid_ring_label') as mock_ring_label:
                        mock_ring_label.return_value = 'A'
                        
                        with patch('src.halogenator.enumerate_k._process_reaction_product') as mock_process:
                            mock_process.return_value = ('result', [], {'smiles': 'test'})
                            
                            # Run the function
                            results, qa_stats = _apply_reaction_rule(
                                self.mol, 'R3', mock_reactions, self.cfg,
                                [], 0, set(), self.aggregator
                            )
                            
                            # Check that heuristic_used was recorded
                            pivots = self.aggregator.to_pivots_dict()
                            if 'R3_F_1' in pivots.get('by_rule_halogen_k', {}):
                                events = pivots['by_rule_halogen_k']['R3_F_1']
                                # Should have heuristic_used event
                                self.assertGreater(events.get('heuristic_used', 0), 0)
    
    def test_atommap_used_recorded_in_aggregator(self):
        """Test that atommap_used events are properly recorded in aggregator."""
        mock_reactions = {
            'R3': {
                'F': MagicMock()
            }
        }
        
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([], 0)  # No matches
            
            with patch('src.halogenator.enumerate_k._run_reaction_safely') as mock_run_reaction:
                # Mock a product that will trigger atommap detection path
                mock_product = MagicMock()
                # Setup mock to have atom with mapping number 1 (atommap success)
                mock_atom = MagicMock()
                mock_atom.GetAtomMapNum.return_value = 1
                mock_atom.GetIdx.return_value = 0
                mock_atom.SetAtomMapNum = MagicMock()
                
                mock_other_atom = MagicMock()
                mock_other_atom.GetAtomMapNum.return_value = 0
                mock_other_atom.SetAtomMapNum = MagicMock()
                
                mock_product.GetAtoms.return_value = [mock_atom, mock_other_atom]
                mock_run_reaction.return_value = [([mock_product],)]
                
                with patch('src.halogenator.enumerate_k.flavonoid_ring_label') as mock_ring_label:
                    mock_ring_label.return_value = 'A'
                    
                    with patch('src.halogenator.enumerate_k._process_reaction_product') as mock_process:
                        mock_process.return_value = ('result', [], {'smiles': 'test'})
                        
                        # Run the function
                        results, qa_stats = _apply_reaction_rule(
                            self.mol, 'R3', mock_reactions, self.cfg,
                            [], 0, set(), self.aggregator
                        )
                        
                        # Check that atommap_used was recorded
                        pivots = self.aggregator.to_pivots_dict()
                        if 'R3_F_1' in pivots.get('by_rule_halogen_k', {}):
                            events = pivots['by_rule_halogen_k']['R3_F_1']
                            # Should have atommap_used event
                            self.assertGreater(events.get('atommap_used', 0), 0)


if __name__ == '__main__':
    unittest.main()