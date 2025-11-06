# -*- coding: ascii -*-
"""Tests for attempt boundary semantics in enumeration."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator, _apply_reaction_rule, EnumConfig


class TestAttemptBoundarySemantics(unittest.TestCase):
    """Test that attempts, products, and no_product_matches follow correct semantics."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
        
        # Create a minimal EnumConfig
        self.cfg = EnumConfig(
            halogens=['F', 'Cl'],
            k_max=2,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
    
    def test_attempt_boundary_zero_products(self):
        """Test attempt with zero products increments no_product_matches."""
        # Mock a scenario where the reaction fails to produce any products
        mock_mol = MagicMock()
        mock_reactions = {
            'R3': {
                'F': MagicMock(),
                'Cl': MagicMock()
            }
        }
        
        # Mock _find_reaction_matches to return no matches
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([], 0)  # No matches, no template_unsupported
            
            # Mock _run_reaction_safely to return no products
            with patch('src.halogenator.enumerate_k._run_reaction_safely') as mock_run_safely:
                mock_run_safely.return_value = []  # No products
                
                # Call _apply_reaction_rule
                results, qa_stats = _apply_reaction_rule(
                    mock_mol, 'R3', mock_reactions, self.cfg, [], 0, set(), self.aggregator
                )
                
                # Verify results
                self.assertEqual(len(results), 0)
                
                # Check aggregator - should have 2 attempts (F and Cl), 0 products, 2 no_product_matches
                pivots = self.aggregator.to_pivots_dict()
                
                # F attempt
                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['attempts'], 1)
                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1'].get('products', 0), 0)
                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['no_product_matches'], 1)
                
                # Cl attempt  
                self.assertEqual(pivots['by_rule_halogen_k']['R3_Cl_1']['attempts'], 1)
                self.assertEqual(pivots['by_rule_halogen_k']['R3_Cl_1'].get('products', 0), 0)
                self.assertEqual(pivots['by_rule_halogen_k']['R3_Cl_1']['no_product_matches'], 1)
    
    def test_attempt_boundary_one_product(self):
        """Test attempt with one product increments products once."""
        mock_mol = MagicMock()
        mock_reactions = {
            'R3': {
                'F': MagicMock()
            }
        }
        
        # Mock _find_reaction_matches to return one match
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([(0, MagicMock())], 0)  # One match
            
            # Mock flavonoid_ring_label
            with patch('src.halogenator.enumerate_k.flavonoid_ring_label') as mock_ring_label:
                mock_ring_label.return_value = 'A'
                
                # Mock the reaction execution to produce one product
                with patch('src.halogenator.enumerate_k.Chem.RWMol') as mock_rwmol:
                    mock_mol_copy = MagicMock()
                    mock_rwmol.return_value = mock_mol_copy
                    mock_atom = MagicMock()
                    mock_mol_copy.GetAtomWithIdx.return_value = mock_atom
                    
                    mock_rxn = MagicMock()
                    mock_reactions['R3']['F'] = mock_rxn
                    
                    # Mock reaction to return one product
                    mock_product = MagicMock()
                    mock_rxn.RunReactants.return_value = [[(mock_product,)]]
                    
                    with patch('src.halogenator.enumerate_k._find_isotope_tagged_site') as mock_find_site:
                        mock_find_site.return_value = 0  # Found tagged site
                        
                        with patch('src.halogenator.enumerate_k._clear_isotope_tags'):
                            with patch('src.halogenator.enumerate_k._process_reaction_product') as mock_process:
                                mock_final_product = (MagicMock(), [], {'smiles': 'test'})
                                mock_process.return_value = mock_final_product
                                
                                # Call _apply_reaction_rule
                                results, qa_stats = _apply_reaction_rule(
                                    mock_mol, 'R3', mock_reactions, self.cfg, [], 0, set(), self.aggregator
                                )
                                
                                # Verify results
                                self.assertEqual(len(results), 1)
                                
                                # Check aggregator - should have 1 attempt, 1 product, 0 no_product_matches
                                pivots = self.aggregator.to_pivots_dict()
                                
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['attempts'], 1)
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['products'], 1)
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1'].get('no_product_matches', 0), 0)
    
    def test_attempt_boundary_multiple_products(self):
        """Test attempt with multiple products still increments products only once."""
        mock_mol = MagicMock()
        mock_reactions = {
            'R3': {
                'F': MagicMock()
            }
        }
        
        # Mock _find_reaction_matches to return two matches (two sites)
        with patch('src.halogenator.enumerate_k._find_reaction_matches') as mock_find_matches:
            mock_find_matches.return_value = ([(0, MagicMock()), (1, MagicMock())], 0)  # Two matches
            
            # Mock flavonoid_ring_label
            with patch('src.halogenator.enumerate_k.flavonoid_ring_label') as mock_ring_label:
                mock_ring_label.return_value = 'A'
                
                # Mock the reaction execution to produce products for both sites
                with patch('src.halogenator.enumerate_k.Chem.RWMol') as mock_rwmol:
                    mock_mol_copy = MagicMock()
                    mock_rwmol.return_value = mock_mol_copy
                    mock_atom = MagicMock()
                    mock_mol_copy.GetAtomWithIdx.return_value = mock_atom
                    
                    mock_rxn = MagicMock()
                    mock_reactions['R3']['F'] = mock_rxn
                    
                    # Mock reaction to return one product each time
                    mock_product1 = MagicMock()
                    mock_product2 = MagicMock()
                    mock_rxn.RunReactants.side_effect = [
                        [[(mock_product1,)]],  # First site produces one product
                        [[(mock_product2,)]]   # Second site produces one product
                    ]
                    
                    with patch('src.halogenator.enumerate_k._find_isotope_tagged_site') as mock_find_site:
                        mock_find_site.return_value = 0  # Found tagged site for both
                        
                        with patch('src.halogenator.enumerate_k._clear_isotope_tags'):
                            with patch('src.halogenator.enumerate_k._process_reaction_product') as mock_process:
                                mock_final_product1 = (MagicMock(), [], {'smiles': 'test1'})
                                mock_final_product2 = (MagicMock(), [], {'smiles': 'test2'})
                                mock_process.side_effect = [mock_final_product1, mock_final_product2]
                                
                                # Call _apply_reaction_rule
                                results, qa_stats = _apply_reaction_rule(
                                    mock_mol, 'R3', mock_reactions, self.cfg, [], 0, set(), self.aggregator
                                )
                                
                                # Verify results - should have 2 products from this attempt
                                self.assertEqual(len(results), 2)
                                
                                # Check aggregator - should have 1 attempt, 1 product (not 2!), 0 no_product_matches
                                # This is the key test: multiple products from one attempt still count as 1 "products" event
                                pivots = self.aggregator.to_pivots_dict()
                                
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['attempts'], 1)
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1']['products'], 1)  # Not 2!
                                self.assertEqual(pivots['by_rule_halogen_k']['R3_F_1'].get('no_product_matches', 0), 0)
    
    def test_attempt_consistency_invariant(self):
        """Test that attempts >= products and attempts == products + no_product_matches."""
        # This test uses the aggregator directly to verify the invariant
        self.aggregator.record_attempt_result('R1', 'F', 1, produced_count=1, qa_events_dict={})
        self.aggregator.record_attempt_result('R1', 'Cl', 1, produced_count=0, qa_events_dict={})
        self.aggregator.record_attempt_result('R2', 'F', 1, produced_count=2, qa_events_dict={})
        
        pivots = self.aggregator.to_pivots_dict()
        
        # Check invariants for each rule-halogen-k combination
        for key, events in pivots['by_rule_halogen_k'].items():
            attempts = events.get('attempts', 0)
            products = events.get('products', 0) 
            no_matches = events.get('no_product_matches', 0)
            
            # Invariant 1: attempts >= products
            self.assertGreaterEqual(attempts, products, f"Failed for {key}: attempts={attempts}, products={products}")
            
            # Invariant 2: attempts == products + no_product_matches
            self.assertEqual(attempts, products + no_matches, 
                           f"Failed for {key}: attempts={attempts}, products={products}, no_matches={no_matches}")


if __name__ == '__main__':
    unittest.main()