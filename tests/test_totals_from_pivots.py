# -*- coding: ascii -*-
"""Tests for totals-from-pivots unification."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator, _compute_totals_from_aggregator, enumerate_with_stats, EnumConfig
from src.halogenator.enumerate_k1 import enumerate_k1_with_stats


class TestTotalsFromPivots(unittest.TestCase):
    """Test that totals are correctly computed from aggregator pivots."""
    
    def test_compute_totals_from_aggregator_basic(self):
        """Test basic totals computation from aggregator."""
        aggregator = QAAggregator()
        
        # Record some events
        aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 2, 'atommap_used': 1})
        aggregator.record_attempt_result('R3', 'Cl', 1, 0, {'isotope_miss': 1, 'heuristic_used': 1})
        aggregator.record('template_unsupported', rule='R4', halogen='Br', k=2, amount=3)
        
        # Compute totals
        totals = _compute_totals_from_aggregator(aggregator)
        
        # Validate totals match pivot sums
        self.assertEqual(totals['attempts'], 2)  # Two record_attempt_result calls
        self.assertEqual(totals['products'], 1)  # One successful attempt (R1_F_1)
        self.assertEqual(totals['no_product_matches'], 1)  # One failed attempt (R3_Cl_1)
        self.assertEqual(totals['template_unsupported'], 3)  # R4 template failures
        
        # Validate qa_paths
        self.assertEqual(totals['qa_paths']['isotope_unavailable'], 2)
        self.assertEqual(totals['qa_paths']['isotope_miss'], 1)
        self.assertEqual(totals['qa_paths']['atommap_used'], 1)
        self.assertEqual(totals['qa_paths']['heuristic_used'], 1)
        
        # Validate invariants hold
        self.assertGreaterEqual(totals['attempts'], totals['products'])
        self.assertEqual(totals['attempts'], totals['products'] + totals['no_product_matches'])
    
    def test_compute_totals_empty_aggregator(self):
        """Test totals computation with empty aggregator."""
        aggregator = QAAggregator()
        totals = _compute_totals_from_aggregator(aggregator)
        
        # Should have all zeros
        self.assertEqual(totals['attempts'], 0)
        self.assertEqual(totals['products'], 0)
        self.assertEqual(totals['no_product_matches'], 0)
        self.assertEqual(totals['template_unsupported'], 0)
        
        for key in totals['qa_paths']:
            self.assertEqual(totals['qa_paths'][key], 0)
    
    def test_multi_match_single_successful_attempt(self):
        """Test scenario with multi-match but only one successful attempt (the key edge case)."""
        aggregator = QAAggregator()
        
        # Simulate scenario: one attempt had multiple matches, some failed, but attempt succeeded overall
        # In old system: multiple failed matches would increment no_product_matches multiple times
        # In new system: attempt semantics means only 1 or 0 per attempt
        aggregator.record_attempt_result('R3', 'F', 1, 1, {})  # Overall successful attempt
        
        totals = _compute_totals_from_aggregator(aggregator)
        
        # Should show attempt-level semantics (not per-match)
        self.assertEqual(totals['attempts'], 1)
        self.assertEqual(totals['products'], 1)  # Attempt succeeded
        self.assertEqual(totals['no_product_matches'], 0)  # Attempt didn't fail
        
        # Validate invariants
        self.assertEqual(totals['attempts'], totals['products'] + totals['no_product_matches'])
    
    def test_enumerate_with_stats_uses_aggregator_totals(self):
        """Test that enumerate_with_stats returns aggregator-derived totals."""
        test_smiles = "c1ccccc1"  # Simple benzene
        cfg = EnumConfig(
            halogens=['F'],
            k_max=1,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        # Mock the enumerate_products generator to return controlled data
        with patch('src.halogenator.enumerate_k.enumerate_products') as mock_enum:
            # Mock return: one product and one QA summary marker
            mock_product = {'smiles': 'Fc1ccccc1', 'rule': 'R1', 'halogen': 'F'}
            mock_qa_final = {
                'no_product_matches': 999,  # This should be IGNORED - old per-match semantics
                'dedup_hits_statesig': 2,
                'dedup_hits_inchi': 1,
                'template_unsupported': 0,
                'qa_paths': {'isotope_unavailable': 888, 'isotope_miss': 0, 'atommap_used': 777, 'heuristic_used': 0}
            }
            
            # Return sequence: (product, snapshot), (marker, final_snapshot)
            mock_enum.return_value = [
                (mock_product, {}),
                ({'is_qa_summary_marker': True}, mock_qa_final)
            ]
            
            # Mock aggregator recording
            with patch.object(QAAggregator, 'record_attempt_result') as mock_record:
                products, qa_stats = enumerate_with_stats(test_smiles, cfg)
                
                # Verify products
                self.assertEqual(len(products), 1)
                self.assertEqual(products[0]['smiles'], 'Fc1ccccc1')
                
                # Key test: totals should come from aggregator, NOT from mock_qa_final
                # Since we didn't actually record anything in aggregator, totals should be 0
                self.assertEqual(qa_stats['no_product_matches'], 0)  # NOT 999
                self.assertEqual(qa_stats['qa_paths']['isotope_unavailable'], 0)  # NOT 888
                self.assertEqual(qa_stats['qa_paths']['atommap_used'], 0)  # NOT 777
                
                # But legacy dedup stats should come from mock_qa_final
                self.assertEqual(qa_stats['dedup_hits_statesig'], 2)
                self.assertEqual(qa_stats['dedup_hits_inchi'], 1)
                
                # Should have version 2 and pivots
                self.assertEqual(qa_stats['version'], '2')
                self.assertIn('pivots', qa_stats)
    
    def test_enumerate_k1_uses_aggregator_totals(self):
        """Test that enumerate_k1_with_stats returns aggregator-derived totals."""
        test_smiles = "c1ccccc1"  # Simple benzene
        cfg = EnumConfig(
            halogens=['F'],
            k_max=1,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
        
        # Mock the internal k1 enumeration
        with patch('src.halogenator.enumerate_k1._enumerate_k1_halogenation_with_stats_tracking') as mock_k1_enum:
            # Return empty (no products, no aggregator recording)
            mock_k1_enum.return_value = []
            
            # Mock molecule parsing
            with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
                mock_mol = MagicMock()
                mock_mol_from_smiles.return_value = mock_mol
                
                products, qa_stats = enumerate_k1_with_stats(test_smiles, cfg)
                
                # Should have empty results but correct structure
                self.assertEqual(len(products), 0)
                
                # Totals should be aggregator-derived (all zeros since no recording)
                self.assertEqual(qa_stats['no_product_matches'], 0)
                self.assertEqual(qa_stats['template_unsupported'], 0)
                for key in qa_stats['qa_paths']:
                    self.assertEqual(qa_stats['qa_paths'][key], 0)
                
                # Should have version 2 and pivots
                self.assertEqual(qa_stats['version'], '2')
                self.assertIn('pivots', qa_stats)
    
    def test_consistency_validation_always_passes_with_unified_totals(self):
        """Test that consistency validation passes when totals come from aggregator."""
        aggregator = QAAggregator()
        
        # Record some complex scenario
        aggregator.record_attempt_result('R1', 'F', 1, 1, {'isotope_unavailable': 1})
        aggregator.record_attempt_result('R1', 'Cl', 1, 0, {'isotope_miss': 2})  
        aggregator.record_attempt_result('R3', 'F', 2, 1, {'atommap_used': 1, 'heuristic_used': 1})
        
        # Get aggregator-derived totals
        totals = _compute_totals_from_aggregator(aggregator)
        
        # Build qa_stats as the enumerate functions would
        qa_stats_dict = {
            'no_product_matches': totals['no_product_matches'],
            'template_unsupported': totals['template_unsupported'], 
            'qa_paths': totals['qa_paths']
        }
        
        # Validation should pass without any errors or warnings
        from src.halogenator.enumerate_k import _validate_totals_pivots_consistency
        
        # Should not raise any exception (even in debug mode)
        import os
        with patch.dict(os.environ, {'HALO_ASSERT_PIVOT_CONSISTENCY': '1'}):
            _validate_totals_pivots_consistency(qa_stats_dict, aggregator)  # Should not raise


if __name__ == '__main__':
    unittest.main()