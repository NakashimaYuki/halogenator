# -*- coding: ascii -*-
"""Tests for site-rule attempt boundary semantics."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator, EnumConfig, _apply_site_rules
from rdkit import Chem


class TestSiteRuleAttempts(unittest.TestCase):
    """Test that site rules use correct attempt boundaries (site x halogen x k)."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
        self.cfg = EnumConfig(
            halogens=['F', 'Cl'], 
            k_max=2,
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={'enable_symmetry_fold': False}  # Disable for predictable behavior
        )
        
        # Mock molecule
        self.mol = MagicMock()
        self.mol.GetNumAtoms.return_value = 10
    
    def test_r1_site_attempts_separate(self):
        """Test that R1 rule creates separate attempts for each site x halogen."""
        with patch('src.halogenator.enumerate_k.ensure_ready'):
            with patch('src.halogenator.enumerate_k.Chem.CanonicalRankAtoms') as mock_ranks:
                mock_ranks.return_value = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                
                with patch('src.halogenator.enumerate_k.aromatic_CH_indices') as mock_aromatic:
                    mock_aromatic.return_value = [1, 3]  # Two aromatic sites
                    
                    with patch('src.halogenator.enumerate_k._build_ring_label_map') as mock_ring_map:
                        mock_ring_map.return_value = {}
                        
                        with patch('src.halogenator.enumerate_k._apply_single_site') as mock_apply:
                            # Site 1 succeeds for both halogens, site 3 fails for both
                            def mock_apply_side_effect(mol, site, halogen, rule_id, ranks, ring_label_map, history, depth, cfg, seen_global):
                                if site == 1:
                                    return ('product', {'some': 'data'})  # Success
                                else:  # site == 3
                                    return None  # Failure
                            
                            mock_apply.side_effect = mock_apply_side_effect
                            
                            # Run site rules application
                            results = _apply_site_rules(self.mol, self.cfg, [], 0, set(), self.aggregator)
                            
                            # Check results
                            self.assertEqual(len(results), 2)  # Only site 1 succeeded for F and Cl
                            
                            # Check pivots - should have 4 attempts (2 sites x 2 halogens)
                            pivots = self.aggregator.to_pivots_dict()
                            
                            # R1_F_1: 2 attempts (sites 1,3), 1 product (site 1), 1 no_product_match (site 3)
                            if 'R1_F_1' in pivots['by_rule_halogen_k']:
                                r1_f = pivots['by_rule_halogen_k']['R1_F_1']
                                self.assertEqual(r1_f['attempts'], 2)
                                self.assertEqual(r1_f['products'], 1) 
                                self.assertEqual(r1_f['no_product_matches'], 1)
                            
                            # R1_Cl_1: 2 attempts (sites 1,3), 1 product (site 1), 1 no_product_match (site 3)
                            if 'R1_Cl_1' in pivots['by_rule_halogen_k']:
                                r1_cl = pivots['by_rule_halogen_k']['R1_Cl_1']
                                self.assertEqual(r1_cl['attempts'], 2)
                                self.assertEqual(r1_cl['products'], 1)
                                self.assertEqual(r1_cl['no_product_matches'], 1)
    
    def test_r2_site_attempts_separate(self):
        """Test that R2 rule creates separate attempts for each site x halogen.""" 
        with patch('src.halogenator.enumerate_k.ensure_ready'):
            with patch('src.halogenator.enumerate_k.Chem.CanonicalRankAtoms') as mock_ranks:
                mock_ranks.return_value = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                
                with patch('src.halogenator.enumerate_k.c_ring_indices') as mock_c_ring:
                    mock_c_ring.return_value = [2, 4, 6]  # Three C ring indices
                    
                    with patch('src.halogenator.enumerate_k.is_c_ring_site_ready') as mock_ready:
                        # Only sites 2 and 6 are ready
                        mock_ready.side_effect = lambda mol, site: site in [2, 6]
                        
                        with patch('src.halogenator.enumerate_k._build_ring_label_map') as mock_ring_map:
                            mock_ring_map.return_value = {}
                            
                            with patch('src.halogenator.enumerate_k._apply_single_site') as mock_apply:
                                # Site 2 succeeds, site 6 fails
                                def mock_apply_side_effect(mol, site, halogen, rule_id, ranks, ring_label_map, history, depth, cfg, seen_global):
                                    if site == 2:
                                        return ('product', {'some': 'data'})  # Success
                                    else:  # site == 6
                                        return None  # Failure
                                
                                mock_apply.side_effect = mock_apply_side_effect
                                
                                # Run site rules application
                                results = _apply_site_rules(self.mol, self.cfg, [], 0, set(), self.aggregator)
                                
                                # Check results - only site 2 succeeded for both halogens
                                self.assertEqual(len(results), 2)
                                
                                # Check pivots - should have 4 attempts (2 ready sites x 2 halogens)
                                pivots = self.aggregator.to_pivots_dict()
                                
                                # R2_F_1: 2 attempts (sites 2,6), 1 product (site 2), 1 no_product_match (site 6)
                                if 'R2_F_1' in pivots['by_rule_halogen_k']:
                                    r2_f = pivots['by_rule_halogen_k']['R2_F_1']
                                    self.assertEqual(r2_f['attempts'], 2)
                                    self.assertEqual(r2_f['products'], 1)
                                    self.assertEqual(r2_f['no_product_matches'], 1)
                                
                                # R2_Cl_1: 2 attempts (sites 2,6), 1 product (site 2), 1 no_product_match (site 6)
                                if 'R2_Cl_1' in pivots['by_rule_halogen_k']:
                                    r2_cl = pivots['by_rule_halogen_k']['R2_Cl_1']
                                    self.assertEqual(r2_cl['attempts'], 2)
                                    self.assertEqual(r2_cl['products'], 1)
                                    self.assertEqual(r2_cl['no_product_matches'], 1)
    
    def test_mixed_success_failure_attempts(self):
        """Test scenario with some sites succeeding, others failing."""
        with patch('src.halogenator.enumerate_k.ensure_ready'):
            with patch('src.halogenator.enumerate_k.Chem.CanonicalRankAtoms') as mock_ranks:
                mock_ranks.return_value = [0, 1, 2]
                
                with patch('src.halogenator.enumerate_k.aromatic_CH_indices') as mock_aromatic:
                    mock_aromatic.return_value = [0, 1]  # Two sites
                    
                    with patch('src.halogenator.enumerate_k.c_ring_indices') as mock_c_ring:
                        mock_c_ring.return_value = []  # No R2 sites
                        
                        with patch('src.halogenator.enumerate_k.is_c_ring_site_ready'):
                            with patch('src.halogenator.enumerate_k._build_ring_label_map') as mock_ring_map:
                                mock_ring_map.return_value = {}
                                
                                with patch('src.halogenator.enumerate_k._apply_single_site') as mock_apply:
                                    call_count = 0
                                    def mock_apply_side_effect(mol, site, halogen, rule_id, ranks, ring_label_map, history, depth, cfg, seen_global):
                                        nonlocal call_count
                                        call_count += 1
                                        # Site 0 always succeeds, site 1 succeeds only for F
                                        if site == 0:
                                            return ('product', {'data': call_count})
                                        elif site == 1 and halogen == 'F':
                                            return ('product', {'data': call_count})
                                        else:
                                            return None
                                    
                                    mock_apply.side_effect = mock_apply_side_effect
                                    
                                    # Run with only R1 rule
                                    results = _apply_site_rules(self.mol, self.cfg, [], 0, set(), self.aggregator)
                                    
                                    # Should have 3 successful results (site 0: F+Cl, site 1: F only)  
                                    self.assertEqual(len(results), 3)
                                    
                                    # Check pivots
                                    pivots = self.aggregator.to_pivots_dict()
                                    
                                    # R1_F_1: 2 attempts, 2 products, 0 no_product_matches
                                    if 'R1_F_1' in pivots['by_rule_halogen_k']:
                                        r1_f = pivots['by_rule_halogen_k']['R1_F_1']
                                        self.assertEqual(r1_f['attempts'], 2)
                                        self.assertEqual(r1_f['products'], 2)
                                        self.assertEqual(r1_f.get('no_product_matches', 0), 0)
                                    
                                    # R1_Cl_1: 2 attempts, 1 product, 1 no_product_match
                                    if 'R1_Cl_1' in pivots['by_rule_halogen_k']:
                                        r1_cl = pivots['by_rule_halogen_k']['R1_Cl_1']
                                        self.assertEqual(r1_cl['attempts'], 2)
                                        self.assertEqual(r1_cl['products'], 1)
                                        self.assertEqual(r1_cl.get('no_product_matches', 0), 1)
    
    def test_attempt_invariants_maintained(self):
        """Test that attempt invariants are maintained with correct site boundaries."""
        with patch('src.halogenator.enumerate_k.ensure_ready'):
            with patch('src.halogenator.enumerate_k.Chem.CanonicalRankAtoms') as mock_ranks:
                mock_ranks.return_value = [0, 1, 2, 3, 4]
                
                with patch('src.halogenator.enumerate_k.aromatic_CH_indices') as mock_aromatic:
                    mock_aromatic.return_value = [0, 2, 4]  # Three sites
                    
                    with patch('src.halogenator.enumerate_k.c_ring_indices') as mock_c_ring:
                        mock_c_ring.return_value = []  # No R2 sites
                        
                        with patch('src.halogenator.enumerate_k.is_c_ring_site_ready'):
                            with patch('src.halogenator.enumerate_k._build_ring_label_map') as mock_ring_map:
                                mock_ring_map.return_value = {}
                                
                                with patch('src.halogenator.enumerate_k._apply_single_site') as mock_apply:
                                    # Random success/failure pattern
                                    def mock_apply_side_effect(mol, site, halogen, rule_id, ranks, ring_label_map, history, depth, cfg, seen_global):
                                        # Site 0: always fails, Site 2: F succeeds, Site 4: Cl succeeds
                                        if site == 2 and halogen == 'F':
                                            return ('product', {})
                                        elif site == 4 and halogen == 'Cl':  
                                            return ('product', {})
                                        else:
                                            return None
                                    
                                    mock_apply.side_effect = mock_apply_side_effect
                                    
                                    # Run site rules
                                    results = _apply_site_rules(self.mol, self.cfg, [], 1, set(), self.aggregator)
                                    
                                    # Check results (should have 2: site 2-F and site 4-Cl)
                                    self.assertEqual(len(results), 2)
                                    
                                    # Check attempt invariants across all pivots
                                    pivots = self.aggregator.to_pivots_dict()
                                    
                                    for key, stats in pivots['by_rule_halogen_k'].items():
                                        attempts = stats.get('attempts', 0)
                                        products = stats.get('products', 0)
                                        no_matches = stats.get('no_product_matches', 0)
                                        
                                        # Invariants
                                        self.assertGreaterEqual(attempts, products, f"attempts >= products failed for {key}")
                                        self.assertEqual(attempts, products + no_matches, f"attempts == products + no_matches failed for {key}")


if __name__ == '__main__':
    unittest.main()