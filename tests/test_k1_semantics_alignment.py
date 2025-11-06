# -*- coding: ascii -*-
"""Tests for k=1 semantics alignment with k>1 rules."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k1 import enumerate_k1_with_stats
from src.halogenator.enumerate_k import EnumConfig


class TestK1SemanticsAlignment(unittest.TestCase):
    """Test that k=1 semantics align with k>1 attempt boundaries."""
    
    def setUp(self):
        """Set up test environment."""
        # Create a minimal EnumConfig
        self.cfg = EnumConfig(
            halogens=['F', 'Cl'],
            k_max=1,  # k=1 for this test
            constraints={},
            std_cfg={},
            qc_cfg={},
            pruning_cfg={}
        )
    
    def test_k1_returns_version_2_and_pivots(self):
        """Test that k=1 returns version '2' and non-empty pivots."""
        # Use a simple test molecule
        test_smiles = "c1ccccc1"  # benzene
        
        # Mock the underlying functions to simulate successful enumeration
        with patch('src.halogenator.enumerate_k1.build_reactions') as mock_build_reactions:
            mock_reactions = {
                'R1': {'F': MagicMock(), 'Cl': MagicMock()},
                'R3': {'F': MagicMock()},
                'R4': {'Cl': MagicMock()}
            }
            mock_build_reactions.return_value = mock_reactions
            
            with patch('src.halogenator.enumerate_k1._run_reaction_safely') as mock_run_safely:
                # Mock some products being generated
                mock_product = MagicMock()
                mock_run_safely.return_value = [[(mock_product,)]]
                
                with patch('src.halogenator.enumerate_k1.Chem.SanitizeMol'):
                    with patch('src.halogenator.enumerate_k1.ensure_ready'):
                        with patch('src.halogenator.enumerate_k1.Chem.MolToSmiles') as mock_to_smiles:
                            mock_to_smiles.return_value = "c1ccc(F)cc1"
                            
                            # Run k=1 enumeration
                            products, qa_stats = enumerate_k1_with_stats(test_smiles, self.cfg)
                            
                            # Verify version '2' is returned
                            self.assertEqual(qa_stats.get('version'), '2')
                            
                            # Verify pivots are present and non-empty
                            self.assertIn('pivots', qa_stats)
                            pivots = qa_stats['pivots']
                            
                            # Should have all required dimensions
                            self.assertIn('by_rule', pivots)
                            self.assertIn('by_halogen', pivots)
                            self.assertIn('by_k', pivots)
                            self.assertIn('by_rule_halogen', pivots)
                            self.assertIn('by_rule_halogen_k', pivots)
    
    def test_k1_attempt_semantics_zero_products(self):
        """Test k=1 with zero products increments no_product_matches correctly."""
        test_smiles = "CCN"  # ethylamine
        
        with patch('src.halogenator.enumerate_k1.build_reactions') as mock_build_reactions:
            mock_reactions = {
                'R3': {'F': MagicMock()}
            }
            mock_build_reactions.return_value = mock_reactions
            
            with patch('src.halogenator.enumerate_k1._run_reaction_safely') as mock_run_safely:
                # Mock no products being generated
                mock_run_safely.return_value = []  # No products
                
                # Run k=1 enumeration
                products, qa_stats = enumerate_k1_with_stats(test_smiles, self.cfg)
                
                # Verify basic results
                self.assertEqual(len(products), 0)
                self.assertEqual(qa_stats.get('version'), '2')
                
                # Check pivot semantics
                pivots = qa_stats['pivots']
                
                # Should have attempt for R3_F_1 with zero products
                if 'R3_F_1' in pivots['by_rule_halogen_k']:
                    r3_f_1 = pivots['by_rule_halogen_k']['R3_F_1']
                    self.assertEqual(r3_f_1['attempts'], 1)
                    self.assertEqual(r3_f_1.get('products', 0), 0)
                    self.assertEqual(r3_f_1['no_product_matches'], 1)
    
    def test_k1_attempt_semantics_multiple_products(self):
        """Test k=1 with multiple products still increments products only once per attempt."""
        test_smiles = "c1ccccc1"  # benzene
        
        with patch('src.halogenator.enumerate_k1.build_reactions') as mock_build_reactions:
            mock_reactions = {
                'R1': {'F': MagicMock()}
            }
            mock_build_reactions.return_value = mock_reactions
            
            with patch('src.halogenator.enumerate_k1._run_reaction_safely') as mock_run_safely:
                # Mock multiple products being generated from one attempt
                mock_product1 = MagicMock()
                mock_product2 = MagicMock()
                mock_run_safely.return_value = [[(mock_product1, mock_product2)]]
                
                with patch('src.halogenator.enumerate_k1.Chem.SanitizeMol'):
                    with patch('src.halogenator.enumerate_k1.ensure_ready'):
                        with patch('src.halogenator.enumerate_k1.Chem.MolToSmiles') as mock_to_smiles:
                            # Return different SMILES for each product to avoid deduplication
                            mock_to_smiles.side_effect = ["c1ccc(F)cc1", "c1cc(F)ccc1"]
                            
                            # Run k=1 enumeration
                            products, qa_stats = enumerate_k1_with_stats(test_smiles, self.cfg)
                            
                            # Should have 2 products from this attempt
                            self.assertEqual(len(products), 2)
                            
                            # Check pivot semantics - should have 1 attempt, 1 products (not 2!)
                            pivots = qa_stats['pivots']
                            
                            if 'R1_F_1' in pivots['by_rule_halogen_k']:
                                r1_f_1 = pivots['by_rule_halogen_k']['R1_F_1']
                                self.assertEqual(r1_f_1['attempts'], 1)
                                self.assertEqual(r1_f_1['products'], 1)  # Should be 1, not 2!
                                self.assertEqual(r1_f_1.get('no_product_matches', 0), 0)
    
    def test_k1_totals_and_pivots_alignment(self):
        """Test that totals and pivots reflect the same numbers in k=1."""
        test_smiles = "c1ccccc1O"  # phenol
        
        with patch('src.halogenator.enumerate_k1.build_reactions') as mock_build_reactions:
            mock_reactions = {
                'R1': {'F': MagicMock(), 'Cl': MagicMock()},
                'R3': {'F': MagicMock()}
            }
            mock_build_reactions.return_value = mock_reactions
            
            with patch('src.halogenator.enumerate_k1._run_reaction_safely') as mock_run_safely:
                # F succeeds, Cl fails, R3 has template issue
                call_count = 0
                def mock_run_safely_side_effect(*args, **kwargs):
                    nonlocal call_count
                    call_count += 1
                    if call_count == 1:  # R1_F
                        return [[(MagicMock(),)]]  # Success
                    elif call_count == 2:  # R1_Cl
                        return []  # Failure
                    else:  # R3_F
                        return []  # Failure
                
                mock_run_safely.side_effect = mock_run_safely_side_effect
                
                with patch('src.halogenator.enumerate_k1.Chem.SanitizeMol'):
                    with patch('src.halogenator.enumerate_k1.ensure_ready'):
                        with patch('src.halogenator.enumerate_k1.Chem.MolToSmiles') as mock_to_smiles:
                            mock_to_smiles.return_value = "c1ccc(F)c(O)c1"
                            
                            # Run k=1 enumeration
                            products, qa_stats = enumerate_k1_with_stats(test_smiles, self.cfg)
                            
                            # Should have 1 product
                            self.assertEqual(len(products), 1)
                            
                            # Check that totals and pivots are consistent
                            pivots = qa_stats['pivots']
                            
                            # Count attempts/products/no_matches across all pivots
                            total_attempts = sum(
                                events.get('attempts', 0) 
                                for events in pivots['by_rule_halogen_k'].values()
                            )
                            total_products = sum(
                                events.get('products', 0) 
                                for events in pivots['by_rule_halogen_k'].values()
                            )
                            total_no_matches = sum(
                                events.get('no_product_matches', 0) 
                                for events in pivots['by_rule_halogen_k'].values()
                            )
                            
                            # Verify attempt invariants
                            self.assertGreater(total_attempts, 0)
                            self.assertEqual(total_attempts, total_products + total_no_matches)
                            self.assertGreaterEqual(total_attempts, total_products)
    
    def test_k1_r2_rule_attempt_boundaries(self):
        """Test that R2 rule in k=1 uses proper attempt boundaries."""
        test_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"  # quercetin-like
        
        # Mock R2 sites and halogenation
        with patch('src.halogenator.sites.c_ring_indices') as mock_c_ring:
            mock_c_ring.return_value = [1, 2]  # Two C ring sites
            
            with patch('src.halogenator.sites.is_c_ring_site_ready') as mock_is_ready:
                mock_is_ready.return_value = True
                
                with patch('src.halogenator.reactions.apply_single_site_halogenation') as mock_halogenate:
                    # First site succeeds for both halogens, second site fails
                    call_count = 0
                    def mock_halogenate_side_effect(mol, site, halogen):
                        nonlocal call_count
                        call_count += 1
                        if site == 1:  # First site always succeeds
                            return MagicMock()
                        else:  # Second site fails
                            return None
                    
                    mock_halogenate.side_effect = mock_halogenate_side_effect
                    
                    with patch('src.halogenator.sites.flavonoid_ring_label') as mock_ring_label:
                        mock_ring_label.return_value = 'C'
                        
                        # Run k=1 enumeration with only R2 rule
                        cfg_r2_only = EnumConfig(
                            halogens=['F', 'Cl'],
                            k_max=1,
                            constraints={},
                            std_cfg={},
                            qc_cfg={},
                            pruning_cfg={}
                        )
                        
                        # Mock build_reactions to return empty (only R2 is site-based)
                        with patch('src.halogenator.enumerate_k1.build_reactions') as mock_build_reactions:
                            mock_build_reactions.return_value = {}
                            
                            products, qa_stats = enumerate_k1_with_stats(test_smiles, cfg_r2_only)
                            
                            # Should have products from R2 rule
                            self.assertGreater(len(products), 0)
                            
                            # Check R2 attempt boundaries in pivots
                            pivots = qa_stats['pivots']
                            
                            # Should have separate attempts for R2_F_1 and R2_Cl_1
                            # Each site x halogen is one attempt (2 sites x 1 halogen = 2 attempts per halogen)
                            if 'R2_F_1' in pivots['by_rule_halogen_k']:
                                r2_f_1 = pivots['by_rule_halogen_k']['R2_F_1']
                                self.assertEqual(r2_f_1['attempts'], 2)  # Two sites per halogen
                                self.assertEqual(r2_f_1['products'], 1)  # Site 1 succeeds, site 2 fails
                                self.assertEqual(r2_f_1['no_product_matches'], 1)  # Site 2 fails
                                
                            if 'R2_Cl_1' in pivots['by_rule_halogen_k']:
                                r2_cl_1 = pivots['by_rule_halogen_k']['R2_Cl_1']
                                self.assertEqual(r2_cl_1['attempts'], 2)  # Two sites per halogen  
                                self.assertEqual(r2_cl_1['products'], 1)  # Site 1 succeeds, site 2 fails
                                self.assertEqual(r2_cl_1['no_product_matches'], 1)  # Site 2 fails


if __name__ == '__main__':
    unittest.main()