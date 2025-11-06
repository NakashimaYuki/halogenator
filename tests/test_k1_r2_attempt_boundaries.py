# -*- coding: ascii -*-
"""Tests specifically for k=1 R2 rule attempt boundaries."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.enumerate_k import QAAggregator, EnumConfig
from src.halogenator.enumerate_k1 import _enumerate_k1_halogenation_with_stats_tracking
from rdkit import Chem


class TestK1R2AttemptBoundaries(unittest.TestCase):
    """Test that k=1 R2 rule uses correct site x halogen attempt boundaries."""
    
    def setUp(self):
        """Set up test environment."""
        self.aggregator = QAAggregator()
        
        # Simple test molecule
        self.mol = Chem.MolFromSmiles('c1ccccc1')  # benzene
        self.parent_smiles = 'c1ccccc1'
        self.parent_inchikey = 'UHOVQNZJYSORNB-UHFFFAOYSA-N'
        
        self.halogens = ['F', 'Cl']
        self.rules = ['R2']
        self.config = {}
        self.stats_dict = {
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': {'isotope_unavailable': 0, 'isotope_miss': 0, 'atommap_used': 0, 'heuristic_used': 0}
        }
    
    def test_r2_attempts_per_site_halogen(self):
        """Test that R2 rule creates one attempt per site x halogen combination."""
        # Mock the site detection to return 2 sites
        with patch('src.halogenator.enumerate_k1.c_ring_indices') as mock_c_ring:
            mock_c_ring.return_value = [0, 1]  # Two sites
            
            with patch('src.halogenator.enumerate_k1.is_c_ring_site_ready') as mock_is_ready:
                mock_is_ready.return_value = True  # Both sites are ready
                
                with patch('src.halogenator.enumerate_k1.apply_single_site_halogenation') as mock_halogenate:
                    # Site 0 succeeds for F, fails for Cl
                    # Site 1 fails for both halogens
                    def mock_halogenate_side_effect(mol, site, halogen):
                        if site == 0 and halogen == 'F':
                            return MagicMock()  # Success
                        else:
                            return None  # Failure
                    
                    mock_halogenate.side_effect = mock_halogenate_side_effect
                    
                    with patch('src.halogenator.enumerate_k1.flavonoid_ring_label') as mock_ring_label:
                        mock_ring_label.return_value = 'A'
                        
                        # Run the k=1 enumeration with stats tracking
                        products = _enumerate_k1_halogenation_with_stats_tracking(
                            self.mol, self.halogens, self.rules, self.config, self.stats_dict, self.aggregator
                        )
                        
                        # Check aggregator state
                        pivots = self.aggregator.to_pivots_dict()
                        
                        # Should have 1 product (site 0 + F)
                        self.assertEqual(len(products), 1)
                        
                        # Verify the product has correct properties
                        product_mol, product_props = products[0]
                        self.assertEqual(product_props['rule'], 'R2')
                        self.assertEqual(product_props['site_index'], 0)
                        self.assertEqual(product_props['halogen'], 'F')
                        self.assertEqual(product_props['ring_tag'], 'A')
                        
                        # Check attempt boundaries in aggregator
                        pivots = self.aggregator.to_pivots_dict()
                        
                        # Should have attempts for each site x halogen combination
                        # 2 sites x 2 halogens = 4 attempts total
                        
                        # R2_F_1: 2 attempts (sites 0,1), 1 product (site 0), 1 no_product_match (site 1)
                        self.assertIn('by_rule_halogen_k', pivots)
                        self.assertIn('R2_F_1', pivots['by_rule_halogen_k'])
                        r2_f_1 = pivots['by_rule_halogen_k']['R2_F_1']
                        self.assertEqual(r2_f_1['attempts'], 2, "R2_F_1 should have 2 attempts")
                        self.assertEqual(r2_f_1['products'], 1, "R2_F_1 should have 1 product") 
                        self.assertEqual(r2_f_1.get('no_product_matches', 0), 1, "R2_F_1 should have 1 no_product_match")
                        
                        # R2_Cl_1: 2 attempts (sites 0,1), 0 products, 2 no_product_matches
                        self.assertIn('R2_Cl_1', pivots['by_rule_halogen_k'])
                        r2_cl_1 = pivots['by_rule_halogen_k']['R2_Cl_1']
                        self.assertEqual(r2_cl_1['attempts'], 2, "R2_Cl_1 should have 2 attempts")
                        self.assertEqual(r2_cl_1.get('products', 0), 0, "R2_Cl_1 should have 0 products")
                        self.assertEqual(r2_cl_1.get('no_product_matches', 0), 2, "R2_Cl_1 should have 2 no_product_matches")
                        
                        # Total verification  
                        self.assertIn('by_rule', pivots)
                        self.assertIn('R2', pivots['by_rule'])
                        r2_totals = pivots['by_rule']['R2']
                        self.assertEqual(r2_totals['attempts'], 4, "R2 should have 4 total attempts")
                        self.assertEqual(r2_totals['products'], 1, "R2 should have 1 total product")
                        self.assertEqual(r2_totals.get('no_product_matches', 0), 3, "R2 should have 3 total no_product_matches")
    
    def test_r2_no_sites_available(self):
        """Test R2 behavior when no sites are available."""
        # Mock no sites available
        with patch('src.halogenator.sites.c_ring_indices') as mock_c_ring:
            mock_c_ring.return_value = []  # No sites
            
            # Run enumeration
            products = _enumerate_k1_halogenation_with_stats_tracking(
                self.mol, self.halogens, self.rules, self.config, self.stats_dict, self.aggregator
            )
            
            # Should have no products
            self.assertEqual(len(products), 0)
            
            # Should have no attempts recorded either (no sites to attempt)
            pivots = self.aggregator.to_pivots_dict()
            by_rule_halogen_k = pivots.get('by_rule_halogen_k', {})
            
            # Should not have any R2 entries
            r2_entries = [key for key in by_rule_halogen_k.keys() if key.startswith('R2_')]
            self.assertEqual(len(r2_entries), 0, "Should have no R2 attempts when no sites available")
    
    def test_r2_exception_handling(self):
        """Test R2 behavior when site processing raises exceptions."""
        # Mock sites that cause exceptions
        with patch('src.halogenator.sites.c_ring_indices') as mock_c_ring:
            mock_c_ring.side_effect = Exception("Site detection failed")
            
            # Run enumeration - should handle exception gracefully
            products = _enumerate_k1_halogenation_with_stats_tracking(
                self.mol, self.halogens, self.rules, self.config, self.stats_dict, self.aggregator
            )
            
            # Should have no products
            self.assertEqual(len(products), 0)
            
            # Should record template_unsupported attempts for each halogen
            pivots = self.aggregator.to_pivots_dict()
            by_rule_halogen_k = pivots.get('by_rule_halogen_k', {})
            
            # Should have template_unsupported for each halogen
            for halogen in self.halogens:
                key = f'R2_{halogen}_1'
                if key in by_rule_halogen_k:
                    self.assertGreater(by_rule_halogen_k[key].get('template_unsupported', 0), 0)


if __name__ == '__main__':
    unittest.main()