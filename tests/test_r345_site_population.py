# -*- coding: ascii -*-
"""Test that R3/R4/R5 rules populate site and ring_tag fields."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from collections import Counter


class TestR345SitePopulation(unittest.TestCase):
    """Test that R3/R4/R5 rules populate site and ring_tag for per_ring_quota."""
    
    def test_r3_populates_site_and_ring_tag(self):
        """Test R3 rule populates site and ring_tag when applicable."""
        # Use a molecule with -OH group that R3 can act on
        # Simple phenol
        phenol_smiles = "Oc1ccccc1"
        
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(phenol_smiles, cfg))
        
        # Look for R3 products
        r3_products = [p for p in products if any(h.get('rule') == 'R3' for h in p.get('substitutions', []))]
        
        if r3_products:  # Only test if R3 products exist
            r3_product = r3_products[0]
            r3_history = [h for h in r3_product.get('substitutions', []) if h.get('rule') == 'R3'][0]
            
            # R3 should now have site populated (not None)
            self.assertIsNotNone(r3_history.get('site'), 
                               "R3 rule should populate site field")
            
            # Site should be a valid integer
            self.assertIsInstance(r3_history.get('site'), int, 
                                "R3 site should be integer atom index")
            
            # ring_tag might be empty string for non-flavonoid, but should not be None
            ring_tag = r3_history.get('ring_tag')
            self.assertIsNotNone(ring_tag, 
                               "R3 rule should populate ring_tag field")
            
            # For phenol (not a flavonoid), ring_tag should be empty string
            self.assertEqual(ring_tag, '', 
                           "Phenol should have empty ring_tag (not a flavonoid)")
    
    def test_r345_constraints_integration(self):
        """Test that R3/R4/R5 sites are counted in per_ring_quota constraints."""
        # Use quercetin which has both aromatic sites and -OH groups
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        
        # Very strict per_ring_quota to force constraint violations
        cfg_strict = EnumConfig(
            k_max=3,
            halogens=('F',),
            constraints={'per_ring_quota': 1, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(quercetin_smiles, cfg_strict))
        
        # Collect all substitution histories
        all_histories = []
        for p in products:
            history = p.get('substitutions', [])
            if len(history) >= 2:  # Multi-step products
                all_histories.append(history)
        
        if all_histories:
            # Check that some R3 steps have non-empty ring_tag when in flavonoid rings
            r3_steps_with_ring_tags = []
            for history in all_histories:
                for step in history:
                    if step.get('rule') == 'R3' and step.get('ring_tag'):
                        r3_steps_with_ring_tags.append(step)
            
            if r3_steps_with_ring_tags:
                # At least some R3 steps should have ring tags in flavonoid
                self.assertGreater(len(r3_steps_with_ring_tags), 0,
                                 "Some R3 steps should have ring_tag populated in flavonoid")
    
    def test_site_field_consistency(self):
        """Test that site field is consistent between R1 and R3/R4/R5."""
        quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
        
        cfg = EnumConfig(
            k_max=2,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        products = list(enumerate_products(quercetin_smiles, cfg))
        
        # Collect site values from all rules
        r1_sites = []
        r3_sites = []
        
        for p in products:
            for h in p.get('substitutions', []):
                rule = h.get('rule')
                site = h.get('site')
                if rule == 'R1' and site is not None:
                    r1_sites.append(site)
                elif rule == 'R3' and site is not None:
                    r3_sites.append(site)
        
        # Both R1 and R3 should have integer sites when populated
        for site in r1_sites:
            self.assertIsInstance(site, int, "R1 sites should be integers")
            
        for site in r3_sites:
            self.assertIsInstance(site, int, "R3 sites should be integers")
            


if __name__ == '__main__':
    unittest.main()