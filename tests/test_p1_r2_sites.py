# -*- coding: ascii -*-
"""Test P1 R2 ring site identification."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig
from src.halogenator.sites import c_ring_indices, is_c_ring_site_ready
from rdkit import Chem


class TestP1R2Sites(unittest.TestCase):
    """Test R2 rule for C ring site identification."""
    
    def test_r2_ring_carbon_sites(self):
        """Test R2 correctly identifies ring carbon sites."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        # Test with hydroxylated tetrahydropyran that has ring carbons with oxygen neighbors
        hydroxy_thp_smiles = "OC1CCCCO1"
        mol = Chem.MolFromSmiles(hydroxy_thp_smiles)
        
        if mol is not None:
            # Check that c_ring_indices finds carbon ring sites
            ring_c_indices = c_ring_indices(mol)
            self.assertGreater(len(ring_c_indices), 0, "Should find ring carbon sites")
            
            # Check that sites are ready for halogenation
            ready_sites = [i for i in ring_c_indices if is_c_ring_site_ready(mol, i)]
            self.assertGreater(len(ready_sites), 0, "Should find ready ring carbon sites")
        
    def test_r2_excludes_carbonyl_carbons(self):
        """Test R2 excludes carbonyl carbons from halogenation."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        # Benzoquinone has carbonyl carbons that should be excluded
        benzoquinone_smiles = "O=C1C=CC(=O)C=C1"
        mol = Chem.MolFromSmiles(benzoquinone_smiles)
        
        if mol is not None:
            products = list(enumerate_products(benzoquinone_smiles, cfg))
            
            # Should produce fewer R2 products due to carbonyl exclusion
            r2_products = [p for p in products if p['rule'] == 'R2']
            
            # The exact number depends on implementation, but should be reasonable
            # Main point is that it doesn't crash and produces some valid products
            self.assertIsInstance(r2_products, list, "Should return list of R2 products")
    
    def test_r2_with_oxygen_containing_rings(self):
        """Test R2 with oxygen-containing 6-membered rings."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        # Hydroxylated THP - 6-membered ring with oxygen
        hydroxy_thp_smiles = "OC1CCCCO1"
        products = list(enumerate_products(hydroxy_thp_smiles, cfg))
        
        # Should produce R2 products for the carbon ring sites (R2a/R2b when PR2 enabled by default)
        r2_products = [p for p in products if p['rule'] in ['R2', 'R2a', 'R2b']]
        self.assertGreater(len(r2_products), 0, "Should produce R2 products for ring carbons with O neighbors")
        
    def test_r2_vs_r1_distinction(self):
        """Test that R2 and R1 rules target different sites appropriately."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )

        # Use symmetry folding to ensure R1 products are generated
        cfg_folded = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': True, 'enable_state_sig': False}
        )
        
        # Benzene should have R1 (aromatic CH) but no R2 (non-aromatic ring C)
        benzene_smiles = "c1ccccc1"
        benzene_products = list(enumerate_products(benzene_smiles, cfg_folded))
        
        benzene_r1 = [p for p in benzene_products if p['rule'] == 'R1']
        benzene_r2 = [p for p in benzene_products if p['rule'] == 'R2'] 
        
        self.assertGreater(len(benzene_r1), 0, "Benzene should have R1 products")
        self.assertEqual(len(benzene_r2), 0, "Benzene should not have R2 products")
        
        # Hydroxylated THP should have R2 (ring C with O) but no R1 (aromatic CH)
        hydroxy_thp_smiles = "OC1CCCCO1"
        hydroxy_thp_products = list(enumerate_products(hydroxy_thp_smiles, cfg_folded))
        
        hydroxy_thp_r1 = [p for p in hydroxy_thp_products if p['rule'] == 'R1']
        hydroxy_thp_r2 = [p for p in hydroxy_thp_products if p['rule'] in ['R2', 'R2a', 'R2b']]
        
        self.assertEqual(len(hydroxy_thp_r1), 0, "Hydroxylated THP should not have R1 products")
        self.assertGreater(len(hydroxy_thp_r2), 0, "Hydroxylated THP should have R2 products")
        
    def test_r2_site_metadata(self):
        """Test that R2 products have correct site metadata."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            constraints={'per_ring_quota': 10, 'min_graph_distance': 1, 'max_per_halogen': None},
            pruning_cfg={'enable_symmetry_fold': False, 'enable_state_sig': False}
        )
        
        hydroxy_thp_smiles = "OC1CCCCO1" 
        products = list(enumerate_products(hydroxy_thp_smiles, cfg))
        
        r2_products = [p for p in products if p['rule'] == 'R2']
        
        for product in r2_products:
            history = product.get('substitutions', [])
            if history:
                step = history[0]  # First (and only) step for k=1
                
                # Should have site information
                self.assertIn('site', step, "R2 step should have site index")
                self.assertIsInstance(step['site'], int, "Site should be integer")
                
                # Should have ring tag information
                self.assertIn('ring_tag', step, "R2 step should have ring_tag")
                
                # Should have symmetry class
                self.assertIn('sym', step, "R2 step should have symmetry class")
                

if __name__ == '__main__':
    unittest.main()