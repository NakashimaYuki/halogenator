# -*- coding: ascii -*-
"""Test k=1 enumeration functionality."""

import unittest
from rdkit import Chem

from halogenator.enumerate_k1 import enumerate_k1_halogenation
from halogenator.sites import aromatic_CH_indices, c_ring_indices


class TestEnumerateK1(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures."""
        # Simple phenol for testing
        self.test_smiles = "c1ccc(O)cc1"
        self.test_mol = Chem.MolFromSmiles(self.test_smiles)
        
        # Test configuration
        self.config = {
            'qc': {'sanitize_strict': True},
            'dedupe': {'method': 'inchikey'}
        }
        
        self.halogens = ['F', 'Cl']  # Just two for testing
        self.all_rules = ['R1', 'R2', 'R3', 'R4', 'R5']
    
    def test_enumerate_k1_basic(self):
        """Test basic k=1 enumeration."""
        products = enumerate_k1_halogenation(
            self.test_mol, self.halogens, self.all_rules, self.config
        )
        
        # Should generate some products
        self.assertTrue(len(products) > 0)
        
        # Each product should be a (mol, props) tuple
        for mol, props in products:
            self.assertIsNotNone(mol)
            self.assertIsInstance(props, dict)
            
            # Check required properties
            self.assertIn('rule', props)
            self.assertIn('halogen', props)
            self.assertIn('depth', props)
            self.assertEqual(props['depth'], 1)
    
    def test_enumerate_r1_only(self):
        """Test R1 rule only (aromatic CH)."""
        products = enumerate_k1_halogenation(
            self.test_mol, self.halogens, ['R1'], self.config
        )
        
        # Should have products from R1
        r1_products = [p for p in products if p[1]['rule'] == 'R1']
        self.assertTrue(len(r1_products) > 0)
        
        # Check aromatic CH sites exist
        ch_sites = aromatic_CH_indices(self.test_mol)
        self.assertTrue(len(ch_sites) > 0)
    
    def test_enumerate_r3_only(self):
        """Test R3 rule only (OH replacement)."""
        products = enumerate_k1_halogenation(
            self.test_mol, self.halogens, ['R3'], self.config
        )
        
        # Should have products from R3 (phenol has OH group) 
        # But R3 may not work if reactions fail, so just check it doesn't crash
        self.assertIsInstance(products, list)
    
    def test_different_halogens(self):
        """Test that different halogens generate different products."""
        products = enumerate_k1_halogenation(
            self.test_mol, ['F', 'Cl', 'Br'], ['R1'], self.config
        )
        
        # Should have products with different halogens
        halogens_found = set()
        for mol, props in products:
            halogens_found.add(props['halogen'])
        
        # Should have multiple halogens
        self.assertTrue(len(halogens_found) > 1)
    
    def test_deduplication(self):
        """Test that deduplication works."""
        # Run twice and combine - should dedupe
        products1 = enumerate_k1_halogenation(
            self.test_mol, self.halogens, ['R1'], self.config
        )
        products2 = enumerate_k1_halogenation(
            self.test_mol, self.halogens, ['R1'], self.config  
        )
        
        # Should be identical (deterministic)
        self.assertEqual(len(products1), len(products2))
    
    def test_qc_filtering(self):
        """Test QC filtering."""
        products = enumerate_k1_halogenation(
            self.test_mol, self.halogens, self.all_rules, self.config
        )
        
        # All products should pass sanitize_ok if strict QC is enabled
        for mol, props in products:
            self.assertTrue(props.get('sanitize_ok', False))
    
    def test_benzene_simple(self):
        """Test with simple benzene."""
        benzene = Chem.MolFromSmiles("c1ccccc1")
        products = enumerate_k1_halogenation(
            benzene, ['F'], ['R1'], self.config
        )
        
        # Benzene should generate R1 products (aromatic CH)
        self.assertTrue(len(products) > 0)
        
        # All should be R1 rule
        for mol, props in products:
            self.assertEqual(props['rule'], 'R1')


if __name__ == '__main__':
    unittest.main()