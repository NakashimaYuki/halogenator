# -*- coding: ascii -*-
"""Test reaction rules functionality."""

import unittest
from rdkit import Chem

from halogenator.rules import build_reactions, get_rule_description


class TestRules(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures."""
        self.reactions = build_reactions()
        
    def test_build_reactions(self):
        """Test reaction building returns correct structure."""
        self.assertIsInstance(self.reactions, dict)
        
        # Should have correct rule keys
        expected_rules = ['R1', 'R3', 'R4', 'R5']
        for rule in expected_rules:
            self.assertIn(rule, self.reactions)
            self.assertIsInstance(self.reactions[rule], dict)
            
            # Each rule should have reactions for each halogen - reactions must be non-None
            expected_halogens = ['F', 'Cl', 'Br', 'I']
            for halogen in expected_halogens:
                self.assertIn(halogen, self.reactions[rule])
                self.assertIsNotNone(self.reactions[rule][halogen], f"{rule}-{halogen} reaction is None")
    
    def test_r3_phenol(self):
        """Test R3 rule with phenol (replace -OH)."""
        phenol = Chem.MolFromSmiles("c1ccc(O)cc1")  # phenol
        
        # Test with Cl - reaction must exist
        reaction = self.reactions['R3']['Cl']
        self.assertIsNotNone(reaction, "R3-Cl reaction must not be None")
        
        products = reaction.RunReactants((phenol,))
        self.assertTrue(len(products) > 0, "R3 should generate products from phenol")
        
        # Check that product contains Cl
        for product_set in products:
            for product in product_set:
                if product is not None:
                    smiles = Chem.MolToSmiles(product)
                    self.assertIn('Cl', smiles, "R3 product should contain Cl")
                    break
    
    def test_r4_ethylamine(self):
        """Test R4 rule with ethylamine (replace -NH2)."""
        ethylamine = Chem.MolFromSmiles("CCN")  # ethylamine
        
        # Test with Br - reaction must exist
        reaction = self.reactions['R4']['Br']
        self.assertIsNotNone(reaction, "R4-Br reaction must not be None")
        
        products = reaction.RunReactants((ethylamine,))
        self.assertTrue(len(products) > 0, "R4 should generate products from ethylamine")
        
        # Check that product contains Br
        for product_set in products:
            for product in product_set:
                if product is not None:
                    smiles = Chem.MolToSmiles(product)
                    self.assertIn('Br', smiles, "R4 product should contain Br")
                    break
    
    def test_r5_benzoic_acid(self):
        """Test R5 rule with benzoic acid (replace -COOH)."""
        benzoic_acid = Chem.MolFromSmiles("c1ccc(C(=O)O)cc1")  # benzoic acid
        
        # Test with F - reaction must exist
        reaction = self.reactions['R5']['F']
        self.assertIsNotNone(reaction, "R5-F reaction must not be None")
        
        products = reaction.RunReactants((benzoic_acid,))
        self.assertTrue(len(products) > 0, "R5 should generate products from benzoic acid")
        
        # Check that product contains F
        for product_set in products:
            for product in product_set:
                if product is not None:
                    smiles = Chem.MolToSmiles(product)
                    self.assertIn('F', smiles, "R5 product should contain F")
                    break
    
    def test_get_rule_description(self):
        """Test rule descriptions."""
        self.assertEqual(get_rule_description('R1'), 'Aromatic CH -> C-halogen')
        self.assertEqual(get_rule_description('R3'), 'Replace -OH (not carboxylic acid)')
        self.assertEqual(get_rule_description('R4'), 'Replace -NHx (not amide)')
        self.assertEqual(get_rule_description('R5'), 'Replace -C(=O)OH (whole carboxyl)')
        self.assertEqual(get_rule_description('UNKNOWN'), 'Unknown rule')


if __name__ == '__main__':
    unittest.main()