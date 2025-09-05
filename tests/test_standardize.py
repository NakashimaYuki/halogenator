# -*- coding: ascii -*-
"""Test standardization functionality."""

import unittest
from rdkit import Chem

from halogenator.standardize import std_from_smiles, to_inchikey


class TestStandardize(unittest.TestCase):
    
    def test_std_from_smiles_valid(self):
        """Test standardization of valid SMILES."""
        # Simple benzene
        mol = std_from_smiles("c1ccccc1")
        self.assertIsNotNone(mol)
        
        # Phenol  
        mol = std_from_smiles("c1ccc(O)cc1")
        self.assertIsNotNone(mol)
        
        # Simple alcohol
        mol = std_from_smiles("CCO")
        self.assertIsNotNone(mol)
    
    def test_std_from_smiles_invalid(self):
        """Test standardization of invalid SMILES."""
        # Invalid SMILES
        mol = std_from_smiles("invalid_smiles")
        self.assertIsNone(mol)
        
        # None input
        mol = std_from_smiles(None)
        self.assertIsNone(mol)
    
    def test_to_inchikey(self):
        """Test InChIKey generation."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        key = to_inchikey(mol)
        self.assertIsInstance(key, str)
        self.assertTrue(len(key) > 0)
        
        # Should be deterministic
        key2 = to_inchikey(mol)
        self.assertEqual(key, key2)
    
    def test_to_inchikey_fallback(self):
        """Test InChIKey fallback to SMILES."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        key = to_inchikey(mol)
        # Should be either InChIKey format or canonical SMILES
        self.assertIsInstance(key, str)
        self.assertTrue(len(key) > 0)


if __name__ == '__main__':
    unittest.main()