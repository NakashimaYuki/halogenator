# -*- coding: ascii -*-
"""Test sugar acceptance via main scoring path."""

import unittest
from rdkit import Chem
from src.halogenator.sugar_mask import get_sugar_mask_status


class TestSugarAcceptMainPath(unittest.TestCase):
    """Test that canonical sugars are accepted via main scoring with strong evidence."""

    def test_glucose_main_accept(self):
        """Test that D-glucopyranose is accepted via main scoring path."""
        # Simple D-glucopyranose; choose a canonical SMILES that RDKit parses robustly
        smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smiles)
        self.assertIsNotNone(mol, "Glucose SMILES should parse correctly")

        mask, accepted, audit = get_sugar_mask_status(mol, mode='heuristic')

        # Must be accepted through main scoring with strong evidence
        self.assertTrue(accepted, "Glucose should be accepted as sugar")

        # accepted_via_score should be True for main-path acceptance
        self.assertTrue(audit.get('accepted_via_score', False),
                        "Glucose should be accepted via main scoring path")

        # Exocyclic O count should be >= 2
        self.assertGreaterEqual(audit.get('exocyclic_O_count_single_bond', 0), 2,
                                "Glucose should have >=2 exocyclic single-bond oxygens")

        self.assertGreater(len(mask), 0, "Glucose ring atoms should be masked")

    def test_common_sugar_acceptance(self):
        """Test that other common sugars are also accepted via main path."""
        # Test additional common sugar structures
        test_sugars = [
            ("fructose", "OC[C@@H]1O[C@@H](CO)[C@H](O)[C@H]1O"),  # D-fructofuranose (5-ring)
            ("galactose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),  # D-galactopyranose (6-ring)
            ("ribose", "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O"),  # D-ribofuranose (5-ring)
        ]

        for sugar_name, smiles in test_sugars:
            with self.subTest(sugar=sugar_name):
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:  # Skip if SMILES doesn't parse
                    mask, accepted, audit = get_sugar_mask_status(mol, mode='heuristic')

                    # These should be accepted via main scoring
                    self.assertTrue(accepted, f"{sugar_name} should be accepted as sugar")
                    self.assertTrue(audit.get('accepted_via_score', False),
                                    f"{sugar_name} should be accepted via main scoring path")


if __name__ == '__main__':
    unittest.main()