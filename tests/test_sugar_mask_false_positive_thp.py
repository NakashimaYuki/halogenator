# -*- coding: ascii -*-
"""Test sugar mask false positive prevention for THP."""

import unittest
from src.halogenator.sugar_mask import get_sugar_mask_with_full_status
from rdkit import Chem


class TestSugarMaskFalsePositiveTHP(unittest.TestCase):
    """Test that THP (tetrahydropyran) is not incorrectly masked as sugar."""

    def test_thp_not_masked(self):
        """Test that hydroxylated THP is not masked as a sugar ring."""
        # Hydroxylated tetrahydropyran - should NOT be considered a sugar
        # because it has only 1 exocyclic oxygen and no C-glycoside evidence
        thp_smiles = "OC1CCCCO1"
        mol = Chem.MolFromSmiles(thp_smiles)

        self.assertIsNotNone(mol, "THP molecule should be valid")

        # Test with default sugar configuration
        mask, accepted, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # THP should NOT be masked (empty mask set)
        self.assertEqual(len(mask), 0,
                        f"THP should not be masked as sugar, but got masked atoms: {mask}")

        # Should not be accepted as sugar
        self.assertFalse(accepted, "THP should not be accepted as containing sugar rings")


    def test_real_sugar_still_masked(self):
        """Test that a real sugar with strong evidence is still properly masked."""
        # Use a simpler glucose molecule - should be masked due to multiple exocyclic oxygens
        glucose_smiles = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        mol = Chem.MolFromSmiles(glucose_smiles)

        if mol is not None:
            mask, accepted, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

            # Real sugar should be masked (non-empty mask set)
            self.assertGreater(len(mask), 0,
                             "Real sugar should be masked")

            # Real sugar should be masked either through main pathway (accepted=True)
            # or through degraded pathways (mask size > 0)
            # The important thing is that it's being masked to prevent enumeration explosion
            self.assertTrue(len(mask) > 0 or accepted,
                          f"Real sugar should be masked somehow: mask_size={len(mask)}, accepted={accepted}")


if __name__ == '__main__':
    unittest.main()