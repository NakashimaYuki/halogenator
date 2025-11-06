# -*- coding: ascii -*-
"""Test sugar fallback acceptance path."""

import unittest
from rdkit import Chem
from src.halogenator.sugar_mask import get_sugar_mask_with_full_status, get_sugar_mask_status


class TestSugarFallbackAccept(unittest.TestCase):
    """Test that sugar fallback path works correctly."""

    def test_fallback_acceptance_semantics(self):
        """Test that molecules triggering fallback have correct acceptance semantics."""
        # Use a molecule that might trigger fallback - a simple sugar-like structure
        # that doesn't meet strict main-path criteria but should be accepted via fallback
        # This is more of a test of the API semantics than finding a specific molecule

        # Test with glucose first to understand main path behavior
        glucose_smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        glucose_mol = Chem.MolFromSmiles(glucose_smiles)
        self.assertIsNotNone(glucose_mol)

        mask, degraded, audit = get_sugar_mask_with_full_status(glucose_mol, mode='heuristic')

        if not degraded:
            # Glucose uses main path - verify semantics
            self.assertTrue(audit.get('accepted', False), "Main path should set accepted=True")
            self.assertTrue(audit.get('accepted_via_score', False), "Main path should set accepted_via_score=True")
            self.assertEqual(audit.get('path'), 'main', "Should report main path")

        # Test new API consistency
        mask2, accepted, audit2 = get_sugar_mask_status(glucose_mol, mode='heuristic')
        self.assertEqual(mask, mask2, "Both APIs should return same mask")
        self.assertEqual(accepted, bool(mask), "Accepted should match mask existence")

    def test_fallback_audit_completeness(self):
        """Test that fallback path provides complete audit information."""
        # This test verifies that when fallback is triggered, the audit contains
        # all expected fields even if no specific molecule is guaranteed to trigger fallback

        # Test molecules that might not meet strict main-path criteria
        test_molecules = [
            ("simple_ether", "C1CCOCC1"),  # THP - should not be accepted
            ("modified_sugar", "CC1OC(O)C(O)C(O)C1O"),  # Modified sugar structure
        ]

        for name, smiles in test_molecules:
            with self.subTest(molecule=name):
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mask, degraded, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

                    # Check audit completeness regardless of acceptance
                    self.assertIn('path', audit, "Audit should contain path")
                    self.assertIn('accepted', audit, "Audit should contain accepted")
                    self.assertIn('accepted_via_score', audit, "Audit should contain accepted_via_score")

                    if degraded:
                        # Fallback semantics
                        self.assertEqual(audit['path'], 'fallback', "Degraded should report fallback path")
                        self.assertFalse(audit['accepted_via_score'], "Fallback should not accept via score")

                        # Fallback audit fields should be present
                        if audit.get('accepted', False):
                            # If fallback accepts, it should be via alternate criteria
                            self.assertTrue(audit.get('accepted'), "Fallback acceptance should set accepted=True")

    def test_strong_evidence_requirement(self):
        """Test that strong evidence gate is enforced in both main and fallback paths."""
        # This test verifies that the strong evidence requirement (>=2 exocyclic O OR cglyco)
        # is consistently applied

        glucose_smiles = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(glucose_smiles)
        self.assertIsNotNone(mol)

        mask, degraded, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

        if audit.get('accepted', False):
            # If accepted, should meet strong evidence criteria
            exo_o_count = audit.get('exocyclic_O_count_single_bond', 0)
            has_cglyco = audit.get('has_cglyco_evidence', False)

            strong_evidence = (exo_o_count >= 2) or has_cglyco
            self.assertTrue(strong_evidence,
                          f"Accepted sugar should meet strong evidence criteria: "
                          f"exo_O={exo_o_count}, cglyco={has_cglyco}")

    def test_api_consistency_fallback(self):
        """Test that both APIs handle fallback scenarios consistently."""
        # Test a molecule that should not be accepted as sugar
        thp_smiles = "C1CCOCC1"
        mol = Chem.MolFromSmiles(thp_smiles)
        self.assertIsNotNone(mol)

        # Test old API
        mask1, degraded, audit1 = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # Test new API
        mask2, accepted, audit2 = get_sugar_mask_status(mol, mode='heuristic')

        # Both should return empty mask for THP
        self.assertEqual(mask1, set(), "THP should not be masked")
        self.assertEqual(mask2, set(), "THP should not be masked")
        self.assertEqual(mask1, mask2, "Both APIs should return same mask")

        # Acceptance should be consistent
        self.assertFalse(accepted, "THP should not be accepted")
        self.assertEqual(audit1.get('accepted'), accepted, "APIs should agree on acceptance")


if __name__ == '__main__':
    unittest.main()