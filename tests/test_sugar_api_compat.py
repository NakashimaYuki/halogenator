# -*- coding: ascii -*-
"""Test backward compatibility of sugar masking API."""

import unittest
from rdkit import Chem
from src.halogenator.sugar_mask import get_sugar_mask_with_full_status, get_sugar_mask_status


class TestSugarAPICompat(unittest.TestCase):
    """Test that sugar masking API maintains backward compatibility."""

    def test_abi_compat_glucose(self):
        """Test API compatibility with glucose example."""
        smi = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(mol)

        # Test old API returns (mask, degraded, audit)
        mask1, degraded, audit1 = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # Test new API returns (mask, accepted, audit)
        mask2, accepted, audit2 = get_sugar_mask_status(mol, mode='heuristic')

        # Both should return the same mask
        self.assertEqual(mask1, mask2, "Both APIs should return identical masks")

        # Acceptance should match mask existence
        self.assertEqual(accepted, bool(mask1), "Accepted should match mask existence")

        # Both audits should contain required fields
        required_fields = ['accepted', 'accepted_via_score', 'path']
        for field in required_fields:
            self.assertIn(field, audit1, f"Old API audit should contain {field}")
            self.assertIn(field, audit2, f"New API audit should contain {field}")

        # Field values should be consistent between APIs
        self.assertEqual(audit1['accepted'], audit2['accepted'])
        self.assertEqual(audit1['accepted_via_score'], audit2['accepted_via_score'])
        self.assertEqual(audit1['path'], audit2['path'])

    def test_abi_compat_thp(self):
        """Test API compatibility with THP (should not be masked)."""
        smi = "C1CCOCC1"  # Simple THP
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(mol)

        # Test both APIs
        mask1, degraded, audit1 = get_sugar_mask_with_full_status(mol, mode='heuristic')
        mask2, accepted, audit2 = get_sugar_mask_status(mol, mode='heuristic')

        # Both should return empty mask
        self.assertEqual(mask1, set(), "THP should not be masked")
        self.assertEqual(mask2, set(), "THP should not be masked")

        # THP should not be accepted as sugar
        self.assertFalse(accepted, "THP should not be accepted as sugar")
        self.assertFalse(audit1['accepted'], "THP should not be accepted in old API")
        self.assertFalse(audit2['accepted'], "THP should not be accepted in new API")

    def test_audit_schema_consistency(self):
        """Test that audit schema is consistent between APIs."""
        smi = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smi)

        _, _, audit1 = get_sugar_mask_with_full_status(mol, mode='heuristic')
        _, _, audit2 = get_sugar_mask_status(mol, mode='heuristic')

        # Core fields should be present in both
        core_fields = ['accepted', 'accepted_via_score', 'path']
        for field in core_fields:
            self.assertIn(field, audit1)
            self.assertIn(field, audit2)
            self.assertEqual(audit1[field], audit2[field])

        # Expected fields for main path glucose
        if audit1.get('path') == 'main':
            expected_main_fields = [
                'ring_size', 'exocyclic_O_count_single_bond',
                'has_cglyco_evidence', 'ring_score'
            ]
            for field in expected_main_fields:
                self.assertIn(field, audit1, f"Main path audit should contain {field}")
                self.assertIn(field, audit2, f"Main path audit should contain {field}")


if __name__ == '__main__':
    unittest.main()