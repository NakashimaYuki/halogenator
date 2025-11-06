# -*- coding: ascii -*-
"""Test sugar audit schema consistency."""

import unittest
from rdkit import Chem
from src.halogenator.sugar_mask import get_sugar_mask_with_full_status, get_sugar_mask_status
from src.halogenator.schema import CORE_SUGAR_AUDIT_KEYS, FALLBACK_SUGAR_AUDIT_KEYS, SUGAR_AUDIT_KEYS


class TestSugarAuditSchema(unittest.TestCase):
    """Test that sugar audit schema is consistent and complete."""

    def test_main_path_audit_schema(self):
        """Test that main path audit contains expected core fields."""
        # Use glucose which should be accepted via main path
        smi = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(mol)

        mask, degraded, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # Should be main path (not degraded)
        self.assertFalse(degraded, "Glucose should use main path")
        self.assertTrue(audit.get('accepted', False), "Glucose should be accepted")
        self.assertEqual(audit.get('path'), 'main', "Should be main path")

        # Check that all core audit keys are present
        for key in CORE_SUGAR_AUDIT_KEYS:
            self.assertIn(key, audit, f"Core audit key '{key}' should be present in main path")

        # Validate key types and values
        self.assertIsInstance(audit['path'], str)
        self.assertIsInstance(audit['accepted'], bool)
        self.assertIsInstance(audit['accepted_via_score'], bool)
        self.assertIsInstance(audit['ring_size'], int)
        self.assertIsInstance(audit['exocyclic_O_count_single_bond'], int)
        self.assertIsInstance(audit['has_cglyco_evidence'], bool)
        self.assertIsInstance(audit['ring_score'], (int, float))

        # Main path specific validations
        self.assertTrue(audit['accepted_via_score'], "Main path should have accepted_via_score=True")
        self.assertGreater(audit['ring_size'], 0, "Should have positive ring size")
        self.assertGreaterEqual(audit['exocyclic_O_count_single_bond'], 2, "Glucose should have >=2 exocyclic O")
        self.assertGreater(audit['ring_score'], 0, "Should have positive ring score")

    def test_thp_no_acceptance_audit_schema(self):
        """Test that non-sugar molecules have consistent audit schema."""
        # Use THP which should not be accepted
        smi = "C1CCOCC1"
        mol = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(mol)

        mask, degraded, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # Should not be accepted
        self.assertFalse(audit.get('accepted', True), "THP should not be accepted")
        self.assertEqual(audit.get('path'), 'main', "Should still report main path attempt")
        self.assertFalse(audit.get('accepted_via_score', True), "Should not be accepted via score")

        # Core fields should still be present with appropriate values
        expected_core_fields = ['path', 'accepted', 'accepted_via_score']
        for key in expected_core_fields:
            self.assertIn(key, audit, f"Core field '{key}' should be present even for non-sugar")

    def test_new_api_audit_consistency(self):
        """Test that new API returns consistent audit schema."""
        smi = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smi)

        mask, accepted, audit = get_sugar_mask_status(mol, mode='heuristic')

        # Should contain core audit keys
        for key in CORE_SUGAR_AUDIT_KEYS:
            if key in ('degraded_reason', 'fallback_ring_count'):
                # These are fallback-specific, may not be present in main path
                continue
            self.assertIn(key, audit, f"New API audit should contain core key '{key}'")

        # Acceptance consistency
        self.assertTrue(accepted, "Glucose should be accepted")
        self.assertEqual(audit['accepted'], accepted, "Audit accepted should match return value")

    def test_audit_key_whitelist_compliance(self):
        """Test that audit only contains whitelisted keys."""
        smi = "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"
        mol = Chem.MolFromSmiles(smi)

        mask, degraded, audit = get_sugar_mask_with_full_status(mol, mode='heuristic')

        # All audit keys should be in the schema whitelist
        audit_keys = set(audit.keys())
        allowed_keys = set(SUGAR_AUDIT_KEYS)

        unexpected_keys = audit_keys - allowed_keys
        self.assertEqual(len(unexpected_keys), 0,
                        f"Audit contains unexpected keys not in schema: {unexpected_keys}")

    def test_schema_constants_non_empty(self):
        """Test that schema constants are properly defined."""
        self.assertGreater(len(CORE_SUGAR_AUDIT_KEYS), 0, "Core audit keys should not be empty")
        self.assertGreater(len(FALLBACK_SUGAR_AUDIT_KEYS), 0, "Fallback audit keys should not be empty")
        self.assertGreater(len(SUGAR_AUDIT_KEYS), 0, "Combined audit keys should not be empty")

        # Check no duplicates
        self.assertEqual(len(SUGAR_AUDIT_KEYS), len(set(SUGAR_AUDIT_KEYS)),
                        "Sugar audit keys should not contain duplicates")

        # Check core keys are included in combined
        for key in CORE_SUGAR_AUDIT_KEYS:
            self.assertIn(key, SUGAR_AUDIT_KEYS, f"Core key '{key}' should be in combined keys")


if __name__ == '__main__':
    unittest.main()