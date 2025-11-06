# -*- coding: ascii -*-
"""Tests for R6 halogen whitelist enforcement."""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.rules_methyl import validate_methyl_halogen, apply_methyl_step, apply_methyl_macro


class TestR6HalogenWhitelist(unittest.TestCase):
    """Test that R6 methyl halogenation only allows F and Cl (no Br/I)."""

    def test_validate_methyl_halogen_allowed(self):
        """Test that F and Cl are allowed by validation function."""
        self.assertTrue(validate_methyl_halogen('F'))
        self.assertTrue(validate_methyl_halogen('Cl'))

    def test_validate_methyl_halogen_blocked(self):
        """Test that Br and I are blocked by validation function."""
        self.assertFalse(validate_methyl_halogen('Br'))
        self.assertFalse(validate_methyl_halogen('I'))

    def test_validate_methyl_halogen_invalid(self):
        """Test that invalid halogen symbols are blocked."""
        self.assertFalse(validate_methyl_halogen('X'))
        self.assertFalse(validate_methyl_halogen(''))
        self.assertFalse(validate_methyl_halogen('At'))
        self.assertFalse(validate_methyl_halogen('Ts'))

    def test_apply_methyl_step_allowed_halogens(self):
        """Test that step halogenation works for F and Cl."""
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit not available")

        mol = Chem.MolFromSmiles('CC')  # Ethane
        if mol is None:
            self.skipTest("Failed to create test molecule")

        # F should work
        result_f = apply_methyl_step(mol, 0, 'F')
        self.assertIsNotNone(result_f, "F halogenation should succeed")

        # Cl should work
        result_cl = apply_methyl_step(mol, 0, 'Cl')
        self.assertIsNotNone(result_cl, "Cl halogenation should succeed")

    def test_apply_methyl_step_blocked_halogens(self):
        """Test that step halogenation is blocked for Br and I."""
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit not available")

        mol = Chem.MolFromSmiles('CC')  # Ethane
        if mol is None:
            self.skipTest("Failed to create test molecule")

        # Br should be blocked
        result_br = apply_methyl_step(mol, 0, 'Br')
        self.assertIsNone(result_br, "Br halogenation should be blocked")

        # I should be blocked
        result_i = apply_methyl_step(mol, 0, 'I')
        self.assertIsNone(result_i, "I halogenation should be blocked")

    def test_apply_methyl_macro_allowed_labels(self):
        """Test that macro halogenation works for CF3 and CCl3."""
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit not available")

        mol = Chem.MolFromSmiles('CC')  # Ethane
        if mol is None:
            self.skipTest("Failed to create test molecule")

        # CF3 should work
        result_cf3 = apply_methyl_macro(mol, 0, 'CF3')
        self.assertIsNotNone(result_cf3, "CF3 macro halogenation should succeed")

        # CCl3 should work
        result_ccl3 = apply_methyl_macro(mol, 0, 'CCl3')
        self.assertIsNotNone(result_ccl3, "CCl3 macro halogenation should succeed")

    def test_apply_methyl_macro_blocked_labels(self):
        """Test that macro halogenation is blocked for Br and I containing labels."""
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit not available")

        mol = Chem.MolFromSmiles('CC')  # Ethane
        if mol is None:
            self.skipTest("Failed to create test molecule")

        # CBr3 should be blocked
        result_cbr3 = apply_methyl_macro(mol, 0, 'CBr3')
        self.assertIsNone(result_cbr3, "CBr3 macro halogenation should be blocked")

        # CI3 should be blocked
        result_ci3 = apply_methyl_macro(mol, 0, 'CI3')
        self.assertIsNone(result_ci3, "CI3 macro halogenation should be blocked")

        # Mixed labels should be blocked
        result_mixed = apply_methyl_macro(mol, 0, 'CFBr2')
        self.assertIsNone(result_mixed, "Mixed halogen macro labels should be blocked")

    def test_whitelist_constants_definition(self):
        """Test that the whitelist constants are defined correctly."""
        from halogenator.rules_methyl import _ALLOWED_STEP, _ALLOWED_MACRO

        # Test step whitelist
        self.assertEqual(_ALLOWED_STEP, {"F", "Cl"})
        self.assertNotIn("Br", _ALLOWED_STEP)
        self.assertNotIn("I", _ALLOWED_STEP)

        # Test macro whitelist
        self.assertEqual(_ALLOWED_MACRO, {"CF3", "CCl3"})
        self.assertNotIn("CBr3", _ALLOWED_MACRO)
        self.assertNotIn("CI3", _ALLOWED_MACRO)


if __name__ == '__main__':
    unittest.main()