# -*- coding: ascii -*-
"""Test configurable beta-to-carbonyl semantics for R2b site detection."""

import unittest
from rdkit import Chem
from src.halogenator.sites import c_ring_sp3_CH2_flavanone_sites, _is_beta_to_carbonyl


class TestBetaCarbonylConfig(unittest.TestCase):
    """Test configurable beta-to-carbonyl distance semantics."""

    def setUp(self):
        """Set up test molecules."""
        # Flavanone - should have beta-to-carbonyl sites
        self.flavanone_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"
        self.flavanone_mol = Chem.MolFromSmiles(self.flavanone_smiles)
        self.assertIsNotNone(self.flavanone_mol)

    def test_strict_beta_semantics_default(self):
        """Test default strict beta=2 bonds semantics."""
        # Default behavior: no rules_cfg means strict beta=2 only
        sites_default = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=None
        )

        # Explicit strict configuration
        rules_cfg_strict = {'R2': {'allow_alpha_as_beta': False}}
        sites_strict = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=rules_cfg_strict
        )

        # Both should yield the same results
        self.assertEqual(sites_default, sites_strict,
                        "Default and explicit strict config should yield same results")

        # With strict beta=2 semantics, flavanone C3 (alpha position) should NOT be detected
        self.assertEqual(len(sites_default), 0,
                        "Flavanone should have no strict beta-to-carbonyl sites (C3 is alpha)")

        # But should find sites with alpha-as-beta compatibility
        rules_cfg_compat = {'R2': {'allow_alpha_as_beta': True}}
        sites_compat = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=rules_cfg_compat
        )
        self.assertGreater(len(sites_compat), 0,
                          "Flavanone should have alpha-to-carbonyl sites with compatibility mode")

    def test_alpha_as_beta_compatibility(self):
        """Test alpha-as-beta compatibility mode."""
        # Enable alpha-as-beta compatibility
        rules_cfg_compat = {'R2': {'allow_alpha_as_beta': True}}
        sites_compat = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=rules_cfg_compat
        )

        # Strict mode
        sites_strict = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=None
        )

        # Compatibility mode should find >= strict mode sites
        self.assertGreaterEqual(len(sites_compat), len(sites_strict),
                               "Alpha-as-beta mode should find at least as many sites as strict mode")

    def test_beta_carbonyl_helper_function(self):
        """Test the _is_beta_to_carbonyl helper function directly."""
        # Find a CH2 carbon in the flavanone
        ch2_carbons = []
        for i in range(self.flavanone_mol.GetNumAtoms()):
            atom = self.flavanone_mol.GetAtomWithIdx(i)
            if (atom.GetSymbol() == 'C' and
                atom.GetHybridization() == Chem.HybridizationType.SP3 and
                atom.GetTotalNumHs() == 2):
                ch2_carbons.append(i)

        self.assertGreater(len(ch2_carbons), 0, "Flavanone should have CH2 carbons")

        # Test strict beta detection
        found_beta_strict = False
        for idx in ch2_carbons:
            if _is_beta_to_carbonyl(self.flavanone_mol, idx, allow_alpha_as_beta=False):
                found_beta_strict = True
                break

        # Test alpha-as-beta detection
        found_beta_compat = False
        for idx in ch2_carbons:
            if _is_beta_to_carbonyl(self.flavanone_mol, idx, allow_alpha_as_beta=True):
                found_beta_compat = True
                break

        self.assertTrue(found_beta_strict or found_beta_compat,
                       "Should find beta-to-carbonyl carbons in flavanone")

    def test_configuration_independence(self):
        """Test that different molecules respond independently to configuration."""
        # Test with a simple non-flavanone molecule
        benzene_smiles = "c1ccccc1"
        benzene_mol = Chem.MolFromSmiles(benzene_smiles)
        self.assertIsNotNone(benzene_mol)

        # Benzene should not have beta-to-carbonyl sites regardless of config
        sites_strict = c_ring_sp3_CH2_flavanone_sites(
            benzene_mol, set(), sugar_cfg=None, rules_cfg=None
        )

        rules_cfg_compat = {'R2': {'allow_alpha_as_beta': True}}
        sites_compat = c_ring_sp3_CH2_flavanone_sites(
            benzene_mol, set(), sugar_cfg=None, rules_cfg=rules_cfg_compat
        )

        # Both should find no sites (benzene has no sp3 CH2 carbons)
        self.assertEqual(len(sites_strict), 0, "Benzene should have no R2b sites")
        self.assertEqual(len(sites_compat), 0, "Benzene should have no R2b sites")
        self.assertEqual(sites_strict, sites_compat, "Config should not affect non-applicable molecules")

    def test_malformed_config_handling(self):
        """Test that malformed configurations are handled gracefully."""
        # Empty rules_cfg
        sites_empty = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg={}
        )

        # Malformed R2 section
        rules_cfg_malformed = {'R2': 'not_a_dict'}
        sites_malformed = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=rules_cfg_malformed
        )

        # Default behavior
        sites_default = c_ring_sp3_CH2_flavanone_sites(
            self.flavanone_mol, set(), sugar_cfg=None, rules_cfg=None
        )

        # All should default to strict behavior
        self.assertEqual(sites_empty, sites_default, "Empty config should default to strict")
        self.assertEqual(sites_malformed, sites_default, "Malformed config should default to strict")


if __name__ == '__main__':
    unittest.main()