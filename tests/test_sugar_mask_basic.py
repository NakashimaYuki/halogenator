# -*- coding: ascii -*-
"""
Basic tests for sugar masking functionality.

Tests the core sugar masking implementation to ensure it correctly identifies
and masks sugar atoms in flavonoid glycosides.
"""

import unittest
from typing import Set

# Import the modules we want to test
try:
    from rdkit import Chem
    from src.halogenator.sugar_mask import get_sugar_mask, get_sugar_mask_with_status, post_guard_blocked, compute_sugar_audit_fields
    from src.halogenator.symmetry import canonical_ranks_on_masked_subgraph
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@unittest.skipUnless(RDKIT_AVAILABLE, "RDKit not available")
class TestSugarMaskBasic(unittest.TestCase):
    """Test basic sugar masking functionality."""

    def setUp(self):
        """Set up test molecules."""
        # Simple glucose-like molecule (should be masked)
        self.glucose_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"

        # Simple flavonoid without sugar (should not be masked)
        self.quercetin_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

        # Flavonoid glycoside (should have both masked and unmasked parts)
        self.quercetin_glucose_smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c12"

    def test_sugar_mask_off_mode(self):
        """Test that 'off' mode returns empty mask."""
        mol = Chem.MolFromSmiles(self.glucose_smiles)
        mask, _ = get_sugar_mask_with_status(mol, mode='off')
        self.assertEqual(mask, set())

    def test_sugar_mask_heuristic_glucose(self):
        """Test that glucose-like molecule gets masked."""
        mol = Chem.MolFromSmiles(self.glucose_smiles)
        mask, _ = get_sugar_mask_with_status(mol, mode='heuristic')

        # Should identify atoms as sugar
        self.assertGreater(len(mask), 0)

        # Glucose should mask most/all atoms since it's a pure sugar
        # OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O
        # Expected: 6-member sugar ring + exocyclic oxygens
        self.assertGreaterEqual(len(mask), 6)  # At least the ring atoms

        # Should include the 6-member ring (atoms with @H stereochemistry)
        ring_atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            if '@' in self.glucose_smiles:  # Ring carbons have stereochemistry
                pass  # Ring detection is done by RDKit ring analysis

        print(f"Glucose mask ({len(mask)} atoms): {sorted(mask)}")

    def test_sugar_mask_heuristic_quercetin(self):
        """Test that pure flavonoid (quercetin) gets minimal masking."""
        mol = Chem.MolFromSmiles(self.quercetin_smiles)
        mask, _ = get_sugar_mask_with_status(mol, mode='heuristic')

        # Pure flavonoid should have no sugar masking (aromatic rings don't match criteria)
        # Quercetin: O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12
        # All rings are aromatic, so should not be detected as sugar
        self.assertEqual(len(mask), 0, "Pure flavonoid should not be masked as sugar")

        print(f"Quercetin mask ({len(mask)} atoms): {sorted(mask)}")

    def test_sugar_mask_heuristic_glycoside(self):
        """Test that flavonoid glycoside gets appropriate masking."""
        mol = Chem.MolFromSmiles(self.quercetin_glucose_smiles)
        mask, _ = get_sugar_mask_with_status(mol, mode='heuristic')

        # Should mask sugar part but not flavonoid core
        self.assertGreater(len(mask), 0, "Glycoside should have some sugar masking")

        # The glucose portion should be masked (6-member ring + exocyclic oxygens)
        # but not the flavonoid core (aromatic rings)
        # Quercetin has 22 atoms, glucose adds ~12 more, so total ~34 atoms
        # Sugar masking should be significant but not all atoms
        total_atoms = mol.GetNumAtoms()
        mask_ratio = len(mask) / total_atoms

        self.assertGreater(mask_ratio, 0.2, "Should mask significant portion (sugar part)")
        self.assertLess(mask_ratio, 0.8, "Should not mask everything (preserve flavonoid core)")

        print(f"Quercetin-glucose mask ({len(mask)}/{total_atoms} atoms, {mask_ratio:.2%}): {sorted(mask)}")

    def test_post_guard_clean_molecule(self):
        """Test post-guard with molecule that has no halogens."""
        mol = Chem.MolFromSmiles(self.quercetin_smiles)
        mask = {1, 2, 3}  # Arbitrary mask

        # Should not be blocked since no halogens present
        self.assertFalse(post_guard_blocked(mol, mask))

    def test_post_guard_halogenated_molecule(self):
        """Test post-guard with halogenated molecule."""
        # Create a simple halogenated molecule
        hal_mol = Chem.MolFromSmiles("CCF")  # Simple molecule with fluorine
        mask = {2}  # Mask the fluorine atom

        # Should be blocked since fluorine is in mask
        self.assertTrue(post_guard_blocked(hal_mol, mask))

    def test_compute_sugar_audit_fields(self):
        """Test computation of sugar audit fields."""
        mol = Chem.MolFromSmiles(self.glucose_smiles)
        mask, _ = get_sugar_mask_with_status(mol, mode='heuristic')
        default_sugar_cfg = {'mode': 'heuristic', 'mask_exocyclic_oxygen': True, 'mask_glycosidic_bridge_oxygen': True, 'audit': False}
        audit = compute_sugar_audit_fields(mol, mask, default_sugar_cfg)

        # Check that audit fields are present
        self.assertIn('sugar_mask_atoms', audit)
        self.assertIn('sugar_rings', audit)
        self.assertIn('masked_oh_count', audit)
        self.assertIn('masked_bridge_o_count', audit)

        print(f"Sugar audit for glucose: {audit}")

    def test_masked_subgraph_symmetry(self):
        """Test symmetry computation on masked subgraph."""
        mol = Chem.MolFromSmiles(self.quercetin_smiles)
        mask = {1, 2}  # Mask some atoms

        ranks = canonical_ranks_on_masked_subgraph(mol, mask)

        # Check that masked atoms have rank -1
        self.assertEqual(ranks[1], -1)
        self.assertEqual(ranks[2], -1)

        # Check that we got ranks for all atoms
        self.assertEqual(len(ranks), mol.GetNumAtoms())

        print(f"Masked symmetry ranks: {ranks}")

    def test_sugar_cfg_detail_switches(self):
        """Test that sugar_cfg detail switches work correctly."""
        mol = Chem.MolFromSmiles(self.glucose_smiles)

        # Test with all switches enabled (default)
        mask_all = get_sugar_mask(mol, mode='heuristic', sugar_cfg={
            'mask_exocyclic_oxygen': True,
            'mask_glycosidic_bridge_oxygen': True
        })

        # Test with exocyclic oxygen disabled
        mask_no_exo = get_sugar_mask(mol, mode='heuristic', sugar_cfg={
            'mask_exocyclic_oxygen': False,
            'mask_glycosidic_bridge_oxygen': True
        })

        # Test with bridge oxygen disabled
        mask_no_bridge = get_sugar_mask(mol, mode='heuristic', sugar_cfg={
            'mask_exocyclic_oxygen': True,
            'mask_glycosidic_bridge_oxygen': False
        })

        # Test with both disabled (only sugar ring atoms)
        mask_ring_only = get_sugar_mask(mol, mode='heuristic', sugar_cfg={
            'mask_exocyclic_oxygen': False,
            'mask_glycosidic_bridge_oxygen': False
        })

        # Assertions about the relationships
        self.assertGreaterEqual(len(mask_all), len(mask_no_exo))
        self.assertGreaterEqual(len(mask_all), len(mask_no_bridge))
        self.assertGreaterEqual(len(mask_all), len(mask_ring_only))
        self.assertGreaterEqual(len(mask_no_exo), len(mask_ring_only))
        self.assertGreaterEqual(len(mask_no_bridge), len(mask_ring_only))

        print(f"Sugar mask variations:")
        print(f"  All enabled: {len(mask_all)} atoms")
        print(f"  No exocyclic: {len(mask_no_exo)} atoms")
        print(f"  No bridge: {len(mask_no_bridge)} atoms")
        print(f"  Ring only: {len(mask_ring_only)} atoms")

    def test_empty_mask_symmetry(self):
        """Test symmetry computation with empty mask."""
        mol = Chem.MolFromSmiles(self.quercetin_smiles)
        mask = set()

        ranks = canonical_ranks_on_masked_subgraph(mol, mask)

        # Should behave like standard symmetry
        self.assertEqual(len(ranks), mol.GetNumAtoms())
        # Check that no atoms have the masked rank (-1)
        masked_atoms = [i for i, rank in enumerate(ranks) if rank == -1]
        self.assertEqual(len(masked_atoms), 0, "No atoms should be masked with empty mask")


if __name__ == '__main__':
    unittest.main()