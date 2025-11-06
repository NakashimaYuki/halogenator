# -*- coding: ascii -*-
"""
Test R3/R4/R5 reaction-level sugar mask filtering.

This test ensures that R3/R4/R5 rules properly filter out reaction matches
that hit sugar-masked atoms and correctly count them in QA statistics.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
from halogenator.sugar_mask import get_sugar_mask_with_status
from halogenator.chem_compat import Chem


class TestR345SugarMaskFiltering(unittest.TestCase):
    """Test R3/R4/R5 reaction matching level sugar mask filtering."""

    def test_r3_sugar_mask_filtering_with_quercetin_glycoside(self):
        """
        Test that R3 (carboxylic acid halogenation) properly filters matches
        that hit sugar-masked atoms and records them in QA statistics.
        """
        # Use quercetin-3-glucoside (sugar attached to quercetin)
        # This molecule has a glucoside sugar group that should be masked
        quercetin_glucoside = "O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

        mol = Chem.MolFromSmiles(quercetin_glucoside)
        self.assertIsNotNone(mol, "Quercetin glycoside SMILES should be valid")

        # Get sugar mask - should mask the glucose portion
        sugar_mask, degraded = get_sugar_mask_with_status(mol, mode='heuristic')
        self.assertFalse(degraded, "Sugar mask should not be degraded for this molecule")
        self.assertGreater(len(sugar_mask), 0, "Sugar mask should identify sugar atoms")

        # Test with R3 only to isolate the effect
        config = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3',),  # Only carboxylic acid halogenation
            sugar_cfg={'mode': 'heuristic', 'audit': True}
        )

        # First, test without sugar masking to establish baseline
        config_no_mask = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R3',),
            sugar_cfg={'mode': 'off'}  # No sugar masking
        )

        products_no_mask, qa_stats_no_mask = enumerate_with_stats(quercetin_glucoside, config_no_mask)
        products_with_mask, qa_stats_with_mask = enumerate_with_stats(quercetin_glucoside, config)

        # There should be fewer products with masking than without
        # (assuming there are R3 sites in the sugar region)
        self.assertGreaterEqual(len(products_no_mask), len(products_with_mask),
                               "Sugar masking should reduce or maintain product count")

        # Check QA statistics for sugar filtering
        qa_paths_with_mask = qa_stats_with_mask.get('qa_paths', {})
        sugar_filtered_count = qa_paths_with_mask.get('sugar_mask_filtered', 0)

        # If there were R3 sites in the sugar region, they should be filtered
        if len(products_no_mask) > len(products_with_mask):
            self.assertGreater(sugar_filtered_count, 0,
                             "Should have sugar_mask_filtered count when products are reduced")

        # Verify that sugar_mask_filtered is correctly tracked
        self.assertIn('sugar_mask_filtered', qa_paths_with_mask,
                     "sugar_mask_filtered should be present in qa_paths")

    def test_r4_sugar_mask_filtering_basic(self):
        """
        Test that R4 (amine halogenation) properly handles sugar mask filtering.
        """
        # Use a simple glucose molecule which should be heavily masked
        glucose = "OCC1OC(O)C(O)C(O)C1O"

        mol = Chem.MolFromSmiles(glucose)
        self.assertIsNotNone(mol, "Glucose SMILES should be valid")

        # Get sugar mask - should mask most/all of the molecule
        sugar_mask, degraded = get_sugar_mask_with_status(mol, mode='heuristic')
        self.assertGreater(len(sugar_mask), 0, "Sugar mask should identify sugar atoms in glucose")

        # Test with R4 only
        config_with_mask = EnumConfig(
            k_max=1,
            halogens=('Cl',),
            rules=('R4',),  # Only amine halogenation
            sugar_cfg={'mode': 'heuristic', 'audit': True}
        )

        config_no_mask = EnumConfig(
            k_max=1,
            halogens=('Cl',),
            rules=('R4',),
            sugar_cfg={'mode': 'off'}
        )

        products_no_mask, qa_stats_no_mask = enumerate_with_stats(glucose, config_no_mask)
        products_with_mask, qa_stats_with_mask = enumerate_with_stats(glucose, config_with_mask)

        # With heavy sugar masking, there should be fewer (possibly zero) products
        self.assertLessEqual(len(products_with_mask), len(products_no_mask),
                            "Sugar masking should reduce product count for glucose")

        # Check that QA tracking is working
        qa_paths = qa_stats_with_mask.get('qa_paths', {})
        self.assertIn('sugar_mask_filtered', qa_paths,
                     "sugar_mask_filtered should be tracked in qa_paths")

    def test_sugar_mask_filtering_counts_properly(self):
        """
        Test that sugar mask filtering counts are accumulated correctly
        across multiple halogens and attempts.
        """
        # Use a molecule with potential R3/R4/R5 sites
        test_molecule = "OCC1OC(O)C(O)C(O)C1O"  # Glucose - should be heavily masked

        config = EnumConfig(
            k_max=1,
            halogens=('F', 'Cl'),  # Multiple halogens
            rules=('R3', 'R4', 'R5'),  # Multiple reaction rules
            sugar_cfg={'mode': 'heuristic', 'audit': True}
        )

        products, qa_stats = enumerate_with_stats(test_molecule, config)

        # Check that QA paths are properly structured
        qa_paths = qa_stats.get('qa_paths', {})

        # Verify expected QA keys are present
        expected_keys = ['sugar_mask_filtered', 'sugar_post_guard_blocked', 'sugar_mask_degraded']
        for key in expected_keys:
            self.assertIn(key, qa_paths, f"QA key '{key}' should be present")
            self.assertIsInstance(qa_paths[key], int, f"QA key '{key}' should be an integer")
            self.assertGreaterEqual(qa_paths[key], 0, f"QA key '{key}' should be non-negative")


if __name__ == '__main__':
    unittest.main(verbosity=2)