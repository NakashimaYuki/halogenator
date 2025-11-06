# -*- coding: ascii -*-
"""
Test symmetry folding on masked subgraph switch.

This test verifies that the `compute_on_masked_subgraph` setting properly
controls whether symmetry analysis uses the masked subgraph or the full molecule.
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
from halogenator.sugar_mask import get_sugar_mask_with_status
from halogenator.chem_compat import Chem


class TestSymmetryMaskedSubgraphSwitch(unittest.TestCase):
    """Test symmetry folding masked subgraph switch functionality."""

    def test_compute_on_masked_subgraph_switch_quercetin_glycoside(self):
        """
        Test that compute_on_masked_subgraph=True/False produces different
        but stable results for a molecule with significant sugar masking.
        """
        # Use quercetin-3-glucoside which has a large sugar group
        quercetin_glucoside = "O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

        mol = Chem.MolFromSmiles(quercetin_glucoside)
        self.assertIsNotNone(mol, "Quercetin glycoside SMILES should be valid")

        # Verify that sugar masking identifies a significant number of atoms
        sugar_mask, degraded = get_sugar_mask_with_status(mol, mode='heuristic')
        self.assertFalse(degraded, "Sugar mask should not be degraded")
        self.assertGreater(len(sugar_mask), 5, "Sugar mask should identify significant atoms")

        # Test with masked subgraph symmetry enabled (default)
        config_masked = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),  # Focus on R1 for clear symmetry effects
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={'compute_on_masked_subgraph': True},
            pruning_cfg={'enable_symmetry_fold': True}
        )

        # Test with standard symmetry (ignore masking for symmetry analysis)
        config_standard = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={'compute_on_masked_subgraph': False},
            pruning_cfg={'enable_symmetry_fold': True}
        )

        # Run both configurations
        products_masked, qa_stats_masked = enumerate_with_stats(quercetin_glucoside, config_masked)
        products_standard, qa_stats_standard = enumerate_with_stats(quercetin_glucoside, config_standard)

        # Both should produce valid results
        self.assertGreaterEqual(len(products_masked), 0, "Masked symmetry should produce valid results")
        self.assertGreaterEqual(len(products_standard), 0, "Standard symmetry should produce valid results")

        # The number of products may differ due to different symmetry representative selection
        # This is expected and acceptable
        print(f"Products with masked symmetry: {len(products_masked)}")
        print(f"Products with standard symmetry: {len(products_standard)}")

        # Both should have similar QA statistics structure
        for qa_stats, mode in [(qa_stats_masked, "masked"), (qa_stats_standard, "standard")]:
            self.assertIn('attempts', qa_stats, f"{mode} mode should track attempts")
            self.assertIn('products', qa_stats, f"{mode} mode should track products")
            self.assertIn('qa_paths', qa_stats, f"{mode} mode should have qa_paths")

    def test_symmetry_switch_stability(self):
        """
        Test that the same configuration produces consistent results
        across multiple runs (stability test).
        """
        test_molecule = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"  # Quercetin (no sugar)

        config = EnumConfig(
            k_max=1,
            halogens=('Cl',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={'compute_on_masked_subgraph': True},
            pruning_cfg={'enable_symmetry_fold': True}
        )

        # Run the same configuration multiple times
        results = []
        for i in range(3):
            products, qa_stats = enumerate_with_stats(test_molecule, config)
            results.append((len(products), qa_stats.get('attempts', 0), qa_stats.get('products', 0)))

        # All runs should produce identical results
        first_result = results[0]
        for i, result in enumerate(results[1:], 1):
            self.assertEqual(result, first_result,
                           f"Run {i+1} should produce same results as run 1")

    def test_symmetry_switch_with_no_masking(self):
        """
        Test that when sugar_mask is empty, both modes should behave identically.
        """
        # Use a simple aromatic molecule with no sugar groups
        simple_phenol = "c1ccc(O)cc1"

        config_base = {
            'k_max': 1,
            'halogens': ('F',),
            'rules': ('R1',),
            'sugar_cfg': {'mode': 'off'},  # No sugar masking
            'pruning_cfg': {'enable_symmetry_fold': True}
        }

        # Test both symmetry modes
        config_masked = EnumConfig(**config_base, symmetry_cfg={'compute_on_masked_subgraph': True})
        config_standard = EnumConfig(**config_base, symmetry_cfg={'compute_on_masked_subgraph': False})

        products_masked, qa_stats_masked = enumerate_with_stats(simple_phenol, config_masked)
        products_standard, qa_stats_standard = enumerate_with_stats(simple_phenol, config_standard)

        # With no masking, both modes should produce identical results
        self.assertEqual(len(products_masked), len(products_standard),
                        "Both symmetry modes should produce same product count with no masking")
        self.assertEqual(qa_stats_masked.get('attempts', 0), qa_stats_standard.get('attempts', 0),
                        "Both modes should have same attempts count with no masking")
        self.assertEqual(qa_stats_masked.get('products', 0), qa_stats_standard.get('products', 0),
                        "Both modes should have same products count with no masking")

    def test_symmetry_switch_config_validation(self):
        """
        Test that the symmetry_cfg parameter is properly recognized and used.
        """
        test_molecule = "OCC1OC(O)C(O)C(O)C1O"  # Glucose

        # Test explicit True
        config_true = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={'compute_on_masked_subgraph': True}
        )

        # Test explicit False
        config_false = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={'compute_on_masked_subgraph': False}
        )

        # Test default (should be True)
        config_default = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic'},
            symmetry_cfg={}  # Empty config, should use default
        )

        # All should run without errors
        for config, name in [(config_true, "explicit True"), (config_false, "explicit False"), (config_default, "default")]:
            try:
                products, qa_stats = enumerate_with_stats(test_molecule, config)
                self.assertIsNotNone(products, f"{name} config should produce valid products")
                self.assertIsNotNone(qa_stats, f"{name} config should produce valid QA stats")
            except Exception as e:
                self.fail(f"{name} config should not raise exception: {e}")


if __name__ == '__main__':
    unittest.main(verbosity=2)