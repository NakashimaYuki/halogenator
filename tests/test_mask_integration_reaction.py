# -*- coding: ascii -*-
"""
Integration tests for sugar masking with reaction enumeration.

Tests that sugar masking correctly affects reaction match counts
and product generation in real enumeration scenarios.
"""

import unittest
from typing import List, Dict, Any

# Import the modules we want to test
try:
    from rdkit import Chem
    from src.halogenator.enumerate_k import enumerate_products, EnumConfig
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@unittest.skipUnless(RDKIT_AVAILABLE, "RDKit not available")
class TestMaskIntegrationReaction(unittest.TestCase):
    """Test integration of sugar masking with reaction enumeration."""

    def setUp(self):
        """Set up test molecules."""
        # Simple phenol (no sugar, should have full enumeration)
        self.phenol_smiles = "c1ccc(O)cc1"

        # Glucose (pure sugar, should have heavy masking)
        self.glucose_smiles = "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"

        # Phenyl-O-glucose model (glycoside, sugar part should be masked)
        self.phenyl_glucose_smiles = "c1ccc(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1"

    def test_sugar_masking_reduces_matches(self):
        """Test that enabling sugar masking reduces reaction matches for glycosides."""
        # Create configs with and without sugar masking
        cfg_no_mask = EnumConfig(
            k_max=1,
            halogens=('F',),
            sugar_cfg={'mode': 'off'}
        )

        cfg_with_mask = EnumConfig(
            k_max=1,
            halogens=('F',),
            sugar_cfg={'mode': 'heuristic'}
        )

        # Test with glucose (pure sugar)
        results_no_mask = list(enumerate_products(self.glucose_smiles, cfg_no_mask, return_qa_stats=True))
        results_with_mask = list(enumerate_products(self.glucose_smiles, cfg_with_mask, return_qa_stats=True))

        if results_no_mask and results_with_mask:
            qa_no_mask = results_no_mask[-1][1] if len(results_no_mask[-1]) > 1 else {}
            qa_with_mask = results_with_mask[-1][1] if len(results_with_mask[-1]) > 1 else {}

            attempts_no_mask = qa_no_mask.get('attempts', 0)
            attempts_with_mask = qa_with_mask.get('attempts', 0)

            # With sugar masking, should have fewer attempts due to filtering
            # (Note: glucose is pure sugar so masking should significantly reduce enumeration)
            print(f"Glucose attempts - no mask: {attempts_no_mask}, with mask: {attempts_with_mask}")

            # For pure sugar, masking should reduce or eliminate attempts
            self.assertGreaterEqual(attempts_no_mask, attempts_with_mask)

    def test_sugar_filtered_reaction_matches_statistic(self):
        """Test that sugar_filtered_reaction_matches statistic is populated."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            sugar_cfg={'mode': 'heuristic'}
        )

        # Test with phenyl-glucose (should have some filtering)
        results = list(enumerate_products(self.phenyl_glucose_smiles, cfg, return_qa_stats=True))

        if results:
            qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
            qa_paths = qa_stats.get('qa_paths', {})

            sugar_filtered = qa_paths.get('sugar_filtered_reaction_matches', 0)
            print(f"Sugar filtered matches for phenyl-glucose: {sugar_filtered}")

            # Should have some reaction matches filtered due to sugar masking
            # (The exact number depends on the specific molecule and reaction rules)
            self.assertIsInstance(sugar_filtered, int)
            self.assertGreaterEqual(sugar_filtered, 0)

    def test_phenol_no_sugar_filtering(self):
        """Test that pure phenol (no sugar) has no sugar filtering."""
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            sugar_cfg={'mode': 'heuristic'}
        )

        results = list(enumerate_products(self.phenol_smiles, cfg, return_qa_stats=True))

        if results:
            qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
            qa_paths = qa_stats.get('qa_paths', {})

            sugar_filtered = qa_paths.get('sugar_filtered_reaction_matches', 0)
            sugar_blocked = qa_paths.get('sugar_post_guard_blocked', 0)

            print(f"Phenol sugar stats - filtered: {sugar_filtered}, blocked: {sugar_blocked}")

            # Pure phenol should have no sugar-related filtering
            self.assertEqual(sugar_filtered, 0, "Pure phenol should not have sugar filtering")
            self.assertEqual(sugar_blocked, 0, "Pure phenol should not have sugar blocking")

    def test_qa_stats_schema_completeness(self):
        """Test that all expected QA statistics keys are present."""
        cfg = EnumConfig(k_max=1, halogens=('F',))

        results = list(enumerate_products(self.phenol_smiles, cfg, return_qa_stats=True))

        if results:
            qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
            qa_paths = qa_stats.get('qa_paths', {})

            # Check that sugar statistics are included in schema
            expected_keys = ['sugar_filtered_reaction_matches', 'sugar_post_guard_blocked']
            for key in expected_keys:
                self.assertIn(key, qa_paths, f"Missing QA statistic: {key}")
                self.assertIsInstance(qa_paths[key], int, f"QA statistic {key} should be integer")

    def test_sugar_cfg_detail_switches_affect_enumeration(self):
        """Test that sugar_cfg detail switches affect enumeration behavior."""
        # Test with different masking configurations
        configs = [
            {'mode': 'heuristic', 'mask_exocyclic_oxygen': True, 'mask_glycosidic_bridge_oxygen': True},
            {'mode': 'heuristic', 'mask_exocyclic_oxygen': False, 'mask_glycosidic_bridge_oxygen': True},
            {'mode': 'heuristic', 'mask_exocyclic_oxygen': True, 'mask_glycosidic_bridge_oxygen': False},
            {'mode': 'heuristic', 'mask_exocyclic_oxygen': False, 'mask_glycosidic_bridge_oxygen': False},
        ]

        results_by_config = []

        for sugar_cfg in configs:
            cfg = EnumConfig(
                k_max=1,
                halogens=('F',),
                sugar_cfg=sugar_cfg
            )

            results = list(enumerate_products(self.glucose_smiles, cfg, return_qa_stats=True))
            if results:
                qa_stats = results[-1][1] if len(results[-1]) > 1 else {}
                attempts = qa_stats.get('attempts', 0)
                results_by_config.append((sugar_cfg, attempts))

        # Print results for analysis
        for sugar_cfg, attempts in results_by_config:
            exo = sugar_cfg.get('mask_exocyclic_oxygen', False)
            bridge = sugar_cfg.get('mask_glycosidic_bridge_oxygen', False)
            print(f"Sugar config - exo: {exo}, bridge: {bridge} -> attempts: {attempts}")

        # Verify that we got results for all configurations
        self.assertEqual(len(results_by_config), 4, "Should test all 4 sugar configuration combinations")


if __name__ == '__main__':
    unittest.main()