# -*- coding: ascii -*-
import unittest, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig, enumerate_with_stats


def sugar_total(qa_paths: dict) -> int:
    """
    Calculate total sugar events using same logic as acceptance script.
    Avoids double counting legacy keys by using fallback logic.
    """
    direct = int(qa_paths.get('sugar_mask_filtered', 0))
    prox = int(qa_paths.get('sugar_proximity_filtered', 0))
    post = int(qa_paths.get('post_guard_blocked',
                   qa_paths.get('sugar_post_guard_blocked', 0)))
    return direct + prox + post


class TestSugarEvents(unittest.TestCase):
    def test_o_glycoside_proximity_guard_triggers(self):
        # phenyl-O-glucoside (simple O-glycoside with aromatic CH sites near sugar bridge O)
        smiles = "Oc1ccccc1OC2OC(CO)C(O)C(O)C2O"
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1','R2'),
            rules_cfg={'R2': {'sp2_CH_in_C_ring': True, 'sp3_CH2_flavanone': False}},
            sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 2}
        )
        products, qa_stats = enumerate_with_stats(smiles, cfg)
        qa_paths = (qa_stats or {}).get('qa_paths', {})
        total = sugar_total(qa_paths)
        self.assertGreater(total, 0, "O-glycoside should produce non-zero sugar_events via proximity guard")

    def test_aglycone_has_zero_sugar_events_even_with_radius(self):
        # toluene (aglycone - no sugar)
        smiles = "Cc1ccccc1"
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 3}
        )
        _, qa_stats = enumerate_with_stats(smiles, cfg)
        qa_paths = (qa_stats or {}).get('qa_paths', {})
        total = sugar_total(qa_paths)
        self.assertEqual(total, 0, "aglycone should have zero sugar_events even with proximity guard enabled")

    def test_glycoside_radius_zero_produces_zero_events(self):
        # Same O-glycoside but with radius=0 should have zero events
        smiles = "Oc1ccccc1OC2OC(CO)C(O)C(O)C2O"
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1','R2'),
            rules_cfg={'R2': {'sp2_CH_in_C_ring': True, 'sp3_CH2_flavanone': False}},
            sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 0}
        )
        _, qa_stats = enumerate_with_stats(smiles, cfg)
        qa_paths = (qa_stats or {}).get('qa_paths', {})
        total = sugar_total(qa_paths)
        self.assertEqual(total, 0, "glycoside with radius=0 should have zero sugar_events")

    def test_legacy_key_backward_compatibility(self):
        # Test that legacy sugar_post_guard_blocked key is properly handled
        from halogenator.enumerate_k import EnumConfig, enumerate_with_stats

        # Create a simple test case
        smiles = "Oc1ccccc1OC2OC(CO)C(O)C(O)C2O"
        cfg = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 2}
        )
        _, qa_stats = enumerate_with_stats(smiles, cfg)
        qa_paths = (qa_stats or {}).get('qa_paths', {})

        # Check that both new and legacy keys exist when post_guard_blocked > 0
        post_guard_blocked = qa_paths.get('post_guard_blocked', 0)
        legacy_post_guard = qa_paths.get('sugar_post_guard_blocked', 0)

        if post_guard_blocked > 0:
            self.assertGreaterEqual(legacy_post_guard, post_guard_blocked,
                "legacy sugar_post_guard_blocked should include post_guard_blocked value")

if __name__ == '__main__':
    unittest.main(verbosity=2)