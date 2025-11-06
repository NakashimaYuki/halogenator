#!/usr/bin/env python
"""Test aggregator paths container functionality for N1 fix."""

import unittest

from src.halogenator.enumerate_k import QAAggregator, _compute_totals_from_aggregator


class TestAggregatorPaths(unittest.TestCase):
    """Test that aggregator paths container works without dimension pollution."""

    def test_aggregator_record_paths_no_rule_noise(self):
        """Test that record_paths doesn't pollute by_rule_halogen_k dimensions and rejects pivot events."""
        agg = QAAggregator()

        # Record non-pivot paths events (should work)
        agg.record_paths({'sugar_proximity_filtered': 1, 'sugar_mask_degraded': 2})

        # Check that non-pivot paths are stored in dedicated container
        self.assertEqual(agg.paths.get('sugar_proximity_filtered', 0), 1)
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 2)

        # Pivot events should be rejected by record_paths (defensive)
        with self.assertLogs(level='WARNING') as log_context:
            agg.record_paths({'sugar_mask_filtered': 3})  # This is a pivot event, should be rejected

        # Verify pivot event was rejected
        self.assertEqual(agg.paths.get('sugar_mask_filtered', 0), 0)
        self.assertTrue(any('ignored pivot event' in log.getMessage() for log in log_context.records))

        # Get pivots
        pivots = agg.to_pivots_dict()

        # by_rule_halogen_k should not have any pollution from paths
        grid = pivots.get('by_rule_halogen_k', {})

        # Check no 'sample_sugar_0' type keys exist
        for key in grid.keys():
            self.assertFalse(
                'sample' in str(key) and 'sugar' in str(key),
                f"Found pollution in by_rule_halogen_k: {key}"
            )

        # Check that non-pivot paths appear in final totals via _compute_totals_from_aggregator
        totals = _compute_totals_from_aggregator(agg)
        qa_paths = totals.get('qa_paths', {})

        self.assertEqual(qa_paths.get('sugar_proximity_filtered', 0), 1)
        self.assertEqual(qa_paths.get('sugar_mask_degraded', 0), 2)
        # Pivot event should NOT appear from paths (should only come from pivot aggregation)
        self.assertEqual(qa_paths.get('sugar_mask_filtered', 0), 0)

    def test_aggregator_merge_paths(self):
        """Test that aggregator merge correctly handles paths containers."""
        agg1 = QAAggregator()
        agg1.record_paths({'sugar_mask_degraded': 3})

        agg2 = QAAggregator()
        agg2.record_paths({'sugar_mask_degraded': 2, 'sugar_proximity_filtered': 1})

        # Merge agg2 into agg1
        agg1.merge(agg2)

        # Check merged paths
        self.assertEqual(agg1.paths.get('sugar_mask_degraded', 0), 5)  # 3 + 2
        self.assertEqual(agg1.paths.get('sugar_proximity_filtered', 0), 1)

        # Verify totals computation works with merged paths
        totals = _compute_totals_from_aggregator(agg1)
        qa_paths = totals.get('qa_paths', {})

        self.assertEqual(qa_paths.get('sugar_mask_degraded', 0), 5)
        self.assertEqual(qa_paths.get('sugar_proximity_filtered', 0), 1)

    def test_aggregator_reset_paths(self):
        """Test that reset clears paths container."""
        agg = QAAggregator()
        agg.record_paths({'sugar_mask_degraded': 5})

        # Verify paths recorded
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 5)

        # Reset and verify cleared
        agg.reset()
        self.assertEqual(len(agg.paths), 0)

    def test_debug_mode_strict_pivot_validation(self):
        """Test that debug mode throws exceptions instead of warnings for pivot events."""
        from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

        # Test with debug mode enabled
        agg_debug = QAAggregator(debug_consistency=True)

        # Should throw ValueError for pivot events in record_paths
        with self.assertRaises(ValueError) as cm:
            agg_debug.record_paths({'sugar_mask_filtered': 1})
        self.assertIn('STRICT MODE', str(cm.exception))
        self.assertIn('sugar_mask_filtered', str(cm.exception))

        # Should throw ValueError for pivot events in record_path
        with self.assertRaises(ValueError) as cm:
            agg_debug.record_path('isotope_unavailable', 2)
        self.assertIn('STRICT MODE', str(cm.exception))
        self.assertIn('isotope_unavailable', str(cm.exception))

        # Non-pivot events should work fine
        agg_debug.record_path('sugar_mask_degraded', 3)
        self.assertEqual(agg_debug.paths.get('sugar_mask_degraded', 0), 3)

        # Test with debug mode disabled (default)
        agg_normal = QAAggregator(debug_consistency=False)

        # Should only warn, not throw, for pivot events
        with self.assertLogs(level='WARNING') as log_context:
            agg_normal.record_paths({'sugar_mask_filtered': 1})
        self.assertTrue(any('ignored pivot event' in log.getMessage() for log in log_context.records))

        # Should not be recorded in paths
        self.assertEqual(agg_normal.paths.get('sugar_mask_filtered', 0), 0)

    def test_pivot_event_sums_validation(self):
        """Test pivot_event_sums method provides correct aggregation for validation."""
        from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

        agg = QAAggregator()

        # Record some pivot events via rule/halogen/k combinations
        agg.record_attempt_result('R3', 'F', 1, 2, {'isotope_unavailable': 1, 'atommap_used': 2})
        agg.record_attempt_result('R3', 'Cl', 1, 1, {'isotope_miss': 1, 'atommap_used': 1})
        agg.record_attempt_result('R4', 'F', 1, 1, {'isotope_unavailable': 1})

        # Get pivot sums
        pivot_sums = agg.pivot_event_sums()

        # Verify expected totals
        self.assertEqual(pivot_sums.get('isotope_unavailable', 0), 2)  # 1 + 1
        self.assertEqual(pivot_sums.get('atommap_used', 0), 3)  # 2 + 1
        self.assertEqual(pivot_sums.get('isotope_miss', 0), 1)  # 1

        # Non-recorded events should not appear
        self.assertEqual(pivot_sums.get('rdkit_error', 0), 0)

        # Compare with _compute_totals_from_aggregator qa_paths for pivot events
        totals = _compute_totals_from_aggregator(agg)
        qa_paths = totals.get('qa_paths', {})

        # Pivot events in qa_paths should match pivot_event_sums exactly
        for pivot_key in PIVOT_EVENT_KEYS:
            pivot_sum = pivot_sums.get(pivot_key, 0)
            qa_paths_count = qa_paths.get(pivot_key, 0)
            self.assertEqual(pivot_sum, qa_paths_count,
                           f"Pivot event '{pivot_key}' mismatch: pivot_sums={pivot_sum}, qa_paths={qa_paths_count}")

    def test_record_path_single_event_api_comprehensive(self):
        """Test comprehensive coverage of record_path() single-event API."""
        from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

        agg = QAAggregator()

        # Test 1: Record non-pivot events successfully
        agg.record_path('sugar_mask_degraded', 3)
        agg.record_path('sugar_proximity_filtered')  # Default amount=1
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 3)
        self.assertEqual(agg.paths.get('sugar_proximity_filtered', 0), 1)

        # Test 2: Multiple calls to same event should accumulate
        agg.record_path('sugar_mask_degraded', 2)
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 5)  # 3 + 2

        # Test 3: Alias handling - input legacy names, get canonical names
        agg.record_path('post_guard_blocked', 1)  # Should become 'sugar_post_guard_blocked'
        # But since this is a pivot event, it should be rejected with warning
        with self.assertLogs(level='WARNING') as log_context:
            agg.record_path('post_guard_blocked', 1)
        self.assertTrue(any('ignored pivot event' in log.getMessage() for log in log_context.records))
        self.assertEqual(agg.paths.get('sugar_post_guard_blocked', 0), 0)

        # Test 4: Rejection of pivot events in non-debug mode
        with self.assertLogs(level='WARNING') as log_context:
            agg.record_path('isotope_unavailable', 5)
        self.assertTrue(any('ignored pivot event' in log.getMessage() for log in log_context.records))
        self.assertEqual(agg.paths.get('isotope_unavailable', 0), 0)

        # Test 5: Zero and negative amounts (edge cases - now restricted to positive values only)
        agg.record_path('sugar_mask_degraded', 0)  # Should not change
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 5)

        agg.record_path('sugar_mask_degraded', -1)  # Should be ignored (T32: positive values only)
        self.assertEqual(agg.paths.get('sugar_mask_degraded', 0), 5)  # Should remain unchanged

        # Test 6: Verify paths container only contains non-pivot events
        contaminated_keys = set(agg.paths.keys()) & PIVOT_EVENT_KEYS
        self.assertEqual(len(contaminated_keys), 0,
                        f"Paths container contaminated with pivot events: {contaminated_keys}")

        # Test 7: Canonical alias processing for non-pivot events
        # (Need to create a test alias that maps to non-pivot event)
        # For now, test direct canonical names work correctly
        expected_non_pivot_events = {'sugar_mask_degraded', 'sugar_proximity_filtered'}
        actual_events = set(agg.paths.keys())
        self.assertTrue(expected_non_pivot_events.issubset(actual_events),
                       f"Expected non-pivot events {expected_non_pivot_events} should be present in {actual_events}")


if __name__ == '__main__':
    unittest.main()