# -*- coding: ascii -*-
"""Tests for QA statistics counting semantics."""

import unittest
from src.halogenator.enumerate_k import enumerate_products, EnumConfig


class TestQACountSemantics(unittest.TestCase):
    """Test QA statistics follow per-attempt semantics correctly."""
    
    def test_qa_count_semantics_when_no_matches(self):
        """
        Test isotope_unavailable is counted per-attempt, not per-product,
        when matches_with_sites is empty.
        """
        # Use a simple substrate that triggers R3 with multiple products but no pre-matches
        cfg = EnumConfig(k_max=1, halogens=('F', 'Cl'), 
                        constraints={'per_ring_quota': 10, 'min_graph_distance': 1})
        
        # Collect all products and final QA stats
        products = list(enumerate_products('c1ccccc1O', cfg, return_qa_stats=True))
        
        if products:
            _, final_qa_stats = products[-1]
            
            # There should be products, but isotope_unavailable should be counted per rule/halogen attempt
            # Not per product count
            isotope_unavailable = final_qa_stats['qa_paths']['isotope_unavailable']
            
            # For R3, R4, R5 rules with F and Cl, should be up to 6 attempts total
            # (3 rules x 2 halogens), regardless of product count
            max_expected = 3 * len(cfg.halogens)  # 3 rules x 2 halogens = 6
            self.assertTrue(isotope_unavailable <= max_expected, 
                          f"isotope_unavailable={isotope_unavailable} should be <={max_expected} (per rule/halogen attempt, not per product)")
            
            # Should be much less than product count (which would indicate per-product counting)
            self.assertTrue(isotope_unavailable < len(products),
                          f"isotope_unavailable={isotope_unavailable} should be much less than products={len(products)}")
            self.assertTrue(len(products) > 0, "Should generate products")
    
    def test_qa_count_semantics_isotope_miss(self):
        """
        Test isotope_miss is counted per failed match, not per product.
        This is harder to trigger naturally, so we test the general structure.
        """
        cfg = EnumConfig(k_max=1, halogens=('F',))
        
        products = list(enumerate_products('c1ccccc1O', cfg, return_qa_stats=True))
        
        if products:
            _, final_qa_stats = products[-1]
            
            # Verify structure contains the new field
            self.assertIn('isotope_miss', final_qa_stats['qa_paths'])
            self.assertIsInstance(final_qa_stats['qa_paths']['isotope_miss'], int)
            
            # In normal operation, isotope_miss should be rare (0 or very low)
            isotope_miss = final_qa_stats['qa_paths']['isotope_miss']
            self.assertTrue(isotope_miss >= 0, "isotope_miss should be non-negative")
    
    def test_qa_stats_structure_matches_semantics(self):
        """Test that QA stats structure matches documented semantics."""
        cfg = EnumConfig(k_max=1, halogens=('F',))
        
        products = list(enumerate_products('c1ccccc1O', cfg, return_qa_stats=True))
        
        if products:
            _, qa_stats = products[-1]
            
            # Verify all expected fields are present
            expected_fields = {
                'no_product_matches', 'dedup_hits_statesig', 'dedup_hits_inchi', 
                'template_unsupported', 'qa_paths'
            }
            self.assertEqual(set(qa_stats.keys()), expected_fields)
            
            # Verify qa_paths substructure contains at least the core fields
            # (may contain additional fields from STANDARD_QA_PATHS_KEYS)
            core_qa_paths = {
                'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'
            }
            actual_qa_paths = set(qa_stats['qa_paths'].keys())
            self.assertTrue(core_qa_paths.issubset(actual_qa_paths),
                           f"qa_paths should contain at least {core_qa_paths}, got {actual_qa_paths}")

            # PIVOT EVENTS VALIDATION:
            # qa_paths should contain core pivot events (legitimate, from pivot aggregation)
            # Detailed contamination testing is covered in test_pivot_events_source_validation
            from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

            # Verify core pivot events exist in qa_paths when expected
            qa_paths_pivot_events = set(qa_stats['qa_paths'].keys()) & PIVOT_EVENT_KEYS
            core_pivot_events = {'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'}
            if core_pivot_events & actual_qa_paths:
                # If any core pivot events exist, they should be properly sourced from pivot aggregation
                # (This is comprehensively tested in test_pivot_events_source_validation)
                pass

            # All counters should be non-negative integers
            for field in expected_fields:
                if field == 'qa_paths':
                    for qa_field, qa_value in qa_stats[field].items():
                        self.assertIsInstance(qa_value, int, f"qa_paths.{qa_field} should be int")
                        self.assertTrue(qa_value >= 0, f"qa_paths.{qa_field} should be non-negative")
                else:
                    self.assertIsInstance(qa_stats[field], int, f"{field} should be int")
                    self.assertTrue(qa_stats[field] >= 0, f"{field} should be non-negative")

    def test_pivot_events_source_validation(self):
        """Test that pivot events appear in qa_paths from pivot aggregation, not from paths contamination."""
        from src.halogenator.enumerate_k import enumerate_with_stats, QAAggregator
        from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

        cfg = EnumConfig(k_max=1, halogens=('F', 'Cl'))

        # Use enumerate_with_stats to get both products and QA stats
        products, qa_stats = enumerate_with_stats('c1ccccc1O', cfg)

        # Create a manual aggregator to test the _compute_totals_from_aggregator logic
        aggregator = QAAggregator()

        # Simulate some pivot events from rule/halogen/k combinations
        aggregator.record_attempt_result('R3', 'F', 1, 2, {'isotope_unavailable': 1, 'atommap_used': 1})
        aggregator.record_attempt_result('R3', 'Cl', 1, 1, {'isotope_miss': 1})

        # Test defensive rejection: try to contaminate paths with pivot events
        with self.assertLogs(level='WARNING') as log_context:
            aggregator.record_paths({'isotope_unavailable': 5})  # Should be rejected

        # Verify contamination was rejected
        self.assertEqual(aggregator.paths.get('isotope_unavailable', 0), 0)
        self.assertTrue(any('ignored pivot event' in log.getMessage() for log in log_context.records))

        # Add some legitimate non-pivot path events
        aggregator.record_paths({'sugar_mask_degraded': 2, 'sugar_proximity_filtered': 1})

        # Verify paths container only has non-pivot events
        contaminated_keys = set(aggregator.paths.keys()) & PIVOT_EVENT_KEYS
        self.assertEqual(len(contaminated_keys), 0,
                        f"aggregator.paths contaminated with pivot events: {contaminated_keys}")

        # Test the totals computation
        from src.halogenator.enumerate_k import _compute_totals_from_aggregator
        totals = _compute_totals_from_aggregator(aggregator)

        # Verify pivot events appear in qa_paths (from pivot aggregation)
        self.assertEqual(totals['qa_paths'].get('isotope_unavailable', 0), 1)
        self.assertEqual(totals['qa_paths'].get('atommap_used', 0), 1)
        self.assertEqual(totals['qa_paths'].get('isotope_miss', 0), 1)

        # Verify non-pivot events also appear in qa_paths (from paths)
        self.assertEqual(totals['qa_paths'].get('sugar_mask_degraded', 0), 2)
        self.assertEqual(totals['qa_paths'].get('sugar_proximity_filtered', 0), 1)

        # Verify that the failed contamination attempt doesn't appear
        # (the isotope_unavailable=5 should have been rejected, only the legitimate=1 should appear)
        self.assertNotEqual(totals['qa_paths'].get('isotope_unavailable', 0), 5)

    def test_streaming_vs_snapshot_pivot_consistency(self):
        """Test that streaming and final snapshot paths produce identical pivot event counts."""
        from src.halogenator.enumerate_k import enumerate_products, enumerate_with_stats
        from src.halogenator.qa_utils import PIVOT_EVENT_KEYS

        cfg = EnumConfig(k_max=1, halogens=('F', 'Cl'))

        # 1) Streaming approach: collect final qa_stats from enumerate_products
        last_qa = None
        for _prod, qa in enumerate_products('c1ccccc1O', cfg, return_qa_stats=True):
            last_qa = qa
        self.assertIsNotNone(last_qa, "Should have received QA stats from streaming")

        # 2) Final snapshot approach: use enumerate_with_stats
        _products, snap_qa = enumerate_with_stats('c1ccccc1O', cfg)

        # 3) Compare pivot event counts between both approaches
        for pivot_key in PIVOT_EVENT_KEYS:
            streaming_count = last_qa['qa_paths'].get(pivot_key, 0)
            snapshot_count = snap_qa['qa_paths'].get(pivot_key, 0)
            self.assertEqual(streaming_count, snapshot_count,
                           f"Pivot event '{pivot_key}' mismatch: streaming={streaming_count}, snapshot={snapshot_count}")

        # 4) Verify that both contain the expected core pivot events
        core_pivot_events = {'isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used'}
        streaming_pivot_keys = set(last_qa['qa_paths'].keys()) & PIVOT_EVENT_KEYS
        snapshot_pivot_keys = set(snap_qa['qa_paths'].keys()) & PIVOT_EVENT_KEYS

        # Both should have the same set of pivot events present
        self.assertEqual(streaming_pivot_keys, snapshot_pivot_keys,
                        "Streaming and snapshot should have identical pivot event key sets")

        # Both should contain the core pivot events (when they occur)
        if core_pivot_events & streaming_pivot_keys:
            self.assertTrue(core_pivot_events.issubset(streaming_pivot_keys | snapshot_pivot_keys),
                           "Core pivot events should be present in both approaches when they occur")


if __name__ == '__main__':
    unittest.main()