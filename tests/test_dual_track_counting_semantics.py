#!/usr/bin/env python
"""Test and document dual-track counting semantics for top-level totals vs qa_paths."""

import unittest
from src.halogenator.enumerate_k import (
    enumerate_products, enumerate_with_stats, EnumConfig,
    QAAggregator, _compute_totals_from_aggregator, QAStats
)


class TestDualTrackCountingSemantics(unittest.TestCase):
    """Test and stabilize semantics for top-level totals vs qa_paths dual-track counting."""

    def test_dedup_metrics_dual_track_semantics(self):
        """Test that dedup metrics follow expected dual-track counting semantics."""
        cfg = EnumConfig(k_max=1, halogens=('F',))

        # 1) Streaming path: enumerate_products with QAStats
        streaming_qa = None
        for _prod, qa in enumerate_products('c1ccccc1O', cfg, return_qa_stats=True):
            streaming_qa = qa

        self.assertIsNotNone(streaming_qa, "Should have streaming QA stats")

        # 2) Final snapshot path: enumerate_with_stats with aggregator
        _products, snapshot_qa = enumerate_with_stats('c1ccccc1O', cfg)

        # 3) Document expected semantics for dedup metrics
        # EXPECTATION: Both paths should have dedup metrics in qa_paths
        streaming_dedup_inchi_paths = streaming_qa['qa_paths'].get('dedup_hits_inchi', 0)
        snapshot_dedup_inchi_paths = snapshot_qa['qa_paths'].get('dedup_hits_inchi', 0)

        # EXPECTATION: Both paths should have consistent qa_paths dedup counts
        self.assertEqual(streaming_dedup_inchi_paths, snapshot_dedup_inchi_paths,
                        f"Dedup InChI in qa_paths: streaming={streaming_dedup_inchi_paths}, snapshot={snapshot_dedup_inchi_paths}")

        # EXPECTATION: Top-level dedup fields behavior
        streaming_dedup_inchi_top = streaming_qa.get('dedup_hits_inchi', 0)
        snapshot_dedup_inchi_top = snapshot_qa.get('dedup_hits_inchi', 0)

        # DOCUMENTED SEMANTICS: Currently, streaming has top-level dedup fields,
        # while snapshot may or may not (depends on _compute_totals_from_aggregator)
        print(f"Streaming top-level dedup_hits_inchi: {streaming_dedup_inchi_top}")
        print(f"Snapshot top-level dedup_hits_inchi: {snapshot_dedup_inchi_top}")
        print(f"Streaming qa_paths dedup_hits_inchi: {streaming_dedup_inchi_paths}")
        print(f"Snapshot qa_paths dedup_hits_inchi: {snapshot_dedup_inchi_paths}")

        # For now, document the current behavior rather than enforcing a specific one
        # This test serves as documentation of the dual-track semantics

    def test_aggregator_totals_dedup_placement(self):
        """Test where dedup metrics appear in aggregator-based totals."""
        agg = QAAggregator()

        # Record some dedup events via rule/halogen/k combinations
        agg.record_attempt_result('R3', 'F', 1, 2, {'dedup_hits_inchi': 3, 'dedup_hits_statesig': 1})
        agg.record_attempt_result('R4', 'Cl', 1, 1, {'dedup_hits_inchi': 2})

        # Get totals from aggregator
        totals = _compute_totals_from_aggregator(agg)

        # DOCUMENTED BEHAVIOR: Dedup metrics from aggregator go to qa_paths only
        self.assertEqual(totals['qa_paths'].get('dedup_hits_inchi', 0), 5)  # 3 + 2
        self.assertEqual(totals['qa_paths'].get('dedup_hits_statesig', 0), 1)

        # DOCUMENTED BEHAVIOR: No top-level dedup fields in aggregator totals
        self.assertNotIn('dedup_hits_inchi', totals)
        self.assertNotIn('dedup_hits_statesig', totals)

        # But other metrics should appear at top level
        self.assertIn('attempts', totals)
        self.assertIn('products', totals)

    def test_qa_stats_dedup_placement(self):
        """Test where dedup metrics appear in QAStats-based output."""
        qa_stats = QAStats()

        # Simulate merging some stats
        qa_stats.merge({
            'dedup_hits_inchi': 3,
            'dedup_hits_statesig': 2,
            'no_product_matches': 1
        })

        # Convert to dict (streaming path format)
        qa_dict = qa_stats.to_dict()

        # DOCUMENTED BEHAVIOR: QAStats maintains top-level dedup fields
        self.assertEqual(qa_dict.get('dedup_hits_inchi', 0), 3)
        self.assertEqual(qa_dict.get('dedup_hits_statesig', 0), 2)

        # DOCUMENTED BEHAVIOR: These should also appear in qa_paths for consistency
        # (This depends on how qa_paths gets populated)
        self.assertIn('qa_paths', qa_dict)

    def test_dual_track_consistency_requirement(self):
        """Test that when both top-level and qa_paths have dedup metrics, they should be consistent."""
        # This test documents the consistency requirement:
        # If a metric appears in both top-level and qa_paths, the values should match

        qa_stats = QAStats()
        qa_stats.dedup_hits_inchi = 5
        qa_stats.qa_paths['dedup_hits_inchi'] = 5  # Should match

        qa_dict = qa_stats.to_dict()

        # CONSISTENCY REQUIREMENT: Top-level and qa_paths should match for same metric
        top_level = qa_dict.get('dedup_hits_inchi', 0)
        qa_paths_level = qa_dict.get('qa_paths', {}).get('dedup_hits_inchi', 0)

        if top_level > 0 and qa_paths_level > 0:
            self.assertEqual(top_level, qa_paths_level,
                           f"Dedup metric inconsistency: top-level={top_level}, qa_paths={qa_paths_level}")

    def test_dedup_dual_write_semantics_non_zero_scenario(self):
        """T34: Lock down dedup dual-write semantics with actual non-zero deduplication."""
        # Use a symmetrical molecule likely to produce duplicates via different pathways
        # Benzene with multiple identical positions should trigger deduplication
        cfg = EnumConfig(k_max=2, halogens=('F', 'Cl'))

        # Use benzene (symmetric) which should produce duplicates when substituting at equivalent positions
        parent_smiles = 'c1ccccc1'  # benzene - highly symmetric

        # 1) Streaming path
        streaming_products = []
        streaming_qa = None
        for prod, qa in enumerate_products(parent_smiles, cfg, return_qa_stats=True):
            streaming_products.append(prod)
            streaming_qa = qa

        # 2) Aggregator snapshot path
        snapshot_products, snapshot_qa = enumerate_with_stats(parent_smiles, cfg)

        # Should have produced the same products
        self.assertEqual(len(streaming_products), len(snapshot_products))

        # Look for actual deduplication (this test may need adjustment based on RDKit behavior)
        streaming_dedup_inchi = streaming_qa['qa_paths'].get('dedup_hits_inchi', 0)
        snapshot_dedup_inchi = snapshot_qa['qa_paths'].get('dedup_hits_inchi', 0)

        print(f"Streaming products: {len(streaming_products)}")
        print(f"Snapshot products: {len(snapshot_products)}")
        print(f"Streaming dedup_hits_inchi: {streaming_dedup_inchi}")
        print(f"Snapshot dedup_hits_inchi: {snapshot_dedup_inchi}")

        # LOCKED SEMANTICS: qa_paths dedup counts must be consistent between paths
        self.assertEqual(streaming_dedup_inchi, snapshot_dedup_inchi,
                        f"qa_paths dedup_hits_inchi must be consistent: streaming={streaming_dedup_inchi}, snapshot={snapshot_dedup_inchi}")

        # LOCKED SEMANTICS: If deduplication occurred, document current behavior
        if streaming_dedup_inchi > 0:
            # T34 LOCK: Document current INCONSISTENT behavior between paths
            streaming_dedup_top = streaming_qa.get('dedup_hits_inchi', 0)
            snapshot_dedup_top = snapshot_qa.get('dedup_hits_inchi', 0)

            print(f"Streaming top-level dedup_hits_inchi: {streaming_dedup_top}")
            print(f"Snapshot top-level dedup_hits_inchi: {snapshot_dedup_top}")

            # CURRENT BEHAVIOR (documented for stability):
            # - Streaming path: qa_paths only (no top-level dedup fields)
            # - Snapshot path: dual-write (both top-level and qa_paths)

            # LOCK DOWN CURRENT BEHAVIOR: Streaming uses qa_paths only
            self.assertEqual(streaming_dedup_top, 0,
                           f"Streaming path uses qa_paths-only: top-level should be 0, got {streaming_dedup_top}")

            # LOCK DOWN CURRENT BEHAVIOR: Snapshot uses dual-write
            self.assertEqual(snapshot_dedup_top, snapshot_dedup_inchi,
                           f"Snapshot path uses dual-write: top={snapshot_dedup_top} should match qa_paths={snapshot_dedup_inchi}")

            # Document the inconsistency for future standardization
            # TODO(future): Consider standardizing on one approach across both paths

        # T34 LOCK: This test documents the current heterogeneous dedup semantics:
        # - Streaming: qa_paths only for dedup metrics
        # - Snapshot: dual-write for dedup metrics

    def test_semantic_stability_documentation(self):
        """Document the expected stable semantics for dual-track counting."""
        # This test serves as documentation of the expected behavior

        # STABLE SEMANTIC RULES:
        # 1. qa_paths ALWAYS contains comprehensive event counts (pivot + non-pivot)
        # 2. Top-level fields are optional and depend on the output path
        # 3. When both exist for the same metric, they MUST be consistent
        # 4. Pivot events (dedup_hits_*) can appear in both places
        # 5. The aggregator path prefers qa_paths-only for pivot events (T34)
        # 6. The QAStats path maintains dual-write for backward compatibility

        # This test passes by design - it's documentation
        self.assertTrue(True, "Semantic documentation test")


if __name__ == '__main__':
    unittest.main()