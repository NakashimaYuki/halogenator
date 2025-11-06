# -*- coding: ascii -*-
"""Contract tests for QAAggregator attempt semantics."""

import unittest

from halogenator.enumerate_k import QAAggregator


class TestAggregatorAttemptsSemantics(unittest.TestCase):
    """Ensure attempts/products/no_match/template_unsupported semantics hold."""

    def test_successful_attempt_counts_once(self):
        agg = QAAggregator()
        # One attempt, multiple products -> products += 1 (not N)
        agg.record_attempt_result('R1', 'F', 1, produced_count=3, qa_events_dict={})
        piv = agg.to_pivots_dict()
        key = 'R1_F_1'
        self.assertEqual(piv['by_rule_halogen_k'][key]['attempts'], 1)
        self.assertEqual(piv['by_rule_halogen_k'][key]['products'], 1)
        self.assertEqual(piv['by_rule_halogen_k'][key].get('no_product_matches', 0), 0)

    def test_failed_attempt_counts_no_match(self):
        agg = QAAggregator()
        agg.record_attempt_result('R1', 'Cl', 1, produced_count=0, qa_events_dict={})
        piv = agg.to_pivots_dict()
        key = 'R1_Cl_1'
        self.assertEqual(piv['by_rule_halogen_k'][key]['attempts'], 1)
        self.assertEqual(piv['by_rule_halogen_k'][key].get('products', 0), 0)
        self.assertEqual(piv['by_rule_halogen_k'][key]['no_product_matches'], 1)

    def test_invariant_holds_in_mixed_batch(self):
        agg = QAAggregator()
        # Success
        agg.record_attempt_result('R1', 'F', 1, produced_count=2, qa_events_dict={})
        # Failure
        agg.record_attempt_result('R1', 'F', 1, produced_count=0, qa_events_dict={})
        # Another success
        agg.record_attempt_result('R2', 'Br', 2, produced_count=1, qa_events_dict={})
        piv = agg.to_pivots_dict()['by_rule_halogen_k']

        # Sum across all keys
        attempts = sum(ev.get('attempts', 0) for ev in piv.values())
        products = sum(ev.get('products', 0) for ev in piv.values())
        no_matches = sum(ev.get('no_product_matches', 0) for ev in piv.values())

        # The contract: attempts == products + no_product_matches (+ template_unsupported if present)
        tmpl_unsupported = sum(ev.get('template_unsupported', 0) for ev in piv.values())
        self.assertEqual(attempts, products + no_matches + tmpl_unsupported)


if __name__ == '__main__':
    unittest.main()

