#!/usr/bin/env python
"""Test alias coverage audit mechanism."""

import unittest
from src.halogenator.qa_utils import (
    canonical_event, audit_alias_coverage, reset_seen_keys,
    get_known_event_keys, SEEN_EVENT_KEYS
)
from src.halogenator.enumerate_k import EnumConfig, enumerate_products


class TestAliasAudit(unittest.TestCase):
    """Test automated alias coverage auditing."""

    def setUp(self):
        """Reset seen keys before each test."""
        reset_seen_keys()

    def test_known_event_keys_coverage(self):
        """Test that get_known_event_keys returns comprehensive set."""
        known_keys = get_known_event_keys()

        # Should include canonical constants
        self.assertIn('isotope_unavailable', known_keys)
        self.assertIn('sugar_mask_filtered', known_keys)
        self.assertIn('rdkit_error', known_keys)

        # Should include alias keys
        self.assertIn('post_guard_blocked', known_keys)
        self.assertIn('pruned_inchikey_dupe', known_keys)

        # Should include legacy dedup keys
        self.assertIn('statesig_hits', known_keys)
        self.assertIn('inchi_hits', known_keys)

        # Should include standard keys
        self.assertIn('template_unsupported', known_keys)
        self.assertIn('no_product_matches', known_keys)

        # Verify minimum expected size
        self.assertGreaterEqual(len(known_keys), 15,
                               f"Expected at least 15 known keys, got {len(known_keys)}: {sorted(known_keys)}")

    def test_alias_audit_with_known_keys(self):
        """Test alias audit with only known keys."""
        # Use some known keys
        canonical_event('isotope_unavailable')  # canonical
        canonical_event('post_guard_blocked')   # alias
        canonical_event('statesig_hits')        # legacy

        audit_result = audit_alias_coverage()

        # Should have no unknown keys
        self.assertEqual(len(audit_result['unknown_keys']), 0,
                        f"Unexpected unknown keys: {audit_result['unknown_keys']}")

        # Coverage should be 100%
        self.assertEqual(audit_result['coverage_stats']['coverage_percentage'], 100.0)
        self.assertEqual(audit_result['coverage_stats']['total_seen'], 3)
        self.assertEqual(audit_result['coverage_stats']['known_keys'], 3)
        self.assertEqual(audit_result['coverage_stats']['unknown_keys'], 0)

    def test_alias_audit_with_unknown_keys(self):
        """Test alias audit with unknown keys that need mapping."""
        # Use some known keys
        canonical_event('isotope_unavailable')
        canonical_event('sugar_mask_filtered')

        # Use some hypothetical unknown keys
        canonical_event('unknown_test_event_1')
        canonical_event('legacy_unmapped_key')

        audit_result = audit_alias_coverage()

        # Should detect unknown keys
        expected_unknown = {'unknown_test_event_1', 'legacy_unmapped_key'}
        self.assertEqual(audit_result['unknown_keys'], expected_unknown)

        # Coverage should be less than 100%
        self.assertEqual(audit_result['coverage_stats']['total_seen'], 4)
        self.assertEqual(audit_result['coverage_stats']['known_keys'], 2)
        self.assertEqual(audit_result['coverage_stats']['unknown_keys'], 2)
        self.assertEqual(audit_result['coverage_stats']['coverage_percentage'], 50.0)

    def test_alias_audit_integration_with_enumeration(self):
        """Test alias audit after running actual enumeration."""
        # Run a minimal enumeration to collect real event keys
        cfg = EnumConfig(k_max=1, halogens=('F',))

        # Collect some products to trigger QA events
        products = list(enumerate_products('c1ccccc1O', cfg, return_qa_stats=True))

        # Should have collected some event keys
        self.assertGreater(len(SEEN_EVENT_KEYS), 0, "Should have seen some event keys from enumeration")

        audit_result = audit_alias_coverage()

        # T33: Event keys from enumeration should be known, but topevel totals (attempts/products) are expected unknowns
        unknown_keys = audit_result['unknown_keys']
        expected_unknown_totals = {'attempts', 'products'}  # These are topevel totals, not event keys

        unexpected_unknown = unknown_keys - expected_unknown_totals
        if unexpected_unknown:
            # If there are unexpected unknown keys, they might be legitimate new keys that need aliases
            print(f"WARNING: Found potentially unmapped event keys: {sorted(unexpected_unknown)}")
            print(f"Consider adding these to ALIAS_EVENT_KEYS if they are legacy aliases")
            print(f"Or add them to qa_utils constants if they are new canonical events")

        # Coverage should be high for actual event keys (allowing for expected totals to be unknown)
        total_seen = audit_result['coverage_stats']['total_seen']
        unknown_totals_count = len(unknown_keys & expected_unknown_totals)
        event_keys_seen = total_seen - unknown_totals_count
        known_event_keys = audit_result['coverage_stats']['known_keys']

        if event_keys_seen > 0:
            event_key_coverage = (known_event_keys / event_keys_seen) * 100
            self.assertGreaterEqual(event_key_coverage, 80.0,
                                   f"Event key coverage too low: {event_key_coverage}% - unexpected unknown: {unexpected_unknown}")

    def test_reset_functionality(self):
        """Test that reset_seen_keys works correctly."""
        canonical_event('test_key')
        self.assertIn('test_key', SEEN_EVENT_KEYS)

        reset_seen_keys()
        self.assertEqual(len(SEEN_EVENT_KEYS), 0)

        audit_result = audit_alias_coverage()
        self.assertEqual(audit_result['coverage_stats']['total_seen'], 0)


if __name__ == '__main__':
    unittest.main()