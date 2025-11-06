# -*- coding: ascii -*-
"""Unit tests for rule family computation."""

import unittest
from halogenator.schema import compute_rule_family


class TestComputeRuleFamily(unittest.TestCase):
    """Test compute_rule_family() function for rule grouping."""

    def test_r2a_maps_to_r2(self):
        """Test R2a sub-rule maps to R2 family."""
        result = compute_rule_family('R2a')
        self.assertEqual(result, 'R2')

    def test_r2b_maps_to_r2(self):
        """Test R2b sub-rule maps to R2 family."""
        result = compute_rule_family('R2b')
        self.assertEqual(result, 'R2')

    def test_r6_methyl_maps_to_r6(self):
        """Test R6_methyl sub-rule maps to R6 family."""
        result = compute_rule_family('R6_methyl')
        self.assertEqual(result, 'R6')

    def test_r1_unchanged(self):
        """Test R1 rule remains unchanged."""
        result = compute_rule_family('R1')
        self.assertEqual(result, 'R1')

    def test_r2_unchanged(self):
        """Test R2 base rule remains unchanged."""
        result = compute_rule_family('R2')
        self.assertEqual(result, 'R2')

    def test_r3_unchanged(self):
        """Test R3 rule remains unchanged."""
        result = compute_rule_family('R3')
        self.assertEqual(result, 'R3')

    def test_r4_unchanged(self):
        """Test R4 rule remains unchanged."""
        result = compute_rule_family('R4')
        self.assertEqual(result, 'R4')

    def test_r5_unchanged(self):
        """Test R5 rule remains unchanged."""
        result = compute_rule_family('R5')
        self.assertEqual(result, 'R5')

    def test_r6_unchanged(self):
        """Test R6 base rule remains unchanged."""
        result = compute_rule_family('R6')
        self.assertEqual(result, 'R6')

    def test_all_r2_variants(self):
        """Test batch R2 sub-rule mappings."""
        r2_variants = ['R2a', 'R2b']
        for variant in r2_variants:
            with self.subTest(variant=variant):
                result = compute_rule_family(variant)
                self.assertEqual(result, 'R2', f"{variant} should map to R2")

    def test_all_r6_variants(self):
        """Test batch R6 sub-rule mappings."""
        r6_variants = ['R6_methyl']
        for variant in r6_variants:
            with self.subTest(variant=variant):
                result = compute_rule_family(variant)
                self.assertEqual(result, 'R6', f"{variant} should map to R6")

    def test_base_rules_unchanged(self):
        """Test all base rules remain unchanged."""
        base_rules = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
        for rule in base_rules:
            with self.subTest(rule=rule):
                result = compute_rule_family(rule)
                self.assertEqual(result, rule, f"{rule} should remain unchanged")

    def test_idempotency_on_base_rules(self):
        """Test applying compute_rule_family twice gives same result for base rules."""
        base_rules = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
        for rule in base_rules:
            with self.subTest(rule=rule):
                first_result = compute_rule_family(rule)
                second_result = compute_rule_family(first_result)
                self.assertEqual(first_result, second_result,
                               f"Idempotency failed for {rule}")

    def test_idempotency_on_sub_rules(self):
        """Test applying compute_rule_family twice gives same result for sub-rules."""
        sub_rules = ['R2a', 'R2b', 'R6_methyl']
        for rule in sub_rules:
            with self.subTest(rule=rule):
                first_result = compute_rule_family(rule)
                second_result = compute_rule_family(first_result)
                self.assertEqual(first_result, second_result,
                               f"Idempotency failed for {rule}")


if __name__ == '__main__':
    unittest.main()
