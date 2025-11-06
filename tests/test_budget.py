# -*- coding: ascii -*-
"""
Unit tests for budget system (PR-3 foundation)
"""

import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.budget import BudgetState


class TestBudgetState(unittest.TestCase):
    """Test budget system for R6 methyl enumeration."""

    def test_ops_mode_basic(self):
        """Test basic ops mode functionality."""
        budget = BudgetState("ops", k_max=2)

        # Initial state
        self.assertEqual(budget.budget_mode, "ops")
        self.assertEqual(budget.k_max, 2)
        self.assertEqual(budget.k_ops, 0)
        self.assertEqual(budget.k_atoms, 0)

        # Charge 1 op, 1 atom - should succeed
        self.assertTrue(budget.charge(1, 1))
        self.assertEqual(budget.k_ops, 1)
        self.assertEqual(budget.k_atoms, 1)

        # Charge 1 more op, 2 more atoms - should succeed (ops=2)
        self.assertTrue(budget.charge(1, 2))
        self.assertEqual(budget.k_ops, 2)
        self.assertEqual(budget.k_atoms, 3)

        # Try to charge 1 more op - should fail (would exceed k_max=2 ops)
        self.assertFalse(budget.charge(1, 1))
        # State should not change on failed charge
        self.assertEqual(budget.k_ops, 2)
        self.assertEqual(budget.k_atoms, 3)

    def test_atoms_mode_basic(self):
        """Test basic atoms mode functionality."""
        budget = BudgetState("atoms", k_max=3)

        # Initial state
        self.assertEqual(budget.budget_mode, "atoms")
        self.assertEqual(budget.k_max, 3)

        # Charge 2 ops, 2 atoms - should succeed
        self.assertTrue(budget.charge(2, 2))
        self.assertEqual(budget.k_ops, 2)
        self.assertEqual(budget.k_atoms, 2)

        # Charge 1 op, 1 atom - should succeed (atoms=3)
        self.assertTrue(budget.charge(1, 1))
        self.assertEqual(budget.k_ops, 3)
        self.assertEqual(budget.k_atoms, 3)

        # Try to charge more atoms - should fail (would exceed k_max=3 atoms)
        self.assertFalse(budget.charge(1, 1))
        # State should not change on failed charge
        self.assertEqual(budget.k_ops, 3)
        self.assertEqual(budget.k_atoms, 3)

    def test_site_tokens_tracking(self):
        """Test site tokens for first-touch tracking."""
        budget = BudgetState("ops", k_max=5)

        # Initially empty
        self.assertEqual(len(budget.site_tokens), 0)

        # Set some site tokens
        budget.site_tokens[10] = 1
        budget.site_tokens[15] = 1

        self.assertEqual(budget.site_tokens[10], 1)
        self.assertEqual(budget.site_tokens[15], 1)
        self.assertEqual(budget.site_tokens.get(20, 0), 0)

    def test_site_state_tracking(self):
        """Test site state for multiset tracking."""
        budget = BudgetState("atoms", k_max=10)

        # Initially empty
        self.assertEqual(len(budget.site_state), 0)

        # Set some site states (sorted multisets)
        budget.site_state[5] = (("F", 1), ("Cl", 1))
        budget.site_state[8] = (("F", 2),)

        self.assertEqual(budget.site_state[5], (("F", 1), ("Cl", 1)))
        self.assertEqual(budget.site_state[8], (("F", 2),))

    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        # k_max=0 ops mode: only ops count matters
        budget = BudgetState("ops", k_max=0)
        self.assertFalse(budget.charge(1, 0))  # 1 op > 0, should fail
        self.assertTrue(budget.charge(0, 1))   # 0 ops <= 0, should succeed

        # k_max=0 atoms mode: only atoms count matters
        budget = BudgetState("atoms", k_max=0)
        self.assertTrue(budget.charge(1, 0))   # 0 atoms <= 0, should succeed
        budget = BudgetState("atoms", k_max=0)  # reset
        self.assertFalse(budget.charge(0, 1))  # 1 atom > 0, should fail

        # Zero charges should succeed
        budget = BudgetState("atoms", k_max=1)
        self.assertTrue(budget.charge(0, 0))
        self.assertEqual(budget.k_ops, 0)
        self.assertEqual(budget.k_atoms, 0)

    def test_invalid_budget_mode(self):
        """Test that invalid budget mode raises assertion error."""
        with self.assertRaises(AssertionError):
            BudgetState("invalid", k_max=2)


if __name__ == '__main__':
    unittest.main(verbosity=2)