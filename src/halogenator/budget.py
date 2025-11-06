# -*- coding: ascii -*-

from collections import Counter
from typing import Dict, Tuple, Optional

class BudgetState:
    def __init__(self, budget_mode: str, k_max: int):
        assert budget_mode in ("ops", "atoms")
        self.budget_mode = budget_mode
        self.k_max = k_max
        self.k_ops = 0
        self.k_atoms = 0
        self.site_tokens: Dict[int, int] = {}   # first-touch flags
        self.site_state: Dict[int, Counter] = {}  # cidx -> Counter({"F":nF, "Cl":nCl})

    def charge(self, cost_ops: int, cost_atoms: int) -> bool:
        next_ops = self.k_ops + cost_ops
        next_atoms = self.k_atoms + cost_atoms
        if self.budget_mode == "ops":
            if next_ops > self.k_max:
                return False
        else:  # atoms mode
            if next_atoms > self.k_max:
                return False
        # commit
        self.k_ops = next_ops
        self.k_atoms = next_atoms
        return True

    def token_cost(self, cidx: int) -> int:
        """
        Calculate the token cost for first touch of a site.
        Returns 1 if site not touched before, 0 if already touched.
        """
        return 1 - int(self.site_tokens.get(cidx, 0))

    def mark_first_touch(self, cidx: int):
        """Mark a site as touched for the first time."""
        self.site_tokens[cidx] = 1

    def get_state_key(self, cidx: int) -> Tuple[Tuple[str, int], ...]:
        """
        Get the state key for local-state dedup.
        Returns tuple sorted for local-state dedup, e.g. (("Cl",1),("F",2))
        """
        counter = self.site_state.get(cidx, Counter())
        return tuple(sorted(counter.items()))

    def bump_state(self, cidx: int, X: str):
        """
        Increment the count of halogen X at site cidx.
        """
        counter = self.site_state.setdefault(cidx, Counter())
        counter[X] += 1

    def copy(self) -> 'BudgetState':
        """Create a deep copy of the budget state."""
        new_budget = BudgetState(self.budget_mode, self.k_max)
        new_budget.k_ops = self.k_ops
        new_budget.k_atoms = self.k_atoms
        new_budget.site_tokens = self.site_tokens.copy()
        new_budget.site_state = {
            cidx: counter.copy()
            for cidx, counter in self.site_state.items()
        }
        return new_budget

    def current_budget(self) -> int:
        """Get the current budget value based on budget_mode."""
        return self.k_ops if self.budget_mode == "ops" else self.k_atoms

    def can_afford(self, cost_ops: int, cost_atoms: int) -> bool:
        """Check if we can afford the given cost without exceeding k_max."""
        if self.budget_mode == "ops":
            return self.k_ops + cost_ops <= self.k_max
        else:
            return self.k_atoms + cost_atoms <= self.k_max

    @classmethod
    def from_payload(cls, payload: dict, k_max: int) -> 'BudgetState':
        """Create BudgetState from serialized payload."""
        obj = cls(payload.get('mode', 'ops'), k_max)
        obj.k_ops = int(payload.get('k_ops', 0))
        obj.k_atoms = int(payload.get('k_atoms', 0))
        obj.site_tokens = dict(payload.get('site_tokens', {}))
        # site_state: {cidx: {X: n}}
        obj.site_state = {
            int(k): Counter(v)
            for k, v in payload.get('site_state', {}).items()
        }
        return obj

    def to_payload(self) -> dict:
        """Serialize BudgetState to lightweight payload for BFS passing."""
        return {
            'mode': self.budget_mode,
            'k_ops': self.k_ops,
            'k_atoms': self.k_atoms,
            'site_tokens': dict(self.site_tokens),
            'site_state': {
                cidx: dict(cnt)
                for cidx, cnt in self.site_state.items()
            },
        }