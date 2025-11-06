# -*- coding: ascii -*-
"""Centralized deduplication utilities for halogenator."""

from typing import Optional, Tuple, Set, Dict
from .standardize import to_inchikey_sanitized
from .qa_utils import qa_mark


def compute_stable_key(mol) -> Optional[str]:
    """
    Compute stable InChIKey for deduplication purposes.

    Uses sanitized InChIKey that strips isotope and atom mapping information
    for consistent deduplication across different enumeration paths.

    Args:
        mol: RDKit molecule object

    Returns:
        Sanitized InChIKey string or None if computation fails
    """
    try:
        return to_inchikey_sanitized(mol)
    except Exception:
        return None


def early_check(
    mol,
    seen: Set[str],
    qa_bus: Optional[Dict[str, int]] = None,
    metric: str = "dedup_hits_inchi",
    policy: str = "stable_key",
    state_sig: Optional[str] = None,
    attempt: Optional[Dict[str, int]] = None,
    enable: bool = True,
) -> Tuple[Optional[str], bool]:
    """
    Check if molecule is a duplicate early in the enumeration process.

    This function computes the dedup key according to the specified policy and
    checks if it's already in the seen set. If it's a duplicate, it updates the
    QA metrics but does NOT add the key to the seen set. The caller must call
    commit() separately after validation passes.

    CRITICAL: Uses qa_mark() to ensure totals and pivots stay synchronized.

    Args:
        mol: RDKit molecule object
        seen: Set of already seen dedup keys
        qa_bus: Optional QA metrics dictionary to update (totals)
        metric: Name of the QA metric to increment for duplicates
        policy: Deduplication policy ("stable_key", "state_sig", "none")
        state_sig: Pre-computed state signature (required if policy="state_sig")
        attempt: Optional per-attempt QA events dictionary (for pivots)
        enable: Master switch for deduplication (set False for raw mode)

    Returns:
        (dedup_key, is_duplicate) tuple
        - dedup_key: Computed dedup key or None if computation failed/policy="none"
        - is_duplicate: True if duplicate or failed to compute key (treated as drop)
    """
    # Master switch: if deduplication disabled, treat as policy="none" (raw mode)
    if not enable or policy == "none":
        return None, False

    key = None
    if policy == "stable_key":
        key = compute_stable_key(mol)
    elif policy == "state_sig":
        key = state_sig

    if key is None:
        qa_mark(qa_bus, attempt, "rdkit_error")
        return None, True  # drop invalid candidate early

    if key in seen:
        qa_mark(qa_bus, attempt, metric)
        return key, True

    return key, False


def commit(key: Optional[str], seen: Set[str]) -> None:
    """
    Commit a dedup key to the seen set after all validations pass.

    This function should only be called after the molecule has passed
    all constraints, post-guards, and other validation checks.

    Args:
        key: Dedup key to add to seen set (None keys are ignored)
        seen: Set of seen dedup keys to update
    """
    if key:
        seen.add(key)