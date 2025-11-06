# -*- coding: ascii -*-
"""
Unified QA counting utilities for halogenator.

This module provides consistent counting mechanisms to ensure QA totals
and pivot sums remain synchronized across all enumeration paths.

QA Event Granularity Definitions:
=================================

PER-ATTEMPT EVENTS (max 1 per attempt):
- rdkit_error: RDKit computation/sanitization failures
- dedup_hits_inchi: Product deduplicated by stable InChIKey
- dedup_hits_statesig: Product deduplicated by state signature
- isotope_unavailable: Isotope tagging strategy unavailable for this attempt
- isotope_miss: Isotope tagged site not found in reaction products
- atommap_used: Site detection used AtomMapNum fallback strategy
- heuristic_used: Site detection used heuristic fallback strategy
- sugar_mask_filtered: Sites filtered out by sugar masking
- sugar_mask_degraded: Sugar mask computation degraded/failed
- sugar_post_guard_blocked: Product blocked by post-enumeration sugar guard
- template_unsupported: Reaction template has unsupported structure

PER-PRODUCT EVENTS:
- Information written to record fields only, not counted as QA events
- Exception: sugar_post_guard_blocked is treated as per-attempt for historical compatibility

CRITICAL REQUIREMENT:
All per-attempt events MUST use qa_mark() to ensure totals and pivots stay synchronized.
Direct modification of qa_bus or attempt dicts is prohibited.
"""

from typing import Dict, Optional


def qa_mark(
    qa_bus: Optional[Dict[str, int]],
    attempt: Optional[Dict[str, int]],
    key: str,
    inc: int = 1,
) -> None:
    """
    Simultaneously write to total accounting (qa_bus) and per-attempt accounting (attempt).

    This function ensures QA totals and pivot sums remain synchronized by writing
    the same increment to both accounting systems atomically. Event keys are
    automatically canonicalized to ensure consistency across the system.

    Args:
        qa_bus: Global QA totals dictionary (may be None)
        attempt: Per-attempt QA events dictionary (may be None)
        key: QA metric key to increment (will be canonicalized)
        inc: Increment amount (default: 1)

    Usage:
        # Instead of:
        # qa_bus['rdkit_error'] += 1
        # attempt_qa['rdkit_error'] += 1

        # Use:
        qa_mark(qa_bus, attempt_qa, 'rdkit_error')
        # Legacy names are automatically converted:
        qa_mark(qa_bus, attempt_qa, 'post_guard_blocked')  # -> 'sugar_post_guard_blocked'
    """
    canonical_key = canonical_event(key)
    if qa_bus is not None:
        qa_bus[canonical_key] = qa_bus.get(canonical_key, 0) + inc
    if attempt is not None:
        attempt[canonical_key] = attempt.get(canonical_key, 0) + inc


def qa_ensure_keys(qa_dict: Dict[str, int], keys: list) -> None:
    """
    Ensure all required QA keys exist in dictionary with zero values.

    Args:
        qa_dict: QA dictionary to initialize
        keys: List of keys that must exist
    """
    for key in keys:
        if key not in qa_dict:
            qa_dict[key] = 0


# Centralized event key constants for maintainability
EV_ISOTOPE_UNAVAILABLE = 'isotope_unavailable'
EV_ISOTOPE_MISS = 'isotope_miss'
EV_ATOMMAP_USED = 'atommap_used'
EV_HEURISTIC_USED = 'heuristic_used'
EV_DEDUP_INCHI = 'dedup_hits_inchi'
EV_DEDUP_STATESIG = 'dedup_hits_statesig'
EV_SUGAR_MASK_FILTERED = 'sugar_mask_filtered'
EV_SUGAR_MASK_DEGRADED = 'sugar_mask_degraded'
EV_SUGAR_POST_GUARD = 'sugar_post_guard_blocked'
EV_SUGAR_PROXIMITY_FILTERED = 'sugar_proximity_filtered'
EV_RDKIT_ERROR = 'rdkit_error'

# Canonical event name mappings (old_key -> new_key)
# This provides single source of truth for event name standardization
ALIAS_EVENT_KEYS = {
    'post_guard_blocked': EV_SUGAR_POST_GUARD,
    'sugar_filtered_reaction_matches': EV_SUGAR_MASK_FILTERED,
    'pruned_inchikey_dupe': EV_DEDUP_INCHI,
}

# Global set to track all event keys seen during runtime for alias audit
SEEN_EVENT_KEYS = set()

def canonical_event(key: str) -> str:
    """
    Convert event key to its canonical form.

    Args:
        key: Event key that may be legacy or canonical

    Returns:
        Canonical event key name
    """
    # Track all keys seen for alias coverage audit
    SEEN_EVENT_KEYS.add(key)
    return ALIAS_EVENT_KEYS.get(key, key)

# Standard QA metrics that should always be present in final outputs
STANDARD_QA_TOTALS_KEYS = [
    'dedup_hits_statesig',
    'dedup_hits_inchi',
    'no_product_matches',
    'template_unsupported'
]

STANDARD_QA_PATHS_KEYS = [
    EV_ISOTOPE_UNAVAILABLE,
    EV_ISOTOPE_MISS,
    EV_ATOMMAP_USED,
    EV_HEURISTIC_USED,
    EV_SUGAR_MASK_FILTERED,
    EV_SUGAR_MASK_DEGRADED,
    EV_SUGAR_POST_GUARD,
    EV_RDKIT_ERROR
]

# Events that belong to pivot-side counting (will be aggregated by rule/halogen/k)
# All event names here are in canonical form
PIVOT_EVENT_KEYS = {
    EV_ISOTOPE_UNAVAILABLE, EV_ISOTOPE_MISS,
    EV_ATOMMAP_USED, EV_HEURISTIC_USED,
    EV_DEDUP_INCHI, EV_DEDUP_STATESIG,
    EV_SUGAR_MASK_FILTERED,
    EV_RDKIT_ERROR,
    EV_SUGAR_POST_GUARD,
}

# Derived sets for different use cases
# NON_PIVOT_PATH_KEYS = set(STANDARD_QA_PATHS_KEYS) - PIVOT_EVENT_KEYS  # Can be computed as needed for merge filtering


def get_known_event_keys() -> set:
    """
    Return the set of all known event keys (canonical constants + alias keys).

    This is used for alias coverage auditing to identify unknown keys.
    """
    canonical_constants = {
        EV_ISOTOPE_UNAVAILABLE, EV_ISOTOPE_MISS, EV_ATOMMAP_USED, EV_HEURISTIC_USED,
        EV_DEDUP_INCHI, EV_DEDUP_STATESIG, EV_SUGAR_MASK_FILTERED, EV_SUGAR_MASK_DEGRADED,
        EV_SUGAR_POST_GUARD, EV_SUGAR_PROXIMITY_FILTERED, EV_RDKIT_ERROR
    }

    alias_keys = set(ALIAS_EVENT_KEYS.keys())
    alias_values = set(ALIAS_EVENT_KEYS.values())

    # Also include legacy dedup keys that are handled specially
    legacy_dedup_keys = {'statesig_hits', 'inchi_hits'}

    # Include other standard event keys (not topevel totals like attempts/products)
    other_event_keys = {'template_unsupported', 'no_product_matches'}

    return canonical_constants | alias_keys | alias_values | legacy_dedup_keys | other_event_keys


def audit_alias_coverage() -> Dict[str, set]:
    """
    Audit alias coverage by comparing seen keys against known keys.

    Returns:
        Dictionary with 'unknown_keys' and 'coverage_stats'
    """
    known_keys = get_known_event_keys()
    unknown_keys = SEEN_EVENT_KEYS - known_keys

    return {
        'unknown_keys': unknown_keys,
        'coverage_stats': {
            'total_seen': len(SEEN_EVENT_KEYS),
            'known_keys': len(SEEN_EVENT_KEYS & known_keys),
            'unknown_keys': len(unknown_keys),
            'coverage_percentage': len(SEEN_EVENT_KEYS & known_keys) / len(SEEN_EVENT_KEYS) * 100 if SEEN_EVENT_KEYS else 100
        }
    }


def reset_seen_keys():
    """Reset the seen keys set for testing purposes."""
    global SEEN_EVENT_KEYS
    SEEN_EVENT_KEYS.clear()


def qa_init_standard_structure() -> Dict[str, any]:
    """
    Initialize a standard QA statistics structure with all required keys.

    Returns:
        Dictionary with standard QA structure and zero values
    """
    qa_stats = {}

    # Initialize top-level counters
    qa_ensure_keys(qa_stats, STANDARD_QA_TOTALS_KEYS)

    # Initialize qa_paths substructure
    qa_stats['qa_paths'] = {}
    qa_ensure_keys(qa_stats['qa_paths'], STANDARD_QA_PATHS_KEYS)

    return qa_stats