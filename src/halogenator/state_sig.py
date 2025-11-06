# -*- coding: ascii -*-
"""Path-independent state signatures for enumeration pruning."""

from typing import List, Dict, Any, Optional, Tuple
import hashlib


def compute_state_signature(history: List[Dict[str, Any]]) -> str:
    """
    Compute path-independent signature for a substitution history.
    
    This creates a signature based on the set of substitutions, independent
    of the order in which they were applied. Useful for detecting equivalent
    states reached via different paths.
    
    Args:
        history: List of substitution history items
    
    Returns:
        Hex string signature
    """
    if not history:
        return ""
    
    # Extract key information from each history item
    substitution_tuples = []
    
    for item in history:
        # Create tuple with key substitution information
        rule = item.get('rule', '')
        sym = item.get('sym', 0)
        ring_tag = item.get('ring_tag')
        halogen = item.get('halogen', '')
        
        # Convert ring_tag to hashable form
        if isinstance(ring_tag, (list, tuple)):
            ring_tag_key = tuple(sorted(ring_tag)) if ring_tag else ()
        else:
            ring_tag_key = ring_tag if ring_tag is not None else ""
        
        # Create substitution signature tuple
        sub_tuple = (rule, sym, ring_tag_key, halogen)
        substitution_tuples.append(sub_tuple)
    
    # Sort tuples to make signature path-independent
    sorted_tuples = tuple(sorted(substitution_tuples))
    
    # Create hash
    signature_str = str(sorted_tuples)
    return hashlib.md5(signature_str.encode('ascii', errors='ignore')).hexdigest()


def compute_canonical_selection_key(rule: str, sym_class: int, ring_tag: Any) -> Tuple:
    """
    Compute canonical key for site selection ordering.
    
    This enables "non-decreasing group key selection" strategy to reduce
    path permutation redundancy by enforcing a canonical ordering of
    site group selection.
    
    Args:
        rule: Rule name (e.g., 'R1', 'R2')
        sym_class: Symmetry class of the site
        ring_tag: Ring tag for the site
    
    Returns:
        Canonical selection key tuple
    """
    # Convert ring_tag to hashable, sortable form
    if isinstance(ring_tag, (list, tuple)):
        ring_key = tuple(sorted(ring_tag)) if ring_tag else ()
    else:
        ring_key = ring_tag if ring_tag is not None else ""
    
    return (rule, sym_class, ring_key)


def should_prune_by_selection_order(current_key: Tuple, previous_key: Optional[Tuple]) -> bool:
    """
    Determine if current selection should be pruned based on canonical ordering.
    
    Returns True if current_key < previous_key (violates canonical order).
    This implements the "non-decreasing group key selection" strategy.
    
    Args:
        current_key: Key for current site group selection
        previous_key: Key from previous selection in path
    
    Returns:
        True if should prune (violates canonical order)
    """
    if previous_key is None:
        return False  # First selection is always allowed
    
    return current_key < previous_key


def is_equivalent_state(sig1: str, sig2: str) -> bool:
    """Check if two state signatures represent equivalent states."""
    return sig1 == sig2 if sig1 and sig2 else False