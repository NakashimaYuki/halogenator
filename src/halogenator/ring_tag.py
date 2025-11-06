# -*- coding: ascii -*-
"""
DEPRECATED: Stable ring labeling system for per-ring quotas.

This module is deprecated and no longer used in the codebase.
All ring tagging functionality has been moved to sites.py using
flavonoid_ring_label() for enhanced A/B/C ring identification.

Use sites.flavonoid_ring_label() and sites.get_ring_tag_for_atom() instead.
"""

from typing import List, Tuple, Optional, Dict


def ring_tags(mol) -> List[Tuple[int, ...]]:
    """
    Generate stable ring tags based on canonical ranks of ring atoms.
    Focus on 6-membered rings for R2 constraint tracking.
    
    Args:
        mol: RDKit molecule with ring info computed
    
    Returns:
        List of ring tags (one per 6-ring), where each tag is a sorted tuple
        of canonical ranks of atoms in that ring
    """
    try:
        from rdkit import Chem
    except Exception:
        return []
    
    try:
        ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
        ri = mol.GetRingInfo()
        tags = []
        for ring in ri.AtomRings():
            if len(ring) == 6:
                tags.append(tuple(sorted(int(ranks[i]) for i in ring)))
        return tags
    except Exception:
        return []


def ring_tag_for_atom(mol, atom_idx: int) -> Optional[Tuple[int, ...]]:
    """
    Get the ring tag for a specific atom (if it's in a ring).
    
    Args:
        mol: RDKit molecule
        atom_idx: Index of the atom
    
    Returns:
        Ring tag tuple if atom is in a ring, None otherwise
    """
    try:
        ring_info = mol.GetRingInfo()
        
        # Find which ring(s) this atom belongs to
        for i, ring_atoms in enumerate(ring_info.AtomRings()):
            if atom_idx in ring_atoms:
                # Generate tag for this ring
                from rdkit import Chem
                ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
                ring_ranks = [ranks[idx] for idx in ring_atoms]
                return tuple(sorted(ring_ranks))
        
        return None
        
    except Exception:
        return None


def ring_tag_hash(ring_tag: Tuple[int, ...]) -> int:
    """Convert ring tag to hash for efficient comparison."""
    return hash(ring_tag) if ring_tag else 0


def atom_to_ring_tags(mol) -> Dict[int, List[Tuple[int, ...]]]:
    """
    Map each atom to its 6-membered ring tags.
    
    Args:
        mol: RDKit molecule
        
    Returns:
        Dict mapping atom index to list of ring tags (for 6-rings only)
    """
    try:
        tags = ring_tags(mol)
        ri = mol.GetRingInfo()
        six_rings = [r for r in ri.AtomRings() if len(r)==6]
        m = {}
        for tag, ring in zip(tags, six_rings):
            s = set(ring)
            for idx in s:
                m.setdefault(idx, []).append(tag)
        return m
    except Exception:
        return {}