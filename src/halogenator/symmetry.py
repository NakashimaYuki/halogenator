# -*- coding: ascii -*-
"""
Symmetry computation on masked subgraphs for enumeration.

This module provides functionality to compute canonical atom ranks while excluding
certain atoms (masked atoms) from the symmetry analysis. This is essential for
proper symmetry folding when sugar atoms are masked from halogenation enumeration.

Algorithm:
==========
1. Create a copy of the molecule with masked atoms removed
2. Compute canonical ranks on the remaining subgraph
3. Map ranks back to original atom indices
4. Masked atoms get rank -1 (excluded from symmetry folding)
"""

import logging
from typing import Set, List, Dict, Any, Optional, Tuple
from .chem_compat import Chem, canonical_rank_atoms_safe

LOG = logging.getLogger(__name__)


def canonical_ranks_on_masked_subgraph(mol, mask: Set[int]) -> List[int]:
    """
    Compute canonical atom ranks considering only non-masked atoms.

    Args:
        mol: RDKit molecule object
        mask: Set of atom indices to exclude from symmetry computation

    Returns:
        List of ranks where result[i] is the canonical rank for atom i.
        Masked atoms get rank -1.
    """
    if not mask:
        # No masking - use standard ranking
        return canonical_rank_atoms_safe(mol, breakTies=False)

    try:
        num_atoms = mol.GetNumAtoms()
        ranks = [-1] * num_atoms  # Initialize all ranks to -1 (masked)

        # Get non-masked atom indices
        non_masked = [i for i in range(num_atoms) if i not in mask]

        if not non_masked:
            # All atoms are masked
            return ranks

        # Create subgraph with only non-masked atoms
        subgraph_mol, old_to_new = _create_subgraph(mol, non_masked)

        if subgraph_mol is None:
            LOG.debug("Failed to create subgraph, using fallback ranking")
            return _fallback_ranking(mol, mask)

        # Compute ranks on subgraph
        subgraph_ranks = canonical_rank_atoms_safe(subgraph_mol, breakTies=False)

        # Map subgraph ranks back to original atom indices
        for original_idx in non_masked:
            if original_idx in old_to_new:
                subgraph_idx = old_to_new[original_idx]
                if subgraph_idx < len(subgraph_ranks):
                    ranks[original_idx] = subgraph_ranks[subgraph_idx]

        return ranks

    except Exception as e:
        LOG.debug(f"Masked subgraph ranking failed: {e}, using fallback")
        return _fallback_ranking(mol, mask)


def _create_subgraph(mol, atom_indices: List[int]) -> Tuple[Optional[Any], Dict[int, int]]:
    """
    Create a subgraph molecule containing only specified atoms.

    Args:
        mol: Original molecule
        atom_indices: List of atom indices to include in subgraph

    Returns:
        (subgraph_molecule, old_to_new_mapping) or (None, {}) if creation fails
    """
    try:
        if not atom_indices:
            return None, {}

        # Create editable molecule
        new_mol = Chem.EditableMol(Chem.Mol())

        # Map original atom indices to new indices
        old_to_new = {}

        # Add atoms
        for new_idx, old_idx in enumerate(atom_indices):
            if old_idx >= mol.GetNumAtoms():
                continue

            old_atom = mol.GetAtomWithIdx(old_idx)
            new_atom = Chem.Atom(old_atom.GetAtomicNum())

            # Copy basic properties
            try:
                new_atom.SetFormalCharge(old_atom.GetFormalCharge())
                new_atom.SetIsAromatic(old_atom.GetIsAromatic())
                new_atom.SetNumExplicitHs(old_atom.GetNumExplicitHs())
            except Exception:
                # If property copying fails, continue with basic atom
                pass

            # Add atom to new molecule
            actual_new_idx = new_mol.AddAtom(new_atom)
            old_to_new[old_idx] = actual_new_idx

        # Add bonds between included atoms
        atom_indices_set = set(atom_indices)
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            # Only add bond if both atoms are in subgraph
            if begin_idx in atom_indices_set and end_idx in atom_indices_set:
                if begin_idx in old_to_new and end_idx in old_to_new:
                    new_begin = old_to_new[begin_idx]
                    new_end = old_to_new[end_idx]

                    try:
                        new_mol.AddBond(new_begin, new_end, bond.GetBondType())
                    except Exception:
                        # If bond addition fails, skip this bond
                        continue

        # Get final molecule
        result_mol = new_mol.GetMol()

        # Try to sanitize the subgraph
        try:
            Chem.SanitizeMol(result_mol)
        except Exception:
            # Sanitization may fail for disconnected subgraphs, but that's OK
            # We can still compute ranks on the structure
            pass

        return result_mol, old_to_new

    except Exception as e:
        LOG.debug(f"Subgraph creation failed: {e}")
        return None, {}


def _fallback_ranking(mol, mask: Set[int]) -> List[int]:
    """
    Fallback ranking strategy when subgraph creation fails.

    Simply computes standard ranks and sets masked atoms to -1.
    """
    try:
        standard_ranks = canonical_rank_atoms_safe(mol, breakTies=False)
        for i in mask:
            if i < len(standard_ranks):
                standard_ranks[i] = -1
        return standard_ranks
    except Exception:
        # Ultimate fallback: sequential numbering
        num_atoms = mol.GetNumAtoms()
        ranks = list(range(num_atoms))
        for i in mask:
            if i < len(ranks):
                ranks[i] = -1
        return ranks


def compute_grouping_key(rule: str, sym_class: int, ring_tag: Optional[str] = None,
                        extra_context: Optional[str] = None) -> Tuple:
    """
    Compute stable grouping key for symmetry folding.

    Different rules use different grouping strategies to ensure proper
    symmetry folding without incorrect merging.

    Args:
        rule: Rule identifier (R1, R2a, R2b, R6, etc.)
        sym_class: Symmetry class from canonical ranking
        ring_tag: Optional ring identifier
        extra_context: Optional extra context for specialized rules

    Returns:
        Tuple that can be used as dictionary key for grouping
    """
    if rule == 'R1':
        # R1: aromatic CH sites - group by symmetry and ring
        return (rule, sym_class, ring_tag)

    elif rule == 'R2a':
        # R2a: sp2-CH in C ring - group by symmetry and ring
        return (rule, sym_class, ring_tag)

    elif rule == 'R2b':
        # R2b: flavanone C3-CH2 - group by symmetry, ring, and CH2 marker
        return (rule, sym_class, ring_tag, 'CH2')

    elif rule == 'R6':
        # R6: methyl sites - group by symmetry, neighbor class, and ring
        return (rule, sym_class, extra_context, ring_tag)

    else:
        # Default: group by rule and symmetry only
        return (rule, sym_class)


def pick_representatives_by_grouping(sites: List[int], ranks: List[int],
                                   rule: str, ring_label_map: Optional[Dict[int, str]] = None,
                                   neighbor_classes: Optional[Dict[int, str]] = None) -> List[int]:
    """
    Pick representative sites for symmetry folding using rule-specific grouping.

    Args:
        sites: List of candidate site atom indices
        ranks: Canonical ranks for all atoms (masked atoms have rank -1)
        rule: Rule identifier for grouping strategy
        ring_label_map: Optional mapping from atom index to ring label
        neighbor_classes: Optional mapping from atom index to neighbor class

    Returns:
        List of representative site indices (one per symmetry group)
    """
    if not sites:
        return []

    # Group sites by their grouping key
    groups = {}

    for site in sites:
        if site >= len(ranks):
            continue

        sym_class = ranks[site]
        if sym_class == -1:
            # This site is masked - skip it
            continue

        ring_tag = ring_label_map.get(site) if ring_label_map else None
        neighbor_class = neighbor_classes.get(site) if neighbor_classes else None

        key = compute_grouping_key(rule, sym_class, ring_tag, neighbor_class)

        if key not in groups:
            groups[key] = []
        groups[key].append(site)

    # Pick first representative from each group
    representatives = []
    for group_sites in groups.values():
        if group_sites:
            # Sort for deterministic selection
            group_sites.sort()
            representatives.append(group_sites[0])

    return representatives


def validate_masked_symmetry(mol, mask: Set[int], original_ranks: List[int]) -> bool:
    """
    Validate that masked symmetry computation is consistent.

    This is a debugging/testing function to ensure masked symmetry
    produces reasonable results.

    Args:
        mol: Molecule object
        mask: Set of masked atom indices
        original_ranks: Ranks computed on full molecule

    Returns:
        True if masked symmetry appears consistent
    """
    try:
        masked_ranks = canonical_ranks_on_masked_subgraph(mol, mask)

        # Basic consistency checks
        if len(masked_ranks) != len(original_ranks):
            return False

        # All masked atoms should have rank -1
        for atom_idx in mask:
            if atom_idx < len(masked_ranks) and masked_ranks[atom_idx] != -1:
                return False

        # Non-masked atoms should have valid ranks >= 0
        non_masked_ranks = [r for i, r in enumerate(masked_ranks)
                           if i not in mask and r >= 0]

        if non_masked_ranks:
            # Should have reasonable rank distribution
            unique_ranks = set(non_masked_ranks)
            if len(unique_ranks) == 0:
                return False

        return True

    except Exception as e:
        LOG.debug(f"Masked symmetry validation failed: {e}")
        return False