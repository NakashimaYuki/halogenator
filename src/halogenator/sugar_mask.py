# -*- coding: ascii -*-
"""
Sugar masking for halogenation enumeration.

This module provides functionality to identify and mask sugar atoms in flavonoid glycosides
to prevent combinatorial explosion during halogenation enumeration.

Heuristic Rules for Sugar Identification:
========================================
1. Sugar Rings: 5/6-member rings that are:
   - Non-aromatic
   - Contain exactly 1 ring oxygen
   - Have no inner carbonyl (C=O where C is in ring)
   - Have >= 2 exocyclic oxygen substituents

2. Glycosidic Bridge Oxygen: Oxygen atoms that are:
   - Degree 2 (exactly 2 bonds)
   - Bridge between sugar ring and non-sugar ring
"""

import logging
from typing import Set, Dict, Any, List, Tuple
from .chem_compat import Chem

LOG = logging.getLogger(__name__)


def _get_default_sugar_cfg() -> Dict[str, Any]:
    """
    Get default sugar masking configuration with parameterized thresholds.

    Returns:
        Dictionary with default configuration values
    """
    return {
        'mask_glycosidic_bridge_oxygen': True,
        'mask_exocyclic_oxygen': True,
        # P0-H1: Parameterized thresholds for fine-tuning
        'sp3_threshold_5_ring': 0.4,  # 5-member rings: >=40% sp3 atoms
        'sp3_threshold_6_ring': 0.5,  # 6-member rings: >=50% sp3 atoms
        'few_masked_atoms_threshold': 5,  # Threshold for "very few" masked atoms
        'allow_one_hop_cglyco': True,  # Enable one-hop C-glycoside topology detection
        'cglyco_minimal_mask_include_intermediate': False,  # Include intermediate carbon in minimal mask
        # P1-1: Evidence-based scoring system for main pathway sugar ring detection
        'sugar_ring_score_threshold': 8.0,  # Minimum score for sugar ring acceptance in main pathway
    }


def get_sugar_mask_with_full_status(mol, mode: str = 'heuristic', sugar_cfg: Dict[str, Any] = None) -> Tuple[Set[int], bool, Dict[str, Any]]:
    """
    Get set of atom indices to mask with comprehensive status information.

    Args:
        mol: RDKit molecule object
        mode: Masking strategy
            - 'off': No masking (empty set)
            - 'heuristic': Rule-based sugar identification
            - 'sru': External SRU service (fallback to heuristic if unavailable)
        sugar_cfg: Sugar configuration dictionary with detail switches

    Returns:
        Tuple of (masked_atom_indices, degraded_flag, status_metadata)
        status_metadata includes: degraded_reason, fallback_ring_count, cglyco_candidates_count
    """
    # Merge user config with defaults
    default_cfg = _get_default_sugar_cfg()
    if sugar_cfg is None:
        sugar_cfg = default_cfg
    else:
        # Merge user config with defaults (user config takes precedence)
        merged_cfg = default_cfg.copy()
        merged_cfg.update(sugar_cfg)
        sugar_cfg = merged_cfg

    if mode == 'off':
        return set(), False, {}

    if mode == 'sru':
        # Try SRU service first, fallback to heuristic
        try:
            mask_result = _get_sugar_mask_sru(mol, sugar_cfg)
            return mask_result, False, {}
        except Exception as e:
            LOG.debug(f"SRU masking failed, falling back to heuristic: {e}")
            mode = 'heuristic'

    if mode == 'heuristic':
        mask, degraded = _get_sugar_mask_heuristic(mol, sugar_cfg)

        # Collect detailed status information for degraded masking
        status_metadata = {}
        if degraded:
            try:
                # Re-run degraded analysis to get detailed info
                _, _, reasons = _perform_degraded_masking(mol, sugar_cfg)
                status_metadata['degraded_reason'] = reasons

                # Count fallback rings
                fallback_rings = _find_sugar_like_oxygen_rings_fallback(mol, sugar_cfg)
                status_metadata['fallback_ring_count'] = len(fallback_rings)

                # Count C-glycoside candidates and collect audit statistics
                cglyco_atoms, cglyco_audit = _detect_c_glycoside_topo(mol, sugar_cfg)
                status_metadata['cglyco_candidates_count'] = len(cglyco_atoms)
                status_metadata['cglyco_audit'] = cglyco_audit

            except Exception as e:
                LOG.debug(f"Failed to collect degraded status metadata: {e}")
                status_metadata = {'degraded_reason': ['unknown'], 'fallback_ring_count': 0, 'cglyco_candidates_count': 0}

        return mask, degraded, status_metadata

    LOG.warning(f"Unknown sugar masking mode: {mode}, using 'off'")
    return set(), False, {}


def get_sugar_mask_with_status(mol, mode: str = 'heuristic', sugar_cfg: Dict[str, Any] = None) -> Tuple[Set[int], bool]:
    """
    Get set of atom indices to mask from halogenation enumeration with status.

    Args:
        mol: RDKit molecule object
        mode: Masking strategy
            - 'off': No masking (empty set)
            - 'heuristic': Rule-based sugar identification
            - 'sru': External SRU service (fallback to heuristic if unavailable)
        sugar_cfg: Sugar configuration dictionary with detail switches

    Returns:
        Tuple of (masked_atom_indices, degraded_flag)
        degraded_flag is True if masking failed and fell back to minimal safe masking
    """
    mask, degraded, _ = get_sugar_mask_with_full_status(mol, mode, sugar_cfg)
    return mask, degraded


def get_sugar_mask(mol, mode: str = 'heuristic', sugar_cfg: Dict[str, Any] = None) -> Set[int]:
    """
    Get set of atom indices to mask from halogenation enumeration.

    Args:
        mol: RDKit molecule object
        mode: Masking strategy ('off', 'heuristic', 'sru')
        sugar_cfg: Sugar configuration dictionary with detail switches

    Returns:
        Set of atom indices to mask
    """
    mask, _ = get_sugar_mask_with_status(mol, mode, sugar_cfg)
    return mask


def post_guard_blocked(product_mol, mask_atoms: Set[int]) -> bool:
    """
    Check if a product molecule contains halogenation in masked regions.

    Args:
        product_mol: RDKit molecule object (enumeration product)
        mask_atoms: Set of atom indices that should have been masked

    Returns:
        True if product violates masking constraints (should be rejected)
    """
    if not mask_atoms:
        return False

    try:
        # Check if any halogen atoms are in masked regions or bonded to masked atoms
        for atom in product_mol.GetAtoms():
            if atom.GetSymbol() in {'F', 'Cl', 'Br', 'I'}:
                atom_idx = atom.GetIdx()

                # Direct violation: halogen in mask
                if atom_idx in mask_atoms:
                    LOG.debug(f"Post-guard violation: halogen at masked atom {atom_idx}")
                    return True

                # Indirect violation: halogen bonded to masked atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in mask_atoms:
                        LOG.debug(f"Post-guard violation: halogen {atom_idx} bonded to masked atom {neighbor.GetIdx()}")
                        return True

        return False

    except Exception as e:
        LOG.debug(f"Post-guard check failed: {e}")
        # Conservative: if we can't check, assume violation
        return True


def compute_sugar_audit_fields(mol, mask_atoms: Set[int], sugar_cfg: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Compute audit fields for sugar masking QA reporting.

    Args:
        mol: RDKit molecule object
        mask_atoms: Set of masked atom indices
        sugar_cfg: Sugar configuration dictionary (optional, defaults to standard config)

    Returns:
        Dictionary with sugar masking audit information including C-glycoside detection
    """
    try:
        if not mask_atoms:
            return {
                'sugar_mask_atoms': [],
                'sugar_rings': [],
                'masked_oh_count': 0,
                'masked_bridge_o_count': 0,
                'is_c_glycoside_like': False,
                # Observation fields
                'accepted_via_score': False,
                'accepted_ring_size': 0,
                'accepted_ring_score': 0.0,
                'exocyclic_O_count_single_bond': 0,
                'has_cglyco_evidence': False
            }

        # Ensure sugar_cfg is available
        if sugar_cfg is None:
            sugar_cfg = _get_default_sugar_cfg()

        # Find sugar rings
        sugar_rings = _find_sugar_rings(mol, sugar_cfg)

        # Detect C-glycoside patterns and compute observation fields on best scoring ring
        is_c_glycoside_like = False
        best_score = 0.0
        best_ring = None
        best_ring_size = 0
        best_ring_exo_o = 0
        best_ring_has_cglyco = False

        # Local scorer replicating main scoring to capture score for audit
        def _score_ring(ring_tuple):
            score = 0.0
            ring_size = len(ring_tuple)
            if ring_size == 6:
                score += 2.0
            elif ring_size == 5:
                score += 1.0

            # Exactly 1 ring oxygen required
            ring_oxygens = 0
            for aidx in ring_tuple:
                a = mol.GetAtomWithIdx(aidx)
                if a.GetSymbol() == 'O':
                    ring_oxygens += 1
            if ring_oxygens != 1:
                return -1.0, 0, False  # disqualify

            # No inner carbonyl
            if _has_inner_carbonyl(mol, ring_tuple):
                return -1.0, 0, False
            else:
                score += 2.0

            # SP3 ratio scoring
            sp3_count = 0
            total_carbons = 0
            for aidx in ring_tuple:
                a = mol.GetAtomWithIdx(aidx)
                if a.GetSymbol() == 'C':
                    total_carbons += 1
                    if a.GetHybridization() == Chem.HybridizationType.SP3:
                        sp3_count += 1
            if total_carbons > 0:
                sp3_ratio = sp3_count / total_carbons
                if ring_size == 5:
                    threshold = sugar_cfg.get('sp3_threshold_5_ring', 0.4)
                    if sp3_ratio >= threshold:
                        score += min(3.0, (sp3_ratio - threshold) / (1.0 - threshold) * 3.0)
                elif ring_size == 6:
                    threshold = sugar_cfg.get('sp3_threshold_6_ring', 0.5)
                    if sp3_ratio >= threshold:
                        score += min(3.0, (sp3_ratio - threshold) / (1.0 - threshold) * 3.0)

            # Exocyclic oxygens (single bonds only)
            exo_o = _count_exocyclic_oxygens_single_bond(mol, ring_tuple)
            score += min(4.0, exo_o * 1.0)

            # C-glycoside topology
            has_cglyco = _has_c_glycoside_pattern(mol, ring_tuple)
            if has_cglyco:
                score += 2.0

            return score, exo_o, has_cglyco

        for ring in sugar_rings:
            ring_score, exo_o_cnt, has_cglyco = _score_ring(ring)
            # Track any cglyco evidence overall
            if has_cglyco:
                is_c_glycoside_like = True
            if ring_score > best_score:
                best_score = ring_score
                best_ring = ring
                best_ring_size = len(ring)
                best_ring_exo_o = exo_o_cnt
                best_ring_has_cglyco = has_cglyco

        # Count exocyclic oxygens in mask
        exocyclic_o_count = 0
        bridge_o_count = 0

        for atom_idx in mask_atoms:
            if atom_idx >= mol.GetNumAtoms():
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                if atom.GetDegree() == 2:
                    # Could be bridge oxygen
                    bridge_o_count += 1
                elif atom.GetDegree() == 1:
                    # Could be OH group
                    exocyclic_o_count += 1

        return {
            'sugar_mask_atoms': sorted(list(mask_atoms)),
            'sugar_rings': [list(ring) for ring in sugar_rings],
            'masked_oh_count': exocyclic_o_count,
            'masked_bridge_o_count': bridge_o_count,
            'is_c_glycoside_like': is_c_glycoside_like,
            # Observation fields
            'accepted_via_score': bool(sugar_rings),
            'accepted_ring_size': int(best_ring_size) if best_ring is not None else 0,
            'accepted_ring_score': float(best_score) if best_ring is not None else 0.0,
            'exocyclic_O_count_single_bond': int(best_ring_exo_o) if best_ring is not None else 0,
            'has_cglyco_evidence': bool(best_ring_has_cglyco)
        }

    except Exception as e:
        LOG.debug(f"Sugar audit field computation failed: {e}")
        return {
            'sugar_mask_atoms': sorted(list(mask_atoms)),
            'sugar_rings': [],
            'masked_oh_count': 0,
            'masked_bridge_o_count': 0,
            'is_c_glycoside_like': False,
            # Observation fields (fallback defaults)
            'accepted_via_score': False,
            'accepted_ring_size': 0,
            'accepted_ring_score': 0.0,
            'exocyclic_O_count_single_bond': 0,
            'has_cglyco_evidence': False
        }


def _get_sugar_mask_heuristic(mol, sugar_cfg: Dict[str, Any]) -> Tuple[Set[int], bool]:
    """
    Heuristic-based sugar masking using ring analysis and connectivity.
    """
    masked_atoms = set()

    try:
        # Ensure molecule has ring info
        mol.GetRingInfo()

        # Step 1: Identify sugar rings
        sugar_rings = _find_sugar_rings(mol, sugar_cfg)
        LOG.debug(f"Found {len(sugar_rings)} sugar rings")

        # Check if we need degraded masking due to zero sugar rings
        if not sugar_rings:
            LOG.debug("No sugar rings found - triggering degraded masking")
            minimal_mask, degraded, reasons = _perform_degraded_masking(mol, sugar_cfg)
            LOG.debug(f"Degraded masking result: degraded={degraded}, reasons={reasons}")
            return minimal_mask, degraded

        # Mask all atoms in sugar rings
        for ring in sugar_rings:
            masked_atoms.update(ring)

        # Step 2: Identify glycosidic bridge oxygens (if enabled)
        if sugar_cfg.get('mask_glycosidic_bridge_oxygen', True):
            bridge_oxygens = _find_glycosidic_bridge_oxygens(mol, sugar_rings)
            LOG.debug(f"Found {len(bridge_oxygens)} bridge oxygens")
            masked_atoms.update(bridge_oxygens)
        else:
            LOG.debug("Glycosidic bridge oxygen masking disabled")

        # Step 3: Mask exocyclic oxygen substituents on sugar rings (if enabled)
        if sugar_cfg.get('mask_exocyclic_oxygen', True):
            exocyclic_oxygens = _find_exocyclic_oxygens(mol, sugar_rings)
            LOG.debug(f"Found {len(exocyclic_oxygens)} exocyclic oxygens")
            masked_atoms.update(exocyclic_oxygens)
        else:
            LOG.debug("Exocyclic oxygen masking disabled")

        # Step 4: Check for potential C-glycoside patterns that standard ring detection might miss
        # If we found sugar rings but very few masked atoms, try C-glycoside topology as supplement
        few_threshold = sugar_cfg.get('few_masked_atoms_threshold', 5)
        if len(masked_atoms) < few_threshold:
            try:
                c_glycoside_atoms, _ = _detect_c_glycoside_topo(mol, sugar_cfg)
                if c_glycoside_atoms:
                    LOG.debug(f"Supplementing standard masking with {len(c_glycoside_atoms)} C-glycoside topology atoms")

                    # If C-glycoside topology found significant evidence, consider using degraded masking instead
                    if len(c_glycoside_atoms) >= len(masked_atoms):
                        LOG.debug("C-glycoside topology found more evidence than standard detection - switching to degraded masking")
                        minimal_mask, degraded, reasons = _perform_degraded_masking(mol, sugar_cfg)
                        LOG.debug(f"Degraded masking via C-glycoside evidence: degraded={degraded}, reasons={reasons}")
                        return minimal_mask, degraded
                    else:
                        # Just supplement existing mask
                        masked_atoms.update(c_glycoside_atoms)
            except Exception as e:
                LOG.debug(f"C-glycoside topology supplementing failed: {e}")

        LOG.debug(f"Sugar masking: {len(masked_atoms)} total atoms masked")

    except Exception as e:
        LOG.debug(f"Heuristic sugar masking failed: {e}")
        # Implement minimal safe masking instead of empty set
        minimal_mask = set()
        degraded = True

        try:
            # Try to identify at least bridge oxygens as minimal safety measure
            mol.GetRingInfo()  # Ensure ring info is available

            # Try to find basic sugar rings first
            try:
                sugar_rings = _find_sugar_rings(mol, sugar_cfg)
                if sugar_rings:
                    # Mask atoms in identified sugar rings
                    for ring in sugar_rings:
                        minimal_mask.update(ring)
                    LOG.debug(f"Minimal masking: found {len(sugar_rings)} sugar rings")
            except Exception:
                pass

            # Try to find bridge oxygens as backup - use independent heuristics
            try:
                # First try context-aware detection if we have any rings
                if minimal_mask:  # If we found some rings, try context-aware detection
                    bridge_oxygens = _find_glycosidic_bridge_oxygens(mol, [])
                    minimal_mask.update(bridge_oxygens)
                    LOG.debug(f"Minimal masking: found {len(bridge_oxygens)} context-aware bridge oxygens")

                # Always try independent detection as additional safety
                independent_bridge_oxygens = _find_independent_bridge_oxygens(mol)
                minimal_mask.update(independent_bridge_oxygens)
                LOG.debug(f"Minimal masking: found {len(independent_bridge_oxygens)} independent bridge oxygens")
            except Exception:
                pass

            # Try to find anomeric carbons for C-glycosides if we have sugar rings
            try:
                if sugar_rings:
                    anomeric_carbons = _find_anomeric_carbons_for_masking(mol, sugar_rings)
                    minimal_mask.update(anomeric_carbons)
                    LOG.debug(f"Minimal masking: found {len(anomeric_carbons)} anomeric carbons for C-glycosides")
            except Exception:
                pass

        except Exception as fallback_e:
            LOG.debug(f"Minimal safe masking also failed: {fallback_e}")
            # Final fallback - truly empty set but marked as degraded
            minimal_mask = set()

        LOG.debug(f"Sugar masking degraded: {len(minimal_mask)} atoms in minimal mask")
        return minimal_mask, degraded

    return masked_atoms, False


def get_sugar_mask(mol, mode: str = 'heuristic', sugar_cfg: Dict[str, Any] = None) -> Set[int]:
    """
    Get set of atom indices to mask from halogenation enumeration.

    Args:
        mol: RDKit molecule object
        mode: Masking strategy
            - 'off': No masking (empty set)
            - 'heuristic': Rule-based sugar identification
            - 'sru': External SRU service (fallback to heuristic if unavailable)
        sugar_cfg: Sugar configuration dictionary with detail switches

    Returns:
        Set of masked atom indices (backwards compatible API)
    """
    mask, _ = get_sugar_mask_with_status(mol, mode=mode, sugar_cfg=sugar_cfg)
    return mask


def _find_sugar_rings(mol, sugar_cfg: Dict[str, Any] = None) -> List[Tuple[int, ...]]:
    """
    Find rings that match sugar criteria using evidence-based scoring system (P1-1).

    Scoring factors:
    - Ring size (5/6-member): +2/+1 points
    - Exactly 1 ring oxygen: +3 points
    - No inner carbonyl: +2 points
    - SP3 ratio (5-ring >=0.4, 6-ring >=0.5): linear scoring
    - Exocyclic oxygens (single bond): points per oxygen
    - Anomeric-like topology: +2 points

    Rings with score >= threshold are accepted as sugar rings.
    """
    # Merge default config
    if sugar_cfg is None:
        sugar_cfg = _get_default_sugar_cfg()

    # P1-1: Configurable evidence score threshold for main pathway acceptance
    evidence_threshold = sugar_cfg.get('sugar_ring_score_threshold', 8.0)

    sugar_rings = []

    try:
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        for ring in atom_rings:
            ring_size = len(ring)
            if ring_size not in (5, 6):
                continue

            # Check if ring is non-aromatic
            if _is_ring_aromatic(mol, ring):
                continue

            # P1-1: Calculate evidence score for this ring
            score = 0.0

            # Factor 1: Ring size preference (6-member pyranose > 5-member furanose)
            if ring_size == 6:
                score += 2.0
            elif ring_size == 5:
                score += 1.0

            # Factor 2: Count ring oxygens
            ring_oxygens = 0
            ring_oxygen_idx = None
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O':
                    ring_oxygens += 1
                    ring_oxygen_idx = atom_idx

            if ring_oxygens == 1:
                score += 3.0  # Exactly 1 ring oxygen is ideal
            else:
                continue  # Must have exactly 1 ring oxygen

            # Factor 3: No inner carbonyl
            if not _has_inner_carbonyl(mol, ring):
                score += 2.0
            else:
                continue  # Inner carbonyl disqualifies

            # Factor 4: SP3 ratio scoring (higher SP3 ratio = more sugar-like)
            sp3_count = 0
            total_carbons = 0
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'C':
                    total_carbons += 1
                    if atom.GetHybridization() == Chem.HybridizationType.SP3:
                        sp3_count += 1

            if total_carbons > 0:
                sp3_ratio = sp3_count / total_carbons
                # Linear scoring based on SP3 ratio thresholds
                if ring_size == 5:
                    threshold = sugar_cfg.get('sp3_threshold_5_ring', 0.4)
                    if sp3_ratio >= threshold:
                        score += min(3.0, (sp3_ratio - threshold) / (1.0 - threshold) * 3.0)
                elif ring_size == 6:
                    threshold = sugar_cfg.get('sp3_threshold_6_ring', 0.5)
                    if sp3_ratio >= threshold:
                        score += min(3.0, (sp3_ratio - threshold) / (1.0 - threshold) * 3.0)

            # Factor 5: Exocyclic oxygens (single bonds only, to avoid counting double-bonded O)
            exocyclic_o_count = _count_exocyclic_oxygens_single_bond(mol, ring)
            score += min(4.0, exocyclic_o_count * 1.0)  # Cap at 4 points

            # Factor 6: Anomeric-like topology (C-glycoside evidence)
            has_cglyco_evidence = _has_c_glycoside_pattern(mol, ring)
            if has_cglyco_evidence:
                score += 2.0

            LOG.debug(f"Ring {ring} evidence score: {score:.1f} (threshold: {evidence_threshold})")

            # P1-1: Accept ring if score meets evidence threshold AND strong sugar evidence
            # Strong evidence requires either:
            # - >=2 exocyclic oxygens OR C-glycoside evidence (original criteria)
            # - Very high score indicating clear sugar pattern (relaxed criteria for common sugars)
            has_strong_evidence = (exocyclic_o_count >= 2) or has_cglyco_evidence
            has_very_high_score = score >= (evidence_threshold + 5.0)  # High confidence threshold

            if score >= evidence_threshold and (has_strong_evidence or has_very_high_score):
                sugar_rings.append(ring)
                evidence_type = "strong" if has_strong_evidence else "high_score"
                LOG.debug(f"Sugar ring accepted via evidence scoring ({evidence_type}): {ring} (score: {score:.1f}, exo_O: {exocyclic_o_count}, cglyco: {has_cglyco_evidence})")
            elif score >= evidence_threshold:
                LOG.debug(f"Sugar ring rejected due to weak evidence: {ring} (score: {score:.1f}, exo_O: {exocyclic_o_count}, cglyco: {has_cglyco_evidence})")

    except Exception as e:
        LOG.debug(f"Sugar ring detection failed: {e}")

    return sugar_rings


def _is_ring_aromatic(mol, ring: Tuple[int, ...]) -> bool:
    """Check if ring is aromatic."""
    try:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsAromatic():
                return True
        return False
    except Exception:
        return False


def _has_inner_carbonyl(mol, ring: Tuple[int, ...]) -> bool:
    """Check if ring contains a carbonyl C=O bond where C is in the ring."""
    try:
        ring_atoms = set(ring)

        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()

                # Check for C=O where C is in ring
                if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O' and
                    begin_atom.GetIdx() in ring_atoms):
                    return True
                elif (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C' and
                      end_atom.GetIdx() in ring_atoms):
                    return True

        return False
    except Exception:
        return False


def _count_exocyclic_oxygens_single_bond(mol, ring: Tuple[int, ...]) -> int:
    """Count exocyclic oxygen substituents connected via single bonds only (P1-1 scoring)."""
    try:
        ring_atoms = set(ring)
        seen_o = set()
        exocyclic_count = 0

        # Pre-fetch ring info for fallback checks
        ring_info = mol.GetRingInfo()

        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)

            # Count oxygen neighbors that are not in the ring and connected by single bond
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() != 'O':
                    continue
                n_idx = neighbor.GetIdx()
                if n_idx in ring_atoms:
                    continue

                # Check bond type - only count single bonds
                bond = mol.GetBondBetweenAtoms(atom_idx, n_idx)
                if not bond or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue

                # Check if this oxygen should be excluded
                # Only exclude ring oxygens that are simple cyclic ethers without glycosidic evidence
                try:
                    in_ring = neighbor.IsInRing()  # Fast path when available
                except Exception:
                    # Fallback to ring_info if IsInRing is unavailable
                    try:
                        in_ring = bool(ring_info.IsAtomInRing(n_idx))
                    except Exception:
                        in_ring = False

                if in_ring:
                    # Check if this is a simple cyclic ether (both neighbors are ring carbons in same ring)
                    # If so, exclude it. If it has glycosidic patterns, keep it as sugar evidence.
                    o_neighbors = list(neighbor.GetNeighbors())
                    if len(o_neighbors) == 2:
                        o_neighbor_idxs = [n.GetIdx() for n in o_neighbors]
                        # Check if both neighbors are carbons and in rings
                        both_carbon_and_ring = all(
                            mol.GetAtomWithIdx(nidx).GetSymbol() == 'C' and
                            mol.GetAtomWithIdx(nidx).IsInRing()
                            for nidx in o_neighbor_idxs
                        )
                        # If both neighbors are ring carbons, this is likely a simple cyclic ether
                        if both_carbon_and_ring:
                            continue  # Exclude simple ring ethers
                    # Otherwise, keep ring oxygens that might be glycosidic links

                # De-duplicate the same oxygen connected to multiple ring atoms
                if n_idx not in seen_o:
                    seen_o.add(n_idx)
                    exocyclic_count += 1

        return exocyclic_count
    except Exception:
        return 0


def _count_exocyclic_oxygens(mol, ring: Tuple[int, ...]) -> int:
    """Count exocyclic oxygen substituents on ring atoms."""
    try:
        ring_atoms = set(ring)
        exocyclic_count = 0

        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)

            # Count oxygen neighbors that are not in the ring
            for neighbor in atom.GetNeighbors():
                if (neighbor.GetSymbol() == 'O' and
                    neighbor.GetIdx() not in ring_atoms and
                    neighbor.GetDegree() <= 2):  # OH or OR groups

                    # Check if this oxygen is truly exocyclic (not in any ring)
                    neighbor_ring_info = mol.GetRingInfo()
                    neighbor_in_any_ring = False
                    for ring_tuple in neighbor_ring_info.AtomRings():
                        if neighbor.GetIdx() in ring_tuple:
                            neighbor_in_any_ring = True
                            break

                    if not neighbor_in_any_ring:
                        exocyclic_count += 1

        return exocyclic_count
    except Exception:
        return 0


def _find_glycosidic_bridge_oxygens(mol, sugar_rings: List[Tuple[int, ...]]) -> List[int]:
    """
    Find glycosidic bridge oxygens that connect sugar rings to non-sugar parts.

    Criteria:
    - Oxygen with exactly 2 bonds (degree 2)
    - Connected to atoms from different rings (one sugar, one non-sugar)
    """
    bridge_oxygens = []

    try:
        sugar_atoms = set()
        for ring in sugar_rings:
            sugar_atoms.update(ring)

        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'O':
                continue

            if atom.GetDegree() != 2:
                continue

            neighbors = list(atom.GetNeighbors())
            if len(neighbors) != 2:
                continue

            neighbor_idxs = [n.GetIdx() for n in neighbors]

            # Check if one neighbor is in sugar ring and other is not
            in_sugar = [idx in sugar_atoms for idx in neighbor_idxs]

            if sum(in_sugar) == 1:  # Exactly one neighbor in sugar
                bridge_oxygens.append(atom.GetIdx())
                LOG.debug(f"Glycosidic bridge oxygen found: {atom.GetIdx()} "
                           f"connecting {neighbor_idxs}")

    except Exception as e:
        LOG.debug(f"Bridge oxygen detection failed: {e}")

    return bridge_oxygens


def _find_independent_bridge_oxygens(mol) -> List[int]:
    """
    Find potential glycosidic bridge oxygens using topology-based heuristics.

    This function doesn't rely on sugar ring context and uses structural patterns
    that are typical of glycosidic bonds:
    - Oxygen with exactly 2 bonds (degree 2)
    - Connected to sp3 carbons (typical of sugar linkages)
    - Part of acetal/ketal-like patterns (C-O-C where carbons have multiple oxygens)
    """
    bridge_oxygens = []

    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'O':
                continue

            if atom.GetDegree() != 2:
                continue

            # Exclude aromatic oxygens (part of aromatic rings like coumarin, chromone)
            if atom.GetIsAromatic():
                continue

            # Exclude ring oxygens that are not part of confirmed sugar rings
            # This prevents false positives from simple ether rings like THP
            if atom.IsInRing():
                continue

            neighbors = list(atom.GetNeighbors())
            if len(neighbors) != 2:
                continue

            # Check if both neighbors are carbons first (this filtering only applies to C-O-C bridges)
            if not all(n.GetSymbol() == 'C' for n in neighbors):
                continue

            # Exclude lactone/ester oxygens: if O is connected to C* that has C*=O double bond to another O
            # This prevents false positives from gamma-butyrolactone and similar non-aromatic lactones/esters
            # Note: This filtering only applies to C-O-C bridges as confirmed above
            is_lactone_ester_oxygen = False
            for carbon_neighbor in neighbors:
                # Check if this carbon has a double bond to another oxygen (C=O pattern)
                for bond in carbon_neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(carbon_neighbor)
                    if (other_atom.GetSymbol() == 'O' and
                        other_atom.GetIdx() != atom.GetIdx() and
                        bond.GetBondType() == Chem.BondType.DOUBLE):
                        is_lactone_ester_oxygen = True
                        break

                if is_lactone_ester_oxygen:
                    break

            if is_lactone_ester_oxygen:
                continue

            # Heuristic: look for acetal/ketal patterns where at least one carbon
            # has multiple oxygen neighbors (typical of sugar linkages)
            has_sugar_like_carbon = False
            for carbon in neighbors:
                oxygen_neighbors = [n for n in carbon.GetNeighbors() if n.GetSymbol() == 'O']
                # If carbon has 2+ oxygen neighbors, it's likely part of sugar structure
                if len(oxygen_neighbors) >= 2:
                    has_sugar_like_carbon = True
                    break

            if has_sugar_like_carbon:
                bridge_oxygens.append(atom.GetIdx())
                LOG.debug(f"Independent bridge oxygen found: {atom.GetIdx()} "
                           f"(acetal/ketal pattern)")

    except Exception as e:
        LOG.debug(f"Independent bridge oxygen detection failed: {e}")

    return bridge_oxygens


def _find_exocyclic_oxygens(mol, sugar_rings: List[Tuple[int, ...]]) -> List[int]:
    """Find exocyclic oxygen atoms attached to sugar rings."""
    exocyclic_oxygens = []

    try:
        ring_atoms = set()
        for ring in sugar_rings:
            ring_atoms.update(ring)

        for ring in sugar_rings:
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)

                # Find oxygen neighbors that are exocyclic
                for neighbor in atom.GetNeighbors():
                    if (neighbor.GetSymbol() == 'O' and
                        neighbor.GetIdx() not in ring_atoms and
                        neighbor.GetDegree() <= 2):  # OH or OR groups

                        # Check if this oxygen is truly exocyclic (not in any ring)
                        neighbor_ring_info = mol.GetRingInfo()
                        neighbor_in_any_ring = False
                        for ring_tuple in neighbor_ring_info.AtomRings():
                            if neighbor.GetIdx() in ring_tuple:
                                neighbor_in_any_ring = True
                                break

                        if not neighbor_in_any_ring:
                            exocyclic_oxygens.append(neighbor.GetIdx())

    except Exception as e:
        LOG.debug(f"Exocyclic oxygen detection failed: {e}")

    return exocyclic_oxygens


def _has_c_glycoside_pattern(mol, ring: Tuple[int, ...]) -> bool:
    """
    Detect C-glycoside pattern by identifying anomeric carbon connections.

    Looks for ring carbons (potential anomeric C) that are:
    1. Adjacent to the ring oxygen
    2. Connected to sp2/aromatic carbons (aglycon)
    3. Have sugar-like substituents
    """
    try:
        ring_atoms = set(ring)

        # Find the ring oxygen
        ring_oxygen_idx = None
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                ring_oxygen_idx = atom_idx
                break

        if ring_oxygen_idx is None:
            return False

        ring_oxygen = mol.GetAtomWithIdx(ring_oxygen_idx)

        # Check carbons adjacent to ring oxygen (potential anomeric carbons)
        for neighbor in ring_oxygen.GetNeighbors():
            if (neighbor.GetSymbol() == 'C' and
                neighbor.GetIdx() in ring_atoms):

                anomeric_c_idx = neighbor.GetIdx()
                anomeric_c = neighbor

                # Check if this anomeric carbon has C-glycoside linkage
                for anomeric_neighbor in anomeric_c.GetNeighbors():
                    if (anomeric_neighbor.GetSymbol() == 'C' and
                        anomeric_neighbor.GetIdx() not in ring_atoms):

                        # Check if connected to sp2/aromatic carbon (aglycon)
                        # Note: SP2AROM not available in all RDKit versions, use aromatic check instead
                        if (anomeric_neighbor.GetHybridization() == Chem.HybridizationType.SP2 or
                            anomeric_neighbor.GetIsAromatic()):

                            LOG.debug(f"C-glycoside pattern detected: anomeric C {anomeric_c_idx} -> sp2/aromatic C {anomeric_neighbor.GetIdx()}")
                            return True

                        # Also check for connection to aromatic rings
                        for arom_neighbor in anomeric_neighbor.GetNeighbors():
                            if (arom_neighbor.GetSymbol() == 'C' and
                                arom_neighbor.GetIsAromatic() and
                                arom_neighbor.GetIdx() not in ring_atoms):

                                LOG.debug(f"C-glycoside pattern detected: anomeric C {anomeric_c_idx} -> aromatic system via {anomeric_neighbor.GetIdx()}")
                                return True

        return False

    except Exception as e:
        LOG.debug(f"C-glycoside pattern detection failed: {e}")
        return False


def _find_anomeric_carbons_for_masking(mol, sugar_rings: List[Tuple[int, ...]]) -> List[int]:
    """
    Find anomeric carbons in C-glycosides for enhanced masking.

    Returns anomeric carbon indices that should be included in mask
    even if traditional bridge oxygen detection fails.
    """
    anomeric_carbons = []

    try:
        for ring in sugar_rings:
            ring_atoms = set(ring)

            # Find ring oxygen
            ring_oxygen_idx = None
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O':
                    ring_oxygen_idx = atom_idx
                    break

            if ring_oxygen_idx is None:
                continue

            ring_oxygen = mol.GetAtomWithIdx(ring_oxygen_idx)

            # Find anomeric carbons (adjacent to ring oxygen)
            for neighbor in ring_oxygen.GetNeighbors():
                if (neighbor.GetSymbol() == 'C' and
                    neighbor.GetIdx() in ring_atoms):

                    anomeric_c_idx = neighbor.GetIdx()

                    # Check for C-glycoside connection
                    has_c_linkage = False
                    for anomeric_neighbor in neighbor.GetNeighbors():
                        if (anomeric_neighbor.GetSymbol() == 'C' and
                            anomeric_neighbor.GetIdx() not in ring_atoms):
                            has_c_linkage = True
                            break

                    if has_c_linkage:
                        anomeric_carbons.append(anomeric_c_idx)
                        LOG.debug(f"Anomeric carbon for masking: {anomeric_c_idx}")

    except Exception as e:
        LOG.debug(f"Anomeric carbon detection failed: {e}")

    return anomeric_carbons


def _detect_c_glycoside_topo(mol, sugar_cfg: Dict[str, Any] = None) -> Tuple[List[int], Dict[str, int]]:
    """
    Detect C-glycoside topology patterns with strict constraints.

    Applies hard constraints to prevent false detection of chromone/flavone cores:
    1. Ring O must be non-aromatic
    2. Anomeric C must be in ring and SP3 hybridized
    3. O-C bond must be single bond (excludes lactone/carbonyl)
    4. Ring must be 5/6-member, non-aromatic, no intracyclic carbonyl
    5. External sp2/aromatic C must be outside the ring

    Returns only minimal evidence subgraph (anomeric C + neighbors), not whole rings.

    Returns:
        tuple: (List of atom indices, audit statistics dict)
        audit_stats contains:
        - 'direct_cglyco_hits': Number of direct C-glycoside connections found
        - 'one_hop_cglyco_hits': Number of one-hop C-glycoside connections found
        - 'intermediate_included_count': Number of intermediate carbons included in mask
        - 'total_anomeric_candidates': Total anomeric carbon candidates examined
    """
    # Merge user config with defaults
    if sugar_cfg is None:
        sugar_cfg = _get_default_sugar_cfg()

    c_glycoside_atoms = []

    # Initialize audit statistics
    audit_stats = {
        'direct_cglyco_hits': 0,
        'one_hop_cglyco_hits': 0,
        'intermediate_included_count': 0,
        'total_anomeric_candidates': 0
    }

    try:
        ring_info = mol.GetRingInfo()

        # Find all qualifying rings first (5/6-member, non-aromatic, exactly one O)
        qualifying_rings = []
        for ring in ring_info.AtomRings():
            if len(ring) not in (5, 6):
                continue

            # Check if ring is non-aromatic
            is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            if is_aromatic:
                continue

            # Count oxygens and check for intracyclic carbonyls
            oxygen_count = 0
            has_intracyclic_carbonyl = False
            ring_oxygen_idx = None

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O':
                    oxygen_count += 1
                    ring_oxygen_idx = atom_idx

                    # Check if this oxygen is part of intracyclic carbonyl
                    for bond in atom.GetBonds():
                        other_atom = bond.GetOtherAtom(atom)
                        if (other_atom.GetIdx() in ring and
                            other_atom.GetSymbol() == 'C' and
                            bond.GetBondType() == Chem.BondType.DOUBLE):
                            has_intracyclic_carbonyl = True
                            break

            if oxygen_count == 1 and not has_intracyclic_carbonyl and ring_oxygen_idx is not None:
                qualifying_rings.append((ring, ring_oxygen_idx))

        # For each qualifying ring, check for C-glycoside patterns
        for ring, ring_oxygen_idx in qualifying_rings:
            ring_atoms = set(ring)
            ring_oxygen = mol.GetAtomWithIdx(ring_oxygen_idx)

            # Check if ring oxygen is non-aromatic
            if ring_oxygen.GetIsAromatic():
                continue

            # Find carbon neighbors of ring oxygen (potential anomeric carbons)
            for neighbor in ring_oxygen.GetNeighbors():
                if neighbor.GetSymbol() != 'C' or neighbor.GetIdx() not in ring_atoms:
                    continue

                anomeric_c = neighbor
                anomeric_idx = anomeric_c.GetIdx()

                # Audit: Count anomeric carbon candidates
                audit_stats['total_anomeric_candidates'] += 1

                # Hard constraint 1: Anomeric C must be SP3 (with redundant criteria for RDKit robustness)
                is_sp3_like = (
                    anomeric_c.GetHybridization() == Chem.HybridizationType.SP3 or
                    (anomeric_c.GetDegree() == 4 and
                     anomeric_c.GetTotalNumHs() in (0, 1) and
                     not any(bond.GetBondType() != Chem.BondType.SINGLE for bond in anomeric_c.GetBonds()))
                )
                if not is_sp3_like:
                    continue

                # Hard constraint 2: O-C bond must be single bond
                o_c_bond = mol.GetBondBetweenAtoms(ring_oxygen_idx, anomeric_idx)
                if o_c_bond is None or o_c_bond.GetBondType() != Chem.BondType.SINGLE:
                    continue

                # Check for C-glycoside linkage to external sp2/aromatic carbon (direct or one-hop)
                has_c_glycoside_linkage = False
                intermediate_carbon_idx = None  # Track intermediate carbon for potential inclusion

                for anomeric_neighbor in anomeric_c.GetNeighbors():
                    if (anomeric_neighbor.GetSymbol() == 'C' and
                        anomeric_neighbor.GetIdx() not in ring_atoms):  # Must be external

                        # Direct connection: anomeric_C -> sp2/aromatic_C
                        # Note: SP2AROM not available in all RDKit versions, use aromatic check instead
                        if (anomeric_neighbor.GetHybridization() == Chem.HybridizationType.SP2 or
                            anomeric_neighbor.GetIsAromatic()):

                            has_c_glycoside_linkage = True
                            # Audit: Count direct C-glycoside hits
                            audit_stats['direct_cglyco_hits'] += 1
                            LOG.debug(f"C-glycoside topology detected (direct): ring O {ring_oxygen_idx} -> anomeric C {anomeric_idx} -> external sp2/aromatic C {anomeric_neighbor.GetIdx()}")
                            break

                        # One-hop connection: anomeric_C -> intermediate_C(SP3/external/single) -> sp2/aromatic_C
                        # Only if enabled in configuration
                        elif (sugar_cfg.get('allow_one_hop_cglyco', True) and
                              anomeric_neighbor.GetHybridization() == Chem.HybridizationType.SP3 and
                              mol.GetBondBetweenAtoms(anomeric_idx, anomeric_neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE):

                            intermediate_c_idx = anomeric_neighbor.GetIdx()

                            for intermediate_neighbor in anomeric_neighbor.GetNeighbors():
                                if (intermediate_neighbor.GetSymbol() == 'C' and
                                    intermediate_neighbor.GetIdx() not in ring_atoms and
                                    intermediate_neighbor.GetIdx() != anomeric_idx):

                                    # Note: SP2AROM not available in all RDKit versions, use aromatic check instead
                                    if (intermediate_neighbor.GetHybridization() == Chem.HybridizationType.SP2 or
                                        intermediate_neighbor.GetIsAromatic()):

                                        has_c_glycoside_linkage = True
                                        # Store intermediate carbon for potential inclusion in mask
                                        intermediate_carbon_idx = intermediate_c_idx
                                        # Audit: Count one-hop C-glycoside hits
                                        audit_stats['one_hop_cglyco_hits'] += 1
                                        LOG.debug(f"C-glycoside topology detected (one-hop): ring O {ring_oxygen_idx} -> anomeric C {anomeric_idx} -> intermediate C {intermediate_c_idx} -> external sp2/aromatic C {intermediate_neighbor.GetIdx()}")
                                        break

                            if has_c_glycoside_linkage:
                                break

                if has_c_glycoside_linkage:
                    # Add only minimal evidence subgraph (anomeric C + immediate neighbors in ring)
                    c_glycoside_atoms.append(anomeric_idx)

                    # Add neighboring ring atoms to anomeric carbon
                    for ring_neighbor in anomeric_c.GetNeighbors():
                        if (ring_neighbor.GetIdx() in ring_atoms and
                            ring_neighbor.GetSymbol() in ['C', 'O']):
                            c_glycoside_atoms.append(ring_neighbor.GetIdx())

                    # P0-H2: Optionally include intermediate carbon in minimal mask
                    if (intermediate_carbon_idx is not None and
                        sugar_cfg.get('cglyco_minimal_mask_include_intermediate', False)):
                        c_glycoside_atoms.append(intermediate_carbon_idx)
                        # Audit: Count intermediate carbons included in mask
                        audit_stats['intermediate_included_count'] += 1
                        LOG.debug(f"C-glycoside minimal evidence subgraph: anomeric C {anomeric_idx} + ring neighbors + intermediate C {intermediate_carbon_idx}")
                    else:
                        LOG.debug(f"C-glycoside minimal evidence subgraph: anomeric C {anomeric_idx} + ring neighbors")

    except Exception as e:
        LOG.debug(f"C-glycoside topology detection failed: {e}")

    # Remove duplicates and return
    unique_atoms = list(set(c_glycoside_atoms))
    LOG.debug(f"C-glycoside topology detection complete: {len(unique_atoms)} atoms identified")
    return unique_atoms, audit_stats


def _perform_degraded_masking(mol, sugar_cfg: Dict[str, Any]) -> Tuple[Set[int], bool, List[str]]:
    """
    Perform degraded masking when sugar ring detection fails.

    This function is triggered when no sugar rings are found, indicating
    either a detection failure or C-glycoside patterns that need alternative
    approaches.

    Returns:
        Tuple of (minimal_mask_atoms, degraded_flag, degraded_reasons)
        degraded_flag is True only when evidence is found AND atoms are masked
    """
    minimal_mask = set()
    reasons = []

    try:
        # Try independent bridge oxygen detection
        try:
            independent_bridge_oxygens = _find_independent_bridge_oxygens(mol)
            if independent_bridge_oxygens:
                minimal_mask.update(independent_bridge_oxygens)
                reasons.append('bridge_O')
                LOG.debug(f"Degraded masking: found {len(independent_bridge_oxygens)} independent bridge oxygens")
        except Exception as e:
            LOG.debug(f"Independent bridge oxygen detection failed: {e}")

        # Try C-glycoside topology detection (enhanced with parameterized constraints)
        try:
            c_glycoside_atoms, _ = _detect_c_glycoside_topo(mol, sugar_cfg)
            if c_glycoside_atoms:
                minimal_mask.update(c_glycoside_atoms)
                reasons.append('cglyco_topo')
                LOG.debug(f"Degraded masking: found {len(c_glycoside_atoms)} C-glycoside atoms")
        except Exception as e:
            LOG.debug(f"C-glycoside topology detection failed: {e}")

        # Try to find sugar-like oxygen rings as fallback (only if C-glycoside topology didn't find evidence)
        try:
            # Only require fallback rings if we don't have C-glycoside topology evidence
            if 'cglyco_topo' not in reasons:
                fallback_rings = _find_oxygen_containing_rings(mol, sugar_cfg)
                if fallback_rings:
                    for ring in fallback_rings:
                        minimal_mask.update(ring)
                    reasons.append('oxygen_ring_fallback')
                    LOG.debug(f"Degraded masking: found {len(fallback_rings)} oxygen-containing rings as fallback")
            else:
                LOG.debug("Skipping oxygen ring fallback - C-glycoside topology already provided evidence")
        except Exception as e:
            LOG.debug(f"Fallback ring detection failed: {e}")

    except Exception as e:
        LOG.debug(f"Degraded masking failed completely: {e}")

    # Only mark as degraded if we found evidence AND have atoms to mask
    degraded = bool(reasons) and bool(minimal_mask)
    LOG.debug(f"Degraded masking complete: {len(minimal_mask)} atoms masked, degraded={degraded}, reasons={reasons}")

    return minimal_mask, degraded, reasons


def _find_sugar_like_oxygen_rings_fallback(mol, sugar_cfg: Dict[str, Any] = None) -> List[Tuple[int, ...]]:
    """
    Find sugar-like oxygen rings with strict constraints for degraded masking fallback.

    Applies strict constraints to exclude chromone/flavone cores:
    1. 5/6-member, non-aromatic rings
    2. Exactly 1 ring oxygen
    3. No intracyclic carbonyl (excludes lactone/aromatic oxygen rings)
    4. Ring not all sp2 (>=50% sp3 atoms)
    5. External oxygen count >= 1

    Returns:
        List of qualifying ring tuples
    """
    # Merge user config with defaults
    if sugar_cfg is None:
        sugar_cfg = _get_default_sugar_cfg()

    sugar_like_rings = []

    try:
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        for ring in atom_rings:
            if len(ring) not in (5, 6):
                continue

            # Constraint 1: Ring must be non-aromatic
            is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            if is_aromatic:
                continue

            # Constraint 2: Count oxygens and check for intracyclic carbonyls
            oxygen_count = 0
            has_intracyclic_carbonyl = False
            ring_oxygen_idx = None

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'O':
                    oxygen_count += 1
                    ring_oxygen_idx = atom_idx

                    # Check if this oxygen is part of intracyclic carbonyl
                    for bond in atom.GetBonds():
                        other_atom = bond.GetOtherAtom(atom)
                        if (other_atom.GetIdx() in ring and
                            other_atom.GetSymbol() == 'C' and
                            bond.GetBondType() == Chem.BondType.DOUBLE):
                            has_intracyclic_carbonyl = True
                            break

            # Must have exactly 1 oxygen and no intracyclic carbonyl
            if oxygen_count != 1 or has_intracyclic_carbonyl:
                continue

            # Constraint 3: Ring should have adequate sp3 ratio (not all sp2)
            # Use parameterized thresholds from configuration
            sp3_count = 0
            total_count = len(ring)
            sp3_threshold = (sugar_cfg['sp3_threshold_5_ring'] if len(ring) == 5
                           else sugar_cfg['sp3_threshold_6_ring'])

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetHybridization() == Chem.HybridizationType.SP3:
                    sp3_count += 1

            if sp3_count < total_count * sp3_threshold:
                continue

            # Constraint 4: External oxygen count >= 1 (single bond, deduplicated)
            external_oxygens = set()
            ring_atoms = set(ring)

            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if (neighbor.GetSymbol() == 'O' and
                        neighbor.GetIdx() not in ring_atoms):

                        # Check if it's a single bond (exclude NO2, SO2, etc.)
                        bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                            external_oxygens.add(neighbor.GetIdx())

            external_oxygen_count = len(external_oxygens)
            if external_oxygen_count < 1:
                continue

            # Apply stronger evidence requirements to prevent false positives
            # Require either: >=2 exocyclic oxygens OR C-glycoside evidence
            has_cglyco_evidence = _has_c_glycoside_pattern(mol, ring)
            has_strong_evidence = (external_oxygen_count >= 2) or has_cglyco_evidence

            if not has_strong_evidence:
                LOG.debug(f"Sugar-like ring rejected due to weak evidence: {ring} (external O: {external_oxygen_count}, cglyco: {has_cglyco_evidence})")
                continue

            # All constraints satisfied - this is a sugar-like ring
            sugar_like_rings.append(ring)
            LOG.debug(f"Sugar-like oxygen ring found: {ring} (sp3: {sp3_count}/{total_count}, external O: {external_oxygen_count}, cglyco: {has_cglyco_evidence})")

    except Exception as e:
        LOG.debug(f"Sugar-like oxygen ring fallback failed: {e}")

    return sugar_like_rings


def _find_oxygen_containing_rings(mol, sugar_cfg: Dict[str, Any] = None) -> List[Tuple[int, ...]]:
    """
    Legacy function - now redirects to strict sugar-like detection.
    """
    return _find_sugar_like_oxygen_rings_fallback(mol, sugar_cfg)


def _get_sugar_mask_sru(mol, sugar_cfg: Dict[str, Any]) -> Set[int]:
    """
    SRU-based sugar masking (placeholder for future implementation).
    """
    # TODO: Implement SRU client in future PR
    raise NotImplementedError("SRU sugar masking not yet implemented")
