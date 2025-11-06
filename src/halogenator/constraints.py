# -*- coding: ascii -*-
"""Constraint engine for k-dimensional enumeration."""

from typing import Dict, Any, List, Tuple, Optional
from collections import Counter


class HistoryItem(dict):
    """
    History item representing one substitution step.
    
    Expected fields:
    - rule: str (e.g., 'R1', 'R2', 'R3')
    - site: int (atom index, None for reaction rules)
    - sym: int (symmetry class)
    - ring_tag: tuple (stable ring identifier)
    - halogen: str (e.g., 'F', 'Cl', 'Br', 'I')
    - depth: int (substitution level)
    """
    pass


def accept(mol, history: List[Dict[str, Any]], cfg: Dict[str, Any]):
    """
    Check if current substitution pattern satisfies all constraints.

    Args:
        mol: Current molecule
        history: List of substitution steps
        cfg: Constraint configuration dict (supports 'enable' master switch)

    Returns:
        (ok: bool, violations: dict)
    """
    # Master switch: if constraints disabled, accept all products (raw mode)
    if not cfg.get('enable', True):
        return (True, {})

    violations = {}
    per_ring_quota = cfg.get('per_ring_quota', None)
    max_per_halogen = cfg.get('max_per_halogen', None)
    d_min = cfg.get('min_graph_distance', None)

    # per_ring_quota - only count non-empty ring tags (flavonoid rings: A/B/C)
    if per_ring_quota is not None:
        cnt = Counter([h.get('ring_tag') for h in history if h.get('ring_tag') and h.get('ring_tag') != ''])
        over = [str(tag) for tag, n in cnt.items() if n > per_ring_quota]
        if over:
            violations['per_ring_quota'] = over

    # max_per_halogen
    if isinstance(max_per_halogen, dict):
        cnt_h = Counter([h.get('halogen') for h in history if h.get('halogen')])
        over_h = [hx for hx, n in cnt_h.items() if max_per_halogen.get(hx) is not None and n > max_per_halogen[hx]]
        if over_h:
            violations['max_per_halogen'] = over_h

    # min_graph_distance (only check site-based rules)
    if d_min and d_min > 0:
        try:
            from .chem_compat import Chem
        except Exception:
            # Skip distance check if RDKit not available
            pass
        else:
            # Collect sites from history (site-based rules only)
            sites = [int(h['site']) for h in history if 'site' in h and isinstance(h['site'], int)]
            if len(sites) >= 2:
                try:
                    D = Chem.GetDistanceMatrix(mol)
                    # Check all pairs for minimum distance
                    bad_pairs = []
                    for i in range(len(sites)):
                        for j in range(i+1, len(sites)):
                            if int(D[sites[i], sites[j]]) < d_min:
                                bad_pairs.append((sites[i], sites[j]))
                    if bad_pairs:
                        violations['min_graph_distance'] = bad_pairs
                except Exception:
                    pass  # Allow if distance calculation fails

    return (len(violations) == 0, violations)


def _check_per_ring_quota(history: List[HistoryItem], per_ring_quota: int) -> Optional[str]:
    """Check per-ring substitution quota."""
    if per_ring_quota <= 0:
        return None
    
    # Count substitutions per ring_tag
    ring_counts = Counter()
    
    for item in history:
        ring_tag = item.get('ring_tag')
        if ring_tag is not None:
            # Convert to hashable representation if needed
            if isinstance(ring_tag, (tuple, list)):
                ring_key = tuple(ring_tag) if isinstance(ring_tag, list) else ring_tag
            else:
                ring_key = ring_tag
            
            ring_counts[ring_key] += 1
    
    # Check if any ring exceeds quota
    for ring_key, count in ring_counts.items():
        if count > per_ring_quota:
            return f"Ring {ring_key} has {count} substitutions (limit: {per_ring_quota})"
    
    return None


def _check_min_distance(mol, history: List[HistoryItem], min_distance: int) -> Optional[str]:
    """Check minimum graph distance between substitution sites."""
    if min_distance <= 0:
        return None
    
    try:
        from .chem_compat import Chem
    except Exception:
        # Skip if RDKit not available
        return None
    
    # Collect site indices from history (skip reaction-based rules)
    sites = []
    for item in history:
        site = item.get('site')
        if site is not None:
            sites.append(int(site))
    
    if len(sites) <= 1:
        return None  # Need at least 2 sites to check distance
    
    try:
        from .chem_compat import Chem
        # Compute distance matrix
        distance_matrix = Chem.GetDistanceMatrix(mol)
        
        # Check all pairs of sites
        for i in range(len(sites)):
            for j in range(i + 1, len(sites)):
                site1, site2 = sites[i], sites[j]
                
                # Ensure indices are within bounds
                if site1 >= len(distance_matrix) or site2 >= len(distance_matrix[0]):
                    continue
                
                actual_distance = int(distance_matrix[site1][site2])
                if actual_distance < min_distance:
                    return f"Sites {site1} and {site2} are {actual_distance} bonds apart (min: {min_distance})"
    
    except Exception as e:
        # Fallback: allow if distance calculation fails
        return None
    
    return None


def _check_max_per_halogen(history: List[HistoryItem], max_per_halogen: Dict[str, int]) -> Optional[str]:
    """Check maximum count per halogen type."""
    if not max_per_halogen:
        return None
    
    # Count halogens in history
    halogen_counts = Counter()
    for item in history:
        halogen = item.get('halogen')
        if halogen:
            halogen_counts[halogen] += 1
    
    # Check limits
    for halogen, count in halogen_counts.items():
        limit = max_per_halogen.get(halogen)
        if limit is not None and count > limit:
            return f"Halogen {halogen} count {count} exceeds limit {limit}"
    
    return None


def default_constraints() -> Dict[str, Any]:
    """Get default constraint configuration."""
    return {
        'enable': True,  # Master switch (set to False for raw mode)
        'per_ring_quota': 2,
        'min_graph_distance': 2,
        'max_per_halogen': None  # No limits by default
    }