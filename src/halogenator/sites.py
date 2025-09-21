# -*- coding: ascii -*-
"""Site identification and symmetry folding.

A/B/C Ring Labeling System:
===========================

This module provides flavonoid-specific A/B/C ring labeling for natural product
halogenation enumeration. The ring labeling follows standard flavonoid nomenclature:

- **C Ring**: 6-membered ring containing oxygen and carbonyl (C=O) - the central 
  heterocyclic ring in flavonoids (pyran-4-one core)
- **A Ring**: Aromatic 6-membered ring fused to C ring (shares 2+ atoms) - typically 
  the resorcinol-like ring in natural flavonoids
- **B Ring**: Aromatic 6-membered ring single-bond connected to C ring - typically 
  the phenolic ring at C-2 position

Scope and Usage:
---------------
- Designed specifically for flavonoid molecules (quercetin, apigenin, kaempferol, etc.)
- Used for per-ring quota constraints (per_ring_quota=1 prevents multiple substitutions on same ring)
- Tracks substitution sites by ring type in enumeration history
- Non-flavonoid molecules return empty labels ('')
- Sugar rings in glycosides are excluded from labeling

Performance Optimization:
------------------------
- flavonoid_ring_label() uses WeakKeyDictionary caching based on molecule objects
- Computes labels for all atoms at once and caches per-molecule to avoid index mismatch
- WeakKeyDictionary automatically cleans up when molecules are garbage collected
- Fallback LRU cache using InChI keys for molecules that can't use weak references

Examples:
--------
Quercetin: O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12
- A ring: positions with resorcinol-like aromatic ring (atoms in fused benzene)
- B ring: positions in phenolic ring (atoms in single-bond connected benzene) 
- C ring: positions in pyran-4-one ring (atoms in O-containing heterocycle)

Benzene: c1ccccc1
- Returns '' for all atoms (not a flavonoid structure)
"""

import functools
import weakref
from typing import List, Dict, Any, Optional, Set, Iterable, Tuple
from .chem_compat import Chem, mol_to_smiles_safe, canonical_rank_atoms_safe

# Sugar ring identification cache to avoid raw+masked duplicate computation
try:
    _sugar_rings_cache = weakref.WeakKeyDictionary()
    _SUGAR_RING_CACHE_MODE = "weakkey"
except TypeError:
    _sugar_rings_cache = None
    _SUGAR_RING_CACHE_MODE = "smiles-lru"

# Cache statistics for debugging
_SUGAR_RING_CACHE_STATS = {"mode": _SUGAR_RING_CACHE_MODE, "hit": 0, "miss": 0}


def _cfg_key_tuple(sugar_cfg: dict) -> tuple:
    """Convert sugar_cfg to a hashable tuple for caching."""
    return tuple(sorted((k, str(v)) for k, v in (sugar_cfg or {}).items()))


@functools.lru_cache(maxsize=1024)
def _find_sugar_rings_by_smiles(smiles: str, cfg_key: tuple) -> tuple:
    """SMILES-based LRU cached version of _find_sugar_rings."""
    # Delayed import to avoid circular imports
    from .sugar_mask import _find_sugar_rings
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ()
    cfg_dict = dict(cfg_key) if cfg_key else {}
    return tuple(_find_sugar_rings(mol, cfg_dict))


def _find_sugar_rings_cached(mol, sugar_cfg):
    """
    Cached version of _find_sugar_rings to avoid duplicate computation in raw+masked calls.
    Uses hybrid strategy: WeakKeyDictionary when available, falls back to SMILES-LRU cache.
    """
    cfg_key = _cfg_key_tuple(sugar_cfg)

    if _SUGAR_RING_CACHE_MODE == "weakkey" and _sugar_rings_cache is not None:
        # Try WeakKeyDictionary approach
        try:
            entry = _sugar_rings_cache.get(mol)
            if entry is not None and cfg_key in entry:
                _SUGAR_RING_CACHE_STATS["hit"] += 1
                return entry[cfg_key]

            # Cache miss - compute and store
            from .sugar_mask import _find_sugar_rings
            rings = list(_find_sugar_rings(mol, sugar_cfg))
            entry = {} if entry is None else dict(entry)
            entry[cfg_key] = rings
            _sugar_rings_cache[mol] = entry
            _SUGAR_RING_CACHE_STATS["miss"] += 1
            return rings
        except (TypeError, AttributeError):
            # WeakKeyDictionary failed, fall back to SMILES mode for this call
            pass

    # SMILES-LRU fallback mode
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False) if mol else ""
    if smiles:
        _SUGAR_RING_CACHE_STATS["miss"] += 1
        return list(_find_sugar_rings_by_smiles(smiles, cfg_key))
    else:
        # Direct computation as last resort
        from .sugar_mask import _find_sugar_rings
        _SUGAR_RING_CACHE_STATS["miss"] += 1
        return list(_find_sugar_rings(mol, sugar_cfg))


# T29: Optional cache hit monitoring for debugging
def get_sugar_ring_cache_stats() -> Dict[str, Any]:
    """Get current sugar ring cache statistics for debugging."""
    global _SUGAR_RING_CACHE_STATS
    stats = dict(_SUGAR_RING_CACHE_STATS)

    # Add cache size if available
    if _SUGAR_RING_CACHE_MODE == "weakkey" and _sugar_rings_cache is not None:
        try:
            stats["cache_size"] = len(_sugar_rings_cache)
        except (TypeError, AttributeError):
            stats["cache_size"] = "unknown"
    else:
        stats["cache_size"] = "n/a"

    # Calculate hit rate
    total = stats["hit"] + stats["miss"]
    stats["hit_rate"] = (stats["hit"] / total * 100) if total > 0 else 0.0

    return stats


def reset_sugar_ring_cache_stats():
    """Reset sugar ring cache statistics for debugging purposes."""
    global _SUGAR_RING_CACHE_STATS
    _SUGAR_RING_CACHE_STATS = {"mode": _SUGAR_RING_CACHE_MODE, "hit": 0, "miss": 0}


def log_sugar_ring_cache_stats(logger=None, prefix="Sugar ring cache"):
    """Log current cache statistics if monitoring is enabled."""
    import os
    import logging

    # Only log if explicitly enabled via environment variable
    if not os.environ.get("HALOGENATOR_CACHE_MONITOR", "").lower() in ("1", "true", "yes"):
        return

    stats = get_sugar_ring_cache_stats()
    log_func = logger.info if logger else logging.getLogger(__name__).info

    log_func(f"{prefix}: mode={stats['mode']}, hits={stats['hit']}, "
             f"misses={stats['miss']}, hit_rate={stats['hit_rate']:.1f}%, "
             f"cache_size={stats['cache_size']}")


def filter_sites_with_mask(sites: List[int], mask: Set[int]) -> List[int]:
    """
    Filter site indices by removing any that intersect with the given mask.

    This provides a unified abstraction for site filtering across all rules (R1, R2, R6).

    Args:
        sites: List of site indices to filter
        mask: Set of atom indices that should be excluded

    Returns:
        List of site indices with masked atoms removed
    """
    return [site for site in sites if site not in mask]


def ensure_ready(mol) -> None:
    """Ensure molecule is ready for ring info and property access."""
    # Check if full RDKit is available by testing Chem functionality
    rdkit_available = hasattr(Chem, 'SanitizeMol') and hasattr(Chem, 'GetSymmSSSR')
    
    if rdkit_available:
        try:
            # Use basic sanitization without specific flags to avoid direct rdkit import
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass
        
        # Force ring info initialization
        try:
            Chem.GetSymmSSSR(mol)
        except Exception:
            pass
    else:
        # Without RDKit, still try basic property access via chem_compat fallbacks
        try:
            # Build canonical ranks using chem_compat fallbacks to ensure atom ordering available
            from .chem_compat import canonical_rank_atoms_safe
            canonical_rank_atoms_safe(mol, breakTies=False)
        except Exception:
            pass


def is_carbonyl_carbon(atom) -> Optional[bool]:
    """Check if atom is carbonyl carbon (has O double bond neighbor).
    
    Returns:
        True if carbonyl carbon, False if not, None if determination not possible (no RDKit)
    """    
    if atom is None:
        return False
    
    try:
        if atom.GetSymbol() != 'C':
            return False
    except AttributeError:
        return False
    
    # Guard RDKit-only branches
    if not hasattr(Chem, 'BondType'):
        return None  # explicit not-available
    
    try:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                bond = atom.GetOwningMol().GetBondBetweenAtoms(
                    atom.GetIdx(), neighbor.GetIdx()
                )
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    return True
    except Exception:
        return None  # Error accessing bond information
    return False


def is_c_ring_site_ready(mol, aidx: int) -> bool:
    """Check if atom is valid C ring site for R2. Assumes mol is already prepared."""
    atom = mol.GetAtomWithIdx(aidx)
    
    # Must be carbon
    if atom.GetSymbol() != 'C':
        return False
    
    # Must be in 6-member ring
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    in_6ring = False
    six_rings = []
    
    for ring in atom_rings:
        if aidx in ring and len(ring) == 6:
            in_6ring = True
            six_rings.append(ring)
    
    if not in_6ring:
        return False
    
    # Must have exactly 1 hydrogen
    try:
        if atom.GetTotalNumHs() != 1:
            return False
    except Exception:
        # Skip if hydrogen count unavailable
        return False
    
    # Must not be carbonyl carbon
    carbonyl_result = is_carbonyl_carbon(atom)
    if carbonyl_result is True:  # Only reject if definitively carbonyl
        return False
    # If None (can't determine) or False, continue processing
    # NOTE: When carbonyl_result is None (RDKit unavailable), we assume non-carbonyl.
    # This is a "relaxed" policy that may include some false positives in site selection
    # but avoids blocking enumeration when RDKit is unavailable.
    
    # Must have O neighbor in the same 6-member ring
    has_ring_oxygen = False
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'O':
            neighbor_idx = neighbor.GetIdx()
            # Check if this O is in any of the same 6-rings as our carbon
            for ring in six_rings:
                if neighbor_idx in ring:
                    has_ring_oxygen = True
                    break
            if has_ring_oxygen:
                break
    
    return has_ring_oxygen


def is_c_ring_site(mol, aidx: int) -> bool:
    """Check if atom is valid C ring site for R2."""
    ensure_ready(mol)
    return is_c_ring_site_ready(mol, aidx)


def aromatic_CH_indices(mol, masked_atoms: set = None) -> List[int]:
    """Get aromatic carbon indices with H > 0."""
    ensure_ready(mol)

    if masked_atoms is None:
        masked_atoms = set()

    indices = []
    for atom in mol.GetAtoms():
        try:
            idx = atom.GetIdx()
            if idx in masked_atoms:
                continue
            if (atom.GetSymbol() == 'C' and
                atom.GetIsAromatic() and
                atom.GetTotalNumHs() > 0):
                indices.append(idx)
        except Exception:
            # Skip atoms that can't provide ring/hydrogen info
            continue
    return indices


def symmetry_groups(mol, indices: List[int], extra_tag: Optional[str] = None) -> Dict[Any, List[int]]:
    """Group indices by symmetry using CanonicalRankAtoms."""
    if not indices:
        return {}
    
    # Ensure molecule is ready for symmetry operations
    ensure_ready(mol)
    
    # Get canonical ranks for symmetry folding using chem_compat
    try:
        # Chem is imported at module level from chem_compat
        pass
    except Exception:
        return {i: [idx] for i, idx in enumerate(indices)}
    
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    
    groups = {}
    for idx in indices:
        rank = ranks[idx]
        
        if extra_tag is not None:
            # For R2, we might want to distinguish by ring membership
            key = (extra_tag, rank)
        else:
            # For R1, just use rank
            key = rank
        
        if key not in groups:
            groups[key] = []
        groups[key].append(idx)
    
    return groups


def c_ring_indices(mol) -> List[int]:
    """Get C ring site indices using is_c_ring_site_ready."""
    ensure_ready(mol)
    
    indices = []
    for i in range(mol.GetNumAtoms()):
        try:
            if is_c_ring_site_ready(mol, i):
                indices.append(i)
        except Exception:
            # Skip atoms that cause ring info errors
            continue
    return indices


def get_ring_tag_for_atom(mol, aidx: int) -> Optional[int]:
    """Get a stable ring tag for atoms in 6-member rings."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    for ring in atom_rings:
        if aidx in ring and len(ring) == 6:
            # Return minimum atom index in ring as stable tag
            return min(ring)
    
    return None


def _is_carbonyl_in_ring(mol, ring_set) -> bool:
    """Check if ring contains a carbonyl C=O bond (C must be in ring, O can be outside)."""
    ring_atoms = set(ring_set)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 2.0:
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            syma, symb = a.GetSymbol(), b.GetSymbol()
            ia, ib = a.GetIdx(), b.GetIdx()
            if ((syma == 'C' and symb == 'O') or (syma == 'O' and symb == 'C')):
                # At least the carbon must be in the ring (oxygen can be outside for lactones/pyranones)
                carbon_idx = ia if syma == 'C' else ib
                if carbon_idx in ring_atoms:
                    return True
    return False


def _ring_is_aromatic_6(mol, ring) -> bool:
    """Check if 6-member ring is fully aromatic carbon."""
    if len(ring) != 6:
        return False
    for i in ring:
        a = mol.GetAtomWithIdx(i)
        if a.GetSymbol() != 'C' or not a.GetIsAromatic():
            return False
    return True


def _shared_atoms(r1, r2) -> int:
    """Count shared atoms between two rings."""
    return len(set(r1) & set(r2))


def _is_b_connected(mol, c_ring, candidate) -> bool:
    """Check if candidate ring is single-bond connected to C ring."""
    ca = set(c_ring)
    cb = set(candidate)
    if ca & cb:  # Share atoms - not single-bond connected
        return False
    edges = 0
    for i in ca:
        ai = mol.GetAtomWithIdx(i)
        for nb in ai.GetNeighbors():
            j = nb.GetIdx()
            if j in cb:
                edges += 1
                if edges > 1:
                    break
        if edges > 1:
            break
    return edges == 1


# Weak reference cache for per-molecule ring labeling with signature-based invalidation
# Maps mol -> (signature, labels) where signature is used for topology change detection
_mol_ring_labels_cache = weakref.WeakKeyDictionary()

# Cache statistics tracking
_cache_stats = {'hits': 0, 'misses': 0, 'invalidations': 0}


def _compute_molecule_signature(mol) -> tuple:
    """
    Compute a lightweight signature for topology change detection.
    Uses non-isomeric canonical SMILES to ignore isotopes but detect topology/aromaticity changes.
    Returns (num_atoms, smiles) tuple.
    """
    try:
        num_atoms = mol.GetNumAtoms()
        smiles = mol_to_smiles_safe(mol, canonical=True, isomericSmiles=False)
        # Deterministic integer signature from topology
        sig_int = (num_atoms << 16)
        for ch in smiles:
            sig_int = (sig_int * 131 + ord(ch)) & 0x7FFFFFFF
        return (num_atoms, sig_int)
    except Exception:
        # Fallback: atom count and sentinel int
        try:
            return (mol.GetNumAtoms(), -1)
        except Exception:
            return (0, -1)


def _compute_all_ring_labels(mol) -> List[str]:
    """
    Compute ring labels for all atoms in the molecule at once.
    Returns a list where result[i] is the ring label for atom i.
    """
    num_atoms = mol.GetNumAtoms()
    labels = []
    for i in range(num_atoms):
        try:
            label = _flavonoid_ring_label_impl(mol, i)
            labels.append(label)
        except Exception:
            labels.append('')
    return labels


def flavonoid_ring_label(mol, atom_idx: int) -> str:
    """
    Return stable ring label for flavonoids: 'A', 'B', or 'C' if applicable; else ''.
    Enhanced implementation with strict C-ring identification and topological A/B assignment.
    Uses WeakKeyDictionary caching with signature-based invalidation when possible.
    On weak reference failures, computes directly on live molecule without caching.
    """
    try:
        # Try WeakKeyDictionary cache with signature-based invalidation
        current_signature = _compute_molecule_signature(mol)
        
        if mol in _mol_ring_labels_cache:
            cached_signature, cached_labels = _mol_ring_labels_cache[mol]
            
            # Check if signature matches (no topology changes)
            if cached_signature == current_signature:
                if atom_idx < len(cached_labels):
                    _cache_stats['hits'] += 1
                    return cached_labels[atom_idx]
            else:
                # Signature changed - invalidate and recompute
                _cache_stats['invalidations'] += 1
                try:
                    del _mol_ring_labels_cache[mol]
                except KeyError:
                    pass  # Already removed by concurrent access or GC
        
        # Cache miss or invalidation - compute labels and cache them
        _cache_stats['misses'] += 1
        labels = _compute_all_ring_labels(mol)
        _mol_ring_labels_cache[mol] = (current_signature, labels)
        
        if atom_idx < len(labels):
            return labels[atom_idx]
        return ''
        
    except (TypeError, ValueError):
        # WeakKeyDictionary failed (object not weak-ref-able)
        # Do NOT cache, just compute on live molecule and return
        try:
            return _flavonoid_ring_label_impl(mol, atom_idx)
        except Exception:
            return ''
    except Exception:
        # Other unexpected errors - fallback to direct computation
        try:
            return _flavonoid_ring_label_impl(mol, atom_idx)
        except Exception:
            return ''


def _flavonoid_ring_label_impl(mol, atom_idx: int) -> str:
    """
    Core implementation of flavonoid ring labeling logic.
    Separated from public interface to enable caching.
    """
    ensure_ready(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    six_rings = [r for r in atom_rings if len(r) == 6]
    
    if not six_rings:
        return ''
    
    # Find atom's 6-member ring
    atom_6ring = None
    for r in six_rings:
        if atom_idx in r:
            atom_6ring = r
            break
    if atom_6ring is None:
        return ''
    
    # C ring candidates: 6-ring with O and carbonyl inside
    c_candidates = []
    for r in six_rings:
        if any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in r):
            if _is_carbonyl_in_ring(mol, r):
                c_candidates.append(r)
    
    if not c_candidates:
        return ''
    
    # Prefer C candidate fused to more aromatic rings (flavonoid core vs sugar rings)
    aro6 = [r for r in six_rings if _ring_is_aromatic_6(mol, r)]
    best = None
    best_score = -1
    for r in c_candidates:
        score = sum(1 for ar in aro6 if _shared_atoms(r, ar) >= 2)
        if score > best_score:
            best = r
            best_score = score
    
    c_ring = best
    if c_ring is None:
        return ''
    
    # If atom is in C-ring, return 'C'
    if atom_idx in c_ring:
        return 'C'
    
    # A rings: aromatic 6-rings fused to C ring (share >= 2 atoms)
    a_rings = [r for r in aro6 if _shared_atoms(r, c_ring) >= 2]
    
    # B rings: aromatic 6-rings single-bond connected to C ring
    b_rings = [r for r in aro6 if _is_b_connected(mol, c_ring, r)]
    
    # Stable ordering by sum of atom indices
    def _ring_key(r):
        return sum(int(x) for x in r)
    
    a_rings.sort(key=_ring_key)
    b_rings.sort(key=_ring_key)
    
    # Check if atom's ring is A or B
    if any(atom_idx in r for r in a_rings[:1]):
        return 'A'
    if any(atom_idx in r for r in b_rings[:1]):
        return 'B'
    
    return ''


def get_flavonoid_ring_label_cache_info():
    """Return cache statistics for flavonoid ring labeling performance monitoring."""
    weak_cache_size = len(_mol_ring_labels_cache)
    
    # Create a cache_info-like object for compatibility
    class CacheInfo:
        def __init__(self, hits, misses, maxsize, currsize, invalidations=0):
            self.hits = hits
            self.misses = misses
            self.maxsize = maxsize
            self.currsize = currsize
            self.invalidations = invalidations
    
    # Return WeakKeyDictionary-only stats (no fallback cache)
    return CacheInfo(
        hits=_cache_stats['hits'],
        misses=_cache_stats['misses'],
        maxsize=None,  # WeakKeyDictionary has no maxsize limit
        currsize=weak_cache_size,
        invalidations=_cache_stats['invalidations']
    )


def reset_ring_cache():
    """Reset ring label cache to clean state."""
    _mol_ring_labels_cache.clear()
    _cache_stats['hits'] = 0
    _cache_stats['misses'] = 0
    _cache_stats['invalidations'] = 0


def invalidate_ring_cache(reason: str):
    """Invalidate ring cache with specified reason for debugging."""
    _mol_ring_labels_cache.clear()
    _cache_stats['invalidations'] += 1


def clear_flavonoid_ring_label_cache():
    """Clear the flavonoid ring labeling cache. Useful for testing or memory management."""
    _mol_ring_labels_cache.clear()
    _cache_stats['hits'] = 0
    _cache_stats['misses'] = 0
    _cache_stats['invalidations'] = 0


# PR2 C-ring Identification Helpers
# ==================================

def _ring_has_oxygen(mol: Chem.Mol, ring_atoms: Iterable[int]) -> bool:
    return any(mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in ring_atoms)

def _soft_c_ring_candidates(mol: Chem.Mol) -> Set[int]:
    """
    Soft fallback for C-ring detection when strict criteria find no matches.
    Identifies rings that are:
    - 5-6 member rings
    - <=1 heteroatom (typically oxygen)
    - Not aromatic
    This provides broader coverage while avoiding 'any ring' expansion.
    """
    atoms: Set[int] = set()
    try:
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            ring_size = len(ring)

            # Size constraint: 5-6 member rings
            if not (5 <= ring_size <= 6):
                continue

            # Heteroatom constraint: <=1 heteroatom (non-carbon, non-hydrogen)
            hetero_count = sum(1 for a in ring_atoms if a.GetAtomicNum() not in (6, 1))
            if hetero_count > 1:
                continue

            # Aromaticity constraint: not aromatic
            is_aromatic = any(a.GetIsAromatic() for a in ring_atoms)
            if is_aromatic:
                continue

            atoms.update(ring)
    except Exception:
        pass
    return atoms

def c_ring_membership_atoms(mol: Chem.Mol) -> Set[int]:
    """
    Return all atom indices in 6-member rings identified as C-rings (independent of H count/hybridization).
    Basic heuristic: 6-member ring AND (contains oxygen) AND (contains in-ring carbonyl).
    If no strict matches found, falls back to broader criteria for better test coverage.
    Note: This replaces historical CH=1-bound C-ring detection, making it reusable for CH2.
    """
    try:
        out: Set[int] = set()
        ring_info = mol.GetRingInfo()

        # First try strict criteria
        for ring in ring_info.AtomRings():
            if len(ring) != 6:
                continue
            if not _ring_has_oxygen(mol, ring):
                continue
            if not _is_carbonyl_in_ring(mol, ring):
                continue
            out.update(ring)

        # If strict criteria found no rings, use soft fallback
        if not out:
            out = _soft_c_ring_candidates(mol)

        return out
    except Exception:
        return set()


# Pure Carbon Ring Helper Functions
# ==================================

def _get_pure_carbon_rings(mol: Chem.Mol) -> List[Tuple[int, ...]]:
    """
    Identify pure carbon rings (carbocyclic rings with no heteroatoms).

    Args:
        mol: RDKit molecule

    Returns:
        List of ring tuples containing only carbon atoms
    """
    try:
        pure_carbon_rings = []
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            # Check if all atoms in ring are carbon
            if all(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in ring):
                pure_carbon_rings.append(ring)
        return pure_carbon_rings
    except Exception:
        return []


def _get_pure_carbon_ring_atoms(mol: Chem.Mol) -> Set[int]:
    """
    Get all atom indices that are part of pure carbon rings.

    Args:
        mol: RDKit molecule

    Returns:
        Set of atom indices in pure carbon rings
    """
    try:
        pure_carbon_atoms = set()
        for ring in _get_pure_carbon_rings(mol):
            pure_carbon_atoms.update(ring)
        return pure_carbon_atoms
    except Exception:
        return set()


def _is_beta_to_carbonyl(mol: Chem.Mol, carbon_idx: int) -> bool:
    """
    Check if a carbon atom is beta (two bonds away) from a carbonyl carbon.

    Args:
        mol: RDKit molecule
        carbon_idx: Index of carbon atom to check

    Returns:
        True if carbon is beta to any carbonyl
    """
    try:
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.GetSymbol() != 'C':
            return False

        # Look for carbonyls that are exactly 2 bonds away
        for neighbor in carbon_atom.GetNeighbors():
            for second_neighbor in neighbor.GetNeighbors():
                if (second_neighbor.GetIdx() != carbon_idx and
                    second_neighbor.GetSymbol() == 'C'):
                    # Check if this carbon is part of a carbonyl
                    for bond in second_neighbor.GetBonds():
                        other_atom = bond.GetOtherAtom(second_neighbor)
                        if (other_atom.GetSymbol() == 'O' and
                            bond.GetBondType() == Chem.BondType.DOUBLE):
                            return True
        return False
    except Exception:
        return False


def _has_ring_oxygen_neighbor(mol: Chem.Mol, carbon_idx: int) -> bool:
    """
    Check if a carbon atom has an oxygen neighbor that is also in a ring.

    Args:
        mol: RDKit molecule
        carbon_idx: Index of carbon atom to check

    Returns:
        True if carbon has a ring oxygen neighbor
    """
    try:
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.GetSymbol() != 'C':
            return False

        for neighbor in carbon_atom.GetNeighbors():
            if (neighbor.GetSymbol() == 'O' and neighbor.IsInRing()):
                return True
        return False
    except Exception:
        return False


# PR2 Site Identification Functions
# =================================

def c_ring_sp2_CH_sites(mol, masked_atoms: set) -> List[int]:
    """R2a: Ring aromatic sp2 CH sites in C-rings.

    Uses c_ring_membership_atoms() for proper C-ring targeting
    instead of accepting any ring atoms.
    """
    try:
        masked = set(masked_atoms or ())
        # Use proper C-ring detection for flavonoid C-ring targeting
        c_ring_atoms = c_ring_membership_atoms(mol)

        out: List[int] = []
        for idx in c_ring_atoms:
            if idx in masked:
                continue
            a = mol.GetAtomWithIdx(idx)
            if (a.GetSymbol() == "C" and
                a.GetIsAromatic() and
                a.GetHybridization() == Chem.HybridizationType.SP2 and
                a.GetTotalNumHs() == 1):
                out.append(idx)
        return out
    except Exception:
        return []



def c_ring_sp3_CH2_flavanone_sites(mol, masked_atoms: set, sugar_cfg: dict | None) -> List[int]:
    """R2b: Ring sp3 CH2 sites with dual condition targeting.

    Targets sp3 CH2 sites that satisfy either:
    1. Has ring oxygen neighbor (THP, oxygen-containing rings)
    2. Is beta to carbonyl (flavanone C3 position)

    Maintains sugar ring exclusion when sugar masking is enabled.
    """
    try:
        masked = set(masked_atoms or ())

        # Get all ring atoms for initial filtering
        ring_atoms = set()
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            ring_atoms.update(ring)

        # Compute sugar_ring_atoms only when sugar masking is actually enabled
        sugar_ring_atoms: Set[int] = set()
        if sugar_cfg is not None and sugar_cfg.get('mode', 'heuristic') not in ('off', 'none'):
            try:
                for ring in _find_sugar_rings_cached(mol, sugar_cfg):
                    sugar_ring_atoms.update(ring)
            except Exception:
                pass

        out: List[int] = []
        for idx in ring_atoms:
            if idx in masked or idx in sugar_ring_atoms:
                continue
            a = mol.GetAtomWithIdx(idx)
            if (a.GetSymbol() == "C" and
                a.GetHybridization() == Chem.HybridizationType.SP3 and
                a.GetTotalNumHs() == 2):
                # Dual condition: ring oxygen neighbor OR beta to carbonyl
                if (_has_ring_oxygen_neighbor(mol, idx) or
                    _is_beta_to_carbonyl(mol, idx)):
                    out.append(idx)
        return out
    except Exception:
        return []
