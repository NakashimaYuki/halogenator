# -*- coding: ascii -*-
"""Site identification and symmetry folding."""

from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import SanitizeFlags


def ensure_ready(mol: Chem.Mol) -> None:
    """Ensure molecule is ready for ring info and property access."""
    try:
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:
        mol.UpdatePropertyCache(strict=False)
    
    # Force ring info initialization
    try:
        Chem.GetSymmSSSR(mol)
    except Exception:
        pass


def is_carbonyl_carbon(atom: Chem.Atom) -> bool:
    """Check if atom is carbonyl carbon (has O double bond neighbor)."""
    if atom.GetSymbol() != 'C':
        return False
    
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() == 'O':
            bond = atom.GetOwningMol().GetBondBetweenAtoms(
                atom.GetIdx(), neighbor.GetIdx()
            )
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                return True
    return False


def is_c_ring_site_ready(mol: Chem.Mol, aidx: int) -> bool:
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
    if is_carbonyl_carbon(atom):
        return False
    
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


def is_c_ring_site(mol: Chem.Mol, aidx: int) -> bool:
    """Check if atom is valid C ring site for R2."""
    ensure_ready(mol)
    return is_c_ring_site_ready(mol, aidx)


def aromatic_CH_indices(mol: Chem.Mol) -> List[int]:
    """Get aromatic carbon indices with H > 0."""
    ensure_ready(mol)
    
    indices = []
    for atom in mol.GetAtoms():
        try:
            if (atom.GetSymbol() == 'C' and 
                atom.GetIsAromatic() and 
                atom.GetTotalNumHs() > 0):
                indices.append(atom.GetIdx())
        except Exception:
            # Skip atoms that can't provide ring/hydrogen info
            continue
    return indices


def symmetry_groups(mol: Chem.Mol, indices: List[int], extra_tag: Optional[str] = None) -> Dict[Any, List[int]]:
    """Group indices by symmetry using CanonicalRankAtoms."""
    if not indices:
        return {}
    
    # Ensure molecule is ready for symmetry operations
    ensure_ready(mol)
    
    # Get canonical ranks for symmetry folding
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


def c_ring_indices(mol: Chem.Mol) -> List[int]:
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


def get_ring_tag_for_atom(mol: Chem.Mol, aidx: int) -> Optional[int]:
    """Get a stable ring tag for atoms in 6-member rings."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    for ring in atom_rings:
        if aidx in ring and len(ring) == 6:
            # Return minimum atom index in ring as stable tag
            return min(ring)
    
    return None