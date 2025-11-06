# -*- coding: ascii -*-

from typing import Optional, Union, Dict, Any
from rdkit import Chem

_ALLOWED_STEP = {"F", "Cl"}
_ALLOWED_MACRO = {"CF3", "CCl3"}

def _get_atomic_num(halogen_symbol: str) -> int:
    """Get atomic number for halogen symbol."""
    halogen_map = {
        'F': 9,
        'Cl': 17,
        'Br': 35,
        'I': 53
    }
    return halogen_map.get(halogen_symbol, 9)  # Default to F if unknown

def apply_methyl_step(mol, site: Union[int, Dict[str, Any]], X: str) -> Optional[Chem.Mol]:
    """
    Replace one H on CH3 carbon with X (F or Cl).

    Args:
        mol: Source molecule
        site: Either carbon index (int) for backward compatibility,
              or site descriptor dict with 'idx' and 'kind' fields
        X: Halogen symbol ('F' or 'Cl')

    Return a new molecule or None if invalid.
    """
    if X not in _ALLOWED_STEP:
        return None

    if mol is None:
        return None

    # Handle both old integer format and new structured format
    if isinstance(site, int):
        cidx = site
        site_kind = "OTHER_CH3"  # Default for backward compatibility
    elif isinstance(site, dict):
        cidx = site['idx']
        site_kind = site.get('kind', 'OTHER_CH3')
    else:
        return None

    try:
        # Create editable molecule
        rwmol = Chem.RWMol(mol)

        # Validate target carbon exists and is carbon
        if cidx >= rwmol.GetNumAtoms():
            return None

        carbon_atom = rwmol.GetAtomWithIdx(cidx)
        if carbon_atom.GetSymbol() != 'C':
            return None

        # Check that this carbon has hydrogens available
        if carbon_atom.GetTotalNumHs() == 0:
            return None

        # Site-specific processing
        if site_kind == "OCH3":
            # Methoxy halogenation: -OCH3 -> -OCH2X
            return _apply_methoxy_halogenation(rwmol, cidx, X)
        elif site_kind == "ALLYLIC_CH3":
            # Allylic methyl halogenation: C=C-CH3 -> C=C-CH2X
            return _apply_allylic_methyl_halogenation(rwmol, cidx, X)
        else:
            # Generic methyl halogenation
            return _apply_generic_methyl_halogenation(rwmol, cidx, X)

    except Exception:
        return None

def _apply_methoxy_halogenation(rwmol, cidx: int, X: str) -> Optional[Chem.Mol]:
    """Apply halogenation to methoxy group: -OCH3 -> -OCH2X"""
    try:
        # Add halogen atom
        halogen_idx = rwmol.AddAtom(Chem.Atom(_get_atomic_num(X)))

        # Add single bond between carbon and halogen
        rwmol.AddBond(cidx, halogen_idx, Chem.BondType.SINGLE)

        # Sanitize
        mol_result = rwmol.GetMol()
        Chem.SanitizeMol(mol_result)
        return mol_result
    except Exception:
        return None

def _apply_allylic_methyl_halogenation(rwmol, cidx: int, X: str) -> Optional[Chem.Mol]:
    """Apply halogenation to allylic methyl: C=C-CH3 -> C=C-CH2X"""
    try:
        # Add halogen atom
        halogen_idx = rwmol.AddAtom(Chem.Atom(_get_atomic_num(X)))

        # Add single bond between carbon and halogen
        rwmol.AddBond(cidx, halogen_idx, Chem.BondType.SINGLE)

        # Sanitize
        mol_result = rwmol.GetMol()
        Chem.SanitizeMol(mol_result)
        return mol_result
    except Exception:
        return None

def _apply_generic_methyl_halogenation(rwmol, cidx: int, X: str) -> Optional[Chem.Mol]:
    """Apply generic methyl halogenation: -CH3 -> -CH2X"""
    try:
        # Add halogen atom
        halogen_idx = rwmol.AddAtom(Chem.Atom(_get_atomic_num(X)))

        # Add single bond between carbon and halogen
        rwmol.AddBond(cidx, halogen_idx, Chem.BondType.SINGLE)

        # Sanitize
        mol_result = rwmol.GetMol()
        Chem.SanitizeMol(mol_result)
        return mol_result
    except Exception:
        return None

def apply_methyl_macro(mol, cidx: int, label: str) -> Optional[Chem.Mol]:
    """
    Macro substitution at methyl: CF3 or CCl3.
    Build by replacing three H with X atoms in a single shot.
    """
    if label not in _ALLOWED_MACRO:
        return None

    if mol is None:
        return None

    X = "F" if label == "CF3" else "Cl"

    try:
        # Create editable molecule
        rwmol = Chem.RWMol(mol)

        # Validate target carbon exists and is carbon
        if cidx >= rwmol.GetNumAtoms():
            return None

        carbon_atom = rwmol.GetAtomWithIdx(cidx)
        if carbon_atom.GetSymbol() != 'C':
            return None

        # Check that this carbon has at least 3 hydrogens available
        if carbon_atom.GetTotalNumHs() < 3:
            return None

        # Add three halogen atoms
        for _ in range(3):
            halogen_idx = rwmol.AddAtom(Chem.Atom(_get_atomic_num(X)))
            rwmol.AddBond(cidx, halogen_idx, Chem.BondType.SINGLE)

        # The hydrogen removal is handled implicitly by RDKit sanitization
        # when we add the new bonds

        # Sanitize
        mol_result = rwmol.GetMol()
        Chem.SanitizeMol(mol_result)
        return mol_result

    except Exception:
        return None

def validate_methyl_halogen(halogen_symbol: str) -> bool:
    """Validate that symbol is a supported halogen for methyl halogenation."""
    return halogen_symbol in _ALLOWED_STEP