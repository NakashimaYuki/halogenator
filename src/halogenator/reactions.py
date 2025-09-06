# -*- coding: ascii -*-
"""Reaction processing utilities."""

from typing import Optional
from rdkit import Chem
from rdkit.Chem import RWMol


def replace_star_with_halogen(prod_mol: Chem.Mol, halogen_symbol: str) -> Optional[Chem.Mol]:
    """Replace * (dummy atom) with halogen in reaction product."""
    if prod_mol is None:
        return None
    
    try:
        # Create editable molecule
        rwmol = RWMol(prod_mol)
        
        # Find dummy atoms (*)
        dummy_indices = []
        for atom in rwmol.GetAtoms():
            if atom.GetSymbol() == '*' or atom.GetAtomicNum() == 0:
                dummy_indices.append(atom.GetIdx())
        
        # Replace dummy atoms with halogen (process in reverse order to avoid index shifts)
        for idx in reversed(dummy_indices):
            atom = rwmol.GetAtomWithIdx(idx)
            atom.SetAtomicNum(_get_atomic_num(halogen_symbol))
            atom.SetFormalCharge(0)
        
        # Sanitize
        mol = rwmol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
        
    except Exception:
        return None


def apply_single_site_halogenation(mol: Chem.Mol, atom_idx: int, halogen_symbol: str) -> Optional[Chem.Mol]:
    """Apply single site halogenation by adding halogen to carbon."""
    if mol is None:
        return None
    
    try:
        # Create editable molecule
        rwmol = RWMol(mol)
        
        # Get the target carbon atom
        if atom_idx >= rwmol.GetNumAtoms():
            return None
        
        carbon_atom = rwmol.GetAtomWithIdx(atom_idx)
        if carbon_atom.GetSymbol() != 'C':
            return None
        
        # Add halogen atom
        halogen_idx = rwmol.AddAtom(Chem.Atom(_get_atomic_num(halogen_symbol)))
        
        # Add single bond between carbon and halogen
        rwmol.AddBond(atom_idx, halogen_idx, Chem.BondType.SINGLE)
        
        # Remove one hydrogen from the carbon (if it has any)
        if carbon_atom.GetTotalNumHs() > 0:
            # The hydrogen removal is handled implicitly by RDKit sanitization
            # when we add the new bond
            pass
        
        # Sanitize
        mol = rwmol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
        
    except Exception:
        return None


def _get_atomic_num(halogen_symbol: str) -> int:
    """Get atomic number for halogen symbol."""
    halogen_map = {
        'F': 9,
        'Cl': 17,
        'Br': 35, 
        'I': 53
    }
    return halogen_map.get(halogen_symbol, 9)  # Default to F if unknown


def validate_halogen(halogen_symbol: str) -> bool:
    """Validate that symbol is a supported halogen."""
    return halogen_symbol in ['F', 'Cl', 'Br', 'I']