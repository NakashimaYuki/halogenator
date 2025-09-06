# -*- coding: ascii -*-
"""Deduplication utilities."""

from typing import List
from rdkit import Chem
from .standardize import to_inchikey


def seen_key(mol: Chem.Mol) -> str:
    """Get deduplication key for molecule (InChIKey preferred, SMILES fallback)."""
    if mol is None:
        return "NONE"
    
    return to_inchikey(mol)


def dedupe_mols(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """Deduplicate molecules using InChIKey/SMILES keys."""
    if not mols:
        return []
    
    seen_keys = set()
    unique_mols = []
    
    for mol in mols:
        if mol is None:
            continue
        
        key = seen_key(mol)
        if key not in seen_keys:
            seen_keys.add(key)
            unique_mols.append(mol)
    
    return unique_mols


def dedupe_with_props(mol_prop_pairs: List[tuple]) -> List[tuple]:
    """Deduplicate (mol, props) pairs, keeping first occurrence."""
    if not mol_prop_pairs:
        return []
    
    seen_keys = set()
    unique_pairs = []
    
    for mol, props in mol_prop_pairs:
        if mol is None:
            continue
        
        key = seen_key(mol)
        if key not in seen_keys:
            seen_keys.add(key)
            unique_pairs.append((mol, props))
    
    return unique_pairs