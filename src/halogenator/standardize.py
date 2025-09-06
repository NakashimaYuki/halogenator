# -*- coding: ascii -*-
"""Molecule standardization utilities."""

from typing import Optional
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize


def std_from_smiles(smi: str, do_tautomer: bool = False) -> Optional[Chem.Mol]:
    """Standardize molecule from SMILES using rdMolStandardize."""
    if smi is None:
        return None
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        
        # Clean up and uncharge
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        
        # Sanitize
        try:
            Chem.SanitizeMol(mol)
        except (Chem.AtomValenceException, Chem.KekulizeException, ValueError):
            return None
        
        # Optional tautomer canonicalization
        if do_tautomer:
            try:
                te = rdMolStandardize.TautomerEnumerator()
                mol = te.Canonicalize(mol)
            except:
                # If tautomer canonicalization fails, continue without it
                pass
        
        return mol
        
    except Exception:
        return None


def to_inchikey(mol: Chem.Mol) -> str:
    """Convert molecule to InChIKey, fallback to canonical SMILES."""
    try:
        # Try InChI first
        inchi_str = inchi.MolToInchi(mol)
        if inchi_str:
            key = inchi.InchiToInchiKey(inchi_str)
            if key:
                return key
    except:
        pass
    
    # Fallback to canonical SMILES
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return "UNKNOWN"