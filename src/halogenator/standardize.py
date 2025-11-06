# -*- coding: ascii -*-
"""Molecule standardization utilities."""

from typing import Optional


def std_from_smiles(smi: str, do_tautomer: bool = False):
    """Standardize molecule from SMILES using rdMolStandardize."""
    if smi is None:
        return None
    
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except Exception:
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
            except Exception:
                # If tautomer canonicalization fails, continue without it
                pass
        
        return mol
        
    except Exception:
        return None


def strip_isotopes_and_props(mol):
    """Create a copy of molecule with isotopes and atom mappings removed for stable InChIKey generation.

    P0 STABILITY FIX: Ensures consistent InChIKey generation by removing isotope labels
    and atom mappings that can vary between equivalent molecules.
    """
    try:
        from rdkit import Chem
        dm = Chem.Mol(mol)
        for atom in dm.GetAtoms():
            if atom.GetIsotope():
                atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
        dm.ClearComputedProps()
        Chem.SanitizeMol(dm)
        return dm
    except Exception:
        return mol

def to_inchikey_sanitized(mol) -> str:
    """Convert molecule to InChIKey with isotopes and mappings stripped for consistency."""
    try:
        from rdkit.Chem import inchi
        cleaned_mol = strip_isotopes_and_props(mol)
        inchi_str = inchi.MolToInchi(cleaned_mol)
        if inchi_str:
            key = inchi.InchiToInchiKey(inchi_str)
            if key:
                return key
    except Exception:
        pass
    return None

def to_inchikey(mol) -> str:
    """Convert molecule to InChIKey, fallback to canonical SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
    except Exception:
        return "UNKNOWN"

    try:
        # Try InChI first
        inchi_str = inchi.MolToInchi(mol)
        if inchi_str:
            key = inchi.InchiToInchiKey(inchi_str)
            if key:
                return key
    except Exception:
        pass
    
    # Fallback to canonical SMILES
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return "UNKNOWN"