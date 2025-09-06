# -*- coding: ascii -*-
"""Quality control utilities."""

from typing import Dict, Any
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors


def sanitize_ok(mol: Chem.Mol) -> bool:
    """Check if molecule sanitizes OK."""
    if mol is None:
        return False
    
    try:
        # Try to sanitize a copy
        test_mol = Chem.Mol(mol)
        Chem.SanitizeMol(test_mol)
        return True
    except (Chem.AtomValenceException, Chem.KekulizeException, ValueError):
        return False
    except Exception:
        return False


def basic_descriptors(mol: Chem.Mol) -> Dict[str, Any]:
    """Calculate basic molecular descriptors."""
    if mol is None:
        return {
            'mw': 0.0,
            'logp': 0.0, 
            'tpsa': 0.0,
            'hba': 0,
            'hbd': 0,
            'aromatic_rings': 0
        }
    
    try:
        # Molecular weight
        mw = Descriptors.MolWt(mol)
        
        # LogP (octanol-water partition coefficient)
        logp = Crippen.MolLogP(mol)
        
        # Topological polar surface area
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        
        # Hydrogen bond acceptors and donors
        hba = Lipinski.NumHBA(mol)
        hbd = Lipinski.NumHBD(mol)
        
        # Count aromatic rings
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        
        return {
            'mw': round(mw, 2),
            'logp': round(logp, 2),
            'tpsa': round(tpsa, 2), 
            'hba': hba,
            'hbd': hbd,
            'aromatic_rings': aromatic_rings
        }
    
    except Exception:
        return {
            'mw': 0.0,
            'logp': 0.0,
            'tpsa': 0.0,
            'hba': 0,
            'hbd': 0,
            'aromatic_rings': 0
        }


def pains_flags(mol: Chem.Mol) -> int:
    """Check PAINS flags. Returns count of PAINS alerts (0 = no alerts)."""
    if mol is None:
        return 0
    
    # TODO: Implement PAINS filtering using RDKit FilterCatalog
    # For P0, return 0 (no PAINS alerts detected)
    # This can be enhanced later with proper PAINS implementation
    return 0


def passes_qc(mol: Chem.Mol, strict: bool = True) -> bool:
    """Check if molecule passes basic QC requirements."""
    if not sanitize_ok(mol):
        return False
    
    if strict:
        descriptors = basic_descriptors(mol)
        
        # Basic drug-like filters (can be adjusted)
        if descriptors['mw'] > 800 or descriptors['mw'] < 100:
            return False
        if descriptors['logp'] > 6 or descriptors['logp'] < -3:
            return False
        if descriptors['tpsa'] > 200:
            return False
    
    return True