# -*- coding: ascii -*-
"""RDKit compatibility shims and patchable Chem alias.

This module exposes a single Chem alias that callers and tests can patch.
Functions provide safe fallbacks when RDKit is unavailable.
"""

import sys
import types

# Create Chem alias with minimal stubs if RDKit is not available
try:
    import rdkit.Chem as _rdkit_Chem  # type: ignore
    Chem = _rdkit_Chem
except Exception:
    Chem = types.ModuleType('Chem')

# Create RDLogger alias with minimal stubs if RDKit is not available
try:
    from rdkit import RDLogger as _rdkit_RDLogger  # type: ignore
    RDLogger = _rdkit_RDLogger
except Exception:
    RDLogger = types.ModuleType('RDLogger')
    
    def _DisableLog_stub(category):
        pass
    
    def _EnableLog_stub(category):
        pass
    
    # Attach stubs
    setattr(RDLogger, 'DisableLog', _DisableLog_stub)
    setattr(RDLogger, 'EnableLog', _EnableLog_stub)

    class _Mol:
        def __init__(self):
            pass

    class _RWMol(_Mol):
        def __init__(self, mol=None):
            super().__init__()

    def _MolToSmiles_stub(mol, canonical=True, isomericSmiles=False):
        return ''

    def _CanonicalRankAtoms_stub(mol, breakTies=False):
        try:
            return list(range(mol.GetNumAtoms()))
        except Exception:
            return []

    # Attach stubs
    setattr(Chem, 'Mol', _Mol)
    setattr(Chem, 'RWMol', _RWMol)
    setattr(Chem, 'MolToSmiles', _MolToSmiles_stub)
    setattr(Chem, 'CanonicalRankAtoms', _CanonicalRankAtoms_stub)

# Register to ease patching in tests
sys.modules[__name__ + '.Chem'] = Chem

def mol_to_smiles_safe(mol, canonical=True, isomericSmiles=False) -> str:
    try:
        return Chem.MolToSmiles(mol, canonical=canonical, isomericSmiles=isomericSmiles)
    except Exception:
        return ''

def canonical_rank_atoms_safe(mol, breakTies=False):
    try:
        return Chem.CanonicalRankAtoms(mol, breakTies=breakTies)
    except Exception:
        try:
            return list(range(mol.GetNumAtoms()))
        except Exception:
            return []

