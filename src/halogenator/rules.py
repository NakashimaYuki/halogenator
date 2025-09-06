# -*- coding: ascii -*-
"""Reaction SMARTS/SMIRKS templates."""

from typing import Dict
from rdkit.Chem import rdChemReactions as RXN

HALOGENS = ("F", "Cl", "Br", "I")

def _r1_smirks(X: str) -> str:
    """R1: Aromatic CH -> C-halogen"""
    return "[cH:1]>>[c:1]" + X

def _r3_smirks(X: str) -> str:
    """R3: Replace -OH (excluding carboxylic acid)"""
    return "[#6;!$(C=O):1]-[O;H1]>>[#6:1]" + X

def _r4_smirks(X: str) -> str:
    """R4: Replace -NHx (excluding amide)"""
    return "[#6;!$(C=O):1]-[N;H2,H1]>>[#6:1]" + X

def _r5_smirks(X: str) -> str:
    """R5: Replace -C(=O)OH (whole carboxyl group)"""
    return "[c,#6:1]-C(=O)[O;H1]>>[#6:1]" + X


def build_reactions() -> Dict[str, Dict[str, RXN.ChemicalReaction]]:
    """Build reaction templates from SMIRKS for each halogen."""
    rxns = {"R1": {}, "R3": {}, "R4": {}, "R5": {}}
    
    for X in HALOGENS:
        try:
            rxns["R1"][X] = RXN.ReactionFromSmarts(_r1_smirks(X))
        except Exception as e:
            print(f"Warning: Failed to create R1 reaction for {X}: {e}")
            
        try:
            rxns["R3"][X] = RXN.ReactionFromSmarts(_r3_smirks(X))
        except Exception as e:
            print(f"Warning: Failed to create R3 reaction for {X}: {e}")
            
        try:
            rxns["R4"][X] = RXN.ReactionFromSmarts(_r4_smirks(X))
        except Exception as e:
            print(f"Warning: Failed to create R4 reaction for {X}: {e}")
            
        try:
            rxns["R5"][X] = RXN.ReactionFromSmarts(_r5_smirks(X))
        except Exception as e:
            print(f"Warning: Failed to create R5 reaction for {X}: {e}")
    
    return rxns


def get_rule_description(rule: str) -> str:
    """Get ASCII description of rule."""
    descriptions = {
        'R1': 'Aromatic CH -> C-halogen',
        'R3': 'Replace -OH (not carboxylic acid)',
        'R4': 'Replace -NHx (not amide)', 
        'R5': 'Replace -C(=O)OH (whole carboxyl)',
    }
    return descriptions.get(rule, 'Unknown rule')
