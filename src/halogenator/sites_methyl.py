# -*- coding: ascii -*-

from typing import List, Set, Dict, Any
from rdkit import Chem

def enumerate_methyl_sites(mol, masked_atoms: Set[int], allow_on_methoxy: bool=False, allow_allylic_methyl: bool=False) -> List[Dict[str, Any]]:
    sites = []

    for atom in mol.GetAtoms():
        if atom.GetIdx() in masked_atoms:
            continue
        if atom.GetSymbol() != "C":
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        # Support multi-step halogenation: CH3 -> CH2X -> CHX2
        if atom.GetTotalNumHs() < 1:
            continue

        # Categorize neighbors: halogens vs non-halogens
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        halogen_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetSymbol() in ('F', 'Cl', 'Br', 'I')]
        non_halogen_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetSymbol() not in ('F', 'Cl', 'Br', 'I')]

        # Must have exactly one non-halogen heavy neighbor (the main structure attachment)
        if len(non_halogen_neighbors) != 1:
            continue

        neighbor = non_halogen_neighbors[0]
        site_type = None

        # Determine site type and apply filtering
        if neighbor.GetSymbol() == "O":
            # This is a methoxy site (O-CH3)
            if allow_on_methoxy:
                site_type = "OCH3"
            else:
                continue  # Skip if methoxy not allowed

        elif (neighbor.GetSymbol() == "C" and
              neighbor.GetHybridization() == Chem.HybridizationType.SP2 and
              not neighbor.GetIsAromatic()):
            # This is an allylic/vinylic methyl site
            if allow_allylic_methyl:
                site_type = "ALLYLIC_CH3"
            else:
                continue  # Skip if allylic methyl not allowed
        else:
            # Other types of methyl sites (aromatic, sp3, etc.)
            site_type = "OTHER_CH3"

        if site_type:
            sites.append({
                'idx': atom.GetIdx(),
                'kind': site_type,
                'neighbor_idx': neighbor.GetIdx(),
                'neighbor_symbol': neighbor.GetSymbol()
            })

    return sites