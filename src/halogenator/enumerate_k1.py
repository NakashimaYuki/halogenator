# -*- coding: ascii -*-
"""k=1 enumeration core logic."""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import SanitizeFlags
from rdkit.Chem.rdChemReactions import ChemicalReaction

from .rules import build_reactions
from .sites import aromatic_CH_indices, c_ring_indices, symmetry_groups, get_ring_tag_for_atom, ensure_ready
from .reactions import apply_single_site_halogenation
from .qc import sanitize_ok, basic_descriptors, pains_flags
from .dedupe import dedupe_with_props
from .standardize import to_inchikey


def enumerate_k1_halogenation(parent_mol: Chem.Mol, halogens: List[str], rules: List[str], config: Dict[str, Any]) -> List[Tuple[Chem.Mol, Dict[str, Any]]]:
    """Enumerate k=1 halogenated products with properties."""
    if parent_mol is None:
        return []
    
    # Ensure molecule is properly sanitized before any operations
    try:
        Chem.SanitizeMol(parent_mol)
    except Exception:
        parent_mol.UpdatePropertyCache(strict=False)
    
    parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)
    parent_inchikey = to_inchikey(parent_mol)
    
    products = []
    reactions = build_reactions()
    
    # Step 1: Apply reaction-based rules (R1, R3, R4, R5)
    for rule in ['R1', 'R3', 'R4', 'R5']:
        if rule in rules and rule in reactions:
            for halogen in halogens:
                if halogen in reactions[rule]:
                    rxn = reactions[rule][halogen]
                    local_seen = set()  # Local deduplication per (rule, halogen)
                    try:
                        reaction_products = rxn.RunReactants((parent_mol,))
                        for product_set in reaction_products:
                            for product_mol in product_set:
                                if product_mol is not None:
                                    # Sanitize product molecule (lightweight first, full if needed)
                                    try:
                                        Chem.SanitizeMol(product_mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES)
                                    except Exception:
                                        try:
                                            Chem.SanitizeMol(product_mol)
                                        except Exception:
                                            product_mol.UpdatePropertyCache(strict=False)
                                    
                                    # Ensure ring info and properties are ready for downstream operations
                                    try:
                                        ensure_ready(product_mol)
                                    except Exception:
                                        pass  # Continue with potentially incomplete initialization
                                    
                                    try:
                                        product_smiles = Chem.MolToSmiles(product_mol, canonical=True)
                                    except Exception:
                                        # fallback: skip local dedup on this product, but keep pipeline alive
                                        product_smiles = None
                                    
                                    if product_smiles is None or product_smiles not in local_seen:
                                        if product_smiles is not None:
                                            local_seen.add(product_smiles)
                                        props = {
                                            'parent_smiles': parent_smiles,
                                            'parent_inchikey': parent_inchikey,
                                            'rule': rule,
                                            'site_index': None,
                                            'sym_class': None,
                                            'ring_tag': None,
                                            'halogen': halogen,
                                            'depth': 1
                                        }
                                        products.append((product_mol, props))
                    except Exception as e:
                        print(f"Warning: Rule {rule} with {halogen} failed: {e}")
    
    # Step 2: Apply site-based rules with symmetry folding (R2)  
    # Note: R1 now uses reaction engine above, may add back symmetry folding in future
    
    if 'R2' in rules:
        r2_products = _apply_r2_rule(
            parent_mol, parent_smiles, parent_inchikey, halogens
        )
        products.extend(r2_products)
    
    # Step 3: QC filtering
    if config.get('qc', {}).get('sanitize_strict', True):
        products = [(mol, props) for mol, props in products if sanitize_ok(mol)]
    
    # Step 4: Deduplication
    if config.get('dedupe', {}).get('method', 'inchikey') == 'inchikey':
        products = dedupe_with_props(products)
    
    # Step 5: Add descriptors and final properties
    final_products = []
    for mol, props in products:
        if mol is None:
            continue
        
        # Add descriptors
        descriptors = basic_descriptors(mol)
        props.update(descriptors)
        
        # Add QC flags
        props['sanitize_ok'] = sanitize_ok(mol)
        props['pains_flags'] = pains_flags(mol)
        
        final_products.append((mol, props))
    
    return final_products






def _apply_r2_rule(parent_mol: Chem.Mol, parent_smiles: str, parent_inchikey: str,
                   halogens: List[str]) -> List[Tuple[Chem.Mol, Dict[str, Any]]]:
    """Apply R2 rule with symmetry folding."""
    products = []
    
    # Ensure parent molecule is ready for ring operations
    ensure_ready(parent_mol)
    
    # Get C ring indices
    cring_indices = c_ring_indices(parent_mol)
    if not cring_indices:
        return products
    
    # Group by symmetry with ring tags
    indexed_sites = []
    for idx in cring_indices:
        ring_tag = get_ring_tag_for_atom(parent_mol, idx)
        indexed_sites.append((idx, ring_tag))
    
    # Group by (rank, ring_tag)
    ranks = Chem.CanonicalRankAtoms(parent_mol, breakTies=False)
    groups = {}
    
    for idx, ring_tag in indexed_sites:
        rank = ranks[idx]
        key = (rank, ring_tag)
        if key not in groups:
            groups[key] = []
        groups[key].append(idx)
    
    # For each group, take representative
    for (rank, ring_tag), indices in groups.items():
        representative_idx = indices[0]
        
        # Apply each halogen
        for halogen in halogens:
            product_mol = apply_single_site_halogenation(parent_mol, representative_idx, halogen)
            if product_mol is not None:
                props = {
                    'parent_smiles': parent_smiles,
                    'parent_inchikey': parent_inchikey,
                    'rule': 'R2',
                    'site_index': representative_idx,
                    'sym_class': str(rank),
                    'ring_tag': str(ring_tag) if ring_tag is not None else None,
                    'halogen': halogen,
                    'depth': 1
                }
                products.append((product_mol, props))
    
    return products