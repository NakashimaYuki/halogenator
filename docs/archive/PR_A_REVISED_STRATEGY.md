# PR-A Revised Strategy: Manual Atom Mapping After Reaction
**Date**: 2025-11-11
**Issue**: RDKit RunReactants() does not preserve atom mapping numbers in products

---

## Problem Discovery

### RDKit Behavior
Running the debug script revealed:
```
Product atom mappings (BEFORE sanitization):
  Atom 0: C mapnum=0
  Atom 1: O mapnum=0   ← Should be mapnum=2, but it's 0!
  ...

Mapped SMILES: c1(OC)ccccc1  ← No :2 in output
```

**Root Cause**: RDKit's `RunReactants()` **strips all atom mapping numbers** from products by default. This is documented RDKit behavior.

---

## Revised Strategy: Post-Reaction Mapping

Instead of relying on SMIRKS to preserve mappings, we **manually set mapping numbers after the reaction**.

### Approach
1. **Keep SMIRKS unchanged** (current versions are fine)
2. **After `rxn.RunReactants()`**:
   - Identify the product atom(s) that should be highlighted
   - Set `atom.SetAtomMapNum(777)` on those atoms
3. **Generate mapped SMILES** with `Chem.MolToSmiles(mol, canonical=False)`
4. **Extract mapnums** and save to data

---

## Implementation Plan

### Step 1: Add Helper Function to Identify Product Atoms

```python
def identify_product_site_atoms(
    prod_mol: Chem.Mol,
    match: Tuple[int, ...],
    xf_label: str
) -> List[int]:
    """
    Identify which product atom(s) should be highlighted for this transformation.

    Args:
        prod_mol: Product molecule from reaction
        match: Tuple of reactant atom indices that matched query
        xf_label: Transformation label (e.g., "OH->OMe")

    Returns:
        List of product atom indices to highlight
    """
    # For OH->OMe and OH->NH2:
    # match = (carbon_idx, oxygen_idx) from reactant
    # In product, the oxygen/nitrogen is at a similar position

    if 'OMe' in xf_label or 'NH2' in xf_label:
        # Find the heteroatom (O or N) in product
        # Strategy: Look for O or N connected to the matched carbon

        carbon_idx = match[0]  # First atom in match is always the aromatic carbon

        # In product, find O or N connected to carbon at position ~carbon_idx
        # Note: Atom indices may shift, so we use chemical logic

        target_symbol = 'O' if 'OMe' in xf_label else 'N'

        for atom in prod_mol.GetAtoms():
            if atom.GetSymbol() != target_symbol:
                continue

            # Check if connected to aromatic carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic():
                    # Additional check for OMe: must be ether oxygen
                    if target_symbol == 'O':
                        if atom.GetDegree() == 2:
                            # Check if other neighbor is alkyl carbon
                            others = [n for n in atom.GetNeighbors() if n.GetIdx() != neighbor.GetIdx()]
                            if others and others[0].GetSymbol() == 'C' and not others[0].GetIsAromatic():
                                return [atom.GetIdx()]

                    # For NH2: primary amine connected to aromatic
                    elif target_symbol == 'N':
                        if atom.GetTotalNumHs() >= 2 and not atom.IsInRing():
                            return [atom.GetIdx()]

    return []
```

### Step 2: Modify apply_to_molecule() in TransformationEngineV2

In `scripts/08_transform_library_v2.py`, around line 477:

```python
try:
    prod_mol = rxn_products[site_idx][0]

    # Sanitize
    Chem.SanitizeMol(prod_mol)

    # NEW: Identify and mark product site atoms with mapping number 777
    site_atom_indices = identify_product_site_atoms(
        prod_mol,
        match,
        self.transform['label']
    )

    # Set mapping number 777 on product site atoms
    for atom_idx in site_atom_indices:
        atom = prod_mol.GetAtomWithIdx(atom_idx)
        atom.SetAtomMapNum(777)

    # Get canonical SMILES (fast)
    prod_smiles_canon = canonical_smiles(prod_mol)
    if not prod_smiles_canon:
        continue

    # NEW: Get mapped SMILES (preserve atom mappings)
    prod_smiles_mapped = Chem.MolToSmiles(prod_mol, canonical=False)

    # InChIKey generation
    prod_inchikey = get_inchikey(prod_mol)

    # Quick properties only
    props = calculate_properties_quick(prod_mol)

    # Calculate halogen fields
    halogen_fields = calculate_halogen_fields(prod_mol, self.transform['label'])

    products.append({
        'smiles': prod_smiles_canon,
        'canonical_smiles': prod_smiles_canon,
        'inchikey': prod_inchikey,
        'xf_success': True,
        'xf_reason': None,
        'xf_label': self.transform['label'],
        'xf_rule_id': self.transform['name'],
        'xf_site_index': site_idx,
        'xf_site_atoms': str(match),
        # NEW FIELDS for PR-A:
        'product_mapped_smiles': prod_smiles_mapped,
        'xf_site_mapnums': json.dumps([777] * len(site_atom_indices)),  # List of 777s
        'source_smiles': smiles,
        **source_data,
        **props,
        **halogen_fields
    })

except Exception as e:
    logger.debug(f"Site {site_idx} failed: {e}")
    continue
```

---

## Advantages of This Approach

1. ✅ **Full control**: We decide exactly which atoms get marked
2. ✅ **No RDKit limitations**: Don't rely on mapping preservation
3. ✅ **Flexible**: Can use different logic for different reaction types
4. ✅ **Debuggable**: Clear code path for atom identification
5. ✅ **Works with site_finder**: The `:777` will be in mapped SMILES

---

## Testing Strategy

### Test 1: Manual Mapping After Reaction

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Reaction
rxn = AllChem.ReactionFromSmarts('[c:1][OX2H:2]>>[c:1][O:2]C')
mol = Chem.MolFromSmiles('c1ccccc1O')
products = rxn.RunReactants((mol,))

prod_mol = products[0][0]
Chem.SanitizeMol(prod_mol)

# Find the oxygen in methoxy group
for atom in prod_mol.GetAtoms():
    if atom.GetSymbol() == 'O':
        # Check if it's an ether oxygen connected to aromatic + alkyl
        neighbors = list(atom.GetNeighbors())
        if len(neighbors) == 2:
            has_aromatic = any(n.GetIsAromatic() for n in neighbors)
            has_alkyl = any(not n.GetIsAromatic() for n in neighbors)
            if has_aromatic and has_alkyl:
                # Mark this atom
                atom.SetAtomMapNum(777)
                print(f"Marked atom {atom.GetIdx()} (O) with :777")
                break

# Get mapped SMILES
mapped_smiles = Chem.MolToSmiles(prod_mol, canonical=False)
print(f"Mapped SMILES: {mapped_smiles}")

# Verify :777 is present
if ':777' in mapped_smiles:
    print("[SUCCESS] :777 found in mapped SMILES!")
else:
    print("[FAIL] :777 not found")
```

---

## Updated File Changes

### File: `scripts/08_transform_library_v2.py`

**New function** (add before TransformationEngineV2 class):
```python
def identify_product_site_atoms(
    prod_mol: Chem.Mol,
    match: Tuple[int, ...],
    xf_label: str
) -> List[int]:
    """Identify product atoms to highlight (see implementation above)."""
    # ... (full implementation)
```

**Modify `apply_to_molecule()` method** (line ~477):
- After sanitization, call `identify_product_site_atoms()`
- Set `atom.SetAtomMapNum(777)` on identified atoms
- Generate `prod_smiles_mapped` with `Chem.MolToSmiles(mol, canonical=False)`
- Add fields to product dictionary

**No changes needed to**:
- `configs/transforms.yaml` (SMIRKS stay the same)
- Reaction queries or patterns

---

## Expected Outcome

With manual mapping after reaction:

```
Product atom mappings (AFTER manual marking):
  Atom 0: C mapnum=0
  Atom 1: O mapnum=777  ← Manually set!
  Atom 2: C mapnum=0
  ...

Mapped SMILES: c1(O:777]C)ccccc1  ← Contains :777
```

Then downstream `site_finder.resolve_sites_from_mapnum()` can extract atoms with mapnum=777.

---

## Implementation Checklist

- [ ] Write `identify_product_site_atoms()` helper function
- [ ] Add unit test for manual mapping
- [ ] Integrate into `TransformationEngineV2.apply_to_molecule()`
- [ ] Update product dictionary with new fields
- [ ] Test on 10-20 molecules
- [ ] Verify mapped SMILES contains `:777`
- [ ] Run regression test with site_finder
- [ ] Achieve 100% accuracy

---

**Status**: Ready to implement revised strategy
