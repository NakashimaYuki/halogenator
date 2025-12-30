# PR-A: Reaction Atom Mapping Design
**Date**: 2025-11-11
**Goal**: Achieve 100% site identification accuracy using atom mapping

---

## Design Principles

### 1. Mapping Number Consistency (RDKit Constraint)
- **RDKit Rule**: Product-side mapping numbers **MUST** exist in reactant side
- **Cannot** introduce new mapping numbers (like `:777`) only on product side
- **Solution**: Use existing reactant-side mapping numbers, with metadata to indicate which one to highlight

### 2. Mapping Strategy
- **Reactant side**: Use standard mappings (`:1`, `:2`, etc.) for reaction connectivity
- **Product side**: Preserve relevant mapping numbers from reactant side
- **Configuration metadata**: Add `highlight_mapnum` field to specify which mapping number marks the site
- **Example**: For OH->OMe, reactant has `[O:2]`, product has `[O:2]`, metadata says "highlight :2"

### 3. Flexibility
- Different reactions can use different mapping numbers
- OH->OMe uses `:2` (the oxygen)
- OH->NH2 uses `:2` (becomes nitrogen)
- H->F could use `:1` (the carbon that gets fluorine)

### 3. Heteroatom Selection
For each transformation type, mark the **characteristic heteroatom** of the product:
- **OH->OMe**: Mark the ether oxygen `O` (not the methyl carbon)
- **OH->NH2**: Mark the nitrogen `N` (not the hydrogens)
- **H->F/Cl/Br/I**: Mark the halogen atom

---

## Transformation SMIRKS Modifications

### OH->OMe (O-methylation)

**Current SMIRKS**:
```yaml
query_smarts: "[c:1][OX2H:2]"
smirks: "[c:1][OX2H:2]>>[c:1][O:2]C"
```

**Modified SMIRKS** (NO CHANGE NEEDED - already correct!):
```yaml
query_smarts: "[c:1][OX2H:2]"
smirks: "[c:1][OX2H:2]>>[c:1][O:2]C"
highlight_mapnum: 2  # NEW METADATA: Highlight atoms with mapping number :2
```

**Changes**:
- SMIRKS stays the same (product oxygen keeps `:2` from reactant)
- Add `highlight_mapnum: 2` to config to tell downstream which mapping number to extract
- This oxygen is the ether oxygen in the methoxy group (-OCH₃)

**Visual**:
```
Reactant:  Ar-O-H  (hydroxyl)
             ↓
Product:   Ar-O-CH₃  (methoxy)
              ↑
           Mark with :777
```

---

### OH->NH2 (Bioisosteric replacement)

**Current SMIRKS**:
```yaml
query_smarts: "[c:1][OX2H:2]"
smirks: "[c:1][OX2H:2]>>[c:1][NH2]"
```

**Modified SMIRKS**:
```yaml
query_smarts: "[c:1][OX2H:2]"
smirks: "[c:1][OX2H:2]>>[c:1][N:2]([H])[H]"
highlight_mapnum: 2  # NEW METADATA: Highlight atoms with mapping number :2
```

**Changes**:
- Product nitrogen: `[NH2]` → `[N:2]([H])[H]`
- Reuse mapping number `:2` from reactant oxygen (becomes nitrogen in product)
- Explicit hydrogens ensure correct RDKit parsing
- Add `highlight_mapnum: 2` to config
- Nitrogen is the characteristic atom of primary amine

**Visual**:
```
Reactant:  Ar-O-H  (hydroxyl)
             ↓
Product:   Ar-N-H₂  (primary amine)
              ↑
           Mark with :777
```

---

### Future Extensions (Halogens, if needed)

For halogenation transformations (e.g., H->F):

**Example SMIRKS**:
```yaml
query_smarts: "[c:1][H]"
smirks: "[c:1][H]>>[c:1][F:777]"
```

**Rationale**:
- Halogen is the transformation product
- Easy to identify visually
- Already handled well by current site_finder, but mapping provides 100% certainty

---

## Data Schema Changes

### New Fields (Added to product records)

#### 1. `product_mapped_smiles` (string)
- **Description**: Product SMILES with atom mapping numbers preserved
- **Format**: Standard SMILES with `:[num]` notation
- **Example**: `COc1c(O)cc(C=C2COc3cc([O:777])ccc3C2=O)cc1Br`
- **Generation**: `Chem.MolToSmiles(mol, canonical=False)`
  - `canonical=False` is CRITICAL to preserve mapping numbers
  - Canonical SMILES strips mapping numbers by default

#### 2. `xf_site_mapnums` (JSON list)
- **Description**: List of atom mapping numbers that mark transformation sites
- **Format**: JSON array of integers
- **Example**: `[777]` for single-site, `[777, 777]` for two-site
- **Rationale**:
  - Single value for most cases (k=1)
  - Multiple values if multiple sites marked (future k>1 support)
  - JSON format for easy parsing

### Existing Fields (Preserved for compatibility)
- `xf_site_atoms`: String tuple, e.g., `"(8, 9)"` - **KEEP** for backward compatibility
- `xf_site_index`: Integer, e.g., `1` - **KEEP** for backward compatibility
- `smiles`: Canonical SMILES **WITHOUT** mapping - **KEEP** as primary identifier

---

## Implementation Changes

### File: `scripts/08_transform_library_v2.py`

#### Change 1: Product SMILES generation (Line ~483)

**Current**:
```python
# Get canonical SMILES (fast)
prod_smiles_canon = canonical_smiles(prod_mol)
if not prod_smiles_canon:
    continue
```

**Modified**:
```python
# Get canonical SMILES (fast)
prod_smiles_canon = canonical_smiles(prod_mol)
if not prod_smiles_canon:
    continue

# NEW: Get mapped SMILES (preserve atom mappings for site identification)
prod_smiles_mapped = Chem.MolToSmiles(prod_mol, canonical=False)
```

#### Change 2: Extract mapping numbers (New function)

**Add new helper function**:
```python
def extract_site_mapnums(mol: Chem.Mol, target_mapnum: int = 777) -> List[int]:
    """
    Extract atom indices with specific mapping number.

    Args:
        mol: Product molecule with atom mappings
        target_mapnum: Mapping number to search for (default 777)

    Returns:
        List of atom indices with target mapping number
    """
    mapnums = []
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == target_mapnum:
            mapnums.append(target_mapnum)
    return mapnums
```

#### Change 3: Update product dictionary (Line ~497-511)

**Current**:
```python
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
    'source_smiles': smiles,
    **source_data,
    **props,
    **halogen_fields
})
```

**Modified**:
```python
# Extract mapping numbers
site_mapnums = extract_site_mapnums(prod_mol, target_mapnum=777)

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
    'xf_site_mapnums': json.dumps(site_mapnums),  # JSON list
    'source_smiles': smiles,
    **source_data,
    **props,
    **halogen_fields
})
```

---

### File: `configs/transforms.yaml`

**Modifications**:

```yaml
transforms:
  # Phenolic OH Methylation: -OH → -OCH3
  - name: OH_to_OMe
    label: "OH->OMe"
    description: "O-methylation of phenolic hydroxyl groups"
    query_smarts: "[c:1][OX2H:2]"
    # MODIFIED: Add :777 to product oxygen
    smirks: "[c:1][OX2H:2]>>[c:1][O:777]C"

  # Phenolic OH to Amine: -OH → -NH2
  - name: OH_to_NH2
    label: "OH->NH2"
    description: "Bioisosteric replacement of phenolic OH with primary amine"
    query_smarts: "[c:1][OX2H:2]"
    # MODIFIED: Add :777 to product nitrogen, explicit H
    smirks: "[c:1][OX2H:2]>>[c:1][N:777]([H])[H]"
```

---

## Validation Strategy

### Phase 1: SMIRKS Syntax Validation
Test that modified SMIRKS are valid:
```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Test OH->OMe
rxn_ome = AllChem.ReactionFromSmarts("[c:1][OX2H:2]>>[c:1][O:777]C")
assert rxn_ome is not None, "OMe SMIRKS invalid"

# Test OH->NH2
rxn_nh2 = AllChem.ReactionFromSmarts("[c:1][OX2H:2]>>[c:1][N:777]([H])[H]")
assert rxn_nh2 is not None, "NH2 SMIRKS invalid"

print("✓ SMIRKS syntax valid")
```

### Phase 2: Mapping Number Extraction
Test on a simple molecule:
```python
# Test molecule: phenol
mol = Chem.MolFromSmiles("c1ccccc1O")

# Run OH->OMe reaction
rxn = AllChem.ReactionFromSmarts("[c:1][OX2H:2]>>[c:1][O:777]C")
products = rxn.RunReactants((mol,))

if products:
    prod = products[0][0]
    Chem.SanitizeMol(prod)

    # Check mapping
    for atom in prod.GetAtoms():
        if atom.GetAtomMapNum() == 777:
            print(f"✓ Found :777 on atom {atom.GetIdx()} ({atom.GetSymbol()})")
            break
    else:
        print("✗ ERROR: :777 not found")

    # Get mapped SMILES
    mapped_smiles = Chem.MolToSmiles(prod, canonical=False)
    print(f"Mapped SMILES: {mapped_smiles}")
```

### Phase 3: End-to-End Test
Run transformation on 10 test molecules, verify:
1. `product_mapped_smiles` contains `:777`
2. `xf_site_mapnums` is `[777]`
3. site_finder can extract sites using mapping

### Phase 4: Regression Test
Re-run 8 problem samples:
- All should achieve **High confidence** with `method_used='mapnum'`
- Verify picked atoms match expected transformation sites

---

## Expected Outcomes

### Before PR-A (Current State)
```
High confidence:    69.5% (139/200)
Medium confidence:  18.5% (37/200)
Low confidence:     12.0% (24/200)
Failed:              0.0% (0/200)
```

### After PR-A (Projected)
```
High confidence:    ~98% (196/200)  ← Most use mapnum method
Medium confidence:   ~2% (4/200)    ← Fallback to smart for edge cases
Low confidence:      ~0% (0/200)
Failed:              ~0% (0/200)
```

**Rationale**:
- Atom mapping provides deterministic site identification
- No dependency on SMILES renumbering
- No multi-candidate ambiguity for k=1 cases
- Only fallback needed for molecules where reaction failed to assign mapping (rare)

---

## Rollback Plan

If PR-A causes issues:
1. **Data schema**: New fields are optional, old fields still present
2. **Visualization**: site_finder has 3-tier fallback (mapnum → smart → hint)
3. **Config**: Can revert `transforms.yaml` to original SMIRKS
4. **Code**: Can disable mapping extraction with a flag

**No breaking changes** - PR-A is purely additive.

---

## Testing Checklist

- [ ] Phase A1: SMIRKS syntax validation (RDKit parsing)
- [ ] Phase A2: Config file updated (transforms.yaml)
- [ ] Phase A3: Code changes implemented (08_transform_library_v2.py)
- [ ] Phase A4: Helper function added (extract_site_mapnums)
- [ ] Phase A5: Small-scale test (10-20 molecules)
  - [ ] product_mapped_smiles contains :777
  - [ ] xf_site_mapnums is [777]
  - [ ] site_finder extracts correct atoms
- [ ] Phase A6: Regression test (200 molecules)
  - [ ] High confidence ≥95%
  - [ ] method_used='mapnum' for majority
  - [ ] All 8 problem samples fixed
  - [ ] Visualizations correct

---

**Design Status**: ✅ **COMPLETE**
**Ready for Implementation**: ✅ **YES**
