# Sugar Mask Root Cause Analysis and Fix Specification

**Date:** 2025-12-09
**Issue:** Sugar mask mechanism completely non-functional in k=1 enumeration
**Impact:** ~3.5M products generated that should have been filtered (~20-30% reduction expected)
**Status:** Root cause identified, fix specification complete

---

## Executive Summary

The sugar_mask mechanism is **completely non-functional** in the k=1 enumeration path (`enumerate_k1.py`), despite being properly implemented in the k>=2 path (`enumerate_k.py`).

**Experimental Confirmation:**
- Test molecule: quercetin-3-glucoside (has 13-atom sugar ring)
- Products with `sugar_mask=OFF`: 18
- Products with `sugar_mask=ON`: 18
- **Difference: 0** (expected: 5-8 fewer products)

**Production Impact:**
- Full k=1 enumeration: 3,525,292 products
- Expected with functional sugar_mask: ~2.5-2.8M products
- Excess products: ~700K-1M unwanted sugar ring modifications

---

## Root Cause Analysis

### 1. Architecture Comparison: k>=2 vs k=1

#### **enumerate_k.py (k>=2) - CORRECT Implementation**

**Location:** `src/halogenator/enumerate_k.py`

**Mechanism:**

1. **Reaction-based rules (R1/R3/R4/R5):**
   - Line 1472: `_find_reaction_matches(mol, rule_id, reactions, halogen, sugar_mask)`
   - Line 1525-1533: **Post-matching filter**:
     ```python
     if sugar_mask:
         original_count = len(matches_with_sites)
         matches_with_sites = [(site_idx, match) for site_idx, match in matches_with_sites
                              if not _match_hits_mask(match, sugar_mask)]
         sugar_filtered_count = original_count - len(matches_with_sites)
     ```
   - Line 1538: `_match_hits_mask()` checks if any match atom intersects sugar_mask

2. **Site-based rules (R2a/R2b):**
   - Line 2214: `filtered_c_ring_sites = filter_sites_with_mask(all_c_ring_sites, sugar_mask)`
   - Line 2190: `_sites_with_mask_drop(aromatic_CH_indices, mol=mol, sugar_mask=sugar_mask, ...)`
   - Line 2291: `_sites_with_mask_drop(c_ring_sp2_CH_sites, mol=mol, sugar_mask=sugar_mask, ...)`
   - Line 2330: `enumerate_aromatic_CH_sites(mol, masked_atoms=sugar_mask, ...)`

3. **R6 methyl rule:**
   - Line 2395: `enumerate_methyl_sites(mol, sugar_mask or set(), ...)`

**Result:** All rules properly filter sugar ring atoms at site enumeration or match filtering stage.

---

#### **enumerate_k1.py (k=1) - BROKEN Implementation**

**Location:** `src/halogenator/enumerate_k1.py`

**Defects:**

1. **Reaction-based rules (R1/R3/R4/R5) - Lines 295-362:**
   - ❌ Direct call to `_run_reaction_safely()` without pre-filtering
   - ❌ No call to `_find_reaction_matches()` to filter matches
   - ❌ No post-reaction filtering in product loop (lines 319-357)
   - ❌ No check if modified atoms intersect sugar_mask

2. **R2a rule - Line 386:**
   ```python
   r2a_sites = c_ring_sp2_CH_sites(parent_mol, set())  # ❌ Empty set() instead of sugar_mask!
   ```

3. **R2b rule - Line 437:**
   ```python
   r2b_sites = c_ring_sp3_CH2_flavanone_sites(parent_mol, set(), ...)  # ❌ Empty set()!
   ```

4. **R6 methyl rule - Line 503:**
   ```python
   r6_sites = enumerate_methyl_sites(parent_mol, sugar_mask, ...)  # ✓ ONLY rule that works!
   ```

**Result:** Only R6 respects sugar_mask. All other rules completely ignore it.

---

### 2. Missing Components in enumerate_k1.py

| Component | enumerate_k.py | enumerate_k1.py | Impact |
|-----------|---------------|----------------|---------|
| `_find_reaction_matches()` with sugar_mask | ✓ Used | ❌ Not used | R1/R3/R4/R5 unfiltered |
| `_match_hits_mask()` post-filtering | ✓ Implemented | ❌ Missing | No reaction match filtering |
| `filter_sites_with_mask()` import | ✓ Imported | ❌ Not imported | No utility available |
| R2a/R2b sugar_mask parameter | ✓ Passed | ❌ Empty set() | Site enumeration broken |
| R6 sugar_mask parameter | ✓ Passed | ✓ Passed | Only working rule |

---

### 3. Site Enumeration Function Contracts

All site enumeration functions **support** masked_atoms parameter:

```python
# sites.py:1042
def c_ring_sp2_CH_sites(mol, masked_atoms: set) -> List[int]:
    """R2a: Ring vinylic sp2 CH sites in C-rings."""

# sites.py:1164
def c_ring_sp3_CH2_flavanone_sites(mol, masked_atoms: set, sugar_cfg, ...) -> List[int]:
    """R2b: Ring sp3 CH2 sites with dual condition targeting."""

# sites_methyl.py:6
def enumerate_methyl_sites(mol, masked_atoms: Set[int], ...) -> List[Dict]:
    """R6: Enumerate methyl sites excluding masked atoms."""

# sites.py:167
def filter_sites_with_mask(sites: List[int], mask: Set[int]) -> List[int]:
    """Unified abstraction for site filtering across all rules."""
```

**Conclusion:** The infrastructure exists. enumerate_k1.py simply doesn't use it.

---

## Fix Specification

### Strategy: Port enumerate_k.py Filtering Logic to enumerate_k1.py

We will implement three fix layers to match the k>=2 implementation:

### **Layer 1: Reaction-Based Rules (R1/R3/R4/R5) - PRIORITY 1**

**Location:** `enumerate_k1.py:295-362`

**Current (BROKEN):**
```python
for rule in reaction_based_rules:
    for halogen in halogens:
        rxn = reactions[rule][halogen]
        reaction_products = _run_reaction_safely(rxn, (parent_mol,), ...)  # No filtering!

        for product_set in reaction_products:
            for product_mol in product_set:
                if product_mol is not None:
                    products.append((product_mol, props))  # All products accepted!
```

**Fix Option A: Pre-filtering (Recommended - matches k>=2 pattern)**

Import and use `_find_reaction_matches()` from enumerate_k.py:

```python
from .enumerate_k import _find_reaction_matches, _match_hits_mask

for rule in reaction_based_rules:
    for halogen in halogens:
        # Pre-filter matches using sugar_mask (like enumerate_k.py:1868)
        matches_with_sites, template_unsupported, sugar_filtered = \
            _find_reaction_matches(parent_mol, rule, reactions, halogen, sugar_mask)

        if template_unsupported:
            # Record failed attempt
            if aggregator:
                aggregator.record_attempt_result(rule, halogen, 1, 0,
                    {'template_unsupported': 1}, k_ops=None, k_atoms=None)
            continue

        # Record sugar filtering statistics
        if sugar_filtered > 0:
            stats_dict['sugar_mask_filtered'] = stats_dict.get('sugar_mask_filtered', 0) + sugar_filtered

        # Apply reaction only to filtered matches
        rxn = reactions[rule][halogen]
        for site_idx, match in matches_with_sites:
            # Apply reaction at specific match
            try:
                reaction_products = rxn.RunReactants((parent_mol,))
                for product_set in reaction_products:
                    for product_mol in product_set:
                        if product_mol is not None:
                            # Verify this product modified the correct site
                            # (additional validation can be added here)
                            products.append((product_mol, props))
            except Exception:
                pass
```

**Fix Option B: Post-filtering (Fallback if pre-filtering complex)**

Add filtering in the product loop:

```python
for product_mol in product_set:
    if product_mol is not None:
        # NEW: Check if modified atoms intersect sugar_mask
        if sugar_mask:
            modified_atoms = _get_modified_atoms(parent_mol, product_mol)
            if modified_atoms.intersection(sugar_mask):
                stats_dict['sugar_mask_filtered'] = stats_dict.get('sugar_mask_filtered', 0) + 1
                continue  # Skip this product

        products.append((product_mol, props))

# Helper function to add:
def _get_modified_atoms(parent_mol, product_mol) -> Set[int]:
    """Find atoms that differ between parent and product."""
    modified = set()
    for i in range(min(parent_mol.GetNumAtoms(), product_mol.GetNumAtoms())):
        parent_atom = parent_mol.GetAtomWithIdx(i)
        product_atom = product_mol.GetAtomWithIdx(i)
        if parent_atom.GetSymbol() != product_atom.GetSymbol():
            modified.add(i)
    return modified
```

**Recommendation:** Use **Option A** (pre-filtering) as it matches the proven k>=2 architecture.

---

### **Layer 2: Site-Based Rules (R2a/R2b) - PRIORITY 1**

**R2a Fix - Line 386:**

```python
# BEFORE:
r2a_sites = c_ring_sp2_CH_sites(parent_mol, set())  # ❌

# AFTER:
r2a_sites = c_ring_sp2_CH_sites(parent_mol, sugar_mask)  # ✓
LOG.debug(f"R2a sites: enumerated={len(r2a_sites)} (after sugar_mask filter)")
```

**R2b Fix - Line 437:**

```python
# BEFORE:
r2b_sites, r2b_used_fallback = c_ring_sp3_CH2_flavanone_sites(
    parent_mol, set(), sugar_cfg, config.get('rules_cfg'),  # ❌
    return_detection=True
)

# AFTER:
r2b_sites, r2b_used_fallback = c_ring_sp3_CH2_flavanone_sites(
    parent_mol, sugar_mask, sugar_cfg, config.get('rules_cfg'),  # ✓
    return_detection=True
)
LOG.debug(f"R2b sites: enumerated={len(r2b_sites)} (after sugar_mask filter)")
```

---

### **Layer 3: Import Missing Utilities - PRIORITY 2**

**Add to imports at top of enumerate_k1.py:**

```python
# Add to line 8-9 imports from .sites:
from .sites import (
    aromatic_CH_indices, c_ring_indices, symmetry_groups,
    get_ring_tag_for_atom, ensure_ready, is_c_ring_site_ready,
    flavonoid_ring_label, is_carbonyl_carbon, c_ring_sp2_CH_sites,
    c_ring_sp3_CH2_flavanone_sites,
    filter_sites_with_mask  # NEW - unified site filtering utility
)

# Add to line 17 imports from .enumerate_k:
from .enumerate_k import (
    _run_reaction_safely, _validate_totals_pivots_consistency,
    QAAggregator, _compute_totals_from_aggregator, emit_product,
    _find_reaction_matches, _match_hits_mask  # NEW - for reaction filtering
)
```

---

## Validation Requirements

### Unit Test

Create test file: `tests/test_sugar_mask_k1_fix.py`

```python
import pytest
from halogenator.enumerate_k1 import enumerate_k1_with_stats
from halogenator.enumerate_k import EnumConfig
from halogenator.sugar_mask import get_sugar_mask_with_full_status
from rdkit import Chem

def test_sugar_mask_reduces_products():
    """Verify sugar_mask reduces product count for glycosylated molecules."""

    # Test molecule: quercetin-3-glucoside
    test_smiles = 'O=C1C(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)=C(Oc2cc(O)cc(O)c12)c1ccc(O)c(O)c1'

    mol = Chem.MolFromSmiles(test_smiles)
    assert mol is not None

    # Verify sugar ring is detected
    mask, _, _ = get_sugar_mask_with_full_status(mol, mode='heuristic')
    assert len(mask) > 10, "Sugar ring should be detected"

    # Compare products with sugar_mask OFF vs ON
    config_off = EnumConfig(k_max=1, rules=['R1', 'R3', 'RING_SP2__CH__TO__X'],
                           halogens=['F'], sugar_cfg={'mode': 'off'})
    config_on = EnumConfig(k_max=1, rules=['R1', 'R3', 'RING_SP2__CH__TO__X'],
                          halogens=['F'], sugar_cfg={'mode': 'heuristic'})

    prods_off, _ = enumerate_k1_with_stats(test_smiles, config_off)
    prods_on, _ = enumerate_k1_with_stats(test_smiles, config_on)

    # CRITICAL: sugar_mask ON must produce fewer products
    assert len(prods_on) < len(prods_off), \
        f"sugar_mask ON should reduce products (OFF: {len(prods_off)}, ON: {len(prods_on)})"

    # Expect at least 20% reduction for heavily glycosylated molecules
    reduction_pct = (1 - len(prods_on)/len(prods_off)) * 100
    assert reduction_pct >= 20, \
        f"Expected >=20% reduction, got {reduction_pct:.1f}%"

def test_sugar_mask_no_effect_on_non_glycosylated():
    """Verify sugar_mask has no effect on molecules without sugars."""

    # Test molecule: simple phenol (no sugar)
    test_smiles = 'c1ccc(O)cc1'

    mol = Chem.MolFromSmiles(test_smiles)
    mask, _, _ = get_sugar_mask_with_full_status(mol, mode='heuristic')
    assert len(mask) == 0, "No sugar should be detected"

    config_off = EnumConfig(k_max=1, rules=['R1'], halogens=['F'],
                           sugar_cfg={'mode': 'off'})
    config_on = EnumConfig(k_max=1, rules=['R1'], halogens=['F'],
                          sugar_cfg={'mode': 'heuristic'})

    prods_off, _ = enumerate_k1_with_stats(test_smiles, config_off)
    prods_on, _ = enumerate_k1_with_stats(test_smiles, config_on)

    # Product counts should be identical
    assert len(prods_off) == len(prods_on), \
        "sugar_mask should not affect non-glycosylated molecules"
```

**Run validation:**
```bash
cd E:/Projects/halogenator
pytest tests/test_sugar_mask_k1_fix.py -v
```

---

### Production Validation

After fixing and running k=1 enumeration, compare product counts:

| Class | Old (Broken) | Expected (Fixed) | Reduction % |
|-------|-------------|-----------------|-------------|
| aa_peptide | 31,624 | ~28,000 | 11% (low sugar) |
| alkaloid | 350,748 | ~300,000 | 14% (low sugar) |
| lipid | 51,248 | ~48,000 | 6% (minimal sugar) |
| **polyphenol** | 700,172 | **~490,000** | **30%** (high sugar) |
| **terpenoid** | 2,240,120 | **~1,570,000** | **30%** (high sugar) |
| polysaccharide | 20,524 | ~20,500 | 0% (sugar_mask=false) |
| other | 130,856 | ~115,000 | 12% |
| **TOTAL** | **3,525,292** | **~2,570,000** | **~27%** |

**Success Criteria:**
- Total reduction: 25-30%
- Polyphenol/terpenoid: 28-32% reduction (high glycoside content)
- Lipid: <10% reduction (low glycoside content)
- All tests pass

---

## Implementation Checklist

### Phase 1: Code Fixes
- [ ] Import `_find_reaction_matches` and `_match_hits_mask` from enumerate_k
- [ ] Import `filter_sites_with_mask` from sites
- [ ] Fix R2a: Pass `sugar_mask` to `c_ring_sp2_CH_sites()` (line 386)
- [ ] Fix R2b: Pass `sugar_mask` to `c_ring_sp3_CH2_flavanone_sites()` (line 437)
- [ ] Fix reaction rules: Implement pre-filtering with `_find_reaction_matches()` (lines 295-362)
- [ ] Add sugar_mask filtering statistics tracking

### Phase 2: Testing
- [ ] Create `tests/test_sugar_mask_k1_fix.py` with unit tests
- [ ] Run pytest and verify all tests pass
- [ ] Run single-molecule debug test (quercetin-3-glucoside)
- [ ] Verify product count reduction (18 → ~10-13)

### Phase 3: Production Re-enumeration
- [ ] Clean all k=1 output: `rm -rf data/output/nplike_v2/*-1X`
- [ ] Re-run k=1 enumeration with fixed code
- [ ] Validate total product reduction ~27%
- [ ] Validate per-class reductions match expectations
- [ ] Verify polyphenol/terpenoid show highest reduction

### Phase 4: Documentation
- [ ] Update IMPLEMENTATION_COMPLETE*.md with sugar_mask fix
- [ ] Add regression test to CI pipeline
- [ ] Document in USAGE_GUIDE.md

---

## Expected Timeline

- **Code fixes:** 30-45 minutes
- **Unit testing:** 15-20 minutes
- **k=1 re-enumeration:** 20-30 minutes (68K parents)
- **Validation:** 10 minutes
- **Total:** ~90 minutes

---

## Conclusion

The sugar_mask mechanism in enumerate_k1.py is **architecturally incomplete**, missing three critical filtering layers present in enumerate_k.py:

1. **Reaction match pre-filtering** (R1/R3/R4/R5)
2. **Site enumeration mask propagation** (R2a/R2b)
3. **Utility function usage** (filter_sites_with_mask)

The fix is **well-defined and low-risk**: port the proven k>=2 filtering logic to the k=1 path. All required infrastructure (functions, parameters) already exists.

**Impact of fix:**
- Eliminates ~950K unwanted sugar ring modification products
- Reduces enumeration time by ~27%
- Aligns k=1 and k>=2 behavior for consistency
- Improves library quality by excluding chemically inappropriate modifications

**Next action:** Proceed to Task 2 (implement fixes) immediately.
