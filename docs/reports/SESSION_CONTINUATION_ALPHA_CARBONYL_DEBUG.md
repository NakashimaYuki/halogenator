# Session Continuation: ALPHA_CARBONYL Bug Debug & k=1/k=2 Enumeration

**Date:** 2025-12-05
**Status:** ALPHA_CARBONYL bug investigation in progress
**Priority:** HIGH - Critical bug blocking complete k=1 enumeration

---

## Part 1: Completed Work Summary

### ✅ Phase 1-6: Classification & Rules Expansion (COMPLETED)

1. **Classification System**
   - Created master parent index (67,983 unique molecules)
   - Removed Flavone overlaps (7,180 molecules reassigned)
   - Exported clean `base_clean.parquet` files (mutually exclusive)
   - Files: `scripts/01_build_master_parent_index.py`, `NP_CLASSIFICATION_OVERLAP_REPORT.md`

2. **Rules Expansion**
   - Enabled R3 (hydroxyl) in all 7 classes
   - Enabled R4 (amine) in alkaloid/aa_peptide/other
   - Added ALPHA_CARBONYL to all classes k=1
   - Removed ALL max_sites_per_parent constraints
   - File: `configs/halogen_rules_by_class.yaml`

3. **CLI Fixes**
   - Added support for "smi" column name (`cli.py:1909-1912`)
   - Modified orchestrator to use base_clean.parquet (`scripts/04_enum_halogen_all_classes.py:71`)

### ✅ k=1 Enumeration Results (R1/R3/R4/R5 only)

**Total: 1,090,256 products from 56,708 parents**

| Class | Parents | Products | Key Rules |
|-------|---------|----------|-----------|
| lipid | 18 | 24 | R3 (83%), R5 (17%) |
| aa_peptide | 3,863 | 90,488 | R1 (48%), R4 (30%), R3 (17%) |
| polyphenol | 2,899 | 78,488 | R1 (66%), R3 (32%) |
| alkaloid | 4,202 | 88,396 | R1 (73%), R3 (13%), R4 (12%) |
| terpenoid | 25,513 | 284,508 | R3 (52%), R1 (43%) |
| glycoside | 20,213 | 548,352 | R3 (77%), R1 (21%) |

**Rule Contributions:**
- R3 (hydroxyl): 620,204 (57%) ← MAJOR SUCCESS
- R1 (aromatic): 418,532 (38%)
- R4 (amine): 38,196 (3.5%)
- R5 (carboxyl): 33,324 (3%)

**Missing:** ALPHA_CARBONYL produces ZERO products despite being in config!

### ⚠️ Critical Bug Discovered: ALPHA_CARBONYL Not Working

**Current Status:** ALPHA_CARBONYL rule is configured and loaded but produces NO products in actual enumeration.

**Files Modified:**
- `src/halogenator/rules.py:33` - Fixed SMIRKS with correct atom mapping
- `src/halogenator/enumerate_k.py:1414-1422` - Added diagnostic logging (not yet tested)

---

## Part 2: ALPHA_CARBONYL Bug Investigation Status

### What We Know (Confirmed Working)

1. **✓ SMIRKS is correct and works in isolation**
   ```python
   [CH2:1][CX3:2](=O)>>[CH:1]([Cl])[CX3:2](=O)
   ```
   - Tested on butanone: produces correct products
   - Tested on terpenoid molecule #6: produces 2 correct products

2. **✓ Rules configuration is correct**
   - ALPHA_CARBONYL in YAML for all 7 classes
   - Config loader correctly maps to 'ALPHA_CARBONYL__CH2__TO__X'
   - Passed to enumerate engine in rules tuple

3. **✓ Reactions dict is built correctly**
   - `build_reactions()` creates ALPHA_CARBONYL reactions for all 4 halogens (F, Cl, Br, I)
   - No errors during reaction building

4. **✓ Matching logic works**
   - `_find_reaction_matches()` correctly finds 2 matches on test molecule #6
   - Atom mapping detection works (finds atoms with map number 1)
   - No exceptions during matching

5. **✓ Isotope tagging works**
   - Site carbon correctly tagged with ISOTOPE_TAG=13
   - Reaction runs on tagged molecule
   - Tagged carbon found in product next to new halogen

6. **✓ NP molecules contain CH2-CO substructure**
   - 21% of terpenoid molecules have α-carbonyl sites
   - Test molecule #6 confirmed to have 2 CH2-CO groups

### What Fails (Confirmed Broken)

1. **✗ Actual enumeration produces ZERO ALPHA_CARBONYL products**
   - Ran full k=1 enumeration: 0 ALPHA_CARBONYL products in 1.09M total
   - Ran test on molecule #6 (has 2 CH2-CO): only R1 products, no ALPHA_CARBONYL
   - POC test on 10 terpenoids: only R1/R3/R5, no ALPHA_CARBONYL

### Investigation Findings

**Code Review Completed:**
- `_find_reaction_matches()` (line 1472-1535): Logic correct, tested working
- `_apply_reaction_rule()` (line 1812-2039): Isotope tagging path (line 1943-1993) looks correct
- Reaction rules list (line 1409-1412): ALPHA_CARBONYL included
- No obvious filtering or blocking conditions found

**Diagnostic Logging Added:**
- Added WARNING logs at line 1415-1416 (before _apply_reaction_rule call)
- Added WARNING logs at line 1421-1422 (after _apply_reaction_rule return)
- **NOT YET TESTED** - need to run enumerate with diagnostic

**Possible Root Causes (Unconfirmed):**
1. Exception silently caught in isotope tagging try-except (line 1989-1993)
2. Product filtered by post-guard or deduplication in _process_reaction_product
3. Issue with sugar_mask interfering despite sugar_mask=False for most classes
4. Reaction execution fails but error not propagated
5. Some conditional check we haven't identified yet

---

## Part 3: Next Steps - Detailed Action Plan

### 🔍 TASK 1: Complete ALPHA_CARBONYL Bug Diagnosis (CRITICAL)

**Objective:** Find exact point where ALPHA_CARBONYL execution fails

**Step 1.1: Run Diagnostic Test**
```bash
cd E:/Projects/halogenator
halogenator enum-parquet \
  --input-parquet data/output/nplike/test_mol6.parquet \
  --outdir data/output/nplike/test_mol6_diagnostic \
  --k 1 \
  --np-class terpenoid \
  2>&1 | tee diagnostic_log.txt
```

**Expected Output:**
```
[DIAGNOSTIC] Applying ALPHA_CARBONYL rule
[DIAGNOSTIC] ALPHA_CARBONYL returned X products
```

**Analysis:**
- If "[DIAGNOSTIC] Applying..." appears: ALPHA_CARBONYL is in loop ✓
- If "returned 0 products": Problem is in _apply_reaction_rule
- If "returned >0 products" but not in final output: Problem is in post-processing

**Step 1.2: Add More Detailed Logging**

If Step 1.1 shows "returned 0 products", add logging inside `_apply_reaction_rule()`:

```python
# Add at line 1868 (after _find_reaction_matches)
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAG] ALPHA matches_with_sites: {len(matches_with_sites)}")

# Add at line 1945 (in isotope tagging loop)
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAG] Processing match: site={site_atom_idx}")

# Add at line 1967 (after finding tagged site)
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAG] Tagged site in product: {tagged_site_in_product}")

# Add at line 1978 (after _process_reaction_product)
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAG] final_prod: {final_prod is not None}")
```

**Step 1.3: Check Exception Paths**

Add logging in exception handlers:

```python
# At line 1989 (isotope tagging exception)
except Exception as e:
    if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
        LOG.error(f"[DIAG] Isotope tagging failed: {type(e).__name__}: {e}")
    isotope_unavailable_occurred = True
    continue

# At line 1995 (top-level exception)
except Exception as e:
    if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
        LOG.error(f"[DIAG] Top-level failure: {type(e).__name__}: {e}")
```

**Step 1.4: Check _process_reaction_product**

If products reach _process_reaction_product but get filtered out, add logging there (line 2042):

```python
def _process_reaction_product(...):
    if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
        LOG.warning(f"[DIAG] _process_reaction_product called for ALPHA_CARBONYL")

    # ... existing code ...

    # Before final return
    if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
        LOG.warning(f"[DIAG] Returning: {result is not None}")
    return result
```

**Step 1.5: Reproduce Issue Minimally**

Create standalone test script:
```python
# test_alpha_carbonyl_enum.py
import sys
sys.path.insert(0, 'E:/Projects/halogenator/src')

from halogenator.enumerate_k import _apply_one_layer
from halogenator.rules import build_reactions
from halogenator.schema import EnumConfig
from rdkit import Chem

# Test molecule #6
test_smi = 'COc1cc2c(cc1OC)[C@@]1(C)CCC(=O)C(C)(C)[C@@H]1CC2=O'
mol = Chem.MolFromSmiles(test_smi)

# Build reactions
reactions = build_reactions()

# Configure with ALPHA_CARBONYL
cfg = EnumConfig(
    k_max=1,
    halogens=('Cl',),
    rules=('ALPHA_CARBONYL__CH2__TO__X',),
    sugar_cfg={'mode': 'off', 'apply_mask': False}
)

# Run enumeration
next_frontier, product_records, qa_stats = _apply_one_layer(
    mol, 0, [], set(), set(), None, reactions, cfg, None, None, None, None, None
)

print(f"Products: {len(product_records)}")
print(f"Frontier: {len(next_frontier)}")
print(f"QA stats: {qa_stats}")

for record in product_records:
    print(f"  Rule: {record['rule']}, SMILES: {record['smiles'][:60]}")
```

Run: `python test_alpha_carbonyl_enum.py`

**Expected:** Should produce 2 ALPHA_CARBONYL products (one for each CH2-CO site)

**If fails:** Problem is in core enumerate logic, not CLI/config layer

---

### 🔧 TASK 2: Fix ALPHA_CARBONYL Bug

**Once root cause identified, apply fix based on findings:**

**Scenario A: Exception Being Thrown**
- Fix the exception cause
- Or: Improve exception handling to not silently fail

**Scenario B: Filtering/Deduplication Issue**
- Adjust filtering logic for ALPHA_CARBONYL
- Or: Fix product structure that causes filter rejection

**Scenario C: SMIRKS Still Incorrect**
- Review product structure from manual test
- Adjust SMIRKS if needed
- Re-test manually before full enum

**Scenario D: Reactions Dict Not Passed Correctly**
- Check enumerate_k.py reaction parameter flow
- Verify reactions dict passed to _apply_reaction_rule contains ALPHA_CARBONYL

**After Fix:**
- Remove diagnostic logging (or set to DEBUG level)
- Test on molecule #6: should produce ALPHA_CARBONYL products
- Test on 10 terpenoids POC: should show ALPHA_CARBONYL in rule distribution

---

### 🚀 TASK 3: Re-run Complete k=1 Enumeration (With ALPHA_CARBONYL)

**Prerequisites:**
- ✅ ALPHA_CARBONYL bug fixed and verified
- ✅ Test on molecule #6 passes
- ✅ POC test (10-100 molecules) shows ALPHA_CARBONYL in output

**Step 3.1: Clean Old Outputs**
```bash
rm -rf E:/Projects/halogenator/data/output/nplike/*-1X
```

**Step 3.2: Run Full Enumeration**
```bash
cd E:/Projects/halogenator
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid terpenoid glycoside \
  --k-values 1 \
  2>&1 | tee logs/k1_full_with_alpha_carbonyl.log
```

**Expected Runtime:** ~60-90 minutes

**Expected Results:**
- Previous: 1,090,256 products (R1/R3/R4/R5 only)
- New: 1,200,000-1,400,000 products (with ALPHA_CARBONYL)
- ALPHA_CARBONYL should contribute ~10-15% additional products

**Step 3.3: Verify Results**
```python
import pandas as pd
from pathlib import Path

classes = ['lipid', 'aa_peptide', 'polyphenol', 'alkaloid', 'terpenoid', 'glycoside']
base_dir = Path('E:/Projects/halogenator/data/output/nplike')

total = 0
alpha_total = 0

for cls in classes:
    df = pd.read_parquet(base_dir / f'{cls}-1X' / 'products.parquet')
    rule_dist = df['rule'].value_counts()

    total += len(df)
    if 'ALPHA_CARBONYL__CH2__TO__X' in rule_dist.index:
        alpha_count = rule_dist['ALPHA_CARBONYL__CH2__TO__X']
        alpha_total += alpha_count
        pct = alpha_count / len(df) * 100
        print(f'{cls}: {alpha_count:,} ALPHA products ({pct:.1f}%)')
    else:
        print(f'{cls}: NO ALPHA_CARBONYL products [WARNING]')

print(f'\nTotal: {total:,} products')
print(f'ALPHA_CARBONYL: {alpha_total:,} ({alpha_total/total*100:.1f}%)')
```

**Success Criteria:**
- ✅ ALPHA_CARBONYL appears in at least 4 classes
- ✅ Total ALPHA products > 50,000
- ✅ Overall products increased by 10-30% vs previous
- ✅ No extreme_site warnings for ALPHA_CARBONYL

**If ALPHA_CARBONYL still shows 0:**
- Bug fix was insufficient, return to TASK 1
- Check if rules were correctly passed in new run
- Verify code changes were saved and environment reloaded

---

### 📊 TASK 4: Generate Final k=1 Completion Report

**Create comprehensive report:**

```markdown
# k=1 Enumeration Final Report

## Summary
- Total products: X,XXX,XXX
- Total parents: 56,708
- Products/parent: XX.X
- Active rules: R1, R3, R4, R5, ALPHA_CARBONYL

## Rule Distribution
| Rule | Products | Percentage | Notes |
|------|----------|------------|-------|
| R3 (hydroxyl) | XXX,XXX | XX% | PRIMARY CONTRIBUTOR |
| R1 (aromatic) | XXX,XXX | XX% | |
| ALPHA_CARBONYL | XXX,XXX | XX% | **NEW - α-carbonyl halogenation** |
| R4 (amine) | XXX,XXX | XX% | |
| R5 (carboxyl) | XXX,XXX | XX% | |

## ALPHA_CARBONYL Validation
- Expected contribution: 10-15%
- Actual contribution: XX%
- Status: [SUCCESS/PARTIAL/FAIL]
- Classes with ALPHA: [list]

## Comparison with Previous
| Metric | Without ALPHA | With ALPHA | Change |
|--------|---------------|------------|--------|
| Total products | 1,090,256 | X,XXX,XXX | +XX% |
| Rules active | 4 | 5 | +1 |
| terpenoid prod/parent | 11.2 | XX.X | +XX% |
| alkaloid prod/parent | 21.0 | XX.X | +XX% |

## Next Steps
- Proceed to k=2 enumeration
- Expected k=2 products: 6M-10M
```

---

### 🚀 TASK 5: k=2 Enumeration

**Prerequisites:**
- ✅ k=1 enumeration complete and validated
- ✅ All 5 rules (R1, R3, R4, R5, ALPHA_CARBONYL) working
- ✅ k=1 results meet expectations (>1.2M products)

**Step 5.1: Fast Batch (Est. 2-4 hours)**
```bash
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid \
  --k-values 2 \
  2>&1 | tee logs/k2_fast_batch.log
```

Expected:
- lipid: 18 → ~50-100 products
- aa_peptide: 3,863 → ~300K-600K products
- polyphenol: 2,899 → ~300K-600K products (high R3 impact)
- alkaloid: 4,202 → ~400K-700K products (R3+R4+ALPHA)

**Step 5.2: Slow Batch (Est. 6-10 hours)**
```bash
# Run terpenoid (largest class)
python scripts/04_enum_halogen_all_classes.py \
  --classes terpenoid \
  --k-values 2 \
  2>&1 | tee logs/k2_terpenoid.log

# Run glycoside (sugar_mask overhead)
python scripts/04_enum_halogen_all_classes.py \
  --classes glycoside \
  --k-values 2 \
  2>&1 | tee logs/k2_glycoside.log
```

Expected:
- terpenoid: 25,513 → ~3M-6M products
- glycoside: 20,213 → ~2M-4M products

**Total Expected k=2:** 6M-12M products

**Step 5.3: Validate k=2 Results**
```python
# Check extreme site guards
import json
for cls in classes:
    with open(f'data/output/nplike/{cls}-2X/SUMMARY.json') as f:
        summary = json.load(f)

    soft = summary.get('qa_paths', {}).get('extreme_site_soft_warning', 0)
    hard = summary.get('qa_paths', {}).get('extreme_site_hard_skip', 0)

    if hard > 0:
        print(f'{cls}: {hard} hard skips [WARNING - may need threshold adjustment]')
    elif soft > 5:
        print(f'{cls}: {soft} soft warnings [INFO]')
```

**Success Criteria:**
- ✅ All classes complete without errors
- ✅ extreme_site_hard_skip < 5% of attempts
- ✅ Products/parent reasonable (<200 for all classes)
- ✅ k=2 products = 5-10x k=1 products

---

### 🔄 TASK 6: Library Integration & QC

**Step 6.1: Merge k=1 + k=2**
```python
for cls in classes:
    k1 = pd.read_parquet(f'data/output/nplike/{cls}-1X/products.parquet')
    k2 = pd.read_parquet(f'data/output/nplike/{cls}-2X/products.parquet')

    combined = pd.concat([k1, k2], ignore_index=True)
    combined_dedup = combined.drop_duplicates(subset='inchikey', keep='first')

    outdir = Path(f'data/output/nplike/{cls}-full')
    outdir.mkdir(exist_ok=True)
    combined_dedup.to_parquet(outdir / 'products.parquet', index=False)

    print(f'{cls}: k1={len(k1):,}, k2={len(k2):,}, merged={len(combined_dedup):,}')
```

**Step 6.2: Cross-Class Deduplication Check**
```python
# Check for product overlap between classes (expected: some overlap)
all_inchikeys = {}
for cls in classes:
    df = pd.read_parquet(f'data/output/nplike/{cls}-full/products.parquet')
    all_inchikeys[cls] = set(df['inchikey'])

overlaps = []
for i, c1 in enumerate(classes):
    for c2 in classes[i+1:]:
        overlap = len(all_inchikeys[c1] & all_inchikeys[c2])
        if overlap > 100:
            overlaps.append((c1, c2, overlap))

print('Significant cross-class overlaps:')
for c1, c2, count in sorted(overlaps, key=lambda x: x[2], reverse=True):
    print(f'  {c1} ∩ {c2}: {count:,} products')
```

**Step 6.3: Generate Statistics**
```bash
python scripts/05_summaries.py \
  --input data/output/nplike/terpenoid-full/products.parquet \
  --output data/output/nplike/terpenoid-full/stats.json

# Repeat for all classes
```

---

### 📈 TASK 7: Descriptor Calculation & Visualization

**Step 7.1: Calculate Descriptors**
```bash
for cls in lipid aa_peptide polyphenol alkaloid terpenoid glycoside; do
  python scripts/08a_fill_descriptors.py \
    --input data/output/nplike/${cls}-full/products.parquet \
    --output data/output/nplike/${cls}-full/products_descriptors.parquet
done
```

**Descriptors to calculate:**
- MW, LogP, HBA, HBD, TPSA
- Aromatic rings, rotatable bonds
- Complexity metrics

**Step 7.2: Generate Visualizations**
```bash
python scripts/09_visualize_library.py html \
  -i data/output/nplike/polyphenol-full/products_descriptors.parquet \
  -o data/output/nplike/polyphenol-full/gallery.html \
  -n 500
```

---

### 📦 TASK 8: Export & Packaging

**Step 8.1: Export for Virtual Screening**
```bash
python scripts/11_export_for_vs.py \
  --input data/output/nplike/terpenoid-full/products_descriptors.parquet \
  --output data/output/vs_libraries/terpenoid_full.smi \
  --format smiles

# Also export SDF format
python scripts/11_export_for_vs.py \
  --input data/output/nplike/terpenoid-full/products_descriptors.parquet \
  --output data/output/vs_libraries/terpenoid_full.sdf \
  --format sdf
```

**Step 8.2: Create Final Package**
```
halogenated_np_library_v1.0/
├── README.md
├── CHANGELOG.md
├── libraries/
│   ├── by_class/
│   │   ├── lipid/
│   │   ├── aa_peptide/
│   │   └── ... (k1, k2, full for each)
│   └── merged/
│       ├── all_k1.parquet
│       ├── all_k2.parquet
│       └── all_full.parquet
├── vs_formats/
│   ├── all_full.smi
│   └── all_full.sdf
└── metadata/
    ├── K1_ENUMERATION_FINAL_REPORT.md
    ├── master_parent_index.parquet
    └── classification_report.md
```

---

## Part 4: Key Context & Reference

### Current File Status

**Data Files:**
- `data/output/nplike/*/base_clean.parquet` - Clean parent libraries (67,983 total)
- `data/output/nplike/*-1X/products.parquet` - k=1 results (1.09M, NO ALPHA_CARBONYL)
- `data/output/nplike/test_mol6.parquet` - Test molecule for ALPHA bug diagnosis

**Source Code:**
- `src/halogenator/rules.py:31-35` - ALPHA_CARBONYL SMIRKS (FIXED with :2 mapping)
- `src/halogenator/enumerate_k.py:1414-1422` - Diagnostic logging (ADDED, not tested)
- `src/halogenator/cli.py:1909-1912` - smi column support (FIXED)
- `src/halogenator/class_config.py` - Config loader (working correctly)

**Configuration:**
- `configs/halogen_rules_by_class.yaml` - All rules enabled, zero max_sites

**Reports:**
- `K1_ENUMERATION_COMPLETION_REPORT.md` - Current k=1 results (without ALPHA)
- `NP_CLASSIFICATION_OVERLAP_REPORT.md` - Overlap analysis
- `NEXT_SESSION_CONTINUATION_GUIDE.md` - Original continuation guide

### Critical Code Locations

**ALPHA_CARBONYL SMIRKS (Fixed):**
```python
# src/halogenator/rules.py:31-35
def _alpha_carbonyl_smirks(X: str) -> str:
    """Alpha to carbonyl CH2 -> CHX (e.g., HVZ-like)"""
    # X must be element symbol (F, Cl, Br, I)
    # Map both CH2 (:1) and carbonyl carbon (:2) to preserve connectivity
    return "[CH2:1][CX3:2](=O)>>[CH:1]([" + X + "])[CX3:2](=O)"
```

**Reaction Rules List:**
```python
# src/halogenator/enumerate_k.py:1409-1412
reaction_rules = [r for r in cfg.rules if r in ('R3', 'R4', 'R5',
                                                'RING_SP3__CH__TO__X',
                                                'ALPHA_CARBONYL__CH2__TO__X',
                                                'PRIMARY_OH__CH2__TO__X')]
```

**Diagnostic Logging Added:**
```python
# src/halogenator/enumerate_k.py:1414-1422
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAGNOSTIC] Applying ALPHA_CARBONYL rule")
reaction_products, rule_qa_stats = _apply_reaction_rule(...)
if rule_id == 'ALPHA_CARBONYL__CH2__TO__X':
    LOG.warning(f"[DIAGNOSTIC] ALPHA_CARBONYL returned {len(reaction_products)} products")
```

### Test Commands

**Quick ALPHA Test:**
```bash
cd E:/Projects/halogenator
halogenator enum-parquet \
  --input-parquet data/output/nplike/test_mol6.parquet \
  --outdir data/output/nplike/test_alpha_debug \
  --k 1 --np-class terpenoid \
  2>&1 | grep -E "DIAGNOSTIC|ALPHA"
```

**Check Results:**
```python
import pandas as pd
df = pd.read_parquet('data/output/nplike/test_alpha_debug/products.parquet')
print(df['rule'].value_counts())
```

### Known Working Tests

These tests ALL pass outside of enumerate:

1. **SMIRKS parsing:** ✓ `RXN.ReactionFromSmarts('[CH2:1][CX3:2](=O)>>[CH:1]([Cl])[CX3:2](=O)')`
2. **Reaction on butanone:** ✓ Produces `CC(=O)C(C)Cl`
3. **Matching logic:** ✓ Finds 2 matches on molecule #6
4. **Isotope tagging:** ✓ Tags site carbon, finds it in product
5. **Config loading:** ✓ ALPHA_CARBONYL in final rules list
6. **Reactions dict:** ✓ ALPHA_CARBONYL reactions exist for F/Cl/Br/I

**The bug is ONLY in the enumerate execution path, not in any individual component.**

---

## Part 5: Success Metrics

### k=1 with ALPHA_CARBONYL (Target)
- ✅ Total products: 1.2M-1.4M
- ✅ ALPHA_CARBONYL: 50K-200K products (5-15%)
- ✅ All 5 rules active in distribution
- ✅ Products/parent: 20-25 average

### k=2 (Target)
- ✅ Total products: 6M-12M
- ✅ k=2/k=1 ratio: 5-10x
- ✅ extreme_site_hard_skip: <5% attempts
- ✅ All classes complete successfully

### Final Library (Target)
- ✅ Combined k=1+k=2: 7M-13M unique products
- ✅ All descriptors calculated
- ✅ VS formats exported (SMILES, SDF)
- ✅ Complete documentation package

---

## Part 6: Troubleshooting Guide

**If ALPHA_CARBONYL still produces 0 products after "fix":**
1. Check diagnostic logs show "[DIAGNOSTIC] Applying ALPHA_CARBONYL rule"
2. If yes but "returned 0 products": Add more detailed logging in _apply_reaction_rule
3. If no diagnostic output: Rules not being iterated, check reaction_rules list construction
4. Run standalone test script (test_alpha_carbonyl_enum.py) to isolate issue

**If k=2 enumeration fails with extreme_site warnings:**
1. Check which class/rule triggering warnings
2. If >10% hard skips: Consider raising EXTREME_SITE_HARD_THRESHOLD to 80
3. If specific molecules problematic: May need to filter them from base

**If descriptor calculation fails:**
1. Check for invalid SMILES in products
2. Run sanitization pass before descriptor calculation
3. Skip molecules that fail, log them separately

---

**Document Status:** ✅ Ready for next session
**Last Updated:** 2025-12-05
**Priority:** Start with TASK 1 (ALPHA_CARBONYL debug) immediately
**Estimated Time to Complete All Tasks:** 15-20 hours (mostly enumeration runtime)
