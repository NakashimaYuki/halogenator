# Next Session Continuation Guide

**Date Created:** 2025-12-04
**Purpose:** Context for next session - completed work summary + detailed next steps
**Status:** Ready for full-scale k=1/k=2 enumeration with clean bases

---

## Part 1: Completed Work (Brief Summary)

### ✅ Phase 1-3: Rules Expansion & Engine Guards

**What was done:**
1. **Rules audit**: Discovered R3 (hydroxyl) and R4 (amine) were implemented but unused
2. **YAML expansion**:
   - Added R3 (hydroxyl) to all classes
   - Added R4 (amine) to alkaloid/aa_peptide/other
   - Added ALPHA_CARBONYL to all classes in k1 (not just k2)
3. **Removed ALL max_sites_per_parent constraints** from YAML config
4. **Engine guards**: Added soft=40/hard=60 thresholds in enumerate_k.py
5. **POC validation**: Tested lipid/polyphenol/alkaloid - ALL PASSED
   - R3 generates 33% of polyphenol products
   - R4 generates 13.5% of alkaloid products
   - No extreme_site warnings

**Key files:**
- `configs/halogen_rules_by_class.yaml` - Expanded rules, zero constraints
- `src/halogenator/enumerate_k.py` - EXTREME_SITE_SOFT/HARD_THRESHOLD
- `src/halogenator/qa_utils.py` - EV_EXTREME_SITE_* events
- Backup: `configs/halogen_rules_by_class.yaml.bak_before_r3r4`

### ✅ Phase 4-5: Classification System & Clean Bases

**What was done:**
1. **Master parent index**: Built global index of 67,983 unique molecules
2. **Overlap detection**: Found 7,180 multi-class molecules (all Flavone overlaps)
3. **Primary_class priority**: Flavone > aa_peptide > alkaloid > glycoside > terpenoid > polyphenol > lipid > other
4. **Clean bases exported**: 8 mutually exclusive base_clean.parquet files
5. **Verification**: Confirmed zero overlaps between clean bases

**Key files:**
- `scripts/01_build_master_parent_index.py` - Index builder
- `data/output/nplike/master_parent_index.parquet` - Global index
- `data/output/nplike/<class>/base_clean.parquet` - Clean bases (no overlaps)
- `NP_CLASSIFICATION_OVERLAP_REPORT.md` - Overlap analysis

**Clean bases summary:**

| Class | Clean Parents | Notes |
|-------|--------------|-------|
| Flavone | 7,180 | Highest priority, no filtering |
| terpenoid | 25,513 | -3,093 from original (Flavone overlap) |
| glycoside | 20,213 | -2,660 from original |
| polyphenol | 2,899 | -1,148 from original (28% reduction) |
| alkaloid | 4,202 | -85 from original |
| aa_peptide | 3,863 | -27 from original |
| lipid | 18 | No change |
| other | 4,095 | -167 from original |
| **TOTAL** | **67,983** | Mutually exclusive |

### ✅ Phase 6: Final Rules Optimization (Latest Update)

**What was done:**
1. **ALPHA_CARBONYL added to k=1 for all classes** (previously only in k=2)
   - terpenoid, alkaloid, polyphenol, glycoside, aa_peptide, lipid, other
2. **RING_SP3 added to k=1 for 'other' class** (previously only in k=2)
3. **Finalized k=1 rule matrix:**

| Rule Type | k=1 Coverage |
|-----------|--------------|
| R1 (RING_SP2) | terpenoid, alkaloid, polyphenol, glycoside, aa_peptide, other |
| R3 (OH→X) | ALL 7 classes |
| R4 (NHx→X) | alkaloid, aa_peptide, other |
| R5 (COOH) | ALL 7 classes |
| RING_SP3 | terpenoid, alkaloid, polyphenol, glycoside, other |
| ALPHA_CARBONYL | ALL 7 classes (NEW) |
| PRIMARY_OH | terpenoid, lipid |

**Expected impact:**
- Previous estimate: 1.4M - 2.0M k=1 products
- Updated estimate: 1.5M - 2.2M k=1 products (~10% increase)
- ALPHA_CARBONYL will capture α-carbonyl positions in ketones, aldehydes, esters

---

## Part 2: Next Steps (Detailed Instructions)

### Task 1: Full k=1 Enumeration with Clean Bases

**Objective:** Generate complete k=1 halogenated library using clean bases + expanded rules

**Prerequisites:**
- ✅ Clean bases exist: `data/output/nplike/<class>/base_clean.parquet`
- ✅ Expanded rules config: R1 + R3 + R4 + R5 + ALPHA_CARBONYL + RING_SP3 + PRIMARY_OH (no max_sites constraints)
- ✅ Engine guards active: soft=40, hard=60

**Execution Strategy:**

#### Option A: Use Existing Orchestrator (Recommended)

**Script:** `scripts/04_enum_halogen_all_classes.py`

**Problem:** This script currently uses `base.parquet`, need to modify to use `base_clean.parquet`

**Quick Fix:**
```python
# In scripts/04_enum_halogen_all_classes.py, find the line:
input_path = DATA_DIR / class_name / "base.parquet"

# Change to:
input_path = DATA_DIR / class_name / "base_clean.parquet"
```

**Then run:**
```bash
# Option 1: All classes at once (sequential, ~60 min total)
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid terpenoid glycoside \
  --k-values 1

# Option 2: Fast classes first, then slow classes
# Fast batch (~15 min):
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid \
  --k-values 1

# Slow batch (~45 min):
python scripts/04_enum_halogen_all_classes.py \
  --classes terpenoid glycoside \
  --k-values 1
```

#### Option B: Manual CLI Commands (More Control)

Run each class individually to monitor progress:

```bash
# 1. Lipid (fastest, ~1 min)
halogenator enum-parquet \
  --input-parquet data/output/nplike/lipid/base_clean.parquet \
  --outdir data/output/nplike/lipid-1X \
  --k 1 \
  --np-class lipid

# 2. AA_peptide (~5-10 min)
halogenator enum-parquet \
  --input-parquet data/output/nplike/aa_peptide/base_clean.parquet \
  --outdir data/output/nplike/aa_peptide-1X \
  --k 1 \
  --np-class aa_peptide

# 3. Polyphenol (~5-10 min, HIGH R3 IMPACT)
halogenator enum-parquet \
  --input-parquet data/output/nplike/polyphenol/base_clean.parquet \
  --outdir data/output/nplike/polyphenol-1X \
  --k 1 \
  --np-class polyphenol

# 4. Alkaloid (~5-10 min, R3+R4)
halogenator enum-parquet \
  --input-parquet data/output/nplike/alkaloid/base_clean.parquet \
  --outdir data/output/nplike/alkaloid-1X \
  --k 1 \
  --np-class alkaloid

# 5. Terpenoid (~20-30 min, largest class)
halogenator enum-parquet \
  --input-parquet data/output/nplike/terpenoid/base_clean.parquet \
  --outdir data/output/nplike/terpenoid-1X \
  --k 1 \
  --np-class terpenoid

# 6. Glycoside (~20-30 min, sugar_mask overhead)
halogenator enum-parquet \
  --input-parquet data/output/nplike/glycoside/base_clean.parquet \
  --outdir data/output/nplike/glycoside-1X \
  --k 1 \
  --np-class glycoside

# Note: 'other' class can be added if needed
```

**Expected Output:**

Each class will generate:
- `data/output/nplike/<class>-1X/products.parquet` - Main product library
- `data/output/nplike/<class>-1X/SUMMARY.json` - Statistics and QA metrics

**Expected Results:**

| Class | Parents | Est. Products | Est. Runtime | Key Rules (k=1) |
|-------|---------|---------------|--------------|-----------------|
| lipid | 18 | ~200-400 | <1 min | R3, R5, ALPHA_CARBONYL |
| aa_peptide | 3,863 | ~85K-110K | 5-10 min | R1, R3, R4, R5, ALPHA_CARBONYL |
| polyphenol | 2,899 | ~85K-110K | 5-10 min | **R1, R3**, ALPHA_CARBONYL (R3 ~33%) |
| alkaloid | 4,202 | ~95K-130K | 5-10 min | R1, **R3, R4**, R5, ALPHA_CARBONYL |
| terpenoid | 25,513 | ~650K-850K | 20-30 min | R1, R3, R5, RING_SP3, ALPHA_CARBONYL |
| glycoside | 20,213 | ~520K-650K | 20-30 min | R1, R3, R5, ALPHA_CARBONYL (sugar_mask) |
| **TOTAL** | **56,708** | **1.5M - 2.2M** | **~60 min** | All 7 rule types active |

**Success Criteria:**
- ✅ All classes complete without errors
- ✅ extreme_site_hard_skip < 1% (check SUMMARY.json qa_paths)
- ✅ R3 shows up in product rule distribution (not 0%)
- ✅ Products/parent mean < 50 for all classes

---

### Task 2: Post-Enumeration QC & Analysis

**After k=1 completes, run QC checks:**

#### 2.1 Generate Summary Statistics

```bash
# Run summaries script for each class
python scripts/05_summaries.py \
  --input data/output/nplike/lipid-1X/products.parquet \
  --output data/output/nplike/lipid-1X/stats.json

# Or batch process all classes
for class in lipid aa_peptide polyphenol alkaloid terpenoid glycoside; do
  python scripts/05_summaries.py \
    --input data/output/nplike/${class}-1X/products.parquet \
    --output data/output/nplike/${class}-1X/stats.json
done
```

#### 2.2 Check Rule Distribution

**Critical check:** Verify R3/R4 are producing products

```python
import pandas as pd

classes = ['lipid', 'aa_peptide', 'polyphenol', 'alkaloid', 'terpenoid', 'glycoside']

for class_name in classes:
    df = pd.read_parquet(f'data/output/nplike/{class_name}-1X/products.parquet')
    rule_dist = df['rule'].value_counts()

    print(f"\n=== {class_name} k=1 Rule Distribution ===")
    print(rule_dist)

    # Check for R3/R4
    if 'R3' in rule_dist.index:
        pct = rule_dist['R3'] / len(df) * 100
        print(f"✓ R3 active: {rule_dist['R3']:,} products ({pct:.1f}%)")
    else:
        print("⚠ R3 not found in products")

    if 'R4' in rule_dist.index:
        pct = rule_dist['R4'] / len(df) * 100
        print(f"✓ R4 active: {rule_dist['R4']:,} products ({pct:.1f}%)")
```

**Expected:**
- Polyphenol: R3 should be ~20-35% of products
- Alkaloid: R3 ~10-15%, R4 ~10-15%
- Terpenoid: R3 ~10-20%

#### 2.3 Check Extreme Site Guards

```python
import json

for class_name in classes:
    with open(f'data/output/nplike/{class_name}-1X/SUMMARY.json') as f:
        summary = json.load(f)

    soft_warnings = summary.get('qa_paths', {}).get('extreme_site_soft_warning', 0)
    hard_skips = summary.get('qa_paths', {}).get('extreme_site_hard_skip', 0)

    print(f"{class_name}: soft_warnings={soft_warnings}, hard_skips={hard_skips}")
```

**Expected:**
- soft_warnings: 0-10 across all classes (very rare)
- hard_skips: 0 (should never trigger with current dataset)

---

### Task 3: Generate k=1 Completion Report

**Create comprehensive report comparing old vs new:**

```bash
# Generate report comparing:
# - Old k=1 results (base.parquet, limited rules)
# - New k=1 results (base_clean.parquet, expanded rules)
```

**Report should include:**
1. Total products by class (old vs new)
2. Rule distribution changes (R1+R5 only → R1+R3+R4+R5)
3. Products/parent statistics
4. Chemical space expansion metrics
5. Extreme site guard statistics

**Template structure:**

```markdown
# k=1 Enumeration Completion Report

## Summary
- Total parents: 56,708 (clean bases)
- Total products: [actual]
- Mean products/parent: [actual]
- Rules active: R1, R3, R4, R5
- Extreme site skips: [count]

## Comparison with Previous k=1

| Class | Old Products | New Products | Increase | Old Rules | New Rules |
|-------|--------------|--------------|----------|-----------|-----------|
| ... | ... | ... | ... | R1+R5 | R1+R3+R4+R5 |

## Rule Distribution

[Table showing % of products from each rule per class]

## Next Steps
- Proceed to k=2 enumeration
- Or: Additional QC/analysis
```

---

### Task 4: k=2 Enumeration (After k=1 Success)

**Objective:** Generate double-halogenation library

**Strategy:** k=2 is ~5-10x more expensive than k=1, run in batches

#### 4.1 Fast Batch (Est. 1-2 hours)

```bash
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid \
  --k-values 2
```

#### 4.2 Slow Batch (Est. 4-8 hours)

```bash
# Terpenoid k=2 (~3-5 hours, largest class)
halogenator enum-parquet \
  --input-parquet data/output/nplike/terpenoid/base_clean.parquet \
  --outdir data/output/nplike/terpenoid-2X \
  --k 2 \
  --np-class terpenoid

# Glycoside k=2 (~2-4 hours, sugar_mask overhead)
halogenator enum-parquet \
  --input-parquet data/output/nplike/glycoside/base_clean.parquet \
  --outdir data/output/nplike/glycoside-2X \
  --k 2 \
  --np-class glycoside
```

**Expected k=2 Results:**

| Class | Parents | Est. k=2 Products | Notes |
|-------|---------|-------------------|-------|
| lipid | 18 | ~50-100 | Very small |
| aa_peptide | 3,863 | ~300K-500K | R3+R4 combinations |
| polyphenol | 2,899 | ~300K-500K | R3 combinations (high) |
| alkaloid | 4,202 | ~400K-600K | R3+R4 combinations |
| terpenoid | 25,513 | ~3M-5M | Largest, may need monitoring |
| glycoside | 20,213 | ~2M-3M | Sugar_mask reduces space |
| **TOTAL** | **56,708** | **~6M-10M** | 3-5x k=1 size |

**k=2 Safety Checks:**
- Monitor extreme_site_hard_skip counts
- If any class shows >5% hard skips, may need to adjust thresholds
- Terpenoid/glycoside most likely to hit guards

---

### Task 5: Post-k=2 QC & Library Integration

#### 5.1 QC Each k=2 Library

Same as Task 2 but for k=2:
- Generate statistics
- Check rule distribution
- Check extreme site guards
- Verify products/parent are reasonable

#### 5.2 Merge k=1 + k=2 Libraries

**For each class:**

```python
import pandas as pd

class_name = 'polyphenol'

# Load k=1 and k=2
k1 = pd.read_parquet(f'data/output/nplike/{class_name}-1X/products.parquet')
k2 = pd.read_parquet(f'data/output/nplike/{class_name}-2X/products.parquet')

# Combine
combined = pd.concat([k1, k2], ignore_index=True)

# Deduplicate by InChIKey
combined_dedup = combined.drop_duplicates(subset='inchikey', keep='first')

# Export
combined_dedup.to_parquet(f'data/output/nplike/{class_name}-full/products.parquet', index=False)

print(f"{class_name}: k=1={len(k1):,}, k=2={len(k2):,}, combined={len(combined_dedup):,}")
```

#### 5.3 Cross-Class Deduplication (Optional)

**Check if products from different classes overlap:**

```python
# Load all class libraries
all_products = {}
for class_name in classes:
    df = pd.read_parquet(f'data/output/nplike/{class_name}-full/products.parquet')
    all_products[class_name] = set(df['inchikey'])

# Check overlaps
for i, c1 in enumerate(classes):
    for c2 in classes[i+1:]:
        overlap = all_products[c1] & all_products[c2]
        if len(overlap) > 0:
            print(f"{c1} ∩ {c2}: {len(overlap):,} products")
```

**Note:** Some product overlap is expected (e.g., simple halogenated aromatics)

---

### Task 6: Descriptor Calculation & Visualization

#### 6.1 Calculate Molecular Descriptors

```bash
# For each class library, calculate descriptors
python scripts/08a_fill_descriptors.py \
  --input data/output/nplike/polyphenol-full/products.parquet \
  --output data/output/nplike/polyphenol-full/products_with_descriptors.parquet
```

**Descriptors to calculate:**
- MW (molecular weight)
- LogP (lipophilicity)
- HBA/HBD (hydrogen bond acceptors/donors)
- TPSA (topological polar surface area)
- Aromatic rings count
- Rotatable bonds
- Complexity metrics

#### 6.2 Generate Visualization Galleries

```bash
# Sample and visualize products
python scripts/09_visualize_library.py html \
  -i data/output/nplike/polyphenol-full/products_with_descriptors.parquet \
  -o data/output/nplike/polyphenol-full/gallery.html \
  -n 500
```

**Visualizations to create:**
- HTML galleries with structure images
- Descriptor distribution plots
- PCA/t-SNE chemical space plots
- Rule distribution pie charts

---

### Task 7: Export for Virtual Screening

#### 7.1 Convert to VS Formats

```bash
# Export to SMILES for docking
python scripts/11_export_for_vs.py \
  --input data/output/nplike/polyphenol-full/products_with_descriptors.parquet \
  --output data/output/vs_libraries/polyphenol_full.smi \
  --format smiles

# Export to SDF for structure-based VS
python scripts/11_export_for_vs.py \
  --input data/output/nplike/polyphenol-full/products_with_descriptors.parquet \
  --output data/output/vs_libraries/polyphenol_full.sdf \
  --format sdf
```

#### 7.2 Apply Filtering (Optional)

**Drug-like filtering:**
- MW: 150-500 Da
- LogP: -0.4 to 5.6
- HBD: ≤5
- HBA: ≤10
- TPSA: ≤140 Å²
- Rotatable bonds: ≤10

```python
# Apply Lipinski/Veber filters
df = pd.read_parquet('products_with_descriptors.parquet')

drug_like = df[
    (df['MW'] >= 150) & (df['MW'] <= 500) &
    (df['LogP'] >= -0.4) & (df['LogP'] <= 5.6) &
    (df['HBD'] <= 5) &
    (df['HBA'] <= 10) &
    (df['TPSA'] <= 140) &
    (df['RotatableBonds'] <= 10)
]

print(f"Total: {len(df):,}, Drug-like: {len(drug_like):,} ({len(drug_like)/len(df)*100:.1f}%)")
```

---

### Task 8: Final Library Package & Documentation

#### 8.1 Create Library Release Package

```
halogenated_np_library_v1.0/
├── README.md                    # Library description, citation, usage
├── CHANGELOG.md                 # Version history
├── statistics_summary.json      # Global statistics
├── libraries/
│   ├── by_class/
│   │   ├── flavone/
│   │   │   ├── k1_products.parquet
│   │   │   ├── k2_products.parquet
│   │   │   └── full_products.parquet
│   │   ├── terpenoid/
│   │   ├── alkaloid/
│   │   └── ... (other classes)
│   ├── merged/
│   │   ├── all_k1.parquet       # All k=1 products
│   │   ├── all_k2.parquet       # All k=2 products
│   │   └── all_full.parquet     # Complete library
│   └── filtered/
│       ├── drug_like.parquet    # Lipinski/Veber filtered
│       └── lead_like.parquet    # Lead-like filtered
├── vs_formats/
│   ├── all_full.smi             # SMILES for VS
│   ├── all_full.sdf             # SDF for docking
│   └── by_class/                # Class-specific exports
└── metadata/
    ├── master_parent_index.parquet
    ├── rule_distribution.json
    ├── descriptor_statistics.json
    └── classification_report.md
```

#### 8.2 Generate Documentation

**README.md should include:**
- Library overview (size, diversity, rules used)
- How to cite
- File format descriptions
- Usage examples (Python, RDKit)
- Descriptor definitions
- Known limitations

**CHANGELOG.md should include:**
- Version 1.0 baseline
- Rules used (R1, R3, R4, R5)
- Classification system (primary_class priority)
- QC metrics

---

## Quick Start Commands for Next Session

```bash
# 1. Verify clean bases exist
ls -lh data/output/nplike/*/base_clean.parquet

# 2. Run k=1 enumeration (all classes)
python scripts/04_enum_halogen_all_classes.py \
  --classes lipid aa_peptide polyphenol alkaloid terpenoid glycoside \
  --k-values 1
# BUT FIRST: Edit script to use base_clean.parquet instead of base.parquet!

# 3. Check results
python -c "
import pandas as pd
for c in ['lipid', 'aa_peptide', 'polyphenol', 'alkaloid', 'terpenoid', 'glycoside']:
    df = pd.read_parquet(f'data/output/nplike/{c}-1X/products.parquet')
    print(f'{c:15} {len(df):>8,} products')
"

# 4. Analyze rule distribution
python -c "
import pandas as pd
for c in ['polyphenol', 'alkaloid']:
    df = pd.read_parquet(f'data/output/nplike/{c}-1X/products.parquet')
    print(f'\n{c}:')
    print(df['rule'].value_counts())
"
```

---

## Important Files & Locations

**Configuration:**
- `configs/halogen_rules_by_class.yaml` - Main config (R3/R4 enabled, no max_sites)
- Backup: `configs/halogen_rules_by_class.yaml.bak_before_r3r4`

**Data:**
- Clean bases: `data/output/nplike/<class>/base_clean.parquet`
- Master index: `data/output/nplike/master_parent_index.parquet`

**Scripts:**
- Enumeration: `scripts/04_enum_halogen_all_classes.py`
- Summaries: `scripts/05_summaries.py`
- Descriptors: `scripts/08a_fill_descriptors.py`
- Visualization: `scripts/09_visualize_library.py`
- Export: `scripts/11_export_for_vs.py`

**Reports:**
- POC validation: `POC_NO_LIMITS_VALIDATION_REPORT.md`
- Overlap analysis: `NP_CLASSIFICATION_OVERLAP_REPORT.md`
- Rules audit: `HALOGEN_RULES_AUDIT_REPORT.md`
- Phase 1-3: `PHASE_1-3_COMPLETION_REPORT.md`

---

## Expected Timeline (If Running Sequentially)

1. **k=1 enumeration (6 classes):** ~60 minutes
2. **k=1 QC & analysis:** ~30 minutes
3. **k=2 enumeration (6 classes):** ~6-8 hours
4. **k=2 QC & analysis:** ~1 hour
5. **Library merging & dedup:** ~30 minutes
6. **Descriptor calculation:** ~1-2 hours
7. **Visualization:** ~1 hour
8. **Export & packaging:** ~30 minutes

**Total estimated time:** ~12-15 hours (can run overnight)

---

## Troubleshooting Guide

**If enumeration fails:**
1. Check input file exists: `base_clean.parquet`
2. Check config is loaded: `--np-class <name>` matches YAML
3. Check disk space (products can be large)
4. Check SUMMARY.json for error messages

**If too many extreme_site_hard_skip:**
1. Identify which molecules trigger guards
2. Consider raising EXTREME_SITE_HARD_THRESHOLD to 80
3. Or: filter out problematic molecules before enumeration

**If R3/R4 not producing products:**
1. Verify YAML config has R3/R4 in include_rules
2. Check class_config.py SEMANTIC_TO_LEGACY_RULES mapping
3. Verify rules.py has R3/R4 SMIRKS defined

**If product count much lower than expected:**
1. Check if using base_clean.parquet (not base.parquet)
2. Verify parent count in clean base
3. Check SUMMARY.json qa_paths for filters

---

## Success Metrics

**k=1 enumeration is successful if:**
- ✅ All 6 classes complete without errors
- ✅ Total products: 1.5M - 2.2M (updated with ALPHA_CARBONYL)
- ✅ R3 appears in polyphenol/alkaloid/terpenoid rule distribution
- ✅ R4 appears in alkaloid/aa_peptide rule distribution
- ✅ ALPHA_CARBONYL appears in at least some classes
- ✅ extreme_site_hard_skip < 1% across all classes
- ✅ No duplicate InChIKeys within each class

**k=2 enumeration is successful if:**
- ✅ All classes complete (may take 6-8 hours)
- ✅ Total products: 6M - 10M
- ✅ Products/parent mean < 200 for all classes
- ✅ extreme_site_hard_skip < 5% across all classes

---

## Contact & Support

**If issues arise:**
1. Check SUMMARY.json for detailed QA statistics
2. Review logs for warnings/errors
3. Compare with POC results to identify deviations
4. Check GitHub issues: https://github.com/anthropics/claude-code/issues

**Key context for debugging:**
- Rules version: "1.0.0" with R3/R4/ALPHA_CARBONYL expansion
- Active rules in k=1: R1, R3, R4, R5, ALPHA_CARBONYL, RING_SP3, PRIMARY_OH
- Strategy: Zero max_sites, engine guards only (soft=40, hard=60)
- Classification: Primary_class priority system
- Base libraries: Clean bases (no overlaps, 67,983 total parents)

---

**Document Status:** ✅ Ready for next session
**Last Updated:** 2025-12-04 (Updated with ALPHA_CARBONYL in k=1)
**Ready to Execute:** Task 1 (Full k=1 Enumeration)
