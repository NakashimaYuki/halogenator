# NP Classification v2.0 - Quick Start Guide

**TL;DR:** New skeleton-based classification that properly handles glycosides. Run tests first, then classify your library.

---

## Quick Test

Verify everything works:
```bash
python scripts/test_np_classification.py
# Expected: 10 passed, 0 failed ✅
```

---

## Basic Usage

### Classify Your Library
```bash
python scripts/02_partition_nplike_by_class.py \
    -i data/input/your_library.parquet \
    -o data/output/classified \
    --split
```

**What this does:**
- Classifies molecules into skeleton classes (aa_peptide, alkaloid, terpenoid, polyphenol, lipid, polysaccharide, other)
- Adds sugar-related tags (glycoside, saponin, cglycoside_like)
- Excludes Flavone library by default (if exists)
- Outputs combined file + per-class shards (with `--split`)

---

## Key Improvements Over v1.0

| What | Before | After |
|------|--------|-------|
| Flavonoid glycosides | "glycoside" ❌ | "polyphenol" + glycoside tag ✅ |
| Saponins | "glycoside" ❌ | "terpenoid" + saponin tag ✅ |
| Classification basis | Full molecule | Sugar-free scaffold ✅ |
| Method | Priority chain | Independent scoring ✅ |

---

## Understanding the Output

### Schema
Each molecule gets 4 new fields:
- `np_primary_class`: Main skeleton class
- `np_tags`: JSON list of all applicable classes/tags
- `np_sugar_fraction`: 0.0-1.0 (fraction of molecule that is sugar)
- `np_num_sugar_rings`: Number of sugar rings detected

### Example Results
```python
# Quercetin-3-glucoside (flavonoid glycoside)
primary_class = "polyphenol"
tags = ["glycoside", "polyphenol", "terpenoid"]
sugar_fraction = 0.394
num_sugar_rings = 1

# Ginsenoside (triterpenoid saponin)
primary_class = "terpenoid"
tags = ["glycoside", "saponin", "terpenoid"]
sugar_fraction = 0.289
num_sugar_rings = 1

# Caffeine (alkaloid, no sugar)
primary_class = "alkaloid"
tags = ["alkaloid"]
sugar_fraction = 0.0
num_sugar_rings = 0
```

---

## Common Queries

### Find All Glycosides
```python
import pandas as pd
df = pd.read_parquet("data/output/classified/nplike_with_classes.parquet")

# Any molecule with sugar
glycosides = df[df['np_tags'].str.contains('glycoside')]
```

### Find Terpenoid Saponins
```python
saponins = df[
    (df['np_primary_class'] == 'terpenoid') &
    (df['np_tags'].str.contains('saponin'))
]
```

### Find Polyphenol Aglycones (No Sugar)
```python
polyphenol_aglycones = df[
    (df['np_primary_class'] == 'polyphenol') &
    (~df['np_tags'].str.contains('glycoside'))
]
```

---

## Advanced Options

### Custom Exclusions
```bash
python scripts/02_partition_nplike_by_class.py \
    -i data/input/library.parquet \
    -o data/output/classified \
    --exclude-parquet path/to/exclude1.parquet path/to/exclude2.parquet
```

### Adjust Batch Size (for memory management)
```bash
python scripts/02_partition_nplike_by_class.py \
    -i data/input/library.parquet \
    -o data/output/classified \
    --batch-size 5000  # default: 2000
```

---

## Architecture (Simplified)

```
Molecule → Remove sugars → Classify skeleton → Add sugar tags → Output
           (sugar_mask)    (scoring system)    (glycoside etc)
```

**Key insight:** Classify the **skeleton** (aglycone), not the decorated molecule.

---

## Troubleshooting

### Tests Fail
```bash
# Check imports
python -c "from halogenator.sugar_mask import analyze_sugars_for_classification"
python -c "from np_classify_utils import make_aglycone_from_mask"
```

### RDKit Warnings
```
[15:45:23] non-ring atom 1 marked aromatic
```
→ Safe to ignore. This is a SMILES parsing warning, not an error.

### Unexpected Classification
Enable debug logging to see scores:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

## Files Created/Modified

### New Files
- `scripts/np_classify_utils.py` - Aglycone generation utility
- `scripts/test_np_classification.py` - Test suite
- `NP_CLASSIFICATION_V2_IMPLEMENTATION.md` - Full documentation
- `NP_CLASSIFICATION_QUICK_START.md` - This guide

### Modified Files
- `scripts/02_partition_nplike_by_class.py` - Complete rewrite with new architecture
- `src/halogenator/sugar_mask.py` - Added SugarFeatures dataclass and analysis function

---

## Next Steps

1. ✅ Run tests to verify installation
2. ✅ Classify a small test set (100 molecules)
3. ✅ Inspect results to validate classification quality
4. ✅ Run on full library
5. ✅ Use classified data for downstream halogenation enumeration

---

## Support

- Full documentation: `NP_CLASSIFICATION_V2_IMPLEMENTATION.md`
- Test suite: `scripts/test_np_classification.py`
- Source code: `scripts/02_partition_nplike_by_class.py`

---

**Status: PRODUCTION READY** 🚀
All tests passing. System ready for full-scale library classification.
