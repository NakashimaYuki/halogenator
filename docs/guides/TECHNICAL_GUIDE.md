# Halogenator Technical Guide

**Last Updated**: 2025-11-04
**Version**: Post-K2-Iteration Complete

---

## Table of Contents

1. [History as Single Source of Truth](#history-as-single-source-of-truth)
2. [JSON Field Type Safety](#json-field-type-safety)
3. [Step Construction Factory](#step-construction-factory)
4. [Folding and Deduplication Effects](#folding-and-deduplication-effects)
5. [Common Phenomena Explained](#common-phenomena-explained)
6. [Best Practices](#best-practices)

---

## History as Single Source of Truth

### Core Principle

**The `substitutions` (history) field is the authoritative source for all k-level metrics.**

All k-level counters (`k`, `k_ops`, `k_atoms`) should be **derived** from the history array, not stored separately or copied from budget state.

### Why This Matters

**Previous Bug**: R6_methyl k-level mislabeling occurred because `emit_product` used `budget_state.k_ops` instead of computing from history. When budget state was incorrect or stale, products were mislabeled (e.g., k=2 products marked as k=1).

**Fix Applied**: All k-level metrics are now computed from history:

```python
def emit_product(product_mol, parent_mol, rule, halogen, history, ...):
    """Create product record with metrics derived from history."""

    # Derive k-level metrics from history (single source of truth)
    k_ops = len(history)  # Number of operations
    k_atoms = sum(step.get('atom_cost', 0) for step in history)  # Number of atoms

    return {
        'smiles': ...,
        'k': k_ops,  # Derived from history length
        'k_ops': k_ops,  # Derived
        'k_atoms': k_atoms,  # Derived from sum of atom_costs
        'substitutions': history or [],  # Typed field (becomes substitutions_json)
        ...
    }
```

### History Structure

Each history item (step) must contain:

**Required Fields**:
- `rule`: Rule identifier ('R1', 'R2a', 'R2b', 'R3', 'R6_methyl')
- `site`: Atom index where substitution occurred
- `halogen`: Halogen symbol ('F', 'Cl', 'Br', 'I')
- `atom_cost`: Number of halogen atoms added (1 for single, 3 for CF3/CCl3)
- `depth`: BFS depth when this step was created

**Optional Fields (R1/R2/R3)**:
- `sym`: Symmetry class from canonical ranking
- `ring_tag`: Ring label ('A', 'B', 'C') for aromatic sites

**Optional Fields (R6_methyl)**:
- `type`: 'step' or 'macro'
- `label`: 'CF3' or 'CCl3' for macro mode

**Legacy Fields (deprecated)**:
- `k_ops`, `k_atoms`, `budget_mode`: Kept for backward compatibility, but should not be used for calculation

---

## JSON Field Type Safety

### Problem Context

The io_utils module serializes typed Python fields (lists/dicts) to JSON strings for storage in Parquet/CSV. Previously, **all missing fields defaulted to `{}`** (empty dict), which caused issues for list-semantic fields like `substitutions`.

**Bug Manifestation**:
```python
# Before fix
record = {'smiles': 'CCO'}  # Missing 'substitutions'
prepared = _prepare_records_for_table([record])
# Result: substitutions_json = '{}'  ← WRONG! Should be '[]'
```

### Fix Applied

Created a **type table** to distinguish list vs dict fields:

```python
# In io_utils.py
_JSON_LIST_FIELDS = frozenset({'substitutions', 'sugar_mask_atoms', 'sugar_rings'})
_JSON_DICT_FIELDS = frozenset({'constraints_violations'})

def _get_json_default(field_name):
    """Get appropriate default value based on semantic type."""
    if field_name in _JSON_LIST_FIELDS:
        return []
    elif field_name in _JSON_DICT_FIELDS:
        return {}
    else:
        return {}  # Fallback
```

Now:
- `substitutions` defaults to `[]` (empty list)
- `constraints_violations` defaults to `{}` (empty dict)

### Validation

```python
# After fix
record = {'smiles': 'CCO'}  # Missing 'substitutions'
prepared = _prepare_records_for_table([record])
# Result: substitutions_json = '[]'  ← CORRECT!
```

---

## Step Construction Factory

### Purpose

Unified factory function `make_history_step()` provides a single, consistent way to create history items across all rules, ensuring all required fields are present and avoiding field omission bugs.

### Usage

**Location**: `src/halogenator/enumerate_k.py`

**Basic Usage (R1/R2/R3)**:
```python
from halogenator.enumerate_k import make_history_step

step = make_history_step(
    rule='R1',
    site=5,
    halogen='F',
    atom_cost=1,
    depth=1,
    sym=2,
    ring_tag='A'
)
```

**R6_methyl Step Mode**:
```python
step = make_history_step(
    rule='R6_methyl',
    site=10,
    halogen='Cl',
    atom_cost=1,
    depth=2,
    step_type='step',
    k_ops=2,  # Legacy field
    k_atoms=2,  # Legacy field
    budget_mode='ops'  # Legacy field
)
```

**R6_methyl Macro Mode**:
```python
step = make_history_step(
    rule='R6_methyl',
    site=10,
    halogen='F',
    atom_cost=3,  # CF3 = 3 atoms
    depth=2,
    step_type='macro',
    macro_label='CF3',
    k_ops=2,
    k_atoms=5,
    budget_mode='atoms'
)
```

### Benefits

1. **Prevents field omission**: All required fields are explicit parameters
2. **Consistent defaults**: Optional fields get appropriate defaults
3. **Self-documenting**: Function signature shows all available fields
4. **Backward compatible**: Existing inline dict construction still works

### Migration Strategy

- **New code**: Use `make_history_step()` for all history item creation
- **Existing code**: Can continue using inline dicts, but consider migrating for better maintainability
- **Mixed approach**: Use factory for new rules, keep existing code unchanged

---

## Folding and Deduplication Effects

### Three Independent Mechanisms

#### 1. Symmetry Folding (`--no-sym-fold` to disable)

**What it does**: At enumeration time, reduce equivalent sites to one representative based on molecular symmetry.

**Effect**: Fewer enumeration paths → fewer products

**Example**: Benzene has 6 equivalent C-H. With folding, only 1 is enumerated. Without folding, all 6 are enumerated (but may produce same structure).

#### 2. InChI Deduplication (`--no-dedup` to disable)

**What it does**: After enumeration, merge products with identical InChI keys (same chemical structure).

**Effect**: Fewer final products (removes structural duplicates)

**Example**: Two enumeration paths produce same molecule → dedup keeps only one

#### 3. Sugar Mask (enabled in strict mode)

**What it does**: Protect glycosidic oxygens and sugar ring -CH₂OH groups from halogenation.

**Effect**: Fewer sites available for R3 (hydroxyl halogenation) → fewer products

**Example**: Glycoside with 5 hydroxyl groups, 2 in sugar → strict mode only allows 3 to react

### Mode Combinations

| Mode | Folding | Dedup | Sugar Mask | Use Case |
|------|---------|-------|------------|----------|
| **Strict** | ✓ | ✓ | ✓ | Production (fewer, high-quality products) |
| **Raw** | ✗ | ✗ | ✗ | Analysis (all enumeration paths) |
| **Balanced** | ✓ | ✓ | ✗ | More products, still deduplicated |

### Impact on Product Counts

**Example: G1 Glycoside**
```
Strict mode:  590 products
Raw mode:    2548 products (4.3x difference)

Breakdown:
- Folding reduction: ~30%
- Dedup reduction: ~60%
- Sugar mask reduction: ~70% (for R3)
```

**Example: 8-PN Prenyl**
```
Strict mode:  340 products
Raw mode:    1008 products (3x difference)

Breakdown:
- Folding reduction: ~40%
- R2b fallback addition: +200%
```

---

## Common Phenomena Explained

### Phenomenon 1: Duplicate Products in Raw Mode (8-PN Example)

**Observation**: 8-prenylnaringenin's prenyl group has two methyl groups. In raw mode, both produce F/Cl products that look identical in PyMOL.

**Explanation**:
- The two methyls are **topologically equivalent** (gem-dimethyl)
- Without symmetry folding (`--no-sym-fold`), both are enumerated separately
- They produce **structurally identical** molecules (same InChI)
- Without dedup (`--no-dedup`), both are retained as separate records
- PyMOL shows identical 3D structures (because they ARE identical)

**Is this a bug?**: ❌ No. This is expected behavior in raw mode.

**How to avoid**:
1. Enable symmetry folding (remove `--no-sym-fold`)
2. Enable dedup (remove `--no-dedup`)
3. Use rule-level micro-folding (advanced option, not implemented)

**Verification**:
```python
import pandas as pd
df = pd.read_parquet('data/output/8pn_raw_v4/products_k2.parquet')
k1 = df[df['k'] == 1]
r6 = k1[k1['rule'] == 'R6_methyl']

# Check for duplicate InChIKeys
for hal in ['F', 'Cl']:
    subset = r6[r6['halogen'] == hal]
    dup_count = len(subset) - subset['inchikey'].nunique()
    print(f"{hal}: {dup_count} duplicates found")

# Expected: F: 1 duplicate, Cl: 1 duplicate
```

---

### Phenomenon 2: Uneven K=2 Counts (G1 Example)

**Observation**: In G1-strict scenario, some k=1 parents produce only 2-3 k=2 children, while others produce 20-30.

**Explanation**: This is due to **rule site consumption asymmetry + sugar mask + dedup**:

#### R1 Parent (芳香C-H halogenation)
```
K=1: R1 consumes 1 aromatic C-H on ring A
↓
K=2 available sites:
- Remaining aromatic C-H on ring A (limited, ~4-5 sites)
- Hydroxyl on ring B (heavily restricted by sugar mask)
↓
Result: Only ~8.8 k=2 children on average
```

#### R3 Parent (hydroxyl halogenation)
```
K=1: R3 consumes 1 hydroxyl
↓
K=2 available sites:
- ALL aromatic C-H on ring A (5 sites × 4 halogens = 20 theoretical)
- Remaining hydroxyl sites
↓
Result: ~23.6 k=2 children on average (2.7x more!)
```

**Data Evidence**:
```
R1 parents (20 total):
  Avg k=2 children: 8.8
  Rule distribution: R1=144 (81.8%), R3=32 (18.2%)

R3 parents (16 total):
  Avg k=2 children: 23.6
  Rule distribution: R1=282 (74.6%), R3=96 (25.4%)

Ratio: 2.68x difference
```

**Extreme Cases**:
- Some R1 parents produce only **4 k=2 children** (F/Cl/Br/I, one each via R3)
- These have exhausted aromatic sites AND only 1 available hydroxyl
- Sugar mask blocks all sugar ring hydroxyls

**Is this a bug?**: ❌ No. This is expected behavior due to:
1. Rule site consumption characteristics
2. Sugar mask protecting glycosidic groups
3. Dedup removing structural duplicates (60% of attempts)

**How to get more uniform counts**:
1. Use raw mode (`--no-sugar-mask --no-dedup`)
2. Or accept that different rule types have different "fertility"

**Verification**:
```python
import pandas as pd
df = pd.read_parquet('data/output/g1_strict_v4/products_k2.parquet')
k1_df = df[df['k'] == 1]
k2_df = df[df['k'] == 2]

# Group k=2 by parent rule
for rule in ['R1', 'R3']:
    parents = k1_df[k1_df['rule'] == rule]['inchikey']
    k2_children = k2_df[k2_df['parent_inchikey'].isin(parents)]
    print(f"{rule} parents: avg {len(k2_children)/len(parents):.1f} k=2 children")

# Expected: R1: ~8.8, R3: ~23.6
```

---

## Best Practices

### 1. Always Use History for K-Level Metrics

❌ **Don't**:
```python
record = {
    'k': budget_state.k_ops,  # May be stale or incorrect
    'k_ops': budget_state.k_ops,
    'k_atoms': budget_state.k_atoms,
    ...
}
```

✅ **Do**:
```python
k_ops = len(history)
k_atoms = sum(step.get('atom_cost', 0) for step in history)

record = {
    'k': k_ops,
    'k_ops': k_ops,
    'k_atoms': k_atoms,
    'substitutions': history or [],  # Typed field
    ...
}
```

### 2. Use Typed Fields for JSON Serialization

❌ **Don't**:
```python
record = {
    'substitutions_json': json.dumps(history)  # Manual serialization
}
```

✅ **Do**:
```python
record = {
    'substitutions': history or []  # Typed field; io_utils handles serialization
}
```

### 3. Use Factory for New History Items

❌ **Don't** (error-prone):
```python
history_item = {
    'rule': 'R1',
    'site': 5,
    'halogen': 'F',
    # Oops, forgot atom_cost!
}
```

✅ **Do**:
```python
from halogenator.enumerate_k import make_history_step

history_item = make_history_step(
    rule='R1',
    site=5,
    halogen='F',
    atom_cost=1,  # Required parameter, can't forget
    depth=1
)
```

### 4. Choose Appropriate Mode for Use Case

**Production / Final Results**:
```bash
# Use strict mode (default)
python -m halogenator.cli enum -c config.yaml
```

**Analysis / Path Exploration**:
```bash
# Use raw mode
python -m halogenator.cli enum -c config.yaml \
  --no-constraints --no-sugar-mask --no-sym-fold --no-dedup --r2-fallback
```

**Balanced**:
```bash
# Keep folding and dedup, disable sugar mask
python -m halogenator.cli enum -c config.yaml --no-sugar-mask
```

### 5. Validate Outputs

Always run validation scripts after enumeration:

```bash
# Field consistency
python scripts/validate_field_consistency.py data/output/*/products_k2.parquet

# Parent chain integrity
python scripts/validate_parent_chain.py data/output/*/products_k2.parquet

# SDF purity
python scripts/check_k1_sdf_purity.py data/output/* --check-k2 --check-structural

# All-in-one
python scripts/ci_validation_gate.py data/output/* --paranoid
```

### 6. Document Non-Obvious Behavior

When observing unexpected product counts or distributions:
1. Check mode configuration (strict vs raw)
2. Analyze QA summary JSON for dedup/filter rates
3. Compare with alternate mode runs
4. Document findings for future reference

---

## Appendix: Testing

### Regression Test Suite

Run complete regression tests:
```bash
python test_regression_suite.py
```

### Individual Component Tests

```bash
# JSON type safety
python tests/test_io_utils_json_type_safety.py

# History factory
python tests/test_history_step_factory.py

# 8-PN phenomenon verification
python tests/test_8pn_equivalent_methyl_verification.py

# G1 distribution analysis
python tests/test_g1_k2_low_count_verification.py
```

---

## References

- **K2_ITERATION_COMPLETE_REPORT.md**: Full validation results for 4 scenarios
- **G1_K2_DISTRIBUTION_ANALYSIS_REPORT.md**: Detailed analysis of k=2 count distribution
- **src/halogenator/enumerate_k.py**: Main enumeration engine
- **src/halogenator/io_utils.py**: I/O and serialization utilities

---

**End of Technical Guide**
