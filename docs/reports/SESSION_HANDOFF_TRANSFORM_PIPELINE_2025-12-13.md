# Session Handoff Report: Transform Pipeline Preparation

**Date:** 2025-12-13
**Session Focus:** Sugar_mask fix validation, K=2 enumeration completion, Transform pipeline preparation
**Status:** Ready for transform pipeline implementation
**Next Session:** Implement StreamingParquetWriter improvements and run transform library generation

---

## Executive Summary

### What's Ready
- ✅ **Complete halogenated library:** 67.4M products (K=1 + K=2, all 7 classes)
- ✅ **Sugar_mask fix:** Production-validated, 0 errors across 35K parents
- ✅ **New transform rules:** 2 rules added to transforms.yaml
- ✅ **Critical bug fixed:** Schema mismatch in parallel_enum.py (788 flushes successful)

### What's Next
- 🔄 **Transform pipeline improvements:** Add flush_interval, memory monitoring, schema enforcement to 08_transform_library_v2.py
- 🔄 **Transform library generation:** Apply all transform rules to 67M halogenated products

---

## PART A: COMPLETED TASKS (Summary)

### 1. Sugar_mask Fix Implementation ✅

**Problem Solved:** Sugar_mask was non-functional in enumerate_k1.py, causing ~975K unwanted products

**Solution:**
- Fixed `src/halogenator/enumerate_k1.py` (4 critical fixes)
- Line 103: Added `sugar_cfg` to config dict
- Lines 343-427: Implemented isotope tagging strategy
- Lines 456, 507: Fixed R2a/R2b to pass sugar_mask

**Validation Results:**
- K=1 reduction: 27.6% (974,576 products eliminated)
- Old total: 3,525,292 → New total: 2,550,716
- All validation criteria passed

**Git Commit:** `2ac15dd` on `fix/pr2-contract-and-sugar-accept`

---

### 2. K=1 + K=2 Enumeration Complete ✅

**Final Library Statistics:**

| Class | K=1 | K=2 | Total | Parents |
|-------|-----|-----|-------|---------|
| aa_peptide | 30,564 | 582,730 | 613,294 | 1,119 |
| alkaloid | 322,280 | 8,461,937 | 8,784,217 | 7,871 |
| lipid | 33,084 | 180,224 | 213,308 | 6,247 |
| polyphenol | 518,852 | 13,790,820 | 14,309,672 | 13,168 |
| terpenoid | 1,528,400 | 40,170,816 | 41,699,216 | 35,131 |
| polysaccharide | 20,524 | 293,616 | 314,140 | 566 |
| other | 97,012 | 1,419,941 | 1,516,953 | 4,146 |
| **TOTAL** | **2,550,716** | **64,900,084** | **67,450,800** | **68,248** |

**Key Achievement:**
- Terpenoid K=2: 40.17M products, 1,143.5 products/parent
- 788 flushes, 0 errors (schema fix successful)

---

### 3. Schema Mismatch Bug Fixed ✅

**Problem:** PyArrow inferred null type for all-null sym/sym_class fields, causing crashes on subsequent flushes

**Solution Applied to `parallel_enum.py`:**

**Lines 433-453:** Cast table to established schema on every flush (not just first)
```python
else:
    # CRITICAL: Cast table to file schema on every flush
    # Subsequent batches may have all-null fields that PyArrow infers as null type
    # Must cast to the established schema to prevent mismatch
    try:
        table = table.cast(schema)
    except Exception as e:
        LOG.warning(f"Schema cast failed, attempting column-by-column cast: {e}")
        # Fallback: column-by-column cast
        import pyarrow as pa
        arrays = []
        for field in schema:
            if field.name in table.column_names:
                try:
                    arrays.append(table.column(field.name).cast(field.type))
                except:
                    arrays.append(table.column(field.name))
            else:
                arrays.append(pa.nulls(len(table), type=field.type))
        table = pa.Table.from_arrays(arrays, schema=schema)
```

**Validation:** 788 consecutive flushes successful on terpenoid K=2

---

### 4. New Transform Rules Added ✅

**File:** `configs/transforms.yaml`

**Rule 1: FG_PHENOL_OH__OH__TO__CH2CH2NH2**
```yaml
- name: FG_PHENOL_OH__OH__TO__CH2CH2NH2
  code: FG_PHENOL_OH__OH__TO__CH2CH2NH2
  legacy_rule_id: null
  label: "OH->CH2CH2NH2"
  family: phenol_mod
  scope_classes: ["polyphenol", "terpenoid"]
  query_smarts: "[c:1][OX2H:2]"
  smirks: "[c:1][OX2H:2]>>[c:1]CCN"
  highlight_mapnums: [1]
  product_highlight_symbol: "C"
```

**Rule 2: FG_PHENOL_OH__OH__TO__CH2CH2COOH**
```yaml
- name: FG_PHENOL_OH__OH__TO__CH2CH2COOH
  code: FG_PHENOL_OH__OH__TO__CH2CH2COOH
  legacy_rule_id: null
  label: "OH->CH2CH2COOH"
  family: phenol_mod
  scope_classes: ["polyphenol", "terpenoid"]
  query_smarts: "[c:1][OX2H:2]"
  smirks: "[c:1][OX2H:2]>>[c:1]CCC(=O)O"
  highlight_mapnums: [1]
  product_highlight_symbol: "C"
```

**Validation:** Both rules produce correct products (phenethylamine, hydrocinnamic acid derivatives)

---

## PART B: PENDING TASKS (Detailed Implementation Guide)

### Task 1: Improve StreamingParquetWriter in 08_transform_library_v2.py

**Goal:** Add flush_interval, memory monitoring, and schema enforcement to prevent OOM and schema errors on 67M product transform job

**Priority:** 🔴 **CRITICAL** - Must complete before running transform pipeline

**Estimated Time:** 1 hour

---

#### Step 1.1: Add Product Buffer and Flush Interval

**File:** `scripts/08_transform_library_v2.py`

**Location:** Lines 301-328 (StreamingParquetWriter class)

**Current Code (Lines 301-328):**
```python
class StreamingParquetWriter:
    """Stream products to Parquet file without accumulating in memory."""

    def __init__(self, output_path: Path, schema: pa.Schema):
        self.output_path = output_path
        self.writer = None
        self.schema = schema
        self.total_written = 0

    def write_batch(self, products: List[Dict]):
        """Write a batch of products."""
        if not products:
            return

        # Convert to pyarrow table
        table = pa.Table.from_pylist(products, schema=self.schema)

        # Open writer on first write
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)

        self.writer.write_table(table)
        self.total_written += len(products)

    def close(self):
        if self.writer:
            self.writer.close()
```

**REPLACE WITH:**
```python
class StreamingParquetWriter:
    """
    Stream products to Parquet file with buffering and memory monitoring.

    Improvements from parallel_enum.py pattern:
    - Product buffer with configurable flush_interval
    - Optional memory monitoring (psutil)
    - Schema enforcement on every flush (prevents null-type inference bugs)
    - Adaptive flushing under memory pressure
    """

    def __init__(
        self,
        output_path: Path,
        schema: pa.Schema,
        flush_interval: int = 10000,
        memory_threshold: float = 0.8
    ):
        self.output_path = output_path
        self.writer = None
        self.schema = schema
        self.total_written = 0
        self.flush_interval = flush_interval
        self.memory_threshold = memory_threshold

        # Product buffer (NEW)
        self.buffer = []
        self.flush_count = 0

        # Memory monitoring (optional, like parallel_enum.py)
        self.memory_monitor_enabled = self._check_psutil()

        import logging
        self.logger = logging.getLogger(__name__)

    def _check_psutil(self) -> bool:
        """Check if psutil available for memory monitoring."""
        try:
            import psutil
            return True
        except ImportError:
            return False

    def _check_memory(self) -> float:
        """Get current memory usage percentage (0.0-1.0)."""
        if not self.memory_monitor_enabled:
            return 0.0
        try:
            import psutil
            return psutil.virtual_memory().percent / 100.0
        except:
            return 0.0

    def write_batch(self, products: List[Dict]):
        """
        Add products to buffer, flush when interval reached or memory high.

        This pattern matches parallel_enum.py for consistency.
        """
        if not products:
            return

        # Add to buffer
        self.buffer.extend(products)

        # Check flush conditions
        mem_usage = self._check_memory()
        should_flush = (
            len(self.buffer) >= self.flush_interval or
            (self.memory_monitor_enabled and mem_usage > self.memory_threshold)
        )

        if should_flush:
            if mem_usage > self.memory_threshold:
                self.logger.warning(
                    f"Memory-triggered flush: {mem_usage:.1%} > {self.memory_threshold:.1%}"
                )
            self._flush()

    def _flush(self):
        """
        Flush buffer to disk with schema enforcement.

        CRITICAL: Cast to schema on every flush to prevent null-type inference issues.
        This matches the fix in parallel_enum.py (lines 433-453).
        """
        if not self.buffer:
            return

        # Convert buffer to PyArrow table
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)

        # Initialize writer on first flush
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
            self.logger.info(f"ParquetWriter initialized (flush_interval={self.flush_interval:,})")
        else:
            # CRITICAL: Cast table to established schema on every flush
            # Prevents null-type inference issues (terpenoid k=2 bug)
            try:
                table = table.cast(self.schema)
            except Exception as e:
                self.logger.warning(f"Schema cast failed, attempting column-by-column: {e}")
                # Fallback: column-by-column cast
                arrays = []
                for field in self.schema:
                    if field.name in table.column_names:
                        try:
                            arrays.append(table.column(field.name).cast(field.type))
                        except:
                            arrays.append(table.column(field.name))
                    else:
                        # Column missing, create null column
                        arrays.append(pa.nulls(len(table), type=field.type))
                table = pa.Table.from_arrays(arrays, schema=self.schema)

        # Write to disk
        self.writer.write_table(table)
        self.total_written += len(self.buffer)
        self.flush_count += 1

        # Log progress
        if self.flush_count % 10 == 0:
            self.logger.info(
                f"Flush #{self.flush_count}: {self.total_written:,} total products written"
            )

        # Clear buffer
        self.buffer.clear()

    def close(self):
        """Flush remaining buffer and close writer."""
        # Final flush
        if self.buffer:
            self.logger.info(f"Final flush: {len(self.buffer):,} products")
            self._flush()

        # Close writer
        if self.writer:
            self.writer.close()
            self.logger.info(
                f"StreamingParquetWriter closed: {self.total_written:,} total products, "
                f"{self.flush_count} flushes"
            )
```

**Why This Matters:**
- Reduces I/O overhead (fewer, larger writes vs. many small writes)
- Prevents OOM via memory monitoring
- Prevents schema mismatch crashes (proven fix from parallel_enum.py)
- Provides progress logging for long-running jobs

---

#### Step 1.2: Add CLI Parameter for flush_interval

**File:** `scripts/08_transform_library_v2.py`

**Location:** Around line 1040 (in argparse section)

**Current Code:**
```python
parser_apply.add_argument('--batch-size', type=int, default=50000, help='Batch size (default: 50000)')
parser_apply.add_argument('--workers', type=int, default=6, help='Parallel workers (default: 6)')
parser_apply.add_argument('--resume', action='store_true', help='Resume from existing output')
```

**ADD AFTER `--workers`:**
```python
parser_apply.add_argument(
    '--flush-interval',
    type=int,
    default=10000,
    help='Products to buffer before flushing to disk (default: 10000). '
         'Larger values = better performance but higher memory usage. '
         'Recommended: 10000-50000 depending on available RAM.'
)
```

---

#### Step 1.3: Pass flush_interval to StreamingParquetWriter

**File:** `scripts/08_transform_library_v2.py`

**Location:** Around line 840 (in cmd_apply function, where writer is created)

**Find This Code (around line 840):**
```python
# Initialize streaming writer
writer = StreamingParquetWriter(output_path, schema)
```

**REPLACE WITH:**
```python
# Initialize streaming writer with flush interval
writer = StreamingParquetWriter(
    output_path,
    schema,
    flush_interval=args.flush_interval  # NEW: Pass flush_interval parameter
)
logger.info(f"StreamingParquetWriter initialized with flush_interval={args.flush_interval:,}")
```

---

#### Step 1.4: Verification Steps

**After making changes, verify:**

1. **Syntax Check:**
```bash
python -m py_compile scripts/08_transform_library_v2.py
```

2. **Help Display:**
```bash
python scripts/08_transform_library_v2.py apply --help
# Should show new --flush-interval parameter
```

3. **Small Test Run:**
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/lipid-1X/products.parquet \
  --output test_transform_output.parquet \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 4 \
  --flush-interval 1000
```

**Expected Output:**
- "ParquetWriter initialized (flush_interval=1,000)"
- "Flush #1: X products written" messages
- "Final flush: Y products"
- "StreamingParquetWriter closed: Z total products, N flushes"
- No errors

4. **Verify Output File:**
```bash
python -c "import pyarrow.parquet as pq; print(f'{pq.read_table(\"test_transform_output.parquet\").num_rows:,} products')"
```

---

### Task 2: Generate Transform Library for All Classes

**Goal:** Apply all transformation rules to the complete halogenated library (67.4M products)

**Priority:** 🟡 **HIGH** - Main deliverable for this phase

**Estimated Time:** Variable (depends on rules and input size, 1-7 days)

---

#### Step 2.1: Identify Transform Rules to Apply

**File:** `configs/transforms.yaml`

**Available Transform Rules:**

**Phenolic OH Modifications (5 rules):**
1. `FG_PHENOL_OH__OH__TO__OMe` - OH → OMe
2. `FG_PHENOL_OH__OH__TO__NH2` - OH → NH2
3. `FG_PHENOL_OH__OH__TO__CH2CH2NH2` - OH → CH2CH2NH2 (NEW)
4. `FG_PHENOL_OH__OH__TO__CH2CH2COOH` - OH → CH2CH2COOH (NEW)

**Aromatic Amine Modifications (8 rules):**
5. `FG_AMINE_AR__NH2__TO__F` - Ar-NH2 → F
6. `FG_AMINE_AR__NH2__TO__Cl` - Ar-NH2 → Cl
7. `FG_AMINE_AR__NH2__TO__Br` - Ar-NH2 → Br
8. `FG_AMINE_AR__NH2__TO__I` - Ar-NH2 → I
9. `FG_AMINE_AR__NH2__TO__OH` - Ar-NH2 → OH
10. `FG_AMINE_AR__NH2__TO__NHMe` - Ar-NH2 → NHMe
11. `FG_AMINE_AR__NH2__TO__NMe2` - Ar-NH2 → NMe2
12. `FG_AMINE_AR__NH2__TO__NHCHO` - Ar-NH2 → NHCHO
13. `FG_AMINE_AR__NH2__TO__NHAc` - Ar-NH2 → NHAc

**Aliphatic Amine Modifications (3 rules):**
14. `FG_AMINE_ALIPH__NH2__TO__NHMe`
15. `FG_AMINE_ALIPH__NH2__TO__NMe2`
16. `FG_AMINE_ALIPH__NH2__TO__NHAc`

**Carboxylic Acid Modifications (3 rules):**
17. `FG_CARBOXYL__COOH__TO__COOMe`
18. `FG_CARBOXYL__COOH__TO__COOEt`
19. `FG_CARBOXYL__COOH__TO__COCl`

**Total:** 19 transformation rules

---

#### Step 2.2: Determine Input-Rule Mapping Strategy

**Strategy:** Apply rules to relevant K=1+K=2 subsets based on scope_classes

**Recommended Approach:**

**Option 1: Per-Class, Per-Rule (RECOMMENDED)**
- Apply each rule to each relevant class separately
- Easier to parallelize and monitor
- Can resume individual jobs if failures occur

**Example for polyphenol class:**
```bash
# For each rule where scope_classes includes "polyphenol"
for rule in FG_PHENOL_OH__OH__TO__OMe FG_PHENOL_OH__OH__TO__NH2 ...; do
  for k in 1X 2X; do
    python scripts/08_transform_library_v2.py apply \
      --input data/output/nplike_v2/polyphenol-${k}/products.parquet \
      --output data/output/transforms/polyphenol-${k}_${rule}.parquet \
      --xf-name ${rule} \
      --workers 16 \
      --flush-interval 50000 \
      --batch-size 50000
  done
done
```

**Option 2: Merged Input, Per-Rule**
- Merge all K=1+K=2 files first
- Apply each rule to merged file
- Simpler workflow but harder to resume

---

#### Step 2.3: Recommended Execution Order

**Batch 1: Small Classes (Fast, ~1-2 hours total)**
- lipid (213K products)
- aa_peptide (613K products)
- polysaccharide (314K products)
- other (1.5M products)

**Batch 2: Medium Classes (~6-12 hours total)**
- alkaloid (8.8M products)
- polyphenol (14.3M products)

**Batch 3: Large Class (~24-48 hours)**
- terpenoid (41.7M products)

**Rationale:**
- Start with small classes to validate pipeline
- Catch any issues before committing to large jobs
- Can parallelize multiple small jobs

---

#### Step 2.4: Sample Command Template

**For a Single Transform Job:**
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/CLASSNAME-KX/products.parquet \
  --output data/output/transforms/CLASSNAME-KX_RULENAME.parquet \
  --xf-name RULENAME \
  --workers 16 \
  --flush-interval 50000 \
  --batch-size 50000 \
  --resume  # Optional: resume if interrupted
```

**Key Parameters:**
- `--workers 16`: Full CPU utilization
- `--flush-interval 50000`: Balanced memory/performance (50K products)
- `--batch-size 50000`: Input chunk size (50K parents)
- `--resume`: Safe to use if job interrupted

---

#### Step 2.5: Automation Script (Optional but Recommended)

**Create:** `scripts/run_all_transforms.sh`

```bash
#!/bin/bash
# Batch transform pipeline orchestrator
# Usage: bash scripts/run_all_transforms.sh

set -e  # Exit on error

WORKERS=16
FLUSH_INTERVAL=50000
BATCH_SIZE=50000

# Define class-rule mappings (based on scope_classes in transforms.yaml)
declare -A CLASS_RULES=(
  ["polyphenol"]="FG_PHENOL_OH__OH__TO__OMe FG_PHENOL_OH__OH__TO__NH2 FG_PHENOL_OH__OH__TO__CH2CH2NH2 FG_PHENOL_OH__OH__TO__CH2CH2COOH FG_AMINE_AR__NH2__TO__F FG_CARBOXYL__COOH__TO__COOMe"
  ["terpenoid"]="FG_PHENOL_OH__OH__TO__OMe FG_PHENOL_OH__OH__TO__NH2 FG_PHENOL_OH__OH__TO__CH2CH2NH2 FG_PHENOL_OH__OH__TO__CH2CH2COOH FG_AMINE_ALIPH__NH2__TO__NHMe FG_CARBOXYL__COOH__TO__COOMe"
  ["alkaloid"]="FG_AMINE_AR__NH2__TO__F FG_AMINE_AR__NH2__TO__OH FG_AMINE_ALIPH__NH2__TO__NHMe"
  ["lipid"]="FG_AMINE_ALIPH__NH2__TO__NHMe FG_CARBOXYL__COOH__TO__COOMe"
  ["aa_peptide"]="FG_AMINE_AR__NH2__TO__F FG_AMINE_ALIPH__NH2__TO__NHMe FG_CARBOXYL__COOH__TO__COOMe"
)

echo "========================================================================"
echo "Transform Pipeline Orchestrator"
echo "========================================================================"

# Iterate over classes and rules
for CLASS in "${!CLASS_RULES[@]}"; do
  for K in 1X 2X; do
    INPUT="data/output/nplike_v2/${CLASS}-${K}/products.parquet"

    # Check if input exists
    if [ ! -f "$INPUT" ]; then
      echo "SKIP: $INPUT not found"
      continue
    fi

    echo ""
    echo "========================================================================"
    echo "Processing: $CLASS-$K"
    echo "========================================================================"

    for RULE in ${CLASS_RULES[$CLASS]}; do
      OUTPUT="data/output/transforms/${CLASS}-${K}_${RULE}.parquet"

      # Skip if already exists (resumable)
      if [ -f "$OUTPUT" ]; then
        echo "SKIP: $OUTPUT already exists"
        continue
      fi

      echo ""
      echo "Transform: $RULE"
      echo "  Input:  $INPUT"
      echo "  Output: $OUTPUT"

      python scripts/08_transform_library_v2.py apply \
        --input "$INPUT" \
        --output "$OUTPUT" \
        --xf-name "$RULE" \
        --workers $WORKERS \
        --flush-interval $FLUSH_INTERVAL \
        --batch-size $BATCH_SIZE

      if [ $? -eq 0 ]; then
        # Count products
        PRODUCT_COUNT=$(python -c "import pyarrow.parquet as pq; print(pq.read_table('$OUTPUT').num_rows)")
        echo "SUCCESS: $PRODUCT_COUNT products generated"
      else
        echo "ERROR: Transform failed for $CLASS-$K + $RULE"
        exit 1
      fi
    done
  done
done

echo ""
echo "========================================================================"
echo "Transform Pipeline Complete"
echo "========================================================================"
```

**Usage:**
```bash
# Make executable
chmod +x scripts/run_all_transforms.sh

# Run (can take 1-7 days depending on parallelization)
nohup bash scripts/run_all_transforms.sh > transform_pipeline.log 2>&1 &

# Monitor
tail -f transform_pipeline.log
```

---

#### Step 2.6: Monitoring and Validation

**During Execution:**

1. **Monitor Progress:**
```bash
tail -f transform_pipeline.log
```

2. **Check Memory Usage:**
```bash
top -b -n 1 | grep python
```

3. **Count Completed Files:**
```bash
ls -lh data/output/transforms/ | wc -l
```

**After Completion:**

4. **Validate All Outputs:**
```bash
python << 'EOF'
import pyarrow.parquet as pq
from pathlib import Path

transform_dir = Path('data/output/transforms')
total_products = 0

print("Transform Library Summary:")
print("="*80)

for file in sorted(transform_dir.glob('*.parquet')):
    count = pq.read_table(file).num_rows
    total_products += count
    print(f"{file.name:<60} {count:>12,}")

print("="*80)
print(f"{'TOTAL TRANSFORM PRODUCTS':<60} {total_products:>12,}")
print("="*80)
EOF
```

5. **Check for Failures:**
```bash
grep -i "error\|failed" transform_pipeline.log
```

---

### Task 3: Post-Transform Quality Control (Optional)

**Goal:** Validate transform library quality

**Steps:**

1. **Check for Duplicates Across Rules:**
```bash
python << 'EOF'
import pyarrow.parquet as pq
from pathlib import Path
from collections import Counter

transform_dir = Path('data/output/transforms')
all_inchikeys = []

for file in transform_dir.glob('*.parquet'):
    df = pq.read_table(file, columns=['inchikey']).to_pandas()
    all_inchikeys.extend(df['inchikey'].tolist())

total = len(all_inchikeys)
unique = len(set(all_inchikeys))
duplicates = total - unique

print(f"Total products: {total:,}")
print(f"Unique InChIKeys: {unique:,}")
print(f"Duplicates: {duplicates:,} ({duplicates/total*100:.2f}%)")
EOF
```

2. **Sample Quality Check:**
```bash
python << 'EOF'
import pyarrow.parquet as pq
import pandas as pd

# Pick a random transform file
file = 'data/output/transforms/polyphenol-1X_FG_PHENOL_OH__OH__TO__OMe.parquet'
df = pq.read_table(file).to_pandas()

print(f"File: {file}")
print(f"Products: {len(df):,}")
print(f"\nSample products:")
print(df[['smiles', 'xf_rule', 'source_smiles']].head(10))

# Check success rate
success_rate = df['xf_success'].sum() / len(df) * 100
print(f"\nTransform success rate: {success_rate:.2f}%")
EOF
```

---

## PART C: CRITICAL FILES REFERENCE

### Modified Files (This Session)

1. **`src/halogenator/enumerate_k1.py`**
   - Lines 103, 343-427, 456, 507: Sugar_mask fix
   - Status: ✅ Committed (2ac15dd)

2. **`src/halogenator/parallel_enum.py`**
   - Lines 433-453: Schema cast on every flush
   - Status: ✅ Production-validated (788 flushes successful)

3. **`configs/transforms.yaml`**
   - Lines 37-57: Added 2 new phenolic OH transform rules
   - Status: ✅ Validated

### Files to Modify (Next Session)

4. **`scripts/08_transform_library_v2.py`**
   - Lines 301-328: StreamingParquetWriter class (REPLACE ENTIRE CLASS)
   - Line ~840: StreamingParquetWriter initialization (ADD flush_interval parameter)
   - Line ~1040: Argparse (ADD --flush-interval argument)
   - Status: ⏳ **PENDING - CRITICAL**

### Input Data Files

5. **Halogenated Library (K=1):**
   - `data/output/nplike_v2/aa_peptide-1X/products.parquet` (30,564)
   - `data/output/nplike_v2/alkaloid-1X/products.parquet` (322,280)
   - `data/output/nplike_v2/lipid-1X/products.parquet` (33,084)
   - `data/output/nplike_v2/polyphenol-1X/products.parquet` (518,852)
   - `data/output/nplike_v2/terpenoid-1X/products.parquet` (1,528,400)
   - `data/output/nplike_v2/polysaccharide-1X/products.parquet` (20,524)
   - `data/output/nplike_v2/other-1X/products.parquet` (97,012)

6. **Halogenated Library (K=2):**
   - `data/output/nplike_v2/aa_peptide-2X/products.parquet` (582,730)
   - `data/output/nplike_v2/alkaloid-2X/products.parquet` (8,461,937)
   - `data/output/nplike_v2/lipid-2X/products.parquet` (180,224)
   - `data/output/nplike_v2/polyphenol-2X/products.parquet` (13,790,820)
   - `data/output/nplike_v2/terpenoid-2X/products.parquet` (40,170,816)
   - `data/output/nplike_v2/polysaccharide-2X/products.parquet` (293,616)
   - `data/output/nplike_v2/other-2X/products.parquet` (1,419,941)

### Output Directory

7. **Transform Library Output:**
   - `data/output/transforms/` (CREATE THIS DIRECTORY)
   - Expected files: ~250+ parquet files (19 rules × 7 classes × 2 K-values)

---

## PART D: TROUBLESHOOTING GUIDE

### Issue 1: OOM During Transform

**Symptoms:**
- Python process killed
- "Killed" message in logs
- No output file created

**Solution:**
1. Reduce `--workers` (try 8 or 4)
2. Reduce `--flush-interval` (try 10000)
3. Reduce `--batch-size` (try 25000)

**Example:**
```bash
# Low-memory mode
python scripts/08_transform_library_v2.py apply \
  --workers 4 \
  --flush-interval 10000 \
  --batch-size 25000 \
  ...
```

---

### Issue 2: Schema Mismatch Error

**Symptoms:**
```
ValueError: Table schema does not match schema used to create file
```

**Solution:**
Verify that StreamingParquetWriter includes the schema cast code:
```python
# In _flush() method, after "if self.writer is None:" block
else:
    table = table.cast(self.schema)  # This line MUST be present
```

If error persists, use the column-by-column fallback (already in code above).

---

### Issue 3: Transform Rule Not Found

**Symptoms:**
```
ERROR: Transform rule 'XYZ' not found in config
```

**Solution:**
1. Check rule name spelling (case-sensitive)
2. Verify rule exists in `configs/transforms.yaml`
3. Check YAML syntax (no tabs, proper indentation)

**Validation:**
```bash
python -c "import yaml; print(yaml.safe_load(open('configs/transforms.yaml'))['transforms'])"
```

---

### Issue 4: No Products Generated

**Symptoms:**
- Job completes successfully
- Output file created but has 0 products

**Possible Causes:**
1. Input SMILES don't match query_smarts
2. Transform rule not applicable to this class
3. All reactions failed (check logs for 'xf_success' rate)

**Debug:**
```bash
# Check input SMILES sample
python -c "import pyarrow.parquet as pq; df = pq.read_table('INPUT.parquet', columns=['smiles']).to_pandas(); print(df.head())"

# Check transform rule
grep -A 10 "name: RULENAME" configs/transforms.yaml
```

---

## PART E: SUCCESS CRITERIA

### Minimum Viable Product (MVP)

✅ **Completed When:**
1. StreamingParquetWriter improvements implemented
2. At least one transform rule applied successfully to one class
3. Output validated (products > 0, no errors)

### Full Success

✅ **Completed When:**
1. All 19 transform rules applied to all applicable classes
2. No schema errors, no OOM crashes
3. Transform library validated (duplicates checked, quality sampled)
4. Total transform products: estimated 500K-2M

---

## PART F: EXPECTED OUTCOMES

### Transform Library Size Estimates

**Conservative Estimate (per rule, per class):**
- Rules match ~1-10% of input molecules
- Each match generates 1 product (no k>1 for transforms)

**Example for polyphenol-2X (13.79M products):**
- FG_PHENOL_OH rules: ~5% match → ~690K products/rule
- 4 phenolic OH rules → ~2.76M transform products

**Total Estimate:**
- Small classes: ~50K products/rule
- Medium classes: ~200K-500K products/rule
- Large classes: ~500K-2M products/rule
- **Grand Total: 5-15M transform products** (depends on hit rates)

---

## PART G: NEXT SESSION STARTUP INSTRUCTIONS

### Quick Start Checklist

1. ✅ **Verify environment:**
```bash
cd E:/Projects/halogenator
conda activate halo-p0  # Or your environment
python --version  # Should be 3.8+
```

2. ✅ **Check input data:**
```bash
ls -lh data/output/nplike_v2/*-1X/products.parquet
ls -lh data/output/nplike_v2/*-2X/products.parquet
# Should see 14 files (7 classes × 2 K-values)
```

3. ✅ **Read this report:**
```bash
cat SESSION_HANDOFF_TRANSFORM_PIPELINE_2025-12-13.md
```

4. ✅ **Implement StreamingParquetWriter improvements** (Task 1)
   - Follow Step 1.1 exactly (copy-paste the new class)
   - Follow Step 1.2 (add argparse parameter)
   - Follow Step 1.3 (pass parameter to writer)
   - Follow Step 1.4 (verification tests)

5. ✅ **Run test transform** (small dataset)
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/lipid-1X/products.parquet \
  --output test_transform_lipid.parquet \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 4 \
  --flush-interval 1000
```

6. ✅ **If test passes, proceed to full pipeline** (Task 2)

---

## PART H: DOCUMENTATION GENERATED THIS SESSION

1. **`SUGAR_MASK_FIX_VALIDATION_REPORT.md`**
   - Complete K=1 validation details
   - Before/after comparison
   - Product distribution by rule

2. **`FINAL_K1_K2_SUMMARY_2025-12-11.md`**
   - Complete K=1+K=2 statistics
   - All 7 classes, 67.4M products
   - Performance metrics

3. **`SCHEMA_FIX_AND_TERPENOID_K2_SUCCESS.md`**
   - Schema bug diagnosis
   - Solution implementation
   - Terpenoid K=2 validation (788 flushes)

4. **`PARALLEL_MECHANISMS_COMPARISON.md`**
   - Detailed comparison of 04 vs 08 mechanisms
   - Implementation recommendations
   - Risk analysis

5. **`SESSION_HANDOFF_TRANSFORM_PIPELINE_2025-12-13.md`** (this file)
   - Complete task list
   - Detailed implementation guide
   - Troubleshooting reference

---

## PART I: CONTEXT FOR AI ASSISTANT

### Key Design Decisions

1. **Why flush_interval instead of just batch_size?**
   - batch_size controls INPUT chunks (how many parents to read)
   - flush_interval controls OUTPUT buffer (how many products to accumulate)
   - Independent tuning allows optimization for different bottlenecks

2. **Why copy parallel_enum.py's pattern?**
   - Proven in production (788 consecutive flushes on 40M products)
   - Already handles schema edge cases
   - Consistent UX across codebase

3. **Why optional psutil?**
   - Not all systems have psutil
   - Graceful degradation (works without it)
   - Only needed for memory-constrained environments

4. **Why schema cast on EVERY flush?**
   - PyArrow infers null type for all-null columns
   - Subsequent batches may have real values
   - Without cast, crashes after 1st or 2nd flush
   - Terpenoid K=2 bug proved this is critical

### Common Pitfalls to Avoid

1. ❌ **Don't skip verification tests** (Step 1.4)
   - Small test catches 90% of bugs
   - Much faster than debugging full pipeline

2. ❌ **Don't use small flush_interval on large jobs**
   - Too many I/O operations = slow
   - Recommended: 10K-50K for production

3. ❌ **Don't forget to call writer.close()**
   - Final flush is in close()
   - Missing close() = lost products

4. ❌ **Don't run transforms before improving 08**
   - High risk of OOM or schema crash
   - 30 min to fix, saves days of reruns

---

## PART J: FINAL CHECKLIST FOR NEXT SESSION

### Before Starting Coding

- [ ] Read Part B (Pending Tasks) completely
- [ ] Review StreamingParquetWriter replacement code (Step 1.1)
- [ ] Check that all input files exist (Part C, items 5-6)
- [ ] Create output directory: `mkdir -p data/output/transforms`

### During Implementation

- [ ] Replace StreamingParquetWriter class (Step 1.1)
- [ ] Add --flush-interval argument (Step 1.2)
- [ ] Pass flush_interval to writer (Step 1.3)
- [ ] Run syntax check (Step 1.4.1)
- [ ] Run help check (Step 1.4.2)
- [ ] Run small test (Step 1.4.3)
- [ ] Verify test output (Step 1.4.4)

### After Implementation

- [ ] Decide on batch strategy (Step 2.2)
- [ ] Choose execution order (Step 2.3)
- [ ] Run small class first (lipid or aa_peptide)
- [ ] Monitor for errors
- [ ] If successful, proceed to larger classes

### Quality Gates

- [ ] No OOM crashes
- [ ] No schema mismatch errors
- [ ] Products > 0 for each transform
- [ ] Log shows "Flush #N" messages
- [ ] Final "StreamingParquetWriter closed" message

---

**End of Handoff Report**

**Session Completed By:** Claude Sonnet 4.5
**Date:** 2025-12-13
**Next Session Task:** Implement StreamingParquetWriter improvements (1 hour) → Run transform pipeline (1-7 days)
**Status:** ✅ Ready for implementation
