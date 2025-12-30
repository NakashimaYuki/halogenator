# Parallel Processing & OOM Protection Comparison

**Scripts:** `04_enum_halogen_all_classes.py` vs `08_transform_library_v2.py`

---

## Current Implementation Analysis

### 04_enum_halogen_all_classes.py (Orchestrator)

**Architecture:**
- Wrapper script that calls `halogenator enum-parquet` CLI
- Actual parallelism handled by `parallel_enum.py`

**Key Parameters:**
```python
--workers 16              # Parallel worker processes
--flush-interval 50000    # Products before disk flush
--batch-size 5000         # Parents per batch
--rdkit-threads 8         # RDKit internal threads
```

**OOM Protection (in parallel_enum.py):**
1. ✅ **Product Buffer + Flush Interval**
   - Accumulates products in memory buffer
   - Flush to disk every N products (configurable)
   - Line 298: `_flush_to_disk()` method

2. ✅ **Memory Monitoring** (optional)
   - psutil-based memory usage tracking
   - Auto-flush when memory > 80% threshold
   - Line ~250: Memory check in enumeration loop

3. ✅ **Adaptive Flushing**
   - Regular interval flushes (50K products)
   - Memory-triggered emergency flushes
   - Prevents OOM on large jobs

4. ✅ **Schema Consistency** (recently fixed)
   - Force cast to established schema on every flush
   - Handles null-type inference issues
   - Lines 405-453: Schema enforcement

---

### 08_transform_library_v2.py (Transform Library)

**Architecture:**
- Direct ProcessPoolExecutor parallelism
- Streaming Parquet writer

**Key Parameters:**
```python
--workers 6               # Parallel worker processes (default)
--batch-size 50000        # Input batch size
# ❌ No flush-interval parameter
# ❌ No memory monitoring
```

**Current Mechanisms:**

1. ✅ **Streaming Input**
   - `iter_batches(batch_size)` avoids loading full dataset
   - Line 879: Batch-by-batch processing
   - Good: Prevents input-side OOM

2. ✅ **Immediate Output**
   - `writer.write_batch()` writes each batch immediately
   - Line 929: Direct write after deduplication
   - Good: No output accumulation

3. ❌ **No Product Buffer**
   - StreamingParquetWriter writes immediately
   - Lines 310-323: No buffering mechanism
   - Risk: Many small I/O operations

4. ❌ **No Memory Monitoring**
   - No psutil checks
   - No adaptive behavior based on memory pressure

5. ❌ **No Flush Control**
   - Fixed batch size, no flush interval
   - Cannot tune output buffer behavior

---

## Key Differences

| Feature | 04 (parallel_enum.py) | 08 (transform_v2.py) | Winner |
|---------|----------------------|---------------------|--------|
| **Product Buffer** | ✅ Yes (flush_interval) | ❌ No (immediate write) | 04 |
| **Memory Monitoring** | ✅ psutil-based | ❌ None | 04 |
| **Adaptive Flushing** | ✅ Memory + interval | ❌ Fixed batch | 04 |
| **Schema Enforcement** | ✅ Cast on every flush | ⚠️ Only on first write | 04 |
| **Streaming Input** | ✅ Batch processing | ✅ iter_batches | Tie |
| **Parallel Architecture** | ✅ multiprocessing (C level) | ✅ ProcessPoolExecutor | Tie |
| **Deduplication** | ✅ InChIKey (enumerate_k) | ✅ Canonical SMILES + SQL | 08 |

---

## Risks in Current 08 Implementation

### 1. Lack of Output Buffering
**Problem:**
- Writes to Parquet on every batch
- Many small write operations → I/O overhead
- No control over memory vs. performance trade-off

**Impact:**
- Slower for large outputs (64M products)
- Cannot tune for performance

### 2. No Memory Monitoring
**Problem:**
- Cannot detect memory pressure
- No emergency flush mechanism
- Workers may OOM if products accumulate

**Impact:**
- Risk of OOM on memory-constrained systems
- No graceful degradation under memory pressure

### 3. Fixed Batch Size
**Problem:**
- batch_size controls INPUT chunks, not output buffering
- No separate control for output flushing

**Impact:**
- Cannot independently tune input/output behavior
- Less flexible than enum pipeline

---

## Recommended Improvements

### Priority 1: Add Flush Interval to StreamingParquetWriter

**Goal:** Buffer products and flush at configurable intervals

**Implementation:**
```python
class StreamingParquetWriter:
    def __init__(self, output_path: Path, schema: pa.Schema, flush_interval: int = 10000):
        self.output_path = output_path
        self.writer = None
        self.schema = schema
        self.total_written = 0
        self.flush_interval = flush_interval
        self.buffer = []  # NEW: Product buffer

    def write_batch(self, products: List[Dict]):
        """Add products to buffer, flush when interval reached."""
        self.buffer.extend(products)

        if len(self.buffer) >= self.flush_interval:
            self._flush()

    def _flush(self):
        """Flush buffer to disk."""
        if not self.buffer:
            return

        table = pa.Table.from_pylist(self.buffer, schema=self.schema)

        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)

        # Cast to schema (like parallel_enum.py fix)
        table = table.cast(self.schema)

        self.writer.write_table(table)
        self.total_written += len(self.buffer)
        self.buffer.clear()

    def close(self):
        """Flush remaining buffer before closing."""
        self._flush()
        if self.writer:
            self.writer.close()
```

**Benefits:**
- Reduces I/O operations (fewer writes)
- Tunable performance (adjust flush_interval)
- Consistent with enum pipeline

---

### Priority 2: Add Memory Monitoring (Optional)

**Goal:** Auto-flush on memory pressure (like parallel_enum.py)

**Implementation:**
```python
import psutil

class StreamingParquetWriter:
    def __init__(self, ..., memory_threshold: float = 0.8):
        # ... existing init ...
        self.memory_threshold = memory_threshold
        self.memory_monitor_enabled = self._check_psutil()

    def _check_psutil(self) -> bool:
        """Check if psutil available for memory monitoring."""
        try:
            import psutil
            return True
        except ImportError:
            logger.warning("psutil not available - memory monitoring disabled")
            return False

    def _check_memory(self) -> float:
        """Get current memory usage percentage."""
        if not self.memory_monitor_enabled:
            return 0.0
        try:
            return psutil.virtual_memory().percent / 100.0
        except:
            return 0.0

    def write_batch(self, products: List[Dict]):
        """Add products to buffer, flush based on interval OR memory."""
        self.buffer.extend(products)

        # Check flush conditions
        mem_usage = self._check_memory()
        should_flush = (
            len(self.buffer) >= self.flush_interval or
            (self.memory_monitor_enabled and mem_usage > self.memory_threshold)
        )

        if should_flush:
            if mem_usage > self.memory_threshold:
                logger.warning(f"Memory-triggered flush: {mem_usage:.1%} > {self.memory_threshold:.1%}")
            self._flush()
```

**Benefits:**
- Prevents OOM on large jobs
- Graceful degradation under memory pressure
- Optional (works without psutil)

---

### Priority 3: Add flush_interval CLI Parameter

**Goal:** Let users tune output buffer size

**Implementation:**
```python
parser_apply.add_argument(
    '--flush-interval',
    type=int,
    default=10000,
    help='Products to buffer before flushing to disk (default: 10000)'
)

# In cmd_apply:
writer = StreamingParquetWriter(
    output_path,
    schema,
    flush_interval=args.flush_interval  # NEW parameter
)
```

**Benefits:**
- User control over memory/performance trade-off
- Consistent with enum pipeline UX
- Easy to tune for different systems

---

### Priority 4: Schema Enforcement on Every Flush

**Goal:** Prevent schema mismatch errors (like terpenoid k=2 bug)

**Implementation:**
```python
def _flush(self):
    """Flush buffer to disk with schema enforcement."""
    if not self.buffer:
        return

    table = pa.Table.from_pylist(self.buffer, schema=self.schema)

    if self.writer is None:
        # First flush: establish schema
        self.writer = pq.ParquetWriter(self.output_path, self.schema)
    else:
        # Subsequent flushes: cast to established schema
        # CRITICAL: Prevents null-type inference issues
        try:
            table = table.cast(self.schema)
        except Exception as e:
            logger.warning(f"Schema cast failed: {e}")
            # Fallback: column-by-column cast
            # (copy logic from parallel_enum.py lines 437-453)

    self.writer.write_table(table)
    self.total_written += len(self.buffer)
    self.buffer.clear()
```

**Benefits:**
- Prevents schema mismatch crashes
- Robust to all-null batches
- Proven solution from enum pipeline

---

## Implementation Plan

### Step 1: Enhance StreamingParquetWriter (CRITICAL)
- Add product buffer
- Add flush_interval parameter
- Add _flush() method
- Update close() to flush remaining buffer

**Estimated Time:** 30 minutes
**Risk:** Low (well-tested pattern from parallel_enum.py)

### Step 2: Add Memory Monitoring (RECOMMENDED)
- Optional psutil import
- _check_memory() method
- Memory-triggered flush logic

**Estimated Time:** 20 minutes
**Risk:** Low (gracefully degrades without psutil)

### Step 3: Add CLI Parameters (EASY)
- --flush-interval argument
- Pass to StreamingParquetWriter

**Estimated Time:** 5 minutes
**Risk:** Minimal

### Step 4: Schema Enforcement (CRITICAL)
- Add cast() on every flush (except first)
- Fallback to column-by-column cast
- Copy logic from parallel_enum.py

**Estimated Time:** 15 minutes
**Risk:** Low (proven solution)

---

## Testing Strategy

### Test 1: Small Dataset (Baseline)
```bash
python 08_transform_library_v2.py apply \
  --input data/output/nplike_v2/lipid-1X/products.parquet \
  --output test_output.parquet \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 4 \
  --flush-interval 1000  # NEW
```

**Expected:** Works like before, but with buffering

### Test 2: Large Dataset (Stress Test)
```bash
python 08_transform_library_v2.py apply \
  --input data/output/nplike_v2/terpenoid-2X/products.parquet \
  --output test_terpenoid_xf.parquet \
  --xf-name FG_PHENOL_OH__OH__TO__CH2CH2NH2 \
  --workers 16 \
  --flush-interval 50000  # Same as enum
```

**Expected:** No OOM, completes successfully

### Test 3: Memory Pressure
```bash
# Artificially constrain memory and run
python 08_transform_library_v2.py apply \
  --input data/output/nplike_v2/polyphenol-2X/products.parquet \
  --output test_polyphenol_xf.parquet \
  --xf-name FG_PHENOL_OH__OH__TO__CH2CH2COOH \
  --workers 8 \
  --flush-interval 100000  # Large buffer
```

**Expected:** Memory-triggered flushes prevent OOM

---

## Conclusion

### Can We Port the Mechanisms?

✅ **YES - Highly Compatible**

The mechanisms from `04_enum_halogen_all_classes.py` (via `parallel_enum.py`) can be successfully ported to `08_transform_library_v2.py`:

1. **Product Buffering + Flush Interval**
   - Direct port from parallel_enum.py
   - Proven, production-tested
   - Low implementation risk

2. **Memory Monitoring**
   - Optional psutil dependency (same as enum)
   - Same thresholds (80%)
   - Graceful degradation

3. **Schema Enforcement**
   - Copy exact logic from parallel_enum.py fix
   - Prevents terpenoid-style crashes
   - No modifications needed

### Recommendations

**CRITICAL (Do Before Running):**
1. Add flush_interval to StreamingParquetWriter ✅
2. Add schema cast on every flush ✅

**RECOMMENDED (Improves Robustness):**
3. Add memory monitoring ✅
4. Add --flush-interval CLI parameter ✅

**Time to Implement:** ~1 hour total
**Risk:** Low (proven patterns)
**Benefit:** High (prevents OOM, improves performance)

---

**Ready to Implement:** These improvements are straightforward ports of proven mechanisms. The code patterns are nearly identical, just need to be adapted to the StreamingParquetWriter class.
