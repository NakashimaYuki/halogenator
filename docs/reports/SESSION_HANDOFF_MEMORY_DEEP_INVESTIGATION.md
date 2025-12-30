# Session Handoff: Transform Pipeline Memory Crisis - Deep Investigation Required

**Date:** 2025-12-26
**Status:** 🔴 CRITICAL - Partial solutions ineffective, requires complete redesign
**Priority:** URGENT - Blocking full pipeline execution
**Estimated Effort:** 8-12 hours deep analysis + implementation

---

## Executive Summary

**Problem:** Transform pipeline for large datasets (polyphenol-2X, terpenoid-2X) consistently reaches 90-95% memory usage and crashes or corrupts outputs, despite multiple attempted fixes.

**Attempts Made:**
1. ❌ **max_in_flight control** - Improved but insufficient
2. ❌ **Bloom filter deduplication** - Helped dedup but final stage still OOMs
3. ❌ Both combined - Task completed but output corrupted, memory hit 94.2%

**Core Issue:** Memory accumulates gradually from multiple sources that compound over time. Need complete data flow analysis to identify ALL accumulation points.

**Next Session Goal:**
1. Trace complete memory data flow from input to output
2. Identify every single point where memory accumulates
3. Design holistic solution addressing all accumulation points
4. Implement and validate on polyphenol-2X

---

## Part A: Background - What We Know (BRIEF)

### Problem History

**Original SQLite Deduplicator Issue:**
- `SqliteDeduplicator.seen_in_memory` grows to 9.9GB for 60M products
- `dedup.db` file reaches 42GB on disk (~15GB in OS cache)
- Combined with other sources → 28-33GB total → OOM on 32GB system

**Fix #1: max_in_flight Control**
- **What:** Limit concurrent futures to workers × 2 (default 32 → 16)
- **File:** `scripts/08_transform_library_v2.py` lines 1170-1249
- **Result:** Reduced futures queue from 7.7GB → 3.8GB
- **Limitation:** Still hit 95%+ memory after 23 hours, dedup.db grew to 42GB

**Fix #2: Bloom Filter**
- **What:** Replace SQLite seen_in_memory with probabilistic Bloom filter
- **File:** `src/halogenator/bloom_dedup.py` (new)
- **Result:** Dedup memory reduced from 9.9GB → 69MB (99.3% saving)
- **Limitation:** Task completed but:
  - Memory peaked at 94.2% in final stage (batch 260-276)
  - Output file corrupted (9GB parquet, no footer magic bytes)
  - Still has critical memory spikes

### Test Results (polyphenol-2X, Batch 260/276, 16h runtime)

**Memory Timeline:**
```
T+0h:     27.8%  (startup)
T+1h:     45.9%  (workers initialized, Bloom filter working)
T+3h:     57.5%  (stable, looked promising)
T+14h:    68.5%  (batch 260/276, still acceptable)
T+14.5h:  72.3%  ("All batches submitted. Processing remaining 15...")
T+15h:    78.8% → 90.3% → 94.2%  (CRITICAL SPIKE!)
T+15-16h: 87-94% fluctuation (corrupted output produced)
```

**Evidence Files:**
- Log: `E:/Projects/halogenator/polyphenol_2x_BLOOM_TEST.log`
- Memory monitor: `E:/Projects/halogenator/bloom_memory_monitor.log`
- Output (corrupted): `data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe/`
  - `products.parquet`: 9.0GB, **CORRUPTED** (no magic bytes in footer)
  - `dedup_bloom.pkl`: 69MB, valid

**Key Observation:** Memory growth in final stage (last 15 batches) was catastrophic: 72% → 94% in 30 minutes.

---

## Part B:待完成任务 - DETAILED ACTION PLAN

### Task 1: Complete Memory Data Flow Tracing

**Objective:** Map every single byte that enters and stays in memory during pipeline execution.

**Why This Matters:**
Previous fixes targeted individual components (dedup, futures queue) but memory still accumulates. We need to find ALL sources, not just the obvious ones.

#### Subtask 1.1: Instrument Code with Memory Profiling

**Goal:** Add detailed memory tracking at every major operation.

**Implementation Steps:**

1. **Install memory profiler:**
   ```bash
   pip install memory_profiler
   pip install pympler  # For detailed object tracking
   ```

2. **Create diagnostic script:** `scripts/diagnose_memory_flow.py`
   ```python
   """
   Memory flow diagnostic tool.

   Tracks memory usage at each step of the transform pipeline.
   """

   import gc
   import psutil
   import tracemalloc
   from pympler import asizeof, muppy, summary

   def get_memory_snapshot(label):
       """Take a memory snapshot and log details."""
       gc.collect()  # Force GC first

       # System memory
       vm = psutil.virtual_memory()

       # Process memory
       process = psutil.Process()
       mem_info = process.memory_info()

       # Python object tracking
       all_objects = muppy.get_objects()
       sum_obj = summary.summarize(all_objects)

       snapshot = {
           'label': label,
           'system_percent': vm.percent,
           'system_used_gb': vm.used / 1024**3,
           'process_rss_gb': mem_info.rss / 1024**3,
           'python_objects': sum_obj[:10],  # Top 10 object types
       }

       return snapshot

   def compare_snapshots(before, after):
       """Compare two snapshots and identify growth."""
       delta_system = after['system_percent'] - before['system_percent']
       delta_process = after['process_rss_gb'] - before['process_rss_gb']

       print(f"\n[MEMORY DELTA] {before['label']} → {after['label']}")
       print(f"  System: {delta_system:+.1f}% ({delta_process:+.2f}GB process)")
       print(f"  Top growing object types:")

       # Compare object summaries
       before_objs = {o[0]: o[1] for o in before['python_objects']}
       after_objs = {o[0]: o[1] for o in after['python_objects']}

       for obj_type in after_objs:
           before_count = before_objs.get(obj_type, 0)
           after_count = after_objs[obj_type]
           if after_count > before_count:
               print(f"    {obj_type}: +{after_count - before_count} bytes")
   ```

3. **Modify `08_transform_library_v2.py` to use profiling:**

   **Location:** After line 1190 (executor initialization)

   ```python
   # Add imports at top
   from diagnose_memory_flow import get_memory_snapshot, compare_snapshots

   # In cmd_apply(), after executor creation:
   snapshots = []

   # Before batch loop
   snapshots.append(get_memory_snapshot("before_batch_loop"))

   # In batch processing loop (every 10 batches)
   if batches_processed % 10 == 0:
       snapshot = get_memory_snapshot(f"batch_{batches_processed}")
       if snapshots:
           compare_snapshots(snapshots[-1], snapshot)
       snapshots.append(snapshot)

   # After all batches submitted
   snapshots.append(get_memory_snapshot("all_batches_submitted"))

   # After completion
   snapshots.append(get_memory_snapshot("pipeline_complete"))

   # Save snapshots to file
   import json
   with open(outdir / "memory_snapshots.json", 'w') as f:
       json.dump(snapshots, f, indent=2, default=str)
   ```

**Expected Output:** Detailed log showing which operations cause memory growth.

**Files to Modify:**
- `scripts/08_transform_library_v2.py` (add profiling hooks)
- `scripts/diagnose_memory_flow.py` (create new)

#### Subtask 1.2: Analyze Worker Process Memory

**Goal:** Track memory usage inside worker processes, not just main process.

**Problem:** Workers run in separate processes, their memory isn't visible in main process profiling.

**Implementation:**

1. **Modify `_process_batch_worker` function** (line 928-969):

   ```python
   def _process_batch_worker(batch_data: Tuple) -> Tuple[List[Dict], Dict]:
       """Worker with memory tracking."""
       import psutil
       import os

       rows = batch_data
       worker_pid = os.getpid()

       # Track worker memory before processing
       process = psutil.Process(worker_pid)
       mem_before = process.memory_info().rss / 1024**2  # MB

       all_products = []
       stats = Counter()

       for i, row in enumerate(rows):
           # ... existing processing code ...

           # Sample memory every 5000 rows
           if i % 5000 == 0 and i > 0:
               mem_current = process.memory_info().rss / 1024**2
               print(f"[Worker {worker_pid}] Processed {i}/{len(rows)}, "
                     f"Memory: {mem_current:.0f}MB (+{mem_current - mem_before:.0f}MB)",
                     flush=True)

       # Track memory after processing
       mem_after = process.memory_info().rss / 1024**2
       stats['_worker_memory_mb'] = int(mem_after)
       stats['_worker_memory_delta_mb'] = int(mem_after - mem_before)

       # CRITICAL: Measure size of products list before return
       import sys
       products_size_mb = sys.getsizeof(all_products) / 1024**2

       # CRITICAL: Estimate deep size of products
       sample_size = sum(sys.getsizeof(p) for p in all_products[:100]) / 100 if all_products else 0
       estimated_total_mb = (sample_size * len(all_products)) / 1024**2

       stats['_products_list_shallow_mb'] = int(products_size_mb)
       stats['_products_estimated_deep_mb'] = int(estimated_total_mb)
       stats['_products_count'] = len(all_products)

       return (all_products, dict(stats))
   ```

2. **Log worker stats in main process** (modify `_process_batch_result`, line 972-1022):

   ```python
   def _process_batch_result(products, stats, deduper, writer, counters, batch_idx):
       # Extract worker memory stats
       worker_mem = stats.get('_worker_memory_mb', 0)
       worker_delta = stats.get('_worker_memory_delta_mb', 0)
       products_shallow = stats.get('_products_list_shallow_mb', 0)
       products_deep = stats.get('_products_estimated_deep_mb', 0)
       products_count = stats.get('_products_count', 0)

       # Log detailed memory info
       logger.info(f"[Batch {batch_idx}] Worker stats: "
                   f"mem={worker_mem}MB (delta={worker_delta:+d}MB), "
                   f"products={products_count} ({products_deep:.0f}MB estimated)")

       # ... rest of existing code ...
   ```

**Expected Output:** Per-worker memory tracking showing which batches consume most memory.

#### Subtask 1.3: Trace Objects in _process_batch_result

**Goal:** Identify memory accumulation during result processing.

**Problem:** Products list gets copied multiple times during dedup and write. Each copy stays in memory temporarily.

**Implementation:**

**Modify `_process_batch_result`** (line 972-1022) to track object lifecycle:

```python
def _process_batch_result(products, stats, deduper, writer, counters, batch_idx):
    import sys
    import gc

    # CHECKPOINT 1: Products received from worker
    mem_1 = psutil.virtual_memory().percent
    products_size_1 = sum(sys.getsizeof(p) for p in products[:100]) * len(products) / 100 / 1024**2
    logger.info(f"[Batch {batch_idx} MEM-1] Received products: {len(products)}, "
                f"~{products_size_1:.1f}MB, sys_mem={mem_1:.1f}%")

    # Update stats
    counters['total_processed'] += len([p for p in products if p.get('source_smiles')])
    counters['total_products'] += len(products)
    counters['total_stats'].update(stats)

    # CHECKPOINT 2: After stats update
    mem_2 = psutil.virtual_memory().percent
    logger.info(f"[Batch {batch_idx} MEM-2] After stats: sys_mem={mem_2:.1f}% (delta={mem_2-mem_1:+.1f}%)")

    # In-batch deduplication
    unique_batch_products = []
    seen_in_batch = set()

    for prod in products:
        if not prod['xf_success']:
            continue
        canon_smi = prod.get('canonical_smiles')
        if not canon_smi or canon_smi in seen_in_batch:
            continue
        seen_in_batch.add(canon_smi)
        unique_batch_products.append(prod)

    # CHECKPOINT 3: After in-batch dedup
    mem_3 = psutil.virtual_memory().percent
    seen_size_mb = sys.getsizeof(seen_in_batch) / 1024**2
    unique_size_mb = sum(sys.getsizeof(p) for p in unique_batch_products[:100]) * len(unique_batch_products) / 100 / 1024**2
    logger.info(f"[Batch {batch_idx} MEM-3] After in-batch dedup: "
                f"{len(unique_batch_products)}/{len(products)} unique, "
                f"seen_set={seen_size_mb:.1f}MB, unique_list=~{unique_size_mb:.1f}MB, "
                f"sys_mem={mem_3:.1f}% (delta={mem_3-mem_2:+.1f}%)")

    # Cross-batch deduplication
    canon_smiles_list = [p['canonical_smiles'] for p in unique_batch_products]

    # CHECKPOINT 4: After creating SMILES list
    mem_4 = psutil.virtual_memory().percent
    smiles_list_mb = sys.getsizeof(canon_smiles_list) / 1024**2
    logger.info(f"[Batch {batch_idx} MEM-4] SMILES list created: "
                f"{len(canon_smiles_list)} items, {smiles_list_mb:.1f}MB, "
                f"sys_mem={mem_4:.1f}% (delta={mem_4-mem_3:+.1f}%)")

    new_keys = deduper.filter_new_keys(canon_smiles_list)

    # CHECKPOINT 5: After deduper
    mem_5 = psutil.virtual_memory().percent
    logger.info(f"[Batch {batch_idx} MEM-5] After deduper.filter_new_keys: "
                f"{len(new_keys)} new, sys_mem={mem_5:.1f}% (delta={mem_5-mem_4:+.1f}%)")

    new_keys_set = set(new_keys)
    final_products = [p for p in unique_batch_products if p['canonical_smiles'] in new_keys_set]

    # CHECKPOINT 6: After creating final list
    mem_6 = psutil.virtual_memory().percent
    final_size_mb = sum(sys.getsizeof(p) for p in final_products[:100]) * len(final_products) / 100 / 1024**2
    logger.info(f"[Batch {batch_idx} MEM-6] Final products: "
                f"{len(final_products)}, ~{final_size_mb:.1f}MB, "
                f"sys_mem={mem_6:.1f}% (delta={mem_6-mem_5:+.1f}%)")

    counters['total_unique'] += len(final_products)

    # Write to parquet
    writer.write_batch(final_products)

    # CHECKPOINT 7: After write
    mem_7 = psutil.virtual_memory().percent
    logger.info(f"[Batch {batch_idx} MEM-7] After write_batch: "
                f"sys_mem={mem_7:.1f}% (delta={mem_7-mem_6:+.1f}%)")

    deduper.commit()

    # CHECKPOINT 8: After deduper commit
    mem_8 = psutil.virtual_memory().percent
    logger.info(f"[Batch {batch_idx} MEM-8] After deduper.commit: "
                f"sys_mem={mem_8:.1f}% (delta={mem_8-mem_7:+.1f}%)")

    # Force cleanup
    del products
    del unique_batch_products
    del canon_smiles_list
    del new_keys
    del new_keys_set
    del final_products
    del seen_in_batch
    gc.collect()

    # CHECKPOINT 9: After explicit cleanup
    mem_9 = psutil.virtual_memory().percent
    logger.info(f"[Batch {batch_idx} MEM-9] After gc.collect: "
                f"sys_mem={mem_9:.1f}% (delta={mem_9-mem_8:+.1f}%)")
    logger.info(f"[Batch {batch_idx} TOTAL] Memory delta: {mem_9-mem_1:+.1f}%")

    return len(final_products)
```

**Expected Output:**
- Identify which step causes most memory retention
- See if `gc.collect()` actually frees memory
- Detect if references are held somewhere preventing cleanup

---

### Task 2: Investigate Final Stage Memory Spike

**Objective:** Understand why memory jumps from 72% to 94% when processing last 15 batches.

**Critical Observation:**
```
Batch 260/276: Memory 68.5%, products 89.3M
"All batches submitted. Processing remaining 15 in-flight batches..."
15 minutes later: Memory 94.2%
```

This suggests the last 15 batches are fundamentally different.

#### Subtask 2.1: Analyze Input Data Distribution

**Goal:** Determine if last batches have more complex molecules.

**Implementation:**

Create script `scripts/analyze_input_complexity.py`:

```python
"""Analyze polyphenol-2X input to find complexity distribution."""

import pyarrow.parquet as pq
from rdkit import Chem

input_file = "data/output/nplike_v2/polyphenol-2X/products.parquet"
table = pq.read_table(input_file)

batch_size = 50000
num_batches = (len(table) + batch_size - 1) // batch_size

print(f"Total rows: {len(table):,}")
print(f"Batches: {num_batches}")
print()

# Analyze each batch
complexity_by_batch = []

for batch_idx in range(num_batches):
    start = batch_idx * batch_size
    end = min(start + batch_size, len(table))
    batch = table.slice(start, end - start)

    smiles_list = batch.column('smiles').to_pylist()

    # Count phenolic OH sites per molecule
    phenolic_oh_counts = []
    for smi in smiles_list[:1000]:  # Sample 1000 per batch
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                # Pattern for phenolic OH: [c][OH]
                pattern = Chem.MolFromSmarts("[c][OX2H]")
                matches = mol.GetSubstructMatches(pattern)
                phenolic_oh_counts.append(len(matches))
        except:
            pass

    if phenolic_oh_counts:
        avg_sites = sum(phenolic_oh_counts) / len(phenolic_oh_counts)
        max_sites = max(phenolic_oh_counts)
        complexity_by_batch.append({
            'batch': batch_idx,
            'avg_phenolic_oh': avg_sites,
            'max_phenolic_oh': max_sites,
            'estimated_products': batch_size * avg_sites
        })

        if batch_idx % 10 == 0 or batch_idx >= num_batches - 20:
            print(f"Batch {batch_idx:3d}: avg_sites={avg_sites:.2f}, "
                  f"max_sites={max_sites}, est_products={batch_size * avg_sites:.0f}")

# Focus on last 20 batches
print("\n" + "="*80)
print("LAST 20 BATCHES ANALYSIS:")
print("="*80)
for batch in complexity_by_batch[-20:]:
    print(f"Batch {batch['batch']:3d}: {batch['avg_phenolic_oh']:.2f} avg sites, "
          f"{batch['estimated_products']:.0f} est products")

# Compare first 50 vs last 20
if len(complexity_by_batch) > 50:
    first_50_avg = sum(b['avg_phenolic_oh'] for b in complexity_by_batch[:50]) / 50
    last_20_avg = sum(b['avg_phenolic_oh'] for b in complexity_by_batch[-20:]) / 20

    print(f"\nFirst 50 batches avg: {first_50_avg:.2f} sites/molecule")
    print(f"Last 20 batches avg: {last_20_avg:.2f} sites/molecule")
    print(f"Complexity ratio: {last_20_avg / first_50_avg:.2f}x")
```

**Run this and analyze output to confirm if last batches are more complex.**

#### Subtask 2.2: Trace as_completed() Behavior

**Goal:** Understand memory behavior when collecting remaining futures.

**Problem Code** (line 1251-1262):
```python
# Process remaining futures
for future in as_completed(futures):
    completed_batch_idx = futures.pop(future)
    products, stats = future.result()  # ← ALL 15 results may already be in memory!
    _process_batch_result(...)
```

**Investigation Steps:**

1. **Add detailed logging before/during/after as_completed loop:**

```python
# After "All batches submitted" log
logger.info(f"Remaining in-flight: {len(futures)}")
logger.info(f"System memory before collecting: {psutil.virtual_memory().percent:.1f}%")

# Track memory per future
remaining_futures_mem = []

for idx, future in enumerate(as_completed(futures)):
    mem_before_result = psutil.virtual_memory().percent

    completed_batch_idx = futures.pop(future)
    products, stats = future.result()

    mem_after_result = psutil.virtual_memory().percent
    mem_after_process = None

    logger.info(f"[REMAINING {idx+1}/{len(futures)+idx+1}] "
                f"Batch {completed_batch_idx}: "
                f"mem_before_result={mem_before_result:.1f}%, "
                f"mem_after_result={mem_after_result:.1f}% "
                f"(delta={mem_after_result - mem_before_result:+.1f}%)")

    _process_batch_result(products, stats, deduper, writer, counters, completed_batch_idx)

    mem_after_process = psutil.virtual_memory().percent
    logger.info(f"[REMAINING {idx+1}/{len(futures)+idx+1}] "
                f"After process: {mem_after_process:.1f}% "
                f"(delta={mem_after_process - mem_after_result:+.1f}%)")

    remaining_futures_mem.append({
        'batch': completed_batch_idx,
        'before_result': mem_before_result,
        'after_result': mem_after_result,
        'after_process': mem_after_process,
        'delta_result': mem_after_result - mem_before_result,
        'delta_process': mem_after_process - mem_after_result
    })

# Log summary
logger.info("="*80)
logger.info("REMAINING FUTURES MEMORY SUMMARY:")
for item in remaining_futures_mem:
    logger.info(f"  Batch {item['batch']}: "
                f"result_delta={item['delta_result']:+.1f}%, "
                f"process_delta={item['delta_process']:+.1f}%")
logger.info("="*80)
```

2. **Check if futures.result() blocks or returns immediately:**

```python
# Before as_completed loop
import time

# Check if futures are already done
done_count = sum(1 for f in futures if f.done())
logger.info(f"Futures already done: {done_count}/{len(futures)}")

# Sample one future to see result size
if futures:
    sample_future = list(futures.keys())[0]
    start_time = time.time()
    sample_result = sample_future.result(timeout=1)
    elapsed = time.time() - start_time

    result_size_mb = sys.getsizeof(sample_result) / 1024**2
    logger.info(f"Sample future.result() took {elapsed:.3f}s, "
                f"returned {result_size_mb:.1f}MB (shallow size)")

    # Put it back (this is just sampling)
    # Note: Can't actually "put back", just log for info
```

**Expected Output:** Identify if memory spike is from:
- All 15 results materialized at once in ProcessPoolExecutor queue
- Individual results being very large
- Processing logic not releasing memory

---

### Task 3: Investigate Parquet Writer Memory

**Objective:** Understand why 9GB parquet file got corrupted and if writer is holding excessive memory.

#### Subtask 3.1: Analyze StreamingParquetWriter Buffer

**File:** `scripts/08_transform_library_v2.py` lines 327-549

**Current Implementation:**
```python
class StreamingParquetWriter:
    def __init__(self, ..., flush_interval=2000, ...):
        self.buffer = []
        self.flush_interval = flush_interval

    def write_batch(self, products):
        self.buffer.extend(products)  # ← Accumulates

        if len(self.buffer) >= self.flush_interval:
            self._flush()

    def _flush(self):
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)
        self.writer.write_table(table)
        self.buffer.clear()
```

**Questions to Answer:**

1. **Is buffer size really controlled?**
   - Multiple batches might call write_batch() before flush triggers
   - If 5 batches each write 2000 products, buffer could reach 10K before flush

2. **Does PyArrow Table.from_pylist() duplicate memory?**
   - Original products dicts in buffer
   - PyArrow table representation
   - Are both in memory simultaneously?

3. **Does writer.close() need extra memory?**
   - Parquet footer contains metadata
   - Large files might need significant memory for final flush

**Investigation Steps:**

1. **Add detailed buffer tracking:**

```python
def write_batch(self, products: List[Dict]):
    if not products:
        return

    # BEFORE extend
    buffer_size_before = len(self.buffer)
    mem_before = psutil.virtual_memory().percent

    self.buffer.extend(products)

    # AFTER extend
    buffer_size_after = len(self.buffer)
    mem_after_extend = psutil.virtual_memory().percent

    logger.info(f"[WRITER] Buffer: {buffer_size_before} → {buffer_size_after} "
                f"(+{len(products)}), mem: {mem_before:.1f}% → {mem_after_extend:.1f}% "
                f"(delta={mem_after_extend - mem_before:+.1f}%)")

    # ... rest of flush logic with detailed tracking ...
```

2. **Profile _flush() method:**

```python
def _flush(self):
    if not self.buffer:
        return

    import time
    flush_start = time.time()
    buffer_size = len(self.buffer)

    mem_start = psutil.virtual_memory().percent
    logger.info(f"[FLUSH START] Buffer: {buffer_size}, mem: {mem_start:.1f}%")

    # Step 1: Convert to table
    table = pa.Table.from_pylist(self.buffer, schema=self.schema)
    mem_after_table = psutil.virtual_memory().percent
    logger.info(f"[FLUSH] After Table.from_pylist: mem={mem_after_table:.1f}% "
                f"(delta={mem_after_table - mem_start:+.1f}%)")

    # Step 2: Write table
    if self.writer is None:
        self.writer = pq.ParquetWriter(self.output_path, self.schema)
    else:
        try:
            table = table.cast(self.schema)
        except Exception as e:
            logger.warning(f"Schema cast failed: {e}")
            # ... fallback logic ...

    self.writer.write_table(table)
    mem_after_write = psutil.virtual_memory().percent
    logger.info(f"[FLUSH] After write_table: mem={mem_after_write:.1f}% "
                f"(delta={mem_after_write - mem_after_table:+.1f}%)")

    # Step 3: Clear buffer
    self.buffer.clear()
    mem_after_clear = psutil.virtual_memory().percent
    logger.info(f"[FLUSH] After buffer.clear: mem={mem_after_clear:.1f}% "
                f"(delta={mem_after_clear - mem_after_write:+.1f}%)")

    # Step 4: GC
    import gc
    gc.collect()
    mem_after_gc = psutil.virtual_memory().percent
    logger.info(f"[FLUSH] After gc.collect: mem={mem_after_gc:.1f}% "
                f"(delta={mem_after_gc - mem_after_clear:+.1f}%)")

    flush_elapsed = time.time() - flush_start
    logger.info(f"[FLUSH COMPLETE] Took {flush_elapsed:.2f}s, "
                f"total mem delta: {mem_after_gc - mem_start:+.1f}%")
```

3. **Investigate writer.close() failure:**

The corrupted file suggests close() failed or was interrupted. Add:

```python
def close(self):
    logger.info(f"[WRITER CLOSE] Starting final flush, buffer: {len(self.buffer)}")
    mem_before_final_flush = psutil.virtual_memory().percent
    logger.info(f"[WRITER CLOSE] Memory before final flush: {mem_before_final_flush:.1f}%")

    if self.buffer:
        logger.info(f"[WRITER CLOSE] Flushing {len(self.buffer)} remaining products")
        try:
            self._flush()
            logger.info(f"[WRITER CLOSE] Final flush succeeded")
        except Exception as e:
            logger.error(f"[WRITER CLOSE] Final flush FAILED: {e}")
            raise

    if self.writer:
        logger.info(f"[WRITER CLOSE] Calling ParquetWriter.close()")
        mem_before_close = psutil.virtual_memory().percent
        logger.info(f"[WRITER CLOSE] Memory before writer.close(): {mem_before_close:.1f}%")

        try:
            self.writer.close()
            logger.info(f"[WRITER CLOSE] ParquetWriter.close() succeeded")
        except Exception as e:
            logger.error(f"[WRITER CLOSE] ParquetWriter.close() FAILED: {e}")
            raise

        mem_after_close = psutil.virtual_memory().percent
        logger.info(f"[WRITER CLOSE] Memory after close: {mem_after_close:.1f}% "
                    f"(delta={mem_after_close - mem_before_close:+.1f}%)")

    logger.info(f"[WRITER CLOSE] Complete. Total written: {self.total_written:,}")
```

---

### Task 4: Design Holistic Solution

**Objective:** Based on findings from Tasks 1-3, design comprehensive fix addressing ALL memory accumulation points.

**Potential Solutions to Consider:**

#### Option A: Reduce Batch Complexity Dynamically

If analysis shows last batches are more complex:

```python
# In batch iteration, estimate complexity and adjust
def estimate_batch_complexity(batch_df):
    """Quick check if batch has complex molecules."""
    sample = batch_df.head(100)
    # Check average molecular weight or SMILES length as proxy
    avg_len = sample['smiles'].str.len().mean()
    return avg_len

if estimate_batch_complexity(batch_df) > threshold:
    # Split batch into smaller sub-batches
    for sub_batch in split_batch(batch_df, smaller_size):
        process(sub_batch)
```

#### Option B: Streaming Product Processing

Instead of returning entire products list from worker:

```python
# Modify worker to yield chunks
def _process_batch_worker_streaming(batch_data):
    """Yield products in chunks instead of all at once."""
    chunk_size = 5000
    chunk = []

    for row in batch_data:
        products = transform(row)
        chunk.extend(products)

        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []

    if chunk:
        yield chunk
```

But this requires significant refactoring of ProcessPoolExecutor usage.

#### Option C: Reduce max_in_flight Further

Currently 16, reduce to 4-8:

```python
max_in_flight = max(4, workers // 4)  # 16 workers → 4 in-flight
```

Trade-off: Slower but safer.

#### Option D: Two-Pass Processing

1. **Pass 1:** Process all batches but only collect statistics, don't write products
2. **Pass 2:** Re-process using statistics to optimize memory (e.g., skip known duplicates)

Slow but guarantees completion.

#### Option E: External Dedup Service

Move deduplication out of process:
- Redis with SET operations
- Separate dedup process with IPC

Complex but isolates memory.

---

### Task 5: Implement and Validate

**After identifying root causes in Tasks 1-3 and designing solution in Task 4:**

1. **Implement chosen solution(s)**
2. **Test on polyphenol-2X with full instrumentation**
3. **Validate:**
   - Memory stays < 75% throughout
   - Output file is valid (not corrupted)
   - Dedup stats match expectations
4. **Run on terpenoid-2X** (even larger, 1.5GB input)
5. **Document final solution**

---

## Part C: File Reference

**Key Files:**

**Main Pipeline:**
- `scripts/08_transform_library_v2.py` - Main transform script (1500+ lines)
  - Line 253-320: SqliteDeduplicator class
  - Line 327-549: StreamingParquetWriter class
  - Line 928-969: _process_batch_worker function
  - Line 972-1022: _process_batch_result function
  - Line 1025-1340: cmd_apply main function
  - Line 1170-1249: Batch processing with max_in_flight control

**Deduplication:**
- `src/halogenator/bloom_dedup.py` - Bloom filter deduplicator (works well, 69MB)

**Test Data:**
- `data/output/nplike_v2/polyphenol-2X/products.parquet` - Input (504MB, 13.79M rows)
- `data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe/` - Failed output
  - `products.parquet` - 9GB, CORRUPTED
  - `dedup_bloom.pkl` - 69MB, valid

**Logs:**
- `polyphenol_2x_BLOOM_TEST.log` - Full task log (16h runtime)
- `bloom_memory_monitor.log` - Memory timeline (60min sample)

**Configuration:**
- `configs/transforms.yaml` - Transform rules

---

## Part D: Critical Unknowns

**Questions Next Session MUST Answer:**

1. **Where exactly does memory accumulate in final stage?**
   - Futures results queue?
   - Products copies during processing?
   - Parquet writer buffer?
   - Bloom filter (unlikely, it's only 69MB)?
   - Python object graph (circular references)?

2. **Why does gc.collect() not free memory?**
   - Are references held somewhere?
   - Are objects in extension modules (RDKit, PyArrow)?
   - Is memory fragmentation the issue?

3. **Why did last 15 batches cause 72% → 94% spike?**
   - Are they more complex molecules?
   - Does as_completed() materialize all results at once?
   - Is there a compounding effect from previous batches?

4. **Why did writer.close() produce corrupted file?**
   - Did it run out of memory during footer write?
   - Was process killed (OOM killer)?
   - Is there a PyArrow bug with large files under memory pressure?

---

## Part E: Success Criteria

**Solution is successful when:**

1. ✅ polyphenol-2X completes without corruption
2. ✅ Memory peak < 70% throughout entire run
3. ✅ Output file is valid and readable
4. ✅ Can process terpenoid-2X (1.5GB input) successfully
5. ✅ Solution is generalizable to all large transform jobs

**Acceptable Trade-offs:**
- Slower speed (e.g., 50% slower) if guaranteed to complete
- Higher disk usage if saves RAM
- Some duplicate products (1-2%) if saves memory

**Unacceptable:**
- Manual intervention required per job
- Hard-coded limits for specific datasets
- Solutions that only work for polyphenol-2X

---

## Part F: Recommended Approach

**Phase 1: Investigation (4-6 hours)**
1. Implement all instrumentation from Task 1
2. Run small test (first 50 batches of polyphenol-2X)
3. Analyze detailed memory logs
4. Identify TOP 3 memory accumulation sources

**Phase 2: Design (1-2 hours)**
1. Based on findings, choose targeted solutions
2. Design implementation plan
3. Identify risks and mitigations

**Phase 3: Implementation (2-3 hours)**
1. Implement chosen solutions
2. Add comprehensive logging
3. Test on full polyphenol-2X

**Phase 4: Validation (2-3 hours)**
1. Verify polyphenol-2X succeeds
2. Test on terpenoid-2X
3. Document solution and update pipeline

**Total: 9-14 hours**

---

## Part G: If You Get Stuck

**Fallback Strategies:**

1. **Split large datasets:**
   - Split polyphenol-2X into 4 chunks
   - Process separately
   - Merge with post-dedup

2. **Use external dedup:**
   - Set up Redis server
   - Use SET data structure for dedup
   - Much larger memory capacity

3. **Reduce parallelism drastically:**
   - workers=4, max_in_flight=2
   - Very slow but safe

4. **Stream to disk intermediate:**
   - Write products to temp files immediately
   - Dedup in second pass
   - Slower but constant memory

---

## Final Notes

**This is a CRITICAL blocker.** Without solving this, cannot complete:
- 18 polyphenol-2X jobs
- 20 terpenoid-2X jobs
- 24 alkaloid-2X jobs
- ~60 total large jobs

**Previous attempts were insufficient because they addressed symptoms, not root causes.**

**This session requires:**
- Deep investigation, not quick fixes
- Complete data flow analysis
- Holistic solution addressing ALL accumulation points

**Success will unlock the entire transform pipeline. Failure means manual splitting or external infrastructure.**

**Good luck. Be thorough. Don't assume previous analysis was complete.**

---

**END OF HANDOFF REPORT**
