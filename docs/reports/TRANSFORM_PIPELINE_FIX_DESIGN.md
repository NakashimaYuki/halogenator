# Transform Pipeline Memory Fix - Architecture Design

## Root Cause Analysis (FINAL)

### Memory Accumulation Points Identified

1. **PRIMARY CAUSE: Submit-All-Then-Collect Pattern**
   ```python
   # Current (BROKEN):
   futures = []
   for batch in all_batches:  # 276 batches
       futures.append(executor.submit(process, batch))  # Submit ALL

   # Results pile up in ProcessPoolExecutor result queue!

   for future in futures:  # THEN collect
       result = future.result()
   ```

   **Problem:** Early-completing batches wait in memory while main process is still submitting later batches.

2. **SECONDARY CAUSE: High Product Multiplier**
   - Input: 1 polyphenol parent with 4 phenolic OH groups
   - Output: 4 products (one per site)
   - 50,000 input × 4 products = 200,000 products per batch
   - Some polyphenols have 6-8 phenolic OHs → up to 400K products/batch

3. **TERTIARY CAUSE: Product Size in Memory**
   - Product dict: ~1.2 KB deep size
   - 500K products × 1.2 KB = 600 MB per batch
   - 50 batches waiting = 30 GB

### Memory Breakdown (from crashed run)

```
Total system memory: 30.2GB / 31.8GB (94.8%)

Main process (PID 23204): 18.5 GB
  - ProcessPoolExecutor result queue: ~15 GB (estimated)
  - Bloom filter: 0.1 GB
  - Parquet writer buffer: 0.5 GB
  - Other: 2.9 GB

Worker processes (8): 6-800 MB each = ~5 GB
  - Processing current batches
  - Temporary Mol objects

Other processes: ~7 GB
```

---

## Solution Design

### Option A: Streaming Collection with In-Flight Limit (RECOMMENDED)

**Principle:** Never have more than N batches in-flight at once.

```python
from concurrent.futures import as_completed, wait, FIRST_COMPLETED

MAX_IN_FLIGHT = workers * 2  # e.g., 8 workers × 2 = 16 batches max

batch_iter = parquet_file.iter_batches(...)
futures = {}  # {future: batch_idx}

for batch_idx, batch in enumerate(batch_iter):
    # Submit new batch
    batch_rows = batch.to_pandas().to_dict('records')
    future = executor.submit(_process_batch_worker, batch_rows)
    futures[future] = batch_idx

    # If at limit, wait for one to complete before submitting next
    while len(futures) >= MAX_IN_FLIGHT:
        done, pending = wait(futures, return_when=FIRST_COMPLETED)
        for future in done:
            batch_idx = futures.pop(future)
            products, stats = future.result()

            # Process result immediately
            process_and_write(products, stats)

# Process remaining futures
for future in as_completed(futures):
    batch_idx = futures.pop(future)
    products, stats = future.result()
    process_and_write(products, stats)
```

**Benefits:**
- Max memory = MAX_IN_FLIGHT × products_per_batch × product_size
- 16 batches × 500K products × 1.2 KB = ~9.6 GB (safe!)
- Results consumed immediately, no queue buildup
- Maintains parallelism (16 in-flight vs 1 sequential)

**Trade-offs:**
- Slightly more complex code
- Order of processing may differ (doesn't matter for us)

### Option B: Chunked Submission

**Principle:** Submit in chunks, collect each chunk before next.

```python
CHUNK_SIZE = workers * 2

all_batches = list(parquet_file.iter_batches(...))

for chunk_start in range(0, len(all_batches), CHUNK_SIZE):
    chunk = all_batches[chunk_start:chunk_start + CHUNK_SIZE]

    # Submit chunk
    futures = []
    for batch in chunk:
        future = executor.submit(process, batch)
        futures.append(future)

    # Collect chunk results
    for future in as_completed(futures):
        result = future.result()
        process_and_write(result)
```

**Benefits:**
- Simple to understand
- Predictable memory usage

**Trade-offs:**
- Requires materializing batch list (minor overhead)
- Less smooth throughput (gaps between chunks)

### Option C: Reduce Batch Size (FALLBACK)

**If Options A/B still have issues:**

```python
batch_size = 10000  # Down from 50000
# Reduces products per batch: 10K × 4 = 40K (vs 200K)
```

**Benefits:**
- Dead simple
- Reduces peak per-batch memory

**Trade-offs:**
- More batches (1380 instead of 276)
- More overhead (Python/parquet I/O)
- Slower overall

---

## Implementation Plan

### Phase 1: Implement Option A (Streaming with Limit)

**Files to modify:**
- `scripts/08_transform_library_v2.py`
  - Function: `cmd_apply()`
  - Lines: 1125-1200 (batch submission and collection loop)

**Changes:**
1. Import `wait, FIRST_COMPLETED` from `concurrent.futures`
2. Replace submit-all-then-collect with streaming pattern
3. Add MAX_IN_FLIGHT parameter (default: workers × 2)
4. Add CLI argument `--max-in-flight` for tuning

**Pseudocode:**
```python
def cmd_apply(args):
    # ... existing setup ...

    MAX_IN_FLIGHT = getattr(args, 'max_in_flight', workers * 2)
    logger.info(f"Max in-flight batches: {MAX_IN_FLIGHT}")

    futures = {}
    batch_iter = parquet_file.iter_batches(...)

    for batch_idx, batch in enumerate(batch_iter):
        # Convert batch
        batch_df = batch.to_pandas()
        batch_df['_source_subset'] = source_subset
        batch_rows = batch_df.to_dict('records')

        # Submit
        future = executor.submit(_process_batch_worker, batch_rows)
        futures[future] = batch_idx

        # If at limit, wait for completions
        while len(futures) >= MAX_IN_FLIGHT:
            done, pending = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                batch_idx = futures.pop(future)

                # Process result
                products, stats = future.result()
                process_products(products, stats, deduper, writer, ...)

    # Process remaining
    for future in as_completed(futures):
        batch_idx = futures.pop(future)
        products, stats = future.result()
        process_products(products, stats, deduper, writer, ...)
```

### Phase 2: Extract Result Processing

**Reason:** DRY - avoid duplicating the result processing logic

```python
def process_batch_result(products, stats, deduper, writer, counters):
    """
    Process batch result: dedup and write.

    Extracted to avoid duplication in streaming loop.
    """
    # Update stats
    counters['total_processed'] += len([p for p in products if p.get('source_smiles')])
    counters['total_products'] += len(products)
    counters['total_stats'].update(stats)

    # In-batch dedup
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

    # Cross-batch dedup
    canon_smiles_list = [p['canonical_smiles'] for p in unique_batch_products]
    new_keys = deduper.filter_new_keys(canon_smiles_list)
    new_keys_set = set(new_keys)
    final_products = [p for p in unique_batch_products if p['canonical_smiles'] in new_keys_set]

    counters['total_unique'] += len(final_products)

    # Write
    writer.write_batch(final_products)
    deduper.commit()
```

### Phase 3: Add CLI Argument

```python
parser_apply.add_argument(
    '--max-in-flight',
    type=int,
    default=None,  # Auto: workers * 2
    help='Maximum in-flight batches (default: workers × 2). '
         'Lower values reduce memory usage. '
         'Recommended: 2-4 for 32GB RAM systems processing large datasets.'
)
```

### Phase 4: Add Monitoring

**Log progress with in-flight count:**
```python
logger.info(f"[Batch {batch_idx+1}/{total_batches}] "
            f"In-flight: {len(futures)}/{MAX_IN_FLIGHT} | "
            f"Processed: {total_processed:,} | "
            f"Products: {total_products:,}")
```

---

## Testing Strategy

### Test 1: Small Dataset (Validation)

**Input:** aa_peptide-1X (30K products)
**Config:**
- workers=2
- max_in_flight=4
- batch_size=5000

**Expected:**
- Memory stable < 40%
- No queue buildup
- Completion successful

### Test 2: Medium Dataset (Stress)

**Input:** polyphenol-1X (9.19M products)
**Config:**
- workers=4
- max_in_flight=8
- batch_size=50000

**Expected:**
- Memory < 60%
- Smooth progress (no long pauses)
- Completion in reasonable time

### Test 3: Full polyphenol-2X (Production)

**Input:** polyphenol-2X (13.79M products)
**Config:**
- workers=8
- max_in_flight=16
- batch_size=50000
- use_bloom_filter=True

**Expected:**
- Memory peak < 75%
- No OOM
- Completion successful
- Valid parquet output

---

## Fallback Strategy

**If Option A still shows memory issues:**

1. **Reduce max_in_flight:**
   - Try max_in_flight=4 (workers/2 instead of workers*2)
   - Trade: slower, but safer

2. **Reduce batch_size:**
   - Try batch_size=20000
   - Reduces peak per-batch memory

3. **Split input file:**
   - Split polyphenol-2X into 5 chunks (2.7M each)
   - Process separately, merge outputs
   - Guaranteed to work but requires manual splitting

---

## Success Metrics

1. ✅ Memory stays < 75% during entire run
2. ✅ No "CRITICAL_SYSTEM_MEMORY" warnings
3. ✅ In-flight batch count never exceeds MAX_IN_FLIGHT
4. ✅ Job completes successfully
5. ✅ Output parquet file is valid
6. ✅ Total execution time reasonable (< 3 hours)

---

## Implementation Checklist

- [ ] Modify cmd_apply() to use streaming pattern
- [ ] Extract process_batch_result() helper function
- [ ] Add --max-in-flight CLI argument
- [ ] Add in-flight count to progress logging
- [ ] Test on aa_peptide-1X (small)
- [ ] Test on polyphenol-1X (medium)
- [ ] Test on polyphenol-2X (full)
- [ ] Document findings
- [ ] Update USAGE_GUIDE.md

---

**END OF DESIGN DOCUMENT**
