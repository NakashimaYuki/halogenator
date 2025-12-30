# Transform Pipeline Memory Crisis - Session Handoff Report
**Date:** 2025-12-22
**Project:** Halogenator Transform Library Generation
**Status:** 🔴 CRITICAL - Memory optimization failed, requires deep investigation
**Next Session Goal:** Root cause analysis of rapid memory growth and design robust solution

---

## 📋 Executive Summary

### Current Situation
- ✅ **Completed:** 60 transform jobs (28 original + 30 aa_peptide + 2 failed polyphenol-2X)
- ❌ **Critical Issue:** Memory still reaches 99-100% despite system-wide monitoring
- ⏸️ **Status:** Pipeline suspended, corrupted outputs from OOM crashes
- 🎯 **Urgent Need:** Deep investigation into why memory grows too fast to flush

### Core Problem
**System memory monitoring is WORKING but INEFFECTIVE:**
- Memory jumps from 56% → 74% → 100% in minutes
- Flush triggers detected in logs
- BUT: Memory growth speed > Flush/GC release speed
- Result: OOM crash, corrupted parquet files

### Key Discovery
**Even with aggressive settings (workers=12, flush=1500, target=70%), memory still hit 100%.**
This suggests a fundamental problem beyond just parameter tuning.

---

## ✅ Part A: Completed Tasks (Summary)

### 1. Memory Monitoring Fix (Session 1)
**Problem Identified:**
```python
# OLD (WRONG) - monitored single process memory
process_mem = self._get_process_memory_mb()  # 800MB for one worker
return (process_mem / 32GB) * 100  # = 2.5% (never triggers!)

# NEW (CORRECT) - monitors system-wide memory
return psutil.virtual_memory().percent  # 65-95% (triggers correctly)
```

**Files Modified:**
- `scripts/08_transform_library_v2.py`: Lines 387-398, 420-442, 499-501
- `scripts/run_transform_pipeline.py`: Lines 378-395

**Changes:**
1. `_get_memory_usage_percent()` now returns system memory (not process)
2. `write_batch()` uses system memory for flush triggers
3. Flush conditions: buffer_full (primary), system_memory_pressure (secondary), critical (85%)
4. Default parameters: flush_interval=2000, target_memory=75%
5. Logs changed to "System memory after flush: XX%"

### 2. Test Validation
**Small dataset test (aa_peptide-1X):**
```
✓ Flush triggered: system_memory_pressure (65.7%)
✓ System memory after flush: 65.6%
✓ Monitoring WORKS correctly
```

### 3. Pipeline Execution Attempts

**Attempt 1:** Full pipeline (workers=16, flush=5000, target=45%)
- Result: Memory reached 95-99%, manual termination

**Attempt 2:** Conservative (workers=12, flush=1500, target=70%)
- Result: Memory reached 100%, OOM crash
- Evidence: memory_monitor.log shows 99-100% repeatedly
- Corrupted: polyphenol-2X jobs (parquet files damaged)

---

## 🔴 Part B: Core Problem - Deep Analysis Required

### Symptoms

**Memory Growth Pattern:**
```
Time    Python Mem    System Mem    Event
+0min   9.3GB (29%)   49%          Pipeline start
+1min   9.9GB (30%)   50%          Stable growth
+?min   18.4GB (57%)  56%          Acceleration begins
+?min   24.2GB (74%)  74%          Rapid jump
+?min   25.5GB (78%)  78%          Warning threshold
+?min   25.2GB (77%)  100%         OOM! System freeze
```

**Critical Observation:** Memory jumped from 56% to 100% in a very short time window (minutes, not hours).

### Evidence Files

**1. Memory Monitor Log**
```
Location: E:\Projects\halogenator\memory_monitor.log
Key findings:
- Peak system memory: 100.0% (multiple occurrences)
- Peak Python memory: 77.4% (25.2GB)
- Warnings triggered: 75-78% range
- Total iterations: 720+ (12 hours runtime)
```

**Key log entries:**
```
System memory used: 18397.4 MB (56.5%)
System memory used: 24193.6 MB (74.3%)  ← Jump!
System memory used: 24891.3 MB (76.4%)
System memory used: 25498.4 MB (78.3%)
...
System memory used: 32557.5 MB (100.0%)  ← OOM
```

**2. Failed Job Output**
```
Location: E:\Projects\halogenator\data\output\transforms\polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe\
Files:
- products.parquet: 2.8GB (CORRUPTED - magic bytes missing)
- dedup.db: 14GB (!)
- dedup.db-wal: 263MB

Error: ArrowInvalid: Parquet magic bytes not found in footer
```

**3. Pipeline Log**
```
Location: E:\Projects\halogenator\transform_pipeline_SAFE.log
Status: Job started but never completed (no completion log)
Last entry: "Starting: polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe"
```

### Root Cause Hypotheses

**Hypothesis 1: Buffer Size Still Too Large**
```
Current: 12 workers × 1500 products/buffer = 18,000 products in memory
If each product ~1-2MB (complex molecule): 18-36GB potential
Problem: All workers fill buffers before ANY flush occurs
```

**Hypothesis 2: Dedup Database Explosion**
```
Evidence: dedup.db reached 14GB for single job
Cause: SQLite WAL mode accumulates changes in memory
Even with checkpoint(), WAL grows faster than truncation
Potential: 14GB dedup.db loaded into memory for queries
```

**Hypothesis 3: Flush Latency**
```
Problem: Flush is triggered BUT takes time to complete:
1. Flush trigger detected (logged)
2. Write 1500 products to parquet (~5-10 seconds with schema cast)
3. GC runs (gc.collect() - can take seconds for 25GB)
4. During this time, OTHER workers continue adding data
5. Net result: Memory still grows during flush
```

**Hypothesis 4: Product Size Variability**
```
polyphenol-2X products are complex molecules:
- Input: 13.79M products
- Each product may have large descriptors, annotations
- 1500 products might be 1.5GB not 15MB
- flush_interval needs to be dynamic based on SIZE not COUNT
```

**Hypothesis 5: Worker Synchronization**
```
Problem: No coordination between workers
- Worker 1 triggers flush at 70% memory
- Workers 2-12 keep adding data
- By the time Worker 1 finishes flush, memory at 90%
- Worker 2 triggers flush, but others keep going
- Memory hits 100% before all flushes complete
```

---

## 📊 Part C: Data for Investigation

### System Configuration
```yaml
Hardware:
  Total RAM: 32GB
  Available for Python: ~28GB (accounting for OS)

Current Job:
  Input: 13.79M products (polyphenol-2X)
  Input size: ~500MB parquet
  Expected output: ~60M products (est.)

Pipeline Config (Conservative):
  workers: 12
  flush_interval: 1500 products
  batch_size: 50000
  target_memory: 70% system
  critical_memory: 85% system
```

### Memory Breakdown (at peak)
```
Total system: 32.56GB
Peak usage: 32.56GB (100%)

Python processes: 25.2GB (77%)
  - Main process: ~7-8GB
  - 12 Workers: ~200-800MB each = ~4-10GB
  - Dedup DB: 14GB

Other: 7.4GB (23%)
  - OS + other processes
```

### Flush Mechanism Current Implementation
```python
# Location: scripts/08_transform_library_v2.py:420-442

def write_batch(self, products: List[Dict]):
    self.buffer.extend(products)
    system_mem_percent = self._get_memory_usage_percent()
    buffer_size = len(self.buffer)

    # Condition 1: Buffer full
    if buffer_size >= self.flush_interval:  # 1500
        should_flush = True

    # Condition 2: System memory
    if system_mem_percent > self.target_memory_percent:  # 70%
        should_flush = True

    # Condition 3: Critical
    if system_mem_percent > 85.0:
        should_flush = True

    if should_flush:
        self._flush()  # Writes to disk + GC
```

### Dedup Checkpoint Implementation
```python
# Location: scripts/08_transform_library_v2.py:1121-1123

# Every 10 batches
if (batch_idx + 1) % 10 == 0:
    deduper.checkpoint()  # PRAGMA wal_checkpoint(TRUNCATE)
```

---

## 🎯 Part D: Tasks for Next Session (DETAILED)

### Task 1: Root Cause Deep Dive (CRITICAL)

**Objective:** Determine exactly WHY memory grows faster than flush can release it.

**Sub-tasks:**

#### 1.1: Measure Product Size Distribution
**Goal:** Understand actual memory footprint per product.

**Method:**
```python
# Create diagnostic script: scripts/measure_product_size.py
import sys
import pyarrow.parquet as pq
from pathlib import Path

def measure_products():
    """
    Read polyphenol-2X input and measure:
    1. Average product size in memory
    2. Size distribution (min, max, median, p95)
    3. Total memory for 1500 products
    """

    input_file = Path("data/output/nplike_v2/polyphenol-2X/products.parquet")
    table = pq.read_table(input_file)

    # Sample 10000 products
    sample = table.slice(0, 10000).to_pylist()

    # Measure size
    import pickle
    sizes = [len(pickle.dumps(p)) for p in sample]

    print(f"Average size: {sum(sizes)/len(sizes)/1024:.1f} KB")
    print(f"Min: {min(sizes)/1024:.1f} KB")
    print(f"Max: {max(sizes)/1024:.1f} KB")
    print(f"For 1500 products: {sum(sizes[:1500])/1024/1024:.1f} MB")
    print(f"For 18000 (12 workers): {sum(sizes[:1500])*12/1024/1024:.1f} MB")

if __name__ == "__main__":
    measure_products()
```

**Run:** `python scripts/measure_product_size.py`

**Expected Output:** Actual memory per product (current assumption: ~1-2MB, might be wrong)

#### 1.2: Profile Flush Timing
**Goal:** Measure how long flush operations take.

**Method:**
```python
# Modify scripts/08_transform_library_v2.py:_flush()
# Add timing instrumentation

def _flush(self):
    import time
    flush_start = time.time()

    # ... existing flush code ...

    flush_elapsed = time.time() - flush_start
    self.logger.info(f"[PROFILING] Flush took {flush_elapsed:.2f}s for {len(self.buffer)} products")
```

**Add to write_batch() before flush:**
```python
if should_flush:
    pre_flush_mem = self._get_memory_usage_percent()
    self.logger.info(f"[PROFILING] Pre-flush memory: {pre_flush_mem:.1f}%")
    self._flush()
    post_flush_mem = self._get_memory_usage_percent()
    self.logger.info(f"[PROFILING] Post-flush memory: {post_flush_mem:.1f}%, delta: {post_flush_mem - pre_flush_mem:.1f}%")
```

**Expected Output:**
- Flush duration (seconds)
- Memory delta (should be negative, but might be positive if workers add faster)

#### 1.3: Monitor Worker Memory Individually
**Goal:** See which workers consume most memory.

**Method:**
Create `scripts/monitor_workers.py`:
```python
import psutil
import time
from collections import defaultdict

def monitor_workers(interval=5, duration=300):
    """Monitor individual worker memory for 5 minutes."""

    worker_history = defaultdict(list)

    for i in range(duration // interval):
        timestamp = time.time()

        for proc in psutil.process_iter(['pid', 'name', 'memory_info', 'cmdline']):
            try:
                if 'python' in proc.info['name'].lower():
                    mem_mb = proc.info['memory_info'].rss / 1024 / 1024
                    pid = proc.info['pid']

                    worker_history[pid].append({
                        'time': timestamp,
                        'memory': mem_mb,
                        'cmdline': ' '.join(proc.info['cmdline'][:3])
                    })
            except:
                pass

        time.sleep(interval)

    # Analyze growth rates
    for pid, history in worker_history.items():
        if len(history) > 1:
            start_mem = history[0]['memory']
            end_mem = history[-1]['memory']
            growth_rate = (end_mem - start_mem) / (len(history) * interval / 60)  # MB/min
            print(f"PID {pid}: {start_mem:.1f}MB -> {end_mem:.1f}MB (growth: {growth_rate:.1f} MB/min)")

if __name__ == "__main__":
    monitor_workers()
```

**Run alongside pipeline:** `python scripts/monitor_workers.py &`

**Expected Output:** Identify if specific workers leak memory or if growth is uniform.

#### 1.4: Dedup Database Analysis
**Goal:** Understand why dedup.db grew to 14GB.

**Method:**
```python
# scripts/analyze_dedup.py
import sqlite3
from pathlib import Path

def analyze_dedup_db(db_path):
    """Analyze dedup database growth."""

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Get record count
    cursor.execute("SELECT COUNT(*) FROM seen_keys")
    count = cursor.fetchone()[0]

    # Get average key size
    cursor.execute("SELECT AVG(LENGTH(key)) FROM seen_keys LIMIT 10000")
    avg_size = cursor.fetchone()[0]

    # Database file size
    db_size = db_path.stat().st_size / 1024 / 1024 / 1024  # GB

    print(f"Records: {count:,}")
    print(f"Avg key size: {avg_size:.1f} bytes")
    print(f"DB size: {db_size:.2f} GB")
    print(f"Expected size: {count * avg_size / 1024 / 1024 / 1024:.2f} GB")
    print(f"Overhead: {db_size / (count * avg_size / 1024 / 1024 / 1024):.1f}x")

    # Check indexes
    cursor.execute("SELECT name, sql FROM sqlite_master WHERE type='index'")
    print(f"\nIndexes: {cursor.fetchall()}")

    conn.close()

if __name__ == "__main__":
    # Use a completed job's dedup.db
    db_path = Path("data/output/transforms/lipid-2X_FG_CARBOXYL__COOH__TO__COOMe/dedup.db")
    if db_path.exists():
        analyze_dedup_db(db_path)
```

**Expected Output:** Identify if dedup.db has excessive overhead or indexing issues.

---

### Task 2: Design Robust Solution

**Based on Task 1 findings, implement ONE of these solutions:**

#### Solution A: Adaptive Flush Based on Memory Size (Not Count)

**Concept:** Monitor actual buffer memory size, not product count.

**Implementation:**
```python
# Modify StreamingParquetWriter.__init__
def __init__(self, ..., max_buffer_mb=100):
    self.max_buffer_mb = max_buffer_mb
    self.buffer = []

# Modify write_batch()
def write_batch(self, products):
    self.buffer.extend(products)

    # Calculate buffer size in MB
    import sys
    buffer_size_mb = sum(sys.getsizeof(p) for p in self.buffer) / 1024 / 1024

    # Flush when buffer exceeds size limit
    if buffer_size_mb > self.max_buffer_mb:
        should_flush = True
        flush_reason = f"buffer_size ({buffer_size_mb:.1f}MB)"
```

**Advantage:** Adapts to product complexity automatically.

#### Solution B: Worker Pool with Memory Budget

**Concept:** Limit total memory budget across all workers.

**Implementation:**
```python
# Create shared memory tracker: scripts/shared_memory_tracker.py
import multiprocessing as mp

class SharedMemoryTracker:
    def __init__(self, max_total_mb=8000):
        self.manager = mp.Manager()
        self.worker_memory = self.manager.dict()
        self.lock = self.manager.Lock()
        self.max_total_mb = max_total_mb

    def can_allocate(self, worker_id, size_mb):
        with self.lock:
            current_total = sum(self.worker_memory.values())
            if current_total + size_mb <= self.max_total_mb:
                self.worker_memory[worker_id] = self.worker_memory.get(worker_id, 0) + size_mb
                return True
            return False

    def release(self, worker_id, size_mb):
        with self.lock:
            self.worker_memory[worker_id] = max(0, self.worker_memory.get(worker_id, 0) - size_mb)
```

**Usage in write_batch():**
```python
if not memory_tracker.can_allocate(worker_id, estimated_size):
    self._flush()  # Force flush to free budget
```

**Advantage:** Global memory coordination across workers.

#### Solution C: Streaming Write (No Buffer)

**Concept:** Write products immediately, no buffering.

**Implementation:**
```python
def write_batch(self, products):
    """Write immediately, no buffering."""
    if not products:
        return

    # Write directly without buffering
    table = pa.Table.from_pylist(products, schema=self.schema)

    if self.writer is None:
        self.writer = pq.ParquetWriter(self.output_path, self.schema)

    self.writer.write_table(table)
    self.total_written += len(products)

    # GC immediately
    import gc
    gc.collect()
```

**Advantage:** Minimizes memory usage, maximum safety.
**Disadvantage:** Slower (more I/O, less batching).

#### Solution D: Reduce Workers Drastically

**Concept:** Use 2-4 workers max, very small buffers.

**Parameters:**
```bash
workers=2
flush_interval=100  # Very small
target_memory=50%
```

**Advantage:** Simplest, most reliable.
**Disadvantage:** Very slow (but safe).

---

### Task 3: Implement and Validate

#### 3.1: Implement Chosen Solution
Based on Task 1 findings, select and implement best solution from Task 2.

#### 3.2: Test with Instrumentation
Run polyphenol-2X job with:
- Profiling enabled
- Memory monitoring (1-minute intervals)
- Worker monitoring
- Flush timing logs

#### 3.3: Success Criteria
- [ ] System memory stays below 75% throughout job
- [ ] No memory growth > 1GB/min sustained
- [ ] Job completes successfully
- [ ] Parquet file valid (readable)
- [ ] Dedup.db < 5GB

---

## 📁 Part E: Key File Locations

### Modified Code Files
```
E:\Projects\halogenator\scripts\08_transform_library_v2.py
  - Lines 387-398: _get_memory_usage_percent() - returns system memory
  - Lines 420-442: write_batch() - flush trigger logic
  - Lines 475-501: _flush() - GC and logging
  - Lines 1241-1256: CLI parameters (flush_interval, target_memory)
  - Lines 293-304: SqliteDeduplicator.checkpoint()

E:\Projects\halogenator\scripts\run_transform_pipeline.py
  - Lines 378-395: CLI parameters
  - Lines 134-164: TransformJob.run() - passes target_memory

E:\Projects\halogenator\scripts\monitor_memory.py
  - Memory monitoring utility (working correctly)
```

### Data Files
```
Input:
E:\Projects\halogenator\data\output\nplike_v2\polyphenol-2X\products.parquet
  - Size: ~500MB
  - Rows: 13,790,820 products

Failed Output:
E:\Projects\halogenator\data\output\transforms\polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe\
  - products.parquet: 2.8GB (CORRUPTED)
  - dedup.db: 14GB (!)
  - Status: Files locked, need cleanup after process termination
```

### Logs
```
E:\Projects\halogenator\transform_pipeline_SAFE.log
  - Last entry: job started but never completed

E:\Projects\halogenator\memory_monitor.log
  - 720+ iterations (12 hours)
  - Peak system memory: 100%
  - Peak Python memory: 77.4% (25.2GB)
```

---

## 🔬 Part F: Recommended Investigation Path

### Priority 1: Measure First (Don't Guess)
1. Run `measure_product_size.py` to get actual memory footprint
2. Add profiling to `_flush()` to measure flush duration
3. Monitor workers individually during test run

### Priority 2: Test Hypothesis
Based on measurements:
- If product size is large (>1MB): Implement Solution A (adaptive flush by size)
- If flush is slow (>10s): Implement Solution B (worker coordination)
- If dedup.db is the bottleneck: Consider alternative dedup strategy
- If unsure: Implement Solution D (reduce workers to 2-4)

### Priority 3: Incremental Testing
1. Test with workers=2, flush=100 first (safest)
2. If successful, gradually increase workers
3. Find optimal point where speed vs safety balance

---

## 🚨 Part G: Critical Warnings for Next Session

### 1. Don't Trust Previous Assumptions
```
WRONG: "flush_interval=1500 is conservative"
RIGHT: Need to measure actual product size first

WRONG: "System memory monitoring will prevent OOM"
RIGHT: Monitoring works but memory grows too fast to react
```

### 2. Dedup Database is a Major Issue
```
14GB for single job is ABNORMAL
Possible causes:
- WAL file not truncating properly
- In-memory cache growing unbounded
- Index overhead

Consider alternatives:
- Bloom filter (probabilistic, much smaller)
- Periodic dedup.db reset (accept some duplicates)
- External dedup service
```

### 3. Worker Parallelism May Be Fundamentally Limited
```
With 32GB RAM:
- OS + overhead: ~4GB
- Available: ~28GB
- Per worker budget: 28GB / N workers

If products are 1-2MB each and workers need 8GB:
Maximum workers = 3-4 (not 12 or 16!)
```

### 4. File Cleanup Required
```bash
# Before starting next session, ensure:
1. No Python processes running
2. Clean corrupted output:
   rm -rf data/output/transforms/polyphenol-2X_*
3. Verify: ls data/output/transforms | wc -l should be 58
```

---

## 📊 Part H: Quick Reference - Current State

### Completed Jobs: 58/118
```
✓ polyphenol-1X: 16 jobs (9.19M products)
✓ lipid-1X: 6 jobs (24K products)
✓ lipid-2X: 6 jobs (292K products)
✓ aa_peptide-1X: 15 jobs (97K products)
✓ aa_peptide-2X: 15 jobs (3.05M products)

✗ polyphenol-2X: 0 jobs (2 failed attempts)
⏳ Remaining: 60 jobs
```

### Failed Jobs (Corrupted)
```
1. polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe (OOM crash, parquet corrupted)
2. polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2 (previous session failure)
```

### System Status
```
Current memory: 42% (healthy, no jobs running)
Python processes: 3 (idle)
Corrupted files: Need cleanup (files locked)
```

---

## 🎯 Part I: Immediate Next Steps

### Step 1: Environment Cleanup
```bash
# Kill any remaining Python processes
tasklist | grep python.exe
# Manually kill if needed

# Clean corrupted outputs
rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe
rm -rf data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__NH2

# Verify clean state
ls data/output/transforms | wc -l  # Should be 58
```

### Step 2: Run Diagnostic Scripts
```bash
# Measure product size
python scripts/measure_product_size.py > product_size_analysis.txt

# Check current code
grep -n "def _flush" scripts/08_transform_library_v2.py
grep -n "flush_interval" scripts/08_transform_library_v2.py
```

### Step 3: Design Solution
Based on diagnostic findings, choose and implement solution from Part D, Task 2.

### Step 4: Conservative Test
```bash
# Test with ultra-safe parameters
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/polyphenol-2X/products.parquet \
  --outdir data/output/transforms/TEST_polyphenol-2X_SAFE \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --workers 2 \
  --flush-interval 100 \
  --batch-size 20000 \
  --target-memory 50.0 \
  > test_ultra_safe.log 2>&1

# Monitor in parallel
python scripts/monitor_memory.py --interval 30 --threshold 60 &
```

### Step 5: If Test Succeeds
Gradually increase workers and flush_interval until finding optimal balance.

### Step 6: If Test Fails
Implement Solution C (streaming write with no buffer) as last resort.

---

## 📞 Part J: Contact Points for Analysis

### When Analyzing, Focus On:
1. **Memory growth rate** - MB/minute during job execution
2. **Flush effectiveness** - Memory before vs after flush
3. **Worker coordination** - Do all workers flush together or independently?
4. **Dedup overhead** - Why 14GB for single job?
5. **Product complexity** - Actual size in memory vs on disk

### Red Flags to Watch:
- Memory growth > 2GB/min sustained
- Flush duration > 30 seconds
- Post-flush memory higher than pre-flush
- Dedup.db > 10GB for single job
- Workers growing memory independently (no coordination)

---

## 🏁 Part K: Success Definition

**The solution is successful when:**
1. polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe job completes
2. System memory peak < 80%
3. Output parquet file is valid (readable)
4. Dedup.db < 5GB
5. Process is repeatable for all 4 polyphenol-2X phenolic OH jobs

**Then:** Apply same solution to remaining 60 jobs and complete full pipeline.

---

## 📝 Part L: Session Completion Checklist

Before ending next session, ensure:
- [ ] Root cause identified with evidence
- [ ] Solution implemented and tested
- [ ] At least 1 polyphenol-2X job completes successfully
- [ ] Documentation updated with findings
- [ ] Parameters finalized for production run
- [ ] Handoff report for final production run (if needed)

---

**End of Handoff Report**

**Priority:** 🔴 CRITICAL
**Complexity:** HIGH - Requires deep investigation, not just parameter tuning
**Estimated Time:** 4-8 hours for investigation + implementation + testing
**Risk:** Memory issues may require fundamental architecture changes

**Key Insight for Next Claude:**
The problem is NOT that monitoring doesn't work - it works perfectly.
The problem is that memory grows FASTER than flush can release it.
You need to find WHY (product size? dedup? worker parallelism? flush latency?)
and design a solution that addresses the root cause, not just symptoms.

**Good luck!** 🚀
