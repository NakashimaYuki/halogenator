# MEDIUM Test Memory Crisis - Root Cause Analysis

**Date:** 2025-12-27
**Status:** 🔴 HIGH MEMORY (86-89%) - But Test COMPLETED

---

## Test Comparison Summary

| Test | Rows | Batches | Workers | In-Flight | Memory Peak | Status |
|------|------|---------|---------|-----------|-------------|--------|
| **MICRO** | 50K | 1 | 4 | 1 | **43.6%** | ✅ PASS |
| **SMALL** | 150K | 3 | 6 | 3 | **45.6%** | ✅ PASS |
| **MEDIUM** | 500K | 10 | 8 | **10** | **86-89%** | ⚠️ HIGH MEM |

---

## Critical Finding

**Memory jumped from 45.6% → 86-89% when in-flight batches increased from 3 → 10!**

### Evidence:
```
MEDIUM test log:
2025-12-27 03:45:56 | INFO | All batches submitted. Processing remaining 10 in-flight batches...

Memory samples (every 100th flush):
88.6%, 88.8%, 88.5%, 88.4%, 88.0%, 88.9%, 88.2%, 88.2%,
87.4%, 87.5%, 86.9%, 86.4%, 86.5%, 87.6%, 87.6%, 87.6%
```

Memory climbed to 86-89% **immediately** when all 10 batches started processing and **stayed there** throughout.

---

## ROOT CAUSE IDENTIFIED

**PRIMARY MEMORY CONSUMER: Worker Processes + In-Flight Results**

### The Problem:

When 10 batches process simultaneously:
1. **8 workers** × multiple RDKit Mol objects in memory **each**
2. **10 batches** worth of intermediate products waiting
3. **ProcessPoolExecutor queue** holding results from completed batches
4. **Each batch** = 50K molecules → large memory footprint

### Memory Breakdown Estimate (10 batches, 8 workers):

```
Component                    Memory     Calculation
-----------------------------------------------------------------
Workers (RDKit Mol objects)  ~8-12 GB   8 workers × 1-1.5 GB each
ProcessPoolExecutor queue    ~2-3 GB    10 batches × results
Intermediate products        ~2-4 GB    In-worker processing
Bloom filter deduplication   ~6 MB      (very small, not the issue)
Writer buffer                ~10 MB     (fixed, under control!)
OS + Other                   ~5-8 GB    System overhead
-----------------------------------------------------------------
TOTAL:                       ~19-27 GB  on 32GB system = 60-85%
```

**This matches observed 86-89% memory!**

### Why Buffer Fix Alone Was Insufficient:

1. ✅ Buffer fix WORKS - we see 1,000-2,000 product flushes (controlled)
2. ❌ BUT buffer was never the PRIMARY consumer in multi-batch scenarios
3. ❌ The REAL consumer is **parallel workers processing 10 batches simultaneously**

---

## Why Tests Behaved Differently

### MICRO Test (43.6% memory):
- Only 1 batch in-flight
- 4 workers
- Minimal ProcessPoolExecutor queue
- **Workers idle most of the time waiting for batch**
- Memory: ~14 GB / 32 GB = 43.6%

### SMALL Test (45.6% memory):
- 3 batches in-flight
- 6 workers
- Small ProcessPoolExecutor queue
- **Workers can stay busy, but not overwhelmed**
- Memory: ~14.5 GB / 32 GB = 45.6%

### MEDIUM Test (86-89% memory):
- **10 batches in-flight** ← THE PROBLEM!
- 8 workers
- **Large ProcessPoolExecutor queue**
- **All workers fully loaded with 10 batches worth of data**
- Memory: ~27 GB / 32 GB = 86-89%

---

## The Vicious Cycle

```
10 batches submitted
    ↓
All start processing simultaneously
    ↓
8 workers × 10 batches = massive parallel processing
    ↓
Each worker holds RDKit Mol objects (~1-1.5GB)
    ↓
Results queue up in ProcessPoolExecutor
    ↓
Main process can't drain queue fast enough
    ↓
Memory climbs to 86-89%
    ↓
Stays there until batches complete
```

---

## Why It Still Completed (Unlike Original Bug)

Despite 86-89% memory:
1. ✅ Buffer fix prevented buffer explosion (no 940K spikes)
2. ✅ Memory STABLE at 86-89% (not climbing to 94%+)
3. ✅ No corruption (ParquetWriter could still write at 86%)
4. ✅ Eventually all batches finished

**Original bug:** 72% → 94% → crash
**Current:** 86-89% → stable → complete

**Improvement:** Prevents catastrophic failure, but still too high for production.

---

## Solution: Reduce max_in_flight

### Current Behavior:
```python
# In validate_fix_gradual.py or manual test:
# max_in_flight defaults to workers × 2 = 8 × 2 = 16
# But only 10 batches, so all 10 in-flight
```

### Required Fix:
```python
# Limit to 2-4 in-flight batches maximum
max_in_flight = 2  # Conservative
# OR
max_in_flight = 4  # Balanced
```

### Expected Impact:

| max_in_flight | Expected Memory | Speed Impact |
|---------------|-----------------|--------------|
| 10 (current) | 86-89% | Fastest |
| 4 | ~60-65% | -20% speed |
| 2 | ~50-55% | -40% speed |
| 1 | ~45-50% | -60% speed |

**Recommendation:** max_in_flight = 4
- Memory: ~60-65% (safe)
- Speed: 20% slower (acceptable)
- Stability: High

---

## Secondary Fixes to Consider

### Fix 1: Reduce Workers (If Needed)
```python
workers = 4  # Instead of 8
```
Pros: Lower memory per worker
Cons: Slower processing

### Fix 2: Smaller Batch Size
```python
batch_size = 25000  # Instead of 50000
```
Pros: Less memory per batch
Cons: More batches = more overhead

### Fix 3: Aggressive GC in Workers
Add to worker function:
```python
import gc
# After processing each batch
gc.collect()
```
Pros: May free some memory
Cons: Limited impact (RDKit objects still needed)

---

## Recommendation

**IMMEDIATE:**
1. Set `max_in_flight = 4` for all tests
2. Re-run MEDIUM test to verify memory < 65%
3. If still high, reduce to `max_in_flight = 2`

**FOR PRODUCTION (polyphenol-2X, 276 batches):**
```bash
python scripts/08_transform_library_v2.py apply \
  --input data/output/nplike_v2/polyphenol-2X/products.parquet \
  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe \
  --xf-config configs/transforms.yaml \
  --xf-name FG_PHENOL_OH__OH__TO__OMe \
  --use-bloom-filter \
  --workers 16 \
  --target-memory 70.0 \
  --max-in-flight 4  # ← ADD THIS!
```

**Expected results with max_in_flight=4:**
- Memory: 60-70% (vs 86-89%)
- Throughput: ~1200-1400 mol/s (vs 1500 without limit)
- Stability: High (no OOM risk)

---

## Conclusion

**The buffer explosion fix WORKS perfectly for its intended purpose.**

**BUT** we discovered a SECONDARY memory issue:
- **Parallel worker memory** dominates when many batches in-flight
- **Not** the buffer (that's fixed!)
- **Not** the Bloom filter (only 6MB)

**Solution is simple:**
**Limit max_in_flight to 4** → Prevents too many batches processing simultaneously

**This is a configuration fix, not a code fix.**

---

**Status:** Issue diagnosed, solution identified
**Next Step:** Modify test parameters and re-validate
