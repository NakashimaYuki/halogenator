# Test Observation Guide - Buffer Fix Validation

## 什么是成功的测试？

### ✅ PASS Criteria

1. **Buffer Sizes**
   - Look for: `Flush triggered: buffer_full (X products)`
   - ✅ GOOD: All X values ≤ 2,000
   - ⚠️ ACCEPTABLE: Some X values between 2,000-3,000 (1.5x limit)
   - ❌ FAIL: Any X value > 5,000

2. **Memory Usage**
   - Look for: `[PROFILING] Pre-flush memory: X%`
   - ✅ GOOD: All X values < 65%
   - ⚠️ ACCEPTABLE: Some X values 65-70%
   - ❌ FAIL: Any X value > 75%

3. **Memory Deltas (After Flush)**
   - Look for: `[PROFILING] Post-flush memory: X%, delta: +Y%`
   - ✅ GOOD: Y is negative or 0 (memory freed)
   - ⚠️ ACCEPTABLE: Y is small positive (<1%)
   - ❌ FAIL: Y is large positive (>2%)

4. **No Critical Warnings**
   - Look for: `[!] Critical system memory:`
   - ✅ GOOD: Zero occurrences
   - ⚠️ ACCEPTABLE: 1-2 occurrences
   - ❌ FAIL: > 3 occurrences

5. **No Errors**
   - Look for: `ERROR` or `FAILED`
   - ✅ GOOD: Zero occurrences
   - ❌ FAIL: Any occurrences

6. **Output Validity**
   - ✅ GOOD: products.parquet file exists and readable
   - ✅ GOOD: SUMMARY.json shows reasonable stats
   - ❌ FAIL: File missing or corrupted

---

## 🔍 Example Output Analysis

### GOOD Example (Should See):
```
2025-12-27 10:15:32 | INFO | Flush triggered: buffer_full (2,000 products)
2025-12-27 10:15:32 | INFO | [PROFILING] Pre-flush memory: 58.3%
2025-12-27 10:15:33 | INFO | System memory after flush: 57.9%
2025-12-27 10:15:33 | INFO | [PROFILING] Post-flush memory: 57.9%, delta: -0.4%
```
✅ Buffer = 2,000 (perfect)
✅ Memory = 58.3% (low)
✅ Delta = -0.4% (freed memory!)

### ACCEPTABLE Example (OK):
```
2025-12-27 10:20:45 | INFO | Flush triggered: buffer_full (2,850 products)
2025-12-27 10:20:45 | INFO | [PROFILING] Pre-flush memory: 62.1%
2025-12-27 10:20:46 | INFO | System memory after flush: 62.5%
2025-12-27 10:20:46 | INFO | [PROFILING] Post-flush memory: 62.5%, delta: +0.4%
```
⚠️ Buffer = 2,850 (1.4x limit, acceptable)
⚠️ Memory = 62.1% (acceptable)
⚠️ Delta = +0.4% (small increase, OK)

### BAD Example (Would See in Old Code):
```
2025-12-26 17:23:16 | INFO | Flush triggered: buffer_full (807,111 products)
2025-12-26 17:23:16 | INFO | [PROFILING] Pre-flush memory: 90.3%
2025-12-26 17:23:16 | WARNING | [!] Critical system memory: 90.3%
2025-12-26 17:23:28 | INFO | System memory after flush: 91.4%
2025-12-26 17:23:28 | INFO | [PROFILING] Post-flush memory: 91.3%, delta: +1.0%
```
❌ Buffer = 807,111 (400x limit!)
❌ Memory = 90.3% (critical!)
❌ Delta = +1.0% (added memory instead of freeing!)
❌ Critical warning present

---

## 📊 Test Levels

### MICRO (Now)
- **Size:** 1 batch, 50K rows
- **Time:** ~2 minutes
- **Goal:** Verify basic functionality
- **Expected buffer flushes:** 5-10
- **Expected max buffer:** ≤ 2,000

### SMALL (Next)
- **Size:** 3 batches, 150K rows
- **Time:** ~5 minutes
- **Goal:** Verify sustained behavior
- **Expected buffer flushes:** 15-30
- **Expected max buffer:** ≤ 2,500

### MEDIUM (After SMALL)
- **Size:** 10 batches, 500K rows
- **Time:** ~15 minutes
- **Goal:** Verify no accumulation over time
- **Expected buffer flushes:** 50-100
- **Expected max buffer:** ≤ 2,500

### FULL (After ALL pass)
- **Size:** 276 batches, 13.79M rows
- **Time:** 16-20 hours
- **Goal:** Production validation
- **Expected buffer flushes:** 1,500-2,000
- **Expected max buffer:** ≤ 2,500

---

## ⚠️ Red Flags - Stop Immediately If You See:

1. **Buffer explosion still happening:**
   ```
   buffer_full (50,000 products)  ← 25x limit!
   ```

2. **Memory climbing steadily:**
   ```
   T+0min: 45%
   T+1min: 52%
   T+2min: 61%  ← Should be stable!
   ```

3. **Large positive deltas:**
   ```
   delta: +2.5%  ← Adding memory, not freeing!
   ```

4. **Frequent critical warnings:**
   ```
   [!] Critical system memory: 85.2%
   [!] Critical system memory: 87.8%
   [!] Critical system memory: 90.1%
   ```

5. **Any ERROR messages:**
   ```
   ERROR: Memory allocation failed
   ERROR: Parquet write failed
   ```

---

## 🎯 Decision Tree

```
Run MICRO test
    │
    ├─► All checks PASS → Run SMALL test
    │
    ├─► Some warnings but acceptable → Run SMALL test with caution
    │
    └─► Any FAIL criteria → STOP
         │
         └─► Investigate what failed
              │
              ├─► Buffer still exploding → Fix logic error
              ├─► Memory still high → Reduce flush_interval or max_in_flight
              └─► Errors → Debug specific error

Run SMALL test
    │
    ├─► All checks PASS → Run MEDIUM test
    │
    └─► Any issues → Stop and analyze

Run MEDIUM test
    │
    ├─► All checks PASS → Ready for FULL test
    │
    └─► Any issues → Stop and analyze

Run FULL test (overnight)
    │
    ├─► Success → Deploy to all transform jobs! 🎉
    │
    └─► Failure → Analyze logs, adjust parameters, retry
```

---

## 📝 Quick Checklist

Use this during/after each test:

```
MICRO TEST CHECKLIST:
[ ] Buffer sizes all ≤ 3,000
[ ] Memory stayed < 70%
[ ] Deltas mostly negative or small positive
[ ] No critical warnings (or < 3)
[ ] No ERROR messages
[ ] Output file valid
[ ] Throughput > 1000 mol/s

If ALL checked → Proceed to SMALL test
If ANY unchecked → Investigate before proceeding
```

---

## 🔬 Advanced Diagnostics

If you need to dig deeper:

1. **Count buffer sizes in log:**
   ```bash
   grep "buffer_full" test.log | grep -oP '\d+(?= products)' | sort -n
   ```

2. **Track memory progression:**
   ```bash
   grep "Pre-flush memory" test.log | grep -oP '\d+\.\d+(?=%)'
   ```

3. **Count warnings:**
   ```bash
   grep -c "Critical system memory" test.log
   ```

4. **Extract all deltas:**
   ```bash
   grep "delta:" test.log | grep -oP '[+-]\d+\.\d+(?=%)'
   ```

---

**Ready to test!** Run `python test_micro_manual.py` and watch carefully.
