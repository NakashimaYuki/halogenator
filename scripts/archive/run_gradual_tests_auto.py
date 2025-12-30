#!/usr/bin/env python3
"""
Automated gradual validation - no interactive prompts.
Runs SMALL and MEDIUM tests sequentially.
"""

import subprocess
import sys
import time
from pathlib import Path
import pyarrow.parquet as pq
import json
import re

def run_test(level, max_rows, workers, timeout_sec, target_memory):
    """Run a single test level."""
    print("\n" + "="*80)
    print(f"TEST LEVEL: {level}")
    print("="*80)

    batch_size = 50000
    expected_batches = (max_rows + batch_size - 1) // batch_size

    print(f"Parameters:")
    print(f"  Input rows: {max_rows:,}")
    print(f"  Expected batches: {expected_batches}")
    print(f"  Workers: {workers}")
    print(f"  Timeout: {timeout_sec//60} minutes")
    print(f"  Target memory: {target_memory}%")
    print()

    # Paths
    input_file = r"E:\Projects\halogenator\data\output\nplike_v2\polyphenol-2X\products.parquet"
    transform_config = r"E:\Projects\halogenator\configs\transforms.yaml"
    transform_name = "FG_PHENOL_OH__OH__TO__OMe"
    script = r"E:\Projects\halogenator\scripts\08_transform_library_v2.py"

    outdir = Path(r"E:\Projects\halogenator\data\output\transforms") / f"TEST_{level}_polyphenol"
    if outdir.exists():
        import shutil
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Create subset
    print(f"Creating subset ({max_rows:,} rows)...")
    table = pq.read_table(input_file)
    subset = table.slice(0, max_rows)
    temp_input = outdir / "input_subset.parquet"
    pq.write_table(subset, temp_input)
    print(f"Subset created: {len(subset):,} rows")
    print()

    # Build command
    cmd = [
        sys.executable,
        script,
        "apply",
        "--input", str(temp_input),
        "--outdir", str(outdir),
        "--xf-config", transform_config,
        "--xf-name", transform_name,
        "--batch-size", "50000",
        "--workers", str(workers),
        "--flush-interval", "2000",
        "--target-memory", str(target_memory),
        "--use-bloom-filter",
        "--bloom-expected-items", str(max_rows * 10),
    ]

    print(f"Starting test at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("-"*80)

    start_time = time.time()
    log_file = outdir / "test.log"

    try:
        with open(log_file, 'w') as log:
            process = subprocess.Popen(
                cmd,
                cwd=r"E:\Projects\halogenator",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            for line in process.stdout:
                print(line, end='')
                log.write(line)
                log.flush()

            process.wait(timeout=timeout_sec)
            elapsed = time.time() - start_time

            print("-"*80)
            print(f"Completed in {elapsed:.1f}s ({elapsed/60:.1f} min)")

            if process.returncode != 0:
                print(f"[FAILED] Exit code: {process.returncode}")
                return False

    except subprocess.TimeoutExpired:
        print(f"\n[FAILED] Timeout after {timeout_sec//60} minutes")
        process.kill()
        return False

    # Validate
    print("\n" + "="*80)
    print("VALIDATION")
    print("="*80)

    products_file = outdir / "products.parquet"
    if not products_file.exists():
        print("[FAILED] Output file not found")
        return False

    try:
        output_table = pq.read_table(products_file)
        file_size_mb = products_file.stat().st_size / 1024**2
        print(f"[OK] Output file: {len(output_table):,} rows, {file_size_mb:.1f} MB")
    except Exception as e:
        print(f"[FAILED] Output corrupted: {e}")
        return False

    # Check summary
    summary_file = outdir / "SUMMARY.json"
    if summary_file.exists():
        with open(summary_file) as f:
            summary = json.load(f)

        processed = summary.get('total_processed', 0)
        products = summary.get('total_products', 0)
        unique = summary.get('unique_products', 0)
        throughput = summary.get('throughput_mol_per_sec', 0)

        print(f"\n[OK] Summary:")
        print(f"  Processed: {processed:,}")
        print(f"  Products: {products:,}")
        print(f"  Unique: {unique:,}")
        print(f"  Throughput: {throughput:.1f} mol/s")

        if throughput < 500:
            print(f"[WARNING] Throughput low (<500 mol/s)")

    # Analyze log
    print(f"\nLog analysis:")
    with open(log_file) as f:
        log_content = f.read()

    # Check buffer sizes
    buffer_sizes = re.findall(r'buffer_full \(([0-9,]+) products\)', log_content)
    if buffer_sizes:
        buffer_nums = [int(b.replace(',', '')) for b in buffer_sizes]
        max_buffer = max(buffer_nums)
        avg_buffer = sum(buffer_nums) / len(buffer_nums)

        print(f"  Buffer flushes: {len(buffer_sizes)}")
        print(f"  Max buffer size: {max_buffer:,} products")
        print(f"  Avg buffer size: {avg_buffer:,.0f} products")

        if max_buffer > 3000:
            print(f"  [FAILED] Buffer exceeded safe limit")
            return False
        elif max_buffer > 2500:
            print(f"  [WARNING] Buffer slightly high (acceptable if <3000)")
        else:
            print(f"  [OK] Buffer within limit")

    # Check memory
    critical_mem = log_content.count('Critical system memory')
    if critical_mem > 0:
        print(f"  [WARNING] {critical_mem} critical memory warnings")
        if critical_mem > 5:
            print(f"  [FAILED] Too many warnings")
            return False
    else:
        print(f"  [OK] No critical memory warnings")

    # Check errors
    errors = log_content.count('ERROR')
    if errors > 0:
        print(f"  [FAILED] {errors} errors in log")
        return False
    else:
        print(f"  [OK] No errors")

    return True

def main():
    print("="*80)
    print("AUTOMATED GRADUAL VALIDATION")
    print("="*80)
    print()
    print("Will run SMALL and MEDIUM tests sequentially.")
    print("MICRO test already passed (43.6% peak memory, 2000 buffer).")
    print()

    tests = [
        ("SMALL",  150000,  6, 900,  65.0),   # 3 batches, 15 min, 6 workers
        ("MEDIUM", 500000,  8, 1800, 68.0),   # 10 batches, 30 min, 8 workers
    ]

    results = {}

    for level, max_rows, workers, timeout, target_mem in tests:
        print(f"\n{'='*80}")
        print(f"Starting {level} test...")
        print(f"{'='*80}\n")

        success = run_test(level, max_rows, workers, timeout, target_mem)
        results[level] = success

        if not success:
            print("\n" + "="*80)
            print(f"TEST {level} FAILED - STOPPING")
            print("="*80)
            return False

        print("\n" + "="*80)
        print(f"TEST {level} PASSED")
        print("="*80)

    # All passed
    print("\n" + "="*80)
    print("ALL TESTS PASSED!")
    print("="*80)
    print("\nValidation Results:")
    print("  [OK] MICRO  (50K rows)")
    for level, success in results.items():
        status = "[OK]" if success else "[FAILED]"
        rows = "150K" if level == "SMALL" else "500K"
        print(f"  {status} {level} ({rows} rows)")

    print("\n" + "="*80)
    print("RECOMMENDATION: PROCEED WITH FULL PRODUCTION RUN")
    print("="*80)
    print("\nNext step:")
    print("python scripts/08_transform_library_v2.py apply \\")
    print("  --input data/output/nplike_v2/polyphenol-2X/products.parquet \\")
    print("  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe \\")
    print("  --xf-config configs/transforms.yaml \\")
    print("  --xf-name FG_PHENOL_OH__OH__TO__OMe \\")
    print("  --use-bloom-filter \\")
    print("  --workers 16 \\")
    print("  --target-memory 70.0")

    return True

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
