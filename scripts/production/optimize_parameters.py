#!/usr/bin/env python3
"""
Parameter optimization for transform pipeline.

Runs grid search over (workers, max-in-flight, batch-size)
to find optimal configuration balancing throughput and memory.
"""

import subprocess
import pandas as pd
import json
from pathlib import Path
import re
import shutil
import time

def run_test(workers, max_in_flight, batch_size, test_name):
    """
    Run single test with given parameters.

    Returns:
        dict: {
            'workers': int,
            'max_in_flight': int,
            'batch_size': int,
            'peak_memory': float,
            'avg_memory': float,
            'throughput': float,
            'time_seconds': float,
            'unique_products': int,
            'success': bool
        }
    """

    # Create output directory
    outdir = Path(f"data/output/transforms/OPT_{test_name}")
    if outdir.exists():
        try:
            shutil.rmtree(outdir)
        except (PermissionError, OSError) as e:
            # Windows file locking issue - try to continue with existing dir
            print(f"  Warning: Could not remove {outdir}: {e}")
            print(f"  Will overwrite files in existing directory")
    outdir.mkdir(parents=True, exist_ok=True)

    # Create 150K subset for quick testing
    import pyarrow.parquet as pq
    input_file = "data/output/nplike_v2/polyphenol-2X/products.parquet"

    print(f"  Creating 150K subset...")
    table = pq.read_table(input_file)
    subset = table.slice(0, 150000)
    subset_file = outdir / "input_150k.parquet"
    pq.write_table(subset, subset_file)

    # Build command
    cmd = [
        "python",
        "scripts/08_transform_library_v2.py",
        "apply",
        "--input", str(subset_file),
        "--outdir", str(outdir),
        "--xf-config", "configs/transforms.yaml",
        "--xf-name", "FG_PHENOL_OH__OH__TO__OMe",
        "--batch-size", str(batch_size),
        "--workers", str(workers),
        "--flush-interval", "2000",
        "--target-memory", "70.0",
        "--use-bloom-filter",
        "--bloom-expected-items", "2000000",
        "--max-in-flight", str(max_in_flight),
    ]

    # Run test with timeout
    log_file = outdir / "test.log"
    start_time = time.time()

    try:
        print(f"  Running test (timeout: 15min)...")
        with open(log_file, 'w', encoding='utf-8') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=900  # 15 min timeout
            )

        elapsed = time.time() - start_time
        success = result.returncode == 0

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"  [X] TIMEOUT after {elapsed:.1f}s")
        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'peak_memory': None,
            'avg_memory': None,
            'throughput': None,
            'time_seconds': elapsed,
            'unique_products': None,
            'success': False,
            'error': 'timeout'
        }

    # Parse results
    if success:
        # Read summary
        summary_file = outdir / "SUMMARY.json"
        if not summary_file.exists():
            print(f"  [X] FAILED: No SUMMARY.json found")
            return {
                'workers': workers,
                'max_in_flight': max_in_flight,
                'batch_size': batch_size,
                'success': False,
                'error': 'no_summary'
            }

        with open(summary_file, encoding='utf-8') as f:
            summary = json.load(f)

        # Parse log for memory
        with open(log_file, encoding='utf-8', errors='replace') as f:
            log_content = f.read()

        mem_values = re.findall(r'Pre-flush memory: ([0-9.]+)%', log_content)
        mem_floats = [float(m) for m in mem_values] if mem_values else []

        peak_mem = max(mem_floats) if mem_floats else None
        avg_mem = sum(mem_floats) / len(mem_floats) if mem_floats else None

        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'peak_memory': peak_mem,
            'avg_memory': avg_mem,
            'throughput': summary.get('throughput_mol_per_sec'),
            'time_seconds': summary.get('elapsed_seconds', elapsed),
            'unique_products': summary.get('unique_products'),
            'success': True
        }
    else:
        print(f"  [X] FAILED: returncode={result.returncode}")
        return {
            'workers': workers,
            'max_in_flight': max_in_flight,
            'batch_size': batch_size,
            'success': False,
            'error': f'returncode_{result.returncode}'
        }

def estimate_memory(workers, max_in_flight, batch_size):
    """
    Estimate peak memory usage using simplified model.

    Parameters from testing:
    - Worker base: ~1.0 GB per worker
    - Queue: ~0.5 GB per in-flight batch
    - Buffer: ~0.01 GB (fixed after chunking)
    - Overhead: ~6 GB (OS + Python + Bloom)
    """
    worker_mem = workers * 1.0  # GB
    queue_mem = max_in_flight * 0.5  # GB
    buffer_mem = 0.01  # GB
    overhead = 6.0  # GB

    total_gb = worker_mem + queue_mem + buffer_mem + overhead
    total_percent = (total_gb / 32.0) * 100  # Assuming 32GB system

    return total_percent

def main():
    print("="*80)
    print("PARAMETER OPTIMIZATION - Grid Search")
    print("="*80)
    print("\nObjective:")
    print("  - Improve throughput from 600-800 mol/s to 1,200-1,500 mol/s")
    print("  - Keep peak memory < 70% (currently 50-55%)")
    print("  - Reduce runtime from 24h to 16-18h")
    print("\nStrategy:")
    print("  - Test priority parameter combinations on 150K subset")
    print("  - Each test takes ~60-90 minutes")
    print("  - Find optimal (workers, max-in-flight, batch-size)")
    print("\n" + "="*80)

    # Priority tests (from document)
    priority_tests = [
        {'workers': 16, 'max_in_flight': 6, 'batch_size': 50000},  # Increase parallelism
        {'workers': 16, 'max_in_flight': 8, 'batch_size': 50000},  # Further increase
        {'workers': 16, 'max_in_flight': 6, 'batch_size': 25000},  # Smaller batch
        {'workers': 12, 'max_in_flight': 8, 'batch_size': 50000},  # Fewer workers, more batches
    ]

    results = []

    for i, params in enumerate(priority_tests, 1):
        print(f"\n{'='*80}")
        print(f"Test {i}/{len(priority_tests)}")
        print(f"{'='*80}")
        print(f"Parameters:")
        print(f"  workers       = {params['workers']}")
        print(f"  max-in-flight = {params['max_in_flight']}")
        print(f"  batch-size    = {params['batch_size']}")

        # Estimate memory
        est_mem = estimate_memory(params['workers'], params['max_in_flight'], params['batch_size'])
        print(f"  Estimated peak memory: {est_mem:.1f}%")

        if est_mem > 75:
            print(f"  [!] WARNING: Estimated memory > 75%, may be unsafe!")

        test_name = f"w{params['workers']}_m{params['max_in_flight']}_b{params['batch_size']}"

        start = time.time()
        result = run_test(**params, test_name=test_name)
        elapsed = time.time() - start

        results.append(result)

        if result['success']:
            print(f"\n  [OK] SUCCESS ({elapsed:.1f}s total)")
            print(f"    Throughput:    {result['throughput']:.1f} mol/s")
            print(f"    Peak Memory:   {result['peak_memory']:.1f}%")
            print(f"    Avg Memory:    {result['avg_memory']:.1f}%")
            print(f"    Runtime:       {result['time_seconds']:.1f}s")
            print(f"    Products:      {result['unique_products']:,}")

            # Calculate improvement vs baseline (600 mol/s)
            improvement = (result['throughput'] / 600.0 - 1) * 100
            print(f"    Improvement:   {improvement:+.1f}% vs baseline (600 mol/s)")
        else:
            print(f"\n  [X] FAILED ({elapsed:.1f}s)")
            print(f"    Error: {result.get('error', 'unknown')}")

    # Create results DataFrame
    df = pd.DataFrame(results)

    # Save results
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    results_file = f'optimization_results_{timestamp}.csv'
    df.to_csv(results_file, index=False)
    print(f"\n\n{'='*80}")
    print(f"Results saved to: {results_file}")
    print(f"{'='*80}")

    # Filter successful tests
    df_success = df[df['success'] == True]

    if len(df_success) > 0:
        # Sort by throughput
        df_sorted = df_success.sort_values('throughput', ascending=False)

        print("\n" + "="*80)
        print("ALL SUCCESSFUL CONFIGURATIONS (sorted by throughput)")
        print("="*80)
        print(df_sorted[['workers', 'max_in_flight', 'batch_size',
                         'throughput', 'peak_memory', 'avg_memory']].to_string(index=False))

        # Find best config with memory < 70%
        df_safe = df_success[df_success['peak_memory'] < 70.0]
        if len(df_safe) > 0:
            best = df_safe.sort_values('throughput', ascending=False).iloc[0]

            improvement = (best['throughput'] / 600.0 - 1) * 100

            print("\n" + "="*80)
            print("*** RECOMMENDED CONFIGURATION")
            print("   (best throughput with peak memory < 70%)")
            print("="*80)
            print(f"  workers:         {int(best['workers'])}")
            print(f"  max-in-flight:   {int(best['max_in_flight'])}")
            print(f"  batch-size:      {int(best['batch_size'])}")
            print(f"\nExpected Performance:")
            print(f"  Throughput:      {best['throughput']:.1f} mol/s")
            print(f"  Peak Memory:     {best['peak_memory']:.1f}%")
            print(f"  Avg Memory:      {best['avg_memory']:.1f}%")
            print(f"  Improvement:     {improvement:+.1f}% vs baseline (600 mol/s)")

            print(f"\nProduction Command:")
            print(f"python scripts/08_transform_library_v2.py apply \\")
            print(f"  --input data/output/nplike_v2/polyphenol-2X/products.parquet \\")
            print(f"  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe_OPT \\")
            print(f"  --xf-config configs/transforms.yaml \\")
            print(f"  --xf-name FG_PHENOL_OH__OH__TO__OMe \\")
            print(f"  --use-bloom-filter \\")
            print(f"  --workers {int(best['workers'])} \\")
            print(f"  --target-memory 70.0 \\")
            print(f"  --max-in-flight {int(best['max_in_flight'])} \\")
            print(f"  --batch-size {int(best['batch_size'])} \\")
            print(f"  --bloom-expected-items 100000000")

            # Achievement level
            print(f"\n{'='*80}")
            if best['throughput'] >= 1200:
                print("[***] LEVEL 3 ACHIEVED: Ideal target (1200+ mol/s, <70% mem)")
            elif best['throughput'] >= 900:
                print("[ + ] LEVEL 2 ACHIEVED: Expected target (900-1000 mol/s, <70% mem)")
            elif best['throughput'] >= 600:
                print("[OK] LEVEL 1 ACHIEVED: Minimum requirement (>600 mol/s, <70% mem)")
            print("="*80)
        else:
            print("\n[!] WARNING: No configuration found with peak memory < 70%")
            print("All successful tests exceeded memory limit.")
            print("RECOMMENDATION: Keep current baseline (workers=16, max-in-flight=4)")
    else:
        print("\n[X] ERROR: All tests failed!")
        print("No successful configurations found.")
        print("RECOMMENDATION: Keep current baseline (workers=16, max-in-flight=4)")
        print("\nCheck individual test logs in data/output/transforms/OPT_*/test.log")

    print(f"\n{'='*80}")
    print("Optimization complete!")
    print(f"{'='*80}\n")

if __name__ == '__main__':
    main()
