#!/usr/bin/env python3
"""
Gradual validation test suite for buffer explosion fix.

Test progression:
1. MICRO:  1 batch   (50K rows,  ~2min)  - Verify basic functionality
2. SMALL:  3 batches (150K rows, ~5min)  - Verify buffer behavior
3. MEDIUM: 10 batches(500K rows, ~15min) - Verify sustained operation
4. FULL:   All data  (13.79M rows, 16-20h) - Production validation

Each level must pass before proceeding to next.
"""

import subprocess
import sys
import time
import json
from pathlib import Path
import pyarrow.parquet as pq

class ValidationSuite:
    def __init__(self):
        self.input_file = r"E:\Projects\halogenator\data\output\nplike_v2\polyphenol-2X\products.parquet"
        self.transform_config = r"E:\Projects\halogenator\configs\transforms.yaml"
        self.transform_name = "FG_PHENOL_OH__OH__TO__OMe"
        self.script = r"E:\Projects\halogenator\scripts\08_transform_library_v2.py"

        # Load input to know total size
        self.input_table = pq.read_table(self.input_file)
        self.total_rows = len(self.input_table)
        print(f"Input dataset: {self.total_rows:,} total rows")

    def run_test(self, level, max_rows, timeout_sec, target_memory):
        """Run a single test level."""
        print("\n" + "="*80)
        print(f"TEST LEVEL: {level}")
        print("="*80)

        # Calculate expected batches
        batch_size = 50000
        expected_batches = (max_rows + batch_size - 1) // batch_size

        print(f"Parameters:")
        print(f"  Input rows: {max_rows:,} / {self.total_rows:,}")
        print(f"  Expected batches: {expected_batches}")
        print(f"  Timeout: {timeout_sec//60} minutes")
        print(f"  Target memory: {target_memory}%")
        print()

        # Create output directory
        outdir = Path(r"E:\Projects\halogenator\data\output\transforms") / f"TEST_{level}_polyphenol"
        if outdir.exists():
            import shutil
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Create temporary input file with limited rows
        if max_rows < self.total_rows:
            print(f"Creating subset input ({max_rows:,} rows)...")
            subset = self.input_table.slice(0, max_rows)
            temp_input = outdir / "input_subset.parquet"
            pq.write_table(subset, temp_input)
            input_to_use = str(temp_input)
            print(f"Subset created: {temp_input}")
        else:
            input_to_use = self.input_file

        # Build command
        cmd = [
            sys.executable,
            self.script,
            "apply",
            "--input", input_to_use,
            "--outdir", str(outdir),
            "--xf-config", self.transform_config,
            "--xf-name", self.transform_name,
            "--batch-size", "50000",
            "--workers", "8",
            "--flush-interval", "2000",
            "--target-memory", str(target_memory),
            "--use-bloom-filter",
            "--bloom-expected-items", str(max_rows * 10),  # Estimate 10x products
        ]

        print(f"\nCommand: {' '.join(cmd)}")
        print(f"\nStarting test at {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("-"*80)

        # Run with timeout
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

                # Stream output
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
        except KeyboardInterrupt:
            print(f"\n[INTERRUPTED] Test cancelled by user")
            process.kill()
            return False

        # Validate output
        print("\n" + "="*80)
        print("VALIDATION")
        print("="*80)

        success = True

        # Check output file
        products_file = outdir / "products.parquet"
        if not products_file.exists():
            print("[FAILED] Output file not found")
            return False

        try:
            output_table = pq.read_table(products_file)
            file_size_mb = products_file.stat().st_size / 1024**2
            print(f"[OK] Output file valid: {len(output_table):,} rows, {file_size_mb:.1f} MB")
        except Exception as e:
            print(f"[FAILED] Output file corrupted: {e}")
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

            print(f"\nSummary:")
            print(f"  Processed: {processed:,}")
            print(f"  Products: {products:,}")
            print(f"  Unique: {unique:,}")
            print(f"  Throughput: {throughput:.1f} mol/s")

            if throughput < 1000:
                print(f"[WARNING] Throughput very low (<1000 mol/s)")
                success = False

        # Analyze log for issues
        print(f"\nLog analysis:")
        with open(log_file) as f:
            log_content = f.read()

        # Check for buffer explosions
        import re
        buffer_sizes = re.findall(r'buffer_full \(([0-9,]+) products\)', log_content)
        if buffer_sizes:
            max_buffer = max(int(b.replace(',', '')) for b in buffer_sizes)
            print(f"  Max buffer size: {max_buffer:,} products")

            if max_buffer > 3000:
                print(f"  [FAILED] Buffer exceeded safe limit (3000)")
                success = False
            elif max_buffer > 2000:
                print(f"  [WARNING] Buffer slightly over limit (acceptable if <3000)")
            else:
                print(f"  [OK] Buffer within limit")

        # Check for memory warnings
        critical_mem = log_content.count('Critical system memory')
        if critical_mem > 0:
            print(f"  [WARNING] {critical_mem} critical memory warnings")
            if critical_mem > 5:
                print(f"  [FAILED] Too many critical memory warnings")
                success = False
        else:
            print(f"  [OK] No critical memory warnings")

        # Check for errors
        errors = log_content.count('ERROR')
        if errors > 0:
            print(f"  [FAILED] {errors} errors in log")
            success = False
        else:
            print(f"  [OK] No errors")

        return success

    def run_all(self):
        """Run all test levels in sequence."""
        print("="*80)
        print("GRADUAL VALIDATION TEST SUITE")
        print("="*80)
        print("\nThis will run 3 progressive tests before recommending full run.")
        print("Each level must pass before proceeding to next.\n")

        input("Press Enter to start (Ctrl+C to cancel)...")

        tests = [
            ("MICRO",  50000,   180,  60.0),  # 1 batch, 3 min, 60% target
            ("SMALL",  150000,  600,  65.0),  # 3 batches, 10 min, 65% target
            ("MEDIUM", 500000,  1200, 68.0),  # 10 batches, 20 min, 68% target
        ]

        results = {}

        for level, max_rows, timeout, target_mem in tests:
            success = self.run_test(level, max_rows, timeout, target_mem)
            results[level] = success

            if not success:
                print("\n" + "="*80)
                print(f"TEST {level} FAILED - STOPPING")
                print("="*80)
                print("\nPlease review the logs and fix issues before proceeding.")
                print(f"Log location: data/output/transforms/TEST_{level}_polyphenol/test.log")
                return False

            print("\n" + "="*80)
            print(f"TEST {level} PASSED")
            print("="*80)

            if level != "MEDIUM":
                print("\nProceed to next test level?")
                response = input("Press Enter to continue (Ctrl+C to stop)...")

        # All tests passed
        print("\n" + "="*80)
        print("ALL TESTS PASSED!")
        print("="*80)
        print("\nValidation Results:")
        for level, success in results.items():
            status = "[OK]" if success else "[FAILED]"
            print(f"  {status} {level}")

        print("\n" + "="*80)
        print("RECOMMENDATION: PROCEED WITH FULL PRODUCTION RUN")
        print("="*80)
        print("\nThe fix has been validated on 500K rows with no issues.")
        print("You can now run the full polyphenol-2X transform job:")
        print()
        print("python scripts/08_transform_library_v2.py apply \\")
        print("  --input data/output/nplike_v2/polyphenol-2X/products.parquet \\")
        print("  --outdir data/output/transforms/polyphenol-2X_FG_PHENOL_OH__OH__TO__OMe \\")
        print("  --xf-config configs/transforms.yaml \\")
        print("  --xf-name FG_PHENOL_OH__OH__TO__OMe \\")
        print("  --use-bloom-filter \\")
        print("  --workers 16 \\")
        print("  --target-memory 70.0")
        print()
        print("Expected runtime: 16-20 hours")
        print("Expected memory: < 70% peak")
        print("Expected output: ~89M unique products")

        return True

if __name__ == '__main__':
    suite = ValidationSuite()
    success = suite.run_all()
    sys.exit(0 if success else 1)
