#!/usr/bin/env python3
"""
Batch Transform Pipeline - Divide and Conquer Strategy

Splits large dataset into manageable chunks, processes each independently,
then merges results. Provides checkpoint functionality and crash resilience.
"""

import subprocess
import pyarrow.parquet as pq
import pyarrow as pa
from pathlib import Path
import json
import time
import sys
from datetime import datetime

class BatchTransformPipeline:
    def __init__(
        self,
        input_file,
        output_dir,
        transform_name,
        workers=16,
        max_in_flight=6,
        batch_size=50000,
        rows_per_chunk=1000000,
    ):
        self.input_file = Path(input_file)
        self.output_dir = Path(output_dir)
        self.transform_name = transform_name
        self.workers = workers
        self.max_in_flight = max_in_flight
        self.batch_size = batch_size
        self.rows_per_chunk = rows_per_chunk

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.chunks_dir = self.output_dir / "chunks"
        self.chunks_dir.mkdir(exist_ok=True)

        # State tracking
        self.state_file = self.output_dir / "pipeline_state.json"
        self.state = self.load_state()

    def load_state(self):
        """Load pipeline state for resume capability."""
        if self.state_file.exists():
            with open(self.state_file) as f:
                return json.load(f)
        return {
            "chunks": [],
            "completed_chunks": [],
            "failed_chunks": [],
            "start_time": None,
            "end_time": None
        }

    def save_state(self):
        """Save pipeline state."""
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2)

    def create_chunks(self):
        """Split input parquet into manageable chunks."""
        print(f"\n{'='*80}")
        print(f"STEP 1: Creating data chunks")
        print(f"{'='*80}")

        # Read input file metadata
        table = pq.read_table(self.input_file)
        total_rows = len(table)

        print(f"Input: {self.input_file}")
        print(f"Total rows: {total_rows:,}")
        print(f"Chunk size: {self.rows_per_chunk:,} rows")

        num_chunks = (total_rows + self.rows_per_chunk - 1) // self.rows_per_chunk
        print(f"Number of chunks: {num_chunks}")

        chunks = []
        for i in range(num_chunks):
            start_row = i * self.rows_per_chunk
            end_row = min((i + 1) * self.rows_per_chunk, total_rows)
            chunk_size = end_row - start_row

            chunk_info = {
                "chunk_id": i,
                "start_row": start_row,
                "end_row": end_row,
                "num_rows": chunk_size,
                "input_file": str(self.chunks_dir / f"chunk_{i:03d}_input.parquet"),
                "output_dir": str(self.chunks_dir / f"chunk_{i:03d}_output"),
                "status": "pending"
            }

            # Create chunk input file if not exists
            chunk_input = Path(chunk_info["input_file"])
            if not chunk_input.exists():
                print(f"  Creating chunk {i}: rows {start_row:,}-{end_row:,} ({chunk_size:,} rows)")
                chunk_table = table.slice(start_row, chunk_size)
                pq.write_table(chunk_table, chunk_input)
            else:
                print(f"  Chunk {i} already exists: {chunk_input}")

            chunks.append(chunk_info)

        self.state["chunks"] = chunks
        self.save_state()

        print(f"\n[OK] Created {num_chunks} chunks")
        return chunks

    def process_chunk(self, chunk_info):
        """Process a single chunk."""
        chunk_id = chunk_info["chunk_id"]

        print(f"\n{'='*80}")
        print(f"Processing Chunk {chunk_id}")
        print(f"{'='*80}")
        print(f"  Input: {chunk_info['input_file']}")
        print(f"  Rows: {chunk_info['num_rows']:,}")
        print(f"  Workers: {self.workers}")
        print(f"  Max-in-flight: {self.max_in_flight}")

        # Create output directory
        output_dir = Path(chunk_info["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)

        # Build command
        cmd = [
            "python",
            "scripts/08_transform_library_v2.py",
            "apply",
            "--input", chunk_info["input_file"],
            "--outdir", str(output_dir),
            "--xf-config", "configs/transforms.yaml",
            "--xf-name", self.transform_name,
            "--workers", str(self.workers),
            "--max-in-flight", str(self.max_in_flight),
            "--batch-size", str(self.batch_size),
            "--flush-interval", "2000",
            "--target-memory", "70.0",
            "--use-bloom-filter",
            "--bloom-expected-items", str(chunk_info['num_rows'] * 10),  # Estimate
        ]

        # Run with logging
        log_file = output_dir / "transform.log"
        start_time = time.time()

        try:
            with open(log_file, 'w', encoding='utf-8') as log:
                result = subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=14400  # 4 hour timeout per chunk
                )

            elapsed = time.time() - start_time

            if result.returncode == 0:
                # Read summary
                summary_file = output_dir / "SUMMARY.json"
                if summary_file.exists():
                    with open(summary_file) as f:
                        summary = json.load(f)

                    print(f"\n  [OK] SUCCESS ({elapsed:.1f}s)")
                    print(f"    Products: {summary.get('unique_products', 0):,}")
                    print(f"    Throughput: {summary.get('throughput_mol_per_sec', 0):.1f} mol/s")

                    chunk_info["status"] = "completed"
                    chunk_info["elapsed_seconds"] = elapsed
                    chunk_info["summary"] = summary
                    return True
                else:
                    print(f"\n  [X] FAILED: No SUMMARY.json")
                    chunk_info["status"] = "failed"
                    chunk_info["error"] = "no_summary"
                    return False
            else:
                print(f"\n  [X] FAILED: returncode={result.returncode}")
                chunk_info["status"] = "failed"
                chunk_info["error"] = f"returncode_{result.returncode}"
                return False

        except subprocess.TimeoutExpired:
            print(f"\n  [X] TIMEOUT after 4 hours")
            chunk_info["status"] = "failed"
            chunk_info["error"] = "timeout"
            return False
        except Exception as e:
            print(f"\n  [X] ERROR: {e}")
            chunk_info["status"] = "failed"
            chunk_info["error"] = str(e)
            return False

    def process_all_chunks(self):
        """Process all chunks sequentially."""
        print(f"\n{'='*80}")
        print(f"STEP 2: Processing chunks")
        print(f"{'='*80}")

        self.state["start_time"] = datetime.now().isoformat()
        self.save_state()

        chunks = self.state["chunks"]
        completed = []
        failed = []

        for i, chunk_info in enumerate(chunks, 1):
            # Skip already completed chunks
            if chunk_info.get("status") == "completed":
                print(f"\n[Chunk {chunk_info['chunk_id']}] Already completed, skipping")
                completed.append(chunk_info)
                continue

            print(f"\n[Progress: {i}/{len(chunks)}]")

            success = self.process_chunk(chunk_info)

            if success:
                completed.append(chunk_info)
                self.state["completed_chunks"] = [c["chunk_id"] for c in completed]
            else:
                failed.append(chunk_info)
                self.state["failed_chunks"] = [c["chunk_id"] for c in failed]

            self.save_state()

            # Summary after each chunk
            print(f"\n  Progress: {len(completed)}/{len(chunks)} completed, {len(failed)} failed")

        self.state["end_time"] = datetime.now().isoformat()
        self.save_state()

        return completed, failed

    def merge_results(self, completed_chunks):
        """Merge all chunk results into final output."""
        print(f"\n{'='*80}")
        print(f"STEP 3: Merging results")
        print(f"{'='*80}")

        if not completed_chunks:
            print("  [X] No completed chunks to merge!")
            return False

        print(f"  Merging {len(completed_chunks)} chunks...")

        # Collect all product parquet files
        all_tables = []
        total_products = 0

        for chunk_info in completed_chunks:
            output_dir = Path(chunk_info["output_dir"])
            products_file = output_dir / "products.parquet"

            if products_file.exists():
                table = pq.read_table(products_file)
                all_tables.append(table)
                total_products += len(table)
                print(f"    Chunk {chunk_info['chunk_id']}: {len(table):,} products")
            else:
                print(f"    [X] Missing: {products_file}")

        if not all_tables:
            print("  [X] No product files found!")
            return False

        # Concatenate all tables
        print(f"\n  Concatenating {len(all_tables)} tables...")
        merged_table = pa.concat_tables(all_tables)

        print(f"  Total products before dedup: {len(merged_table):,}")

        # Write merged output
        final_output = self.output_dir / "products_merged.parquet"
        pq.write_table(merged_table, final_output)

        print(f"\n  [OK] Merged output: {final_output}")
        print(f"    Total products: {len(merged_table):,}")
        print(f"    File size: {final_output.stat().st_size / (1024**3):.2f} GB")

        # Create final summary
        total_elapsed = sum(c.get("elapsed_seconds", 0) for c in completed_chunks)

        summary = {
            "total_chunks": len(completed_chunks),
            "total_products": len(merged_table),
            "total_elapsed_seconds": total_elapsed,
            "average_throughput": len(merged_table) / total_elapsed if total_elapsed > 0 else 0,
            "output_file": str(final_output),
            "timestamp": datetime.now().isoformat()
        }

        summary_file = self.output_dir / "FINAL_SUMMARY.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        print(f"\n  [OK] Summary: {summary_file}")

        return True

    def run(self):
        """Run the complete pipeline."""
        print(f"\n{'='*80}")
        print(f"BATCH TRANSFORM PIPELINE")
        print(f"{'='*80}")
        print(f"Transform: {self.transform_name}")
        print(f"Input: {self.input_file}")
        print(f"Output: {self.output_dir}")
        print(f"Workers: {self.workers}")
        print(f"Max-in-flight: {self.max_in_flight}")
        print(f"Chunk size: {self.rows_per_chunk:,} rows")

        # Step 1: Create chunks
        chunks = self.create_chunks()

        # Step 2: Process chunks
        completed, failed = self.process_all_chunks()

        # Step 3: Merge results
        if completed:
            self.merge_results(completed)

        # Final report
        print(f"\n{'='*80}")
        print(f"PIPELINE COMPLETE")
        print(f"{'='*80}")
        print(f"Total chunks: {len(chunks)}")
        print(f"Completed: {len(completed)}")
        print(f"Failed: {len(failed)}")

        if failed:
            print(f"\nFailed chunks:")
            for chunk in failed:
                print(f"  - Chunk {chunk['chunk_id']}: {chunk.get('error', 'unknown')}")
            print(f"\nTo retry failed chunks, run this script again.")

        return len(failed) == 0


def main():
    pipeline = BatchTransformPipeline(
        input_file="data/output/nplike_v2/polyphenol-2X/products.parquet",
        output_dir="data/output/transforms/polyphenol-2X_BATCHED",
        transform_name="FG_PHENOL_OH__OH__TO__OMe",
        workers=16,
        max_in_flight=6,
        batch_size=50000,
        rows_per_chunk=1000000,  # 1M rows per chunk
    )

    success = pipeline.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
