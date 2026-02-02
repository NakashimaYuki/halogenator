#!/usr/bin/env python3
"""
Shuffled Chunks Transform Pipeline
Processes pre-created shuffled chunks to avoid timeout issues.

Uses existing 28 shuffled chunks (500K rows each) from:
  data/output/transforms/polyphenol-2X_SHUFFLED_500K/

Based on shuffle-rechunking solution that reduced complexity variance by 94.7%
"""

import subprocess
import pyarrow.parquet as pq
import pyarrow as pa
from pathlib import Path
import json
import time
import sys
from datetime import datetime


class ShuffledChunksPipeline:
    def __init__(
        self,
        chunks_dir,
        output_dir,
        transform_name,
        workers=16,
        max_in_flight=6,
        batch_size=50000,
        timeout=14400,  # 4 hours
    ):
        self.chunks_dir = Path(chunks_dir)
        self.output_dir = Path(output_dir)
        self.transform_name = transform_name
        self.workers = workers
        self.max_in_flight = max_in_flight
        self.batch_size = batch_size
        self.timeout = timeout

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

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

    def discover_chunks(self):
        """Discover all pre-created shuffled chunks."""
        print(f"\n{'='*80}")
        print(f"STEP 1: Discovering shuffled chunks")
        print(f"{'='*80}")

        # Find all chunk files
        chunk_files = sorted(self.chunks_dir.glob("chunk_*_input.parquet"))

        if not chunk_files:
            print(f"ERROR: No chunks found in {self.chunks_dir}")
            sys.exit(1)

        print(f"Chunks directory: {self.chunks_dir}")
        print(f"Found {len(chunk_files)} chunk files")

        # Load metadata if available
        metadata_file = self.chunks_dir / "chunks_metadata.json"
        if metadata_file.exists():
            with open(metadata_file) as f:
                metadata = json.load(f)
            print(f"\nMetadata:")
            print(f"  Total rows: {metadata['total_rows']:,}")
            print(f"  Chunk size: {metadata['chunk_size']:,}")
            print(f"  Shuffle seed: {metadata['shuffle_seed']}")

        # Create chunk info list
        chunks = []
        total_rows = 0

        for chunk_file in chunk_files:
            # Extract chunk ID from filename (chunk_000_input.parquet -> 0)
            chunk_id = int(chunk_file.stem.split('_')[1])

            # Read chunk to get row count
            table = pq.read_table(chunk_file)
            num_rows = len(table)
            total_rows += num_rows

            chunk_info = {
                "chunk_id": chunk_id,
                "num_rows": num_rows,
                "input_file": str(chunk_file),
                "output_dir": str(self.output_dir / f"chunk_{chunk_id:03d}_output"),
                "status": "pending"
            }

            chunks.append(chunk_info)

            print(f"  Chunk {chunk_id:03d}: {num_rows:,} rows ({chunk_file.stat().st_size / (1024**2):.1f} MB)")

        print(f"\nTotal rows across all chunks: {total_rows:,}")

        # Update state
        self.state["chunks"] = chunks
        self.save_state()

        return chunks

    def process_chunk(self, chunk_info):
        """Process a single chunk."""
        chunk_id = chunk_info["chunk_id"]

        print(f"\n{'='*80}")
        print(f"Processing Chunk {chunk_id:03d}")
        print(f"{'='*80}")
        print(f"  Input: {chunk_info['input_file']}")
        print(f"  Rows: {chunk_info['num_rows']:,}")
        print(f"  Workers: {self.workers}")
        print(f"  Max-in-flight: {self.max_in_flight}")
        print(f"  Timeout: {self.timeout}s ({self.timeout/3600:.1f}h)")

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
            "--bloom-expected-items", "5000000",  # Per chunk dedup
        ]

        # Run with logging
        log_file = output_dir / "transform.log"
        start_time = time.time()

        print(f"\n  Starting transform... (log: {log_file})")

        try:
            with open(log_file, 'w', encoding='utf-8') as log:
                result = subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout
                )

            elapsed = time.time() - start_time

            if result.returncode == 0:
                # Read summary
                summary_file = output_dir / "SUMMARY.json"
                if summary_file.exists():
                    with open(summary_file) as f:
                        summary = json.load(f)

                    print(f"\n  [OK] SUCCESS ({elapsed:.1f}s = {elapsed/3600:.2f}h)")
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
                # Read last 20 lines of log for diagnostics
                if log_file.exists():
                    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                        lines = f.readlines()
                        print("\n  Last 20 lines of log:")
                        for line in lines[-20:]:
                            print(f"    {line.rstrip()}")

                chunk_info["status"] = "failed"
                chunk_info["error"] = f"returncode_{result.returncode}"
                return False

        except subprocess.TimeoutExpired:
            elapsed = time.time() - start_time
            print(f"\n  [X] TIMEOUT after {elapsed:.1f}s ({elapsed/3600:.2f}h)")
            print(f"    This should NOT happen with shuffled chunks!")
            print(f"    Check log: {log_file}")

            chunk_info["status"] = "failed"
            chunk_info["error"] = "timeout"
            chunk_info["elapsed_seconds"] = elapsed
            return False

        except Exception as e:
            elapsed = time.time() - start_time
            print(f"\n  [X] ERROR: {e}")
            chunk_info["status"] = "failed"
            chunk_info["error"] = str(e)
            chunk_info["elapsed_seconds"] = elapsed
            return False

    def process_all_chunks(self):
        """Process all chunks sequentially."""
        print(f"\n{'='*80}")
        print(f"STEP 2: Processing all chunks")
        print(f"{'='*80}")

        self.state["start_time"] = datetime.now().isoformat()
        self.save_state()

        chunks = self.state["chunks"]
        completed = []
        failed = []

        # Calculate initial stats from state
        for chunk in chunks:
            if chunk.get("status") == "completed":
                completed.append(chunk)
            elif chunk.get("status") == "failed":
                failed.append(chunk)

        print(f"\nInitial state:")
        print(f"  Completed: {len(completed)}/{len(chunks)}")
        print(f"  Failed: {len(failed)}/{len(chunks)}")
        print(f"  Remaining: {len(chunks) - len(completed) - len(failed)}/{len(chunks)}")

        for i, chunk_info in enumerate(chunks, 1):
            # Skip already completed chunks
            if chunk_info.get("status") == "completed":
                print(f"\n[Progress: {i}/{len(chunks)}] Chunk {chunk_info['chunk_id']:03d}: Already completed, skipping")
                continue

            # Ask user if they want to retry failed chunks
            if chunk_info.get("status") == "failed":
                print(f"\n[Progress: {i}/{len(chunks)}] Chunk {chunk_info['chunk_id']:03d}: Previously failed")
                print(f"  Retrying...")

            print(f"\n[Progress: {i}/{len(chunks)}]")

            success = self.process_chunk(chunk_info)

            # Update tracking lists
            if success:
                if chunk_info not in completed:
                    completed.append(chunk_info)
                if chunk_info in failed:
                    failed.remove(chunk_info)
                self.state["completed_chunks"] = [c["chunk_id"] for c in completed]
            else:
                if chunk_info not in failed:
                    failed.append(chunk_info)
                if chunk_info in completed:
                    completed.remove(chunk_info)
                self.state["failed_chunks"] = [c["chunk_id"] for c in failed]

            self.save_state()

            # Summary after each chunk
            remaining = len(chunks) - len(completed) - len(failed)
            print(f"\n  Progress Summary:")
            print(f"    Completed: {len(completed)}/{len(chunks)}")
            print(f"    Failed: {len(failed)}/{len(chunks)}")
            print(f"    Remaining: {remaining}/{len(chunks)}")

            if completed:
                avg_time = sum(c.get("elapsed_seconds", 0) for c in completed) / len(completed)
                estimated_remaining = avg_time * remaining / 3600
                print(f"    Avg time per chunk: {avg_time:.1f}s ({avg_time/3600:.2f}h)")
                print(f"    Estimated time remaining: {estimated_remaining:.1f}h ({estimated_remaining/24:.1f}d)")

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

        for chunk_info in sorted(completed_chunks, key=lambda x: x["chunk_id"]):
            output_dir = Path(chunk_info["output_dir"])
            products_file = output_dir / "products.parquet"

            if products_file.exists():
                table = pq.read_table(products_file)
                all_tables.append(table)
                total_products += len(table)
                print(f"    Chunk {chunk_info['chunk_id']:03d}: {len(table):,} products")
            else:
                print(f"    [X] Missing: {products_file}")

        if not all_tables:
            print("  [X] No product files found!")
            return False

        # Concatenate all tables
        print(f"\n  Concatenating {len(all_tables)} tables...")
        merged_table = pa.concat_tables(all_tables)

        print(f"  Total products: {len(merged_table):,}")

        # Write merged output
        final_output = self.output_dir / "polyphenol-2X_FINAL_MERGED.parquet"
        print(f"  Writing merged parquet...")
        pq.write_table(merged_table, final_output)

        file_size_gb = final_output.stat().st_size / (1024**3)

        print(f"\n  [OK] Merged output: {final_output}")
        print(f"    Total products: {len(merged_table):,}")
        print(f"    File size: {file_size_gb:.2f} GB")

        # Create final summary
        total_elapsed = sum(c.get("elapsed_seconds", 0) for c in completed_chunks)
        total_hours = total_elapsed / 3600

        summary = {
            "total_chunks": len(completed_chunks),
            "total_products": len(merged_table),
            "total_elapsed_seconds": total_elapsed,
            "total_elapsed_hours": total_hours,
            "average_throughput_mol_per_sec": len(merged_table) / total_elapsed if total_elapsed > 0 else 0,
            "average_time_per_chunk_seconds": total_elapsed / len(completed_chunks),
            "output_file": str(final_output),
            "timestamp": datetime.now().isoformat()
        }

        summary_file = self.output_dir / "FINAL_SUMMARY.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        print(f"\n  [OK] Summary: {summary_file}")
        print(f"    Total processing time: {total_hours:.1f}h ({total_hours/24:.1f}d)")
        print(f"    Avg per chunk: {summary['average_time_per_chunk_seconds']/3600:.2f}h")
        print(f"    Avg throughput: {summary['average_throughput_mol_per_sec']:.1f} mol/s")

        return True

    def run(self):
        """Run the complete pipeline."""
        print(f"\n{'='*80}")
        print(f"SHUFFLED CHUNKS TRANSFORM PIPELINE")
        print(f"{'='*80}")
        print(f"Transform: {self.transform_name}")
        print(f"Chunks directory: {self.chunks_dir}")
        print(f"Output directory: {self.output_dir}")
        print(f"Workers: {self.workers}")
        print(f"Max-in-flight: {self.max_in_flight}")
        print(f"Batch size: {self.batch_size:,}")
        print(f"Timeout per chunk: {self.timeout}s ({self.timeout/3600:.1f}h)")

        # Step 1: Discover chunks
        chunks = self.discover_chunks()

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

        if completed:
            total_time = sum(c.get("elapsed_seconds", 0) for c in completed) / 3600
            print(f"\nTotal processing time: {total_time:.1f}h ({total_time/24:.1f}d)")

        if failed:
            print(f"\nFailed chunks:")
            for chunk in failed:
                print(f"  - Chunk {chunk['chunk_id']:03d}: {chunk.get('error', 'unknown')}")
            print(f"\nTo retry failed chunks, run this script again (checkpoint will resume).")

        return len(failed) == 0


def main():
    """Run pipeline on shuffled chunks."""
    pipeline = ShuffledChunksPipeline(
        chunks_dir="data/output/transforms/polyphenol-2X_SHUFFLED_500K",
        output_dir="data/output/transforms/SHUFFLED_RESULTS",
        transform_name="FG_PHENOL_OH__OH__TO__OMe",
        workers=16,
        max_in_flight=6,
        batch_size=50000,
        timeout=14400,  # 4 hours - should be sufficient with shuffled chunks
    )

    success = pipeline.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
