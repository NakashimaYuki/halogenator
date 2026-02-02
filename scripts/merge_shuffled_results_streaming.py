#!/usr/bin/env python3
"""
Ultra Low Memory Streaming Merge

Uses record batch iteration to minimize memory usage.
Memory usage: ~100-200MB peak (instead of 10GB+).
"""

import pyarrow.parquet as pq
import pyarrow as pa
from pathlib import Path
import json
from datetime import datetime
import sys
import gc


def merge_chunks_ultra_low_memory(chunks_dir, output_file, num_chunks=28, batch_rows=50000):
    """
    Ultra low memory merge using record batch iteration.

    Args:
        chunks_dir: Directory containing chunk outputs
        output_file: Final merged parquet file
        num_chunks: Total number of chunks
        batch_rows: Rows per batch (lower = less memory)
    """

    print(f"\n{'='*80}")
    print(f"ULTRA LOW MEMORY STREAMING MERGE")
    print(f"{'='*80}")
    print(f"Chunks directory: {chunks_dir}")
    print(f"Expected chunks: {num_chunks}")
    print(f"Batch size: {batch_rows:,} rows")
    print(f"Output file: {output_file}")
    print(f"Memory strategy: Record batch iteration (~100-200MB peak)")

    chunks_dir = Path(chunks_dir)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Verify all chunks exist and get total count
    print(f"\nStep 1: Verifying chunks...")
    chunk_files = []
    total_expected = 0

    for i in range(num_chunks):
        chunk_output_dir = chunks_dir / f"chunk_{i:03d}_output"
        products_file = chunk_output_dir / "products.parquet"

        if not products_file.exists():
            print(f"  [X] Chunk {i:03d}: Missing products.parquet")
            return False

        # Get row count from metadata (no loading)
        pf = pq.ParquetFile(products_file)
        num_rows = pf.metadata.num_rows
        total_expected += num_rows

        chunk_files.append({
            'chunk_id': i,
            'file': products_file,
            'num_rows': num_rows
        })

        print(f"  [OK] Chunk {i:03d}: {num_rows:,} rows")

    print(f"\nTotal expected rows: {total_expected:,}")

    # Step 2: Get schema from first chunk (minimal read)
    print(f"\nStep 2: Reading schema...")
    pf = pq.ParquetFile(chunk_files[0]['file'])
    schema = pf.schema_arrow
    print(f"  Schema: {len(schema)} columns")

    # Step 3: Stream write using record batches
    print(f"\nStep 3: Streaming write using record batches...")
    print(f"  (Each chunk read in {batch_rows:,}-row batches)")

    total_written = 0

    with pq.ParquetWriter(output_file, schema, compression='snappy') as writer:
        for idx, chunk_info in enumerate(chunk_files):
            chunk_id = chunk_info['chunk_id']
            chunk_rows = chunk_info['num_rows']

            print(f"\n  [{idx+1}/{len(chunk_files)}] Chunk {chunk_id:03d} ({chunk_rows:,} rows)...")

            # Open parquet file for batch iteration
            pf = pq.ParquetFile(chunk_info['file'])

            # Iterate through row groups
            batch_count = 0
            chunk_written = 0

            for batch in pf.iter_batches(batch_size=batch_rows):
                # Convert to table and write
                batch_table = pa.Table.from_batches([batch], schema=schema)
                writer.write_table(batch_table)

                chunk_written += len(batch)
                batch_count += 1

                # Cleanup immediately
                del batch_table
                del batch

            total_written += chunk_written
            print(f"    Written: {chunk_written:,} rows in {batch_count} batches")

            # Force garbage collection after each chunk
            gc.collect()

    file_size_gb = output_file.stat().st_size / (1024**3)

    print(f"\n  [OK] Streaming write complete!")
    print(f"    Total rows written: {total_written:,}")
    print(f"    File size: {file_size_gb:.2f} GB")

    # Step 4: Create summary
    print(f"\nStep 4: Creating summary...")

    merge_summary = {
        "timestamp": datetime.now().isoformat(),
        "num_chunks_merged": len(chunk_files),
        "total_rows_written": total_written,
        "expected_rows": total_expected,
        "output_file": str(output_file),
        "file_size_gb": file_size_gb,
        "batch_size": batch_rows
    }

    summary_file = output_file.parent / "FINAL_MERGE_SUMMARY.json"
    with open(summary_file, 'w') as f:
        json.dump(merge_summary, f, indent=2)

    print(f"  Summary saved: {summary_file}")

    # Final summary
    print(f"\n{'='*80}")
    print(f"MERGE COMPLETE")
    print(f"{'='*80}")
    print(f"Chunks merged: {len(chunk_files)}/{num_chunks}")
    print(f"Total rows: {total_written:,}")
    print(f"Output: {output_file}")
    print(f"Size: {file_size_gb:.2f} GB")

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Ultra low memory streaming merge")
    parser.add_argument(
        "--chunks-dir",
        default="data/output/transforms/SHUFFLED_RESULTS",
        help="Directory containing chunk outputs"
    )
    parser.add_argument(
        "--output",
        default="data/output/transforms/polyphenol-2X_FINAL_MERGED.parquet",
        help="Output merged parquet file"
    )
    parser.add_argument(
        "--num-chunks",
        type=int,
        default=28,
        help="Expected number of chunks"
    )
    parser.add_argument(
        "--batch-rows",
        type=int,
        default=50000,
        help="Rows per batch (lower = less memory, default 50000)"
    )

    args = parser.parse_args()

    success = merge_chunks_ultra_low_memory(
        chunks_dir=args.chunks_dir,
        output_file=args.output,
        num_chunks=args.num_chunks,
        batch_rows=args.batch_rows
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
