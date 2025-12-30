#!/usr/bin/env python3
"""
Simple shuffle-based re-chunking.

Strategy: Randomly shuffle the dataset, then split into equal-sized chunks.
This distributes complex and simple molecules evenly across all chunks.
"""

import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np
import json
from pathlib import Path


def shuffle_and_rechunk(input_file, output_dir, chunk_size=1_000_000, dry_run=False):
    """
    Shuffle dataset and create fixed-size chunks.

    Args:
        input_file: Path to input parquet file
        output_dir: Directory to save new chunks
        chunk_size: Number of rows per chunk (default: 1,000,000)
        dry_run: If True, only show plan without creating files
    """
    print("="*80)
    print("SHUFFLE-BASED RE-CHUNKING")
    print("="*80)

    # Read dataset
    print(f"Reading dataset: {input_file}")
    table = pq.read_table(input_file)
    total_rows = len(table)

    print(f"Total rows: {total_rows:,}")
    print(f"Chunk size: {chunk_size:,} rows")

    num_chunks = int(np.ceil(total_rows / chunk_size))
    print(f"Number of chunks: {num_chunks}")

    # Generate random permutation
    print("\nGenerating random shuffle permutation...")
    np.random.seed(42)  # For reproducibility
    shuffle_indices = np.random.permutation(total_rows)

    print("Shuffle indices generated")

    # Plan chunks
    chunks_metadata = []
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, total_rows)
        num_rows = end_idx - start_idx

        chunks_metadata.append({
            'chunk_id': i,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'num_rows': num_rows
        })

    # Print plan
    print("\n" + "="*80)
    print("CHUNKING PLAN")
    print("="*80)
    for chunk in chunks_metadata:
        print(f"  Chunk {chunk['chunk_id']:2d}: {chunk['num_rows']:8,} rows")

    avg_size = total_rows / num_chunks
    print(f"\nAverage chunk size: {avg_size:,.0f} rows")
    print(f"Expected complexity variance: ~{100/np.sqrt(num_chunks):.1f}% (by Central Limit Theorem)")

    if dry_run:
        print("\n[DRY RUN] Skipping chunk creation")
        return

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nCreating chunks in: {output_dir}")

    # Create chunks
    for chunk_meta in chunks_metadata:
        chunk_id = chunk_meta['chunk_id']
        start = chunk_meta['start_idx']
        end = chunk_meta['end_idx']

        # Get shuffled indices for this chunk
        chunk_indices = shuffle_indices[start:end]

        # Take rows from original table using shuffled indices
        chunk_table = table.take(chunk_indices)

        # Write chunk
        output_file = output_dir / f"chunk_{chunk_id:03d}_input.parquet"
        pq.write_table(chunk_table, output_file)

        file_size_mb = output_file.stat().st_size / (1024**2)
        print(f"  Chunk {chunk_id:2d}: {len(chunk_table):8,} rows -> {output_file.name} ({file_size_mb:.1f} MB)")

    # Save metadata
    metadata_file = output_dir / "chunks_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump({
            'input_file': str(input_file),
            'total_rows': total_rows,
            'chunk_size': chunk_size,
            'num_chunks': num_chunks,
            'shuffle_seed': 42,
            'chunks': chunks_metadata
        }, f, indent=2)

    print(f"\n[OK] Metadata saved to: {metadata_file}")
    print("\n[OK] Shuffle-based re-chunking complete!")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Shuffle-based re-chunking')
    parser.add_argument('--input', required=True, help='Input parquet file')
    parser.add_argument('--output-dir', required=True, help='Output directory for new chunks')
    parser.add_argument('--chunk-size', type=int, default=1_000_000,
                       help='Number of rows per chunk (default: 1,000,000)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show plan without creating files')

    args = parser.parse_args()

    shuffle_and_rechunk(
        args.input,
        args.output_dir,
        chunk_size=args.chunk_size,
        dry_run=args.dry_run
    )


if __name__ == '__main__':
    main()
