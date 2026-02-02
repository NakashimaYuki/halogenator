#!/usr/bin/env python3
"""
Adaptive re-chunking based on complexity analysis.

This script re-chunks the polyphenol-2X dataset to create chunks with
balanced computational complexity rather than equal row counts.
"""

import pyarrow.parquet as pq
import pyarrow as pa
import json
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors


def estimate_row_complexity(smiles):
    """Estimate complexity of a single molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0

    # Count phenolic OH groups
    pattern = Chem.MolFromSmarts('[OH][c]')
    phenolic_oh = len(mol.GetSubstructMatches(pattern)) if pattern else 0

    # Complexity approximation: phenolic_OH^2
    return phenolic_oh ** 2


def analyze_full_dataset(input_file, sample_rate=0.01):
    """
    Analyze the full dataset with sampling to estimate complexity distribution.

    Args:
        input_file: Path to input parquet file
        sample_rate: Fraction of rows to sample for analysis (default 1%)

    Returns:
        List of complexity scores for sampled rows
    """
    print(f"Analyzing dataset: {input_file}")
    print(f"Sample rate: {sample_rate*100:.1f}%")

    table = pq.read_table(input_file)
    total_rows = len(table)

    # Sample rows
    sample_size = int(total_rows * sample_rate)
    sample_indices = np.random.choice(total_rows, sample_size, replace=False)
    sample_indices = sorted(sample_indices)  # Sort for sequential access

    print(f"Total rows: {total_rows:,}")
    print(f"Sampling {sample_size:,} rows for complexity estimation...")

    complexities = []
    row_positions = []

    print("Analyzing samples...")
    for i, idx in enumerate(sample_indices):
        if (i + 1) % 1000 == 0 or (i + 1) == len(sample_indices):
            print(f"  Progress: {i+1:,}/{len(sample_indices):,} ({100*(i+1)/len(sample_indices):.1f}%)")
        smiles = table['smiles'][idx].as_py()
        complexity = estimate_row_complexity(smiles)
        complexities.append(complexity)
        row_positions.append(idx)

    return row_positions, complexities, total_rows


def calculate_adaptive_splits(row_positions, complexities, total_rows, target_complexity):
    """
    Calculate optimal split points based on cumulative complexity.

    Args:
        row_positions: List of sampled row indices
        complexities: List of complexity scores for sampled rows
        total_rows: Total number of rows in dataset
        target_complexity: Target complexity per chunk

    Returns:
        List of (start_row, end_row) tuples for each chunk
    """
    print(f"\nCalculating adaptive splits...")
    print(f"Target complexity per chunk: {target_complexity:.1f}")

    # Estimate complexity density (complexity per row)
    # Use rolling average to smooth out variations
    window_size = max(100, len(complexities) // 100)
    complexity_density = np.convolve(complexities, np.ones(window_size)/window_size, mode='valid')

    # Extend density estimate to all rows
    density_positions = row_positions[window_size//2:len(complexity_density)+window_size//2]

    # Interpolate to get density for all positions
    all_positions = np.arange(total_rows)
    estimated_density = np.interp(all_positions, density_positions, complexity_density)

    # Calculate cumulative complexity
    cumulative_complexity = np.cumsum(estimated_density)
    total_complexity = cumulative_complexity[-1]

    print(f"Total estimated complexity: {total_complexity:.1f}")
    print(f"Estimated number of chunks: {int(np.ceil(total_complexity / target_complexity))}")

    # Find split points
    splits = [0]  # Start with row 0
    current_target = target_complexity

    while current_target < total_complexity:
        # Find row where cumulative complexity exceeds current target
        split_idx = np.searchsorted(cumulative_complexity, current_target)
        if split_idx >= total_rows:
            break
        splits.append(split_idx)
        current_target += target_complexity

    splits.append(total_rows)  # End with last row

    # Create chunk ranges
    chunks = []
    for i in range(len(splits) - 1):
        start = splits[i]
        end = splits[i + 1]
        chunk_complexity = cumulative_complexity[end-1] - (cumulative_complexity[start-1] if start > 0 else 0)
        chunks.append({
            'chunk_id': i,
            'start_row': start,
            'end_row': end,
            'num_rows': end - start,
            'estimated_complexity': chunk_complexity
        })

    return chunks


def create_adaptive_chunks(input_file, output_dir, chunks_metadata):
    """
    Create new chunk files based on calculated splits.

    Args:
        input_file: Path to input parquet file
        output_dir: Directory to save new chunks
        chunks_metadata: List of chunk metadata from calculate_adaptive_splits
    """
    print(f"\nCreating adaptive chunks...")
    print(f"Output directory: {output_dir}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read full table
    print("Reading full dataset...")
    table = pq.read_table(input_file)

    # Create chunks
    print("Writing chunks...")
    for chunk_meta in chunks_metadata:
        chunk_id = chunk_meta['chunk_id']
        start = chunk_meta['start_row']
        end = chunk_meta['end_row']

        # Slice table
        chunk_table = table.slice(start, end - start)

        # Write chunk
        output_file = output_dir / f"chunk_{chunk_id:03d}_input.parquet"
        pq.write_table(chunk_table, output_file)

        print(f"  Chunk {chunk_id}: rows {start:,}-{end:,} ({end-start:,} rows, "
              f"complexity ~{chunk_meta['estimated_complexity']:.1f}) -> {output_file.name}")

    # Save metadata
    metadata_file = output_dir / "chunks_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump({
            'input_file': str(input_file),
            'total_rows': len(table),
            'num_chunks': len(chunks_metadata),
            'chunks': chunks_metadata
        }, f, indent=2)

    print(f"\n[OK] Metadata saved to: {metadata_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Adaptive re-chunking based on complexity')
    parser.add_argument('--input', required=True, help='Input parquet file')
    parser.add_argument('--output-dir', required=True, help='Output directory for new chunks')
    parser.add_argument('--target-complexity', type=float, default=150.0,
                       help='Target complexity per chunk (default: 150)')
    parser.add_argument('--sample-rate', type=float, default=0.01,
                       help='Sample rate for complexity estimation (default: 0.01 = 1%%)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Calculate splits without creating files')

    args = parser.parse_args()

    print("="*80)
    print("ADAPTIVE RE-CHUNKING")
    print("="*80)

    # Step 1: Analyze dataset
    row_positions, complexities, total_rows = analyze_full_dataset(
        args.input,
        sample_rate=args.sample_rate
    )

    # Step 2: Calculate splits
    chunks_metadata = calculate_adaptive_splits(
        row_positions,
        complexities,
        total_rows,
        args.target_complexity
    )

    # Print summary
    print("\n" + "="*80)
    print("PROPOSED CHUNKING PLAN")
    print("="*80)
    print(f"Number of chunks: {len(chunks_metadata)}")
    print(f"\nChunk details:")
    for chunk in chunks_metadata:
        print(f"  Chunk {chunk['chunk_id']:2d}: "
              f"{chunk['num_rows']:8,} rows, "
              f"complexity ~{chunk['estimated_complexity']:6.1f}")

    avg_complexity = np.mean([c['estimated_complexity'] for c in chunks_metadata])
    std_complexity = np.std([c['estimated_complexity'] for c in chunks_metadata])
    print(f"\nComplexity statistics:")
    print(f"  Mean: {avg_complexity:.1f}")
    print(f"  Std:  {std_complexity:.1f}")
    print(f"  CV:   {std_complexity/avg_complexity:.2%}")

    # Step 3: Create chunks (if not dry-run)
    if args.dry_run:
        print("\n[DRY RUN] Skipping chunk creation")
    else:
        create_adaptive_chunks(args.input, args.output_dir, chunks_metadata)
        print("\n[OK] Adaptive re-chunking complete!")


if __name__ == '__main__':
    main()
