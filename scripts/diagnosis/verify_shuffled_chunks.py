#!/usr/bin/env python3
"""
Verify that shuffled chunks have balanced complexity distribution.
"""

import pyarrow.parquet as pq
from rdkit import Chem
import numpy as np
from pathlib import Path
import json


def estimate_chunk_complexity(chunk_file, sample_size=1000):
    """Estimate complexity of a chunk by sampling."""

    table = pq.read_table(chunk_file)
    total_rows = len(table)

    # Sample rows
    sample_indices = np.random.choice(total_rows, min(sample_size, total_rows), replace=False)
    smiles_list = [table['smiles'][i].as_py() for i in sample_indices]

    complexities = []
    phenolic_oh_counts = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        # Count phenolic OH groups
        pattern = Chem.MolFromSmarts('[OH][c]')
        phenolic_oh = len(mol.GetSubstructMatches(pattern)) if pattern else 0

        complexities.append(phenolic_oh ** 2)
        phenolic_oh_counts.append(phenolic_oh)

    return {
        'total_rows': total_rows,
        'sample_size': len(complexities),
        'avg_complexity': float(np.mean(complexities)),
        'std_complexity': float(np.std(complexities)),
        'avg_phenolic_oh': float(np.mean(phenolic_oh_counts)),
        'max_phenolic_oh': int(np.max(phenolic_oh_counts)) if phenolic_oh_counts else 0
    }


def verify_shuffled_chunks(chunks_dir):
    """Verify complexity distribution across shuffled chunks."""

    chunks_dir = Path(chunks_dir)

    print("="*80)
    print("SHUFFLED CHUNKS COMPLEXITY VERIFICATION")
    print("="*80)
    print(f"Chunks directory: {chunks_dir}\n")

    chunk_files = sorted(chunks_dir.glob("chunk_*_input.parquet"))

    if not chunk_files:
        print(f"No chunk files found in {chunks_dir}")
        return

    print(f"Found {len(chunk_files)} chunks\n")

    results = []

    for chunk_file in chunk_files:
        chunk_id = int(chunk_file.stem.split('_')[1])
        print(f"Analyzing Chunk {chunk_id}...")

        stats = estimate_chunk_complexity(chunk_file)
        stats['chunk_id'] = chunk_id

        results.append(stats)

        print(f"  Rows: {stats['total_rows']:,}")
        print(f"  Avg complexity: {stats['avg_complexity']:.1f}")
        print(f"  Avg phenolic OH: {stats['avg_phenolic_oh']:.2f}")
        print(f"  Max phenolic OH: {stats['max_phenolic_oh']}")

    # Calculate statistics across all chunks
    complexities = [r['avg_complexity'] for r in results]
    phenolic_ohs = [r['avg_phenolic_oh'] for r in results]

    print("\n" + "="*80)
    print("CROSS-CHUNK STATISTICS")
    print("="*80)

    print(f"\nComplexity Distribution:")
    print(f"  Mean: {np.mean(complexities):.1f}")
    print(f"  Std:  {np.std(complexities):.1f}")
    print(f"  Min:  {np.min(complexities):.1f} (Chunk {np.argmin(complexities)})")
    print(f"  Max:  {np.max(complexities):.1f} (Chunk {np.argmax(complexities)})")
    print(f"  Range: {np.max(complexities) - np.min(complexities):.1f}")
    print(f"  CV (coefficient of variation): {np.std(complexities) / np.mean(complexities):.2%}")

    print(f"\nPhenolic OH Distribution:")
    print(f"  Mean: {np.mean(phenolic_ohs):.2f}")
    print(f"  Std:  {np.std(phenolic_ohs):.2f}")
    print(f"  Range: {np.max(phenolic_ohs) - np.min(phenolic_ohs):.2f}")

    # Compare with original non-shuffled chunks
    print(f"\n" + "="*80)
    print("COMPARISON WITH ORIGINAL (NON-SHUFFLED) CHUNKS")
    print("="*80)

    # Load original analysis if available
    original_file = Path('chunk_complexity_analysis.json')
    if original_file.exists():
        with open(original_file, 'r') as f:
            original_data = json.load(f)

        original_complexities = [r['complexity_score']['mean'] for r in original_data]

        print(f"\nOriginal chunks:")
        print(f"  Complexity range: {np.min(original_complexities):.1f} - {np.max(original_complexities):.1f}")
        print(f"  Std: {np.std(original_complexities):.1f}")
        print(f"  CV: {np.std(original_complexities) / np.mean(original_complexities):.2%}")

        print(f"\nShuffled chunks:")
        print(f"  Complexity range: {np.min(complexities):.1f} - {np.max(complexities):.1f}")
        print(f"  Std: {np.std(complexities):.1f}")
        print(f"  CV: {np.std(complexities) / np.mean(complexities):.2%}")

        improvement = (np.std(original_complexities) - np.std(complexities)) / np.std(original_complexities)
        print(f"\n[OK] Variance reduction: {improvement:.1%}")

        print(f"\nOriginal: {len([c for c in original_complexities if c > 200])} chunks with complexity > 200")
        print(f"Shuffled: {len([c for c in complexities if c > 200])} chunks with complexity > 200")

    else:
        print("Original analysis not found, skipping comparison")

    # Save results
    output_file = chunks_dir / 'shuffled_complexity_verification.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n[OK] Verification results saved to: {output_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Verify shuffled chunks complexity')
    parser.add_argument('--chunks-dir', required=True, help='Directory containing shuffled chunks')

    args = parser.parse_args()

    verify_shuffled_chunks(args.chunks_dir)


if __name__ == '__main__':
    main()
