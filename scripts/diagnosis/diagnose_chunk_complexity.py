#!/usr/bin/env python3
"""
Diagnose chunk complexity distribution to understand timeout issues.
"""

import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np
from pathlib import Path

def analyze_chunk_complexity(chunk_id, sample_size=1000):
    """Analyze molecular complexity in a chunk."""

    chunk_file = Path(f"data/output/transforms/polyphenol-2X_BATCHED/chunks/chunk_{chunk_id:03d}_input.parquet")

    if not chunk_file.exists():
        print(f"Chunk {chunk_id} input file not found")
        return None

    table = pq.read_table(chunk_file)
    total_rows = len(table)

    # Sample molecules
    sample_indices = np.random.choice(total_rows, min(sample_size, total_rows), replace=False)
    smiles_list = [table['smiles'][i].as_py() for i in sample_indices]

    metrics = {
        'chunk_id': chunk_id,
        'total_rows': total_rows,
        'phenolic_oh_counts': [],
        'mol_weights': [],
        'num_aromatic_rings': [],
        'num_rotatable_bonds': [],
        'complexity_scores': []
    }

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        # Count phenolic OH groups
        pattern = Chem.MolFromSmarts('[OH][c]')
        phenolic_oh = len(mol.GetSubstructMatches(pattern)) if pattern else 0

        # Molecular descriptors
        mw = Descriptors.MolWt(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)

        # Complexity score (custom metric)
        # Assume each phenolic OH can generate multiple products (per_site mode)
        # Complexity ~ phenolic_oh^k where k depends on transformation
        complexity = phenolic_oh ** 2  # Rough approximation

        metrics['phenolic_oh_counts'].append(phenolic_oh)
        metrics['mol_weights'].append(mw)
        metrics['num_aromatic_rings'].append(aromatic_rings)
        metrics['num_rotatable_bonds'].append(rotatable_bonds)
        metrics['complexity_scores'].append(complexity)

    # Statistics (convert numpy types to Python native for JSON serialization)
    results = {
        'chunk_id': int(chunk_id),
        'total_rows': int(total_rows),
        'sample_size': int(len(metrics['phenolic_oh_counts'])),
        'phenolic_oh': {
            'mean': float(np.mean(metrics['phenolic_oh_counts'])),
            'median': float(np.median(metrics['phenolic_oh_counts'])),
            'std': float(np.std(metrics['phenolic_oh_counts'])),
            'max': int(np.max(metrics['phenolic_oh_counts'])),
            'min': int(np.min(metrics['phenolic_oh_counts']))
        },
        'mol_weight': {
            'mean': float(np.mean(metrics['mol_weights'])),
            'std': float(np.std(metrics['mol_weights']))
        },
        'complexity_score': {
            'mean': float(np.mean(metrics['complexity_scores'])),
            'median': float(np.median(metrics['complexity_scores'])),
            'max': float(np.max(metrics['complexity_scores']))
        },
        'predicted_products_per_mol': float(np.mean(metrics['complexity_scores']) * 3)  # Rough estimate
    }

    return results

def main():
    print("="*80)
    print("CHUNK COMPLEXITY ANALYSIS")
    print("="*80)

    # Analyze chunks 0-2 (and optionally more)
    chunks_to_analyze = [0, 1, 2, 3, 4]  # First 5 chunks

    all_results = []
    for chunk_id in chunks_to_analyze:
        print(f"\nAnalyzing Chunk {chunk_id}...")
        result = analyze_chunk_complexity(chunk_id)
        if result:
            all_results.append(result)

            print(f"  Total rows: {result['total_rows']:,}")
            print(f"  Phenolic OH (mean): {result['phenolic_oh']['mean']:.2f}")
            print(f"  Phenolic OH (max): {result['phenolic_oh']['max']}")
            print(f"  Complexity score (mean): {result['complexity_score']['mean']:.1f}")
            print(f"  Predicted products/mol: {result['predicted_products_per_mol']:.1f}")

    # Comparison
    print("\n" + "="*80)
    print("COMPARATIVE ANALYSIS")
    print("="*80)

    if len(all_results) >= 2:
        baseline = all_results[0]
        for result in all_results[1:]:
            complexity_ratio = result['complexity_score']['mean'] / baseline['complexity_score']['mean']
            print(f"\nChunk {result['chunk_id']} vs Chunk 0:")
            print(f"  Complexity ratio: {complexity_ratio:.2f}x")
            print(f"  Phenolic OH ratio: {result['phenolic_oh']['mean'] / baseline['phenolic_oh']['mean']:.2f}x")

            if complexity_ratio > 2.0:
                print(f"  WARNING: Chunk {result['chunk_id']} is {complexity_ratio:.1f}x more complex!")
                print(f"     Expected processing time: ~{44 * complexity_ratio:.0f} minutes")

    # Save results
    import json
    with open('chunk_complexity_analysis.json', 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"\n[OK] Results saved to: chunk_complexity_analysis.json")

if __name__ == '__main__':
    main()
