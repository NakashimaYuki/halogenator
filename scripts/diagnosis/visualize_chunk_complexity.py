#!/usr/bin/env python3
"""
Visualize chunk complexity distribution and create analysis report.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def load_complexity_data():
    """Load complexity analysis results."""
    with open('chunk_complexity_analysis.json', 'r') as f:
        return json.load(f)

def create_visualizations(data):
    """Create comprehensive visualization plots."""

    chunk_ids = [r['chunk_id'] for r in data]
    complexities = [r['complexity_score']['mean'] for r in data]
    phenolic_oh = [r['phenolic_oh']['mean'] for r in data]
    predicted_products = [r['predicted_products_per_mol'] for r in data]
    total_rows = [r['total_rows'] for r in data]

    # Calculate expected processing time (baseline = 44 min for chunk 0)
    baseline_complexity = complexities[0]
    baseline_time = 44  # minutes
    expected_times = [baseline_time * (c / baseline_complexity) for c in complexities]

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Chunk Complexity Analysis - polyphenol-2X Dataset', fontsize=16, fontweight='bold')

    # Plot 1: Complexity Score Distribution
    ax1 = axes[0, 0]
    bars1 = ax1.bar(chunk_ids, complexities, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axhline(y=baseline_complexity, color='green', linestyle='--', linewidth=2, label='Baseline (Chunk 0)')
    ax1.axhline(y=baseline_complexity * 2, color='orange', linestyle='--', linewidth=1, label='2x Baseline')
    ax1.axhline(y=baseline_complexity * 10, color='red', linestyle='--', linewidth=1, label='10x Baseline')

    # Highlight extreme chunks
    for i, (chunk_id, complexity) in enumerate(zip(chunk_ids, complexities)):
        if complexity > baseline_complexity * 10:
            bars1[i].set_color('red')
            bars1[i].set_alpha(1.0)
            ax1.text(chunk_id, complexity + 20, f'{complexity:.0f}\n({complexity/baseline_complexity:.1f}x)',
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax1.set_xlabel('Chunk ID', fontsize=12)
    ax1.set_ylabel('Complexity Score', fontsize=12)
    ax1.set_title('Complexity Score by Chunk', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # Plot 2: Expected Processing Time
    ax2 = axes[0, 1]
    colors = ['green' if t <= 120 else 'orange' if t <= 240 else 'red' for t in expected_times]
    bars2 = ax2.bar(chunk_ids, [t/60 for t in expected_times], color=colors, alpha=0.7, edgecolor='black')
    ax2.axhline(y=4, color='red', linestyle='--', linewidth=2, label='4h Timeout Limit')

    for i, (chunk_id, time) in enumerate(zip(chunk_ids, expected_times)):
        if time > 240:  # > 4 hours
            ax2.text(chunk_id, time/60 + 0.5, f'{time/60:.1f}h',
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax2.set_xlabel('Chunk ID', fontsize=12)
    ax2.set_ylabel('Expected Time (hours)', fontsize=12)
    ax2.set_title('Estimated Processing Time by Chunk', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    # Plot 3: Phenolic OH Distribution
    ax3 = axes[1, 0]
    ax3.bar(chunk_ids, phenolic_oh, color='forestgreen', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Chunk ID', fontsize=12)
    ax3.set_ylabel('Avg Phenolic OH Count', fontsize=12)
    ax3.set_title('Average Phenolic OH Groups by Chunk', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # Plot 4: Products per Molecule
    ax4 = axes[1, 1]
    ax4.bar(chunk_ids, predicted_products, color='purple', alpha=0.7, edgecolor='black')
    ax4.set_xlabel('Chunk ID', fontsize=12)
    ax4.set_ylabel('Predicted Products/Molecule', fontsize=12)
    ax4.set_title('Product Generation Potential by Chunk', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    # Save figure
    output_path = Path('chunk_complexity_visualization.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f'[OK] Visualization saved to: {output_path}')

    return expected_times

def generate_summary_report(data, expected_times):
    """Generate detailed summary report."""

    print("\n" + "="*80)
    print("FULL CHUNK COMPLEXITY ANALYSIS SUMMARY")
    print("="*80)

    baseline_complexity = data[0]['complexity_score']['mean']

    # Calculate statistics
    total_expected_time = sum(expected_times)
    chunks_over_4h = sum(1 for t in expected_times if t > 240)
    chunks_over_8h = sum(1 for t in expected_times if t > 480)
    max_time_chunk = max(range(len(expected_times)), key=lambda i: expected_times[i])

    print(f"\nDataset: polyphenol-2X")
    print(f"Total chunks: {len(data)}")
    print(f"Total rows: {sum(r['total_rows'] for r in data):,}")

    print(f"\nComplexity Distribution:")
    print(f"  Min complexity: {min(r['complexity_score']['mean'] for r in data):.1f} (Chunk {min(range(len(data)), key=lambda i: data[i]['complexity_score']['mean'])})")
    print(f"  Max complexity: {max(r['complexity_score']['mean'] for r in data):.1f} (Chunk {max(range(len(data)), key=lambda i: data[i]['complexity_score']['mean'])})")
    print(f"  Mean complexity: {np.mean([r['complexity_score']['mean'] for r in data]):.1f}")
    print(f"  Std deviation: {np.std([r['complexity_score']['mean'] for r in data]):.1f}")

    print(f"\nProcessing Time Estimates:")
    print(f"  Total estimated time: {total_expected_time/60:.1f} hours ({total_expected_time/60/24:.1f} days)")
    print(f"  Chunks exceeding 4h timeout: {chunks_over_4h}/{len(data)}")
    print(f"  Chunks exceeding 8h: {chunks_over_8h}/{len(data)}")
    print(f"  Longest chunk: Chunk {max_time_chunk} ({expected_times[max_time_chunk]/60:.1f} hours)")

    print(f"\nTop 5 Most Complex Chunks:")
    sorted_chunks = sorted(enumerate(data), key=lambda x: x[1]['complexity_score']['mean'], reverse=True)
    for i, (chunk_idx, chunk_data) in enumerate(sorted_chunks[:5], 1):
        complexity = chunk_data['complexity_score']['mean']
        ratio = complexity / baseline_complexity
        time_h = expected_times[chunk_idx] / 60
        print(f"  {i}. Chunk {chunk_idx}: complexity {complexity:.1f} ({ratio:.1f}x) -> {time_h:.1f}h")

    print(f"\nRecommendations:")
    if chunks_over_4h > 0:
        print(f"  [CRITICAL] {chunks_over_4h} chunks will exceed 4h timeout")
        print(f"  [ACTION] Implement adaptive chunking immediately")
    if chunks_over_8h > 0:
        print(f"  [WARNING] {chunks_over_8h} chunks would need >8h even with doubled timeout")

    print("\n" + "="*80)

def main():
    print("Loading complexity analysis data...")
    data = load_complexity_data()

    print("Creating visualizations...")
    expected_times = create_visualizations(data)

    generate_summary_report(data, expected_times)

if __name__ == '__main__':
    main()
