#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Re-shard a shard file into multiple sub-shards based on complexity analysis
"""

import os
import sys
import argparse
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def reshard_by_complexity(shard_file, complexity_csv, output_dir, n_shards=6):
    """
    Split a shard file into multiple sub-shards based on complexity scores

    Args:
        shard_file: Original .smi shard file
        complexity_csv: CSV file with complexity analysis (from analyze_shard_complexity.py)
        output_dir: Output directory for sub-shard files
        n_shards: Number of sub-shards to create
    """

    logger.info(f"Re-sharding: {shard_file}")
    logger.info(f"Using complexity analysis: {complexity_csv}")
    logger.info(f"Target: {n_shards} sub-shards")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load complexity analysis
    if not os.path.exists(complexity_csv):
        logger.error(f"Complexity CSV not found: {complexity_csv}")
        return 1

    df_complexity = pd.read_csv(complexity_csv)
    logger.info(f"Loaded complexity data for {len(df_complexity)} molecules")

    # Check required columns
    required_cols = ['index', 'name', 'smiles', 'complexity_score']
    missing_cols = [col for col in required_cols if col not in df_complexity.columns]
    if missing_cols:
        logger.error(f"Missing required columns in CSV: {missing_cols}")
        return 1

    # Sort by complexity score
    df_sorted = df_complexity.sort_values('complexity_score').reset_index(drop=True)
    logger.info(f"Sorted molecules by complexity score")

    # Calculate molecules per shard
    molecules_per_shard = len(df_sorted) // n_shards
    remainder = len(df_sorted) % n_shards

    logger.info(f"Molecules per sub-shard: {molecules_per_shard} (with {remainder} remainder)")

    # Determine base name for output files
    shard_base = os.path.basename(shard_file)
    shard_name = os.path.splitext(shard_base)[0]  # Remove .smi extension

    # Split into sub-shards
    for shard_id in range(n_shards):
        # Calculate start and end indices for this sub-shard
        start_idx = shard_id * molecules_per_shard

        # Distribute remainder molecules to first few shards
        if shard_id < remainder:
            start_idx += shard_id
            end_idx = start_idx + molecules_per_shard + 1
        else:
            start_idx += remainder
            end_idx = start_idx + molecules_per_shard

        # Extract molecules for this sub-shard
        shard_df = df_sorted.iloc[start_idx:end_idx]

        # Generate output filename (e.g., flav_shard_0002a.smi)
        # Use letters a-f for sub-shards
        sub_shard_suffix = chr(ord('a') + shard_id)
        output_file = os.path.join(output_dir, f"{shard_name}{sub_shard_suffix}.smi")

        # Write sub-shard file
        with open(output_file, 'w', encoding='utf-8') as f:
            for _, row in shard_df.iterrows():
                f.write(f"{row['smiles']}\t{row['name']}\n")

        logger.info(f"  Sub-shard {shard_id} ({sub_shard_suffix}): {len(shard_df)} molecules → {output_file}")
        logger.info(f"    Complexity: {shard_df['complexity_score'].min():.3f} - {shard_df['complexity_score'].max():.3f}")
        logger.info(f"    SMILES length: {shard_df['smiles_len'].mean():.1f} ± {shard_df['smiles_len'].std():.1f}")
        logger.info(f"    O atoms: {shard_df['o_count'].mean():.1f} ± {shard_df['o_count'].std():.1f}")

    logger.info(f"Re-sharding complete! Created {n_shards} sub-shards in {output_dir}")

    return 0


def main():
    parser = argparse.ArgumentParser(description='Re-shard a shard file by complexity')
    parser.add_argument('shard_file', help='Original .smi shard file')
    parser.add_argument('complexity_csv', help='CSV file with complexity analysis')
    parser.add_argument('--output-dir', default='data/work/shards_resharded',
                        help='Output directory for sub-shard files')
    parser.add_argument('--n-shards', type=int, default=6,
                        help='Number of sub-shards to create (default: 6)')

    args = parser.parse_args()

    return reshard_by_complexity(
        args.shard_file,
        args.complexity_csv,
        args.output_dir,
        args.n_shards
    )


if __name__ == '__main__':
    sys.exit(main())
