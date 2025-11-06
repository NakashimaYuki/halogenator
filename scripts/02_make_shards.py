#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate .smi Shard Files from Flavonoid Database

Takes flavonoids_final.parquet and splits into manageable .smi files
for parallel batch enumeration.

Shard size: 2000-5000 compounds (configurable)
Output format: <SMILES> <PARENT_NAME>
"""

import os
import sys
import math
import argparse
import pandas as pd
import logging
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def sanitize_name(name: str, max_len: int = 60) -> str:
    """
    Create safe filename-compatible parent name
    - Remove/replace problematic characters
    - Limit length
    - Ensure uniqueness will be handled by appending index if needed
    """
    if not name or not isinstance(name, str):
        return "UNKNOWN"

    # Replace spaces and problematic characters
    safe_name = name.replace(" ", "_")
    safe_name = safe_name.replace("/", "_")
    safe_name = safe_name.replace("\\", "_")
    safe_name = safe_name.replace(":", "_")
    safe_name = safe_name.replace("*", "_")
    safe_name = safe_name.replace("?", "_")
    safe_name = safe_name.replace('"', "_")
    safe_name = safe_name.replace("<", "_")
    safe_name = safe_name.replace(">", "_")
    safe_name = safe_name.replace("|", "_")
    safe_name = safe_name.replace("\t", "_")
    safe_name = safe_name.replace("\n", "_")
    safe_name = safe_name.replace("\r", "_")

    # Remove any non-ASCII characters (optional - keep if you want Chinese names)
    # For safety in cross-platform file handling, converting to ASCII
    try:
        safe_name = safe_name.encode('ascii', errors='ignore').decode('ascii')
    except:
        safe_name = "UNKNOWN"

    # Truncate to max length
    if len(safe_name) > max_len:
        safe_name = safe_name[:max_len]

    # Fallback if empty after sanitization
    if not safe_name or safe_name.strip() == "":
        safe_name = "UNKNOWN"

    return safe_name


def generate_parent_name(row: pd.Series, idx: int) -> str:
    """
    Generate unique parent name with priority:
    1. Pubchem-ID (if available)
    2. ID field
    3. Name field (sanitized)
    4. InChIKey
    5. Index-based fallback

    Append index to ensure uniqueness
    """
    # Try different fields in priority order
    for field in ["Pubchem-ID", "ID", "Name"]:
        if field in row and pd.notna(row[field]):
            base_name = str(row[field])
            safe_name = sanitize_name(base_name)
            if safe_name != "UNKNOWN":
                return f"{safe_name}_{idx:06d}"

    # Fallback to InChIKey
    if "InChIKey" in row and pd.notna(row["InChIKey"]):
        inchikey = str(row["InChIKey"])
        return f"{inchikey[:14]}_{idx:06d}"

    # Ultimate fallback
    return f"FLAV_{idx:06d}"


def make_shards(
    input_parquet: str,
    output_dir: str,
    shard_size: int = 3000,
    prefix: str = "flav_shard"
) -> None:
    """
    Split flavonoid database into .smi shard files

    Args:
        input_parquet: Path to flavonoids_final.parquet
        output_dir: Directory to write shard files
        shard_size: Number of compounds per shard
        prefix: Prefix for shard filenames
    """

    # Load data
    logger.info(f"Loading flavonoid database from {input_parquet}")
    if not os.path.exists(input_parquet):
        raise FileNotFoundError(f"Input file not found: {input_parquet}")

    df = pd.read_parquet(input_parquet)
    logger.info(f"Loaded {len(df)} flavonoid compounds")

    # Verify required columns
    if "Smiles_clean" not in df.columns:
        raise ValueError("Missing 'Smiles_clean' column in input data")

    # Reset index to ensure clean sequential indices
    df = df.reset_index(drop=True)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Calculate number of shards
    num_shards = math.ceil(len(df) / shard_size)
    logger.info(f"Generating {num_shards} shards with max {shard_size} compounds each")

    # Generate shards
    shard_manifest = []

    for shard_idx in range(num_shards):
        # Calculate slice boundaries
        start_idx = shard_idx * shard_size
        end_idx = min((shard_idx + 1) * shard_size, len(df))
        shard_df = df.iloc[start_idx:end_idx]

        # Shard filename (1-indexed for human readability)
        shard_filename = f"{prefix}_{shard_idx + 1:04d}.smi"
        shard_path = os.path.join(output_dir, shard_filename)

        # Write .smi file
        logger.info(f"Writing shard {shard_idx + 1}/{num_shards}: {shard_filename} ({len(shard_df)} compounds)")

        with open(shard_path, "w", encoding="utf-8") as f:
            for local_idx, (global_idx, row) in enumerate(shard_df.iterrows()):
                smiles = row["Smiles_clean"]
                parent_name = generate_parent_name(row, global_idx)

                # Write SMILES and name (tab-separated)
                f.write(f"{smiles}\t{parent_name}\n")

        # Track manifest
        shard_manifest.append({
            "shard_id": shard_idx + 1,
            "filename": shard_filename,
            "path": shard_path,
            "num_compounds": len(shard_df),
            "start_idx": start_idx,
            "end_idx": end_idx
        })

    # Write manifest
    manifest_path = os.path.join(output_dir, "shard_manifest.csv")
    manifest_df = pd.DataFrame(shard_manifest)
    manifest_df.to_csv(manifest_path, index=False)
    logger.info(f"Shard manifest written to {manifest_path}")

    # Summary
    logger.info("=" * 80)
    logger.info("SHARD GENERATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total compounds: {len(df)}")
    logger.info(f"Number of shards: {num_shards}")
    logger.info(f"Shard size: {shard_size} (target)")
    logger.info(f"Average compounds per shard: {len(df) / num_shards:.1f}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Manifest: {manifest_path}")

    # Shard size distribution
    shard_sizes = manifest_df["num_compounds"].values
    logger.info(f"\nShard size statistics:")
    logger.info(f"  Min: {shard_sizes.min()}")
    logger.info(f"  Max: {shard_sizes.max()}")
    logger.info(f"  Mean: {shard_sizes.mean():.1f}")
    logger.info(f"  Median: {int(pd.Series(shard_sizes).median())}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate .smi shard files from flavonoid database",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        default="data/work/flavonoids_final.parquet",
        help="Input parquet file (flavonoids_final.parquet)"
    )

    parser.add_argument(
        "-o", "--output",
        default="data/work/shards",
        help="Output directory for shard files"
    )

    parser.add_argument(
        "-s", "--shard-size",
        type=int,
        default=3000,
        help="Number of compounds per shard (2000-5000 recommended)"
    )

    parser.add_argument(
        "-p", "--prefix",
        default="flav_shard",
        help="Prefix for shard filenames"
    )

    args = parser.parse_args()

    # Validate shard size
    if args.shard_size < 100:
        logger.warning("Shard size < 100 may create too many small files")
    if args.shard_size > 10000:
        logger.warning("Shard size > 10000 may cause memory issues during enumeration")

    # Run
    make_shards(
        input_parquet=args.input,
        output_dir=args.output,
        shard_size=args.shard_size,
        prefix=args.prefix
    )


if __name__ == "__main__":
    main()
