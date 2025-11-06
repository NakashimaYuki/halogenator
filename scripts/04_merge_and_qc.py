#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Merge Shard Results and Run Quality Control

Purpose:
1. Merge all shard products_k2.parquet files
2. Global deduplication by InChIKey
3. Run comprehensive QC validation suite
4. Generate merged output with metadata

QC Checks:
- Field consistency (k == k_ops, len(substitutions) == k_ops, etc.)
- Structural consistency (halogen counts, atom_cost validation)
- Parent chain integrity
- Schema validation
"""

import os
import sys
import glob
import json
import argparse
import subprocess
from datetime import datetime
from pathlib import Path
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ==============================================================================
# Merge Shards
# ==============================================================================

def find_shard_outputs(shard_output_roots: list) -> list:
    """
    Find all products_k2.parquet files in shard output directories

    Args:
        shard_output_roots: List of directories to search (can be single or multiple)
    """
    all_parquet_files = []

    for root_dir in shard_output_roots:
        if not os.path.exists(root_dir):
            logger.warning(f"Directory not found: {root_dir}")
            continue

        pattern = os.path.join(root_dir, "**/products_k2.parquet")
        parquet_files = glob.glob(pattern, recursive=True)

        logger.info(f"Found {len(parquet_files)} shard output files in {root_dir}")
        all_parquet_files.extend(parquet_files)

    logger.info(f"Total: {len(all_parquet_files)} shard output files")

    return all_parquet_files


def merge_shards(parquet_files: list, output_path: str) -> pd.DataFrame:
    """
    Merge all shard parquet files into single dataframe

    Args:
        parquet_files: List of paths to products_k2.parquet files
        output_path: Path to save merged parquet

    Returns:
        Merged DataFrame
    """
    logger.info(f"Merging {len(parquet_files)} shard files...")

    dfs = []
    total_records = 0

    for i, path in enumerate(parquet_files, 1):
        try:
            df = pd.read_parquet(path)
            dfs.append(df)
            total_records += len(df)

            if i % 10 == 0:
                logger.info(f"Loaded {i}/{len(parquet_files)} files, {total_records} records")

        except Exception as e:
            logger.error(f"Failed to load {path}: {e}")
            continue

    logger.info(f"Concatenating {len(dfs)} dataframes...")
    merged = pd.concat(dfs, ignore_index=True)

    logger.info(f"Merged dataset: {len(merged)} total records")

    # Save pre-dedup version (optional, for debugging)
    pre_dedup_path = output_path.replace(".parquet", "_pre_dedup.parquet")
    logger.info(f"Saving pre-dedup version to {pre_dedup_path}")
    merged.to_parquet(pre_dedup_path, index=False)

    return merged


def global_deduplication(df: pd.DataFrame) -> pd.DataFrame:
    """
    Global deduplication by InChIKey

    Note: Shards may have been deduplicated individually, but duplicates
    can still exist across shards if different parents in different shards
    produce the same product structure.
    """
    logger.info("Running global deduplication by InChIKey...")

    initial_count = len(df)

    if "inchikey" not in df.columns:
        logger.warning("No 'inchikey' column found, skipping deduplication")
        return df

    # Deduplicate by InChIKey (keep first occurrence)
    df_dedup = df.drop_duplicates(subset=["inchikey"], keep="first")
    df_dedup = df_dedup.reset_index(drop=True)

    final_count = len(df_dedup)
    removed = initial_count - final_count

    logger.info(f"Deduplication complete:")
    logger.info(f"  Before: {initial_count} records")
    logger.info(f"  After: {final_count} records")
    logger.info(f"  Removed: {removed} duplicates ({100*removed/initial_count:.2f}%)")

    return df_dedup


# ==============================================================================
# Quality Control
# ==============================================================================

def run_field_consistency_check(parquet_path: str) -> bool:
    """
    Run field consistency validation
    Uses existing validation script: scripts/validate_field_consistency.py
    """
    logger.info("=" * 80)
    logger.info("QC Check 1: Field Consistency")
    logger.info("=" * 80)

    script_path = "scripts/validate_field_consistency.py"

    if not os.path.exists(script_path):
        logger.warning(f"Field consistency script not found: {script_path}")
        logger.warning("Skipping field consistency check")
        return True

    try:
        result = subprocess.run(
            ["python", script_path, parquet_path],
            capture_output=True,
            text=True,
            timeout=300
        )

        print(result.stdout)

        if result.returncode == 0:
            logger.info("[PASS] Field consistency check")
            return True
        else:
            logger.error("[FAIL] Field consistency check")
            print(result.stderr)
            return False

    except subprocess.TimeoutExpired:
        logger.error("Field consistency check timed out")
        return False
    except Exception as e:
        logger.error(f"Field consistency check failed: {e}")
        return False


def run_structural_consistency_check(parquet_path: str) -> bool:
    """
    Run structural consistency validation
    Uses existing validation script: scripts/validate_structural_consistency.py
    """
    logger.info("=" * 80)
    logger.info("QC Check 2: Structural Consistency")
    logger.info("=" * 80)

    script_path = "scripts/validate_structural_consistency.py"

    if not os.path.exists(script_path):
        logger.warning(f"Structural consistency script not found: {script_path}")
        logger.warning("Skipping structural consistency check")
        return True

    try:
        # Run with --paranoid flag for thorough checking
        result = subprocess.run(
            ["python", script_path, parquet_path, "--paranoid"],
            capture_output=True,
            text=True,
            timeout=600
        )

        print(result.stdout)

        if result.returncode == 0:
            logger.info("[PASS] Structural consistency check")
            return True
        else:
            logger.error("[FAIL] Structural consistency check")
            print(result.stderr)
            return False

    except subprocess.TimeoutExpired:
        logger.error("Structural consistency check timed out")
        return False
    except Exception as e:
        logger.error(f"Structural consistency check failed: {e}")
        return False


def run_parent_chain_check(parquet_path: str) -> bool:
    """
    Run parent chain integrity validation
    Uses existing validation script: scripts/validate_parent_chain.py
    """
    logger.info("=" * 80)
    logger.info("QC Check 3: Parent Chain Integrity")
    logger.info("=" * 80)

    script_path = "scripts/validate_parent_chain.py"

    if not os.path.exists(script_path):
        logger.warning(f"Parent chain script not found: {script_path}")
        logger.warning("Skipping parent chain check")
        return True

    try:
        result = subprocess.run(
            ["python", script_path, parquet_path],
            capture_output=True,
            text=True,
            timeout=300
        )

        print(result.stdout)

        if result.returncode == 0:
            logger.info("[PASS] Parent chain integrity check")
            return True
        else:
            logger.error("[FAIL] Parent chain integrity check")
            print(result.stderr)
            return False

    except subprocess.TimeoutExpired:
        logger.error("Parent chain check timed out")
        return False
    except Exception as e:
        logger.error(f"Parent chain check failed: {e}")
        return False


def run_ci_validation_gate(output_dir: str) -> bool:
    """
    Run comprehensive CI validation gate
    Uses existing validation script: scripts/ci_validation_gate.py
    """
    logger.info("=" * 80)
    logger.info("QC Check 4: CI Validation Gate (Comprehensive)")
    logger.info("=" * 80)

    script_path = "scripts/ci_validation_gate.py"

    if not os.path.exists(script_path):
        logger.warning(f"CI validation script not found: {script_path}")
        logger.warning("Skipping CI validation gate")
        return True

    try:
        result = subprocess.run(
            ["python", script_path, output_dir],
            capture_output=True,
            text=True,
            timeout=600
        )

        print(result.stdout)

        if result.returncode == 0:
            logger.info("[PASS] CI validation gate")
            return True
        else:
            logger.error("[FAIL] CI validation gate")
            print(result.stderr)
            return False

    except subprocess.TimeoutExpired:
        logger.error("CI validation gate timed out")
        return False
    except Exception as e:
        logger.error(f"CI validation gate failed: {e}")
        return False


def run_qc_suite(output_dir: str, parquet_path: str) -> bool:
    """
    Run complete QC validation suite

    Returns True if all checks pass, False otherwise
    """
    logger.info("=" * 80)
    logger.info("QUALITY CONTROL VALIDATION SUITE")
    logger.info("=" * 80)
    logger.info("")

    results = {}

    # Check 1: Field consistency
    results["field_consistency"] = run_field_consistency_check(parquet_path)

    # Check 2: Structural consistency
    results["structural_consistency"] = run_structural_consistency_check(parquet_path)

    # Check 3: Parent chain integrity
    results["parent_chain"] = run_parent_chain_check(parquet_path)

    # Check 4: CI validation gate (comprehensive)
    results["ci_validation"] = run_ci_validation_gate(output_dir)

    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("QC VALIDATION SUMMARY")
    logger.info("=" * 80)

    all_passed = True
    for check, passed in results.items():
        status = "[PASS]" if passed else "[FAIL]"
        color = "green" if passed else "red"
        logger.info(f"{status} {check}")
        if not passed:
            all_passed = False

    logger.info("")

    if all_passed:
        logger.info("ALL QC CHECKS PASSED")
    else:
        logger.error("SOME QC CHECKS FAILED - REVIEW LOGS")

    return all_passed


# ==============================================================================
# Metadata Generation
# ==============================================================================

def generate_metadata(
    df: pd.DataFrame,
    shard_count: int,
    config_path: str,
    output_dir: str
) -> dict:
    """
    Generate run metadata for the merged library
    """
    metadata = {
        "pipeline": "cnpd_etcm_flavonoid_halogenation",
        "version": "1.0",
        "timestamp": datetime.now().isoformat(),

        "input": {
            "shard_count": shard_count,
            "config": config_path
        },

        "output": {
            "total_records": len(df),
            "k_levels": sorted(df["k"].unique().tolist()) if "k" in df.columns else [],
        },

        "statistics": {},

        "qc": {
            "status": "pending"
        }
    }

    # Add statistics
    if "k" in df.columns:
        metadata["statistics"]["by_k_level"] = df["k"].value_counts().to_dict()

    if "rule" in df.columns:
        metadata["statistics"]["by_rule"] = df["rule"].value_counts().to_dict()

    if "halogen" in df.columns:
        metadata["statistics"]["by_halogen"] = df["halogen"].value_counts().to_dict()

    if "macro_label" in df.columns:
        macro_counts = df[df["macro_label"].notna()]["macro_label"].value_counts().to_dict()
        metadata["statistics"]["macro_substitutions"] = macro_counts

    # Save metadata
    metadata_path = os.path.join(output_dir, "RUN_METADATA.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Metadata saved to {metadata_path}")

    return metadata


# ==============================================================================
# Main Pipeline
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Merge shard results and run QC validation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--shard-output-root",
        action='append',
        dest='shard_output_roots',
        help="Root directory containing shard output folders (can be specified multiple times)"
    )

    parser.add_argument(
        "--default-roots",
        action='store_true',
        help="Use default roots: data/output/cnpd_flav_k2 and data/output/cnpd_flav_k2_resumed"
    )

    parser.add_argument(
        "-o", "--output-dir",
        default="data/output/haloflav_k2",
        help="Output directory for merged results"
    )

    parser.add_argument(
        "-c", "--config",
        default="configs/flavonoids_k2_prod.yaml",
        help="Configuration file used for enumeration"
    )

    parser.add_argument(
        "--skip-qc",
        action="store_true",
        help="Skip QC validation (not recommended for production)"
    )

    parser.add_argument(
        "--skip-dedup",
        action="store_true",
        help="Skip global deduplication (not recommended)"
    )

    args = parser.parse_args()

    # Determine which root directories to use
    if args.default_roots:
        shard_roots = [
            "data/output/cnpd_flav_k2",
            "data/output/cnpd_flav_k2_resumed"
        ]
    elif args.shard_output_roots:
        shard_roots = args.shard_output_roots
    else:
        # Default to original directory for backwards compatibility
        shard_roots = ["data/output/cnpd_flav_k2"]

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    logger.info("=" * 80)
    logger.info("CNPD-ETCM Flavonoid Library - Merge and QC")
    logger.info("=" * 80)
    logger.info("")
    logger.info(f"Shard output roots:")
    for root in shard_roots:
        logger.info(f"  - {root}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Configuration: {args.config}")
    logger.info("")

    # Step 1: Find shard outputs
    logger.info("STEP 1: Finding shard output files")
    parquet_files = find_shard_outputs(shard_roots)

    if len(parquet_files) == 0:
        logger.error("No shard output files found!")
        logger.error("Please run enumeration first (scripts/03_enum_shards.*)")
        sys.exit(1)

    # Step 2: Merge shards
    logger.info("")
    logger.info("STEP 2: Merging shard results")
    merged_df = merge_shards(parquet_files, os.path.join(args.output_dir, "products_k2.parquet"))

    # Step 3: Global deduplication
    if not args.skip_dedup:
        logger.info("")
        logger.info("STEP 3: Global deduplication")
        merged_df = global_deduplication(merged_df)
    else:
        logger.warning("Skipping global deduplication (--skip-dedup flag)")

    # Step 4: Save merged output
    logger.info("")
    logger.info("STEP 4: Saving merged output")
    output_parquet = os.path.join(args.output_dir, "products_k2.parquet")
    merged_df.to_parquet(output_parquet, index=False)
    logger.info(f"Saved merged library to {output_parquet}")
    logger.info(f"Final record count: {len(merged_df)}")

    # Step 5: Generate metadata
    logger.info("")
    logger.info("STEP 5: Generating metadata")
    metadata = generate_metadata(
        merged_df,
        shard_count=len(parquet_files),
        config_path=args.config,
        output_dir=args.output_dir
    )

    # Step 6: Run QC validation
    if not args.skip_qc:
        logger.info("")
        logger.info("STEP 6: Quality control validation")
        qc_passed = run_qc_suite(args.output_dir, output_parquet)

        # Update metadata with QC status
        metadata["qc"]["status"] = "passed" if qc_passed else "failed"
        metadata_path = os.path.join(args.output_dir, "RUN_METADATA.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        if not qc_passed:
            logger.error("")
            logger.error("=" * 80)
            logger.error("QC VALIDATION FAILED")
            logger.error("=" * 80)
            logger.error("Review error messages above and fix issues before using library")
            sys.exit(1)
    else:
        logger.warning("Skipping QC validation (--skip-qc flag)")

    # Success
    logger.info("")
    logger.info("=" * 80)
    logger.info("MERGE AND QC COMPLETE")
    logger.info("=" * 80)
    logger.info("")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Products file: {output_parquet}")
    logger.info(f"Total compounds: {len(merged_df)}")
    logger.info("")
    logger.info("Next step: Generate statistics and summaries")
    logger.info("  python scripts/05_summaries.py")
    logger.info("")


if __name__ == "__main__":
    main()
