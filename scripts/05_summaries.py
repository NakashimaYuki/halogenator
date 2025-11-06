#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate Statistics and Summaries for Halogenated Flavonoid Library

Generates:
1. by_rule.csv - Product counts by rule family, halogen, and k-level
2. parents_coverage.csv - Per-parent analysis (k=1/k=2 counts, rules triggered)
3. macro_summary.csv - Macro substitution (CF3/CCl3) statistics
4. overall_stats.json - Global statistics and metrics
5. Text summary report

Analyzes:
- Product distribution by rule, halogen, k-level
- Parent coverage and productivity
- Macro vs step substitution ratios
- Deduplication effectiveness
- Average/median products per parent
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ==============================================================================
# Statistics Generation
# ==============================================================================

def generate_by_rule_csv(df: pd.DataFrame, output_path: str) -> pd.DataFrame:
    """
    Generate by_rule.csv: Product counts by rule_family × halogen × k
    """
    logger.info("Generating by_rule.csv...")

    if not all(col in df.columns for col in ["rule", "halogen", "k"]):
        logger.warning("Missing required columns for by_rule analysis")
        return pd.DataFrame()

    # Group by rule, halogen, k
    by_rule = df.groupby(["rule", "halogen", "k"]).size().reset_index(name="count")

    # Sort for readability
    by_rule = by_rule.sort_values(["rule", "halogen", "k"]).reset_index(drop=True)

    # Save
    by_rule.to_csv(output_path, index=False)
    logger.info(f"Saved by_rule.csv to {output_path}")

    return by_rule


def generate_parents_coverage_csv(df: pd.DataFrame, output_path: str) -> pd.DataFrame:
    """
    Generate parents_coverage.csv: Per-parent productivity analysis

    For each unique parent (k=0 or inferred from k=1 parent_inchikey):
    - Number of k=1 products
    - Number of k=2 products
    - Set of rules triggered
    - Set of halogens used
    """
    logger.info("Generating parents_coverage.csv...")

    if "k" not in df.columns:
        logger.warning("Missing 'k' column for parent coverage analysis")
        return pd.DataFrame()

    # Get k=1 products (direct children of parents)
    k1_df = df[df["k"] == 1].copy()

    # Get k=2 products
    k2_df = df[df["k"] == 2].copy()

    # Build parent coverage map
    parent_stats = defaultdict(lambda: {
        "k1_count": 0,
        "k2_count": 0,
        "rules": set(),
        "halogens": set()
    })

    # Process k=1 products
    for _, row in k1_df.iterrows():
        parent_ik = row.get("parent_inchikey", "")
        if not parent_ik:
            parent_ik = row.get("parent_name", "UNKNOWN")

        parent_stats[parent_ik]["k1_count"] += 1
        if "rule" in row:
            parent_stats[parent_ik]["rules"].add(row["rule"])
        if "halogen" in row:
            parent_stats[parent_ik]["halogens"].add(row["halogen"])

    # Process k=2 products (link back to k=1 parent)
    for _, row in k2_df.iterrows():
        k1_parent_ik = row.get("parent_inchikey", "")

        # Find the corresponding k=0 parent by looking up k1_parent_ik in k1_df
        if k1_parent_ik:
            k1_matches = k1_df[k1_df["inchikey"] == k1_parent_ik]
            if len(k1_matches) > 0:
                k0_parent_ik = k1_matches.iloc[0].get("parent_inchikey", "")
                if k0_parent_ik:
                    parent_stats[k0_parent_ik]["k2_count"] += 1
                    if "rule" in row:
                        parent_stats[k0_parent_ik]["rules"].add(row["rule"])
                    if "halogen" in row:
                        parent_stats[k0_parent_ik]["halogens"].add(row["halogen"])

    # Convert to DataFrame
    rows = []
    for parent_ik, stats in parent_stats.items():
        rows.append({
            "parent_inchikey": parent_ik,
            "k1_products": stats["k1_count"],
            "k2_products": stats["k2_count"],
            "total_products": stats["k1_count"] + stats["k2_count"],
            "rules_triggered": "|".join(sorted(stats["rules"])),
            "halogens_used": "|".join(sorted(stats["halogens"]))
        })

    coverage_df = pd.DataFrame(rows)

    # Sort by total products (descending)
    coverage_df = coverage_df.sort_values("total_products", ascending=False).reset_index(drop=True)

    # Save
    coverage_df.to_csv(output_path, index=False)
    logger.info(f"Saved parents_coverage.csv to {output_path}")

    return coverage_df


def generate_macro_summary_csv(df: pd.DataFrame, output_path: str) -> pd.DataFrame:
    """
    Generate macro_summary.csv: Macro substitution (CF3/CCl3) statistics
    """
    logger.info("Generating macro_summary.csv...")

    if "macro_label" not in df.columns:
        logger.warning("No 'macro_label' column found, skipping macro summary")
        return pd.DataFrame()

    # Filter for macro substitutions
    macro_df = df[df["macro_label"].notna()].copy()

    if len(macro_df) == 0:
        logger.warning("No macro substitutions found")
        return pd.DataFrame()

    # Group by macro_label, k, rule
    macro_summary = macro_df.groupby(["macro_label", "k", "rule"]).size().reset_index(name="count")

    # Sort
    macro_summary = macro_summary.sort_values(["macro_label", "k", "rule"]).reset_index(drop=True)

    # Add total row
    totals = macro_summary.groupby("macro_label")["count"].sum().reset_index(name="total")

    # Save
    macro_summary.to_csv(output_path, index=False)
    logger.info(f"Saved macro_summary.csv to {output_path}")

    logger.info("Macro substitution totals:")
    for _, row in totals.iterrows():
        logger.info(f"  {row['macro_label']}: {row['total']} products")

    return macro_summary


def generate_overall_stats(
    df: pd.DataFrame,
    by_rule_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    output_path: str
) -> dict:
    """
    Generate overall_stats.json: Global statistics and metrics
    """
    logger.info("Generating overall_stats.json...")

    stats = {
        "library": {
            "total_products": len(df),
            "unique_structures": len(df["inchikey"].unique()) if "inchikey" in df.columns else len(df)
        },

        "k_distribution": {},
        "rule_distribution": {},
        "halogen_distribution": {},

        "parent_productivity": {},

        "macro_substitutions": {},

        "quality_metrics": {}
    }

    # K-level distribution
    if "k" in df.columns:
        stats["k_distribution"] = df["k"].value_counts().to_dict()

    # Rule distribution
    if "rule" in df.columns:
        stats["rule_distribution"] = df["rule"].value_counts().to_dict()

    # Halogen distribution
    if "halogen" in df.columns:
        stats["halogen_distribution"] = df["halogen"].value_counts().to_dict()

    # Parent productivity (if coverage data available)
    if len(coverage_df) > 0:
        stats["parent_productivity"] = {
            "num_parents": len(coverage_df),
            "avg_k1_per_parent": float(coverage_df["k1_products"].mean()),
            "median_k1_per_parent": int(coverage_df["k1_products"].median()),
            "avg_k2_per_parent": float(coverage_df["k2_products"].mean()),
            "median_k2_per_parent": int(coverage_df["k2_products"].median()),
            "avg_total_per_parent": float(coverage_df["total_products"].mean()),
            "median_total_per_parent": int(coverage_df["total_products"].median()),
            "max_products_from_single_parent": int(coverage_df["total_products"].max())
        }

    # Macro substitutions
    if "macro_label" in df.columns:
        macro_df = df[df["macro_label"].notna()]
        stats["macro_substitutions"] = {
            "total_macro_products": len(macro_df),
            "macro_ratio": float(len(macro_df) / len(df)) if len(df) > 0 else 0,
            "by_type": macro_df["macro_label"].value_counts().to_dict()
        }

    # Quality metrics
    if "inchikey" in df.columns:
        stats["quality_metrics"]["structural_diversity"] = len(df["inchikey"].unique())

    # Save
    with open(output_path, "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"Saved overall_stats.json to {output_path}")

    return stats


def generate_text_summary(
    df: pd.DataFrame,
    stats: dict,
    output_path: str
) -> None:
    """
    Generate human-readable text summary report
    """
    logger.info("Generating summary report...")

    lines = []
    lines.append("=" * 80)
    lines.append("CNPD-ETCM HALOGENATED FLAVONOID LIBRARY - SUMMARY REPORT")
    lines.append("=" * 80)
    lines.append("")

    # Library overview
    lines.append("LIBRARY OVERVIEW")
    lines.append("-" * 80)
    lines.append(f"Total products: {stats['library']['total_products']:,}")
    lines.append(f"Unique structures: {stats['library']['unique_structures']:,}")
    lines.append("")

    # K-level distribution
    if stats.get("k_distribution"):
        lines.append("K-LEVEL DISTRIBUTION")
        lines.append("-" * 80)
        for k, count in sorted(stats["k_distribution"].items()):
            pct = 100 * count / stats['library']['total_products']
            lines.append(f"  k={k}: {count:,} products ({pct:.1f}%)")
        lines.append("")

    # Rule distribution
    if stats.get("rule_distribution"):
        lines.append("RULE FAMILY DISTRIBUTION")
        lines.append("-" * 80)
        for rule, count in sorted(stats["rule_distribution"].items(), key=lambda x: -x[1]):
            pct = 100 * count / stats['library']['total_products']
            lines.append(f"  {rule}: {count:,} products ({pct:.1f}%)")
        lines.append("")

    # Halogen distribution
    if stats.get("halogen_distribution"):
        lines.append("HALOGEN DISTRIBUTION")
        lines.append("-" * 80)
        for hal, count in sorted(stats["halogen_distribution"].items()):
            pct = 100 * count / stats['library']['total_products']
            lines.append(f"  {hal}: {count:,} products ({pct:.1f}%)")
        lines.append("")

    # Parent productivity
    if stats.get("parent_productivity"):
        pp = stats["parent_productivity"]
        lines.append("PARENT PRODUCTIVITY")
        lines.append("-" * 80)
        lines.append(f"Number of parent compounds: {pp['num_parents']:,}")
        lines.append("")
        lines.append("K=1 products per parent:")
        lines.append(f"  Average: {pp['avg_k1_per_parent']:.1f}")
        lines.append(f"  Median: {pp['median_k1_per_parent']}")
        lines.append("")
        lines.append("K=2 products per parent:")
        lines.append(f"  Average: {pp['avg_k2_per_parent']:.1f}")
        lines.append(f"  Median: {pp['median_k2_per_parent']}")
        lines.append("")
        lines.append("Total products per parent:")
        lines.append(f"  Average: {pp['avg_total_per_parent']:.1f}")
        lines.append(f"  Median: {pp['median_total_per_parent']}")
        lines.append(f"  Maximum: {pp['max_products_from_single_parent']}")
        lines.append("")

    # Macro substitutions
    if stats.get("macro_substitutions") and stats["macro_substitutions"].get("total_macro_products", 0) > 0:
        ms = stats["macro_substitutions"]
        lines.append("MACRO SUBSTITUTIONS (CF3/CCl3)")
        lines.append("-" * 80)
        lines.append(f"Total macro products: {ms['total_macro_products']:,}")
        lines.append(f"Macro ratio: {100*ms['macro_ratio']:.2f}%")
        lines.append("")
        lines.append("By type:")
        for mtype, count in sorted(ms["by_type"].items()):
            lines.append(f"  {mtype}: {count:,}")
        lines.append("")

    # Footer
    lines.append("=" * 80)
    lines.append("END OF REPORT")
    lines.append("=" * 80)

    # Write to file
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    logger.info(f"Saved summary report to {output_path}")

    # Also print to console
    print("\n" + "\n".join(lines))


# ==============================================================================
# Main Pipeline
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate statistics and summaries for halogenated flavonoid library",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--input",
        default="data/output/haloflav_k2/products_k2.parquet",
        help="Input parquet file (merged library)"
    )

    parser.add_argument(
        "-o", "--output-dir",
        default="data/output/haloflav_k2",
        help="Output directory for statistics files"
    )

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        logger.error("Please run merge script first (scripts/04_merge_and_qc.py)")
        sys.exit(1)

    # Load data
    logger.info("=" * 80)
    logger.info("FLAVONOID LIBRARY STATISTICS GENERATION")
    logger.info("=" * 80)
    logger.info("")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info("")

    logger.info("Loading library data...")
    df = pd.read_parquet(args.input)
    logger.info(f"Loaded {len(df)} products")
    logger.info("")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate statistics
    logger.info("Generating statistics files...")
    logger.info("")

    # 1. by_rule.csv
    by_rule_path = os.path.join(args.output_dir, "by_rule.csv")
    by_rule_df = generate_by_rule_csv(df, by_rule_path)

    # 2. parents_coverage.csv
    coverage_path = os.path.join(args.output_dir, "parents_coverage.csv")
    coverage_df = generate_parents_coverage_csv(df, coverage_path)

    # 3. macro_summary.csv
    macro_path = os.path.join(args.output_dir, "macro_summary.csv")
    macro_df = generate_macro_summary_csv(df, macro_path)

    # 4. overall_stats.json
    stats_path = os.path.join(args.output_dir, "overall_stats.json")
    stats = generate_overall_stats(df, by_rule_df, coverage_df, stats_path)

    # 5. Text summary report
    summary_path = os.path.join(args.output_dir, "SUMMARY_REPORT.txt")
    generate_text_summary(df, stats, summary_path)

    # Done
    logger.info("")
    logger.info("=" * 80)
    logger.info("STATISTICS GENERATION COMPLETE")
    logger.info("=" * 80)
    logger.info("")
    logger.info("Generated files:")
    logger.info(f"  - {by_rule_path}")
    logger.info(f"  - {coverage_path}")
    logger.info(f"  - {macro_path}")
    logger.info(f"  - {stats_path}")
    logger.info(f"  - {summary_path}")
    logger.info("")


if __name__ == "__main__":
    main()
