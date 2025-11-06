#!/usr/bin/env python
"""
Field verification script for scenario parquet outputs.
Checks presence of sub_rule, detection, rule_family fields and validates distributions.

Usage:
  python scripts/verify_scenario_fields.py <parquet_path>
"""

import pandas as pd
import sys
from pathlib import Path


def verify_fields(parquet_path: str) -> dict:
    """Verify critical fields in parquet output."""

    if not Path(parquet_path).exists():
        return {"error": f"File not found: {parquet_path}"}

    df = pd.read_parquet(parquet_path)

    results = {
        "path": parquet_path,
        "total_products": len(df),
        "has_sub_rule": "sub_rule" in df.columns,
        "has_detection": "detection" in df.columns,
        "has_rule_family": "rule_family" in df.columns,
        "has_macro_label": "macro_label" in df.columns,
        "k_equals_k_ops": (df["k"] == df["k_ops"]).all() if "k" in df.columns and "k_ops" in df.columns else None,
    }

    # Family distribution
    if "rule_family" in df.columns:
        family_counts = df.groupby("rule_family")["smiles"].count().to_dict()
        results["family_counts"] = family_counts

    # R2 family analysis
    if "rule_family" in df.columns:
        r2_df = df[df["rule_family"] == "R2"]
        if len(r2_df) > 0:
            results["r2_count"] = len(r2_df)
            if "sub_rule" in df.columns:
                results["r2_sub_rule_dist"] = r2_df["sub_rule"].value_counts(dropna=False).to_dict()
            if "detection" in df.columns:
                results["r2_detection_dist"] = r2_df["detection"].value_counts(dropna=False).to_dict()
        else:
            results["r2_count"] = 0

    # R6 family analysis
    if "rule_family" in df.columns:
        r6_df = df[df["rule_family"] == "R6"]
        if len(r6_df) > 0:
            results["r6_count"] = len(r6_df)
        else:
            results["r6_count"] = 0

    # Macro label analysis
    if "macro_label" in df.columns:
        macro_df = df[df["macro_label"].notna()]
        results["macro_count"] = len(macro_df)
        if len(macro_df) > 0:
            results["macro_labels"] = macro_df["macro_label"].value_counts().to_dict()
            # Sample macro products k_ops/k_atoms
            if "k_ops" in df.columns and "k_atoms" in df.columns:
                sample = macro_df[["macro_label", "k_ops", "k_atoms"]].head(5).to_dict("records")
                results["macro_sample"] = sample
    else:
        results["macro_count"] = 0

    return results


def print_results(results: dict):
    """Pretty print verification results."""

    if "error" in results:
        print(f"ERROR: {results['error']}")
        return

    print(f"\n{'='*60}")
    print(f"Field Verification Report")
    print(f"{'='*60}")
    print(f"Path: {results['path']}")
    print(f"Total products: {results['total_products']}")
    print(f"\nCore Fields:")
    print(f"  [OK] Has sub_rule: {results['has_sub_rule']}")
    print(f"  [OK] Has detection: {results['has_detection']}")
    print(f"  [OK] Has rule_family: {results['has_rule_family']}")
    print(f"  [OK] Has macro_label: {results['has_macro_label']}")
    print(f"  [OK] k == k_ops: {results['k_equals_k_ops']}")

    if results.get("family_counts"):
        print(f"\nRule Family Distribution:")
        for family, count in sorted(results["family_counts"].items(), key=lambda x: -x[1]):
            print(f"  {family}: {count}")

    if results.get("r2_count", 0) > 0:
        print(f"\nR2 Family Analysis ({results['r2_count']} products):")
        if results.get("r2_sub_rule_dist"):
            print(f"  sub_rule: {results['r2_sub_rule_dist']}")
        if results.get("r2_detection_dist"):
            print(f"  detection: {results['r2_detection_dist']}")

    if results.get("r6_count", 0) > 0:
        print(f"\nR6 Family: {results['r6_count']} products")

    if results.get("macro_count", 0) > 0:
        print(f"\nMacro Products ({results['macro_count']} total):")
        if results.get("macro_labels"):
            print(f"  Labels: {results['macro_labels']}")
        if results.get("macro_sample"):
            print(f"  Sample:")
            for row in results["macro_sample"]:
                print(f"    {row}")

    print(f"{'='*60}\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    parquet_path = sys.argv[1]
    results = verify_fields(parquet_path)
    print_results(results)
