#!/usr/bin/env python3
"""
Cross-platform batch statistics script.
Unified replacement for inline Python code in batch processing scripts.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any


def check_pandas_availability() -> bool:
    """Check if pandas and pyarrow are available for advanced statistics."""
    try:
        import pandas
        import pyarrow
        return True
    except ImportError:
        return False


def print_basic_file_stats(file_path: str, label: str) -> None:
    """Print basic file information without pandas."""
    # Try to find the best available products file
    actual_path = find_products_file(file_path)
    if actual_path != file_path:
        print(f"[INFO] Using products file: {actual_path}")

    path = Path(actual_path)
    if path.exists():
        size_bytes = path.stat().st_size
        size_mb = size_bytes / (1024 * 1024)
        print(f"[INFO] {label} file exists: {actual_path}")
        print(f"[INFO] File size: {size_bytes:,} bytes ({size_mb:.2f} MB)")
    else:
        print(f"[WARN] {label} file not found: {actual_path}")


def find_products_file(base_path: str) -> Optional[str]:
    """Find products file, preferring standardized name but falling back to k2/k3 naming."""
    from pathlib import Path

    base_dir = Path(base_path).parent

    # First preference: standardized name
    std_path = base_dir / "products.parquet"
    if std_path.exists():
        return str(std_path)

    # Fall back to the provided path (likely with k2/k3 naming)
    if Path(base_path).exists():
        return base_path

    # Try to infer k2/k3 variants if original doesn't exist
    for variant in ["products_k2.parquet", "products_k3.parquet"]:
        variant_path = base_dir / variant
        if variant_path.exists():
            return str(variant_path)

    return base_path  # Return original path even if not found for error reporting


def print_pandas_stats(products_path: str, label: str) -> None:
    """Print detailed statistics using pandas."""
    try:
        import pandas as pd

        # Try to find the best available products file
        actual_path = find_products_file(products_path)
        if actual_path != products_path:
            print(f"[INFO] Using products file: {actual_path}")

        df = pd.read_parquet(actual_path)
        print(f"[INFO] {label} Products generated: {len(df):,}")

        # Parent molecules
        if "parent_key" in df.columns:
            unique_parents = df["parent_key"].nunique()
            print(f"[INFO] Unique parent molecules: {unique_parents:,}")
        else:
            print(f"[INFO] Unique parent molecules: N/A (no parent_key column)")

        # Rule distribution
        if "rule" in df.columns:
            rules_dist = dict(df["rule"].value_counts())
            print(f"[INFO] Rules distribution: {rules_dist}")
        else:
            print(f"[INFO] Rules distribution: N/A (no rule column)")

        # Halogen distribution
        if "halogen" in df.columns:
            halogens_dist = dict(df["halogen"].value_counts())
            print(f"[INFO] Halogens distribution: {halogens_dist}")
        else:
            print(f"[INFO] Halogens distribution: N/A (no halogen column)")

    except Exception as e:
        print(f"[ERROR] Failed to analyze {label} results: {e}")
        # Fallback to basic stats
        print_basic_file_stats(products_path, label)


def print_qa_summary(qa_path: str) -> None:
    """Print QA warnings summary from JSON file."""
    if not os.path.exists(qa_path):
        print(f"[WARN] QA summary not found: {qa_path}")
        return

    try:
        with open(qa_path, 'r') as f:
            qa = json.load(f)

        metadata = qa.get('metadata', {})
        if 'warnings_count' in metadata:
            warnings_count = metadata['warnings_count']
            warnings_returned = metadata.get('warnings_returned', 0)
            print(f"[INFO] QA warnings: {warnings_count} (returned: {warnings_returned})")
        else:
            print(f"[INFO] QA metadata available, no warnings count found")

    except Exception as e:
        print(f"[ERROR] Failed to read QA summary: {e}")


def generate_audit_sample(products_path: str, output_dir: str, sample_size: int = 100) -> None:
    """Generate audit sample CSV file."""
    try:
        import pandas as pd
        import random

        # Try to find the best available products file
        actual_path = find_products_file(products_path)
        if actual_path != products_path:
            print(f"[INFO] Using products file for audit: {actual_path}")

        df = pd.read_parquet(actual_path)

        # Find parent key column with fallback priority
        key_col = None
        for candidate in ("parent_key", "unified_parent_key", "parent_inchikey", "parent_name"):
            if candidate in df.columns:
                key_col = candidate
                print(f"[INFO] Using parent identifier column: {key_col}")
                break

        if key_col is None:
            print(f"[WARN] Cannot generate audit sample: no parent identifier column found")
            return

        parents = list(df[key_col].dropna().unique())
        audit_parents = random.sample(parents, min(sample_size, len(parents)))

        audit_products = []
        for p in audit_parents:
            parent_products = df[df[key_col] == p]
            if len(parent_products) > 0:
                audit_products.append(parent_products.iloc[0])

        if audit_products:
            audit_df = pd.DataFrame(audit_products)
            audit_path = os.path.join(output_dir, "audit", f"sample_{len(audit_products)}.csv")
            os.makedirs(os.path.dirname(audit_path), exist_ok=True)
            audit_df.to_csv(audit_path, index=False)
            print(f"[INFO] Generated audit sample: {len(audit_products)} products from {len(audit_parents)} parents")
            print(f"[INFO] Audit sample saved to: {audit_path}")
        else:
            print(f"[WARN] No audit products generated")

    except ImportError:
        print(f"[WARN] Pandas not available - skipping audit sample generation")
    except Exception as e:
        print(f"[ERROR] Failed to generate audit sample: {e}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Cross-platform batch statistics")
    parser.add_argument("--products", help="Path to products parquet file")
    parser.add_argument("--qa", help="Path to QA JSON file")
    parser.add_argument("--label", default="k=?", help="Label for output (e.g., 'k=2', 'k=3')")
    parser.add_argument("--audit", help="Generate audit sample to this directory")
    parser.add_argument("--audit-size", type=int, default=100, help="Audit sample size")
    parser.add_argument("--check-deps", action="store_true", help="Check pandas/pyarrow availability")

    args = parser.parse_args()

    # Check dependencies if requested
    if args.check_deps:
        if check_pandas_availability():
            print("Pandas and PyArrow available")
            sys.exit(0)
        else:
            print("Pandas and/or PyArrow not available")
            sys.exit(1)

    # Print products statistics
    if args.products:
        if check_pandas_availability():
            print_pandas_stats(args.products, args.label)
        else:
            print("[WARN] Pandas not available - using basic file-based statistics")
            print_basic_file_stats(args.products, args.label)

    # Print QA summary
    if args.qa:
        print_qa_summary(args.qa)

    # Generate audit sample
    if args.audit and args.products:
        generate_audit_sample(args.products, args.audit, args.audit_size)


if __name__ == "__main__":
    main()