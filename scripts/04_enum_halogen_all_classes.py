#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Orchestrator for full-scale per-class halogenation library generation.

Usage:
    python 04_enum_halogen_all_classes.py --classes terpenoid polyphenol --k-values 1 2
    python 04_enum_halogen_all_classes.py --all --k-values 1 2

This script:
1. Loops over specified NP classes
2. For each class and k value, runs halogenator enum-parquet
3. Generates summary statistics
4. Outputs consolidated report
"""

import argparse
import json
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
LOG = logging.getLogger(__name__)

# All available NP classes
ALL_CLASSES = [
    "terpenoid",
    "alkaloid",
    "polyphenol",
    "glycoside",
    "lipid",
    "aa_peptide",
    "other"
]


def run_enum_for_class(
    np_class: str,
    k: int,
    batch_size: int = 5000,
    rdkit_threads: int = 8,
    max_parents: int = 0
) -> Dict[str, Any]:
    """
    Run full enumeration for a single class and k value.

    Args:
        np_class: NP class name
        k: Halogenation depth
        batch_size: Batch size
        rdkit_threads: Number of RDKit threads
        max_parents: Max parents to process (0 = all)

    Returns:
        Statistics dict
    """
    LOG.info(f"=" * 70)
    LOG.info(f"ENUMERATING: {np_class} k={k}")
    LOG.info(f"=" * 70)

    # Setup paths
    data_dir = Path("E:/Projects/halogenator/data/output/nplike")
    input_parquet = data_dir / np_class / "base.parquet"
    output_dir = data_dir / f"{np_class}-{k}X"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_parquet.exists():
        LOG.error(f"Input not found: {input_parquet}")
        return {"status": "error", "message": "Input not found"}

    # Build command
    cmd = [
        "halogenator",
        "enum-parquet",
        "--input-parquet", str(input_parquet),
        "--outdir", str(output_dir),
        "--k", str(k),
        "--np-class", np_class,
        "--batch-size", str(batch_size),
        "--rdkit-threads", str(rdkit_threads)
    ]

    if max_parents > 0:
        cmd.extend(["--max-parents", str(max_parents)])

    LOG.info(f"Command: {' '.join(cmd)}")

    # Run enumeration
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        LOG.error(f"Enumeration failed for {np_class} k={k}")
        LOG.error(f"stderr: {result.stderr}")
        return {
            "status": "error",
            "class": np_class,
            "k": k,
            "message": result.stderr[:500]
        }

    # Collect statistics
    products_path = output_dir / "products.parquet"
    if not products_path.exists():
        return {
            "status": "success_no_products",
            "class": np_class,
            "k": k,
            "message": "No products generated"
        }

    df_products = pd.read_parquet(products_path)
    num_products = len(df_products)
    num_unique = df_products['smiles'].nunique() if 'smiles' in df_products.columns else num_products

    df_parents = pd.read_parquet(input_parquet)
    num_parents = len(df_parents)

    LOG.info(f"✓ SUCCESS: {num_products} products from {num_parents} parents")

    return {
        "status": "success",
        "class": np_class,
        "k": k,
        "num_parents": num_parents,
        "num_products": num_products,
        "num_unique_products": num_unique,
        "products_per_parent": num_products / num_parents if num_parents > 0 else 0,
        "output_path": str(products_path)
    }


def generate_summary_report(
    results: List[Dict[str, Any]],
    output_path: Path
) -> None:
    """
    Generate markdown summary report.

    Args:
        results: List of enumeration results
        output_path: Path to output .md file
    """
    LOG.info(f"Generating summary report: {output_path}")

    lines = [
        "# Full-Scale NP Halogenation Library Generation Summary",
        "",
        f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Overview",
        "",
        f"- **Classes Processed:** {len(set(r['class'] for r in results if 'class' in r))}",
        f"- **Total Enumeration Jobs:** {len(results)}",
        f"- **Successful Jobs:** {sum(1 for r in results if r.get('status') == 'success')}",
        f"- **Failed Jobs:** {sum(1 for r in results if r.get('status') == 'error')}",
        "",
        "## Results by Class and K",
        "",
        "| Class | k | Parents | Products | Unique | Prod/Parent | Status |",
        "|-------|---|---------|----------|--------|-------------|--------|"
    ]

    for result in sorted(results, key=lambda x: (x.get('class', ''), x.get('k', 0))):
        cls = result.get('class', 'N/A')
        k = result.get('k', 'N/A')
        status = result.get('status', 'unknown')

        if status == 'success':
            parents = result.get('num_parents', 0)
            products = result.get('num_products', 0)
            unique = result.get('num_unique_products', 0)
            ppp = result.get('products_per_parent', 0)
            lines.append(f"| {cls} | {k} | {parents:,} | {products:,} | {unique:,} | {ppp:.1f} | ✅ |")
        else:
            msg = result.get('message', 'Unknown error')[:30]
            lines.append(f"| {cls} | {k} | - | - | - | - | ❌ {msg} |")

    lines.extend([
        "",
        "## Total Statistics",
        ""
    ])

    successful = [r for r in results if r.get('status') == 'success']
    if successful:
        total_parents = sum(r.get('num_parents', 0) for r in successful)
        total_products = sum(r.get('num_products', 0) for r in successful)
        total_unique = sum(r.get('num_unique_products', 0) for r in successful)

        lines.extend([
            f"- **Total Parents Processed:** {total_parents:,}",
            f"- **Total Products Generated:** {total_products:,}",
            f"- **Total Unique Products:** {total_unique:,}",
            f"- **Overall Products/Parent:** {total_products/total_parents:.1f}" if total_parents > 0 else "- **Overall Products/Parent:** N/A",
            ""
        ])

    lines.extend([
        "## Output Locations",
        "",
        "Generated libraries are stored in: `data/output/nplike/<class>-<k>X/products.parquet`",
        "",
        "## Next Steps",
        "",
        "1. Perform QC on generated libraries using 09_visualize_library.py",
        "2. Run 05_summaries.py to generate descriptor statistics",
        "3. Export to virtual screening formats using 11_export_for_vs.py",
        ""
    ])

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

    LOG.info(f"Summary report written: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Full-scale per-class halogenation library generation"
    )
    parser.add_argument(
        "--classes",
        nargs="+",
        choices=ALL_CLASSES,
        help="NP classes to enumerate (default: all enabled classes)"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Enumerate all available classes"
    )
    parser.add_argument(
        "--k-values",
        nargs="+",
        type=int,
        choices=[1, 2],
        default=[1, 2],
        help="K values to enumerate (default: 1 2)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5000,
        help="Batch size (default: 5000)"
    )
    parser.add_argument(
        "--rdkit-threads",
        type=int,
        default=8,
        help="RDKit threads (default: 8)"
    )
    parser.add_argument(
        "--max-parents",
        type=int,
        default=0,
        help="Max parents per class (0 = all; for testing)"
    )

    args = parser.parse_args()

    # Determine classes to process
    if args.all:
        classes = ALL_CLASSES
    elif args.classes:
        classes = args.classes
    else:
        # Default: main classes (exclude polysaccharide)
        classes = ["terpenoid", "alkaloid", "polyphenol", "glycoside", "lipid", "aa_peptide"]

    LOG.info(f"Classes to process: {classes}")
    LOG.info(f"K values: {args.k_values}")

    # Run enumerations
    results = []
    for np_class in classes:
        for k in args.k_values:
            result = run_enum_for_class(
                np_class,
                k,
                args.batch_size,
                args.rdkit_threads,
                args.max_parents
            )
            results.append(result)

    # Generate summary report
    report_path = Path("E:/Projects/halogenator/FULL_ENUM_SUMMARY.md")
    generate_summary_report(results, report_path)

    LOG.info("=" * 70)
    LOG.info("ORCHESTRATION COMPLETE")
    LOG.info(f"Summary: {report_path}")
    LOG.info("=" * 70)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
