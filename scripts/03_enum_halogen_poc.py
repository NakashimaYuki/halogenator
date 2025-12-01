#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
POC (Proof of Concept) + QC Script for Per-Class Halogenation.

Usage:
    python 03_enum_halogen_poc.py --class terpenoid --k 1 --max-parents 1000
    python 03_enum_halogen_poc.py --class polyphenol --k 2 --max-parents 500

This script:
1. Samples N parents from the class base.parquet
2. Runs enum-parquet with class-specific configuration
3. Generates statistics (products/parent, unique products, rule breakdown)
4. Runs 09_visualize_library.py to create HTML gallery
5. Outputs a PocReport_<class>_k<k>.md
"""

import argparse
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any

import pandas as pd
import pyarrow.parquet as pq

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
LOG = logging.getLogger(__name__)


def sample_parents(
    input_parquet: Path,
    output_parquet: Path,
    max_parents: int,
    strategy: str = "random"
) -> int:
    """
    Sample parents from input parquet and write to output parquet.

    Args:
        input_parquet: Path to class base.parquet
        output_parquet: Path to output sampled parents
        max_parents: Number of parents to sample
        strategy: Sampling strategy ('random', 'first', 'diverse')

    Returns:
        Number of parents sampled
    """
    LOG.info(f"Sampling {max_parents} parents from {input_parquet}")

    df = pd.read_parquet(input_parquet)
    total_available = len(df)

    if max_parents >= total_available:
        LOG.warning(f"Requested {max_parents} but only {total_available} available, using all")
        sampled = df
    else:
        if strategy == "random":
            sampled = df.sample(n=max_parents, random_state=42)
        elif strategy == "first":
            sampled = df.head(max_parents)
        else:
            # TODO: implement diverse sampling (e.g., by MW bins)
            LOG.warning(f"Strategy '{strategy}' not implemented, using random")
            sampled = df.sample(n=max_parents, random_state=42)

    sampled.to_parquet(output_parquet, index=False)
    LOG.info(f"Sampled {len(sampled)} parents -> {output_parquet}")
    return len(sampled)


def run_enum_parquet(
    input_parquet: Path,
    outdir: Path,
    np_class: str,
    k: int,
    batch_size: int = 2000,
    rdkit_threads: int = 8
) -> Dict[str, Any]:
    """
    Run halogenator enum-parquet command.

    Returns:
        Statistics dict
    """
    LOG.info(f"Running enum-parquet for class='{np_class}' k={k}")

    cmd = [
        "halogenator",
        "enum-parquet",
        "--input-parquet", str(input_parquet),
        "--outdir", str(outdir),
        "--k", str(k),
        "--np-class", np_class,
        "--batch-size", str(batch_size),
        "--rdkit-threads", str(rdkit_threads)
    ]

    LOG.info(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        LOG.error(f"enum-parquet failed with code {result.returncode}")
        LOG.error(f"stderr: {result.stderr}")
        raise RuntimeError("enum-parquet failed")

    LOG.info("enum-parquet completed successfully")

    # Parse output for statistics
    products_path = outdir / "products.parquet"
    if not products_path.exists():
        raise FileNotFoundError(f"Products file not found: {products_path}")

    df_products = pd.read_parquet(products_path)
    total_products = len(df_products)
    unique_smiles = df_products['smiles'].nunique() if 'smiles' in df_products.columns else total_products

    stats = {
        'total_products': total_products,
        'unique_products': unique_smiles,
        'stdout': result.stdout,
        'stderr': result.stderr
    }

    LOG.info(f"Products generated: {total_products} (unique: {unique_smiles})")
    return stats


def compute_statistics(
    parents_path: Path,
    products_path: Path
) -> Dict[str, Any]:
    """
    Compute detailed statistics on enumeration results.

    Returns:
        Statistics dict with:
        - num_parents
        - num_products
        - num_unique_products
        - products_per_parent_mean
        - products_per_parent_median
        - products_per_parent_max
        - (TODO) per_rule_breakdown
    """
    LOG.info("Computing statistics...")

    df_parents = pd.read_parquet(parents_path)
    df_products = pd.read_parquet(products_path)

    num_parents = len(df_parents)
    num_products = len(df_products)
    num_unique = df_products['smiles'].nunique() if 'smiles' in df_products.columns else num_products

    # Compute products per parent
    if 'parent_smiles' in df_products.columns:
        products_per_parent = df_products.groupby('parent_smiles').size()
        ppp_mean = products_per_parent.mean()
        ppp_median = products_per_parent.median()
        ppp_max = products_per_parent.max()
    else:
        ppp_mean = num_products / num_parents if num_parents > 0 else 0
        ppp_median = ppp_mean
        ppp_max = ppp_mean

    stats = {
        'num_parents': num_parents,
        'num_products': num_products,
        'num_unique_products': num_unique,
        'products_per_parent_mean': float(ppp_mean),
        'products_per_parent_median': float(ppp_median),
        'products_per_parent_max': int(ppp_max) if ppp_max == ppp_max else 0,  # NaN check
    }

    LOG.info(f"Statistics: {json.dumps(stats, indent=2)}")
    return stats


def run_visualization(
    products_path: Path,
    output_html: Path,
    sample_size: int = 500
) -> None:
    """
    Run 09_visualize_library.py to create HTML gallery.

    Args:
        products_path: Path to products.parquet
        output_html: Path to output HTML gallery
        sample_size: Number of products to visualize
    """
    LOG.info(f"Generating HTML gallery with {sample_size} samples")

    # Use 09_visualize_library.py
    viz_script = Path(__file__).parent / "09_visualize_library.py"
    if not viz_script.exists():
        LOG.warning(f"Visualization script not found: {viz_script}, skipping visualization")
        return

    cmd = [
        "python",
        str(viz_script),
        "html",
        "-i", str(products_path),
        "-o", str(output_html),
        "--sample", str(sample_size)
    ]

    LOG.info(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        LOG.warning(f"Visualization failed: {result.stderr}")
    else:
        LOG.info(f"HTML gallery generated: {output_html}")


def generate_poc_report(
    np_class: str,
    k: int,
    stats: Dict[str, Any],
    output_report: Path
) -> None:
    """
    Generate markdown POC report.

    Args:
        np_class: NP class name
        k: Halogenation depth
        stats: Statistics dict
        output_report: Path to output .md file
    """
    LOG.info(f"Generating POC report: {output_report}")

    report_lines = [
        f"# POC Report: {np_class} k={k}",
        "",
        f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Summary",
        "",
        f"- **NP Class:** {np_class}",
        f"- **Halogenation Depth:** k={k}",
        f"- **Parents Processed:** {stats.get('num_parents', 'N/A')}",
        f"- **Total Products:** {stats.get('num_products', 'N/A')}",
        f"- **Unique Products:** {stats.get('num_unique_products', 'N/A')}",
        "",
        "## Products per Parent",
        "",
        f"- **Mean:** {stats.get('products_per_parent_mean', 0):.2f}",
        f"- **Median:** {stats.get('products_per_parent_median', 0):.2f}",
        f"- **Max:** {stats.get('products_per_parent_max', 'N/A')}",
        "",
        "## Verdict",
        "",
    ]

    # Add verdict based on products/parent ratio
    ppp_mean = stats.get('products_per_parent_mean', 0)
    if ppp_mean < 50:
        verdict = "✅ **PASS** - Products/parent ratio is reasonable"
    elif ppp_mean < 100:
        verdict = "⚠️ **WARNING** - Products/parent ratio is high, consider tightening rules"
    else:
        verdict = "❌ **FAIL** - Products/parent ratio is too high, rules too aggressive"

    report_lines.append(verdict)
    report_lines.extend([
        "",
        "## Next Steps",
        "",
        "1. Review HTML gallery for visual QC",
        "2. Adjust `configs/halogen_rules_by_class.yaml` if needed",
        "3. Re-run POC with adjusted parameters",
        "4. Proceed to full-scale enumeration if results are acceptable",
        ""
    ])

    with open(output_report, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))

    LOG.info(f"POC report written: {output_report}")


def main():
    parser = argparse.ArgumentParser(
        description="POC+QC script for per-class halogenation"
    )
    parser.add_argument(
        "--class",
        dest="np_class",
        required=True,
        choices=["terpenoid", "alkaloid", "polyphenol", "glycoside", "lipid", "aa_peptide", "other"],
        help="NP class to test"
    )
    parser.add_argument(
        "--k",
        type=int,
        required=True,
        choices=[1, 2],
        help="Halogenation depth"
    )
    parser.add_argument(
        "--max-parents",
        type=int,
        default=1000,
        help="Number of parents to sample (default: 1000)"
    )
    parser.add_argument(
        "--sample-strategy",
        default="random",
        choices=["random", "first", "diverse"],
        help="Sampling strategy (default: random)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=2000,
        help="Batch size for enumeration (default: 2000)"
    )
    parser.add_argument(
        "--rdkit-threads",
        type=int,
        default=8,
        help="RDKit threads (default: 8)"
    )
    parser.add_argument(
        "--viz-samples",
        type=int,
        default=500,
        help="Number of products to visualize (default: 500)"
    )

    args = parser.parse_args()

    # Setup paths
    base_dir = Path("E:/Projects/halogenator")
    data_dir = base_dir / "data" / "output" / "nplike"

    input_base_parquet = data_dir / args.np_class / "base.parquet"
    poc_outdir = data_dir / "poc" / f"{args.np_class}-k{args.k}"
    poc_outdir.mkdir(parents=True, exist_ok=True)

    sampled_parquet = poc_outdir / "sampled_parents.parquet"
    enum_outdir = poc_outdir / "enum"
    enum_outdir.mkdir(exist_ok=True)
    products_path = enum_outdir / "products.parquet"
    html_gallery = poc_outdir / "gallery.html"
    report_path = base_dir / f"PocReport_{args.np_class}_k{args.k}.md"

    # Check input exists
    if not input_base_parquet.exists():
        LOG.error(f"Input parquet not found: {input_base_parquet}")
        return 1

    try:
        # Step 1: Sample parents
        num_sampled = sample_parents(
            input_base_parquet,
            sampled_parquet,
            args.max_parents,
            args.sample_strategy
        )

        # Step 2: Run enumeration
        enum_stats = run_enum_parquet(
            sampled_parquet,
            enum_outdir,
            args.np_class,
            args.k,
            args.batch_size,
            args.rdkit_threads
        )

        # Step 3: Compute statistics
        stats = compute_statistics(sampled_parquet, products_path)

        # Step 4: Generate visualization
        run_visualization(products_path, html_gallery, args.viz_samples)

        # Step 5: Generate report
        generate_poc_report(args.np_class, args.k, stats, report_path)

        LOG.info("=" * 60)
        LOG.info("POC+QC COMPLETE")
        LOG.info(f"Report: {report_path}")
        LOG.info(f"Gallery: {html_gallery}")
        LOG.info(f"Products: {products_path}")
        LOG.info("=" * 60)

        return 0

    except Exception as e:
        LOG.error(f"POC failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
