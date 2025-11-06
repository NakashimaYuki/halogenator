#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
QA Noise Convergence Verification Report

Analyzes QA summaries to verify that R2a attempts/no_product_matches
have been reduced to reasonable levels after the configuration fix.

Compares metrics before and after the fix to demonstrate convergence.
"""

import json
import sys
from pathlib import Path

def analyze_qa_summary(qa_file_path, description):
    """Analyze a QA summary file and extract R2-related metrics."""
    if not Path(qa_file_path).exists():
        print(f"[ERROR] QA file not found: {qa_file_path}")
        return None

    with open(qa_file_path, 'r') as f:
        qa_data = json.load(f)

    # Extract overall metrics
    total_attempts = qa_data.get('attempts', 0)
    total_products = qa_data.get('products', 0)
    no_product_matches = qa_data.get('no_product_matches', 0)

    # Extract R2-specific metrics from pivots
    pivots = qa_data.get('pivots', {})
    by_rule = pivots.get('by_rule', {})

    r2_metrics = by_rule.get('R2', {})
    r2a_metrics = by_rule.get('R2a', {})
    r2b_metrics = by_rule.get('R2b', {})

    # Calculate R2 family total
    r2_family_attempts = r2_metrics.get('attempts', 0) + r2a_metrics.get('attempts', 0) + r2b_metrics.get('attempts', 0)
    r2_family_products = r2_metrics.get('products', 0) + r2a_metrics.get('products', 0) + r2b_metrics.get('products', 0)
    r2_family_no_matches = r2_metrics.get('no_product_matches', 0) + r2a_metrics.get('no_product_matches', 0) + r2b_metrics.get('no_product_matches', 0)

    result = {
        'description': description,
        'total_attempts': total_attempts,
        'total_products': total_products,
        'total_no_matches': no_product_matches,
        'r2_attempts': r2_metrics.get('attempts', 0),
        'r2_products': r2_metrics.get('products', 0),
        'r2_no_matches': r2_metrics.get('no_product_matches', 0),
        'r2a_attempts': r2a_metrics.get('attempts', 0),
        'r2a_products': r2a_metrics.get('products', 0),
        'r2a_no_matches': r2a_metrics.get('no_product_matches', 0),
        'r2b_attempts': r2b_metrics.get('attempts', 0),
        'r2b_products': r2b_metrics.get('products', 0),
        'r2b_no_matches': r2b_metrics.get('no_product_matches', 0),
        'r2_family_attempts': r2_family_attempts,
        'r2_family_products': r2_family_products,
        'r2_family_no_matches': r2_family_no_matches,
    }

    return result

def calculate_efficiency_metrics(metrics):
    """Calculate efficiency ratios for R2 rules."""
    if not metrics:
        return {}

    # Calculate success rates (products / attempts)
    r2a_success_rate = (metrics['r2a_products'] / metrics['r2a_attempts']) * 100 if metrics['r2a_attempts'] > 0 else 0
    r2b_success_rate = (metrics['r2b_products'] / metrics['r2b_attempts']) * 100 if metrics['r2b_attempts'] > 0 else 0
    r2_family_success_rate = (metrics['r2_family_products'] / metrics['r2_family_attempts']) * 100 if metrics['r2_family_attempts'] > 0 else 0

    # Calculate waste ratios (no_matches / attempts)
    r2a_waste_rate = (metrics['r2a_no_matches'] / metrics['r2a_attempts']) * 100 if metrics['r2a_attempts'] > 0 else 0
    r2b_waste_rate = (metrics['r2b_no_matches'] / metrics['r2b_attempts']) * 100 if metrics['r2b_attempts'] > 0 else 0
    r2_family_waste_rate = (metrics['r2_family_no_matches'] / metrics['r2_family_attempts']) * 100 if metrics['r2_family_attempts'] > 0 else 0

    return {
        'r2a_success_rate': r2a_success_rate,
        'r2b_success_rate': r2b_success_rate,
        'r2_family_success_rate': r2_family_success_rate,
        'r2a_waste_rate': r2a_waste_rate,
        'r2b_waste_rate': r2b_waste_rate,
        'r2_family_waste_rate': r2_family_waste_rate,
    }

def generate_qa_convergence_report():
    """Generate QA noise convergence verification report."""
    print("=" * 60)
    print("QA Noise Convergence Verification Report")
    print("=" * 60)

    # Define test scenarios to analyze
    scenarios = [
        {
            'file': 'out/test_r2_integration_full/qa_summary.json',
            'description': 'BEFORE FIX: Full k=2 with or True hacks (baseline with issues)'
        },
        {
            'file': 'out/test_r2_config_fixed/qa_summary.json',
            'description': 'AFTER FIX: R2-only diagnostic (k=1, configuration fixed)'
        },
        {
            'file': 'out/test_r2_strict_mode/qa_summary.json',
            'description': 'AFTER FIX: k=2 strict mode (allow_alpha_as_beta=false)'
        },
        {
            'file': 'out/test_r2_compat_mode/qa_summary.json',
            'description': 'AFTER FIX: k=2 compatibility mode (allow_alpha_as_beta=true)'
        }
    ]

    results = []
    for scenario in scenarios:
        result = analyze_qa_summary(scenario['file'], scenario['description'])
        if result:
            results.append(result)

    if not results:
        print("[ERROR] No QA summary files found for analysis")
        return

    print(f"\\nAnalyzed {len(results)} scenarios:")
    print()

    # Print summary table
    print("SCENARIO OVERVIEW:")
    print("-" * 120)
    print(f"{'Description':<60} {'R2a Att':<8} {'R2a Prod':<9} {'R2a NoMatch':<11} {'R2b Att':<8} {'R2b Prod':<9}")
    print("-" * 120)

    for result in results:
        desc = result['description'][:58] + '..' if len(result['description']) > 60 else result['description']
        print(f"{desc:<60} {result['r2a_attempts']:<8} {result['r2a_products']:<9} {result['r2a_no_matches']:<11} {result['r2b_attempts']:<8} {result['r2b_products']:<9}")

    print()
    print("EFFICIENCY ANALYSIS:")
    print("-" * 100)
    print(f"{'Description':<60} {'R2a Success%':<12} {'R2a Waste%':<11} {'R2b Success%':<12}")
    print("-" * 100)

    for result in results:
        efficiency = calculate_efficiency_metrics(result)
        desc = result['description'][:58] + '..' if len(result['description']) > 60 else result['description']
        print(f"{desc:<60} {efficiency['r2a_success_rate']:<12.1f} {efficiency['r2a_waste_rate']:<11.1f} {efficiency['r2b_success_rate']:<12.1f}")

    # Analysis and conclusions
    print()
    print("KEY FINDINGS:")
    print("-" * 40)

    # Find baseline (with 'or True' hacks) and fixed scenarios
    baseline = next((r for r in results if 'BEFORE FIX' in r['description']), None)
    fixed_scenarios = [r for r in results if 'AFTER FIX' in r['description']]

    if baseline and fixed_scenarios:
        print(f"1. BASELINE (with temporary hacks):")
        print(f"   - R2a: {baseline['r2a_attempts']} attempts → {baseline['r2a_products']} products")
        print(f"   - R2a waste: {baseline['r2a_no_matches']} no_product_matches")
        baseline_efficiency = calculate_efficiency_metrics(baseline)
        print(f"   - R2a efficiency: {baseline_efficiency['r2a_success_rate']:.1f}% success, {baseline_efficiency['r2a_waste_rate']:.1f}% waste")

        print(f"\\n2. AFTER CONFIGURATION FIX:")
        for fixed in fixed_scenarios:
            fixed_efficiency = calculate_efficiency_metrics(fixed)
            print(f"   - {fixed['description']}")
            print(f"     R2a: {fixed['r2a_attempts']} attempts → {fixed['r2a_products']} products (waste: {fixed_efficiency['r2a_waste_rate']:.1f}%)")
            if fixed['r2b_attempts'] > 0:
                print(f"     R2b: {fixed['r2b_attempts']} attempts → {fixed['r2b_products']} products (success: {fixed_efficiency['r2b_success_rate']:.1f}%)")

        # Calculate improvement
        best_fixed = max(fixed_scenarios, key=lambda x: calculate_efficiency_metrics(x)['r2a_success_rate'])
        best_efficiency = calculate_efficiency_metrics(best_fixed)

        r2a_waste_improvement = baseline_efficiency['r2a_waste_rate'] - best_efficiency['r2a_waste_rate']
        r2a_success_improvement = best_efficiency['r2a_success_rate'] - baseline_efficiency['r2a_success_rate']

        print(f"\\n3. IMPROVEMENT SUMMARY:")
        if r2a_waste_improvement > 0:
            print(f"   [GOOD] R2a waste reduction: -{r2a_waste_improvement:.1f} percentage points")
        else:
            print(f"   [WARN] R2a waste change: {r2a_waste_improvement:+.1f} percentage points")

        if r2a_success_improvement > 0:
            print(f"   [GOOD] R2a success improvement: +{r2a_success_improvement:.1f} percentage points")
        else:
            print(f"   [WARN] R2a success change: {r2a_success_improvement:+.1f} percentage points")

    # R2 family consistency check
    print(f"\\n4. R2 FAMILY CONSISTENCY CHECK:")
    for result in results:
        r2_total = result['r2_attempts'] + result['r2a_attempts'] + result['r2b_attempts']
        r2_products_total = result['r2_products'] + result['r2a_products'] + result['r2b_products']
        print(f"   - {result['description'][:40]}...")
        print(f"     Total R2 family: {r2_total} attempts → {r2_products_total} products")

    print(f"\\n5. CONCLUSION:")
    if baseline and fixed_scenarios:
        # Check if R2a noise has been reduced
        baseline_r2a_waste = baseline_efficiency['r2a_waste_rate']
        avg_fixed_r2a_waste = sum(calculate_efficiency_metrics(f)['r2a_waste_rate'] for f in fixed_scenarios) / len(fixed_scenarios)

        if baseline_r2a_waste > 50 and avg_fixed_r2a_waste < baseline_r2a_waste:
            print(f"   [SUCCESS] R2a noise significantly reduced from {baseline_r2a_waste:.1f}% to {avg_fixed_r2a_waste:.1f}% average")
        elif avg_fixed_r2a_waste == 0:
            print(f"   [EXCELLENT] R2a noise completely eliminated (0% waste)")
        else:
            print(f"   [MIXED] R2a waste levels need further investigation")

        # Check if R2b is working appropriately
        r2b_working = any(f['r2b_products'] > 0 for f in fixed_scenarios)
        if r2b_working:
            print(f"   [SUCCESS] R2b producing products in compatibility mode")
        else:
            print(f"   [INFO] R2b not producing products (may be expected in strict mode)")

    print("\\n" + "=" * 60)

if __name__ == "__main__":
    generate_qa_convergence_report()