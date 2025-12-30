"""Verify k=1 enumeration results and check ALPHA_CARBONYL presence"""
import pandas as pd
from pathlib import Path

classes = ['lipid', 'aa_peptide', 'polyphenol', 'alkaloid', 'terpenoid', 'glycoside']
base_dir = Path('E:/Projects/halogenator/data/output/nplike')

print("=" * 80)
print("k=1 ENUMERATION VERIFICATION - RULE DISTRIBUTIONS")
print("=" * 80)

total_products = 0
total_alpha = 0
all_rules = {}

for cls in classes:
    products_file = base_dir / f'{cls}-1X' / 'products.parquet'

    if not products_file.exists():
        print(f"\n{cls}: FILE NOT FOUND - {products_file}")
        continue

    df = pd.read_parquet(products_file)
    rule_dist = df['rule'].value_counts()

    total_products += len(df)

    print(f"\n{cls.upper()}: {len(df):,} products")
    print("-" * 60)

    for rule, count in rule_dist.items():
        pct = count / len(df) * 100
        print(f"  {rule:30s}: {count:8,} ({pct:5.1f}%)")

        # Aggregate across classes
        if rule not in all_rules:
            all_rules[rule] = 0
        all_rules[rule] += count

        if rule == 'ALPHA_CARBONYL__CH2__TO__X':
            total_alpha += count

    # Check if ALPHA_CARBONYL is present
    if 'ALPHA_CARBONYL__CH2__TO__X' in rule_dist.index:
        alpha_count = rule_dist['ALPHA_CARBONYL__CH2__TO__X']
        print(f"  {'':30s}  [OK] ALPHA_CARBONYL present: {alpha_count:,}")
    else:
        print(f"  {'':30s}  [MISSING] ALPHA_CARBONYL NOT FOUND")

print("\n" + "=" * 80)
print("OVERALL SUMMARY")
print("=" * 80)
print(f"Total products: {total_products:,}")
print(f"\nRule totals across all classes:")
for rule in sorted(all_rules.keys(), key=lambda x: all_rules[x], reverse=True):
    count = all_rules[rule]
    pct = count / total_products * 100
    marker = " ← ALPHA_CARBONYL" if rule == 'ALPHA_CARBONYL__CH2__TO__X' else ""
    print(f"  {rule:30s}: {count:10,} ({pct:5.1f}%){marker}")

print("\n" + "=" * 80)
print("ALPHA_CARBONYL VALIDATION")
print("=" * 80)
if total_alpha > 0:
    print(f"[SUCCESS] ALPHA_CARBONYL produced {total_alpha:,} products ({total_alpha/total_products*100:.2f}%)")
    print(f"  Expected contribution: 5-15%")
    print(f"  Actual contribution: {total_alpha/total_products*100:.2f}%")
    if total_alpha/total_products*100 >= 5:
        print(f"  Status: EXCELLENT")
    else:
        print(f"  Status: LOW (needs investigation)")
else:
    print(f"[FAILURE] ALPHA_CARBONYL produced 0 products")
    print(f"  BUG STILL PRESENT OR REGRESSION")

print("=" * 80)
