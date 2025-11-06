#!/usr/bin/env python
# -*- coding: ascii -*-
"""
Verify parquet schema fields for PR2 compliance.
Usage: python scripts/verify_parquet_fields.py <parquet_path>
"""
import pandas as pd
import sys

if len(sys.argv) < 2:
    print("Usage: python scripts/verify_parquet_fields.py <parquet_path>")
    sys.exit(1)

path = sys.argv[1]
df = pd.read_parquet(path)
cols = set(df.columns)

print("==== Parquet Field Verification ====")
print(f"File: {path}")
print(f"Total products: {len(df)}")
print()

print("Field checks:")
print(f"  Has sub_rule: {'YES' if 'sub_rule' in cols else 'NO'}")
print(f"  Has detection: {'YES' if 'detection' in cols else 'NO'}")
print(f"  Has rule_family: {'YES' if 'rule_family' in cols else 'NO'}")
print()

print("Family counts:")
if 'rule_family' in cols:
    print(df.groupby('rule_family')['smiles'].count().sort_values(ascending=False))
else:
    print("  (rule_family field not found)")
print()

r2 = df[df['rule_family']=='R2'] if 'rule_family' in cols else pd.DataFrame()
print(f"R2 family products: {len(r2)}")
if len(r2) > 0:
    if 'sub_rule' in r2.columns:
        print("  R2 sub_rule distribution:")
        print("   ", r2['sub_rule'].value_counts(dropna=False).to_dict())
    if 'detection' in r2.columns:
        print("  R2 detection distribution:")
        print("   ", r2['detection'].value_counts(dropna=False).to_dict())
print()

print(f"k == k_ops consistency: {'YES' if (df['k']==df['k_ops']).all() else 'NO'}")
if 'k' in cols and 'k_ops' in cols:
    mismatches = df[df['k'] != df['k_ops']]
    if len(mismatches) > 0:
        print(f"  Mismatches: {len(mismatches)}")
        print(f"  Sample: k={mismatches['k'].iloc[0]}, k_ops={mismatches['k_ops'].iloc[0]}")
print()

print("==== Verification Complete ====")
