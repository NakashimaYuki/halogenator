#!/usr/bin/env python
"""抽样验证萜类分类的准确性"""

import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import Descriptors
import random

def analyze_sample(row):
    """详细分析单个分子"""
    smiles = row.get("smiles") or row.get("Smiles")
    mol_id = row.get("id") or row.get("ID") or "unknown"

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    num_C = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    num_N = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
    num_O = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
    mw = Descriptors.MolWt(mol)
    num_rings = mol.GetRingInfo().NumRings()

    # 检查是否是C5倍数（允许±3误差，考虑氧化/降解）
    c5_multiples = [10, 15, 20, 25, 30, 40]
    closest_c5 = min(c5_multiples, key=lambda x: abs(num_C - x))
    c5_deviation = abs(num_C - closest_c5)

    return {
        "id": mol_id,
        "smiles": smiles[:60] + "..." if len(smiles) > 60 else smiles,
        "C": num_C,
        "N": num_N,
        "O": num_O,
        "MW": f"{mw:.1f}",
        "rings": num_rings,
        "C5_fit": f"{closest_c5}±{c5_deviation}",
        "C5_ok": "OK" if c5_deviation <= 3 else "NO"
    }

def main():
    path = "E:/Projects/halogenator/data/output/nplike/nplike_with_classes.parquet"
    pq_file = pq.ParquetFile(path)

    terpenoid_rows = []
    for batch in pq_file.iter_batches(batch_size=1000):
        rows = batch.to_pylist()
        for row in rows:
            tags = row.get("np_tags", "[]")
            if "terpenoid" in tags or row.get("np_primary_class") == "terpenoid":
                terpenoid_rows.append(row)

    # 随机抽样30个
    samples = random.sample(terpenoid_rows, min(30, len(terpenoid_rows)))

    print("=== 萜类分子抽样验证（随机30个） ===\n")

    results = []
    for row in samples:
        result = analyze_sample(row)
        if result:
            results.append(result)

    # 统计C5倍数符合度
    c5_ok_count = sum(1 for r in results if r["C5_ok"] == "OK")

    print(f"{'ID':<15} {'C':<4} {'N':<3} {'O':<3} {'MW':<8} {'Rings':<6} {'C5匹配':<10} {'符合'}")
    print("-" * 80)
    for r in results:
        print(f"{r['id']:<15} {r['C']:<4} {r['N']:<3} {r['O']:<3} {r['MW']:<8} {r['rings']:<6} {r['C5_fit']:<10} {r['C5_ok']}")

    print("\n" + "=" * 80)
    print(f"C5倍数符合率: {c5_ok_count}/{len(results)} ({c5_ok_count/len(results)*100:.1f}%)")
    print(f"含氮比例: {sum(1 for r in results if int(r['N']) > 0)}/{len(results)} ({sum(1 for r in results if int(r['N']) > 0)/len(results)*100:.1f}%)")

if __name__ == "__main__":
    main()
