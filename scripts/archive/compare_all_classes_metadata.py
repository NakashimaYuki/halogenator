#!/usr/bin/env python
"""对比所有类别的元数据 vs 结构识别依赖度"""

import pyarrow.parquet as pq
from rdkit import Chem

# 复制SMARTS模式
SUGAR_SMARTS = [
    Chem.MolFromSmarts("C1OCCC1"),
    Chem.MolFromSmarts("C1OCCCC1"),
]
PEPTIDE_BOND = Chem.MolFromSmarts("C(=O)N")
ALPHA_AA = Chem.MolFromSmarts("[NX3;H2,H1]-[CH1](-[CX4H2])[C](=O)[OX1H,OX2H1]")
BASIC_N = Chem.MolFromSmarts("[NX3;H2,H1;!$([N]-C(=O));!$([N+])]")
AROMATIC_N = Chem.MolFromSmarts("n")
ISOPROPYL = Chem.MolFromSmarts("C(C)(C)")
C12_CHAIN = Chem.MolFromSmarts("CCCCCCCCCCCC")
PHENOL_OH = Chem.MolFromSmarts("c[OX2H]")

def _lower_meta(row):
    keys = ("np_class", "np_subclass", "class", "subclass", "kingdom", "category")
    values = []
    for k in keys:
        v = row.get(k)
        if isinstance(v, str):
            values.append(v.lower())
    return values

def _has_keyword(meta_values, keywords):
    return any(any(kw in m for kw in keywords) for m in meta_values)

def analyze_class(class_name, rows):
    """分析单个类别的元数据覆盖率"""
    metadata_only = 0
    structure_only = 0
    both = 0

    for row in rows:
        meta_values = _lower_meta(row)
        smiles = row.get("smiles") or row.get("Smiles")
        mol = Chem.MolFromSmiles(smiles) if smiles else None

        if not mol:
            continue

        # 检查元数据和结构匹配
        has_meta = False
        has_structure = False

        if class_name == "glycoside":
            has_meta = _has_keyword(meta_values, ["glycoside", "saponin"])
            num_sugar_rings = sum(len(mol.GetSubstructMatches(p)) for p in SUGAR_SMARTS if p)
            has_structure = num_sugar_rings >= 1

        elif class_name == "aa_peptide":
            has_meta = _has_keyword(meta_values, ["peptide", "protein", "amino acid"])
            has_structure = (mol.HasSubstructMatch(ALPHA_AA) or
                           len(mol.GetSubstructMatches(PEPTIDE_BOND)) >= 1)

        elif class_name == "alkaloid":
            has_meta = _has_keyword(meta_values, ["alkaloid"])
            num_rings = mol.GetRingInfo().NumRings()
            has_structure = ((mol.HasSubstructMatch(BASIC_N) or
                            mol.HasSubstructMatch(AROMATIC_N)) and num_rings >= 1)

        elif class_name == "terpenoid":
            has_meta = _has_keyword(meta_values, ["terpen"])
            num_C = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
            num_N = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
            has_structure = (num_C >= 10 and num_N <= 1 and
                           mol.HasSubstructMatch(ISOPROPYL))

        elif class_name == "polyphenol":
            has_meta = _has_keyword(meta_values, ["phenylpropanoid", "flavonoid",
                                                  "isoflavonoid", "coumarin", "lignan"])
            from rdkit.Chem import rdMolDescriptors
            num_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            num_phenolic_OH = len(mol.GetSubstructMatches(PHENOL_OH))
            has_structure = (num_arom_rings >= 1 and num_phenolic_OH >= 1)

        elif class_name == "lipid":
            has_meta = _has_keyword(meta_values, ["lipid", "fatty", "glyceride", "phospho"])
            num_rings = mol.GetRingInfo().NumRings()
            has_structure = (mol.HasSubstructMatch(C12_CHAIN) and num_rings <= 2)

        if has_meta and has_structure:
            both += 1
        elif has_meta:
            metadata_only += 1
        elif has_structure:
            structure_only += 1

    total = metadata_only + structure_only + both
    return {
        "total": total,
        "metadata_only": metadata_only,
        "structure_only": structure_only,
        "both": both,
        "meta_coverage": (metadata_only + both) / total * 100 if total > 0 else 0,
        "struct_coverage": (structure_only + both) / total * 100 if total > 0 else 0
    }

def main():
    path = "E:/Projects/halogenator/data/output/nplike/nplike_with_classes.parquet"
    pq_file = pq.ParquetFile(path)

    # 收集各类别的分子
    classes_data = {
        "glycoside": [],
        "aa_peptide": [],
        "alkaloid": [],
        "terpenoid": [],
        "polyphenol": [],
        "lipid": []
    }

    for batch in pq_file.iter_batches(batch_size=1000):
        rows = batch.to_pylist()
        for row in rows:
            primary = row.get("np_primary_class", "")
            if primary in classes_data:
                classes_data[primary].append(row)

    print("=== 各类别元数据 vs 结构识别依赖度对比 ===\n")
    print(f"{'类别':<15} {'总数':<8} {'仅元数据':<10} {'仅结构':<10} {'两者':<8} {'元数据%':<10} {'结构%'}")
    print("-" * 85)

    for class_name, rows in classes_data.items():
        if not rows:
            continue
        stats = analyze_class(class_name, rows[:5000])  # 限制5000个样本以加速
        print(f"{class_name:<15} {stats['total']:<8} "
              f"{stats['metadata_only']:<10} {stats['structure_only']:<10} "
              f"{stats['both']:<8} {stats['meta_coverage']:<10.1f} {stats['struct_coverage']:.1f}")

if __name__ == "__main__":
    main()
