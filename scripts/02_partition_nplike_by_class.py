#!/usr/bin/env python
# -*- coding: ascii -*-
"""
Partition CNPD-ETCM merged library into major natural product classes.

Outputs:
- A unified parquet with `np_primary_class` and `np_tags` columns added.
- Optional per-class parquet shards under the output directory.

Classification priority (first hit becomes primary; all hits become tags):
1. polysaccharide / glycoside
2. aa_peptide (amino acid / peptide / protein)
3. alkaloid
4. terpenoid
5. polyphenol / phenylpropanoid
6. lipid
7. other

Heuristics combine source metadata (class/subclass fields when present)
with structural cues (SMARTS-based).
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
LOG = logging.getLogger("partition_nplike")

# Precompiled SMARTS patterns
SUGAR_SMARTS = [
    Chem.MolFromSmarts("C1OCCC1"),  # furanose skeleton
    Chem.MolFromSmarts("C1OCCCC1"),  # pyranose skeleton
    Chem.MolFromSmarts("[C;R1]1[O;R1][C;R1][C;R1][C;R1][C;R1]1"),
    Chem.MolFromSmarts("[C;R1]1[O;R1][C;R1][C;R1][C;R1]1"),
]

PEPTIDE_BOND = Chem.MolFromSmarts("C(=O)N")
ALPHA_AA = Chem.MolFromSmarts("[NX3;H2,H1]-[CH1](-[CX4H2])[C](=O)[OX1H,OX2H1]")
BASIC_N = Chem.MolFromSmarts("[NX3;H2,H1;!$([N]-C(=O));!$([N+])]")
AROMATIC_N = Chem.MolFromSmarts("n")
ISOPROPYL = Chem.MolFromSmarts("C(C)(C)")
C12_CHAIN = Chem.MolFromSmarts("CCCCCCCCCCCC")
PHENOL_OH = Chem.MolFromSmarts("c[OX2H]")


def _lower_meta(row: Dict[str, any]) -> List[str]:
    """Collect lowercased metadata strings to aid keyword checks."""
    keys = ("np_class", "np_subclass", "class", "subclass", "kingdom", "category")
    values: List[str] = []
    for k in keys:
        v = row.get(k)
        if isinstance(v, str):
            values.append(v.lower())
    return values


def _has_keyword(meta_values: List[str], keywords: List[str]) -> bool:
    return any(any(kw in m for kw in keywords) for m in meta_values)


def _sugar_atoms(mol: Chem.Mol) -> Set[int]:
    atoms: Set[int] = set()
    for patt in SUGAR_SMARTS:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            atoms.update(match)
    return atoms


def _count_matches(mol: Chem.Mol, patt: Chem.Mol) -> int:
    try:
        return len(mol.GetSubstructMatches(patt))
    except Exception:
        return 0


def classify_molecule(row: Dict[str, any]) -> Tuple[str, List[str]]:
    """Return (primary_class, tags)."""
    tags: List[str] = []
    primary: Optional[str] = None

    smiles = row.get("smiles") or row.get("Smiles")
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return "invalid", ["invalid_smiles"]

    meta_values = _lower_meta(row)

    # Common descriptors
    num_rings = mol.GetRingInfo().NumRings()
    num_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    num_C = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    num_N = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")

    sugar_atoms = _sugar_atoms(mol)
    num_sugar_atoms = len(sugar_atoms)
    num_sugar_rings = _count_matches(mol, SUGAR_SMARTS[1]) + _count_matches(mol, SUGAR_SMARTS[0])
    aglycone_atoms = mol.GetNumHeavyAtoms() - num_sugar_atoms

    num_pept_bond = _count_matches(mol, PEPTIDE_BOND)
    has_alpha_aa = mol.HasSubstructMatch(ALPHA_AA)
    has_basic_N = mol.HasSubstructMatch(BASIC_N)
    has_aromatic_N = mol.HasSubstructMatch(AROMATIC_N)
    has_isopropyl = mol.HasSubstructMatch(ISOPROPYL)
    has_C12_chain = mol.HasSubstructMatch(C12_CHAIN)
    num_phenolic_OH = _count_matches(mol, PHENOL_OH)

    # 1) polysaccharide / glycoside
    if (
        _has_keyword(meta_values, ["glycoside", "saponin", "oligosaccharide", "polysaccharide"])
        or num_sugar_rings >= 1
    ):
        if num_sugar_rings >= 2 and aglycone_atoms <= 10:
            tags.append("polysaccharide")
            primary = primary or "polysaccharide"
        elif num_sugar_rings >= 1:
            tags.append("glycoside")
            primary = primary or "glycoside"

    # 2) amino acid / peptide
    if (
        _has_keyword(meta_values, ["peptide", "protein", "amino acid"])
        or has_alpha_aa
        or num_pept_bond >= 1
    ):
        tags.append("aa_peptide")
        primary = primary or "aa_peptide"

    # 3) alkaloid
    if (
        _has_keyword(meta_values, ["alkaloid"])
        or ((has_basic_N or has_aromatic_N) and num_rings >= 1)
    ):
        tags.append("alkaloid")
        primary = primary or "alkaloid"

    # 4) terpenoid
    if (
        _has_keyword(meta_values, ["terpen"])
        or (num_C >= 10 and num_N <= 1 and has_isopropyl)
    ):
        tags.append("terpenoid")
        primary = primary or "terpenoid"

    # 5) polyphenol / phenylpropanoid
    if (
        _has_keyword(meta_values, ["phenylpropanoid", "flavonoid", "isoflavonoid", "coumarin", "lignan"])
        or (num_arom_rings >= 1 and num_phenolic_OH >= 1)
    ):
        tags.append("polyphenol")
        primary = primary or "polyphenol"

    # 6) lipid
    if (
        _has_keyword(meta_values, ["lipid", "fatty", "glyceride", "phospho"])
        or (has_C12_chain and num_rings <= 2)
    ):
        tags.append("lipid")
        primary = primary or "lipid"

    if primary is None:
        primary = "other"
        tags.append("other")

    return primary, tags


def write_batch(rows: List[Dict[str, any]], writer: pq.ParquetWriter) -> None:
    if not rows:
        return
    table = pa.Table.from_pylist(rows)
    # Align schema to writer to avoid null/string mismatches across batches
    try:
        table = table.cast(writer.schema)
    except Exception:
        pass
    writer.write_table(table)
    rows.clear()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Partition CNPD-ETCM merged library into NP classes"
    )
    parser.add_argument("-i", "--input", required=True, help="Merged parquet file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--split", action="store_true", help="Also write per-class parquet files")
    parser.add_argument("--batch-size", type=int, default=2000, help="Rows per batch")
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    combined_path = outdir / "nplike_with_classes.parquet"
    per_class_writers: Dict[str, pq.ParquetWriter] = {}

    pq_file = pq.ParquetFile(str(input_path))
    schema = pq_file.schema_arrow
    combined_schema = schema.append(pa.field("np_primary_class", pa.string())).append(
        pa.field("np_tags", pa.string())
    )
    combined_writer = pq.ParquetWriter(str(combined_path), combined_schema)

    try:
        for batch_idx, batch in enumerate(pq_file.iter_batches(batch_size=args.batch_size)):
            rows = batch.to_pylist()
            combined_rows: List[Dict[str, any]] = []
            per_class_rows: Dict[str, List[Dict[str, any]]] = {}

            for row in rows:
                primary, tags = classify_molecule(row)
                if primary == "invalid":
                    continue  # drop invalid SMILES entirely

                row["np_primary_class"] = primary
                row["np_tags"] = json.dumps(tags, ensure_ascii=False)
                combined_rows.append(row)

                if args.split:
                    per_class_rows.setdefault(primary, []).append(row)

            write_batch(combined_rows, combined_writer)

            if args.split:
                for cls, cls_rows in per_class_rows.items():
                    if cls not in per_class_writers:
                        per_path = outdir / cls
                        per_path.mkdir(parents=True, exist_ok=True)
                        per_writer = pq.ParquetWriter(
                            str(per_path / "base.parquet"), combined_schema
                        )
                        per_class_writers[cls] = per_writer
                    write_batch(cls_rows, per_class_writers[cls])

            LOG.info(f"Processed batch {batch_idx+1}")
    finally:
        combined_writer.close()
        for writer in per_class_writers.values():
            writer.close()

    LOG.info(f"Done. Combined output: {combined_path}")
    if args.split:
        LOG.info(f"Per-class outputs under: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
