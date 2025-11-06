#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export derivatives per parent molecule into SDF or MOL2 (one file per parent).
- Input: products.parquet or products.csv from k=2 enumeration
- Group key: parent_key > parent_inchikey (fallback)
- Optional: 3D coordinates via RDKit ETKDG + UFF
- One file per parent: {safe_parent_label}.sdf/.mol2
"""

import argparse, os, sys, re, pathlib
from typing import List, Optional, Iterable

# --- deps check ---
try:
    import pandas as pd
except Exception as e:
    print("[ERROR] pandas not available. pip install pandas pyarrow", file=sys.stderr); sys.exit(2)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception as e:
    print("[ERROR] RDKit not available. Install RDKit (e.g. conda install -c conda-forge rdkit).", file=sys.stderr); sys.exit(2)


def load_products(path: str) -> "pd.DataFrame":
    path = pathlib.Path(path)
    if not path.exists():
        print(f"[ERROR] products not found: {path}", file=sys.stderr); sys.exit(2)

    if path.suffix.lower() in [".parquet", ".pq"]:
        try:
            return pd.read_parquet(path)
        except Exception as e:
            print(f"[WARN] failed to read parquet: {e}; trying csv...", file=sys.stderr)

    # CSV fallback
    return pd.read_csv(path)


def sanitize_filename(name: str) -> str:
    """ASCII-only + safe chars; cut overly long names."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    if len(s) > 140:  # avoid crazy-long filenames
        s = s[:140]
    return s or "parent"


def pick_group_key(df: "pd.DataFrame") -> str:
    for key in ["parent_key", "parent_inchikey", "parent_smiles"]:
        if key in df.columns:
            return key
    raise RuntimeError("No parent identifier column found (expected one of: parent_key, parent_inchikey, parent_smiles).")


def smiles_to_mol(smiles: str, add3d: bool = True, seed: int = 1337) -> Optional["Chem.Mol"]:
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        m = Chem.AddHs(m)
        if add3d:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            if AllChem.EmbedMolecule(m, params) != 0:
                # Embedding failed: try without constraints
                if AllChem.EmbedMolecule(m, AllChem.ETKDG()) != 0:
                    return None
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
        return m
    except Exception:
        return None


def set_props(m: "Chem.Mol", row: "pd.Series", title_fmt: str) -> None:
    # Title
    title = title_fmt.format(**{c: str(row.get(c, "")) for c in row.index})
    m.SetProp("_Name", title)

    # Write key columns as SDF data fields
    for col in [
        "inchikey", "smiles", "parent_inchikey", "parent_smiles", "parent_name",
        "k", "rule", "halogen", "constraints_ok", "mw", "logp", "tpsa", "hba", "hbd",
        "aromatic_rings", "sanitize_ok", "pains_flags"
    ]:
        if col in row and pd.notna(row[col]):
            m.SetProp(col, str(row[col]))


def write_group_sdf(rows: "pd.DataFrame", out_path: str, add3d: bool, dedupe: bool, title_fmt: str, max_n: int, parent_smiles: str = None, parent_name: str = None):
    writer = Chem.SDWriter(out_path)
    writer.SetKekulize(False)
    seen = set()
    n = 0

    # First, write the parent molecule if provided
    if parent_smiles and parent_name:
        parent_mol = smiles_to_mol(parent_smiles, add3d=add3d)
        if parent_mol is not None:
            parent_mol.SetProp("_Name", f"{parent_name}_parent")
            parent_mol.SetProp("smiles", parent_smiles)
            parent_mol.SetProp("parent_name", parent_name)
            parent_mol.SetProp("molecule_type", "parent")
            if dedupe:
                parent_key = Chem.MolToSmiles(Chem.MolFromSmiles(parent_smiles), isomericSmiles=True, canonical=True)
                seen.add(parent_key)
            writer.write(parent_mol)
            n += 1

    # Then write all derivatives for this parent
    for _, r in rows.iterrows():
        if max_n and n >= max_n: break
        smiles = r.get("smiles", None)
        if not isinstance(smiles, str): continue
        if dedupe:
            from rdkit.Chem.rdmolfiles import MolToSmiles
            # fast dedupe by canonical smiles string
            key = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True, canonical=True)
            if key in seen: continue
            seen.add(key)

        m = smiles_to_mol(smiles, add3d=add3d)
        if m is None: continue
        set_props(m, r, title_fmt)
        # Mark as derivative
        m.SetProp("molecule_type", "derivative")
        writer.write(m); n += 1
    writer.close()
    return n


def write_group_mol2(rows: "pd.DataFrame", out_path: str, add3d: bool, dedupe: bool, title_fmt: str, max_n: int):
    with open(out_path, "w", encoding="utf-8") as f:
        seen = set(); n = 0
        for _, r in rows.iterrows():
            if max_n and n >= max_n: break
            smiles = r.get("smiles", None)
            if not isinstance(smiles, str): continue
            if dedupe:
                key = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True, canonical=True)
                if key in seen: continue
                seen.add(key)

            m = smiles_to_mol(smiles, add3d=add3d)
            if m is None: continue
            set_props(m, r, title_fmt)
            block = Chem.MolToMol2Block(m)
            f.write(block + "\n")
            n += 1
    return n


def read_parent_list(arg: Optional[str]) -> Optional[List[str]]:
    if not arg: return None
    if os.path.exists(arg):
        with open(arg, "r", encoding="utf-8") as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    # comma separated list
    return [x.strip() for x in arg.split(",") if x.strip()]


def read_parent_smiles(smiles_file: str) -> dict:
    """Read parent SMILES from tab-separated file (SMILES\tname format)."""
    parents = {}
    if not smiles_file or not os.path.exists(smiles_file):
        return parents

    with open(smiles_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                smiles, name = parts[0], parts[1]
                parents[name] = smiles
    return parents


def main():
    ap = argparse.ArgumentParser(description="Export k=2 derivatives per parent into SDF/MOL2.")
    ap.add_argument("--products", required=True, help="Path to products.parquet (or products.csv).")
    ap.add_argument("--outdir", required=True, help="Output directory, per-parent files will be created here.")
    ap.add_argument("--format", choices=["sdf","mol2"], default="sdf")
    ap.add_argument("--parents", help="Comma-separated parent ids or a file path (one id per line).")
    ap.add_argument("--parent-smiles", help="Parent SMILES file (tab-separated: SMILES\\tname) to include parent molecules first.")
    ap.add_argument("--name-by", choices=["parent_name","parent_inchikey","parent_key","auto"], default="auto",
                    help="How to name each file (auto: choose the best available).")
    ap.add_argument("--title", default="{rule}_{halogen}_k{k}", help="Per-molecule title in SDF/MOL2.")
    ap.add_argument("--max-per-parent", type=int, default=0, help="Limit exported molecules per parent (0=no limit).")
    ap.add_argument("--no3d", action="store_true", help="Do NOT generate 3D coordinates.")
    ap.add_argument("--dedupe", action="store_true", help="Dedupe by canonical SMILES per parent.")
    args = ap.parse_args()

    df = load_products(args.products)
    gkey = pick_group_key(df)

    # Read parent SMILES if provided
    parent_molecules = {}
    if args.parent_smiles:
        parent_molecules = read_parent_smiles(args.parent_smiles)
        print(f"[INFO] loaded {len(parent_molecules)} parent molecules from {args.parent_smiles}")

    # Filter to requested parents (optional)
    selected = read_parent_list(args.parents)
    if selected:
        df = df[df[gkey].isin(set(selected))].copy()
        if df.empty:
            print(f"[ERROR] No rows matched requested parents via key={gkey}", file=sys.stderr); sys.exit(3)

    # Choose filename label column
    name_col = None
    if args.name_by == "auto":
        for c in ["parent_name", "parent_inchikey", "parent_key"]:
            if c in df.columns: name_col = c; break
        if name_col is None: name_col = gkey
    else:
        name_col = args.name_by if args.name_by in df.columns else gkey

    outdir = pathlib.Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    add3d = not args.no3d
    total_files, total_mols = 0, 0
    for pid, rows in df.groupby(gkey):
        label = rows.iloc[0].get(name_col, pid)
        fname = sanitize_filename(str(label))
        fpath = outdir / f"{fname}.{args.format}"

        # Get parent information for this group
        parent_smiles = parent_molecules.get(str(label), None)
        parent_name = str(label) if parent_smiles else None

        if args.format == "sdf":
            n = write_group_sdf(rows, str(fpath), add3d=add3d, dedupe=args.dedupe,
                                title_fmt=args.title, max_n=args.max_per_parent,
                                parent_smiles=parent_smiles, parent_name=parent_name)
        else:
            n = write_group_mol2(rows, str(fpath), add3d=add3d, dedupe=args.dedupe,
                                 title_fmt=args.title, max_n=args.max_per_parent)

        print(f"[INFO] wrote {n} molecules -> {fpath.name}")
        if n > 0: total_files += 1; total_mols += n

    print(f"[DONE] {total_files} files, {total_mols} molecules exported.")


if __name__ == "__main__":
    main()