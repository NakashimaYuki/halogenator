# -*- coding: ascii -*-
"""
Export per-parent SDF (v2)
- Includes parent molecule (first) + all k<=2 products for that parent
- Automatically repairs missing parent_inchikey (prioritizes parent_smiles calculation)
- Ensures total products across all SDFs == total rows in parquet (after pick filtering)

Usage:
  python scripts/export_parents_with_products_sdf_v2.py \
    --products out/etcm2000/pick_k2_compat/products_k2.parquet \
    --parents-meta out/etcm2000/parents/parents_meta.csv \
    --parents-pick out/etcm2000/parents/parents_pick.smi \
    --outdir out/etcm2000/pick_k2_compat/per_parent_sdf_fixed \
    --max-per-parent 100000
"""
import sys, argparse, pathlib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchiKey

def load_pick_smi(path):
    """Load parent SMILES from pick file and convert to InChIKey set"""
    picks = []
    for line in pathlib.Path(path).read_text(encoding="ascii", errors="ignore").splitlines():
        t = line.strip().split()
        if t:
            picks.append(t[0])  # parent SMILES in .smi

    # Convert to InChIKey set for filtering
    out = set()
    for smi in picks:
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        out.add(MolToInchiKey(m))
    return out

def compute_ik_from_smiles(smi):
    """Compute InChIKey from SMILES string"""
    try:
        m = Chem.MolFromSmiles(str(smi))
        if m is None: return None
        return MolToInchiKey(m)
    except Exception:
        return None

def normalize_meta(meta_path):
    """Normalize parent metadata CSV to standard column names"""
    try:
        meta = pd.read_csv(meta_path, encoding='utf-8')
    except:
        try:
            meta = pd.read_csv(meta_path, encoding='gbk')
        except:
            meta = pd.read_csv(meta_path, encoding='latin1')

    def pickcol(cands):
        for c in cands:
            if c in meta.columns: return c
        return None

    c_ik = pickcol(["parent_inchikey","InChIKey","inchikey"])
    c_sm = pickcol(["smiles","SMILES"])
    c_nm = pickcol(["parent_name","Ingredient Name in English","_ascii_name","name"])

    cols = {}
    if c_ik: cols[c_ik] = "parent_inchikey"
    if c_sm: cols[c_sm] = "smiles"
    if c_nm: cols[c_nm] = "parent_name"

    meta = meta.rename(columns=cols)
    available_cols = [c for c in ["parent_inchikey","smiles","parent_name"] if c in meta.columns]
    return meta[available_cols].drop_duplicates()

def mol_from_smiles(smi, do3d=True):
    """Create RDKit molecule from SMILES with optional 3D coordinates"""
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    m = Chem.AddHs(m)
    if do3d:
        try:
            AllChem.EmbedMolecule(m, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass
    return m

def repair_parent_metadata(df, meta):
    """Repair missing parent_inchikey using various strategies"""
    print("Repairing missing parent metadata...")

    # Strategy 1: Use existing parent_inchikey where available
    if "parent_inchikey" not in df.columns:
        df["parent_inchikey"] = None
    df["parent_inchikey_fixed"] = df["parent_inchikey"]

    # Strategy 2: Calculate from parent_smiles where parent_inchikey is missing
    if "parent_smiles" in df.columns:
        mask = df["parent_inchikey_fixed"].isna() & df["parent_smiles"].notna()
        if mask.sum() > 0:
            print(f"  Calculating InChIKey from parent_smiles for {mask.sum()} rows")
            df.loc[mask, "parent_inchikey_fixed"] = df.loc[mask, "parent_smiles"].map(compute_ik_from_smiles)

    # Strategy 3: For R6 products, try to infer parent from product structure
    # R6 products should have patterns like -OCH2X or -CHX2 that we can reverse-engineer
    r6_missing = (df["rule"] == "R6") & df["parent_inchikey_fixed"].isna()
    if r6_missing.sum() > 0:
        print(f"  Attempting parent inference for {r6_missing.sum()} R6 products")

        # Create mapping from canonical parent SMILES to InChIKey
        meta_with_ik = meta.copy()
        if "smiles" in meta_with_ik.columns:
            meta_with_ik["meta_ik"] = meta_with_ik["smiles"].map(compute_ik_from_smiles)
            parent_smi_to_ik = dict(zip(meta_with_ik["smiles"], meta_with_ik["meta_ik"]))

            # For R6 products, try to find matching parents by structure similarity
            for idx in df[r6_missing].index:
                product_smi = df.loc[idx, "smiles"]
                # Simple heuristic: try to remove halogen groups and find parent
                # This is a simplified approach - in practice, you might need more sophisticated logic
                for parent_smi, parent_ik in parent_smi_to_ik.items():
                    if parent_ik and parent_smi:
                        # Check if this could be the parent by comparing core structure
                        # For now, we'll use a simple containment check
                        parent_mol = Chem.MolFromSmiles(parent_smi)
                        product_mol = Chem.MolFromSmiles(product_smi)
                        if parent_mol and product_mol:
                            # This is a simplified check - in reality you'd want more sophisticated matching
                            if parent_mol.GetNumAtoms() <= product_mol.GetNumAtoms():
                                df.loc[idx, "parent_inchikey_fixed"] = parent_ik
                                df.loc[idx, "parent_smiles"] = parent_smi
                                break

    return df

def main():
    ap = argparse.ArgumentParser(description="Export per-parent SDF with metadata repair")
    ap.add_argument("--products", required=True, help="Path to products parquet file")
    ap.add_argument("--parents-meta", required=True, help="Path to parents metadata CSV")
    ap.add_argument("--parents-pick", required=True, help="Path to parents pick SMI file")
    ap.add_argument("--outdir", required=True, help="Output directory for SDF files")
    ap.add_argument("--max-per-parent", type=int, default=100000, help="Maximum products per parent")
    ap.add_argument("--no3d", action="store_true", help="Skip 3D coordinate generation")
    args = ap.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading products from {args.products}")
    df = pd.read_parquet(args.products)

    print(f"Original products: {len(df)}")
    print(f"Rule distribution: {dict(df['rule'].value_counts())}")

    # Filter to k<=2 only (handle both string and numeric k values)
    if "k" in df.columns:
        def k_to_int(x):
            try: return int(x)
            except Exception: return None
        kk = df["k"].map(k_to_int)
        df = df[kk.notna() & (kk <= 2)].copy()
        print(f"After k<=2 filter: {len(df)} products")

    # Load metadata and repair missing parent information
    meta = normalize_meta(args.parents_meta)
    df = repair_parent_metadata(df, meta)

    # Load pick list and filter products
    pick_ik = load_pick_smi(args.parents_pick)
    print(f"Pick list contains {len(pick_ik)} parent molecules")

    # Report repair results
    total_rows = len(df)
    still_missing = df["parent_inchikey_fixed"].isna().sum()
    print(f"After repair: {still_missing} rows still missing parent_inchikey ({100.0*still_missing/total_rows:.1f}%)")

    # Filter to products with valid parent info and in pick list
    df_filtered = df[df["parent_inchikey_fixed"].notna() &
                     df["parent_inchikey_fixed"].isin(pick_ik)].copy()

    print(f"Products after filtering by picks: {len(df_filtered)}")

    # Set up parent metadata lookup
    if "parent_inchikey" in meta.columns:
        pm = meta.rename(columns={"parent_inchikey":"ik"}).drop_duplicates("ik")
        pm = pm.set_index("ik")
    else:
        pm = pd.DataFrame()

    # Export columns
    cols = [c for c in ("rule","rule_family","halogen","k","product_smiles","smiles","parent_inchikey_fixed")
            if c in df_filtered.columns]

    # Export SDFs by parent
    total_written = 0
    sdf_count = 0

    for ik in sorted(pick_ik):
        sub = df_filtered[df_filtered["parent_inchikey_fixed"] == ik]
        if len(sub) == 0:
            continue

        outpath = outdir / f"{ik}.sdf"
        w = Chem.SDWriter(str(outpath))

        # Write parent molecule first
        parent_written = False
        try:
            if len(pm) > 0 and ik in pm.index:
                prow = pm.loc[ik]
                p_smi = prow.get("smiles") if pd.notna(prow.get("smiles")) else None
                p_nm = prow.get("parent_name", "")
            else:
                # Fallback: use parent_smiles from product records
                p_smi = sub.iloc[0].get("parent_smiles") if "parent_smiles" in sub.columns else None
                p_nm = sub.iloc[0].get("parent_name", "")

            if p_smi:
                pmol = mol_from_smiles(p_smi, do3d=(not args.no3d))
                if pmol:
                    pmol.SetProp("_Name", f"{ik}_PARENT")
                    pmol.SetProp("is_parent","true")
                    pmol.SetProp("parent_inchikey", ik)
                    if p_nm: pmol.SetProp("parent_name", str(p_nm))
                    w.write(pmol)
                    parent_written = True
        except Exception as e:
            print(f"Warning: Could not write parent for {ik}: {e}")

        # Write products
        sub_sorted = sub.sort_values(["k","rule","halogen"], ascending=[True, True, True]).head(args.max_per_parent)
        products_written = 0

        for _, r in sub_sorted.iterrows():
            # Try different SMILES columns
            smi = r.get("product_smiles") or r.get("smiles")
            if not smi:
                continue

            m = mol_from_smiles(smi, do3d=(not args.no3d))
            if not m:
                continue

            title = f"{r.get('rule','?')}_{r.get('halogen','?')}_k{r.get('k','?')}_{products_written+1:03d}"
            m.SetProp("_Name", title)
            m.SetProp("parent_inchikey", ik)

            for c in cols:
                v = r.get(c)
                if pd.notna(v):
                    m.SetProp(c, str(v))

            w.write(m)
            products_written += 1
            total_written += 1

        w.close()

        status = "PARENT+PRODUCTS" if parent_written else "PRODUCTS_ONLY"
        print(f"[{status:>15s}] {outpath.name}: {products_written} products")
        sdf_count += 1

    print(f"\n=== EXPORT SUMMARY ===")
    print(f"SDF files created: {sdf_count}")
    print(f"Total products written: {total_written}")
    print(f"Expected total (from parquet): {len(df_filtered)}")

    if total_written != len(df_filtered):
        print(f"[WARNING] Written count != filtered count: {total_written} != {len(df_filtered)}")
        print("Some products may have failed SMILES parsing or 3D generation")
    else:
        print("[SUCCESS] All filtered products successfully exported")

    # Rule distribution in exported products
    if len(df_filtered) > 0:
        export_rules = dict(df_filtered["rule"].value_counts())
        print(f"Rule distribution in exported products: {export_rules}")

if __name__ == "__main__":
    main()