# -*- coding: ascii -*-
import sys
import os
import re
import unicodedata
import pandas as pd

def ascii_slug(s):
    if s is None:
        return "unknown"
    s = unicodedata.normalize("NFKD", str(s))
    s = s.encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        return "unknown"
    if len(s) > 80:
        s = s[:80]
    return s

def _norm_cols(cols):
    out = []
    for c in cols:
        c = str(c).strip()
        c_low = c.lower()
        out.append((c, c_low))
    return out

def _pick_col(df, wanted):
    # wanted: list of lower-case candidate names
    cols = _norm_cols(df.columns)
    for (raw, low) in cols:
        if low in wanted:
            return raw
    return None

def load_table(path, force_csv=False):
    ext = os.path.splitext(path)[1].lower()
    if force_csv or ext == ".csv":
        return pd.read_csv(path)
    if ext == ".xlsx":
        try:
            return pd.read_excel(path, engine="openpyxl")
        except Exception as e:
            raise SystemExit(
                "Reading .xlsx requires openpyxl. Please install openpyxl or convert the file to CSV.\n"
                "Example: pip install openpyxl\n"
                "Or: save as CSV and rerun with --csv"
            )
    raise SystemExit("Unsupported input format. Use .xlsx or .csv")

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Convert ETCM flavonoid Excel/CSV to .smi and meta CSV (ASCII-safe)."
    )
    ap.add_argument("-i", "--input", required=True, help="Input .xlsx or .csv")
    ap.add_argument("-o", "--out-smi", default="data/input/etcm_flavonoids.smi", help="Output .smi")
    ap.add_argument("--out-meta", default="data/input/etcm_flavonoids_meta.csv", help="Output meta CSV")
    ap.add_argument("--csv", action="store_true", help="Force read as CSV")
    ap.add_argument("--sample", type=int, default=0, help="Take first N rows after filtering")
    args = ap.parse_args()

    df = load_table(args.input, force_csv=args.csv)

    # Column picking with English and Chinese names using unicode escapes
    id_candidates = {"ingredient id"}
    name_en_candidates = {"ingredient name in english"}
    name_zh_candidates = {"\u6210\u5206\u4e2d\u6587\u540d"}  # Unicode for Chinese name
    smiles_candidates = {"\u6807\u51c6smiles", "smiles"}      # Unicode for standard SMILES
    flag_candidates = {"is_flavonoid"}

    id_col = _pick_col(df, id_candidates)
    name_en_col = _pick_col(df, name_en_candidates)
    name_zh_col = _pick_col(df, name_zh_candidates)
    smiles_col = _pick_col(df, smiles_candidates)
    flag_col = _pick_col(df, flag_candidates)

    if not smiles_col:
        raise SystemExit("Cannot find SMILES column (looked for standard SMILES or SMILES).")

    # Filter flavonoids if column exists
    if flag_col and flag_col in df.columns:
        df = df[df[flag_col] == True]

    # Down-sample if required
    if args.sample and args.sample > 0:
        df = df.head(args.sample)

    # Deduplicate by SMILES
    df = df.dropna(subset=[smiles_col])
    df[smiles_col] = df[smiles_col].astype(str).str.strip()
    df = df[df[smiles_col] != ""]
    df = df.drop_duplicates(subset=[smiles_col], keep="first")

    # Build ascii names
    ascii_names = []
    for idx, row in df.iterrows():
        ingr_id = str(row[id_col]).strip() if id_col else ""
        name_en = str(row[name_en_col]).strip() if name_en_col else ""
        slug = ascii_slug(name_en) if name_en else "unknown"
        if ingr_id:
            name = "{}__{}".format(ascii_slug(ingr_id), slug)
        else:
            name = slug
        ascii_names.append(name)

    df["_ascii_name"] = ascii_names

    # Write .smi
    os.makedirs(os.path.dirname(args.out_smi), exist_ok=True)
    with open(args.out_smi, "w", encoding="ascii", errors="ignore") as f:
        for _, row in df.iterrows():
            f.write("{}\t{}\n".format(row[smiles_col], row["_ascii_name"]))

    # Write meta CSV
    os.makedirs(os.path.dirname(args.out_meta), exist_ok=True)
    keep_cols = []
    for c in [id_col, name_en_col, name_zh_col, smiles_col, flag_col]:
        if c and c not in keep_cols:
            keep_cols.append(c)
    df_meta = df[keep_cols + ["_ascii_name"]].copy()
    df_meta.to_csv(args.out_meta, index=False)

    print("Wrote {} rows".format(len(df)))
    print("SMI:  {}".format(args.out_smi))
    print("META: {}".format(args.out_meta))

if __name__ == "__main__":
    main()