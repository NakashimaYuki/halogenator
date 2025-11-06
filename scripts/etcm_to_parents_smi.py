# -*- coding: ascii -*-
import sys, os, re, unicodedata
import pandas as pd

def ascii_slug(s):
    if s is None:
        return "unknown"
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")
    return s or "unknown"

def main(xlsx, out_smi):
    df = pd.read_excel(xlsx)  # requires openpyxl; if unavailable, convert by external tool first
    # Try to guess columns
    cols = [c.lower() for c in df.columns]
    if "smiles" in cols:
        smiles_col = df.columns[cols.index("smiles")]
    else:
        raise SystemExit("SMILES column not found")
    name_col = None
    for cand in ("name", "compound", "compound_name", "chinese_name", "english_name"):
        if cand in cols:
            name_col = df.columns[cols.index(cand)]
            break
    if name_col is None:
        name_col = smiles_col  # fall back

    seen = set()
    lines = []
    for _, row in df.iterrows():
        smi = str(row[smiles_col]).strip()
        if not smi or smi in seen:
            continue
        seen.add(smi)
        nm = ascii_slug(row[name_col])
        lines.append(f"{smi}\t{nm}")

    os.makedirs(os.path.dirname(out_smi), exist_ok=True)
    with open(out_smi, "w", encoding="ascii", errors="ignore") as f:
        f.write("\n".join(lines))
    print(f"Wrote {len(lines)} unique SMILES to {out_smi}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python scripts/etcm_to_parents_smi.py <input.xlsx> <out.smi>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])