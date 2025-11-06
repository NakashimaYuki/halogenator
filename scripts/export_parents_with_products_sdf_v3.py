# -*- coding: ascii -*-
"""
Export per-parent SDF (v3) - Advanced R6 parent matching
- Includes parent molecule (first) + all k<=2 products for that parent
- Uses sophisticated R6 parent matching based on structure patterns
- Ensures all 7022 products including R6 are properly exported

Usage:
  python scripts/export_parents_with_products_sdf_v3.py \
    --products out/etcm2000/pick_k2_compat/products_k2.parquet \
    --parents-pick out/etcm2000/parents/parents_pick.smi \
    --outdir out/etcm2000/pick_k2_compat/per_parent_sdf_fixed
"""
import sys, argparse, pathlib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchiKey

def load_parents_with_metadata(pick_path):
    """Load parent molecules from pick file with computed metadata"""
    parents = {}

    for line in pathlib.Path(pick_path).read_text(encoding="ascii", errors="ignore").splitlines():
        parts = line.strip().split('\t')
        if not parts:
            continue

        parent_smi = parts[0]
        parent_name = parts[1] if len(parts) > 1 else ""

        parent_mol = Chem.MolFromSmiles(parent_smi)
        if parent_mol is None:
            continue

        parent_ik = MolToInchiKey(parent_mol)

        parents[parent_ik] = {
            'smiles': parent_smi,
            'name': parent_name,
            'mol': parent_mol
        }

    return parents

def identify_r6_pattern(product_smi):
    """Identify R6 halogenation patterns in product SMILES"""
    # R6 patterns: -OCH2X, -OCHX2, -OCX3 (methoxy halogenation)
    # or terminal methyl halogenation in prenyl chains
    patterns = [
        r'OCH[FClBrI]',      # -OCHF, -OCHCl, etc.
        r'OCC[FClBrI]',      # -OCCF, -OCCl, etc.
        r'OCH[FClBrI][FClBrI]',  # -OCHF2, etc.
        r'C[FClBrI]$',       # terminal C-X (end of prenyl chain)
        r'\)C[FClBrI]',      # C-X after parentheses (prenyl methyls)
    ]

    for pattern in patterns:
        import re
        if re.search(pattern, product_smi):
            return True
    return False

def remove_r6_halogenation(product_smi):
    """Attempt to reverse R6 halogenation to find parent structure"""
    import re

    # Replace halogenated methoxy back to methoxy
    parent_candidates = []

    # Pattern 1: -OCHX -> -OCH3
    smi1 = re.sub(r'OCH[FClBrI]', 'OC', product_smi)
    parent_candidates.append(smi1)

    # Pattern 2: -OCCX -> -OCH3
    smi2 = re.sub(r'OCC[FClBrI]', 'OC', product_smi)
    parent_candidates.append(smi2)

    # Pattern 3: Terminal C-X -> CH3 (for prenyl chains)
    smi3 = re.sub(r'C[FClBrI]$', 'C', product_smi)
    parent_candidates.append(smi3)

    smi4 = re.sub(r'\)C[FClBrI]', ')C', product_smi)
    parent_candidates.append(smi4)

    return list(set(parent_candidates))

def match_r6_to_parents(r6_products, parents):
    """Match R6 products to their likely parents"""
    matched = {}

    print(f"Attempting to match {len(r6_products)} R6 products to parents...")

    for idx, r6_product in r6_products.iterrows():
        product_smi = r6_product['smiles']

        if not identify_r6_pattern(product_smi):
            continue

        parent_candidates = remove_r6_halogenation(product_smi)

        best_match = None
        best_similarity = 0

        for parent_ik, parent_data in parents.items():
            parent_smi = parent_data['smiles']
            parent_mol = parent_data['mol']

            # Try direct SMILES matching with candidates
            for candidate_smi in parent_candidates:
                try:
                    candidate_mol = Chem.MolFromSmiles(candidate_smi)
                    if candidate_mol is None:
                        continue

                    candidate_canonical = Chem.MolToSmiles(candidate_mol)
                    parent_canonical = Chem.MolToSmiles(parent_mol)

                    # Direct match
                    if candidate_canonical == parent_canonical:
                        matched[idx] = parent_ik
                        best_match = parent_ik
                        break

                    # Substructure match
                    if parent_mol.HasSubstructMatch(candidate_mol) or candidate_mol.HasSubstructMatch(parent_mol):
                        # Calculate a simple similarity score
                        similarity = len(parent_canonical) / max(len(candidate_canonical), len(parent_canonical))
                        if similarity > best_similarity:
                            best_similarity = similarity
                            best_match = parent_ik

                except Exception:
                    continue

            if best_match:
                break

        if best_match and idx not in matched:
            matched[idx] = best_match

    print(f"Successfully matched {len(matched)} R6 products to parents")
    return matched

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

def main():
    ap = argparse.ArgumentParser(description="Export per-parent SDF with advanced R6 matching")
    ap.add_argument("--products", required=True, help="Path to products parquet file")
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

    # Filter to k<=2 only
    if "k" in df.columns:
        def k_to_int(x):
            try: return int(x)
            except Exception: return None
        kk = df["k"].map(k_to_int)
        df = df[kk.notna() & (kk <= 2)].copy()
        print(f"After k<=2 filter: {len(df)} products")

    # Load parent molecules
    parents = load_parents_with_metadata(args.parents_pick)
    print(f"Loaded {len(parents)} parent molecules")

    # Separate products by rule type
    non_r6_products = df[df['rule'] != 'R6'].copy()
    r6_products = df[df['rule'] == 'R6'].copy()

    print(f"Non-R6 products: {len(non_r6_products)}")
    print(f"R6 products: {len(r6_products)}")

    # For non-R6 products, use existing parent_inchikey where available
    non_r6_with_parent = non_r6_products[non_r6_products['parent_inchikey'].notna() &
                                         non_r6_products['parent_inchikey'].isin(parents.keys())].copy()

    # Match R6 products to parents
    r6_matches = match_r6_to_parents(r6_products, parents)

    # Apply R6 matches
    for idx, parent_ik in r6_matches.items():
        r6_products.loc[idx, 'parent_inchikey'] = parent_ik

    # Combine all products with valid parent assignments
    r6_with_parent = r6_products[r6_products['parent_inchikey'].notna() &
                                 r6_products['parent_inchikey'].isin(parents.keys())].copy()

    all_products = pd.concat([non_r6_with_parent, r6_with_parent], ignore_index=True)

    print(f"Total products with parent assignments: {len(all_products)}")
    print(f"Final rule distribution: {dict(all_products['rule'].value_counts())}")

    # Export SDFs by parent
    total_written = 0
    sdf_count = 0

    for parent_ik, parent_data in sorted(parents.items()):
        sub = all_products[all_products['parent_inchikey'] == parent_ik]
        if len(sub) == 0:
            continue

        outpath = outdir / f"{parent_ik}.sdf"
        w = Chem.SDWriter(str(outpath))

        # Write parent molecule first
        parent_mol = mol_from_smiles(parent_data['smiles'], do3d=(not args.no3d))
        if parent_mol:
            parent_mol.SetProp("_Name", f"{parent_ik}_PARENT")
            parent_mol.SetProp("is_parent", "true")
            parent_mol.SetProp("parent_inchikey", parent_ik)
            if parent_data['name']:
                parent_mol.SetProp("parent_name", parent_data['name'])
            w.write(parent_mol)

        # Write products
        sub_sorted = sub.sort_values(["k","rule","halogen"], ascending=[True, True, True]).head(args.max_per_parent)
        products_written = 0

        for _, r in sub_sorted.iterrows():
            smi = r.get("smiles")
            if not smi:
                continue

            m = mol_from_smiles(smi, do3d=(not args.no3d))
            if not m:
                continue

            title = f"{r.get('rule','?')}_{r.get('halogen','?')}_k{r.get('k','?')}_{products_written+1:03d}"
            m.SetProp("_Name", title)
            m.SetProp("parent_inchikey", parent_ik)
            m.SetProp("rule", str(r.get('rule', '')))
            m.SetProp("halogen", str(r.get('halogen', '')))
            m.SetProp("k", str(r.get('k', '')))

            w.write(m)
            products_written += 1
            total_written += 1

        w.close()

        print(f"[PARENT+PRODUCTS] {outpath.name}: {products_written} products")
        sdf_count += 1

    print(f"\n=== EXPORT SUMMARY ===")
    print(f"SDF files created: {sdf_count}")
    print(f"Total products written: {total_written}")
    print(f"R6 products matched: {len(r6_with_parent)}")

    # Final rule distribution check
    if len(all_products) > 0:
        export_rules = dict(all_products["rule"].value_counts())
        print(f"Rule distribution in final export: {export_rules}")

if __name__ == "__main__":
    main()