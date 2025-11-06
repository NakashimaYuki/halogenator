# -*- coding: ascii -*-
"""
Export per-parent SDF (v4) - Direct R6 product mapping
- Uses direct structural comparison to map R6 products to parents
- Ensures all products including R6 are properly exported
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

def compute_r6_products_for_parent(parent_mol, parent_ik):
    """Generate expected R6 products for a given parent molecule"""
    import sys
    sys.path.append('src')
    try:
        from halogenator.sites_methyl import enumerate_methyl_sites
        from halogenator.rules_methyl import apply_methyl_step
        from halogenator.sugar_mask import get_sugar_mask
    except ImportError:
        return []

    expected_products = []

    try:
        # Get R6 sites
        masked_atoms = get_sugar_mask(parent_mol, mode='heuristic')
        r6_sites = enumerate_methyl_sites(parent_mol, masked_atoms,
                                         allow_on_methoxy=True,
                                         allow_allylic_methyl=True)

        # Generate R6 products
        for site in r6_sites:
            for halogen in ['F', 'Cl']:  # R6 only supports F and Cl
                product_mol = apply_methyl_step(parent_mol, site, halogen)
                if product_mol:
                    product_smi = Chem.MolToSmiles(product_mol)
                    expected_products.append(product_smi)

    except Exception as e:
        print(f"Warning: Could not compute R6 products for {parent_ik}: {e}")

    return expected_products

def match_r6_products_direct(r6_products, parents):
    """Match R6 products to parents using direct structural comparison"""
    matched = {}

    print(f"Computing expected R6 products for {len(parents)} parents...")

    # Create mapping of expected products to parents
    expected_to_parent = {}
    for parent_ik, parent_data in parents.items():
        parent_mol = parent_data['mol']
        expected_products = compute_r6_products_for_parent(parent_mol, parent_ik)

        for expected_smi in expected_products:
            # Normalize SMILES
            try:
                mol = Chem.MolFromSmiles(expected_smi)
                canonical_smi = Chem.MolToSmiles(mol) if mol else expected_smi
                expected_to_parent[canonical_smi] = parent_ik
            except:
                continue

    print(f"Generated {len(expected_to_parent)} expected R6 products")

    # Match actual products to expected
    matched_count = 0
    for idx, r6_product in r6_products.iterrows():
        product_smi = r6_product['smiles']

        # Direct match
        if product_smi in expected_to_parent:
            matched[idx] = expected_to_parent[product_smi]
            matched_count += 1
        else:
            # Try normalizing the product SMILES
            try:
                mol = Chem.MolFromSmiles(product_smi)
                if mol:
                    canonical_smi = Chem.MolToSmiles(mol)
                    if canonical_smi in expected_to_parent:
                        matched[idx] = expected_to_parent[canonical_smi]
                        matched_count += 1
            except:
                continue

    print(f"Successfully matched {matched_count} R6 products to parents")
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
    ap = argparse.ArgumentParser(description="Export per-parent SDF with direct R6 matching")
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

    print(f"Non-R6 products with valid parents: {len(non_r6_with_parent)}")

    # Match R6 products to parents using direct computation
    r6_matches = match_r6_products_direct(r6_products, parents)

    # Apply R6 matches
    for idx, parent_ik in r6_matches.items():
        r6_products.loc[idx, 'parent_inchikey'] = parent_ik

    # Get R6 products with valid parent assignments
    r6_with_parent = r6_products[r6_products['parent_inchikey'].notna() &
                                 r6_products['parent_inchikey'].isin(parents.keys())].copy()

    print(f"R6 products with valid parents: {len(r6_with_parent)}")

    # Combine all products with valid parent assignments
    all_products = pd.concat([non_r6_with_parent, r6_with_parent], ignore_index=True)

    print(f"Total products with parent assignments: {len(all_products)}")
    print(f"Final rule distribution: {dict(all_products['rule'].value_counts())}")

    # Export SDFs by parent
    total_written = 0
    sdf_count = 0

    for parent_ik, parent_data in sorted(parents.items()):
        sub = all_products[all_products['parent_inchikey'] == parent_ik]
        if len(sub) == 0:
            print(f"No products for parent {parent_ik}")
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

        # Show rule distribution for this parent
        parent_rules = dict(sub["rule"].value_counts())
        print(f"[PARENT+PRODUCTS] {outpath.name}: {products_written} products {parent_rules}")
        sdf_count += 1

    print(f"\n=== EXPORT SUMMARY ===")
    print(f"SDF files created: {sdf_count}")
    print(f"Total products written: {total_written}")
    print(f"R6 products successfully matched and exported: {len(r6_with_parent)}")

    # Final rule distribution check
    if len(all_products) > 0:
        export_rules = dict(all_products["rule"].value_counts())
        print(f"Final rule distribution in export: {export_rules}")

if __name__ == "__main__":
    main()