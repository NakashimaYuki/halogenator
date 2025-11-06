# -*- coding: ascii -*-
"""
Export per-parent SDF (v5) - Complete k<=2 genealogy tracking
- Traces k=2 products back to their root parents through k=1 intermediates
- Ensures all k<=2 products (both direct k=1 and k=1->k=2 chains) are included
- Includes parent molecule (first) + all k=1 and k=2 descendants
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

def build_genealogy_mapping(df):
    """Build complete genealogy: root_parent -> k=1 -> k=2 products"""
    print("Building product genealogy mapping...")

    # Step 1: Map k=1 products to their root parents
    k1_products = df[df['k'] == 1].copy()
    k1_to_root = {}

    for _, row in k1_products.iterrows():
        k1_inchikey = row['inchikey']
        root_parent_ik = row['parent_inchikey']
        k1_to_root[k1_inchikey] = root_parent_ik

    print(f"Mapped {len(k1_to_root)} k=1 products to root parents")

    # Step 2: Map k=2 products to their root parents via k=1 intermediates
    k2_products = df[df['k'] == 2].copy()
    k2_root_assignments = {}

    for idx, row in k2_products.iterrows():
        k1_parent_ik = row['parent_inchikey']  # This points to k=1 product

        if k1_parent_ik in k1_to_root:
            root_parent_ik = k1_to_root[k1_parent_ik]
            k2_root_assignments[idx] = root_parent_ik
        else:
            print(f"Warning: k=2 product {idx} has orphaned parent {k1_parent_ik}")

    print(f"Traced {len(k2_root_assignments)} k=2 products back to root parents")

    return k1_to_root, k2_root_assignments

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
        masked_atoms = get_sugar_mask(parent_mol, mode='heuristic')
        r6_sites = enumerate_methyl_sites(parent_mol, masked_atoms,
                                         allow_on_methoxy=True,
                                         allow_allylic_methyl=True)

        for site in r6_sites:
            for halogen in ['F', 'Cl']:
                product_mol = apply_methyl_step(parent_mol, site, halogen)
                if product_mol:
                    product_smi = Chem.MolToSmiles(product_mol)
                    expected_products.append(product_smi)

    except Exception as e:
        print(f"Warning: Could not compute R6 products for {parent_ik}: {e}")

    return expected_products

def assign_r6_products_to_roots(r6_products, parents):
    """Assign orphaned R6 products to their root parents using structural matching"""
    print(f"Assigning {len(r6_products)} orphaned R6 products to root parents...")

    # Generate expected R6 products for each root parent
    expected_to_parent = {}
    for parent_ik, parent_data in parents.items():
        parent_mol = parent_data['mol']
        expected_products = compute_r6_products_for_parent(parent_mol, parent_ik)

        for expected_smi in expected_products:
            try:
                mol = Chem.MolFromSmiles(expected_smi)
                canonical_smi = Chem.MolToSmiles(mol) if mol else expected_smi
                expected_to_parent[canonical_smi] = parent_ik
            except:
                continue

    print(f"Generated {len(expected_to_parent)} expected R6 products")

    # Match actual R6 products
    r6_assignments = {}
    for idx, r6_product in r6_products.iterrows():
        product_smi = r6_product['smiles']

        # Direct match
        if product_smi in expected_to_parent:
            r6_assignments[idx] = expected_to_parent[product_smi]
        else:
            # Try canonical form
            try:
                mol = Chem.MolFromSmiles(product_smi)
                if mol:
                    canonical_smi = Chem.MolToSmiles(mol)
                    if canonical_smi in expected_to_parent:
                        r6_assignments[idx] = expected_to_parent[canonical_smi]
            except:
                continue

    print(f"Successfully assigned {len(r6_assignments)} R6 products to root parents")
    return r6_assignments

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
    ap = argparse.ArgumentParser(description="Export per-parent SDF with complete k<=2 genealogy")
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
    print(f"k-value distribution: {dict(df['k'].value_counts())}")

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

    # Build complete genealogy mapping
    k1_to_root, k2_root_assignments = build_genealogy_mapping(df)

    # Separate products by type and k-value
    k1_products = df[df['k'] == 1].copy()
    k2_products = df[df['k'] == 2].copy()
    r6_products = df[df['rule'] == 'R6'].copy()

    print(f"k=1 products: {len(k1_products)}")
    print(f"k=2 products: {len(k2_products)}")
    print(f"R6 products (all k): {len(r6_products)}")

    # Get k=1 products with valid root parents
    k1_with_root = k1_products[k1_products['parent_inchikey'].notna() &
                              k1_products['parent_inchikey'].isin(parents.keys())].copy()
    print(f"k=1 products with valid root parents: {len(k1_with_root)}")

    # Assign k=2 products to root parents
    for idx, root_parent_ik in k2_root_assignments.items():
        if root_parent_ik in parents:
            k2_products.loc[idx, 'root_parent_inchikey'] = root_parent_ik

    k2_with_root = k2_products[k2_products.get('root_parent_inchikey', pd.Series()).notna() &
                              k2_products['root_parent_inchikey'].isin(parents.keys())].copy()
    print(f"k=2 products with valid root parents: {len(k2_with_root)}")

    # Handle orphaned R6 products
    r6_orphaned = r6_products[r6_products['parent_inchikey'].isna()]
    r6_assignments = assign_r6_products_to_roots(r6_orphaned, parents)

    for idx, root_parent_ik in r6_assignments.items():
        r6_products.loc[idx, 'root_parent_inchikey'] = root_parent_ik

    r6_with_root = r6_products[r6_products.get('root_parent_inchikey', pd.Series()).notna() &
                              r6_products['root_parent_inchikey'].isin(parents.keys())].copy()
    print(f"R6 products with valid root parents: {len(r6_with_root)}")

    # Combine all products with proper root parent assignments
    # For k=1, use existing parent_inchikey as root_parent_inchikey
    k1_with_root['root_parent_inchikey'] = k1_with_root['parent_inchikey']

    all_products = pd.concat([
        k1_with_root[['smiles', 'inchikey', 'rule', 'rule_family', 'halogen', 'k', 'root_parent_inchikey']],
        k2_with_root[['smiles', 'inchikey', 'rule', 'rule_family', 'halogen', 'k', 'root_parent_inchikey']],
        r6_with_root[['smiles', 'inchikey', 'rule', 'rule_family', 'halogen', 'k', 'root_parent_inchikey']]
    ], ignore_index=True)

    print(f"Total products with root parent assignments: {len(all_products)}")
    print(f"Final rule distribution: {dict(all_products['rule'].value_counts())}")
    print(f"Final k-value distribution: {dict(all_products['k'].value_counts())}")

    # Export SDFs by root parent
    total_written = 0
    sdf_count = 0

    for parent_ik, parent_data in sorted(parents.items()):
        sub = all_products[all_products['root_parent_inchikey'] == parent_ik]
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

        # Write products sorted by k, then rule, then halogen
        sub_sorted = sub.sort_values(["k", "rule", "halogen"], ascending=[True, True, True]).head(args.max_per_parent)
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
            m.SetProp("root_parent_inchikey", parent_ik)
            m.SetProp("rule", str(r.get('rule', '')))
            m.SetProp("halogen", str(r.get('halogen', '')))
            m.SetProp("k", str(r.get('k', '')))

            w.write(m)
            products_written += 1
            total_written += 1

        w.close()

        # Show detailed breakdown for this parent
        parent_rules = dict(sub["rule"].value_counts())
        parent_k_values = dict(sub["k"].value_counts())
        print(f"[COMPLETE] {outpath.name}: {products_written} products")
        print(f"  Rules: {parent_rules}")
        print(f"  k-values: {parent_k_values}")
        sdf_count += 1

    print(f"\n=== FINAL EXPORT SUMMARY ===")
    print(f"SDF files created: {sdf_count}")
    print(f"Total products written: {total_written}")
    print(f"Products by k-value: {dict(all_products['k'].value_counts())}")
    print(f"Products by rule: {dict(all_products['rule'].value_counts())}")

if __name__ == "__main__":
    main()