# -*- coding: utf-8 -*-
"""
Append target parent molecules (Naringenin + 8-Prenylnaringenin) to parents_pick.smi
Since these specific molecules are not in parents_meta.csv, we'll add them directly
with their known SMILES to enhance R2/R6 testing coverage
"""
import pandas as pd
import pathlib
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchiKey

def main():
    pick_smi = pathlib.Path("out/etcm2000/parents/parents_pick.smi")

    # Target molecules with their canonical SMILES
    # These are specifically chosen for R2/R6 testing
    target_molecules = [
        {
            "smiles": "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12",
            "name": "Naringenin",
            "type": "flavanone"
        },
        {
            "smiles": "COc1cc(O)c2c(=O)cc(-c3ccc(O)cc3)oc2c1CC=C(C)C",
            "name": "8-Prenylnaringenin",
            "type": "prenylated_flavonoid"
        }
    ]

    print("Target molecules for R2/R6 testing enhancement:")
    for mol in target_molecules:
        print(f"  {mol['name']} ({mol['type']})")
        print(f"    SMILES: {mol['smiles']}")

        # Validate SMILES
        rdkit_mol = Chem.MolFromSmiles(mol['smiles'])
        if rdkit_mol is None:
            raise ValueError(f"Invalid SMILES for {mol['name']}: {mol['smiles']}")
        print(f"    Valid: Yes, MW = {Descriptors.MolWt(rdkit_mol):.1f}")
        print()

    # Read existing InChIKeys from pick.smi for deduplication
    existing_ik = set()
    if pick_smi.exists():
        try:
            with pick_smi.open("r", encoding="utf-8") as f:
                content = f.read()
        except UnicodeDecodeError:
            with pick_smi.open("r", encoding="ascii", errors="ignore") as f:
                content = f.read()

        for line in content.splitlines():
            smi = line.strip().split()[0] if line.strip() else ""
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    existing_ik.add(MolToInchiKey(mol))

    # Check and prepare molecules to add
    to_add = []
    for mol_data in target_molecules:
        mol = Chem.MolFromSmiles(mol_data["smiles"])
        if mol is None:
            raise ValueError(f"SMILES parse failed for: {mol_data['name']}")

        ik = MolToInchiKey(mol)
        if ik in existing_ik:
            print(f"[SKIP] Already in pick.smi: {mol_data['name']}")
            continue

        # Format for .smi file
        to_add.append(f"{mol_data['smiles']}\t{mol_data['name']}")

    # Append to pick.smi
    pick_smi.parent.mkdir(parents=True, exist_ok=True)
    if to_add:
        with pick_smi.open("a", encoding="utf-8") as f:
            for line in to_add:
                f.write(line + "\n")

    print(f"\nAdded {len(to_add)} new records to {pick_smi}")
    if to_add:
        for mol_data in target_molecules:
            if any(mol_data['name'] in line for line in to_add):
                print(f"  + {mol_data['name']} ({mol_data['type']})")

    # Show final count
    if pick_smi.exists():
        with pick_smi.open("r", encoding="utf-8") as f:
            total_lines = len([l for l in f.readlines() if l.strip()])
        print(f"\nTotal parents in pick.smi: {total_lines}")
        print("Enhanced dataset ready for R2/R6 testing with flavanone and prenylated structures")

if __name__ == "__main__":
    main()