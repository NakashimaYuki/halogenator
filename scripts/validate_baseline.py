# -*- coding: utf-8 -*-
"""Validate 10-flavonoids baseline file for SMILES parsing."""

from rdkit import Chem
import sys
import os

def main():
    baseline_file = "examples/input/parents_flavonoids_10.smi"
    
    if not os.path.exists(baseline_file):
        print(f"[ERROR] Baseline file not found: {baseline_file}")
        sys.exit(1)
    
    ok, total = 0, 0
    with open(baseline_file, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # Check TAB separation
            parts = line.split("\t")
            if len(parts) != 2:
                print(f"[ERROR] Line {ln}: not TAB-separated -> {line}")
                sys.exit(1)
            
            smi, name = parts
            if not smi.strip() or not name.strip():
                print(f"[ERROR] Line {ln}: empty SMILES or name -> {line}")
                sys.exit(1)
            
            # Check RDKit parsing
            mol = Chem.MolFromSmiles(smi)
            total += 1
            if mol is None:
                print(f"[ERROR] Line {ln}: RDKit failed to parse SMILES for '{name}': {smi}")
                sys.exit(1)
            ok += 1

    print(f"[OK] Validated {ok}/{total} baseline molecules successfully.")
    if total != 10:
        print(f"[WARNING] Expected 10 molecules, found {total}")
        sys.exit(1)

if __name__ == "__main__":
    main()