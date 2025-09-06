# -*- coding: ascii -*-
"""Input/output utilities."""

import os
import time
from typing import List, Tuple, Dict, Any
import urllib.request
import urllib.parse
import json

from rdkit import Chem
from rdkit.Chem import SanitizeMol
import pandas as pd


def load_names(path: str) -> List[str]:
    """Load compound names from file, removing duplicates and empty lines."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Names file not found: {path}")
    
    names = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            name = line.strip()
            if name and name not in names:
                names.append(name)
    return names


def resolve_names_to_smiles(names: List[str]) -> List[Tuple[str, str]]:
    """Resolve compound names to SMILES using PubChem PUG REST API."""
    results = []
    
    for name in names:
        print(f"Resolving: {name}")
        try:
            # PubChem PUG REST API to get canonical SMILES
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/property/CanonicalSMILES/JSON"
            with urllib.request.urlopen(url, timeout=30) as response:
                data = json.loads(response.read().decode('utf-8'))
                smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                results.append((smiles, name))
                print(f"  -> {smiles}")
                time.sleep(0.2)  # Be nice to PubChem
        except Exception as e:
            print(f"  -> Failed: {e}")
            continue
    
    if not results:
        raise RuntimeError(
            "No compounds resolved. Please provide parents.smi manually with format: "
            "SMILES<TAB>Name"
        )
    
    return results


def write_smi(pairs: List[Tuple[str, str]], path: str) -> None:
    """Write SMILES file in format: SMILES<TAB>Name."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        for smiles, name in pairs:
            f.write(f"{smiles}\t{name}\n")


def read_smi(path: str) -> List[Tuple[str, str]]:
    """Read SMILES file with robust separator handling: SMILES<TAB>Name or SMILES<2+spaces>Name."""
    pairs = []
    if not os.path.exists(path):
        return pairs
    
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if '\t' in line:
                # Tab separated - preferred format
                parts = line.split('\t', 1)
                if len(parts) == 2:
                    pairs.append((parts[0], parts[1]))
            else:
                # Split by multiple spaces (2+)
                import re
                parts = re.split(r'  +', line)  # 2 or more spaces
                if len(parts) >= 2:
                    # SMILES should be the first token, name is the last token
                    smiles = parts[0]       # SMILES should be the first token
                    name = parts[-1]        # Name is the last token
                    pairs.append((smiles, name))
    return pairs


def write_sdf_with_props(mols: List[Chem.Mol], path: str) -> None:
    """Write SDF file with properties."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    writer = Chem.SDWriter(path)
    try:
        for mol in mols:
            if mol is not None:
                writer.write(mol)
    finally:
        writer.close()


def write_table(records: List[Dict[str, Any]], path: str) -> None:
    """Write table file (parquet/csv)."""
    if not records:
        return
    
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df = pd.DataFrame(records)
    
    if path.endswith('.parquet'):
        df.to_parquet(path, index=False)
    elif path.endswith('.csv'):
        df.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported file format: {path}")


def read_table(path: str) -> List[Dict[str, Any]]:
    """Read table file (parquet/csv)."""
    if not os.path.exists(path):
        return []
    
    if path.endswith('.parquet'):
        df = pd.read_parquet(path)
    elif path.endswith('.csv'):
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported file format: {path}")
    
    return df.to_dict('records')