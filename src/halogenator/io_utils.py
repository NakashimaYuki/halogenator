# -*- coding: ascii -*-
"""Input/output utilities."""

import os
import time
import re
from typing import List, Tuple, Dict, Any, Iterator
import urllib.request
import urllib.parse
import json
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
                parts = re.split(r'  +', line)  # 2 or more spaces
                if len(parts) >= 2:
                    # SMILES should be the first token, name is the last token
                    smiles = parts[0]       # SMILES should be the first token
                    name = parts[-1]        # Name is the last token
                    pairs.append((smiles, name))
                else:
                    # SMILES-only line - generate default name
                    smiles = line.strip()
                    if smiles:  # Make sure it's not empty
                        name = f"mol_{len(pairs) + 1}"
                        pairs.append((smiles, name))
    return pairs


def write_sdf_with_props(mols: List, path: str) -> None:
    """Write SDF file with properties."""
    try:
        from rdkit import Chem
    except ImportError as e:
        raise ImportError("RDKit is required for SDF file writing") from e
    
    os.makedirs(os.path.dirname(path), exist_ok=True)
    writer = Chem.SDWriter(path)
    try:
        for mol in mols:
            if mol is not None:
                writer.write(mol)
    finally:
        writer.close()


_JSON_FIELDS = ("substitutions", "constraints_violations", "sugar_mask_atoms", "sugar_rings")

# Type table for JSON fields: distinguish list vs dict semantics
# List fields default to [], dict fields default to {}
_JSON_LIST_FIELDS = frozenset({"substitutions", "sugar_mask_atoms", "sugar_rings"})
_JSON_DICT_FIELDS = frozenset({"constraints_violations"})

def _get_json_default(field_name):
    """Get the appropriate default value for a JSON field based on its semantic type."""
    if field_name in _JSON_LIST_FIELDS:
        return []
    elif field_name in _JSON_DICT_FIELDS:
        return {}
    else:
        # Fallback for unknown fields (should not happen if _JSON_FIELDS is complete)
        return {}

def _prepare_records_for_table(records):
    out = []
    for r in records:
        r2 = dict(r)
        # Create parallel *_json string columns, keep typed fields in memory
        for f in _JSON_FIELDS:
            v = r2.get(f, None)
            try:
                default_val = _get_json_default(f)
                r2[f + "_json"] = json.dumps(v if v is not None else default_val)
            except Exception:
                # Fallback on serialization error: use appropriate default
                default_val = _get_json_default(f)
                r2[f + "_json"] = json.dumps(default_val)
        out.append(r2)
    return out

def write_table(records: List[Dict[str, Any]], path: str) -> None:
    """Write table file (parquet/csv)."""
    if not records:
        return
    
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = _prepare_records_for_table(records)
    df = pd.DataFrame(rows)
    # Drop raw typed fields to make parquet robust
    for f in _JSON_FIELDS:
        if f in df.columns:
            df.drop(columns=[f], inplace=True)
    
    if path.endswith('.parquet'):
        df.to_parquet(path, index=False)
    elif path.endswith('.csv'):
        df.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported file format: {path}")


def read_table(path: str) -> List[Dict[str, Any]]:
    """Read table file (parquet/csv) and rebuild typed fields, removing JSON columns from memory."""
    if not os.path.exists(path):
        return []
    
    if path.endswith('.parquet'):
        df = pd.read_parquet(path)
    elif path.endswith('.csv'):
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported file format: {path}")
    
    # Rebuild typed fields from *_json string columns
    for f in _JSON_FIELDS:
        col = f + "_json"
        if col in df.columns:
            try:
                default_val = _get_json_default(f)
                df[f] = df[col].apply(lambda s: json.loads(s) if isinstance(s, str) and len(s) else default_val)
            except Exception:
                default_val = _get_json_default(f)
                df[f] = [default_val for _ in range(len(df))]
            # Remove the *_json column from memory objects
            df.drop(columns=[col], inplace=True)
    
    return df.to_dict('records')


def iter_table_records(path: str, csv_chunksize: int = 200_000) -> Iterator[Dict[str, Any]]:
    """Yield records incrementally from parquet/csv without loading all into memory."""
    if not os.path.exists(path):
        return
        
    if path.endswith(".parquet"):
        try:
            import pyarrow.parquet as pq
            import pyarrow as pa
        except ImportError as e:
            raise RuntimeError("pyarrow is required for parquet streaming") from e
        
        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches():
            tbl = pa.Table.from_batches([batch])
            for row in tbl.to_pylist():
                # Rebuild typed fields from *_json columns like read_table does
                for f in _JSON_FIELDS:
                    col = f + "_json"
                    if col in row and row[col] is not None:
                        try:
                            row[f] = json.loads(row[col])
                        except Exception:
                            row[f] = _get_json_default(f)
                        row.pop(col, None)
                yield row
                
    elif path.endswith(".csv"):
        for chunk in pd.read_csv(path, chunksize=csv_chunksize):
            # Rebuild typed fields from *_json columns
            for f in _JSON_FIELDS:
                col = f + "_json"
                if col in chunk.columns:
                    try:
                        default_val = _get_json_default(f)
                        chunk[f] = chunk[col].apply(lambda s: json.loads(s) if isinstance(s, str) and len(s) else default_val)
                    except Exception:
                        default_val = _get_json_default(f)
                        chunk[f] = [default_val for _ in range(len(chunk))]
                    chunk.drop(columns=[col], inplace=True)
            
            for rec in chunk.to_dict("records"):
                yield rec
    else:
        raise ValueError(f"Unsupported file format: {path}")