# -*- coding: ascii -*-
"""
Hierarchical SDF output writer for halogenator products.

This module implements a structured directory layout for enumeration products:
- k=1: Organized by halogen type (X/k1/F/, X/k1/Cl/, etc.)
- k=2: Organized by first halogen, then by k=1 product, then by second halogen

Each SDF file contains:
- First entry: Parent molecule for that group
- Subsequent entries: All products derived from that parent

Index and summary files provide navigation and statistics.
"""

import os
import json
import csv
import hashlib
from typing import Dict, List, Any, Optional, Set
from collections import defaultdict
import logging

LOG = logging.getLogger(__name__)


def _split_by_k_strict(records: List[Dict[str, Any]]):
    """
    Split records by k value with strict type checking and validation.

    This function ensures:
    1. Uses k_ops (preferred) or k (fallback) as the authoritative k column
    2. Coerces values to int to avoid string '1' != integer 1 issues
    3. Filters strictly by k==1 and k==2
    4. Returns empty lists for missing k values (no fallback mixing)

    Args:
        records: List of product record dicts

    Returns:
        Tuple of (k1_records, k2_records, kcol) where kcol is the column used

    Raises:
        ValueError: If neither k_ops nor k exists in records
    """
    if not records:
        return [], [], 'k_ops'

    # Determine which column to use (prefer k_ops)
    sample = records[0]
    if 'k_ops' in sample:
        kcol = 'k_ops'
    elif 'k' in sample:
        kcol = 'k'
    else:
        raise ValueError("Records missing both 'k_ops' and 'k' columns - cannot determine k-level")

    # Filter records with valid k values and coerce to int
    k1_records = []
    k2_records = []

    for rec in records:
        if kcol not in rec:
            LOG.warning(f"Record missing {kcol} column, skipping: {rec.get('inchikey', 'UNKNOWN')}")
            continue

        k_value = rec[kcol]
        if k_value is None:
            LOG.warning(f"Record has None {kcol} value, skipping: {rec.get('inchikey', 'UNKNOWN')}")
            continue

        # Coerce to int (handles both '1' strings and 1 integers)
        try:
            k_int = int(k_value)
        except (ValueError, TypeError) as e:
            LOG.warning(f"Record has invalid {kcol}={k_value} (cannot convert to int), skipping: {rec.get('inchikey', 'UNKNOWN')}")
            continue

        # Strictly filter by k value
        if k_int == 1:
            k1_records.append(rec)
        elif k_int == 2:
            k2_records.append(rec)
        else:
            LOG.warning(f"Record has unexpected {kcol}={k_int} (not 1 or 2), skipping: {rec.get('inchikey', 'UNKNOWN')}")

    LOG.info(f"[_split_by_k_strict] Split {len(records)} records -> k=1: {len(k1_records)}, k=2: {len(k2_records)} (using column '{kcol}')")

    return k1_records, k2_records, kcol


def _write_sdf_strict(records: List[Dict[str, Any]], k_expected: int, kcol: str,
                      filepath: str, parent_mol=None):
    """
    Write records to SDF file with strict k-level validation.

    This function enforces:
    1. ALL records must have kcol == k_expected (no mixing allowed)
    2. Creates a new SDWriter for each call (no handle reuse)
    3. Closes writer immediately after use (no leakage)
    4. Writes k property to SDF for verification

    Args:
        records: List of product records to write
        k_expected: Expected k value (1 or 2)
        kcol: Column name used for k ('k_ops' or 'k')
        filepath: Output SDF file path
        parent_mol: Optional parent molecule to write as first entry

    Raises:
        AssertionError: If any record has k != k_expected
        ValueError: If kcol missing from records
    """
    from .chem_compat import Chem

    if not records:
        LOG.debug(f"[_write_sdf_strict] No records to write for k={k_expected}: {filepath}")
        return

    # STRONG ASSERTION: All records must have the expected k value
    bad_records = []
    for rec in records:
        if kcol not in rec:
            raise ValueError(f"Record missing {kcol} column in k={k_expected} writer: {rec.get('inchikey', 'UNKNOWN')}")

        try:
            k_value = int(rec[kcol])
        except (ValueError, TypeError):
            raise ValueError(f"Record has non-integer {kcol}={rec[kcol]} in k={k_expected} writer: {rec.get('inchikey', 'UNKNOWN')}")

        if k_value != k_expected:
            bad_records.append((rec.get('inchikey', 'UNKNOWN'), k_value))

    if bad_records:
        bad_summary = ', '.join([f"{ik}(k={kv})" for ik, kv in bad_records[:5]])
        if len(bad_records) > 5:
            bad_summary += f" ... and {len(bad_records)-5} more"
        raise AssertionError(
            f"k={k_expected} SDF writer received {len(bad_records)} records with wrong k values: {bad_summary} -> {filepath}"
        )

    # Create output directory
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # Create NEW writer (per-call isolation, no handle reuse)
    LOG.info(f"[_write_sdf_strict] Writing {len(records)} molecules -> {filepath} (k={k_expected}, using {kcol})")

    writer = Chem.SDWriter(filepath)
    try:
        # Write parent molecule first (if provided)
        if parent_mol:
            try:
                parent_mol.SetProp('parent', 'TRUE')
                writer.write(parent_mol)
            except Exception as e:
                LOG.warning(f"Failed to write parent mol to {filepath}: {e}")

        # Write product molecules with metadata
        for rec in records:
            mol = None

            # Try to get mol from record (if cached)
            if '_rdkit_mol' in rec:
                mol = rec['_rdkit_mol']

            # Otherwise parse from SMILES
            if mol is None:
                smiles = rec.get('smiles', '')
                if smiles:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                    except Exception as e:
                        LOG.warning(f"Failed to parse SMILES '{smiles}': {e}")

            if mol is None:
                LOG.warning(f"Skipping record with no valid molecule: {rec.get('inchikey', 'UNKNOWN')}")
                continue

            # Write REQUIRED k property for verification
            mol.SetProp('k', str(int(rec[kcol])))

            # Write other metadata properties
            for prop_name in ['inchikey', 'rule', 'rule_family', 'halogen', 'sub_rule', 'detection',
                              'substitutions_json', 'macro_label', 'k_ops']:
                if prop_name in rec and rec[prop_name] is not None:
                    mol.SetProp(prop_name, str(rec[prop_name]))

            try:
                writer.write(mol)
            except Exception as e:
                LOG.warning(f"Failed to write mol {rec.get('inchikey', 'UNKNOWN')} to {filepath}: {e}")

    finally:
        # ALWAYS close and delete writer to ensure proper cleanup
        writer.flush()
        writer.close()
        del writer

    LOG.debug(f"[_write_sdf_strict] Successfully wrote k={k_expected} SDF: {filepath}")


def _safe_filename(text: str, max_len: int = 50) -> str:
    """
    Create a safe filename from text (InChIKey, name, etc.).

    Args:
        text: Input text
        max_len: Maximum length for filename component

    Returns:
        Safe filename string
    """
    if not text:
        return "unknown"

    # Replace problematic characters
    safe = text.replace('/', '_').replace('\\', '_').replace(':', '_')
    safe = safe.replace(' ', '_').replace('"', '').replace("'", '')

    # Truncate and add hash if too long
    if len(safe) > max_len:
        hash_suffix = hashlib.md5(text.encode()).hexdigest()[:8]
        safe = safe[:max_len-9] + '_' + hash_suffix

    return safe


def _group_by_k(records: List[Dict[str, Any]]) -> Dict[int, List[Dict]]:
    """Group records by k value (substitution depth)."""
    by_k = defaultdict(list)
    for rec in records:
        k = rec.get('k', rec.get('depth', 1))
        by_k[int(k)].append(rec)
    return dict(by_k)


def _extract_history_halogens(record: Dict[str, Any]) -> List[str]:
    """
    Extract halogen sequence from substitutions_json or history.

    Returns:
        List of halogens in order of application (e.g., ['F', 'Cl'])
    """
    # Try substitutions_json first (new format)
    if 'substitutions_json' in record:
        try:
            import json as json_lib
            history = json_lib.loads(record['substitutions_json'])
            return [step.get('halogen', 'X') for step in history]
        except Exception:
            pass

    # Fallback: if k=1, use record's halogen field
    k = record.get('k', record.get('depth', 1))
    if k == 1:
        return [record.get('halogen', 'X')]
    elif k == 2:
        # For k=2 without history, try to infer from parent tracking
        # This is a best-effort fallback
        return [record.get('halogen', 'X')] * 2

    return []


def write_group_sdf(filepath: str, parent_mol, product_mols: List, product_records: List[Dict]):
    """
    Write an SDF file with parent as first entry, followed by products.

    Args:
        filepath: Output SDF file path
        parent_mol: RDKit molecule object for parent (or None)
        product_mols: List of RDKit molecule objects for products
        product_records: List of product record dicts (for metadata)
    """
    from .chem_compat import Chem

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, 'w') as f:
        # Write parent molecule first
        if parent_mol:
            try:
                sdf_block = Chem.MolToMolBlock(parent_mol)
                f.write(sdf_block)
                f.write("> <parent>\nTRUE\n\n")
                f.write("$$$$\n")
            except Exception as e:
                LOG.warning(f"Failed to write parent mol to {filepath}: {e}")

        # Write product molecules
        for mol, rec in zip(product_mols, product_records):
            if mol is None:
                continue

            try:
                sdf_block = Chem.MolToMolBlock(mol)
                f.write(sdf_block)

                # Add metadata fields
                if 'inchikey' in rec:
                    f.write(f"> <inchikey>\n{rec['inchikey']}\n\n")
                if 'k' in rec:
                    f.write(f"> <k>\n{rec['k']}\n\n")
                if 'rule' in rec:
                    f.write(f"> <rule>\n{rec['rule']}\n\n")
                if 'halogen' in rec:
                    f.write(f"> <halogen>\n{rec['halogen']}\n\n")
                if 'substitutions_json' in rec:
                    f.write(f"> <substitutions_json>\n{rec['substitutions_json']}\n\n")

                f.write("$$$$\n")
            except Exception as e:
                LOG.warning(f"Failed to write product mol to {filepath}: {e}")


def write_hierarchical_k1(parent_mol, k1_records: List[Dict], outdir: str,
                          parent_name: str, halogens_order: List[str], kcol: str = 'k_ops'):
    """
    Write k=1 products in hierarchical structure with strict k-level validation.

    Directory structure:
        outdir/<parent_name>/k1/F/<parent_name>_F.sdf
        outdir/<parent_name>/k1/Cl/<parent_name>_Cl.sdf
        ...

    Args:
        parent_mol: RDKit molecule for parent
        k1_records: List of k=1 product records (must all have k==1)
        outdir: Base output directory
        parent_name: Parent molecule name
        halogens_order: Order for halogen subdirectories
        kcol: Column name for k value ('k_ops' or 'k')

    Returns:
        List of file metadata dicts for index

    Raises:
        AssertionError: If any k1_records have k != 1
    """
    files_meta = []
    parent_dir = os.path.join(outdir, parent_name, 'k1')

    # Group by halogen (keep k=1 records separate per halogen)
    by_halogen = defaultdict(list)
    for rec in k1_records:
        hal = rec.get('halogen', 'X')
        by_halogen[hal].append(rec)

    # Write one SDF per halogen in specified order using STRICT writer
    for halogen in halogens_order:
        if halogen not in by_halogen:
            continue

        recs = by_halogen[halogen]
        halogen_dir = os.path.join(parent_dir, halogen)
        filename = f"{parent_name}_{halogen}.sdf"
        filepath = os.path.join(halogen_dir, filename)

        # Use strict writer with k=1 validation
        _write_sdf_strict(recs, k_expected=1, kcol=kcol, filepath=filepath, parent_mol=parent_mol)

        # Record metadata
        files_meta.append({
            'path': os.path.relpath(filepath, outdir),
            'k': 1,
            'halogen_first': halogen,
            'halogen_second': None,
            'parent_inchikey': recs[0].get('parent_inchikey', 'UNKNOWN') if recs else 'UNKNOWN',
            'product_count': len(recs)
        })

    return files_meta


def write_hierarchical_k2(parent_mol, k1_records: List[Dict], k2_records: List[Dict],
                          outdir: str, parent_name: str, halogens_order: List[str], kcol: str = 'k_ops'):
    """
    Write k=2 products in hierarchical structure with strict k-level validation.

    Directory structure:
        outdir/<parent_name>/k2/F/<k1_product_key>/<k1_product_key>_F.sdf
        outdir/<parent_name>/k2/F/<k1_product_key>/<k1_product_key>_Cl.sdf
        ...

    Args:
        parent_mol: RDKit molecule for root parent
        k1_records: List of k=1 product records (for parent lookup)
        k2_records: List of k=2 product records (must all have k==2)
        outdir: Base output directory
        parent_name: Root parent molecule name
        halogens_order: Order for halogen subdirectories
        kcol: Column name for k value ('k_ops' or 'k')

    Returns:
        List of file metadata dicts for index

    Raises:
        AssertionError: If any k2_records have k != 2
    """
    from .chem_compat import Chem

    files_meta = []
    parent_dir = os.path.join(outdir, parent_name, 'k2')

    # Build mapping: k1_product_inchikey -> k1_record
    k1_by_inchikey = {}
    for rec in k1_records:
        key = rec.get('inchikey', rec.get('product_inchikey', ''))
        if key:
            k1_by_inchikey[key] = rec

    # Group k=2 products by (first_halogen, parent_inchikey, second_halogen)
    # parent_inchikey here refers to the k=1 product that is the immediate parent
    by_first_hal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for rec in k2_records:
        # Extract halogen sequence
        halogens_seq = _extract_history_halogens(rec)
        if len(halogens_seq) < 2:
            # Fallback: use record's halogen as second, guess first
            h2 = rec.get('halogen', 'X')
            h1 = h2  # Best guess
        else:
            h1, h2 = halogens_seq[0], halogens_seq[1]

        # Get immediate parent (k=1 product)
        parent_key = rec.get('parent_inchikey', 'UNKNOWN')

        by_first_hal[h1][parent_key][h2].append(rec)

    # Write SDFs organized by first halogen -> k1 parent -> second halogen using STRICT writer
    for h1 in halogens_order:
        if h1 not in by_first_hal:
            continue

        for parent_key, by_h2 in by_first_hal[h1].items():
            # Get k=1 parent molecule for this group
            k1_parent_rec = k1_by_inchikey.get(parent_key)
            k1_parent_mol = None
            if k1_parent_rec:
                try:
                    k1_parent_mol = Chem.MolFromSmiles(k1_parent_rec['smiles'])
                except Exception:
                    pass

            # Create safe directory name for k1 parent
            parent_dir_name = _safe_filename(parent_key, max_len=30)
            k1_parent_dir = os.path.join(parent_dir, h1, parent_dir_name)

            # Write one SDF per second halogen using STRICT writer
            for h2 in halogens_order:
                if h2 not in by_h2:
                    continue

                recs = by_h2[h2]
                filename = f"{parent_dir_name}_{h2}.sdf"
                filepath = os.path.join(k1_parent_dir, filename)

                # Use strict writer with k=2 validation
                _write_sdf_strict(recs, k_expected=2, kcol=kcol, filepath=filepath, parent_mol=k1_parent_mol)

                # Record metadata
                files_meta.append({
                    'path': os.path.relpath(filepath, outdir),
                    'k': 2,
                    'halogen_first': h1,
                    'halogen_second': h2,
                    'parent_inchikey': parent_key,
                    'root_parent_inchikey': recs[0].get('root_parent_inchikey', 'UNKNOWN') if recs else 'UNKNOWN',
                    'product_count': len(recs)
                })

    return files_meta


def write_hierarchical_outputs(parent_record: Dict[str, Any],
                               all_records: List[Dict[str, Any]],
                               outdir: str,
                               halogens_order: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Write enumeration products in hierarchical directory structure with strict k-level validation.

    This function ensures:
    1. All k=1 products go ONLY to k1/ directories
    2. All k=2 products go ONLY to k2/ directories
    3. No cross-contamination between k levels
    4. Assertions will fail-fast if k-level mixing is detected

    Args:
        parent_record: Dict with parent molecule info (smiles, inchikey, name)
        all_records: List of all product records (k=1 and k=2)
        outdir: Base output directory
        halogens_order: Order for halogen subdirectories (default: ['F','Cl','Br','I'])

    Returns:
        Summary dict with paths and statistics

    Raises:
        AssertionError: If k-level mixing is detected during write
        ValueError: If records lack k information
    """
    from .chem_compat import Chem

    if halogens_order is None:
        halogens_order = ['F', 'Cl', 'Br', 'I']

    # Extract parent info
    parent_smiles = parent_record.get('smiles', parent_record.get('parent_smiles', ''))
    parent_inchikey = parent_record.get('inchikey', parent_record.get('parent_inchikey', 'UNKNOWN'))
    parent_name = parent_record.get('name', _safe_filename(parent_inchikey, max_len=40))

    # Convert parent SMILES to mol
    parent_mol = None
    try:
        parent_mol = Chem.MolFromSmiles(parent_smiles)
    except Exception as e:
        LOG.warning(f"Failed to parse parent SMILES: {e}")

    # Use STRICT k-level splitting with type coercion and validation
    k1_records, k2_records, kcol = _split_by_k_strict(all_records)

    LOG.info(f"[write_hierarchical_outputs] Parent {parent_name}: k1={len(k1_records)}, k2={len(k2_records)}, using '{kcol}'")

    # Write k=1 products with strict validation
    k1_files = write_hierarchical_k1(parent_mol, k1_records, outdir, parent_name, halogens_order, kcol=kcol)

    # Write k=2 products with strict validation
    k2_files = write_hierarchical_k2(parent_mol, k1_records, k2_records, outdir, parent_name, halogens_order, kcol=kcol)

    # Generate index.json
    index_data = {
        'parent': {
            'name': parent_name,
            'smiles': parent_smiles,
            'inchikey': parent_inchikey
        },
        'files': k1_files + k2_files,
        'summary': {
            'total_k1_products': len(k1_records),
            'total_k2_products': len(k2_records),
            'total_k1_files': len(k1_files),
            'total_k2_files': len(k2_files)
        }
    }

    index_path = os.path.join(outdir, parent_name, 'index.json')
    os.makedirs(os.path.dirname(index_path), exist_ok=True)
    with open(index_path, 'w') as f:
        json.dump(index_data, f, indent=2)

    # Generate k1_summary.csv
    if k1_files:
        k1_csv_path = os.path.join(outdir, parent_name, 'k1_summary.csv')
        with open(k1_csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['path', 'halogen', 'product_count', 'parent_inchikey'])
            writer.writeheader()
            for fmeta in k1_files:
                writer.writerow({
                    'path': fmeta['path'],
                    'halogen': fmeta['halogen_first'],
                    'product_count': fmeta['product_count'],
                    'parent_inchikey': fmeta['parent_inchikey']
                })

    # Generate k2_summary.csv
    if k2_files:
        k2_csv_path = os.path.join(outdir, parent_name, 'k2_summary.csv')
        with open(k2_csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['path', 'halogen_first', 'halogen_second',
                                                   'product_count', 'parent_inchikey', 'root_parent_inchikey'])
            writer.writeheader()
            for fmeta in k2_files:
                writer.writerow({
                    'path': fmeta['path'],
                    'halogen_first': fmeta['halogen_first'],
                    'halogen_second': fmeta['halogen_second'],
                    'product_count': fmeta['product_count'],
                    'parent_inchikey': fmeta['parent_inchikey'],
                    'root_parent_inchikey': fmeta.get('root_parent_inchikey', '')
                })

    return {
        'parent_name': parent_name,
        'index_path': os.path.relpath(index_path, outdir),
        'total_files': len(k1_files) + len(k2_files),
        'total_products': len(k1_records) + len(k2_records)
    }
