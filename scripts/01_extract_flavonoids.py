#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flavonoid Extraction from CNPD-ETCM Database

Two-stage robust filtering:
1. Stage I (loose): SMARTS OR morphology OR name keyword → candidates
2. Stage II (strict): SMARTS AND morphology≥2 → final

Uses RDKit standardization: largest fragment, metal disconnect, uncharge, tautomer canonical
"""

import os
import sys
import pandas as pd
import datetime
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize as std
from rdkit.Chem import rdMolDescriptors
import logging
from typing import List, Tuple, Dict, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# SMARTS Library for Flavonoid Classes
# =============================================================================

FLAVONOID_SMARTS = {
    # Flavone / Flavonol (2-phenylchromen-4-one)
    # Multiple patterns to reduce false negatives
    "flavone_v1": "O=C1C=C(Oc2ccccc12)c3ccccc3",
    "flavone_v2": "O=C1C=C(c2ccccc2)Oc2ccccc12",
    "flavone_v3": "[O;D1]=,:[C;R1]1[C;R1]=,:[C;R1]([c,C;R1]2[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]2)[O;R1][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]12",

    # Flavonol (3-hydroxyflavone)
    "flavonol": "O=C1C(O)=C(Oc2ccccc12)c3ccccc3",

    # Isoflavone (3-phenylchromen-4-one, B-ring at C3)
    "isoflavone_v1": "O=C1COc2ccccc2C1=Cc3ccccc3",
    "isoflavone_v2": "O=C1C(=Cc2ccccc2)Oc2ccccc12",
    "isoflavone_v3": "[O;D1]=,:[C;R1]1[C;R1][O;R1][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]2[C;R1]1=,:[C;R0][c,C;R1]3[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]3",

    # Flavanone (2-phenylchroman-4-one, C2-C3 saturated)
    "flavanone_v1": "O=C1CC(Oc2ccccc12)c3ccccc3",
    "flavanone_v2": "[O;D1]=,:[C;R1]1[C;R0][C;R1]([c,C;R1]2[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]2)[O;R1][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]12",

    # Flavanonol (dihydroflavonol, 3-hydroxyflavanone)
    "flavanonol_v1": "O=C1C(O)C(Oc2ccccc12)c3ccccc3",
    "flavanonol_v2": "[O;D1]=,:[C;R1]1[C;R0]([O;D1])[C;R1]([c,C;R1]2[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]2)[O;R1][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]12",

    # Flavan-3-ol (catechin/epicatechin, C4 no ketone)
    "flavan3ol_v1": "OC1Cc2ccccc2OC1c3ccccc3",
    "flavan3ol_v2": "[O;D1][C;R1]1[C;R0][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]2[O;R1][C;R1]1[c,C;R1]3[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]3",

    # Chalcone (1,3-diphenyl-2-propen-1-one, open-chain precursor)
    "chalcone_v1": "O=C(C=Cc1ccccc1)c2ccccc2",
    "chalcone_v2": "[c,C;R1]1[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]1[C;R0](=,:[O;D1])[C;R0]=,:[C;R0][c,C;R1]2[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]2",

    # Aurone (2-benzylidene-benzofuran-3(2H)-one)
    "aurone_v1": "O=C1Oc2ccccc2C1=Cc3ccccc3",
    "aurone_v2": "[O;D1]=,:[C;R1]1[O;R1][c;R1]2[c;R1][c;R1][c;R1][c;R1][c;R1]2[C;R1]1=,:[C;R0][c,C;R1]3[c,C;R1][c,C;R1][c,C;R1][c,C;R1][c,C;R1]3",

    # Anthocyanidin (flavylium cation)
    "anthocyanidin": "[O+]=c1cc(-c2ccccc2)oc2ccccc12",

    # Biflavone (dimer)
    "biflavone": "O=C1C=C(Oc2ccc(cc12)-c3cc4c(O)cc(cc4[o+]c3-c5ccccc5)c6ccccc6)c7ccccc7",
}

# Compile SMARTS patterns
COMPILED_SMARTS = {}
for name, smarts in FLAVONOID_SMARTS.items():
    try:
        patt = Chem.MolFromSmarts(smarts)
        if patt:
            COMPILED_SMARTS[name] = patt
        else:
            logger.warning(f"Failed to compile SMARTS for {name}: {smarts}")
    except Exception as e:
        logger.warning(f"Exception compiling SMARTS {name}: {e}")

logger.info(f"Compiled {len(COMPILED_SMARTS)}/{len(FLAVONOID_SMARTS)} SMARTS patterns")


# Name keywords (case-insensitive, multi-language)
FLAVONOID_KEYWORDS = [
    'flavone', 'flavonol', 'flavanone', 'flavanonol',
    'flavan-3-ol', 'flavan', 'isoflavone', 'isoflavonoid',
    'chalcone', 'aurone', 'anthocyanidin', 'anthocyanin',
    'biflavone', 'biflavonoid', 'prenylflavonoid',
    'flavonolignan', 'chromone', 'chromen',
    # Chinese
    '黄酮', '异黄酮', '二氢黄酮', '黄烷', '查尔酮', '奥罗酮',
    '花色素', '双黄酮', '色原酮'
]


# =============================================================================
# Molecule Standardization
# =============================================================================

def sanitize_molecule(smiles: str) -> Optional[Chem.Mol]:
    """
    Standardize molecule using RDKit best practices:
    1. Parse SMILES
    2. Largest fragment (remove salts/solvents)
    3. Metal disconnect
    4. Uncharge
    5. Tautomer canonicalization (optional but recommended)
    6. Sanitize

    Returns canonical mol or None if invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None

        # Largest fragment
        lfc = std.LargestFragmentChooser()
        mol = lfc.choose(mol)

        # Metal disconnect
        md = std.MetalDisconnector()
        mol = md.Disconnect(mol)

        # Uncharge
        uc = std.Uncharger()
        mol = uc.uncharge(mol)

        # Tautomer canonicalization (choose canonical tautomer)
        try:
            enumerator = std.TautomerEnumerator()
            mol = enumerator.Canonicalize(mol)
        except Exception as e:
            logger.debug(f"Tautomer canonicalization failed: {e}")
            pass

        # Final sanitization
        Chem.SanitizeMol(mol)

        return mol

    except Exception as e:
        logger.debug(f"Sanitization failed for '{smiles[:50]}...': {e}")
        return None


# =============================================================================
# Flavonoid Detection
# =============================================================================

def check_smarts_match(mol: Chem.Mol) -> List[str]:
    """Check which SMARTS patterns match the molecule"""
    matches = []
    for name, patt in COMPILED_SMARTS.items():
        try:
            if mol.HasSubstructMatch(patt):
                matches.append(f"smarts:{name}")
        except Exception as e:
            logger.debug(f"SMARTS match error for {name}: {e}")
    return matches


def check_morphology(mol: Chem.Mol) -> List[str]:
    """
    Morphological heuristics for flavonoid detection:
    1. ≥2 aromatic rings
    2. C6-C3-C6 backbone (shortest path 2-4 atoms between aromatic rings, contains O or C=O)
    3. Benzopyran/chromone ring system (6-member O-heterocycle fused with benzene)
    """
    hints = []

    try:
        # 1. Aromatic ring count
        ri = mol.GetRingInfo()
        aromatic_rings = []
        for ring in ri.AtomRings():
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                aromatic_rings.append(ring)

        if len(aromatic_rings) >= 2:
            hints.append("morph:>=2_aromatic_rings")

        # 2. C6-C3-C6 backbone check
        # Find pairs of aromatic carbon atoms in different rings
        # Check shortest path length and presence of O or C=O
        aromatic_carbons_by_ring = []
        for ring in aromatic_rings:
            ring_aromatic_c = [i for i in ring
                              if mol.GetAtomWithIdx(i).GetSymbol() == 'C'
                              and mol.GetAtomWithIdx(i).GetIsAromatic()]
            if ring_aromatic_c:
                aromatic_carbons_by_ring.append(ring_aromatic_c)

        if len(aromatic_carbons_by_ring) >= 2:
            # Check shortest path between first two aromatic rings
            ring1_atoms = aromatic_carbons_by_ring[0]
            ring2_atoms = aromatic_carbons_by_ring[1]

            min_path_len = float('inf')
            best_path = None

            for a1 in ring1_atoms[:3]:  # Sample first 3 to avoid excessive computation
                for a2 in ring2_atoms[:3]:
                    try:
                        path = Chem.GetShortestPath(mol, a1, a2)
                        if len(path) < min_path_len:
                            min_path_len = len(path)
                            best_path = path
                    except:
                        pass

            # C6-C3-C6 means 3 atoms between rings (path length ~3-5)
            if best_path and 2 <= len(best_path) <= 6:
                path_atoms = [mol.GetAtomWithIdx(i) for i in best_path]
                has_oxygen = any(a.GetSymbol() == 'O' for a in path_atoms)
                has_carbonyl = False
                for atom in path_atoms:
                    if atom.GetSymbol() == 'C':
                        for bond in atom.GetBonds():
                            if (bond.GetBondType() == Chem.BondType.DOUBLE and
                                bond.GetOtherAtom(atom).GetSymbol() == 'O'):
                                has_carbonyl = True
                                break

                if has_oxygen or has_carbonyl:
                    hints.append("morph:C6-C3-C6_with_O")

        # 3. Benzopyran/chromone ring system
        # Look for 6-member ring with O fused to aromatic 6-member ring
        for ring in ri.AtomRings():
            if len(ring) == 6:
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                has_oxygen = any(a.GetSymbol() == 'O' for a in ring_atoms)

                if has_oxygen:
                    # Check if fused with aromatic ring
                    for other_ring in aromatic_rings:
                        if len(set(ring) & set(other_ring)) == 2:  # Fused (share 2 atoms)
                            hints.append("morph:benzopyran_system")
                            break

    except Exception as e:
        logger.debug(f"Morphology check error: {e}")

    return hints


def check_name_keywords(name: str) -> List[str]:
    """Check if compound name contains flavonoid keywords"""
    if not name or not isinstance(name, str):
        return []

    name_lower = name.lower()
    matches = []
    for kw in FLAVONOID_KEYWORDS:
        if kw.lower() in name_lower:
            matches.append(f"name:{kw}")

    return matches


def is_flavonoid_stage1(mol: Chem.Mol, name: str = "") -> Tuple[bool, List[str]]:
    """
    Stage I (loose): SMARTS OR morphology OR name keyword
    Returns (is_candidate, reasons)
    """
    reasons = []

    # SMARTS
    smarts_hits = check_smarts_match(mol)
    reasons.extend(smarts_hits)

    # Morphology
    morph_hits = check_morphology(mol)
    reasons.extend(morph_hits)

    # Name
    name_hits = check_name_keywords(name)
    reasons.extend(name_hits)

    is_candidate = len(reasons) > 0

    return is_candidate, reasons


def is_flavonoid_stage2(mol: Chem.Mol, name: str = "") -> Tuple[bool, List[str]]:
    """
    Stage II (strict): SMARTS AND morphology≥2
    (Or strong SMARTS match with ≥1 morphology)
    Returns (is_flavonoid, reasons)
    """
    reasons = []

    # SMARTS
    smarts_hits = check_smarts_match(mol)
    reasons.extend(smarts_hits)

    # Morphology
    morph_hits = check_morphology(mol)
    reasons.extend(morph_hits)

    # Strict criteria:
    # - At least one SMARTS match AND at least 2 morphology hits
    # OR at least 2 SMARTS matches AND at least 1 morphology hit
    # OR at least 3 morphology hits (strong morphological evidence)
    has_smarts = len(smarts_hits) >= 1
    has_morph_strong = len(morph_hits) >= 2
    has_morph_very_strong = len(morph_hits) >= 3
    has_morph_weak = len(morph_hits) >= 1
    has_smarts_strong = len(smarts_hits) >= 2

    is_flavonoid = (
        (has_smarts and has_morph_strong) or
        (has_smarts_strong and has_morph_weak) or
        has_morph_very_strong
    )

    return is_flavonoid, reasons


# =============================================================================
# Main Processing
# =============================================================================

def load_cnpd_excel(excel_path: str) -> pd.DataFrame:
    """Load CNPD-ETCM Excel file"""
    logger.info(f"Loading Excel file: {excel_path}")

    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel file not found: {excel_path}")

    df = pd.read_excel(excel_path, engine="openpyxl")
    logger.info(f"Loaded {len(df)} rows from Excel")

    # Strip whitespace from column names
    df = df.rename(columns=str.strip)

    # Filter out rows with missing or invalid SMILES
    initial_count = len(df)
    df = df[df.get("Smiles", pd.Series()).notna()].copy()

    # Filter out common invalid SMILES values
    invalid_smiles = ['SDF', 'NA', 'N/A', '', 'NULL', 'null', 'NaN']
    df['Smiles'] = df['Smiles'].astype(str).str.strip()
    df = df[~df['Smiles'].isin(invalid_smiles)].copy()
    df = df[df['Smiles'].str.len() > 2].copy()  # Require at least 3 characters

    filtered_count = initial_count - len(df)
    logger.info(f"Retained {len(df)} rows with valid SMILES (filtered out {filtered_count} invalid/missing)")

    return df


def clean_dataframe_for_parquet(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean DataFrame to ensure compatibility with Parquet format
    Converts datetime objects and other problematic types to strings
    """
    df_clean = df.copy()

    for col in df_clean.columns:
        # Convert datetime objects to strings
        if df_clean[col].dtype == 'object':
            try:
                # Check if column contains datetime objects
                if df_clean[col].apply(lambda x: isinstance(x, (pd.Timestamp, datetime.datetime))).any():
                    df_clean[col] = df_clean[col].astype(str)
            except:
                pass

    return df_clean


def standardize_and_deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize molecules and deduplicate by InChIKey
    """
    logger.info("Starting standardization and deduplication...")

    rows = []
    seen_inchikeys = set()
    failed_count = 0

    for idx, row in df.iterrows():
        smiles = row.get("Smiles", "")

        # Sanitize
        mol = sanitize_molecule(smiles)
        if not mol:
            failed_count += 1
            continue

        # Generate canonical SMILES and InChIKey
        try:
            canonical_smi = Chem.MolToSmiles(mol, canonical=True)
            inchikey = Chem.inchi.MolToInchiKey(mol)

            if not inchikey:
                failed_count += 1
                continue

            # Deduplicate by InChIKey
            if inchikey in seen_inchikeys:
                continue
            seen_inchikeys.add(inchikey)

            # Add standardized data - convert to simple types
            row_dict = {}
            for k, v in row.items():
                # Convert datetime/timestamp to string
                if isinstance(v, (pd.Timestamp, datetime.datetime)):
                    row_dict[k] = str(v)
                else:
                    row_dict[k] = v

            row_dict["Smiles_clean"] = canonical_smi
            row_dict["InChIKey"] = inchikey
            row_dict["MW"] = float(Descriptors.MolWt(mol))
            row_dict["Formula"] = str(rdMolDescriptors.CalcMolFormula(mol))

            rows.append(row_dict)

        except Exception as e:
            logger.debug(f"Failed to process molecule: {e}")
            failed_count += 1
            continue

        if (idx + 1) % 5000 == 0:
            logger.info(f"Processed {idx + 1}/{len(df)} rows, {len(rows)} valid, {failed_count} failed")

    logger.info(f"Standardization complete: {len(rows)} unique molecules, {failed_count} failed")

    return pd.DataFrame(rows)


def extract_flavonoids_stage1(df: pd.DataFrame) -> pd.DataFrame:
    """
    Stage I: Loose filtering (SMARTS OR morphology OR name)
    """
    logger.info("Stage I: Loose flavonoid filtering...")

    candidates = []

    for idx, row in df.iterrows():
        smiles = row["Smiles_clean"]
        name = row.get("Name", "") or row.get("ID", "") or ""

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue

        is_candidate, reasons = is_flavonoid_stage1(mol, name)

        if is_candidate:
            row_dict = row.to_dict()
            row_dict["flav_match_reasons"] = "|".join(reasons)
            candidates.append(row_dict)

        if (idx + 1) % 5000 == 0:
            logger.info(f"Stage I: Processed {idx + 1}/{len(df)} rows, {len(candidates)} candidates")

    logger.info(f"Stage I complete: {len(candidates)} candidates from {len(df)} molecules")

    return pd.DataFrame(candidates)


def extract_flavonoids_stage2(df: pd.DataFrame) -> pd.DataFrame:
    """
    Stage II: Strict filtering (SMARTS AND morphology≥2)
    """
    logger.info("Stage II: Strict flavonoid filtering...")

    flavonoids = []

    for idx, row in df.iterrows():
        smiles = row["Smiles_clean"]
        name = row.get("Name", "") or row.get("ID", "") or ""

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue

        is_flav, reasons = is_flavonoid_stage2(mol, name)

        if is_flav:
            row_dict = row.to_dict()
            row_dict["flav_match_reasons"] = "|".join(reasons)
            flavonoids.append(row_dict)

        if (idx + 1) % 2000 == 0:
            logger.info(f"Stage II: Processed {idx + 1}/{len(df)} rows, {len(flavonoids)} flavonoids")

    logger.info(f"Stage II complete: {len(flavonoids)} flavonoids from {len(df)} candidates")

    return pd.DataFrame(flavonoids)


def main():
    """Main execution"""
    # Paths
    excel_path = "data/input/CNPD-ETCM-合并去重.xlsx"
    work_dir = "data/work"
    os.makedirs(work_dir, exist_ok=True)

    raw_path = os.path.join(work_dir, "cnpd_raw.parquet")
    candidates_path = os.path.join(work_dir, "flavonoids_candidates.parquet")
    final_path = os.path.join(work_dir, "flavonoids_final.parquet")

    # Step 1: Load and standardize
    if not os.path.exists(raw_path):
        logger.info("=" * 80)
        logger.info("STEP 1: Load and standardize CNPD-ETCM database")
        logger.info("=" * 80)

        df = load_cnpd_excel(excel_path)
        df_clean = standardize_and_deduplicate(df)

        logger.info(f"Saving standardized data to {raw_path}")
        df_clean.to_parquet(raw_path, index=False)
    else:
        logger.info(f"Loading existing standardized data from {raw_path}")
        df_clean = pd.read_parquet(raw_path)

    logger.info(f"Standardized database: {len(df_clean)} unique molecules")

    # Step 2: Stage I filtering
    if not os.path.exists(candidates_path):
        logger.info("=" * 80)
        logger.info("STEP 2: Stage I - Loose flavonoid filtering")
        logger.info("=" * 80)

        df_candidates = extract_flavonoids_stage1(df_clean)

        logger.info(f"Saving candidates to {candidates_path}")
        df_candidates.to_parquet(candidates_path, index=False)
    else:
        logger.info(f"Loading existing candidates from {candidates_path}")
        df_candidates = pd.read_parquet(candidates_path)

    logger.info(f"Stage I candidates: {len(df_candidates)} molecules")

    # Step 3: Stage II filtering
    logger.info("=" * 80)
    logger.info("STEP 3: Stage II - Strict flavonoid filtering")
    logger.info("=" * 80)

    df_final = extract_flavonoids_stage2(df_candidates)

    logger.info(f"Saving final flavonoids to {final_path}")
    df_final.to_parquet(final_path, index=False)

    # Summary statistics
    logger.info("=" * 80)
    logger.info("EXTRACTION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total compounds processed: {len(df_clean)}")
    logger.info(f"Stage I candidates: {len(df_candidates)} ({100*len(df_candidates)/len(df_clean):.2f}%)")
    logger.info(f"Stage II flavonoids: {len(df_final)} ({100*len(df_final)/len(df_clean):.2f}%)")
    logger.info(f"Selectivity (Stage II / Stage I): {100*len(df_final)/max(1,len(df_candidates)):.2f}%")

    # Show match reason distribution
    if len(df_final) > 0:
        logger.info("\nMatch reason distribution:")
        all_reasons = []
        for reasons_str in df_final["flav_match_reasons"]:
            all_reasons.extend(reasons_str.split("|"))

        from collections import Counter
        reason_counts = Counter(all_reasons)
        for reason, count in reason_counts.most_common(10):
            logger.info(f"  {reason}: {count}")


if __name__ == "__main__":
    main()
