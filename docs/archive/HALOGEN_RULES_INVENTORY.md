# Halogenation Rules Inventory

**Generated:** 2025-12-04
**Purpose:** Complete reference of all halogenation rules available in the halogenator codebase

---

## Quick Reference Table

| Rule ID | Legacy | Semantic Name | Current Usage | Implementation | SMIRKS |
|---------|--------|---------------|---------------|----------------|--------|
| R1 | R1 | RING_SP2__CH__TO__X | ✅ ACTIVE | ✅ Complete | `[cH:1]>>[c:1]X` |
| R3 | R3 | - | ❌ NOT IN YAML | ✅ Complete | `[#6;!$(C=O):1]-[O;H1]>>[#6:1]X` |
| R4 | R4 | - | ❌ NOT IN YAML | ✅ Complete | `[#6;!$(C=O):1]-[N;H2,H1]>>[#6:1]X` |
| R5 | R5 | COOH__TO__CX | ✅ ACTIVE | ✅ Complete | `[c,#6:1]-C(=O)[O;H1]>>[#6:1]X` |
| R6 | R6_methyl | - | ❌ NOT IN YAML | ⚠️ Partial (config-gated) | Not implemented |
| - | - | RING_SP3__CH__TO__X | ⚠️ IN YAML, 0 PRODUCTS | ✅ Complete | `[C;R;H1,H2:1]>>[C;R:1]X` |
| - | - | ALPHA_CARBONYL__CH2__TO__X | ⚠️ IN YAML, 0 PRODUCTS | ✅ Complete | `[CH2:1][CX3](=O)>>[CH:1]([X])[CX3](=O)` |
| - | - | PRIMARY_OH__CH2OH__TO__X | ⚠️ IN YAML, 0 PRODUCTS | ✅ Complete | `[CH2:1][OX2H]>>[CH2:1]X` |
| R2 | R2 | - | ⚠️ FLAVONE-SPECIFIC | ✅ Complete | Multiple variants |

---

## Detailed Rule Specifications

### R1: Aromatic/Heteroaromatic C-H Halogenation

**Legacy Name:** R1
**Semantic Name:** RING_SP2__CH__TO__X
**Status:** ✅ ACTIVE - Produces 90-98% of current k=1 products

**SMIRKS:**
```
[cH:1]>>[c:1]X
```

**Description:**
Halogenates aromatic or heteroaromatic C-H positions. Matches any aromatic carbon with a hydrogen atom.

**Chemical Scope:**
- Benzene rings
- Pyridine, pyrrole, furan, thiophene rings
- Polycyclic aromatics (naphthalene, anthracene, etc.)
- Indole, quinoline, isoquinoline scaffolds

**Typical Products/Parent:**
- Polyphenol: 15-20 (multiple aromatic rings)
- Alkaloid: 15-20 (heteroaromatic scaffolds)
- Terpenoid: 5-10 (fewer aromatic positions)

**Mapping:** RING_SP2__CH__TO__X → R1

---

### R3: Hydroxyl Halogenation

**Legacy Name:** R3
**Semantic Name:** None (use legacy R3)
**Status:** ❌ NOT IN YAML CONFIG (CRITICAL MISSING RULE)

**SMIRKS:**
```
[#6;!$(C=O):1]-[O;H1]>>[#6:1]X
```

**Description:**
Replaces hydroxyl (-OH) with halogen, excluding carboxylic acid OH (handled by R5).

**Chemical Scope:**
- Phenolic -OH (tyrosine, flavonoids, tannins)
- Aliphatic alcohols (terpenols, glycerol)
- Secondary/tertiary alcohols
- Excludes: Carboxylic acid -OH

**Expected Impact if Enabled:**
- **Polyphenol:** 2-5x increase (phenolic OH abundant)
- **Terpenoid:** 2-3x increase (terpene alcohols common)
- **Glycoside:** Moderate (aglycone OH only, sugar masked)
- **Lipid:** Small (glycerol backbone)

**Recommendation:** **HIGH PRIORITY - ADD TO YAML CONFIG**

**Mapping:** R3 → R3

---

### R4: Amine Halogenation

**Legacy Name:** R4
**Semantic Name:** None (use legacy R4)
**Status:** ❌ NOT IN YAML CONFIG

**SMIRKS:**
```
[#6;!$(C=O):1]-[N;H2,H1]>>[#6:1]X
```

**Description:**
Replaces primary or secondary amine with halogen, excluding amide nitrogen.

**Chemical Scope:**
- Primary amines (-NH₂): amino acids, aminoglycosides
- Secondary amines (-NH-): alkaloids, peptides
- Excludes: Amide nitrogen, tertiary amines

**Expected Impact if Enabled:**
- **Alkaloid:** Moderate-High (N-containing scaffolds)
- **AA_peptide:** Moderate (side chains only, NOT backbone)
- **Other classes:** Low (rare)

**Recommendation:** MEDIUM PRIORITY - ADD selectively (alkaloid, aa_peptide)

**Warning:** Be conservative with aa_peptide - avoid backbone amine halogenation

**Mapping:** R4 → R4

---

### R5: Carboxyl Halogenation

**Legacy Name:** R5
**Semantic Name:** COOH__TO__CX
**Status:** ✅ ACTIVE - Produces 2-10% of current k=1 products

**SMIRKS:**
```
[c,#6:1]-C(=O)[O;H1]>>[#6:1]X
```

**Description:**
Replaces entire carboxyl group (-C(=O)OH) with halogen attached to the α-carbon.

**Chemical Scope:**
- Aromatic acids (benzoic, phenylacetic)
- Aliphatic acids (fatty acids, amino acid C-terminus)
- Excludes: Esters, amides (no -OH)

**Typical Products/Parent:**
- Lipid: High (fatty acids)
- AA_peptide: Low-Moderate (C-terminus)
- Other classes: Low (carboxylic acids less common)

**Mapping:** COOH__TO__CX → R5

---

### RING_SP3__CH__TO__X: Ring sp³ C-H Halogenation

**Legacy Name:** None (semantic only)
**Semantic Name:** RING_SP3__CH__TO__X
**Status:** ⚠️ IN YAML CONFIG - Produced 0 products in k=1

**SMIRKS:**
```
[C;R;H1,H2:1]>>[C;R:1]X
```

**Description:**
Halogenates sp³ hybridized carbons in ring systems (saturated or partially saturated rings).

**Chemical Scope:**
- Cyclohexane, cyclopentane rings
- Partially saturated heterocycles (piperidine, pyrrolidine)
- Terpene ring systems (decalin, steroid-like)

**Why 0 Products in k=1:**
- Natural products in dataset are predominantly aromatic
- sp³ rings less common than sp² aromatics
- May be too strict max_sites_per_parent (6 for terpenoid)

**Recommendation:** Keep in config, may produce products in:
- k=2 enumeration
- Specific terpenoid/steroid-rich subsets

**Mapping:** RING_SP3__CH__TO__X → RING_SP3__CH__TO__X (pass-through)

---

### ALPHA_CARBONYL__CH2__TO__X: α-Carbonyl Halogenation

**Legacy Name:** None (semantic only)
**Semantic Name:** ALPHA_CARBONYL__CH2__TO__X
**Status:** ⚠️ IN YAML CONFIG (mostly k=2) - Produced 0 products in k=1

**SMIRKS:**
```
[CH2:1][CX3](=O)>>[CH:1]([X])[CX3](=O)
```

**Description:**
Halogenates α-position (CH₂) adjacent to carbonyl groups. Similar to Hell-Volhard-Zelinsky reaction.

**Chemical Scope:**
- Ketones with α-CH₂
- Aldehydes with α-CH₂
- Esters/amides with α-CH₂
- Excludes: α-CH (quaternary) or no α-hydrogen

**Why 0 Products in k=1:**
- Only enabled in k=2 for most classes
- Motif may be rare in dataset

**Recommendation:** Keep in k=2 configs, potentially add to k=1 for terpenoid/alkaloid

**Mapping:** ALPHA_CARBONYL__CH2__TO__X → ALPHA_CARBONYL__CH2__TO__X (pass-through)

---

### PRIMARY_OH__CH2OH__TO__X: Primary Alcohol Halogenation

**Legacy Name:** None (semantic only)
**Semantic Name:** PRIMARY_OH__CH2OH__TO__X
**Status:** ⚠️ IN YAML CONFIG (limited) - Produced 0 products in k=1

**SMIRKS:**
```
[CH2:1][OX2H]>>[CH2:1]X
```

**Description:**
Converts primary alcohols (-CH₂OH) to alkyl halides (-CH₂X).

**Chemical Scope:**
- Primary alcohols: ethanol-like, fatty alcohols
- Glycerol -CH₂OH
- Excludes: Phenolic OH (handled by R3), secondary/tertiary alcohols

**Why 0 Products in k=1:**
- Disabled in many classes (alkaloid, aa_peptide, polyphenol)
- Enabled only in lipid/terpenoid but rare motif
- May overlap with R3 (which is more general)

**Recommendation:**
- Consider removing in favor of R3 (more general -OH halogenation)
- Or keep for specific primary alcohol targeting

**Mapping:** PRIMARY_OH__CH2OH__TO__X → PRIMARY_OH__CH2OH__TO__X (pass-through)

---

### R2: Flavone C-Ring Halogenation

**Legacy Name:** R2 (with extensions R2a, R2b)
**Semantic Name:** None
**Status:** ⚠️ FLAVONE-SPECIFIC (not used in per-class configs)

**Description:**
Specialized rules for flavonoid C-ring halogenation:
- **R2a:** sp² C-H in C-ring (α to carbonyl)
- **R2b:** sp³ C-H₂ in flavanone C-ring

**Chemical Scope:**
- Flavones, flavonols, flavanones
- C-ring positions adjacent to or in the heterocyclic ring

**Status:** Not included in current per-class configs (flavonoids likely in polyphenol class)

**Recommendation:** Only enable if processing flavonoid-specific library

---

### R6_methyl: Methyl Halogenation

**Legacy Name:** R6_methyl
**Semantic Name:** None
**Status:** ⚠️ MENTIONED BUT NOT IMPLEMENTED

**Description:**
Intended for benzylic/allylic methyl halogenation (not yet implemented).

**Implementation Status:**
- Mentioned in `rules.build_reactions()` dict
- No SMIRKS function defined
- Gated by config: `cfg.rules_cfg.get('R6_methyl', {}).get('enable', False)`

**Recommendation:** Low priority - not needed for initial expansion

---

## Rule Processing Architecture

### In enumerate_k.py

**Reaction-type rules** (via `_apply_reaction_rule`):
- R3, R4, R5
- RING_SP3__CH__TO__X
- ALPHA_CARBONYL__CH2__TO__X
- PRIMARY_OH__CH2OH__TO__X

**Site-type rules** (via `_apply_site_rules`):
- R1 / RING_SP2__CH__TO__X
- R2 (flavone-specific)
- R6_methyl (config-gated)

---

## Name Mapping Reference

### Semantic → Legacy

| Semantic Name | Legacy Name | Use Case |
|---------------|-------------|----------|
| RING_SP2__CH__TO__X | R1 | Preferred in configs, auto-mapped |
| COOH__TO__CX | R5 | Preferred in configs, auto-mapped |
| RING_SP3__CH__TO__X | RING_SP3__CH__TO__X | Use as-is (no legacy alias) |
| ALPHA_CARBONYL__CH2__TO__X | ALPHA_CARBONYL__CH2__TO__X | Use as-is (no legacy alias) |
| PRIMARY_OH__CH2OH__TO__X | PRIMARY_OH__CH2OH__TO__X | Use as-is (no legacy alias) |

### Legacy → Legacy

| Legacy Name | Use Case |
|-------------|----------|
| R1 | Accepted in configs, maps to self |
| R3 | **ADD TO YAML** - hydroxyl halogenation |
| R4 | **ADD TO YAML** - amine halogenation |
| R5 | Accepted in configs, maps to self |

---

## Priority Action Items

### Immediate (Phase 2)

1. **Add R3 to YAML configs:**
   - polyphenol: max_sites 6-8
   - terpenoid: max_sites 4-6
   - glycoside: max_sites 4 (with sugar_mask)
   - alkaloid: max_sites 3-4
   - lipid: max_sites 2-3
   - other: max_sites 4

2. **Add R4 to YAML configs (selective):**
   - alkaloid: max_sites 3-4
   - aa_peptide: max_sites 1-2 (side chains only)
   - other: max_sites 2-3

### Future Investigation

1. **RING_SP3__CH__TO__X:** Why 0 products? Debug or accept as rare motif
2. **PRIMARY_OH vs R3:** Consider removing PRIMARY_OH in favor of more general R3
3. **R6_methyl:** Implement if benzylic/allylic halogenation is desired

---

**Document Status:** ✅ Complete
**Next Steps:** Phase 2.1 - Rewrite halogen_rules_by_class.yaml with R3/R4
