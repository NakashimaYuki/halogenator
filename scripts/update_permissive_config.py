#!/usr/bin/env python3
"""
Update halogen_rules_by_class.yaml with permissive enumeration strategy.
Strategy: Maximum enumeration with safety nets only for extreme cases.
"""

import yaml
from pathlib import Path

# Define the new permissive configuration
config = {
    'defaults': {
        'halogens': ['F', 'Cl', 'Br', 'I'],
        'batch_size': 5000,
        'rdkit_threads': 8
    },
    'classes': {
        # Terpenoids - No global max_sites; rule-level safety nets
        'terpenoid': {
            'description': 'Isoprene-derived polycyclic natural products',
            'molecule_count': 28606,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX', 'PRIMARY_OH__CH2OH__TO__X'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX', 'ALPHA_CARBONYL__CH2__TO__X', 'PRIMARY_OH__CH2OH__TO__X'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'RING_SP3__CH__TO__X': {'max_sites_per_parent': 6},
                'ALPHA_CARBONYL__CH2__TO__X': {'max_sites_per_parent': 3},
                'COOH__TO__CX': {'max_sites_per_parent': 2},
                'PRIMARY_OH__CH2OH__TO__X': {'max_sites_per_parent': 3}
            }
        },
        # Alkaloids - Very permissive aromatic halogenation
        'alkaloid': {
            'description': 'N-containing heterocyclic natural products',
            'molecule_count': 4287,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX', 'RING_SP3__CH__TO__X'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX', 'RING_SP3__CH__TO__X', 'ALPHA_CARBONYL__CH2__TO__X'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'RING_SP2__CH__TO__X': {'max_sites_per_parent': 10},
                'RING_SP3__CH__TO__X': {'max_sites_per_parent': 4},
                'ALPHA_CARBONYL__CH2__TO__X': {'max_sites_per_parent': 4},
                'COOH__TO__CX': {'max_sites_per_parent': 4},
                'PRIMARY_OH__CH2OH__TO__X': {'enabled': False}
            }
        },
        # Polyphenols - Very high aromatic limit (12)
        'polyphenol': {
            'description': 'Flavonoids, lignans, phenylpropanoids, coumarins',
            'molecule_count': 4047,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'RING_SP2__CH__TO__X': {'max_sites_per_parent': 12},
                'COOH__TO__CX': {'max_sites_per_parent': 3},
                'PRIMARY_OH__CH2OH__TO__X': {'enabled': False}
            }
        },
        # Glycosides - Focus on aglycone with sugar exclusion
        'glycoside': {
            'description': 'Glycosides with aglycone + sugar moieties',
            'molecule_count': 22873,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': True
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'RING_SP3__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': True
            },
            'per_rule_overrides': {
                'RING_SP2__CH__TO__X': {'exclude_sugar_like': True},
                'RING_SP3__CH__TO__X': {'exclude_sugar_like': True},
                'COOH__TO__CX': {'exclude_sugar_like': True},
                'PRIMARY_OH__CH2OH__TO__X': {'enabled': False}
            }
        },
        # Polysaccharides - DISABLED
        'polysaccharide': {
            'description': 'Oligosaccharides and polysaccharides',
            'molecule_count': 265,
            'enabled': False,
            'k1': {'include_rules': [], 'sugar_mask': False},
            'k2': {'include_rules': [], 'sugar_mask': False}
        },
        # Lipids - Rule-level safety nets at 3
        'lipid': {
            'description': 'Fatty acids, glycerides, phospholipids',
            'molecule_count': 18,
            'enabled': True,
            'k1': {
                'include_rules': ['COOH__TO__CX', 'ALPHA_CARBONYL__CH2__TO__X', 'PRIMARY_OH__CH2OH__TO__X'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['COOH__TO__CX', 'ALPHA_CARBONYL__CH2__TO__X', 'PRIMARY_OH__CH2OH__TO__X'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'COOH__TO__CX': {'max_sites_per_parent': 3},
                'ALPHA_CARBONYL__CH2__TO__X': {'max_sites_per_parent': 3},
                'PRIMARY_OH__CH2OH__TO__X': {'max_sites_per_parent': 3},
                'RING_SP2__CH__TO__X': {'enabled': False},
                'RING_SP3__CH__TO__X': {'enabled': False}
            }
        },
        # AA/Peptides - Soft constraints at 2
        'aa_peptide': {
            'description': 'Amino acids, peptides, small proteins',
            'molecule_count': 3890,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX', 'ALPHA_CARBONYL__CH2__TO__X'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'RING_SP2__CH__TO__X': {'max_sites_per_parent': 2},
                'COOH__TO__CX': {'max_sites_per_parent': 2},
                'ALPHA_CARBONYL__CH2__TO__X': {'max_sites_per_parent': 2},
                'RING_SP3__CH__TO__X': {'enabled': False},
                'PRIMARY_OH__CH2OH__TO__X': {'enabled': False}
            }
        },
        # Other - Moderate safety nets
        'other': {
            'description': 'Unclassified natural products',
            'molecule_count': 4262,
            'enabled': True,
            'k1': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX'],
                'sugar_mask': False
            },
            'k2': {
                'include_rules': ['RING_SP2__CH__TO__X', 'COOH__TO__CX', 'RING_SP3__CH__TO__X'],
                'sugar_mask': False
            },
            'per_rule_overrides': {
                'RING_SP2__CH__TO__X': {'max_sites_per_parent': 8},
                'COOH__TO__CX': {'max_sites_per_parent': 3}
            }
        }
    }
}

def main():
    output_path = Path('E:/Projects/halogenator/configs/halogen_rules_by_class.yaml')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('# Per-class halogenation rule configurations\n')
        f.write('# STRATEGY: Maximum enumeration with safety nets only for extreme cases\n')
        f.write('# Updated: 2025-12-01 - Permissive enumeration strategy\n')
        f.write('#\n')
        f.write('# Key changes from conservative strategy:\n')
        f.write('# - Removed global max_sites_per_parent limits\n')
        f.write('# - Rule-level safety nets only (3-12 depending on rule)\n')
        f.write('# - Polysaccharide: disabled entirely\n')
        f.write('# - Glycoside: exclude_sugar_like on all rules\n')
        f.write('# - Terpenoid: enabled PRIMARY_OH exploration\n')
        f.write('# - Alkaloid: aromatic sites up to 10\n')
        f.write('# - Polyphenol: aromatic sites up to 12\n')
        f.write('\n')
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    print(f'[OK] Successfully updated {output_path}')
    print('\nKey changes:')
    print('  - Terpenoid: No global max_sites, rule-level safety nets (2-6)')
    print('  - Alkaloid: Aromatic up to 10 sites')
    print('  - Polyphenol: Aromatic up to 12 sites')
    print('  - Glycoside: exclude_sugar_like on all rules')
    print('  - Polysaccharide: DISABLED')
    print('  - Lipid: All rules capped at 3')
    print('  - AA_peptide: Soft limits at 2')

if __name__ == '__main__':
    main()
