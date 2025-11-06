#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test CLI configuration pipeline to prevent R2/R6 regression issues.

These tests ensure that user configuration values from YAML are preserved
through the CLI configuration merging pipeline and actually affect enumeration.
"""

import tempfile
import yaml
import os
import sys
from pathlib import Path

import pytest

# Add the src directory to Python path for imports
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

from halogenator.cli import load_config, cmd_enum


class TestCLIConfigPipeline:
    """Test configuration pipeline from YAML through CLI to enumeration."""

    def test_r2b_naringenin_alpha_as_beta(self):
        """
        Test R2b configuration: naringenin with allow_alpha_as_beta=true should
        produce products through full CLI path.
        """
        config_content = {
            'rules': ['R2'],
            'rules_cfg': {
                'R2': {
                    'sp3_CH2_flavanone': True,
                    'allow_alpha_as_beta': True
                }
            },
            'k_max': 1,
            'halogens': ['F'],
            'io': {
                'parents_file': None,  # Will be set dynamically
                'output_dir': None     # Will be set dynamically
            }
        }

        # Naringenin SMILES (flavanone with CH2 at C3 position)
        naringenin_smiles = "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create parent SMILES file
            parents_file = tmpdir_path / "parents.smi"
            with open(parents_file, 'w') as f:
                f.write(f"{naringenin_smiles} naringenin\n")

            # Create config file
            config_file = tmpdir_path / "config.yaml"
            config_content['io']['parents_file'] = str(parents_file)
            config_content['io']['output_dir'] = str(tmpdir_path / "output")

            with open(config_file, 'w') as f:
                yaml.dump(config_content, f)

            # Load and run enumeration
            config = load_config(str(config_file))
            cmd_enum(config)

            # Check that products were generated
            output_dir = Path(config_content['io']['output_dir'])
            qa_summary_file = output_dir / "qa_summary.json"

            assert qa_summary_file.exists(), "QA summary file should be generated"

            import json
            with open(qa_summary_file) as f:
                qa_data = json.load(f)

            # Should have at least 1 product from R2b
            assert qa_data['products'] > 0, f"R2b should produce products, got {qa_data['products']}"

            # Check rule-specific statistics
            by_rule = qa_data.get('pivots', {}).get('by_rule', {})
            assert 'R2' in by_rule, "R2 rule should be present in statistics"
            assert by_rule['R2']['attempts'] > 0, "R2 should have attempts > 0"

    def test_r6_toluene_basic(self):
        """
        Test R6 configuration: toluene with R6_methyl enabled should
        produce exactly 1 product through full CLI path.
        """
        config_content = {
            'rules': ['R6'],  # This should get mapped to R6_methyl
            'rules_cfg': {
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F']
                }
            },
            'k_max': 1,
            'halogens': ['F'],
            'io': {
                'parents_file': None,  # Will be set dynamically
                'output_dir': None     # Will be set dynamically
            }
        }

        # Toluene SMILES
        toluene_smiles = "Cc1ccccc1"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create parent SMILES file
            parents_file = tmpdir_path / "parents.smi"
            with open(parents_file, 'w') as f:
                f.write(f"{toluene_smiles} toluene\n")

            # Create config file
            config_file = tmpdir_path / "config.yaml"
            config_content['io']['parents_file'] = str(parents_file)
            config_content['io']['output_dir'] = str(tmpdir_path / "output")

            with open(config_file, 'w') as f:
                yaml.dump(config_content, f)

            # Load and run enumeration
            config = load_config(str(config_file))
            cmd_enum(config)

            # Check that products were generated
            output_dir = Path(config_content['io']['output_dir'])
            qa_summary_file = output_dir / "qa_summary.json"

            assert qa_summary_file.exists(), "QA summary file should be generated"

            import json
            with open(qa_summary_file) as f:
                qa_data = json.load(f)

            # Should have exactly 1 product from R6 (toluene has 1 methyl)
            assert qa_data['products'] == 1, f"R6 on toluene should produce 1 product, got {qa_data['products']}"

            # Check that R6 attempts were made
            assert qa_data['attempts'] > 0, "R6 should have attempts > 0"

    def test_r6_prenyl_allylic_control(self):
        """
        Test R6 allylic control: prenyl benzene should produce 0 products with
        allow_allylic_methyl=false, but >=1 with allow_allylic_methyl=true.
        """
        # Prenyl benzene (benzene with prenyl side chain)
        prenyl_benzene_smiles = "CC(C)=CCc1ccccc1"

        # Test 1: allow_allylic_methyl=false (default) - should produce 0 products
        config_content_false = {
            'rules': ['R6'],
            'rules_cfg': {
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_allylic_methyl': False
                }
            },
            'k_max': 1,
            'halogens': ['F'],
            'io': {
                'parents_file': None,
                'output_dir': None
            }
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create parent SMILES file
            parents_file = tmpdir_path / "parents.smi"
            with open(parents_file, 'w') as f:
                f.write(f"{prenyl_benzene_smiles} prenyl_benzene\n")

            # Test with allow_allylic_methyl=false
            config_file = tmpdir_path / "config_false.yaml"
            config_content_false['io']['parents_file'] = str(parents_file)
            config_content_false['io']['output_dir'] = str(tmpdir_path / "output_false")

            with open(config_file, 'w') as f:
                yaml.dump(config_content_false, f)

            config = load_config(str(config_file))
            cmd_enum(config)

            # Check results for false case
            output_dir_false = Path(config_content_false['io']['output_dir'])
            qa_summary_file_false = output_dir_false / "qa_summary.json"

            assert qa_summary_file_false.exists(), "QA summary file should be generated"

            import json
            with open(qa_summary_file_false) as f:
                qa_data_false = json.load(f)

            # Should have 0 products when allylic methyls are disabled
            assert qa_data_false['products'] == 0, f"R6 on prenyl with allow_allylic_methyl=false should produce 0 products, got {qa_data_false['products']}"

        # Test 2: allow_allylic_methyl=true - should produce >=1 products
        config_content_true = {
            'rules': ['R6'],
            'rules_cfg': {
                'R6_methyl': {
                    'enable': True,
                    'allowed': ['F'],
                    'allow_allylic_methyl': True
                }
            },
            'k_max': 1,
            'halogens': ['F'],
            'io': {
                'parents_file': None,
                'output_dir': None
            }
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create parent SMILES file
            parents_file = tmpdir_path / "parents.smi"
            with open(parents_file, 'w') as f:
                f.write(f"{prenyl_benzene_smiles} prenyl_benzene\n")

            # Test with allow_allylic_methyl=true
            config_file = tmpdir_path / "config_true.yaml"
            config_content_true['io']['parents_file'] = str(parents_file)
            config_content_true['io']['output_dir'] = str(tmpdir_path / "output_true")

            with open(config_file, 'w') as f:
                yaml.dump(config_content_true, f)

            config = load_config(str(config_file))
            cmd_enum(config)

            # Check results for true case
            output_dir_true = Path(config_content_true['io']['output_dir'])
            qa_summary_file_true = output_dir_true / "qa_summary.json"

            assert qa_summary_file_true.exists(), "QA summary file should be generated"

            with open(qa_summary_file_true) as f:
                qa_data_true = json.load(f)

            # Should have >=1 products when allylic methyls are enabled
            assert qa_data_true['products'] > 0, f"R6 on prenyl with allow_allylic_methyl=true should produce >=1 products, got {qa_data_true['products']}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])