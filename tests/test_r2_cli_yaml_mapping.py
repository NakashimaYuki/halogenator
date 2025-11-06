# -*- coding: ascii -*-
"""
Test R2 CLI and YAML configuration mapping.

This test ensures that R2a/R2b settings can be properly controlled
through both CLI arguments and YAML configuration files.
"""

import tempfile
import unittest
import sys
import yaml
from pathlib import Path
from unittest.mock import patch

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import EnumConfig


class TestR2CliYamlMapping(unittest.TestCase):
    """Test R2 CLI and YAML configuration mapping."""

    def test_default_r2_config(self):
        """Test that R2 rules are disabled by default."""
        config = EnumConfig()

        # Check default R2 configuration
        self.assertIn('R2', config.rules_cfg)
        self.assertEqual(config.rules_cfg['R2']['sp2_CH_in_C_ring'], False)
        self.assertEqual(config.rules_cfg['R2']['sp3_CH2_flavanone'], False)

    def test_yaml_r2_config_mapping(self):
        """Test that YAML configuration properly maps to EnumConfig."""
        # Create a test YAML configuration
        yaml_config = {
            'rules_cfg': {
                'R2': {
                    'sp2_CH_in_C_ring': True,
                    'sp3_CH2_flavanone': True
                }
            }
        }

        # Create EnumConfig with YAML-style rules_cfg
        config = EnumConfig(rules_cfg=yaml_config['rules_cfg'])

        # Check that YAML config is properly applied
        self.assertEqual(config.rules_cfg['R2']['sp2_CH_in_C_ring'], True)
        self.assertEqual(config.rules_cfg['R2']['sp3_CH2_flavanone'], True)

    def test_partial_yaml_r2_config(self):
        """Test that partial YAML configuration works correctly."""
        # Test enabling only R2a
        yaml_config_r2a = {
            'rules_cfg': {
                'R2': {
                    'sp2_CH_in_C_ring': True,
                    'sp3_CH2_flavanone': False
                }
            }
        }

        config_r2a = EnumConfig(rules_cfg=yaml_config_r2a['rules_cfg'])
        self.assertEqual(config_r2a.rules_cfg['R2']['sp2_CH_in_C_ring'], True)
        self.assertEqual(config_r2a.rules_cfg['R2']['sp3_CH2_flavanone'], False)

        # Test enabling only R2b
        yaml_config_r2b = {
            'rules_cfg': {
                'R2': {
                    'sp2_CH_in_C_ring': False,
                    'sp3_CH2_flavanone': True
                }
            }
        }

        config_r2b = EnumConfig(rules_cfg=yaml_config_r2b['rules_cfg'])
        self.assertEqual(config_r2b.rules_cfg['R2']['sp2_CH_in_C_ring'], False)
        self.assertEqual(config_r2b.rules_cfg['R2']['sp3_CH2_flavanone'], True)

    def test_yaml_file_loading_with_r2_config(self):
        """Test that R2 configuration can be loaded from actual YAML files."""
        # Create a temporary YAML file with R2 configuration
        yaml_content = """
        k_max: 2
        halogens: [F, Cl]
        rules: [R1, R2]
        rules_cfg:
          R2:
            sp2_CH_in_C_ring: true
            sp3_CH2_flavanone: false
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(yaml_content)
            yaml_path = f.name

        try:
            # Load the YAML configuration
            with open(yaml_path, 'r') as f:
                config_dict = yaml.safe_load(f)

            # Create EnumConfig from loaded configuration
            enum_config = EnumConfig(
                k_max=config_dict['k_max'],
                halogens=tuple(config_dict['halogens']),
                rules=tuple(config_dict['rules']),
                rules_cfg=config_dict['rules_cfg']
            )

            # Verify the configuration
            self.assertEqual(enum_config.k_max, 2)
            self.assertEqual(enum_config.halogens, ('F', 'Cl'))
            self.assertIn('R2', enum_config.rules)
            self.assertEqual(enum_config.rules_cfg['R2']['sp2_CH_in_C_ring'], True)
            self.assertEqual(enum_config.rules_cfg['R2']['sp3_CH2_flavanone'], False)

        finally:
            Path(yaml_path).unlink(missing_ok=True)

    def test_missing_r2_config_defaults(self):
        """Test that missing R2 configuration gets reasonable defaults."""
        # Test with no rules_cfg at all
        config_no_rules_cfg = EnumConfig()
        self.assertIn('R2', config_no_rules_cfg.rules_cfg)
        self.assertEqual(config_no_rules_cfg.rules_cfg['R2']['sp2_CH_in_C_ring'], False)
        self.assertEqual(config_no_rules_cfg.rules_cfg['R2']['sp3_CH2_flavanone'], False)

        # Test with rules_cfg present but no R2 section
        empty_rules_cfg = {}
        config_empty_rules = EnumConfig(rules_cfg=empty_rules_cfg)
        # Should get defaults from EnumConfig.__init__
        self.assertIn('R2', config_empty_rules.rules_cfg)

    def test_cli_args_simulation(self):
        """
        Test that CLI-style arguments can be properly processed into rules_cfg.

        This simulates how the CLI parsing would work.
        """
        # Simulate CLI args processing for --enable-R2a
        base_rules_cfg = {'R2': {'sp2_CH_in_C_ring': False, 'sp3_CH2_flavanone': False}}

        # Simulate CLI processing logic
        args_r2a_enabled = type('Args', (), {'enable_r2a': True, 'enable_r2b': None})()

        rules_cfg_after_cli = dict(base_rules_cfg)
        if getattr(args_r2a_enabled, 'enable_r2a', None) is not None:
            if 'R2' not in rules_cfg_after_cli:
                rules_cfg_after_cli['R2'] = {}
            rules_cfg_after_cli['R2']['sp2_CH_in_C_ring'] = bool(args_r2a_enabled.enable_r2a)

        if getattr(args_r2a_enabled, 'enable_r2b', None) is not None:
            if 'R2' not in rules_cfg_after_cli:
                rules_cfg_after_cli['R2'] = {}
            rules_cfg_after_cli['R2']['sp3_CH2_flavanone'] = bool(args_r2a_enabled.enable_r2b)

        # Create config with CLI-processed rules_cfg
        config = EnumConfig(rules_cfg=rules_cfg_after_cli)

        # Should have R2a enabled, R2b disabled
        self.assertEqual(config.rules_cfg['R2']['sp2_CH_in_C_ring'], True)
        self.assertEqual(config.rules_cfg['R2']['sp3_CH2_flavanone'], False)

        # Test the opposite case: --enable-R2b only
        args_r2b_enabled = type('Args', (), {'enable_r2a': None, 'enable_r2b': True})()

        # Start with fresh base config for R2b test
        rules_cfg_r2b = {'R2': {'sp2_CH_in_C_ring': False, 'sp3_CH2_flavanone': False}}
        if getattr(args_r2b_enabled, 'enable_r2a', None) is not None:
            rules_cfg_r2b['R2']['sp2_CH_in_C_ring'] = bool(args_r2b_enabled.enable_r2a)

        if getattr(args_r2b_enabled, 'enable_r2b', None) is not None:
            rules_cfg_r2b['R2']['sp3_CH2_flavanone'] = bool(args_r2b_enabled.enable_r2b)

        config_r2b = EnumConfig(rules_cfg=rules_cfg_r2b)

        # Should have R2a disabled, R2b enabled
        self.assertEqual(config_r2b.rules_cfg['R2']['sp2_CH_in_C_ring'], False)
        self.assertEqual(config_r2b.rules_cfg['R2']['sp3_CH2_flavanone'], True)


if __name__ == '__main__':
    unittest.main(verbosity=2)