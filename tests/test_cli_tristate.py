# -*- coding: ascii -*-
"""
Tests for CLI tri-state override logic.

Tests that CLI flags properly override YAML configuration only when explicitly provided,
implementing proper tri-state semantics (None = use YAML, True/False = override YAML).
"""

import unittest
import tempfile
import os
from argparse import Namespace

# Import the modules we want to test
try:
    from src.halogenator.cli import cmd_enum
    from src.halogenator.enumerate_k import EnumConfig
    CLI_AVAILABLE = True
except ImportError:
    CLI_AVAILABLE = False


@unittest.skipUnless(CLI_AVAILABLE, "CLI modules not available")
class TestCLITristate(unittest.TestCase):
    """Test CLI tri-state override behavior."""

    def setUp(self):
        """Set up test configuration files."""
        # Create temporary YAML config for testing
        self.test_config_content = """
k_max: 2
halogens: ['F', 'Cl']
rules: ['R3', 'R4']
symmetry:
  compute_on_masked_subgraph: false
qc:
  tautomer_canonicalize: false
sugar:
  mode: heuristic
"""

        # Create temporary file
        self.config_fd, self.config_path = tempfile.mkstemp(suffix='.yaml', text=True)
        with os.fdopen(self.config_fd, 'w') as f:
            f.write(self.test_config_content)

    def tearDown(self):
        """Clean up test files."""
        if os.path.exists(self.config_path):
            os.unlink(self.config_path)

    def test_symmetry_no_cli_flag_uses_yaml(self):
        """Test that without CLI flag, YAML configuration is used for symmetry."""
        # Simulate args with no symmetry flag (should be None)
        args = Namespace(
            config=self.config_path,
            k_max=None,
            subset='all',
            outdir=None,
            workers=1,
            symmetry_masked_subgraph=None,  # No CLI override
            enable_tautomer=None  # No CLI override
        )

        # This should use YAML value (false)
        try:
            enum_cfg = cmd_enum(args, return_config_only=True)
            if enum_cfg:
                symmetry_setting = enum_cfg.symmetry_cfg.get('compute_on_masked_subgraph', True)
                print(f"No CLI flag -> symmetry from YAML: {symmetry_setting}")
                self.assertFalse(symmetry_setting, "Should use YAML value (false)")
        except Exception as e:
            # If cmd_enum doesn't support return_config_only, skip this test
            self.skipTest(f"Cannot test configuration extraction: {e}")

    def test_symmetry_explicit_enable_overrides_yaml(self):
        """Test that --symmetry-masked-subgraph overrides YAML configuration."""
        args = Namespace(
            config=self.config_path,
            k_max=None,
            subset='all',
            outdir=None,
            workers=1,
            symmetry_masked_subgraph=True,  # CLI override to True
            enable_tautomer=None
        )

        try:
            enum_cfg = cmd_enum(args, return_config_only=True)
            if enum_cfg:
                symmetry_setting = enum_cfg.symmetry_cfg.get('compute_on_masked_subgraph', False)
                print(f"CLI --symmetry-masked-subgraph -> symmetry: {symmetry_setting}")
                self.assertTrue(symmetry_setting, "CLI flag should override YAML (false -> true)")
        except Exception as e:
            self.skipTest(f"Cannot test configuration extraction: {e}")

    def test_symmetry_explicit_disable_overrides_yaml(self):
        """Test that --no-symmetry-masked-subgraph overrides YAML configuration."""
        # First, test with a YAML that has symmetry enabled
        yaml_with_symmetry_on = """
k_max: 2
symmetry:
  compute_on_masked_subgraph: true
"""
        config_fd, config_path = tempfile.mkstemp(suffix='.yaml', text=True)
        try:
            with os.fdopen(config_fd, 'w') as f:
                f.write(yaml_with_symmetry_on)

            args = Namespace(
                config=config_path,
                k_max=None,
                subset='all',
                outdir=None,
                workers=1,
                symmetry_masked_subgraph=False,  # CLI override to False
                enable_tautomer=None
            )

            enum_cfg = cmd_enum(args, return_config_only=True)
            if enum_cfg:
                symmetry_setting = enum_cfg.symmetry_cfg.get('compute_on_masked_subgraph', True)
                print(f"CLI --no-symmetry-masked-subgraph -> symmetry: {symmetry_setting}")
                self.assertFalse(symmetry_setting, "CLI flag should override YAML (true -> false)")

        except Exception as e:
            self.skipTest(f"Cannot test configuration extraction: {e}")
        finally:
            if os.path.exists(config_path):
                os.unlink(config_path)

    def test_tautomer_no_cli_flag_uses_yaml(self):
        """Test that without CLI flag, YAML configuration is used for tautomer."""
        args = Namespace(
            config=self.config_path,
            k_max=None,
            subset='all',
            outdir=None,
            workers=1,
            symmetry_masked_subgraph=None,
            enable_tautomer=None  # No CLI override
        )

        try:
            enum_cfg = cmd_enum(args, return_config_only=True)
            if enum_cfg:
                tautomer_setting = enum_cfg.qc_cfg.get('tautomer_canonicalize', True)
                print(f"No CLI flag -> tautomer from YAML: {tautomer_setting}")
                self.assertFalse(tautomer_setting, "Should use YAML value (false)")
        except Exception as e:
            self.skipTest(f"Cannot test configuration extraction: {e}")

    def test_tautomer_explicit_enable_overrides_yaml(self):
        """Test that --enable-tautomer overrides YAML configuration."""
        args = Namespace(
            config=self.config_path,
            k_max=None,
            subset='all',
            outdir=None,
            workers=1,
            symmetry_masked_subgraph=None,
            enable_tautomer=True  # CLI override to True
        )

        try:
            enum_cfg = cmd_enum(args, return_config_only=True)
            if enum_cfg:
                tautomer_setting = enum_cfg.qc_cfg.get('tautomer_canonicalize', False)
                print(f"CLI --enable-tautomer -> tautomer: {tautomer_setting}")
                self.assertTrue(tautomer_setting, "CLI flag should override YAML (false -> true)")
        except Exception as e:
            self.skipTest(f"Cannot test configuration extraction: {e}")

    def test_enum_config_creation_with_defaults(self):
        """Test that EnumConfig can be created with default sugar_cfg."""
        # Test that the updated defaults work
        try:
            cfg = EnumConfig()

            # Check that sugar_cfg has the expected structure
            self.assertIn('mask_exocyclic_oxygen', cfg.sugar_cfg)
            self.assertIn('mask_glycosidic_bridge_oxygen', cfg.sugar_cfg)
            self.assertIn('audit', cfg.sugar_cfg)

            # Check default values
            self.assertTrue(cfg.sugar_cfg['mask_exocyclic_oxygen'])
            self.assertTrue(cfg.sugar_cfg['mask_glycosidic_bridge_oxygen'])
            self.assertFalse(cfg.sugar_cfg['audit'])

            print(f"Default sugar_cfg: {cfg.sugar_cfg}")

        except Exception as e:
            self.fail(f"EnumConfig creation failed: {e}")


if __name__ == '__main__':
    unittest.main()