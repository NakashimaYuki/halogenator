"""
Test pre-commit config parsability and format compliance.

This module tests that .pre-commit-config.yaml is parsable and follows
single-line ASCII requirements for all entry commands.
"""
import unittest
import os
import yaml


class TestPrecommitConfigParsable(unittest.TestCase):
    """Test pre-commit configuration parsability and format compliance."""

    def setUp(self):
        """Set up test fixtures."""
        self.config_path = '.pre-commit-config.yaml'
        self.assertFileExists(self.config_path, "Pre-commit config file must exist")

    def assertFileExists(self, file_path, message=None):
        """Assert that a file exists."""
        if not os.path.exists(file_path):
            self.fail(message or f"File {file_path} does not exist")

    def test_precommit_config_yaml_is_parsable(self):
        """Test that .pre-commit-config.yaml is valid YAML."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        # Should be a dictionary with 'repos' key
        self.assertIsInstance(config, dict)
        self.assertIn('repos', config)
        self.assertIsInstance(config['repos'], list)

    def test_all_entries_are_single_line_ascii_strings(self):
        """Test that all entry commands are single-line ASCII strings."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        # Traverse all hooks and check entries
        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    if 'entry' in hook:
                        entry = hook['entry']
                        hook_id = hook.get('id', 'unknown')

                        # Entry must be a string
                        self.assertIsInstance(entry, str,
                                            f"Entry for hook '{hook_id}' must be a string")

                        # Entry must be single line (no newlines)
                        self.assertNotIn('\n', entry,
                                       f"Entry for hook '{hook_id}' must be single line")
                        self.assertNotIn('\r', entry,
                                       f"Entry for hook '{hook_id}' must be single line")

                        # Entry must be ASCII-only
                        try:
                            entry.encode('ascii')
                        except UnicodeEncodeError:
                            self.fail(f"Entry for hook '{hook_id}' contains non-ASCII characters")

    def test_specific_hooks_are_present(self):
        """Test that expected hooks are present."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        # Collect all hook IDs
        hook_ids = []
        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    if 'id' in hook:
                        hook_ids.append(hook['id'])

        # Expected hooks should be present
        expected_hooks = ['ascii-check']
        for expected_hook in expected_hooks:
            self.assertIn(expected_hook, hook_ids,
                         f"Expected hook '{expected_hook}' not found in config")

    def test_ascii_check_hook_configuration(self):
        """Test that ascii-check hook is properly configured."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        # Find ascii-check hook
        ascii_check_hook = None
        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    if hook.get('id') == 'ascii-check':
                        ascii_check_hook = hook
                        break

        self.assertIsNotNone(ascii_check_hook, "ascii-check hook must be present")

        # Check hook configuration
        self.assertEqual(ascii_check_hook['entry'], 'python scripts/check_ascii_portable.py')
        self.assertEqual(ascii_check_hook['language'], 'system')
        self.assertEqual(ascii_check_hook['pass_filenames'], False)

    def test_python_entry_commands_are_valid(self):
        """Test that Python entry commands are syntactically valid."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    entry = hook.get('entry', '')
                    hook_id = hook.get('id', 'unknown')

                    # If entry starts with "python -c", check it's a valid Python command
                    if entry.startswith('python -c '):
                        # Extract the Python code part (after 'python -c ')
                        python_code = entry[len('python -c '):]

                        # Remove quotes if present
                        if python_code.startswith('"') and python_code.endswith('"'):
                            python_code = python_code[1:-1]
                        elif python_code.startswith("'") and python_code.endswith("'"):
                            python_code = python_code[1:-1]

                        # Try to compile the Python code to check syntax
                        try:
                            compile(python_code, f'<hook:{hook_id}>', 'exec')
                        except SyntaxError as e:
                            self.fail(f"Hook '{hook_id}' has invalid Python syntax: {e}")

    def test_all_hooks_have_valid_language_field(self):
        """Test that all hooks have a non-empty language field."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    hook_id = hook.get('id', 'unknown')

                    # Language field must exist and have a value
                    self.assertIn('language', hook, f"Hook '{hook_id}' must have a 'language' field")
                    language = hook['language']
                    self.assertIsInstance(language, str, f"Hook '{hook_id}' language must be a string")
                    self.assertTrue(language.strip(), f"Hook '{hook_id}' language must be non-empty")

    def test_all_hooks_have_valid_stages_field(self):
        """Test that all hooks have valid stages field with known stage values."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        # Known valid stages for pre-commit
        known_stages = {
            'commit', 'merge-commit', 'prepare-commit-msg', 'commit-msg',
            'post-commit', 'manual', 'post-checkout', 'post-merge',
            'push', 'post-rewrite'
        }

        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    hook_id = hook.get('id', 'unknown')

                    # If stages field exists, validate it
                    if 'stages' in hook:
                        stages = hook['stages']
                        self.assertIsInstance(stages, list, f"Hook '{hook_id}' stages must be a list")

                        for stage in stages:
                            self.assertIsInstance(stage, str, f"Hook '{hook_id}' stage values must be strings")
                            self.assertIn(stage, known_stages,
                                        f"Hook '{hook_id}' has unknown stage '{stage}'. Known stages: {sorted(known_stages)}")

    def test_all_hooks_have_valid_pass_filenames_field(self):
        """Test that all hooks with pass_filenames field have boolean values."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    hook_id = hook.get('id', 'unknown')

                    # If pass_filenames field exists, validate it's boolean
                    if 'pass_filenames' in hook:
                        pass_filenames = hook['pass_filenames']
                        self.assertIsInstance(pass_filenames, bool,
                                            f"Hook '{hook_id}' pass_filenames must be a boolean, got {type(pass_filenames).__name__}")

    def test_hook_structural_completeness(self):
        """Test comprehensive structural validation of all hook fields."""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        for repo in config['repos']:
            if 'hooks' in repo:
                for hook in repo['hooks']:
                    hook_id = hook.get('id', 'unknown')

                    # Required fields
                    required_fields = ['id', 'entry', 'language']
                    for field in required_fields:
                        self.assertIn(field, hook, f"Hook '{hook_id}' must have required field '{field}'")
                        self.assertTrue(hook[field], f"Hook '{hook_id}' field '{field}' must have a value")

                    # Validate field types
                    self.assertIsInstance(hook['id'], str, f"Hook '{hook_id}' id must be string")
                    self.assertIsInstance(hook['entry'], str, f"Hook '{hook_id}' entry must be string")
                    self.assertIsInstance(hook['language'], str, f"Hook '{hook_id}' language must be string")


if __name__ == '__main__':
    unittest.main()