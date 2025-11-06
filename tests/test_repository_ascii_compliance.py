# -*- coding: ascii -*-
"""Test that all repository text files contain only ASCII characters."""

import unittest
import pathlib
import os


class TestRepositoryASCIICompliance(unittest.TestCase):
    """Test that all text files in the repository contain only ASCII characters."""

    def test_repository_has_no_non_ascii_text_files(self):
        """Test that all text files contain only ASCII characters."""
        project_root = pathlib.Path(__file__).parent.parent

        # File extensions to check
        text_extensions = {'.py', '.md', '.yml', '.yaml', '.sh', '.txt', '.toml'}

        # Whitelist paths to skip (binary files, licenses, etc.)
        skip_paths = {
            '.git',
            '__pycache__',
            '.pytest_cache',
            'node_modules',
            '.venv',
            'venv',
            'env',
            'data',
            'benchmarks',
            '.claude',
            'out',
            'tmp',
            'artifacts',
            'dist',
            'build',
            '.eggs',
            'htmlcov',
            '.coverage',
            '.cache'
        }

        non_ascii_files = []

        def should_skip_path(path):
            """Check if a path should be skipped."""
            parts = path.parts
            return any(skip_part in parts for skip_part in skip_paths)

        # Check all text files
        for file_path in project_root.rglob('*'):
            if file_path.is_file() and not should_skip_path(file_path):
                if file_path.suffix in text_extensions or file_path.name in {'README', 'LICENSE', 'CHANGELOG'}:
                    try:
                        # Try to read as bytes and decode as ASCII
                        content = file_path.read_bytes()
                        content.decode('ascii')
                    except UnicodeDecodeError as e:
                        relative_path = file_path.relative_to(project_root)
                        non_ascii_files.append({
                            'file': str(relative_path),
                            'error': str(e)
                        })
                    except (OSError, IOError):
                        # Skip files that can't be read (permissions, etc.)
                        continue

        # Report any non-ASCII files found
        if non_ascii_files:
            error_msg = "Non-ASCII characters found in the following files:\n"
            for file_info in non_ascii_files:
                error_msg += f"  - {file_info['file']}: {file_info['error']}\n"
            error_msg += "\nAll text files must contain only ASCII characters. Please remove or replace non-ASCII characters."
            self.fail(error_msg)

    def test_ascii_checking_script_exists(self):
        """Test that ASCII checking script exists and is executable."""
        project_root = pathlib.Path(__file__).parent.parent

        # Check for pre-push script
        pre_push_script = project_root / 'scripts' / 'pre-push.sh'
        self.assertTrue(pre_push_script.exists(), "pre-push.sh script should exist")

        # Check that it contains ASCII checking logic
        if pre_push_script.exists():
            content = pre_push_script.read_text(encoding='utf-8')
            self.assertIn('ASCII', content, "pre-push.sh should contain ASCII checking logic")

    def test_pre_commit_config_has_ascii_check(self):
        """Test that pre-commit configuration includes ASCII checking."""
        project_root = pathlib.Path(__file__).parent.parent

        # Check for pre-commit config
        pre_commit_config = project_root / '.pre-commit-config.yaml'
        if pre_commit_config.exists():
            content = pre_commit_config.read_text(encoding='utf-8')
            self.assertIn('ascii', content.lower(), ".pre-commit-config.yaml should include ASCII checking")


if __name__ == '__main__':
    unittest.main()