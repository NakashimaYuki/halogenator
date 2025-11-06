# -*- coding: ascii -*-
"""Test that pre-push ASCII checker is portable across platforms."""

import unittest
import tempfile
import os
import subprocess
import pathlib
from unittest.mock import patch


class TestPrePushASCIIPortable(unittest.TestCase):
    """Test that the portable ASCII checker in pre-push.sh works correctly."""
    
    def test_ascii_checker_detects_non_ascii(self):
        """Test that the ASCII checker correctly detects non-ASCII content."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a temporary scripts directory structure
            scripts_dir = os.path.join(temp_dir, 'scripts')
            os.makedirs(scripts_dir)
            
            # Create a file with non-ASCII content
            bad_script = os.path.join(scripts_dir, 'bad.sh')
            with open(bad_script, 'wb') as f:
                # Write content with non-ASCII Chinese characters using bytes
                f.write(b'#!/bin/bash\necho "Hello \xe4\xb8\x96\xe7\x95\x8c"\n')
            
            # Run the Python ASCII checker inline
            check_code = '''
import sys, pathlib
bad = []
for p in pathlib.Path('scripts').rglob('*.sh'):
    try:
        content = p.read_bytes()
        content.decode('ascii')
    except UnicodeDecodeError:
        bad.append(str(p))
if bad:
    print("Non-ASCII found in:", ", ".join(bad))
    sys.exit(1)
'''
            
            # Change to temp directory and run the check
            old_cwd = os.getcwd()
            try:
                os.chdir(temp_dir)
                result = subprocess.run([
                    'python', '-c', check_code
                ], capture_output=True, text=True)
                
                # Should exit with code 1 and report the bad file
                self.assertEqual(result.returncode, 1)
                self.assertIn('Non-ASCII found in:', result.stdout)
                self.assertIn('bad.sh', result.stdout)
            finally:
                os.chdir(old_cwd)
    
    def test_ascii_checker_passes_clean_files(self):
        """Test that the ASCII checker passes files with only ASCII content."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a temporary scripts directory structure
            scripts_dir = os.path.join(temp_dir, 'scripts')
            os.makedirs(scripts_dir)
            
            # Create a file with only ASCII content
            good_script = os.path.join(scripts_dir, 'good.sh')
            with open(good_script, 'w', encoding='ascii') as f:
                f.write('#!/bin/bash\necho "Hello World"\n')
            
            # Run the Python ASCII checker inline
            check_code = '''
import sys, pathlib
bad = []
for p in pathlib.Path('scripts').rglob('*.sh'):
    try:
        content = p.read_bytes()
        content.decode('ascii')
    except UnicodeDecodeError:
        bad.append(str(p))
if bad:
    print("Non-ASCII found in:", ", ".join(bad))
    sys.exit(1)
'''
            
            # Change to temp directory and run the check
            old_cwd = os.getcwd()
            try:
                os.chdir(temp_dir)
                result = subprocess.run([
                    'python', '-c', check_code
                ], capture_output=True, text=True)
                
                # Should exit with code 0 (success)
                self.assertEqual(result.returncode, 0)
                self.assertEqual(result.stdout.strip(), '')
            finally:
                os.chdir(old_cwd)
    
    def test_pre_push_script_uses_portable_solution(self):
        """Test that pre-push.sh uses the portable Python solution instead of grep -P."""
        # Read the pre-push.sh script
        pre_push_path = os.path.join(os.getcwd(), 'scripts', 'pre-push.sh')
        
        if os.path.exists(pre_push_path):
            with open(pre_push_path, 'r') as f:
                content = f.read()
            
            # Should NOT contain grep -P (non-portable)
            self.assertNotIn('grep -P', content, "pre-push.sh should not use non-portable grep -P")
            
            # Should contain the portable Python solution
            self.assertIn('python -', content, "pre-push.sh should use inline Python")
            self.assertIn('pathlib.Path', content, "pre-push.sh should use pathlib for portability")
            self.assertIn('decode(\'ascii\')', content, "pre-push.sh should check ASCII encoding")
    
    def test_python_ascii_checker_is_cross_platform(self):
        """Test that the Python ASCII checker logic works on all platforms."""
        # Test the core logic directly
        import pathlib
        
        # Create temporary content for testing
        with tempfile.NamedTemporaryFile(mode='w+b', suffix='.sh', delete=False) as f:
            # Write non-ASCII content using bytes literal to avoid source file encoding issues
            f.write(b'#!/bin/bash\necho "Test \xe6\xb5\x8b\xe8\xaf\x95"\n')
            temp_path = f.name
        
        try:
            # Test the detection logic
            try:
                content = pathlib.Path(temp_path).read_bytes()
                content.decode('ascii')
                found_non_ascii = False
            except UnicodeDecodeError:
                found_non_ascii = True
            
            self.assertTrue(found_non_ascii, "Should detect non-ASCII content")
            
        finally:
            os.unlink(temp_path)
        
        # Test with ASCII content
        with tempfile.NamedTemporaryFile(mode='w+b', suffix='.sh', delete=False) as f:
            # Write ASCII-only content
            f.write('#!/bin/bash\necho "Test"\n'.encode('ascii'))
            temp_path = f.name
        
        try:
            # Test the detection logic
            try:
                content = pathlib.Path(temp_path).read_bytes()
                content.decode('ascii')
                found_non_ascii = False
            except UnicodeDecodeError:
                found_non_ascii = True
            
            self.assertFalse(found_non_ascii, "Should not detect non-ASCII in ASCII content")
            
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    unittest.main()