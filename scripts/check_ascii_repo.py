# -*- coding: ascii -*-
"""
Repository-wide ASCII compliance checker.

Scans all text files in the repository to ensure they contain only ASCII characters.
This helps maintain consistency and avoid encoding issues across different platforms.

Usage:
    python scripts/check_ascii_repo.py [root_directory]

By default, checks the current directory recursively.
"""

import sys
import os

# Directories to exclude from scanning
EXCLUDE_DIRS = {
    ".git", ".venv", "venv", "__pycache__", ".mypy_cache", ".pytest_cache",
    ".idea", ".vscode", "node_modules", ".tox", "dist", "build", ".eggs",
    "htmlcov", ".coverage", ".cache", "artifacts", "out", "tmp", "data",
    "halogenator.egg-info", "src/halogenator.egg-info"
}

# File extensions that are binary and should be skipped
BINARY_EXT = {
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".bin", ".pkl", ".gz", ".bz2",
    ".zip", ".7z", ".tar", ".rar", ".exe", ".dll", ".so", ".dylib", ".ico",
    ".woff", ".woff2", ".ttf", ".otf", ".eot", ".mp3", ".mp4", ".avi", ".mov",
    ".xlsx", ".xls", ".doc", ".docx", ".ppt", ".pptx", ".parquet", ".feather",
    ".h5", ".hdf5", ".sqlite", ".db"
}

# File names to exclude entirely
EXCLUDE_FILES = {
    ".coverage", "coverage.xml", "PKG-INFO", "METADATA", "RECORD", "WHEEL",
    "top_level.txt", "dependency_links.txt", "requires.txt"
}

def is_binary(path):
    """Check if a file is binary based on its extension."""
    ext = os.path.splitext(path)[1].lower()
    return ext in BINARY_EXT

def check_ascii(path):
    """Check if a file contains only ASCII characters."""
    try:
        with open(path, "rb") as f:
            data = f.read()
        try:
            text = data.decode("ascii")
            return True, None
        except UnicodeDecodeError as e:
            return False, str(e)
    except (IOError, OSError) as e:
        # Skip files that can't be read (permissions, etc.)
        return True, f"Skipped (read error): {e}"

def main(root="."):
    """Main function to scan repository for non-ASCII files."""
    bad_files = []
    total_checked = 0

    print(f"Scanning {os.path.abspath(root)} for ASCII compliance...")

    for dirpath, dirnames, filenames in os.walk(root):
        # Filter out excluded directories
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]

        for filename in filenames:
            filepath = os.path.join(dirpath, filename)

            # Skip binary files
            if is_binary(filepath):
                continue

            # Skip excluded filenames
            if filename in EXCLUDE_FILES:
                continue

            # Check ASCII compliance
            total_checked += 1
            is_ascii, error = check_ascii(filepath)

            if not is_ascii:
                # Convert absolute path to relative for cleaner output
                rel_path = os.path.relpath(filepath, root)
                bad_files.append((rel_path, error))

    # Report results
    print(f"Checked {total_checked} text files.")

    if bad_files:
        print(f"\nFound {len(bad_files)} files with non-ASCII characters:")
        for file_path, error in bad_files:
            print(f"  - {file_path}: {error}")
        print("\nPlease fix these files to contain only ASCII characters.")
        sys.exit(1)
    else:
        print("ASCII compliance check PASSED - all text files contain only ASCII characters.")

if __name__ == "__main__":
    root_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    main(root_dir)