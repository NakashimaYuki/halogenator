#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Cross-platform ASCII compliance checker for root documentation files.

Checks specific root files for non-ASCII bytes to ensure they remain
ASCII-only as required by the project standards.
"""

import sys
import pathlib
from typing import List, Tuple

# Root files that must remain ASCII-only
ROOT_FILES = [
    'README.md',
    'CHANGELOG.md',
    'CONTRIBUTING.md',
    '.pre-commit-config.yaml'
]


def first_non_ascii_line(path: pathlib.Path) -> int:
    """Return the 1-based line number of the first non-ASCII byte, or 0 if none."""
    try:
        data = path.read_bytes()
    except OSError as exc:
        print(f"Failed to read {path}: {exc}", file=sys.stderr)
        return -1

    try:
        data.decode('ascii')
        return 0
    except UnicodeDecodeError:
        with path.open('rb') as handle:
            for idx, raw_line in enumerate(handle, start=1):
                try:
                    raw_line.decode('ascii')
                except UnicodeDecodeError:
                    return idx
    return -1


def main() -> int:
    """Check root files for ASCII compliance."""
    failures: List[Tuple[pathlib.Path, int]] = []

    for filename in ROOT_FILES:
        path = pathlib.Path(filename)

        # Skip files that don't exist (optional files)
        if not path.exists():
            continue

        line_no = first_non_ascii_line(path)
        if line_no > 0:
            failures.append((path, line_no))
        elif line_no < 0:
            return 2  # Early exit on IO failures

    if failures:
        print("Non-ASCII found in root files:", file=sys.stderr)
        for path, line_no in failures:
            print(f"{path}:{line_no}: non-ascii byte found", file=sys.stderr)
        return 1

    print('ASCII check passed for root files')
    return 0


if __name__ == '__main__':
    sys.exit(main())