#!/usr/bin/env python3
# -*- coding: ascii -*-
"""Portable ASCII compliance checker for halogenator project.

Checks src/, scripts/, and docs/ for non-ASCII bytes. Tests are ignored because
historical fixtures contain intentional non-ASCII characters. Use this script
in CI to guard newly written code and tooling.
"""

import sys
import pathlib
from typing import Iterable

CHECK_DIRECTORIES = ('src', 'scripts', 'docs')
BINARY_SUFFIXES = {
    '.png', '.jpg', '.jpeg', '.gif', '.bmp', '.ico', '.pdf', '.zip', '.gz', '.bz2',
    '.xz', '.7z', '.rar', '.tar', '.pyc', '.pyo', '.exe', '.dll', '.so', '.dylib',
    '.xlsx', '.xls', '.ppt', '.pptx', '.doc', '.docx'
}
IGNORED_NAMES = {'.git', '__pycache__'}


def iter_files() -> Iterable[pathlib.Path]:
    """Yield candidate files to scan for ASCII compliance."""
    for directory in CHECK_DIRECTORIES:
        base = pathlib.Path(directory)
        if not base.exists():
            continue
        for path in base.rglob('*'):
            if path.is_dir():
                continue
            if any((part in IGNORED_NAMES) or part.endswith('.egg-info') for part in path.parts):
                continue
            if path.suffix.lower() in BINARY_SUFFIXES:
                continue
            yield path


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
    failures = []

    for path in iter_files():
        line_no = first_non_ascii_line(path)
        if line_no > 0:
            failures.append((path, line_no))
        elif line_no < 0:
            return 2  # Early exit on IO failures

    if failures:
        for path, line_no in failures:
            print(f"{path}:{line_no}: non-ascii byte found")
        return 1

    print('ASCII check passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
