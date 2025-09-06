#!/usr/bin/env bash
# -*- coding: ascii -*-
set -euo pipefail
# Fail if any non-ASCII byte is found in tracked source files.

# Check source code files in src/, tests/, scripts/, and configs/
echo "Checking ASCII compliance in source files..."

found_non_ascii=false

while read -r file; do
  if [ -f "$file" ]; then
    if LC_ALL=C grep -q $'[\x80-\xFF]' "$file" 2>/dev/null; then
      echo "ERROR: Non-ASCII characters found in $file"
      found_non_ascii=true
    fi
  fi
done < <(git ls-files | grep -E '^(src|tests|scripts|configs)/.*\.(py|sh|yaml)$')

if [ "$found_non_ascii" = true ]; then
  echo "ASCII check FAILED: Non-ASCII characters found in source files"
  exit 1
else
  echo "ASCII check PASSED: All source files contain only ASCII characters"
fi