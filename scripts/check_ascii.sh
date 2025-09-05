#!/usr/bin/env bash
# -*- coding: ascii -*-
set -euo pipefail
# Fail if any non-ASCII byte is found in tracked source files.
git ls-files | grep -E '^(src|tests|scripts)/.*\.(py|sh)$' \
| xargs -I{} sh -c "LC_ALL=C grep -n -P '[^\x00-\x7F]' '{}' && echo 'Non-ASCII in {}' && exit 1" \
|| exit 0