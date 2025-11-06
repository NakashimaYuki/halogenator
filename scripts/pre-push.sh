#!/usr/bin/env bash
set -euo pipefail
echo "[pre-push] 1/5 Syntax check"
git ls-files '*.py' | xargs -I{} python -m py_compile {}
echo "[pre-push] 2/5 ASCII-only check (scripts, src, tests, README)"
python scripts/check_ascii_portable.py
echo "[pre-push] 3/5 Greps"
find src/ -name "*.py" -exec grep -n "from rdkit\|rdkit\.Chem" {} + | cat
echo "[pre-push] 4/5 Smoke tests"
python -m unittest tests.test_granular_qa_json -v
echo "[pre-push] 5/5 Core tests"
python -m unittest tests.test_site_rule_attempts -v