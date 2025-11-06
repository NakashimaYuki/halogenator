=========

# Halogenator Agents Guide (AGENTS.md)

This document is the long-lived operating manual for AI agents working on this
repo. It encodes build/test commands, contracts, code style, safety rules,
and PR workflow. Subpackages may include their own AGENTS.md to locally
override sections (nearest file wins).

All code, filenames, logs, and test output must be ASCII-only.

---

1. Project overview and goals

---

* Purpose: enumerate halogenation products and produce QA stats and reports.
* Key components:

  * src/halogenator/enumerate\_k.py
  * src/halogenator/report.py
  * src/halogenator/cli.py
  * src/halogenator/sites.py
  * src/halogenator/chem\_compat.py (RDKit compatibility shim)
  * tests/\* (unittest suite)
* Primary outputs:

  * products tables
  * QA JSON: v2 with pivots if available; otherwise v1 builder composes
    a stable single-write file with totals + granular slices.

---

2. Golden contracts (do not break)

---

2.1 QA loader acceptance (report.\_load\_qa\_stats\_with\_fallback)

* Strict by default: accept only schemas that contain "attempts".
* Degrade mode: if config\["qa\_loader\_degrade"] == True and only "products"
  is present, accept with attempts=0. Mark qa\_source.degrade=True and
  qa\_source.has\_attempts=False. three\_way\_mutex is skipped in this case.
* Canonicalize to 4 integer keys: attempts, products, no\_product\_matches,
  template\_unsupported.
* Return shape: (qa\_stats: dict, qa\_source: dict) and never mutate inputs.

2.2 Consistency checks in QA report

* Key: "global\_product\_conservation"

  * Sub-buckets:

    * "global\_products": 3D sum over (rule, halogen, k) must equal product\_count
    * "product\_count\_vs\_enum": product\_count must equal enumeration QA "products"
    * "product\_vs\_halogen\_sum": product\_count must equal sum(halogen\_counts)
* Key: "three\_way\_mutex" (separate)

  * Invariant: attempts == products + no\_product\_matches + template\_unsupported
  * Present only when enumeration QA stats include attempts and not in degrade.
  * Expose "qa\_source" and boolean "three\_way\_mutex\_checked" in result.

2.3 QA JSON write contract (single write)

* If payload is v2 (version == "2" or has "pivots"): write as-is, single write.
* If payload is totals-only: writer composes a v1 object with:
  { "version": "1", "total": {4 counters}, "by\_rule": {}, "by\_halogen": {},
  "by\_rule\_halogen": {} } and writes once.
* No double-write hacks.

2.4 Streaming vs stable interface

* enumerate\_with\_stats returns a v2 snapshot (pivots + totals + dedup counters).
* enumerate\_products with return\_qa\_stats:

  * Default per-record is legacy QA dict.
  * Final marker always emits a v2 snapshot.
  * Opt-in v2 per-record via stream\_shape="v2".
* Attempts semantics:
  attempts is incremented once per attempt;
  products increments by 1 when N>0 products for an attempt;
  0 products -> no\_product\_matches += 1;
  unsupported template -> template\_unsupported += 1.

2.5 RDKit compatibility

* All RDKit access goes through src/halogenator/chem\_compat.py:
  from .chem\_compat import Chem, mol\_to\_smiles\_safe, canonical\_rank\_atoms\_safe
* Tests may patch halogenator.chem\_compat.Chem.
* Never import rdkit.Chem directly in other modules.

2.6 Caches and test-order independence

* Ring-label cache must not leak between tests.
* Provide a reset function and call it at top-level enumeration entry points.

---

3. Build, test, and verification

---

3.1 Environment

* Python >= 3.8
* Optional: RDKit. The code must run without RDKit via chem\_compat shims.
* Windows and Linux supported. Prefer PowerShell on Windows and bash on Linux.

3.2 Install

* Use your standard venv. No special package manager is required.

3.3 Test commands

* Full suite:
  python -m unittest -v
* Targeted modules (examples):
  python -m unittest tests.test\_consistency\_qa\_loader\_e2e -v
  python -m unittest tests.test\_cli\_end\_to\_end\_consistency -v

3.4 ASCII check (if script exists)

* ./scripts/check\_ascii.sh
* All code and test content must be ASCII-only.

3.5 Syntax preflight

* Windows:
  git ls-files \*.py | % { python -m py\_compile $\_ }
* Linux:
  python - <<'PY'
  import pathlib, py\_compile
  \[py\_compile.compile(str(p), doraise=True) for p in pathlib.Path('.').rglob('\*.py')]
  print("Syntax OK")
  PY

3.6 PEP 585 generics guard

* No usage of list\[...], dict\[...], set\[...], tuple\[...].
* Grep:
  rg -n "(tuple|list|dict|set)$[^$]+]" src tests || true

3.7 Single-write QA JSON guard

* Grep:
  rg -n "qa\_summary.json" src tests | rg -n "(open(|with\s+open)" || true
* Only the writer function should perform the single write.

---

4. File map and key pointers

---

* enumerate\_k.py

  * QAAggregator, enumerate\_products, enumerate\_with\_stats
  * \_apply\_reaction\_rule, \_iter\_reaction\_mols (must exist as a single helper)
* report.py

  * \_load\_qa\_stats\_with\_fallback (returns qa\_stats, qa\_source)
  * \_validate\_rule\_halogen\_k\_consistency
  * write\_qa\_summary\_json (single-write logic)
* cli.py

  * cmd\_enum orchestrates enumeration and QA JSON writing
* sites.py

  * ring label mapping, cache, compute signature

---

5. Code style and safety

---

* ASCII-only. No Unicode in code, comments, logs, or tests.
* Do not run unknown binaries. Network access is disabled unless explicitly allowed.
* On Windows, avoid "sed" and GNU-only tools; prefer PowerShell equivalents.
* Keep diffs minimal; do not refactor unrelated code in the same patch.

---

6. PR rules and Done definition

---

6.1 PR title format

* <scope>: <short summary>
* Examples:
  enumerate\_k: centralize reaction flatten helper
  report: strict QA loader, add qa\_source flags

6.2 PR content must include

* What changed and why.
* Before vs after behavior.
* Links to tests added/updated.

6.3 Done definition (all must pass)

* python -m py\_compile for all \*.py
* Full test suite green (python -m unittest -v)
* No new skips added
* ASCII-only checks pass (if script exists)
* PEP 585 grep returns no matches
* No double-write of qa\_summary.json
* Contracts in section 2 remain true

---

7. How to work with Codex (prompt templates)

---

Mode: Code
Task: Edit \[file/function] to \[goal]
Pointers: \[exact file paths, function names, grep anchors, short code snippet]
Constraints: ASCII-only; do not use sed on Windows; do not add double-write
Steps:

1. Print a short plan (subtasks and order)
2. Implement minimal changes
3. Run: py\_compile, targeted tests, full suite
4. If failing, apply minimal fix and rerun
   Verification:

* Proceed only if all checks pass
* Output: change summary and diff-ready patch
  Example:
* Locate src/halogenator/enumerate\_k.py:\_apply\_reaction\_rule and replace
  inlined nested loops with \_iter\_reaction\_mols(...) helper.

Mode: Ask

* Use to clarify ambiguity, architecture, or contracts. Keep questions minimal.

---

8. Pre-push checks (recommended script)

---

Linux (tools/prepush\_checks.sh):
\#!/usr/bin/env bash
set -euo pipefail

python - <<'PY'
import pathlib, py\_compile
for p in pathlib.Path('.').rglob('\*.py'):
py\_compile.compile(str(p), doraise=True)
print('Syntax OK')
PY

if rg -n "(tuple|list|dict|set)$[^$]+]" src tests; then
echo "Found PEP 585 generics"; exit 1; fi

rg -n "qa\_summary.json" src tests | rg -n "(open(|with\s+open)" || true

python -m unittest -v

Windows (PowerShell):
git ls-files \*.py | % { python -m py\_compile $\_ }
rg -n "(tuple|list|dict|set)$[^$]+]" src tests
rg -n "qa\_summary.json" src tests | rg -n "(open(|with\s+open)"
python -m unittest -v

---

9. Local overrides

---

Place AGENTS.md in a subpackage directory to override sections locally when
package-specific rules are needed. Nearest AGENTS.md takes precedence.