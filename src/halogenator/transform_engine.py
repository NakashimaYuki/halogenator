# -*- coding: ascii -*-
"""
Shared transformation engine facade.

This module provides a stable Python entrypoint that wraps the existing
`scripts/08_transform_library_v2.py` apply command. It gives downstream code
(including future CLI thin wrappers) a single place to invoke the library
transformation pipeline without duplicating orchestration logic.

Design choices:
- Keep orchestration here, reuse the proven script implementation via a
  subprocess call. This avoids code drift while we gradually migrate the
  heavy logic into `src/halogenator`.
- Return the SUMMARY path to allow callers to chain reporting/visualization.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional


class TransformEngineError(RuntimeError):
    """Raised when the transformation subprocess fails."""


def run_transform_library(
    input_path: Path,
    output_dir: Path,
    xf_name: str,
    xf_config: Path,
    *,
    workers: int = 8,
    batch_size: int = 50_000,
    resume: bool = False,
    python_executable: Optional[str] = None,
    extra_args: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Run the transformation pipeline (08) as a library call.

    Args:
        input_path: Source parquet file.
        output_dir: Destination directory (created if missing).
        xf_name: Transformation name (must exist in xf_config).
        xf_config: YAML configuration describing the transformation.
        workers: Parallel workers for RDKit/processing.
        batch_size: Batch size for streaming parquet.
        resume: Whether to resume a previous run (reuses dedup.db).
        python_executable: Override Python interpreter (defaults to sys.executable).
        extra_args: Additional CLI flags as a dict, e.g. {"--max-sites": "4"}.

    Returns:
        Path to the generated SUMMARY.json for downstream use.

    Raises:
        TransformEngineError: when the subprocess exits non-zero.
    """

    # Locate the repo-level scripts directory (one level above "src")
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "scripts" / "08_transform_library_v2.py"
    if not script_path.exists():
        raise TransformEngineError(f"Script not found: {script_path}")

    python_bin = python_executable or sys.executable

    cmd = [
        python_bin,
        str(script_path),
        "apply",
        "-i",
        str(input_path),
        "-o",
        str(output_dir),
        "--xf-config",
        str(xf_config),
        "--xf-name",
        xf_name,
        "--workers",
        str(workers),
        "--batch-size",
        str(batch_size),
    ]

    if resume:
        cmd.append("--resume")

    if extra_args:
        for key, value in extra_args.items():
            cmd.append(str(key))
            if value is not None and value is not True:
                cmd.append(str(value))

    output_dir.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise TransformEngineError(
            f"Transform failed (code {proc.returncode}).\n"
            f"stdout:\n{proc.stdout[:4000]}\n\nstderr:\n{proc.stderr[:4000]}"
        )

    summary_path = output_dir / "SUMMARY.json"
    if not summary_path.exists():
        # Try to surface more context to callers.
        raise TransformEngineError(
            f"SUMMARY.json not found at {summary_path}. "
            f"stdout tail:\n{proc.stdout[-2000:]}"
        )

    # Quick sanity check: ensure SUMMARY is readable JSON.
    try:
        with summary_path.open("r", encoding="utf-8") as f:
            json.load(f)
    except Exception as exc:  # pragma: no cover - defensive
        raise TransformEngineError(f"SUMMARY.json is not valid JSON: {exc}") from exc

    return summary_path


__all__ = ["run_transform_library", "TransformEngineError"]
