"""Run-provenance helpers.

CLAUDE.md mandates that every run output carry the full configuration plus the
git commit hash. This module produces a compact JSON manifest for that purpose.
"""

from __future__ import annotations

import json
import platform
import socket
import subprocess
import sys
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def _git_commit(repo_root: Path) -> str | None:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def _git_dirty(repo_root: Path) -> bool | None:
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=repo_root,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return bool(out.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def build_manifest(
    config: Mapping[str, Any],
    *,
    repo_root: Path | str = ".",
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble a run-provenance dict ready for JSON serialization."""
    repo_root = Path(repo_root).resolve()
    manifest: dict[str, Any] = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "git_commit": _git_commit(repo_root),
        "git_dirty": _git_dirty(repo_root),
        "config": dict(config),
    }
    if extra:
        manifest["extra"] = dict(extra)
    return manifest


def write_manifest(manifest: Mapping[str, Any], path: Path | str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return path
