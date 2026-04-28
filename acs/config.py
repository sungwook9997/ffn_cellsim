"""YAML configuration loader.

Each Tier/Layer is independently togglable via config flags so stages can be
activated incrementally — see CLAUDE.md "Code Conventions".
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: Path | str) -> dict[str, Any]:
    """Load a YAML config file into a plain dict."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"Config root must be a mapping, got {type(data).__name__} in {path}")
    return data
