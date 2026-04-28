"""Structured logging — one helper, called once per run.

Two handlers: stderr (human-readable) and a per-run rotating file under `logs/`.
Production runs are headless so the file handler is the durable record;
the stderr handler is for live tailing over SSH.
"""

from __future__ import annotations

import logging
import logging.handlers
import sys
from datetime import datetime
from pathlib import Path

_DEFAULT_FMT = "%(asctime)s %(levelname)-7s %(name)s | %(message)s"


def setup_logging(
    log_dir: Path | str = "logs",
    *,
    run_name: str | None = None,
    level: int = logging.INFO,
    fmt: str = _DEFAULT_FMT,
) -> Path:
    """Configure the root logger and return the path of the per-run log file."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = f"_{run_name}" if run_name else ""
    log_path = log_dir / f"{stamp}{suffix}.log"

    formatter = logging.Formatter(fmt)

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)

    root = logging.getLogger()
    root.setLevel(level)
    # Replace handlers idempotently so re-configuring (e.g., in tests) is safe.
    for h in list(root.handlers):
        root.removeHandler(h)
    root.addHandler(file_handler)
    root.addHandler(stream_handler)

    logging.getLogger(__name__).info("Logging initialized; file=%s", log_path)
    return log_path
