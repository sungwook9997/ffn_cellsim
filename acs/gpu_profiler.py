"""nvidia-smi sampler — runs in a daemon thread and writes per-run gpu_log.csv.

CLI tooling like `nvidia-smi --query-gpu=utilization.gpu,memory.used,power.draw
--format=csv` is the cheapest way to capture utilization that survives across
GPU vendors / driver upgrades. We sample once a minute by default (configurable)
to avoid measurement overhead — the simulation runs at 1.5 ms/step, so even a
1 Hz sampler would be invisible.

Aggregate statistics (avg / p50 / p95 utilization, peak VRAM, average power)
are computed at shutdown by `summarise()` and rendered into the gate report's
"Performance Profiling" section.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import threading
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def _nvidia_smi_available() -> bool:
    return shutil.which("nvidia-smi") is not None


def _query_once() -> tuple[float, float, float] | None:
    """Single nvidia-smi read. Returns (gpu_util_pct, mem_used_MiB, power_W) for
    GPU 0, or None on failure / when nvidia-smi is missing.
    """
    try:
        out = subprocess.check_output(
            [
                "nvidia-smi",
                "--query-gpu=utilization.gpu,memory.used,power.draw",
                "--format=csv,noheader,nounits",
            ],
            text=True,
            timeout=2,
        )
        first = out.strip().splitlines()[0]
        parts = [p.strip() for p in first.split(",")]
        util = float(parts[0])
        mem = float(parts[1])
        power = float(parts[2]) if parts[2] not in ("", "[N/A]") else float("nan")
        return util, mem, power
    except Exception as exc:
        logger.debug("nvidia-smi query failed: %s", exc)
        return None


class GpuProfiler:
    """Background sampler. Use as a context manager.

    Example:
        with GpuProfiler(out_dir / "gpu_log.csv", interval_s=60) as prof:
            run_simulation()
        stats = prof.summarise()
    """

    def __init__(self, csv_path: Path | str, *, interval_s: float = 60.0):
        self.csv_path = Path(csv_path)
        self.interval_s = float(interval_s)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._available = _nvidia_smi_available()
        self._samples: list[tuple[str, float, float, float]] = []

    def _loop(self) -> None:
        self.csv_path.parent.mkdir(parents=True, exist_ok=True)
        with self.csv_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["timestamp_utc", "gpu_util_pct", "mem_used_MiB", "power_W"])
            # Sample once at startup so a short run still produces ≥ 1 row.
            sample = _query_once()
            if sample is not None:
                ts = datetime.now(UTC).isoformat()
                self._samples.append((ts, *sample))
                w.writerow([ts, *sample])
                fh.flush()
            while not self._stop.wait(self.interval_s):
                sample = _query_once()
                if sample is None:
                    continue
                ts = datetime.now(UTC).isoformat()
                self._samples.append((ts, *sample))
                w.writerow([ts, *sample])
                fh.flush()

    def __enter__(self) -> GpuProfiler:
        if not self._available:
            logger.info("nvidia-smi not on PATH — GPU profiling disabled.")
            return self
        self._thread = threading.Thread(target=self._loop, daemon=True, name="gpu-profiler")
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval_s + 5.0)

    def summarise(self) -> dict[str, Any]:
        if not self._samples:
            return {"available": False, "n_samples": 0}
        utils = sorted(s[1] for s in self._samples)
        mems = [s[2] for s in self._samples]
        powers = [s[3] for s in self._samples if s[3] == s[3]]  # filter NaN

        def percentile(xs: list[float], p: float) -> float:
            if not xs:
                return float("nan")
            k = max(0, min(len(xs) - 1, int(round((p / 100.0) * (len(xs) - 1)))))
            return xs[k]

        return {
            "available": True,
            "n_samples": len(self._samples),
            "avg_gpu_util_pct": sum(utils) / len(utils),
            "gpu_util_p50": percentile(utils, 50),
            "gpu_util_p95": percentile(utils, 95),
            "peak_vram_GB": max(mems) / 1024.0,
            "avg_power_W": (sum(powers) / len(powers)) if powers else float("nan"),
        }
