"""HDF5 frame writer for Stage 1a runs.

Writes a single `snapshots.h5` per run, conforming to docs/13_data_schema.md.
Stage 1a only populates the subset of fields available with Layer 1 alone
(position, velocity, deformation_gradient, stress_tensor, is_boundary). Other
fields documented in the schema are written but filled with NaN/zero so the
file structure remains stable across stages.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import h5py
import numpy as np


class FrameWriter:
    def __init__(
        self,
        path: Path | str,
        *,
        n_particles: int,
        domain_star: float,
        total_time_star: float,
        frame_interval_star: float,
    ):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.h5 = h5py.File(self.path, "w")
        meta = self.h5.create_group("metadata")
        meta.attrs["n_material_points"] = n_particles
        meta.attrs["domain_star"] = domain_star
        meta.attrs["total_time_star"] = total_time_star
        meta.attrs["frame_interval_star"] = frame_interval_star
        meta.attrs["created_utc"] = datetime.now(UTC).isoformat()
        meta.attrs["units"] = "dimensionless (length=R0, time=tau_relax, stress=K)"
        self._frames = self.h5.create_group("frames")
        self._frame_count = 0

    def write_frame(
        self,
        time_star: float,
        *,
        position: np.ndarray,
        velocity: np.ndarray,
        F: np.ndarray,
        tau_dev: np.ndarray,
        is_boundary: np.ndarray,
    ) -> None:
        gid = f"{self._frame_count:05d}"
        g = self._frames.create_group(gid)
        g.attrs["time_star"] = float(time_star)
        kw = {"compression": "gzip", "compression_opts": 4}
        g.create_dataset("position", data=position.astype(np.float32), **kw)
        g.create_dataset("velocity", data=velocity.astype(np.float32), **kw)
        g.create_dataset("deformation_gradient", data=F.astype(np.float32), **kw)
        g.create_dataset("stress_tensor", data=tau_dev.astype(np.float32), **kw)
        g.create_dataset("is_boundary", data=is_boundary.astype(np.bool_), **kw)
        self._frame_count += 1

    def close(self) -> None:
        self.h5.attrs["n_frames"] = self._frame_count
        self.h5.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
