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
        # Stage 1b.b / 1c / 1d.b per-particle state fields (PI viz
        # directive 2026-04-29). All optional — when None, the dataset
        # is omitted, preserving backwards compat with older readers.
        phi: np.ndarray | None = None,
        phi_memory: np.ndarray | None = None,
        c_act: np.ndarray | None = None,
        rho_osm: np.ndarray | None = None,
        gamma_p: np.ndarray | None = None,
        # Stage 1d.c per-particle ECM/protrusion + stress-channel
        # disambiguation (PI directive 2026-04-30). The legacy
        # `stress_tensor` name has been a misleading alias for tau_dev
        # only — going forward we save tau_dev under its true name plus
        # sigma_vol, sigma_active, sigma_total, pressure, dev_norm.
        sigma_vol: np.ndarray | None = None,
        sigma_active: np.ndarray | None = None,
        sigma_total: np.ndarray | None = None,
        pressure: np.ndarray | None = None,
        dev_norm: np.ndarray | None = None,
        traction_ecm: np.ndarray | None = None,
        fa_strength: np.ndarray | None = None,
        protrusion_state: np.ndarray | None = None,
        ecm_signal: np.ndarray | None = None,
        polarity: np.ndarray | None = None,
    ) -> None:
        gid = f"{self._frame_count:05d}"
        g = self._frames.create_group(gid)
        g.attrs["time_star"] = float(time_star)
        kw = {"compression": "gzip", "compression_opts": 4}
        g.create_dataset("position", data=position.astype(np.float32), **kw)
        g.create_dataset("velocity", data=velocity.astype(np.float32), **kw)
        g.create_dataset("deformation_gradient", data=F.astype(np.float32), **kw)
        # PI directive 2026-04-30: save tau_dev under its true name. The
        # legacy `stress_tensor` alias is kept (same data) for backwards
        # read-compat with older notebooks; new code should read `tau_dev`.
        g.create_dataset("tau_dev", data=tau_dev.astype(np.float32), **kw)
        g.create_dataset("stress_tensor", data=tau_dev.astype(np.float32), **kw)
        g.create_dataset("is_boundary", data=is_boundary.astype(np.bool_), **kw)
        # Optional per-particle state fields enable state_overlays viz.
        if phi is not None:
            g.create_dataset("phi", data=phi.astype(np.float32), **kw)
        if phi_memory is not None:
            g.create_dataset("phi_memory", data=phi_memory.astype(np.float32), **kw)
        if c_act is not None:
            g.create_dataset("c_act", data=c_act.astype(np.float32), **kw)
        if rho_osm is not None:
            g.create_dataset("rho_osm", data=rho_osm.astype(np.float32), **kw)
        if gamma_p is not None:
            g.create_dataset("gamma_p", data=gamma_p.astype(np.float32), **kw)
        # Stage 1d.c stress-channel disambiguation.
        if sigma_vol is not None:
            g.create_dataset("sigma_vol", data=sigma_vol.astype(np.float32), **kw)
        if sigma_active is not None:
            g.create_dataset("sigma_active", data=sigma_active.astype(np.float32), **kw)
        if sigma_total is not None:
            g.create_dataset("sigma_total", data=sigma_total.astype(np.float32), **kw)
        if pressure is not None:
            g.create_dataset("pressure", data=pressure.astype(np.float32), **kw)
        if dev_norm is not None:
            g.create_dataset("dev_norm", data=dev_norm.astype(np.float32), **kw)
        # Stage 1d.c ECM + protrusion fields.
        if traction_ecm is not None:
            g.create_dataset("traction_ecm", data=traction_ecm.astype(np.float32), **kw)
        if fa_strength is not None:
            g.create_dataset("fa_strength", data=fa_strength.astype(np.float32), **kw)
        if protrusion_state is not None:
            g.create_dataset("protrusion_state", data=protrusion_state.astype(np.int32), **kw)
        if ecm_signal is not None:
            g.create_dataset("ecm_signal", data=ecm_signal.astype(np.float32), **kw)
        if polarity is not None:
            g.create_dataset("polarity", data=polarity.astype(np.float32), **kw)
        self._frame_count += 1

    def close(self) -> None:
        self.h5.attrs["n_frames"] = self._frame_count
        self.h5.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
