"""GPU-resident EnclosedVolumePressure / turgor (compartment-force GPU-main port).

Device-resident twin of ``cortex.enclosed_volume.EnclosedVolumePressure`` (option
b, 2026-06-07 handoff). IDENTICAL physics — the sphere-equivalent osmotic turgor +
bulk-elastic Young-Laplace pressure of enclosed_volume.py — recomputed in cupy via
gpu_local_snapshot + gpu_local_force_arrays, fully masked-vectorised so NO host-sync
occurs in the hot loop. Falls back to the numpy original on a CPU device. NOT wired
into any build — the Lead evaluates + wires (cell-physics lane).
"""

from __future__ import annotations

import math

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cortex.enclosed_volume import (
    ResolvedEnclosedVolume, EnclosedVolumePressure,
)


class EnclosedVolumePressureGPU(md.force.Custom):
    """Device-resident twin of :class:`EnclosedVolumePressure`."""

    def __init__(
        self,
        p: ResolvedEnclosedVolume,
        shell_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(shell_tag_range[0])
        self.tag_end = int(shell_tag_range[1])
        self.last_volume = float("nan")
        self.last_pressure = float("nan")
        self.last_R_mean = float("nan")
        self._cp = None

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        if not isinstance(self._state._simulation.device, hoomd.device.GPU):
            EnclosedVolumePressure.set_forces(self, timestep)
            return
        if self._cp is None:
            import cupy as cp
            self._cp = cp
        cp = self._cp
        p = self.p
        with self._state.gpu_local_snapshot as snap:
            tag = cp.asarray(snap.particles.tag)
            pos = cp.asarray(snap.particles.position, dtype=cp.float64)
            mask = ((tag >= self.tag_start) & (tag < self.tag_end)).astype(cp.float64)
            mcount = cp.maximum(mask.sum(), 1.0)
            centroid = (pos * mask[:, None]).sum(axis=0) / mcount
            dx = pos - centroid
            radii = cp.sqrt((dx * dx).sum(axis=1))
            R_mean = (radii * mask).sum() / mcount
            V = (4.0 / 3.0) * math.pi * R_mean ** 3
            S = 4.0 * math.pi * R_mean * R_mean
            dP = p.turgor_dP0 - p.K_vol * (V - p.V0) / p.V0
            r_safe = cp.where(radii > 0.0, radii, 1.0)
            n_hat = (dx / r_safe[:, None]) * (radii > 0.0)[:, None]
            A_i = S / mcount
            F_vec = ((dP * A_i) * n_hat) * mask[:, None]            # outward turgor
            U_per = ((0.5 * (p.K_vol / p.V0) * (V - p.V0) ** 2) / mcount) * mask
        with self.gpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per
        # no host-sync: self.last_* left stale (off-hot-path accessor if needed)
