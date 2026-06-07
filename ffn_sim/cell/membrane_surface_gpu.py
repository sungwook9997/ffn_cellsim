"""GPU-resident MembraneSurfaceTension (compartment-force GPU-main port).

Device-resident twin of ``cell.membrane_surface.MembraneSurfaceTension`` (option
b, 2026-06-07 handoff). IDENTICAL physics — the sphere-equivalent Young-Laplace
surface tension of membrane_surface.py — recomputed in cupy via gpu_local_snapshot
+ gpu_local_force_arrays, fully masked-vectorised so NO host-sync occurs in the
hot loop (the per-step cpu_local_snapshot of the original costs ~1.4 ms/force at
the full-cell N — see h7_native_nucleus_force_parity). Falls back to the numpy
original on a CPU device. NOT wired into any build — the Lead evaluates + wires.
"""

from __future__ import annotations

import math

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cell.membrane_surface import (
    ResolvedMembraneSurface, MembraneSurfaceTension,
)


class MembraneSurfaceTensionGPU(md.force.Custom):
    """Device-resident twin of :class:`MembraneSurfaceTension`."""

    def __init__(
        self,
        p: ResolvedMembraneSurface,
        shell_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(shell_tag_range[0])
        self.tag_end = int(shell_tag_range[1])
        self.last_area = float("nan")
        self.last_tension = float("nan")
        self.last_pressure = float("nan")
        self.last_R_mean = float("nan")
        self._cp = None

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        if not isinstance(self._state._simulation.device, hoomd.device.GPU):
            MembraneSurfaceTension.set_forces(self, timestep)
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
            S = 4.0 * math.pi * R_mean * R_mean
            gamma_tot = p.gamma_mem + p.K_A * (S - p.A0) / p.A0
            dP = 2.0 * gamma_tot / cp.maximum(R_mean, 1e-300)
            r_safe = cp.where(radii > 0.0, radii, 1.0)
            n_hat = (dx / r_safe[:, None]) * (radii > 0.0)[:, None]
            A_i = S / mcount
            F_vec = ((-(dP * A_i)) * n_hat) * mask[:, None]          # inward (compressive)
            U_const = cp.maximum(p.gamma_mem * (S - p.A0), 0.0)
            U_area = 0.5 * p.K_A * (S - p.A0) ** 2 / p.A0
            U_per = ((U_const + U_area) / mcount) * mask
        with self.gpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per
        # no host-sync: self.last_* left stale (off-hot-path accessor if needed)
