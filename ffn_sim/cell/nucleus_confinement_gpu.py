"""GPU-resident NucleusConfinement (compartment-force GPU-main port PROTOTYPE).

Sub-session prototype (option b of the 2026-06-07 handoff): a device-resident
drop-in for ``cell.nucleus.NucleusConfinement`` whose ``set_forces`` stays on the
GPU (``gpu_local_snapshot`` + cupy + ``gpu_local_force_arrays``) instead of the
per-step ``cpu_local_snapshot`` host-sync that the Lead's full-cell profile found
dominating the real-op-point step (the 2.43× ceiling). IDENTICAL physics — the
two-regime radial bilinear confinement law of nucleus.py, recomputed in cupy.

Fully masked-vectorised (no boolean fancy-index, no ``int(mask.sum())``) so NOT
ONE host-sync occurs in the hot loop. Falls back to the numpy path on a CPU
device (bit-for-bit the original) so it is safe to wire unconditionally.

NOT wired into any build — the Lead evaluates + wires (their cell-physics lane).
Validated by ``scripts/h7_native_nucleus_force_parity.py`` against the original.
"""

from __future__ import annotations

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cell.nucleus import ResolvedNucleus, NucleusConfinement


class NucleusConfinementGPU(md.force.Custom):
    """Device-resident twin of :class:`cell.nucleus.NucleusConfinement`.

    Same constructor + same bilinear law; ``set_forces`` runs entirely on the
    GPU when the device is a ``hoomd.device.GPU`` (cupy), else falls back to the
    original numpy host path.
    """

    def __init__(
        self,
        p: ResolvedNucleus,
        nucleus_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(nucleus_tag_range[0])
        self.tag_end = int(nucleus_tag_range[1])
        if self.tag_end < self.tag_start:
            raise ValueError("nucleus_tag_range must satisfy end >= start")
        self.last_n: int = 0
        self.last_max_strain: float = float("nan")
        self.last_cross_section_area: float = float("nan")
        self._cp = None  # cupy module, resolved on first GPU step

    def _bilinear(self, xp, d):
        """Vectorised two-regime radial law (xp = numpy or cupy). Returns F_mag, U."""
        p = self.p
        ad = xp.abs(d)
        sign = xp.sign(d)
        interior = ad <= p.d_knee
        F_mag = xp.where(
            interior,
            -p.k_chrom * d,
            -sign * (p.F_knee + (p.k_chrom + p.k_lamin) * (ad - p.d_knee)),
        )
        u_knee = 0.5 * p.k_chrom * p.d_knee * p.d_knee
        e = ad - p.d_knee
        U = xp.where(
            interior,
            0.5 * p.k_chrom * d * d,
            u_knee + p.F_knee * e + 0.5 * (p.k_chrom + p.k_lamin) * e * e,
        )
        return F_mag, U

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        on_gpu = isinstance(self._state._simulation.device, hoomd.device.GPU)
        if not on_gpu:
            # CPU fallback: defer to the original (bit-for-bit) to stay safe.
            NucleusConfinement.set_forces(self, timestep)
            return
        if self._cp is None:
            import cupy as cp
            self._cp = cp
        cp = self._cp

        with self._state.gpu_local_snapshot as snap:
            tag = cp.asarray(snap.particles.tag)
            pos = cp.asarray(snap.particles.position, dtype=cp.float64)
            mask = ((tag >= self.tag_start) & (tag < self.tag_end)).astype(cp.float64)
            mcount = mask.sum()
            # masked centroid (device-only; mcount>0 in any real nucleus build)
            centroid = (pos * mask[:, None]).sum(axis=0) / cp.maximum(mcount, 1.0)
            dx = pos - centroid
            r = cp.sqrt((dx * dx).sum(axis=1))
            r_safe = cp.where(r > 0.0, r, 1.0)
            r_hat = dx / r_safe[:, None]
            r_hat = r_hat * (r > 0.0)[:, None]
            d = r - self.p.R_nuc
            F_mag, U = self._bilinear(cp, d)
            F_vec = (F_mag[:, None] * r_hat) * mask[:, None]
            U_per = U * mask

        with self.gpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per
        # NOTE: deliberately NO host-sync here (no int(mcount.get())) — the hot
        # loop stays fully on-device. self.last_* diagnostics are left stale; a
        # consumer that needs them calls a separate (off-hot-path) accessor.
