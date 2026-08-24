r"""Microtubule (MT) compartment for the Active Cell — the tensegrity compression strut (Warp-CUDA).

The `ac/` engine's rule is to REUSE audited `ff/` physics, not re-derive it. Microtubule mechanics already
live, verified, in :mod:`ff.microtubule` (KAPPA_MT = 20 pN·µm², the sole force constant; Nédélec & Foethke
2007 / Gittes 1993 EI ≈ 22 pN·µm²; DOI-sourcing 2026-07-22 confirms 20–22 in-band). This module wraps that
builder into a §1.4 force primitive with the SAME ``accumulate(pos, f)`` contract as
:class:`ac.cell.compartments.MembraneCompartment` / ``NucleusCompartment``, so an MT aster can be composed
into the assembled cell alongside cortex / membrane / nucleus.

An MT is a STIFF hollow tube = an FF fiber whose only changed constant is the bending rigidity κ = KAPPA_MT.
``n_mt`` tubes radiate from one MTOC (centrosome) along Fibonacci-sphere directions; compression load-bearing
is EMERGENT from the stiff angle term resisting Euler buckling — never a lumped strut. The bending force is
the reused ``ff.forces_warp.cytosim_bending_kernel`` with the DERIVED per-triple α = κ/seg³
(``_per_triple_alpha``); inextensibility is the NF2007 reshape the outer solver already applies to every fiber.

Sanity Gate (reuses the ff/ analytic ground truth):
    * Dimensions: κ [pN·µm²]; α = κ/seg³ [pN/µm]; force [pN].
    * Buckling: a pinned-pinned MT of length L buckles at the Euler load F_crit = π²·EI/L²
      (:func:`ff.microtubule.euler_buckling_load`; 7.90 pN at L=5 µm, κ=20) — the discrete bending Hessian
      reproduces it (ff Sanity Gate VERIFIED 2026-07-07, err linear in ℓ₀). ⚠ In the FULL cell, MTs buckle at
      ~100× the free-Euler load with wavelength λ≈3 µm because the surrounding actomyosin reinforces them
      (Brangwynne 2006, DOI 10.1083/jcb.200601060) — the elastic-foundation reinforcement is EMERGENT from the
      composed cortex, not added here; bare Euler is the isolated-MT verification gate only.
    * Persistence: L_p = EI/k_BT ≈ 4.9 mm ∈ [1,8] mm (:func:`ff.microtubule.persistence_length_um`).
    * Boundary: a straight MT has zero discrete curvature ⇒ zero bending force; sign restores toward straight.
    * Residency: triples/alpha/positions are CUDA arrays; ``accumulate`` launches one kernel, no host readback.

GAP (surface to PI, never tuned): the MT COUNT per MCF7 cell is MCF7-absent — proxy ~250–600 MTs/epithelial
cell (order-of-magnitude); ``n_mt`` here is a provisional geometry, not a sourced production count. Length-
dependent softening of short MTs (Pampaloni 2006) and dynein/kinesin point forces (Leidel 2013, ~1–7 pN) are
deferred kinetic layers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from aleph.laws.microtubule import (
    MT_RADIUS_UM,
    build_microtubule_aster,
    euler_buckling_load,
    persistence_length_um,
)
from aleph.laws.units import KAPPA_MT

__all__ = [
    "MicrotubuleCompartment",
    "build_microtubule_compartment",
    "euler_buckling_load",
    "persistence_length_um",
    "KAPPA_MT",
    "MT_RADIUS_UM",
]


@dataclass
class MicrotubuleCompartment:
    """MT-aster force primitive — ``accumulate(pos, f)`` adds the reused Cytosim bending force (§1.4).

    ``tri_d`` are GLOBAL node indices (offset by ``node_off`` at build) so the compartment composes into the
    combined ``[actin | myosin | nucleus | membrane | MT]`` node array; for an isolated test ``node_off=0`` and
    ``pos0`` holds the standalone aster geometry.
    """

    device: str
    node_off: int
    n_nodes: int
    n_mt: int
    L_mt_um: float
    mtoc: npt.NDArray[np.float64]
    pos0: npt.NDArray[np.float64]      # (n_nodes,3) reference aster geometry [µm] (LOCAL, for isolated runs)
    tri_d: wp.array                    # (T,3) int32 GLOBAL bending triples
    alpha_d: wp.array                  # (T,) float64 α = κ/seg³ [pN/µm]
    seg_rest_d: wp.array               # (S,) float64 segment rest lengths [µm] (for NF2007 reshape)
    fiber_off: npt.NDArray[np.int32]   # (n_mt+1,) local fiber offsets
    kappa_pn_um2: float = KAPPA_MT
    radius_um: float = MT_RADIUS_UM
    n_tri: int = field(init=False)

    def __post_init__(self) -> None:
        self.n_tri = int(self.tri_d.shape[0])

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        """Add the reused Cytosim bending force of every MT triple to the global node force."""
        if self.n_tri:
            wp.launch(cytosim_bending_kernel, dim=self.n_tri,
                      inputs=[pos, self.tri_d, self.alpha_d], outputs=[f], device=self.device)

    def euler_buckling_load_pn(self, L_um: float | None = None) -> float:
        """Analytic pinned-pinned Euler load F_crit = π²·EI/L² [pN] (isolated-MT verification gate)."""
        return euler_buckling_load(self.kappa_pn_um2, float(L_um if L_um is not None else self.L_mt_um))

    def persistence_length_um(self) -> float:
        """L_p = EI/k_BT [µm] (thermal-fluctuation gate; MT band 1–8 mm)."""
        return persistence_length_um(self.kappa_pn_um2)


def build_microtubule_compartment(
    *, centre=(0.0, 0.0, 0.0), n_mt: int = 40, reach_R_um: float | None = 7.4, seg_um: float = 0.5,
    node_off: int = 0, device: str | None = None,
) -> MicrotubuleCompartment:
    """Build an MT-aster compartment: ``n_mt`` stiff tubes from the MTOC, tips reaching ``reach_R_um``.

    ``reach_R_um`` (≈ the cortex radius) makes the tubes load-bearing compression struts (tips press on the
    cortex); ``None`` keeps the ff default arm length. ``node_off`` offsets the triples/segments to global
    indices for composition; the CUDA device is required (I0-A).
    """
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError("MicrotubuleCompartment is Warp-CUDA-only (I0-A)")
    dev_s = str(dev)
    aster = build_microtubule_aster(centre=centre, n_mt=n_mt, reach_R_um=reach_R_um, seg_um=seg_um)
    net = aster.net
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    tri = np.ascontiguousarray(net.bend_triples + int(node_off), np.int32)
    # ff FiberNetwork stores SI (metres); the engine runtime is µm. Reference geometry is built in µm by the
    # aster builder, so pos0/seg_rest are taken from the µm-space fibers via fiber offsets + spacing.
    pos0 = np.ascontiguousarray(_aster_positions_um(aster), np.float64)
    seg_rest = np.ascontiguousarray(_aster_seg_rest_um(aster), np.float64)
    with wp.ScopedDevice(dev_s):
        tri_d = wp.array(tri, dtype=wp.int32, device=dev_s)
        alpha_d = wp.array(alpha, dtype=wp.float64, device=dev_s)
        seg_rest_d = wp.array(seg_rest, dtype=wp.float64, device=dev_s)
    return MicrotubuleCompartment(
        device=dev_s, node_off=int(node_off), n_nodes=int(net.n_nodes), n_mt=int(aster.n_mt),
        L_mt_um=float(aster.L_mt_um), mtoc=np.asarray(aster.mtoc, np.float64), pos0=pos0,
        tri_d=tri_d, alpha_d=alpha_d, seg_rest_d=seg_rest_d,
        fiber_off=np.asarray(net.fiber_offsets, np.int32), radius_um=float(aster.radius_um),
    )


def _aster_positions_um(aster) -> npt.NDArray[np.float64]:
    """Reconstruct the aster node positions [µm] from the Fibonacci arms (builder works in µm)."""
    from aleph.laws.microtubule import _fibonacci_directions
    c = np.asarray(aster.mtoc, np.float64)
    dirs = _fibonacci_directions(aster.n_mt)
    off = np.asarray(aster.net.fiber_offsets, np.int64)
    seg = aster.L_mt_um / max(1, (off[1] - off[0]))            # per-arm spacing
    pts = []
    for a in range(aster.n_mt):
        nb = int(off[a + 1] - off[a])
        pts.append(c[None, :] + (1 + np.arange(nb))[:, None] * seg * dirs[a][None, :])
    return np.concatenate(pts, axis=0)


def _aster_seg_rest_um(aster) -> npt.NDArray[np.float64]:
    """Per-segment rest length [µm] for every MT arm (uniform spacing)."""
    off = np.asarray(aster.net.fiber_offsets, np.int64)
    seg = aster.L_mt_um / max(1, (off[1] - off[0]))
    rests = []
    for a in range(aster.n_mt):
        nb = int(off[a + 1] - off[a])
        rests.append(np.full(max(0, nb - 1), seg, np.float64))
    return np.concatenate(rests, axis=0) if rests else np.zeros(0, np.float64)
