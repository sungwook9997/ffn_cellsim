r"""Microtubule aster + MTOC (FF, Warp, µm·pN·s) — the second cytoskeletal network.

Ported from the archived HOOMD module (``archive/hoomd_legacy/cell/microtubules.py``, 1032 lines) onto the live
FF fiber infrastructure: a microtubule is a STIFF hollow tube = an FF fiber whose only changed constant is the
bending rigidity ``EI = KAPPA_MT = 20 pN·µm²`` (Nédélec & Foethke 2007 NJP 9:427 p9; cross-check Gittes 1993
EI ≈ 2.2e-23 N·m² ≈ 300× actin). ``n_mt`` tubes radiate from ONE MTOC (centrosome) along even Fibonacci-sphere
directions; each is a linear bead chain with the FF bending triples (κ=KAPPA_MT) + reshape inextensibility.

Compression load-bearing is EMERGENT from the stiff angle term resisting Euler buckling — NOT a lumped strut.
This is the tensegrity compression element that balances actomyosin tension, positions the nucleus, and bears
the load in confined migration; the audit found it ABSENT from the live FF engine.

Mechanistic, no magic numbers: the ONLY force constant is ``EI = k_B·T·L_p`` (KB-anchored); the per-triple
bending stiffness is DERIVED (``ff/forces_warp._per_triple_alpha``), inextensibility is the NF2007 §5.3 reshape
(no invented stretch modulus). Dynamic instability (growth/catastrophe/rescue, Mitchison-Kirschner 1984) is a
DEFAULT-OFF kinetic layer — a stable aster is the physiological baseline — and is deferred (needs KB DI rates).

Verification gate (analytic ground truth, always-on): a pinned-pinned MT of length L buckles at the Euler load
``F_crit = π²·EI/L²`` (:func:`euler_buckling_load`); at KAPPA_MT and L=5 µm that is 7.90 pN. The persistence
length recovered from thermal fluctuations must return ``L_p = EI/k_BT`` (≈ 4.86 mm, band 1–8 mm).

Sanity Gate (VERIFIED 2026-07-07): the DISCRETE buckling load, computed as the smallest generalized eigenvalue
of the real ``cytosim_bending_kernel`` Hessian (finite-differenced) vs the axial geometric stiffness, CONVERGES
to the analytic Euler load — L=5 µm: err 5.0% (seg 0.25) → 2.5% (0.125) → 1.3% (0.0625), i.e. linear in ℓ₀:
the MT reproduces π²EI/L² as an EMERGENT compression strut, not a lumped one. L_p = 4.67 mm ∈ [1,8] mm.
Aster builds: 40 tubes / 520 nodes / κ=20 pN·µm². (Dynamic instability + MTOC-as-rigid-hub merge deferred.)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ffn_sim.ff.fiber_network import build_fiber_network, FiberNetwork
from ffn_sim.ff.units import KAPPA_MT, KBT


@dataclass
class MicrotubuleAster:
    """A centrosomal MT aster as an FF :class:`FiberNetwork` (nodes = MTOC-adjacent arm beads)."""
    net: FiberNetwork            # the fiber network (κ = KAPPA_MT) — reuse cytosim_bending_kernel + reshape
    mtoc: np.ndarray             # (3,) MTOC / centrosome position [µm]
    n_mt: int                    # number of tubes
    L_mt_um: float               # arm contour length [µm]
    arm_base: np.ndarray         # (n_mt,) node index of each arm's innermost (MTOC-side) bead


def _fibonacci_directions(n: int) -> np.ndarray:
    """``n`` even directions on the unit sphere (deterministic Fibonacci spiral) — the aster arm axes."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)                     # polar
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k                 # golden-angle azimuth
    return np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)


def build_microtubule_aster(centre=(0.0, 0.0, 0.0), *, n_mt: int = 40, L_mt_um: float = 6.0,
                            seg_um: float = 0.5) -> MicrotubuleAster:
    """Build an MT aster: ``n_mt`` stiff tubes radiating from the MTOC at ``centre``.

    Each tube is a chain of beads at spacing ``seg_um`` from just outside the MTOC out to contour length
    ``L_mt_um`` (so it can reach the R≈7.5 µm cortex). Arms are independent FF fibers sharing the MTOC point
    (coincident inner beads = the hub); the MTOC is the ENDPOINT of each arm's first bending triple, so arms
    bend independently. κ = KAPPA_MT on every fiber.
    """
    if n_mt < 1 or L_mt_um <= 0 or seg_um <= 0:
        raise ValueError(f"n_mt={n_mt}, L_mt_um={L_mt_um}, seg_um={seg_um} must be positive")
    c = np.asarray(centre, np.float64)
    n_beads = max(2, int(round(L_mt_um / seg_um)) + 1)     # includes the MTOC-side base bead
    dirs = _fibonacci_directions(n_mt)
    fibers = [c[None, :] + np.arange(n_beads)[:, None] * seg_um * dirs[a][None, :] for a in range(n_mt)]
    net = build_fiber_network(fibers, kappa=KAPPA_MT)
    arm_base = (np.arange(n_mt) * n_beads).astype(np.int64)  # first bead index of each arm
    return MicrotubuleAster(net=net, mtoc=c, n_mt=n_mt, L_mt_um=float(L_mt_um), arm_base=arm_base)


def euler_buckling_load(kappa_pn_um2: float = KAPPA_MT, L_um: float = 5.0) -> float:
    """Pinned-pinned Euler critical buckling load ``F_crit = π²·EI/L²`` [pN] — the analytic MT gate.

    A stiff rod under axial compression is stable below F_crit and buckles above it; the FF aster must reproduce
    this from its emergent bending resistance alone (no strut force). EI = ``kappa_pn_um2`` [pN·µm²]."""
    return float(np.pi ** 2 * kappa_pn_um2 / (L_um ** 2))


def persistence_length_um(kappa_pn_um2: float = KAPPA_MT, kbt_pn_um: float = KBT) -> float:
    """``L_p = EI / k_BT`` [µm] — the thermal-fluctuation gate (MT band ~1–8 mm)."""
    return float(kappa_pn_um2 / kbt_pn_um)


__all__ = ["MicrotubuleAster", "build_microtubule_aster", "euler_buckling_load", "persistence_length_um",
           "_fibonacci_directions"]
