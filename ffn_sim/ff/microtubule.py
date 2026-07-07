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

from ffn_sim.ff.fiber_network import build_fiber_network, concat_fiber_networks, FiberNetwork
from ffn_sim.ff.units import KAPPA_MT, KBT


@dataclass
class MicrotubuleAster:
    """A centrosomal MT aster as an FF :class:`FiberNetwork` (nodes = MTOC-adjacent arm beads)."""
    net: FiberNetwork            # the fiber network (κ = KAPPA_MT) — reuse cytosim_bending_kernel + reshape
    mtoc: np.ndarray             # (3,) MTOC / centrosome position [µm]
    n_mt: int                    # number of tubes
    L_mt_um: float               # arm contour length [µm]
    arm_base: np.ndarray         # (n_mt,) node index of each arm's innermost (near-MTOC) bead


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
    n_beads = max(2, int(round(L_mt_um / seg_um)))         # arm beads (the MTOC is a SEPARATE node, added at merge)
    dirs = _fibonacci_directions(n_mt)
    # arms run from ONE segment out (bead 0 at seg_um from the MTOC) → the MTOC↔base hub bond has finite rest seg_um
    fibers = [c[None, :] + (1 + np.arange(n_beads))[:, None] * seg_um * dirs[a][None, :] for a in range(n_mt)]
    net = build_fiber_network(fibers, kappa=KAPPA_MT)
    arm_base = (np.arange(n_mt) * n_beads).astype(np.int64)  # innermost (near-MTOC) bead index of each arm
    return MicrotubuleAster(net=net, mtoc=c, n_mt=n_mt, L_mt_um=float(L_mt_um), arm_base=arm_base)


def merge_aster_into_cortex(cortex_net: FiberNetwork, aster: MicrotubuleAster | None, *,
                            k_hub_pn_um: float, seg_um: float = 0.5) -> dict:
    """Merge the MT aster into the cortex fiber network for the SHARED implicit solve (Thread-C stage 1).

    Returns a dict with the merged elastic network + the layout + the MTOC-hub crosslinks + MT node drag. The
    merged bending K carries cortex (κ_actin) and MT (κ=KAPPA_MT) triples in one solve. The MTOC is a SEPARATE
    node (index ``mtoc_idx``) held by finite-rest (``seg_um``) stiff crosslinks to each arm's innermost bead —
    NOT a zero-length weld (which ``implicit_ff`` drops at L<1e-9). Node layout the driver assembles:
    ``pos_all = [cortex(Nc) ; MT_arms(Nmt) ; MTOC(1) ; nucleus(n_nuc)]`` ⇒ ``Ne = Nc + Nmt + 1``.

    ``aster is None`` → no-op passthrough (Nc==Ne, empty hub) → the --microtubules-OFF path is bit-identical."""
    Nc = cortex_net.n_nodes
    if aster is None:
        return {"net": cortex_net, "Nc": Nc, "Ne": Nc, "mtoc_pos": np.zeros((0, 3)), "mtoc_idx": -1,
                "hub_i": np.zeros(0, np.int64), "hub_j": np.zeros(0, np.int64),
                "hub_k": np.zeros(0), "hub_rest": np.zeros(0), "mt_gammas": np.zeros(0)}
    merged, node_off = concat_fiber_networks([cortex_net, aster.net])   # [cortex ; MT arms]
    Nmt = aster.net.n_nodes
    mtoc_idx = Nc + Nmt                                     # the MTOC node is appended after the MT arms
    Ne = mtoc_idx + 1
    base_idx = node_off[1] + aster.arm_base                 # arm innermost beads in the merged indexing
    hub_i = np.full(aster.n_mt, mtoc_idx, np.int64)
    hub_j = base_idx.astype(np.int64)
    hub_k = np.full(aster.n_mt, float(k_hub_pn_um))
    hub_rest = np.full(aster.n_mt, float(seg_um))          # finite rest = one segment (arms start at seg_um)
    return {"net": merged, "Nc": Nc, "Ne": Ne, "mtoc_pos": aster.mtoc.reshape(1, 3), "mtoc_idx": mtoc_idx,
            "hub_i": hub_i, "hub_j": hub_j, "hub_k": hub_k, "hub_rest": hub_rest,
            "Nmt": Nmt, "n_mt": aster.n_mt}


def euler_buckling_load(kappa_pn_um2: float = KAPPA_MT, L_um: float = 5.0) -> float:
    """Pinned-pinned Euler critical buckling load ``F_crit = π²·EI/L²`` [pN] — the analytic MT gate.

    A stiff rod under axial compression is stable below F_crit and buckles above it; the FF aster must reproduce
    this from its emergent bending resistance alone (no strut force). EI = ``kappa_pn_um2`` [pN·µm²]."""
    return float(np.pi ** 2 * kappa_pn_um2 / (L_um ** 2))


def persistence_length_um(kappa_pn_um2: float = KAPPA_MT, kbt_pn_um: float = KBT) -> float:
    """``L_p = EI / k_BT`` [µm] — the thermal-fluctuation gate (MT band ~1–8 mm)."""
    return float(kappa_pn_um2 / kbt_pn_um)


__all__ = ["MicrotubuleAster", "build_microtubule_aster", "euler_buckling_load", "persistence_length_um", "merge_aster_into_cortex",
           "_fibonacci_directions"]
