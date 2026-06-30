"""Focal-adhesion / traction ANCHOR layer (increment 3) — integrin molecular clutch on a woven network.

Traction is NOT a clean filament weaving (it has no `weave()` manifold); it is an integrin CLUTCH layer
in series with the SF/cortex end nodes (per the increment-3 spec). Each focal adhesion is a population
of integrin–ECM clutches: a spring ``k_int`` from a stress-fiber/cortex actin END node to a FIXED ECM
anchor point on the substrate, bonded by an integrin α5β1–fibronectin CATCH-SLIP bond (Kong 2009,
``hand_kmc.INTEGRIN_A5B1``) — so de-adhesion is EMERGENT from load-dependent bond rupture, never a
binary latch (``feedback-junction-switch-fine-grained``).

Mechanistic reuse, no new force law: the clutch spring is a Hookean link (``network_warp.link_spring``
form); the bond kinetics are the existing two-pathway ``hand_kmc`` catch-slip off-rate. This module is
host-side (numpy) for the prototype; a GPU-resident clutch (fixed-anchor support in relax_on_device) is
the follow-up. Constants (k_int, k_on) are PI-gated Phase defaults (KB-2.4/2.18); the catch-slip rates
are KB-2.5 (Kong 2009) — see the INTEGRIN_A5B1 ⚠️ note (recorded params give F*≈7 pN, KB claim says 30).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ffn_sim.ff.hand_kmc import INTEGRIN_A5B1, HandParams, detach_probability


def fa_end_nodes(cortex, *, axis: int = 2, frac: float = 0.12) -> np.ndarray:
    """The FA-anchored END nodes of a bundle: nodes within ``frac`` of each axial extreme (the two FA
    ends of a stress fiber; barbed-end clusters). Returns node indices."""
    z = cortex.net.pos[:, axis]
    lo, hi = z.min(), z.max()
    span = hi - lo
    return np.where((z < lo + frac * span) | (z > hi - frac * span))[0]


@dataclass(slots=True)
class FAClutchPopulation:
    """A population of integrin clutches: actin END nodes ↔ fixed ECM anchor points (catch-slip bonded)."""

    actin_nodes: np.ndarray          # (M,) cortex/SF node indices (the clutch actin side)
    ecm_anchors: np.ndarray          # (M, 3) FIXED ECM anchor positions [µm] (substrate)
    bound: np.ndarray                # (M,) bool engaged state
    k_int: float                     # clutch spring stiffness [pN/µm]
    rest_um: float                   # clutch rest length (FA bond length l ≈ 0.03 µm)
    integrin: HandParams = field(default_factory=lambda: INTEGRIN_A5B1)


def attach_fa_clutches(cortex, *, end_nodes=None, substrate_offset_um: float = 0.03,
                       k_int: float | None = None, integrin: HandParams = INTEGRIN_A5B1,
                       rng: np.random.Generator | None = None) -> FAClutchPopulation:
    """Decorate ``cortex`` with one integrin clutch per FA end node → a fixed ECM anchor directly
    'below' it (offset ``substrate_offset_um`` along −axis = the FA bond length l ≈ 30 nm). All clutches
    start engaged. ``k_int`` defaults to the INTEGRIN_A5B1 link_k (1e3 pN/µm; PI-gated)."""
    if rng is None:
        rng = np.random.default_rng(0)
    if end_nodes is None:
        end_nodes = fa_end_nodes(cortex)
    pos = cortex.net.pos[end_nodes]
    anchors = pos.copy()
    anchors[:, 2] -= substrate_offset_um                       # ECM anchor 30 nm 'below' the actin node
    return FAClutchPopulation(
        actin_nodes=np.asarray(end_nodes, np.int64), ecm_anchors=anchors,
        bound=np.ones(len(end_nodes), bool), k_int=float(k_int if k_int is not None else integrin.link_k),
        rest_um=substrate_offset_um, integrin=integrin)


def clutch_forces(pop: FAClutchPopulation, pos: np.ndarray) -> np.ndarray:
    """Per-clutch tension magnitude [pN] = k_int·(L − rest) for the BOUND clutches (0 for detached)."""
    d = pos[pop.actin_nodes] - pop.ecm_anchors
    L = np.linalg.norm(d, axis=1)
    f = pop.k_int * np.abs(L - pop.rest_um)
    return np.where(pop.bound, f, 0.0)


def traction_force_nN(pop: FAClutchPopulation, pos: np.ndarray) -> float:
    """Total FA traction = Σ engaged clutch tension [nN] (1 nN = 1000 pN). Literature overlay vs
    Gil-Redondo 2023 ~102 nN per MCF-7 cell (NOT a model-validation gate — H.7 SF-array is HALTED)."""
    return float(clutch_forces(pop, pos).sum() / 1000.0)


def step_clutch_turnover(pop: FAClutchPopulation, pos: np.ndarray, tau: float,
                         rng: np.random.Generator) -> None:
    """One KMC tick (in place): each BOUND clutch detaches with the integrin CATCH-SLIP probability
    1−exp(−τ·k_off(F)) on its current load F; each DETACHED clutch re-attaches with 1−exp(−τ·k_on).
    De-adhesion is emergent from the load-dependent catch-slip rupture (no binary latch)."""
    f = clutch_forces(pop, pos)
    p = pop.integrin
    # bound → detach (load-dependent off-rate)
    bd = pop.bound
    if bd.any():
        p_det = detach_probability(tau, p.off_rate(f[bd]))
        pop.bound[bd] = ~(rng.random(int(bd.sum())) < np.atleast_1d(p_det))
    # detached → re-attach (Poisson k_on)
    ub = ~pop.bound
    if ub.any():
        p_att = 1.0 - np.exp(-tau * p.k_on)
        pop.bound[ub] = rng.random(int(ub.sum())) < p_att
