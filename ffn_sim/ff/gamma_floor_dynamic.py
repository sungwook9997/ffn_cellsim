"""γ-floor with a DYNAMIC Hand kinetic layer (Stage 6e) — regime-B event-driven on the solve.

Stage 6d measured γ on a *static* link set. This wires the `hand_kmc` layer in dynamically: links
(crosslinkers + myosin) attach / detach over KMC ticks interleaved with the mechanical equilibrium
solve — the regime-B "event-driven on the mechanical solve" of ENGINE.md §1. Each tick:

  1. settle the mechanical equilibrium at the current bound-link set (bending + active links + turgor);
  2. read the load (tension) on each bound link from the mechanics;
  3. KMC: bound links detach at their force-dependent off-rate (Bell for α-actinin/NMIIA, Pereverzev
     catch–slip for filamin); free candidate links (cross-fiber pairs currently within the mesoscale
     reach) attach at rate k_on;
  4. measure γ.

A link is modelled as a single Hand toggled bound/unbound over a precomputed candidate-pair pool
(the geometric site search done once; HandPopulation carries the bound state + per-link off-rate).

Validation target (ENGINE.md §4 / the archived finding): dynamic turnover should give ≈ the static-
link γ — i.e. **binding kinetics is NOT the lever for the γ-floor** (the floor is force-magnitude /
transmission). Reproducing that here closes the loop on "is the MD-free γ kinetics-sensitive?".
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ffn_sim.ff.constraints import project_constraint_forces, reshape, segment_axial_tension
from ffn_sim.ff.cortex_assembly import CortexParams, build_cortex_network
from ffn_sim.ff.fiber_network import FiberNetwork
from ffn_sim.ff.forces_warp import make_bending_force_fn
from ffn_sim.ff.gamma_estimator import method_of_planes_gamma
from ffn_sim.ff.gamma_floor import (
    NMIIA_MINIFIL_STALL_PN,
    TURGOR_DP0,
    _link_spring_force,
    _myosin_force,
    _cross_fiber_pairs,
    gamma_passive_young_laplace,
    mesoscale_reach,
)
from ffn_sim.ff.hand_kmc import ALPHA_ACTININ, FILAMIN, NMIIA_MYOSIN, HandPopulation


@dataclass(slots=True)
class DynamicCortex:
    """A cortex with dynamic crosslinker + myosin Hand populations over candidate cross-fiber pairs."""

    net: FiberNetwork
    R_um: float
    reach: float
    # crosslinker candidates (pairs) + per-candidate type/stiffness/rest + Hand population
    xl_pairs: np.ndarray            # (Cxl, 2) candidate node pairs
    xl_is_alpha: np.ndarray         # (Cxl,) bool — α-actinin vs filamin
    xl_k: np.ndarray                # (Cxl,) link stiffness
    xl_rest: np.ndarray             # (Cxl,) rest length (formation distance)
    xl_pop: HandPopulation
    # myosin candidates
    myo_pairs: np.ndarray           # (Cmyo, 2)
    myo_pop: HandPopulation


def build_dynamic_cortex(params: CortexParams | None = None, *, n_filaments: int = 100,
                         n_xl_cand: int = 800, n_myo_cand: int = 400,
                         alpha_fraction: float = 0.30,
                         rng: np.random.Generator | None = None) -> DynamicCortex:
    """Assemble a cortex + candidate-pair pools for the dynamic Hand layer (all start UNBOUND)."""
    if params is None:
        params = CortexParams()
    if rng is None:
        rng = np.random.default_rng(0)
    net, _ = build_cortex_network(params, rng=rng, n_filaments=n_filaments)
    reach = mesoscale_reach(params.R_um, n_filaments)

    xl_pairs = _cross_fiber_pairs(net, reach, n_xl_cand, rng)
    cxl = xl_pairs.shape[0]
    is_alpha = rng.random(cxl) < alpha_fraction
    xl_k = np.where(is_alpha, ALPHA_ACTININ.link_k, FILAMIN.link_k)
    xl_rest = (np.linalg.norm(net.pos[xl_pairs[:, 1]] - net.pos[xl_pairs[:, 0]], axis=1)
               if cxl else np.zeros(0))

    used = {(int(a), int(b)) for a, b in xl_pairs}
    myo_pairs = _cross_fiber_pairs(net, reach, n_myo_cand, rng, exclude=used)

    # one Hand per candidate; α-actinin params stand in for the crosslinker population off-rate
    # dispatch is per-candidate (handled in the tick), so the population's params are a placeholder.
    xl_pop = HandPopulation.empty(cxl, ALPHA_ACTININ)
    myo_pop = HandPopulation.empty(myo_pairs.shape[0], NMIIA_MYOSIN)
    return DynamicCortex(net=net, R_um=params.R_um, reach=reach, xl_pairs=xl_pairs,
                         xl_is_alpha=is_alpha, xl_k=xl_k, xl_rest=xl_rest, xl_pop=xl_pop,
                         myo_pairs=myo_pairs, myo_pop=myo_pop)


def _active_links(dc: DynamicCortex):
    """Bound crosslinker + myosin link arrays (endpoints, k, rest) for the current state."""
    xb = dc.xl_pop.bound
    mb = dc.myo_pop.bound
    return (dc.xl_pairs[xb], dc.xl_k[xb], dc.xl_rest[xb], dc.myo_pairs[mb])


def _external_force_dynamic(dc: DynamicCortex, pos: np.ndarray, bending_fn, f_myo: float,
                            *, turgor: bool) -> np.ndarray:
    from ffn_sim.ff.gamma_floor import _turgor_force, turgor_pressure
    N = dc.net.n_nodes
    F = np.asarray(bending_fn(pos.reshape(-1)), dtype=np.float64).reshape(N, 3)
    xl_p, xl_k, xl_rest, myo_p = _active_links(dc)
    if xl_p.shape[0]:
        _link_spring_force(pos, xl_p[:, 0], xl_p[:, 1], xl_k, xl_rest, F)
    if myo_p.shape[0]:
        _myosin_force(pos, myo_p[:, 0], myo_p[:, 1], f_myo, F)
    if turgor:
        dP, R_mean = turgor_pressure(dc, pos)
        _turgor_force(dc, pos, dP, R_mean, F)
    return F


def _settle(dc: DynamicCortex, f_myo: float, *, n_steps: int, turgor: bool) -> None:
    """Projected overdamped settle of the current bound-link configuration (in place on net.pos)."""
    net = dc.net
    bending_fn = make_bending_force_fn(net)
    seg = float(net.seg_rest.mean())
    k_max = max(float(net.kappa.max()) / seg**3, float(dc.xl_k.max()) if dc.xl_k.size else 0.0)
    dt_mu = 0.1 / k_max
    x = net.pos.reshape(-1, 3).copy()
    for step in range(n_steps):
        net.pos = x
        F = _external_force_dynamic(dc, x, bending_fn, f_myo, turgor=turgor)
        x = x + dt_mu * project_constraint_forces(net, F)
        if (step + 1) % 25 == 0:
            net.pos = x
            x = reshape(net, n_iter=2)
    net.pos = reshape(net, n_iter=4)


def _link_loads(dc: DynamicCortex, f_myo: float):
    """Tension magnitude on each BOUND crosslinker / myosin link at the current geometry [pN]."""
    pos = dc.net.pos
    xb = dc.xl_pop.bound
    if xb.any():
        p = dc.xl_pairs[xb]
        L = np.linalg.norm(pos[p[:, 1]] - pos[p[:, 0]], axis=1)
        xl_load = np.abs(dc.xl_k[xb] * (L - dc.xl_rest[xb]))
    else:
        xl_load = np.zeros(0)
    myo_load = np.full(int(dc.myo_pop.bound.sum()), abs(f_myo))   # myosin carries ~its prestress
    return xl_load, myo_load


def _kmc_tick(dc: DynamicCortex, f_myo: float, tau: float, rng: np.random.Generator) -> None:
    """One KMC sub-step: detach bound links (force-dependent), attach free in-range candidates."""
    from ffn_sim.ff.hand_kmc import attach_probability, detach_probability

    pos = dc.net.pos
    # ---- crosslinkers ----
    cxl = dc.xl_pairs.shape[0]
    if cxl:
        d = np.linalg.norm(pos[dc.xl_pairs[:, 1]] - pos[dc.xl_pairs[:, 0]], axis=1)
        # detach bound: per-candidate off-rate (α-actinin Bell / filamin catch–slip) at its load
        bound = dc.xl_pop.bound
        if bound.any():
            L = d[bound]
            load = np.abs(dc.xl_k[bound] * (L - dc.xl_rest[bound]))
            is_a = dc.xl_is_alpha[bound]
            p_off = np.where(is_a, ALPHA_ACTININ.off_rate(load), FILAMIN.off_rate(load))
            det = rng.random(int(bound.sum())) < detach_probability(tau, p_off)
            idx = np.where(bound)[0][det]
            dc.xl_pop.bound[idx] = False
        # attach free in-range (rest reset to current distance at (re)binding → force-free)
        free = (~dc.xl_pop.bound) & (d <= dc.reach)
        roll = rng.random(cxl) < attach_probability(tau, ALPHA_ACTININ.k_on)
        new = free & roll
        dc.xl_pop.bound[new] = True
        dc.xl_rest[new] = d[new]
    # ---- myosin ----
    cmyo = dc.myo_pairs.shape[0]
    if cmyo:
        dm = np.linalg.norm(pos[dc.myo_pairs[:, 1]] - pos[dc.myo_pairs[:, 0]], axis=1)
        bound = dc.myo_pop.bound
        if bound.any():
            p_off = NMIIA_MYOSIN.off_rate(np.full(int(bound.sum()), abs(f_myo)))
            det = rng.random(int(bound.sum())) < detach_probability(tau, p_off)
            idx = np.where(bound)[0][det]
            dc.myo_pop.bound[idx] = False
        free = (~dc.myo_pop.bound) & (dm <= dc.reach)
        roll = rng.random(cmyo) < attach_probability(tau, NMIIA_MYOSIN.k_on)
        dc.myo_pop.bound[free & roll] = True


def measure_gamma_dynamic(dc: DynamicCortex, f_myo: float, *, n_planes: int = 50) -> dict:
    """γ over the currently-bound links (actin axial + crosslinker + myosin) + passive turgor."""
    net = dc.net
    pos = net.pos
    bending_fn = make_bending_force_fn(net)
    F = _external_force_dynamic(dc, pos, bending_fn, f_myo, turgor=False)
    seg_tau = segment_axial_tension(net, F)
    seg = net.segments
    rA, rB, T = [pos[seg[:, 0]]], [pos[seg[:, 1]]], [seg_tau]
    xl_p, xl_k, xl_rest, myo_p = _active_links(dc)
    if xl_p.shape[0]:
        L = np.linalg.norm(pos[xl_p[:, 1]] - pos[xl_p[:, 0]], axis=1)
        rA.append(pos[xl_p[:, 0]]); rB.append(pos[xl_p[:, 1]]); T.append(xl_k * (L - xl_rest))
    if myo_p.shape[0]:
        rA.append(pos[myo_p[:, 0]]); rB.append(pos[myo_p[:, 1]])
        T.append(np.full(myo_p.shape[0], f_myo))
    centre = pos.mean(axis=0)
    g = method_of_planes_gamma(np.concatenate(rA), np.concatenate(rB), np.concatenate(T),
                               dc.R_um, n_planes=n_planes, centre=centre)
    return {"gamma_active": g, "gamma_total": g,
            "gamma_passive": gamma_passive_young_laplace(TURGOR_DP0, dc.R_um),
            "n_xl_bound": int(dc.xl_pop.bound.sum()), "n_myo_bound": int(dc.myo_pop.bound.sum())}


def run_dynamic(f_myo: float, *, n_filaments: int = 100, n_xl_cand: int = 800, n_myo_cand: int = 400,
                n_ticks: int = 30, tau: float = 0.05, settle_steps: int = 120, seed: int = 0,
                turgor: bool = False, burn_in: int = 10) -> dict:
    """Run the dynamic KMC γ loop; return the time-averaged steady-state γ (after burn-in) + history.

    ``tau`` = KMC tick interval [s]. With α-actinin k_on=10/s, attach prob/tick = 1−exp(−0.5)≈0.39.
    """
    rng = np.random.default_rng(seed)
    dc = build_dynamic_cortex(n_filaments=n_filaments, n_xl_cand=n_xl_cand, n_myo_cand=n_myo_cand,
                              rng=rng)
    _settle(dc, 0.0, n_steps=settle_steps, turgor=False)     # settle resting passive shell first
    hist = []
    for t in range(n_ticks):
        _kmc_tick(dc, f_myo, tau, rng)
        _settle(dc, f_myo, n_steps=settle_steps, turgor=turgor)
        m = measure_gamma_dynamic(dc, f_myo)
        hist.append(m)
    g = np.array([h["gamma_active"] for h in hist[burn_in:]])
    xlb = np.array([h["n_xl_bound"] for h in hist[burn_in:]])
    return {"f_myo": f_myo, "gamma_active_mean": float(g.mean()), "gamma_active_std": float(g.std()),
            "gamma_passive": hist[-1]["gamma_passive"], "n_xl_bound_mean": float(xlb.mean()),
            "n_xl_cand": int(dc.xl_pairs.shape[0]), "n_myo_cand": int(dc.myo_pairs.shape[0]),
            "history": hist}
