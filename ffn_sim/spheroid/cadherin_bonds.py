"""L2.5 force-dependent E-cadherin catch-bond cohesion for the Layer-2 CBM (replaces Morse).

The L2.4 static Morse well cannot hold a *growing* spheroid together (proliferation tension
fragments it). This module makes cohesion FORCE-DEPENDENT via the faithful Rakshit-2012
sliding-rebinding catch bond (``validation.cadherin_sliding_rebinding``): the per-contact
cohesive force STRENGTHENS under tension up to f0≈29 pN — exactly the regime proliferation
drives — so cohesion self-resists fragmentation, the mechanistic fix the static well lacked.

Why an effective force law (adiabatic elimination), not explicit stochastic bonds
--------------------------------------------------------------------------------
Catch-bond lifetimes are ~0.02–0.2 s (k_off ~ 5–60 s⁻¹) while the overdamped CBM relaxation
timestep is dt ≈ 1.7 s and a cell-cycle is hours — a >4-order timescale gap. Resolving
individual stochastic bind/rupture events would need dt ~ 1 ms (~10⁸ steps for hours of
growth): infeasible. But that same separation (bond kinetics ≪ cell motion, τ_relax≈17 s)
makes **adiabatic elimination** valid and standard: the fast bond ensemble reaches its
steady-state occupancy between cell moves, so we use that steady state as a force-dependent
effective cohesion. The faithful sliding-rebinding k_off(f) enters directly — only the
sub-second stochasticity (which the CBM cannot and need not resolve) is averaged out.

Effective contact cohesion (adiabatic)
--------------------------------------
A cell-cell contact = an ensemble of N_cad cadherins (scale bridge):

    N_cad = F_detach / f0   ≈ 6.5 nN / 29.2 pN ≈ 223     (Iturri 2020 de-adhesion / Rakshit f0)

so the ENSEMBLE catch peak (N_cad·f0) equals the MEASURED MCF7-MCF7 de-adhesion force. The
contact is an elastic spring k_bond (rest length r0); at extension ext = d − r0 the elastic
tension is F_el = k_bond·ext, the per-cadherin force f = F_el/N_cad, and

    k_bond = N_cad · f0 / contact_zone_width            (f reaches f0 at ext = one contact zone)

The contact is bound with steady-state probability (occupancy)

    φ(f) = k_on / (k_on + k_off(f)),   k_on ≡ 1/τ(0) = k_off(0)   (rest-symmetric reformation)

where k_off(f) = 1/τ(f) is the faithful sliding-rebinding effective off-rate. The mean
cohesive (attractive) force the contact transmits is

    F_coh(ext) = φ(f) · F_el = φ(k_bond·ext/N_cad) · k_bond·ext          (ext ≥ 0; 0 if ext<0)

φ is HIGH near f≈f0 (catch: k_off minimal → long-lived → occupied) and COLLAPSES in the slip
regime (f≫f0: k_off large → rupture), so F_coh has a catch hump (peak holding force ~ a large
fraction of F_detach) then decays — the force-strengthening that resists growth fragmentation.
Excluded volume for ext<0 (overlap) is a separate repulsive-only WCA core.

Pool-ready: applied as a ``cell``-``cell`` tabulated pair force; parked ``void`` particles
(typeid 1, the L2.4b growth pool) have r_cut 0 → no interaction.

Sanity Gate
-----------
- Dimensional: ext,d,σ [m]; F_el,F_coh [N]; f [N]; φ [–]; U [J].
- Boundary: ext≤0 → F_coh=0 (WCA repels); φ∈[0,1]; F_coh→0 as ext→∞ (slip rupture).
- Sign/sense: dF_coh/dext>0 in the catch regime (strengthening), <0 past the peak (slip).
- Conservation: a conservative tabulated pair (U,F with F=−dU/dr); BAOAB integrates as before.
- Measurement: peak F_coh is the effective holding force; reported against F_detach (Iturri).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import gsd.hoomd
import hoomd
import numpy as np
import numpy.typing as npt
from hoomd import md

from ffn_sim.integrator.baoab_device import make_baoab_updater_for_device
from ffn_sim.spheroid.cbm import get_positions, make_blob_positions, pool_cluster_radius
from ffn_sim.spheroid.params import ResolvedL2
from ffn_sim.validation.cadherin_sliding_rebinding import (
    RAKSHIT_W2A,
    CadherinCatchParams,
    effective_k_off,
    mean_lifetime,
)

__all__ = [
    "ResolvedCadherin",
    "resolve_cadherin",
    "occupancy",
    "catch_cohesion_force",
    "catch_cohesion_energy",
    "build_cbm_catch",
    "run_g1_catch",
]


@dataclass(frozen=True)
class ResolvedCadherin:
    """Resolved adiabatic catch-bond cohesion parameters (SI)."""

    catch: CadherinCatchParams   # single-cadherin sliding-rebinding kinetics (Rakshit 2012)
    n_cad: float                 # cadherins per contact = F_detach / f0
    k_bond: float                # N/m  contact elastic stiffness = n_cad·f0/contact_zone
    k_on: float                  # s⁻¹  contact reformation rate = k_off(0) (rest-symmetric)
    r0: float                    # m    rest length (cell diameter)
    r_cut: float                 # m    cohesion cutoff (slip-ruptured beyond)
    wca_epsilon: float           # J    repulsive-core energy (excluded volume)
    wca_sigma: float             # m    WCA sigma (= r0/2^(1/6) → repulsion onset at r0)


def resolve_cadherin(
    resolved: ResolvedL2,
    *,
    catch: CadherinCatchParams = RAKSHIT_W2A,
    r_cut_zones: float = 3.0,
) -> ResolvedCadherin:
    """Derive the adiabatic catch-bond cohesion from measured anchors (no tuned constants).

    Args:
        resolved: resolved CBM params (F_detach, r0, contact_zone).
        catch: single-cadherin sliding-rebinding kinetics (Rakshit 2012 Table S1 default).
        r_cut_zones: cohesion cutoff in contact-zone widths beyond r0 (slip rupture tail).
    """
    n_cad = resolved.deadhesion_force_mature / catch.f0
    k_bond = n_cad * catch.f0 / resolved.contact_zone_width
    k_on = 1.0 / mean_lifetime(0.0, catch)                  # = k_off(0), rest-symmetric
    sigma = resolved.morse_r0 / (2.0 ** (1.0 / 6.0))
    wca_eps = resolved.deadhesion_force_mature * resolved.contact_zone_width
    return ResolvedCadherin(
        catch=catch, n_cad=float(n_cad), k_bond=float(k_bond), k_on=float(k_on),
        r0=float(resolved.morse_r0),
        r_cut=float(resolved.morse_r0 + r_cut_zones * resolved.contact_zone_width),
        wca_epsilon=float(wca_eps), wca_sigma=float(sigma),
    )


def occupancy(f: float | npt.NDArray[np.float64], cad: ResolvedCadherin) -> npt.NDArray[np.float64]:
    """Steady-state bound occupancy φ(f) = k_on/(k_on + k_off(f)) of a catch contact."""
    f_arr = np.atleast_1d(np.asarray(f, dtype=np.float64))
    koff = np.array([effective_k_off(float(x), cad.catch) for x in f_arr])
    return cad.k_on / (cad.k_on + koff)


def catch_cohesion_force(
    d: float | npt.NDArray[np.float64], cad: ResolvedCadherin
) -> npt.NDArray[np.float64]:
    """Effective attractive cohesion force magnitude F_coh(d) [N] (0 for d≤r0 or d≥r_cut).

    F_coh = φ(f)·k_bond·(d−r0), f = k_bond·(d−r0)/n_cad. Returns the attractive magnitude
    (the pair force pulls cells together); excluded volume (d<r0) is the separate WCA.
    """
    d_arr = np.atleast_1d(np.asarray(d, dtype=np.float64))
    ext = d_arr - cad.r0
    F_el = cad.k_bond * np.maximum(ext, 0.0)
    f_pc = F_el / cad.n_cad
    phi = occupancy(f_pc, cad)
    F = phi * F_el
    F[(ext <= 0.0) | (d_arr >= cad.r_cut)] = 0.0
    return F


def catch_cohesion_energy(
    r_grid: npt.NDArray[np.float64], cad: ResolvedCadherin
) -> npt.NDArray[np.float64]:
    """Cohesion potential U(r) on a grid s.t. F=−dU/dr (for md.pair.Table; U(r_cut)=0)."""
    F = catch_cohesion_force(r_grid, cad)            # attractive magnitude (force toward smaller r)
    # F_r (radial component) = +F_coh for attraction means force is −F (pulls inward): U' = +F_coh
    # integrate from r_cut inward: U(r) = ∫_r^{r_cut} F_coh dr' (so U(r_cut)=0, U<... attractive well)
    dr = np.gradient(r_grid)
    # cumulative integral from the outer cutoff inward
    U = np.zeros_like(r_grid)
    # trapezoid cumulative from the right
    for i in range(len(r_grid) - 2, -1, -1):
        U[i] = U[i + 1] - 0.5 * (F[i] + F[i + 1]) * (r_grid[i + 1] - r_grid[i])
    return U


def _table_pair(cad: ResolvedCadherin, nlist: "md.nlist.NeighborList", n_grid: int = 200):
    """Build the cell-cell tabulated pair: WCA repulsion (r<r0) + catch cohesion (r0..r_cut)."""
    r_min = 0.5 * cad.r0
    r = np.linspace(r_min, cad.r_cut, n_grid)
    # WCA repulsion (purely repulsive LJ truncated at the minimum r0)
    sig6 = (cad.wca_sigma / r) ** 6
    U_wca = 4.0 * cad.wca_epsilon * (sig6 ** 2 - sig6) + cad.wca_epsilon
    F_wca = 24.0 * cad.wca_epsilon * (2.0 * sig6 ** 2 - sig6) / r
    rep = r < cad.r0
    U_wca = np.where(rep, U_wca, 0.0)
    F_wca = np.where(rep, F_wca, 0.0)
    # catch cohesion (attraction) for r>=r0
    F_coh = catch_cohesion_force(r, cad)             # attractive magnitude
    U_coh = catch_cohesion_energy(r, cad)
    U = U_wca + U_coh
    F = F_wca - F_coh                                 # radial force: + repulsive, − attractive
    table = md.pair.Table(nlist=nlist, default_r_cut=0.0)
    table.params[("cell", "cell")] = dict(r_min=r_min, U=U, F=F)
    table.r_cut[("cell", "cell")] = cad.r_cut
    for pair in (("cell", "void"), ("void", "void")):
        table.params[pair] = dict(r_min=r_min, U=np.zeros(n_grid), F=np.zeros(n_grid))
        table.r_cut[pair] = 0.0
    return table


def build_cbm_catch(
    resolved: ResolvedL2,
    cad: ResolvedCadherin,
    n_cells: int,
    *,
    device: hoomd.device.Device | None = None,
    init_spacing_factor: float = 1.05,
    seed: int | None = None,
    positions: npt.NDArray[np.float64] | None = None,
    n_max: int | None = None,
) -> tuple[hoomd.Simulation, float]:
    """Build a CBM sim with the adiabatic catch-bond cohesion (tabulated pair) + BAOAB.

    If ``positions``/``n_max`` are given, builds a pre-allocated pool (cells = given positions,
    rest parked as ``void``) for leak-free growth (L2.4b style). Returns ``(sim, r_cut)``.
    """
    seed = resolved.seed if seed is None else int(seed)
    rng = np.random.default_rng(seed)
    r0 = cad.r0

    if positions is not None:
        act = np.asarray(positions, dtype=np.float64)
        act = act - act.mean(axis=0)
        n_active = act.shape[0]
        N = n_active if n_max is None else int(n_max)
        if N < n_active:
            raise ValueError("n_max must be >= number of active positions.")
        # worst-case active-cluster radius (3D pack vs substrate-wetting 2D disk + spread
        # safety); park voids beyond it so a spreading spheroid never leaves the box (B1 fix).
        r_cluster_max = pool_cluster_radius(r0, N)
        n_void = N - n_active
        if n_void > 0:
            m = int(np.ceil(n_void ** (1.0 / 3.0)))
            g = (np.arange(m) - (m - 1) / 2.0) * (2.0 * r0)
            xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
            grid = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])[:n_void]
            # Numerical-policy (void containment), NOT physics: park the void block's centre
            # one full worst-case cluster radius BEYOND the worst-case active edge — i.e. at
            # 2·r_cluster_max from the active COM (was a fixed +40·r0, tiny vs a ~400 µm
            # cluster at native N). The void block half-extent is (m−1)·r0 ≈ n_void^(1/3)·r0,
            # and r_cluster_max ≥ 2·r0·N^(1/3) ≥ 2·r0·n_void^(1/3) > that half-extent, so the
            # nearest void stays ≥ (r_cluster_max − block_halfwidth) > 0 beyond the worst-case
            # active edge — always ≫ r_cut (r_cut < 2·r0, r_cluster_max grows with N). Voids
            # have ×1e6 drag (frozen below) so this offset only guarantees they never enter
            # r_cut of an active cell; it changes no force or measurement.
            void_offset = 2.0 * r_cluster_max
            void = grid + np.array([void_offset, 0.0, 0.0])
            pos = np.vstack([act, void])
        else:
            pos = act
        typeid = np.zeros(N, dtype=np.uint32)
        typeid[n_active:] = 1
    else:
        pos = make_blob_positions(n_cells, r0 * init_spacing_factor, rng=rng)
        N = n_cells
        typeid = np.zeros(N, dtype=np.uint32)

    # Box L (HOOMD nlist/PBC only — never a force or a measurement): 2·extent covers every
    # particle (active cells AND the parked void block, since extent is the max norm over all
    # of pos), and the margin is the per-axis minimum-image gap. PBC safety needs that gap
    # ≥ r_cut so no particle ever sees its own image inside the cohesion cutoff. We keep the
    # historical 20·r0 floor (unchanged small-N behaviour) and additionally guarantee the gap
    # spans r_cut even if r_cut_zones is raised — numerical-policy (containment), not physics.
    extent = float(np.linalg.norm(pos, axis=1).max())
    pbc_margin = max(20.0 * r0, 2.0 * cad.r_cut)
    L = 2.0 * extent + pbc_margin

    snap = gsd.hoomd.Frame()
    snap.particles.N = N
    snap.particles.types = ["cell", "void"]
    snap.particles.typeid = typeid
    snap.particles.position = pos.astype(np.float64)
    snap.particles.mass = np.ones(N, dtype=np.float64)
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]

    sim = hoomd.Simulation(device=device or hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(snap)

    nlist = md.nlist.Tree(buffer=resolved.contact_zone_width)
    table = _table_pair(cad, nlist)
    ig = md.Integrator(dt=resolved.dt_cfl)
    ig.forces.append(table)
    sim.operations.integrator = ig

    # Device-aware BAOAB: frozen numpy on CPU (byte-stable), device-resident cupy on GPU
    # (the B2 native-N unlock — removes the per-step host sync). See integrator/baoab_device.py.
    _action, baoab = make_baoab_updater_for_device(
        device or sim.device,
        kT=resolved.kT, gamma={"cell": resolved.gamma_cell, "void": resolved.gamma_cell * 1e6},
        dt=resolved.dt_cfl, seed=seed,
    )
    sim.operations.updaters.append(baoab)
    sim._baoab_action = _action  # noqa: SLF001 lifetime anchor
    return sim, cad.r_cut


def run_g1_catch(
    resolved: ResolvedL2,
    cad: ResolvedCadherin,
    n_cells: int = 200,
    *,
    settle_steps: int = 20_000,
    measure_steps: int = 10_000,
    device: hoomd.device.Device | None = None,
    seed: int | None = None,
) -> dict[str, Any]:
    """G1-with-catch-cohesion: a loose blob must settle into a stable cohesive aggregate."""
    from ffn_sim.spheroid.observables import (
        detached_fraction,
        nearest_neighbor_stats,
        radius_of_gyration,
    )

    sim, _rc = build_cbm_catch(resolved, cad, n_cells, device=device, seed=seed)
    sim.run(0)
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos_settled = get_positions(sim)
    sim.run(measure_steps)
    pos_final = get_positions(sim)
    nn = nearest_neighbor_stats(pos_final)
    det = detached_fraction(pos_final, d_crit=3.0 * resolved.morse_r0, neighbor_radius=1.5 * resolved.morse_r0)
    return {
        "n_cells": n_cells,
        "nn_median_over_r0": nn["median"] / resolved.morse_r0,
        "detached_fraction": det,
        "rg_settled": radius_of_gyration(pos_settled),
        "rg_final": radius_of_gyration(pos_final),
        "rg_growth_factor": radius_of_gyration(pos_final) / radius_of_gyration(pos_settled),
        "pos_init": pos_init, "pos_settled": pos_settled, "pos_final": pos_final,
    }
