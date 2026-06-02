"""Center-based multicellular spheroid builder (L2.1): N cells = N Morse particles on HOOMD.

Each cell is ONE particle of type ``cell``. A single ``md.pair.Morse`` provides cell-cell
adhesion (attractive well, depth D_e at the rest separation r0 = cell diameter) AND
excluded volume (repulsive core for r < r0). Dynamics is the frozen overdamped
Leimkuhler-Matthews BAOAB integrator (``ffn_sim.integrator.baoab``), reused verbatim.

This is the MVP physics for the Layer-2 line. The STATIC Morse well is the documented MVP
scaffold for the mechanistic KU-4.2 cadherin catch-bond upgrade (phase L2.5); D_e is
derived from that catch-bond chain (see ``params.py`` / ``docs/LAYER2_ANCHORS_2026-06-02``).

Sanity Gate
-----------
- Dimensional: position [m]; Morse D0 [J], alpha [1/m], r0 [m], r_cut [m]; gamma [N·s/m];
  dt [s]; kT [J]. All consumed from a resolved ``ResolvedL2`` (SI, validated).
- CFL: dt = safety·gamma/k_spring (resolve_layer2); overdamped, sub-relaxation-time.
- Boundary (no explosive overlap): cells are seeded at >= r0 separation (init_spacing_factor
  >= 1), so they start in the ATTRACTIVE branch — never inside the stiff repulsive core that
  would blow up the timestep.
- Boundary (no PBC self-interaction): the cubic box L is set to 2·(aggregate extent) plus
  several r_cut of margin, so no particle sees its periodic image within r_cut.
- Conservation: G1 has no active force => aggregate count conserved; COM drifts only by
  thermal noise (no net propulsion).
- Sign-sense: D0 > 0 => attractive well at r0 pulls a loose blob inward; r < r0 repulsion
  resists collapse => a stable aggregate (the G1 contract).
- Integrator reuse: ``md.Integrator`` carries the Morse force only (methods=[] default);
  the BAOAB Action advances positions (the frozen-integrator contract).

r_cut and the nlist skin buffer are builder-level NUMERICAL-POLICY choices (like the
Mikado builder's ``buffer = 0.5·sigma``): r_cut = r0 + 5/alpha truncates the Morse tail at
e^-5 (~0.7% of well depth); ``mode="shift"`` makes the potential continuous at r_cut.
"""

from __future__ import annotations

from typing import Any

import gsd.hoomd
import hoomd
import numpy as np
import numpy.typing as npt
from hoomd import md

from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.spheroid.observables import (
    detached_fraction,
    nearest_neighbor_stats,
    radius_of_gyration,
)
from ffn_sim.spheroid.params import ResolvedL2

# Numerical-policy constant: truncate the Morse tail this many decay lengths (1/alpha)
# beyond the rest separation. e^-5 ~ 0.7% of the well depth -> negligible truncation.
_CUTOFF_N_RANGES = 5.0


def make_blob_positions(
    n_cells: int,
    spacing: float,
    *,
    rng: np.random.Generator,
    jitter_frac: float = 0.05,
) -> npt.NDArray[np.float64]:
    """Generate a compact, roughly-spherical initial blob of ``n_cells`` cell centers.

    A jittered simple-cubic lattice at the given ``spacing`` is generated and the
    ``n_cells`` points closest to the origin are kept (=> approximately spherical), then
    perturbed by Gaussian jitter. Centered at the origin.

    Args:
        n_cells: number of cells (particles).
        spacing: lattice spacing (metres); use >= r0 to start in the attractive branch.
        rng: numpy Generator (seeded by the caller for reproducibility).
        jitter_frac: Gaussian jitter as a fraction of ``spacing``.

    Returns:
        (n_cells, 3) positions (metres), centered on the centroid.
    """
    if n_cells < 1:
        raise ValueError("n_cells must be >= 1.")
    # SC fill fraction ~ pi/6; oversize the lattice so the inner sphere holds n_cells.
    m = int(np.ceil((n_cells * 6.0 / np.pi) ** (1.0 / 3.0))) + 2
    g = (np.arange(m) - (m - 1) / 2.0) * spacing
    xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
    pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    keep = np.argsort(np.linalg.norm(pts, axis=1))[:n_cells]
    pts = pts[keep]
    pts = pts + rng.normal(0.0, jitter_frac * spacing, pts.shape)
    return pts - pts.mean(axis=0)


def build_cbm_simulation(
    resolved: ResolvedL2,
    n_cells: int,
    *,
    device: hoomd.device.Device | None = None,
    init_spacing_factor: float = 1.1,
    box_margin_ranges: float = 2.0,
    rng: np.random.Generator | None = None,
    seed: int | None = None,
    positions: npt.NDArray[np.float64] | None = None,
) -> tuple[hoomd.Simulation, Any, Any, float]:
    """Construct and wire a HOOMD center-based spheroid simulation.

    Args:
        resolved: fully-resolved Layer-2 parameters (``resolve_layer2``).
        n_cells: number of cells (1 particle each). Ignored when ``positions`` is given
            (the count is taken from ``positions``).
        device: HOOMD device (defaults to CPU).
        init_spacing_factor: initial lattice spacing as a multiple of r0 (>= 1 keeps the
            seed in the attractive branch; default 1.1 = slightly loose, adhesion settles).
        box_margin_ranges: extra box half-margin in units of r_cut beyond the blob extent.
        rng: optional numpy Generator (defaults to one seeded by ``resolved.seed``).
        positions: optional explicit (N, 3) cell centers (metres). When given, the blob
            generator is bypassed and these centers are used verbatim (centered on their
            centroid) — the entry point the proliferation epoch loop uses to re-seed a
            grown population into a fresh HOOMD state. Must already be in the attractive
            branch (centre separations >= r0) to keep the timestep stable.

    Returns:
        ``(sim, baoab_action, baoab_updater, r_cut)``.

    Raises:
        ValueError: on non-positive init_spacing_factor or n_cells < 1.
    """
    if init_spacing_factor < 1.0:
        raise ValueError("init_spacing_factor must be >= 1 (else cells start overlapping).")
    seed = resolved.seed if seed is None else int(seed)  # realization key (ensemble)
    rng = rng if rng is not None else np.random.default_rng(seed)

    r_cut = resolved.morse_r0 + _CUTOFF_N_RANGES / resolved.morse_alpha
    if positions is not None:
        pos = np.asarray(positions, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 3 or pos.shape[0] < 1:
            raise ValueError("positions must be a non-empty (N, 3) array.")
        pos = pos - pos.mean(axis=0)
        n_cells = pos.shape[0]
    else:
        pos = make_blob_positions(
            n_cells, resolved.morse_r0 * init_spacing_factor, rng=rng
        )
    extent = float(np.linalg.norm(pos, axis=1).max())
    L = 2.0 * extent + 2.0 * box_margin_ranges * r_cut + resolved.morse_r0

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_cells
    snap.particles.types = ["cell"]
    snap.particles.typeid = np.zeros(n_cells, dtype=np.uint32)
    snap.particles.position = pos.astype(np.float64)
    snap.particles.mass = np.ones(n_cells, dtype=np.float64)
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]

    sim = hoomd.Simulation(device=device or hoomd.device.CPU(), seed=seed)
    sim.create_state_from_snapshot(snap)

    nlist = md.nlist.Tree(buffer=resolved.contact_zone_width)
    morse = md.pair.Morse(nlist=nlist, default_r_cut=0.0)
    morse.params[("cell", "cell")] = dict(
        D0=resolved.D_e, alpha=resolved.morse_alpha, r0=resolved.morse_r0
    )
    morse.r_cut[("cell", "cell")] = r_cut
    morse.mode = "shift"

    ig = md.Integrator(dt=resolved.dt_cfl)  # methods=[] (default) — L-M Action contract
    ig.forces.append(morse)
    sim.operations.integrator = ig

    action, updater = make_baoab_updater(
        kT=resolved.kT,
        gamma={"cell": resolved.gamma_cell},
        dt=resolved.dt_cfl,
        seed=seed,
    )
    sim.operations.updaters.append(updater)
    return sim, action, updater, r_cut


def get_positions(sim: hoomd.Simulation) -> npt.NDArray[np.float64]:
    """Return a copy of current particle positions (N x 3, metres)."""
    snap = sim.state.get_snapshot()
    return np.array(snap.particles.position, dtype=np.float64, copy=True)


def run_g1(
    resolved: ResolvedL2,
    n_cells: int = 200,
    *,
    settle_steps: int = 20_000,
    measure_steps: int = 10_000,
    device: hoomd.device.Device | None = None,
    detached_d_crit_over_r0: float = 3.0,
    detached_neighbor_over_r0: float = 1.5,
) -> dict[str, Any]:
    """Run the G1 stable-aggregate test and return metrics (no pass/fail — that's the gate).

    Builds a loose blob, settles it, then measures over a window. Returns the metrics the
    G1 acceptance bands (oracle config ``layer2_cbm.yaml``) are checked against. This
    function stays oracle-free (runtime); the band comparison lives in the smoke script.

    Args:
        resolved: resolved Layer-2 parameters.
        n_cells: number of cells.
        settle_steps: steps to equilibrate the blob before measuring.
        measure_steps: steps between the settled and final snapshots (dispersal check).
        device: HOOMD device (CPU default).
        detached_d_crit_over_r0, detached_neighbor_over_r0: detached-fraction thresholds as
            multiples of r0 (from the G1 gate contract).

    Returns:
        Dict of metrics + the three position snapshots (init/settled/final).
    """
    sim, _action, _updater, r_cut = build_cbm_simulation(
        resolved, n_cells, device=device
    )
    sim.run(0)  # seed net_force for the L-M Action before it reads it
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos_settled = get_positions(sim)
    sim.run(measure_steps)
    pos_final = get_positions(sim)

    nn = nearest_neighbor_stats(pos_final)
    det = detached_fraction(
        pos_final,
        d_crit=detached_d_crit_over_r0 * resolved.morse_r0,
        neighbor_radius=detached_neighbor_over_r0 * resolved.morse_r0,
    )
    rg_settled = radius_of_gyration(pos_settled)
    rg_final = radius_of_gyration(pos_final)

    return {
        "n_cells": n_cells,
        "r_cut": r_cut,
        "nn_median": nn["median"],
        "nn_median_over_r0": nn["median"] / resolved.morse_r0,
        "detached_fraction": det,
        "rg_settled": rg_settled,
        "rg_final": rg_final,
        "rg_growth_factor": rg_final / rg_settled if rg_settled > 0 else float("inf"),
        "pos_init": pos_init,
        "pos_settled": pos_settled,
        "pos_final": pos_final,
    }
