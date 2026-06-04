"""L2.4 contact-inhibited proliferation: the mechanistic, size-dependent spreading driver.

The minimal CBM (cohesion + edge-traction, ``cbm.py`` + ``spreading.py``) is COHESION-LOCKED
— A/A0 ≈ 1 at the measured MCF7 scales and the PI law ``A/A0 = a + b/R + c/R²`` does NOT
emerge (a model-limit finding, L2.3 REPORT). The experiment runs over **days**; spreading is
partly **proliferation-driven**. This module adds proliferation as a fine-grained MECHANISM,
not a fitted term:

* Each cell carries a stochastic **cell-cycle timer** (mean = MCF7 doubling time, per-cell CV).
* A cell **divides** when its timer elapses **and** it sits on a **free surface** — i.e. its
  first-coordination shell is not full (fewer than the sphere kissing number Z=12 neighbours).
  This is the canonical Drasdo-Höhme "free-space" division rule: a buried bulk cell (≈12
  contacts) is **contact-inhibited / quiescent**; a rim cell (a free face) proliferates.
* On division the cell **buds** a daughter into its least-crowded direction at the rest
  separation r0 (no repulsive-core blow-up, no gap).

Because the proliferating cells form a **rim of ~constant thickness**, the proliferating
fraction ∝ surface/volume ∝ **1/R**. That geometric surface-to-volume ratio is the EMERGENT
origin of the law's 1/R (and curvature 1/R²) terms — small spheroids grow more (relatively)
than large ones. Nothing is hard-coded to 1/R; it falls out of where division is allowed.

HOOMD has a fixed particle count per ``State``, so growth runs in **epochs**: relax the
current population for ``epoch_steps``, read positions back, apply the divisions that became
due, then rebuild a fresh HOOMD state (``cbm.build_cbm_simulation(positions=…)``) for the
next epoch. Overdamped BAOAB-limit dynamics have no persistent velocity, so rebuilding loses
no physical state. The frozen integrator + single ``md.pair.Morse`` are reused verbatim.

Sanity Gate
-----------
- Dimensional: positions [m]; ages/targets/times [s]; shell cutoff & split [m]; Z, counts [–].
- Boundary (proliferation OFF): with an infinite cycle time no division fires ⇒ this reduces
  exactly to the G1 stable aggregate (count conserved). Verified in tests.
- Boundary (dilute limit): an isolated cell (0 neighbours < Z) divides every cycle ⇒
  N(t) → 2^(t/τ) exponential. Verified in tests.
- Sign-sense: more crowding (higher first-shell count) ⇒ FEWER cells eligible ⇒ slower growth;
  division is biased to the RIM (above-median radial position), never the buried core.
- Conservation: each division is +1 cell (one parent → two daughters); COM budding shift is
  ≤ r0/2 (one daughter inherits the parent site, the other buds at the rest separation).
- Measurement-protocol: A(t)/A0 uses the same ``projected_area`` convex hull as G1/G3.
"""

from __future__ import annotations

import gc
from typing import Any

import gsd.hoomd
import hoomd
import numpy as np
import numpy.typing as npt
from hoomd import md
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from ffn_sim.integrator.baoab_device import make_baoab_updater_for_device
from ffn_sim.spheroid.cbm import (
    _CUTOFF_N_RANGES,
    build_cbm_simulation,
    get_positions,
    make_blob_positions,
    pool_cluster_radius,
)
from ffn_sim.spheroid.observables import (
    core_projected_area,
    projected_area,
    radius_of_gyration,
)

# Single-linkage cluster threshold for the fragmentation-robust core area, in units of r0:
# just past the rest separation (1st shell) so cohesive/touching cells link but a detached
# fragment does not. A measurement-policy choice (between the 1st and 2nd coordination shell).
_CORE_LINK_FACTOR = 1.6
from ffn_sim.spheroid.params import ResolvedL2, ResolvedProliferation

__all__ = [
    "first_shell_counts",
    "candidate_directions",
    "best_bud_direction",
    "sample_cycle_targets",
    "apply_divisions",
    "run_growth",
    "build_pool_simulation",
    "run_growth_pooled",
]

# Number of candidate bud directions probed for free space (a NUMERICAL sampling count, not
# physics). A Fibonacci sphere of 12 points ~ evenly tiles the 4π solid angle at the kissing
# scale — enough to find a free face on a rim cell.
_N_BUD_CANDIDATES = 12


def first_shell_counts(
    positions: npt.NDArray[np.float64], shell_cutoff: float
) -> npt.NDArray[np.int64]:
    """First-coordination-shell neighbour count per cell (excludes self).

    Args:
        positions: (N, 3) cell centers (m).
        shell_cutoff: first-shell radius (m); cells within this distance are "contacts".

    Returns:
        (N,) integer neighbour counts (number of OTHER cells within ``shell_cutoff``).
    """
    if shell_cutoff <= 0.0:
        raise ValueError("shell_cutoff must be strictly positive.")
    p = np.asarray(positions, dtype=np.float64)
    if p.shape[0] < 2:
        return np.zeros(p.shape[0], dtype=np.int64)
    tree = cKDTree(p)
    # query_ball_point with count_only includes self → subtract 1.
    counts = tree.query_ball_point(p, shell_cutoff, return_length=True) - 1
    return counts.astype(np.int64)


def candidate_directions(n: int = _N_BUD_CANDIDATES) -> npt.NDArray[np.float64]:
    """``n`` ~evenly distributed unit vectors on the sphere (deterministic Fibonacci sphere)."""
    if n < 1:
        raise ValueError("n must be >= 1.")
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)              # polar angle, area-uniform
    golden = np.pi * (1.0 + 5.0**0.5)               # golden-angle azimuth
    theta = golden * i
    return np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )


def best_bud_direction(
    positions: npt.NDArray[np.float64],
    i: int,
    neighbor_idx: npt.NDArray[np.int64],
    split: float,
    *,
    candidates: npt.NDArray[np.float64] | None = None,
) -> tuple[npt.NDArray[np.float64], float]:
    """Freest bud direction for cell ``i`` and the free-space gap it achieves.

    Probes a fixed set of candidate directions; for each, the daughter would sit at
    ``positions[i] + split·dir``. Returns the direction MAXIMISING the distance from that
    daughter site to the nearest existing neighbour, together with that distance (the
    "free-space gap"). A buried bulk cell has a small gap (no room); a rim cell a large one.
    With no neighbours the gap is ``+inf`` (anywhere is free) and the first candidate is used.

    Args:
        positions: (N, 3) all cell centers (m).
        i: index of the dividing cell.
        neighbor_idx: indices of ``i``'s nearby cells (the local set to clear; self excluded).
        split: daughter offset from the parent (m), = rest separation r0.
        candidates: optional precomputed candidate unit vectors (defaults to a 12-pt sphere).

    Returns:
        ``(nhat, gap)`` — chosen unit direction and the daughter's nearest-neighbour distance.
    """
    dirs = candidate_directions() if candidates is None else candidates
    if neighbor_idx.size == 0:
        return dirs[0], float("inf")
    cand_sites = positions[i] + split * dirs                 # (K, 3)
    nb = positions[neighbor_idx]                             # (M, 3)
    dmin = cdist(cand_sites, nb).min(axis=1)                 # (K,)
    k = int(np.argmax(dmin))
    return dirs[k], float(dmin[k])


def sample_cycle_targets(
    n: int, prolif: ResolvedProliferation, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Sample ``n`` per-cell cycle-completion times (Gaussian: mean·(1+CV·N(0,1)), floored)."""
    t = prolif.cycle_time_mean * (1.0 + prolif.cycle_time_cv * rng.standard_normal(n))
    # Floor at 10% of the mean so a rare large negative draw cannot give a non-positive target.
    return np.maximum(t, 0.1 * prolif.cycle_time_mean)


def apply_divisions(
    positions: npt.NDArray[np.float64],
    ages: npt.NDArray[np.float64],
    targets: npt.NDArray[np.float64],
    prolif: ResolvedProliferation,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Apply all due, non-inhibited divisions to the population (one pass).

    A timer-elapsed cell (``age >= target``) divides iff it passes the contact-inhibition
    gate: (PRIMARY) the Drasdo-Höhme FREE-SPACE rule — its freest bud direction leaves the
    daughter at least ``prolif.min_gap`` (the Morse repulsive-core onset) from every existing
    cell, i.e. there is room; AND (SECONDARY) its first shell is below the kissing-number
    ceiling. Each dividing cell becomes two: daughter 1 inherits the parent site, daughter 2
    buds at the rest separation into that freest direction. Both daughters reset age 0 and
    draw fresh cycle targets. Quiescent or not-yet-due cells pass through unchanged.

    The free-space gate (not a neighbour count) is what makes division RIM-localised: the
    settled liquid-like packing has bulk free-gap ≈ 0.7·r0 < min_gap (no room) but rim
    free-gap ≈ 1.3·r0 (room) — see the L2.4 REPORT diagnostic.

    Returns a dict with the new ``positions``/``ages``/``targets`` arrays, the parent indices
    that divided (``divided_idx``), their pre-division radial distances from the population
    COM (``divided_radial``), and the population ``median_radial`` (rim-localisation gate).
    """
    p = np.asarray(positions, dtype=np.float64)
    n = p.shape[0]
    counts = first_shell_counts(p, prolif.shell_cutoff)
    com = p.mean(axis=0)
    radial_all = np.linalg.norm(p - com, axis=1)

    due = np.where((ages >= targets) & (counts < prolif.kissing_number))[0]
    if due.size == 0:
        return {
            "positions": p, "ages": ages, "targets": targets,
            "divided_idx": due, "divided_radial": radial_all[due],
            "median_radial": float(np.median(radial_all)) if n else 0.0,
        }

    cands = candidate_directions()
    tree = cKDTree(p)
    # search radius for the free-space test: a daughter sits at `split`; clearing min_gap
    # around it means any neighbour within split + min_gap could matter.
    search_r = prolif.split_distance + prolif.min_gap
    keep = np.ones(n, dtype=bool)
    divided, new_pos, new_ages, new_targets = [], [], [], []
    for i in due:
        nbr = np.array(
            [j for j in tree.query_ball_point(p[i], search_r) if j != i], dtype=np.int64
        )
        nhat, gap = best_bud_direction(p, int(i), nbr, prolif.split_distance, candidates=cands)
        if gap < prolif.min_gap:
            continue  # no room — contact-inhibited (quiescent); timer keeps running
        new_pos.append(p[i])                                   # daughter 1: parent site
        new_pos.append(p[i] + prolif.split_distance * nhat)    # daughter 2: budded
        new_ages.extend([0.0, 0.0])
        new_targets.extend(sample_cycle_targets(2, prolif, rng).tolist())
        keep[i] = False
        divided.append(int(i))

    divided_idx = np.array(divided, dtype=np.int64)
    if divided_idx.size == 0:
        return {
            "positions": p, "ages": ages, "targets": targets,
            "divided_idx": divided_idx, "divided_radial": radial_all[divided_idx],
            "median_radial": float(np.median(radial_all)),
        }

    pos2 = np.vstack([p[keep], np.array(new_pos, dtype=np.float64)])
    ages2 = np.concatenate([ages[keep], np.array(new_ages, dtype=np.float64)])
    targets2 = np.concatenate([targets[keep], np.array(new_targets, dtype=np.float64)])
    return {
        "positions": pos2, "ages": ages2, "targets": targets2,
        "divided_idx": divided_idx, "divided_radial": radial_all[divided_idx],
        "median_radial": float(np.median(radial_all)),
    }


def build_pool_simulation(
    resolved: ResolvedL2,
    n_active_init: int,
    n_max: int,
    *,
    device: hoomd.device.Device,
    seed: int,
) -> tuple[hoomd.Simulation, npt.NDArray[np.bool_], float]:
    """Build ONE leak-free Simulation with a pre-allocated particle pool (L2.4b).

    The State holds ``n_max`` particles for the whole run (N never changes → no Simulation
    rebuild → no HOOMD per-rebuild leak). Two types: ``cell`` (active, interacting via the
    Morse cohesion exactly as ``cbm.py``) and ``void`` (parked, non-interacting — every pair
    involving ``void`` has ``r_cut=0``, and ``void`` is frozen with a large drag so it does not
    drift). Division ACTIVATES a parked ``void`` into a ``cell`` in place via ``set_snapshot``
    (no new Simulation). The first ``n_active_init`` particles are active (a settled-ready
    blob at the origin); the rest are parked far away.

    Returns ``(sim, active_mask, link_unused)`` — ``active_mask`` (n_max,) marks the initial
    active cells; the caller mutates it as voids are activated.
    """
    if not (1 <= n_active_init <= n_max):
        raise ValueError("require 1 <= n_active_init <= n_max.")
    rng = np.random.default_rng(seed)
    r0 = resolved.morse_r0
    r_cut = r0 + _CUTOFF_N_RANGES / resolved.morse_alpha

    # active blob at origin (slightly loose so adhesion settles it, like G1)
    act = make_blob_positions(n_active_init, 1.1 * r0, rng=rng)

    # worst-case active-cluster radius (3D pack vs substrate-wetting 2D disk + spread safety);
    # park voids well beyond it so a spreading spheroid never reaches the box edge (B1 fix).
    r_cluster_max = pool_cluster_radius(r0, n_max)
    n_void = n_max - n_active_init
    park_center = np.array([r_cluster_max + 40.0 * r0, 0.0, 0.0])
    if n_void > 0:
        m = int(np.ceil(n_void ** (1.0 / 3.0)))
        g = (np.arange(m) - (m - 1) / 2.0) * (2.0 * r0)
        xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
        grid = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])[:n_void]
        void = grid + park_center
    else:
        void = np.zeros((0, 3))

    pos = np.vstack([act, void]).astype(np.float64)
    typeid = np.zeros(n_max, dtype=np.uint32)
    typeid[n_active_init:] = 1  # void
    # box must contain the grown cluster AND the parking region with margin (no PBC contact)
    extent = float(np.linalg.norm(pos, axis=1).max())
    L = 2.0 * (extent + r_cluster_max) + 20.0 * r_cut

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_max
    snap.particles.types = ["cell", "void"]
    snap.particles.typeid = typeid
    snap.particles.position = pos
    snap.particles.mass = np.ones(n_max, dtype=np.float64)
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]

    sim = hoomd.Simulation(device=device, seed=seed)
    sim.create_state_from_snapshot(snap)

    nlist = md.nlist.Tree(buffer=resolved.contact_zone_width)
    morse = md.pair.Morse(nlist=nlist, default_r_cut=0.0)
    morse.params[("cell", "cell")] = dict(
        D0=resolved.D_e, alpha=resolved.morse_alpha, r0=resolved.morse_r0
    )
    morse.r_cut[("cell", "cell")] = r_cut
    # void interacts with NOTHING (parked ghost): r_cut 0 for any pair involving void
    for pair in (("cell", "void"), ("void", "void")):
        morse.params[pair] = dict(D0=0.0, alpha=resolved.morse_alpha, r0=resolved.morse_r0)
        morse.r_cut[pair] = 0.0
    morse.mode = "shift"

    ig = md.Integrator(dt=resolved.dt_cfl)
    ig.forces.append(morse)
    sim.operations.integrator = ig

    # void frozen via a large drag (1e6×) so thermal kicks don't move the parked ghosts.
    # Device-aware BAOAB: frozen numpy on CPU, device-resident cupy on GPU (B2 native-N unlock).
    _action, updater = make_baoab_updater_for_device(
        device,
        kT=resolved.kT,
        gamma={"cell": resolved.gamma_cell, "void": resolved.gamma_cell * 1.0e6},
        dt=resolved.dt_cfl,
        seed=seed,
    )
    sim.operations.updaters.append(updater)
    # keep the action alive for the sim's lifetime by stashing it on the sim object
    sim._baoab_action = _action  # noqa: SLF001 — intentional lifetime anchor

    active = np.zeros(n_max, dtype=bool)
    active[:n_active_init] = True
    return sim, active, r_cut


def run_growth_pooled(
    resolved: ResolvedL2,
    prolif: ResolvedProliferation,
    n_cells_init: int = 120,
    *,
    total_time: float,
    epoch_steps: int = 1_200,
    settle_steps: int = 1_000,
    max_cells: int = 4_000,
    device: hoomd.device.Device | None = None,
    seed: int | None = None,
    cohesion: str = "morse",
    cad: "Any" = None,
    substrate: "Any" = None,
    f_traction: float = 0.0,
    Lp: float = 11.0e-6,
) -> dict[str, Any]:
    """Leak-free contact-inhibited growth (L2.4b): ONE Simulation + pre-allocated pool.

    Same physics as ``run_growth`` (contact-inhibited free-space division, MCF7 cycle timer)
    but a parked-``void`` pool replaces the per-epoch Simulation rebuild, so memory stays flat
    (one State of ``max_cells`` particles for the whole run — see the HOOMD-leak note on
    ``run_growth``). Division activates a parked void via ``set_snapshot`` (no rebuild).
    Returns the same record dict as ``run_growth``.

    ``cohesion``: ``"morse"`` (L2.4 static Morse well) or ``"catch"`` (L2.5 force-dependent
    E-cadherin catch bond, ``spheroid.cadherin_bonds``; pass the resolved ``cad``). The catch
    cohesion is a tabulated ``cell``-``cell`` pair, so an activated void daughter feels it
    automatically — no bond bookkeeping. The contact-inhibition / division logic is identical.
    """
    seed = resolved.seed if seed is None else int(seed)
    rng = np.random.default_rng(seed)
    dt = resolved.dt_cfl
    r0 = resolved.morse_r0
    link_r = _CORE_LINK_FACTOR * r0
    device = device or hoomd.device.CPU(notice_level=0)

    if cohesion == "catch":
        if cad is None:
            from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
            cad = resolve_cadherin(resolved)
        from ffn_sim.spheroid.cadherin_bonds import build_cbm_catch
        init_pos = make_blob_positions(n_cells_init, 1.1 * r0, rng=np.random.default_rng(seed))
        sim, _r_cut = build_cbm_catch(
            resolved, cad, n_cells_init, device=device, seed=seed,
            positions=init_pos, n_max=max_cells,
        )
        active = np.asarray(sim.state.get_snapshot().particles.typeid) == 0
    elif cohesion == "morse":
        sim, active, _r_cut = build_pool_simulation(
            resolved, n_cells_init, max_cells, device=device, seed=seed
        )
    else:
        raise ValueError(f"cohesion must be 'morse' or 'catch'; got {cohesion!r}.")

    if substrate is not None:
        # quasi-2D wetting: lift the active cells above the z=0 wall so they settle ONTO the
        # substrate (cells can only live on the +z side), then add the adhesive Morse wall.
        from ffn_sim.spheroid.substrate import add_substrate_wall
        snap0 = sim.state.get_snapshot()
        pos0 = np.array(snap0.particles.position, dtype=np.float64, copy=True)
        zmin = pos0[active, 2].min()
        pos0[active, 2] += (substrate.r0_sub - zmin) + 2.0 * resolved.contact_zone_width
        snap0.particles.position[:] = pos0
        sim.state.set_snapshot(snap0)
        add_substrate_wall(sim, substrate)

    # optional active edge-directed traction (L2.2 active wetting) — the ligand-modulated
    # spreading driver. Per-cell force needs tag==index, so disable the particle sorter.
    edge_force = None
    if f_traction > 0.0:
        from ffn_sim.spheroid.spreading import SettableForce, edge_outward_forces
        sim.operations.tuners.clear()
        edge_force = SettableForce(max_cells)
        sim.operations.integrator.forces.append(edge_force)
        _edge_outward = edge_outward_forces  # local alias

    sim.run(0)
    snap = sim.state.get_snapshot()
    pos_all = np.array(snap.particles.position, dtype=np.float64, copy=True)
    pos_init = pos_all[active].copy()

    sim.run(settle_steps)

    n_max = max_cells
    ages = np.zeros(n_max, dtype=np.float64)
    targets = np.full(n_max, np.inf, dtype=np.float64)  # voids never divide (inf target)
    n0 = int(active.sum())
    targets[active] = sample_cycle_targets(n0, prolif, rng)
    ages[active] = rng.uniform(0.0, 1.0, n0) * targets[active]  # random cycle phase

    def active_positions() -> npt.NDArray[np.float64]:
        s = sim.state.get_snapshot()
        p = np.array(s.particles.position, dtype=np.float64, copy=True)
        return p[active]

    a0 = projected_area(active_positions())
    a0_core = core_projected_area(active_positions(), link_r)

    t = 0.0
    ts, ns, areas, rgs, areas_core = [0.0], [n0], [a0], [radius_of_gyration(active_positions())], [a0_core]
    rim_frac_mean = []
    capped = False
    ejected = False
    cands = candidate_directions()
    search_r = prolif.split_distance + prolif.min_gap
    # box-containment guard (B1): if an active cell wanders past this radius it is about to
    # wrap across the periodic boundary (PBC) and corrupt the run — stop cleanly and flag it
    # rather than letting HOOMD throw a C++ "particle out of bounds" mid-epoch. 0.45·L leaves
    # a 5% margin before the L/2 wrap. A legitimately wetting/growing cluster stays well inside
    # (the box is sized for the worst-case spread via pool_cluster_radius); exceeding this is a
    # runaway (e.g. traction ≳ cohesion → genuine detachment, beyond the model's valid regime).
    eject_radius = 0.45 * float(sim.state.box.Lx)

    while t < total_time:
        if edge_force is not None:
            # recompute edge-directed outward traction for active cells (cheap, per epoch)
            ap = active_positions()
            ef_active = _edge_outward(ap, f_traction=f_traction, Lp=Lp)
            full = np.zeros((max_cells, 3), dtype=np.float64)
            full[np.where(active)[0]] = ef_active
            edge_force.set_vectors(full)
        try:
            sim.run(epoch_steps)
        except RuntimeError as exc:
            # HOOMD "particle out of bounds" from a within-epoch force runaway (traction ≳
            # cohesion). Stop gracefully on the last good state instead of crashing the run.
            ejected = True
            print(f"[run_growth_pooled] ejection at t={t/3600:.1f} h (HOOMD: {exc}); "
                  f"stopping cleanly (ejected=True).")
            break
        dt_epoch = epoch_steps * dt
        t += dt_epoch
        ages[active] += dt_epoch

        snap = sim.state.get_snapshot()
        pos_all = np.array(snap.particles.position, dtype=np.float64, copy=True)
        # post-epoch containment / finiteness guard (catches a slow drift before it wraps)
        act_pos_now = pos_all[np.where(active)[0]]
        if (not np.all(np.isfinite(act_pos_now))) or (
            np.abs(act_pos_now).max() > eject_radius
        ):
            ejected = True
            print(f"[run_growth_pooled] ejection at t={t/3600:.1f} h "
                  f"(cell beyond {eject_radius*1e6:.0f} µm box-guard); stopping cleanly.")
            break
        act_idx = np.where(active)[0]
        pos_act = pos_all[act_idx]
        com = pos_act.mean(axis=0)
        radial = np.linalg.norm(pos_act - com, axis=1)
        median_radial = float(np.median(radial))
        counts = first_shell_counts(pos_act, prolif.shell_cutoff)
        tree = cKDTree(pos_act)

        free_voids = list(np.where(~active)[0])
        due_local = np.where(
            (ages[act_idx] >= targets[act_idx]) & (counts < prolif.kissing_number)
        )[0]
        divided_radial = []
        new_typeid = np.array(snap.particles.typeid, dtype=np.uint32, copy=True)
        changed = False
        for li in due_local:
            if not free_voids:
                capped = True
                break
            nbr = np.array(
                [j for j in tree.query_ball_point(pos_act[li], search_r) if j != li],
                dtype=np.int64,
            )
            nhat, gap = best_bud_direction(pos_act, int(li), nbr, prolif.split_distance, candidates=cands)
            if gap < prolif.min_gap:
                continue  # contact-inhibited (no room)
            v = free_voids.pop()                     # activate a parked void as the daughter
            pos_all[v] = pos_act[li] + prolif.split_distance * nhat
            new_typeid[v] = 0                        # void -> cell
            active[v] = True
            ages[v] = 0.0
            targets[v] = float(sample_cycle_targets(1, prolif, rng)[0])
            gi = act_idx[li]
            ages[gi] = 0.0
            targets[gi] = float(sample_cycle_targets(1, prolif, rng)[0])
            divided_radial.append(radial[li])
            changed = True

        if changed:
            snap.particles.position[:] = pos_all
            snap.particles.typeid[:] = new_typeid
            sim.state.set_snapshot(snap)
        if divided_radial:
            rim_frac_mean.append(float(np.mean(np.array(divided_radial) >= median_radial)))

        ap = active_positions()
        ts.append(t); ns.append(int(active.sum()))
        areas.append(projected_area(ap)); rgs.append(radius_of_gyration(ap))
        areas_core.append(core_projected_area(ap, link_r))
        if capped:
            break

    areas = np.asarray(areas); areas_core = np.asarray(areas_core)
    try:
        pos_final = active_positions()
        if not np.all(np.isfinite(pos_final)):
            raise ValueError("non-finite final positions")
    except (RuntimeError, ValueError):
        pos_final = pos_init  # ejected run: the final snapshot is corrupt; report the seed
    del sim
    gc.collect()
    return {
        "n_cells_init": n_cells_init,
        "a0": a0, "a0_core": a0_core,
        "t": np.asarray(ts), "n_cells": np.asarray(ns),
        "area": areas, "area_over_a0": areas / a0 if a0 > 0 else areas,
        "area_core": areas_core,
        "area_core_over_a0": areas_core / a0_core if a0_core > 0 else areas_core,
        "rg": np.asarray(rgs),
        "rim_fraction_mean": float(np.mean(rim_frac_mean)) if rim_frac_mean else float("nan"),
        "n_division_epochs": len(rim_frac_mean),
        "growth_factor": ns[-1] / ns[0] if ns[0] else float("nan"),
        "capped_at_max_cells": capped,
        "ejected": ejected,
        "f_traction": 0.0,
        "pos_init": pos_init, "pos_final": pos_final,
    }


def run_growth(
    resolved: ResolvedL2,
    prolif: ResolvedProliferation,
    n_cells_init: int = 120,
    *,
    total_time: float,
    epoch_steps: int = 2_000,
    settle_steps: int = 1_500,
    max_cells: int = 4_000,
    max_epochs: int = 80,
    device: hoomd.device.Device | None = None,
    seed: int | None = None,
    f_traction: float = 0.0,
    Lp: float = 11.0e-6,
) -> dict[str, Any]:
    """Grow a CBM spheroid by contact-inhibited proliferation; record A(t), N(t).

    Settles a loose blob (G1), then advances the population in epochs: relax → read back →
    divide → rebuild. Optionally superposes edge-directed active-wetting traction
    (``f_traction`` > 0, recomputed per epoch) so growth and motility act together.

    ⚠️ MEMORY (known limitation): HOOMD's ``State`` has a fixed particle count, so growth
    rebuilds the ``Simulation`` each epoch. HOOMD leaks ~tens of MB per ``Simulation`` rebuild
    on the C++ side (NOT freeable from Python via ``del``/``gc``), so a long run (many epochs)
    can balloon to multi-GB. ``max_epochs`` is a HARD GUARD that stops the run before it can
    run away (it raises if the biological ``total_time`` would need more than ``max_epochs``
    rebuilds at the given ``epoch_steps``). The leak-free fix is a single-Simulation
    pre-allocated particle pool (activate parked particles on division via ``set_snapshot``) —
    the planned L2.4b refactor; until then keep runs SMALL and never launch many concurrently.

    Args:
        resolved: resolved CBM parameters.
        prolif: resolved proliferation parameters.
        n_cells_init: starting cell count.
        total_time: biological time to simulate (s).
        epoch_steps: BAOAB steps between division checks (≈ epoch_steps·dt seconds/epoch).
        settle_steps: pre-growth settle (defines A0).
        max_cells: hard cap (stops the run; logged, not silently truncated).
        device: HOOMD device (CPU default).
        seed: realization key.
        f_traction: optional edge-traction magnitude (N); 0 = proliferation only.
        Lp: edge-traction screening length (m) when ``f_traction`` > 0.

    Returns:
        Dict with time/count/area/Rg series, A0, A/A0 series, rim-localisation diagnostics,
        and init/settled/final positions.
    """
    seed = resolved.seed if seed is None else int(seed)
    rng = np.random.default_rng(seed)
    dt = resolved.dt_cfl

    # HARD memory guard: cap the number of epoch rebuilds (HOOMD leaks per rebuild — see
    # docstring). Refuse up front rather than ballooning to multi-GB mid-run.
    needed_epochs = int(np.ceil(total_time / (epoch_steps * dt)))
    if needed_epochs > max_epochs:
        raise ValueError(
            f"run_growth would need {needed_epochs} epoch rebuilds "
            f"(total_time/{epoch_steps}·dt) but max_epochs={max_epochs}. HOOMD leaks per "
            f"rebuild — raise epoch_steps, lower total_time, or wait for the single-Simulation "
            f"L2.4b refactor. (Guard against the multi-GB runaway.)"
        )

    # Reuse ONE device across all epoch rebuilds. HOOMD State has a fixed particle count, so
    # growth must rebuild the Simulation each epoch — but creating a fresh device every epoch
    # (and not releasing the old Simulation) leaks GBs over ~100 epochs. One device + an
    # explicit per-epoch release (below) keeps the footprint at a single sim's worth.
    device = device or hoomd.device.CPU(notice_level=0)

    # --- initial blob + settle (defines A0) ---
    sim, _a, _u, _rc = build_cbm_simulation(resolved, n_cells_init, device=device, seed=seed)
    sim.run(0)
    pos_init = get_positions(sim)
    sim.run(settle_steps)
    pos = get_positions(sim)
    del sim, _a, _u, _rc
    gc.collect()
    link_r = _CORE_LINK_FACTOR * resolved.morse_r0
    a0 = projected_area(pos)
    a0_core = core_projected_area(pos, link_r)

    ages = np.zeros(pos.shape[0], dtype=np.float64)
    targets = sample_cycle_targets(pos.shape[0], prolif, rng)
    # Desynchronise: start cells at a random phase of their cycle (else all divide at once).
    ages = rng.uniform(0.0, 1.0, pos.shape[0]) * targets

    t = 0.0
    epoch = 0
    ts, ns, areas, rgs = [0.0], [pos.shape[0]], [a0], [radius_of_gyration(pos)]
    areas_core = [a0_core]
    rim_frac_mean = []  # fraction of divisions that were above-median-radial (rim) per epoch
    capped = False

    while t < total_time:
        if pos.shape[0] >= max_cells:
            capped = True
            break
        epoch += 1
        # rebuild a fresh state from the current (possibly grown) population
        sim, _a, _u, _rc = build_cbm_simulation(
            resolved, pos.shape[0], device=device, seed=seed + epoch, positions=pos
        )
        if f_traction > 0.0:
            from ffn_sim.spheroid.spreading import SettableForce, edge_outward_forces

            sim.operations.tuners.clear()
            edge = SettableForce(pos.shape[0])
            edge.set_vectors(edge_outward_forces(pos, f_traction=f_traction, Lp=Lp))
            sim.operations.integrator.forces.append(edge)
        sim.run(0)
        sim.run(epoch_steps)
        pos = get_positions(sim)
        # release this epoch's Simulation/forces before building the next (leak fix)
        del sim, _a, _u, _rc
        if f_traction > 0.0:
            del edge
        gc.collect()

        dt_epoch = epoch_steps * dt
        ages = ages + dt_epoch
        t += dt_epoch

        div = apply_divisions(pos, ages, targets, prolif, rng)
        pos, ages, targets = div["positions"], div["ages"], div["targets"]
        if div["divided_idx"].size > 0:
            rim = float(np.mean(div["divided_radial"] >= div["median_radial"]))
            rim_frac_mean.append(rim)

        ts.append(t); ns.append(pos.shape[0])
        areas.append(projected_area(pos)); rgs.append(radius_of_gyration(pos))
        areas_core.append(core_projected_area(pos, link_r))

    areas = np.asarray(areas)
    areas_core = np.asarray(areas_core)
    return {
        "n_cells_init": n_cells_init,
        "a0": a0,
        "a0_core": a0_core,
        "t": np.asarray(ts),
        "n_cells": np.asarray(ns),
        "area": areas,
        "area_over_a0": areas / a0 if a0 > 0 else areas,
        "area_core": areas_core,
        "area_core_over_a0": areas_core / a0_core if a0_core > 0 else areas_core,
        "rg": np.asarray(rgs),
        "rim_fraction_mean": float(np.mean(rim_frac_mean)) if rim_frac_mean else float("nan"),
        "n_division_epochs": len(rim_frac_mean),
        "growth_factor": ns[-1] / ns[0] if ns[0] else float("nan"),
        "capped_at_max_cells": capped,
        "f_traction": f_traction,
        "pos_init": pos_init,
        "pos_final": pos,
    }
