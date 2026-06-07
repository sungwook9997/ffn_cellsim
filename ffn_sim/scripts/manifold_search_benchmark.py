"""Patch-local vs global crosslinker candidate-search benchmark (H.7 Step-0).

The smallest H.7 surface-manifold prototype
(``docs/v2_audit/H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md`` §Smallest
prototype spec; ``H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`` §VG-1/
VG-6): a STANDALONE benchmark that compares crosslinker candidate generation by

  * **GLOBAL** — the current approach: one ``scipy.spatial.cKDTree`` over all
    cortex beads, ``query_pairs`` at the physical reach (exactly what
    ``crosslinkers.py`` / ``myosin.py`` rebuild every batch tick), and

  * **PATCH-LOCAL** — map each bead to its manifold triangle
    (``surface_manifold.nearest_patch``), then search only within the geodesic
    k-ring patch neighbourhood (``surface_manifold.patch_kring``), where the
    k-ring hop count is DERIVED from ``reach / mean_edge_length`` (never tuned).

It proves the broad-phase win **without touching the simulator or any physics**:
no HOOMD, no Updater, no integrator, no γ. The cortex bead cloud is the REAL
bimodal finite-thickness shell (``generate_bimodal_cortex_layout``) at the REAL
``ℓ₀ = L_filament/(bpf−1)``, ``beads_per_filament``, ``max_bind_dist = 60 nm``
from ``configs/phase1_h3.yaml`` — exactly like ``scripts/crosslink_2d_test.py``.

HARD guardrails reported as asserts/checks (per the task + the design gates):

* **BROAD-PHASE NEUTRALITY (VG-6).** At the same physical reach the patch-local
  candidate *set* must equal the global candidate set (``C_mesh ⊇ C_true`` with
  the k-ring covering the reach, then the SAME physical reach test → identity).
  If they differ, the k-ring is too small → REPORTED, never silently accepted.
* **RESOLUTION INDEPENDENCE (VG-1).** Candidate count and connectivity (z,
  giant-component fraction) must not drift with manifold resolution. If they
  drift, it is FLAGGED (mesh resolution controlling a physical observable = the
  manifold doing physics = halt).
* **NO magic numbers.** The k-ring size is DERIVED from
  (reach / mean triangle edge length); the manifold supplies geometry only, the
  reach comes from the resolved config.

The candidate physics is faithful to ``crosslink_2d_test.py`` /
``connected_mesh``: a candidate is a bead pair of DIFFERENT filaments within the
reach; the filament-pair graph gives mean coordination ``z`` and the giant-
component fraction. The two reaches reproduce the fragment-vs-percolate result
(60 nm → fragmented; ``√(A/n)`` → percolated), numerically demonstrating §8 of
the design: the manifold reorganises/accelerates the search, it does NOT change
the physical bead spacing.

Usage::

    python -m ffn_sim.scripts.manifold_search_benchmark
    python -m ffn_sim.scripts.manifold_search_benchmark --n-filaments 1000 --seed 42
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from ffn_sim.cell.manifest import load_manifest  # noqa: E402
from ffn_sim.cortex.cortex import (  # noqa: E402
    resolve_h3_derived,
    generate_bimodal_cortex_layout,
)
from ffn_sim.cortex.surface_manifold import SurfaceManifold  # noqa: E402

_OUT = (
    Path(__file__).resolve().parents[1]
    / "outputs" / "h7" / "figs" / "manifold_search_benchmark.png"
)

# MCF7 cell radius (Wagner 2011), per the §Smallest-prototype spec. The cortex
# config ships R_cell = 10 µm (generic); the prototype is specified at MCF7.
_R_CELL_MCF7 = 7.5e-6

# Canonical resolution-invariance reference grid: icosphere subdivisions {2,3,4}
# → n_tri {320, 1280, 5120} (VG-1 reference grid), plus level 1 (n_tri 80) for
# breadth → ≥3 resolutions. Each entry is an icosphere subdivision level.
_SUBDIVISIONS = (1, 2, 3, 4)


# ---------------------------------------------------------------------------
# Candidate generation — GLOBAL vs PATCH-LOCAL (both return the SAME object:
# a set of canonical (i<j) DIFFERENT-filament bead-index pairs within `reach`).
# ---------------------------------------------------------------------------
def _candidate_pairs_global(
    beads: np.ndarray, fil_id: np.ndarray, reach: float
) -> tuple[set[tuple[int, int]], int]:
    """Current approach: one global cKDTree (built+queried), ``query_pairs``.

    Mirrors ``crosslinkers.py`` / ``myosin.py``: build a global ``cKDTree`` over
    every bead (rebuilt and discarded every batch tick), find all within-``reach``
    pairs, keep only DIFFERENT-filament pairs (the bridge-different-filament
    narrow-phase). Returns the candidate set as canonical ``(i, j)`` bead-index
    pairs with ``i < j``, plus the within-reach pair count (its narrow-phase
    broad-phase output, for the work-reduction metric).

    Args:
        beads: Bead positions, shape ``(n_bead, 3)`` [m].
        fil_id: Per-bead filament index, shape ``(n_bead,)``.
        reach: Physical candidate-search radius [m].

    Returns:
        ``(candidate_set, n_within_reach_pairs)``.
    """
    tree = cKDTree(beads)
    pairs = tree.query_pairs(reach, output_type="ndarray")
    if len(pairs) == 0:
        return set(), 0
    n_within = int(len(pairs))
    diff = fil_id[pairs[:, 0]] != fil_id[pairs[:, 1]]
    pairs = pairs[diff]
    lo = np.minimum(pairs[:, 0], pairs[:, 1])
    hi = np.maximum(pairs[:, 0], pairs[:, 1])
    return set(zip(lo.tolist(), hi.tolist())), n_within


@dataclass(slots=True)
class PatchIndex:
    """Amortizable bead→patch broad-phase index (built once, reused per query).

    Mirrors the production design: ``bead_tri`` (home patch per bead) + the
    per-patch k-ring pool are STATIC geometry, refreshed lazily every K batch
    ticks — NOT rebuilt per candidate query. So the fair per-query comparison is
    against the global cKDTree, which production rebuilds-and-discards every tick;
    here the patch index is the amortized counterpart.

    Attributes:
        home: Per-bead home-patch id, shape ``(n_bead,)``.
        patch_pool: For each patch ``t`` with ≥1 home bead, the int64 array of
            bead indices in ``t``'s k-ring union (the broad-phase candidate pool
            for any bead homed in ``t``). Empty patches are absent.
        k: The DERIVED k-ring hop count this index was built for.
    """

    home: np.ndarray
    patch_pool: dict[int, np.ndarray]
    k: int


def build_patch_index(
    beads: np.ndarray, manifold: SurfaceManifold, k: int
) -> PatchIndex:
    """Build the amortizable bead→patch k-ring index (the once-per-K-batches cost).

    Args:
        beads: Bead positions, shape ``(n_bead, 3)`` [m].
        manifold: The triangulated surface manifold.
        k: Geodesic k-ring hop count (DERIVED via ``kring_for_reach``).

    Returns:
        A :class:`PatchIndex`.
    """
    home = manifold.nearest_patch(beads)
    # Bucket bead indices by home patch (vectorised via argsort). `groups` holds
    # ORIGINAL bead indices; each group's patch id is home[<any member>].
    order = np.argsort(home, kind="stable")
    home_sorted = home[order]
    boundaries = np.flatnonzero(np.diff(home_sorted)) + 1
    groups = np.split(order, boundaries)
    patch_to_beads: dict[int, np.ndarray] = {
        int(home[g[0]]): g for g in groups if g.size
    }
    occupied = list(patch_to_beads.keys())

    patch_pool: dict[int, np.ndarray] = {}
    ring_cache: dict[int, np.ndarray] = {}
    for t in occupied:
        ring = ring_cache.get(t)
        if ring is None:
            ring = manifold.patch_kring(t, k)
            ring_cache[t] = ring
        pools = [patch_to_beads[int(rt)] for rt in ring.tolist()
                 if int(rt) in patch_to_beads]
        patch_pool[t] = (
            np.concatenate(pools) if pools else np.empty((0,), dtype=np.int64)
        )
    return PatchIndex(home=home, patch_pool=patch_pool, k=k)


def _candidate_pairs_patch_local(
    beads: np.ndarray,
    fil_id: np.ndarray,
    reach: float,
    index: PatchIndex,
) -> tuple[set[tuple[int, int]], int]:
    """Patch-local candidate gather using a PREBUILT patch index (the per-query work).

    For each centre bead, its broad-phase candidate pool is the prebuilt
    ``index.patch_pool[home]`` (beads in the home patch's k-ring); the SAME
    physical narrow-phase then applies (within ``reach`` by a direct vectorised
    distance test — no per-patch tree rebuild — then DIFFERENT filament). The
    candidate *set* equals the global set when the k-ring covers the reach
    (broad-phase neutrality / VG-6).

    Args:
        beads: Bead positions, shape ``(n_bead, 3)`` [m].
        fil_id: Per-bead filament index, shape ``(n_bead,)``.
        reach: Physical candidate-search radius [m].
        index: A :class:`PatchIndex` from :func:`build_patch_index`.

    Returns:
        ``(candidate_set, n_broadphase_pairs_examined)`` — the candidate set as
        canonical ``(i, j)`` bead-index pairs, and the number of broad-phase
        bead pairs the narrow-phase distance-tested (the local-work proxy, which
        is index-quality- and resolution-dependent — the real algorithmic signal).
    """
    home = index.home
    candidates: set[tuple[int, int]] = set()
    n_examined = 0
    # For each occupied patch, restrict the search to its prebuilt k-ring pool:
    # one LOCAL cKDTree over only that pool (spatially pruned, not a dense
    # centre×pool matrix), queried by the patch's home beads at the reach.
    for t, pool in index.patch_pool.items():
        home_mask = home[pool] == t        # beads homed in t (each bead once)
        if not home_mask.any() or pool.size < 2:
            continue
        pool_xyz = beads[pool]             # (P, 3) the broad-phase pool
        local_tree = cKDTree(pool_xyz)
        # local positions of the centre beads within the pool array.
        centre_local = np.flatnonzero(home_mask)
        neigh = local_tree.query_ball_point(pool_xyz[centre_local], reach)
        for cl, nb_list in zip(centre_local.tolist(), neigh):
            bi = int(pool[cl])
            for nj in nb_list:
                n_examined += 1                # within-reach pair examined
                if nj == cl:
                    continue
                bj = int(pool[nj])
                if fil_id[bi] == fil_id[bj]:
                    continue
                candidates.add((bi, bj) if bi < bj else (bj, bi))
    return candidates, n_examined


# ---------------------------------------------------------------------------
# Connectivity (faithful to crosslink_2d_test / connected_mesh diagnostics)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class Connectivity:
    """Filament-graph connectivity from a candidate bead-pair set."""

    n_edges: int          # unique different-filament FILAMENT-pair edges
    z: float              # mean coordination = 2·n_edges / n_filaments
    giant: float          # giant-component fraction
    n_components: int


def _connectivity(
    candidate_pairs: set[tuple[int, int]],
    fil_id: np.ndarray,
    n_fil: int,
) -> Connectivity:
    """Build the filament-pair graph and report z / giant-component fraction.

    Args:
        candidate_pairs: DIFFERENT-filament bead-index pairs (``i < j``).
        fil_id: Per-bead filament index.
        n_fil: Number of filaments.

    Returns:
        A :class:`Connectivity`.
    """
    if not candidate_pairs:
        return Connectivity(n_edges=0, z=0.0, giant=1.0 / max(1, n_fil),
                            n_components=n_fil)
    arr = np.fromiter(
        (c for pair in candidate_pairs for c in pair),
        dtype=np.int64, count=2 * len(candidate_pairs),
    ).reshape(-1, 2)
    f0 = fil_id[arr[:, 0]]
    f1 = fil_id[arr[:, 1]]
    edges = np.unique(np.sort(np.stack([f0, f1], axis=1), axis=1), axis=0)
    data = np.ones(len(edges) * 2)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    adj = coo_matrix((data, (rows, cols)), shape=(n_fil, n_fil)).tocsr()
    n_comp, labels = connected_components(adj, directed=False)
    sizes = np.bincount(labels)
    return Connectivity(
        n_edges=int(len(edges)),
        z=float(2.0 * len(edges) / n_fil),
        giant=float(sizes.max() / n_fil),
        n_components=int(n_comp),
    )


# ---------------------------------------------------------------------------
# Benchmark driver
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ReachResult:
    """Global-vs-patch result for one (reach, resolution) point."""

    reach: float
    subdivisions: int
    n_tri: int
    mean_edge: float
    k: int
    n_global: int           # global candidate count
    n_patch: int            # patch-local candidate count
    sets_equal: bool        # broad-phase neutrality
    t_global: float         # global per-query wall-time [s] (median)
    t_patch: float          # patch-local per-query wall-time [s] (median)
    t_index: float          # patch-index build wall-time [s] (amortized/K)
    n_examined_global: int  # broad-phase pairs the global tree examines
    n_examined_patch: int   # broad-phase pairs the patch gather examines
    conn_global: Connectivity
    conn_patch: Connectivity


def _time_median(fn, repeats: int) -> tuple[object, float]:
    """Run ``fn`` ``repeats`` times; return (last result, median wall-time)."""
    ts = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        ts.append(time.perf_counter() - t0)
    return result, float(np.median(ts))


def run_benchmark(
    n_filaments: int,
    seed: int,
    repeats: int,
    subdivisions: tuple[int, ...] = _SUBDIVISIONS,
    R_cell: float = _R_CELL_MCF7,
) -> tuple[list[ReachResult], dict]:
    """Run the full global-vs-patch benchmark over reaches × resolutions.

    Args:
        n_filaments: Number of cortex filaments (×40 mesoscale; 1000 = MCF7).
        seed: RNG seed for the bimodal layout.
        repeats: Wall-time timing repeats (median reported).
        subdivisions: Icosphere subdivision levels to sweep.
        R_cell: Cell radius [m].

    Returns:
        ``(results, meta)`` — the per-(reach, resolution) results and a metadata
        dict (cloud size, reaches, ℓ₀, etc.).
    """
    # --- resolve the REAL cortex config; override R_cell to MCF7 ---
    cfg = load_manifest("phase1_h3.yaml")
    p = resolve_h3_derived(cfg)
    from dataclasses import replace as _replace
    p = _replace(p, R_cell=R_cell)

    ell0 = p.rest_length                      # = L_filament/(bpf−1)
    max_bind = p.xl_max_bind_dist             # 60 nm dynamic reach (faithful)
    area = 4.0 * math.pi * R_cell ** 2
    bridge_reach = math.sqrt(area / n_filaments)   # √(A/n) connected-mesh recipe

    rng = np.random.default_rng(seed)
    layout = generate_bimodal_cortex_layout(
        p, n_filaments=n_filaments, project_to_shell=True,
        arp_branch_fraction=0.7, rng=rng,
    )
    beads = layout.positions_flat
    fil_id = layout.filament_idx
    n_bead = int(beads.shape[0])
    n_fil_present = int(np.unique(fil_id).size)

    reaches = {"60nm": max_bind, "bridge": bridge_reach}

    results: list[ReachResult] = []
    for subdiv in subdivisions:
        manifold = SurfaceManifold.icosphere(subdiv, R_cell)
        for _label, reach in reaches.items():
            k = manifold.kring_for_reach(reach)   # DERIVED, never tuned

            (global_set, n_ex_global), t_global = _time_median(
                lambda r=reach: _candidate_pairs_global(beads, fil_id, r),
                repeats,
            )
            # Patch index is amortized (lazy, every K batch ticks) — timed apart
            # from the per-query gather, which is the fair counterpart of the
            # per-tick global cKDTree build+query.
            index, t_index = _time_median(
                lambda kk=k: build_patch_index(beads, manifold, kk),
                repeats,
            )
            (patch_set, n_ex_patch), t_patch = _time_median(
                lambda r=reach, idx=index: _candidate_pairs_patch_local(
                    beads, fil_id, r, idx
                ),
                repeats,
            )

            results.append(
                ReachResult(
                    reach=reach, subdivisions=subdiv, n_tri=manifold.n_tri,
                    mean_edge=manifold.mean_edge_length, k=k,
                    n_global=len(global_set), n_patch=len(patch_set),
                    sets_equal=(global_set == patch_set),
                    t_global=t_global, t_patch=t_patch, t_index=t_index,
                    n_examined_global=n_ex_global, n_examined_patch=n_ex_patch,
                    conn_global=_connectivity(global_set, fil_id, n_filaments),
                    conn_patch=_connectivity(patch_set, fil_id, n_filaments),
                )
            )

    meta = {
        "n_filaments": n_filaments,
        "n_filaments_present": n_fil_present,
        "n_bead": n_bead,
        "R_cell": R_cell,
        "ell0": ell0,
        "max_bind": max_bind,
        "bridge_reach": bridge_reach,
        "seed": seed,
        "repeats": repeats,
        "area": area,
    }
    return results, meta


# ---------------------------------------------------------------------------
# Density crossover sweep — WHERE does the patch-local win materialise?
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class DensityPoint:
    """Per-query cost at one bead density (fixed reach, fixed fine resolution)."""

    n_filaments: int
    n_bead: int
    t_global: float          # global per-query wall-time [s]
    t_patch: float           # patch per-query wall-time [s]
    n_examined_global: int   # global within-reach pairs (broad-phase output)
    n_examined_patch: int    # patch within-reach pairs examined
    mean_pool: float         # mean k-ring pool size per occupied patch


def run_density_sweep(
    n_fil_list: tuple[int, ...],
    seed: int,
    repeats: int,
    subdivisions: int,
    R_cell: float = _R_CELL_MCF7,
) -> tuple[list[DensityPoint], dict]:
    """Sweep bead density at FIXED reach + FIXED resolution → cost scaling.

    The design's load-bearing claim (companion §3): the patch-local win is
    *marginal at the ×40 mesoscale and real at native/ECM density*. This holds
    the cell radius and manifold resolution fixed and raises ``n_filaments`` (so
    the shell gets denser), measuring how global vs patch per-query cost scale.
    The physical reach is fixed at the 60 nm dynamic reach (the per-tick binding
    reach), so the patch pool per centre stays bounded by the physical
    neighbourhood while the global tree must (re)build over ALL beads each tick.

    Args:
        n_fil_list: Filament counts to sweep (rising density).
        seed: RNG seed.
        repeats: Wall-time timing repeats (median).
        subdivisions: Fixed icosphere subdivision level for the index.
        R_cell: Cell radius [m].

    Returns:
        ``(points, meta)``.
    """
    cfg = load_manifest("phase1_h3.yaml")
    p = resolve_h3_derived(cfg)
    from dataclasses import replace as _replace
    p = _replace(p, R_cell=R_cell)
    reach = p.xl_max_bind_dist                    # 60 nm dynamic reach (fixed)

    manifold = SurfaceManifold.icosphere(subdivisions, R_cell)
    k = manifold.kring_for_reach(reach)           # DERIVED, fixed reach

    points: list[DensityPoint] = []
    for n_fil in n_fil_list:
        rng = np.random.default_rng(seed)
        layout = generate_bimodal_cortex_layout(
            p, n_filaments=n_fil, project_to_shell=True,
            arp_branch_fraction=0.7, rng=rng,
        )
        beads = layout.positions_flat
        fil_id = layout.filament_idx

        (_gset, n_ex_g), t_global = _time_median(
            lambda: _candidate_pairs_global(beads, fil_id, reach), repeats
        )
        index = build_patch_index(beads, manifold, k)
        (_pset, n_ex_p), t_patch = _time_median(
            lambda idx=index: _candidate_pairs_patch_local(
                beads, fil_id, reach, idx
            ),
            repeats,
        )
        pools = [pool.size for pool in index.patch_pool.values()]
        points.append(
            DensityPoint(
                n_filaments=n_fil, n_bead=int(beads.shape[0]),
                t_global=t_global, t_patch=t_patch,
                n_examined_global=n_ex_g, n_examined_patch=n_ex_p,
                mean_pool=float(np.mean(pools)) if pools else 0.0,
            )
        )
    meta = {"reach": reach, "subdivisions": subdivisions, "n_tri": manifold.n_tri,
            "k": k, "R_cell": R_cell}
    return points, meta


# ---------------------------------------------------------------------------
# Reporting + figure
# ---------------------------------------------------------------------------
def _print_report(results: list[ReachResult], meta: dict) -> dict:
    """Print the benchmark report; return a gate-summary dict."""
    print("=" * 78, flush=True)
    print("MANIFOLD SEARCH BENCHMARK — patch-local k-ring vs global cKDTree", flush=True)
    print(
        f"  cortex: {meta['n_filaments']} filaments "
        f"({meta['n_filaments_present']} placed), {meta['n_bead']} beads, "
        f"R_cell {meta['R_cell']*1e6:.1f} µm, ℓ₀ {meta['ell0']*1e9:.0f} nm",
        flush=True,
    )
    print(
        f"  reaches: 60 nm (dynamic α-actinin/filamin) | "
        f"√(A/n) bridge {meta['bridge_reach']*1e9:.0f} nm | "
        f"timing median of {meta['repeats']} repeats",
        flush=True,
    )
    print("-" * 78, flush=True)

    neutrality_ok = all(r.sets_equal for r in results)

    # ---- broad-phase neutrality (VG-6) ----
    print("BROAD-PHASE NEUTRALITY (VG-6): patch-local set == global set?", flush=True)
    for r in results:
        rl = "60nm" if abs(r.reach - meta["max_bind"]) < 1e-15 else "bridge"
        flag = "PASS" if r.sets_equal else "*** FAIL (k-ring too small) ***"
        print(
            f"  s={r.subdivisions} n_tri={r.n_tri:5d} reach={rl:6s} k={r.k} "
            f"| global={r.n_global:7d} patch={r.n_patch:7d}  {flag}",
            flush=True,
        )
    print(f"  => NEUTRALITY GATE: {'PASS' if neutrality_ok else 'FAIL'}", flush=True)
    print("-" * 78, flush=True)

    # ---- search-cost win (HONEST: conditional, not realized in CPU-Python) ----
    print("SEARCH COST: per-query global cKDTree (rebuilt+discarded each tick) vs", flush=True)
    print("  patch-local k-ring gather over a PREBUILT (amortized/K) patch index.", flush=True)
    print("  HONEST FINDING: scipy's C cKDTree is very fast; at the ×40 mesoscale on", flush=True)
    print("  CPU the patch Python loop is SLOWER, and the coarse-vs-reach patches", flush=True)
    print("  OVER-select (examined ratio <1). The design's win (avoid per-tick rebuild;", flush=True)
    print("  reach-sized patches) is GPU-resident / native-density — see density sweep.", flush=True)
    speedups_60 = []
    speedups_br = []
    work_red_60 = []
    work_red_br = []
    for r in results:
        rl = "60nm" if abs(r.reach - meta["max_bind"]) < 1e-15 else "bridge"
        speed = r.t_global / r.t_patch if r.t_patch > 0 else float("inf")
        work = (r.n_examined_global / r.n_examined_patch
                if r.n_examined_patch > 0 else float("inf"))
        if rl == "60nm":
            speedups_60.append(speed)
            work_red_60.append(work)
        else:
            speedups_br.append(speed)
            work_red_br.append(work)
        print(
            f"  s={r.subdivisions} n_tri={r.n_tri:5d} reach={rl:6s} k={r.k} "
            f"| t_query[g/p]={r.t_global*1e3:6.2f}/{r.t_patch*1e3:6.2f} ms "
            f"(wall ×{speed:4.2f})  t_index={r.t_index*1e3:6.2f} ms "
            f"| broad-phase examined[g/p]={r.n_examined_global}/{r.n_examined_patch} "
            f"(ratio {work:4.2f}; <1 = patch over-selects)",
            flush=True,
        )
    print("-" * 78, flush=True)

    # ---- connectivity (fragment vs percolate) ----
    print("CONNECTIVITY (filament graph): fragment (60nm) vs percolate (√A/n)", flush=True)
    for r in results:
        rl = "60nm" if abs(r.reach - meta["max_bind"]) < 1e-15 else "bridge"
        cg, cp = r.conn_global, r.conn_patch
        verdict = "PERCOLATED" if cg.giant > 0.5 else "FRAGMENTED"
        match = "OK" if (abs(cg.z - cp.z) < 1e-9 and abs(cg.giant - cp.giant) < 1e-9) else "*** MISMATCH ***"
        print(
            f"  s={r.subdivisions} n_tri={r.n_tri:5d} reach={rl:6s} "
            f"| GLOBAL z={cg.z:5.2f} giant={cg.giant*100:5.1f}%  "
            f"PATCH z={cp.z:5.2f} giant={cp.giant*100:5.1f}%  "
            f"-> {verdict}  [{match}]",
            flush=True,
        )
    print("-" * 78, flush=True)

    # ---- resolution independence (VG-1) ----
    print("RESOLUTION INDEPENDENCE (VG-1): observables invariant under n_tri?", flush=True)
    res_ok = True
    for rl, reach in (("60nm", meta["max_bind"]), ("bridge", meta["bridge_reach"])):
        rows = [r for r in results if abs(r.reach - reach) < 1e-15]
        ns = {r.n_patch for r in rows}
        zs = {round(r.conn_patch.z, 9) for r in rows}
        gs = {round(r.conn_patch.giant, 9) for r in rows}
        invariant = (len(ns) == 1 and len(zs) == 1 and len(gs) == 1)
        res_ok = res_ok and invariant
        cand = next(iter(ns))
        z_v = rows[0].conn_patch.z
        g_v = rows[0].conn_patch.giant
        flag = "INVARIANT" if invariant else "*** DRIFTS WITH n_tri — HALT ***"
        print(
            f"  reach={rl:6s}: candidate_count={sorted(ns)} z={sorted(zs)} "
            f"giant={sorted(gs)} -> {flag}",
            flush=True,
        )
        if invariant:
            print(
                f"            (plateau: candidates={cand}, z={z_v:.2f}, "
                f"giant={g_v*100:.1f}% across n_tri ∈ "
                f"{sorted({r.n_tri for r in rows})})",
                flush=True,
            )
    print(f"  => RESOLUTION-INVARIANCE GATE (VG-1): {'PASS' if res_ok else 'FAIL'}", flush=True)
    print("=" * 78, flush=True)

    med_speed_60 = float(np.median(speedups_60)) if speedups_60 else float("nan")
    med_speed_br = float(np.median(speedups_br)) if speedups_br else float("nan")
    med_work_60 = float(np.median(work_red_60)) if work_red_60 else float("nan")
    med_work_br = float(np.median(work_red_br)) if work_red_br else float("nan")
    return {
        "neutrality_ok": neutrality_ok,
        "resolution_ok": res_ok,
        "median_speedup_60nm": med_speed_60,
        "median_speedup_bridge": med_speed_br,
        "median_workreduction_60nm": med_work_60,
        "median_workreduction_bridge": med_work_br,
    }


def _make_figure(
    results: list[ReachResult],
    meta: dict,
    gates: dict,
    out: Path,
    density: list[DensityPoint] | None = None,
    density_meta: dict | None = None,
) -> None:
    """Emit the global-vs-patch comparison figure (viz-integrity compliant)."""
    fig, axes = plt.subplots(2, 3, figsize=(21, 11), constrained_layout=True)
    subdivs = sorted({r.subdivisions for r in results})
    n_tri_by_s = {r.subdivisions: r.n_tri for r in results}
    xs = [n_tri_by_s[s] for s in subdivs]

    def _rows(reach_key):
        reach = meta["max_bind"] if reach_key == "60nm" else meta["bridge_reach"]
        return [next(r for r in results if r.subdivisions == s
                     and abs(r.reach - reach) < 1e-15) for s in subdivs]

    rows60 = _rows("60nm")
    rowsbr = _rows("bridge")

    # -- (a) candidate-set agreement: global vs patch counts vs n_tri --
    ax = axes[0, 0]
    ax.plot(xs, [r.n_global for r in rows60], "o-", color="C0", lw=2,
            label="GLOBAL cKDTree, 60 nm")
    ax.plot(xs, [r.n_patch for r in rows60], "x--", color="C1", ms=10, mew=2,
            label="PATCH k-ring, 60 nm")
    ax.plot(xs, [r.n_global for r in rowsbr], "s-", color="C2", lw=2,
            label="GLOBAL cKDTree, √(A/n)")
    ax.plot(xs, [r.n_patch for r in rowsbr], "+--", color="C3", ms=12, mew=2,
            label="PATCH k-ring, √(A/n)")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("manifold resolution  n_tri  (icosphere faces)")
    ax.set_ylabel("different-filament candidate pairs  [count]")
    neutral = "PASS" if gates["neutrality_ok"] else "FAIL"
    ax.set_title(f"(a) Candidate-set agreement — broad-phase neutrality {neutral}\n"
                 "patch markers must lie ON global lines (sets identical)")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.3)

    # -- (b) per-query wall-time global vs patch vs n_tri (honest: CPU/scale) --
    ax = axes[0, 1]
    ax.plot(xs, [r.t_global * 1e3 for r in rows60], "o-", color="C0", lw=2,
            label="GLOBAL build+query, 60 nm")
    ax.plot(xs, [r.t_patch * 1e3 for r in rows60], "x--", color="C1", ms=10, mew=2,
            label="PATCH k-ring query, 60 nm")
    ax.plot(xs, [r.t_global * 1e3 for r in rowsbr], "s-", color="C2", lw=2,
            label="GLOBAL build+query, √(A/n)")
    ax.plot(xs, [r.t_patch * 1e3 for r in rowsbr], "+--", color="C3", ms=12, mew=2,
            label="PATCH k-ring query, √(A/n)")
    ax.plot(xs, [r.t_index * 1e3 for r in rows60], "v:", color="0.5", lw=1.5,
            label="PATCH index build (amortized/K)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("manifold resolution  n_tri  (icosphere faces)")
    ax.set_ylabel("per-query candidate wall-time  [ms]  (log)")
    sp60 = gates["median_speedup_60nm"]
    ax.set_title(
        "(b) Per-query wall-time (CPU, log-y) — global cKDTree is fast in C;\n"
        f"patch Python loop is ×{1/sp60:.0f} SLOWER here (win is GPU-resident/native-density)"
    )
    ax.legend(fontsize=7, loc="best")
    ax.grid(True, which="both", alpha=0.3)

    # -- (c) connectivity z + giant vs n_tri (resolution invariance) --
    ax = axes[1, 0]
    ax.plot(xs, [r.conn_patch.z for r in rows60], "o-", color="C1", lw=2,
            label="z (patch), 60 nm")
    ax.plot(xs, [r.conn_patch.z for r in rowsbr], "s-", color="C3", lw=2,
            label="z (patch), √(A/n)")
    ax.axhspan(3.0, 3.5, color="0.85", alpha=0.6, label="connected-mesh band z∈[3.0,3.5]")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("manifold resolution  n_tri  (icosphere faces)")
    ax.set_ylabel("mean coordination  z  [neighbours/filament]")
    ax2 = ax.twinx()
    ax2.plot(xs, [r.conn_patch.giant * 100 for r in rows60], "o:", color="C4",
             label="giant % (patch), 60 nm")
    ax2.plot(xs, [r.conn_patch.giant * 100 for r in rowsbr], "s:", color="C5",
             label="giant % (patch), √(A/n)")
    ax2.set_ylabel("giant-component fraction  [%]")
    ax2.set_ylim(0, 105)
    res = "PASS" if gates["resolution_ok"] else "FAIL — HALT"
    ax.set_title(f"(c) Connectivity vs resolution — VG-1 invariance {res}\n"
                 "flat lines = observables do NOT track n_tri")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="center right")
    ax.grid(True, which="both", alpha=0.3)

    # -- (d) the §8 honesty claim: reach vs spacing (fragment vs percolate) --
    ax = axes[1, 1]
    finest = max(subdivs)
    r60 = next(r for r in results if r.subdivisions == finest
               and abs(r.reach - meta["max_bind"]) < 1e-15)
    rbr = next(r for r in results if r.subdivisions == finest
               and abs(r.reach - meta["bridge_reach"]) < 1e-15)
    labels = ["60 nm\n(dynamic reach)", f"√(A/n)={meta['bridge_reach']*1e9:.0f} nm\n(bridge reach)"]
    z_vals = [r60.conn_patch.z, rbr.conn_patch.z]
    g_vals = [r60.conn_patch.giant * 100, rbr.conn_patch.giant * 100]
    x = np.arange(2)
    w = 0.35
    b1 = ax.bar(x - w / 2, z_vals, w, color="C1", label="z (coordination)")
    ax.axhspan(3.0, 3.5, color="0.85", alpha=0.6)
    ax.set_ylabel("mean coordination  z")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    axb = ax.twinx()
    b2 = axb.bar(x + w / 2, g_vals, w, color="C4", label="giant %")
    axb.set_ylabel("giant-component fraction  [%]")
    axb.set_ylim(0, 105)
    for xi, zi, gi in zip(x, z_vals, g_vals):
        verdict = "PERCOLATED" if gi > 50 else "FRAGMENTED"
        ax.text(xi, max(z_vals) * 1.02 + 0.1, verdict, ha="center",
                fontsize=9, fontweight="bold")
    ax.set_title("(d) §8 honesty: reach—not the index—sets spacing\n"
                 "manifold returns the SAME sets the cKDTree does (does NOT fix √A/n)")
    ax.legend([b1, b2], ["z (coordination)", "giant %"], fontsize=8, loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)

    # -- (e) density crossover: per-query wall-time vs bead count (fixed reach) --
    ax = axes[0, 2]
    if density:
        nb = [d.n_bead for d in density]
        ax.plot(nb, [d.t_global * 1e3 for d in density], "o-", color="C0", lw=2,
                label="GLOBAL build+query")
        ax.plot(nb, [d.t_patch * 1e3 for d in density], "s-", color="C1", lw=2,
                label="PATCH k-ring query")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("cortex bead count  N  (rising shell density, fixed R_cell)")
        ax.set_ylabel("per-query wall-time  [ms]  (log)")
        ax.set_title(
            "(e) Density scaling at fixed reach (60 nm) + fixed mesh\n"
            f"(n_tri={density_meta['n_tri']}, k={density_meta['k']}): both ~linear; "
            "patch pool ∝ N at fixed mesh"
        )
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.axis("off")

    # -- (f) density: broad-phase pool size + within-reach work vs N --
    ax = axes[1, 2]
    if density:
        nb = [d.n_bead for d in density]
        ax.plot(nb, [d.mean_pool for d in density], "o-", color="C1", lw=2,
                label="mean k-ring pool size / patch (patch broad-phase)")
        ax.plot(nb, [d.n_examined_global for d in density], "s-", color="C0", lw=2,
                label="global within-reach pairs (its broad-phase output)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("cortex bead count  N  (rising shell density)")
        ax.set_ylabel("broad-phase selectivity  [count]  (log)")
        ax.set_title(
            "(f) Why the win needs a reach-sized patch: at FIXED mesh the\n"
            "patch pool grows ∝ N (the patch must shrink with density for O(1))"
        )
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.axis("off")

    fig.suptitle(
        f"H.7 manifold search benchmark — patch-local k-ring vs global cKDTree  |  "
        f"{meta['n_filaments']} filaments, {meta['n_bead']} beads @ ×40 mesoscale, "
        f"R_cell {meta['R_cell']*1e6:.1f} µm (MCF7)\n"
        f"GEOMETRY-only manifold (no physics): broad-phase NEUTRALITY "
        f"{'PASS' if gates['neutrality_ok'] else 'FAIL'}  |  "
        f"resolution-invariance VG-1 {'PASS' if gates['resolution_ok'] else 'FAIL'}  |  "
        f"k-ring DERIVED from reach/mean_edge (no magic number)  |  "
        f"search-cost win is GPU/native-density-conditional, NOT realized in CPU-Python here "
        f"(matches design §3/§13)",
        fontweight="bold", fontsize=10,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=1000,
                    help="cortex filament count (×40 mesoscale; 1000 = MCF7)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--repeats", type=int, default=5,
                    help="wall-time timing repeats (median reported)")
    ap.add_argument("--no-density-sweep", action="store_true",
                    help="skip the density crossover sweep (faster)")
    ap.add_argument("--out", default=str(_OUT))
    args = ap.parse_args()

    results, meta = run_benchmark(args.n_filaments, args.seed, args.repeats)
    gates = _print_report(results, meta)

    density = density_meta = None
    if not args.no_density_sweep:
        print("DENSITY CROSSOVER SWEEP (fixed 60 nm reach, fixed finest mesh):", flush=True)
        density, density_meta = run_density_sweep(
            (250, 500, 1000, 2000, 4000), seed=args.seed,
            repeats=max(2, args.repeats // 2), subdivisions=max(_SUBDIVISIONS),
        )
        for d in density:
            sp = d.t_global / d.t_patch if d.t_patch > 0 else float("inf")
            print(
                f"  N={d.n_bead:6d} | t_query[g/p]={d.t_global*1e3:6.2f}/"
                f"{d.t_patch*1e3:7.2f} ms (×{sp:4.2f}) | global within-reach="
                f"{d.n_examined_global:6d}  mean k-ring pool={d.mean_pool:6.1f}",
                flush=True,
            )
        print("  => at FIXED mesh the patch pool grows ∝ N (mean pool rises); the "
              "win needs reach-sized patches / GPU residency, per design §3.", flush=True)
        print("=" * 78, flush=True)

    _make_figure(results, meta, gates, Path(args.out), density, density_meta)
    print(f"  fig: {args.out}", flush=True)

    # Hard guardrails as a final assert block (report, do not silently accept).
    if not gates["neutrality_ok"]:
        print("GUARDRAIL VIOLATION: broad-phase neutrality FAILED — patch and "
              "global candidate sets differ (k-ring too small for the reach).",
              flush=True)
    if not gates["resolution_ok"]:
        print("GUARDRAIL VIOLATION: resolution-independence FAILED — a physical "
              "observable tracks manifold resolution (VG-1 halt).", flush=True)
    return 0 if (gates["neutrality_ok"] and gates["resolution_ok"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
