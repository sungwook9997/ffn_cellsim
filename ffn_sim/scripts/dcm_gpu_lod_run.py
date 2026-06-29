"""GPU-friendly DCM spheroid + activity-LOD run-script (gbook RTX A5000 / CPU dev).

Builds the GPU-FRIENDLY DCM spheroid (``cell/dcm_gpu_build.py``: custom K1-turgor,
single-bond-type edge springs, tent contact, adhesive substrate, optional active
rim traction — NO native md.mesh, NO per-cell mesh types) at a parametrised N,
runs it on the GPU (auto-selected; CPU fallback on a CUDA-less Mac), and measures:

  * A/A₀ — basal footprint area (convex hull of basal nodes) / initial.
  * effective R — radius of gyration of all nodes (spheroid spread).
  * active-rim / necrotic / inert fractions from the activity-LOD classifier.
  * wall-time — and, with ``--ab``, the wall-time WITHOUT LOD at the same N (the
    LOD speedup number).

Writes a JSON (metrics) + a PNG figure (footprint + LOD fractions + wall-time bar)
into ``ffn_sim/outputs/h_dcm_gpu_lod/``.

THIS is what runs on the gbook A5000 — there it picks ``hoomd.device.GPU`` and the
GPU-friendly stack runs device-resident (no many-mesh-types stall). On this Mac it
runs the SAME code on the CPU fallback to validate the finite gate + LOD speedup.

Run from repo root:
    PYTHONPATH=. python scripts/dcm_gpu_lod_run.py --n-cells 40 --steps 8000 --ab
    # gbook A5000 (GPU auto-selected):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python scripts/dcm_gpu_lod_run.py \
        --n-cells 200 --steps 20000 --active
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd

import hoomd as _hoomd

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import (
    ResolvedGpuDCM,
    build_gpu_dcm_simulation,
    build_gpu_spheroid_prolif,
    GpuProliferationUpdater,
    attach_activity_lod_prolif,
)
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_lod import (
    ResolvedLOD,
    attach_activity_lod,
)
from ffn_sim.common.sim_realtime import map_realtime, division_probability

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod")


def _positions(sim) -> np.ndarray:
    """Tag-ordered (N,3) positions from the active device's snapshot."""
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _footprint_area(pos: np.ndarray, z0: float, band: float) -> float:
    """Convex-hull area of basal nodes (z < z0 + band) — the spread footprint."""
    basal = pos[pos[:, 2] < z0 + band]
    if basal.shape[0] < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(basal[:, :2]).volume)  # 2D hull "volume" = area
    except Exception:
        # fallback: bounding-box area
        xy = basal[:, :2]
        return float((xy[:, 0].ptp()) * (xy[:, 1].ptp()))


def _live_node_mask(h: dict, n: int) -> np.ndarray:
    """Boolean node mask of cells that are LIVE in the spheroid (not dormant).

    A pool cell is live iff its nodes carry a non-negative ``cell_of_node`` — the
    proliferation updater sets that ON for activated daughters and leaves the
    parked reserve at −1. This is the authoritative liveness signal: the LOD
    classifier may flip the (shared) traction ``active`` mask for isolated parked
    cells, but it NEVER touches ``cell_of_node``/``face_cell`` (the turgor + tent
    governing arrays), so counting via ``cell_of_node`` excludes the dormant pool.
    """
    return np.asarray(h["cell_of_node"]) >= 0


def _footprint_active(h: dict, pos: np.ndarray, z0: float, band: float) -> float:
    """Footprint of LIVE pool cells only (excludes parked dormant nodes)."""
    return _footprint_area(pos[_live_node_mask(h, pos.shape[0])], z0, band)


def _active_pos(h: dict, pos: np.ndarray) -> np.ndarray:
    """Positions of LIVE pool-cell nodes only (parked reserve dropped)."""
    return pos[_live_node_mask(h, pos.shape[0])]


def _radius_of_gyration(pos: np.ndarray) -> float:
    c = pos.mean(0)
    return float(np.sqrt(((pos - c) ** 2).sum(1).mean()))


def _live_necrotic_count(h: dict, pos: np.ndarray, necrotic_depth_um: float) -> int:
    """Necrotic count over LIVE pool cells ONLY (parked dormant cells excluded).

    The LOD classifier (``dcm_gpu_lod.py``, read-only) computes its depth /
    centroid over ALL n_max pool cells, so the parked dormant cells (z≈60·R, far
    off) corrupt the surface radius and make every LIVE cell read "deep" →
    spuriously necrotic. Here we restrict to LIVE cells (``cell_of_node ≥ 0``),
    compute their OWN centroid + surface radius, and count cells deeper than the
    necrosis onset. At R_surface ≪ 150 µm this is 0, as physics requires.
    """
    live = _live_node_mask(h, pos.shape[0])
    ranges = h["ranges"]
    con = np.asarray(h["cell_of_node"])
    live_cells = [c for c in range(len(ranges))
                  if con[ranges[c][0]] >= 0]
    if len(live_cells) == 0:
        return 0
    cents = np.array([pos[ranges[c][0]:ranges[c][1]].mean(0) for c in live_cells])
    cluster_cen = cents.mean(0)
    r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
    depth = float(r_cell.max()) - r_cell
    return int((depth > necrotic_depth_um * 1.0e-6).sum())


def _run_one(p: ResolvedGpuDCM, n_cells: int, steps: int, *, active: bool,
             use_lod: bool, lod_cfg: ResolvedLOD, device=None,
             prolif: bool = False, p_div: float = 0.04, div_every: int = 4000,
             n_max: int | None = None, auto_pdiv: bool = False,
             cell_cycle_h: float = 20.0):
    """Build + run one configuration; return (metrics dict, sim).

    With ``prolif=True`` the GPU pool build is used: ``n_cells`` cells start
    active and a dormant reserve up to ``n_max`` is parked for division. A
    ``GpuProliferationUpdater`` fires every ``div_every`` steps (necrosis-gated by
    the LOD ``cell_inert`` mask when LOD is on).

    With ``auto_pdiv=True`` the division probability is DERIVED from physics
    (``common/sim_realtime.division_probability``) instead of the passed ``p_div``:
    p_div = S·dt·div_every / T_cycle, so the realised divisions reproduce the
    cell-cycle fraction the spreading episode occupies (a PREDICTION, not a fit).
    """
    band = lod_cfg.contact_band * p.R_cell

    if prolif:
        n_max = n_max if n_max is not None else int(np.ceil(n_cells * 1.6))
        h = build_gpu_spheroid_prolif(p, n_cells, n_max, device=device,
                                      active=active)
        n_start = n_cells
    else:
        h = build_gpu_dcm_simulation(p, n_cells, device=device, active=active)
        n_start = n_cells
    sim = h["sim"]

    if use_lod:
        # PROLIF builds use the LIVE-cell-only classifier (Problem 2 fix): the
        # parked dormant pool would otherwise mis-mark live cells necrotic/inert
        # and collapse the active rim. Plain builds use the frozen attach.
        if prolif:
            attach_activity_lod_prolif(h, lod_cfg)
        else:
            attach_activity_lod(h, lod_cfg)

    # PHYSICAL division-rate anchoring (Problem 1, PI 2026-06-11). Derive p_div
    # from the real cell-cycle time vs the run's real-time-equivalent duration,
    # NOT a free knob tuned to an A/A₀ band. p_div = S·dt·div_every / T_cycle.
    div_anchor = None
    if prolif and auto_pdiv:
        # rim estimate for the EXPECTED-TOTAL print only (does not enter p_div).
        # A small spheroid is nearly all-rim → n_start is a fair rim estimate.
        div_anchor = division_probability(
            steps, p.dt, t_cycle_h=cell_cycle_h, div_every=div_every,
            n_rim_est=n_start)
        p_div = div_anchor.p_div

    prolif_upd = None
    if prolif:
        prolif_upd = GpuProliferationUpdater(
            handles=h, p_div=p_div, gap_factor=0.4,
            cell_inert=h.get("cell_inert"),
            necrotic_depth_um=lod_cfg.necrotic_depth_um, rng_seed=p.seed + 3)
        sim.operations.updaters.append(_hoomd.update.CustomUpdater(
            action=prolif_upd, trigger=_hoomd.trigger.Periodic(div_every)))

    # finite gate FIRST
    sim.run(0)
    pos0 = _positions(sim)
    # footprint / Rg of the ACTIVE subset only (parked dormant pool cells excluded).
    A0 = _footprint_active(h, pos0, p.z_substrate, band) if prolif \
        else _footprint_area(pos0, p.z_substrate, band)
    Rg0 = _radius_of_gyration(_active_pos(h, pos0) if prolif else pos0)

    t0 = time.perf_counter()
    sim.run(steps)
    sim.device  # ensure the device has finished (host sync via snapshot below)
    pos1 = _positions(sim)
    wall = time.perf_counter() - t0

    finite = bool(np.all(np.isfinite(pos1)))
    A1 = _footprint_active(h, pos1, p.z_substrate, band) if prolif \
        else _footprint_area(pos1, p.z_substrate, band)
    Rg1 = _radius_of_gyration(_active_pos(h, pos1) if prolif else pos1)

    lod = h.get("lod_updater")
    # LIVE cell count = cells with a non-negative cell_of_node (daughters that the
    # proliferation updater activated). NOT h["active"].sum() — the LOD classifier
    # shares + overwrites that mask, so it is not the liveness authority.
    n_live = int(len(np.unique(np.asarray(h["cell_of_node"])[
        np.asarray(h["cell_of_node"]) >= 0]))) if prolif else n_start
    n_active_final = n_live
    m = dict(
        n_cells=n_cells, n_start=n_start, steps=steps, active=active,
        use_lod=use_lod, prolif=prolif,
        device=type(sim.device).__name__, N_particles=int(sim.state.N_particles),
        finite=finite, wall_s=wall,
        A0=A0, A1=A1, A_over_A0=(A1 / A0 if A0 > 0 else float("nan")),
        Rg0=Rg0, Rg1=Rg1, Rg_ratio=(Rg1 / Rg0 if Rg0 > 0 else float("nan")),
        steps_per_s=(steps / wall if wall > 0 else float("nan")),
    )
    # REAL-TIME readout (PI 2026-06-11): the raw integrator time t_sim = dt·steps
    # is only ~µs (accelerated kinetics); t_real is the real-time equivalent via
    # the spreading-front acceleration factor S. MAPPING, not native real time.
    rt = map_realtime(steps, p.dt)
    m.update(rt.to_dict())
    if prolif:
        # Necrosis over LIVE cells ONLY (the parked dormant pool corrupts the
        # pool-wide LOD depth/centroid → spurious all-necrotic; see helper).
        n_necrotic_live = _live_necrotic_count(h, pos1, lod_cfg.necrotic_depth_um)
        m.update(dict(
            n_active_final=n_active_final,
            n_divisions=(prolif_upd.n_divisions if prolif_upd else 0),
            n_necrotic_live=n_necrotic_live,
            n_max=h["n_max"], p_div=p_div, div_every=div_every,
            auto_pdiv=auto_pdiv))
        if div_anchor is not None:
            m.update(div_anchor.to_dict())
    if lod is not None:
        if prolif:
            # Report the LIVE-cell necrosis (parked dormant pool excluded), not the
            # LOD classifier's pool-wide count (which mislabels every live cell).
            n_live = n_active_final
            m.update(dict(
                n_active=lod.n_active, n_inert=max(0, n_live - n_necrotic_live),
                n_necrotic=n_necrotic_live,
                active_frac=lod.active_frac,
                necrotic_frac=(n_necrotic_live / max(1, n_live)),
                k_freeze=h.get("lod_kfreeze")))
        else:
            m.update(dict(
                n_active=lod.n_active, n_inert=lod.n_inert,
                n_necrotic=lod.n_necrotic,
                active_frac=lod.active_frac, necrotic_frac=lod.necrotic_frac,
                k_freeze=h.get("lod_kfreeze")))
    elif prolif:
        # No LOD wired but prolif on → still report the live-cell necrosis count.
        m["n_necrotic"] = n_necrotic_live
    return m, sim


def _figure(metrics: dict, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lod = metrics["lod"]
    no = metrics.get("no_lod")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # (1) LOD cell fractions
    ax = axes[0]
    fr = [lod.get("active_frac", 1.0), lod.get("n_inert", 0) / max(1, lod["n_cells"]),
          lod.get("necrotic_frac", 0.0)]
    ax.bar(["active rim", "inert", "necrotic"], fr,
           color=["#2c7fb8", "#999999", "#d95f0e"])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("cell fraction")
    ax.set_title(f"activity-LOD classification (N={lod['n_cells']})")

    # (2) spread metrics
    ax = axes[1]
    ax.bar(["A/A0", "Rg ratio"], [lod["A_over_A0"], lod["Rg_ratio"]],
           color=["#31a354", "#756bb1"])
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_ylabel("ratio (final / initial)")
    ax.set_title("spheroid spread")

    # (3) wall-time A/B
    ax = axes[2]
    if no is not None:
        ax.bar(["no LOD", "LOD"], [no["wall_s"], lod["wall_s"]],
               color=["#bdbdbd", "#2c7fb8"])
        sp = no["wall_s"] / lod["wall_s"] if lod["wall_s"] > 0 else float("nan")
        ax.set_title(f"wall-time ({metrics['steps']} steps)  speedup {sp:.2f}x")
    else:
        ax.bar(["LOD"], [lod["wall_s"]], color=["#2c7fb8"])
        ax.set_title(f"wall-time ({metrics['steps']} steps)")
    ax.set_ylabel("wall-time [s]")

    # REAL-TIME readout: raw integrator time (~µs) AND the real-time equivalent
    # via the spreading-front accel factor (a MAPPING, not native real time).
    rt_line = (f"t = {lod.get('t_sim_human', '?')} sim  "
               f"≈ {lod.get('t_real_human', '?')} real  "
               f"(accel {lod.get('accel_factor', float('nan')):.0e}, spreading-front)")
    fig.suptitle(f"GPU-friendly DCM + activity-LOD  ({lod['device']}, "
                 f"N={lod['n_cells']}, {lod['N_particles']} nodes)\n{rt_line}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=40)
    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--subdivisions", type=int, default=1)
    ap.add_argument("--dt", type=float, default=1.0e-9)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--active", action="store_true",
                    help="wire the active rim traction")
    ap.add_argument("--ab", action="store_true",
                    help="also run WITHOUT LOD at the same N (measure the speedup)")
    ap.add_argument("--cpu", action="store_true", help="force CPU device")
    ap.add_argument("--cadence", type=int, default=500)
    ap.add_argument("--prolif", action="store_true",
                    help="enable live proliferation (pre-allocated pool + division)")
    ap.add_argument("--p-div", type=float, default=0.04,
                    help="per-eligible-rim-cell division prob per cadence (slow). "
                         "Ignored when --auto-pdiv is set (then DERIVED from physics)")
    ap.add_argument("--auto-pdiv", action="store_true",
                    help="DERIVE p_div from physics (real cell cycle vs the run's "
                         "real-time-equivalent via sim_realtime); not a free knob")
    ap.add_argument("--cell-cycle-h", type=float, default=20.0,
                    help="real cell-cycle / division time [hours] for --auto-pdiv "
                         "(MCF7 ~18-24 h; default 20)")
    ap.add_argument("--div-every", type=int, default=4000,
                    help="division-check period in steps (large = slow vs spreading)")
    ap.add_argument("--n-max", type=int, default=None,
                    help="proliferation pool size (default ceil(1.6*n_cells))")
    ap.add_argument("--r-cell-um", type=float, default=7.5,
                    help="per-cell radius [um]; >7.5 = coarse-grained tissue patch "
                         "(reaches R>150um necrosis onset at fewer cells)")
    args = ap.parse_args()

    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)

    p = ResolvedGpuDCM(subdivisions=args.subdivisions, dt=args.dt, seed=args.seed,
                       R_cell=args.r_cell_um * 1e-6)
    lod_cfg = ResolvedLOD(cadence=args.cadence)
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None

    print(f"[lod-run] N={args.n_cells} steps={args.steps} active={args.active} "
          f"ab={args.ab} prolif={args.prolif}")

    m_lod, _ = _run_one(p, args.n_cells, args.steps, active=args.active,
                        use_lod=True, lod_cfg=lod_cfg, device=device,
                        prolif=args.prolif, p_div=args.p_div,
                        div_every=args.div_every, n_max=args.n_max,
                        auto_pdiv=args.auto_pdiv, cell_cycle_h=args.cell_cycle_h)
    print(f"[LOD]    finite={m_lod['finite']} wall={m_lod['wall_s']:.2f}s "
          f"A/A0={m_lod['A_over_A0']:.3f} active={m_lod.get('n_active')}/"
          f"{args.n_cells} inert={m_lod.get('n_inert')} "
          f"necrotic={m_lod.get('n_necrotic')}")
    print(f"[TIME]   t_sim={m_lod['t_sim_human']} (raw integrator, dt={args.dt:g}s) "
          f"≈ {m_lod['t_real_human']} real-equivalent "
          f"(accel {m_lod['accel_factor']:.0e}, {m_lod['realtime_basis']})")
    if args.prolif:
        if args.auto_pdiv and ("auto_p_div" in m_lod):
            # PHYSICAL division-rate derivation header (Problem 1): show t_real,
            # T_cycle, expected divisions, and the DERIVED p_div — a prediction.
            print(f"[AUTO-PDIV] DERIVED p_div={m_lod['auto_p_div']:.3e} "
                  f"(= S·dt·div_every/T_cycle) — NOT tuned to A/A0")
            print(f"[AUTO-PDIV]   t_real={m_lod['div_t_real_human']} "
                  f"(=S·dt·steps)  T_cycle={m_lod['cell_cycle_human']} "
                  f"(={args.cell_cycle_h:g} h)  S={m_lod['div_accel_factor']:.0e}")
            print(f"[AUTO-PDIV]   div/cell={m_lod['div_per_cell']:.3e}  "
                  f"expected total divisions≈{m_lod['expected_divisions']:.3g} "
                  f"over n_rim≈{m_lod.get('n_start')} "
                  f"({m_lod['div_n_checks']} checks)")
        print(f"[PROLIF] divisions={m_lod.get('n_divisions')} "
              f"n_active: {m_lod.get('n_start')} -> {m_lod.get('n_active_final')} "
              f"(pool n_max={m_lod.get('n_max')}, p_div={m_lod.get('p_div'):.3e}, "
              f"div_every={args.div_every})")

    out = dict(steps=args.steps, n_cells=args.n_cells, active=args.active,
               prolif=args.prolif, lod=m_lod)

    if args.ab:
        dev2 = hoomd.device.CPU(notice_level=0) if args.cpu else None
        m_no, _ = _run_one(p, args.n_cells, args.steps, active=args.active,
                           use_lod=False, lod_cfg=lod_cfg, device=dev2,
                           prolif=args.prolif, p_div=args.p_div,
                           div_every=args.div_every, n_max=args.n_max,
                           auto_pdiv=args.auto_pdiv,
                           cell_cycle_h=args.cell_cycle_h)
        speedup = m_no["wall_s"] / m_lod["wall_s"] if m_lod["wall_s"] > 0 else float("nan")
        out["no_lod"] = m_no
        out["lod_speedup"] = speedup
        print(f"[no-LOD] finite={m_no['finite']} wall={m_no['wall_s']:.2f}s")
        print(f"[SPEEDUP] LOD vs no-LOD = {speedup:.2f}x "
              f"(N={args.n_cells}, {args.steps} steps, {m_lod['device']})")

    tag = f"n{args.n_cells}_s{args.steps}_{m_lod['device']}"
    (_OUT / f"metrics_{tag}.json").write_text(json.dumps(out, indent=2))
    fig_png = _OUT / "figs" / f"dcm_gpu_lod_{tag}.png"
    try:
        _figure(out, fig_png)
        print(f"[fig] {fig_png}")
    except Exception as e:  # matplotlib optional on the gbook headless
        print(f"[fig] skipped ({e})")
    print(f"[json] {_OUT / f'metrics_{tag}.json'}")


if __name__ == "__main__":
    main()
