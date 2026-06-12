"""TWO-STAGE DCM spheroid production (PI 2026-06-12, v2 — BIOLOGICAL aggregation).

PI rejected v1: spreading must NOT start from a just-placed lattice, AND the v1
"aggregation" (AggregationDrive central pull) was a forced hack — cells were dragged
to a point, not aggregated. Reference: SimuCell3D organoid cell aggregation. v2 makes
the aggregation EMERGENT (active-matter coalescence, search-and-capture):

  STAGE 1 · AGGREGATION (free-float hanging drop, all cells full physics):
    Loosely-seeded cells are SELF-PROPELLED (DcmActiveMotilitySPP: each cell carries a
    polarity reoriented with persistence tau_p by DcmPolarityUpdater) inside a SOFT
    SPHERICAL DROP (DcmDropConfinement: force-free interior, inward only at r>R_drop —
    the hanging-drop meniscus, NOT a central pull). Motile cells collide, cadherin-
    adhere (existing tent), and COALESCE into ONE aggregate that ROUNDS by the emergent
    adhesion-cortex-turgor surface tension. Run TO CONVERGENCE (Rg plateau + contact
    fraction > 0.7 + asphericity → 0), not a fixed step budget. NO central force, NO
    activity-LOD, NO rim/interior split — every cell wanders + has full physics.

  STAGE 2 · SPREADING: the CONVERGED aggregate (stage-1 final positions, same topology)
    is placed on the substrate via build init_pos and run with substrate + settling +
    active rim traction (drop + motility OFF). Hand-off verified (stage-2 frame 0 ==
    stage-1 final).

Diagnostics per frame: per-cell V/V0 (turgor held?), Rg, asphericity, maxZ, basal
footprint, AND contact-node fraction. gbook computes; the MacBook renders. Real-time
labels (S, t_sim, t_real) reported honestly.

Run (gbook, one size, to convergence):
    PYTHONPATH=. python ffn_sim/scripts/dcm_two_stage_production.py --n 200 \
        --agg-max-steps 10000000 --spread-steps 30000
Run (CPU dev smoke — confirm FORM only):
    PYTHONPATH=. python ffn_sim/scripts/dcm_two_stage_production.py --n 30 \
        --agg-max-steps 300000 --spread-steps 3000 --frames 6
"""

from __future__ import annotations

import argparse
import dataclasses
import pickle
import time
from pathlib import Path

import numpy as np
import hoomd

from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.cell.dcm_gpu_build import (
    DcmDropConfinement, ResolvedGpuDCM, build_gpu_dcm_simulation, pick_device)
from ffn_sim.cell.dcm_gpu_forces import DcmActiveMotilitySPP, DcmPolarityUpdater
from ffn_sim.cell.dcm_confluence import capture_positions
from ffn_sim.common.sim_realtime import map_realtime
from ffn_sim.scripts.dcm_native_capstone import _footprint_area, _effective_radius
from ffn_sim.scripts.dcm_cohesion_check import _contact_fraction

OUT = Path("ffn_sim/outputs/h_dcm_two_stage")
OUT.mkdir(parents=True, exist_ok=True)
UM = 1.0e6
S_ACCEL = 6.0e5            # accelerated-kinetics factor (common/sim_realtime basis)


def cell_volume(verts, tris):
    t = verts[tris]
    return abs(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()) / 6.0


def centroids(pos, ranges):
    return np.array([pos[lo:hi].mean(0) for lo, hi in ranges])


def _topdown_area(pos):
    """TOP-DOWN projected silhouette area A (xy convex hull of ALL nodes).

    This is the EXPERIMENTAL observable (PI's spheroid spreading assay is imaged
    from above): a spheroid projects its full shadow (≈ π R_spheroid² at t=0), and
    A grows ONLY when cells migrate OUT beyond that silhouette → A/A₀ ≥ 1 from then
    on. It is the correct denominator for A/A₀ = a + b/R + c/R².

    Distinct from ``_footprint_area`` (basal CONTACT footprint), which counts only
    substrate-touching nodes and gives a tiny t=0 denominator (the bottom contact
    cap), inflating A/A₀ ~10–20× as the contact patch grows toward the silhouette —
    a measurement artifact, NOT spreading. Never report basal A/A₀ as the assay value.
    """
    from scipy.spatial import ConvexHull
    xy = pos[:, :2]
    if xy.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(xy).volume)  # 2D hull "volume" == area
    except Exception:  # noqa: BLE001
        return float(np.ptp(xy[:, 0]) * np.ptp(xy[:, 1]))


def diagnostics(pos, ranges, cell_of_node, tris0, *, R, z0, V0, c_adh):
    """V/V0, Rg, asphericity, maxZ, top-down silhouette + basal footprint, contact."""
    vols = np.array([cell_volume(pos[lo:hi], tris0) for lo, hi in ranges])
    vr = vols / V0
    cen = centroids(pos, ranges)
    cc = cen.mean(0)
    r = np.linalg.norm(cen - cc, axis=1)
    Rg = float(np.sqrt((r ** 2).mean())) if len(cen) else 0.0
    if len(cen) >= 4:
        d = cen - cc
        G = (d[:, :, None] * d[:, None, :]).mean(0)
        ev = np.sort(np.linalg.eigvalsh(G))[::-1]
        asph = float((ev[0] - 0.5 * (ev[1] + ev[2])) / max(ev.sum(), 1e-30))
    else:
        asph = 0.0
    maxZ = float(pos[:, 2].max())
    topd = _topdown_area(pos)              # A — experimental top-down silhouette
    foot = _footprint_area(pos, z0, R)     # basal contact (NOT the assay A)
    cfrac = float(_contact_fraction(pos, cell_of_node, c_adh))
    return dict(VV0_mean=float(vr.mean()), VV0_min=float(vr.min()),
                VV0_max=float(vr.max()), Rg_um=Rg * UM, asphericity=asph,
                maxZ_um=maxZ * UM, topdown_um2=topd * UM * UM,
                footprint_um2=foot * UM * UM, contact_frac=cfrac)


def _rt(steps, dt):
    rt = map_realtime(int(steps), float(dt))
    return f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real"


# ---------------------------------------------------------------------------
# STAGE 1 — biological aggregation to CONVERGENCE
# ---------------------------------------------------------------------------
def aggregate(p, n_cells, *, dev, f_active, tau_p_min, reorient_every, R_drop_factor,
              k_wall, max_steps, frames, R, z0, V0, tris0, seed):
    h = build_gpu_dcm_simulation(p, n_cells, device=dev, active=False,
                                 with_substrate=False, settle_force=0.0)
    sim, ranges = h["sim"], h["ranges"]
    nbuilt = h["n_cells"]
    cell_of_node = h["cell_of_node"]
    # drop centre (fixed; interior is force-free so the exact value barely matters)
    snap0 = sim.state.get_snapshot()
    center = np.asarray(snap0.particles.position).mean(0)
    R_drop = (nbuilt ** (1.0 / 3.0) + R_drop_factor) * R
    # SPP motility (all live cells) + soft drop wall + polarity reorientation updater
    mot = DcmActiveMotilitySPP(cell_of_node=cell_of_node, ranges=ranges,
                               n_cells=nbuilt, f_active=f_active, ramp_steps=4000,
                               seed=seed)
    sim.operations.integrator.forces.append(mot)
    sim.operations.integrator.forces.append(
        DcmDropConfinement(center=center, R_drop=R_drop, k_wall=k_wall))
    dt_reorient = reorient_every * float(p.dt)
    D_r_eff = S_ACCEL / (tau_p_min * 60.0)           # accelerated rotational diffusion
    pol = DcmPolarityUpdater(motility=mot, D_r_eff=D_r_eff, dt_reorient=dt_reorient,
                             seed=seed)
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=pol, trigger=hoomd.trigger.Periodic(reorient_every)))

    print(f"STAGE 1 · AGGREGATION (n_cells={nbuilt}, f_active={f_active:.2e} N/cell, "
          f"tau_p={tau_p_min}min, R_drop={R_drop*UM:.0f}µm, drop+motility, ALL cells):",
          flush=True)
    sim.run(0)
    chunk = max(2000, max_steps // max(1, frames - 1))
    Frames, diags, steps = [], [], []

    def snap_diag():
        pos = capture_positions(sim)
        return pos, diagnostics(pos, ranges, cell_of_node, tris0,
                                R=R, z0=z0, V0=V0, c_adh=p.c_adh)

    pos, dg = snap_diag()
    Frames.append(pos.copy()); diags.append(dg); steps.append(int(sim.timestep))
    print(f"  [agg] f0 step 0: Rg={dg['Rg_um']:.1f}µm asph={dg['asphericity']:.3f} "
          f"contact={dg['contact_frac']:.2f} V/V0={dg['VV0_mean']:.3f}", flush=True)
    t0 = time.time()
    converged = False
    rg_hist = [dg["Rg_um"]]
    while int(sim.timestep) < max_steps:
        sim.run(chunk)
        pos, dg = snap_diag()
        if not np.isfinite(pos).all():
            print("  [agg] NON-FINITE — truncating", flush=True); break
        Frames.append(pos.copy()); diags.append(dg); steps.append(int(sim.timestep))
        rg_hist.append(dg["Rg_um"])
        print(f"  [agg] f{len(Frames)-1} step {sim.timestep}: Rg={dg['Rg_um']:.1f}µm "
              f"asph={dg['asphericity']:.3f} contact={dg['contact_frac']:.2f} "
              f"V/V0={dg['VV0_mean']:.3f}  ({_rt(sim.timestep, p.dt)})", flush=True)
        # CONVERGENCE for a TURGID-ROUNDED-CELL aggregate (not confluent/deformed):
        # the cohesive-ball equilibrium is contact≈0.45 (rounded cells touch ~half
        # their surface), NOT 0.7 (which needs deformed polyhedral cells). Converge
        # when the contact fraction has RISEN substantially AND plateaued (cells
        # reached their cohesive equilibrium), Rg stable, aggregate round.
        if len(rg_hist) >= 5:
            cf = [g["contact_frac"] for g in diags[-3:]]
            cf_plateau = (max(cf) - min(cf)) < 0.02
            rgw = rg_hist[-3:]
            rg_plateau = (max(rgw) - min(rgw)) / max(np.mean(rgw), 1e-9) < 0.01
            if cf_plateau and rg_plateau and dg["contact_frac"] > 0.35 \
                    and dg["asphericity"] < 0.10:
                converged = True
                print(f"  [agg] CONVERGED at step {sim.timestep} "
                      f"(contact {dg['contact_frac']:.2f} plateau, Rg plateau, "
                      f"asph {dg['asphericity']:.3f})", flush=True)
                break
    wall = time.time() - t0
    return dict(frames=Frames, diags=diags, steps=steps, wall_s=round(wall, 1),
                converged=converged, R_drop_um=R_drop * UM, f_active=f_active,
                tau_p_min=tau_p_min), h, ranges, Frames[-1].copy(), cell_of_node


# ---------------------------------------------------------------------------
# STAGE 2 — spreading from the converged aggregate
# ---------------------------------------------------------------------------
def _remesh_aware_diag(pos, cell_of_node, turgor, V0, *, ranges=None, z0=0.0, R=7.5e-6):
    """Remesh-safe diagnostics: TOP-DOWN A over LIVE nodes only (parked pool nodes
    excluded so the hull is not inflated), V/V0 from the CURRENT turgor faces grouped
    by face_cell, and the closed-manifold check on the current mesh. Topology-aware so
    it stays correct as SPLIT/COLLAPSE change the per-cell face set. Also carries the
    Rg/asphericity/footprint keys the surface renderer (dcm_two_stage_viz) reads,
    computed over LIVE cells."""
    from ffn_sim.cell.dcm_remesh import mesh_edges, enclosed_volume
    cof = np.asarray(cell_of_node)
    live = cof >= 0
    faces = np.asarray(turgor.faces)
    fc = np.asarray(turgor.face_cell)
    vr = []
    for c in np.unique(fc):
        cf = faces[fc == c]
        if cf.shape[0]:
            vr.append(enclosed_volume(pos, cf) / V0)
    vr = np.array(vr) if vr else np.array([1.0])
    _e, ef, _a, bnd = mesh_edges(faces)
    manifold = bool((not bnd.any()) and (ef >= 0).all())
    # per-LIVE-cell centroids → Rg + asphericity (renderer keys); footprint over live
    Rg = asph = 0.0
    if ranges is not None:
        cents = np.array([pos[a:b].mean(0) for (a, b) in ranges
                          if cof[a] >= 0])
        if len(cents):
            cc = cents.mean(0); rr = np.linalg.norm(cents - cc, axis=1)
            Rg = float(np.sqrt((rr ** 2).mean()))
            if len(cents) >= 4:
                dd = cents - cc
                G = (dd[:, :, None] * dd[:, None, :]).mean(0)
                ev = np.sort(np.linalg.eigvalsh(G))[::-1]
                asph = float((ev[0] - 0.5 * (ev[1] + ev[2])) / max(ev.sum(), 1e-30))
    foot = _footprint_area(pos[live], z0, R) if live.any() else 0.0
    return dict(topdown_um2=_topdown_area(pos[live]) * UM * UM,
                footprint_um2=foot * UM * UM, Rg_um=Rg * UM, asphericity=asph,
                VV0_mean=float(vr.mean()), VV0_min=float(vr.min()),
                VV0_max=float(vr.max()), maxZ_um=float(pos[live, 2].max()) * UM,
                manifold=manifold, n_faces=int(faces.shape[0]),
                n_live=int(live.sum()))


def spread(p, n_cells, *, dev, init_pos, steps, frames, R, z0, V0, tris0,
           node_face_contact=False, n_pool=0, remesh=False, remesh_period=2000,
           remesh_max_ops=24):
    h = build_gpu_dcm_simulation(p, n_cells, device=dev, active=True,
                                 with_substrate=True, init_pos=init_pos,
                                 node_face_contact=node_face_contact, n_pool=n_pool)
    sim, ranges = h["sim"], h["ranges"]
    cell_of_node = h["cell_of_node"]
    turgor = h["turgor"]
    remesh_action = None
    if remesh:
        from ffn_sim.cell.dcm_remesh_updater import attach_remesh_updater
        remesh_action, _ = attach_remesh_updater(
            h, period=remesh_period, max_ops=remesh_max_ops)
    print(f"STAGE 2 · SPREADING (from converged aggregate, active rim traction; "
          f"contact={'node-FACE' if node_face_contact else 'node-node'}, "
          f"remesh={'ON p=%d pool=%d' % (remesh_period, n_pool) if remesh else 'OFF'}):",
          flush=True)
    sim.run(0)
    spf = max(1, steps // max(1, frames - 1))
    Frames, diags, st = [], [], []
    # remesh changes the per-cell topology, so V/V0 + top-down come from the
    # remesh-aware diagnostic (current faces, live nodes); the static-topology
    # diagnostics() (Rg/asph/contact from ranges) is overlaid for the non-remesh path.
    use_rd = remesh or node_face_contact or n_pool > 0

    def snap_diag():
        pos = capture_positions(sim)
        if use_rd:
            return pos, _remesh_aware_diag(pos, cell_of_node, turgor, V0,
                                           ranges=ranges, z0=z0, R=R)
        return pos, diagnostics(pos, ranges, cell_of_node, tris0,
                                R=R, z0=z0, V0=V0, c_adh=p.c_adh)

    pos, dg = snap_diag()
    Frames.append(pos.copy()); diags.append(dg); st.append(int(sim.timestep))
    A0_top = max(dg["topdown_um2"], 1e-9)   # A₀ = spheroid top-down shadow at t=0
    print(f"  [spread] f0 step 0: A/A0=1.00 (A0_topdown={A0_top:.0f}µm²) "
          f"maxZ={dg['maxZ_um']:.1f}µm V/V0={dg['VV0_mean']:.3f} "
          f"faces={dg.get('n_faces','-')}", flush=True)
    t0 = time.time()
    for f in range(1, frames):
        sim.run(spf)
        pos, dg = snap_diag()
        if not np.isfinite(pos).all():
            print("  [spread] NON-FINITE — truncating", flush=True); break
        Frames.append(pos.copy()); diags.append(dg); st.append(int(sim.timestep))
        print(f"  [spread] f{f} step {sim.timestep}: "
              f"A/A0={dg['topdown_um2']/A0_top:.3f} (top-down, the assay obs) "
              f"maxZ={dg['maxZ_um']:.1f}µm V/V0={dg['VV0_mean']:.3f} "
              f"manifold={dg.get('manifold','-')} faces={dg.get('n_faces','-')}",
              flush=True)
    if remesh_action is not None:
        print(f"  [spread] remesh totals: {remesh_action.totals} "
              f"({remesh_action.n_acts} acts)", flush=True)
    return dict(frames=Frames, diags=diags, steps=st, wall_s=round(time.time() - t0, 1),
                A0_topdown_um2=A0_top), h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--agg-max-steps", type=int, default=10_000_000,
                    help="aggregation safety cap; runs to CONVERGENCE before this")
    ap.add_argument("--spread-steps", type=int, default=30000)
    ap.add_argument("--frames", type=int, default=16)
    ap.add_argument("--agg-spacing", type=float, default=2.6,
                    help="loose seed spacing ·R (within cadherin reach; motility coalesces)")
    ap.add_argument("--f-active", type=float, default=1.6e-10,
                    help="per-cell active self-propulsion [N] (accelerated kinetics)")
    ap.add_argument("--tau-p-min", type=float, default=20.0, help="polarity persistence [min]")
    ap.add_argument("--reorient-every", type=int, default=200)
    ap.add_argument("--agg-dt", type=float, default=1.0e-8,
                    help="aggregation timestep [s]. Larger than the spreading dt (1e-9) "
                         "is stable while cells are loose (stiff contact inactive) and "
                         "covers the slow aggregation in feasible step counts.")
    ap.add_argument("--agg-k-edge", type=float, default=None,
                    help="cortex edge-spring stiffness during aggregation [N/m]. "
                         "Default keeps the validated stiff band (1e-3); a softer value "
                         "(~5e-5, the confluent regime) lets cells DEFORM into dense "
                         "polygonal contact (higher contact fraction, rounding) as a "
                         "real organoid does — SimuCell3D deformable-cell aggregation.")
    ap.add_argument("--r-cell-um", type=float, default=7.5)
    ap.add_argument("--seed", type=int, default=7)
    # Phase 2 step 2: node-face contact + remeshing on the SPREAD stage.
    ap.add_argument("--node-face-contact", action="store_true",
                    help="STAGE 2 uses the SimuCell3D node-vs-face penalty contact "
                         "(no shell interpenetration) instead of the node-node tent.")
    ap.add_argument("--remesh", action="store_true",
                    help="attach the low-cadence DcmRemeshUpdater during spreading "
                         "(keeps edges in [l_min,3·l_min]; prevents sliver blow-up).")
    ap.add_argument("--remesh-period", type=int, default=2000,
                    help="steps between remesh passes (low cadence).")
    ap.add_argument("--n-pool", type=int, default=0,
                    help="dormant node-pool size per build (SPLIT headroom). Default "
                         "0; set ~ (frames·max_ops·2) headroom for a long spread.")
    ap.add_argument("--remesh-max-ops", type=int, default=24,
                    help="max mutations per remesh pass.")
    ap.add_argument("--spread-dt", type=float, default=1.0e-3,
                    help="STAGE-2 spreading timestep [s]. The physiological-γ "
                         "foundation raised the overdamped CFL ceiling (~0.25 s) so a "
                         "much larger dt than the legacy 1e-9 is stable; default 1e-3 "
                         "(3e-3 ceiling-validated, 1e-3 is a safe margin) so the spread "
                         "actually progresses in feasible step counts.")
    # node-face contact stiffness CALIBRATION overrides (Phase-2 re-tune; the
    # ResolvedGpuDCM defaults 2e8/8e8 are node-NODE patch_area values, ~5x too stiff /
    # ~10x over the mcf7_p0.xml dimensionless ratio for the A_face law). Anchor:
    # mcf7 xi≈4e7 Pa/m at our K, adhesion via cadherin (mcf7 omega=0) — these let the
    # calibration sweep the lit-anchored ratio without editing the frozen defaults.
    ap.add_argument("--rep-strength", type=float, default=None,
                    help="override node-face repulsion xi [Pa/m] (default: p value).")
    ap.add_argument("--adh-strength", type=float, default=None,
                    help="override node-face adhesion omega [Pa/m] (default: p value).")
    args = ap.parse_args()

    R = args.r_cell_um * 1e-6
    V0 = (4.0 / 3.0) * np.pi * R ** 3
    z0 = 0.0
    _v, _e, tris0 = icosphere_mesh(R, 1)
    dev = pick_device(None)
    is_gpu = type(dev).__name__.endswith("GPU")
    print(f"[two-stage v2 BIOLOGICAL] N={args.n} R_cell={args.r_cell_um}µm "
          f"device={'GPU' if is_gpu else 'CPU'} S={S_ACCEL:.0e}\n", flush=True)

    agg_over = dict(R_cell=R, spacing_factor=args.agg_spacing, dt=args.agg_dt)
    if args.agg_k_edge is not None:
        agg_over["k_edge"] = args.agg_k_edge       # soft cortex → deformable cells
    p_agg = dataclasses.replace(ResolvedGpuDCM(seed=args.seed), **agg_over)
    s1, h1, ranges, agg_pos, _con = aggregate(
        p_agg, args.n, dev=dev, f_active=args.f_active, tau_p_min=args.tau_p_min,
        reorient_every=args.reorient_every, R_drop_factor=2.0, k_wall=1.0e-3,
        max_steps=args.agg_max_steps, frames=args.frames, R=R, z0=z0, V0=V0,
        tris0=tris0, seed=args.seed)
    # place the converged aggregate onto the dish (bottom at z0, centred in xy)
    agg_pos[:, 2] -= (agg_pos[:, 2].min() - z0)
    cxy = centroids(agg_pos, ranges).mean(0)[:2]
    agg_pos[:, 0] -= cxy[0]; agg_pos[:, 1] -= cxy[1]
    print(f"  -> aggregated: Rg {s1['diags'][0]['Rg_um']:.1f}→{s1['diags'][-1]['Rg_um']:.1f}µm, "
          f"asph {s1['diags'][0]['asphericity']:.3f}→{s1['diags'][-1]['asphericity']:.3f}, "
          f"contact {s1['diags'][0]['contact_frac']:.2f}→{s1['diags'][-1]['contact_frac']:.2f}, "
          f"converged={s1['converged']}\n", flush=True)

    # stage 2 uses the cohesive contact bands (ResolvedGpuDCM defaults) but the SAME
    # spacing_factor as stage 1 so _cluster_centers yields the SAME cell count/topology
    # → the converged aggregate's init_pos lines up (the lattice positions it places
    # are overridden by init_pos anyway, so spacing here only fixes the topology).
    p = dataclasses.replace(ResolvedGpuDCM(seed=args.seed), R_cell=R,
                            spacing_factor=args.agg_spacing, dt=args.spread_dt)
    # node-face contact stiffness calibration overrides (Phase-2 re-tune)
    _contact_over = {}
    if args.rep_strength is not None:
        _contact_over["rep_strength"] = args.rep_strength
    if args.adh_strength is not None:
        _contact_over["adh_strength"] = args.adh_strength
    if _contact_over:
        p = dataclasses.replace(p, **_contact_over)
        print(f"[contact re-tune] {_contact_over}", flush=True)
    s2, h2 = spread(p, args.n, dev=dev, init_pos=agg_pos, steps=args.spread_steps,
                    frames=args.frames, R=R, z0=z0, V0=V0, tris0=tris0,
                    node_face_contact=args.node_face_contact, n_pool=args.n_pool,
                    remesh=args.remesh, remesh_period=args.remesh_period,
                    remesh_max_ops=args.remesh_max_ops)
    tr = h2["traction"]

    out = dict(
        n_cells=int(h2["n_cells"]), R_cell=R, V0=V0, z0=z0, S=S_ACCEL,
        tris0=tris0, ranges=h2["ranges"], device="GPU" if is_gpu else "CPU",
        agg=s1, spread=s2, p_agg=p_agg, p=p,
        rim_params=dict(r_neigh=float(tr.r_neigh), max_neigh=int(tr.max_neigh),
                        contact_band=float(tr.contact_band), R=R, z0=z0),
        realtime=dict(agg=_rt(s1["steps"][-1], p_agg.dt),
                      spread=_rt(s2["steps"][-1], p.dt), accel=S_ACCEL))
    path = OUT / f"two_stage_n{args.n}.pkl"
    with open(path, "wb") as fh:
        pickle.dump(out, fh)
    print(f"\nwrote {path} (agg {s1['wall_s']}s converged={s1['converged']} + "
          f"spread {s2['wall_s']}s)", flush=True)
    print(f"AGGREGATION biological check: contact {s1['diags'][-1]['contact_frac']:.2f} "
          f"(>0.7?), asph {s1['diags'][-1]['asphericity']:.3f} (→0?), "
          f"Rg {s1['diags'][0]['Rg_um']:.0f}→{s1['diags'][-1]['Rg_um']:.0f}µm; "
          f"VOLUME held V/V0={s1['diags'][-1]['VV0_mean']:.3f}", flush=True)
    A0t = s2["A0_topdown_um2"]
    print(f"SPREADING A/A0 (TOP-DOWN silhouette = the experimental observable): "
          f"{s2['diags'][0]['topdown_um2']/A0t:.3f} → "
          f"{s2['diags'][-1]['topdown_um2']/A0t:.3f}  "
          f"(A0={A0t:.0f}µm²; basal-contact ratio is NOT the assay value)", flush=True)


if __name__ == "__main__":
    main()
