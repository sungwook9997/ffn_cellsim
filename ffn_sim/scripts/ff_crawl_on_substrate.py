"""FF crawl — a polarized whole cell CRAWLS across a substrate in real physical time (piece-4 + B).

Answers PI 2026-07-06: the cell must be seen DEFORMING and MOVING, membrane + filaments + nucleus co-moving —
not a thin filopodial appendage, and not equilibrate-between-ticks quasi-statics. Here:

  - **Leading-edge protrusion** (Mogilner-Oster ratchet as a lamellipodial stress on the FRONT cortex/membrane,
    throttled by the emergent counter-load) pushes the front of the CELL BODY forward → the front deforms out.
  - **Basal FA integrin catch-slip clutches** grip the substrate → the protrusion is converted to net forward
    COM translocation by TRACTION (KU-2.4/2.5). Rear clutches rupture under load, front ones re-form (nascent
    adhesion) → a clutch treadmill.
  - **Real η-dynamics** (γ = units.fiber_point_drag, η = 65.9 Pa·s) → physical-time trajectory; the crawl SPEED
    is a measurable observable (validate vs retrograde flow 10–100 nm/s, KU-3.12).

Decisive adversarial audit (``--audit``): run again with clutches OFF. Polymerization alone is an INTERNAL
force (Newton) → the no-clutch net COM drift must be ≈ 0 (front deforms in place, cell does not translocate);
only WITH traction does the cell move. If the no-clutch run also translocates, the "crawl" is an artifact.

    python -m ffn_sim.scripts.ff_crawl_on_substrate --device cuda:0 --cortex-fil 900 --steps 8000 --audit
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp
from scipy.spatial import ConvexHull, cKDTree

from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0,
                                    VMIN_FRAC, NMIIA_MINIFIL_STALL_PN)
from ffn_sim.ff.network_warp import (_zero, link_spring_kernel, myosin_kernel, turgor_kernel,
                                     reshape_kernel, nucleus_shell_kernel, _seed_nucleus_cloud,
                                     substrate_plane_kernel, simulate_whole_cell_compression_on_device)
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.microtubule import build_microtubule_aster, merge_aster_into_cortex
from ffn_sim.ff.ff_virial_stress import cortex_node_stress, cortex_node_areal_strain, face_areas
from ffn_sim.ff.fa_clutch_warp import (clutch_spring_kernel, clutch_catchslip_kmc_kernel, resolve_clutch)
from ffn_sim.ff.motility_warp import (axpy_physical_kernel, leading_edge_push_kernel, protrusion_reaction_kernel,
                                      spreading_push_kernel, spreading_reaction_kernel, gravity_kernel,
                                      cortex_volume_kernel, xl_turnover_kernel, actin_assembly_kernel,
                                      barbed_end_growth_kernel, directed_front_growth_kernel,
                                      pointed_end_depoly_kernel, fiber_treadmill_kernel,
                                      clutch_slip_accumulate_kernel, clutch_slip_traction_kernel,
                                      sum_pos_kernel, sum_radius_kernel,
                                      volume_gradient, physical_node_gammas, crawl_cfl_dt)
from ffn_sim.ff.units import ETA_CYTOPLASM
from ffn_sim.ff.implicit_ff import implicit_step_current
from ffn_sim.ff.polymerization_warp import resolve_polymerization
from ffn_sim.ff.myosin_linear import minifilament_kernel, resolve_myosin
from ffn_sim.ff.substrate import resolve_substrate, substrate_anchor_equilibrium_kernel
from ffn_sim.ff.fa_maturation import (MaturationParams, talin_vinculin_kernel, fa_growth_kernel,
                                      clutch_load_kernel, fa_disassemble_kernel)
from ffn_sim.ff.piezo import resolve_piezo, p_open
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane

# leading-edge protrusion magnitude — DERIVED, not tuned (KU-3.6 / KU-3.18):
FIL_AREAL_DENSITY_PER_UM2 = 100.0      # cortical/lamellipodial actin areal density [1/µm²] (KB-3.18)
F_STALL_ACTIN_PN = 4.0                  # per-filament Brownian-ratchet stall force [pN] (KB-3.6: 2–5 pN)


def _fps_subsample(pts: np.ndarray, k: int) -> np.ndarray:
    """Farthest-point subsample: pick ``k`` well-spread points from ``pts`` (M,3). Deterministic (no RNG).
    Used to place DISCRETE focal-adhesion sites — real FAs are ~10² discrete integrin clusters, NOT one clutch
    per cortex node; at native Nc the per-node model dilutes per-clutch load ~100× → clutches never turn over →
    symmetric pinning → no crawl. Capping to ~physiological FA count restores the per-clutch load that drives
    the catch-slip treadmill. Returns local indices into ``pts``."""
    n = pts.shape[0]
    if k >= n:
        return np.arange(n)
    sel = [n // 2]                                            # deterministic seed point
    d = np.full(n, np.inf)
    for _ in range(1, k):
        d = np.minimum(d, np.sum((pts - pts[sel[-1]]) ** 2, axis=1))
        sel.append(int(np.argmax(d)))
    return np.asarray(sel, dtype=np.int64)


def build(n_cortex_fil=900, seed=7, contact_h=0.6, front_frac=0.5, phat=(1.0, 0.0, 0.0),
          length_dist="mono", microtubules=False, n_mt=40, L_mt_um=6.0,
          from_resting=False, relax_steps=6000, relax_device="cpu",
          n_myo_ratio=160, f_excess=0.0, f_myo=NMIIA_MINIFIL_STALL_PN, n_fa=0,
          ecm=False, ecm_fibers=1500, ecm_depth=4.0, ecm_capture=1.0, ecm_lp_um=20.0):
    """Polarized cell on a substrate: cortex + nucleus + membrane, basal FA clutches on the contact cap, and a
    FRONT cap (nodes with (x−com)·phat > front_frac·R) that carries the leading-edge protrusion.

    ``length_dist`` — cortex filament contour-length model: ``"mono"`` (every filament = L_filament_um, the
    γ-validated default) or ``"exponential"`` (KB-3.18-distributed 1–10 µm, mean-preserving; real cortical
    F-actin is length-distributed, not identical).

    ``microtubules`` (Thread-C Stage 1) — merge an MT aster (``n_mt`` stiff κ=KAPPA_MT tubes of contour length
    ``L_mt_um`` radiating from one MTOC at the centroid) into the cortex fiber network, so MT bending rides the
    SAME implicit K as the cortex at 300× the bending stiffness. Node layout becomes
    ``[cortex(Nc) ; MT_arms(Nmt) ; MTOC(1) ; nucleus(n_nuc)]`` with ``Ne = Nc+Nmt+1`` elastic nodes; the cortex
    block ``[:Nc]`` keeps ALL surface physics. ``microtubules=False`` → no-op passthrough (Ne==Nc, empty hub) so
    the OFF path is bit-identical (needs the ``--implicit`` GPU solver; explicit is dev-only).

    A3 (2026-07-08) — ``from_resting`` starts the ADHERENT run from the VALIDATED resting checkpoint cell rather
    than a fresh un-relaxed geometry (the physiological-baseline HARD rule: adhere FROM the resting state, don't
    let a null baseline emergently align). It matches the ``ff_resting_full_compartment.py`` recipe
    (``n_myo_ratio``=10 dense myosin, membrane ``f_excess``=0.25 reservoir, ``turgor_every``=50) and PRE-RELAXES
    to the resting turgor set-point (ΔP=TURGOR_DP0=40 Pa → γ=ΔP·R/2, biphasic drained solid K_drained=300 Pa) via
    the SAME ``simulate_whole_cell_compression_on_device`` the checkpoint uses, seeding the nucleus into the
    IDENTICAL ``[cortex ; MT_arms ; MTOC ; nucleus]`` layout, so the relaxed positions map back with no
    reconciliation. Substrate contact (z_sub / basal cap / anchors / front) is then derived from the RELAXED
    cortex geometry.

    The MT aster is included only when ``microtubules`` is set (as in the full checkpoint); at the resting
    set-point the MT arms are mechanically DECOUPLED from the cortex (only MTOC↔arm hub crosslinks, no MT↔cortex
    bond, and the tip↔cortex contact is inactive since L_mt<R), so the resting cortical-mechanics baseline (ΔP,
    γ, R_eq, V/V0) is MT-independent — ``--from-resting`` alone reproduces that baseline; the FULL checkpoint
    compartment adds the aster. Production is native (``--cortex-fil 38000 --from-resting --microtubules
    --implicit`` on the A5000); the CPU path is the small-N dev demo."""
    rng = np.random.default_rng(seed)
    n_myo = max(1, n_cortex_fil // n_myo_ratio)               # from_resting → dense (//10) checkpoint myosin
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_cortex_fil, n_xl=n_cortex_fil,
                                  n_myo=n_myo, length_dist=length_dist, rng=rng)
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    R = cx.R0_mean
    Nc = cx.net.n_nodes
    phat = np.asarray(phat, np.float64); phat = phat / np.linalg.norm(phat)
    nuc = resolve_nucleus(R_nuc_um=0.70 * R, n_beads=3000)     # 0.70R: MCF7 nucleus (Moore2016 0.68-0.77, PI-ratified 2026-07-07;
    #                                                            Ø~12µm ≈ 0.8R, N:C 1.9 ~50% cell vol) — was 0.25R (~3× too small)
    mem = resolve_membrane(f_excess=f_excess)                 # from_resting → 0.25 reservoir (checkpoint); else plateau
    c0 = cx.net.pos.mean(0)
    nuc_pos = _seed_nucleus_cloud(c0, nuc.R_nuc_um, nuc.n_beads, np.random.default_rng(seed + 2))
    # --- Thread-C Stage 1: MT aster merges into the cortex elastic network (one shared implicit K) ---
    aster = build_microtubule_aster(centre=c0, n_mt=n_mt, L_mt_um=L_mt_um) if microtubules else None
    m = merge_aster_into_cortex(cx.net, aster, k_hub_pn_um=float(cx.xl_k.max()))
    merged = m["net"]                                          # [cortex(Nc) ; MT_arms(Nmt)]; MTOC appended below
    Ne = m["Ne"]                                               # elastic node count Nc+Nmt+1 (MTOC); == Nc when OFF
    # node layout [cortex(Nc) ; MT_arms(Nmt) ; MTOC(1) ; nucleus(n_nuc)] — mtoc_pos empty ⇒ OFF pos_all unchanged
    pos_all = np.concatenate([merged.pos, m["mtoc_pos"], nuc_pos], 0)
    resting = None
    if from_resting:
        # PRE-RELAX to the validated resting checkpoint (same builder the checkpoint uses → same node layout).
        pos_relaxed, mrelax = simulate_whole_cell_compression_on_device(
            cx, f_myo, strain=0.0, nucleus=nuc, membrane=mem, microtubule=aster,
            pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=300.0, n_steps=relax_steps,
            turgor_every=50, device=relax_device)                 # turgor_every=50 = exact checkpoint recipe
        pos_all = np.ascontiguousarray(pos_relaxed, np.float64)
        merged.pos = pos_all[:merged.pos.shape[0]].copy()     # keep net.pos consistent (run() reads it for V0 + faces)
        resting = dict(dP_Pa=float(mrelax["dP_turgor_Pa"]), gamma_mN_m=float(mrelax["gamma_apparent_mN_m"]),
                       V_over_V0=float(mrelax["V_over_V0"]), R_eq_um=float(mrelax["R_eq_um"]))
    # contact + polarity from the (possibly relaxed) CORTEX geometry pos_all[:Nc]
    cortex_pos = pos_all[:Nc]
    c = cortex_pos.mean(0)
    z_sub = float(cortex_pos[:, 2].min())
    basal = np.where(cortex_pos[:, 2] < z_sub + contact_h)[0]
    if n_fa and 0 < n_fa < basal.size:                       # DISCRETE FA sites (physiological ~10² integrin clusters,
        basal = basal[_fps_subsample(cortex_pos[basal], n_fa)]   # not one-per-node) → physical per-clutch load → turnover
    anchors = cortex_pos[basal].copy(); anchors[:, 2] = z_sub
    proj = (cortex_pos - c) @ phat
    front = np.where(proj > front_frac * R)[0]
    node_area = 4.0 * np.pi * R**2 / Nc                       # mean cortex area per node [µm²]
    f_pro = FIL_AREAL_DENSITY_PER_UM2 * node_area * F_STALL_ACTIN_PN   # per front-node protrusive force [pN]
    # S4 (Phase B) — build a collagen-I (Mikado) ECM slab UNDER the cell + attach the basal clutches to fiber nodes.
    # BUILD-side only (the FA↔ECM force coupling in run() is the next, PI-overseen step — see
    # FF_S4_INTEGRATION_DESIGN_2026-07-09); with --ecm alone the crawl still runs on the fixed-substrate clutch, so
    # the validated path is UNTOUCHED. This just makes the matrix + the nascent bonds available in S.
    ecm_net = ecm_node = None
    if ecm:
        from ffn_sim.ff import units as U
        from ffn_sim.ff.ecm_mikado import build_mikado_network
        from ffn_sim.ff.fa_ecm import attach_clutches_to_ecm
        xy = cortex_pos[:, :2]
        lo = np.array([xy[:, 0].min() - 2.0, xy[:, 1].min() - 2.0, z_sub - ecm_depth])
        hi = np.array([xy[:, 0].max() + 2.0, xy[:, 1].max() + 2.0, z_sub + 0.5])  # slab from below z_sub to the basal cap
        ecm_net = build_mikado_network(lo, hi, n_fibers=ecm_fibers, pin_face="z_lo",  # bulk collagen below = pinned
                                       kappa=U.KBT * float(ecm_lp_um),                 # κ=kBT·Lp; Lp is the biphasic-sweep stiffness knob (PI-gated)
                                       rng=np.random.default_rng(seed + 11))
        ecm_node = attach_clutches_to_ecm(cortex_pos[basal.astype(np.int64)], ecm_net.net.pos, capture_um=ecm_capture)
    return dict(cx=cx, net=merged, Nc=Nc, Ne=Ne, n_nuc=nuc_pos.shape[0], pos_all=pos_all, nuc=nuc,
                mem=mem, basal=basal.astype(np.int32), anchors=anchors, z_sub=z_sub, R=R, c=c,
                phat=phat, front=front.astype(np.int32), front_frac=front_frac, f_pro=f_pro,
                mtoc_idx=m["mtoc_idx"], hub_i=m["hub_i"], hub_j=m["hub_j"], hub_k=m["hub_k"],
                hub_rest=m["hub_rest"], n_mt=int(m.get("n_mt", 0)), resting=resting,
                ecm=ecm_net, ecm_node=ecm_node)


def ecm_remodel_metrics(pos0, posf, ecm_node, basal, cortex_final, foff, near_um=8.0):
    """S6 collagen-remodel decomposition — separate COHERENT physiological remodeling from passive settling.

    The raw gripped-node displacement magnitude (the old "266 nm") is undirected and dominated by the cell
    settling onto the compliant slab. This decomposes it into physiological components so the *coherent inward*
    remodel is isolated from settling, and measures fiber reorientation toward the footprint.

    Args:
        pos0, posf: (Ne,3) collagen node positions, frame-0 and final [µm].
        ecm_node: (M,) clutch→collagen node index (<0 = clutch not attached to any fiber).
        basal: (M,) cortex (actin) node index each clutch grips through.
        cortex_final: (N,3) final cortex node positions [µm].
        foff: (F+1,) collagen fiber node offsets.
        near_um: only fibers whose midpoint is within this horizontal distance of the footprint centroid count
            toward the alignment index (a measurement window ≈ the cell radius, NOT a tuned parameter).

    Returns:
        dict of nm-scale components (recruit toward the bound actin, footprint-radial densification, vertical
        settling, undirected magnitude, coherence∈[0,1]) + the radial fiber-alignment index (frame0/final/Δ).
    """
    ecm_node = np.asarray(ecm_node); bnd = ecm_node >= 0
    out = {"n_grip": int(bnd.sum())}
    if bnd.any():
        gr = ecm_node[bnd]
        p0 = pos0[gr]; pf = posf[gr]; d = pf - p0                         # gripped-node displacement (Ng,3) [µm]
        act = cortex_final[np.asarray(basal)[bnd]]                        # final actin node each clutch grips through
        tc = act - p0; tc /= (np.linalg.norm(tc, axis=1, keepdims=True) + 1e-12)
        out["recruit_nm"] = float((d * tc).sum(1).mean()) * 1e3           # projection onto the toward-actin direction (f93683a-correct)
        fc = p0[:, :2].mean(0)                                            # footprint centroid (xy)
        tr = fc[None, :] - p0[:, :2]; tr /= (np.linalg.norm(tr, axis=1, keepdims=True) + 1e-12)
        out["dens_nm"] = float((d[:, :2] * tr).sum(1).mean()) * 1e3       # centripetal (footprint-radial-inward) densification
        out["dz_nm"] = float(d[:, 2].mean()) * 1e3                        # vertical settling component (down = −)
        out["mag_nm"] = float(np.linalg.norm(d, axis=1).mean()) * 1e3     # undirected magnitude (the old headline number)
        out["coh"] = float(np.linalg.norm(d.mean(0)) / (np.linalg.norm(d, axis=1).mean() + 1e-12))   # coherence [0,1]
    else:
        out.update(recruit_nm=0.0, dens_nm=0.0, dz_nm=0.0, mag_nm=0.0, coh=0.0)
        fc = pos0[:, :2].mean(0)
    foff = np.asarray(foff)
    def _rai(pos):                                                       # radial-alignment index |t̂_xy·r̂| near the footprint
        vals = []
        for f in range(len(foff) - 1):
            a, b = int(foff[f]), int(foff[f + 1])
            if b - a < 2:
                continue
            t = pos[b - 1] - pos[a]; n = np.linalg.norm(t[:2])           # end-to-end tangent (xy)
            mid = pos[(a + b) // 2, :2]; r = fc - mid; rn = np.linalg.norm(r)
            if n < 1e-9 or rn < 1e-9 or rn > near_um:
                continue
            vals.append(abs(float((t[:2] @ r) / (n * rn))))
        return float(np.mean(vals)) if vals else 0.0
    out["rai0"] = _rai(pos0); out["raif"] = _rai(posf); out["drai"] = out["raif"] - out["rai0"]
    return out


def run(S, *, steps=600000, dt=None, safety=0.1, f_myo=NMIIA_MINIFIL_STALL_PN, clutches=True, protrude=True,
        spread=False, rupture=True, gravity=True, delta_rho=55.0, koff_xl=0.4, implicit=False, dt_impl=1.0e-2,
        assembly=False, k_assembly=0.4, growth=False, myosin_linear=False, substrate_E=0.0, fa_maturation=False,
        treadmill=False, rear_depoly=False, flow=False, cortex_treadmill=True, v_retro_um_s=0.03, bulk_drag=False,
        com_drag=False, ecm_regrip=False, refresh_every=50, reshape_every=20, kmc_every=2000, xl_turn_every=50,
        assembly_every=20, record_every=2500, device="cpu"):
    """PHYSICAL-TIME crawl via a single EXPLICIT overdamped loop (CFL-stable — cannot diverge) + cortical
    crosslink turnover. Every force ticks at the same physical ``dt`` (= safety·γ_min/kmax, ~5.5 µs — set by
    the stiff α-actinin crosslinks) and every node (cortex + membrane law + nucleus) co-moves in real time:

      x_i ← x_i + (dt/γ_i)·F_i,   F = bending + crosslink + myosin + turgor + membrane(inward) + nucleus
                                       + substrate + clutch + leading-edge protrusion

    γ_i = physical per-node cytoplasm drag (units.fiber_point_drag, η=65.9 Pa·s), so trajectory time and crawl
    speed are physical. **Crosslink turnover** (α-actinin k_off≈``koff_xl``/s, on-device Maxwell rest-length
    relaxation) fluidizes the cortex so the leading-edge protrusion can deform the cell body (a crawling cell's
    cortex flows). ``clutches=False`` = the no-traction audit (protrusion is internal → net COM drift ≈ 0).

    Explicit + tiny dt ⇒ many steps ⇒ built for the GPU (A5000). Returns frames, COM trajectory, crawl metrics."""
    cx = S["cx"]; net = S["net"]; Nc = S["Nc"]; Ne = S["Ne"]; n_nuc = S["n_nuc"]; N = Ne + n_nuc; d = device
    cp = resolve_clutch(); poly = resolve_polymerization()
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(net.fiber_offsets) - 1
    soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    n_cortex_fib = int(cx.net.fiber_offsets.shape[0] - 1)      # cortex fibers only (assembly/growth must skip MT arms)
    seg = float(net.seg_rest.mean()); V0 = float(ConvexHull(net.pos[:Nc]).volume); vmin = VMIN_FRAC * V0
    k_plane = 1.0e3                                             # substrate excluded-volume stiffness [pN/µm]
    n_fib_nodes = int(net.fiber_offsets[-1])                   # Nc (MT off) or Nc+Nmt (MT on); MTOC+nucleus are beads
    gammas = physical_node_gammas(net, n_fib_nodes, N - n_fib_nodes)   # per-node drag from η (fiber log-drag + Stokes beads)
    if bulk_drag:
        # DIAGNOSIS (2026-07-09): the crawl speed is grid-dependent — the NF2007 single-fiber log-drag gives Σγ ∝ Nc,
        # so at native Nc the COM drag is ~100-300× the physical whole-cell Stokes drag 6πηR → v ∝ 1/Nc (native crawls
        # ~0 while the SAME physics crawls physiologically at coarse Nc; mac diag 60→13 nm/s as fil 120→500). The
        # physically-correct drag is Σγ_cortex = 6πηR (η=65.9 Pa·s, the AFM-validated bulk-η). ⚠️ BUT setting it here
        # DESTABILIZES the implicit solver: γ/dt (≈0.7 at native) ≪ K (crosslink stiffness ~1e6) leaves the cortex
        # rigid-body modes unregularized → NaN. So grid-consistent crawl drag is an OPEN solver-side item (regularize
        # rigid modes or add inertia), NOT a one-liner. Flag kept to document the diagnosis; do not use for production.
        gammas[:Nc] = 6.0 * np.pi * ETA_CYTOPLASM * S["R"] / Nc
    # hub crosslinks (MTOC↔arm-base, finite rest seg_um) appended to cortex crosslinks → same link force + same K
    xl_i = np.concatenate([cx.xl_i, S["hub_i"]]).astype(np.int64)
    xl_j = np.concatenate([cx.xl_j, S["hub_j"]]).astype(np.int64)
    xl_k_all = np.concatenate([cx.xl_k, S["hub_k"]]).astype(np.float64)
    xl_rest_all = np.concatenate([cx.xl_rest, S["hub_rest"]]).astype(np.float64)
    n_xl = int(xl_i.size)                                       # cortex + hub bonds (link force + K assembly)
    n_xl_cortex = int(cx.xl_i.size)                             # cortex-only (crosslink turnover EXCLUDES the stiff hub)
    kmax = max(float(net.kappa.max()) / seg**3, float(xl_k_all.max()),
               S["nuc"].k_chrom + S["nuc"].k_lamin, cp.k_int, k_plane)
    if dt is None:
        dt = crawl_cfl_dt(gammas, kmax, safety=safety)         # explicit CFL-stable (stiff crosslink sets it)
    if implicit:
        dt = dt_impl                                           # NF2007 implicit: dt bounded by accuracy, not the CFL
    gamma_rep = float(np.median(gammas[:Nc]))                  # cortex-representative scalar drag (implicit solver)
    bend_triples_np = np.ascontiguousarray(net.bend_triples, np.int64)
    xl_ij_np = np.ascontiguousarray(np.stack([xl_i, xl_j], 1), np.int64)
    kxl_np = np.ascontiguousarray(xl_k_all, np.float64)
    mem = S["mem"]; z_sub = S["z_sub"]; adh_h = 0.4; phat = S["phat"]
    front_cos_R = S["front_frac"] * S["R"]
    ph = wp.vec3d(float(phat[0]), float(phat[1]), float(phat[2])); kT = float(cp.kT)
    _tri = ConvexHull(net.pos[:Nc]).simplices.copy()           # fixed CORTEX surface triangulation (per-step volume)
    _v = net.pos[:Nc]; _cen = _v.mean(0)                       # orient every face OUTWARD (scipy hull isn't consistent)
    _n = np.cross(_v[_tri[:, 1]] - _v[_tri[:, 0]], _v[_tri[:, 2]] - _v[_tri[:, 0]])
    _flip = (_n * (_v[_tri[:, 0]] - _cen)).sum(1) < 0
    _tri[_flip] = _tri[_flip][:, [0, 2, 1]]                    # swap last two verts → consistent outward normal
    faces = np.ascontiguousarray(_tri, np.int32)
    # net sedimentation body force (gravity − buoyancy), DERIVED not tuned: −Δρ·g·v_node, Δρ=55 kg/m³ (≈1 pN/cell)
    fz_node = -(delta_rho * 9.80665 * (V0 * 1e-18) / Nc) * 1e12 if gravity else 0.0

    pos_d = wp.array(S["pos_all"].copy(), dtype=wp.vec3d, device=d); f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    total_d = wp.zeros(1, dtype=wp.float64, device=d)          # leading-edge push total (for the retrograde reaction)
    spread_total_d = wp.zeros(2, dtype=wp.float64, device=d)   # spreading push (x,y) total for the retrograde reaction
    vol_d = wp.zeros(1, dtype=wp.float64, device=d)            # per-step enclosed volume (no refresh lag)
    csum_d = wp.zeros(3, dtype=wp.float64, device=d)           # device centroid reduction (GPU-only, no host ConvexHull)
    rsum_d = wp.zeros(1, dtype=wp.float64, device=d)           # device mean-radius reduction
    faces_d = wp.array(faces, dtype=wp.int32, ndim=2, device=d)
    gamma_d = wp.array(gammas, dtype=wp.float64, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_ij = np.ascontiguousarray(np.stack([xl_i, xl_j], 1), np.int32)
    xl_d = wp.array(xl_ij, dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(xl_k_all, dtype=wp.float64, device=d)
    r0_d = wp.array(np.ascontiguousarray(xl_rest_all, np.float64), dtype=wp.float64, device=d)
    has_myo = cx.myo_i.size > 0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cx.myo_i, cx.myo_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    # Myosin-LINEAR (PI-ratified): explicit Stam-Hocky minifilament, DERIVED stall + linear force-velocity,
    # replacing the swept constant f_myo. Quasi-static ⇒ v_slide≈0 ⇒ full stall prestress (the γ source).
    myo_lin = resolve_myosin() if myosin_linear else None
    if has_myo and myosin_linear:
        vslide_d = wp.zeros(cx.myo_i.size, dtype=wp.float64, device=d)   # inter-anchor sliding rate (0 in quasi-static)
        f_myo = myo_lin.f_stall_pn                                       # 60 pN derived stall (KB-3.18), not swept
    ac_d = wp.array(S["basal"], dtype=wp.int32, device=d)
    anch_d = wp.array(S["anchors"], dtype=wp.vec3d, device=d)
    bd_d = wp.array(np.ones(S["basal"].size, np.int32), dtype=wp.int32, device=d)
    slip_d = wp.zeros(S["basal"].size, dtype=wp.float64, device=d)   # per-clutch retrograde slip (molecular-clutch traction; anchor FIXED)
    # COMPLIANT SUBSTRATE (PI experimental axis): the FA anchors become movable Winkler-spring DOFs (k_sub∝E_sub)
    # instead of fixed pins → traction becomes E-dependent (Bangasser-Odde). E_sub=0 keeps the rigid-pin path.
    substrate = resolve_substrate(E_pa=substrate_E) if substrate_E > 0 else None
    if substrate is not None:
        anch_rest_d = wp.array(np.ascontiguousarray(S["anchors"], np.float64), dtype=wp.vec3d, device=d)  # fixed dish points
    # FA MATURATION (KB-2.x): per-clutch talin unfolding → vinculin → force-gated FA growth/disassembly (mechanosensor)
    mp = MaturationParams() if fa_maturation else None
    if fa_maturation:
        Mb = S["basal"].size
        pu_d = wp.zeros(Mb, dtype=wp.float64, device=d); nv_d = wp.zeros(Mb, dtype=wp.float64, device=d)
        area_d = wp.array(np.ones(Mb), dtype=wp.float64, device=d); load_d = wp.zeros(Mb, dtype=wp.float64, device=d)
    M = S["basal"].size; nuc = S["nuc"]
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    cx_c, cy_c = float(S["c"][0]), float(S["c"][1])            # cell xy-centre (for radial spreading)
    h_basal = 0.6                                              # basal-cap height for the spreading push [µm]
    dP_area = dP_mem_area = 0.0
    frac_xl = 1.0 - np.exp(-koff_xl * dt * xl_turn_every)      # crosslink turnover fraction per turnover tick
    frac_asm = 1.0 - np.exp(-k_assembly * dt * assembly_every) if assembly else 0.0   # actin-assembly area growth
    # Per-filament BARBED-END growth (KB-3.6/Pollard ratchet): tip segment elongates at v0·dt·e^{−load·δ/kT}.
    # RATE-LIMIT the per-tick increment to ≤0.1·ℓ₀: at large implicit dt the raw v0·dt (≈0.6 µm) would jump many
    # tips to the cap in ONE step and shock the hard volume constraint (apex over-corrects → collapse). Capping the
    # increment approaches seg_max GRADUALLY over several ticks — same bounded total growth, numerically stable
    # ("growth CFL", not a physical rate change). Load-gating (ratchet) still stalls the loaded (interior) tips.
    _seg0 = CortexParams().seg_um
    v0_dt = min(poly.v0_um_s * dt * assembly_every, 0.1 * _seg0)   # ≤0.05 µm elongation per growth tick
    seg_max = 2.0 * _seg0                                      # tip-segment cap [µm]; sustained growth past 2ℓ₀ ⇒ bead insertion (staged)
    seg_min = 0.3 * _seg0                                      # pointed-end depoly floor [µm]; never shrink a segment below this (no inversion)
    grown_d = wp.zeros(1, dtype=wp.float64, device=d)          # Σ Δlength this tick (G-actin-pool budget diagnostic, KB-3.21)
    shrunk_d = wp.zeros(1, dtype=wp.float64, device=d)         # Σ Δlength REMOVED at rear pointed ends (treadmill mass balance vs grown_d, gate G5)

    # ---- S6: ECM Mikado collagen CO-SIMULATION (2026-07-09) — the cell's basal clutches grip collagen fibers and
    #      the myosin-contracted cortex pulls them → the compliant matrix REMODELS (MCF7's real mechanobiology on
    #      collagen-I). Implicit cell (large dt) + explicit-substepped collagen (stiff → small CFL); clutch_ecm_spring
    #      is the two-sided coupling. Replaces the fixed-dish clutch when --ecm. ----
    ecm_on = S.get("ecm") is not None
    if ecm_on:
        from ffn_sim.ff.fa_ecm import clutch_ecm_spring_kernel, mask_ecm_by_bound_kernel
        from ffn_sim.scripts.ff_ecm_remodel_demo import _segment_pairs, K_SEG   # _per_triple_alpha already imported at module level
        _emk = S["ecm"]; _enet = _emk.net
        _Ep0 = np.ascontiguousarray(_enet.pos, np.float64); _En = _Ep0.shape[0]
        _Etri = np.ascontiguousarray(_enet.bend_triples, np.int32); _EnT = _Etri.shape[0]
        _Eal = np.ascontiguousarray(_per_triple_alpha(_enet), np.float64)
        _Exl = np.ascontiguousarray(np.stack([_emk.xl_i, _emk.xl_j], 1), np.int32) if _emk.xl_i.size else np.zeros((0, 2), np.int32)
        _Eseg = _segment_pairs(_enet.fiber_offsets)
        _Esegr = np.linalg.norm(_Ep0[_Eseg[:, 0]] - _Ep0[_Eseg[:, 1]], axis=1) if _Eseg.shape[0] else np.zeros(0)
        _Egam = np.where(_emk.pinned, 1.0e18, 1.0).astype(np.float64)      # pinned boundary immovable; bulk drag 1.0 (demo units)
        _Ekmax = max(float(_enet.kappa.max()) / (0.5 ** 3) if _enet.kappa.size else 1.0, K_SEG,
                     float(_emk.xl_k.max()) if _emk.xl_k.size else 1.0, float(cp.k_int))
        _Edt = 0.1 * 1.0 / _Ekmax                                          # collagen explicit CFL-stable substep [s]
        _Ensub = 20                                                        # collagen substeps per cell step (partial relax; remodel accumulates over steps)
        Ep_d = wp.array(_Ep0.copy(), dtype=wp.vec3d, device=d); Ef_d = wp.zeros(_En, dtype=wp.vec3d, device=d)
        Etri_d = wp.array(_Etri, dtype=wp.int32, ndim=2, device=d); Eal_d = wp.array(_Eal, dtype=wp.float64, device=d)
        Exl_d = wp.array(_Exl, dtype=wp.int32, ndim=2, device=d)
        Exlk_d = wp.array(np.ascontiguousarray(_emk.xl_k, np.float64), dtype=wp.float64, device=d)
        Exlr_d = wp.array(np.ascontiguousarray(_emk.xl_rest, np.float64), dtype=wp.float64, device=d)
        Eseg_d = wp.array(np.ascontiguousarray(_Eseg, np.int32), dtype=wp.int32, ndim=2, device=d)
        Esegk_d = wp.array(np.full(_Eseg.shape[0], K_SEG), dtype=wp.float64, device=d)
        Esegr_d = wp.array(np.ascontiguousarray(_Esegr, np.float64), dtype=wp.float64, device=d)
        Egam_d = wp.array(_Egam, dtype=wp.float64, device=d)
        en_base_d = wp.array(np.ascontiguousarray(S["ecm_node"], np.int32), dtype=wp.int32, device=d)   # clutch → collagen node (frame-0 attach)
        en_d = wp.zeros(S["basal"].size, dtype=wp.int32, device=d)         # per-step EFFECTIVE attach (masked by bound state below)
        Edummy = wp.zeros(N, dtype=wp.vec3d, device=d)                     # discarded cell-side out during collagen substeps
        Edummy_e = wp.zeros(_En, dtype=wp.vec3d, device=d)                 # discarded collagen-side out during the cell force-eval
        Ep0_host = _Ep0                                                    # frame-0 collagen (remodel reference)

    vs = {"dP": 0.0, "g": np.zeros((Nc, 3))}                   # osmotic ΔP + exact volume gradient (set by the loop)

    def full_force(x_np):
        """Full FF force (elastic + active) at x — the implicit solver's RHS. The cortex osmotic/turgor force is
        the EXACT volume-gradient force ``ΔP·g`` (added host-side from ``vs``), so its stiffness is the clean
        rank-1 ``k_vol·g·gᵀ`` the implicit solve treats implicitly; the membrane inward tension stays a kernel."""
        pos_d.assign(np.ascontiguousarray(x_np, np.float64).reshape(N, 3))
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            if myosin_linear:
                wp.launch(minifilament_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, vslide_d, wp.float64(f_myo), wp.float64(myo_lin.v0_um_s), f_d], device=d)
            else:
                wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if fz_node != 0.0:
            wp.launch(gravity_kernel, dim=Nc, inputs=[wp.float64(fz_node), f_d], device=d)
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        if clutches:
            if ecm_on:                                         # S6: collagen-I IS the substrate → the clutch pulls the LIVE fiber node (two-way; the cell FEELS the matrix)
                wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),
                          wp.float64(cp.rest_um), f_d, Edummy_e], device=d)                  # +f traction on the cell; the −f reaction is applied in the collagen substep (partitioned/staggered)
            else:
                wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                          wp.float64(cp.rest_um), f_d], device=d)
                if flow:                                       # molecular-clutch retrograde-flow traction (anchor FIXED)
                    wp.launch(clutch_slip_traction_kernel, dim=M, inputs=[f_d, ac_d, bd_d, slip_d, ph, wp.float64(cp.k_int)], device=d)
        if spread:
            spread_total_d.zero_()
            wp.launch(spreading_push_kernel, dim=Nc, inputs=[pos_d, wp.float64(cx_c), wp.float64(cy_c),
                      wp.float64(z_sub), wp.float64(h_basal), wp.float64(0.1), wp.float64(S["f_pro"]),
                      wp.float64(poly.delta_um), wp.float64(kT), f_d, spread_total_d], device=d)
            wp.launch(spreading_reaction_kernel, dim=Nc, inputs=[spread_total_d, wp.float64(Nc), f_d], device=d)
        # NOTE (2026-07-09): the leading-edge PROTRUSION is no longer a body force here — it is real DIRECTED
        # front barbed-end polymerization applied at the growth tick (directed_front_growth_kernel, post-step),
        # advanced by reshape. The retired leading_edge_push/protrusion_reaction body-force pair tore the cell.
        wp.synchronize_device(d)
        F = f_d.numpy().reshape(-1)
        # EXACT cortex osmotic/turgor force f = ΔP·g, computed from the CURRENT x (consistent across Newton iters)
        pcx = np.ascontiguousarray(x_np, np.float64).reshape(N, 3)[:Nc]; cen = pcx.mean(0)
        aa = pcx[faces[:, 0]] - cen; bb = pcx[faces[:, 1]] - cen; ccf = pcx[faces[:, 2]] - cen
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        dPv = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dPv = min(max(dPv, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        F[:3 * Nc] += (dPv * volume_gradient(pcx, faces, cen)).reshape(-1)
        return F

    # ---- GPU-RESIDENT implicit path (cupy) — the GPU-only fast solver for native runs ----
    gpu_impl = implicit and str(device).startswith("cuda")
    if gpu_impl:
        import cupy as cpx
        from ffn_sim.ff.implicit_ff import ff_implicit_step_gpu
        bt_cp = cpx.asarray(np.ascontiguousarray(net.bend_triples, np.int64))
        al_cp = cpx.asarray(alpha); xlij_cp = cpx.asarray(xl_ij_np); kxl_cp = cpx.asarray(kxl_np)
        faces_cp = cpx.asarray(faces.astype(np.int64))
        bas_cp = cpx.asarray(S["basal"].astype(np.int64))     # basal clutch/contact node indices (for diag_extra)
        # MT-node drag correction (R-γ / R7): the solver's diagonal is a scalar gamma_rep/dt, so the added MT arms +
        # MTOC (nodes [Nc,Ne)) would get the wrong overdamped timescale. Fix ONLY those nodes via the K diagonal:
        # diag[3i(+1,+2)] += (γ_i − γ_rep)/dt (constant → precomputed). Cortex + nucleus are untouched ⇒ their
        # bit-identity holds and the pre-existing scalar-γ nucleus mismatch is left for a separate PI decision.
        _diag_mt = np.zeros(3 * N)
        if Ne > Nc:
            _mt = np.arange(Nc, Ne)
            _c = (gammas[_mt] - gamma_rep) / dt
            _diag_mt[3 * _mt] = _c; _diag_mt[3 * _mt + 1] = _c; _diag_mt[3 * _mt + 2] = _c
        diag_base_cp = cpx.asarray(_diag_mt)                   # zeros when MT OFF ⇒ diag == the pre-change cpx.zeros(3N)

        def gpu_force_fn(x_cp):
            """Full FF force at x as a device-resident cupy array (no host round-trip). Warp force kernels write
            f_d; the exact osmotic force ΔP·g is added in cupy."""
            pos_d.assign(wp.array(cpx.ascontiguousarray(x_cp.reshape(N, 3)), dtype=wp.vec3d, device=d))
            wp.launch(_zero, dim=N, inputs=[f_d], device=d)
            wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
            wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
            if has_myo:
                if myosin_linear:
                    wp.launch(minifilament_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, vslide_d, wp.float64(f_myo), wp.float64(myo_lin.v0_um_s), f_d], device=d)
                else:
                    wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
            wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
            if fz_node != 0.0:
                wp.launch(gravity_kernel, dim=Nc, inputs=[wp.float64(fz_node), f_d], device=d)
            if n_nuc:
                wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), centre,
                          wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                          wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
            wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
            if clutches:
                if ecm_on:                                       # S6: collagen-I IS the substrate → the clutch pulls the LIVE fiber node (two-way; native production path)
                    wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),
                              wp.float64(cp.rest_um), f_d, Edummy_e], device=d)              # +f traction on the cell; −f reaction applied in the collagen substep
                else:
                    wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                              wp.float64(cp.rest_um), f_d], device=d)
                    if flow:                                     # molecular-clutch retrograde-flow traction (anchor FIXED)
                        wp.launch(clutch_slip_traction_kernel, dim=M, inputs=[f_d, ac_d, bd_d, slip_d, ph, wp.float64(cp.k_int)], device=d)
            if spread:
                spread_total_d.zero_()
                wp.launch(spreading_push_kernel, dim=Nc, inputs=[pos_d, wp.float64(cx_c), wp.float64(cy_c),
                          wp.float64(z_sub), wp.float64(h_basal), wp.float64(0.1), wp.float64(S["f_pro"]),
                          wp.float64(poly.delta_um), wp.float64(kT), f_d, spread_total_d], device=d)
                wp.launch(spreading_reaction_kernel, dim=Nc, inputs=[spread_total_d, wp.float64(Nc), f_d], device=d)
            # protrusion is directed front polymerization (post-step tick), NOT a body force here — see full_force note
            wp.synchronize_device(d)                              # mirrors full_force / the explicit loop
            F = cpx.asarray(f_d).reshape(-1).copy()
            # exact osmotic force ΔP·g in cupy (g = ∂V/∂x from the fixed face triangulation)
            pc = x_cp.reshape(N, 3)[:Nc]; ce = pc.mean(0)
            a3 = pc[faces_cp[:, 0]] - ce; b3 = pc[faces_cp[:, 1]] - ce; c3 = pc[faces_cp[:, 2]] - ce
            Vg = float(cpx.abs((a3 * cpx.cross(b3, c3)).sum() / 6.0))
            g = cpx.zeros((Nc, 3))
            cpx.add.at(g, faces_cp[:, 0], cpx.cross(b3, c3) / 6.0)
            cpx.add.at(g, faces_cp[:, 1], cpx.cross(c3, a3) / 6.0)
            cpx.add.at(g, faces_cp[:, 2], cpx.cross(a3, b3) / 6.0)
            dPg = TURGOR_PI_IN0 * (V0 - vmin) / max(Vg - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
            dPg = min(max(dPg, -TURGOR_PI_IN0), TURGOR_PI_IN0)
            F[:3 * Nc] += (dPg * g).reshape(-1)
            return F

    frames, com_traj, times, vol_traj, diverged = [], [], [], [], False
    # FEM-field (per recorded frame) + FA-junction recording (cortex stress topology = cortex xl + myosin only)
    _area0 = face_areas(net.pos[:Nc], faces)                  # rest cortex face areas (frame-0)
    _xl_st = np.stack([cx.xl_i, cx.xl_j], 1).astype(np.int64)
    _kxl_st = np.ascontiguousarray(cx.xl_k, np.float64); _r0_st = np.ascontiguousarray(cx.xl_rest, np.float64)
    _myo_st = np.stack([cx.myo_i, cx.myo_j], 1).astype(np.int64) if cx.myo_i.size else np.zeros((0, 2), np.int64)
    svm_frames, strain_frames, bound_frames, anch_frames = [], [], [], []
    ecm_frames = []                                            # S6: collagen positions over time (the remodel/recruitment PROCESS)
    t0 = time.time()
    for step in range(steps):
        if ecm_on and clutches:                                # S6: refresh en_d (bound-masked clutch→collagen map) at the TOP so the
            wp.launch(mask_ecm_by_bound_kernel, dim=M, inputs=[en_base_d, bd_d, en_d], device=d)   # two-way force sites read the CURRENT bound state (not the wp.zeros init = node-0)
        if step % refresh_every == 0:                          # slow modes via DEVICE reductions (GPU-only, no host ConvexHull)
            csum_d.zero_(); wp.launch(sum_pos_kernel, dim=Nc, inputs=[pos_d, csum_d], device=d)
            cc = csum_d.numpy() / Nc                            # centroid — 3 scalars cross the bus
            centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2])); cx_c, cy_c = float(cc[0]), float(cc[1])
            rsum_d.zero_(); wp.launch(sum_radius_kernel, dim=Nc, inputs=[pos_d, centre, rsum_d], device=d)
            Rm = float(rsum_d.numpy()[0]) / Nc                  # mean radius — 1 scalar; area ≈ 4πR² (no ConvexHull)
            area = 4.0 * np.pi * Rm**2
            dP_mem_area = -2.0 * min(mem.gamma_mem, mem.tau_lysis) / max(Rm, 1e-9) * area / Nc
            if not np.isfinite(Rm) or Rm > 50.0 * S["R"]:
                print(f"  [!] divergence at step {step} — truncating"); diverged = True; break
        # PER-STEP enclosed volume (on-device, NO refresh lag) → BIDIRECTIONAL osmotic ΔP (incompressible cytoplasm:
        # outward when V<V0, INWARD when V>V0). No lag ⇒ the stiff osmotic constraint is stable (no balloon, no collapse).
        vol_d.zero_()
        wp.launch(cortex_volume_kernel, dim=faces.shape[0], inputs=[pos_d, faces_d, centre, vol_d], device=d)
        vol = abs(float(vol_d.numpy()[0]))
        dP = TURGOR_PI_IN0 * (V0 - vmin) / max(vol - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dP = min(max(dP, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        dP_area = dP * area / Nc
        if flow and clutches:                                  # RETROGRADE FLOW (anchor FIXED): accumulate per-clutch slip ONCE
            wp.launch(clutch_slip_accumulate_kernel, dim=M,    # per step (the traction force k·slip·phat is applied per force-eval
                      inputs=[bd_d, slip_d, wp.float64(v_retro_um_s * dt)], device=d)        # with clutch_spring, below)
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            if myosin_linear:
                wp.launch(minifilament_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, vslide_d, wp.float64(f_myo), wp.float64(myo_lin.v0_um_s), f_d], device=d)
            else:
                wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_area), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if fz_node != 0.0:                                     # gravity − buoyancy (real body force, ≈1 pN/cell)
            wp.launch(gravity_kernel, dim=Nc, inputs=[wp.float64(fz_node), f_d], device=d)
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        if clutches:
            if ecm_on:                                         # S6: collagen-I IS the substrate → the clutch pulls the LIVE fiber node (two-way; also the tip-load the ratchet reads)
                wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),
                          wp.float64(cp.rest_um), f_d, Edummy_e], device=d)                  # +f traction on the cell; −f reaction applied in the collagen substep
            else:
                wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                          wp.float64(cp.rest_um), f_d], device=d)
                if flow:                                       # molecular-clutch retrograde-flow traction (anchor FIXED): forward k·slip·phat
                    wp.launch(clutch_slip_traction_kernel, dim=M, inputs=[f_d, ac_d, bd_d, slip_d, ph, wp.float64(cp.k_int)], device=d)
        if spread:                                             # EMERGENT spreading: peripheral radial-outward ratchet
            spread_total_d.zero_()
            wp.launch(spreading_push_kernel, dim=Nc, inputs=[pos_d, wp.float64(cx_c), wp.float64(cy_c),
                      wp.float64(z_sub), wp.float64(h_basal), wp.float64(0.1), wp.float64(S["f_pro"]),
                      wp.float64(poly.delta_um), wp.float64(kT), f_d, spread_total_d], device=d)
            wp.launch(spreading_reaction_kernel, dim=Nc, inputs=[spread_total_d, wp.float64(Nc), f_d], device=d)
        # protrusion is directed front polymerization (post-step tick), NOT a body force here — see full_force note
        if gpu_impl:                                           # GPU-RESIDENT implicit (cupy) — GPU-only, native-scale
            x_cp = cpx.asarray(pos_d).reshape(-1).copy()
            pc = x_cp.reshape(N, 3)[:Nc]; ce = pc.mean(0)
            a3 = pc[faces_cp[:, 0]] - ce; b3 = pc[faces_cp[:, 1]] - ce; c3 = pc[faces_cp[:, 2]] - ce
            Vc = float(cpx.abs((a3 * cpx.cross(b3, c3)).sum() / 6.0))
            if not np.isfinite(Vc) or Vc <= 1.02 * vmin:       # volume collapsed toward vmin (implicit overshoot at too-large dt)
                print(f"  [!] volume collapse at step {step} (Vc={Vc:.1f} ≤ vmin={vmin:.1f}) — truncating"); diverged = True; break
            g = cpx.zeros((Nc, 3))
            cpx.add.at(g, faces_cp[:, 0], cpx.cross(b3, c3) / 6.0)
            cpx.add.at(g, faces_cp[:, 1], cpx.cross(c3, a3) / 6.0)
            cpx.add.at(g, faces_cp[:, 2], cpx.cross(a3, b3) / 6.0)
            k_vol = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-3 * V0) ** 2   # floor 1e-3·V0 (not 1e-9→overflow when squared)
            vol_g = cpx.zeros(3 * N); vol_g[:3 * Nc] = g.reshape(-1)
            # CLUTCH + SUBSTRATE stiffness → implicit K diagonal (so dt is not capped by their explicit CFL);
            # start from the constant MT-γ correction (zeros when MT OFF ⇒ bit-identical to the pre-change zeros)
            diag = diag_base_cp.copy()
            if clutches:
                bset = bas_cp[cpx.asarray(bd_d) > 0]           # bound basal actin nodes → k_int·I (spring-to-anchor ≈ pin)
                diag[3 * bset] += cp.k_int; diag[3 * bset + 1] += cp.k_int; diag[3 * bset + 2] += cp.k_int
            diag[3 * bas_cp + 2] += k_plane                    # substrate excluded-volume on basal z
            x_cp = ff_implicit_step_gpu(x_cp, gpu_force_fn, bt_cp, al_cp, xlij_cp, kxl_cp,
                                        gamma=gamma_rep, dt=dt, vol_g=vol_g, k_vol=k_vol, diag_extra=diag,
                                        com_gamma=(6.0 * np.pi * ETA_CYTOPLASM * S["R"]) if com_drag else None,
                                        n_cortex=Nc)
            pos_d.assign(wp.array(cpx.ascontiguousarray(x_cp.reshape(N, 3)), dtype=wp.vec3d, device=d))
        elif implicit:                                         # host implicit (CPU dev fallback)
            xv = pos_d.numpy().reshape(-1)
            pcx2 = xv.reshape(N, 3)[:Nc]; cen2 = pcx2.mean(0)
            aa = pcx2[faces[:, 0]] - cen2; bb = pcx2[faces[:, 1]] - cen2; ccf = pcx2[faces[:, 2]] - cen2
            Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
            gN = volume_gradient(pcx2, faces, cen2)            # ∂V/∂x (Nc,3) — exact osmotic force direction
            dP = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
            vs["dP"] = min(max(dP, -TURGOR_PI_IN0), TURGOR_PI_IN0); vs["g"] = gN
            k_vol = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-3 * V0) ** 2    # floor 1e-3·V0 (osmotic stiffness; not 1e-9→overflow)
            vol_g = np.zeros(3 * N); vol_g[:3 * Nc] = gN.reshape(-1)
            _de = _cg = None
            if com_drag:                                       # MODAL drag: rigid-COM feels the physical 6πηR; the bound
                _cg = 6.0 * np.pi * ETA_CYTOPLASM * S["R"]      # clutches (added to K) regularise the now-soft rigid mode
                _bd = bd_d.numpy(); _bb = S["basal"][_bd > 0]
                _de = np.zeros(3 * N)
                _de[3 * _bb] = cp.k_int; _de[3 * _bb + 1] = cp.k_int; _de[3 * _bb + 2] = cp.k_int
            xv, _info = implicit_step_current(xv, full_force, bend_triples_np, alpha, xl_ij_np, kxl_np,
                                              gamma=gamma_rep, dt=dt, n_newton=1, vol_g=vol_g, k_vol=k_vol,
                                              diag_extra=_de, com_gamma=_cg, n_cortex=Nc)
            # HARD incompressibility: project the cortex to V=V0 exactly (Newton on the volume constraint along
            # g=∂V/∂x). The basal cap is held by the clutches, so the inward correction drops the free APEX ⇒ as
            # the base area grows (assembly) the cell FLATTENS at constant volume instead of inflating.
            xr = xv.reshape(N, 3); pcx3 = xr[:Nc].copy(); cen3 = pcx3.mean(0)
            for _ in range(3):
                a3 = pcx3[faces[:, 0]] - cen3; b3 = pcx3[faces[:, 1]] - cen3; c3 = pcx3[faces[:, 2]] - cen3
                Vp = abs(float((a3 * np.cross(b3, c3)).sum() / 6.0))
                gp = volume_gradient(pcx3, faces, cen3); den = float((gp * gp).sum())
                if den > 1e-9:
                    pcx3 = pcx3 - ((Vp - V0) / den) * gp
            xr[:Nc] = pcx3; xv = xr.reshape(-1)
            pos_d.assign(np.ascontiguousarray(xv, np.float64).reshape(N, 3))
        else:                                                  # explicit physical-γ step (CFL-bound)
            wp.launch(axpy_physical_kernel, dim=N, inputs=[pos_d, wp.float64(dt), gamma_d, f_d], device=d)
        if substrate is not None and clutches:                 # COMPLIANT SUBSTRATE: movable FA anchors at the clutch↔substrate series equilibrium (stable)
            wp.launch(substrate_anchor_equilibrium_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, anch_rest_d, bd_d,
                      wp.float64(cp.k_int), wp.float64(substrate.k_sub)], device=d)
        if step % reshape_every == 0:                          # inextensibility (NF2007 §5.3)
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        if step % xl_turn_every == 0 and step > 0:             # crosslink turnover (on-device Maxwell relax; cortex xl only, NOT the hub)
            wp.launch(xl_turnover_kernel, dim=n_xl_cortex, inputs=[pos_d, xl_d, r0_d, wp.float64(frac_xl)], device=d)
        if assembly and step % assembly_every == 0 and step > 0:   # dynamic cortex AREA GROWTH (actin assembly, tension-gated)
            wp.launch(actin_assembly_kernel, dim=n_cortex_fib, inputs=[pos_d, foff_d, soff_d, sr_d, wp.float64(frac_asm)], device=d)
        if growth and step % assembly_every == 0 and step > 0:   # per-filament BARBED-END polymerization (KB-3.6 ratchet; reads tip load in f_d)
            wp.launch(barbed_end_growth_kernel, dim=n_cortex_fib, inputs=[pos_d, foff_d, soff_d, sr_d,
                      wp.float64(v0_dt), wp.float64(poly.delta_um), wp.float64(kT), f_d, wp.float64(seg_max), grown_d], device=d)
        if flow and cortex_treadmill and step % assembly_every == 0 and step > 0:  # PER-FIBER ACTIN TREADMILL (COM-conserving): each forward
            wp.launch(fiber_treadmill_kernel, dim=n_cortex_fib,    # fiber grows its barbed end = shrinks its pointed end → material
                      inputs=[pos_d, ph, wp.float64(0.0), foff_d, soff_d, sr_d, wp.float64(v0_dt),   # flows barbed→pointed (retrograde),
                      wp.float64(poly.delta_um), wp.float64(kT), f_d, wp.float64(seg_max),           # fiber length+COG conserved ⇒ ONLY
                      wp.float64(seg_min), grown_d], device=d)                                       # the clutch traction translocates (G4)
        if protrude and not flow and step % assembly_every == 0 and step > 0:  # (non-treadmill) DIRECTED FRONT polymerization only
            _pcm = pos_d.numpy()[:Nc].mean(0)
            _cmv = wp.vec3d(float(_pcm[0]), float(_pcm[1]), float(_pcm[2]))
            wp.launch(directed_front_growth_kernel, dim=n_cortex_fib, inputs=[pos_d, _cmv, ph,
                      wp.float64(front_cos_R), wp.float64(0.0), foff_d, soff_d, sr_d,
                      wp.float64(v0_dt), wp.float64(poly.delta_um), wp.float64(kT), f_d,
                      wp.float64(seg_max), grown_d], device=d)
        if rear_depoly and not flow and step % assembly_every == 0 and step > 0:  # (non-treadmill) rear depoly only
            _pcm = pos_d.numpy()[:Nc].mean(0)
            _cmv = wp.vec3d(float(_pcm[0]), float(_pcm[1]), float(_pcm[2]))
            wp.launch(pointed_end_depoly_kernel, dim=n_cortex_fib, inputs=[pos_d, _cmv, ph, wp.float64(front_cos_R),
                      foff_d, soff_d, sr_d, wp.float64(v0_dt), wp.float64(seg_min), shrunk_d], device=d)
        if clutches and rupture and step % kmc_every == 0 and step > 0:   # catch-slip turnover + nascent-adhesion rebind
            if ecm_on:                                         # S6: seat the catch-slip anchor at the LIVE collagen node so rupture reads the ECM-clutch
                _eh = Ep_d.numpy(); _enb = np.asarray(S["ecm_node"]); _bh = bd_d.numpy()   # extension (actin↔collagen), NOT a stale dish pin
                _sel = (_enb >= 0) & (_bh > 0)
                if _sel.any():
                    _ah = anch_d.numpy(); _ah[_sel] = _eh[_enb[_sel]]; anch_d.assign(_ah)
            wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                      wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                      wp.float64(dt * kmc_every), wp.int32(step)], device=d)
            p = pos_d.numpy(); bd = bd_d.numpy(); anchors = anch_d.numpy(); zb = p[S["basal"], 2]
            can = (bd == 0) & (zb < z_sub + adh_h)
            p_on = 1.0 - np.exp(-dt * kmc_every * cp.k_on)
            if treadmill:
                # C2 EXPLICIT FRONT-BIAS — REJECTED DEAD-END (2026-07-09). Idea: nascent adhesions re-form only in the
                # LEADING half (s>0), rear releases. BUT the directed crawl already EMERGES from protrusion + uniform
                # clutch turnover (below) at physiological ~60 nm/s, traction-driven; and this explicit bias FAILS —
                # as the COM advances, formerly-front clutches fall BEHIND it (s→<0) and stop re-forming, so the bound
                # population runs down to ~0 → traction collapses → LESS crawl (mac diag: 46 vs 61 nm/s, bound 0.00).
                # Left in, off by default, as a documented negative result. The working C2 is the emergent path.
                com = p[:Nc].mean(0)
                s = (p[S["basal"]] - com) @ phat                  # signed leading(+)/trailing(−) position [µm]
                w = (s > 0.0).astype(np.float64)
                reb = can & (np.random.default_rng(step).random(M) < p_on * w)
            else:
                reb = can & (np.random.default_rng(step).random(M) < p_on)
            if reb.any():
                if ecm_on and ecm_regrip:                      # S6 MIGRATION: the rebinding clutch RE-GRIPS the nearest CURRENT collagen
                    _ep = Ep_d.numpy(); _ri = np.where(reb)[0]  # node under it (KB-1.23 FA captures fibres within R_FA≈1.5µm) — as the
                    _dq, _jq = cKDTree(_ep).query(p[S["basal"][reb]], k=1)   # cell protrudes forward, front clutches grab NEW fibres ahead,
                    _ok = _dq <= 1.5                            # rear clutches (no fibre in reach) stay released → the grip rolls forward
                    _enb = en_base_d.numpy()
                    _enb[_ri[_ok]] = _jq[_ok].astype(np.int32)  # re-map the clutch→collagen node (drives en_d at the force sites next step)
                    en_base_d.assign(_enb)
                    anchors[_ri[_ok]] = _ep[_jq[_ok]]           # seat the catch-slip anchor at the new node
                    bd[_ri[_ok]] = 1                            # only clutches that found collagen bind
                else:
                    anchors[reb, 0] = p[S["basal"][reb], 0]; anchors[reb, 1] = p[S["basal"][reb], 1]
                    anchors[reb, 2] = z_sub; bd[reb] = 1
                anch_d.assign(anchors); bd_d.assign(bd)
        if fa_maturation and clutches and step % kmc_every == 0 and step > 0:   # FA MATURATION mechanosensor (KB-2.x)
            wp.launch(clutch_load_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int), wp.float64(cp.rest_um), load_d], device=d)
            wp.launch(talin_vinculin_kernel, dim=M, inputs=[load_d, pu_d, nv_d, bd_d, wp.float64(mp.ku0), wp.float64(mp.dx_um),
                      wp.float64(mp.kT), wp.float64(mp.k_refold), wp.float64(mp.k_rec), wp.float64(mp.k_diss), wp.float64(mp.n_max), wp.float64(dt * kmc_every)], device=d)
            wp.launch(fa_growth_kernel, dim=M, inputs=[load_d, area_d, wp.float64(mp.kg0), wp.float64(mp.kd), wp.float64(mp.n_hill), wp.float64(mp.fth), wp.float64(dt * kmc_every)], device=d)
            wp.launch(fa_disassemble_kernel, dim=M, inputs=[area_d, bd_d, wp.float64(0.1)], device=d)   # sub-threshold FAs unbind (force-gated adhesion)
        if ecm_on:                                             # S6: SUBSTEP the collagen under the cell's clutch traction → REMODEL
            if clutches:                                       # only ENGAGED clutches grip the collagen (catch-slip release stops pulling;
                wp.launch(mask_ecm_by_bound_kernel, dim=M, inputs=[en_base_d, bd_d, en_d], device=d)   # OFF control → no grip → remodel is traction-driven)
            for _es in range(_Ensub):                          # (cell pos_d held; collagen is stiff → many small CFL substeps per cell step)
                wp.launch(_zero, dim=_En, inputs=[Ef_d], device=d)
                wp.launch(cytosim_bending_kernel, dim=_EnT, inputs=[Ep_d, Etri_d, Eal_d, Ef_d], device=d)
                if _Exl.shape[0]:
                    wp.launch(link_spring_kernel, dim=_Exl.shape[0], inputs=[Ep_d, Exl_d, Exlk_d, Exlr_d, Ef_d], device=d)
                if _Eseg.shape[0]:
                    wp.launch(link_spring_kernel, dim=_Eseg.shape[0], inputs=[Ep_d, Eseg_d, Esegk_d, Esegr_d, Ef_d], device=d)
                if clutches:
                    wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),   # BOUND cell basal
                              wp.float64(cp.rest_um), Edummy, Ef_d], device=d)                                          # clutch pulls the fiber
                wp.launch(axpy_physical_kernel, dim=_En, inputs=[Ep_d, wp.float64(_Edt), Egam_d, Ef_d], device=d)   # pinned bulk BC via huge γ
        if step % record_every == 0:
            p = pos_d.numpy(); pcxr = p[:Nc]; cc = pcxr.mean(0)      # frames need the host copy (rare)
            frames.append(p.astype(np.float32)); com_traj.append(cc.copy()); times.append(step * dt)
            if ecm_on:                                           # S6: record the collagen positions this frame (the recruitment PROCESS)
                ecm_frames.append(Ep_d.numpy().astype(np.float32))
            ar = pcxr[faces[:, 0]] - cc; br = pcxr[faces[:, 1]] - cc; cr = pcxr[faces[:, 2]] - cc
            vol_traj.append(abs(float((ar * np.cross(br, cr)).sum() / 6.0)) / V0)   # face-based volume (no ConvexHull)
            # FEM fields on the cortex (virial σ_vm + areal strain) + FA-junction snapshot (bound clutch↔anchor)
            svmf, _pf = cortex_node_stress(pcxr, np.full(Nc, V0 / Nc), xl_ij=_xl_st, k_xl=_kxl_st, r0_xl=_r0_st,
                                           myo_ij=_myo_st, f_myo=f_myo, dP=dP)
            svm_frames.append(svmf.astype(np.float32))
            strain_frames.append(cortex_node_areal_strain(pcxr, faces, _area0).astype(np.float32))
            bound_frames.append(bd_d.numpy().astype(np.int8)); anch_frames.append(anch_d.numpy().astype(np.float32))
    wp.synchronize_device(d)
    p = pos_d.numpy(); cc_final = p[:Nc].mean(0)
    frames.append(p.astype(np.float32)); com_traj.append(cc_final.copy()); times.append(step * dt)
    svmf, _pf = cortex_node_stress(p[:Nc], np.full(Nc, V0 / Nc), xl_ij=_xl_st, k_xl=_kxl_st, r0_xl=_r0_st,
                                   myo_ij=_myo_st, f_myo=f_myo, dP=dP)
    svm_frames.append(svmf.astype(np.float32))
    strain_frames.append(cortex_node_areal_strain(p[:Nc], faces, _area0).astype(np.float32))
    bound_frames.append(bd_d.numpy().astype(np.int8)); anch_frames.append(anch_d.numpy().astype(np.float32))
    xp = p
    com_traj = np.array(com_traj); times = np.array(times)
    disp_along = float((com_traj[-1] - com_traj[0]) @ phat)        # net COM displacement along crawl axis [µm]
    disp_perp = float(np.linalg.norm((com_traj[-1] - com_traj[0]) - disp_along * phat))
    T = times[-1] if times[-1] > 0 else 1.0
    v_crawl_nm_s = disp_along / T * 1e3                            # nm/s
    bd = bd_d.numpy(); anchors = anch_d.numpy()
    ecm_pos0 = Ep0_host if ecm_on else None                    # S6: frame-0 collagen (remodel reference)
    ecm_posf = Ep_d.numpy() if ecm_on else None                # S6: remodeled collagen
    if ecm_on:                                                 # S6: traction = the ENGAGED clutch extension against the LIVE collagen (Option B: collagen is the substrate)
        _enb = np.asarray(S["ecm_node"]); _be = (_enb >= 0) & (bd > 0)
        _Le = np.linalg.norm(xp[S["basal"]][_be] - ecm_posf[_enb[_be]], axis=1) if _be.any() else np.zeros(0)
        Fclutch = cp.k_int * np.maximum(_Le - cp.rest_um, 0.0)
    else:
        L = np.linalg.norm(xp[S["basal"]] - anchors, axis=1)
        Fclutch = cp.k_int * np.maximum(L - cp.rest_um, 0.0) * bd
    ecm_metrics = (ecm_remodel_metrics(ecm_pos0, ecm_posf, S["ecm_node"], S["basal"], xp,
                                       S["ecm"].net.fiber_offsets) if ecm_on else None)   # S6: coherent-remodel decomposition
    # ---- CONTACT / ADHESION diagnostics (foundation state) ----
    pcx = xp[:Nc]                                              # cortex nodes only
    zc = pcx[:, 2] - z_sub                                     # cortex node height above the substrate plane
    in_contact = zc < 0.25                                     # nodes within 0.25 µm of the substrate
    contact_r = float(np.hypot(pcx[in_contact, 0] - cc_final[0], pcx[in_contact, 1] - cc_final[1]).max()) if in_contact.any() else 0.0
    basal_gap = float(xp[S["basal"], 2].mean() - z_sub)        # mean basal-node lift-off (0 = flat contact)
    cell_h = float(zc.max())                                   # apical height (spread dome should be < R)
    ecm_bound = (np.asarray(S["ecm_node"]) >= 0) if ecm_on else None
    return dict(frames=frames, com=com_traj, times=times, vol=np.array(vol_traj), dt=dt, wall_s=time.time() - t0,
                ecm_pos0=ecm_pos0, ecm_posf=ecm_posf, ecm_node=(np.asarray(S["ecm_node"]) if ecm_on else None),
                ecm_frames=(np.array(ecm_frames) if ecm_on and ecm_frames else None), ecm_metrics=ecm_metrics,
                disp_along_um=disp_along, disp_perp_um=disp_perp, v_crawl_nm_s=v_crawl_nm_s,
                traction_nN=float(Fclutch.sum() / 1e3), bound_frac=float(bd.mean()), n_clutch=M,
                n_contact=int(in_contact.sum()), contact_radius_um=contact_r, basal_gap_um=basal_gap,
                cell_height_um=cell_h, R_um=S["R"], vol_final=float(vol_traj[-1]) if vol_traj else np.nan,
                f_pro_pN=S["f_pro"], gamma_min=float(gammas.min()), gamma_max=float(gammas.max()),
                F_star_pN=cp.F_star_pN, kmax=float(kmax),
                grown_um=float(grown_d.numpy()[0]), shrunk_um=float(shrunk_d.numpy()[0]),   # treadmill mass balance (gate G5)
                svm=np.array(svm_frames), cstrain=np.array(strain_frames), faces=faces.astype(np.int32),
                bound_frames=np.array(bound_frames), anch_frames=np.array(anch_frames), basal=S["basal"],
                k_int=float(cp.k_int), clutch_rest_um=float(cp.rest_um))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cortex-fil", type=int, default=900)
    ap.add_argument("--seed", type=int, default=7, help="cortex/nucleus RNG seed (cross-seed ensemble: is the crawl drift seed-robust, not a single-seed artifact?)")
    ap.add_argument("--bulk-drag", action="store_true", help="EXPERIMENTAL — grid-invariant whole-cell drag Σγ_cortex=6πηR (the physical "
                    "Stokes drag; the crawl speed is otherwise ∝1/Nc, native crawls ~0). ⚠️ DESTABILIZES the current implicit solver: at native "
                    "Nc the per-node γ=6πηR/Nc≈0.035 makes the solver diagonal γ/dt ≪ the crosslink stiffness K, so the cortex rigid-body modes go "
                    "unregularized → NaN. The grid-consistent crawl drag needs a solver-side fix (regularize rigid modes / add an inertial term) — "
                    "OPEN ITEM for PI, not a one-liner. Left as a flag documenting the diagnosis.")
    ap.add_argument("--com-drag", action="store_true", help="MODAL DRAG (2026-07-09, the WORKING native-crawl fix): the crawl COM must feel the "
                    "physical whole-cell Stokes drag 6πηR, but a uniform small γ makes the CG ill-conditioned (γ/dt≪K → NaN). Instead keep the LARGE "
                    "per-node γ (deformation modes stay well-conditioned) and lower ONLY the rigid-COM-translation mode's drag to 6πηR via a symmetric "
                    "rank-3 correction inside the implicit solve (regularized by the bound clutches in K). Stable + far less Nc-dependent than the "
                    "single-fiber drag. (Earlier --bulk-drag = uniform-small-γ NaN; the old post-hoc override runaway — both superseded by this.)")
    ap.add_argument("--ecm", action="store_true", help="S4 (Phase B): build a collagen-I (Mikado) ECM slab under the cell + attach basal "
                    "clutches to fiber nodes (BUILD-side; the FA↔ECM run() coupling is the next PI-overseen step — FF_S4_INTEGRATION_DESIGN). "
                    "Validated path untouched: with --ecm alone the crawl still uses the fixed-substrate clutch.")
    ap.add_argument("--ecm-fibers", type=int, default=1500, help="collagen-I fiber count for the --ecm Mikado slab")
    ap.add_argument("--ecm-regrip", action="store_true", help="S6 MIGRATION: rebinding basal clutches RE-GRIP the nearest "
                    "CURRENT collagen fibre (KB-1.23 FA-captures-fibres, R_FA≈1.5µm) instead of staying tethered to their "
                    "frame-0 node — so as the cell protrudes, front clutches grab NEW fibres ahead + rear releases → the cell "
                    "WALKS across the collagen (needs --ecm; pair with --com-drag). Default off = validated remodel-in-place path.")
    ap.add_argument("--ecm-lp-um", type=float, default=20.0, help="collagen fiber persistence length Lp [µm] → κ=kBT·Lp "
                    "(the BIPHASIC-sweep stiffness knob; default 20 = LP_COLLAGEN_UM; thin-fibril→thick-bundle range 2–2e5, PI-gated). "
                    "Chan-Odde/Bangasser cross-check: clutch traction should peak at an intermediate effective stiffness.")
    ap.add_argument("--n-fa", type=int, default=0, help="cap basal clutches to N DISCRETE spatially-spread focal-adhesion sites "
                    "(0=one-per-cortex-node). Real FAs are ~10² integrin clusters; at native Nc the per-node model dilutes per-clutch "
                    "load ~100× so clutches never turn over (no crawl). ~150-300 restores physical per-clutch load → the catch-slip treadmill.")
    ap.add_argument("--kmc-every", type=int, default=2000, help="steps between clutch catch-slip turnover + nascent rebind + FA-maturation ticks. "
                    "MUST be << steps or the clutch treadmill never fires (default 2000 suits the explicit tiny-dt path; for the implicit "
                    "large-dt path set so dt·kmc_every ≈ 1-2 s, e.g. ~25 at dt=0.05, to resolve the ~1 s clutch lifetime)")
    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--record-every", type=int, default=100)
    ap.add_argument("--audit", action="store_true", help="also run clutches-OFF control (must give ~0 net drift)")
    ap.add_argument("--static", action="store_true", help="FOUNDATION mode: no protrusion — just adhere + contact")
    ap.add_argument("--spread", action="store_true", help="EMERGENT spreading: peripheral polymerization + clutch (no wetting)")
    ap.add_argument("--mature", action="store_true", help="mature (stable, non-rupturing) basal FA — firm adhesion")
    ap.add_argument("--implicit", action="store_true", help="NF2007 implicit stepping (large dt, unconditionally stable)")
    ap.add_argument("--dt-impl", type=float, default=1.0e-2, help="implicit timestep [s] (clutch-in-K allows up to ~0.2)")
    ap.add_argument("--assembly", action="store_true", help="dynamic cortex area growth (actin assembly) → cell can flatten")
    ap.add_argument("--growth", action="store_true", help="per-filament barbed-end polymerization (KB-3.6 ratchet) → filaments elongate individually")
    ap.add_argument("--myosin-linear", action="store_true", help="explicit Stam-Hocky minifilament (derived 60pN stall + linear force-velocity) instead of swept f_myo")
    ap.add_argument("--substrate-E", type=float, default=0.0, help="compliant substrate Young's modulus [Pa] (movable FA anchors, k_sub∝E; 0=rigid pins; KB-1.5 5000)")
    ap.add_argument("--fa-maturation", action="store_true", help="per-clutch talin unfolding → vinculin → force-gated FA growth/disassembly (KB-2.x mechanosensor)")
    ap.add_argument("--piezo", action="store_true", help="Piezo1 tension-gated open-probability reporter (KB-3.10; diagnostic)")
    ap.add_argument("--treadmill", action="store_true", help="EXPERIMENTAL/REJECTED (2026-07-09): explicit front-bias nascent rebind. "
                    "The directed crawl already EMERGES from protrusion + uniform clutch turnover (set --kmc-every so turnover fires) at "
                    "physiological ~60 nm/s, traction-driven (OFF-audit PASS). This explicit front-bias instead OVER-de-adheres — as the COM "
                    "advances, clutches fall behind it and stop re-forming → bound→0 → traction collapses → LESS crawl. Kept only as a documented dead-end.")
    ap.add_argument("--rear-depoly", action="store_true", help="TREADMILL S1: pointed-end depolymerization at the rear cap "
                    "(mirror of directed front growth; mass-matched) — the rear half of the actin treadmill (mass-conserving).")
    ap.add_argument("--flow", action="store_true", help="TREADMILL S2/S3: molecular-clutch RETROGRADE FLOW — bound-clutch anchors "
                    "drift forward at v_retro in the mesh frame (= actin flows rearward), so the clutch builds forward load from FLOW "
                    "→ traction → compact translocation; catch-slip releases at F*, nascent-rebind reseeds at the front.")
    ap.add_argument("--v-retro", type=float, default=0.03, help="retrograde actin flow speed [µm/s] (Chan-Odde 0.01–0.10; default 0.03 mid-band)")
    ap.add_argument("--fil-length-dist", default="mono", choices=["mono", "exponential"],
                    help="cortex filament length model: mono (identical L) or exponential (KB-3.18 distributed 1–10µm)")
    ap.add_argument("--microtubules", action="store_true",
                    help="Thread-C: merge an MT aster (κ=KAPPA_MT) into the shared implicit K — needs --implicit")
    ap.add_argument("--n-mt", type=int, default=40, help="number of MT tubes radiating from the MTOC (aster size — flag to PI)")
    ap.add_argument("--l-mt", type=float, default=6.0, help="MT arm contour length [µm] (reach toward the R≈7.5µm cortex — flag to PI)")
    ap.add_argument("--from-resting", action="store_true",
                    help="A3: start the adherent run from the VALIDATED resting checkpoint cell (dense myosin //10 + "
                         "membrane reservoir f_excess=0.25 + pre-relax to the ΔP=40 Pa turgor set-point), not a fresh "
                         "un-relaxed geometry — the PI's MCF7-on-substrate physiological baseline")
    ap.add_argument("--relax-steps", type=int, default=6000, help="pre-relaxation steps to the resting set-point (--from-resting)")
    ap.add_argument("--tag", default="crawl")
    ap.add_argument("--out", default="ffn_sim/outputs/ff")
    args = ap.parse_args()
    wp.init(); t0 = time.time()
    if args.microtubules and not args.implicit:
        print("[!] --microtubules needs the implicit solver (MT bending rides K); add --implicit for production.")
    S = build(n_cortex_fil=args.cortex_fil, seed=args.seed, n_fa=args.n_fa, length_dist=args.fil_length_dist,
              microtubules=args.microtubules, n_mt=args.n_mt, L_mt_um=args.l_mt,
              from_resting=args.from_resting, relax_steps=args.relax_steps, relax_device=args.device,
              n_myo_ratio=(10 if args.from_resting else 160), f_excess=(0.25 if args.from_resting else 0.0),
              ecm=args.ecm, ecm_fibers=args.ecm_fibers, ecm_lp_um=args.ecm_lp_um)
    if S.get("ecm") is not None:
        _en = S["ecm_node"]
        print(f"[ecm] collagen-I Mikado slab: {len(S['ecm'].net.fiber_offsets)-1} fibers, {S['ecm'].net.pos.shape[0]} nodes, "
              f"mesh ξ={S['ecm'].mesh_size_um:.2f} µm; basal clutches attached to fibers: {int((_en>=0).sum())}/{_en.size} "
              f"(BUILD-side; run() FA↔ECM coupling = next PI-overseen step)")
    mt_note = (f"; MT aster {S['n_mt']} tubes → Ne {S['Ne']} (Nmt+MTOC {S['Ne']-S['Nc']})" if args.microtubules else "")
    if S.get("resting") is not None:
        rr = S["resting"]
        print(f"[from-resting] pre-relaxed to checkpoint: ΔP={rr['dP_Pa']:.1f} Pa  γ={rr['gamma_mN_m']:.3f} mN/m  "
              f"V/V0={rr['V_over_V0']:.3f}  R_eq={rr['R_eq_um']:.2f} µm  (myosin //10, membrane f_excess=0.25)")
    print(f"[build] cortex {S['Nc']} + nucleus {S['n_nuc']}{mt_note}; basal FA clutches {S['basal'].size}; "
          f"front-cap nodes {S['front'].size}; f_pro {S['f_pro']:.1f} pN/node; z_sub {S['z_sub']:.2f}  "
          f"({time.time()-t0:.0f}s)")
    r = run(S, steps=args.steps, record_every=args.record_every, clutches=True,
            protrude=(not args.static and not args.spread), spread=args.spread,
            rupture=not args.mature, implicit=args.implicit, dt_impl=args.dt_impl,
            assembly=args.assembly, growth=args.growth, myosin_linear=args.myosin_linear, substrate_E=args.substrate_E,
            fa_maturation=args.fa_maturation, treadmill=args.treadmill, bulk_drag=args.bulk_drag,
            rear_depoly=args.rear_depoly, flow=args.flow, v_retro_um_s=args.v_retro,
            com_drag=args.com_drag, ecm_regrip=args.ecm_regrip, kmc_every=args.kmc_every, device=args.device)
    tag_mode = "SPREAD" if args.spread else ("STATIC adhere" if args.static else "CRAWL clutch ON")
    print(f"[{tag_mode}] dt={r['dt']*1e3:.3g} ms  T={r['times'][-1]:.1f} s  "
          f"disp∥={r['disp_along_um']:+.3f} µm  v_crawl={r['v_crawl_nm_s']:+.2f} nm/s  "
          f"traction={r['traction_nN']:.2f} nN  bound={r['bound_frac']:.2f}  wall {r['wall_s']:.0f}s")
    print(f"[CONTACT] V/V0={r['vol_final']:.3f}  basal-gap={r['basal_gap_um']:+.3f} µm (0=flat contact)  "
          f"contact-radius={r['contact_radius_um']:.2f} µm  n-contact={r['n_contact']}  "
          f"cell-height={r['cell_height_um']:.2f} µm (R={r['R_um']:.1f})  → "
          f"{'STABLE ADHERED' if 0.8 < r['vol_final'] < 1.25 and r['basal_gap_um'] < 0.4 and r['bound_frac'] > 0.5 else 'NOT STABLE/ADHERED'}")
    if r.get("ecm_metrics") is not None:                      # S6: collagen matrix REMODELING under the cell's clutch traction (DECOMPOSED)
        _m = r["ecm_metrics"]; _en = np.asarray(r["ecm_node"])
        print(f"[ECM-REMODEL] gripped-node |disp|={_m['mag_nm']:.1f} nm  →  COHERENT-inward: recruit(toward-actin)={_m['recruit_nm']:+.1f} nm  "
              f"densification(footprint-radial)={_m['dens_nm']:+.1f} nm  settling(dz)={_m['dz_nm']:+.1f} nm  coherence={_m['coh']:.2f}")
        print(f"[ECM-ALIGN]   fiber radial-alignment index (near footprint) {_m['rai0']:.3f} → {_m['raif']:.3f}  "
              f"(Δ={_m['drai']:+.3f}; >0 = fibers reorient toward the cell)  attached={_m['n_grip']}/{_en.size}")
    if args.piezo:                                             # Piezo1 tension reporter (KB-3.10, diagnostic; feedback OFF)
        pz = resolve_piezo(); g_mem = S["mem"].gamma_mem
        print(f"[Piezo] membrane tension {g_mem:.1f} pN/µm → P_open={float(p_open(g_mem, pz)):.4f} (rest≈closed; opens as tension→γ_half=5000; feedback gain=0)")
    out = {"clutch_on": {k: (r[k] if not isinstance(r[k], np.ndarray) else None)
                          for k in ("dt", "disp_along_um", "disp_perp_um", "v_crawl_nm_s", "traction_nN",
                                    "bound_frac", "n_clutch", "f_pro_pN", "gamma_min", "gamma_max",
                                    "F_star_pN", "kmax", "wall_s")}}
    if r.get("ecm_metrics") is not None:                      # S6: the remodel decomposition (recruit/densification/alignment) for the biphasic sweep
        out["ecm_metrics"] = {k: float(v) for k, v in r["ecm_metrics"].items()}
        out["ecm_lp_um"] = float(args.ecm_lp_um)
    _npz = dict(frames=np.array(r["frames"]), com=r["com"], times=r["times"], vol=r["vol"], Nc=S["Nc"], Ne=S["Ne"],
                basal=S["basal"], z_sub=S["z_sub"], R=S["R"], phat=S["phat"], foff=S["net"].fiber_offsets,
                n_nuc=S["n_nuc"], n_mt=S["n_mt"], mtoc_idx=S["mtoc_idx"], svm=r["svm"], cstrain=r["cstrain"],
                faces=r["faces"], bound_frames=r["bound_frames"], anch_frames=r["anch_frames"],
                k_int=r["k_int"], clutch_rest_um=r["clutch_rest_um"])
    if r.get("ecm_posf") is not None:                          # S6: collagen frame-0 + remodeled + attachment (for the remodel viewer)
        _npz.update(ecm_pos0=r["ecm_pos0"], ecm_posf=r["ecm_posf"], ecm_node=np.asarray(r["ecm_node"]),
                    ecm_foff=S["ecm"].net.fiber_offsets)
        if r.get("ecm_frames") is not None:
            _npz.update(ecm_frames=r["ecm_frames"])            # the collagen recruitment PROCESS (per-frame positions)
    np.savez_compressed(f"{args.out}/figs/{args.tag}_on.npz", **_npz)
    if args.audit:
        rc = run(S, steps=args.steps, record_every=args.record_every, clutches=False, device=args.device)
        print(f"[AUDIT clutch OFF] disp∥={rc['disp_along_um']:+.3f} µm  disp⊥={rc['disp_perp_um']:.3f} µm  "
              f"v_crawl={rc['v_crawl_nm_s']:+.2f} nm/s  (must be ≈0: protrusion alone is internal)")
        if rc.get("ecm_metrics") is not None:                  # S6: OFF control — no clutch grip → the collagen must NOT remodel (traction-driven proof)
            _mo = rc["ecm_metrics"]; _mn = r["ecm_metrics"]
            print(f"[ECM-REMODEL AUDIT OFF] recruit={_mo['recruit_nm']:+.1f} nm  densification={_mo['dens_nm']:+.1f} nm  "
                  f"|disp|={_mo['mag_nm']:.1f} nm  (must be ≈0)")
            _rr = abs(_mn['recruit_nm']) / max(abs(_mo['recruit_nm']), 1e-6)
            _ev = "PASS" if abs(_mo['recruit_nm']) < 0.2 * abs(_mn['recruit_nm']) + 1.0 else "FAIL"
            print(f"[ECM-REMODEL AUDIT verdict] traction-driven recruit ratio |ON/OFF|={_rr:.1f}×  → {_ev}")
            out["ecm_remodel_verdict"] = _ev
        ratio = abs(r["disp_along_um"]) / max(abs(rc["disp_along_um"]), 1e-6)
        verdict = "PASS" if abs(rc["disp_along_um"]) < 0.2 * abs(r["disp_along_um"]) + 0.02 else "FAIL"
        print(f"[AUDIT verdict] traction-driven ratio |on/off|={ratio:.1f}×  → {verdict}")
        out["clutch_off"] = {k: rc[k] for k in ("disp_along_um", "disp_perp_um", "v_crawl_nm_s")}
        out["audit_verdict"] = verdict
        _offnpz = dict(frames=np.array(rc["frames"]), com=rc["com"], times=rc["times"], Nc=S["Nc"], Ne=S["Ne"],
                       z_sub=S["z_sub"], phat=S["phat"], foff=S["net"].fiber_offsets, n_nuc=S["n_nuc"],
                       n_mt=S["n_mt"], mtoc_idx=S["mtoc_idx"])
        if rc.get("ecm_posf") is not None:                     # S6: OFF-control collagen (should equal frame-0 ⇒ visual proof of traction-driven)
            _offnpz.update(ecm_pos0=rc["ecm_pos0"], ecm_posf=rc["ecm_posf"], ecm_node=np.asarray(rc["ecm_node"]),
                           ecm_foff=S["ecm"].net.fiber_offsets)
        np.savez_compressed(f"{args.out}/figs/{args.tag}_off.npz", **_offnpz)
    json.dump(out, open(f"{args.out}/figs/{args.tag}.json", "w"), indent=2)
    print(f"wrote {args.tag}_on.npz + {args.tag}.json  (total {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
