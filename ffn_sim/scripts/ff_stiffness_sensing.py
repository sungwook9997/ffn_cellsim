"""FF substrate-stiffness sensing — a resting full-compartment MCF7 cell reads substrate rigidity E and the
engaged-clutch TRACTION traces the Bangasser-Odde motor-clutch BIPHASIC optimum (quasi-static, PATH-b Winkler).

PI 2026-07-10: the cell must SENSE how stiff the substrate is. A rigid dish (``substrate_plane_kernel`` +
FIXED FA anchors) makes stiffness structurally impossible — every anchor is infinitely rigid, so traction can
only rise with contraction, never read E. Here stiffness enters through **ONE knob**: each focal-adhesion
anchor is a movable Winkler-spring DOF tied to the dish by ``k_sub = 2·E·a/(1−ν²)`` (flat-punch on an elastic
half-space; :func:`ff.substrate.resolve_substrate`). When the myosin-contracted cortex pulls a basal clutch,
its anchor gives way by F/k_sub, so the WHOLE molecular-clutch traction becomes E-dependent — soft → the anchor
slips (low force per bond, bonds live long), stiff → the anchor holds (high force per bond, bonds load-and-fail
via the Pereverzev catch-slip KMC). Total traction = (bound fraction)·(force per clutch) therefore peaks at an
INTERMEDIATE stiffness: the Bangasser-Odde biphasic optimum, EMERGENT from the fine-grained clutches (no lumped
traction law). This is quasi-static: the cortex + turgor + membrane + nucleus + substrate relax to their
adhered steady state at each E and we read the engaged traction, bound fraction, and top-down spread vs E.

PATH (b) — the compliant Winkler substrate is the robust single-knob E axis: stiffness enters ONLY through
``k_sub`` on the FA anchors, NOT through building an ECM network (that would double-count compliance). We do NOT
build an :mod:`ff.ecm_library` fiber network here. Sweeping E across the library's PA-gel ladder
(:data:`ff.ecm_library.PA_FORMULATIONS`, 150 Pa brain → 40 kPa stiff-osteoid) keeps the axis library-connected.
PATH (a) — ``sense_ecm_network`` (``--mode ecm-network``): the cell instead grips the actual ``ecm_library``
ECM network built UNDER it (a soft PA-gel E=150 Pa vs a stiff one E=40 kPa, collagen, …), so the substrate
stiffness is the ECM's OWN calibrated stiffness and traction emerges from the TWO-sided cell↔ECM clutch coupling
(:func:`ff.fa_ecm.clutch_ecm_spring_kernel`) plus the ECM's own deformation (its bending + segment/crosslink
springs relax against the pulled clutches each step). This is the deepest "cell senses the real library ECM"
demonstration; the Winkler path (b) stays the robust single-knob analytic axis.

Recipe (extracted from ``ff_crawl_on_substrate.build``/``run`` — same cortex force-loop + readout):
  1. build a crosslinked cortex, seed a nucleus + membrane, RE-RELAX to the resting turgor set-point
     (ΔP=TURGOR_DP0=40 Pa, biphasic drained-solid K=300 Pa) via
     :func:`network_warp.simulate_whole_cell_compression_on_device` — that relaxed cell IS the resting baseline
     (physiological-baseline HARD rule: adhere FROM the resting state). MT OFF → layout [cortex(Nc); nucleus].
  2. derive the basal FA cap + fixed dish anchors from the RELAXED cortex geometry; cap to ~200 discrete FAs.
  3. per E: run a quasi-static implicit loop (bending + crosslink + myosin + osmotic turgor + membrane + nucleus
     + substrate plane + clutch spring), moving the anchors under the Winkler compliance
     (:func:`ff.substrate.substrate_anchor_equilibrium_kernel`) and turning the clutches over with the
     Pereverzev catch-slip KMC + nascent-adhesion rebind. Read traction / bound-fraction / silhouette.

No parameter is tuned to an outcome: **E is the input, k_sub is the DERIVED Winkler relation, and the biphasic
optimum is MEASURED** (never a swept target). Full compartments are ON at their physiological set-points
(turgor / membrane / nucleus). The ABSOLUTE traction carries the adhesion-patch radius ``a`` (KB-PIV-5, the same
missing-datum family as ρ_L) — surfaced, not fitted; only the RELATIVE optimum (which E maximises traction) is
``a``-invariant and claimed here.

⚠️ NATIVE-scale rule (CLAUDE.md HARD): a coarse cortex (~800–2000 filaments) LACKS structural stability and is a
dev-only smoke — every finding MUST be reconfirmed at NATIVE ``--nf 38000`` on the A5000. The driver PRINTS this
whenever ``--nf`` < 38000.

Sanity Gate:
  - Dimensional: k_sub = 2·E·a/(1−ν²) → [Pa·µm] = [pN/µm] (1 Pa≡1 pN/µm²); traction Σ k_int·ΔL [pN] → nN (÷1e3).
  - Boundary: E→∞ ⇒ k_sub→∞ ⇒ anchor→fixed dish point (rigid-pin path recovered); E→0 ⇒ k_sub→0 ⇒ anchor follows
    the actin ⇒ zero clutch extension ⇒ zero traction. Both are limits of the same series spring.
  - Conservation: cytoplasm is incompressible — the cortex volume is HARD-projected to V0 every step (turgor).
  - Sign: the clutch spring pulls the actin node TOWARD its anchor (traction on the cell); myosin CONTRACTS the
    cortex, which is what loads the clutches. A rigid substrate ⇒ monotone traction(E); the biphasic optimum
    requires the catch-slip KMC (bonds rupture faster past F*), so the KMC MUST fire (dt·kmc_every ≈ 1–2 s).
  - Measurement: traction = engaged-clutch extension × k_int (identical readout to ``ff_crawl_on_substrate.run``).

Smoke (dev): ``python -m ffn_sim.scripts.ff_stiffness_sensing --nf 800 --steps 200 --relax-steps 400 \
    --E-list 150,8000 --device cpu``
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import warp as wp
from scipy.spatial import ConvexHull

from ffn_sim.ff.gamma_floor import (build_crosslinked_cortex, CortexParams, TURGOR_DP0, TURGOR_PI_IN0,
                                    VMIN_FRAC, NMIIA_MINIFIL_STALL_PN)
from ffn_sim.ff.network_warp import (_zero, link_spring_kernel, myosin_kernel, turgor_kernel, reshape_kernel,
                                     nucleus_shell_kernel, substrate_plane_kernel,
                                     simulate_whole_cell_compression_on_device)
from ffn_sim.ff.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from ffn_sim.ff.fa_clutch_warp import clutch_spring_kernel, clutch_catchslip_kmc_kernel, resolve_clutch
from ffn_sim.ff.motility_warp import (axpy_physical_kernel, cortex_volume_kernel, sum_pos_kernel,
                                      sum_radius_kernel, xl_turnover_kernel, physical_node_gammas,
                                      volume_gradient, crawl_cfl_dt)
from ffn_sim.ff.implicit_ff import implicit_step_current
from ffn_sim.ff.substrate import resolve_substrate, substrate_anchor_equilibrium_kernel
from ffn_sim.common.compartments import resolve_nucleus, resolve_membrane
from ffn_sim.ff.ecm_library import PA_FORMULATIONS, get_spec, build_ecm
from ffn_sim.ff.fa_ecm import (attach_clutches_to_ecm, clutch_ecm_spring_kernel,
                               mask_ecm_by_bound_kernel, traction_director_anisotropy, _inplane_unit)
from ffn_sim.ff.ecm_mechanics import ecm_material_stress

NATIVE_NF = 38000                 # native cortex filament count (the ONLY authoritative basis; CLAUDE.md HARD)
DEFAULT_E_LIST = (150.0, 500.0, 2000.0, 8000.0, 40000.0)   # Pa: brain-PA → stiff muscle/osteoid-PA (library ladder)
ECM_CFL_SAFETY = 0.01             # ecm-network explicit-substep CFL safety fraction of 0.1/kmax. 0.01 (10× under
                                  # the ff_ecm_remodel_demo's 0.1) leaves headroom for the k·(L−r0)/L spring's 1/L
                                  # geometric stiffening when a Delaunay bond transiently compresses (L<r0) during
                                  # relaxation — a numerical-stability margin, grid-invariant, physics-answer-invariant.


def _fps_subsample(pts: np.ndarray, k: int) -> np.ndarray:
    """Farthest-point subsample: pick ``k`` well-spread points from ``pts`` (M,3). Deterministic (no RNG).

    Places DISCRETE focal-adhesion sites — real FAs are ~10² discrete integrin clusters, not one clutch per
    cortex node; at native Nc the per-node model dilutes the per-clutch load ~100× so clutches never turn over.
    Returns local indices into ``pts`` (copy of ``ff_crawl_on_substrate._fps_subsample``)."""
    n = pts.shape[0]
    if k >= n:
        return np.arange(n)
    sel = [n // 2]
    d = np.full(n, np.inf)
    for _ in range(1, k):
        d = np.minimum(d, np.sum((pts - pts[sel[-1]]) ** 2, axis=1))
        sel.append(int(np.argmax(d)))
    return np.asarray(sel, dtype=np.int64)


def _oriented_faces(v: np.ndarray) -> np.ndarray:
    """Outward-oriented ConvexHull triangulation of cortex nodes ``v`` (Nc,3) → (T,3) int32.

    scipy's hull simplices are not consistently wound; flip any face whose normal points inward so the signed
    tetra-volume sum (``cortex_volume_kernel`` / :func:`volume_gradient`) is a positive enclosed volume. Copy of
    the orientation block in ``ff_crawl_on_substrate.run``."""
    tri = ConvexHull(v).simplices.copy()
    cen = v.mean(0)
    nrm = np.cross(v[tri[:, 1]] - v[tri[:, 0]], v[tri[:, 2]] - v[tri[:, 0]])
    flip = (nrm * (v[tri[:, 0]] - cen)).sum(1) < 0
    tri[flip] = tri[flip][:, [0, 2, 1]]
    return np.ascontiguousarray(tri, np.int32)


def _pa_label(E_pa: float) -> tuple[str, str]:
    """Nearest PA-gel formulation (log-distance) in :data:`ff.ecm_library.PA_FORMULATIONS` → (recipe, tissue).

    Keeps the E axis library-connected: each swept E is annotated with the polyacrylamide recipe (acrylamide/bis)
    that casts that Young's modulus + the tissue that stiffness mimics (DOI-verified, KB-1.V.4.1)."""
    best = min(PA_FORMULATIONS, key=lambda r: abs(np.log(r[3]) - np.log(max(E_pa, 1e-9))))
    return best[0], best[4]


def build_resting_cell(nf: int = 2000, *, relax_steps: int = 4000, n_fa: int = 200, seed: int = 1,
                       contact_h: float = 0.6, device: str = "cpu") -> dict:
    """Build the resting full-compartment cell (E-INDEPENDENT) once, to be re-used across every substrate E.

    Follows the ``ff_crawl_on_substrate.build(from_resting=True)`` recipe with microtubules OFF: a crosslinked
    cortex + dense myosin (//10), a stiff nucleus (0.70·R, 3000 beads) and a reservoir-buffered membrane
    (f_excess=0.25), RE-RELAXED to the resting turgor set-point (ΔP=TURGOR_DP0, biphasic drained-solid K=300 Pa)
    via :func:`network_warp.simulate_whole_cell_compression_on_device` (which seeds the nucleus internally →
    ``pos_all`` layout ``[cortex(Nc); nucleus(n_nuc)]``, Ne==Nc). The relaxed positions ARE the resting cell; the
    basal FA cap + fixed dish anchors are then derived from that RELAXED cortex geometry.

    Args:
        nf: cortex filament count (native = 38000; coarse = dev smoke only, non-authoritative).
        relax_steps: pre-relaxation steps to the resting set-point.
        n_fa: cap the basal cap to N discrete focal-adhesion sites (0 = one clutch per basal node).
        seed: cortex/crosslinker RNG seed.
        contact_h: basal-cap height above the substrate plane [µm].
        device: ``"cpu"`` or ``"cuda:0"``.

    Returns:
        dict with the resolved cortex ``cx``, relaxed ``pos_all``, ``Nc``/``n_nuc``, ``nuc``/``mem`` params, the
        basal FA indices + fixed dish ``anchors``, ``z_sub``, ``R``, centroid ``c``, and the resting metrics.
    """
    rng = np.random.default_rng(seed)
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=max(1, nf // 10),
                                  length_dist="mono", rng=rng)
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    R = cx.R0_mean
    Nc = cx.net.n_nodes
    nuc = resolve_nucleus(R_nuc_um=0.70 * R, n_beads=3000)      # MCF7 nucleus (0.70·R, Moore2016; E_nuc 399 Pa)
    mem = resolve_membrane(f_excess=0.25)                      # reservoir-buffered baseline tension (checkpoint)
    # RE-RELAX to the validated resting checkpoint — this IS the resting cell (nucleus seeded internally,
    # nucleus_seed=3; MT OFF ⇒ Ne==Nc ⇒ pos_all = [cortex(Nc) ; nucleus(n_nuc)]).
    pos_all, mr = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=0.0, nucleus=nuc, membrane=mem, microtubule=None,
        pressure_setpoint=float(TURGOR_DP0), K_drained_Pa=300.0, n_steps=relax_steps,
        turgor_every=50, nucleus_seed=3, device=device)
    pos_all = np.ascontiguousarray(pos_all, np.float64)
    n_nuc = nuc.n_beads
    cx.net.pos = pos_all[:Nc].copy()                          # keep net.pos consistent (V0 + faces read it)
    # basal FA cap + fixed dish anchors from the RELAXED cortex geometry
    cortex_pos = pos_all[:Nc]
    z_sub = float(cortex_pos[:, 2].min())
    basal = np.where(cortex_pos[:, 2] < z_sub + contact_h)[0]
    if n_fa and 0 < n_fa < basal.size:                        # discrete FA sites → physical per-clutch load
        basal = basal[_fps_subsample(cortex_pos[basal], n_fa)]
    anchors = cortex_pos[basal].copy(); anchors[:, 2] = z_sub  # fixed dish rest points (projected to the plane)
    return dict(cx=cx, Nc=Nc, n_nuc=n_nuc, pos_all=pos_all, nuc=nuc, mem=mem,
                basal=basal.astype(np.int32), anchors=anchors, z_sub=z_sub, R=R,
                c=cortex_pos.mean(0), resting=dict(dP_Pa=float(mr["dP_turgor_Pa"]),
                gamma_mN_m=float(mr["gamma_apparent_mN_m"]), V_over_V0=float(mr["V_over_V0"]),
                R_eq_um=float(mr["R_eq_um"])))


def sense(S: dict, E_pa: float, *, steps: int = 1500, dt: float = 0.05, kmc_every: int = 0,
          a_um: float = 0.05, nu: float = 0.45, f_myo: float = NMIIA_MINIFIL_STALL_PN,
          koff_xl: float = 0.4, refresh_every: int = 50, reshape_every: int = 20, xl_turn_every: int = 50,
          adh_h: float = 0.4, device: str = "cpu") -> dict:
    """Quasi-static stiffness-sensing loop on a Winkler substrate of modulus ``E_pa`` → engaged traction & co.

    Copies the cortex force accumulation + integration of ``ff_crawl_on_substrate.run`` (implicit host path),
    stripped to the adhesion-sensing forces (NO protrusion / spreading / treadmill / ECM). Stiffness enters ONLY
    through ``k_sub = 2·E·a/(1−ν²)`` moving the FA anchors to the clutch↔substrate SERIES equilibrium each step
    (:func:`ff.substrate.substrate_anchor_equilibrium_kernel`). The Pereverzev catch-slip KMC ruptures loaded
    bonds and nascent adhesions re-form, so the engaged population reaches an E-dependent steady state.

    Args:
        S: resting cell from :func:`build_resting_cell` (re-used across all E).
        E_pa: substrate Young's modulus [Pa] (= pN/µm²).
        steps: quasi-static implicit steps.
        dt: implicit timestep [s] (large-dt stable via clutch-in-K + rank-1 osmotic).
        kmc_every: steps between catch-slip KMC + nascent-rebind ticks (0 = auto → dt·kmc_every ≈ 1.5 s).
        a_um: adhesion-patch radius [µm] (surfaced magnitude lever; the relative optimum is a-invariant).
        nu: substrate Poisson ratio.
        f_myo: myosin per-minifilament stall force [pN] (loads the clutches).
        koff_xl: α-actinin crosslink off-rate [1/s] (cortex fluidisation / stress relaxation).
        refresh_every / reshape_every / xl_turn_every: slow-mode / inextensibility / turnover cadences.
        adh_h: nascent-adhesion capture height above the plane [µm].
        device: ``"cpu"`` or ``"cuda:0"``.

    Returns:
        dict(E_pa, k_sub, traction_nN, bound_frac, force_per_clutch_pN, silhouette_um2, n_clutch, n_bound,
             kmc_every, dt, wall_s, diverged).
    """
    d = device
    cx = S["cx"]; net = cx.net; Nc = S["Nc"]; n_nuc = S["n_nuc"]; nuc = S["nuc"]; mem = S["mem"]
    Ne = Nc; N = Nc + n_nuc                                    # MT OFF ⇒ Ne==Nc; nucleus beads follow the cortex
    cp = resolve_clutch()
    substrate = resolve_substrate(E_pa=E_pa, nu=nu, a_um=a_um)  # the single stiffness knob: k_sub ∝ E
    if kmc_every <= 0:
        kmc_every = max(1, int(round(1.5 / dt)))              # dt·kmc_every ≈ 1.5 s → resolve the ~1 s clutch life

    # ---- cortex elastic topology (cortex only; no MT hub) ----
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(net.fiber_offsets) - 1
    soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean())
    xl_ij_np = np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int64)
    kxl_np = np.ascontiguousarray(cx.xl_k, np.float64)
    xl_rest_np = np.ascontiguousarray(cx.xl_rest, np.float64)
    n_xl = xl_ij_np.shape[0]
    bend_triples_np = np.ascontiguousarray(net.bend_triples, np.int64)

    # ---- turgor volume reference + oriented faces from the RELAXED cortex (run() convention: hull volume) ----
    pos0 = S["pos_all"].copy()
    faces = _oriented_faces(pos0[:Nc])
    V0 = float(ConvexHull(pos0[:Nc]).volume); vmin = VMIN_FRAC * V0
    k_plane = 1.0e3                                            # substrate excluded-volume stiffness [pN/µm]

    # ---- per-node drag + implicit-solver scalars ----
    gammas = physical_node_gammas(net, Nc, n_nuc)             # [cortex log-drag ; nucleus Stokes beads]
    gamma_rep = float(np.median(gammas[:Nc]))
    kmax = max(float(net.kappa.max()) / seg**3, float(kxl_np.max()) if n_xl else 1.0,
               nuc.k_chrom + nuc.k_lamin, cp.k_int, k_plane)
    dt_cfl = crawl_cfl_dt(gammas, kmax, safety=0.1)           # explicit CFL (logged; the implicit step uses dt)

    # ---- device arrays ----
    pos_d = wp.array(pos0, dtype=wp.vec3d, device=d); f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_d = wp.array(np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(kxl_np, dtype=wp.float64, device=d)
    r0_d = wp.array(xl_rest_np.copy(), dtype=wp.float64, device=d)   # crosslink rest length (Maxwell-relaxed in place)
    faces_d = wp.array(faces, dtype=wp.int32, ndim=2, device=d)
    csum_d = wp.zeros(3, dtype=wp.float64, device=d); rsum_d = wp.zeros(1, dtype=wp.float64, device=d)
    has_myo = cx.myo_i.size > 0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cx.myo_i, cx.myo_j], 1), np.int32),
                         dtype=wp.int32, ndim=2, device=d)
    # FA clutches: actin-end index, movable anchor, its fixed dish rest point, bound state
    M = S["basal"].size
    ac_d = wp.array(S["basal"], dtype=wp.int32, device=d)
    anch_d = wp.array(S["anchors"], dtype=wp.vec3d, device=d)
    anch_rest_d = wp.array(np.ascontiguousarray(S["anchors"], np.float64), dtype=wp.vec3d, device=d)
    bd_d = wp.array(np.ones(M, np.int32), dtype=wp.int32, device=d)
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    z_sub = S["z_sub"]
    dP_mem_area = 0.0                                          # membrane inward tension per node (set at refresh)
    frac_xl = 1.0 - np.exp(-koff_xl * dt * xl_turn_every)      # crosslink turnover fraction per turnover tick

    def full_force(x_np):
        """Full FF force (elastic + active) at x — the implicit solver RHS. The osmotic/turgor force is the EXACT
        volume-gradient ΔP·g added host-side (clean rank-1 stiffness for the implicit solve); the membrane inward
        tension stays a kernel. Mirrors ``ff_crawl_on_substrate.run.full_force`` minus protrusion/spreading."""
        pos_d.assign(np.ascontiguousarray(x_np, np.float64).reshape(N, 3))
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)  # membrane inward
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        wp.launch(clutch_spring_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                  wp.float64(cp.rest_um), f_d], device=d)
        wp.synchronize_device(d)
        F = f_d.numpy().reshape(-1)
        pcx = np.ascontiguousarray(x_np, np.float64).reshape(N, 3)[:Nc]; cen = pcx.mean(0)
        aa = pcx[faces[:, 0]] - cen; bb = pcx[faces[:, 1]] - cen; ccf = pcx[faces[:, 2]] - cen
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        dPv = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dPv = min(max(dPv, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        F[:3 * Nc] += (dPv * volume_gradient(pcx, faces, cen)).reshape(-1)
        return F

    t0 = time.time(); diverged = False
    for step in range(steps):
        if step % refresh_every == 0:                          # slow modes via device reductions (centroid, R_mean)
            csum_d.zero_(); wp.launch(sum_pos_kernel, dim=Nc, inputs=[pos_d, csum_d], device=d)
            cc = csum_d.numpy() / Nc
            centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2]))
            rsum_d.zero_(); wp.launch(sum_radius_kernel, dim=Nc, inputs=[pos_d, centre, rsum_d], device=d)
            Rm = float(rsum_d.numpy()[0]) / Nc
            area = 4.0 * np.pi * Rm**2
            dP_mem_area = -2.0 * min(mem.gamma_mem, mem.tau_lysis) / max(Rm, 1e-9) * area / Nc
            if not np.isfinite(Rm) or Rm > 50.0 * S["R"]:
                print(f"    [!] divergence at step {step} (Rm={Rm:.2f}) — truncating"); diverged = True; break
        # ---- implicit host step (unconditionally stable; large dt) ----
        xv = pos_d.numpy().reshape(-1)
        pcx2 = xv.reshape(N, 3)[:Nc]; cen2 = pcx2.mean(0)
        aa = pcx2[faces[:, 0]] - cen2; bb = pcx2[faces[:, 1]] - cen2; ccf = pcx2[faces[:, 2]] - cen2
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        if not np.isfinite(Vc) or Vc <= 1.02 * vmin:
            print(f"    [!] volume collapse at step {step} (Vc={Vc:.1f} ≤ vmin={vmin:.1f}) — truncating")
            diverged = True; break
        gN = volume_gradient(pcx2, faces, cen2)
        k_vol = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-3 * V0) ** 2
        vol_g = np.zeros(3 * N); vol_g[:3 * Nc] = gN.reshape(-1)
        _de = np.zeros(3 * N)                                  # clutch + substrate-plane stiffness → implicit K diagonal
        _bd = bd_d.numpy(); _bb = S["basal"][_bd > 0]
        _de[3 * _bb] = cp.k_int; _de[3 * _bb + 1] = cp.k_int; _de[3 * _bb + 2] = cp.k_int
        _de[3 * S["basal"] + 2] += k_plane
        xv, _ = implicit_step_current(xv, full_force, bend_triples_np, alpha, xl_ij_np, kxl_np,
                                      gamma=gamma_rep, dt=dt, n_newton=1, vol_g=vol_g, k_vol=k_vol, diag_extra=_de)
        # HARD incompressibility: project the cortex back to V=V0 (Newton on the volume constraint along g=∂V/∂x)
        xr = xv.reshape(N, 3); pcx3 = xr[:Nc].copy(); cen3 = pcx3.mean(0)
        for _ in range(3):
            a3 = pcx3[faces[:, 0]] - cen3; b3 = pcx3[faces[:, 1]] - cen3; c3 = pcx3[faces[:, 2]] - cen3
            Vp = abs(float((a3 * np.cross(b3, c3)).sum() / 6.0))
            gp = volume_gradient(pcx3, faces, cen3); den = float((gp * gp).sum())
            if den > 1e-9:
                pcx3 = pcx3 - ((Vp - V0) / den) * gp
        xr[:Nc] = pcx3; xv = xr.reshape(-1)
        pos_d.assign(np.ascontiguousarray(xv, np.float64).reshape(N, 3))
        # ---- COMPLIANT SUBSTRATE: move each movable FA anchor to the clutch↔substrate SERIES equilibrium ----
        wp.launch(substrate_anchor_equilibrium_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, anch_rest_d, bd_d,
                  wp.float64(cp.k_int), wp.float64(substrate.k_sub)], device=d)
        if step % reshape_every == 0:                          # NF2007 §5.3 inextensibility
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        if step % xl_turn_every == 0 and step > 0:             # α-actinin crosslink turnover (fluidise → stress relax)
            wp.launch(xl_turnover_kernel, dim=n_xl, inputs=[pos_d, xl_d, r0_d, wp.float64(frac_xl)], device=d)
        if step % kmc_every == 0 and step > 0:                 # catch-slip rupture + nascent-adhesion rebind
            wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                      wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                      wp.float64(dt * kmc_every), wp.int32(step)], device=d)
            p = pos_d.numpy(); bd = bd_d.numpy(); anchors = anch_d.numpy(); rest = anch_rest_d.numpy()
            can = (bd == 0) & (p[S["basal"], 2] < z_sub + adh_h)
            p_on = 1.0 - np.exp(-dt * kmc_every * cp.k_on)
            reb = can & (np.random.default_rng(step).random(M) < p_on)
            if reb.any():                                      # a re-forming clutch adheres at its CURRENT basal xy,
                anchors[reb, 0] = p[S["basal"][reb], 0]; anchors[reb, 1] = p[S["basal"][reb], 1]
                anchors[reb, 2] = z_sub
                rest[reb] = anchors[reb]                        # its NEW dish rest point (correct Winkler reference)
                bd[reb] = 1
                anch_d.assign(anchors); anch_rest_d.assign(rest); bd_d.assign(bd)

    wp.synchronize_device(d)
    xp = pos_d.numpy(); bd = bd_d.numpy(); anchors = anch_d.numpy()
    L = np.linalg.norm(xp[S["basal"]] - anchors, axis=1)
    Fclutch = cp.k_int * np.maximum(L - cp.rest_um, 0.0) * bd  # engaged-clutch extension force [pN]
    n_bound = int(bd.sum())
    sil = float(ConvexHull(xp[:Nc][:, :2]).volume)             # top-down xy silhouette area [µm²] (2-D hull)
    return dict(E_pa=float(E_pa), k_sub=float(substrate.k_sub), traction_nN=float(Fclutch.sum() / 1e3),
                bound_frac=float(bd.mean()), force_per_clutch_pN=float(Fclutch[bd > 0].mean()) if n_bound else 0.0,
                silhouette_um2=sil, n_clutch=int(M), n_bound=n_bound, kmc_every=int(kmc_every),
                dt=float(dt), dt_cfl_s=float(dt_cfl), wall_s=float(time.time() - t0), diverged=bool(diverged))


def sense_ecm_network(S: dict, material: str, *, E_gel_Pa: float | None = None, conc: float | None = None,
                      steps: int = 1500, dt: float = 0.05, kmc_every: int = 0,
                      f_myo: float = NMIIA_MINIFIL_STALL_PN, koff_xl: float = 0.4, refresh_every: int = 50,
                      reshape_every: int = 20, xl_turn_every: int = 50, adh_h: float = 0.4,
                      capture_um: float = 1.0, ecm_substeps: int = 20, target_z: float = 3.2,
                      alignment_S: float = 0.0, director: tuple = (1.0, 0.0, 0.0),
                      seed: int = 7, device: str = "cpu") -> dict:
    """PATH (a): quasi-static stiffness sensing against a LIVE :mod:`ff.ecm_library` ECM network under the cell.

    Same resting full-compartment cell + same cortex force-loop as :func:`sense`, but the FA clutches grip an
    ACTUAL ECM network (a soft PA-gel E=150 Pa vs a stiff one E=40 kPa, or collagen) built directly under the
    basal cap, instead of a Winkler-spring anchor. Stiffness is now the ECM's OWN calibrated stiffness — for a
    continuum gel the Delaunay bond stiffness is the analytic inverse of ``E_gel_Pa``
    (:func:`ff.ecm_library.build_continuum_ecm`); for a fibrillar matrix it EMERGES from the fiber microstructure
    at ``conc``. Traction emerges from the TWO-sided cell↔ECM coupling (:func:`ff.fa_ecm.clutch_ecm_spring_kernel`,
    +f on the cell actin, −f on the fiber node) and the ECM's deformation: each cell step the cortex is stepped
    implicitly (clutch pulling toward the live ECM node), then the ECM relaxes for ``ecm_substeps`` explicit
    substeps under those same clutch reactions + its own bending/segment/crosslink springs + a pinned far face.
    A stiff ECM barely gives → the clutch stays extended → high sustained traction; a soft ECM deforms → the
    clutch relaxes → lower traction. No parameter is tuned to an outcome: **E (or conc) is the input, the ECM
    bond stiffness is the DERIVED constitutive value, and traction is MEASURED.**

    Args:
        S: resting cell from :func:`build_resting_cell` (re-used across all E).
        material: :mod:`ff.ecm_library` registry key/alias (continuum: ``pa_gel``/``matrigel``/``hyaluronic_acid``;
            fibrillar: ``collagen_I``/``fibrin``).
        E_gel_Pa: continuum-gel Young's modulus [Pa] (the swept stiffness knob for continuum materials).
        conc: fibrillar concentration (mg/mL) — only for fibrillar materials (density → emergent modulus).
        steps: quasi-static implicit cell steps.
        dt: implicit cell timestep [s].
        kmc_every: steps between catch-slip KMC + nascent-rebind ticks (0 = auto → dt·kmc_every ≈ 1.5 s).
        f_myo: myosin per-minifilament stall force [pN] (loads the clutches).
        koff_xl: α-actinin crosslink off-rate [1/s].
        refresh_every / reshape_every / xl_turn_every: slow-mode / inextensibility / turnover cadences.
        adh_h: nascent-adhesion capture height above the ECM top [µm].
        capture_um: FA↔ECM-node capture radius [µm].
        ecm_substeps: explicit ECM relaxation substeps per cell step.
        target_z: fibrillar crosslink connectivity ⟨z⟩ target (ignored for continuum gels).
        seed: ECM-build RNG seed.
        device: ``"cpu"`` or ``"cuda:0"``.

    Returns:
        dict(material, E_pa, E_or_conc, is_fibrillar, traction_nN, bound_frac, force_per_clutch_pN,
             silhouette_um2, n_clutch, n_bound, n_ecm_nodes, mesh_um, connectivity_z, ecm_mean_disp_um,
             ecm_max_disp_um, kmc_every, dt, dt_cfl_s, ecm_dt_s, ecm_substeps, wall_s, diverged) plus the ECM
             deformation arrays (``_ecm_pos``/``_ecm_pos0``) for a viewer (stripped before JSON).
    """
    d = device
    cx = S["cx"]; net = cx.net; Nc = S["Nc"]; n_nuc = S["n_nuc"]; nuc = S["nuc"]; mem = S["mem"]
    Ne = Nc; N = Nc + n_nuc                                    # MT OFF ⇒ Ne==Nc; nucleus beads follow the cortex
    cp = resolve_clutch()
    if kmc_every <= 0:
        kmc_every = max(1, int(round(1.5 / dt)))              # dt·kmc_every ≈ 1.5 s → resolve the ~1 s clutch life

    # ---- cortex elastic topology (identical to sense) ----
    tri = np.ascontiguousarray(net.bend_triples, np.int32); nT = tri.shape[0]
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    foff = np.ascontiguousarray(net.fiber_offsets, np.int32)
    seg_per = np.diff(net.fiber_offsets) - 1
    soff = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    seg = float(net.seg_rest.mean())
    xl_ij_np = np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int64)
    kxl_np = np.ascontiguousarray(cx.xl_k, np.float64)
    xl_rest_np = np.ascontiguousarray(cx.xl_rest, np.float64)
    n_xl = xl_ij_np.shape[0]
    bend_triples_np = np.ascontiguousarray(net.bend_triples, np.int64)

    # ---- turgor volume reference + oriented faces from the RELAXED cortex (identical to sense) ----
    pos0 = S["pos_all"].copy()
    faces = _oriented_faces(pos0[:Nc])
    V0 = float(ConvexHull(pos0[:Nc]).volume); vmin = VMIN_FRAC * V0
    k_plane = 1.0e3                                            # substrate excluded-volume floor [pN/µm]

    # ---- BUILD THE ECM UNDER THE CELL (the stiffness lives HERE, not in a Winkler spring) ----
    cortex_pos = pos0[:Nc]
    basal = S["basal"]                                         # int32 basal actin-node indices
    z_sub = S["z_sub"]
    xy = cortex_pos[:, :2]
    lo = np.array([xy[:, 0].min() - 2.0, xy[:, 1].min() - 2.0, z_sub - 4.0])
    hi = np.array([xy[:, 0].max() + 2.0, xy[:, 1].max() + 2.0, z_sub + 0.5])
    spec = get_spec(material)
    ecm_rng = np.random.default_rng(seed)
    if spec.is_fibrillar:                                      # emergent modulus from the fiber microstructure
        ecm = build_ecm(material, lo, hi, pin_faces=("z_lo",), rng=ecm_rng, concentration=conc, target_z=target_z,
                        alignment_S=alignment_S, director=director)   # NEAR #6: aligned collagen (contact guidance)
        E_or_conc = float(conc) if conc is not None else float(spec.ref_conc)
    else:                                                      # continuum gel: bond k = analytic inverse of E_gel
        E_gel = spec.E_gel_Pa if E_gel_Pa is None else E_gel_Pa
        ecm = build_ecm(material, lo, hi, pin_faces=("z_lo",), rng=ecm_rng, E_gel_Pa=E_gel)
        E_or_conc = float(E_gel)
    ep0 = np.ascontiguousarray(ecm.net.pos, np.float64)
    n_em = ep0.shape[0]
    M = basal.size
    ecm_node = attach_clutches_to_ecm(cortex_pos[basal], ep0, capture_um=capture_um)   # int64[M]; −1 = no fiber

    # ---- ECM elastic topology + drag + CFL ----
    Etri = np.ascontiguousarray(ecm.net.bend_triples, np.int32); nT_e = Etri.shape[0]
    Eal = np.ascontiguousarray(_per_triple_alpha(ecm.net), np.float64) if nT_e else np.zeros(0)
    Eb = np.ascontiguousarray(ecm.links, np.int32); nB_e = Eb.shape[0]     # xl + seg (Delaunay) bonds stacked
    Ebk = np.ascontiguousarray(ecm.link_k, np.float64)
    Ebr = np.ascontiguousarray(ecm.link_rest, np.float64)
    Egam = np.where(ecm.pinned, 1.0e18, 1.0).astype(np.float64)            # pin the z_lo far face (immovable)
    seg_e = float(ecm.net.seg_rest.mean()) if (nT_e and ecm.net.seg_rest.size) else 1.0
    bend_term = (float(ecm.net.kappa.max()) / seg_e ** 3) if nT_e else 1.0  # continuum: κ≡0 → no bending mode
    Ekmax = max(bend_term, float(Ebk.max()) if Ebk.size else 1.0, cp.k_int)
    Edt = ECM_CFL_SAFETY / Ekmax                                           # explicit overdamped CFL (γ_free=1)
    Ensub = int(ecm_substeps)

    # ---- per-node drag + implicit-solver scalars (cell — identical to sense) ----
    gammas = physical_node_gammas(net, Nc, n_nuc)
    gamma_rep = float(np.median(gammas[:Nc]))
    kmax = max(float(net.kappa.max()) / seg**3, float(kxl_np.max()) if n_xl else 1.0,
               nuc.k_chrom + nuc.k_lamin, cp.k_int, k_plane)
    dt_cfl = crawl_cfl_dt(gammas, kmax, safety=0.1)

    # ---- device arrays: cell (identical to sense) ----
    pos_d = wp.array(pos0, dtype=wp.vec3d, device=d); f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d); alpha_d = wp.array(alpha, dtype=wp.float64, device=d)
    foff_d = wp.array(foff, dtype=wp.int32, device=d); soff_d = wp.array(soff, dtype=wp.int32, device=d)
    sr_d = wp.array(np.ascontiguousarray(net.seg_rest, np.float64), dtype=wp.float64, device=d)
    xl_d = wp.array(np.ascontiguousarray(np.stack([cx.xl_i, cx.xl_j], 1), np.int32), dtype=wp.int32, ndim=2, device=d)
    kxl_d = wp.array(kxl_np, dtype=wp.float64, device=d)
    r0_d = wp.array(xl_rest_np.copy(), dtype=wp.float64, device=d)
    faces_d = wp.array(faces, dtype=wp.int32, ndim=2, device=d)
    csum_d = wp.zeros(3, dtype=wp.float64, device=d); rsum_d = wp.zeros(1, dtype=wp.float64, device=d)
    has_myo = cx.myo_i.size > 0
    if has_myo:
        myo_d = wp.array(np.ascontiguousarray(np.stack([cx.myo_i, cx.myo_j], 1), np.int32),
                         dtype=wp.int32, ndim=2, device=d)
    ac_d = wp.array(basal.astype(np.int32), dtype=wp.int32, device=d)      # actin-end node per clutch (into cell)
    bd_d = wp.array((ecm_node >= 0).astype(np.int32), dtype=wp.int32, device=d)   # bound only where a fiber was found
    centre = wp.vec3d(float(S["c"][0]), float(S["c"][1]), float(S["c"][2]))
    dP_mem_area = 0.0
    frac_xl = 1.0 - np.exp(-koff_xl * dt * xl_turn_every)

    # ---- device arrays: ECM (its own position/force arrays; the clutch is the only bridge) ----
    Ep_d = wp.array(ep0, dtype=wp.vec3d, device=d); Ef_d = wp.zeros(n_em, dtype=wp.vec3d, device=d)
    Etri_d = wp.array(Etri, dtype=wp.int32, device=d) if nT_e else None
    Eal_d = wp.array(Eal, dtype=wp.float64, device=d) if nT_e else None
    Eb_d = wp.array(Eb, dtype=wp.int32, ndim=2, device=d) if nB_e else None
    Ebk_d = wp.array(Ebk, dtype=wp.float64, device=d) if nB_e else None
    Ebr_d = wp.array(Ebr, dtype=wp.float64, device=d) if nB_e else None
    Egam_d = wp.array(Egam, dtype=wp.float64, device=d)
    en_base_d = wp.array(ecm_node.astype(np.int32), dtype=wp.int32, device=d)   # per-clutch bound ECM node
    en_d = wp.zeros(M, dtype=wp.int32, device=d)              # en_base masked by bound state (only engaged grip)
    Edummy_e = wp.zeros(n_em, dtype=wp.vec3d, device=d)       # discarded ECM reaction in the cell force loop
    Edummy_c = wp.zeros(N, dtype=wp.vec3d, device=d)          # discarded cell reaction in the ECM substeps
    kmc_anch_d = wp.zeros(M, dtype=wp.vec3d, device=d)        # per-clutch live ECM-node pos for the KMC load

    def full_force(x_np):
        """Full FF cortex force at x — the implicit solver RHS (mirrors :func:`sense.full_force` but the clutch
        is the two-sided :func:`clutch_ecm_spring_kernel` pulling the actin toward its LIVE ECM node; the ECM
        reaction is discarded here and applied in the ECM substep loop instead). The osmotic/turgor force is the
        exact volume-gradient ΔP·g added host-side; the membrane inward tension stays a kernel."""
        pos_d.assign(np.ascontiguousarray(x_np, np.float64).reshape(N, 3))
        wp.launch(_zero, dim=N, inputs=[f_d], device=d)
        wp.launch(cytosim_bending_kernel, dim=nT, inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
        wp.launch(link_spring_kernel, dim=n_xl, inputs=[pos_d, xl_d, kxl_d, r0_d, f_d], device=d)
        if has_myo:
            wp.launch(myosin_kernel, dim=cx.myo_i.size, inputs=[pos_d, myo_d, wp.float64(f_myo), f_d], device=d)
        wp.launch(turgor_kernel, dim=Nc, inputs=[pos_d, centre, wp.float64(dP_mem_area), f_d], device=d)
        if n_nuc:
            wp.launch(nucleus_shell_kernel, dim=n_nuc, inputs=[pos_d, wp.int32(Ne), centre,
                      wp.float64(nuc.R_nuc_um), wp.float64(nuc.k_chrom), wp.float64(nuc.k_lamin),
                      wp.float64(nuc.d_knee_um), wp.float64(nuc.F_knee_pN), f_d], device=d)
        wp.launch(substrate_plane_kernel, dim=Nc, inputs=[pos_d, wp.float64(z_sub), wp.float64(k_plane), f_d], device=d)
        wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),
                  wp.float64(cp.rest_um), f_d, Edummy_e], device=d)      # +f on the cell (traction), reaction discarded
        wp.synchronize_device(d)
        F = f_d.numpy().reshape(-1)
        pcx = np.ascontiguousarray(x_np, np.float64).reshape(N, 3)[:Nc]; cen = pcx.mean(0)
        aa = pcx[faces[:, 0]] - cen; bb = pcx[faces[:, 1]] - cen; ccf = pcx[faces[:, 2]] - cen
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        dPv = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-12 * V0) - (TURGOR_PI_IN0 - TURGOR_DP0)
        dPv = min(max(dPv, -TURGOR_PI_IN0), TURGOR_PI_IN0)
        F[:3 * Nc] += (dPv * volume_gradient(pcx, faces, cen)).reshape(-1)
        return F

    t0 = time.time(); diverged = False
    for step in range(steps):
        wp.launch(mask_ecm_by_bound_kernel, dim=M, inputs=[en_base_d, bd_d, en_d], device=d)  # only ENGAGED clutches grip
        if step % refresh_every == 0:                          # slow modes via device reductions (centroid, R_mean)
            csum_d.zero_(); wp.launch(sum_pos_kernel, dim=Nc, inputs=[pos_d, csum_d], device=d)
            cc = csum_d.numpy() / Nc
            centre = wp.vec3d(float(cc[0]), float(cc[1]), float(cc[2]))
            rsum_d.zero_(); wp.launch(sum_radius_kernel, dim=Nc, inputs=[pos_d, centre, rsum_d], device=d)
            Rm = float(rsum_d.numpy()[0]) / Nc
            area = 4.0 * np.pi * Rm**2
            dP_mem_area = -2.0 * min(mem.gamma_mem, mem.tau_lysis) / max(Rm, 1e-9) * area / Nc
            if not np.isfinite(Rm) or Rm > 50.0 * S["R"]:
                print(f"    [!] divergence at step {step} (Rm={Rm:.2f}) — truncating"); diverged = True; break
        # ---- implicit host step (unconditionally stable; large dt) — identical machinery to sense ----
        xv = pos_d.numpy().reshape(-1)
        pcx2 = xv.reshape(N, 3)[:Nc]; cen2 = pcx2.mean(0)
        aa = pcx2[faces[:, 0]] - cen2; bb = pcx2[faces[:, 1]] - cen2; ccf = pcx2[faces[:, 2]] - cen2
        Vc = abs(float((aa * np.cross(bb, ccf)).sum() / 6.0))
        if not np.isfinite(Vc) or Vc <= 1.02 * vmin:
            print(f"    [!] volume collapse at step {step} (Vc={Vc:.1f} ≤ vmin={vmin:.1f}) — truncating")
            diverged = True; break
        gN = volume_gradient(pcx2, faces, cen2)
        k_vol = TURGOR_PI_IN0 * (V0 - vmin) / max(Vc - vmin, 1e-3 * V0) ** 2
        vol_g = np.zeros(3 * N); vol_g[:3 * Nc] = gN.reshape(-1)
        _de = np.zeros(3 * N)                                  # clutch + substrate-plane stiffness → implicit K diagonal
        _bd = bd_d.numpy(); _bb = basal[_bd > 0]
        _de[3 * _bb] = cp.k_int; _de[3 * _bb + 1] = cp.k_int; _de[3 * _bb + 2] = cp.k_int
        _de[3 * basal + 2] += k_plane
        xv, _ = implicit_step_current(xv, full_force, bend_triples_np, alpha, xl_ij_np, kxl_np,
                                      gamma=gamma_rep, dt=dt, n_newton=1, vol_g=vol_g, k_vol=k_vol, diag_extra=_de)
        # HARD incompressibility: project the cortex back to V=V0
        xr = xv.reshape(N, 3); pcx3 = xr[:Nc].copy(); cen3 = pcx3.mean(0)
        for _ in range(3):
            a3 = pcx3[faces[:, 0]] - cen3; b3 = pcx3[faces[:, 1]] - cen3; c3 = pcx3[faces[:, 2]] - cen3
            Vp = abs(float((a3 * np.cross(b3, c3)).sum() / 6.0))
            gp = volume_gradient(pcx3, faces, cen3); den = float((gp * gp).sum())
            if den > 1e-9:
                pcx3 = pcx3 - ((Vp - V0) / den) * gp
        xr[:Nc] = pcx3; xv = xr.reshape(-1)
        pos_d.assign(np.ascontiguousarray(xv, np.float64).reshape(N, 3))
        if step % reshape_every == 0:                          # NF2007 §5.3 inextensibility
            wp.launch(reshape_kernel, dim=foff.shape[0] - 1, inputs=[pos_d, foff_d, soff_d, sr_d, wp.int32(2)], device=d)
        if step % xl_turn_every == 0 and step > 0:             # α-actinin crosslink turnover (fluidise → stress relax)
            wp.launch(xl_turnover_kernel, dim=n_xl, inputs=[pos_d, xl_d, r0_d, wp.float64(frac_xl)], device=d)
        if step % kmc_every == 0 and step > 0:                 # catch-slip rupture (against the LIVE ECM node) + rebind
            ep = Ep_d.numpy(); en_host = en_base_d.numpy()
            anch = np.full((M, 3), 1.0e6, np.float64)          # far sentinel for node-less clutches → no spurious load
            hn = en_host >= 0
            if hn.any():
                anch[hn] = ep[en_host[hn]]
            kmc_anch_d.assign(np.ascontiguousarray(anch, np.float64))
            wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, kmc_anch_d, bd_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                      wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(0.0),
                      wp.float64(dt * kmc_every), wp.int32(step)], device=d)     # k_on=0 → host does the rebind
            p = pos_d.numpy(); bd = bd_d.numpy()
            can = (bd == 0) & (p[basal, 2] < z_sub + adh_h)    # detached + near the ECM top → can re-adhere
            p_on = 1.0 - np.exp(-dt * kmc_every * cp.k_on)
            reb = can & (np.random.default_rng(step).random(M) < p_on)
            if reb.any():                                      # nascent adhesion: grab the NEAREST ECM node underneath
                newn = attach_clutches_to_ecm(p[basal[reb]], ep, capture_um=capture_um)
                ok = newn >= 0
                if ok.any():
                    idx = np.where(reb)[0][ok]
                    en_host[idx] = newn[ok]; bd[idx] = 1
                    en_base_d.assign(np.ascontiguousarray(en_host.astype(np.int32)))
                    bd_d.assign(np.ascontiguousarray(bd.astype(np.int32)))
        # ---- LIVE ECM relaxes under the clutch reactions + its own elasticity (explicit overdamped substeps) ----
        wp.launch(mask_ecm_by_bound_kernel, dim=M, inputs=[en_base_d, bd_d, en_d], device=d)   # post-KMC bound state
        for _ in range(Ensub):
            wp.launch(_zero, dim=n_em, inputs=[Ef_d], device=d)
            if nT_e:
                wp.launch(cytosim_bending_kernel, dim=nT_e, inputs=[Ep_d, Etri_d, Eal_d, Ef_d], device=d)
            if nB_e:
                wp.launch(link_spring_kernel, dim=nB_e, inputs=[Ep_d, Eb_d, Ebk_d, Ebr_d, Ef_d], device=d)
            wp.launch(clutch_ecm_spring_kernel, dim=M, inputs=[pos_d, ac_d, Ep_d, en_d, wp.float64(cp.k_int),
                      wp.float64(cp.rest_um), Edummy_c, Ef_d], device=d)         # −f on the fiber (matrix remodel)
            wp.launch(axpy_physical_kernel, dim=n_em, inputs=[Ep_d, wp.float64(Edt), Egam_d, Ef_d], device=d)

    wp.synchronize_device(d)
    xp = pos_d.numpy(); bd = bd_d.numpy(); ep = Ep_d.numpy(); en_host = en_base_d.numpy()
    eng = (en_host >= 0) & (bd > 0)                            # engaged clutches: bound AND holding a fiber node
    if eng.any():
        Le = np.linalg.norm(xp[basal][eng] - ep[en_host[eng]], axis=1)
        Fclutch = cp.k_int * np.maximum(Le - cp.rest_um, 0.0)  # engaged-clutch extension force [pN]
    else:
        Fclutch = np.zeros(0)
    n_bound = int(bd.sum())
    sil = float(ConvexHull(xp[:Nc][:, :2]).volume)             # top-down xy silhouette area [µm²]
    ecm_disp = np.linalg.norm(ep - ep0, axis=1)                # ECM node displacement (the matrix deformation)

    # ---- NEAR #6 contact guidance: director-decomposition of the engaged-clutch traction + ECM virial ratio ----
    dvec = np.asarray(ecm.director, float)
    aniso = traction_director_anisotropy(
        xp[basal][eng] if eng.any() else np.zeros((0, 3)),
        ep[en_host[eng]] if eng.any() else np.zeros((0, 3)),
        cp.k_int, cp.rest_um, dvec)                            # A_F=F∥/F⊥ (Ray-2017 force anisotropy)
    dhat = _inplane_unit(dvec); ehat = np.array([-dhat[1], dhat[0], 0.0])
    try:                                                       # grid-invariant cross-check: virial σ on the director
        sig = ecm_material_stress(ecm, ep)                     # σ[3,3] [Pa] (clutch count-independent)
        sig_par = float(dhat @ sig @ dhat); sig_perp = float(ehat @ sig @ ehat)
        R_sigma = float(sig_par / sig_perp) if abs(sig_perp) > 1e-30 else float("nan")
    except Exception:                                          # virial degenerate (e.g. no bonds) → nan, not a crash
        sig_par = sig_perp = R_sigma = float("nan")
    return dict(material=str(material), E_pa=float(E_or_conc), E_or_conc=float(E_or_conc),
                alignment_S=float(alignment_S), S_measured=float(getattr(ecm, "S_measured", 0.0)),
                director=[float(x) for x in dvec],
                A_F=float(aniso["A_F"]), A_F_rms=float(aniso["A_F_rms"]),
                F_par_pN=float(aniso["F_par_pN"]), F_perp_pN=float(aniso["F_perp_pN"]),
                S_bound=float(aniso["S_bound"]), n_engaged=int(aniso["n_engaged"]),
                R_sigma=R_sigma, sigma_par_Pa=sig_par, sigma_perp_Pa=sig_perp,
                is_fibrillar=bool(spec.is_fibrillar), traction_nN=float(Fclutch.sum() / 1e3),
                bound_frac=float(bd.mean()), force_per_clutch_pN=float(Fclutch.mean()) if Fclutch.size else 0.0,
                silhouette_um2=sil, n_clutch=int(M), n_bound=n_bound, n_ecm_nodes=int(n_em),
                mesh_um=float(ecm.mesh_size_um), connectivity_z=float(ecm.connectivity_z),
                ecm_mean_disp_um=float(ecm_disp.mean()), ecm_max_disp_um=float(ecm_disp.max()),
                kmc_every=int(kmc_every), dt=float(dt), dt_cfl_s=float(dt_cfl), ecm_dt_s=float(Edt),
                ecm_substeps=int(Ensub), wall_s=float(time.time() - t0), diverged=bool(diverged),
                _ecm_pos=ep, _ecm_pos0=ep0)


def plot(rows: list[dict], meta: dict, path: str) -> None:
    """Traction (nN) + bound-fraction vs substrate E (log-x) with the biphasic optimum marked → PNG.

    Twin-axis: left = engaged-clutch traction, right = bound fraction. The Bangasser-Odde optimum (E at max
    traction) is marked with a vertical line + annotation. Units annotated; no axis truncation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    E = np.array([r["E_pa"] for r in rows]); tr = np.array([r["traction_nN"] for r in rows])
    bf = np.array([r["bound_frac"] for r in rows])
    opt = int(np.argmax(tr))
    fig, ax1 = plt.subplots(figsize=(7.4, 5.0))
    c_tr, c_bf = "#1f77b4", "#d62728"
    ax1.plot(E, tr, "-o", color=c_tr, lw=2, ms=7, label="traction")
    ax1.set_xscale("log")
    ecm_mode = meta.get("mode") == "ecm-network"
    ax1.set_xlabel("ECM network Young's modulus  E_gel  [Pa]" if ecm_mode
                   else "substrate Young's modulus  E  [Pa]")
    ax1.set_ylabel("engaged-clutch traction  [nN]", color=c_tr)
    ax1.tick_params(axis="y", labelcolor=c_tr)
    ax1.set_ylim(bottom=0.0)
    ax2 = ax1.twinx()
    ax2.plot(E, bf, "--s", color=c_bf, lw=1.8, ms=6, label="bound fraction")
    ax2.set_ylabel("bound fraction  [—]", color=c_bf)
    ax2.tick_params(axis="y", labelcolor=c_bf)
    ax2.set_ylim(0.0, 1.05)
    ax1.axvline(E[opt], color="0.4", ls=":", lw=1.5)
    ax1.annotate(f"biphasic optimum\nE*={E[opt]:.0f} Pa\n({_pa_label(E[opt])[1].split('(')[0].strip()})",
                 xy=(E[opt], tr[opt]), xytext=(0.02, 0.82), textcoords="axes fraction",
                 fontsize=9, color="0.2", arrowprops=dict(arrowstyle="->", color="0.4"))
    note = "NATIVE" if meta.get("nf", 0) >= NATIVE_NF else f"COARSE nf={meta.get('nf')} — DEV SMOKE (non-authoritative)"
    sub = (f"live ecm_library '{meta.get('ecm_material')}' network under the cell (two-sided clutch↔ECM)  ·  "
           f"Nc={meta.get('Nc')}  ·  tag={meta.get('tag')}" if meta.get("mode") == "ecm-network" else
           f"Winkler k_sub=2Ea/(1−ν²), a={meta.get('a_um')} µm, ν={meta.get('nu')}  ·  "
           f"Nc={meta.get('Nc')}  ·  tag={meta.get('tag')}")
    ax1.set_title(f"FF substrate-stiffness sensing  ·  motor-clutch biphasic  ·  {note}\n{sub}", fontsize=9)
    lines = ax1.get_lines()[:1] + ax2.get_lines()[:1]
    ax1.legend(lines, [ln.get_label() for ln in lines], loc="upper right", fontsize=9, framealpha=0.9)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--E-list", default=",".join(f"{e:g}" for e in DEFAULT_E_LIST),
                    help="comma-separated substrate Young's moduli [Pa] (default = the PA-gel library ladder "
                         "150 brain → 40000 stiff-osteoid)")
    ap.add_argument("--nf", type=int, default=2000,
                    help="cortex filament count. Default 2000 = a COARSE CPU smoke — NON-AUTHORITATIVE; the "
                         f"native cell (--nf {NATIVE_NF}) on the A5000 is the only valid basis for a conclusion.")
    ap.add_argument("--steps", type=int, default=1500, help="quasi-static implicit steps per E")
    ap.add_argument("--relax-steps", type=int, default=4000, help="pre-relaxation steps to the resting set-point")
    ap.add_argument("--n-fa", type=int, default=200, help="discrete focal-adhesion sites (0 = one per basal node)")
    ap.add_argument("--dt-impl", type=float, default=0.05, help="implicit timestep [s]")
    ap.add_argument("--kmc-every", type=int, default=0, help="steps between catch-slip KMC ticks (0 = auto → dt·kmc≈1.5 s)")
    ap.add_argument("--a-um", type=float, default=0.05, help="adhesion-patch radius [µm] (surfaced magnitude lever)")
    ap.add_argument("--nu", type=float, default=0.45, help="substrate Poisson ratio")
    ap.add_argument("--seed", type=int, default=1, help="cortex/crosslinker RNG seed")
    ap.add_argument("--device", default="cpu", help="cpu or cuda:0")
    ap.add_argument("--out", default="ffn_sim/outputs/ff/ecm_lib", help="output directory")
    ap.add_argument("--tag", default="sensing")
    ap.add_argument("--mode", choices=("winkler", "ecm-network"), default="winkler",
                    help="winkler (PATH-b, default): stiffness via k_sub=2Ea/(1−ν²) on Winkler FA anchors. "
                         "ecm-network (PATH-a): grip a LIVE ff.ecm_library network built under the cell — the "
                         "substrate stiffness is the ECM's OWN calibrated modulus, traction two-sided-emergent.")
    ap.add_argument("--ecm-material", default="pa_gel",
                    help="ecm-network mode: ff.ecm_library registry key. Continuum gels (pa_gel/matrigel/"
                         "hyaluronic_acid) sweep --E-list as E_gel; fibrillar (collagen_I/fibrin) use --ecm-conc.")
    ap.add_argument("--ecm-conc", type=float, default=None,
                    help="ecm-network mode: fibrillar concentration (mg/mL) for a fibrillar --ecm-material "
                         "(default = the material's reference concentration). Ignored for continuum gels.")
    ap.add_argument("--ecm-substeps", type=int, default=20, help="ecm-network mode: explicit ECM relax substeps per cell step")
    args = ap.parse_args()

    E_list = [float(x) for x in args.E_list.split(",") if x.strip()]
    wp.init(); t0 = time.time()
    coarse = args.nf < NATIVE_NF
    if coarse:
        print(f"[!] COARSE nf={args.nf} < native {NATIVE_NF} — this is a DEV SMOKE, NON-AUTHORITATIVE. "
              f"A coarse cortex lacks structural stability; reconfirm every finding at NATIVE --nf {NATIVE_NF} "
              f"--device cuda:0 on the A5000 before reporting (CLAUDE.md HARD rule).")
    print(f"[build] resting cell: nf={args.nf}, relax_steps={args.relax_steps}, n_fa={args.n_fa} ...")
    S = build_resting_cell(args.nf, relax_steps=args.relax_steps, n_fa=args.n_fa, seed=args.seed, device=args.device)
    rr = S["resting"]
    print(f"[from-resting] Nc={S['Nc']} nucleus={S['n_nuc']} basal-FA={S['basal'].size}  R={S['R']:.2f} µm  "
          f"ΔP={rr['dP_Pa']:.1f} Pa  γ={rr['gamma_mN_m']:.3f} mN/m  V/V0={rr['V_over_V0']:.3f}  ({time.time()-t0:.0f}s)")

    os.makedirs(args.out, exist_ok=True)
    os.makedirs(f"{args.out}/figs", exist_ok=True)

    if args.mode == "winkler":                                     # ── PATH (b): Winkler compliant substrate ──
        rows = []
        for E in E_list:
            r = sense(S, E, steps=args.steps, dt=args.dt_impl, kmc_every=args.kmc_every, a_um=args.a_um,
                      nu=args.nu, device=args.device)
            recipe, tissue = _pa_label(E)
            r["pa_recipe"] = recipe; r["pa_tissue"] = tissue
            rows.append(r)
            flag = "  [DIVERGED]" if r["diverged"] else ""
            print(f"[E={E:>7.0f} Pa | {recipe:>8s} {tissue.split('(')[0].strip():<22.22s}] "
                  f"k_sub={r['k_sub']:.1f} pN/µm  traction={r['traction_nN']:.3f} nN  "
                  f"bound={r['bound_frac']:.2f}  f/clutch={r['force_per_clutch_pN']:.2f} pN  "
                  f"silhouette={r['silhouette_um2']:.1f} µm²  ({r['wall_s']:.0f}s){flag}")

        opt = int(np.argmax([r["traction_nN"] for r in rows]))
        E_opt = rows[opt]["E_pa"]
        interior = 0 < opt < len(rows) - 1
        verdict = (f"BIPHASIC optimum at E*={E_opt:.0f} Pa ({rows[opt]['pa_tissue'].split('(')[0].strip()})"
                   if interior else
                   f"traction is MONOTONIC over this E range (peak at the {'stiff' if opt else 'soft'} endpoint "
                   f"E={E_opt:.0f} Pa) — widen --E-list or --steps to resolve an interior optimum")
        print(f"[BIPHASIC] {verdict}")

        meta = dict(tag=args.tag, nf=args.nf, Nc=S["Nc"], n_nuc=S["n_nuc"], n_fa=int(S["basal"].size), R_um=S["R"],
                    a_um=args.a_um, nu=args.nu, dt_impl=args.dt_impl, steps=args.steps, relax_steps=args.relax_steps,
                    seed=args.seed, device=args.device, resting=rr, E_list_pa=E_list, E_opt_pa=E_opt,
                    biphasic_interior=bool(interior), coarse_nonauthoritative=bool(coarse),
                    native_nf=NATIVE_NF, path_note="PATH-b Winkler compliant substrate; PATH-a (grip the "
                    "ecm_library fiber network) is a documented future increment",
                    a_note="absolute traction carries the adhesion-patch radius a (KB-PIV-5, surfaced); the "
                    "RELATIVE optimum is a-invariant")
        json.dump({"meta": meta, "rows": rows}, open(f"{args.out}/ecm_stiffness_sensing.json", "w"), indent=2)
        plot(rows, meta, f"{args.out}/figs/stiffness_sensing.png")
        print(f"wrote {args.out}/ecm_stiffness_sensing.json + figs/stiffness_sensing.png  (total {time.time()-t0:.0f}s)")
        return

    # ── PATH (a): grip a LIVE ff.ecm_library ECM network built under the cell ──
    spec = get_spec(args.ecm_material)
    print(f"[ecm-network] material='{args.ecm_material}' ({'fibrillar' if spec.is_fibrillar else 'continuum gel'}) — "
          f"stiffness is the ECM's OWN calibrated modulus; traction two-sided-emergent (clutch↔ECM).")
    sweep = [None] if spec.is_fibrillar else list(E_list)          # continuum sweeps E_gel; fibrillar = 1 conc
    rows = []
    for E in sweep:
        r = sense_ecm_network(S, args.ecm_material, E_gel_Pa=E, conc=args.ecm_conc, steps=args.steps,
                              dt=args.dt_impl, kmc_every=args.kmc_every, ecm_substeps=args.ecm_substeps,
                              seed=args.seed + 6, device=args.device)
        ep = r.pop("_ecm_pos"); ep0 = r.pop("_ecm_pos0")          # strip big arrays out of the JSON row
        if spec.is_fibrillar:
            r["pa_recipe"], r["pa_tissue"] = "—", f"{spec.name} @ {r['E_or_conc']:g} mg/mL"
            label = f"conc={r['E_or_conc']:>6.2f} mg/mL"
        else:
            recipe, tissue = _pa_label(r["E_or_conc"]); r["pa_recipe"], r["pa_tissue"] = recipe, tissue
            label = f"E={r['E_or_conc']:>7.0f} Pa | {recipe:>8s} {tissue.split('(')[0].strip():<20.20s}"
        rows.append(r)
        tag = f"{args.ecm_material}_{('conc%g' % r['E_or_conc']) if spec.is_fibrillar else ('E%g' % r['E_or_conc'])}"
        np.savez_compressed(f"{args.out}/ecm_deform_{tag}.npz", ecm_pos0=ep0.astype(np.float32),
                            ecm_pos=ep.astype(np.float32))
        flag = "  [DIVERGED]" if r["diverged"] else ""
        print(f"[{label}] traction={r['traction_nN']:.3f} nN  bound={r['bound_frac']:.2f} "
              f"({r['n_bound']}/{r['n_clutch']})  f/clutch={r['force_per_clutch_pN']:.2f} pN  "
              f"ECM-nodes={r['n_ecm_nodes']} mesh={r['mesh_um']:.2f}µm ⟨z⟩={r['connectivity_z']:.1f}  "
              f"ECM-disp⟨{r['ecm_mean_disp_um']*1e3:.0f}⟩/max{r['ecm_max_disp_um']*1e3:.0f}nm  "
              f"silhouette={r['silhouette_um2']:.1f} µm²  ({r['wall_s']:.0f}s){flag}")

    opt = int(np.argmax([r["traction_nN"] for r in rows]))
    E_opt = rows[opt]["E_or_conc"]
    if len(rows) >= 2 and not spec.is_fibrillar:
        stiff_hi = rows[-1]["traction_nN"] >= rows[0]["traction_nN"]
        interior = 0 < opt < len(rows) - 1
        verdict = (f"BIPHASIC optimum at E*={E_opt:.0f} Pa" if interior else
                   f"traction is MONOTONIC — {'STIFFER ECM → MORE traction (durotaxis-consistent)' if stiff_hi else 'softer ECM → more traction'} "
                   f"(peak at E={E_opt:.0f} Pa); widen --E-list/--steps to resolve an interior optimum")
        print(f"[ECM-SENSING] {verdict}")
    else:
        interior = False
        print(f"[ECM-SENSING] single {args.ecm_material} network — traction={rows[0]['traction_nN']:.3f} nN")

    meta = dict(tag=args.tag, mode="ecm-network", ecm_material=args.ecm_material, is_fibrillar=bool(spec.is_fibrillar),
                nf=args.nf, Nc=S["Nc"], n_nuc=S["n_nuc"], n_fa=int(S["basal"].size), R_um=S["R"],
                dt_impl=args.dt_impl, steps=args.steps, relax_steps=args.relax_steps, ecm_substeps=args.ecm_substeps,
                seed=args.seed, device=args.device, resting=rr, E_list_pa=E_list, E_opt_pa=E_opt,
                biphasic_interior=bool(interior), coarse_nonauthoritative=bool(coarse), native_nf=NATIVE_NF,
                path_note="PATH-a: the cell grips a LIVE ff.ecm_library network built under it; the substrate "
                "stiffness is the ECM's OWN calibrated modulus (continuum: k_bond=inverse of E_gel; fibrillar: "
                "emergent from microstructure). Traction = two-sided clutch↔ECM coupling + ECM deformation.",
                a_note="absolute traction carries the clutch/adhesion patch magnitude; the RELATIVE stiffness "
                "ranking (which ECM sustains more traction) is the claim.")
    fname = f"{args.out}/ecm_stiffness_ecmnet.json"
    json.dump({"meta": meta, "rows": rows}, open(fname, "w"), indent=2)
    if len(rows) >= 2 and not spec.is_fibrillar:
        plot(rows, meta, f"{args.out}/figs/stiffness_sensing_ecmnet.png")
        print(f"wrote {fname} + figs/stiffness_sensing_ecmnet.png + ecm_deform_*.npz  (total {time.time()-t0:.0f}s)")
    else:
        print(f"wrote {fname} + ecm_deform_*.npz  (total {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
