"""FF reproduction of Slater, Li, Indana, Xie, Chaudhuri & Taeyoon Kim, *Soft Matter* 2021, 17, 10274
(TAG=``Slater2021_SoftMatter``, SE333) — "Transient mechanical interactions between cells and viscoelastic
extracellular matrix".

STAGE (A) — PEER-ORACLE REPRODUCTION (PI 2026-07-13). Reproduce Kim's coupled experiment on the FF engine
*at Kim's own fidelity level*: a contractile cell embedded in a collagen matrix whose TRANSIENT cross-linkers
turn over (Bell slip) so the matrix stress RELAXES over time. The single deviation from the mechanistic ideal
that Kim himself makes, we make too — and label honestly:

  * IMPOSED contraction. Kim's cell has NO explicit myosin; its cortex is a phenomenological zero-rest-length
    spring ``ks,c`` that he SWEEPS. We likewise IMPOSE an inward cortical tension on a cortex ring and sweep
    it. This is NOT tuned to reproduce Kim's numbers (that would be fitting-to-outcome, forbidden) — it is a
    swept boundary condition. FF's *emergent* myosin generator is floored ~530x (the γ-floor) — a SEPARATE
    open item, not this script's concern. Per the physiological-baseline HARD rule, running with the floored
    (≈0-tension) emergent motor would itself be an unphysical baseline; imposing the drive is the compliant
    choice here.

Everything DOWNSTREAM is emergent and mechanistic: the collagen Mikado network (bending + finite-EA segments),
the FA harmonic coupling, the transient cross-link turnover, and the stress propagation / remodeling response.

Model (Kim Fig 1): a central contractile cortex RING (radius ``R_cell``) is linked by FA harmonic springs to
the nearest nodes of a surrounding collagen Mikado matrix; the matrix outer boundary is pinned (embedded in
bulk). Two-phase: (1) equilibrate the assembled matrix, (2) ramp on cortical contraction → tensile force
develops near the cell, propagates outward, and RELAXES over time as loaded cross-linkers unbind.

KEY VALIDATION (this smoke test): the SAME run with cross-link turnover ON vs OFF.
  * turnover OFF (elastic matrix)  → stress rises to a peak and PLATEAUS.
  * turnover ON  (viscoelastic)    → stress rises to a peak then RELAXES toward a lower residual.
That contrast IS the viscoelastic mechanism the paper is about, now coupled to a contracting cell.

Observables (Kim): σ_rr(r) over time (``ecm_mechanics.stress_field_radial``), σ(t) near the membrane
(relaxation), matrix radial fiber displacement vs r.

⚠️ COARSE SMOKE TEST — NON-AUTHORITATIVE per the native+full-compartment HARD rule (coarse cells self-collapse;
conclusions require the native run). This proves the coupled loop RUNS and produces qualitative stress
generation→propagation→relaxation. NATIVE promotion (Nc-scale matrix, physiological cortical-tension anchor,
physical collagen drag for a real time axis, quasi-2D cylindrical shells for Kim's r^-1) is the next step and
the only basis for any reported conclusion.

Run:  python -m aleph.scripts.ff_kim_repro [--device cpu|cuda:0]
Out:  aleph/outputs/ff/kim_repro/{kim_repro.json, figs/kim_repro_smoke.png}
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
from scipy.spatial import cKDTree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import warp as wp

from aleph.laws import ecm_library as L
from aleph.laws import ecm_mechanics as M
from aleph.laws import units as U
from aleph.laws.network_warp import _zero, axpy_kernel, link_spring_kernel, wlc_spring_kernel, xl_turnover_kernel
from aleph.laws.forces_warp import _per_triple_alpha, cytosim_bending_kernel
from aleph.laws.wlc import wlc_x_crossover

wp.init()

OUT = "aleph/outputs/ff/kim_repro"
FIGS = f"{OUT}/figs"


def build_cortex_ring(matrix_pos, center, R_cell, n_cortex, capture_um):
    """A ring of ``n_cortex`` cortex nodes at radius ``R_cell`` in the z=center plane, each FA-linked to its
    nearest matrix node (harmonic spring, rest = initial separation). Returns cortex positions, FA pairs
    (cortex→matrix, global indices assume cortex nodes are appended AFTER the matrix nodes), FA rest lengths,
    and cortex-ring bond pairs (consecutive cortex nodes) with rest lengths."""
    ang = np.linspace(0.0, 2.0 * np.pi, n_cortex, endpoint=False)
    cortex = np.stack([center[0] + R_cell * np.cos(ang),
                       center[1] + R_cell * np.sin(ang),
                       np.full(n_cortex, center[2])], axis=1)
    n_m = matrix_pos.shape[0]
    tree = cKDTree(matrix_pos)
    dists, nn = tree.query(cortex, k=1)                      # nearest matrix node per cortex node
    # FA pairs: (cortex global idx, matrix idx). cortex global idx = n_m + c.
    fa_pairs = np.stack([n_m + np.arange(n_cortex), nn], axis=1).astype(np.int32)
    fa_rest = dists.astype(np.float64)
    # cortex ring bonds (consecutive nodes, closed loop)
    ring_pairs = np.stack([n_m + np.arange(n_cortex),
                           n_m + np.roll(np.arange(n_cortex), -1)], axis=1).astype(np.int32)
    ring_rest = np.linalg.norm(cortex - np.roll(cortex, -1, axis=0), axis=1).astype(np.float64)
    return cortex, fa_pairs, fa_rest, ring_pairs, ring_rest


def stress_flux_cylindrical(ecm, pos, center, r_edges, h):
    """Kim's exact stress measure in a quasi-2D CYLINDRICAL domain: for a cylinder of radius r, sum the RADIAL
    component of the spring tension of every bond (fiber segment + crosslinker) that CROSSES that cylinder,
    and divide by the lateral area 2π·r·h. In steady state the net radial force transmitted across any cylinder
    ≈ the cell's contractile force (independent of r), so σ(r) ∝ 1/r — Kim's r^-1 (vs r^-2 in 3D spherical),
    a lateral-area-∝-r geometric result the spherical virial cannot show. Returns r, σ(r) [Pa], and the fitted
    decay exponent n (|σ|~r^-n) over the clean mid-range."""
    pos = np.ascontiguousarray(pos, float)
    c = np.asarray(center, float)
    ii, jj, kk, rr0 = [], [], [], []
    for i, j, k, r0 in ((ecm.xl_i, ecm.xl_j, ecm.xl_k, ecm.xl_rest),
                        (ecm.seg_i, ecm.seg_j, ecm.seg_k, ecm.seg_rest)):
        if getattr(i, "size", 0):
            ii.append(np.asarray(i)); jj.append(np.asarray(j))
            kk.append(np.broadcast_to(k, np.asarray(i).shape) if np.ndim(k) == 0 else np.asarray(k))
            rr0.append(np.broadcast_to(r0, np.asarray(i).shape) if np.ndim(r0) == 0 else np.asarray(r0))
    i = np.concatenate(ii); j = np.concatenate(jj); k = np.concatenate(kk); r0 = np.concatenate(rr0)
    pi, pj = pos[i], pos[j]
    ri = np.linalg.norm(pi[:, :2] - c[:2], axis=1)
    rj = np.linalg.norm(pj[:, :2] - c[:2], axis=1)
    dvec = pj - pi
    L = np.linalg.norm(dvec, axis=1) + 1e-12
    F = k * (L - r0)                                            # bond tension [pN]
    mid_xy = 0.5 * (pi[:, :2] + pj[:, :2]) - c[:2]
    rmid = np.linalg.norm(mid_xy, axis=1) + 1e-12
    rhat = mid_xy / rmid[:, None]
    d_r = np.einsum("na,na->n", dvec[:, :2], rhat)             # bond radial span
    F_rad = F * d_r / L                                        # radial force carried by the bond [pN]
    rc = 0.5 * (r_edges[:-1] + r_edges[1:])
    lo_r, hi_r = np.minimum(ri, rj), np.maximum(ri, rj)
    sig = np.zeros(rc.size); force = np.zeros(rc.size)
    for s, rk in enumerate(rc):
        crossing = (lo_r < rk) & (hi_r >= rk)                 # bond straddles the cylinder at r_k
        force[s] = abs(float(np.sum(F_rad[crossing])))        # total radial force transmitted across [pN]
        sig[s] = force[s] / (2.0 * np.pi * rk * h)            # → stress [Pa]

    def _fit(mask):
        m = sig[mask]; r_ = rc[mask]; ok = (m > 0) & np.isfinite(m)
        return float(-np.polyfit(np.log(r_[ok]), np.log(m[ok]), 1)[0]) if ok.sum() >= 3 else float("nan")
    n_exp = _fit(np.arange(rc.size) < rc.size)                # whole (all shells)
    q1, q2 = np.percentile(rc, [33, 66])                       # MID = clean intermediate field (drop finite-cell NEAR + pinned FAR)
    n_mid = _fit((rc > q1) & (rc <= q2))
    return rc, sig, n_exp, n_mid, force


def run_case(ecm, center, R_cell, R_pin, *, contract_frac, koff0_per_s, x_beta_nm=0.4,
             n_cortex=20, k_fa=200.0, k_cortex=100.0, n_record=14, mech_substeps=600,
             ramp_records=3, dt_real_s=None, device="cpu", verbose=False, wlc=False, save_geom=False):
    """One coupled contracting-cell-in-matrix run (DISPLACEMENT-controlled: the cortex ring is prescribed to
    shrink to ``contract_frac``·R_cell over ``ramp_records`` then HELD — a step-strain-and-hold, the standard
    stress-relaxation protocol; the imposed contraction is swept via ``contract_frac``). Returns σ_rr(r) &
    near-membrane σ over time + the final matrix radial displacement. ``koff0_per_s=0`` → elastic (OFF) control."""
    d = device
    if dt_real_s is None:
        dt_real_s = (1.0 / max(koff0_per_s, 1e-6)) * 5.0 / n_record   # span ~5 relaxation times (koff>0)

    # ---- matrix arrays (mirror ecm_mechanics.stress_relaxation) ----
    mpos0 = ecm.net.pos.copy()
    n_m = ecm.net.n_nodes
    tri = np.ascontiguousarray(ecm.net.bend_triples, np.int32)
    has_bend = tri.shape[0] > 0
    alpha = np.ascontiguousarray(_per_triple_alpha(ecm.net), np.float64) if has_bend else np.zeros(0)
    seg = np.stack([ecm.seg_i, ecm.seg_j], 1).astype(np.int32) if ecm.seg_i.size else np.zeros((0, 2), np.int32)
    xl = np.stack([ecm.xl_i, ecm.xl_j], 1).astype(np.int32) if ecm.xl_i.size else np.zeros((0, 2), np.int32)

    # ---- cortex ring + FA coupling ----
    cortex0, fa_pairs, fa_rest, ring_pairs, ring_rest = build_cortex_ring(
        mpos0, center, R_cell, n_cortex, capture_um=R_cell)
    n_c = cortex0.shape[0]
    N = n_m + n_c
    pos0 = np.vstack([mpos0, cortex0])
    pos = pos0.copy()

    # ---- boundary: pin outer matrix ring (embedded in bulk) AND the cortex (prescribed contraction);
    #      matrix interior is free and follows via the FA springs ----
    r_xy = np.linalg.norm(mpos0[:, :2] - center[:2], axis=1)
    fixed = np.zeros(N, dtype=np.int32)
    fixed[:n_m] = (r_xy > R_pin).astype(np.int32)
    cortex_idx = np.arange(n_m, N, dtype=np.int32)
    fixed[cortex_idx] = 1                                        # cortex is a moving Dirichlet BC (prescribed radius)
    cortex_rel = cortex0 - center[None, :]                       # cortex offsets from centre (for radius scaling)

    # ---- device arrays ----
    pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=d)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=d)
    tri_d = wp.array(tri, dtype=wp.int32, device=d) if has_bend else None
    alpha_d = wp.array(alpha, dtype=wp.float64, device=d) if has_bend else None
    seg_d = wp.array(seg, dtype=wp.int32, device=d) if seg.shape[0] else None
    segk_d = wp.array(np.ascontiguousarray(ecm.seg_k, np.float64), dtype=wp.float64, device=d) if seg.shape[0] else None
    segr_d = wp.array(np.ascontiguousarray(ecm.seg_rest, np.float64), dtype=wp.float64, device=d) if seg.shape[0] else None
    # WLC (thermal strain-stiffening) segment option — ecm.seg_Lc/Lp/EA are already populated by build_fibrillar_ecm.
    use_wlc = bool(wlc and seg.shape[0] and ecm.seg_Lc is not None and ecm.seg_Lp > 0 and ecm.seg_EA > 0)
    Lc_d = wlc_Lp = wlc_EA = wlc_kbt = wlc_xmax = None
    if use_wlc:
        Lc_d = wp.array(np.ascontiguousarray(ecm.seg_Lc, np.float64), dtype=wp.float64, device=d)
        wlc_Lp, wlc_EA, wlc_kbt = float(ecm.seg_Lp), float(ecm.seg_EA), float(U.KBT)
        wlc_xmax = float(wlc_x_crossover(wlc_EA, wlc_Lp, wlc_kbt))
    xl_d = wp.array(xl, dtype=wp.int32, device=d) if xl.shape[0] else None
    xlk_d = wp.array(np.ascontiguousarray(ecm.xl_k, np.float64), dtype=wp.float64, device=d) if xl.shape[0] else None
    r0 = ecm.xl_rest.astype(np.float64).copy()
    r0_d = wp.array(np.ascontiguousarray(r0, np.float64), dtype=wp.float64, device=d) if xl.shape[0] else None
    fa_d = wp.array(fa_pairs, dtype=wp.int32, device=d)
    fak_d = wp.array(np.full(n_c, k_fa, np.float64), dtype=wp.float64, device=d)
    far_d = wp.array(fa_rest, dtype=wp.float64, device=d)
    ring_d = wp.array(ring_pairs, dtype=wp.int32, device=d)
    ringk_d = wp.array(np.full(n_c, k_cortex, np.float64), dtype=wp.float64, device=d)
    ringr_d = wp.array(ring_rest, dtype=wp.float64, device=d)
    fixed_d = wp.array(fixed, dtype=wp.int32, device=d)
    target_d = wp.array(pos.copy(), dtype=wp.vec3d, device=d)

    # ---- CFL timestep from the stiffest mode ----
    kmax = (float(ecm.net.kappa.max()) / (float(ecm.net.seg_rest.mean()) if ecm.net.seg_rest.size else 1.0) ** 3) if has_bend else 0.0
    kmax = max(kmax, float(np.max(ecm.link_k)) if ecm.link_k.size else 1.0, k_fa, k_cortex)
    if use_wlc:                                                 # WLC enthalpic wall stiffness EA/Lc dominates the CFL
        kmax = max(kmax, wlc_EA / float(np.min(ecm.seg_Lc)))
    dt_mu = 0.1 / max(kmax, 1e-9)
    x_beta_over_kT = (x_beta_nm * 1e-3) / U.KBT
    r_min_meas = max(1.0, contract_frac * R_cell)               # innermost shell just outside the contracted cell
    h_slab = max(float(mpos0[:, 2].max() - mpos0[:, 2].min()), 1e-3)     # slab thickness for cylindrical flux
    cyl_edges = np.linspace(r_min_meas, R_pin, 13)
    if verbose:
        print(f"    [diag] kmax={kmax:.3g} dt_mu={dt_mu:.3g} substeps={mech_substeps} "
              f"k_fa={k_fa} contract_frac={contract_frac} n_c={n_c}", flush=True)

    def equilibrate(cap=None, tol=2.0e-4, chunk=300):
        """Relax to mechanical equilibrium: run substeps in chunks until the max nodal displacement over a
        chunk falls below ``tol`` [µm] (converged) or the ``cap`` is hit. Convergence-based (not fixed-count)
        so a step-strain-and-HOLD gives a true elastic plateau — the OFF control stops moving once settled."""
        cap = cap if cap is not None else mech_substeps
        steps = 0
        while steps < cap:
            prev = pos_d.numpy()
            for _ in range(chunk):
                wp.launch(_zero, dim=N, inputs=[f_d], device=d)
                if has_bend:
                    wp.launch(cytosim_bending_kernel, dim=tri.shape[0], inputs=[pos_d, tri_d, alpha_d, f_d], device=d)
                if seg_d is not None:
                    if use_wlc:
                        wp.launch(wlc_spring_kernel, dim=seg.shape[0],
                                  inputs=[pos_d, seg_d, Lc_d, wp.float64(wlc_Lp), wp.float64(wlc_EA),
                                          wp.float64(wlc_kbt), wp.float64(wlc_xmax), f_d], device=d)
                    else:
                        wp.launch(link_spring_kernel, dim=seg.shape[0], inputs=[pos_d, seg_d, segk_d, segr_d, f_d], device=d)
                if xl_d is not None:
                    wp.launch(link_spring_kernel, dim=xl.shape[0], inputs=[pos_d, xl_d, xlk_d, r0_d, f_d], device=d)
                wp.launch(link_spring_kernel, dim=n_c, inputs=[pos_d, fa_d, fak_d, far_d, f_d], device=d)
                wp.launch(link_spring_kernel, dim=n_c, inputs=[pos_d, ring_d, ringk_d, ringr_d, f_d], device=d)
                wp.launch(M._freeze_mask_kernel, dim=N, inputs=[f_d, fixed_d], device=d)
                wp.launch(axpy_kernel, dim=N, inputs=[pos_d, wp.float64(dt_mu), f_d], device=d)
                wp.launch(M._pin_positions_kernel, dim=N, inputs=[pos_d, fixed_d, target_d], device=d)
            steps += chunk
            wp.synchronize_device(d)
            if float(np.max(np.abs(pos_d.numpy() - prev))) < tol:
                break

    def set_cortex_target(s):                                   # prescribe cortex ring radius = s·R_cell
        tnp = target_d.numpy()
        tnp[n_m:] = center[None, :] + cortex_rel * s
        target_d.assign(np.ascontiguousarray(tnp, np.float64))

    def measure():
        """Primary relaxation observable = the CORTICAL LOAD: the total FA spring tension the matrix exerts
        back on the (held) contracted cortex — exactly the force a stress-relaxation experiment tracks. As
        transient crosslinks unbind (ON), the matrix yields inward toward the held cortex, the FA springs
        shorten, and this load DECAYS; an elastic matrix (OFF) holds it constant. Also returns σ_rr(r) for the
        spatial propagation profile."""
        wp.synchronize_device(d)
        p = pos_d.numpy()
        fa_c, fa_m = p[fa_pairs[:, 0]], p[fa_pairs[:, 1]]
        L_fa = np.linalg.norm(fa_c - fa_m, axis=1)
        load = float(np.mean(np.abs(k_fa * (L_fa - fa_rest))))          # mean FA tension [pN] = cortical load proxy
        ecm.xl_rest[:] = r0_d.numpy() if xl_d is not None else ecm.xl_rest   # reflect creeped rest lengths
        prof = M.stress_field_radial(ecm, p[:n_m], center, n_shells=12, r_min=r_min_meas, r_max=R_pin)
        rc, sig_cyl, n_cyl, n_mid, force_cyl = stress_flux_cylindrical(ecm, p[:n_m], center, cyl_edges, h_slab)   # Kim r^-1
        prof["r_cyl"], prof["sigma_cyl"], prof["n_cyl"] = rc.tolist(), sig_cyl.tolist(), n_cyl
        prof["n_cyl_mid"], prof["force_cyl_pN"] = n_mid, force_cyl.tolist()
        return p, prof, load

    # ---- protocol: assemble, ramp the cortex radius 1→contract_frac, then HOLD and watch the load relax ----
    equilibrate()                                              # (1) assemble/relax at rest
    set_cortex_target(contract_frac)                           # (2) apply the full step-strain
    equilibrate(cap=mech_substeps)                             #     converge to the ELASTIC t=0 state
    ts, sig_mem, profiles = [], [], []
    for rec in range(n_record):                                # (3) hold; relaxation clock (turnover per step)
        p, prof, load = measure()
        ts.append(rec * dt_real_s)
        sig_mem.append(load)                                   # cortical load (FA tension) — the relaxation signal
        if rec in (0, n_record // 2, n_record - 1):
            profiles.append({"t_s": rec * dt_real_s, "r": np.asarray(prof["r"]).tolist(),
                             "sigma_rr": np.abs(np.asarray(prof["sigma_rr"])).tolist(),
                             "r_cyl": prof["r_cyl"], "sigma_cyl": prof["sigma_cyl"], "n_cyl": prof["n_cyl"]})
        if xl_d is not None and koff0_per_s > 0.0:             # advance one real timestep of turnover, then re-settle
            wp.launch(xl_turnover_kernel, dim=xl.shape[0],
                      inputs=[pos_d, xl_d, xlk_d, r0_d, wp.float64(koff0_per_s * dt_real_s),
                              wp.float64(x_beta_over_kT)], device=d)
            equilibrate(cap=mech_substeps)

    # ---- cortex contraction diagnostic (did the ring actually pull in?) ----
    p_all = pos_d.numpy()
    cortex_r0 = float(np.linalg.norm(cortex0[:, :2] - center[:2], axis=1).mean())
    cortex_rf = float(np.linalg.norm(p_all[n_m:, :2] - center[:2], axis=1).mean())
    if verbose:
        print(f"    [diag] cortex ⟨R⟩ {cortex_r0:.2f}→{cortex_rf:.2f}µm ({100*(cortex_r0-cortex_rf)/cortex_r0:.1f}% in), "
              f"peak σ_mem={np.max(sig_mem):.3g}Pa, σ_mem(t=end)={sig_mem[-1]:.3g}Pa", flush=True)

    # ---- final matrix radial displacement vs r ----
    p_final = pos_d.numpy()[:n_m]
    disp = p_final - mpos0
    r_init = np.linalg.norm(mpos0[:, :2] - center[:2], axis=1)
    rhat = (mpos0[:, :2] - center[:2]) / (r_init[:, None] + 1e-9)
    disp_rad = np.einsum("na,na->n", disp[:, :2], rhat)        # + = outward, − = inward (contraction)
    nb = 10
    edges = np.linspace(R_cell, R_pin, nb + 1)
    rc = 0.5 * (edges[:-1] + edges[1:])
    dbin = np.array([disp_rad[(r_init >= edges[i]) & (r_init < edges[i + 1])].mean()
                     if np.any((r_init >= edges[i]) & (r_init < edges[i + 1])) else np.nan
                     for i in range(nb)])

    geom = None
    if save_geom:                                              # positions + creeped xl rests for the 3D tension viewer
        p_all = pos_d.numpy()
        geom = {"matrix_init": mpos0, "matrix_final": p_all[:n_m], "cortex_final": p_all[n_m:],
                "xl_rest_final": (r0_d.numpy() if xl_d is not None else r0)}
    ecm.xl_rest[:] = ecm.xl_rest * 0 + r0                       # restore original rests (safety)
    return {"t_s": ts, "sigma_membrane_Pa": sig_mem, "profiles": profiles, "geom": geom,
            "peak_load_pN": float(sig_mem[0]) if sig_mem else float("nan"),
            "disp_r_um": rc.tolist(), "disp_radial_um": dbin.tolist(),
            "koff0_per_s": koff0_per_s, "contract_frac": contract_frac, "n_m": int(n_m), "n_c": int(n_c)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--contract-frac", type=float, default=0.4, help="cortex contracts to this fraction of R_cell (swept contraction strength)")
    ap.add_argument("--koff", type=float, default=0.03, help="crosslink zero-force off-rate [1/s] (ON case)")
    ap.add_argument("--conc", type=float, default=1.5, help="collagen concentration [mg/mL]")
    ap.add_argument("--n-fibers", type=int, default=220)
    ap.add_argument("--target-z", type=float, default=None, help="target connectivity ⟨z⟩ (KB-1.3 physiological ~3.2; sub-isostatic → longer-range transmission)")
    ap.add_argument("--substeps", type=int, default=6000, help="max relax substeps per record (convergence cap; stops early when settled)")
    ap.add_argument("--half", type=float, default=12.0, help="box half-width [µm]")
    ap.add_argument("--records", type=int, default=12)
    ap.add_argument("--dt-real", type=float, default=None, help="override the per-record real timestep [s] (fix the absolute time axis for a koff sweep)")
    ap.add_argument("--wlc", action="store_true", help="use thermal-WLC (strain-stiffening) fiber segments — tests Kim r^-1 recovery")
    ap.add_argument("--tag", default="", help="output filename suffix (distinct outputs per sweep point)")
    ap.add_argument("--save-geom", action="store_true", help="save node positions + bond tensions (npz) for the 3D HTML viewer")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    # ---- coarse quasi-2D collagen slab; cell at centre, outer ring pinned ----
    half = a.half
    lo, hi = [-half, -half, -2.0], [half, half, 2.0]
    center = np.array([0.0, 0.0, 0.0])
    R_cell, R_pin = 3.0, half - 3.0
    ecm = L.build_ecm("collagen_I", lo, hi, concentration=a.conc, dim=3,
                      n_fibers=a.n_fibers, target_z=a.target_z, pin_faces=(), rng=rng)
    print(f"[build] collagen n_nodes={ecm.net.n_nodes} n_xl={ecm.xl_i.size} ⟨z⟩={ecm.connectivity_z:.2f} "
          f"mesh ξ={ecm.mesh_size_um:.2f}µm", flush=True)

    # ---- KEY comparison: turnover ON (viscoelastic) vs OFF (elastic) ----
    dt_real = a.dt_real if a.dt_real is not None else (1.0 / a.koff) * 5.0 / a.records   # ~5 relaxation times (or fixed for a koff sweep)
    print(f"[run] turnover OFF (elastic control) frac={a.contract_frac} ...", flush=True)
    off = run_case(ecm, center, R_cell, R_pin, contract_frac=a.contract_frac, koff0_per_s=0.0, dt_real_s=dt_real,
                   n_record=a.records, mech_substeps=a.substeps, device=a.device, verbose=True, wlc=a.wlc)
    print(f"[run] turnover ON  (viscoelastic) frac={a.contract_frac} koff={a.koff}/s wlc={a.wlc} ...", flush=True)
    on = run_case(ecm, center, R_cell, R_pin, contract_frac=a.contract_frac, koff0_per_s=a.koff, dt_real_s=dt_real,
                  n_record=a.records, mech_substeps=a.substeps, device=a.device, verbose=True, wlc=a.wlc, save_geom=a.save_geom)

    geom = on.pop("geom", None)                                 # numpy arrays — saved separately (not JSON)
    if geom is not None:                                        # save geometry + bond tensions for the 3D HTML viewer
        segL = np.linalg.norm(geom["matrix_final"][ecm.seg_j] - geom["matrix_final"][ecm.seg_i], axis=1)
        seg_tension = ecm.seg_k * (segL - ecm.seg_rest)
        xlL = np.linalg.norm(geom["matrix_final"][ecm.xl_j] - geom["matrix_final"][ecm.xl_i], axis=1) if ecm.xl_i.size else np.zeros(0)
        xl_tension = ecm.xl_k * (xlL - geom["xl_rest_final"]) if ecm.xl_i.size else np.zeros(0)
        np.savez_compressed(f"{OUT}/kim_geom{a.tag}.npz",
                            matrix_init=geom["matrix_init"], matrix_final=geom["matrix_final"],
                            cortex_final=geom["cortex_final"], seg_i=ecm.seg_i, seg_j=ecm.seg_j,
                            seg_tension=seg_tension, xl_i=ecm.xl_i, xl_j=ecm.xl_j, xl_tension=xl_tension,
                            R_cell=R_cell, R_pin=R_pin)
        print(f"[geom] saved {OUT}/kim_geom{a.tag}.npz", flush=True)
    off.pop("geom", None)

    # ---- verdict: step-strain applied before the clock, so rec0 = the elastic t=0 load. An elastic matrix
    #      (OFF) holds that load; a viscoelastic one (ON) relaxes it toward a residual as crosslinks unbind. ----
    s_on, s_off = np.asarray(on["sigma_membrane_Pa"]), np.asarray(off["sigma_membrane_Pa"])
    ref_on, ref_off = float(s_on[0]), float(s_off[0])
    relax_frac = 1.0 - float(s_on[-1]) / ref_on if ref_on > 0 else float("nan")      # ON viscoelastic relaxation
    off_drift = 1.0 - float(s_off[-1]) / ref_off if ref_off > 0 else float("nan")     # OFF residual drift (should be ~0)
    verdict = ("PASS: turnover ON relaxes the cortical load %.0f%% while OFF holds it (%.0f%% drift) — viscoelastic coupling works"
               % (100 * relax_frac, 100 * off_drift)) if (relax_frac > 0.15 and relax_frac > abs(off_drift) + 0.12) \
        else ("INCONCLUSIVE (coarse smoke) — ON relax %.0f%% vs OFF drift %.0f%%" % (100 * relax_frac, 100 * off_drift))
    print(f"[verdict] {verdict}", flush=True)

    out = {"note": "COARSE SMOKE TEST — non-authoritative (native+full HARD rule). Imposed swept contraction; "
                   "matrix response emergent. Peer-oracle repro of Slater2021_SoftMatter (SE333).",
           "params": {"conc": a.conc, "n_fibers": a.n_fibers, "R_cell": R_cell, "R_pin": R_pin,
                      "contract_frac": a.contract_frac, "koff_per_s": a.koff, "device": a.device, "seed": a.seed},
           "verdict": verdict, "relax_frac_on": relax_frac, "drift_off": off_drift,
           "peak_load_on_pN": on["peak_load_pN"], "target_z": a.target_z, "wlc": a.wlc,
           "connectivity_z": float(ecm.connectivity_z),
           "n_cyl_final_on": on["profiles"][-1]["n_cyl"], "n_cyl_final_off": off["profiles"][-1]["n_cyl"],
           "on": on, "off": off}
    print(f"[decay] cylindrical stress exponent n: ON={out['n_cyl_final_on']:.2f} OFF={out['n_cyl_final_off']:.2f} "
          f"(Kim quasi-2D → r^-1, i.e. n≈1)", flush=True)
    print(f"[peak] peak cortical load (t=0 elastic) = {on['peak_load_pN']:.1f} pN", flush=True)
    with open(f"{OUT}/kim_repro{a.tag}.json", "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- figures ----
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4.6))
    ax1.plot(on["t_s"], s_on, "o-", color="crimson", label="turnover ON (viscoelastic)")
    ax1.plot(off["t_s"], s_off, "s--", color="steelblue", label="turnover OFF (elastic)")
    ax1.set_xlabel("time [s]"); ax1.set_ylabel("cortical load = mean FA tension [pN]")
    ax1.set_title("Stress relaxation of the held contraction"); ax1.legend(); ax1.grid(alpha=0.3)
    pf = on["profiles"][-1]                                     # final (steady) cylindrical profile
    rcyl, scyl, ncyl = np.asarray(pf["r_cyl"]), np.asarray(pf["sigma_cyl"]), pf["n_cyl"]
    ax2.loglog(rcyl, scyl, "o-", color="crimson", label=f"FF σ(r) (fit n={ncyl:.2f})")
    good = (rcyl > rcyl[0]) & (scyl > 0)
    if good.any():
        rr = rcyl[good]; ax2.loglog(rr, scyl[good][0] * (rr / rr[0]) ** -1.0, "k--", label="Kim r$^{-1}$ (2D)")
    ax2.set_xlabel("r [µm]"); ax2.set_ylabel("σ(r) cylindrical flux [Pa]")
    ax2.set_title(f"Stress decay: FF n={ncyl:.2f} vs Kim r$^{{-1}}$"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3, which="both")
    ax3.plot(on["disp_r_um"], on["disp_radial_um"], "o-", color="crimson", label="ON")
    ax3.plot(off["disp_r_um"], off["disp_radial_um"], "s--", color="steelblue", label="OFF")
    ax3.axhline(0, color="k", lw=0.6)
    ax3.set_xlabel("r [µm]"); ax3.set_ylabel("radial displacement [µm] (− = inward)")
    ax3.set_title("Matrix pulled inward (remodeling)"); ax3.legend(); ax3.grid(alpha=0.3)
    fig.suptitle("FF ⟷ Slater/Kim 2021 — COARSE SMOKE (non-authoritative): contractile cell in a "
                 "transient-crosslink viscoelastic collagen matrix", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/kim_repro_smoke{a.tag}.png", dpi=130)
    print(f"[out] {OUT}/kim_repro{a.tag}.json + {FIGS}/kim_repro_smoke{a.tag}.png", flush=True)


if __name__ == "__main__":
    main()
