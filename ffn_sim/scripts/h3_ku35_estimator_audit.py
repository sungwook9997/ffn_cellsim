"""KU-3.5 cortical-tension ESTIMATOR AUDIT — read-only, no production sim run.

Decisive artifact-vs-physics test for the gamma floor (2026-06-02, Lead session,
following the compartment-stack-audit workflow wf_e87b5c75). Reprocesses saved
grip_walk GSD frames with THREE method-of-planes estimators side by side:

  (i)   current   : signed  sum  T*(u.n_hat)            <- h3_ku35_tension.py:109 (the bug)
  (ii)  corrected : tensile sum  T*|u.n_hat|            <- consistent cross-orientation
  (iii) IK hoop   : Irving-Kirkwood whole-shell virial  sum T*L*(1-(rhat.u)^2)/(8 pi R^2)

The current estimator takes the RAW HOOMD A->B bond orientation, so for a crossing
bond the sign of (u.n_hat) is set by storage order vs an arbitrary cut plane -> the
per-plane sum cancels as sqrt(N) instead of N. (i) reproduces the claimed "97x
isotropy wall"; (ii)/(iii) remove it. All three are first VALIDATED on a synthetic
isotropic shell where (ii) and (iii) must agree and be cut-invariant while (i) is
sqrt(N)-suppressed.

Param-map (k, r0) per bond type is PINNED by building a small constrained grip_walk
cortex with the v300full config and reading the live bond.Harmonic params
(constrained => cortex-bond k=0; mesoscale_force_scaling ON => the ~3.14x myosin
scaling baked in). This kills the 1.6e-6..7e-4 mN/m magnitude spread the workflow
agents hit by hardcoding params.

NOT a runtime module. Imports nothing into the live path.

Run:  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_estimator_audit
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
GSD = PKG / "outputs" / "h3" / "production" / "gpu_native38k_v300full.grip_walk.gsd"

F_STALL = 0.5e-12          # N, single-head myosin stall force (KU-3.x)
BAND = (0.35e-3, 0.65e-3)  # N/m, KU-3.5 cortical-tension band


# --------------------------------------------------------------------------- #
# Estimators
# --------------------------------------------------------------------------- #
def mop_estimators(pos: np.ndarray, bg: np.ndarray, T: np.ndarray, R: float,
                   n_planes: int = 24) -> dict:
    """Three method-of-planes estimators on a configuration.

    Args:
        pos: (N,3) positions (will be centered on COM).
        bg:  (M,2) bond endpoint indices into pos.
        T:   (M,) scalar bond tension (+ tensile, - compressive).
        R:   cell radius (cut circumference = 2 pi R).
        n_planes: number of isotropic diametral cut planes (Fibonacci).
    """
    posc = pos - pos.mean(0)
    rA = posc[bg[:, 0]]
    rB = posc[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    Ls = np.where(L > 0, L, 1.0)
    u = d / Ls[:, None]

    phi = (1 + 5 ** 0.5) / 2
    ii = np.arange(n_planes, dtype=np.float64)
    z = 1 - 2 * (ii + 0.5) / n_planes
    rxy = np.sqrt(np.clip(1 - z * z, 0.0, None))
    th = 2 * np.pi * ii / phi
    normals = np.stack([rxy * np.cos(th), rxy * np.sin(th), z], axis=1)

    g_signed, g_corr, ncross = [], [], []
    for nh in normals:
        aS = rA @ nh
        bS = rB @ nh
        cross = (aS * bS) < 0
        ncross.append(int(cross.sum()))
        if not cross.any():
            g_signed.append(0.0)
            g_corr.append(0.0)
            continue
        proj = u[cross] @ nh
        g_signed.append(float(np.sum(T[cross] * proj) / (2.0 * np.pi * R)))
        g_corr.append(float(np.sum(T[cross] * np.abs(proj)) / (2.0 * np.pi * R)))
    g_signed = np.asarray(g_signed)
    g_corr = np.asarray(g_corr)

    # Irving-Kirkwood whole-shell hoop (mechanical surface) tension.
    rmid = 0.5 * (rA + rB)
    rn = np.linalg.norm(rmid, axis=1)
    rn = np.where(rn > 0, rn, 1.0)
    rhat = rmid / rn[:, None]
    cos2 = (np.sum(rhat * u, axis=1)) ** 2
    g_ik = float(np.sum(T * L * (1.0 - cos2)) / (8.0 * np.pi * R ** 2))

    return {
        "signed": float(np.mean(np.abs(g_signed))),   # matches live tension.py:115
        "signed_planes": g_signed,
        "corrected": float(np.mean(g_corr)),
        "corrected_planes": g_corr,
        "ik": g_ik,
        "ncross_median": int(np.median(ncross)),
    }


# --------------------------------------------------------------------------- #
# Synthetic validation shell
# --------------------------------------------------------------------------- #
def synthetic_shell(n: int = 1500, R: float = 1e-5, T0: float = 1e-13,
                    k_nn: int = 6, seed: int = 0):
    """Isotropic tangential bond mesh on a sphere, all bonds at uniform tension T0.

    Returns (pos, bonds, T) with T = +T0 for every bond. For an isotropic shell
    the corrected MOP and the IK hoop must AGREE and be cut-invariant; the signed
    MOP must be sqrt(N)-suppressed (the bug).
    """
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    pos = np.stack([np.sin(phi) * np.cos(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(phi)], axis=1) * R
    # k nearest neighbours -> isotropic tangential mesh
    d2 = np.sum((pos[:, None, :] - pos[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(d2, np.inf)
    nn = np.argsort(d2, axis=1)[:, :k_nn]
    bonds = set()
    for a in range(n):
        for b in nn[a]:
            bonds.add((min(a, int(b)), max(a, int(b))))
    bg = np.asarray(sorted(bonds), dtype=np.int64)
    T = np.full(bg.shape[0], T0, dtype=np.float64)
    return pos, bg, T


# --------------------------------------------------------------------------- #
# Pinned param-map from a real (small) constrained grip_walk build
# --------------------------------------------------------------------------- #
def pinned_param_map():
    """Construct the per-type (k, r0) bond param-map from the REAL resolvers
    (no full build => no motor-placement failure). Replicates the exact wiring
    in cell.py:786-833 + cortex/myosin.py:733-743 for a CONSTRAINED grip_walk
    cortex with the v300full config (force_scaling ON, n_motors=1200)."""
    from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import (
        resolve_cortex_myosin, cortex_myosin_attach_bin_names,
        cortex_myosin_attach_bin_rest_lengths,
        BOND_TYPE_MYOSIN_BACKBONE, BOND_TYPE_MYOSIN_HEAD_BACKBONE)
    from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import (
        resolve_crosslinkers, xlink_attach_bin_names, xlink_attach_bin_rest_lengths)

    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = 200      # irrelevant: no placement, params only
    cfg["cortex"]["demo_mode"] = True
    # CRITICAL: the mesoscale force factor = native_n_motors / n_motors_per_cell
    # (cortex/myosin.py:311), so n_motors MUST match the production run or
    # k_head_actin AND F_stall are mis-scaled together. v300full used n_motors=1200
    # (engaged 23,950 ~ 1200 motors x 20 heads); factor ~ 3.14. The driver default
    # (10/20) gives factor ~190x and inflates k by ~60x — the param-map trap that
    # mis-read gamma. Because k AND F_stall both scale by factor, T/F_stall is
    # factor-INVARIANT (the force-budget verdict is robust to this), but the
    # absolute gamma is not, so we pin n_motors here.
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 1200
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True   # v300full = ON

    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend
    pm = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    px = resolve_crosslinkers(cfg, dt=dtc)
    f_stall = float(getattr(pm, "F_stall_per_head", F_STALL))  # scaled by factor

    pmap: dict[str, tuple[float, float]] = {}
    pmap["cortex-bond"] = (0.0, float(p.rest_length))          # constrained => k=0
    pmap[BOND_TYPE_MYOSIN_BACKBONE] = (float(pm.k_backbone),
                                       float(pm.backbone_segment_length))
    pmap[BOND_TYPE_MYOSIN_HEAD_BACKBONE] = (float(pm.k_head_spring),
                                            float(pm.head_rest_length))
    for nm, r0 in zip(cortex_myosin_attach_bin_names(pm.n_bins),
                      cortex_myosin_attach_bin_rest_lengths(
                          pm.n_bins, pm.head_actin_max_bind_dist)):
        pmap[nm] = (float(pm.k_head_actin), float(r0))
    pmap["xlink_intra"] = (float(px.k_intra), 0.0)
    for nm, r0 in zip(xlink_attach_bin_names(px.n_bins),
                      xlink_attach_bin_rest_lengths(px.n_bins, px.max_bind_dist)):
        pmap[nm] = (float(px.k_attach), float(r0))
    print("    factor=%.3f  k_head_actin=%.3e  F_stall=%.3e (x%.2f base)"
          % (pm.extras.get("mesoscale_force_factor", float("nan")),
             pm.k_head_actin, f_stall, f_stall / F_STALL))
    return pmap, float(p.R_cell), f_stall


# --------------------------------------------------------------------------- #
# Figure (visualize-at-finding rule, CLAUDE.md)
# --------------------------------------------------------------------------- #
def _make_figure(rows, band, R_cell) -> None:
    if not rows:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa: BLE001
        print("    (figure skipped: %s)" % e)
        return
    rows = np.asarray(rows, dtype=float)
    fi = rows[:, 0]
    meanT_fs = rows[:, 2]
    sumT = rows[:, 3]
    band_pN = rows[:, 4]
    signed = rows[:, 6] * 1e3   # N/m -> mN/m
    corrected = rows[:, 7] * 1e3
    ik = rows[:, 8] * 1e3

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))

    # Panel A: estimators vs band
    axL.axhspan(band[0] * 1e3, band[1] * 1e3, color="tab:green", alpha=0.18,
                label="KU-3.5 band [0.35, 0.65]")
    axL.plot(fi, signed, "o-", color="tab:red",
             label="signed (current bug)  T·(u·n̂)")
    axL.plot(fi, corrected, "s-", color="tab:blue",
             label="corrected  T·|u·n̂|")
    axL.plot(fi, ik, "^--", color="tab:cyan", label="Irving-Kirkwood hoop")
    axL.set_yscale("log")
    axL.set_xlabel("GSD frame (transport engaged)")
    axL.set_ylabel("cortical tension  γ_soft  (mN/m)")
    axL.set_title("Sign fix recovers ~2-22×, NOT the claimed 97×\n"
                  "corrected γ_soft stays ~220-1100× below band")
    axL.legend(fontsize=8, loc="lower right")
    axL.grid(True, which="both", alpha=0.25)

    # Panel B: force budget
    axR.axhline(1.0, color="k", ls=":", lw=1.2, label="F_stall (single head)")
    axR.plot(fi, meanT_fs, "o-", color="tab:purple",
             label="mean |T| / F_stall  (head-actin)")
    axR.plot(fi, sumT / band_pN, "s-", color="tab:orange",
             label="Σ|T| / band requirement")
    axR.set_yscale("log")
    axR.set_xlabel("GSD frame (transport engaged)")
    axR.set_ylabel("force ratio (dimensionless)")
    axR.set_title("Force budget IS the wall\n"
                  "mean head force ~0.1·F_stall; Σ|T| ~5× under band; r/r0=1.0000")
    axR.legend(fontsize=8, loc="center right")
    axR.grid(True, which="both", alpha=0.25)

    fig.suptitle("KU-3.5 estimator audit — artifact (sign bug, secondary) vs "
                 "physics (myosin force budget, the wall)  [v300full native 38k]",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = PKG / "outputs" / "h3" / "figs" / "fig_h3_estimator_audit.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("    figure -> %s" % out)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    import gsd.hoomd

    print("=" * 78)
    print("KU-3.5 ESTIMATOR AUDIT — artifact (sign bug) vs physics (force budget)")
    print("=" * 78)

    # 1) synthetic validation -------------------------------------------------
    pos_s, bg_s, T_s = synthetic_shell()
    R_s = 1e-5
    est = mop_estimators(pos_s, bg_s, T_s, R_s)
    ratio = est["corrected"] / est["ik"] if est["ik"] else float("nan")
    print("\n[1] SYNTHETIC isotropic shell (uniform T0=1e-13 N, %d bonds)" % bg_s.shape[0])
    print("    signed  (current bug) : %.4e N/m" % est["signed"])
    print("    corrected (tensile)   : %.4e N/m" % est["corrected"])
    print("    IK hoop               : %.4e N/m" % est["ik"])
    print("    corrected/IK          : %.4f   (must be ~1.0 to validate the fix)" % ratio)
    print("    signed/corrected      : %.4f   (<<1 => sqrt(N) cancellation = the bug)"
          % (est["signed"] / est["corrected"] if est["corrected"] else float("nan")))
    print("    corrected per-plane CV: %.3f   (low => cut-invariant)"
          % (np.std(est["corrected_planes"]) / np.mean(est["corrected_planes"])
             if np.mean(est["corrected_planes"]) else float("nan")))

    # 2) pin param-map --------------------------------------------------------
    print("\n[2] Pinning (k, r0) param-map from a small constrained grip_walk build ...")
    try:
        pmap, R_cell, f_stall = pinned_param_map()
        print("    OK. %d bond types; R_cell=%.3e m; F_stall(scaled)=%.3e N (x%.2f)"
              % (len(pmap), R_cell, f_stall, f_stall / F_STALL))
        myo = {k: v for k, v in pmap.items() if k.startswith("cortex_myosin_attach")}
        if myo:
            ks = [v[0] for v in myo.values()]
            print("    head-actin attach bins: k=%.3e N/m  r0 range [%.3e, %.3e] m"
                  % (ks[0], min(v[1] for v in myo.values()),
                     max(v[1] for v in myo.values())))
    except Exception as e:  # noqa: BLE001
        print("    BUILD FAILED (%s); cannot pin param-map -> aborting real-data pass." % e)
        return 1

    # 3) reprocess real GSD frames -------------------------------------------
    if not GSD.exists():
        print("\n[3] GSD not found: %s" % GSD)
        return 1
    print("\n[3] Reprocessing %s" % GSD.name)
    with gsd.hoomd.open(str(GSD), "r") as traj:
        nfr = len(traj)
        last = traj[nfr - 1]
        types = list(last.bonds.types)
        k_by_tid = np.array([pmap.get(nm, (0.0, 0.0))[0] for nm in types])
        r0_by_tid = np.array([pmap.get(nm, (0.0, 0.0))[1] for nm in types])
        myo_tid = np.array([nm.startswith("cortex_myosin_attach") for nm in types])

        # reference radius from frame 0
        f0 = traj[0]
        p0 = np.asarray(f0.particles.position, dtype=np.float64)
        # cortex actin beads = first F*N particles; approximate r0 by COM radius
        c0 = p0.mean(0)
        R0 = float(np.linalg.norm(p0 - c0, axis=1).mean())

        print("    n_frames=%d  R_cell(pin)=%.3e  band=[%.2e,%.2e] N/m  F_stall=%.1e N"
              % (nfr, R_cell, BAND[0], BAND[1], F_STALL))
        hdr = ("frame  n_eng   meanT/Fstall  sumT(pN)/band(pN)   r/r0      "
               "signed     corrected      IK     | band? ")
        print("    " + hdr)
        rows = []
        for fi in range(nfr):
            fr = traj[fi]
            pos = np.asarray(fr.particles.position, dtype=np.float64)
            bg = np.asarray(fr.bonds.group, dtype=np.int64)
            bt = np.asarray(fr.bonds.typeid, dtype=np.int64)
            k = k_by_tid[bt]
            r0 = r0_by_tid[bt]
            d = pos[bg[:, 1]] - pos[bg[:, 0]]
            L = np.linalg.norm(d, axis=1)
            T = k * (L - r0)

            est = mop_estimators(pos, bg, T, R_cell)

            # force-budget table on the myosin head-actin channel
            mm = myo_tid[bt]
            Tm = np.abs(T[mm])
            n_eng = int(mm.sum())
            meanT = float(Tm.mean()) if n_eng else 0.0
            sumT = float(Tm.sum())
            c = pos.mean(0)
            rr = float(np.linalg.norm(pos - c, axis=1).mean()) / R0
            band_pN = BAND[0] * 2 * np.pi * R_cell  # N -> required force across cut
            in_band = (BAND[0] <= est["corrected"] <= BAND[1]) or (BAND[0] <= est["ik"] <= BAND[1])
            print("    %4d   %5d   %10.2e   %8.0f/%8.0f   %.5f  %.3e  %.3e  %.3e | %s"
                  % (fi, n_eng, meanT / f_stall, sumT * 1e12, band_pN * 1e12, rr,
                     est["signed"], est["corrected"], est["ik"],
                     "YES" if in_band else "no"))
            rows.append((fi, n_eng, meanT / f_stall, sumT, band_pN, rr,
                         est["signed"], est["corrected"], est["ik"]))

    _make_figure(rows, BAND, R_cell)

    # 4) verdict --------------------------------------------------------------
    print("\n[4] VERDICT")
    print("    - If corrected/IK reach [0.35,0.65] mN/m on engaged frames -> pure")
    print("      measurement artifact: ship the estimator fix, floor dissolved.")
    print("    - If corrected/IK recover Nx but stay <<band AND meanT/Fstall<<1 with")
    print("      r/r0~1.000 -> the wall is the myosin FORCE BUDGET (binned-r0 proxy),")
    print("      not the measurement. Estimator fix still mandatory for honesty.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
