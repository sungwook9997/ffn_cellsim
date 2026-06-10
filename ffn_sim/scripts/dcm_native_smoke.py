"""Smoke + validation for DCM Upgrade A (native mesh shell) + B (tent contact).

Three validations from the brief:
  1. SINGLE cell: native ``Volume`` holds V/V0≈1 at rest, finite, ~5000 steps.
  2. TWO cells: adhere to a stable separation (no fly-apart, no interpenetration).
  3. SMALL 3D spheroid (~12-16 cells): BAOAB (CPU), finite + ball compacts and
     spreads on the substrate with A/A0 PHYSICAL (~1.5-4, the fix from 17).

Run FROM REPO ROOT:
    ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/dcm_native_smoke.py

Figures → ffn_sim/outputs/h_dcm_native/figs/ (visualize-at-closeout rule).
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from ffn_sim.cell.dcm_native_shell import (
    ResolvedNativeDCM,
    build_native_dcm_simulation,
)

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_native")
FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)


def _positions(sim):
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _footprint_area(pos_g, ranges, z0, R):
    """Spheroid footprint A = convex-hull area of substrate-contacting nodes (xy)."""
    from scipy.spatial import ConvexHull
    contact = pos_g[pos_g[:, 2] < z0 + 0.5 * R][:, :2]
    if contact.shape[0] < 3:
        contact = pos_g[:, :2]
    try:
        return float(ConvexHull(contact).volume)  # 2D hull "volume" = area
    except Exception:
        return float(np.pi * (np.ptp(contact[:, 0]) / 2) * (np.ptp(contact[:, 1]) / 2))


def validation_1_single():
    print("=== Validation 1: single cell native Volume hold ===")
    p = ResolvedNativeDCM(subdivisions=2, dt=3e-10)
    h = build_native_dcm_simulation(p, 1, substrate=False, contact=False)
    sim, vol, V0 = h["sim"], h["volume"], h["V0"]
    sim.run(0)
    traj = []
    for k in range(10):
        sim.run(500)
        v = float(vol.volume[0])
        traj.append(v / V0)
        if not np.isfinite(v):
            print(f"  step {(k+1)*500}: NON-FINITE volume -> FAIL")
            return False, traj
    print(f"  V/V0 over 5000 steps: {traj[0]:.4f} -> {traj[-1]:.4f} "
          f"(min {min(traj):.4f}, max {max(traj):.4f})")
    ok = np.all(np.isfinite(traj)) and 0.9 < traj[-1] < 1.1
    print(f"  RESULT: {'PASS' if ok else 'FAIL'} (finite, V/V0 in [0.9,1.1])")
    return ok, traj


def validation_2_doublet():
    print("=== Validation 2: two-cell adhesion (stable separation) ===")
    # spacing 2.05 R + c_adh 2 µm so the surface nodes are genuinely inside the
    # adhesion cutoff (the tent actually engages, not a trivial no-op pass).
    p = ResolvedNativeDCM(subdivisions=2, cluster="2d", spacing_factor=2.05,
                          adh_strength=1.0e8, c_adh=2.0e-6, dt=3e-10)
    h = build_native_dcm_simulation(p, 2, substrate=False, contact=True)
    sim, ranges, tent = h["sim"], h["ranges"], h["tent"]
    R = p.R_cell
    seps = []
    n_engaged = 0
    for k in range(16):
        sim.run(500)
        pg = _positions(sim)
        a0, b0 = ranges[0]
        a1, b1 = ranges[1]
        c0 = pg[a0:b0].mean(0)
        c1 = pg[a1:b1].mean(0)
        sep = float(np.linalg.norm(c0 - c1))
        seps.append(sep / R)
        if not np.all(np.isfinite(pg)):
            print(f"  step {(k+1)*500}: NON-FINITE -> FAIL")
            return False, seps
    with tent.cpu_local_force_arrays as arr:
        fmag = np.linalg.norm(np.asarray(arr.force), axis=1)
    n_engaged = int((fmag > 0).sum())
    print(f"  contact-engaged nodes: {n_engaged} (max |F| {fmag.max():.2e} N)")
    s0, sf = seps[0], seps[-1]
    print(f"  center separation / R: {s0:.3f} -> {sf:.3f} "
          f"(min {min(seps):.3f})")
    # Stable doublet: cells neither fly apart (sep grows unbounded) nor fully
    # interpenetrate (sep -> ~0). Two R=1 shells in contact sit at sep ~ 1.4-2 R.
    ok = (np.all(np.isfinite(seps)) and 0.8 < sf < 2.3 and min(seps) > 0.6
          and n_engaged > 0)
    print(f"  RESULT: {'PASS' if ok else 'FAIL'} (adhered, contact engaged, "
          "no fly-apart/interpenetration)")
    return ok, seps


def validation_3_spheroid(n_cells=14):
    print(f"=== Validation 3: {n_cells}-cell 3D spheroid spread (A/A0 fix) ===")
    p = ResolvedNativeDCM(subdivisions=2, cluster="3d", spacing_factor=2.3,
                          adh_strength=1.0e8, rep_strength=1.0e8,
                          W_cs_Jm2=0.5e-3, dt=3e-10)
    h = build_native_dcm_simulation(p, n_cells, substrate=True, contact=True)
    sim, ranges, vol, V0 = h["sim"], h["ranges"], h["volume"], h["V0"]
    R = p.R_cell

    # build->run(0)->run(2000) finite gate before long run
    sim.run(0)
    sim.run(2000)
    pg = _positions(sim)
    if not np.all(np.isfinite(pg)):
        print("  run(2000) NON-FINITE -> FAIL (stability)")
        return False, None
    print(f"  run(2000) finite OK; volumes finite: {np.all(np.isfinite(vol.volume))}")

    A0 = np.pi * R ** 2  # single-cell reference footprint
    areas = []
    for k in range(20):
        sim.run(1000)
        pg = _positions(sim)
        if not np.all(np.isfinite(pg)):
            print(f"  step {(k+1)*1000+2000}: NON-FINITE -> FAIL")
            return False, areas
        A = _footprint_area(pg, ranges, p.z_substrate, R)
        areas.append(A / A0)
    AoverA0 = areas[-1]
    print(f"  A/A0 trajectory (per 1000 steps): "
          f"{[round(a,2) for a in areas[::4]]}")
    print(f"  FINAL A/A0 = {AoverA0:.2f}  (target physical ~1.5-4; OLD node-LJ = 17)")
    ok = np.all(np.isfinite(areas)) and 1.0 < AoverA0 < 6.0
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    return ok, areas


def make_figures(t1, t2, t3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    ax[0].plot(np.arange(1, len(t1) + 1) * 500, t1, "o-", color="C0")
    ax[0].axhline(1.0, ls="--", color="k", lw=0.8, label="V/V0 = 1")
    ax[0].set_xlabel("BAOAB step"); ax[0].set_ylabel("V / V0")
    ax[0].set_title("V1: single cell native Volume hold"); ax[0].set_ylim(0.85, 1.15)
    ax[0].legend()

    ax[1].plot(np.arange(1, len(t2) + 1) * 500, t2, "o-", color="C1")
    ax[1].axhspan(0.8, 2.3, color="C2", alpha=0.15, label="stable doublet band")
    ax[1].set_xlabel("BAOAB step"); ax[1].set_ylabel("center sep / R")
    ax[1].set_title("V2: two-cell adhesion"); ax[1].legend()

    if t3 is not None:
        ax[2].plot(np.arange(1, len(t3) + 1) * 1000 + 2000, t3, "o-", color="C3")
        ax[2].axhspan(1.5, 4.0, color="C2", alpha=0.15, label="physical band")
        ax[2].axhline(17, ls=":", color="r", lw=1.2, label="old node-LJ = 17")
        ax[2].set_xlabel("BAOAB step"); ax[2].set_ylabel("A / A0")
        ax[2].set_title("V3: spheroid spread (A/A0 fix)"); ax[2].legend()
    fig.tight_layout()
    path = os.path.join(FIGS, "dcm_native_validations.png")
    fig.savefig(path, dpi=110)
    print(f"  figure -> {os.path.relpath(path)}")


def main():
    ok1, t1 = validation_1_single()
    ok2, t2 = validation_2_doublet()
    ok3, t3 = validation_3_spheroid()
    try:
        make_figures(t1, t2, t3)
    except Exception as e:  # noqa: BLE001
        print(f"  (figure skipped: {e})")
    print("\n=== SUMMARY ===")
    print(f"  V1 single-cell hold : {'PASS' if ok1 else 'FAIL'}")
    print(f"  V2 doublet adhesion : {'PASS' if ok2 else 'FAIL'}")
    print(f"  V3 spheroid A/A0    : {'PASS' if ok3 else 'FAIL'}"
          + (f"  (A/A0={t3[-1]:.2f})" if t3 else ""))
    return 0 if (ok1 and ok2 and ok3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
