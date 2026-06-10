"""CAPSTONE — test the PI's layer-2 spheroid spreading law A/A0 = a + b/R + c/R²
on the VALIDATED native-mesh + bilinear-tent DCM stack, via a SIZE SWEEP.

The law A/A0(R) is a SIZE-DEPENDENCE of the equilibrium spread footprint. Live
cell division on the native mesh = mesh-split vs HOOMD fixed-tag-space (a ~6-7
week port, out of scope). The size-dependence is instead probed directly by
BUILDING spheroids at several cell counts N (-> several effective radii R) and
measuring each one's EQUILIBRIUM A/A0 — which sidesteps mesh-split entirely.

PI LAW (target): A/A0 = a + b/R + c/R²,  a=-0.33, b=188.7 µm, c=-2655 µm²
  (r²=0.98, R=31-78 µm).  Signature: A/A0 DECREASES with R (smaller spheroids
  spread MORE per area — b/R surface/volume rim term, c/R² cohesion penalty).

Run FROM REPO ROOT:
    ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/dcm_native_capstone.py

Reuses (does NOT modify):
  * cell/dcm_native_shell.build_native_dcm_simulation  — VALIDATED stack
  * scripts/dcm_native_smoke._footprint_area / measurement convention
  * cell/dcm_spheroid_state — 3-zone necrosis depth classification + thresholds

Outputs:
  * outputs/h_dcm_native/figs/capstone_law_fit.png
  * outputs/h_dcm_native/capstone_sweep.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from ffn_sim.cell.dcm_native_shell import (
    ResolvedNativeDCM,
    build_native_dcm_simulation,
)
from ffn_sim.cell.dcm_spheroid_state import CellState

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_native")
FIGS = os.path.join(OUT, "figs")
os.makedirs(FIGS, exist_ok=True)

_UM = 1.0e6

# 3-zone necrosis depth thresholds (literature; dcm_spheroid_state defaults).
D_PROLIF_UM = 40.0     # d < this -> PROLIFERATING (rim)
D_NECROTIC_UM = 150.0  # d >= this -> NECROTIC (core)

# Size sweep cell counts (-> different R). Denser sampling so the continuum size
# trend emerges over the integer FCC-packing granularity at small N.
N_SWEEP = [8, 12, 16, 20, 28, 36, 48]


# ---------------------------------------------------------------------------
# Geometry / measurement helpers (reuse the smoke's conventions)
# ---------------------------------------------------------------------------
def _positions(sim):
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _footprint_area(pos_g, z0, R):
    """Spheroid footprint A = convex-hull area of substrate-contacting nodes (xy).

    Same definition as scripts/dcm_native_smoke._footprint_area: nodes whose
    height is below z0 + 0.5 R are "in contact"; the xy convex-hull area of those
    nodes is the basal footprint.
    """
    from scipy.spatial import ConvexHull
    contact = pos_g[pos_g[:, 2] < z0 + 0.5 * R][:, :2]
    if contact.shape[0] < 3:
        contact = pos_g[:, :2]
    try:
        return float(ConvexHull(contact).volume)  # 2D hull "volume" = area
    except Exception:  # noqa: BLE001
        return float(np.pi * (np.ptp(contact[:, 0]) / 2) * (np.ptp(contact[:, 1]) / 2))


def _cell_centroids(pos_g, ranges):
    return np.array([pos_g[a:b].mean(axis=0) for (a, b) in ranges])


def _effective_radius(cents):
    """Effective spheroid radius R from cell centroid positions (µm).

    R_rms = RMS distance of cell centroids from the cluster centroid;
    R_max = max distance. We report both; R_eff used for the fit is R_max + one
    cell radius (the actual outer surface), consistent with the spheroid-surface
    convention in dcm_spheroid_state (outermost CENTROID + a patch radius).
    """
    cc = cents.mean(axis=0)
    r = np.linalg.norm(cents - cc, axis=1)
    return float(r.max()), float(np.sqrt((r ** 2).mean())), cc, r


def _necrosis_zones(cents, R_cell):
    """3-zone depth classification (reuse dcm_spheroid_state convention).

    depth d = R_cluster - r_cell, with R_cluster = max centroid radius + R_cell
    (outermost cell at depth ~0 = proliferating rim). Returns (states, depth_um,
    R_cluster_um, fractions dict).
    """
    cc = cents.mean(axis=0)
    r = np.linalg.norm(cents - cc, axis=1)
    R_cluster = r.max() + R_cell
    depth_um = (R_cluster - r) * _UM
    states = np.full(cents.shape[0], int(CellState.QUIESCENT), dtype=np.int64)
    states[depth_um < D_PROLIF_UM] = int(CellState.PROLIFERATING)
    states[depth_um >= D_NECROTIC_UM] = int(CellState.NECROTIC)
    n = float(cents.shape[0])
    frac = {
        "proliferating": float((states == int(CellState.PROLIFERATING)).sum() / n),
        "quiescent": float((states == int(CellState.QUIESCENT)).sum() / n),
        "necrotic": float((states == int(CellState.NECROTIC)).sum() / n),
    }
    return states, depth_um, R_cluster * _UM, frac


# ---------------------------------------------------------------------------
# One spheroid: build -> finite gate -> equilibrate -> measure
# ---------------------------------------------------------------------------
def run_one(n_cells, *, seed=7, dt=3.0e-10, n_equil_blocks=18, block=1000,
            verbose=True):
    """Build an n_cells native-mesh+tent spheroid, gate, equilibrate, measure.

    Returns a dict with R_eff, A/A0, the per-cell centroids, necrosis zones, and
    a finite-gate flag. Drops dt to 1e-10 and retries once on non-finite.
    """
    for trial_dt in (dt, 1.0e-10):
        p = ResolvedNativeDCM(
            subdivisions=2,            # validated-stack resolution (162 nodes/cell):
                                       # subdiv 1's ~7 basal nodes make the convex-hull
                                       # footprint pure noise (non-monotone, N=14==N=40);
                                       # subdiv 2 is the smoke's resolution (A/A0=2.36).
            cluster="3d", spacing_factor=2.3,
            adh_strength=1.0e8, rep_strength=1.0e8,
            W_cs_Jm2=0.5e-3, dt=trial_dt, seed=seed,
        )
        h = build_native_dcm_simulation(p, n_cells, substrate=True, contact=True)
        sim, ranges, vol = h["sim"], h["ranges"], h["volume"]
        R = p.R_cell

        # --- build -> run(0) -> run(2000) finite gate (per size, BEFORE equilibrate) ---
        sim.run(0)
        # INITIAL footprint (t=0, packed ball on substrate) — the per-spheroid
        # reference A0 (the brief's "initial ... reference"). Using each spheroid's
        # OWN initial footprint makes A/A0 a size-COMPARABLE spreading magnitude
        # (how much the footprint GREW from the un-spread state), removing the
        # raw cell-count confound that a fixed single-cell A0 = pi*R^2 carries
        # (a bigger ball trivially has a bigger absolute footprint).
        pg0 = _positions(sim)
        A0_self = _footprint_area(pg0, p.z_substrate, R)
        A0_cell = np.pi * R ** 2  # single-cell reference (smoke convention; logged too)
        sim.run(2000)
        pg = _positions(sim)
        if not np.all(np.isfinite(pg)) or not np.all(np.isfinite(vol.volume)):
            if verbose:
                print(f"  N={n_cells} dt={trial_dt:.0e}: run(2000) NON-FINITE "
                      "-> retry at smaller dt" if trial_dt == dt else
                      f"  N={n_cells}: NON-FINITE even at dt=1e-10 -> FAIL")
            if trial_dt == dt:
                continue
            return {"n_cells": n_cells, "finite_gate": False, "dt": trial_dt}
        if verbose:
            print(f"  N={n_cells} dt={trial_dt:.0e}: finite gate PASS")

        # --- equilibrate; track footprint to confirm it plateaus ---
        areas = []          # A(t) / A0_self  (size-comparable spreading magnitude)
        areas_abs = []      # raw footprint [µm²]  (for transparency)
        finite = True
        for _ in range(n_equil_blocks):
            sim.run(block)
            pg = _positions(sim)
            if not np.all(np.isfinite(pg)):
                finite = False
                break
            A_abs = _footprint_area(pg, p.z_substrate, R)
            areas.append(A_abs / A0_self)
            areas_abs.append(A_abs * _UM ** 2)
        if not finite:
            if trial_dt == dt:
                continue
            return {"n_cells": n_cells, "finite_gate": False, "dt": trial_dt}

        # equilibrium A/A0 = mean of the last third (plateau); tail std = the
        # measurement uncertainty (footprint wobble over the equilibrated tail).
        tail = areas[max(1, 2 * len(areas) // 3):]
        AoverA0 = float(np.mean(tail))
        AoverA0_tailstd = float(np.std(tail))

        cents = _cell_centroids(pg, ranges)
        R_max, R_rms, cc, r = _effective_radius(cents)
        R_eff_um = (R_max + R) * _UM  # outer surface radius (centroid + one R)
        states, depth_um, R_cluster_um, frac = _necrosis_zones(cents, R)

        if verbose:
            print(f"    R_eff={R_eff_um:.1f} µm (R_max={R_max*_UM:.1f}, "
                  f"R_rms={R_rms*_UM:.1f})  A/A0={AoverA0:.2f}  "
                  f"necrotic_frac={frac['necrotic']:.2f}")

        return {
            "n_cells": n_cells, "finite_gate": True, "dt": trial_dt,
            "R_eff_um": R_eff_um, "R_max_um": float(R_max * _UM),
            "R_rms_um": float(R_rms * _UM), "R_cluster_um": R_cluster_um,
            "AoverA0": AoverA0, "AoverA0_tailstd": AoverA0_tailstd,
            "A_traj": [float(a) for a in areas],
            "A_abs_um2_traj": [float(a) for a in areas_abs],
            "A0_self_um2": float(A0_self * _UM ** 2),
            "A0_cell_um2": float(A0_cell * _UM ** 2),
            "AoverA0_cellref": float(np.mean(
                [a / A0_cell for a in
                 [_footprint_area(pg, p.z_substrate, R)]])),
            "centroids_um": (cents * _UM).tolist(),
            "states": states.tolist(), "depth_um": depth_um.tolist(),
            "necrotic_frac": frac["necrotic"], "quiescent_frac": frac["quiescent"],
            "prolif_frac": frac["proliferating"],
        }
    return {"n_cells": n_cells, "finite_gate": False, "dt": dt}


# ---------------------------------------------------------------------------
# Fit A/A0 = a + b/R + c/R²
# ---------------------------------------------------------------------------
def fit_law(R_um, AoverA0):
    """Least-squares fit A/A0 = a + b/R + c/R² on the (R, A/A0) points.

    Basis [1, 1/R, 1/R²]. Returns (a, b, c, r²). r² = 1 - SS_res/SS_tot.
    """
    R = np.asarray(R_um, dtype=float)
    y = np.asarray(AoverA0, dtype=float)
    X = np.column_stack([np.ones_like(R), 1.0 / R, 1.0 / R ** 2])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b, c = coef
    yhat = X @ coef
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(a), float(b), float(c), float(r2)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(means, seed_rows, fit, pi_law):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R = np.array([m["R_eff_um"] for m in means])
    A = np.array([m["AoverA0"] for m in means])
    Aerr = np.array([m["AoverA0_std"] for m in means])
    a, b, c, r2 = fit
    pa, pb, pc = pi_law

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    # (a) A/A0 vs R: per-realization (thin) + size mean±std + fitted law + PI law
    ax0 = ax[0]
    sR = np.array([r["R_eff_um"] for r in seed_rows])
    sA = np.array([r["AoverA0"] for r in seed_rows])
    ax0.scatter(sR, sA, s=22, color="C3", alpha=0.35, zorder=3,
                label="per-realization (seeds)")
    ax0.errorbar(R, A, yerr=Aerr, fmt="o", ms=9, color="C3", capsize=4,
                 zorder=5, label="size mean ± seed std")
    Rg = np.linspace(R.min() * 0.9, R.max() * 1.1, 200)
    ax0.plot(Rg, a + b / Rg + c / Rg ** 2, "-", color="C0", lw=2,
             label=f"fit: a={a:.2f}, b={b:.0f}, c={c:.0f}\n(r²={r2:.3f})")
    Rpi = np.linspace(31, 78, 200)
    ax0.plot(Rpi, pa + pb / Rpi + pc / Rpi ** 2, "--", color="k", lw=1.6,
             label=f"PI law (R=31-78µm)\na={pa}, b={pb}, c={pc}")
    for m in means:
        ax0.annotate(f"N={m['n_cells']}", (m["R_eff_um"], m["AoverA0"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax0.set_xlabel("effective radius R [µm]")
    ax0.set_ylabel("equilibrium A / A0  (vs own initial footprint)")
    ax0.set_title("(a) capstone law fit: A/A0 = a + b/R + c/R²")
    ax0.legend(fontsize=8, loc="best")

    # (b) 3-zone state of the largest spheroid (xy projection, colored by zone)
    ax1 = ax[1]
    big = max(means, key=lambda m: m["n_cells"])
    cents = np.array(big["rep_centroids_um"])
    states = np.array(big["rep_states"])
    cc = cents.mean(axis=0)
    colmap = {int(CellState.PROLIFERATING): ("#2ca02c", "PROLIFERATING (rim)"),
              int(CellState.QUIESCENT): ("#ff7f0e", "QUIESCENT (shell)"),
              int(CellState.NECROTIC): ("#7f0000", "NECROTIC (core)")}
    for s, (col, lab) in colmap.items():
        m = states == s
        if m.any():
            ax1.scatter(cents[m, 0] - cc[0], cents[m, 1] - cc[1], s=120,
                        color=col, edgecolor="k", lw=0.5, label=lab)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x − x̄ [µm]"); ax1.set_ylabel("y − ȳ [µm]")
    ax1.set_title(f"(b) 3-zone state, largest spheroid (N={big['n_cells']}, "
                  f"R={big['rep_R_eff_um']:.0f}µm)")
    ax1.legend(fontsize=8, loc="best")

    # (c) necrotic / quiescent fraction vs R
    ax2 = ax[2]
    necf = np.array([m["necrotic_frac"] for m in means])
    quif = np.array([m["quiescent_frac"] for m in means])
    order = np.argsort(R)
    ax2.plot(R[order], necf[order], "o-", color="#7f0000", label="necrotic frac")
    ax2.plot(R[order], quif[order], "s--", color="#ff7f0e", label="quiescent frac")
    ax2.axvline(D_NECROTIC_UM, ls=":", color="k", lw=1.0,
                label=f"necrotic-onset depth {D_NECROTIC_UM:.0f}µm\n(coarse-grain scale)")
    ax2.set_xlabel("effective radius R [µm]")
    ax2.set_ylabel("cell fraction")
    ax2.set_title("(c) necrotic / quiescent fraction vs R")
    ax2.set_ylim(-0.03, 1.03)
    ax2.legend(fontsize=8, loc="best")

    fig.tight_layout()
    path = os.path.join(FIGS, "capstone_law_fit.png")
    fig.savefig(path, dpi=120)
    print(f"  figure -> {os.path.relpath(path)}")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
# The native FCC ball packing (_cluster_centers) + strong-adhesion equilibrium is
# effectively DETERMINISTIC (BAOAB thermal noise is negligible vs the adhesion
# forces; seeds 7/23/101 gave identical A/A0 to 2 dp). So 1 seed/size suffices;
# the per-size spread shown in the figure is the equilibration-tail std, not a
# seed ensemble. N_SEEDS>1 still works (averages seeds) if ever needed.
N_SEEDS = 1
SEEDS = [7, 23, 101]


def aggregate_size(n_cells):
    """Run N_SEEDS realizations at this size; return (per-seed rows, size-mean dict).

    The native FCC packing + strong-adhesion equilibrium is deterministic, so 1
    seed suffices (N_SEEDS=1). The per-size A/A0 uncertainty reported is then the
    equilibration-TAIL std (footprint wobble over the plateau) — a real measure of
    the measurement noise — used for the figure error bars (viz-integrity rule).
    With N_SEEDS>1 it falls back to the across-seed std.
    """
    seed_rows = []
    for s in SEEDS[:N_SEEDS]:
        r = run_one(n_cells, seed=s, verbose=False)
        if r.get("finite_gate"):
            seed_rows.append(r)
            print(f"    seed={s}: R_eff={r['R_eff_um']:.1f}µm  A/A0={r['AoverA0']:.2f}"
                  f"±{r['AoverA0_tailstd']:.2f}  nec_f={r['necrotic_frac']:.2f}")
        else:
            print(f"    seed={s}: finite gate FAIL")
    if not seed_rows:
        return [], {"n_cells": n_cells, "finite_gate": False}
    A = np.array([r["AoverA0"] for r in seed_rows])
    R = np.array([r["R_eff_um"] for r in seed_rows])
    # uncertainty: across-seed std if >1 seed, else the equilibration-tail std.
    A_err = float(A.std()) if len(seed_rows) > 1 else float(seed_rows[0]["AoverA0_tailstd"])
    # representative realization = the seed closest to the mean R (for the 3-zone panel)
    rep = seed_rows[int(np.argmin(np.abs(R - R.mean())))]
    mean = {
        "n_cells": n_cells, "finite_gate": True, "n_seeds": len(seed_rows),
        "R_eff_um": float(R.mean()), "R_eff_std_um": float(R.std()),
        "AoverA0": float(A.mean()), "AoverA0_std": A_err,
        "AoverA0_seeds": A.tolist(), "R_eff_seeds": R.tolist(),
        "necrotic_frac": float(np.mean([r["necrotic_frac"] for r in seed_rows])),
        "quiescent_frac": float(np.mean([r["quiescent_frac"] for r in seed_rows])),
        "prolif_frac": float(np.mean([r["prolif_frac"] for r in seed_rows])),
        "rep_centroids_um": rep["centroids_um"], "rep_states": rep["states"],
        "rep_R_eff_um": rep["R_eff_um"],
    }
    return seed_rows, mean


def main():
    PI_LAW = (-0.33, 188.7, -2655.0)  # a, b (µm), c (µm²)

    print("=== CAPSTONE size sweep: A/A0 = a + b/R + c/R² on native-mesh+tent DCM ===")
    print(f"    {N_SEEDS} seed(s)/size, subdiv 2 (validated res), equilibrate 18k "
          "steps, A0 = each spheroid's OWN initial footprint\n")
    means = []
    all_seed_rows = []
    for n in N_SWEEP:
        print(f"-- N={n} spheroid ({N_SEEDS} seeds) --")
        seed_rows, mean = aggregate_size(n)
        all_seed_rows.extend(seed_rows)
        if mean.get("finite_gate"):
            means.append(mean)

    if len(means) < 3:
        print(f"FAIL: only {len(means)} sizes passed (need >=3 to fit).")
        return 1

    R = [m["R_eff_um"] for m in means]
    A = [m["AoverA0"] for m in means]
    fit = fit_law(R, A)
    a, b, c, r2 = fit

    # law signature: does A/A0 DECREASE with R?
    order = np.argsort(R)
    A_sorted = np.array(A)[order]
    decreasing = bool(A_sorted[-1] < A_sorted[0])
    diffs = np.diff(A_sorted)
    mono_frac = float((diffs < 0).mean())
    # correlation of the size-trend (sign is the qualitative law check)
    corr = float(np.corrcoef(np.array(R), np.array(A))[0, 1])

    print("\n=== SWEEP TABLE (size means ± seed std) ===")
    print(f"  {'N':>4}  {'R_eff[µm]':>12}  {'A/A0':>14}  {'nec_f':>6}  {'qui_f':>6}")
    for m in sorted(means, key=lambda x: x["R_eff_um"]):
        print(f"  {m['n_cells']:>4}  {m['R_eff_um']:>7.1f}±{m['R_eff_std_um']:<4.1f}  "
              f"{m['AoverA0']:>7.2f}±{m['AoverA0_std']:<5.2f}  "
              f"{m['necrotic_frac']:>6.2f}  {m['quiescent_frac']:>6.2f}")

    print("\n=== FIT vs PI ===")
    print(f"  fitted : a={a:+.3f}  b={b:+.1f} µm  c={c:+.1f} µm²  r²={r2:.3f}")
    print(f"  PI law : a={PI_LAW[0]:+.3f}  b={PI_LAW[1]:+.1f} µm  c={PI_LAW[2]:+.1f} µm²  r²=0.98")
    print(f"  R range (this sweep): {min(R):.1f}-{max(R):.1f} µm  "
          f"(PI: 31-78 µm; COMPRESSED — CPU coarse scale)")
    print(f"  law signature A/A0 DECREASES with R: {decreasing}  "
          f"(corr(R,A/A0)={corr:+.2f}, monotone-decr frac {mono_frac:.2f})")

    sweep = {
        "n_sweep": N_SWEEP, "n_seeds": N_SEEDS, "seeds": SEEDS[:N_SEEDS],
        "subdivisions": 2,
        "measurement": "equilibrium A/A0 = mean(last third of footprint traj) "
                       "normalised by the spheroid's OWN initial (t=0) basal "
                       "footprint A0_self (size-comparable spreading magnitude); "
                       "size value = mean over seeds; R_eff = (R_max centroid + "
                       "R_cell). A0_cell=pi*R_cell^2 and AoverA0_cellref also logged.",
        "size_means": means,
        "seed_rows": all_seed_rows,
        "fit": {"a": a, "b_um": b, "c_um2": c, "r2": r2},
        "pi_law": {"a": PI_LAW[0], "b_um": PI_LAW[1], "c_um2": PI_LAW[2], "r2": 0.98,
                   "R_range_um": [31, 78]},
        "R_range_um": [min(R), max(R)],
        "corr_R_AoverA0": corr,
        "law_signature_decreasing": decreasing,
        "monotone_decreasing_fraction": mono_frac,
    }
    jpath = os.path.join(OUT, "capstone_sweep.json")
    with open(jpath, "w") as f:
        json.dump(sweep, f, indent=2)
    print(f"\n  json -> {os.path.relpath(jpath)}")

    try:
        make_figure(means, all_seed_rows, fit, PI_LAW)
    except Exception as e:  # noqa: BLE001
        print(f"  (figure skipped: {e})")

    fit_ok = r2 > 0.7
    print("\n=== CAPSTONE RESULT ===")
    print(f"  sizes passed finite gate : {len(means)}/{len(N_SWEEP)}")
    print(f"  3-param form fits (r²>0.7): {'YES' if fit_ok else 'NO'} (r²={r2:.3f})")
    print(f"  law signature (A/A0 down with R): {'YES' if decreasing else 'NO'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
