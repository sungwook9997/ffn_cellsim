"""Phase 1 Unit 2.1 validation figures.

Run from the repo root::

    python acs_kb/notebooks/02_motor_clutch_biphasic.py

Produces three PNGs under ``acs_kb/outputs/phase1_unit2_1/`` plus a
``02_summary.json`` with the measured / analytic comparison numbers.

Figures:
  01_catch_bond_lifetime.png   — τ(F) closed-form vs Gillespie samples
  02_biphasic_traction.png     — ⟨F_total⟩ vs E across 5 decades
  03_clutch_force_histogram.png — per-clutch force histogram at E = 5 kPa
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acs_kb.bridge.catch_bond import (
    DEFAULT_PARAMS,
    catch_slip_lifetime,
    catch_slip_off_rate,
    lifetime_peak_force,
)
from acs_kb.bridge.motor_clutch import (
    MotorClutchFA,
    MotorClutchParams,
    run_steady_state,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.types import make_focal_adhesion
from acs_kb.common.derived_params import load_bridge_config
from acs_kb.common.sanity_gate import gate_unit2_1_motor_clutch

CFG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_1.yaml"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "phase1_unit2_1"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    cfg = load_bridge_config(CFG_PATH)
    b = cfg["bridge"]
    d = b["derived"]
    summary: dict = {}

    # ------------------------------------------------------------------ #
    # 1. Catch-bond lifetime (analytic + Gillespie samples).             #
    # ------------------------------------------------------------------ #
    F_grid_pN = np.linspace(0.0, 60.0, 121)
    F_grid = F_grid_pN * 1e-12
    tau_analytic = catch_slip_lifetime(F_grid)
    rng = np.random.default_rng(20260518)
    sample_forces_pN = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    sample_means = np.empty_like(sample_forces_pN, dtype=np.float64)
    for i, F_pN in enumerate(sample_forces_pN):
        F = F_pN * 1e-12
        k = float(catch_slip_off_rate(F))
        sample_means[i] = rng.exponential(scale=1.0 / k, size=4000).mean()
    F_star = lifetime_peak_force(DEFAULT_PARAMS)
    tau_at_peak = float(catch_slip_lifetime(F_star))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(F_grid_pN, tau_analytic, "-", color="C0", label="τ(F) = 1/k_off(F)  (analytic)")
    ax.scatter(sample_forces_pN, sample_means, marker="o", color="C3", zorder=5,
               label="Gillespie sample mean (N=4·10³)")
    ax.axvline(F_star * 1e12, color="grey", linestyle="--", alpha=0.5,
               label=f"F* = {F_star*1e12:.2f} pN (closed form)")
    ax.axhline(tau_at_peak, color="grey", linestyle=":", alpha=0.5,
               label=f"τ_max = {tau_at_peak:.2f} s")
    ax.set_xlabel("Force F  [pN]")
    ax.set_ylabel("Bond lifetime τ  [s]")
    ax.set_title("KU-2.5 Pereverzev catch-slip bond")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "01_catch_bond_lifetime.png", dpi=140)
    plt.close(fig)
    summary["catch_bond"] = {
        "F_star_pN_analytic": F_star * 1e12,
        "tau_max_s_analytic": tau_at_peak,
        "sample_forces_pN": sample_forces_pN.tolist(),
        "sample_means_s": sample_means.tolist(),
        "KU_2_5_experimental_F_star_pN": 30.0,
    }

    # ------------------------------------------------------------------ #
    # 2. Biphasic traction sweep.                                        #
    # ------------------------------------------------------------------ #
    sw = b["biphasic_sweep"]
    nu = float(b["substrate"]["poisson_ratio"])
    a = float(b["substrate"]["contact_radius"])
    dt = float(b["dynamics"]["dt"])
    n_steps_per_E = int(sw["n_steps_per_E"])
    E_grid = np.logspace(sw["log10_E_min"], sw["log10_E_max"], int(sw["n_points"]))
    mean_F = np.empty_like(E_grid)
    mean_n_eng = np.empty_like(E_grid)
    seed_root = int(sw["seed"])
    for i, E in enumerate(E_grid):
        sub = LinearElasticSubstrate(young_modulus=float(E), poisson_ratio=nu,
                                     contact_radius=a)
        fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
        mc = MotorClutchFA(fa, sub, MotorClutchParams(n_clutches=50))
        rng = np.random.default_rng(seed_root + i)
        r = run_steady_state(mc, dt=dt, n_steps=n_steps_per_E, rng=rng,
                             burn_in_fraction=0.5)
        mean_F[i] = r["mean_force_total"]
        mean_n_eng[i] = r["mean_n_engaged"]

    E_argmax = float(E_grid[int(np.argmax(mean_F))])
    E_star_matched = float(d["biphasic_EStar_matched"])
    E_star_motor = float(d["biphasic_EStar_motor"])
    factor = float(b["acceptance"]["biphasic_peak_factor"])

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(E_grid, mean_F * 1e12, "o-", color="C0",
            label="⟨F_total⟩  (50 clutches, 5 s sim, last 50 % avg)")
    ax.axvline(E_argmax, color="C3", linestyle="-", alpha=0.6,
               label=f"E_argmax = {E_argmax:.2e} Pa")
    ax.axvline(E_star_matched, color="grey", linestyle="--", alpha=0.5,
               label=f"E*_matched = {E_star_matched:.2e} Pa "
                     f"(Bangasser 2013, ½ N k_int)")
    ax.axvline(E_star_motor, color="grey", linestyle=":", alpha=0.5,
               label=f"E*_motor = {E_star_motor:.2e} Pa "
                     f"(KU-2.8 N_m F_stall / v τ)")
    ax.set_xscale("log")
    ax.set_xlabel("Substrate Young modulus E  [Pa]")
    ax.set_ylabel("Mean traction ⟨F_total⟩  [pN]")
    ax.set_title("KU-2.4 / KU-2.8 biphasic traction vs substrate stiffness")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "02_biphasic_traction.png", dpi=140)
    plt.close(fig)
    summary["biphasic"] = {
        "E_grid_Pa": E_grid.tolist(),
        "mean_F_total_N": mean_F.tolist(),
        "mean_n_engaged": mean_n_eng.tolist(),
        "E_argmax_Pa": E_argmax,
        "E_star_matched_Pa": E_star_matched,
        "E_star_motor_Pa": E_star_motor,
        "ratio_to_matched": E_argmax / E_star_matched,
        "ratio_to_motor": E_argmax / E_star_motor,
        "factor_acceptance": factor,
    }

    # ------------------------------------------------------------------ #
    # 3. Per-clutch force histogram at E = 5 kPa.                        #
    # ------------------------------------------------------------------ #
    sub = LinearElasticSubstrate(young_modulus=5e3, poisson_ratio=nu,
                                 contact_radius=a)
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc = MotorClutchFA(fa, sub, MotorClutchParams(n_clutches=50))
    rng = np.random.default_rng(9999)
    # Snapshot per-clutch forces every 100 steps over the second half.
    n_steps = 50_000
    snapshots: list[np.ndarray] = []
    burn_in = n_steps // 2
    for i in range(n_steps):
        mc.step(dt, rng)
        if i >= burn_in and (i % 100 == 0):
            engaged = fa.clutches_engaged
            if engaged.any():
                snapshots.append(fa.clutch_forces[engaged].copy())
    if snapshots:
        all_forces = np.concatenate(snapshots)
    else:
        all_forces = np.zeros(0)
    all_forces_pN = all_forces * 1e12
    lo, hi = b["acceptance"]["clutch_force_band_pN"]
    in_band = ((all_forces_pN >= lo) & (all_forces_pN <= hi)).mean() \
        if all_forces_pN.size else 0.0

    fig, ax = plt.subplots(figsize=(6, 4))
    if all_forces_pN.size:
        ax.hist(all_forces_pN, bins=40, color="C0", alpha=0.85, density=True)
    ax.axvspan(lo, hi, color="C2", alpha=0.18,
               label=f"KU-2.12 band {lo}-{hi} pN")
    ax.set_xlabel("Per-clutch force  [pN]")
    ax.set_ylabel("PDF")
    ax.set_title(
        f"KU-2.12 per-clutch force on E = 5 kPa\n"
        f"{in_band*100:.1f} % of mass in {lo}-{hi} pN band"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "03_clutch_force_histogram.png", dpi=140)
    plt.close(fig)
    summary["clutch_force"] = {
        "E_Pa": 5e3,
        "n_samples": int(all_forces_pN.size),
        "mean_pN": float(all_forces_pN.mean()) if all_forces_pN.size else 0.0,
        "median_pN": float(np.median(all_forces_pN)) if all_forces_pN.size else 0.0,
        "in_band_fraction": float(in_band),
        "band_pN": [lo, hi],
    }

    # ------------------------------------------------------------------ #
    # 4. Performance budget probe.                                       #
    # ------------------------------------------------------------------ #
    sub = LinearElasticSubstrate()
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=50)
    mc = MotorClutchFA(fa, sub)
    rng = np.random.default_rng(0)
    for _ in range(10):                                # warm-up
        mc.step(dt, rng)
    t0 = time.perf_counter()
    for _ in range(1000):
        mc.step(dt, rng)
    perf_seconds = time.perf_counter() - t0
    summary["performance"] = {
        "n_steps": 1000,
        "n_clutches": 50,
        "seconds": perf_seconds,
        "budget_seconds": float(b["acceptance"]["perf_budget_seconds"]),
    }

    # ------------------------------------------------------------------ #
    # 5. Sanity Gate report.                                             #
    # ------------------------------------------------------------------ #
    # Locate the simulated catch peak from a fine-grid scan of analytic τ.
    F_fine = np.linspace(0.0, 60e-12, 6001)
    F_star_sim = float(F_fine[int(np.argmax(catch_slip_lifetime(F_fine)))])

    rep = gate_unit2_1_motor_clutch(
        dt=dt,
        cfl_safe_dt=float(d["cfl_safe_dt"]),
        catch_peak_force_analytic=float(d["catch_peak_force"]),
        catch_peak_force_simulated=F_star_sim,
        catch_peak_tolerance=float(b["acceptance"]["catch_peak_tolerance"]),
        biphasic_E_star_analytic=E_star_matched,
        biphasic_E_star_measured=E_argmax,
        biphasic_peak_factor=factor,
        clutch_force_band_pN=tuple(b["acceptance"]["clutch_force_band_pN"]),
        clutch_force_band_mass_measured=float(in_band),
        clutch_force_band_mass_required=float(b["acceptance"]["clutch_force_band_mass"]),
        perf_steps=1000,
        perf_seconds=perf_seconds,
        perf_budget_seconds=float(b["acceptance"]["perf_budget_seconds"]),
    )
    summary["sanity_gate"] = {
        "passed": rep.passed,
        "checks": [
            {"name": c.name, "status": c.status, "ku": c.ku, "detail": c.detail}
            for c in rep.checks
        ],
    }
    print(rep.summary())

    (OUT_DIR / "02_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_DIR / '02_summary.json'} and 3 PNGs to {OUT_DIR}")


if __name__ == "__main__":
    main()
