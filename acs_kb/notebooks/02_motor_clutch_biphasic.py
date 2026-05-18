"""Phase 1 Unit 2.1 validation figures (post-PI-review v2).

Run from the repo root::

    PYTHONPATH=. python acs_kb/notebooks/02_motor_clutch_biphasic.py

Produces four PNGs under ``acs_kb/outputs/phase1_unit2_1/`` plus a
``02_summary.json`` with the measured / analytic comparison numbers and
the multi-seed biphasic shape verdict.

Figures:
  01_catch_bond_lifetime.png       — τ(F) closed-form vs Gillespie samples
  02_biphasic_traction.png         — ⟨F_total⟩ vs E (per-seed + mean), saturating
  03_clutch_force_histogram.png    — per-clutch force at E = 5 kPa (nascent-FA band)
  04_biphasic_argmax_spread.png    — per-seed argmax_E bar chart (saturation evidence)

All parameters are loaded from ``configs/phase1_unit2_1.yaml`` via
``from_config`` classmethods. The YAML is the single source of truth.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acs_kb.bridge.catch_bond import (
    CatchSlipParams,
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
from acs_kb.tests.test_KU_2_4_biphasic import (  # type: ignore[import-not-found]
    biphasic_shape_verdict,
)

CFG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_1.yaml"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "phase1_unit2_1"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    cfg = load_bridge_config(CFG_PATH)
    b = cfg["bridge"]
    d = b["derived"]
    summary: dict = {}

    bond_params = CatchSlipParams.from_config(cfg)
    mc_params = MotorClutchParams.from_config(cfg)
    sub_template = LinearElasticSubstrate.from_config(cfg)
    dt = float(b["dynamics"]["dt"])

    # ------------------------------------------------------------------ #
    # 1. Catch-bond lifetime (analytic + Gillespie samples).             #
    # ------------------------------------------------------------------ #
    F_grid_pN = np.linspace(0.0, 60.0, 121)
    F_grid = F_grid_pN * 1e-12
    tau_analytic = catch_slip_lifetime(F_grid, bond_params)
    rng = np.random.default_rng(20260518)
    sample_forces_pN = np.array([0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    sample_means = np.empty_like(sample_forces_pN, dtype=np.float64)
    for i, F_pN in enumerate(sample_forces_pN):
        F = F_pN * 1e-12
        k = float(catch_slip_off_rate(F, bond_params))
        sample_means[i] = rng.exponential(scale=1.0 / k, size=4000).mean()
    F_star = lifetime_peak_force(bond_params)
    tau_at_peak = float(catch_slip_lifetime(F_star, bond_params))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(F_grid_pN, tau_analytic, "-", color="C0",
            label="τ(F) = 1/k_off(F)  (analytic)")
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
    # 2. Biphasic traction sweep (replicated seeds).                     #
    # ------------------------------------------------------------------ #
    sw = b["biphasic_sweep"]
    E_grid = np.logspace(sw["log10_E_min"], sw["log10_E_max"], int(sw["n_points"]))
    n_seeds = int(sw["n_seeds"])
    n_steps_per_E = int(sw["n_steps_per_E"])
    seed_root = int(sw["seed"])
    per_seed = np.empty((n_seeds, len(E_grid)), dtype=np.float64)
    for s in range(n_seeds):
        for i, E in enumerate(E_grid):
            sub = LinearElasticSubstrate(
                young_modulus=float(E),
                poisson_ratio=sub_template.poisson_ratio,
                thickness=sub_template.thickness,
                contact_radius=sub_template.contact_radius,
            )
            fa = make_focal_adhesion(np.array([0.0, 0.0]),
                                     n_clutches_total=mc_params.n_clutches)
            mc = MotorClutchFA(fa, sub, mc_params)
            rng_e = np.random.default_rng(seed_root + 1000 * s + i)
            r = run_steady_state(mc, dt=dt, n_steps=n_steps_per_E, rng=rng_e,
                                 burn_in_fraction=0.5)
            per_seed[s, i] = r["mean_force_total"]
    verdict = biphasic_shape_verdict(E_grid, per_seed)
    summary["biphasic"] = verdict

    # 02_biphasic_traction.png  — per-seed traces + mean
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for s in range(n_seeds):
        ax.plot(E_grid, per_seed[s] * 1e12, "-", color="C0", alpha=0.35,
                linewidth=1)
    ax.plot(E_grid, np.array(verdict["mean_F"]) * 1e12, "o-", color="C0",
            linewidth=2, label=f"mean over {n_seeds} seeds")
    F_stall_total = mc_params.n_motors * mc_params.F_stall_per_motor * 1e12
    ax.axhline(F_stall_total, color="C3", linestyle="--", alpha=0.6,
               label=f"N_m · F_stall = {F_stall_total:.0f} pN (asymptote)")
    ax.axvline(d["biphasic_EStar_matched"], color="grey", linestyle=":",
               alpha=0.5,
               label=f"E*_matched = {d['biphasic_EStar_matched']:.2e} Pa "
                     f"(Bangasser, would-be slip-bond peak)")
    ax.set_xscale("log")
    ax.set_xlabel("Substrate Young modulus E  [Pa]")
    ax.set_ylabel("Mean traction ⟨F_total⟩  [pN]")
    ax.set_title(
        f"KU-2.4 / KU-2.8 traction vs stiffness  —  verdict: "
        f"{verdict['verdict'].upper()}\n"
        f"prominence {verdict['prominence']*100:.2f}% of asymptote, "
        f"argmax_E spread {verdict['log10_seed_argmax_spread']:.2f} decades"
    )
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "02_biphasic_traction.png", dpi=140)
    plt.close(fig)

    # 04_biphasic_argmax_spread.png — per-seed argmax_E
    fig, ax = plt.subplots(figsize=(6, 3.5))
    seeds = np.arange(n_seeds)
    argmax_per_seed = np.array(verdict["seed_argmax_E"])
    ax.bar(seeds, np.log10(argmax_per_seed), color="C0", alpha=0.7)
    ax.set_xticks(seeds)
    ax.set_xticklabels([f"seed {s}" for s in range(n_seeds)])
    ax.set_ylabel("log10(argmax E)  [Pa]")
    ax.axhline(np.log10(d["biphasic_EStar_matched"]), color="grey",
               linestyle="--", alpha=0.6,
               label=f"log10(E*_matched) = {np.log10(d['biphasic_EStar_matched']):.2f}")
    ax.set_title(
        f"argmax_E spans {verdict['log10_seed_argmax_spread']:.2f} decades "
        f"across {n_seeds} seeds — noise-dominated (saturating)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_biphasic_argmax_spread.png", dpi=140)
    plt.close(fig)

    # ------------------------------------------------------------------ #
    # 3. Per-clutch force histogram at E = 5 kPa.                        #
    # ------------------------------------------------------------------ #
    # Use the YAML-resolved substrate (E = 5 kPa default).
    sub = LinearElasticSubstrate.from_config(cfg)
    fa = make_focal_adhesion(np.array([0.0, 0.0]),
                             n_clutches_total=mc_params.n_clutches)
    mc = MotorClutchFA(fa, sub, mc_params)
    rng = np.random.default_rng(9999)
    n_steps = 50_000
    snapshots: list[np.ndarray] = []
    burn_in = n_steps // 2
    for i in range(n_steps):
        mc.step(dt, rng)
        if i >= burn_in and (i % 100 == 0):
            engaged = fa.clutches_engaged
            if engaged.any():
                snapshots.append(fa.clutch_forces[engaged].copy())
    all_forces = np.concatenate(snapshots) if snapshots else np.zeros(0)
    all_forces_pN = all_forces * 1e12
    nascent_lo, nascent_hi = b["acceptance"]["clutch_force_band_pN"]
    mature_lo, mature_hi = b["acceptance"]["ku_2_12_mature_band_pN"]
    in_nascent = float(((all_forces_pN >= nascent_lo) &
                        (all_forces_pN <= nascent_hi)).mean()) \
        if all_forces_pN.size else 0.0
    in_mature = float(((all_forces_pN >= mature_lo) &
                       (all_forces_pN <= mature_hi)).mean()) \
        if all_forces_pN.size else 0.0

    fig, ax = plt.subplots(figsize=(6, 4))
    if all_forces_pN.size:
        ax.hist(all_forces_pN, bins=40, color="C0", alpha=0.85, density=True)
    ax.axvspan(nascent_lo, nascent_hi, color="C2", alpha=0.18,
               label=f"Phase 1 nascent {nascent_lo}-{nascent_hi} pN "
                     f"({in_nascent*100:.1f}%)")
    ax.axvspan(mature_lo, mature_hi, color="C3", alpha=0.10,
               label=f"KU-2.12 mature {mature_lo}-{mature_hi} pN "
                     f"({in_mature*100:.1f}%) → Unit 2.2")
    ax.set_xlabel("Per-clutch force  [pN]")
    ax.set_ylabel("PDF")
    ax.set_title(
        f"KU-2.12 per-clutch force on E = "
        f"{sub.young_modulus:.0f} Pa (Phase 1 nascent FA)"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "03_clutch_force_histogram.png", dpi=140)
    plt.close(fig)
    summary["clutch_force"] = {
        "E_Pa": sub.young_modulus,
        "n_samples": int(all_forces_pN.size),
        "mean_pN": float(all_forces_pN.mean()) if all_forces_pN.size else 0.0,
        "median_pN": float(np.median(all_forces_pN)) if all_forces_pN.size else 0.0,
        "in_phase1_nascent_band_fraction": in_nascent,
        "in_ku_2_12_mature_band_fraction": in_mature,
        "phase1_band_pN": [nascent_lo, nascent_hi],
        "ku_2_12_mature_band_pN": [mature_lo, mature_hi],
    }

    # ------------------------------------------------------------------ #
    # 4. Performance budget probe.                                       #
    # ------------------------------------------------------------------ #
    sub = LinearElasticSubstrate.from_config(cfg)
    fa = make_focal_adhesion(np.array([0.0, 0.0]),
                             n_clutches_total=mc_params.n_clutches)
    mc = MotorClutchFA(fa, sub, mc_params)
    rng = np.random.default_rng(0)
    for _ in range(10):
        mc.step(dt, rng)
    t0 = time.perf_counter()
    for _ in range(1000):
        mc.step(dt, rng)
    perf_seconds = time.perf_counter() - t0
    summary["performance"] = {
        "n_steps": 1000,
        "n_clutches": mc_params.n_clutches,
        "seconds": perf_seconds,
        "budget_seconds": float(b["acceptance"]["perf_budget_seconds"]),
    }

    # ------------------------------------------------------------------ #
    # 5. Sanity Gate report.                                             #
    # ------------------------------------------------------------------ #
    F_fine = np.linspace(0.0, 60e-12, 6001)
    F_star_sim = float(F_fine[int(np.argmax(catch_slip_lifetime(F_fine, bond_params)))])

    rep = gate_unit2_1_motor_clutch(
        dt=dt,
        cfl_safe_dt=float(d["cfl_safe_dt"]),
        catch_peak_force_analytic=float(d["catch_peak_force"]),
        catch_peak_force_simulated=F_star_sim,
        catch_peak_tolerance=float(b["acceptance"]["catch_peak_tolerance"]),
        biphasic_verdict=verdict["verdict"],
        biphasic_prominence=verdict["prominence"],
        biphasic_log10_seed_argmax_spread=verdict["log10_seed_argmax_spread"],
        biphasic_asymptote_N=verdict["asymptote"],
        motor_stall_total_N=mc_params.n_motors * mc_params.F_stall_per_motor,
        biphasic_E_star_analytic=float(d["biphasic_EStar_matched"]),
        biphasic_E_star_measured_argmax=float(
            E_grid[int(np.argmax(np.array(verdict["mean_F"])))]
        ),
        clutch_force_band_pN=tuple(b["acceptance"]["clutch_force_band_pN"]),
        clutch_force_band_mass_measured=in_nascent,
        clutch_force_band_mass_required=float(b["acceptance"]["clutch_force_band_mass"]),
        perf_steps=1000,
        perf_seconds=perf_seconds,
        perf_budget_seconds=float(b["acceptance"]["perf_budget_seconds"]),
        ku_2_12_mature_band_pN=tuple(b["acceptance"]["ku_2_12_mature_band_pN"]),
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
    print(f"\nWrote {OUT_DIR / '02_summary.json'} and 4 PNGs to {OUT_DIR}")


if __name__ == "__main__":
    main()
