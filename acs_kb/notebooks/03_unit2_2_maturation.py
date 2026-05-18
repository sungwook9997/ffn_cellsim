"""Phase 1 Unit 2.2 validation figures.

Run from the repo root::

    PYTHONPATH=. python acs_kb/notebooks/03_unit2_2_maturation.py

Produces four PNGs under ``acs_kb/outputs/phase1_unit2_2/`` plus a
``03_summary.json``.

Figures:
  01_talin_unfold_rate.png         — Bell-Evans rate vs force
  02_vinculin_recruitment.png      — N_vin(t) trajectory at n_unfolded = 1
  03_mature_force_histogram.png    — per-clutch force histogram, mature regime
  04_fa_growth_response.png        — FA area trajectory under varying F_total
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acs_kb.bridge.fa_growth import (
    DEFAULT_PARAMS as DEFAULT_FA_GROWTH,
    FAGrowthParams,
    fa_growth_step,
    hill_growth_rate,
)
from acs_kb.bridge.motor_clutch import (
    MotorClutchFA,
    MotorClutchParams,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.talin import (
    DEFAULT_PARAMS as DEFAULT_TALIN,
    TalinParams,
    talin_unfold_rate,
)
from acs_kb.bridge.types import make_focal_adhesion
from acs_kb.bridge.vinculin import (
    DEFAULT_PARAMS as DEFAULT_VIN,
    VinculinParams,
    vinculin_recruit_step,
)
from acs_kb.common.derived_params import load_bridge_22_config

CFG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_2.yaml"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "phase1_unit2_2"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    cfg = load_bridge_22_config(CFG_PATH)
    b = cfg["bridge"]
    d = b["derived"]
    summary: dict = {}

    talin_params = TalinParams.from_config(cfg)
    vin_params = VinculinParams.from_config(cfg)
    growth_params = FAGrowthParams.from_config(cfg)
    mc_params = MotorClutchParams.from_config(cfg)
    sub = LinearElasticSubstrate.from_config(cfg)
    dt = float(b["dynamics"]["dt"])

    # ------------------------------------------------------------------ #
    # 1. Talin Bell-Evans rate                                           #
    # ------------------------------------------------------------------ #
    F_pN = np.linspace(0.0, 30.0, 121)
    k_u = np.asarray(talin_unfold_rate(F_pN * 1e-12, talin_params))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(F_pN, k_u, "-", color="C0",
            label="k_u(F) = k_{u,0}·exp(F Δx*/kT)")
    ax.axhline(talin_params.k_u0, color="grey", linestyle="--", alpha=0.5,
               label=f"k_u0 = {talin_params.k_u0:.3f} s⁻¹")
    ax.set_yscale("log")
    ax.set_xlabel("Force F  [pN]")
    ax.set_ylabel("Talin unfold rate k_u  [s⁻¹]  (log)")
    ax.set_title(
        f"KU-2.6 talin Bell-Evans (Δx* = {talin_params.dx_star*1e9:.2f} nm)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "01_talin_unfold_rate.png", dpi=140)
    plt.close(fig)
    summary["talin"] = {
        "k_u0_s_inv": talin_params.k_u0,
        "dx_star_m": talin_params.dx_star,
        "k_unfold_at_5pN_s_inv": float(talin_unfold_rate(5e-12, talin_params)),
        "k_unfold_at_15pN_s_inv": float(talin_unfold_rate(15e-12, talin_params)),
    }

    # ------------------------------------------------------------------ #
    # 2. Vinculin recruitment trajectory                                 #
    # ------------------------------------------------------------------ #
    fa = make_focal_adhesion(np.array([0.0, 0.0]), n_clutches_total=20)
    fa.talin_unfolded_domains = 1
    rng = np.random.default_rng(42)
    n_steps = 30000
    dt_v = 1e-2
    trace = np.empty(n_steps, dtype=np.int32)
    for i in range(n_steps):
        vinculin_recruit_step(fa, dt_v, vin_params, rng)
        trace[i] = fa.vinculin_count
    t = np.arange(n_steps) * dt_v
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, trace, "-", color="C0", linewidth=1)
    ss = vin_params.steady_state_count(1)
    ax.axhline(ss, color="C3", linestyle="--", alpha=0.6,
               label=f"N_ss = k_rec·1·N_free/k_diss = {ss:.0f}")
    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("Vinculin count")
    ax.set_title("KU-2.7 vinculin recruitment (Poisson tau-leap)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "02_vinculin_recruitment.png", dpi=140)
    plt.close(fig)
    summary["vinculin"] = {
        "steady_state_count_analytic": ss,
        "measured_late_mean": float(trace[n_steps // 2:].mean()),
        "measured_late_std": float(trace[n_steps // 2:].std()),
    }

    # ------------------------------------------------------------------ #
    # 3. Mature regime per-clutch force histogram                        #
    # ------------------------------------------------------------------ #
    n_clutches_mature = int(b["acceptance"]["mature_test_n_clutches"])
    fa = make_focal_adhesion(np.array([0.0, 0.0]),
                             n_clutches_total=n_clutches_mature)
    fa.talin_unfolded_domains = 1
    mc = MotorClutchFA(fa, sub, mc_params,
                       talin_params=talin_params, vinculin_params=vin_params,
                       fa_growth_params=None)
    rng = np.random.default_rng(int(b["dynamics"]["seed"]))
    n_steps_m = int(b["dynamics"]["n_steps_mature"])
    snapshots: list[np.ndarray] = []
    burn_in = n_steps_m // 2
    for i in range(n_steps_m):
        mc.step(dt, rng)
        if i >= burn_in and (i % 200 == 0):
            engaged = fa.clutches_engaged
            if engaged.any():
                snapshots.append(fa.clutch_forces[engaged].copy())
    forces_pN = (np.concatenate(snapshots) * 1e12) if snapshots else np.zeros(0)
    lo, hi = b["acceptance"]["mature_band_pN"]
    in_band = ((forces_pN >= lo) & (forces_pN <= hi)).mean() \
        if forces_pN.size else 0.0
    fig, ax = plt.subplots(figsize=(6, 4))
    if forces_pN.size:
        ax.hist(forces_pN, bins=40, color="C0", alpha=0.85, density=True)
    ax.axvspan(lo, hi, color="C2", alpha=0.18,
               label=f"KU-2.12 mature {lo}-{hi} pN ({in_band*100:.1f}%)")
    ax.set_xlabel("Per-clutch force  [pN]")
    ax.set_ylabel("PDF")
    ax.set_title(
        f"KU-2.12 mature regime  —  N_clutches={n_clutches_mature}, "
        f"vinculin active\n"
        f"mean {forces_pN.mean():.2f} pN, median {float(np.median(forces_pN)):.2f} pN"
    )
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "03_mature_force_histogram.png", dpi=140)
    plt.close(fig)
    summary["mature_regime"] = {
        "n_clutches": n_clutches_mature,
        "n_samples": int(forces_pN.size),
        "mean_pN": float(forces_pN.mean()) if forces_pN.size else 0.0,
        "median_pN": float(np.median(forces_pN)) if forces_pN.size else 0.0,
        "in_band_fraction": float(in_band),
        "mature_band_pN": [lo, hi],
    }

    # ------------------------------------------------------------------ #
    # 4. FA growth response — Hill kinetics                              #
    # ------------------------------------------------------------------ #
    forces_test = np.array([0.0, 25e-12, 50e-12, 75e-12, 100e-12, 200e-12])
    times = np.arange(0.0, 60.0 + dt, 1e-2)
    fig, ax = plt.subplots(figsize=(6, 4))
    for F_const in forces_test:
        fa = make_focal_adhesion(np.array([0.0, 0.0]),
                                 n_clutches_total=growth_params.n_clutches_nascent,
                                 area=growth_params.A_nascent)
        traj = np.empty(len(times), dtype=np.float64)
        for j, _ in enumerate(times):
            traj[j] = fa.area
            fa_growth_step(fa, F_const, 1e-2, growth_params)
        ax.plot(times, traj / growth_params.A_nascent,
                label=f"F = {F_const*1e12:.0f} pN")
    ax.set_xlabel("Time  [s]")
    ax.set_ylabel("FA area / A_nascent")
    ax.axhline(1.0, color="grey", linestyle=":", alpha=0.4)
    ax.set_title(
        f"KU-2.17 Hill-function FA growth  "
        f"(F_th = {growth_params.F_th*1e12:.0f} pN, n = {growth_params.n_hill})"
    )
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "04_fa_growth_response.png", dpi=140)
    plt.close(fig)
    summary["fa_growth"] = {
        "F_th_pN": growth_params.F_th * 1e12,
        "k_g0_s_inv": growth_params.k_g0,
        "k_d_s_inv": growth_params.k_d,
        "kg_at_F_th_s_inv": hill_growth_rate(growth_params.F_th, growth_params),
        "kg_at_2_F_th_s_inv": hill_growth_rate(2 * growth_params.F_th, growth_params),
    }

    summary["derived"] = {k: v for k, v in d.items() if v is not None}

    (OUT_DIR / "03_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"Wrote 4 PNGs and 03_summary.json to {OUT_DIR}")


if __name__ == "__main__":
    main()
