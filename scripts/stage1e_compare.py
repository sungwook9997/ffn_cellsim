"""Stage 1e Sim A vs Sim B comparison driver.

Reads the existing Stage 1d 3D pilot result (Sim B), runs the analytical
Sim A 1D radial ODE on the same time grid, fits A/A₀ = a + b/R + c/R² to
both, computes comparison metrics, classifies into Bucket E1/E2/E3/E4 per
`docs/outcomes_stage1e.md`, and writes a report.

Usage:

    python scripts/stage1e_compare.py
        [--sim_b_dir results/stage1d_pilot]
        [--sim_b_config configs/stage1d_pilot.yaml]
        [--out_dir results/stage1e_comparison]

Default values match the Stage 1e default scope (single Pre carrier).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import yaml

from acs.analysis.radial_reduction import (
    SimAParams,
    compare_sim_a_vs_b,
    fit_radial_ansatz,
    simulate_sim_a,
)


def _load_sim_b_trajectory(metrics_csv: Path) -> dict:
    """Load Sim B (Stage 1d pilot) trajectory from metrics.csv.

    Reconciles A/A₀ measurement to match Sim A's circular-area
    convention (per docs/stage1e_sanity.md check 6): use
    R_eff_xy = √(A_hull / π) as the effective spreading radius.
    """
    rows = []
    with metrics_csv.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)

    t = np.array([float(r["time_star"]) for r in rows])
    A_hull = np.array([
        float(r.get("contact_area_xy_hull", "nan") or "nan") for r in rows
    ])

    # R_eff_xy from contact area (per Stage 1e sanity-md §6 reconciliation).
    # Where A_hull is missing/NaN, fall back to NaN.
    valid = np.isfinite(A_hull) & (A_hull > 0)
    R_eff_xy = np.full_like(A_hull, np.nan)
    R_eff_xy[valid] = np.sqrt(A_hull[valid] / np.pi)

    # A/A₀ normalised to first valid frame.
    A_over_A0 = np.full_like(A_hull, np.nan)
    A_init = next((a for a in A_hull if np.isfinite(a) and a > 0), np.nan)
    if np.isfinite(A_init) and A_init > 0:
        A_over_A0[valid] = A_hull[valid] / A_init

    return {
        "t": t,
        "R": R_eff_xy,
        "A_over_A0": A_over_A0,
        "A_hull": A_hull,
    }


def _params_from_yaml(cfg_path: Path) -> SimAParams:
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    nd = cfg["nondim"]
    layer3 = cfg.get("layer3", {})
    layer4 = cfg.get("layer4", {})
    return SimAParams(
        R0_star=float(nd["radius_star"]),
        gamma_star=(
            float(nd["capillary_number"]) * float(nd["K_star"]) * float(nd["radius_star"])
        ),
        xi_star=float(nd["drag_xi_star"]),
        zeta_min=float(layer3.get("zeta_min", 0.0)),
        zeta_max=float(layer3.get("zeta_max", 0.0)),
        k_plus_star=float(layer3.get("k_plus_star", 0.0)),
        k_minus_star=float(layer3.get("k_minus_star", 0.0)),
        phi_initial=float(layer3.get("phi_initial", 0.05)),
        layer3_spatial_S=bool(layer4.get("layer3_spatial_S", False)),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Stage 1e Sim A vs Sim B comparison")
    ap.add_argument("--sim_b_dir", default="results/stage1d_pilot")
    ap.add_argument("--sim_b_config", default="configs/stage1d_pilot.yaml")
    ap.add_argument("--out_dir", default="results/stage1e_comparison")
    args = ap.parse_args()

    sim_b_dir = Path(args.sim_b_dir)
    cfg_path = Path(args.sim_b_config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_csv = sim_b_dir / "metrics.csv"
    if not metrics_csv.exists():
        print(f"ERROR: Sim B metrics CSV not found: {metrics_csv}")
        return 2
    if not cfg_path.exists():
        print(f"ERROR: Sim B config not found: {cfg_path}")
        return 2

    print(f"=== Stage 1e: Sim A (1D radial ODE) vs Sim B ({sim_b_dir.name}) ===")

    sim_b = _load_sim_b_trajectory(metrics_csv)
    print(f"Sim B loaded: {len(sim_b['t'])} frames, t ∈ [{sim_b['t'][0]:.1f}, {sim_b['t'][-1]:.1f}]")

    # Sim A on the same time grid.
    params = _params_from_yaml(cfg_path)
    print(
        f"Sim A params: R0={params.R0_star:.3f}, γ={params.gamma_star:.4e}, "
        f"ξ={params.xi_star:.3f}, ζ∈[{params.zeta_min:.3f}, {params.zeta_max:.3f}], "
        f"φ_init={params.phi_initial:.3f}, k_+={params.k_plus_star:.4e}, "
        f"k_-={params.k_minus_star:.4e}"
    )

    try:
        sim_a = simulate_sim_a(params, sim_b["t"])
    except Exception as exc:
        print(f"BUCKET E4 — Sim A integration failed: {exc}")
        return 1
    print(
        f"Sim A integrated: R(0)={sim_a['R'][0]:.4f}, R(end)={sim_a['R'][-1]:.4f}, "
        f"A/A0(end)={sim_a['A_over_A0'][-1]:.4f}"
    )

    # Fit A/A₀ = a + b/R + c/R² to BOTH.
    try:
        fit_a = fit_radial_ansatz(sim_a["R"], sim_a["A_over_A0"])
        print(
            f"Sim A fit: a={fit_a['a']:.4f}, b={fit_a['b']:.4f}, c={fit_a['c']:.4f}, "
            f"R²={fit_a['r_squared']:.4f}, RMSE={fit_a['rmse']:.4f}"
        )
    except Exception as exc:
        print(f"BUCKET E4 — Sim A radial fit failed: {exc}")
        return 1

    try:
        fit_b = fit_radial_ansatz(sim_b["R"], sim_b["A_over_A0"])
        print(
            f"Sim B fit: a={fit_b['a']:.4f}, b={fit_b['b']:.4f}, c={fit_b['c']:.4f}, "
            f"R²={fit_b['r_squared']:.4f}, RMSE={fit_b['rmse']:.4f}"
        )
    except Exception as exc:
        print(f"BUCKET E4 — Sim B radial fit failed: {exc}")
        return 1

    # Comparison.
    cmp = compare_sim_a_vs_b(sim_a, sim_b, fit_a, fit_b)
    print(
        f"Comparison: bucket={cmp['bucket']}, "
        f"a_dev={cmp['a_rel_dev']:.3f}, b_dev={cmp['b_rel_dev']:.3f}, "
        f"c_dev={cmp['c_rel_dev']:.3f}, RMS_norm={cmp['trajectory_rms_normalised']:.3f}, "
        f"Pearson={cmp['pearson_correlation']:.3f}"
    )

    # Write fits.json
    fits = {
        "sim_a_params": {
            "R0_star": params.R0_star, "gamma_star": params.gamma_star,
            "xi_star": params.xi_star,
            "zeta_min": params.zeta_min, "zeta_max": params.zeta_max,
            "k_plus_star": params.k_plus_star, "k_minus_star": params.k_minus_star,
            "phi_initial": params.phi_initial,
        },
        "sim_a_fit": fit_a,
        "sim_b_fit": fit_b,
        "comparison": cmp,
        "sim_a_trajectory": {
            "t": sim_a["t"].tolist(),
            "R": sim_a["R"].tolist(),
            "A_over_A0": sim_a["A_over_A0"].tolist(),
            "phi": sim_a["phi"].tolist(),
            "zeta_eff": sim_a["zeta_eff"].tolist(),
        },
        "sim_b_trajectory": {
            "t": sim_b["t"].tolist(),
            "R": [float(x) if np.isfinite(x) else None for x in sim_b["R"]],
            "A_over_A0": [float(x) if np.isfinite(x) else None for x in sim_b["A_over_A0"]],
        },
    }
    (out_dir / "fits.json").write_text(json.dumps(fits, indent=2), encoding="utf-8")

    # Write report.md
    report = [
        f"# Stage 1e — Sim A vs Sim B comparison report",
        "",
        f"- Sim B source: `{sim_b_dir.name}` ({len(sim_b['t'])} frames, t ∈ [0, {sim_b['t'][-1]:.1f}]·τ_relax)",
        f"- Sim A: 1D radial ODE reduction on the same time grid",
        "",
        f"## Bucket classification: **{cmp['bucket']}**",
        "",
        "Per `docs/outcomes_stage1e.md`:",
        "- E1: Sim A reproduces Sim B (radial reduction valid) — all (a,b,c) ≤ 25% deviation, RMS ≤ 0.10, Pearson ≥ 0.90",
        "- E2: partial agreement — at least one of (a,b,c) within 25% but at least one > 50%, OR RMS in [0.10, 0.30]",
        "- E3: large difference (anisotropy / 3D-only physics important) — RMS > 0.30 OR Pearson < 0.50 OR all (a,b,c) > 50%",
        "- E4: Sim A or Sim B failed (NaN, fit divergent)",
        "",
        "## Radial-fit (a, b, c) parameters",
        "",
        f"| Coeff | Sim A     | Sim B     | |A−B|/|B| |",
        f"|---|---|---|---|",
        f"| a | {fit_a['a']:+.4f} | {fit_b['a']:+.4f} | {cmp['a_rel_dev']:.3f} |",
        f"| b | {fit_a['b']:+.4f} | {fit_b['b']:+.4f} | {cmp['b_rel_dev']:.3f} |",
        f"| c | {fit_a['c']:+.4f} | {fit_b['c']:+.4f} | {cmp['c_rel_dev']:.3f} |",
        "",
        f"- Sim A fit R² = {fit_a['r_squared']:.4f}, RMSE = {fit_a['rmse']:.4f}, n={fit_a['n_samples']}",
        f"- Sim B fit R² = {fit_b['r_squared']:.4f}, RMSE = {fit_b['rmse']:.4f}, n={fit_b['n_samples']}",
        "",
        "## A/A₀ trajectory comparison",
        "",
        f"- RMS deviation (Sim A − Sim B) = {cmp['trajectory_rms']:.4f}",
        f"- Normalised RMS = {cmp['trajectory_rms_normalised']:.4f} (RMS / mean A_B)",
        f"- Pearson correlation = {cmp['pearson_correlation']:.4f}",
        f"- n valid samples = {cmp['n_valid_samples']}",
        "",
        "## Mechanism interpretation",
        "",
        "Sim A captures: surface-tension recovery (κ=2/R), Layer 2/3 active stress",
        "ζ_eff(t) via Layer 3 single-particle ODE limit, effective drag ξ·R²,",
        "Layer 5 ρ_osm baseline. Sim A does NOT capture: Path C effective gravity,",
        "Layer 4 Marangoni tangential surface flow, anisotropic spreading, spatial",
        "heterogeneity beyond mean-field, full 3D substrate-anchoring dynamics.",
        "",
        "Bucket-E classification quantifies the validity gap.",
        "",
        "## Files",
        "",
        f"- `fits.json` — all (a,b,c), trajectory data, comparison metrics",
        f"- `report.md` — this report",
        "",
        "Stop conditions per `docs/outcomes_stage1e.md`. Stage 2 / 1e.b / inherited",
        "fixes auto-entry FORBIDDEN; PI input required.",
    ]
    (out_dir / "report.md").write_text("\n".join(report), encoding="utf-8")

    print(f"\nReport written to: {out_dir / 'report.md'}")
    print(f"Fits JSON written to: {out_dir / 'fits.json'}")
    print(f"Stage 1e Bucket: {cmp['bucket']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
