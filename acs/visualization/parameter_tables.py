"""acs.visualization.parameter_tables — production parameter summary.

Per PI visualization directive 2026-04-29: each run gets a parameter
summary table in BOTH CSV and figure form. Includes the alias
column (poster_label, code_name, value, unit, source) for cross-
reference with PI poster/manuscript terminology.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml


# (poster_label, code_name dotted-path, default unit, source)
PARAM_TABLE_SPEC = [
    ("Cell line", "_const:MCF7", "—", "PI experimental setup"),
    ("Formation condition", "_phenotype", "—", "config layer3.phi_initial → Bare/Pre/Lam4"),
    ("Spreading substrate", "_const:Col1-coated dish", "—", "PI experimental setup"),
    ("Initial spheroid radius", "initialization.radius_um", "μm", "config"),
    ("Domain size", "simulation.domain_size_um", "μm", "config"),
    ("Grid resolution dx", "simulation.dx_um", "μm", "config"),
    ("Time step dt", "simulation.dt_s", "s", "config"),
    ("Total sim time", "simulation.total_sim_time_s", "s", "config"),
    ("Frame interval", "simulation.frame_interval_s", "s", "config"),
    ("n_particles", "simulation.n_material_points", "#", "config"),
    ("Density ρ", "physics.density_kg_m3", "kg/m³", "config"),
    ("Cortical bulk K", "physics.cortical_K_Pa", "Pa", "Fischer-Friedrich 2014, Krieg 2008"),
    ("Shear modulus μ", "physics.shear_mu_Pa", "Pa", "K, μ/K dimensional consistency"),
    ("Maxwell τ", "physics.maxwell_tau_s", "s", "Moeendarbary 2013 *Nat Mater* IF 47"),
    ("Cohesive γ_cc", "physics.cohesive_gamma_J_m2", "J/m²", "Foty & Steinberg 2005 *Dev Biol*"),
    ("Substrate enabled", "substrate.enabled", "—", "config"),
    ("γ_sub_α (Option β)", "substrate.gamma_sub_alpha", "—", "Stage 1a+ Option β PI sweep"),
    ("ζ_min (Layer 3)", "layer3.zeta_min", "—", "stage1a_plus_plus_layer2_sanity.md"),
    ("ζ_max (Layer 3)", "layer3.zeta_max", "—", "Phase 2.2 sweep first-PASS"),
    ("φ_initial", "layer3.phi_initial", "—", "Cho 2020 + phenotype mapping"),
    ("k_+ (Layer 3)", "layer3.k_plus_star", "1/τ", "Cho 2020 Fig. 4 slope"),
    ("k_- (Layer 3)", "layer3.k_minus_star", "1/τ", "Cho 2020 reverse rate"),
    ("Layer 3 split", "layer3.split", "—", "Stage 1b.b (PI 2026-04-29)"),
    ("κ contact-act amplifier", "layer3.kappa_act", "—", "Stage 1b.b PI default"),
    ("ε memory decay", "layer3.memory_eps_star", "1/τ", "Stage 1b.b PI default 0"),
    ("γ_max (Marangoni)", "layer4.gamma_max_star", "—", "Pajic-Lijakovic 2022, Maître 2012"),
    ("γ_min (Marangoni)", "layer4.gamma_min_star", "—", "Pajic-Lijakovic 2022, Maître 2012"),
    ("Dynamic γ (Mech A/F)", "layer4.dynamic_gamma", "—", "Stage 1d.b"),
    ("τ_γ relaxation", "layer4.tau_gamma_star", "—", "Yadav 2022 *PRF* L031101"),
    ("α_A reaccumulation", "layer4.alpha_A_star", "—", "Stage 1d.b PI default"),
    ("α_F osmotic coupling", "layer4.alpha_F_star", "—", "Yadav 2022 + Guo 2017"),
    ("Layer 2.b stochastic events", "layer2_b.enabled", "—", "Stage 1a++.b"),
    ("λ_lam (event rate)", "layer2_b.lambda_lam_star", "1/τ", "Mattila & Lappalainen 2008 *NRMCB*"),
    ("Lam impulse", "layer2_b.impulse_lam_star", "v*", "Stage 1a++.b dimensional"),
    ("ρ_osm initial", "layer5.rho_osm_initial", "—", "Guo 2017 *PNAS*"),
    ("α_osm", "layer5.alpha_osm_star", "1/τ", "Guo 2017"),
    ("β_osm", "layer5.beta_osm_star", "1/τ", "Guo 2017"),
    ("ρ_osm range", "layer5.rho_osm_min/max", "—", "Guo 2017 + sanity_md"),
    ("MMP α (Layer 6)", "layer6.alpha_mmp_star", "1/τ", "Egeblad-Werb 2002 *NRC* IF 78"),
    ("ECM degradation β", "layer6.beta_deg_star", "1/τ", "Egeblad-Werb 2002"),
    ("ecm_strength range", "layer6.ecm_strength_initial/min", "—", "Stage 2 sanity_md"),
    ("Gravity (Path C)", "gravity.gravity_star", "—", "Path C dimensional"),
    ("GPU backend", "gpu.backend", "—", "config"),
]


def _resolve(cfg: Dict[str, Any], path: str) -> Any:
    """Resolve a dotted path with support for special markers."""
    if path.startswith("_const:"):
        return path[len("_const:"):]
    if path == "_phenotype":
        phi_init = cfg.get("layer3", {}).get("phi_initial", None)
        if phi_init is None:
            return "?"
        if phi_init <= 0.40:
            return "Bare"
        elif phi_init <= 0.65:
            return "Pre"
        else:
            return "Lam4"
    parts = path.split(".")
    # Support "a.b/c" notation: take parent path, then return a/b joined values.
    if "/" in parts[-1]:
        tail = parts[-1]
        keys = tail.split("/")
        node = cfg
        for p in parts[:-1]:
            if not isinstance(node, dict):
                return "?"
            node = node.get(p, {})
        if not isinstance(node, dict):
            return "?"
        vals = [node.get(k, "?") for k in keys]
        return " / ".join(str(v) for v in vals)
    node: Any = cfg
    for p in parts:
        if not isinstance(node, dict):
            return "?"
        node = node.get(p, "?")
    return node


def write_parameter_table_csv(run_dir: Path, out_path: Optional[Path] = None) -> Path:
    """Write a CSV row-per-parameter table for the run."""
    run_dir = Path(run_dir)
    cfg_path = run_dir / "config.yaml"
    cfg: Dict[str, Any] = {}
    if cfg_path.exists():
        with open(cfg_path, encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
    if out_path is None:
        out_path = run_dir / "figures" / "parameter_tables" / "run_parameters.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Optional run identity rows.
    run_name = cfg.get("run", {}).get("name", run_dir.name)
    seed = cfg.get("run", {}).get("seed", "?")
    rows = [
        {"poster_label": "Run name", "code_name": "run.name",
         "value": run_name, "unit": "—", "source": "config"},
        {"poster_label": "Random seed", "code_name": "run.seed",
         "value": str(seed), "unit": "—", "source": "config"},
    ]
    for poster_label, code_name, unit, source in PARAM_TABLE_SPEC:
        val = _resolve(cfg, code_name)
        rows.append({
            "poster_label": poster_label,
            "code_name": code_name,
            "value": str(val),
            "unit": unit,
            "source": source,
        })
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "poster_label", "code_name", "value", "unit", "source"
        ])
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def render_parameter_table_figure(run_dir: Path, out_path: Optional[Path] = None) -> Path:
    """Render the parameter table as a matplotlib figure (poster-ready)."""
    run_dir = Path(run_dir)
    csv_path = write_parameter_table_csv(run_dir)
    rows: List[Dict[str, str]] = []
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if out_path is None:
        out_path = run_dir / "figures" / "parameter_tables" / "run_parameters.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig_h = max(4.0, 0.30 * len(rows) + 1.0)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")
    table_data = [[r["poster_label"], r["code_name"], r["value"],
                   r["unit"], r["source"]] for r in rows]
    table = ax.table(
        cellText=table_data,
        colLabels=["poster_label", "code_name", "value", "unit", "source"],
        cellLoc="left",
        loc="upper left",
        colWidths=[0.22, 0.22, 0.18, 0.10, 0.28],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.2)
    # Header styling
    for i, key in enumerate(["poster_label", "code_name", "value", "unit", "source"]):
        cell = table[(0, i)]
        cell.set_text_props(fontweight="bold", color="white")
        cell.set_facecolor("#0d4ec9")
    fig.suptitle(f"Run parameters — {run_dir.name}", fontsize=12, y=0.99)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return out_path
