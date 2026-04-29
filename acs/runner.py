"""Stage 1a runner — drives a single staircase step end-to-end.

Reads a YAML config (stage1a_*.yaml), constructs the solver in dimensionless
units, runs the time loop, writes HDF5 frames + a metrics CSV, and emits a
PASS/FAIL gate report at `results/{run.name}/gate_report.md`.

Designed to be importable (so tests can drive it) and runnable as a CLI:

    python -m acs.runner configs/stage1a_pilot.yaml
"""

from __future__ import annotations

import csv
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from acs.analysis.shape_metrics import (
    shape_metrics,
    shell_density_profile,
    top_down_projection_area,
)
from acs.config import load_config
from acs.gpu import init_taichi
from acs.gpu_profiler import GpuProfiler
from acs.io.hdf5_writer import FrameWriter
from acs.logging_setup import setup_logging
from acs.physics.mlsmpm import MLSMPMSolver, SolverConfig
from acs.provenance import build_manifest, write_manifest

logger = logging.getLogger(__name__)


@dataclass
class GateResult:
    name: str
    passed: bool
    detail: str

    def render(self) -> str:
        marker = "PASS" if self.passed else "FAIL"
        return f"- **[{marker}]** {self.name}: {self.detail}"


def _solver_cfg_from_yaml(cfg: dict) -> SolverConfig:
    nd = cfg["nondim"]
    nm = cfg["numerics"]
    sim = cfg["simulation"]
    sub = cfg.get("substrate", {})
    layer2 = cfg.get("layer2", {})
    layer2_b = cfg.get("layer2_b", {})
    layer7 = cfg.get("layer7", {})
    layer3 = cfg.get("layer3", {})
    gravity = cfg.get("gravity", {})
    layer4 = cfg.get("layer4", {})
    layer5 = cfg.get("layer5", {})
    layer6 = cfg.get("layer6", {})
    # Stage 2 Layer 6 chemistry/ECM remodeling — minimal scope per
    # docs/stage2_sanity.md. PI full authorisation 2026-04-29.
    # PARTIAL Magic-Number Block analogous to ζ_star Option α' / α_osm /
    # gravity_star precedent (Egeblad-Werb 2002 IF 70 + Lu 2011 IF 113
    # framework anchored, dimensionless rates derived from cited
    # timescales rescaled to overdamped τ_relax = 60 s calibration).
    layer6_enabled = bool(layer6.get("enabled", False))
    alpha_mmp_star = float(layer6.get("alpha_mmp_star", 0.0))
    beta_deg_star = float(layer6.get("beta_deg_star", 0.0))
    ecm_strength_initial = float(layer6.get("ecm_strength_initial", 1.0))
    ecm_strength_min = float(layer6.get("ecm_strength_min", 0.1))
    # Stage 1d Layer 4 cellular Marangoni (per
    # docs/07_internal_flow_dynamics.md framework citing Pajic-Lijakovic
    # & Milivojevic Eur Biophys J 2022 + Fütterer Phys Rev Fluids 2022).
    # γ_max/γ_min anchored to Maître Science 2012 IF 47 cell-cell
    # adhesion energy range. PI full authorisation 2026-04-29 covers
    # PARTIAL Magic-Number Block per ζ_star Option α' precedent.
    layer4_enabled = bool(layer4.get("enabled", False))
    gamma_max_star = float(layer4.get("gamma_max_star", 0.0))
    gamma_min_star = float(layer4.get("gamma_min_star", 0.0))
    # Stage 1d.b Marangoni Mechanism A + F (per
    # docs/stage1d_b_marangoni_sanity.md). Default OFF (alpha_A = 0,
    # alpha_F = 0) → recovers legacy γ(φ) static map even with
    # dynamic_gamma=True.
    layer4_dynamic_gamma = bool(layer4.get("dynamic_gamma", False))
    tau_gamma_star = float(layer4.get("tau_gamma_star", 1.0))
    alpha_A_star = float(layer4.get("alpha_A_star", 0.0))
    alpha_F_star = float(layer4.get("alpha_F_star", 0.0))
    # Stage 1a++.b stochastic boundary events (per
    # docs/stage1a_pp_b_stochastic_sanity.md).
    layer2_b_enabled = bool(layer2_b.get("enabled", False))
    lambda_lam_star = float(layer2_b.get("lambda_lam_star", 0.0))
    impulse_lam_star = float(layer2_b.get("impulse_lam_star", 0.0))
    # Stage 1d.c ECM communication + protrusion state machine (per
    # docs/stage1d_c_ecm_communication_sanity.md).
    layer7_enabled = bool(layer7.get("enabled", False))
    eta_ecm_star = float(layer7.get("eta_ecm_star", 5.0))
    G_ecm_star = float(layer7.get("G_ecm_star", 1.0))
    lambda_ecm_star = float(layer7.get("lambda_ecm_star", 4.0))
    R_dep_star = float(layer7.get("R_dep_star", 1.0))
    T0_star = float(layer7.get("T0_star", 0.1))
    tau_lam_star = float(layer7.get("tau_lam_star", 10.0))
    tau_release_star = float(layer7.get("tau_release_star", 5.0))
    lambda_filo_star = float(layer7.get("lambda_filo_star", 0.05))
    r_form_star = float(layer7.get("r_form_star", 0.1))
    r_mature_star = float(layer7.get("r_mature_star", 0.05))
    alpha_edge_star = float(layer7.get("alpha_edge_star", 1.0))
    beta_ecm_star = float(layer7.get("beta_ecm_star", 0.5))
    beta_traction_star = float(layer7.get("beta_traction_star", 0.5))
    # Stage 1b.b (PI directive 2026-04-29 per docs/layer3_phi_audit.md):
    # `layer3_spatial_S` is DEPRECATED in favour of the φ_memory + c_act
    # split (intrinsically encodes contact-band gating in c_act ODE
    # while preserving formation memory). The flag is read for backwards
    # compat (silent ignore) but `layer3_split` (default True) is the
    # new control.
    layer3_spatial_S_legacy = bool(layer4.get("layer3_spatial_S", False))
    layer3_split = bool(layer3.get("split", True))
    layer3_kappa_act = float(layer3.get("kappa_act", 1.0))
    layer3_memory_eps_star = float(layer3.get("memory_eps_star", 0.0))
    if layer3_spatial_S_legacy and layer3_split:
        # Deprecated flag set with new path active — ignore the flag and
        # log a one-time warning (the c_act ODE already gates spatially).
        import logging
        logging.getLogger("acs.runner").warning(
            "layer4.layer3_spatial_S is DEPRECATED under layer3_split=True "
            "(Stage 1b.b φ_memory + c_act split). The flag is ignored; "
            "the c_act ODE intrinsically encodes contact-band spatial gating."
        )
    layer3_spatial_S = layer3_spatial_S_legacy
    # Stage 1c Layer 5 mechano-osmotic Tier 2 (per
    # docs/08_mechano_osmotic.md framework citing Guo PNAS 2017 IF 12 +
    # Venkova eLife 2022). PI full authorisation 2026-04-29 covers PARTIAL
    # Magic-Number Block on the ODE constants per docs/stage1c_sanity.md.
    layer5_enabled = bool(layer5.get("enabled", False))
    rho_osm_initial = float(layer5.get("rho_osm_initial", 1.0))
    alpha_osm_star = float(layer5.get("alpha_osm_star", 0.0))
    beta_osm_star = float(layer5.get("beta_osm_star", 0.0))
    rho_osm_min = float(layer5.get("rho_osm_min", 0.5))
    rho_osm_max = float(layer5.get("rho_osm_max", 1.6))
    # Path C effective gravity — body-force coefficient framing
    # (PARTIAL Magic-Number Block per ζ_star Option α' precedent;
    # framework anchored to Stewart Nature 2011 IF 65 cell density,
    # value framed as a body-force coefficient under the overdamped
    # solver's ξ_star = 1 calibration). See docs/path_c_sanity.md.
    gravity_star = float(gravity.get("gravity_star", 0.0))
    # Stage 1a++ Layer 2 active stress (Option α' resolution 2026-04-29):
    # ζ/K dimensionless ratio, no Pa claim. K is anchored to Fischer-
    # Friedrich Nat Cell Biol 2014 IF 30 (already cited in
    # docs/02_force_models.md §1.1). Marchetti Rev Mod Phys 2013 IF 50
    # retained as framework reference.
    zeta_star = float(layer2.get("zeta_star", 0.0))
    # Stage 1b Layer 3 φ-ODE (PI full authorisation 2026-04-29):
    # Cho et al. 2020 mechanism, Halbleib & Nelson 2006, Hynes 2002 framework.
    # See docs/stage1b_layer3_sanity.md and docs/03_adhesion_dynamics.md.
    layer3_enabled = bool(layer3.get("enabled", False))
    phi_initial = float(layer3.get("phi_initial", 0.05))
    k_plus_star = float(layer3.get("k_plus_star", 0.0))
    k_minus_star = float(layer3.get("k_minus_star", 0.0))
    zeta_min = float(layer3.get("zeta_min", zeta_star))
    zeta_max = float(layer3.get("zeta_max", zeta_star))
    # Stage 1a+ Option β: substrate adhesion energy is anchored to Ca_cc
    # via the sweep multiplier α (Maître Science 2012, IF 47, anchors
    # γ_cc; α is parameter-free at result level — see
    # docs/stage1a_plus_substrate_sanity.md §"Option β addendum"
    # Magic-Number Block). The dimensionless substrate adhesion is then
    #   γ_sub_star = α · γ_cc_star = α · Ca_cc · K_star · radius_star.
    # Default α = 0 recovers Option α (mechanical anchor only, γ_sub = 0).
    Ca_cc = float(nd["capillary_number"])
    K_star = float(nd["K_star"])
    R_star = float(nd["radius_star"])
    alpha = float(sub.get("gamma_sub_alpha", 0.0))
    gamma_sub_star = alpha * Ca_cc * K_star * R_star
    return SolverConfig(
        n_particles=int(sim["n_material_points"]),
        grid_n=int(nm["background_grid_resolution"]),
        domain_star=float(nd["domain_star"]),
        radius_star=R_star,
        dt_star=float(nd["dt_star"]),
        K_star=K_star,
        mu_star=float(nd["mu_star"]),
        tau_star=float(nd["tau_star"]),
        capillary_number=Ca_cc,
        drag_xi_star=float(nd["drag_xi_star"]),
        density_star=float(nd["density_star"]),
        free_surface_threshold=float(nm["free_surface_density_threshold"]),
        seed=int(cfg["run"]["seed"]),
        substrate_enabled=bool(sub.get("enabled", False)),
        n_contact_band=int(sub.get("n_contact_band", 3)),
        gamma_sub_star=gamma_sub_star,
        zeta_star=zeta_star,
        layer3_enabled=layer3_enabled,
        phi_initial=phi_initial,
        k_plus_star=k_plus_star,
        k_minus_star=k_minus_star,
        zeta_min=zeta_min,
        zeta_max=zeta_max,
        gravity_star=gravity_star,
        layer5_enabled=layer5_enabled,
        rho_osm_initial=rho_osm_initial,
        alpha_osm_star=alpha_osm_star,
        beta_osm_star=beta_osm_star,
        rho_osm_min=rho_osm_min,
        rho_osm_max=rho_osm_max,
        layer4_enabled=layer4_enabled,
        gamma_max_star=gamma_max_star,
        gamma_min_star=gamma_min_star,
        layer3_spatial_S=layer3_spatial_S,
        layer3_split=layer3_split,
        layer3_kappa_act=layer3_kappa_act,
        layer3_memory_eps_star=layer3_memory_eps_star,
        layer4_dynamic_gamma=layer4_dynamic_gamma,
        tau_gamma_star=tau_gamma_star,
        alpha_A_star=alpha_A_star,
        alpha_F_star=alpha_F_star,
        layer2_b_enabled=layer2_b_enabled,
        lambda_lam_star=lambda_lam_star,
        impulse_lam_star=impulse_lam_star,
        layer7_enabled=layer7_enabled,
        eta_ecm_star=eta_ecm_star,
        G_ecm_star=G_ecm_star,
        lambda_ecm_star=lambda_ecm_star,
        R_dep_star=R_dep_star,
        T0_star=T0_star,
        tau_lam_star=tau_lam_star,
        tau_release_star=tau_release_star,
        lambda_filo_star=lambda_filo_star,
        r_form_star=r_form_star,
        r_mature_star=r_mature_star,
        alpha_edge_star=alpha_edge_star,
        beta_ecm_star=beta_ecm_star,
        beta_traction_star=beta_traction_star,
        layer6_enabled=layer6_enabled,
        alpha_mmp_star=alpha_mmp_star,
        beta_deg_star=beta_deg_star,
        ecm_strength_initial=ecm_strength_initial,
        ecm_strength_min=ecm_strength_min,
    )


def _build_stage1dc_kwargs(solver, solver_cfg) -> dict:
    """Compute the optional Stage 1d.c HDF5 fields for write_frame.

    Returns kwargs for fields that are populated only when the relevant
    solver layers are active. Stress-channel disambiguation
    (sigma_vol/sigma_active/sigma_total/pressure/dev_norm) is also
    computed here per PI directive 2026-04-30 (HDF5 schema fix).
    """
    out: dict = {}
    # Stress-channel disambiguation. tau_dev is already computed and
    # passed; here we derive sigma_vol, sigma_active, sigma_total,
    # pressure, dev_norm host-side from the solver's per-particle
    # state.
    K = float(solver_cfg.K_star)
    rho_ref = float(solver._rho_ref_kernel[None])
    rho_p = solver._rho_kernel_p.to_numpy().astype(np.float64)
    tau_dev_np = solver.tau_dev.to_numpy().astype(np.float64)
    n = rho_p.shape[0]
    # K_eff per particle (Stage 1c K(ρ_osm) coupling).
    if solver_cfg.layer5_enabled:
        rho_osm_np = solver.rho_osm_p.to_numpy().astype(np.float64)
        K_eff = K * rho_osm_np
    else:
        K_eff = np.full(n, K)
    # σ_vol = K_eff (ρ_ref/ρ − 1) I (v15 (k.3) form)
    vol_strain = rho_ref / np.clip(rho_p, 1e-30, None) - 1.0
    sigma_vol_diag = K_eff * vol_strain  # scalar diagonal value
    sigma_vol = np.zeros((n, 3, 3), dtype=np.float32)
    for d in range(3):
        sigma_vol[:, d, d] = sigma_vol_diag.astype(np.float32)
    # σ_active = -ζ_eff·K I, only on boundary particles.
    is_b = solver.is_boundary.to_numpy().astype(bool)
    sigma_active = np.zeros((n, 3, 3), dtype=np.float32)
    if solver_cfg.layer3_enabled:
        phi_arr = solver.phi_p.to_numpy().astype(np.float64)
        zeta_eff = solver_cfg.zeta_min * (1.0 - phi_arr) + solver_cfg.zeta_max * phi_arr
    else:
        zeta_eff = np.full(n, float(solver_cfg.zeta_star))
    sigma_active_diag = -(zeta_eff * K)
    for d in range(3):
        sigma_active[is_b, d, d] = sigma_active_diag[is_b].astype(np.float32)
    sigma_total = (tau_dev_np + sigma_vol + sigma_active).astype(np.float32)
    pressure = (-(np.trace(sigma_total, axis1=1, axis2=2) / 3.0)).astype(np.float32)
    dev_norm = np.sqrt((tau_dev_np ** 2).sum(axis=(1, 2))).astype(np.float32)
    out["sigma_vol"] = sigma_vol
    out["sigma_active"] = sigma_active
    out["sigma_total"] = sigma_total
    out["pressure"] = pressure
    out["dev_norm"] = dev_norm
    # Stage 1d.c per-particle ECM/protrusion fields.
    if solver_cfg.layer7_enabled:
        out["traction_ecm"] = solver.traction_p.to_numpy().astype(np.float32)
        out["fa_strength"] = solver.fa_strength_p.to_numpy().astype(np.float32)
        out["protrusion_state"] = solver.protrusion_state_p.to_numpy().astype(np.int32)
        out["ecm_signal"] = solver.ecm_signal_p.to_numpy().astype(np.float32)
        out["polarity"] = solver.polarity_p.to_numpy().astype(np.float32)
    return out


def run_stage1a(config_path: Path | str) -> Path:
    """Run a Stage 1a staircase step. Returns path to the gate report."""
    cfg_path = Path(config_path)
    cfg = load_config(cfg_path)

    out_dir = Path(cfg["run"]["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(log_dir=out_dir, run_name=cfg["run"]["name"])
    logger.info("Starting %s ← %s", cfg["run"]["name"], cfg_path)

    # Provenance: config snapshot + git hash + host.
    manifest = build_manifest(cfg, repo_root=".", extra={"config_path": str(cfg_path)})
    write_manifest(manifest, out_dir / "run_manifest.json")
    (out_dir / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    # GPU init (env var ACS_GPU_BACKEND wins if set).
    backend = cfg["gpu"].get("backend", "auto")
    arch = init_taichi(
        backend=backend,
        device_memory_GB=float(cfg["gpu"].get("device_memory_GB", 4)),
    )
    logger.info("Taichi arch=%s", arch)

    # Build solver.
    solver_cfg = _solver_cfg_from_yaml(cfg)
    solver = MLSMPMSolver(solver_cfg)
    if solver_cfg.substrate_enabled:
        # Stage 1a+ Option α: place spheroid in contact with the rigid
        # substrate at z = 0. Spheroid centre at z* = R₀ → bottommost
        # particle at z* = 0. No free-fall transient (no gravity at this
        # stage; see docs/stage1a_plus_substrate_sanity.md scope §).
        centre = np.array(
            [solver_cfg.domain_star * 0.5,
             solver_cfg.domain_star * 0.5,
             solver_cfg.radius_star],
            dtype=np.float32,
        )
        if solver_cfg.gamma_sub_star > 0.0:
            alpha = solver_cfg.gamma_sub_star / max(
                solver_cfg.capillary_number * solver_cfg.K_star * solver_cfg.radius_star,
                1e-30,
            )
            logger.info(
                "Substrate enabled (Stage 1a+ Option β): n_contact_band=%d, "
                "γ_sub*=%.4e (α=γ_sub/γ_cc=%.3f, γ_cc*=%.4e), spheroid centre "
                "at z* = R₀ = %.4f",
                solver_cfg.n_contact_band,
                solver_cfg.gamma_sub_star,
                alpha,
                solver_cfg.capillary_number * solver_cfg.K_star * solver_cfg.radius_star,
                solver_cfg.radius_star,
            )
        else:
            logger.info(
                "Substrate enabled (Stage 1a+ Option α): n_contact_band=%d, "
                "γ_sub*=0 (mechanical anchor only), spheroid centre at z* = R₀ = %.4f",
                solver_cfg.n_contact_band, solver_cfg.radius_star,
            )
    if solver_cfg.zeta_star > 0.0:
        # Stage 1a++ Layer 2 active stress (Option α' framing): ζ/K
        # dimensionless ratio, no Pa claim. K anchored to Fischer-Friedrich
        # Nat Cell Biol 2014 IF 30; Marchetti 2013 IF 50 framework
        # reference. Result reported as response curve in ζ/K.
        logger.info(
            "Stage 1a++ Layer 2 active stress: ζ/K = ζ_star = %.4f "
            "(σ_act = -ζ·K·I on boundary particles; contractile cortex). "
            "Energy-monotone gate SUSPENDED per Cousin-Rule contract change "
            "(active stress injects energy by construction).",
            solver_cfg.zeta_star,
        )
    if solver_cfg.layer3_enabled:
        phi_eq_pred = solver_cfg.k_plus_star / max(
            solver_cfg.k_plus_star + solver_cfg.k_minus_star, 1e-30,
        )
        logger.info(
            "Stage 1b Layer 3 φ-ODE active: phi_initial=%.3f, k_+_star=%.4e, "
            "k_-_star=%.4e, φ_eq=%.3f. ζ(φ) coupling: ζ ∈ [%.3f, %.3f] "
            "(Cho 2020 framework; PI full authorization 2026-04-29).",
            solver_cfg.phi_initial, solver_cfg.k_plus_star,
            solver_cfg.k_minus_star, phi_eq_pred,
            solver_cfg.zeta_min, solver_cfg.zeta_max,
        )
    if solver_cfg.gravity_star > 0.0:
        # Path C: effective gravity body force.
        logger.info(
            "Path C effective gravity active: gravity_star = %.4e "
            "(framework anchored to Stewart Nature 2011 IF 65 cell density "
            "≈ 1.05 g/cm³, ρ_medium ≈ 1.00 g/cm³; specific value framed as "
            "body-force coefficient under overdamped ξ_star = 1 calibration "
            "per docs/path_c_sanity.md Magic-Number Block PARTIAL "
            "resolution; PI full authorisation 2026-04-29).",
            solver_cfg.gravity_star,
        )
    if solver_cfg.layer6_enabled:
        logger.info(
            "Stage 2 Layer 6 chemistry/ECM remodeling active: "
            "α_MMP_star=%.4e, β_deg_star=%.4e, ecm_strength_initial=%.3f, "
            "ecm_strength_min=%.3f (Egeblad-Werb 2002 Nat Rev Cancer "
            "IF 70 + Lu 2011 Nat Rev Mol Cell Biol IF 113 framework "
            "anchors per docs/stage2_sanity.md PARTIAL Magic-Number "
            "Block; PI full authorisation 2026-04-29).",
            solver_cfg.alpha_mmp_star, solver_cfg.beta_deg_star,
            solver_cfg.ecm_strength_initial, solver_cfg.ecm_strength_min,
        )
    if solver_cfg.layer4_enabled:
        logger.info(
            "Stage 1d Layer 4 cellular Marangoni active: "
            "γ_max_star=%.4f, γ_min_star=%.4f, γ_1=γ_min−γ_max=%.4f "
            "(Maître Science 2012 IF 47 anchor for γ range; Pajic-Lijakovic "
            "& Milivojevic Eur Biophys J 2022 framework reference per "
            "docs/07_internal_flow_dynamics.md; PI full authorisation "
            "2026-04-29). Layer 3 spatial S_p extension: %s.",
            solver_cfg.gamma_max_star, solver_cfg.gamma_min_star,
            solver_cfg.gamma_min_star - solver_cfg.gamma_max_star,
            solver_cfg.layer3_spatial_S,
        )
    if solver_cfg.layer5_enabled:
        logger.info(
            "Stage 1c Layer 5 mechano-osmotic Tier 2 active: "
            "rho_osm_initial=%.3f, α_osm_star=%.4e, β_osm_star=%.4e, "
            "ρ_osm ∈ [%.2f, %.2f]. K(ρ_osm) coupling: K_eff = K · ρ_osm "
            "(Guo PNAS 2017 IF 12 + Venkova eLife 2022 framework per "
            "docs/08_mechano_osmotic.md; PI full authorisation 2026-04-29).",
            solver_cfg.rho_osm_initial, solver_cfg.alpha_osm_star,
            solver_cfg.beta_osm_star,
            solver_cfg.rho_osm_min, solver_cfg.rho_osm_max,
        )
    else:
        centre = np.full(3, solver_cfg.domain_star * 0.5, dtype=np.float32)
    solver.initialize_sphere(centre)

    # Reference-state calibration (Hu et al. 2018 §4.3, Adami-Hu-Adams 2010
    # population partition).
    calib = solver.calibrate_reference_state()
    logger.info(
        "Reference calibration: ρ_ref(harmonic)=%.4f vs arith=%.4f, "
        "ρ_actual∈[%.4f, %.4f], F_scale∈[%.4f, %.4f]",
        calib["rho_ref_harmonic"], calib["rho_arith_mean"],
        calib["rho_actual_min"], calib["rho_actual_max"],
        calib["F_scale_min"], calib["F_scale_max"],
    )
    logger.info(
        "Calibration <J> diagnostic: well_resolved=%.8f (n=%d), "
        "boundary_subset=%.8f (n=%d), all=%.8f",
        calib["J_mean_well_resolved"], calib["n_well_resolved"],
        calib["J_mean_boundary_subset"], calib["n_boundary_subset"],
        calib["J_mean_all"],
    )
    logger.info(
        "Calibration consistency (task-7): n_W_well_resolved=%d vs "
        "n_rho_well_resolved=%d; ρ↔W Jaccard = %.4f",
        calib["n_W_well_resolved"], calib["n_well_resolved"],
        calib["rho_W_jaccard"],
    )
    (out_dir / "reference_calibration.json").write_text(
        json.dumps(calib, indent=2), encoding="utf-8"
    )

    # Static-sphere curvature validation (Brackbill 1992 implementation check).
    # For a sphere of radius R*, analytical κ = 2/R*. We assert the measured
    # surface κ matches within 10% before the time loop starts. This is a
    # SCHEME-CORRECTNESS gate — it does not depend on Ca, K, μ, or τ.
    curv = solver.measure_surface_curvature()
    if curv.get("valid"):
        rel_err = abs(curv["kappa_measured_mean"] / curv["kappa_analytical"] - 1.0)
        logger.info(
            "Curvature validation: κ_measured=%.4f vs analytical 2/R=%.4f "
            "(rel err %.2f%%, n=%d surface cells)",
            curv["kappa_measured_mean"], curv["kappa_analytical"],
            rel_err * 100.0, curv["n_surface_cells"],
        )
        curv["relative_error"] = rel_err
        logger.info(
            "CSF localisation (Adami-Hu-Adams 2010): "
            "bulk |∇c| mean = %.4f, surface |∇c| peak = %.4f, "
            "ratio = %.4f (gate < 0.10)",
            curv["bulk_grad_c_mean"], curv["surface_grad_c_peak"],
            curv["bulk_to_surface_grad_ratio"],
        )
    else:
        logger.warning("Curvature validation invalid: %s", curv.get("reason"))
        curv["relative_error"] = float("nan")
    (out_dir / "curvature_validation.json").write_text(
        json.dumps(curv, indent=2), encoding="utf-8"
    )

    # Initial-state diagnostics & R₀ baseline.
    inv0 = solver.invariants()
    x0 = solver.x.to_numpy()
    metrics0 = shape_metrics(x0)
    R0 = metrics0["effective_radius"]
    logger.info(
        "Initial: m=%.6e mom=%s KE=%.3e U=%.3e R0=%.4f ψ=%.3f",
        inv0["mass_star"], inv0["momentum_star"], inv0["kinetic_energy_star"],
        inv0["strain_energy_star"], R0, metrics0["wadell_sphericity"],
    )

    # Time loop with periodic frame writes + diagnostic snapshots.
    total_t = float(cfg["nondim"]["total_time_star"])
    frame_dt = float(cfg["nondim"]["frame_interval_star"])
    diag_every = int(cfg["simulation"]["diagnostic_interval_steps"])
    dt = solver_cfg.dt_star
    n_steps = int(np.ceil(total_t / dt))

    frames_per_save = max(1, int(round(frame_dt / dt)))
    metrics_rows: list[dict[str, Any]] = []
    # shell_rows is the long-format witness for the v15 (k.3) bulk-transmission
    # mechanism (Sanity-Gate check 6, measurement-protocol consistency). One row
    # per (frame, radial bin); read by analysis to verify the equilibrium
    # ρ_kernel(r/R₀) profile is flat across the bulk shell rather than
    # surface-only. Specification: docs/stage1a_interior_pressure_sanity.md
    # §6 (e) and docs/outcomes_v15.md.
    shell_rows: list[dict[str, Any]] = []
    SHELL_N_BINS = 10
    SHELL_R_MAX_FRAC = 1.2
    # Stage 1a+ Option α (γ_sub = 0): per-frame substrate diagnostics
    # (anchor force balance, contact-band ρ_kernel, contact area, apparent
    # contact angle). Written to contact_metrics.csv. Specification:
    # docs/stage1a_plus_substrate_sanity.md and docs/outcomes_stage1a_plus.md.
    contact_rows: list[dict[str, Any]] = []
    step_times: list[float] = []
    halted = False
    halt_reason = ""

    snap_path = out_dir / "snapshots.h5"
    gpu_log_path = out_dir / "gpu_log.csv"
    with GpuProfiler(gpu_log_path, interval_s=60.0) as gpu_prof, FrameWriter(
        snap_path,
        n_particles=solver_cfg.n_particles,
        domain_star=solver_cfg.domain_star,
        total_time_star=total_t,
        frame_interval_star=frame_dt,
    ) as writer:
        # Frame 0 (initial state). Per PI viz directive 2026-04-29 +
        # Stage 1d.c HDF5 schema rationalisation 2026-04-30: save
        # per-particle state fields for state_overlays.py + stress
        # channels properly disambiguated.
        _stage1dc_fields = _build_stage1dc_kwargs(solver, solver_cfg)
        writer.write_frame(
            0.0,
            position=solver.x.to_numpy(),
            velocity=solver.v.to_numpy(),
            F=solver.F.to_numpy(),
            tau_dev=solver.tau_dev.to_numpy(),
            is_boundary=solver.is_boundary.to_numpy(),
            phi=(solver.phi_p.to_numpy() if solver_cfg.layer3_enabled else None),
            phi_memory=(solver.phi_memory_p.to_numpy()
                        if (solver_cfg.layer3_enabled and solver_cfg.layer3_split) else None),
            c_act=(solver.c_act_p.to_numpy()
                   if (solver_cfg.layer3_enabled and solver_cfg.layer3_split) else None),
            rho_osm=(solver.rho_osm_p.to_numpy() if solver_cfg.layer5_enabled else None),
            gamma_p=(solver.gamma_p_state.to_numpy()
                     if (solver_cfg.layer4_enabled and solver_cfg.layer4_dynamic_gamma) else None),
            **_stage1dc_fields,
        )
        # v15 shell-density witness — frame 0.
        x0_np = solver.x.to_numpy()
        rho0_np = solver._rho_kernel_p.to_numpy()
        shell0 = shell_density_profile(
            x0_np, rho0_np,
            centre=np.asarray(metrics0["centroid"], dtype=np.float64),
            R0=R0, n_bins=SHELL_N_BINS, r_max_frac=SHELL_R_MAX_FRAC,
        )
        for b in range(SHELL_N_BINS):
            shell_rows.append({
                "frame_index": 0,
                "time_star": 0.0,
                "bin_index": b,
                "r_lo_over_R0": float(shell0["bin_lo_over_R0"][b]),
                "r_hi_over_R0": float(shell0["bin_hi_over_R0"][b]),
                "count": int(shell0["count_per_bin"][b]),
                "mean_rho_kernel": (
                    float(shell0["mean_rho_per_bin"][b])
                    if not np.isnan(shell0["mean_rho_per_bin"][b])
                    else float("nan")
                ),
            })
        # Phase 1.2 top-down projection area (PI Outstanding-Issue
        # Resolution 2026-04-29): xy-plane convex hull of ALL particles,
        # z-independent — simulation analog of PI's experimental top-down
        # microscope projection (the `Area_um2` column in
        # `data/experimental/260313_*.csv`). This is the meaningful
        # comparison metric vs PI experimental A/A₀, distinct from the
        # substrate-contact area which depopulates on lift-off.
        A_topdown_init = top_down_projection_area(x0_np)
        metrics_row0 = {
            "frame_index": 0,
            "time_star": 0.0,
            **{k: v for k, v in inv0.items() if not isinstance(v, np.ndarray)},
            "momentum_x": float(inv0["momentum_star"][0]),
            "momentum_y": float(inv0["momentum_star"][1]),
            "momentum_z": float(inv0["momentum_star"][2]),
            "wadell_sphericity": metrics0["wadell_sphericity"],
            "effective_radius": metrics0["effective_radius"],
            "radius_of_gyration": metrics0["radius_of_gyration"],
            "shell_bulk_mean_rho": shell0["bulk_mean_rho"],
            "shell_bulk_std_rho": shell0["bulk_std_rho"],
            "shell_bulk_n_particles": shell0["bulk_n_particles"],
            "A_topdown_star": A_topdown_init,
            "A_over_A0_topdown": 1.0,
            # active_power_star and n_boundary already in inv0 dict — see
            # MLSMPMSolver.invariants(). Keep them explicit for clarity.
        }
        if solver_cfg.substrate_enabled:
            # Frame 0: no step has run, so the impulse accumulator is at its
            # constructor-default 0.0 (and would be 0 anyway as no penetration
            # has been clamped). Substrate diagnostics still computable from
            # the initial-state particle positions.
            sub_diag0 = solver.substrate_diagnostics()
            if sub_diag0.get("valid"):
                metrics_row0.update({
                    "n_contact_band_particles": sub_diag0["n_contact_band_particles"],
                    "rho_kernel_contact_over_ref": sub_diag0["rho_kernel_contact_over_ref"],
                    "F_substrate_per_step": sub_diag0["F_substrate_per_step"],
                    "F_pressure_down": sub_diag0["F_pressure_down"],
                    "anchor_force_balance_rel_err": sub_diag0["anchor_force_balance_rel_err"],
                    "contact_area_xy_hull": sub_diag0["contact_area_xy_hull"],
                    "apparent_contact_angle_deg": sub_diag0["apparent_contact_angle_deg"],
                })
                contact_rows.append({
                    "frame_index": 0,
                    "time_star": 0.0,
                    **sub_diag0,
                })
        metrics_rows.append(metrics_row0)

        wall_start = time.perf_counter()
        last_snap_step = 0
        for step in range(1, n_steps + 1):
            t0 = time.perf_counter()
            solver.step()
            step_times.append(time.perf_counter() - t0)

            if step % diag_every == 0:
                inv = solver.invariants()
                if inv["nan_count"] > 0:
                    halted = True
                    halt_reason = f"NaN/Inf detected at step {step}"
                    logger.error(halt_reason)
                    break

            if step - last_snap_step >= frames_per_save:
                last_snap_step = step
                inv = solver.invariants()
                xs = solver.x.to_numpy()
                m = shape_metrics(xs)
                # v15 shell-density witness — runtime frame.
                rho_runtime = solver._rho_kernel_p.to_numpy()
                shell = shell_density_profile(
                    xs, rho_runtime,
                    centre=np.asarray(m["centroid"], dtype=np.float64),
                    R0=R0, n_bins=SHELL_N_BINS, r_max_frac=SHELL_R_MAX_FRAC,
                )
                _stage1dc_runtime = _build_stage1dc_kwargs(solver, solver_cfg)
                writer.write_frame(
                    step * dt,
                    position=xs,
                    velocity=solver.v.to_numpy(),
                    F=solver.F.to_numpy(),
                    tau_dev=solver.tau_dev.to_numpy(),
                    is_boundary=solver.is_boundary.to_numpy(),
                    phi=(solver.phi_p.to_numpy() if solver_cfg.layer3_enabled else None),
                    phi_memory=(solver.phi_memory_p.to_numpy()
                                if (solver_cfg.layer3_enabled and solver_cfg.layer3_split) else None),
                    c_act=(solver.c_act_p.to_numpy()
                           if (solver_cfg.layer3_enabled and solver_cfg.layer3_split) else None),
                    rho_osm=(solver.rho_osm_p.to_numpy() if solver_cfg.layer5_enabled else None),
                    gamma_p=(solver.gamma_p_state.to_numpy()
                             if (solver_cfg.layer4_enabled and solver_cfg.layer4_dynamic_gamma) else None),
                    **_stage1dc_runtime,
                )
                frame_idx = writer._frame_count - 1
                for b in range(SHELL_N_BINS):
                    shell_rows.append({
                        "frame_index": frame_idx,
                        "time_star": step * dt,
                        "bin_index": b,
                        "r_lo_over_R0": float(shell["bin_lo_over_R0"][b]),
                        "r_hi_over_R0": float(shell["bin_hi_over_R0"][b]),
                        "count": int(shell["count_per_bin"][b]),
                        "mean_rho_kernel": (
                            float(shell["mean_rho_per_bin"][b])
                            if not np.isnan(shell["mean_rho_per_bin"][b])
                            else float("nan")
                        ),
                    })
                # Phase 1.2 top-down projection area per frame.
                A_topdown_now = top_down_projection_area(xs)
                A_over_A0_topdown_now = (
                    A_topdown_now / A_topdown_init
                    if (A_topdown_init and np.isfinite(A_topdown_init) and A_topdown_init > 0)
                    else float("nan")
                )
                metrics_row = {
                    "frame_index": frame_idx,
                    "time_star": step * dt,
                    **{k: v for k, v in inv.items() if not isinstance(v, np.ndarray)},
                    "momentum_x": float(inv["momentum_star"][0]),
                    "momentum_y": float(inv["momentum_star"][1]),
                    "momentum_z": float(inv["momentum_star"][2]),
                    "wadell_sphericity": m["wadell_sphericity"],
                    "effective_radius": m["effective_radius"],
                    "radius_of_gyration": m["radius_of_gyration"],
                    "shell_bulk_mean_rho": shell["bulk_mean_rho"],
                    "shell_bulk_std_rho": shell["bulk_std_rho"],
                    "shell_bulk_n_particles": shell["bulk_n_particles"],
                    "A_topdown_star": A_topdown_now,
                    "A_over_A0_topdown": A_over_A0_topdown_now,
                }
                if solver_cfg.substrate_enabled:
                    sub_diag = solver.substrate_diagnostics()
                    if sub_diag.get("valid"):
                        metrics_row.update({
                            "n_contact_band_particles": sub_diag["n_contact_band_particles"],
                            "rho_kernel_contact_over_ref": sub_diag["rho_kernel_contact_over_ref"],
                            "F_substrate_per_step": sub_diag["F_substrate_per_step"],
                            "F_pressure_down": sub_diag["F_pressure_down"],
                            "anchor_force_balance_rel_err": sub_diag["anchor_force_balance_rel_err"],
                            "contact_area_xy_hull": sub_diag["contact_area_xy_hull"],
                            "apparent_contact_angle_deg": sub_diag["apparent_contact_angle_deg"],
                        })
                        contact_rows.append({
                            "frame_index": frame_idx,
                            "time_star": step * dt,
                            **sub_diag,
                        })
                metrics_rows.append(metrics_row)
                logger.info(
                    "step=%d t*=%.3f KE=%.3e U=%.3e R/R0=%.3f ψ=%.3f vmax=%.3e "
                    "ρ_bulk=%.4f±%.4f (n=%d)",
                    step, step * dt,
                    inv["kinetic_energy_star"],
                    inv["strain_energy_star"],
                    m["effective_radius"] / R0,
                    m["wadell_sphericity"],
                    inv["max_speed_star"],
                    shell["bulk_mean_rho"],
                    shell["bulk_std_rho"],
                    shell["bulk_n_particles"],
                )

        wall_total = time.perf_counter() - wall_start

    # Metrics CSV. Use the UNION of all rows' keys as fieldnames (some
    # rows may have substrate / Layer 4 / Layer 6 fields only when
    # solver.substrate_diagnostics() returns valid=True; the union covers
    # all populated fields). `restval=""` writes empty string for missing
    # keys (e.g. frame 0 has no diag if substrate band empty at t=0).
    metrics_path = out_dir / "metrics.csv"
    if metrics_rows:
        all_keys = []
        seen = set()
        for r in metrics_rows:
            for k in r.keys():
                if k not in seen:
                    seen.add(k)
                    all_keys.append(k)
        with metrics_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=all_keys, restval="")
            w.writeheader()
            w.writerows(metrics_rows)

    # v15 shell-density witness CSV (long format: one row per (frame, bin)).
    # Read by analysis to plot ρ_kernel(r/R₀) profile over time and verify
    # the v15 bulk-transmission mechanism. See
    # docs/stage1a_interior_pressure_sanity.md §6 (e).
    shell_path = out_dir / "shell_profile.csv"
    if shell_rows:
        keys = list(shell_rows[0].keys())
        with shell_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(shell_rows)

    # Stage 1a+ Option α substrate diagnostics CSV (one row per frame).
    # Specification: docs/outcomes_stage1a_plus.md §"Files of record".
    if contact_rows:
        contact_path = out_dir / "contact_metrics.csv"
        keys = list(contact_rows[0].keys())
        with contact_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(contact_rows)

    # ---------- Gate evaluation ----------
    g = cfg["gate"]
    results: list[GateResult] = []

    inv_final = solver.invariants()

    # Reference-calibration scheme-correctness gate. By construction (Adami-
    # Hu-Adams 2010 population-aware partition), <J>_well_resolved must equal
    # 1.0 to floating-point round-off. A violation indicates a population
    # mismatch in calibrate_reference_state — exactly the v10 failure mode.
    # Tolerance 1e-6 is f32-roundoff scale (np.cbrt + clamp); not a physics
    # tolerance and not a fittable parameter.
    J_GATE_TOL = 1.0e-6
    j_resolved = float(calib["J_mean_well_resolved"])
    j_resolved_err = abs(j_resolved - 1.0)
    results.append(GateResult(
        "calibration <J>_well_resolved == 1",
        j_resolved_err <= J_GATE_TOL,
        f"<J>_well_resolved = {j_resolved:.8f} "
        f"(|err| = {j_resolved_err:.2e}, limit {J_GATE_TOL:.0e}); "
        f"<J>_boundary_subset = {calib['J_mean_boundary_subset']:.8f}, "
        f"<J>_all = {calib['J_mean_all']:.8f}",
    ))

    # Curvature validation — a scheme-correctness gate independent of
    # integration tolerances.
    if curv.get("valid"):
        results.append(GateResult(
            "curvature operator (κ vs 2/R)",
            curv["relative_error"] <= 0.10,
            f"κ_measured={curv['kappa_measured_mean']:.4f}, "
            f"analytical=2/R={curv['kappa_analytical']:.4f}, "
            f"rel err {curv['relative_error'] * 100:.1f}% (limit 10%)",
        ))
        # CSF localisation gate: with Adami-Hu-Adams 2010 reproducing-kernel
        # normalisation, c ≡ 1 in any cell whose kernel sees particles, so
        # |∇c| should be confined to the surface band. The bulk-shell
        # mean / surface-peak ratio quantifies the residual interior
        # penetration. Limit 0.10 is set by analogy to the 10% curvature
        # tolerance — both are scheme-correctness contracts at floating-point
        # scale rather than physics tolerances. v11 baseline (no AHA): 0.61
        # at n_smoothing_passes=2 (FAIL), 0.50 at n_passes=4 (still FAIL).
        BULK_GRAD_LIMIT = 0.10
        ratio = curv["bulk_to_surface_grad_ratio"]
        results.append(GateResult(
            "CSF localisation (bulk/surface |∇c| ratio)",
            ratio <= BULK_GRAD_LIMIT,
            f"bulk |∇c| mean = {curv['bulk_grad_c_mean']:.4f}, "
            f"surface |∇c| peak = {curv['surface_grad_c_peak']:.4f}, "
            f"ratio = {ratio:.4f} (limit {BULK_GRAD_LIMIT:.2f})",
        ))
    else:
        results.append(GateResult(
            "curvature operator (κ vs 2/R)",
            False,
            f"validation invalid: {curv.get('reason')}",
        ))
        results.append(GateResult(
            "CSF localisation (bulk/surface |∇c| ratio)",
            False,
            f"curvature validation invalid: {curv.get('reason')}",
        ))

    mass_drift = abs(inv_final["mass_star"] - inv0["mass_star"]) / max(inv0["mass_star"], 1e-30)
    results.append(GateResult(
        "mass conservation",
        mass_drift <= float(g["mass_drift_rel_max"]),
        f"|Δm/m₀| = {mass_drift:.2e} (limit {g['mass_drift_rel_max']:.0e})",
    ))

    # Momentum drift gate — Option F Week 2 contract change (per
    # docs/horizontal_momentum_drift_investigation.md, PI-authorised
    # 2026-04-29). The previous denominator m · max(v_rms, V_FLOOR=1e-3)
    # collapses to the V_FLOOR scale in overdamped equilibrium, putting the
    # required |Δp| below the f32 grain noise floor for ~5e3 particles ×
    # 5e5 steps. Across 4 runs (Production Lam4 + 3 Phase 4-v2 pilots) the
    # absolute |Δp_xy| was uniformly bounded ~1e-3, dominated by Poisson-
    # disk pack asymmetry (~1.4% N⁻¹/²) + accumulated atomic-op f32 noise
    # — i.e. NOT a directional solver bug.
    #
    # New gate: absolute drift ≤ a regime-floor calibrated empirically
    # (default 2.0e-3, configurable via gate.momentum_drift_abs_max).
    # The legacy ratio is still reported for backwards-comparable
    # information.
    #
    # Stage 1a+ Option α: −z reflective BC is a substrate momentum-leak
    # channel by design, so vertical momentum is excluded under
    # substrate_enabled; gate measures only (x, y) horizontal drift.
    p_init = np.array(inv0["momentum_star"])
    p_final = np.array(inv_final["momentum_star"])
    v_rms = float(np.sqrt(2.0 * inv_final["kinetic_energy_star"] / max(inv_final["mass_star"], 1e-30)))
    V_FLOOR = 1.0e-3  # legacy denominator floor (kept for the info-only ratio)
    p_norm_scale = inv_final["mass_star"] * max(v_rms, V_FLOOR)
    p_abs_limit = float(g.get("momentum_drift_abs_max", 2.0e-3))
    if solver_cfg.substrate_enabled:
        delta_p_horiz = (p_final - p_init)[:2]
        delta_p_norm = float(np.linalg.norm(delta_p_horiz))
        p_drift_ratio = delta_p_norm / max(p_norm_scale, 1e-30)
        results.append(GateResult(
            "momentum drift (horizontal abs, substrate absorbs vertical)",
            delta_p_norm <= p_abs_limit,
            f"|Δp_xy| = {delta_p_norm:.2e} (limit {p_abs_limit:.1e}); "
            f"info: |Δp_xy|/(m·max(v_rms,{V_FLOOR:.0e})) = {p_drift_ratio:.2e}, "
            f"|Δp_z|={abs(float((p_final - p_init)[2])):.2e} (substrate-leak, not gated), "
            f"v_rms={v_rms:.2e}",
        ))
    else:
        delta_p = p_final - p_init
        delta_p_norm = float(np.linalg.norm(delta_p))
        p_drift_ratio = delta_p_norm / max(p_norm_scale, 1e-30)
        results.append(GateResult(
            "momentum drift (abs)",
            delta_p_norm <= p_abs_limit,
            f"|Δp| = {delta_p_norm:.2e} (limit {p_abs_limit:.1e}); "
            f"info: |Δp|/(m·max(v_rms,{V_FLOOR:.0e})) = {p_drift_ratio:.2e}, "
            f"v_rms={v_rms:.2e}",
        ))

    # Stage 1a++ contract change extended to L3/L4/L5/L6 (Cousin-Rule
    # extension finalised in Phase 1.1, 2026-04-29). When any active-matter
    # / chemistry channel is on (Layer 2 ζ, Layer 3 ζ(φ), Layer 4 Marangoni
    # surface work, Layer 5 osmotic K(ρ), Layer 6 chemistry/ECM remodeling),
    # energy is injected by construction; the energy-monotone gate is
    # suspended and replaced with a finiteness check on the cumulative
    # active power. The Stage 1d sanity-md flagged this extension; this
    # commit implements it. PI full authorisation 2026-04-29.
    energy_gate_suspended = (
        solver_cfg.zeta_star > 0.0
        or solver_cfg.layer3_enabled
        or solver_cfg.layer4_enabled
        or solver_cfg.layer5_enabled
        or solver_cfg.layer6_enabled
    )
    if energy_gate_suspended:
        # Suspended; report as informational.
        active_powers = np.array(
            [r.get("active_power_star", 0.0) for r in metrics_rows],
            dtype=float,
        )
        max_abs_power = float(np.abs(active_powers).max() if len(active_powers) else 0.0)
        # Strain-energy reference for the "active power finite vs strain
        # energy" sanity check.
        U_strain_max = float(max(
            (r.get("strain_energy_star", 0.0) for r in metrics_rows),
            default=1e-30,
        ))
        # Active-work-finite gate: per-frame active power bounded by 10×
        # strain-energy scale (rough sanity check that active power doesn't
        # blow up unboundedly).
        results.append(GateResult(
            "active power finite (Stage 1a++ Cousin-Rule replacement for energy-monotone)",
            np.isfinite(max_abs_power) and max_abs_power <= 10.0 * max(U_strain_max, 1e-30),
            f"max |P_act| = {max_abs_power:.3e}, max U_strain = {U_strain_max:.3e}, "
            f"ratio = {max_abs_power / max(U_strain_max, 1e-30):.2f} (limit 10.0)",
        ))
    else:
        # Path C extension: when gravity is active, include U_grav in the
        # energy-monotone sum (per docs/path_c_sanity.md check 3 contract
        # extension). When gravity_star == 0, the term is identically zero
        # and the sum is identical to the original Stage 1a baseline.
        energies = np.array([
            r["kinetic_energy_star"] + r["strain_energy_star"]
            + r.get("grav_pe_star", 0.0)
            for r in metrics_rows
        ], dtype=float)
        if len(energies) >= 2:
            E_max = float(energies.max())
            increases = np.diff(energies)
            max_inc = float(increases.max() if len(increases) else 0.0)
            E_tol = float(g["energy_increase_tol_rel"]) * max(E_max, 1e-30)
            results.append(GateResult(
                "energy monotone (KE + strain)",
                max_inc <= E_tol,
                f"max(ΔE) = {max_inc:.3e} ≤ tol {E_tol:.3e}",
            ))
        else:
            results.append(GateResult("energy monotone (KE + strain)", False, "insufficient frames"))

    nan_ok = inv_final["nan_count"] == 0 and not halted
    results.append(GateResult(
        "no NaN/Inf",
        nan_ok,
        f"halt={halted} reason='{halt_reason}' nan_count={inv_final['nan_count']}",
    ))

    R_check_t = float(g["radius_drift_check_after_s"]) / float(cfg["physics"]["maxwell_tau_s"])
    R_after = [r for r in metrics_rows if r["time_star"] >= R_check_t]
    if R_after:
        R_drift = max(abs(r["effective_radius"] / R0 - 1.0) for r in R_after)
        results.append(GateResult(
            "radius drift |R/R₀ − 1|",
            R_drift <= float(g["radius_drift_rel_max"]),
            f"max drift = {R_drift:.3f} (limit {g['radius_drift_rel_max']:.3f})",
        ))
    else:
        R_drift = float("nan")
        results.append(GateResult("radius drift |R/R₀ − 1|", False, "no frames after check_after time"))

    # Stage 1a+ Option α additional gates (substrate enabled).
    if solver_cfg.substrate_enabled:
        # (i) R drift improvement vs Stage 1a v15 baseline (24.4%). This is a
        # *relative-improvement* gate, not an absolute tolerance: the substrate
        # must do *some* mechanical work. See
        # docs/outcomes_stage1a_plus.md §"Bounded outcomes".
        v15_baseline = float(g.get("v15_R_drift_baseline", 0.244))
        if not np.isnan(R_drift):
            results.append(GateResult(
                "R drift improvement vs v15 baseline",
                R_drift < v15_baseline,
                f"R_drift_1aplus = {R_drift:.3f} vs v15 baseline {v15_baseline:.3f} "
                f"(strict-less requirement; bucketing → docs/outcomes_stage1a_plus.md)",
            ))

        # (ii) Anchor force balance: substrate reaction = compressive bulk
        # pressure + gravity (Option F Week 2 contract change, per
        # docs/anchor_force_balance_investigation.md).
        # F_substrate ≈ F_compressive + F_gravity_band; tolerance ≤ 0.20.
        rel_errs_after = [
            r.get("anchor_force_balance_rel_err", float("nan"))
            for r in metrics_rows if r["time_star"] >= R_check_t
        ]
        rel_errs_clean = [e for e in rel_errs_after if not (e is None or np.isnan(e))]
        if rel_errs_clean:
            balance_med = float(np.median(rel_errs_clean))
            results.append(GateResult(
                "anchor force balance |F_sub − (F_compr + F_grav)| / |F_sub|",
                balance_med <= 0.20,
                f"median over post-transient frames = {balance_med:.3f} (limit 0.20, "
                f"n={len(rel_errs_clean)} frames)",
            ))
        else:
            results.append(GateResult(
                "anchor force balance |F_sub − (F_compr + F_grav)| / |F_sub|",
                False,
                "no valid anchor-force-balance samples",
            ))

        # (iii) Contact-band ρ_kernel / ρ_ref ∈ [0.65, 1.15] (Option F
        # Week 2 contract change, per docs/anchor_force_balance_
        # investigation.md). The lower bound 0.65 reflects the AHA §3
        # kernel truncation envelope at the −z reflective boundary
        # (~30% under-densification expected); the upper bound 1.15 is
        # unchanged (over-densification would still indicate a packing
        # artefact).
        RHO_LO, RHO_HI = 0.65, 1.15
        rho_after = [
            r.get("rho_kernel_contact_over_ref", float("nan"))
            for r in metrics_rows if r["time_star"] >= R_check_t
        ]
        rho_clean = [e for e in rho_after if not (e is None or np.isnan(e))]
        if rho_clean:
            rho_med = float(np.median(rho_clean))
            results.append(GateResult(
                f"contact-band ρ_kernel / ρ_ref ∈ [{RHO_LO}, {RHO_HI}]",
                RHO_LO <= rho_med <= RHO_HI,
                f"median over post-transient frames = {rho_med:.3f} "
                f"(window [{RHO_LO}, {RHO_HI}], n={len(rho_clean)} frames)",
            ))
        else:
            results.append(GateResult(
                f"contact-band ρ_kernel / ρ_ref ∈ [{RHO_LO}, {RHO_HI}]",
                False,
                "no valid contact-band samples",
            ))

    # Stage 1a++ Layer 2 additional gates (zeta_star > 0).
    if solver_cfg.zeta_star > 0.0:
        # (iv) R drift improvement vs Stage 1a+ Option β α=1.0 baseline
        # (= 0.247, the best-anchored Layer-1+substrate result). Stage
        # 1a++ must improve on this to demonstrate Layer 2 contribution.
        # See docs/outcomes_stage1a_plus_plus.md §"Mechanism question".
        beta_alpha1_baseline = float(g.get("stage1a_plus_beta_alpha1_R_drift_baseline", 0.247))
        if not np.isnan(R_drift):
            results.append(GateResult(
                "R drift improvement vs Stage 1a+ Option β α=1.0 baseline",
                R_drift < beta_alpha1_baseline,
                f"R_drift_1aplusplus = {R_drift:.3f} vs β α=1.0 baseline "
                f"{beta_alpha1_baseline:.3f} (strict-less; bucketing → "
                f"docs/outcomes_stage1a_plus_plus.md)",
            ))

        # (v) Boundary-tag stability: |Δn_boundary / n_boundary| per frame
        # ≤ 0.10. Numerical hygiene — if the boundary tag flickers
        # pathologically, the active-stress measurement is meaningless.
        n_b_series = [int(r.get("n_boundary", 0)) for r in metrics_rows]
        if len(n_b_series) >= 2:
            denom = max(np.mean(n_b_series), 1.0)
            flickers = np.abs(np.diff(n_b_series)) / denom
            max_flicker = float(flickers.max())
            results.append(GateResult(
                "boundary-tag stability |Δn_boundary| / <n_boundary> per frame",
                max_flicker <= 0.10,
                f"max frame-to-frame flicker = {max_flicker:.3f} "
                f"(<n_boundary>={denom:.1f}, limit 0.10, n_frames={len(n_b_series)})",
            ))
        else:
            results.append(GateResult(
                "boundary-tag stability",
                False,
                "insufficient frames for flicker computation",
            ))

    # Stage 1b Layer 3 additional gates (layer3_enabled).
    if solver_cfg.layer3_enabled:
        # (i) φ ∈ [0, 1] per-particle invariant.
        phi_min_series = [r.get("phi_min", float("nan")) for r in metrics_rows]
        phi_max_series = [r.get("phi_max", float("nan")) for r in metrics_rows]
        phi_min_overall = float(np.nanmin(phi_min_series)) if phi_min_series else float("nan")
        phi_max_overall = float(np.nanmax(phi_max_series)) if phi_max_series else float("nan")
        results.append(GateResult(
            "φ ∈ [0, 1] per-particle invariant",
            phi_min_overall >= 0.0 and phi_max_overall <= 1.0,
            f"min(φ) over all frames = {phi_min_overall:.4f}, "
            f"max(φ) over all frames = {phi_max_overall:.4f}",
        ))

        n_p = float(solver_cfg.n_particles)
        if solver_cfg.layer3_split:
            # Stage 1b.b (PI directive 2026-04-29 per
            # docs/layer3_phi_audit.md §5): the legacy F9 "φ trajectory
            # toward predicted φ_eq" is DEPRECATED — it compared global
            # <φ> against the boundary-only φ_eq (category error) AND
            # the v11 spatial S=0 interior decay erased formation
            # phenotype memory on the Cho 2020 transition timescale
            # (semantic conflation). F9 is replaced by:
            #   5a φ_memory preservation: |<φ_memory>(end) − phi_init|
            #      ≤ memory_drift_max (default 0.01)
            #   5b boundary c_act trajectory: <c_act>_band(end) ∈
            #      [c_eq − tol, c_eq + tol] where c_eq = k_+/(k_++k_-)
            phi_memory_means = [
                float(r.get("phi_memory_sum", float("nan"))) / n_p
                for r in metrics_rows
                if not np.isnan(r.get("phi_memory_sum", float("nan")))
            ]
            if phi_memory_means:
                phi_mem_end = float(phi_memory_means[-1])
                phi_mem_drift = abs(phi_mem_end - solver_cfg.phi_initial)
                drift_max = float(g.get("layer3_memory_drift_max", 0.01))
                results.append(GateResult(
                    "5a φ_memory preservation |<φ_memory> − phi_init|",
                    phi_mem_drift <= drift_max,
                    f"<φ_memory>(end) = {phi_mem_end:.4f}, phi_init = "
                    f"{solver_cfg.phi_initial:.4f}, |drift| = {phi_mem_drift:.4f} "
                    f"(limit {drift_max:.4f})",
                ))
            else:
                results.append(GateResult(
                    "5a φ_memory preservation |<φ_memory> − phi_init|",
                    False,
                    "no φ_memory samples (Stage 1b.b split inactive?)",
                ))

            c_eq_pred = solver_cfg.k_plus_star / max(
                solver_cfg.k_plus_star + solver_cfg.k_minus_star, 1e-30,
            )
            c_act_band_means = []
            for r in metrics_rows:
                count = int(r.get("c_act_band_count", 0) or 0)
                if count > 0:
                    band_sum = float(r.get("c_act_band_sum", 0.0) or 0.0)
                    c_act_band_means.append(band_sum / count)
            if c_act_band_means:
                c_act_band_end = float(c_act_band_means[-1])
                c_eq_err = abs(c_act_band_end - c_eq_pred)
                c_eq_tol = float(g.get("layer3_c_act_band_tol", 0.10))
                results.append(GateResult(
                    "5b c_act boundary trajectory toward c_eq",
                    c_eq_err <= c_eq_tol,
                    f"<c_act>_band(end) = {c_act_band_end:.4f}, c_eq = "
                    f"{c_eq_pred:.4f}, |err| = {c_eq_err:.4f} (limit {c_eq_tol:.4f})",
                ))
            else:
                results.append(GateResult(
                    "5b c_act boundary trajectory toward c_eq",
                    False,
                    "no boundary-band c_act samples",
                ))
        else:
            # Legacy F9 retained when layer3_split=False (backwards compat).
            phi_means = [
                float(r.get("phi_sum", 0.0)) / n_p
                for r in metrics_rows
            ]
            phi_eq_pred = solver_cfg.k_plus_star / max(
                solver_cfg.k_plus_star + solver_cfg.k_minus_star, 1e-30,
            )
            if phi_means:
                phi_end = float(phi_means[-1])
                phi_eq_err = abs(phi_end - phi_eq_pred)
                results.append(GateResult(
                    "φ trajectory toward predicted φ_eq (legacy F9)",
                    phi_eq_err <= max(0.5, abs(solver_cfg.phi_initial - phi_eq_pred)),
                    f"<φ>(end) = {phi_end:.4f}, predicted φ_eq = {phi_eq_pred:.4f}, "
                    f"|err| = {phi_eq_err:.4f} (tolerance: ≤ max(0.5, "
                    f"|φ_init − φ_eq|) = {max(0.5, abs(solver_cfg.phi_initial - phi_eq_pred)):.4f})",
                ))
            else:
                results.append(GateResult(
                    "φ trajectory toward predicted φ_eq (legacy F9)",
                    False,
                    "no φ samples",
                ))

        # Stage 1c Layer 5 additional gates (layer5_enabled).
        if solver_cfg.layer5_enabled:
            # (Layer 5 i) ρ_osm ∈ [ρ_osm_min, ρ_osm_max] per-particle invariant.
            rho_osm_min_series = [r.get("rho_osm_min", float("nan")) for r in metrics_rows]
            rho_osm_max_series = [r.get("rho_osm_max", float("nan")) for r in metrics_rows]
            rho_osm_min_overall = float(np.nanmin(rho_osm_min_series)) if rho_osm_min_series else float("nan")
            rho_osm_max_overall = float(np.nanmax(rho_osm_max_series)) if rho_osm_max_series else float("nan")
            results.append(GateResult(
                f"ρ_osm ∈ [{solver_cfg.rho_osm_min}, {solver_cfg.rho_osm_max}] per-particle invariant",
                (
                    rho_osm_min_overall >= solver_cfg.rho_osm_min - 1e-6
                    and rho_osm_max_overall <= solver_cfg.rho_osm_max + 1e-6
                ),
                f"min(ρ_osm) over all frames = {rho_osm_min_overall:.4f}, "
                f"max(ρ_osm) over all frames = {rho_osm_max_overall:.4f}",
            ))

            # (Layer 5 ii) <ρ_osm> trajectory finite & non-pathological.
            n_p_l5 = float(solver_cfg.n_particles)
            rho_osm_means = [
                float(r.get("rho_osm_sum", 0.0)) / n_p_l5
                for r in metrics_rows
            ]
            if rho_osm_means:
                rho_osm_end = float(rho_osm_means[-1])
                # Allow slight relaxation below 1.0 but no catastrophic loss
                # (Guo 2017 mechanism: spreading raises ρ_osm above 1.0).
                results.append(GateResult(
                    "<ρ_osm> trajectory finite & non-pathological",
                    (
                        np.isfinite(rho_osm_end)
                        and 0.95 <= rho_osm_end <= solver_cfg.rho_osm_max + 1e-6
                    ),
                    f"<ρ_osm>(end) = {rho_osm_end:.4f} (window [0.95, "
                    f"{solver_cfg.rho_osm_max}], rho_osm_initial = "
                    f"{solver_cfg.rho_osm_initial:.3f})",
                ))
            else:
                results.append(GateResult(
                    "<ρ_osm> trajectory finite & non-pathological",
                    False,
                    "no ρ_osm samples",
                ))

        # Stage 2 Layer 6 additional gates (layer6_enabled).
        if solver_cfg.layer6_enabled:
            # (Layer 6 i) ecm_strength ∈ [ecm_strength_min, 1.0] invariant.
            ecm_series = [r.get("ecm_strength", float("nan")) for r in metrics_rows]
            ecm_series_clean = [e for e in ecm_series if not (e is None or np.isnan(e))]
            if ecm_series_clean:
                ecm_min_overall = float(min(ecm_series_clean))
                ecm_max_overall = float(max(ecm_series_clean))
                results.append(GateResult(
                    f"ecm_strength ∈ [{solver_cfg.ecm_strength_min}, 1.0] invariant",
                    (
                        ecm_min_overall >= solver_cfg.ecm_strength_min - 1e-6
                        and ecm_max_overall <= 1.0 + 1e-6
                    ),
                    f"min(ecm_strength) over all frames = {ecm_min_overall:.4f}, "
                    f"max(ecm_strength) over all frames = {ecm_max_overall:.4f}",
                ))
            else:
                results.append(GateResult(
                    f"ecm_strength ∈ [{solver_cfg.ecm_strength_min}, 1.0] invariant",
                    False,
                    "no ecm_strength samples",
                ))

            # (Layer 6 ii) mmp_total finite & non-decreasing.
            mmp_series = [r.get("mmp_total", float("nan")) for r in metrics_rows]
            mmp_series_clean = [m for m in mmp_series if not (m is None or np.isnan(m))]
            if len(mmp_series_clean) >= 2:
                mmp_end = float(mmp_series_clean[-1])
                mmp_diffs = np.diff(np.array(mmp_series_clean, dtype=float))
                # Allow tiny negative jitter from f32 round-off (1e-9 tolerance).
                non_decreasing = bool((mmp_diffs >= -1e-9).all())
                finite = np.isfinite(mmp_end) and abs(mmp_end) <= 1e6
                results.append(GateResult(
                    "mmp_total finite & non-decreasing",
                    finite and non_decreasing,
                    f"mmp_total(end) = {mmp_end:.4e} (finite={finite}, "
                    f"min(Δmmp) = {float(mmp_diffs.min()):.3e}, n_frames={len(mmp_series_clean)})",
                ))
            else:
                results.append(GateResult(
                    "mmp_total finite & non-decreasing",
                    False,
                    "insufficient mmp_total samples",
                ))

        # (iii) A/A₀_topdown trajectory finite & non-pathological.
        # Option F Week 2 contract change (per CLAUDE.md Hard Rule 11 +
        # docs/gate_fail_taxonomy.md F8): the previous gate measured
        # contact_area_xy_hull, which is the *substrate-contact patch*
        # (depopulates on lift-off → artefactual A/A₀ ≈ 0 even when the
        # spheroid is intact). The PI experimental measurement is a top-
        # down microscope projection; the matching simulation metric is
        # A_over_A0_topdown (xy-plane convex hull of ALL particles), already
        # computed per frame and reported in metrics.csv. The legacy
        # contact-hull series is kept for diagnostic purposes only and is
        # not gated.
        A_series = [
            r.get("A_over_A0_topdown", float("nan"))
            for r in metrics_rows
        ]
        A_clean = [a for a in A_series if not (a is None or np.isnan(a))]
        if A_clean:
            min_ratio = float(min(A_clean))
            max_ratio = float(max(A_clean))
            results.append(GateResult(
                "A/A₀_topdown trajectory finite & non-pathological",
                np.isfinite(min_ratio) and np.isfinite(max_ratio) and min_ratio >= 0.5,
                f"A/A₀_topdown ∈ [{min_ratio:.3f}, {max_ratio:.3f}] "
                f"(limit min ≥ 0.5; n_frames = {len(A_clean)})",
            ))
        else:
            results.append(GateResult(
                "A/A₀_topdown trajectory finite & non-pathological",
                False,
                "no valid A_over_A0_topdown samples",
            ))

    # Wadell sphericity gate — Option F Week 2 contract change. The Stage 1a
    # threshold 0.95 was correct for free-floating relaxation but
    # incompatible with substrate spreading: under Stage 1a+ and beyond,
    # the spheroid *must* deform (Codex review item 6 in
    # docs/codex_review_synthesis.md; F6 in docs/gate_fail_taxonomy.md).
    # Threshold is now stage-aware:
    #   - Stage 1a (no substrate)            : sphericity_min (default 0.95)
    #   - Stage 1a+ and beyond (substrate)   : sphericity_min_post_spread
    #                                          (default 0.70)
    # The post-spread threshold is set so that catastrophic deformation
    # (ψ < 0.7 indicates the spheroid has lost its spheroidal identity)
    # is still caught while normal spreading-induced flattening passes.
    psi_check_t = float(g["sphericity_check_after_s"]) / float(cfg["physics"]["maxwell_tau_s"])
    psi_after = [r["wadell_sphericity"] for r in metrics_rows if r["time_star"] >= psi_check_t]
    if solver_cfg.substrate_enabled:
        psi_threshold = float(g.get("sphericity_min_post_spread", 0.70))
        psi_label = "Wadell sphericity ψ (post-spread)"
    else:
        psi_threshold = float(g["sphericity_min"])
        psi_label = "Wadell sphericity ψ (free-floating)"
    if psi_after:
        psi_min_seen = float(np.nanmin(psi_after))
        results.append(GateResult(
            psi_label,
            psi_min_seen >= psi_threshold,
            f"min ψ after equilibration = {psi_min_seen:.3f} (limit {psi_threshold:.3f})",
        ))
    else:
        results.append(GateResult(psi_label, False, "no frames after check_after time"))

    speeds_max = float(inv_final["max_speed_star"])
    results.append(GateResult(
        "max speed bound",
        speeds_max <= float(g["max_speed_over_vrms"]) * max(v_rms, 1e-12) or v_rms == 0,
        f"max v* = {speeds_max:.3e}, v_rms = {v_rms:.3e}",
    ))

    gpu_stats = gpu_prof.summarise()
    peak_vram = gpu_stats.get("peak_vram_GB", 0.0) if gpu_stats.get("available") else 0.0
    if peak_vram > 0:
        results.append(GateResult(
            "peak VRAM ≤ ceiling",
            peak_vram <= float(g["vram_peak_gb_max"]),
            f"peak {peak_vram:.2f} GB (limit {g['vram_peak_gb_max']:.1f} GB)",
        ))

    mean_step_ms = (sum(step_times) / len(step_times) * 1000.0) if step_times else float("nan")
    wall_clock_min = wall_total / 60.0 if step_times else 0.0

    overall_pass = all(r.passed for r in results)

    # ---------- Gate report ----------
    curv_lines: list[str] = ["", "## Curvature validation"]
    if curv.get("valid"):
        curv_lines += [
            f"- measured κ at surface: mean {curv['kappa_measured_mean']:.4f}, "
            f"median {curv['kappa_measured_median']:.4f}, "
            f"std {curv['kappa_measured_std']:.4f}",
            f"- analytical 2/R = {curv['kappa_analytical']:.4f}",
            f"- relative error: {curv['relative_error'] * 100:.2f}% (limit 10%)",
            f"- surface band: {curv['n_surface_cells']} cells "
            f"(selection: {curv.get('selection', 'unknown')})",
            "",
            "### CSF localisation (Adami-Hu-Adams 2010 §3)",
            f"- bulk |∇c| mean (r/R₀ ∈ [0.2, 0.7], n={curv['n_bulk_cells_sampled']} cells) = "
            f"{curv['bulk_grad_c_mean']:.4f}",
            f"- surface |∇c| peak = {curv['surface_grad_c_peak']:.4f}",
            f"- ratio bulk/surface = {curv['bulk_to_surface_grad_ratio']:.4f} "
            f"(limit 0.10)",
        ]
    else:
        curv_lines.append(f"- validation invalid: {curv.get('reason')}")

    perf_lines: list[str] = ["", "## Performance Profiling"]
    if gpu_stats.get("available"):
        perf_lines += [
            f"- nvidia-smi samples: {gpu_stats['n_samples']} (1/min)",
            f"- avg GPU util: {gpu_stats['avg_gpu_util_pct']:.1f}%  "
            f"(p50 {gpu_stats['gpu_util_p50']:.1f}%, p95 {gpu_stats['gpu_util_p95']:.1f}%)",
            f"- peak VRAM: {gpu_stats['peak_vram_GB']:.2f} GB",
            f"- avg power: {gpu_stats['avg_power_W']:.1f} W",
        ]
    else:
        perf_lines.append("- nvidia-smi sampler unavailable on this host.")

    report = [
        f"# Gate report — {cfg['run']['name']}",
        "",
        f"- **Overall**: {'PASS' if overall_pass else 'FAIL'}",
        f"- Wall-clock: {wall_clock_min:.2f} min",
        f"- Steps: {len(step_times)} (mean {mean_step_ms:.3f} ms/step)",
        f"- Frames: {len(metrics_rows)}",
        f"- Peak VRAM: {peak_vram:.2f} GB" if peak_vram > 0 else "- Peak VRAM: n/a (nvidia-smi unavailable)",
        f"- Initial R₀\\* = {R0:.4f}",
        f"- Reference calibration: ρ_ref(harmonic)={calib['rho_ref_harmonic']:.4f}, "
        f"<J>_well_resolved={calib['J_mean_well_resolved']:.8f} "
        f"(n={calib['n_well_resolved']}), "
        f"<J>_boundary_subset={calib['J_mean_boundary_subset']:.8f} "
        f"(n={calib['n_boundary_subset']}), "
        f"<J>_all={calib['J_mean_all']:.8f}, "
        f"F_scale∈[{calib['F_scale_min']:.4f}, {calib['F_scale_max']:.4f}]",
        f"- ρ↔W well-resolved Jaccard (task-7 consistency): "
        f"{calib['rho_W_jaccard']:.4f}; "
        f"n_W_well_resolved={calib['n_W_well_resolved']} vs "
        f"n_ρ_well_resolved={calib['n_well_resolved']}",
        "",
        "## Checks",
        *(r.render() for r in results),
        *curv_lines,
        *perf_lines,
        "",
        "## Final invariants",
        "```json",
        json.dumps({k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in inv_final.items()}, indent=2),
        "```",
    ]
    report_path = out_dir / "gate_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    logger.info("Gate report → %s", report_path)
    logger.info("Overall: %s", "PASS" if overall_pass else "FAIL")
    return report_path


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python -m acs.runner <config.yaml>", file=sys.stderr)
        raise SystemExit(2)
    run_stage1a(sys.argv[1])
