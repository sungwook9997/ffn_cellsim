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

from acs.analysis.shape_metrics import shape_metrics, shell_density_profile
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
    layer3 = cfg.get("layer3", {})
    gravity = cfg.get("gravity", {})
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
    )


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
        # Frame 0 (initial state).
        writer.write_frame(
            0.0,
            position=solver.x.to_numpy(),
            velocity=solver.v.to_numpy(),
            F=solver.F.to_numpy(),
            tau_dev=solver.tau_dev.to_numpy(),
            is_boundary=solver.is_boundary.to_numpy(),
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
                writer.write_frame(
                    step * dt,
                    position=xs,
                    velocity=solver.v.to_numpy(),
                    F=solver.F.to_numpy(),
                    tau_dev=solver.tau_dev.to_numpy(),
                    is_boundary=solver.is_boundary.to_numpy(),
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

    # Metrics CSV.
    metrics_path = out_dir / "metrics.csv"
    if metrics_rows:
        keys = list(metrics_rows[0].keys())
        with metrics_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
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

    # Momentum drift: normalised by m · max(v_rms, v_floor) so the ratio does
    # not blow up as v_rms → 0 in a near-equilibrium run. The v_floor is
    # 0.001 in dimensionless units, i.e. one part per thousand of a ballistic
    # R₀/τ_relax; below that the system is at rest by any reasonable physical
    # standard and the *absolute* |Δp| is what matters.
    #
    # Stage 1a+ Option α contract change (per
    # docs/stage1a_plus_substrate_sanity.md §3): the −z reflective BC is a
    # deliberate momentum-leak channel for the substrate reaction force, so
    # vertical momentum is allowed to drift. Under substrate_enabled, the
    # gate measures only horizontal (x, y) momentum drift. This is a Cousin-
    # Rule contract change demanded by the new physical scenario, surfaced
    # to PI in the sanity-md and approved 2026-04-29.
    p_init = np.array(inv0["momentum_star"])
    p_final = np.array(inv_final["momentum_star"])
    v_rms = float(np.sqrt(2.0 * inv_final["kinetic_energy_star"] / max(inv_final["mass_star"], 1e-30)))
    V_FLOOR = 1.0e-3
    p_norm_scale = inv_final["mass_star"] * max(v_rms, V_FLOOR)
    if solver_cfg.substrate_enabled:
        delta_p_horiz = (p_final - p_init)[:2]
        p_drift = float(np.linalg.norm(delta_p_horiz) / max(p_norm_scale, 1e-30))
        results.append(GateResult(
            "momentum drift (horizontal only — substrate absorbs vertical)",
            p_drift <= float(g["momentum_drift_rel_max"]),
            f"|Δp_xy|/(m·max(v_rms,{V_FLOOR:.0e})) = {p_drift:.2e}, "
            f"|Δp_z|={abs(float((p_final - p_init)[2])):.2e} (substrate-leak, not gated), "
            f"v_rms={v_rms:.2e} (limit {g['momentum_drift_rel_max']:.0e})",
        ))
    else:
        p_drift = float(np.linalg.norm(p_final - p_init) / max(p_norm_scale, 1e-30))
        results.append(GateResult(
            "momentum drift",
            p_drift <= float(g["momentum_drift_rel_max"]),
            f"|Δp|/(m·max(v_rms,{V_FLOOR:.0e})) = {p_drift:.2e}, "
            f"v_rms={v_rms:.2e} (limit {g['momentum_drift_rel_max']:.0e})",
        ))

    # Stage 1a++ contract change: when ζ_star > 0 the active stress injects
    # energy by construction (active matter is non-equilibrium); the energy-
    # monotone gate is suspended and replaced with a finiteness check on
    # the cumulative active power. PI-approved 2026-04-29 in
    # docs/stage1a_plus_plus_layer2_sanity.md §3.
    if solver_cfg.zeta_star > 0.0:
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

        # (ii) Anchor force balance: substrate reaction = bulk pressure
        # transmitted into the contact band (Newton's 3rd, finite-Maxwell
        # tolerance ≤ 0.20). Read the median of the per-frame values from
        # frames after the same R-check time to skip the start transient.
        rel_errs_after = [
            r.get("anchor_force_balance_rel_err", float("nan"))
            for r in metrics_rows if r["time_star"] >= R_check_t
        ]
        rel_errs_clean = [e for e in rel_errs_after if not (e is None or np.isnan(e))]
        if rel_errs_clean:
            balance_med = float(np.median(rel_errs_clean))
            results.append(GateResult(
                "anchor force balance |F_sub − ∫σ_zz dA| / |F_sub|",
                balance_med <= 0.20,
                f"median over post-transient frames = {balance_med:.3f} (limit 0.20, "
                f"n={len(rel_errs_clean)} frames)",
            ))
        else:
            results.append(GateResult(
                "anchor force balance |F_sub − ∫σ_zz dA| / |F_sub|",
                False,
                "no valid anchor-force-balance samples",
            ))

        # (iii) Contact-band ρ_kernel / ρ_ref ∈ [0.85, 1.15]. Witnesses
        # substrate-anchored kernel (vs free-surface ~0.5·ρ_ref).
        rho_after = [
            r.get("rho_kernel_contact_over_ref", float("nan"))
            for r in metrics_rows if r["time_star"] >= R_check_t
        ]
        rho_clean = [e for e in rho_after if not (e is None or np.isnan(e))]
        if rho_clean:
            rho_med = float(np.median(rho_clean))
            results.append(GateResult(
                "contact-band ρ_kernel / ρ_ref ∈ [0.85, 1.15]",
                0.85 <= rho_med <= 1.15,
                f"median over post-transient frames = {rho_med:.3f} "
                f"(window [0.85, 1.15], n={len(rho_clean)} frames)",
            ))
        else:
            results.append(GateResult(
                "contact-band ρ_kernel / ρ_ref ∈ [0.85, 1.15]",
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

        # (ii) φ trajectory monotone toward predicted φ_eq.
        n_p = float(solver_cfg.n_particles)
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
            # Pilot duration is ~0.33 of φ relaxation time; allow up to
            # the full equilibrium offset as a generous tolerance (the gate
            # mainly catches divergence, not slow-relaxation behaviour).
            results.append(GateResult(
                "φ trajectory toward predicted φ_eq",
                phi_eq_err <= max(0.5, abs(solver_cfg.phi_initial - phi_eq_pred)),
                f"<φ>(end) = {phi_end:.4f}, predicted φ_eq = {phi_eq_pred:.4f}, "
                f"|err| = {phi_eq_err:.4f} (tolerance: ≤ max(0.5, "
                f"|φ_init − φ_eq|) = {max(0.5, abs(solver_cfg.phi_initial - phi_eq_pred)):.4f})",
            ))
        else:
            results.append(GateResult(
                "φ trajectory toward predicted φ_eq",
                False,
                "no φ samples",
            ))

        # (iii) A/A₀ trajectory finite & non-pathological.
        # A_contact_xy_hull is the substrate contact area projected to xy
        # (already computed by substrate_diagnostics). Pathology: NaN/Inf,
        # or A/A₀ contracts below 0.5 (spheroid disintegrating; spreading
        # context expects A/A₀ ≥ 1 monotone-ish, allowing small
        # equilibration transients).
        A_series = [
            r.get("contact_area_xy_hull", float("nan"))
            for r in metrics_rows
        ]
        A0 = next((a for a in A_series if not np.isnan(a) and a > 0), float("nan"))
        A_clean = [a for a in A_series if not (a is None or np.isnan(a))]
        if A_clean and A0 > 0:
            A_ratios = [a / A0 for a in A_clean]
            min_ratio = float(min(A_ratios))
            max_ratio = float(max(A_ratios))
            results.append(GateResult(
                "A/A₀ trajectory finite & non-pathological",
                np.isfinite(min_ratio) and np.isfinite(max_ratio) and min_ratio >= 0.5,
                f"A/A₀ ∈ [{min_ratio:.3f}, {max_ratio:.3f}] (limit min ≥ 0.5; "
                f"A₀ = {A0:.4f}, n_frames = {len(A_ratios)})",
            ))
        else:
            results.append(GateResult(
                "A/A₀ trajectory finite & non-pathological",
                False,
                "no valid A_contact_xy_hull samples",
            ))

    psi_check_t = float(g["sphericity_check_after_s"]) / float(cfg["physics"]["maxwell_tau_s"])
    psi_after = [r["wadell_sphericity"] for r in metrics_rows if r["time_star"] >= psi_check_t]
    if psi_after:
        psi_min_seen = float(np.nanmin(psi_after))
        results.append(GateResult(
            "Wadell sphericity ψ",
            psi_min_seen >= float(g["sphericity_min"]),
            f"min ψ after equilibration = {psi_min_seen:.3f} (limit {g['sphericity_min']:.3f})",
        ))
    else:
        results.append(GateResult("Wadell sphericity ψ", False, "no frames after check_after time"))

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
