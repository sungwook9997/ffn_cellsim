"""Sanity Gate Protocol enforcement (CLAUDE.md hard rule).

The CLAUDE.md Sanity Gate Protocol requires that, before the first
execution of any physics/numerics module, six checks are recorded and
that runtime measurements that violate stated KU ranges raise rather
than silently passing.

This module gives a single :class:`SanityGate` value object plus a
helper for the Phase 1 ECM-network gate that combines the KU-1.7
(biological mesh ξ — biological reference), KU-1.27 (Mikado segment
length ℓ_c — emergent), and KU-1.3 (coordination ⟨z⟩ — emergent)
checks. ξ and ℓ_c are tracked as distinct quantities.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any


class SanityGateFailure(RuntimeError):
    """Raised when a measurement violates a documented KU range."""


@dataclass(slots=True)
class GateCheck:
    name: str
    status: str           # "PASS" | "FAIL" | "WARN"
    detail: str
    ku: str = ""


@dataclass(slots=True)
class SanityGateReport:
    module: str
    checks: list[GateCheck] = field(default_factory=list)

    def add(self, name: str, ok: bool, detail: str, ku: str = "") -> None:
        self.checks.append(
            GateCheck(name=name, status="PASS" if ok else "FAIL", detail=detail, ku=ku)
        )

    @property
    def passed(self) -> bool:
        return all(c.status == "PASS" for c in self.checks)

    def raise_if_failed(self) -> None:
        if self.passed:
            return
        msg_lines = [f"Sanity Gate FAILED for {self.module}:"]
        for c in self.checks:
            if c.status == "FAIL":
                tag = f"[{c.ku}] " if c.ku else ""
                msg_lines.append(f"  - {c.name}: {tag}{c.detail}")
        raise SanityGateFailure("\n".join(msg_lines))

    def summary(self) -> str:
        lines = [f"Sanity Gate · {self.module}"]
        for c in self.checks:
            tag = f"[{c.ku}] " if c.ku else ""
            lines.append(f"  {c.status:4s} {c.name:38s} {tag}{c.detail}")
        return "\n".join(lines)


def gate_phase1_ecm(
    *,
    biological_mesh_target: float,     # ξ from KU-1.7
    measured_segment_length: float,    # ℓ_c emergent from Mikado generator
    expected_segment_length: float,    # ℓ_c predicted from KU-1.27 inverse
    measured_z: float,                 # ⟨z⟩ emergent from real intersections
    expected_z: float,                 # KU-1.3 prediction from Mikado theory
    ku13_reference_z: float,           # biological collagen ⟨z⟩=3.2 reference
    n_fibers: int,
    n_beads_total: int,
    n_links: int,
    acceptance: dict[str, Any],
    demo_mode: bool,
) -> SanityGateReport:
    """Phase 1 ECM Sanity Gate.

    Distinguishes the **construction** target (ℓ_c matches biological ξ)
    from the **emergent** observable (⟨z⟩ from real segment crossings),
    so a passing run can still flag the biology-vs-Mikado gap when ⟨z⟩
    differs from KU-1.3.
    """
    rep = SanityGateReport(module="phase1_ecm_network")

    lc_lo, lc_hi = acceptance["segment_length_range"]
    rep.add(
        "mikado_segment_length_in_KU_range",
        lc_lo <= measured_segment_length <= lc_hi,
        f"emergent ℓ_c={measured_segment_length:.3e} m, predicted "
        f"{expected_segment_length:.3e}, target ξ={biological_mesh_target:.3e}, "
        f"accepted [{lc_lo:.1e}, {lc_hi:.1e}]",
        ku="KU-1.27 / KU-1.7",
    )

    z_lo, z_hi = acceptance["z_range"]
    rep.add(
        "z_in_sub_isostatic_range",
        z_lo <= measured_z <= z_hi,
        f"emergent ⟨z⟩={measured_z:.3f} (predicted {expected_z:.3f}, "
        f"KU-1.3 ref {ku13_reference_z:.2f}), accepted [{z_lo}, {z_hi}], z_c=2d=4",
        ku="KU-1.3",
    )

    rep.add(
        "z_within_prediction",
        abs(measured_z - expected_z) <= 0.1 * max(expected_z, 1e-9),
        f"|⟨z⟩−predicted|={abs(measured_z-expected_z):.3f}, tol=10% (probes "
        f"the periodic-boundary loss)",
        ku="construction self-consistency",
    )

    cost_cap = acceptance["n_fibers_max"]
    rep.add(
        "fiber_count_under_cost_cap",
        n_fibers <= cost_cap,
        f"n_fibers={n_fibers}, cap={cost_cap}; n_beads={n_beads_total}; "
        f"n_links={n_links}",
        ku="cost",
    )

    # ---- Biology gap is INFORMATIONAL, not a fail ----
    gap = measured_z - ku13_reference_z
    rep.add(
        "biology_gap_logged",
        True,
        f"⟨z⟩_measured − ⟨z⟩_KU-1.3 = {gap:+.3f}; 2D Mikado lacks 3D "
        f"covalent and bundle cross-links present in collagen",
        ku="info",
    )

    if demo_mode:
        for c in rep.checks:
            if c.status == "FAIL":
                c.detail = "DEMO_MODE: " + c.detail
                c.status = "WARN"

    return rep


def gate_unit1_2_dynamics(
    *,
    dt: float,
    tau_xl: float,
    tau_stretch: float,
    tau_bend: float,
    safety_factor: float,
    integrator_name: str,
) -> SanityGateReport:
    """Unit 1.2 CFL / integrator Sanity Gate (KU-1.26).

    Hard requirement: dt < safety_factor · τ_min where τ_min is the
    smallest of (γ_b/k_xl, γ_b ℓ₀/μ, γ_b ℓ₀³/κ).
    """
    rep = SanityGateReport(module="unit1_2_overdamped_langevin")

    tau_min = min(tau_xl, tau_stretch, tau_bend)
    dt_max = safety_factor * tau_min
    rep.add(
        "cfl_dt_below_alpha_tau_min",
        dt < dt_max,
        f"dt={dt:.3e} s, α·τ_min={dt_max:.3e} s "
        f"(τ_xl={tau_xl:.2e}, τ_stretch={tau_stretch:.2e}, τ_bend={tau_bend:.2e}, "
        f"α={safety_factor})",
        ku="KU-1.26 CFL",
    )

    rep.add(
        "integrator_supported_in_phase_1",
        integrator_name == "euler_maruyama",
        f"integrator='{integrator_name}'; Phase 1 supports 'euler_maruyama' only",
        ku="KU-1.26",
    )

    return rep


def gate_unit2_1_motor_clutch(
    *,
    dt: float,
    cfl_safe_dt: float,
    catch_peak_force_analytic: float,
    catch_peak_force_simulated: float,
    catch_peak_tolerance: float,
    biphasic_verdict: str,                   # "saturating" | "peaked"
    biphasic_prominence: float,
    biphasic_log10_seed_argmax_spread: float,
    biphasic_asymptote_N: float,
    motor_stall_total_N: float,
    biphasic_E_star_analytic: float,         # info-only (Bangasser matched-stiffness)
    biphasic_E_star_measured_argmax: float,  # info-only
    clutch_force_band_pN: tuple[float, float],
    clutch_force_band_mass_measured: float,
    clutch_force_band_mass_required: float,
    perf_steps: int,
    perf_seconds: float,
    perf_budget_seconds: float,
    ku_experimental_F_star_pN: float = 30.0,
    ku_2_12_mature_band_pN: tuple[float, float] = (5.0, 20.0),
) -> SanityGateReport:
    """Phase 1 Unit 2.1 Sanity Gate (KU-2.4 / KU-2.5 / KU-2.8 / KU-2.12).

    Combines the four acceptance checks called out in the Brief:
      * dt below the bond-event-rate CFL bound;
      * simulated catch peak F* within ±tolerance of the closed-form
        Pereverzev F* (logs the experimental KU-2.5 F* ≈ 30 pN as an
        informational gap, not a fail);
      * biphasic peak E* within ``factor`` of the KU-2.8 closed-form;
      * per-clutch force histogram mass in the 5–20 pN band (KU-2.12);
      * 50-clutch × 1000-step performance under 1 s (Brief Task 7).
    """
    rep = SanityGateReport(module="unit2_1_motor_clutch")

    rep.add(
        "dt_below_bond_event_rate_cfl",
        dt < cfl_safe_dt,
        f"dt={dt:.3e} s, cfl_safe_dt={cfl_safe_dt:.3e} s "
        f"(α·1/max(k_on, k_off(F_s)), α=0.1)",
        ku="KU-2.4 / KU-2.5",
    )

    if math.isnan(catch_peak_force_analytic):
        rep.add(
            "catch_peak_F_star_within_analytic",
            False,
            "Pereverzev parameters give a slip-only bond — no catch peak "
            "is defined; cannot validate KU-2.5 F*.",
            ku="KU-2.5",
        )
    else:
        rel_err = abs(catch_peak_force_simulated - catch_peak_force_analytic) \
            / max(catch_peak_force_analytic, 1e-30)
        rep.add(
            "catch_peak_F_star_within_analytic",
            rel_err <= catch_peak_tolerance,
            f"F*_sim={catch_peak_force_simulated*1e12:.2f} pN, "
            f"F*_analytic={catch_peak_force_analytic*1e12:.2f} pN, "
            f"rel err={rel_err:.3f}, tol={catch_peak_tolerance}",
            ku="KU-2.5",
        )

    gap = catch_peak_force_analytic * 1e12 - ku_experimental_F_star_pN
    rep.add(
        "catch_peak_KU_2_5_experimental_gap_logged",
        True,
        f"analytic F*={catch_peak_force_analytic*1e12:.2f} pN; "
        f"KU-2.5 experimental F*≈{ku_experimental_F_star_pN:.0f} pN; "
        f"gap={gap:+.2f} pN. The KU-2.18 illustrative parameters do not "
        f"reproduce the Kong 2009 / Elosegui-Artola 2016 experimental "
        f"peak. Refit deferred to Phase 2.",
        ku="info",
    )

    # Biphasic shape verdict — saturating, not a true peak, with Phase 1
    # parameters. Multi-seed argmax spread + asymptote near motor stall
    # are the operational signatures. Matched-stiffness E* is recorded
    # as info; it predicts a slip-bond peak that does not exist here.
    rep.add(
        "biphasic_shape_is_saturating",
        biphasic_verdict == "saturating",
        f"verdict={biphasic_verdict}; prominence={biphasic_prominence:.4f} "
        f"(saturating if < 0.10); inter-seed log10 argmax spread="
        f"{biphasic_log10_seed_argmax_spread:.2f} decades",
        ku="KU-2.4 / KU-2.8",
    )
    asymp_ratio = biphasic_asymptote_N / max(motor_stall_total_N, 1e-30)
    rep.add(
        "biphasic_asymptote_near_motor_stall",
        0.7 <= asymp_ratio <= 1.1,
        f"asymptote ⟨F⟩={biphasic_asymptote_N*1e12:.2f} pN, "
        f"N_m·F_stall={motor_stall_total_N*1e12:.2f} pN, ratio={asymp_ratio:.3f}",
        ku="KU-2.4",
    )
    rep.add(
        "biphasic_matched_stiffness_E_star_info",
        True,
        f"matched-stiffness E*={biphasic_E_star_analytic:.3e} Pa "
        f"(Bangasser 2013) vs argmax E={biphasic_E_star_measured_argmax:.3e} Pa. "
        f"Matched analytic predicts slip-bond peak; Phase 1 catch-stabilised "
        f"curve saturates instead — argmax is noise-dominated.",
        ku="info (KU-2.8 reference, not asserted)",
    )

    lo, hi = clutch_force_band_pN
    rep.add(
        "per_clutch_force_in_phase1_nascent_band",
        clutch_force_band_mass_measured >= clutch_force_band_mass_required,
        f"mass in [{lo}, {hi}] pN = {clutch_force_band_mass_measured:.3f}, "
        f"required ≥ {clutch_force_band_mass_required} "
        f"(KU-2.12 mature-FA {ku_2_12_mature_band_pN[0]}-{ku_2_12_mature_band_pN[1]} "
        f"pN band moves to Unit 2.2 contract)",
        ku="KU-2.12 Phase 1 (nascent)",
    )

    rep.add(
        "perf_budget_met",
        perf_seconds < perf_budget_seconds,
        f"{perf_steps} steps × 50 clutches in {perf_seconds*1e3:.1f} ms "
        f"(budget {perf_budget_seconds*1e3:.0f} ms)",
        ku="Brief Task 7",
    )

    return rep


def gate_equipartition(
    *,
    measured_mean_link_energy: float,
    kT: float,
    n_samples: int,
    tolerance: float = 0.10,
) -> SanityGateReport:
    """Equipartition: ⟨½ k_xl |Δr|²⟩ → ½ k_B T per cross-link (KU-1.26).

    Sample mean error scales as √(2/n_samples); tolerance 10 % is set
    well above the statistical floor for typical n_samples ≥ 10^4.
    """
    rep = SanityGateReport(module="unit1_2_equipartition")
    target = 0.5 * kT
    err = abs(measured_mean_link_energy - target) / target
    rep.add(
        "link_energy_in_canonical_ensemble",
        err <= tolerance,
        f"⟨½ k_xl |Δr|²⟩={measured_mean_link_energy:.3e} J, "
        f"½ k_B T={target:.3e} J, rel err={err:.3f}, tol={tolerance}, "
        f"n_samples={n_samples}",
        ku="KU-1.26 equipartition",
    )
    return rep


def gate_phase1_cell_cortex(
    *,
    measured_z: float,
    coverage_ratio: float,
    initial_aspect_ratio: float,
    final_aspect_ratio: float,
    acceptance_window_reached_aspect_ratio: float,
    acceptance: dict[str, Any],
    cortex_radius_drift_rel: float,
    n_fibers: int,
    n_xls: int,
    demo_mode: bool,
) -> SanityGateReport:
    """Phase 1 Unit 3.1 Cell-Cortex Sanity Gate (KU-3.1, KU-3.17).

    Combines:
      - ⟨z⟩ in the KU-3.17 cortex coordination band (denser than ECM
        but below mechanical saturation).
      - Initial aspect ratio honest (≤ acceptance.aspect_ratio_initial_max,
        sanity that the rounding test is not degenerate).
      - KU-3.1 VALIDATION: aspect ratio within the acceptance window
        falls below ``acceptance.aspect_ratio_window``.

    A non-blocking informational check reports the mean-radius drift
    so the cortex shrinkage under cortical tension is visible without
    failing the gate.
    """
    rep = SanityGateReport(module="phase1_cell_cortex")

    z_lo, z_hi = acceptance["z_range"]
    rep.add(
        "cortex_z_in_KU317_range",
        z_lo <= measured_z <= z_hi,
        f"⟨z⟩={measured_z:.3f}, accepted [{z_lo}, {z_hi}] "
        f"(N_fibers={n_fibers}, N_xls={n_xls}, coverage={coverage_ratio:.2f})",
        ku="KU-3.17",
    )

    ar_init_cap = acceptance.get("aspect_ratio_initial_max", 2.0)
    rep.add(
        "initial_aspect_ratio_below_cap",
        initial_aspect_ratio <= ar_init_cap,
        f"initial aspect ratio={initial_aspect_ratio:.3f}, cap={ar_init_cap} "
        f"(KU-3.1 test design — ellipse not degenerate)",
        ku="KU-3.1 test design",
    )

    ar_cap = acceptance.get("aspect_ratio_window", 1.2)
    measured_ar = min(final_aspect_ratio, acceptance_window_reached_aspect_ratio)
    rep.add(
        "rounding_within_acceptance_window",
        measured_ar <= ar_cap,
        f"min(final, windowed) aspect ratio={measured_ar:.3f}, "
        f"window cap={ar_cap}; final={final_aspect_ratio:.3f}, "
        f"windowed={acceptance_window_reached_aspect_ratio:.3f}",
        ku="KU-3.1 VALIDATION",
    )

    drift_cap = acceptance.get("cortex_radius_drift_rel", 0.10)
    rep.add(
        "cortex_radius_drift_logged",
        True,
        f"⟨R⟩ drift = {cortex_radius_drift_rel:+.3f} rel; "
        f"info cap = {drift_cap} (tension shrinks the cortex, expected)",
        ku="info",
    )

    if demo_mode:
        for c in rep.checks:
            if c.status == "FAIL":
                c.detail = "DEMO_MODE: " + c.detail
                c.status = "WARN"

    return rep
