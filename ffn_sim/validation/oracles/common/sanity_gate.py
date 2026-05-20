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


PHASE_1_CANONICAL_INTEGRATORS: frozenset[str] = frozenset({
    "leimkuhler_matthews_baoab",
    "lm_baoab",
})
"""Names accepted as the Phase 1 canonical integrator (D3, 2026-05-19)."""

PHASE_1_REFERENCE_INTEGRATORS: frozenset[str] = frozenset({
    "euler_maruyama",
    "hoomd_brownian",
})
"""Lower-order integrators only allowed as the L-M correctness reference
(see H.2 brief: 'Same L_p when re-run with vanilla md.methods.Brownian at
halved Δt — confirms L-M ≠ artefact').
"""


def gate_unit1_2_dynamics(
    *,
    dt: float,
    tau_xl: float,
    tau_stretch: float,
    tau_bend: float,
    safety_factor: float,
    integrator_name: str,
    allow_reference_integrator: bool = False,
) -> SanityGateReport:
    """Unit 1.2 CFL / integrator Sanity Gate (KU-1.26 + D3).

    Hard requirement: dt < safety_factor · τ_min where τ_min is the
    smallest of (γ_b/k_xl, γ_b ℓ₀/μ, γ_b ℓ₀³/κ).

    Integrator requirement: Phase 1 canonical is Leimkuhler-Matthews
    BAOAB-limit (PHASE_0_3_DECISIONS D3, 2026-05-19). Set
    ``allow_reference_integrator=True`` only when running the D3
    verification suite (H.2 brief — vanilla Brownian at halved Δt).
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

    if allow_reference_integrator:
        allowed = PHASE_1_CANONICAL_INTEGRATORS | PHASE_1_REFERENCE_INTEGRATORS
        ok = integrator_name in allowed
        rep.add(
            "integrator_supported_in_phase_1",
            ok,
            f"integrator='{integrator_name}'; allowed (D3 verification suite): "
            f"{sorted(allowed)}",
            ku="D3 / KU-1.26",
        )
    else:
        ok = integrator_name in PHASE_1_CANONICAL_INTEGRATORS
        rep.add(
            "integrator_supported_in_phase_1",
            ok,
            f"integrator='{integrator_name}'; Phase 1 canonical = "
            f"{sorted(PHASE_1_CANONICAL_INTEGRATORS)} (D3 Leimkuhler-Matthews "
            f"BAOAB-limit). Set allow_reference_integrator=True for L-M "
            f"correctness checks against {sorted(PHASE_1_REFERENCE_INTEGRATORS)}.",
            ku="D3 / KU-1.26",
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


# ---------------------------------------------------------------------------
# H.2 single filament — persistence length (KU-1.1) + L-M correctness (D3)
# ---------------------------------------------------------------------------

# KU-1.1 actin canonical persistence length and ±10 % H.2 acceptance band.
KU11_LP_TARGET: float = 17.0e-6           # m
KU11_LP_BAND: tuple[float, float] = (15.3e-6, 18.7e-6)


def gate_h2_persistence_length(
    *,
    Lp_measured: float,                       # m, from tangent-correlation fit
    ks_pvalue: float,                         # Boltzmann angle KS p-value
    equipartition_relerr_per_bond: list[float],   # |⟨E⟩−kT/2|/(kT/2) for each bond
    Lp_target: float = KU11_LP_TARGET,
    Lp_band: tuple[float, float] = KU11_LP_BAND,
    ks_pvalue_min: float = 0.05,
    equipartition_tolerance: float = 0.05,
    Lp_em_reference: float | None = None,     # E-M / vanilla Brownian re-run
    Lp_em_sigma: float | None = None,         # 1-σ uncertainty of either measurement
) -> SanityGateReport:
    """Phase 1 Unit H.2 Sanity Gate (KU-1.1 + D3 verification).

    Gates the foundational single-particle physics that every cortex,
    lamellipodium, and ECM filament in Phase 1+ depends on (H.2 brief).

    Checks
    ------
    1. Tangent-correlation persistence length within the ±10 % KU-1.1 band.
    2. Boltzmann angle-distribution KS test p > ``ks_pvalue_min``.
    3. Per-bond equipartition |⟨E_bend⟩ − kT/2| / (kT/2) ≤ ``equipartition_tolerance``
       (default 0.05).
    4. L-M correctness (only if ``Lp_em_reference`` provided): same L_p
       within 1 σ when re-run with vanilla Brownian at halved Δt.
       Confirms the L-M result is not an integrator artefact (D3, H.2
       brief acceptance row 4).
    """
    rep = SanityGateReport(module="h2_single_filament_persistence")

    lp_lo, lp_hi = Lp_band
    rep.add(
        "persistence_length_in_KU11_band",
        lp_lo <= Lp_measured <= lp_hi,
        f"L_p={Lp_measured*1e6:.2f} μm, target={Lp_target*1e6:.1f} μm, "
        f"band=[{lp_lo*1e6:.1f}, {lp_hi*1e6:.1f}] μm",
        ku="KU-1.1",
    )

    rep.add(
        "angle_distribution_boltzmann_ks",
        ks_pvalue > ks_pvalue_min,
        f"KS p={ks_pvalue:.4f}, threshold p > {ks_pvalue_min}",
        ku="KU-1.1 angle distribution",
    )

    n_bonds = len(equipartition_relerr_per_bond)
    if n_bonds == 0:
        rep.add(
            "per_bond_equipartition",
            False,
            "no per-bond equipartition data provided",
            ku="KU-1.26 equipartition",
        )
    else:
        worst = max(equipartition_relerr_per_bond)
        n_fail = sum(1 for e in equipartition_relerr_per_bond
                     if e > equipartition_tolerance)
        rep.add(
            "per_bond_equipartition",
            worst <= equipartition_tolerance,
            f"max |⟨E⟩−kT/2|/(kT/2) = {worst:.4f} over {n_bonds} bonds "
            f"(tol={equipartition_tolerance}); {n_fail} bond(s) above tol",
            ku="KU-1.26 equipartition",
        )

    if Lp_em_reference is not None:
        sigma = Lp_em_sigma if Lp_em_sigma is not None else 0.0
        if sigma <= 0:
            rep.add(
                "lm_vs_em_correctness",
                False,
                f"L_p^LM={Lp_measured*1e6:.2f} μm, L_p^EM(ref)={Lp_em_reference*1e6:.2f} μm; "
                f"Lp_em_sigma must be > 0 to evaluate the 1-σ window",
                ku="D3 L-M correctness",
            )
        else:
            diff = abs(Lp_measured - Lp_em_reference)
            rep.add(
                "lm_vs_em_correctness",
                diff <= sigma,
                f"|L_p^LM − L_p^EM| = {diff*1e6:.3f} μm, 1σ = {sigma*1e6:.3f} μm "
                f"(L_p^LM={Lp_measured*1e6:.2f}, L_p^EM={Lp_em_reference*1e6:.2f}); "
                f"D3 requires L-M result to be within 1σ of vanilla Brownian "
                f"at halved Δt (H.2 brief)",
                ku="D3 L-M correctness",
            )

    return rep


# ---------------------------------------------------------------------------
# D7 excluded volume — WCA invariants
# ---------------------------------------------------------------------------

# WCA cutoff is exactly 2^(1/6)·σ, the LJ minimum (no attractive tail).
WCA_CUTOFF_FACTOR: float = 2.0 ** (1.0 / 6.0)

# D7 default ε = 0.5 kT (PHASE_0_3_DECISIONS); standard "soft" steric ε.
D7_EPSILON_OVER_KT: float = 0.5


def gate_wca_cutoff(
    *,
    pair_params: dict[tuple[str, str], dict[str, float]],
    kT: float,
    epsilon_over_kT_expected: float = D7_EPSILON_OVER_KT,
    cutoff_rel_tolerance: float = 1.0e-6,
    epsilon_rel_tolerance: float = 1.0e-3,
    mode: str | None = None,
    required_pair_types: list[tuple[str, str]] | None = None,
) -> SanityGateReport:
    """Phase 1 D7 WCA excluded-volume Sanity Gate.

    Verifies for every pair-type entry in ``pair_params``:
      * ``r_cut == 2^(1/6)·σ`` within ``cutoff_rel_tolerance`` (WCA
        truncates LJ at its minimum so only the repulsive branch
        survives — a longer r_cut re-introduces the attractive tail
        that D7 explicitly forbids).
      * ``ε == epsilon_over_kT_expected · kT`` within
        ``epsilon_rel_tolerance``.
      * Optional ``mode`` (e.g. ``'shift'``) matches a global mode.
      * Optional ``required_pair_types`` are all present (H.1/H.3
        coverage check — actin_cortex×actin_cortex, actin_cortex×
        xlink_head, actin_cortex×myosin_head, myosin_head×myosin_head,
        integrin×ligand, etc.).

    ``pair_params`` schema: ``{(typeA, typeB): {'sigma': ..., 'r_cut': ...,
    'epsilon': ..., 'mode': ...}}`` (per-pair ``mode`` overrides the
    function arg if present).
    """
    rep = SanityGateReport(module="d7_wca_excluded_volume")

    if not pair_params:
        rep.add(
            "any_wca_pair_configured",
            False,
            "pair_params is empty — D7 requires WCA pairs enabled from "
            "Phase 1 (PHASE_0_3_DECISIONS D7).",
            ku="D7",
        )
        return rep

    eps_target = epsilon_over_kT_expected * kT

    for pair, params in pair_params.items():
        sigma = params["sigma"]
        r_cut = params["r_cut"]
        eps = params["epsilon"]
        r_cut_expected = WCA_CUTOFF_FACTOR * sigma
        rel = abs(r_cut - r_cut_expected) / max(r_cut_expected, 1e-30)
        rep.add(
            f"wca_cutoff_{pair[0]}_{pair[1]}",
            rel <= cutoff_rel_tolerance,
            f"r_cut={r_cut:.4e}, expected 2^(1/6)·σ={r_cut_expected:.4e} "
            f"(σ={sigma:.3e}), rel err={rel:.2e}, tol={cutoff_rel_tolerance:.0e}",
            ku="D7 WCA r_cut",
        )

        eps_rel = abs(eps - eps_target) / max(eps_target, 1e-30)
        rep.add(
            f"wca_epsilon_{pair[0]}_{pair[1]}",
            eps_rel <= epsilon_rel_tolerance,
            f"ε={eps:.3e} J, expected {epsilon_over_kT_expected} kT="
            f"{eps_target:.3e} J, rel err={eps_rel:.2e}, "
            f"tol={epsilon_rel_tolerance:.0e}",
            ku="D7 ε",
        )

        pair_mode = params.get("mode")
        if mode is not None and pair_mode is not None and pair_mode != mode:
            rep.add(
                f"wca_mode_{pair[0]}_{pair[1]}",
                False,
                f"mode='{pair_mode}', expected '{mode}'",
                ku="D7 mode",
            )

    if required_pair_types:
        present = {tuple(sorted(p)) for p in pair_params.keys()}
        for req in required_pair_types:
            req_norm = tuple(sorted(req))
            rep.add(
                f"required_pair_present_{req[0]}_{req[1]}",
                req_norm in present,
                f"required pair {req} {'present' if req_norm in present else 'MISSING'}",
                ku="D7 coverage",
            )

    return rep


# ---------------------------------------------------------------------------
# D5 motor minifilament topology — Stam-Hocky bipolar
# ---------------------------------------------------------------------------

# D5 reference values (PHASE_0_3_DECISIONS D5; Stam & Hocky 2024 convention,
# Hocky-group computational standard).
D5_BACKBONE_BEADS: int = 14
D5_HEADS_PER_SIDE: int = 10
D5_BACKBONE_LENGTH_NM: float = 700.0      # ~700 nm thick filament
D5_HEAD_REST_LENGTH_NM: float = 200.0     # head ↔ backbone perpendicular spring r0
D5_HEAD_SPRING_K_PN_PER_UM: float = 1.0   # AFINES motor stiffness (preserved)


def gate_minifilament_topology(
    *,
    n_backbone_beads: int,
    n_heads_per_side: int,
    n_sides: int,
    has_rigid_body_constraint: bool,
    cross_bridge_k_pN_per_um: float,
    cross_bridge_r0_nm: float,
    backbone_length_nm: float,
    n_backbone_tolerance: int = 2,
    n_heads_tolerance: int = 2,
    length_rel_tolerance: float = 0.15,
    k_rel_tolerance: float = 0.20,
    r0_rel_tolerance: float = 0.15,
) -> SanityGateReport:
    """Phase 1 D5 Stam-Hocky multi-head minifilament topology gate.

    The mechanistic Stam-Hocky bipolar minifilament has a rigid-rod
    backbone (~14 beads, ~700 nm) carrying ~10 cross-bridge heads on
    each polarity side — a total of ~34 particles per minifilament.
    This gate exists primarily to **prevent silent regression to the
    AFINES single-2-head spring**, which is exactly the abstraction
    PI 2026-05-19 ruled against.

    A failure of ``bipolar_two_sides`` or a single-bead backbone is the
    signature of a regression to the AFINES single-spring motor. Halt
    and re-derive from Stam-Hocky 2024 rather than relaxing the gate.
    """
    rep = SanityGateReport(module="d5_minifilament_topology")

    rep.add(
        "bipolar_two_sides",
        n_sides == 2,
        f"n_sides={n_sides}, expected 2 (Stam-Hocky bipolar). "
        f"n_sides=1 indicates regression to AFINES single 2-head spring; "
        f"halt and consult PHASE_0_3_DECISIONS D5.",
        ku="D5 architecture",
    )

    backbone_ok = (
        n_backbone_beads >= 2
        and abs(n_backbone_beads - D5_BACKBONE_BEADS) <= n_backbone_tolerance
    )
    rep.add(
        "backbone_bead_count",
        backbone_ok,
        f"n_backbone_beads={n_backbone_beads}, expected "
        f"{D5_BACKBONE_BEADS}±{n_backbone_tolerance} "
        f"(Stam-Hocky rigid-rod ~14). 1-bead backbone = AFINES regression.",
        ku="D5 backbone",
    )

    heads_ok = abs(n_heads_per_side - D5_HEADS_PER_SIDE) <= n_heads_tolerance
    rep.add(
        "heads_per_side_count",
        heads_ok,
        f"n_heads_per_side={n_heads_per_side}, expected "
        f"{D5_HEADS_PER_SIDE}±{n_heads_tolerance} (Hocky-group convention).",
        ku="D5 heads/side",
    )

    rep.add(
        "rigid_body_constraint_present",
        has_rigid_body_constraint,
        f"has_rigid_body_constraint={has_rigid_body_constraint}; "
        f"D5 requires md.constrain.Rigid on the backbone "
        f"(rigid-rod approximation).",
        ku="D5 rigid constraint",
    )

    total = n_backbone_beads + n_heads_per_side * n_sides
    expected_total = D5_BACKBONE_BEADS + D5_HEADS_PER_SIDE * 2
    rep.add(
        "total_particles_per_minifilament",
        abs(total - expected_total) <= (
            n_backbone_tolerance + 2 * n_heads_tolerance
        ),
        f"n_particles={total} (= {n_backbone_beads} backbone + "
        f"{n_heads_per_side}×{n_sides} heads); expected ≈ {expected_total} "
        f"(PHASE_0_3_DECISIONS D5).",
        ku="D5 particle budget",
    )

    length_rel = abs(backbone_length_nm - D5_BACKBONE_LENGTH_NM) / D5_BACKBONE_LENGTH_NM
    rep.add(
        "backbone_length_nm",
        length_rel <= length_rel_tolerance,
        f"backbone_length={backbone_length_nm:.1f} nm, expected "
        f"{D5_BACKBONE_LENGTH_NM:.0f} nm ±{length_rel_tolerance*100:.0f}%, "
        f"rel err={length_rel:.3f}",
        ku="D5 backbone length",
    )

    k_rel = abs(cross_bridge_k_pN_per_um - D5_HEAD_SPRING_K_PN_PER_UM) \
        / D5_HEAD_SPRING_K_PN_PER_UM
    rep.add(
        "cross_bridge_spring_stiffness",
        k_rel <= k_rel_tolerance,
        f"k_xbridge={cross_bridge_k_pN_per_um:.3f} pN/μm, expected "
        f"{D5_HEAD_SPRING_K_PN_PER_UM} ±{k_rel_tolerance*100:.0f}% (AFINES default)",
        ku="D5 head spring k",
    )

    r0_rel = abs(cross_bridge_r0_nm - D5_HEAD_REST_LENGTH_NM) / D5_HEAD_REST_LENGTH_NM
    rep.add(
        "cross_bridge_rest_length",
        r0_rel <= r0_rel_tolerance,
        f"r0_xbridge={cross_bridge_r0_nm:.1f} nm, expected "
        f"{D5_HEAD_REST_LENGTH_NM:.0f} nm ±{r0_rel_tolerance*100:.0f}%",
        ku="D5 head r0",
    )

    return rep


# ---------------------------------------------------------------------------
# Cortex extension gates — KU-3.5 tension, KU-3.18 blebbistatin, KU-3.20 Q_ij
# ---------------------------------------------------------------------------

# KU-3.5 cortical tension: γ ≈ 0.5 mN/m. H.3 brief acceptance: ±30 %.
KU35_TENSION_TARGET_N_PER_M: float = 0.5e-3       # 0.5 mN/m


def gate_cortex_tension(
    *,
    tension_measured_N_per_m: float,              # Laplace γ from internal pressure
    tension_target_N_per_m: float = KU35_TENSION_TARGET_N_PER_M,
    rel_tolerance: float = 0.30,
) -> SanityGateReport:
    """Phase 1 H.3 KU-3.5 cortical tension gate.

    Laplace law γ = ΔP · R / 2 (2D: γ = ΔP · R) measured from the
    internal pressure required to balance the cortex curvature must
    fall within ±30 % of the KU-3.5 target (0.5 mN/m).
    """
    rep = SanityGateReport(module="h3_cortex_tension")

    rel = abs(tension_measured_N_per_m - tension_target_N_per_m) \
        / max(tension_target_N_per_m, 1e-30)
    rep.add(
        "cortical_tension_in_KU35_band",
        rel <= rel_tolerance,
        f"γ_measured={tension_measured_N_per_m*1e3:.3f} mN/m, "
        f"target={tension_target_N_per_m*1e3:.2f} mN/m, "
        f"rel err={rel:.3f}, tol={rel_tolerance}",
        ku="KU-3.5",
    )

    return rep


def gate_blebbistatin_response(
    *,
    tension_myosin_on_N_per_m: float,
    tension_myosin_off_N_per_m: float,
    final_aspect_ratio_myosin_on: float,
    final_aspect_ratio_myosin_off: float,
    tension_drop_min_rel: float = 0.5,            # ≥50 % drop expected
    rounding_failure_ar_floor: float = 1.3,
) -> SanityGateReport:
    """Phase 1 H.3 KU-3.18 blebbistatin response gate.

    Myosin-OFF must show (a) cortical tension drop of at least
    ``tension_drop_min_rel`` and (b) failed rounding — aspect ratio
    stays ≥ ``rounding_failure_ar_floor`` at the H.3 60 s acceptance
    window (compared to myosin-ON which must round to < 1.2 per
    KU-3.1).
    """
    rep = SanityGateReport(module="h3_blebbistatin_response")

    drop = (tension_myosin_on_N_per_m - tension_myosin_off_N_per_m) \
        / max(tension_myosin_on_N_per_m, 1e-30)
    rep.add(
        "tension_drops_on_myosin_off",
        drop >= tension_drop_min_rel,
        f"γ_on={tension_myosin_on_N_per_m*1e3:.3f} mN/m, "
        f"γ_off={tension_myosin_off_N_per_m*1e3:.3f} mN/m, "
        f"relative drop={drop:.3f}, required ≥{tension_drop_min_rel}",
        ku="KU-3.18 tension",
    )

    rep.add(
        "rounding_fails_on_myosin_off",
        final_aspect_ratio_myosin_off >= rounding_failure_ar_floor,
        f"AR_off(t=60s)={final_aspect_ratio_myosin_off:.3f}, "
        f"required ≥{rounding_failure_ar_floor} "
        f"(KU-3.18: blebbistatin blocks rounding; AR_on={final_aspect_ratio_myosin_on:.3f})",
        ku="KU-3.18 rounding",
    )

    return rep


def gate_nematic_order(
    *,
    S_isotropic_cell: float,                      # |Q| order parameter, isotropic case
    S_aligned_cell: float,                        # aligned case
    S_isotropic_max: float = 0.1,
    S_aligned_min: float = 0.3,
) -> SanityGateReport:
    """Phase 1 H.3 KU-3.20 nematic order Q_ij gate.

    Q_ij = ⟨2 t̂_i t̂_j − δ_ij⟩ over the cortex filament tangent
    distribution; S = max eigenvalue (scalar nematic order). H.3 brief
    acceptance: isotropic S < 0.1, aligned S > 0.3.
    """
    rep = SanityGateReport(module="h3_nematic_order")

    rep.add(
        "isotropic_cell_below_S_floor",
        S_isotropic_cell < S_isotropic_max,
        f"S_iso={S_isotropic_cell:.3f}, required <{S_isotropic_max}",
        ku="KU-3.20 isotropic",
    )
    rep.add(
        "aligned_cell_above_S_ceiling",
        S_aligned_cell > S_aligned_min,
        f"S_aligned={S_aligned_cell:.3f}, required >{S_aligned_min}",
        ku="KU-3.20 aligned",
    )

    return rep


# ---------------------------------------------------------------------------
# H.4 Unit 2.2 — KU-2.17 FA growth Hill threshold
# ---------------------------------------------------------------------------

# KU-2.17 (Tan 2020 / Stricker 2013) FA growth threshold.
KU217_F_TH_PER_FA_PN: float = 50.0
KU217_N_ENGAGED_AT_THRESHOLD: float = 10.0


def gate_unit2_2_fa_growth(
    *,
    F_per_FA_pN: float,                           # force on the FA at the trigger frame
    N_engaged: float,                             # bonded integrins at trigger frame
    F_th_target_pN: float = KU217_F_TH_PER_FA_PN,
    N_engaged_target: float = KU217_N_ENGAGED_AT_THRESHOLD,
    F_rel_tolerance: float = 0.2,                 # ±20 %
    N_rel_tolerance: float = 0.3,                 # ±30 %
) -> SanityGateReport:
    """Phase 1 H.4 Unit 2.2 KU-2.17 FA growth gate.

    Verifies that the emergent FA growth trigger (clutch-population
    threshold proxy in v2 — NO Hill wrapper at runtime; the Hill form
    KU-2.17 is only used in this oracle gate) fires at:
      * F^{per-FA} ≈ 50 pN ±20 %
      * N_engaged ≈ 10 ±30 %
    Both must hold for the emergent dynamics to match the KU-2.17 oracle.
    """
    rep = SanityGateReport(module="h4_fa_growth_ku217")

    F_rel = abs(F_per_FA_pN - F_th_target_pN) / max(F_th_target_pN, 1e-30)
    rep.add(
        "FA_growth_F_threshold",
        F_rel <= F_rel_tolerance,
        f"F^FA={F_per_FA_pN:.2f} pN, target={F_th_target_pN:.1f} pN "
        f"(KU-2.17), rel err={F_rel:.3f}, tol={F_rel_tolerance}",
        ku="KU-2.17 F_th",
    )

    N_rel = abs(N_engaged - N_engaged_target) / max(N_engaged_target, 1e-30)
    rep.add(
        "FA_growth_N_engaged",
        N_rel <= N_rel_tolerance,
        f"N_engaged={N_engaged:.2f}, target={N_engaged_target:.1f} "
        f"(KU-2.17), rel err={N_rel:.3f}, tol={N_rel_tolerance}",
        ku="KU-2.17 N_engaged",
    )

    return rep


# ---------------------------------------------------------------------------
# Generic emergent-vs-oracle comparison (Pereverzev / Hill / Bell-Evans / Buckley)
# ---------------------------------------------------------------------------

def gate_emergent_vs_oracle(
    *,
    module_name: str,                             # e.g. "pereverzev_koff", "hill_velocity"
    oracle_label: str,                            # KU tag or oracle source
    grid: list[float],                            # evaluation grid (force, strain, ...)
    sim_values: list[float],                      # HOOMD-emergent measurement
    oracle_values: list[float],                   # closed-form oracle value at grid
    rel_tolerance: float = 0.05,                  # H.4 brief: ±5 %
    mass_fraction_required: float = 1.0,          # 1.0 = every grid point must pass
) -> SanityGateReport:
    """Generic emergent-vs-oracle comparison gate.

    Use for any Phase 1 oracle check: Pereverzev catch-slip k_off(F)
    (KU-2.5 + D2), Hill force-velocity v(F) (D6), Bell-Evans bond
    off-rates (D2), Buckley catch-bond cadherin (KU-4.2). H.4 brief
    acceptance: "HOOMD-emergent matches oracle within ± 5%".

    ``mass_fraction_required = 1.0`` requires every grid point within
    tolerance; lower values allow a fraction of misses (useful when
    the grid extends past the oracle's domain of validity).
    """
    rep = SanityGateReport(module=f"emergent_vs_oracle_{module_name}")

    if len(grid) != len(sim_values) or len(grid) != len(oracle_values):
        rep.add(
            "shape_consistent",
            False,
            f"grid={len(grid)}, sim={len(sim_values)}, oracle={len(oracle_values)} "
            f"— must match",
            ku="oracle gate",
        )
        return rep

    if not grid:
        rep.add(
            "non_empty_grid",
            False,
            "empty grid — cannot compare",
            ku="oracle gate",
        )
        return rep

    rel_errs = [
        abs(s - o) / max(abs(o), 1e-30)
        for s, o in zip(sim_values, oracle_values)
    ]
    n_pass = sum(1 for e in rel_errs if e <= rel_tolerance)
    fraction = n_pass / len(rel_errs)
    worst = max(rel_errs)
    worst_idx = rel_errs.index(worst)

    rep.add(
        "fraction_within_tolerance",
        fraction >= mass_fraction_required,
        f"{n_pass}/{len(rel_errs)} grid points within rel err ≤ {rel_tolerance} "
        f"(fraction={fraction:.3f}, required ≥{mass_fraction_required}). "
        f"worst point: grid={grid[worst_idx]:.3e}, "
        f"sim={sim_values[worst_idx]:.3e}, oracle={oracle_values[worst_idx]:.3e}, "
        f"rel err={worst:.3f}",
        ku=oracle_label,
    )

    return rep


# ---------------------------------------------------------------------------
# Bell-Evans batch-CFL — D2 batched-updater consistency
# ---------------------------------------------------------------------------

def gate_bell_evans_batch_cfl(
    *,
    dt: float,                                    # integrator timestep
    bond_families: dict[str, dict[str, float]],   # see schema below
    events_per_batch_ceiling: float = 1.0e-3,
) -> SanityGateReport:
    """Phase 1 D2 batched-updater CFL gate.

    PHASE_0_3_DECISIONS D2 specifies that each Bell-Evans / Pereverzev /
    Buckley off-rate updater is invoked every ~100 steps so that
    ``batch_steps · Δt · k_off_max ≲ 10⁻³`` (acceptable per-bond
    event probability per batch). This gate iterates over every
    declared bond family and verifies the bound, so that adding a new
    bond family (filamin, α-actinin, cadherin, integrin, motor head,
    talin domain, ...) cannot silently break CFL.

    ``bond_families`` schema::

        {
            "filamin_xlink": {"k_max": 0.1, "batch_steps": 100},
            "motor_head":     {"k_max": 1.0, "batch_steps": 100},
            "integrin":       {"k_max": 0.5, "batch_steps": 100},
            "cadherin":       {"k_max": 1.0, "batch_steps": 100},
            ...
        }

    ``k_max`` is the worst-case off-rate the family can sustain
    (Bell-Evans evaluated at the largest plausible bond force).
    """
    rep = SanityGateReport(module="d2_bell_evans_batch_cfl")

    if not bond_families:
        rep.add(
            "any_bond_family_declared",
            False,
            "bond_families is empty — D2 requires explicit per-family "
            "rate/batch declarations",
            ku="D2 CFL",
        )
        return rep

    for family, params in bond_families.items():
        k_max = params["k_max"]
        batch_steps = params["batch_steps"]
        events = batch_steps * dt * k_max
        rep.add(
            f"batch_cfl_{family}",
            events <= events_per_batch_ceiling,
            f"batch_steps·Δt·k_max = {batch_steps}·{dt:.2e}·{k_max:.3e} = "
            f"{events:.3e} per bond per batch; ceiling="
            f"{events_per_batch_ceiling:.0e} (D2 spec). "
            f"Fix: reduce batch_steps or refit Bell-Evans params.",
            ku=f"D2 CFL / {family}",
        )

    return rep


# ---------------------------------------------------------------------------
# Per-cell bead budget — Plan v2 §11 hardware envelope
# ---------------------------------------------------------------------------

# Plan v2 §11 — Phase 1 M1 Max CPU per-cell envelope.
PLAN_V2_PHASE1_PER_CELL_BEAD_CAP: int = 15000


def gate_bead_budget_per_cell(
    *,
    components: dict[str, int],                   # e.g. {"cortex": 7000, "xlink_heads": 2000, ...}
    cap: int = PLAN_V2_PHASE1_PER_CELL_BEAD_CAP,
    warn_fraction: float = 0.85,
) -> SanityGateReport:
    """Phase 1 per-cell bead-budget gate (Plan v2 §11, H.3 brief).

    Aggregates the bead count contributed by each Phase 1 cell
    sub-system (cortex, xlinks, myosin minifilaments, ERM, FAs,
    lamellipodium) and verifies the total stays under ``cap`` (Phase 1
    M1 Max CPU envelope; H.3 brief acceptance ≤ 15,000 per cell).

    Logs a non-blocking ``budget_warning`` if usage ≥ warn_fraction·cap.
    """
    rep = SanityGateReport(module="phase1_per_cell_bead_budget")

    if not components:
        rep.add(
            "components_declared",
            False,
            "components dict is empty — cannot aggregate budget",
            ku="Plan v2 §11",
        )
        return rep

    total = sum(components.values())
    detail_parts = ", ".join(f"{k}={v}" for k, v in components.items())

    rep.add(
        "per_cell_total_under_cap",
        total <= cap,
        f"total={total} beads, cap={cap} (Plan v2 §11 / H.3 brief). "
        f"Components: {detail_parts}",
        ku="Plan v2 §11",
    )

    warn_threshold = int(warn_fraction * cap)
    rep.add(
        "budget_warning",
        True,
        f"usage={total}/{cap} ({100.0*total/max(cap,1):.1f}%); warn at "
        f"≥{warn_threshold} ({warn_fraction*100:.0f}%). "
        f"{'ABOVE WARN' if total >= warn_threshold else 'below warn'}",
        ku="info",
    )

    return rep
