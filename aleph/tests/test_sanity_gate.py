"""Sanity gate unit tests (added 2026-05-20, pre-Phase-1-dispatch).

Each new gate is exercised with one PASS configuration and one FAIL
configuration so that future edits cannot silently relax the contract.
"""

from __future__ import annotations

import math

import pytest

from aleph.validation.oracles.common.sanity_gate import (
    D5_BACKBONE_BEADS,
    D5_BACKBONE_LENGTH_NM,
    D5_HEADS_PER_SIDE,
    D5_HEAD_REST_LENGTH_NM,
    D5_HEAD_SPRING_K_PN_PER_UM,
    D7_EPSILON_OVER_KT,
    KU11_LP_TARGET,
    KU217_F_TH_PER_FA_PN,
    KU217_N_ENGAGED_AT_THRESHOLD,
    KU35_TENSION_TARGET_N_PER_M,
    PHASE_1_CANONICAL_INTEGRATORS,
    PHASE_1_REFERENCE_INTEGRATORS,
    PLAN_V2_PHASE1_PER_CELL_BEAD_CAP,
    WCA_CUTOFF_FACTOR,
    gate_bead_budget_per_cell,
    gate_bell_evans_batch_cfl,
    gate_blebbistatin_response,
    gate_cortex_tension,
    gate_emergent_vs_oracle,
    gate_h2_persistence_length,
    gate_minifilament_topology,
    gate_nematic_order,
    gate_unit1_2_dynamics,
    gate_unit2_2_fa_growth,
    gate_wca_cutoff,
)


# ---------- D3 integrator name --------------------------------------------

_DYNAMICS_DEFAULTS = dict(
    dt=1.0e-5,
    tau_xl=1.0e-3,
    tau_stretch=1.0e-3,
    tau_bend=1.0e-3,
    safety_factor=0.1,
)


def test_dynamics_gate_accepts_lm_baoab():
    rep = gate_unit1_2_dynamics(
        integrator_name="leimkuhler_matthews_baoab",
        **_DYNAMICS_DEFAULTS,
    )
    assert rep.passed, rep.summary()


def test_dynamics_gate_rejects_em_when_not_reference_mode():
    rep = gate_unit1_2_dynamics(
        integrator_name="euler_maruyama",
        **_DYNAMICS_DEFAULTS,
    )
    assert not rep.passed
    # CFL passes; integrator check fails.
    failing = [c for c in rep.checks if c.status == "FAIL"]
    assert any("integrator" in c.name for c in failing)


def test_dynamics_gate_allows_em_in_reference_mode():
    rep = gate_unit1_2_dynamics(
        integrator_name="euler_maruyama",
        allow_reference_integrator=True,
        **_DYNAMICS_DEFAULTS,
    )
    assert rep.passed, rep.summary()


def test_canonical_and_reference_sets_disjoint():
    assert PHASE_1_CANONICAL_INTEGRATORS.isdisjoint(PHASE_1_REFERENCE_INTEGRATORS)


# ---------- H.2 persistence length ----------------------------------------


def _h2_pass_kwargs(**overrides):
    base = dict(
        Lp_measured=17.0e-6,
        ks_pvalue=0.4,
        equipartition_relerr_per_bond=[0.01, 0.02, 0.015, 0.018, 0.04],
    )
    base.update(overrides)
    return base


def test_h2_pass():
    rep = gate_h2_persistence_length(**_h2_pass_kwargs())
    assert rep.passed, rep.summary()


def test_h2_fail_lp_out_of_band():
    rep = gate_h2_persistence_length(**_h2_pass_kwargs(Lp_measured=12.0e-6))
    assert not rep.passed


def test_h2_fail_ks():
    rep = gate_h2_persistence_length(**_h2_pass_kwargs(ks_pvalue=0.01))
    assert not rep.passed


def test_h2_fail_equipartition():
    rep = gate_h2_persistence_length(
        **_h2_pass_kwargs(equipartition_relerr_per_bond=[0.01, 0.4])
    )
    assert not rep.passed


def test_h2_lm_correctness_pass():
    rep = gate_h2_persistence_length(
        **_h2_pass_kwargs(
            Lp_measured=17.0e-6,
            Lp_em_reference=17.05e-6,
            Lp_em_sigma=0.2e-6,
        )
    )
    assert rep.passed, rep.summary()


def test_h2_lm_correctness_fail_far_from_em():
    rep = gate_h2_persistence_length(
        **_h2_pass_kwargs(
            Lp_measured=17.0e-6,
            Lp_em_reference=16.0e-6,
            Lp_em_sigma=0.2e-6,
        )
    )
    assert not rep.passed


def test_h2_lm_correctness_requires_positive_sigma():
    rep = gate_h2_persistence_length(
        **_h2_pass_kwargs(
            Lp_em_reference=17.05e-6,
            Lp_em_sigma=0.0,
        )
    )
    assert not rep.passed


def test_h2_constants():
    assert KU11_LP_TARGET == 17.0e-6


# ---------- D7 WCA cutoff -------------------------------------------------


def _wca_pair(sigma: float, kT: float) -> dict[str, float]:
    return dict(
        sigma=sigma,
        r_cut=WCA_CUTOFF_FACTOR * sigma,
        epsilon=D7_EPSILON_OVER_KT * kT,
    )


def test_wca_pass():
    kT = 4.28e-21
    sigma = 60.0e-9
    rep = gate_wca_cutoff(
        pair_params={("actin", "actin"): _wca_pair(sigma, kT)},
        kT=kT,
    )
    assert rep.passed, rep.summary()


def test_wca_fail_cutoff_too_long():
    kT = 4.28e-21
    sigma = 60.0e-9
    bad = _wca_pair(sigma, kT)
    bad["r_cut"] = 2.5 * sigma  # vanilla LJ-style cutoff: attractive tail leaks in
    rep = gate_wca_cutoff(pair_params={("a", "a"): bad}, kT=kT)
    assert not rep.passed


def test_wca_fail_wrong_epsilon():
    kT = 4.28e-21
    sigma = 60.0e-9
    bad = _wca_pair(sigma, kT)
    bad["epsilon"] = 5.0 * kT
    rep = gate_wca_cutoff(pair_params={("a", "a"): bad}, kT=kT)
    assert not rep.passed


def test_wca_required_pair_missing():
    kT = 4.28e-21
    rep = gate_wca_cutoff(
        pair_params={("actin", "actin"): _wca_pair(60e-9, kT)},
        kT=kT,
        required_pair_types=[("actin", "actin"), ("actin", "myosin_head")],
    )
    assert not rep.passed


def test_wca_empty_pair_params_fails():
    rep = gate_wca_cutoff(pair_params={}, kT=4.28e-21)
    assert not rep.passed


# ---------- D5 minifilament topology --------------------------------------


def _minifilament_pass_kwargs(**overrides):
    base = dict(
        n_backbone_beads=D5_BACKBONE_BEADS,
        n_heads_per_side=D5_HEADS_PER_SIDE,
        n_sides=2,
        has_rigid_body_constraint=True,
        cross_bridge_k_pN_per_um=D5_HEAD_SPRING_K_PN_PER_UM,
        cross_bridge_r0_nm=D5_HEAD_REST_LENGTH_NM,
        backbone_length_nm=D5_BACKBONE_LENGTH_NM,
    )
    base.update(overrides)
    return base


def test_minifilament_pass():
    rep = gate_minifilament_topology(**_minifilament_pass_kwargs())
    assert rep.passed, rep.summary()


def test_minifilament_single_2head_regression_fails():
    """The exact AFINES single-2-head spring regression PI is guarding against."""
    rep = gate_minifilament_topology(
        **_minifilament_pass_kwargs(
            n_backbone_beads=1,
            n_heads_per_side=1,
            n_sides=2,
            has_rigid_body_constraint=False,
        )
    )
    assert not rep.passed
    fails = {c.name for c in rep.checks if c.status == "FAIL"}
    assert "backbone_bead_count" in fails
    assert "heads_per_side_count" in fails
    assert "rigid_body_constraint_present" in fails


def test_minifilament_unipolar_fails_bipolar_check():
    rep = gate_minifilament_topology(**_minifilament_pass_kwargs(n_sides=1))
    assert not rep.passed
    assert any(
        c.name == "bipolar_two_sides" and c.status == "FAIL" for c in rep.checks
    )


def test_minifilament_wrong_spring_stiffness_fails():
    rep = gate_minifilament_topology(
        **_minifilament_pass_kwargs(cross_bridge_k_pN_per_um=10.0)
    )
    assert not rep.passed


# ---------- Cortex tension / blebbistatin / nematic -----------------------


def test_cortex_tension_pass():
    rep = gate_cortex_tension(
        tension_measured_N_per_m=KU35_TENSION_TARGET_N_PER_M * 1.1
    )
    assert rep.passed, rep.summary()


def test_cortex_tension_fail():
    rep = gate_cortex_tension(
        tension_measured_N_per_m=KU35_TENSION_TARGET_N_PER_M * 2.0
    )
    assert not rep.passed


def test_blebbistatin_response_pass():
    rep = gate_blebbistatin_response(
        tension_myosin_on_N_per_m=0.5e-3,
        tension_myosin_off_N_per_m=0.1e-3,
        final_aspect_ratio_myosin_on=1.1,
        final_aspect_ratio_myosin_off=1.5,
    )
    assert rep.passed, rep.summary()


def test_blebbistatin_response_fails_when_off_still_rounds():
    rep = gate_blebbistatin_response(
        tension_myosin_on_N_per_m=0.5e-3,
        tension_myosin_off_N_per_m=0.1e-3,
        final_aspect_ratio_myosin_on=1.1,
        final_aspect_ratio_myosin_off=1.15,
    )
    assert not rep.passed


def test_blebbistatin_response_fails_when_tension_does_not_drop():
    rep = gate_blebbistatin_response(
        tension_myosin_on_N_per_m=0.5e-3,
        tension_myosin_off_N_per_m=0.4e-3,
        final_aspect_ratio_myosin_on=1.1,
        final_aspect_ratio_myosin_off=1.5,
    )
    assert not rep.passed


def test_nematic_order_pass():
    rep = gate_nematic_order(S_isotropic_cell=0.05, S_aligned_cell=0.4)
    assert rep.passed, rep.summary()


def test_nematic_order_fail_isotropic_too_ordered():
    rep = gate_nematic_order(S_isotropic_cell=0.2, S_aligned_cell=0.4)
    assert not rep.passed


def test_nematic_order_fail_aligned_too_disordered():
    rep = gate_nematic_order(S_isotropic_cell=0.05, S_aligned_cell=0.2)
    assert not rep.passed


# ---------- KU-2.17 FA growth ---------------------------------------------


def test_fa_growth_pass():
    rep = gate_unit2_2_fa_growth(
        F_per_FA_pN=KU217_F_TH_PER_FA_PN,
        N_engaged=KU217_N_ENGAGED_AT_THRESHOLD,
    )
    assert rep.passed, rep.summary()


def test_fa_growth_fail_force():
    rep = gate_unit2_2_fa_growth(
        F_per_FA_pN=5.0,
        N_engaged=KU217_N_ENGAGED_AT_THRESHOLD,
    )
    assert not rep.passed


def test_fa_growth_fail_N_engaged():
    rep = gate_unit2_2_fa_growth(
        F_per_FA_pN=KU217_F_TH_PER_FA_PN,
        N_engaged=2.0,
    )
    assert not rep.passed


# ---------- Emergent vs oracle --------------------------------------------


def test_emergent_vs_oracle_pass_hill_like():
    grid = [0.0, 1.0, 2.0, 3.0]
    oracle = [10.0, 6.0, 3.0, 1.0]
    sim = [10.05, 5.95, 3.06, 0.98]
    rep = gate_emergent_vs_oracle(
        module_name="hill_velocity",
        oracle_label="D6 Hill",
        grid=grid,
        sim_values=sim,
        oracle_values=oracle,
    )
    assert rep.passed, rep.summary()


def test_emergent_vs_oracle_fail():
    rep = gate_emergent_vs_oracle(
        module_name="hill_velocity",
        oracle_label="D6 Hill",
        grid=[0.0, 1.0],
        sim_values=[10.0, 100.0],
        oracle_values=[10.0, 6.0],
    )
    assert not rep.passed


def test_emergent_vs_oracle_shape_mismatch_fails():
    rep = gate_emergent_vs_oracle(
        module_name="x",
        oracle_label="x",
        grid=[0.0, 1.0],
        sim_values=[10.0],
        oracle_values=[10.0, 6.0],
    )
    assert not rep.passed


def test_emergent_vs_oracle_empty_grid_fails():
    rep = gate_emergent_vs_oracle(
        module_name="x",
        oracle_label="x",
        grid=[],
        sim_values=[],
        oracle_values=[],
    )
    assert not rep.passed


# ---------- D2 Bell-Evans batch CFL ---------------------------------------


def test_batch_cfl_pass():
    rep = gate_bell_evans_batch_cfl(
        dt=2.0e-5,
        bond_families={
            "filamin_xlink": {"k_max": 0.1, "batch_steps": 100},
            "motor_head":    {"k_max": 1.0, "batch_steps": 100},
            "integrin":      {"k_max": 0.5, "batch_steps": 100},
        },
    )
    # 100 · 2e-5 · 1.0 = 2e-3 — fails because motor_head trips ceiling
    # Use a stricter case below; here verify we caught it.
    assert not rep.passed


def test_batch_cfl_pass_safe_params():
    rep = gate_bell_evans_batch_cfl(
        dt=2.0e-5,
        bond_families={
            "filamin_xlink": {"k_max": 0.1, "batch_steps": 100},   # 2e-4 ✓
            "alpha_actinin": {"k_max": 0.1, "batch_steps": 100},
        },
    )
    assert rep.passed, rep.summary()


def test_batch_cfl_fail_high_rate():
    rep = gate_bell_evans_batch_cfl(
        dt=2.0e-5,
        bond_families={
            "cadherin": {"k_max": 10.0, "batch_steps": 100},       # 2e-2 ✗
        },
    )
    assert not rep.passed


def test_batch_cfl_empty_fails():
    rep = gate_bell_evans_batch_cfl(dt=2.0e-5, bond_families={})
    assert not rep.passed


# ---------- Bead budget ---------------------------------------------------


def test_bead_budget_pass():
    rep = gate_bead_budget_per_cell(
        components={
            "cortex_actin": 7000,
            "xlink_heads":  2000,
            "myosin":       3400,
            "erm":           500,
        },
    )
    assert rep.passed, rep.summary()


def test_bead_budget_fail_over_cap():
    rep = gate_bead_budget_per_cell(
        components={
            "cortex_actin": 7000,
            "xlink_heads":  2000,
            "myosin":       3400,
            "erm":           500,
            "lamellipodium": 5000,  # pushes over 15000
        },
    )
    assert not rep.passed


def test_bead_budget_empty_fails():
    rep = gate_bead_budget_per_cell(components={})
    assert not rep.passed


def test_bead_budget_cap_constant():
    assert PLAN_V2_PHASE1_PER_CELL_BEAD_CAP == 15000


# ---------- WCA cutoff factor is the LJ minimum ---------------------------


def test_wca_cutoff_factor_is_lj_minimum():
    """WCA cutoff is exactly the LJ minimum (force-zero point)."""
    assert math.isclose(WCA_CUTOFF_FACTOR, 2.0 ** (1.0 / 6.0))


if __name__ == "__main__":  # pragma: no cover
    pytest.main([__file__, "-q"])
