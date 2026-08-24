"""Host gate for the two-filament sarcomere DECISION reference (pure NumPy — no Warp/CUDA).

Certifies the analytic reference (``ac.motor.two_filament_reference``) that the A5000 Warp runner
(``ac.motor.native_gates.ng1_two_filament_stall``) is gated against:

  * it REPRODUCES the ensemble-stall oracle (``ensemble_stall_analytic``) — no reinvented stall math;
  * ``F_side ≈ 4.96 pN`` at claim_a (N_side=10, f_stall=0.5, the AFINES values in ac.cell.assemble);
  * ``k_eff ≈ 90.9 pN/µm`` (the per-side crossbridge+arm series compliance);
  * the ideal collinear transmission ratio = 1 (PROVEN by a 1-D series-spring equilibrium, not asserted);
  * the two clamp reactions are equal-and-opposite ±F_side;
  * the γ ceiling at the sourced Nie-2015 density is FLOORED far below the MCF7 active band (the density
    finding — the reason a faithful low per-side force does NOT need more heads);
  * the reference's sourced constants do NOT drift from ``ac.cell.assemble`` / ``ff.gamma_floor``.

It also prints the (N_side, f_stall) sweep table + k_eff (run ``pytest -s`` to see it).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor import two_filament_reference as ref
from aleph.components.motor.bell_kinetics_analytic import bell_f0_from_x_beta, engaged_fraction_steady
from aleph.components.motor.ensemble_stall_analytic import ensemble_stall_meanfield


# ── reproduces the reused oracle (no reinvented Hill/Bell/stall math) ─────────────────────────────────────
def test_reference_reproduces_ensemble_stall_oracle() -> None:
    """per_side_stall_force's F_side IS ensemble_stall_meanfield's f_ensemble; φ_b matches the Bell oracle."""
    for n_side in (10, 28):
        for f_stall in (0.5, 2.0):
            r = ref.per_side_stall_force(n_side, f_stall)
            oracle = ensemble_stall_meanfield(n_side, f_stall, ref.K_ON, ref.K_OFF0, ref.F0_PN)
            assert r["f_side"] == pytest.approx(oracle["f_ensemble"], rel=1e-12)
            assert r["oracle_f_ensemble"] == pytest.approx(oracle["f_ensemble"], rel=1e-12)
            # φ_b is the EMERGENT Bell engaged fraction, cross-checked independently
            phi = float(engaged_fraction_steady(ref.K_ON, f_stall, ref.K_OFF0, ref.F0_PN))
            assert r["phi_b"] == pytest.approx(phi, rel=1e-12)
            assert r["phi_b_crosscheck"] == pytest.approx(phi, rel=1e-12)
            # emergent, NOT the naive product (P2): strictly below N_side·f_stall since φ_b < 1
            assert r["f_side"] < r["f_naive"]
            assert r["f_side"] == pytest.approx(n_side * phi * f_stall, rel=1e-12)


def test_f0_matches_bell_x_beta() -> None:
    """The reference f0 equals kBT/x_beta (bell_kinetics_analytic), i.e. assemble.NMII_F0 ≈ 7.13 pN."""
    assert ref.F0_PN == pytest.approx(bell_f0_from_x_beta(0.6e-3), rel=1e-12)
    assert ref.F0_PN == pytest.approx(7.1333333333, rel=1e-6)


# ── the headline claim_a number ──────────────────────────────────────────────────────────────────────────
def test_f_side_claim_a_is_4p96() -> None:
    """F_side ≈ 4.96 pN at claim_a (N_side=10, f_stall=0.5) — the emergent per-side isometric reaction."""
    r = ref.per_side_stall_force(ref.N_SIDE_CLAIM_A, ref.F_STALL_CLAIM_A)
    assert r["f_side"] == pytest.approx(4.9627, abs=0.01)
    assert r["phi_b"] == pytest.approx(0.99255, abs=1e-4)
    assert r["f_naive"] == pytest.approx(5.0, rel=1e-12)


# ── k_eff (series compliance) ────────────────────────────────────────────────────────────────────────────
def test_k_eff_series_stiffness() -> None:
    """k_eff = (1/k_xb + 1/k_head_arm)^-1 ≈ 90.9 pN/µm (1000 & 100 in series)."""
    assert ref.series_stiffness() == pytest.approx(90.90909090909, rel=1e-9)
    assert ref.series_stiffness(1000.0, 100.0) == pytest.approx(1.0 / (1.0 / 1000.0 + 1.0 / 100.0), rel=1e-12)
    # full clamp-to-clamp path: crossbridge·arm·(N_BB−1 backbone bonds)·arm·crossbridge — the backbone is a
    # bond CHAIN (per-bond k), so its end-to-end term is (n_bb−1)/k_backbone ⇒ soft ≈ 28.6 pN/µm.
    assert ref.full_transmission_path_stiffness() == pytest.approx(
        1.0 / (2 / 1000 + 2 / 100 + (ref.N_BB - 1) / 1000), rel=1e-12)
    assert ref.full_transmission_path_stiffness() < ref.series_stiffness()


# ── ideal transmission = 1 (proven by a static series-spring equilibrium) ─────────────────────────────────
def test_ideal_transmission_is_one() -> None:
    """A lossless collinear static chain transmits ratio = 1 (real loss is geometric, exposed by the native gate)."""
    assert ref.ideal_transmission_ratio() == 1.0
    proof = ref.collinear_transmission_proof(f_generated=ref.per_side_stall_force()["f_side"])
    assert proof["ratio"] == pytest.approx(1.0, abs=1e-9)          # uniform tension ⇒ far == near
    assert proof["tension_spread_pN"] < 1e-9                        # every series bond carries the same tension
    assert proof["k_series_pN_per_um"] == pytest.approx(ref.full_transmission_path_stiffness(), rel=1e-12)
    # the common tension equals the generated force (the reaction the far clamp receives)
    f_side = ref.per_side_stall_force()["f_side"]
    assert proof["tension_uniform_pN"] == pytest.approx(f_side, rel=1e-9)


def test_transmission_ratio_independent_of_active_spring_and_force() -> None:
    """Ratio = 1 regardless of which spring is active or how large the generated force is (pure series)."""
    for active in range(5):
        for fg in (0.5, 4.96, 55.0):
            proof = ref.collinear_transmission_proof(f_generated=fg, active_index=active)
            assert proof["ratio"] == pytest.approx(1.0, abs=1e-9)


# ── two-filament isometric bundle: equal-and-opposite clamp reactions ─────────────────────────────────────
def test_isometric_clamp_reactions_equal_and_opposite() -> None:
    """The two clamps hold ±F_side (the minifilament pulls both actins inward with the same magnitude)."""
    b = ref.two_filament_isometric(ref.N_SIDE_CLAIM_A, ref.F_STALL_CLAIM_A)
    f_side = b["f_side"]
    assert b["clamp_reaction_A_pN"] == pytest.approx(+f_side, rel=1e-12)
    assert b["clamp_reaction_B_pN"] == pytest.approx(-f_side, rel=1e-12)
    assert b["clamp_reaction_A_pN"] + b["clamp_reaction_B_pN"] == pytest.approx(0.0, abs=1e-12)
    assert b["ideal_transmission"] == 1.0
    assert b["k_eff_pN_per_um"] == pytest.approx(90.909, rel=1e-4)
    # the compliance the native gate must probe: chain stretch = F_side/k_eff
    assert b["chain_stretch_at_load_nm"] == pytest.approx(f_side / ref.series_stiffness() * 1e3, rel=1e-9)


# ── γ ceiling: the density floor at the sourced Nie-2015 density ──────────────────────────────────────────
def test_gamma_ceiling_floored_below_active_band() -> None:
    """Even the most generous GAP claim leaves the R=7.5µm γ ceiling far below the MCF7 active band."""
    band_lo, _ = ref.MCF7_ACTIVE_BAND_PN_UM
    rows = ref.sweep_gap_claims()
    for row in rows:
        # dipole areal stress γ = ρ · F_side · L_bb / 2 (ff.network_contractility convention)
        expect = ref.NIE2015_DENSITY_UM2 * row["f_side_pN"] * ref.L_BB_UM / 2.0
        assert row["gamma_pn_um"] == pytest.approx(expect, rel=1e-12)
        assert row["gamma_pn_um"] < band_lo                       # floored under the band (density finding)
        assert row["floor_factor_lo"] > 1.0
    # even the max claim (N_side=28, f_stall=2.0) is many-fold under the band
    top = max(rows, key=lambda r: r["gamma_pn_um"])
    assert top["n_side"] == 28 and top["f_stall_pN"] == 2.0
    assert top["floor_factor_lo"] > 10.0                          # ≫ 10× under even at the structural upper end


def test_gamma_ceiling_matches_gamma_floor_finding() -> None:
    """The Nie density (0.625) is ~1-2 orders below the ρ needed for the band — the FF γ-floor finding."""
    r = ref.gamma_ceiling(ref.per_side_stall_force()["f_side"])
    assert r["rho_needed_band_mid_um2"] > 20.0 * ref.NIE2015_DENSITY_UM2   # need ≫ the measured density
    assert r["gamma_mN_per_m"] == pytest.approx(r["gamma_pn_um"] * 1e-3, rel=1e-12)


# ── sweep monotonicity + independent recomputation ───────────────────────────────────────────────────────
def test_sweep_monotone_and_independent() -> None:
    """F_side and γ rise with both N_side and f_stall; each row recomputes independently."""
    rows = {(r["n_side"], r["f_stall_pN"]): r for r in ref.sweep_gap_claims()}
    assert set(rows) == {(10, 0.5), (10, 2.0), (28, 0.5), (28, 2.0)}
    assert rows[(10, 0.5)]["f_side_pN"] < rows[(10, 2.0)]["f_side_pN"] < rows[(28, 2.0)]["f_side_pN"]
    assert rows[(10, 0.5)]["f_side_pN"] < rows[(28, 0.5)]["f_side_pN"] < rows[(28, 2.0)]["f_side_pN"]
    for (n_side, f_stall), r in rows.items():
        phi = float(engaged_fraction_steady(ref.K_ON, f_stall, ref.K_OFF0, ref.F0_PN))
        assert r["f_side_pN"] == pytest.approx(n_side * phi * f_stall, rel=1e-12)


# ── DRIFT GUARD: the reference's sourced constants must equal ac.cell.assemble / ff.gamma_floor ───────────
def test_no_drift_from_assemble_constants() -> None:
    """The reference mirrors ac.cell.assemble NMII_* exactly (no silent divergence of sourced values)."""
    from aleph.components.incumbent import assemble as A

    assert ref.N_SIDE_CLAIM_A == A.NMII_N_SIDE
    assert ref.F_STALL_CLAIM_A == A.NMII_F_STALL_TEST
    assert ref.K_XB == A.NMII_K_XB_TEST
    assert ref.K_HEAD_ARM == A.NMII_K_HEAD_ARM_TEST
    assert ref.K_BACKBONE == A.NMII_K_BACKBONE_TEST
    assert ref.K_ON == A.NMII_KON_TEST
    assert ref.K_OFF0 == A.NMII_KOFF0
    assert ref.F0_PN == pytest.approx(A.NMII_F0, rel=1e-12)
    assert ref.L_BB_UM == A.NMII_L_BB_UM
    assert ref.R_CELL_UM == A.R_CELL_UM


def test_no_drift_from_gamma_floor_density() -> None:
    """The minifilament areal density mirrors ff.gamma_floor's Nie-2015 datum (the only direct measurement)."""
    from aleph.laws import gamma_floor as G

    assert ref.NIE2015_DENSITY_UM2 == G.NIE2015_DENSITY_UM2
    assert ref.NIE2015_DENSITY_RANGE_UM2 == G.NIE2015_DENSITY_RANGE
    # the MCF7 active band mirrors ff.gamma_estimator.active_band_pn_um(SALBREUX_BAND_PN_UM)
    from aleph.laws.gamma_estimator import SALBREUX_BAND_PN_UM, active_band_pn_um

    band = active_band_pn_um(SALBREUX_BAND_PN_UM)
    assert ref.MCF7_ACTIVE_BAND_PN_UM[0] == pytest.approx(band[0], rel=1e-9)
    assert ref.MCF7_ACTIVE_BAND_PN_UM[1] == pytest.approx(band[1], rel=1e-9)


# ── print the decision table (visible with `pytest -s`) ──────────────────────────────────────────────────
def test_print_sweep_table_and_k_eff(capsys: pytest.CaptureFixture) -> None:
    """Emit the (N_side, f_stall) sweep table + k_eff + ideal transmission for the record."""
    table = ref.format_sweep_table()
    with capsys.disabled():
        print("\n" + table)
    assert "k_eff" in table and "gamma" in table
    assert "90.909" in table                     # the k_eff headline is in the rendered table
    # claim_a row is present with F_side ≈ 4.96
    rows = ref.sweep_gap_claims()
    claim_a = next(r for r in rows if r["n_side"] == 10 and r["f_stall_pN"] == 0.5)
    assert claim_a["f_side_pN"] == pytest.approx(4.9627, abs=0.01)
