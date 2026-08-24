r"""No-drift bridge: the density force-budget ledger uses the SAME per-side force the NG-1 gate validated.

Round-2 closeout. The A5000 NG-1 two-filament isometric-stall gate PASSED 6/6 (transmission 0.61 → 1.0056,
per-head load → f_stall, reactions ≈ ±F_side) — the motor is faithful, so the "low per-minifilament force" is
REAL (not a motor bug). Any downstream DENSITY claim must therefore rest on the EXACT force the gate validated
(:mod:`aleph.components.motor.two_filament_reference`, the analytic reference the native gate is scored against). This
file (a) locks that per-side force to the density ledger (:mod:`aleph.components.motor.force_budget_ledger`), and
(b) pins the HONEST density bracket — the loose upper bound and the force-dipole lower bound differ ~100×, so
whether physiological density closes the active band is NOT settled by the ledger arithmetic (it needs the
assembled-cell method-of-planes with the fixed motor).

This gate locks that identity so a future edit to either the ledger's ``per_side_force_pn`` or the reference's
``per_side_stall_force`` (φ_b, f_stall, N_side, or the Bell constants behind them) cannot silently DIVERGE the
density conclusion from the motor the A5000 actually ran. If they drift, the density-floor claim is no longer
grounded in the validated motor — surface to PI (do NOT re-tune either to re-agree).
"""

from __future__ import annotations

import pytest

from aleph.components.motor import force_budget_ledger as ledger
from aleph.components.motor import two_filament_reference as ref

# claim_a (AFINES; ac.cell.assemble t0) — the NG-1 gate topology magnitudes
N_SIDE_A, F_STALL_A = 10, 0.5
F_SIDE_A_EXPECTED = 4.9627384744700995   # N_side·φ_b·f_stall with φ_b = k_on/(k_on+k_off0 e^{f/f0})


def test_ledger_per_side_equals_ng1_reference() -> None:
    """The density ledger's per-side force == the two-filament reference == the NG-1 isometric F_side."""
    lf = ledger.per_side_force_pn(N_SIDE_A, F_STALL_A)
    rf = ref.per_side_stall_force(N_SIDE_A, F_STALL_A)["f_side"]
    ri = ref.two_filament_isometric(N_SIDE_A, F_STALL_A)["f_side"]
    assert lf == pytest.approx(rf, rel=1e-12)
    assert lf == pytest.approx(ri, rel=1e-12)
    assert lf == pytest.approx(F_SIDE_A_EXPECTED, rel=1e-9)


def test_engaged_fraction_agrees() -> None:
    """Both code paths derive the SAME emergent Bell engaged fraction φ_b (not an imposed duty)."""
    phi_ledger = ledger.engaged_fraction(F_STALL_A)
    phi_ref = ref.per_side_stall_force(N_SIDE_A, F_STALL_A)["phi_b"]
    assert phi_ledger == pytest.approx(phi_ref, rel=1e-12)
    assert phi_ledger == pytest.approx(0.99254769489402, rel=1e-9)


def test_bell_constants_do_not_drift() -> None:
    """The Bell constants behind F_side (k_on, k_off0, f0) match between the ledger and the reference."""
    assert ledger.K_ON == pytest.approx(ref.K_ON)
    assert ledger.K_OFF0 == pytest.approx(ref.K_OFF0)
    assert ledger.F0_PN == pytest.approx(ref.F0_PN, rel=1e-9)


def test_density_band_crossing_is_bracketed_not_settled() -> None:
    r"""HONEST bracket: the two density estimators disagree ~100×, so the band-crossing is NOT settled.

    The Round-2 closeout must NOT claim "physiological density closes the band" — that only holds for the loose
    UPPER bound. This gate pins the honest picture so a future edit cannot quietly re-assert the optimistic
    single-estimator claim:

      * loose bound (``force_budget_ledger``, all minifilaments on one cut) — an explicit UPPER bound — sits
        ABOVE the active band at physiological density; while
      * the force-dipole / method-of-planes lower bound (``two_filament_reference.gamma_ceiling`` = ρ·F·L_bb/2,
        the convention the NG-1 gate's reference uses) is still WELL BELOW the band at the same density;
      * the two differ by ~2 orders of magnitude (the isotropic-Kirkwood / network-propagation gap).

    The real γ lies between them and requires the assembled-cell method-of-planes with the FIXED motor —
    ``force_budget_ledger``'s loose-bound "above band" is NOT the physical answer.
    """
    from aleph.laws.gamma_estimator import active_band_pn_um
    band_lo, band_hi = active_band_pn_um()
    rho = 21.0                                                       # upper physiological (gamma-floor 2026-06-09)
    loose = ledger.gamma_at_density(rho, N_SIDE_A, F_STALL_A)[1]     # pN/µm — upper bound
    dipole = ref.gamma_ceiling(F_SIDE_A_EXPECTED, density_um2=rho)["gamma_pn_um"]   # pN/µm — physical lower bound
    assert loose > band_hi                     # loose UPPER bound clears the band ...
    assert dipole < band_lo                    # ... but the physical dipole lower bound does NOT
    assert loose / dipole > 50.0               # ~100× apart ⇒ band-crossing is unresolved, not "legitimate"


def test_dipole_needs_superphysiological_density_for_claim_a() -> None:
    """On the physical (force-dipole) estimator, claim_a needs ρ ≫ physiological to reach the band midpoint.

    Guards against the loose-bound optimism: the sourced physiological density band (16-21 /µm²) does NOT reach
    the required density on the physical estimator, so density alone is not shown to close the gap for claim_a.
    """
    gc = ref.gamma_ceiling(F_SIDE_A_EXPECTED)
    rho_need = gc["rho_needed_band_mid_um2"]
    assert rho_need > 21.0                     # far above the upper physiological 21 /µm² (it is ~469)
