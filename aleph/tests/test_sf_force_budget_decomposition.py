"""B5(ii) — SF generation-limit force-budget DECOMPOSITION consistency.

Locks the labeled gap decomposition (`h7_basal_sf_force_budget.decompose_generation_gap`)
so the "generation-bound, N_filaments-INDEPENDENT" conclusion stays internally consistent
and the cross-line link to the cortical γ-floor (same ½·n·f·ℓ budget) is exact.

These are ARITHMETIC consistency gates on the synthesis, not new physics: the raw budget
is the brief-literal minifilament; the decomposition asserts that gap = (per-minifilament
fidelity) × (cross-sectional NMII count) at the upper engagement, and that the literature
reference minifilament matches the cortical γ-floor's 56 pN dipole.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


# ── Hidden by a collection gap until 2026-08-09 ───────────────────────────────────────────────
# This module sat in `aleph/tests/` while `testpaths` listed only its sibling directories, so it was
# **never collected**. The merge widened `testpaths` to the whole tree (+410 tests) and these came
# back red. They are pre-existing failures, not regressions, and they are marked rather than fixed:
# the charter forbids editing an assertion to make a gate pass. `strict=True` means the mark itself
# is a ratchet — if one starts passing, the suite goes red and somebody has to say why.
#
# `scripts/h7_basal_sf_force_budget.py` was deleted by `00858821` ("finish the HOOMD removal — 167
# more files") on 2026-07-29. This module has been loading it by path ever since. Restoring the
# script is a physics decision — its numbers are h7-era — so it is reported, not resurrected.
pytestmark = pytest.mark.xfail(
    strict=True,
    reason="scripts/h7_basal_sf_force_budget.py deleted by 00858821 (2026-07-29); PI triage",
)

_FB_PATH = Path(__file__).resolve().parents[1] / "scripts" / "h7_basal_sf_force_budget.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("h7_basal_sf_force_budget", _FB_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def fb():
    return _load_module()


def _synthetic_budget(fb):
    """Brief-literal budget at the values the production config carries."""
    phi = 0.993                      # k_on/(k_on+k_off0) = 50/(50+0.35)
    H = 10                           # n_heads_per_side brief literal
    F_stall = 0.5e-12                # 0.5 pN brief literal
    f_mini = phi * H * F_stall
    return {
        "phi_engaged": phi,
        "n_heads_per_side": H,
        "F_stall_per_head_N": F_stall,
        "f_per_minifilament_N": f_mini,
        "T_cable_raw_N": f_mini,
        "kumar_band_N": list(fb.KUMAR_BAND_N),
        "in_kumar_band_raw": False,
        "gap_to_band_lo": fb.KUMAR_BAND_N[0] / f_mini,
        "mesoscale_force_factor_needed_to_band_centre": 0.5 * sum(fb.KUMAR_BAND_N) / f_mini,
        "n_motors": 40,
    }


def test_lit_minifilament_matches_cortical_gamma_floor(fb):
    """The reference dipole is the cortical γ-floor's 56 pN (28 heads × 2 pN)."""
    assert fb.LIT_HEADS_PER_SIDE == 28
    assert fb.LIT_F_PER_HEAD_N == pytest.approx(2.0e-12)
    assert fb.LIT_F_MINIFILAMENT_N == pytest.approx(56.0e-12, rel=1e-9)


def test_decomposition_product_recovers_raw_gap(fb):
    """gap = per-minifilament fidelity × cross-sectional count (to band centre)."""
    budget = _synthetic_budget(fb)
    d = fb.decompose_generation_gap(budget)
    raw_factor_to_centre = budget["mesoscale_force_factor_needed_to_band_centre"]
    # product check field must equal the raw factor-to-centre
    assert d["decomposition_check_product"] == pytest.approx(raw_factor_to_centre, rel=1e-9)
    assert (
        d["per_minifilament_fidelity_factor"]
        * d["cross_sectional_nmii_count_to_band_centre"]
        == pytest.approx(raw_factor_to_centre, rel=1e-9)
    )


def test_per_minifilament_fidelity_is_a_parameter_factor(fb):
    """Brief 5 pN → lit 56 pN ≈ 11× = (28/10 heads) × (2/0.5 pN/head) = 2.8× × 4×.

    Both are PARAMETER choices (the budget uses F_stall directly, already assuming the §9
    delivery fix) — same order as, but a distinct mechanism from, the §9 per-head recovery.
    """
    budget = _synthetic_budget(fb)
    d = fb.decompose_generation_gap(budget)
    assert d["per_minifilament_fidelity_factor"] == pytest.approx(11.2, abs=0.3)
    # the factor is exactly (lit heads × lit per-head) / (brief heads × brief per-head)
    heads_ratio = fb.LIT_HEADS_PER_SIDE / budget["n_heads_per_side"]
    perhead_ratio = fb.LIT_F_PER_HEAD_N / budget["F_stall_per_head_N"]
    assert d["per_minifilament_fidelity_factor"] == pytest.approx(heads_ratio * perhead_ratio, rel=1e-9)


def test_cross_sectional_count_is_order_hundreds(fb):
    """A Kumar SF needs ~180 (floor) to ~360 (centre) lit-minifilaments / cross-section."""
    budget = _synthetic_budget(fb)
    d = fb.decompose_generation_gap(budget)
    assert 150.0 <= d["cross_sectional_nmii_count_to_band_lo"] <= 220.0
    assert 300.0 <= d["cross_sectional_nmii_count_to_band_centre"] <= 420.0


def test_engagement_realism_is_a_penalty_not_a_credit(fb):
    """Realistic ~11% engagement makes the bound ~9× WORSE, not better."""
    budget = _synthetic_budget(fb)
    d = fb.decompose_generation_gap(budget)
    assert d["engagement_realism_penalty_factor"] > 1.0
    assert d["engagement_realism_penalty_factor"] == pytest.approx(8.7, abs=1.0)
