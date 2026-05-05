"""Tests for V2 Adhesion-Network Layer (cadherin/integrin program dynamics).

16 tests per locked §4 of
``docs/v2_adhesion_network_layer_locked.md`` (commit ``7acf0fd``).

PI ``id=2003`` directive forward guards (Y15 + Y16):
- coexistence under symmetric high signals (test 13)
- sum constraint NOT enforced (test 5)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import acs.v2 as v2_pkg
from acs.v2 import dynamics as v2_dynamics
from acs.v2.adhesion_network_state import (
    AdhesionNetworkRates,
    AdhesionNetworkSignals,
    AdhesionNetworkStepDiagnostics,
    AdhesionNetworkStepResult,
    CellAdhesionNetworkState,
)
from acs.v2.dynamics.adhesion_network_dynamics import (
    compute_fa_engagement_signal,
    compute_junction_engagement_signal,
    step_adhesion_network_state,
)
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.junction import JunctionState


def _make_state(
    *, cell_id: str = "cell-an", cadherin: float = 0.5, integrin: float = 0.5
) -> CellAdhesionNetworkState:
    return CellAdhesionNetworkState(
        cell_id=cell_id,
        cadherin_program=cadherin,
        integrin_program=integrin,
    )


def _make_rates(value: float = 1e-3) -> AdhesionNetworkRates:
    return AdhesionNetworkRates(
        k_int_on=value,
        k_int_off=value,
        k_cad_on=value,
        k_cad_off=value,
    )


# ---------------------------------------------------------------------------
# State validation (5)
# ---------------------------------------------------------------------------


def test_cell_adhesion_network_state_cadherin_in_unit_interval():
    """Test 1: cadherin_program ∉ [0, 1] -> ValueError."""
    with pytest.raises(ValueError, match="cadherin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=-0.1, integrin_program=0.5)
    with pytest.raises(ValueError, match="cadherin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=1.5, integrin_program=0.5)


def test_cell_adhesion_network_state_integrin_in_unit_interval():
    """Test 2: integrin_program ∉ [0, 1] -> ValueError."""
    with pytest.raises(ValueError, match="integrin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=0.5, integrin_program=-0.1)
    with pytest.raises(ValueError, match="integrin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=0.5, integrin_program=1.5)


def test_cell_adhesion_network_state_bool_program_rejected():
    """Test 3: Python bool ⊂ int trap; both programs reject True/False."""
    with pytest.raises(ValueError, match="cadherin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=True, integrin_program=0.5)
    with pytest.raises(ValueError, match="integrin_program"):
        CellAdhesionNetworkState(cell_id="c", cadherin_program=0.5, integrin_program=False)


def test_cell_adhesion_network_state_empty_cell_id_raises():
    """Test 4: empty cell_id -> ValueError."""
    with pytest.raises(ValueError, match="cell_id"):
        CellAdhesionNetworkState(cell_id="", cadherin_program=0.5, integrin_program=0.5)


def test_cell_adhesion_network_state_sum_constraint_not_enforced():
    """Test 5 (Y16 PI directive forward guard): cadherin=0.7, integrin=0.7
    (sum=1.4) is valid; __post_init__ does NOT reject. Coexistence
    representable per PI id=2003."""
    state = CellAdhesionNetworkState(
        cell_id="c", cadherin_program=0.7, integrin_program=0.7
    )
    assert state.cadherin_program == 0.7
    assert state.integrin_program == 0.7


# ---------------------------------------------------------------------------
# Derived φ diagnostic (2)
# ---------------------------------------------------------------------------


def test_phi_integrin_dominance_both_zero_returns_none():
    """Test 6 (Y11): both-zero -> phi_integrin_dominance is None
    (NOT 0.0; would falsely imply E-cad dominance)."""
    state = CellAdhesionNetworkState(
        cell_id="c", cadherin_program=0.0, integrin_program=0.0
    )
    assert state.phi_integrin_dominance is None


def test_phi_integrin_dominance_general_returns_ratio():
    """Test 7: cadherin=0.3, integrin=0.7 -> phi == 0.7 / 1.0 == 0.7."""
    state = CellAdhesionNetworkState(
        cell_id="c", cadherin_program=0.3, integrin_program=0.7
    )
    assert state.phi_integrin_dominance == pytest.approx(0.7, rel=1e-12)


# ---------------------------------------------------------------------------
# Signals + Rates validation (2)
# ---------------------------------------------------------------------------


def test_adhesion_network_signals_validation():
    """Test 8 (Y14): bool rejected, finite, [0, 1] for both fields."""
    AdhesionNetworkSignals(fa_engagement_signal=0.5, junction_engagement_signal=0.3)
    with pytest.raises(ValueError, match="fa_engagement_signal"):
        AdhesionNetworkSignals(fa_engagement_signal=True, junction_engagement_signal=0.5)
    with pytest.raises(ValueError, match="junction_engagement_signal"):
        AdhesionNetworkSignals(
            fa_engagement_signal=0.5, junction_engagement_signal=float("nan")
        )
    with pytest.raises(ValueError, match="fa_engagement_signal"):
        AdhesionNetworkSignals(fa_engagement_signal=-0.1, junction_engagement_signal=0.5)
    with pytest.raises(ValueError, match="junction_engagement_signal"):
        AdhesionNetworkSignals(fa_engagement_signal=0.5, junction_engagement_signal=1.1)


def test_adhesion_network_rates_validation():
    """Test 9 (Y6): all 4 rates required, finite positive, bool rejected."""
    AdhesionNetworkRates(k_int_on=1e-3, k_int_off=1e-3, k_cad_on=1e-3, k_cad_off=1e-3)
    with pytest.raises(ValueError, match="k_int_on"):
        AdhesionNetworkRates(k_int_on=True, k_int_off=1e-3, k_cad_on=1e-3, k_cad_off=1e-3)
    with pytest.raises(ValueError, match="k_int_off"):
        AdhesionNetworkRates(k_int_on=1e-3, k_int_off=0.0, k_cad_on=1e-3, k_cad_off=1e-3)
    with pytest.raises(ValueError, match="k_cad_on"):
        AdhesionNetworkRates(k_int_on=1e-3, k_int_off=1e-3, k_cad_on=-1.0, k_cad_off=1e-3)
    with pytest.raises(ValueError, match="k_cad_off"):
        AdhesionNetworkRates(
            k_int_on=1e-3, k_int_off=1e-3, k_cad_on=1e-3, k_cad_off=float("inf")
        )


# ---------------------------------------------------------------------------
# ODE evolution (5)
# ---------------------------------------------------------------------------


def test_step_adhesion_network_state_both_signals_zero_no_op():
    """Test 10 (Y9 boundary): signals = (0, 0) -> state unchanged exactly;
    diagnostics steady states None + relaxation rates 0."""
    state = _make_state(cadherin=0.4, integrin=0.6)
    signals = AdhesionNetworkSignals(fa_engagement_signal=0.0, junction_engagement_signal=0.0)
    rates = _make_rates(1e-3)
    result = step_adhesion_network_state(state, signals, rates, dt_s=10.0)

    assert result.state.cadherin_program == 0.4
    assert result.state.integrin_program == 0.6
    assert result.diagnostics.integrin_steady_state is None
    assert result.diagnostics.cadherin_steady_state is None
    assert result.diagnostics.integrin_relaxation_rate_per_s == 0.0
    assert result.diagnostics.cadherin_relaxation_rate_per_s == 0.0


def test_step_adhesion_network_state_fa_only_integrin_up_cadherin_down():
    """Test 11 (Codex id=2034 fixture correction): initial cadherin=0.7,
    integrin=0.0; S_fa=1.0, S_junction=0.0; equal rates; assert integrin
    increases AND cadherin decreases (cross term active)."""
    state = _make_state(cadherin=0.7, integrin=0.0)
    signals = AdhesionNetworkSignals(fa_engagement_signal=1.0, junction_engagement_signal=0.0)
    rates = _make_rates(1.0)
    result = step_adhesion_network_state(state, signals, rates, dt_s=0.5)

    assert result.state.integrin_program > 0.0
    assert result.state.cadherin_program < 0.7


def test_step_adhesion_network_state_junction_only_cadherin_up_integrin_down():
    """Test 12 (Codex id=2034 fixture correction): initial integrin=0.7,
    cadherin=0.0; S_fa=0.0, S_junction=1.0; equal rates; assert cadherin
    increases AND integrin decreases."""
    state = _make_state(cadherin=0.0, integrin=0.7)
    signals = AdhesionNetworkSignals(fa_engagement_signal=0.0, junction_engagement_signal=1.0)
    rates = _make_rates(1.0)
    result = step_adhesion_network_state(state, signals, rates, dt_s=0.5)

    assert result.state.cadherin_program > 0.0
    assert result.state.integrin_program < 0.7


def test_step_adhesion_network_state_coexistence_under_symmetric_high_signals():
    """Test 13 (Y15 — PI directive forward guard): both signals=1.0,
    all 4 rates equal, run 100 steps; assert both programs converge to
    ≈0.5 (coexistence, NOT winner-take-all)."""
    state = _make_state(cadherin=0.0, integrin=0.0)
    signals = AdhesionNetworkSignals(fa_engagement_signal=1.0, junction_engagement_signal=1.0)
    rates = _make_rates(1.0)
    for _ in range(100):
        result = step_adhesion_network_state(state, signals, rates, dt_s=0.5)
        state = result.state

    assert state.cadherin_program == pytest.approx(0.5, abs=1e-6)
    assert state.integrin_program == pytest.approx(0.5, abs=1e-6)


def test_step_adhesion_network_state_large_dt_stays_bounded():
    """Test 14 (Y17): dt_s=1e6 with arbitrary rates and signals; assert
    both programs ∈ [0, 1] after step (analytical boundedness, not
    approximate)."""
    state = _make_state(cadherin=0.3, integrin=0.7)
    signals = AdhesionNetworkSignals(fa_engagement_signal=0.6, junction_engagement_signal=0.4)
    rates = AdhesionNetworkRates(
        k_int_on=2.0, k_int_off=1.0, k_cad_on=1.5, k_cad_off=2.5
    )
    result = step_adhesion_network_state(state, signals, rates, dt_s=1e6)

    assert 0.0 <= result.state.cadherin_program <= 1.0
    assert 0.0 <= result.state.integrin_program <= 1.0


# ---------------------------------------------------------------------------
# Junction signal fail-closed + dt validation (2)
# ---------------------------------------------------------------------------


def test_compute_junction_engagement_signal_none_cadherin_proxy_raises():
    """Test 15 (Y13): junction with cadherin_proxy=None -> ValueError
    (no silent zero coercion)."""
    j = JunctionState(
        junction_id="j-1",
        cell_id_a="cell-a",
        cell_id_b="cell-b",
        contact_length_um=1.0,
        age_s=0.0,
        maturity=0.5,
        cadherin_proxy=None,
    )
    with pytest.raises(ValueError, match="cadherin_proxy"):
        compute_junction_engagement_signal((j,))


def test_step_adhesion_network_state_dt_validation():
    """Test 16 (Y14): bool dt_s rejected, dt_s=-1.0 rejected, dt_s=inf
    rejected, dt_s=0.0 valid no-op."""
    state = _make_state()
    signals = AdhesionNetworkSignals(fa_engagement_signal=0.5, junction_engagement_signal=0.5)
    rates = _make_rates(1e-3)

    with pytest.raises(ValueError, match="dt_s"):
        step_adhesion_network_state(state, signals, rates, dt_s=True)
    with pytest.raises(ValueError, match="dt_s"):
        step_adhesion_network_state(state, signals, rates, dt_s=-1.0)
    with pytest.raises(ValueError, match="dt_s"):
        step_adhesion_network_state(state, signals, rates, dt_s=float("inf"))

    # dt_s=0.0 valid no-op
    result = step_adhesion_network_state(state, signals, rates, dt_s=0.0)
    assert result.state.cadherin_program == state.cadherin_program
    assert result.state.integrin_program == state.integrin_program


# ---------------------------------------------------------------------------
# Helpful additional tests
# ---------------------------------------------------------------------------


def test_compute_fa_engagement_signal_empty_returns_zero():
    """Empty FA list -> 0.0 (no integrin engagement signal)."""
    assert compute_fa_engagement_signal(()) == 0.0


def test_compute_junction_engagement_signal_empty_returns_zero():
    """Empty junction list -> 0.0 (no cadherin engagement signal)."""
    assert compute_junction_engagement_signal(()) == 0.0


def test_compute_fa_engagement_signal_invalid_fa_raises(monkeypatch):
    """Codex `id=2079` P0 BLOCKER fix regression: invalid FA schema
    (maturity > 1) must raise ValueError at the helper boundary
    instead of producing an out-of-range signal (e.g., 2.0)."""
    # Build a FocalAdhesionState that bypasses normal construction
    # validation, then verify the helper enforces validate() before
    # consuming its fields.
    fa = FocalAdhesionState(
        adhesion_id="fa-bad",
        cell_id="cell-x",
        position_um_xy=(0.0, 0.0),
        age_s=0.0,
        maturity=1.0,  # construct valid
        bound_fraction=1.0,
        state="mature",
        traction_force_nN_xy=(0.0, 0.0),
    )
    # Mutate maturity past schema bound by replacing dataclass field.
    object.__setattr__(fa, "maturity", 2.0)
    with pytest.raises(ValueError, match="maturity"):
        compute_fa_engagement_signal((fa,))


def test_compute_junction_engagement_signal_invalid_junction_raises():
    """Codex `id=2079` P0 BLOCKER fix regression: invalid Junction
    schema (cadherin_proxy > 1) must raise ValueError at the helper
    boundary."""
    j = JunctionState(
        junction_id="j-bad",
        cell_id_a="cell-a",
        cell_id_b="cell-b",
        contact_length_um=1.0,
        age_s=0.0,
        maturity=0.5,
        cadherin_proxy=0.5,
    )
    object.__setattr__(j, "cadherin_proxy", 2.0)
    with pytest.raises(ValueError, match="cadherin_proxy"):
        compute_junction_engagement_signal((j,))


def test_exports_through_both_init():
    """All public symbols accessible through both surfaces."""
    expected_state_classes = {
        "CellAdhesionNetworkState",
        "AdhesionNetworkSignals",
        "AdhesionNetworkRates",
        "AdhesionNetworkStepDiagnostics",
        "AdhesionNetworkStepResult",
    }
    expected_dynamics_funcs = {
        "step_adhesion_network_state",
        "compute_fa_engagement_signal",
        "compute_junction_engagement_signal",
    }
    for sym in expected_state_classes:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym}"
        assert sym in v2_pkg.__all__
    for sym in expected_dynamics_funcs:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym}"
        assert hasattr(v2_dynamics, sym), f"acs.v2.dynamics missing {sym}"
        assert sym in v2_pkg.__all__
        assert sym in v2_dynamics.__all__
