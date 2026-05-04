"""Tests for 6.3b protrusion-coupled FA dynamics
(``acs.v2.dynamics.protrusion_coupled_focal_adhesion``).

Test catalog mirrors ``docs/v2_63b_protrusion_coupled_fa_sanity_gate.md`` §9.
Each test maps to one Sanity Gate item or one forbidden behavior from
the locked design (``docs/v2_63b_protrusion_coupling_locked.md``).
"""

from __future__ import annotations

import numpy as np
import pytest

from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    step_focal_adhesions_static,
)
from acs.v2.dynamics.protrusion_coupled_focal_adhesion import (
    ProtrusionStateMultipliers,
    step_protrusion_coupled_focal_adhesions,
)
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.protrusion import ProtrusionEvent

CENTROID = (0.0, 0.0)


def _make_fa(
    adhesion_id: str,
    *,
    position=(1.0, 0.0),
    state="mature",
    maturity=0.4,
    bound_fraction=0.5,
    linked_protrusion_id=None,
    cell_id="cell-A",
):
    return FocalAdhesionState(
        adhesion_id=adhesion_id,
        cell_id=cell_id,
        position_um_xy=position,
        age_s=0.0,
        maturity=maturity,
        bound_fraction=bound_fraction,
        state=state,
        linked_protrusion_id=linked_protrusion_id,
    )


def _make_protrusion(
    event_id: str,
    *,
    state="growing",
    associated_adhesion_ids=(),
    cell_id="cell-A",
):
    return ProtrusionEvent(
        cell_id=cell_id,
        start_time_s=0.0,
        boundary_angle_rad=0.0,
        length_um=1.0,
        event_id=event_id,
        state=state,
        associated_adhesion_ids=associated_adhesion_ids,
    )


def _base_params(*, k_maturity=0.1, k_bind=0.05, k_unbind=None):
    return FocalAdhesionDynamicsParameters(
        dt_fa_s=0.5,
        traction_scale_nN=1.0,
        k_maturity_per_s=k_maturity,
        k_bind_per_s=k_bind,
        k_unbind_per_s=k_unbind,
    )


# ---------------------------------------------------------------------------
# Test 1 (gate §3 determinism + §5 multiplier 1.0 = neutral)
# ---------------------------------------------------------------------------


def test_neutral_multipliers_match_6_3a_baseline():
    """6.3b with all-1.0 multiplier table reproduces 6.3a output bit-for-bit."""

    fa1 = _make_fa("fa-1", linked_protrusion_id="p-1")
    fa2 = _make_fa("fa-2", position=(0.0, 1.0), linked_protrusion_id="p-1")
    params = _base_params()
    registry = {"p-1": _make_protrusion("p-1", state="growing")}
    neutral = ProtrusionStateMultipliers(
        table={
            "growing": {
                "k_maturity_per_s": 1.0,
                "k_bind_per_s": 1.0,
                "k_unbind_per_s": 1.0,
            }
        }
    )

    coupled = step_protrusion_coupled_focal_adhesions(
        (fa1, fa2), CENTROID, params, registry, neutral
    )
    baseline = step_focal_adhesions_static((fa1, fa2), CENTROID, params)

    np.testing.assert_array_equal(
        coupled.cell_force_nN_xy, baseline.cell_force_nN_xy
    )
    np.testing.assert_array_equal(
        coupled.substrate_reaction_nN_xy,
        baseline.substrate_reaction_nN_xy,
    )
    for c, b in zip(coupled.updated_adhesions, baseline.updated_adhesions):
        assert c.maturity == b.maturity
        assert c.bound_fraction == b.bound_fraction
        assert c.state == b.state


# ---------------------------------------------------------------------------
# Test 2 (gate §1 unit chain + §5 multiplier > 1 = boost)
# ---------------------------------------------------------------------------


def test_growing_multiplier_boosts_maturity_and_bind():
    """k_maturity multiplier 2.0 produces 2× the maturity increment vs 1.0."""

    fa = _make_fa("fa-1", linked_protrusion_id="p-1")
    params = _base_params(k_maturity=0.1, k_bind=0.1, k_unbind=None)
    registry = {"p-1": _make_protrusion("p-1", state="growing")}
    boosted = ProtrusionStateMultipliers(
        table={"growing": {"k_maturity_per_s": 2.0, "k_bind_per_s": 2.0}}
    )

    boosted_run = step_protrusion_coupled_focal_adhesions(
        (fa,), CENTROID, params, registry, boosted
    )
    base_run = step_protrusion_coupled_focal_adhesions(
        (fa,),
        CENTROID,
        params,
        registry,
        ProtrusionStateMultipliers(table={}),  # all neutral
    )

    boosted_dmaturity = boosted_run.updated_adhesions[0].maturity - fa.maturity
    base_dmaturity = base_run.updated_adhesions[0].maturity - fa.maturity
    boosted_dbound = (
        boosted_run.updated_adhesions[0].bound_fraction - fa.bound_fraction
    )
    base_dbound = (
        base_run.updated_adhesions[0].bound_fraction - fa.bound_fraction
    )

    assert boosted_dmaturity == pytest.approx(2.0 * base_dmaturity)
    assert boosted_dbound == pytest.approx(2.0 * base_dbound)


# ---------------------------------------------------------------------------
# Test 3 (gate §5 multiplier 1.0 = neutral on a non-growing state)
# ---------------------------------------------------------------------------


def test_retracting_neutral_no_boost():
    """Retracting protrusion with neutral table = base 6.3a behavior for FA."""

    fa = _make_fa("fa-1", linked_protrusion_id="p-1")
    params = _base_params()
    registry = {"p-1": _make_protrusion("p-1", state="retracting")}
    table = ProtrusionStateMultipliers(
        table={"retracting": {"k_maturity_per_s": 1.0}}
    )

    coupled = step_protrusion_coupled_focal_adhesions(
        (fa,), CENTROID, params, registry, table
    )
    baseline = step_focal_adhesions_static((fa,), CENTROID, params)

    assert (
        coupled.updated_adhesions[0].maturity
        == baseline.updated_adhesions[0].maturity
    )
    assert (
        coupled.updated_adhesions[0].bound_fraction
        == baseline.updated_adhesions[0].bound_fraction
    )


# ---------------------------------------------------------------------------
# Test 4 (gate §2 typed key validation: unknown state)
# ---------------------------------------------------------------------------


def test_unknown_state_key_raises_multiplier_table_unknown_key():
    table = ProtrusionStateMultipliers(
        table={"sprinting": {"k_maturity_per_s": 2.0}}  # not a known state
    )
    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        table.validate()
    assert excinfo.value.failure_kind == "multiplier_table_unknown_key"
    assert "sprinting" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Test 5 (gate §2 typed key validation: unknown rate)
# ---------------------------------------------------------------------------


def test_unknown_rate_key_raises_multiplier_table_unknown_key():
    table = ProtrusionStateMultipliers(
        table={"growing": {"k_bind": 2.0}}  # missing _per_s suffix
    )
    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        table.validate()
    assert excinfo.value.failure_kind == "multiplier_table_unknown_key"
    assert "k_bind" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Test 6 (gate §2 multiplier ≥ 0)
# ---------------------------------------------------------------------------


def test_negative_multiplier_raises_at_validate():
    table = ProtrusionStateMultipliers(
        table={"growing": {"k_maturity_per_s": -0.5}}
    )
    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        table.validate()
    assert excinfo.value.failure_kind == "multiplier_table_unknown_key"
    assert "non-negative" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Test 7 (gate §2 bool not allowed)
# ---------------------------------------------------------------------------


def test_bool_multiplier_value_raises_at_validate():
    table = ProtrusionStateMultipliers(
        table={"growing": {"k_maturity_per_s": True}}  # type: ignore[dict-item]
    )
    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        table.validate()
    assert excinfo.value.failure_kind == "multiplier_table_unknown_key"


# ---------------------------------------------------------------------------
# Test 8 (gate §2 linked_protrusion_missing)
# ---------------------------------------------------------------------------


def test_missing_linked_protrusion_id_raises_linked_protrusion_missing():
    fa = _make_fa("fa-1", linked_protrusion_id="p-ghost")
    params = _base_params()
    registry: dict[str, ProtrusionEvent] = {}  # FA links to ghost
    table = ProtrusionStateMultipliers(table={})

    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        step_protrusion_coupled_focal_adhesions(
            (fa,), CENTROID, params, registry, table
        )
    assert excinfo.value.failure_kind == "linked_protrusion_missing"
    assert "p-ghost" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Test 9 (gate §2 / §6 asymmetric link consistency)
# ---------------------------------------------------------------------------


def test_reciprocal_missing_records_diagnostic_not_raise():
    """FA links to protrusion, but protrusion's associated_adhesion_ids does
    not list the FA → diagnostic counter, no failure."""

    fa = _make_fa("fa-1", linked_protrusion_id="p-1")
    params = _base_params()
    registry = {
        "p-1": _make_protrusion(
            "p-1", state="growing", associated_adhesion_ids=()
        ),  # missing back-link
    }
    table = ProtrusionStateMultipliers(table={})

    result = step_protrusion_coupled_focal_adhesions(
        (fa,), CENTROID, params, registry, table
    )
    assert result.diagnostics["reciprocal_missing"] == 1
    assert result.diagnostics["linked_missing"] == 0


# ---------------------------------------------------------------------------
# Test 10 (gate §4 pre-step gate over all FAs/rate names)
# ---------------------------------------------------------------------------


def test_dt_rate_violation_uses_max_effective_rate_over_all_fas():
    """One FA's high-multiplier protrusion pushes effective rate past the
    global gate even if another FA's rate is well within it."""

    fa1 = _make_fa("fa-1", linked_protrusion_id="p-low")
    fa2 = _make_fa("fa-2", position=(0.0, 1.0), linked_protrusion_id="p-high")
    params = _base_params(k_maturity=0.5, k_bind=0.5, k_unbind=None)
    # base_rate 0.5/s, dt 0.5s → base dt·rate = 0.25 (well within 0.5)
    # multiplier 4 → effective dt·rate = 1.0 > 0.5 → must raise
    registry = {
        "p-low": _make_protrusion("p-low", state="growing"),
        "p-high": _make_protrusion("p-high", state="stalled"),
    }
    table = ProtrusionStateMultipliers(
        table={
            "growing": {"k_maturity_per_s": 1.0, "k_bind_per_s": 1.0},
            "stalled": {"k_maturity_per_s": 4.0, "k_bind_per_s": 1.0},
        }
    )

    with pytest.raises(FocalAdhesionDynamicsError) as excinfo:
        step_protrusion_coupled_focal_adhesions(
            (fa1, fa2), CENTROID, params, registry, table
        )
    assert excinfo.value.failure_kind == "dt_rate_violation"


# ---------------------------------------------------------------------------
# Test 11 (gate §6 input order preservation)
# ---------------------------------------------------------------------------


def test_per_fa_input_order_preserved_in_output():
    fa1 = _make_fa(
        "alpha", position=(1.0, 0.0), linked_protrusion_id="p-1"
    )
    fa2 = _make_fa(
        "beta", position=(0.0, 1.0), linked_protrusion_id="p-2"
    )
    fa3 = _make_fa(
        "gamma", position=(-1.0, 0.0), linked_protrusion_id="p-1"
    )
    params = _base_params()
    registry = {
        "p-1": _make_protrusion("p-1", state="growing"),
        "p-2": _make_protrusion("p-2", state="retracting"),
    }
    table = ProtrusionStateMultipliers(
        table={"growing": {"k_maturity_per_s": 2.0}}
    )

    result = step_protrusion_coupled_focal_adhesions(
        (fa1, fa2, fa3), CENTROID, params, registry, table
    )
    assert [fa.adhesion_id for fa in result.updated_adhesions] == [
        "alpha",
        "beta",
        "gamma",
    ]


# ---------------------------------------------------------------------------
# Test 12 (§0 forbidden + §5 multiplier 0 ≠ state change)
# ---------------------------------------------------------------------------


def test_no_state_label_change_after_step():
    """Even multiplier 0 (no rate-driven update) does not change the FA
    state label. The state label is only set by the caller."""

    fa = _make_fa(
        "fa-1", state="mature", linked_protrusion_id="p-1"
    )
    params = _base_params()
    registry = {"p-1": _make_protrusion("p-1", state="growing")}
    table = ProtrusionStateMultipliers(
        table={
            "growing": {
                "k_maturity_per_s": 0.0,
                "k_bind_per_s": 0.0,
                "k_unbind_per_s": 0.0,
            }
        }
    )

    result = step_protrusion_coupled_focal_adhesions(
        (fa,), CENTROID, params, registry, table
    )
    # State label preserved
    assert result.updated_adhesions[0].state == "mature"
    # Maturity / bound_fraction unchanged because multiplier 0 zeros
    # the effective rate
    assert result.updated_adhesions[0].maturity == fa.maturity
    assert result.updated_adhesions[0].bound_fraction == fa.bound_fraction
