"""Structural gates for the membrane effective-rest-area law (PI `AC_EXECUTION_PLAN_2026-07-25.md` §9 D2 option B, 2026-07-28 — a
DIFFERENT D2 from `COMPARTMENT_VALIDATION_TRACKS` §6 D2, which is the force-assembly omit mask).

Every gate here is a SHAPE or a LIMIT — continuity, monotonicity, the two reservoir limits, sign — and
none is an absolute magnitude, which is what makes them evaluable under the 2026-07-25 gate.
"""

from __future__ import annotations

import math

import pytest

from aleph.engine.membrane_area import (
    MembraneAreaCard,
    area_stiffness,
    area_tension,
    effective_rest_area,
)

A0 = 100.0                     # um^2 reference area; any positive value works, none is asserted


def _card(capacity: float) -> MembraneAreaCard:
    return MembraneAreaCard(
        gamma_mem_pn_per_um=10.0,      # KB-3.B1.1 plateau; the gates never assert this value
        k_area_pn_per_um=2.4e5,        # Rawicz K_A scale
        tau_lysis_pn_per_um=1.0e4,
        capacity=capacity,
        capacity_provenance="sweep point, test",
    )


def test_infinite_capacity_reproduces_the_incumbent_plateau_exactly() -> None:
    """The change must be a strict GENERALISATION, not a different model.

    With an unlimited reservoir the law is the constant-tension plateau the `ac/` path passes today, so
    no existing result can move until a caller deliberately declares a finite capacity.
    """
    card = _card(math.inf)
    for area in (A0, A0 * 1.05, A0 * 2.0, A0 * 100.0):
        assert area_tension(area, A0, card) == card.gamma_mem_pn_per_um
        assert area_stiffness(area, A0, card) == 0.0     # dsigma/dA identically zero — the defect


def test_zero_capacity_recovers_pure_elasticity_from_the_rest_area() -> None:
    card = _card(0.0)
    assert area_tension(A0, A0, card) == card.gamma_mem_pn_per_um     # unstretched: plateau
    strain = 0.01
    expected = card.gamma_mem_pn_per_um + card.k_area_pn_per_um * strain
    assert area_tension(A0 * (1.0 + strain), A0, card) == pytest.approx(
        min(expected, card.tau_lysis_pn_per_um)
    )


def test_the_law_is_continuous_at_the_capacity_boundary() -> None:
    """A jump at the boundary would be a spurious force, not a stiffening membrane."""
    capacity = 0.35
    card = _card(capacity)
    edge = A0 * (1.0 + capacity)
    below = area_tension(edge * (1.0 - 1e-12), A0, card)
    at = area_tension(edge, A0, card)
    above = area_tension(edge * (1.0 + 1e-12), A0, card)
    assert below == pytest.approx(card.gamma_mem_pn_per_um, rel=1e-9)
    assert at == pytest.approx(card.gamma_mem_pn_per_um, rel=1e-9)
    assert above == pytest.approx(card.gamma_mem_pn_per_um, rel=1e-6)


def test_tension_is_monotone_and_gains_real_stiffness_past_capacity() -> None:
    """The point of D2: dsigma/dA is EXACTLY zero inside the reservoir and strictly positive outside."""
    capacity = 0.20
    card = _card(capacity)
    areas = [A0 * f for f in (1.0, 1.1, 1.2, 1.2001, 1.21, 1.3)]
    tensions = [area_tension(a, A0, card) for a in areas]
    assert tensions == sorted(tensions)                       # non-decreasing everywhere

    assert area_stiffness(A0 * 1.1, A0, card) == 0.0          # inside the reservoir
    assert area_stiffness(A0 * 1.2001, A0, card) > 0.0        # outside it — the missing stiffness
    # and the plateau really is flat inside, not merely slowly varying
    assert area_tension(A0 * 1.05, A0, card) == area_tension(A0 * 1.19, A0, card)


def test_a_slack_membrane_is_not_in_compression() -> None:
    """Below A0 a bilayer goes slack; it cannot carry a negative tension."""
    card = _card(0.2)
    assert effective_rest_area(A0 * 0.5, A0, card.capacity) == A0
    assert area_tension(A0 * 0.5, A0, card) == card.gamma_mem_pn_per_um


def test_tension_saturates_at_lysis_and_stops_stiffening() -> None:
    card = _card(0.0)
    huge = A0 * 10.0
    assert area_tension(huge, A0, card) == card.tau_lysis_pn_per_um
    assert area_stiffness(huge, A0, card) == 0.0


def test_capacity_carries_no_default_and_no_bare_number() -> None:
    """Supplying a default would recreate RESERVOIR_STRAIN=0.60 under a new name."""
    import inspect

    params = inspect.signature(MembraneAreaCard).parameters
    assert params["capacity"].default is inspect.Parameter.empty
    assert params["capacity_provenance"].default is inspect.Parameter.empty

    with pytest.raises(ValueError, match="under a new name"):
        MembraneAreaCard(10.0, 2.4e5, 1.0e4, 0.6, "   ")
    with pytest.raises(ValueError, match="cannot be negative"):
        MembraneAreaCard(10.0, 2.4e5, 1.0e4, -0.1, "negative reservoir")
    with pytest.raises(ValueError, match="already lysed"):
        MembraneAreaCard(10.0, 2.4e5, 5.0, 0.2, "lysis below plateau")


def test_capacity_is_a_sweep_axis_that_orders_the_response() -> None:
    """Recovering capacity as an OUTPUT requires the observable to depend on it monotonically."""
    area = A0 * 1.30
    tensions = [area_tension(area, A0, _card(c)) for c in (0.0, 0.1, 0.2, 0.3, 0.5)]
    assert tensions == sorted(tensions, reverse=True)     # more reservoir -> less tension at fixed area
    assert tensions[-1] == _card(0.5).gamma_mem_pn_per_um  # capacity beyond the strain: still plateau
