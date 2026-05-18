"""Unit tests for Phase 1 Unit 4.1 Young-equation contact angle.

These tests exercise the formula in isolation: each test fixes
``gamma_J`` explicitly via the params override and verifies that the
returned θ matches the brief's

.. math::
    \\cos\\theta = (\\gamma_{c,1} + \\gamma_{c,2} - 2\\gamma_J) / (2\\gamma_c).

Coupling to the n_engaged-driven γ_J interpolation is verified
separately in :mod:`acs_kb.tests.test_two_cell_pair`.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from acs_kb.cell.cell import Cell
from acs_kb.cell.cortex import generate_cortex
from acs_kb.junction.contact_angle import compute_contact_angle, junction_tension
from acs_kb.junction.types import make_ecadherin_junction


def _pair(*, gamma_a: float = 5.0e-4, gamma_b: float = 5.0e-4):
    """Build a deterministic two-cell pair for contact-angle tests."""
    ctx_a = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([-9e-6, 0.0]), seed=1,
        params={"gamma_cortex": gamma_a},
    )
    ctx_b = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([+9e-6, 0.0]), seed=2,
        params={"gamma_cortex": gamma_b},
    )
    cell_a = Cell.from_cortex(0, ctx_a)
    cell_b = Cell.from_cortex(1, ctx_b)
    j = make_ecadherin_junction(cell_a, cell_b, n_bonds_total=100)
    return cell_a, cell_b, j


def test_symmetric_gamma_j_equal_gamma_c_gives_pi_over_two() -> None:
    """γ_J = γ_c (symmetric mature contact) ⇒ θ = π/2 (KU-4.4 nominal)."""
    cell_a, cell_b, j = _pair()
    theta = compute_contact_angle(
        cell_a, cell_b, j,
        params={"gamma_c_a": 5e-4, "gamma_c_b": 5e-4, "gamma_J": 5e-4},
    )
    assert theta == pytest.approx(math.pi / 2.0, abs=1e-12)


def test_full_adhesion_limit_zero_gamma_j() -> None:
    """γ_J = 0 ⇒ θ = 0 (brief's full-adhesion limit)."""
    cell_a, cell_b, j = _pair()
    theta = compute_contact_angle(
        cell_a, cell_b, j,
        params={"gamma_c_a": 5e-4, "gamma_c_b": 5e-4, "gamma_J": 0.0},
    )
    assert theta == pytest.approx(0.0, abs=1e-12)


def test_no_adhesion_limit_gamma_j_equals_sum_of_cortices() -> None:
    """γ_J = γ_c1 + γ_c2 ⇒ θ = π (brief's no-adhesion limit)."""
    cell_a, cell_b, j = _pair()
    theta = compute_contact_angle(
        cell_a, cell_b, j,
        params={"gamma_c_a": 5e-4, "gamma_c_b": 5e-4, "gamma_J": 1e-3},
    )
    assert theta == pytest.approx(math.pi, abs=1e-12)


def test_intermediate_gamma_j_matches_arccos_formula() -> None:
    """Spot-check the formula at a non-trivial γ_J value."""
    cell_a, cell_b, j = _pair()
    gJ = 0.3e-3
    gc = 5e-4
    theta = compute_contact_angle(
        cell_a, cell_b, j,
        params={"gamma_c_a": gc, "gamma_c_b": gc, "gamma_J": gJ},
    )
    expected = math.acos((2 * gc - 2 * gJ) / (2 * gc))
    assert theta == pytest.approx(expected, abs=1e-12)


def test_asymmetric_cortical_tensions() -> None:
    """The formula still holds for asymmetric γ_c (Maître 2012 Fig 2g)."""
    cell_a, cell_b, j = _pair(gamma_a=4e-4, gamma_b=6e-4)
    gJ = 5e-4
    theta = compute_contact_angle(
        cell_a, cell_b, j,
        params={"gamma_c_a": 4e-4, "gamma_c_b": 6e-4, "gamma_J": gJ},
    )
    # Mean cortex tension is 5e-4; (4 + 6 - 2·5)/(2·5) = 0 ⇒ θ = π/2.
    assert theta == pytest.approx(math.pi / 2.0, abs=1e-12)


def test_cortex_params_fallback_picks_up_gamma_cortex() -> None:
    """If params omits γ_c_a / γ_c_b they are read from cortex.params."""
    cell_a, cell_b, j = _pair(gamma_a=5e-4, gamma_b=5e-4)
    # No gamma_c_a / gamma_c_b in params; the function should fall back
    # to cortex.params['gamma_cortex'] (= 5e-4 by the _pair fixture).
    theta = compute_contact_angle(
        cell_a, cell_b, j, params={"gamma_J": 5e-4},
    )
    assert theta == pytest.approx(math.pi / 2.0, abs=1e-12)


def test_missing_cortex_tension_raises_keyerror() -> None:
    """If both params and cortex.params are missing γ_c, raise KeyError."""
    # Build a cortex without gamma_cortex in params.
    ctx = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([-9e-6, 0.0]), seed=1, params={},
    )
    cell_a = Cell.from_cortex(0, ctx)
    ctx_b = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([+9e-6, 0.0]), seed=2, params={},
    )
    cell_b = Cell.from_cortex(1, ctx_b)
    j = make_ecadherin_junction(cell_a, cell_b, n_bonds_total=100)
    with pytest.raises(KeyError):
        compute_contact_angle(cell_a, cell_b, j, params={"gamma_J": 5e-4})


def test_junction_cell_pair_mismatch_raises() -> None:
    cell_a, cell_b, j = _pair()
    # Build a third cell with id=99 and try to evaluate the junction on
    # the wrong pair.
    ctx_c = generate_cortex(
        R_cell=10e-6, n_cortex_fibers=20, L_cortex_fiber=3e-6, beads_per_fiber=5,
        cell_center=np.array([0.0, 0.0]), seed=3, params={"gamma_cortex": 5e-4},
    )
    cell_c = Cell.from_cortex(99, ctx_c)
    with pytest.raises(ValueError, match="junction connects cells"):
        compute_contact_angle(
            cell_a, cell_c, j,
            params={"gamma_c_a": 5e-4, "gamma_c_b": 5e-4, "gamma_J": 5e-4},
        )


def test_zero_cortex_tension_raises() -> None:
    cell_a, cell_b, j = _pair()
    with pytest.raises(ValueError, match="cortical tensions must be positive"):
        compute_contact_angle(
            cell_a, cell_b, j,
            params={"gamma_c_a": 0.0, "gamma_c_b": 5e-4, "gamma_J": 5e-4},
        )


def test_junction_tension_endpoints() -> None:
    """``junction_tension`` returns the linear-interpolation endpoints."""
    assert junction_tension(0, 100, 5e-4, 5e-4) == pytest.approx(1.0e-3)
    assert junction_tension(50, 100, 5e-4, 5e-4) == pytest.approx(5.0e-4)
    assert junction_tension(100, 100, 5e-4, 5e-4) == pytest.approx(0.0, abs=1e-12)


def test_junction_tension_rejects_zero_n_total() -> None:
    with pytest.raises(ValueError, match="n_total must be positive"):
        junction_tension(0, 0, 5e-4, 5e-4)
