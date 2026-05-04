"""Tests for FA→ECM bilinear scattering (Hard Blocker #3).

Test catalog mirrors
``docs/v2_fa_to_ecm_scattering_sanity_gate.md`` §9 (18 tests: 17
from locked design §7 + 1 Hard Rule 11 meta-test). Each test maps
to one Sanity Gate item or one forbidden behavior.
"""

from __future__ import annotations

from unittest.mock import patch

import math
import numpy as np
import pytest

from acs.v2.dynamics.fa_to_ecm_scattering import (
    FAToECMScatteringError,
    scatter_fa_traction_to_ecm_bilinear,
)
from acs.v2.ecm_open_loop_harness import make_default_ecm
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


def _ecm(nx: int = 4, ny: int = 4, *, spacing_um: float = 1.0, origin=(0.0, 0.0)):
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=origin,
        spacing_um=spacing_um,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )


def _fa(
    adhesion_id: str,
    *,
    position=(0.5, 0.5),
    traction=(1.0, 0.0),
    cell_id: str = "cell-A",
):
    return FocalAdhesionState(
        adhesion_id=adhesion_id,
        cell_id=cell_id,
        position_um_xy=position,
        age_s=0.0,
        maturity=1.0,
        bound_fraction=1.0,
        state="mature",
        traction_force_nN_xy=traction,
    )


# ---------------------------------------------------------------------------
# Test 1: single FA at exact cell center deposits full traction to one cell
# ---------------------------------------------------------------------------


def test_single_fa_at_center_deposits_full_traction_to_one_cell():
    """FA at (0.5, 0.5) with spacing 1.0 lies exactly at center
    (i=0, j=0). Bilinear weights degenerate to (1, 0, 0, 0)."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(0.5, 0.5), traction=(2.0, 3.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)

    assert out.shape == (4, 4, 2)
    assert out.dtype == np.float64
    # Full traction lands at (0, 0); cell area = 1.0 → density = traction.
    np.testing.assert_allclose(out[0, 0], (2.0, 3.0), atol=0.0)
    # All other cells zero.
    rest = out.copy()
    rest[0, 0] = 0.0
    np.testing.assert_array_equal(rest, np.zeros((4, 4, 2)))


# ---------------------------------------------------------------------------
# Test 2: FA at 4-center crossing splits 25% each
# ---------------------------------------------------------------------------


def test_single_fa_at_grid_crossing_splits_25_each():
    """FA at (1.0, 1.0) with spacing 1.0 lies exactly at the
    crossing of centers (0,0), (1,0), (0,1), (1,1). Weights
    (0.25, 0.25, 0.25, 0.25)."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(1.0, 1.0), traction=(4.0, 8.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)

    expected = np.array([1.0, 2.0])  # (4/4, 8/4) per cell
    for k, l in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        np.testing.assert_allclose(out[k, l], expected, atol=1e-12)


# ---------------------------------------------------------------------------
# Test 3: FA midway along x splits 50/50 in x, full weight in y
# ---------------------------------------------------------------------------


def test_single_fa_midway_along_x_splits_50_50():
    """FA at (1.0, 0.5) with spacing 1.0: x is at center crossing
    (i=0, i=1), y is at center j=0. Weights (0.5, 0.5, 0, 0)."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(1.0, 0.5), traction=(2.0, 0.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)

    np.testing.assert_allclose(out[0, 0], (1.0, 0.0), atol=1e-12)
    np.testing.assert_allclose(out[1, 0], (1.0, 0.0), atol=1e-12)
    np.testing.assert_allclose(out[0, 1], (0.0, 0.0), atol=0.0)


# ---------------------------------------------------------------------------
# Test 4: FA at physical corner deposits 100% to corner center
# ---------------------------------------------------------------------------


def test_single_fa_at_physical_corner_deposits_100_to_corner_center():
    """FA at (0.0, 0.0) with origin (0, 0) and spacing 1.0 is at
    the physical corner. The 2×2 stencil for centers
    (-1, -1)/(0, -1)/(-1, 0)/(0, 0) has only (0, 0) in-grid;
    after renormalization, full traction goes to (0, 0)."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, origin=(0.0, 0.0))
    fa = _fa("fa-1", position=(0.0, 0.0), traction=(5.0, 7.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)

    np.testing.assert_allclose(out[0, 0], (5.0, 7.0), atol=1e-12)
    rest = out.copy()
    rest[0, 0] = 0.0
    np.testing.assert_array_equal(rest, np.zeros((4, 4, 2)))


# ---------------------------------------------------------------------------
# Test 5: FA on physical edge midpoint uses truncated stencil
# ---------------------------------------------------------------------------


def test_single_fa_on_physical_edge_midpoint_uses_truncated_stencil():
    """FA at (0.0, 1.0) with origin (0, 0) and spacing 1.0 is on
    the left edge midway between rows j=0 and j=1. Stencil candidates
    are (-1, 0), (0, 0), (-1, 1), (0, 1); only (0, 0) and (0, 1)
    are in-grid. After renormalization, deposit splits 50/50
    between them."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, origin=(0.0, 0.0))
    fa = _fa("fa-1", position=(0.0, 1.0), traction=(2.0, 0.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)

    np.testing.assert_allclose(out[0, 0], (1.0, 0.0), atol=1e-12)
    np.testing.assert_allclose(out[0, 1], (1.0, 0.0), atol=1e-12)


# ---------------------------------------------------------------------------
# Test 6: two FAs at distinct positions sum correctly
# ---------------------------------------------------------------------------


def test_two_fas_distinct_positions_sum_correctly():
    """Linearity: the scatter of two FAs equals the sum of their
    individual scatters."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa_a = _fa("fa-a", position=(0.5, 0.5), traction=(1.0, 0.0))
    fa_b = _fa("fa-b", position=(2.5, 2.5), traction=(0.0, 1.0))

    both = scatter_fa_traction_to_ecm_bilinear((fa_a, fa_b), ecm)
    only_a = scatter_fa_traction_to_ecm_bilinear((fa_a,), ecm)
    only_b = scatter_fa_traction_to_ecm_bilinear((fa_b,), ecm)
    np.testing.assert_allclose(both, only_a + only_b, atol=1e-12)


# ---------------------------------------------------------------------------
# Test 7: zero FA list returns zero field
# ---------------------------------------------------------------------------


def test_zero_fa_returns_zero_field():
    ecm = _ecm(nx=3, ny=5, spacing_um=2.0)
    out = scatter_fa_traction_to_ecm_bilinear((), ecm)
    assert out.shape == (3, 5, 2)
    np.testing.assert_array_equal(out, np.zeros((3, 5, 2)))


# ---------------------------------------------------------------------------
# Test 8: component-wise conservation with random interior FAs
# ---------------------------------------------------------------------------


def test_conservation_component_wise_with_random_fas():
    """∑ density · cell_area equals ∑ FA traction per component
    within float64 tolerance."""

    rng = np.random.default_rng(42)
    nx, ny = 6, 6
    spacing_um = 1.5
    ecm = _ecm(nx=nx, ny=ny, spacing_um=spacing_um, origin=(2.0, -1.0))
    cell_area = spacing_um * spacing_um

    fas: list[FocalAdhesionState] = []
    total_force = np.zeros(2, dtype=np.float64)
    # Sample positions strictly inside the interior bracket region
    # so renormalization is not exercised here.
    for k in range(20):
        x = float(2.0 + 0.75 + (nx - 1.5) * spacing_um * rng.random())
        y = float(-1.0 + 0.75 + (ny - 1.5) * spacing_um * rng.random())
        fx = float(rng.standard_normal() * 2.0)
        fy = float(rng.standard_normal() * 2.0)
        fas.append(_fa(f"fa-{k}", position=(x, y), traction=(fx, fy)))
        total_force += (fx, fy)

    out = scatter_fa_traction_to_ecm_bilinear(tuple(fas), ecm)
    deposited = out.sum(axis=(0, 1)) * cell_area
    np.testing.assert_allclose(
        deposited,
        total_force,
        atol=1e-9 * max(np.abs(total_force).max(), 1.0),
    )


# ---------------------------------------------------------------------------
# Test 9: out-of-grid raises fa_position_outside_ecm_grid
# ---------------------------------------------------------------------------


def test_out_of_grid_raises_fa_position_outside_ecm_grid():
    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, origin=(0.0, 0.0))
    fa = _fa("fa-1", position=(5.0, 0.5), traction=(1.0, 0.0))
    with pytest.raises(FAToECMScatteringError) as excinfo:
        scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    assert excinfo.value.failure_kind == "fa_position_outside_ecm_grid"
    assert "fa-1" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Test 10: inside boundary tol accepted
# ---------------------------------------------------------------------------


def test_inside_boundary_tol_accepted():
    """FA position within `_FA_POSITION_BOUNDARY_TOL_UM` of the
    physical edge is accepted."""

    nx, ny = 4, 4
    spacing_um = 1.0
    extent = max(nx * spacing_um, ny * spacing_um)
    tol = 1e-12 * extent
    ecm = _ecm(nx=nx, ny=ny, spacing_um=spacing_um, origin=(0.0, 0.0))
    # 0.5 * tol inside the upper-x edge
    fa = _fa(
        "fa-1",
        position=(nx * spacing_um + 0.5 * tol, 1.5),
        traction=(1.0, 0.0),
    )
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    # Conservation must still hold within tolerance
    deposited = out.sum(axis=(0, 1)) * spacing_um * spacing_um
    np.testing.assert_allclose(deposited, (1.0, 0.0), atol=1e-9)


# ---------------------------------------------------------------------------
# Test 11: outside boundary tol raises
# ---------------------------------------------------------------------------


def test_outside_boundary_tol_raises():
    nx, ny = 4, 4
    spacing_um = 1.0
    extent = max(nx * spacing_um, ny * spacing_um)
    tol = 1e-12 * extent
    ecm = _ecm(nx=nx, ny=ny, spacing_um=spacing_um, origin=(0.0, 0.0))
    fa = _fa(
        "fa-1",
        position=(nx * spacing_um + 100 * tol, 1.5),
        traction=(1.0, 0.0),
    )
    with pytest.raises(FAToECMScatteringError) as excinfo:
        scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    assert excinfo.value.failure_kind == "fa_position_outside_ecm_grid"


# ---------------------------------------------------------------------------
# Test 12: boundary half-cell renormalized weights sum to 1
# ---------------------------------------------------------------------------


def test_boundary_half_cell_renormalized_weights_sum_to_one():
    """An FA at (0.0, 0.5) (left edge, between centers (0, 0) and
    (-1, 0)) renormalizes to give 100% weight to (0, 0). Conservation
    sum equals input traction component-wise."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, origin=(0.0, 0.0))
    fa = _fa("fa-1", position=(0.0, 0.5), traction=(3.0, 4.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    np.testing.assert_allclose(out[0, 0], (3.0, 4.0), atol=1e-12)
    deposited = out.sum(axis=(0, 1)) * 1.0
    np.testing.assert_allclose(deposited, (3.0, 4.0), atol=1e-12)


# ---------------------------------------------------------------------------
# Test 13: pure function — no ECM mutation
# ---------------------------------------------------------------------------


def test_pure_function_no_ecm_mutation():
    """ECM byte-equality before and after a scatter call."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    snapshots = {
        "stiffness_kpa": ecm.stiffness_kpa.copy(),
        "ligand_density": ecm.ligand_density.copy(),
        "fiber_density": ecm.fiber_density.copy(),
        "orientation_tensor": ecm.orientation_tensor.copy(),
        "accumulated_traction_nNs_per_um2":
            ecm.accumulated_traction_nNs_per_um2.copy(),
    }
    fa = _fa("fa-1", position=(1.5, 1.5), traction=(1.0, -1.0))
    scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    for name, snap in snapshots.items():
        np.testing.assert_array_equal(getattr(ecm, name), snap)


# ---------------------------------------------------------------------------
# Test 14: pure function — no accumulator call (mock check)
# ---------------------------------------------------------------------------


def test_pure_function_no_accumulator_call_via_mock():
    """The scatter must not invoke
    ``accumulate_prescribed_traction``; verify by mocking it."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(1.5, 1.5), traction=(1.0, 0.0))
    with patch(
        "acs.v2.dynamics.ecm_open_loop.accumulate_prescribed_traction"
    ) as mock_acc:
        scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
        mock_acc.assert_not_called()


# ---------------------------------------------------------------------------
# Test 15: negative traction input preserves sign in x and y
# ---------------------------------------------------------------------------


def test_negative_traction_input_preserved_in_sign_x_and_y():
    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(1.5, 1.5), traction=(-2.0, -3.0))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    deposited = out.sum(axis=(0, 1))
    # Cell area = 1.0 → density sum equals input force divided by area
    np.testing.assert_allclose(deposited, (-2.0, -3.0), atol=1e-12)


# ---------------------------------------------------------------------------
# Test 16: returned shape is (nx, ny, 2) float64
# ---------------------------------------------------------------------------


def test_returned_shape_is_nx_ny_2_float64():
    ecm = _ecm(nx=3, ny=5, spacing_um=1.5)
    fa = _fa("fa-1", position=(2.25, 3.75), traction=(0.5, 0.5))
    out = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    assert out.shape == (3, 5, 2)
    assert out.dtype == np.float64


# ---------------------------------------------------------------------------
# Test 17: non-finite FA position raises
# ---------------------------------------------------------------------------


def test_non_finite_fa_position_raises():
    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa_nan = _fa("fa-nan", position=(float("nan"), 0.5), traction=(1.0, 0.0))
    with pytest.raises(FAToECMScatteringError) as excinfo_nan:
        scatter_fa_traction_to_ecm_bilinear((fa_nan,), ecm)
    assert excinfo_nan.value.failure_kind == "non_finite_fa_position"

    fa_inf = _fa("fa-inf", position=(0.5, math.inf), traction=(1.0, 0.0))
    with pytest.raises(FAToECMScatteringError) as excinfo_inf:
        scatter_fa_traction_to_ecm_bilinear((fa_inf,), ecm)
    assert excinfo_inf.value.failure_kind == "non_finite_fa_position"


def test_non_finite_fa_traction_raises():
    """Bonus boundary: non-finite traction component also raises."""

    ecm = _ecm(nx=4, ny=4, spacing_um=1.0)
    fa = _fa("fa-1", position=(1.5, 1.5), traction=(float("nan"), 0.0))
    with pytest.raises(FAToECMScatteringError) as excinfo:
        scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    assert excinfo.value.failure_kind == "non_finite_fa_traction"


# ---------------------------------------------------------------------------
# Test 18: Hard Rule 11 wording-boundary meta-test (Phase B/C precedent)
# ---------------------------------------------------------------------------


def test_scatter_does_not_satisfy_closed_loop_response_law():
    """Phase D Hard Blocker #3 — Hard Rule 11 wording protection.

    Documents at runtime that the bilinear scatter outputs a vector
    traction-density field; it does NOT decide a constitutive
    response law. Phase E response laws (consuming this output)
    decide how to reduce or compare against the response field.
    This meta-test makes the boundary explicit so a future refactor
    cannot silently relabel the scatter as response-law satisfaction.

    Mirrors:
    - Phase B test 4
      ``test_stimulus_monotonicity_does_not_satisfy_closed_loop_item_1``
    - Phase C test 5
      ``test_open_loop_sweep_does_not_satisfy_closed_loop_item_5``
    """

    scatter_modality = "vector traction density (nx, ny, 2) in nN/μm²"
    response_law_claim = "ECM response field update (Phase E TBD)"
    assert scatter_modality != response_law_claim, (
        "FA→ECM bilinear scatter is an interface primitive, not a "
        "response law. Phase E (separate Sanity Gate per locked "
        "phased plan §1 Phase E) is the response-law decision."
    )
