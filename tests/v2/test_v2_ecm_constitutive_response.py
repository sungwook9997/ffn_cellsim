"""Tests for HB#1+#2 ECM orientation constitutive response.

Test catalog mirrors
``docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`` §4
(commits 7b1d3cf + c18a5dc + 42d34e7 + 7154342) and the impl-work
Sanity Gate
``docs/v2/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md``
§7 (commit 186fb75 + same amendments). 20 tests; count not
capped (extras allowed if regression class surfaces).

Per Codex review id=1488: tests 3-9 assert
ECMConstitutiveResponseError + specific failure_kind. Per Codex
id=1486: test 20 asserts plain schema ValueError for the
inherited ecm.validate() rejection (NOT wrapped).
"""

from __future__ import annotations

import dataclasses
import inspect
import math

import numpy as np
import pytest

from acs.v2 import dynamics as v2_dynamics
import acs.v2 as v2_pkg
from acs.v2.dynamics import ecm_constitutive_response as ecr
from acs.v2.dynamics.ecm_constitutive_response import (
    ECMConstitutiveResponseError,
    ECMOrientationResponseDiagnostics,
    ECMOrientationResponseResult,
    K_ORIENT_PER_S,
    TAU_ALIGN_RANGE_S,
    TRACTION_REF_NN_PER_UM2,
    step_ecm_orientation_response,
)
from acs.v2.ecm_substrate import ECMSubstrateState


def _ecm(
    nx: int = 4,
    ny: int = 4,
    *,
    spacing_um: float = 1.0,
    origin=(0.0, 0.0),
    stiffness_value: float = 0.0,
    fiber_value: float = 0.0,
    ligand_value: float = 0.0,
    accumulated_value: float = 0.0,
    orientation_value: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 1.0),
):
    """Build a minimal valid ECMSubstrateState for tests.

    `orientation_value` is `(T_00, T_01, T_10, T_11)`. Default
    `(1, 0, 0, 1)` is the identity tensor (isotropic) — valid
    componentwise but not the rank-1 target of any nonzero
    traction.
    """
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = orientation_value[0]
    orientation[..., 0, 1] = orientation_value[1]
    orientation[..., 1, 0] = orientation_value[2]
    orientation[..., 1, 1] = orientation_value[3]
    return ECMSubstrateState(
        origin_um_xy=origin,
        spacing_um=spacing_um,
        stiffness_kpa=np.full((nx, ny), stiffness_value, dtype=np.float64),
        ligand_density=np.full((nx, ny), ligand_value, dtype=np.float64),
        fiber_density=np.full((nx, ny), fiber_value, dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.full(
            (nx, ny), accumulated_value, dtype=np.float64
        ),
        source="synthetic",
    )


def _zero_traction(nx: int = 4, ny: int = 4) -> np.ndarray:
    return np.zeros((nx, ny, 2), dtype=np.float64)


def _uniform_traction(
    nx: int = 4, ny: int = 4, *, magnitude: float = 1.0, direction=(1.0, 0.0)
) -> np.ndarray:
    """Uniform traction across the grid, magnitude × direction."""

    arr = np.zeros((nx, ny, 2), dtype=np.float64)
    arr[..., 0] = magnitude * direction[0]
    arr[..., 1] = magnitude * direction[1]
    return arr


# ---------------------------------------------------------------------------
# Test 1: zero traction is identity
# ---------------------------------------------------------------------------
def test_orientation_response_zero_traction_is_identity():
    ecm = _ecm()
    traction = _zero_traction()
    result = step_ecm_orientation_response(ecm, traction, dt_s=60.0)
    np.testing.assert_array_equal(
        result.updated_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert result.diagnostics.n_nonzero_cells == 0
    assert result.diagnostics.max_convex_weight == 0.0
    assert result.diagnostics.max_orientation_delta_frobenius == 0.0


# ---------------------------------------------------------------------------
# Test 2: dt=0 valid no-op (Y15)
# ---------------------------------------------------------------------------
def test_orientation_response_dt_zero_is_identity_no_op():
    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0)  # nonzero stimulus
    result = step_ecm_orientation_response(ecm, traction, dt_s=0.0)
    np.testing.assert_array_equal(
        result.updated_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert result.diagnostics.max_convex_weight == 0.0


# ---------------------------------------------------------------------------
# Tests 3-5: dt invalid (negative / bool / nonfinite) → dt_invalid
# ---------------------------------------------------------------------------
def test_orientation_response_dt_negative_raises():
    ecm = _ecm()
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(ecm, _zero_traction(), dt_s=-0.1)
    assert exc.value.failure_kind == "dt_invalid"


def test_orientation_response_dt_bool_raises():
    ecm = _ecm()
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(ecm, _zero_traction(), dt_s=True)
    assert exc.value.failure_kind == "dt_invalid"


def test_orientation_response_dt_nonfinite_raises():
    ecm = _ecm()
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(
            ecm, _zero_traction(), dt_s=float("nan")
        )
    assert exc.value.failure_kind == "dt_invalid"


# ---------------------------------------------------------------------------
# Test 6: traction_ref nonpositive → traction_ref_invalid
# ---------------------------------------------------------------------------
def test_orientation_response_traction_ref_nonpositive_raises():
    ecm = _ecm()
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(
            ecm, _zero_traction(), dt_s=60.0,
            traction_ref_nN_per_um2=0.0,
        )
    assert exc.value.failure_kind == "traction_ref_invalid"


# ---------------------------------------------------------------------------
# Test 7: k_orient nonpositive (strict, no zero no-op) → k_orient_invalid
# ---------------------------------------------------------------------------
def test_orientation_response_k_orient_nonpositive_raises():
    ecm = _ecm()
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(
            ecm, _zero_traction(), dt_s=60.0,
            k_orient_per_s=0.0,
        )
    assert exc.value.failure_kind == "k_orient_invalid"


# ---------------------------------------------------------------------------
# Test 8: traction shape mismatch
# ---------------------------------------------------------------------------
def test_orientation_response_traction_shape_mismatch_raises():
    ecm = _ecm(nx=4, ny=4)
    bad_traction = np.zeros((4, 5, 2), dtype=np.float64)  # wrong ny
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(ecm, bad_traction, dt_s=60.0)
    assert exc.value.failure_kind == "traction_shape_mismatch"


# ---------------------------------------------------------------------------
# Test 9: nonfinite traction → non_finite_traction
# ---------------------------------------------------------------------------
def test_orientation_response_traction_nonfinite_raises():
    ecm = _ecm()
    bad = _zero_traction()
    bad[0, 0, 0] = float("inf")
    with pytest.raises(ECMConstitutiveResponseError) as exc:
        step_ecm_orientation_response(ecm, bad, dt_s=60.0)
    assert exc.value.failure_kind == "non_finite_traction"


# ---------------------------------------------------------------------------
# Test 10: no-aliasing of all 5 ECM array fields (Y12)
# ---------------------------------------------------------------------------
def test_orientation_response_no_aliasing_all_arrays():
    ecm = _ecm(stiffness_value=2.0, fiber_value=0.5, ligand_value=0.3)
    traction = _uniform_traction(magnitude=2.0)
    result = step_ecm_orientation_response(ecm, traction, dt_s=60.0)
    new_ecm = result.updated_ecm

    assert not np.shares_memory(new_ecm.stiffness_kpa, ecm.stiffness_kpa)
    assert not np.shares_memory(new_ecm.ligand_density, ecm.ligand_density)
    assert not np.shares_memory(new_ecm.fiber_density, ecm.fiber_density)
    assert not np.shares_memory(
        new_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert not np.shares_memory(
        new_ecm.accumulated_traction_nNs_per_um2,
        ecm.accumulated_traction_nNs_per_um2,
    )


# ---------------------------------------------------------------------------
# Test 11: unchanged fields bytewise equal
# ---------------------------------------------------------------------------
def test_orientation_response_unchanged_fields_bytewise_equal():
    ecm = _ecm(
        stiffness_value=2.0, fiber_value=0.5, ligand_value=0.3,
        accumulated_value=1.5,
    )
    traction = _uniform_traction(magnitude=2.0)
    result = step_ecm_orientation_response(ecm, traction, dt_s=60.0)
    new_ecm = result.updated_ecm

    np.testing.assert_array_equal(new_ecm.stiffness_kpa, ecm.stiffness_kpa)
    np.testing.assert_array_equal(new_ecm.ligand_density, ecm.ligand_density)
    np.testing.assert_array_equal(new_ecm.fiber_density, ecm.fiber_density)
    np.testing.assert_array_equal(
        new_ecm.accumulated_traction_nNs_per_um2,
        ecm.accumulated_traction_nNs_per_um2,
    )


# ---------------------------------------------------------------------------
# Test 12: componentwise bound preserved (|T_ij| <= 1)
# ---------------------------------------------------------------------------
def test_orientation_response_componentwise_bound_preserved():
    ecm = _ecm()
    traction = _uniform_traction(magnitude=10.0)  # high stimulus
    result = step_ecm_orientation_response(ecm, traction, dt_s=3600.0 * 10)
    assert np.abs(result.updated_ecm.orientation_tensor).max() <= 1.0 + 1e-12
    # validate() called inside step; if we got here it passed.


# ---------------------------------------------------------------------------
# Test 13: symmetry preserved
# ---------------------------------------------------------------------------
def test_orientation_response_symmetry_preserved():
    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0)
    result = step_ecm_orientation_response(ecm, traction, dt_s=60.0)
    T = result.updated_ecm.orientation_tensor
    np.testing.assert_allclose(T, np.swapaxes(T, -1, -2), atol=1e-12)


# ---------------------------------------------------------------------------
# Test 14: distance-to-target monotone per step under fixed target (Y2)
# ---------------------------------------------------------------------------
def test_orientation_response_distance_to_target_monotone_per_step():
    ecm = _ecm()  # identity orientation
    traction = _uniform_traction(magnitude=2.0, direction=(1.0, 0.0))
    # T_target = (1,0;0,0) since n=(1,0); T_old = (1,0;0,1)
    # dist_old = ‖(0,0;0,1)‖_F = 1
    T_target = np.zeros((2, 2))
    T_target[0, 0] = 1.0
    T_old = ecm.orientation_tensor[0, 0]
    dist_old = np.linalg.norm(T_old - T_target, ord="fro")

    result = step_ecm_orientation_response(ecm, traction, dt_s=3600.0)
    T_new = result.updated_ecm.orientation_tensor[0, 0]
    dist_new = np.linalg.norm(T_new - T_target, ord="fro")
    assert dist_new <= dist_old + 1e-12


# ---------------------------------------------------------------------------
# Test 15: uniform traction converges to target (Y2 corollary)
# ---------------------------------------------------------------------------
def test_orientation_response_uniform_traction_converges_to_target():
    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0, direction=(1.0, 0.0))
    state = ecm
    for _ in range(2000):
        state = step_ecm_orientation_response(
            state, traction, dt_s=3600.0
        ).updated_ecm
    # T_target = n⊗n with n=(1,0) → ((1,0),(0,0))
    T_final = state.orientation_tensor[0, 0]
    np.testing.assert_allclose(T_final, np.array([[1.0, 0.0], [0.0, 0.0]]), atol=1e-6)


# ---------------------------------------------------------------------------
# Test 16: typed dataclass + __getitem__ raises TypeError (B1 sister-pattern)
# ---------------------------------------------------------------------------
def test_orientation_response_diagnostics_typed_dataclass():
    ecm = _ecm()
    result = step_ecm_orientation_response(ecm, _uniform_traction(), dt_s=60.0)
    assert isinstance(result.diagnostics, ECMOrientationResponseDiagnostics)
    with pytest.raises(TypeError):
        _ = result.diagnostics["max_convex_weight"]


# ---------------------------------------------------------------------------
# Test 17: diagnostics records input parameters (Y13 sweep traceback)
# ---------------------------------------------------------------------------
def test_orientation_response_diagnostics_records_input_params():
    ecm = _ecm()
    custom_traction_ref = 2.5
    custom_k_orient = 5e-4
    result = step_ecm_orientation_response(
        ecm, _uniform_traction(magnitude=3.0), dt_s=60.0,
        traction_ref_nN_per_um2=custom_traction_ref,
        k_orient_per_s=custom_k_orient,
    )
    assert result.diagnostics.traction_ref_nN_per_um2 == custom_traction_ref
    assert result.diagnostics.k_orient_per_s == custom_k_orient


# ---------------------------------------------------------------------------
# Test 18: 7 symbols exported through both __init__ surfaces (Codex id=1488)
# ---------------------------------------------------------------------------
def test_orientation_response_exports_through_both_init():
    expected_symbols = {
        "ECMConstitutiveResponseError",
        "ECMOrientationResponseDiagnostics",
        "ECMOrientationResponseResult",
        "K_ORIENT_PER_S",
        "TAU_ALIGN_RANGE_S",
        "TRACTION_REF_NN_PER_UM2",
        "step_ecm_orientation_response",
    }

    # Module-level attribute presence + identity match.
    for sym in expected_symbols:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym!r}"
        assert hasattr(v2_dynamics, sym), f"acs.v2.dynamics missing {sym!r}"
        assert getattr(v2_pkg, sym) is getattr(ecr, sym), (
            f"acs.v2.{sym} is not acs.v2.dynamics.ecm_constitutive_response.{sym}"
        )
        assert getattr(v2_dynamics, sym) is getattr(ecr, sym), (
            f"acs.v2.dynamics.{sym} is not "
            f"acs.v2.dynamics.ecm_constitutive_response.{sym}"
        )
        assert sym in v2_pkg.__all__, f"acs.v2.__all__ missing {sym!r}"
        assert sym in v2_dynamics.__all__, (
            f"acs.v2.dynamics.__all__ missing {sym!r}"
        )


# ---------------------------------------------------------------------------
# Test 19: Y12-guard static check — uses current schema, validates BOTH ways
# ---------------------------------------------------------------------------
def test_orientation_response_uses_current_schema_not_stale_naming():
    """Per Codex id=1486 + id=1488: implementation must use
    spacing_um (not stale dx_um), must NOT pass grid_shape to
    ECMSubstrateState constructor (it's a derived property), and
    must call BOTH ecm.validate() at function entry AND
    updated_ecm.validate() before return."""

    src = inspect.getsource(step_ecm_orientation_response)

    # Stale field names absent.
    assert "dx_um" not in src, (
        "step_ecm_orientation_response references stale dx_um field"
    )

    # No grid_shape arg to ECMSubstrateState constructor (derived
    # property only).
    # Allow grid_shape attribute reads (e.g., ecm.grid_shape) but not
    # construction kwarg.
    constructor_call_ok = (
        "ECMSubstrateState(\n" in src or "ECMSubstrateState(" in src
    )
    assert constructor_call_ok
    # Grep ensures no `grid_shape=` kwarg in any constructor call.
    assert "grid_shape=" not in src, (
        "step_ecm_orientation_response passes grid_shape to "
        "ECMSubstrateState constructor"
    )

    # Both validation calls present.
    assert "ecm.validate()" in src, (
        "step_ecm_orientation_response missing ecm.validate() at entry"
    )
    assert "updated_ecm.validate()" in src, (
        "step_ecm_orientation_response missing updated_ecm.validate() before return"
    )


# ---------------------------------------------------------------------------
# Test 20: invalid input ECM rejected before update (Codex id=1486)
# ---------------------------------------------------------------------------
def test_orientation_response_invalid_input_ecm_rejected_before_update():
    """Per Codex id=1486 silent-heal-path closure: passing an
    ECMSubstrateState whose orientation_tensor violates validate()
    must raise plain ValueError from ecm.validate() at function
    entry, BEFORE any algebra runs. Per Codex id=1488: inherited
    schema ValueError is NOT wrapped in
    ECMConstitutiveResponseError — schema-level errors stay
    schema-level."""

    # Build an ECM with healable invalid orientation (component > 1
    # which the convex update could "heal" toward T_target under
    # high traction × dt). We bypass the validate() in the
    # constructor by mutating the array post-construction.
    ecm = _ecm()
    ecm.orientation_tensor[0, 0, 0, 0] = 5.0  # violates |T_ij| <= 1
    traction = _uniform_traction(magnitude=10.0)

    with pytest.raises(ValueError) as exc:
        step_ecm_orientation_response(ecm, traction, dt_s=3600.0 * 100)

    # Per Codex id=1488: not wrapped in ECMConstitutiveResponseError.
    assert not isinstance(exc.value, ECMConstitutiveResponseError), (
        "schema ValueError must NOT be wrapped in ECMConstitutiveResponseError; "
        "schema-level errors stay schema-level per Codex id=1488"
    )
