"""Tests for HB#5 Lyapunov-like ECM orientation metric.

Test catalog mirrors
``docs/v2_hard_blocker_5_lyapunov_metric_locked.md`` §4 (commit
340c353) and the impl-work Sanity Gate
``docs/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`` §8
(commit 1c9bb0a). 15 active tests after locked §4 14/15 merge;
count not capped (extras allowed if regression class surfaces).

Per Codex `id=1518` Sanity Gate PASS, this catalog covers:
- Y1 active-cells-only V + passive diagnostic separate (test 12)
- Y2 positive-semidefinite (no PSD claim — implicit)
- Y3 two-channel evidence: contraction (test 5) + boundedness (test 6)
- Y4 ECM-half scope (implicit; tests do not assert FA loop closure)
- Y7 branch-defined active mask + ecm.validate() at entry (tests 8-10)
- Y8 9 → 4.5 derivation chain (test 4 dual + test 16)
- Y10 dual assertion in test 4
- Y11 AST-walk meta-test scoped to function body (test 14)
"""

from __future__ import annotations

import ast
import inspect
import math
import textwrap

import numpy as np
import pytest

from acs.v2 import dynamics as v2_dynamics
import acs.v2 as v2_pkg
from acs.v2.dynamics import ecm_lyapunov_metric as ecm_lm
from acs.v2.dynamics.ecm_constitutive_response import (
    step_ecm_orientation_response,
)
from acs.v2.dynamics.ecm_lyapunov_metric import (
    ECMOrientationLyapunovMetricDiagnostics,
    ECMOrientationLyapunovMetricResult,
    MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL,
    MAX_SQ_FROBENIUS_DIFF_PER_CELL,
    compute_ecm_orientation_lyapunov_metric,
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
    """Build a minimal valid ECMSubstrateState for tests."""
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


def test_lyapunov_metric_zero_traction_returns_zero():
    """Test 1: zero traction returns V_active = 0; passive diagnostic
    = full grid (per locked §4 + Sanity Gate §2.1)."""

    ecm = _ecm()
    traction = _zero_traction()
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    assert result.v_active_um2 == 0.0
    assert result.diagnostics.n_active_cells == 0
    expected_passive = 0.5 * 1.0 * 16 * 2.0
    assert math.isclose(
        result.diagnostics.passive_orientation_magnitude_um2,
        expected_passive,
    )


def test_lyapunov_metric_no_active_cells_no_passive_penalty_in_v():
    """Test 2: V_active does NOT include passive cells (Y1 forward
    guard per locked §4 + Sanity Gate §6.1)."""

    ecm = _ecm(orientation_value=(1.0, 0.0, 0.0, -1.0))
    traction = _zero_traction()
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    assert result.v_active_um2 == 0.0
    assert result.diagnostics.passive_orientation_magnitude_um2 > 0.0


def test_lyapunov_metric_at_target_returns_zero():
    """Test 3: at-target orientation returns V_active = 0 (per locked §4
    + Sanity Gate §2.3)."""

    nx, ny = 4, 4
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    traction = _uniform_traction(nx, ny, magnitude=2.0, direction=(1.0, 0.0))
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    assert result.v_active_um2 == 0.0


def test_lyapunov_metric_v_active_within_bound_4_5_times_active_area():
    """Test 4: adversarial-saturation dual assertion (Y10): per-cell
    raw squared Frobenius diff = 9.0 AND v_active = 4.5 · active_area
    (per locked §4 + Y10 + Sanity Gate §2.4 + §6.2)."""

    nx, ny = 4, 4
    orientation = np.full((nx, ny, 2, 2), -1.0, dtype=np.float64)
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    sqrt_2_inv = 1.0 / math.sqrt(2.0)
    traction = _uniform_traction(
        nx, ny, magnitude=1.0, direction=(sqrt_2_inv, sqrt_2_inv)
    )
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)

    n_xy = traction[0, 0] / np.linalg.norm(traction[0, 0])
    T_target = np.outer(n_xy, n_xy)
    diff = orientation[0, 0] - T_target
    raw_sq_frobenius = float(np.sum(diff ** 2))
    assert math.isclose(raw_sq_frobenius, MAX_SQ_FROBENIUS_DIFF_PER_CELL, abs_tol=1e-12)

    expected_v = MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL * (nx * ny * 1.0)
    assert math.isclose(result.v_active_um2, expected_v, abs_tol=1e-9)
    assert math.isclose(
        result.v_active_um2,
        result.diagnostics.v_active_max_bound_um2,
        abs_tol=1e-9,
    )


def test_lyapunov_metric_fixed_target_contraction_under_hb12_step():
    """Test 5: fixed-target per-step contraction under HB#1+#2 step
    (Y3 channel A per locked §4 + Sanity Gate §3.2)."""

    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0, direction=(1.0, 0.0))
    v_before = compute_ecm_orientation_lyapunov_metric(ecm, traction).v_active_um2
    ecm_after = step_ecm_orientation_response(ecm, traction, dt_s=3600.0).updated_ecm
    v_after = compute_ecm_orientation_lyapunov_metric(ecm_after, traction).v_active_um2
    assert v_after <= v_before + 1e-12


def test_lyapunov_metric_varying_target_boundedness():
    """Test 6: varying-target boundedness (Y3 channel B per locked §4
    + Sanity Gate §3.3)."""

    ecm = _ecm()
    state = ecm
    for k in range(20):
        theta = 2 * math.pi * k / 20
        traction = _uniform_traction(
            magnitude=2.0, direction=(math.cos(theta), math.sin(theta))
        )
        result = compute_ecm_orientation_lyapunov_metric(state, traction)
        assert result.v_active_um2 >= 0.0
        assert result.v_active_um2 <= result.diagnostics.v_active_max_bound_um2 + 1e-9
        state = step_ecm_orientation_response(state, traction, dt_s=600.0).updated_ecm


def test_lyapunov_metric_pure_no_mutation_of_ecm():
    """Test 7: pure read-only no mutation (per locked §4 + Sanity Gate
    §3.1)."""

    ecm = _ecm(stiffness_value=2.0, fiber_value=0.5, ligand_value=0.3, accumulated_value=1.5)
    traction = _uniform_traction(magnitude=2.0)
    stiffness_before = ecm.stiffness_kpa.copy()
    fiber_before = ecm.fiber_density.copy()
    ligand_before = ecm.ligand_density.copy()
    orientation_before = ecm.orientation_tensor.copy()
    accumulated_before = ecm.accumulated_traction_nNs_per_um2.copy()

    _ = compute_ecm_orientation_lyapunov_metric(ecm, traction)

    np.testing.assert_array_equal(ecm.stiffness_kpa, stiffness_before)
    np.testing.assert_array_equal(ecm.fiber_density, fiber_before)
    np.testing.assert_array_equal(ecm.ligand_density, ligand_before)
    np.testing.assert_array_equal(ecm.orientation_tensor, orientation_before)
    np.testing.assert_array_equal(
        ecm.accumulated_traction_nNs_per_um2, accumulated_before
    )


def test_lyapunov_metric_traction_nonfinite_raises():
    """Test 8: nonfinite traction raises plain ValueError (per locked §5
    no error class; Sanity Gate §2.7)."""

    ecm = _ecm()
    bad = _zero_traction()
    bad[0, 0, 0] = float("inf")
    with pytest.raises(ValueError):
        compute_ecm_orientation_lyapunov_metric(ecm, bad)


def test_lyapunov_metric_traction_shape_mismatch_raises():
    """Test 9: traction shape mismatch raises plain ValueError."""

    ecm = _ecm(nx=4, ny=4)
    bad = np.zeros((4, 5, 2), dtype=np.float64)
    with pytest.raises(ValueError):
        compute_ecm_orientation_lyapunov_metric(ecm, bad)


def test_lyapunov_metric_invalid_ecm_raises_at_validate():
    """Test 10: invalid input ECM raises plain ValueError at
    ecm.validate() entry (HB#1+#2 id=1486 sister-precedent;
    Sanity Gate §2.6)."""

    ecm = _ecm()
    ecm.orientation_tensor[0, 0, 0, 0] = 5.0
    traction = _uniform_traction(magnitude=2.0)
    with pytest.raises(ValueError):
        compute_ecm_orientation_lyapunov_metric(ecm, traction)


def test_lyapunov_metric_diagnostics_typed_dataclass():
    """Test 11: typed dataclass + __getitem__ raises TypeError
    (B1 sister-pattern meta-test per locked §4)."""

    ecm = _ecm()
    result = compute_ecm_orientation_lyapunov_metric(ecm, _uniform_traction())
    assert isinstance(result.diagnostics, ECMOrientationLyapunovMetricDiagnostics)
    with pytest.raises(TypeError):
        _ = result.diagnostics["v_active_max_bound_um2"]


def test_lyapunov_metric_passive_diagnostic_separate_from_v_active():
    """Test 12: passive diagnostic separate from V_active (Y1 forward
    guard per locked §4 + Sanity Gate §6.1)."""

    nx, ny = 4, 4
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    orientation[0:2, :, 0, 1] = 1.0
    orientation[0:2, :, 1, 0] = 1.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    traction = np.zeros((nx, ny, 2), dtype=np.float64)
    traction[2:4, :, 0] = 2.0

    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    assert result.diagnostics.n_active_cells == 2 * ny
    assert result.diagnostics.passive_orientation_magnitude_um2 > 0.0

    orientation_small_passive = orientation.copy()
    orientation_small_passive[0:2, :, 0, 1] = 0.0
    orientation_small_passive[0:2, :, 1, 0] = 0.0
    ecm_small_passive = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation_small_passive,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    result_small_passive = compute_ecm_orientation_lyapunov_metric(
        ecm_small_passive, traction
    )

    assert math.isclose(
        result.v_active_um2,
        result_small_passive.v_active_um2,
        abs_tol=1e-12,
    )
    assert (
        result.diagnostics.passive_orientation_magnitude_um2
        > result_small_passive.diagnostics.passive_orientation_magnitude_um2
    )


def test_lyapunov_metric_exports_through_both_init():
    """Test 13: 5 symbols exported through both __init__ surfaces
    (per locked §4 + Sanity Gate §6.7 sister-gate-mirror API surface)."""

    expected_symbols = {
        "ECMOrientationLyapunovMetricDiagnostics",
        "ECMOrientationLyapunovMetricResult",
        "MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL",
        "MAX_SQ_FROBENIUS_DIFF_PER_CELL",
        "compute_ecm_orientation_lyapunov_metric",
    }

    for sym in expected_symbols:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym!r}"
        assert hasattr(v2_dynamics, sym), f"acs.v2.dynamics missing {sym!r}"
        assert getattr(v2_pkg, sym) is getattr(ecm_lm, sym), (
            f"acs.v2.{sym} is not acs.v2.dynamics.ecm_lyapunov_metric.{sym}"
        )
        assert getattr(v2_dynamics, sym) is getattr(ecm_lm, sym), (
            f"acs.v2.dynamics.{sym} is not "
            f"acs.v2.dynamics.ecm_lyapunov_metric.{sym}"
        )
        assert sym in v2_pkg.__all__, f"acs.v2.__all__ missing {sym!r}"
        assert sym in v2_dynamics.__all__, (
            f"acs.v2.dynamics.__all__ missing {sym!r}"
        )


def test_lyapunov_metric_does_not_reference_other_ecm_fields():
    """Test 14 (locked §4 14/15 merged): AST-walk meta-test for 4
    forbidden ECM field names (Y11 + Codex id=1506 review focus 5).
    Walks textwrap.dedent(inspect.getsource(...)) of the function
    body; asserts none of `accumulated_traction_nNs_per_um2`,
    `stiffness_kpa`, `fiber_density`, `ligand_density` appear as
    `ast.Attribute.attr` or `ast.Name.id`."""

    src = textwrap.dedent(inspect.getsource(compute_ecm_orientation_lyapunov_metric))
    tree = ast.parse(src)
    referenced = (
        {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
        | {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    )
    forbidden = {
        "accumulated_traction_nNs_per_um2",
        "stiffness_kpa",
        "fiber_density",
        "ligand_density",
    }
    overlap = forbidden & referenced
    assert not overlap, (
        f"compute_ecm_orientation_lyapunov_metric must not reference "
        f"{overlap}; HB#5 is orientation-only metric per locked §1 Y7 "
        f"+ id=1506 review focus"
    )


def test_lyapunov_metric_active_area_consistent_with_n_active():
    """Test 15 (locked §4 test 16): active_area consistent with
    n_active_cells (sanity)."""

    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0)
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    expected_active_area = (
        result.diagnostics.n_active_cells * result.diagnostics.cell_area_um2
    )
    assert math.isclose(
        result.diagnostics.active_area_um2, expected_active_area, abs_tol=1e-12
    )


def test_lyapunov_metric_max_bound_consistent_with_active_area():
    """Test 16 (locked §4 test 17): Y8 derivation chain check —
    v_active_max_bound_um2 == 4.5 · active_area_um2."""

    ecm = _ecm()
    traction = _uniform_traction(magnitude=2.0)
    result = compute_ecm_orientation_lyapunov_metric(ecm, traction)
    expected_bound = (
        MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL
        * result.diagnostics.active_area_um2
    )
    assert math.isclose(
        result.diagnostics.v_active_max_bound_um2,
        expected_bound,
        abs_tol=1e-12,
    )
    assert MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL == 0.5 * MAX_SQ_FROBENIUS_DIFF_PER_CELL
    assert MAX_SQ_FROBENIUS_DIFF_PER_CELL == 9.0
    assert MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL == 4.5
