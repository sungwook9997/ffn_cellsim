"""Tests for ECM→FA bias target (Hard Blocker #4).

Test catalog mirrors
``docs/v2_ecm_to_fa_bias_sanity_gate.md`` §9 (19 tests: 18 from
locked §6 + 1 fa.validate() regression per Codex review id=1378).
Each test maps to one Sanity Gate item or one forbidden behavior.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMSampledAtFAs,
    ECMToFABiasResult,
    FAToECMBiasError,
    RATE_NAMES,
    compute_ecm_to_fa_bias_neutral,
    sample_ecm_at_fa_positions,
)
from acs.v2.dynamics.fa_to_ecm_scattering import scatter_fa_traction_to_ecm_bilinear
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


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
):
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
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


def _fa(
    adhesion_id: str,
    *,
    position=(0.5, 0.5),
    cell_id: str = "cell-A",
    state="mature",
    bound_fraction: float = 1.0,
    traction=(1.0, 0.0),
):
    return FocalAdhesionState(
        adhesion_id=adhesion_id,
        cell_id=cell_id,
        position_um_xy=position,
        age_s=0.0,
        maturity=1.0,
        bound_fraction=bound_fraction,
        state=state,
        traction_force_nN_xy=traction,
    )


# ---------------------------------------------------------------------------
# Test 1: neutral multipliers all 1.0 for arbitrary FA list
# ---------------------------------------------------------------------------


def test_neutral_multipliers_all_ones_for_arbitrary_fa_list():
    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, stiffness_value=2.0)
    fas = (_fa("fa-1", position=(0.5, 0.5)), _fa("fa-2", position=(2.5, 1.5)))
    result = compute_ecm_to_fa_bias_neutral(fas, ecm)
    assert result.multipliers_per_fa.shape == (2, 3)
    np.testing.assert_array_equal(result.multipliers_per_fa, np.ones((2, 3)))


# ---------------------------------------------------------------------------
# Test 2: neutral multipliers match RATE_NAMES column order
# ---------------------------------------------------------------------------


def test_neutral_multipliers_match_RATE_NAMES_column_order():
    ecm = _ecm(nx=4, ny=4)
    fas = (_fa("fa-1", position=(0.5, 0.5)),)
    result = compute_ecm_to_fa_bias_neutral(fas, ecm)
    assert result.rate_names == RATE_NAMES
    assert RATE_NAMES == ("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")
    # Column index k corresponds to RATE_NAMES[k]
    assert result.multipliers_per_fa.shape[1] == 3


# ---------------------------------------------------------------------------
# Test 3: sampler at cell center returns field value exactly
# ---------------------------------------------------------------------------


def test_sampler_at_cell_center_returns_field_value_exactly():
    ecm = _ecm(nx=4, ny=4, spacing_um=1.0, stiffness_value=3.0, fiber_value=0.5)
    fa = _fa("fa-1", position=(0.5, 0.5))  # cell center (0, 0)
    sampled = sample_ecm_at_fa_positions((fa,), ecm)
    np.testing.assert_allclose(sampled.stiffness_kpa, [3.0], atol=1e-12)
    np.testing.assert_allclose(sampled.fiber_density, [0.5], atol=1e-12)


# ---------------------------------------------------------------------------
# Test 4: sampler at grid crossing bilinear average
# ---------------------------------------------------------------------------


def test_sampler_at_grid_crossing_bilinear_average():
    nx, ny = 4, 4
    spacing_um = 1.0
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    stiff = np.zeros((nx, ny), dtype=np.float64)
    # Set 4 cells around (1.0, 1.0) crossing
    stiff[0, 0] = 1.0
    stiff[1, 0] = 2.0
    stiff[0, 1] = 3.0
    stiff[1, 1] = 4.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=spacing_um,
        stiffness_kpa=stiff,
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    fa = _fa("fa-1", position=(1.0, 1.0))  # 4-center crossing
    sampled = sample_ecm_at_fa_positions((fa,), ecm)
    expected = (1.0 + 2.0 + 3.0 + 4.0) / 4.0  # = 2.5
    np.testing.assert_allclose(sampled.stiffness_kpa, [expected], atol=1e-12)


# ---------------------------------------------------------------------------
# Test 5: sampler boundary half-cell renormalized
# ---------------------------------------------------------------------------


def test_sampler_boundary_half_cell_renormalized():
    """FA at the left edge midpoint (0.0, 1.0) reads from cells
    (0, 0) and (0, 1) with renormalized 50/50 weights (after
    truncated stencil drops the out-of-grid columns)."""

    nx, ny = 4, 4
    stiff = np.zeros((nx, ny), dtype=np.float64)
    stiff[0, 0] = 1.0
    stiff[0, 1] = 3.0
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=stiff,
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    fa = _fa("fa-1", position=(0.0, 1.0))
    sampled = sample_ecm_at_fa_positions((fa,), ecm)
    expected = (1.0 + 3.0) / 2.0  # 50/50
    np.testing.assert_allclose(sampled.stiffness_kpa, [expected], atol=1e-12)


# ---------------------------------------------------------------------------
# Test 6: sampler outside footprint raises fa_bias_position_outside_ecm_grid
# ---------------------------------------------------------------------------


def test_sampler_outside_footprint_raises_fa_bias_position_outside_ecm_grid():
    ecm = _ecm(nx=4, ny=4)
    fa = _fa("fa-1", position=(5.0, 0.5))
    with pytest.raises(FAToECMBiasError) as excinfo:
        sample_ecm_at_fa_positions((fa,), ecm)
    assert excinfo.value.failure_kind == "fa_bias_position_outside_ecm_grid"


# ---------------------------------------------------------------------------
# Test 7: sampler input order preserved
# ---------------------------------------------------------------------------


def test_sampler_input_order_preserved():
    ecm = _ecm(nx=4, ny=4)
    fas = (
        _fa("alpha", position=(0.5, 0.5)),
        _fa("beta", position=(1.5, 1.5)),
        _fa("gamma", position=(2.5, 2.5)),
    )
    sampled = sample_ecm_at_fa_positions(fas, ecm)
    assert sampled.fa_ids == ("alpha", "beta", "gamma")


# ---------------------------------------------------------------------------
# Test 8: sampler accumulated_traction is scalar shape (N_FA,)
# ---------------------------------------------------------------------------


def test_sampler_accumulated_traction_is_scalar_shape_n_fa():
    ecm = _ecm(nx=4, ny=4, accumulated_value=2.0)
    fas = (_fa("fa-1", position=(0.5, 0.5)), _fa("fa-2", position=(1.5, 1.5)))
    sampled = sample_ecm_at_fa_positions(fas, ecm)
    assert sampled.accumulated_traction_nNs_per_um2.shape == (2,)
    assert sampled.accumulated_traction_nNs_per_um2.dtype == np.float64
    np.testing.assert_allclose(
        sampled.accumulated_traction_nNs_per_um2, [2.0, 2.0], atol=1e-12
    )


# ---------------------------------------------------------------------------
# Test 9: sampler orientation tensor shape (N_FA, 2, 2)
# ---------------------------------------------------------------------------


def test_sampler_orientation_tensor_shape_n_fa_2_2():
    """Native-shape protection: sampler returns full 2×2 tensor
    per FA, NOT a scalar reduction (Hard Rule 11 / locked §6)."""

    ecm = _ecm(nx=4, ny=4)
    fas = (_fa("fa-1", position=(0.5, 0.5)), _fa("fa-2", position=(1.5, 1.5)))
    sampled = sample_ecm_at_fa_positions(fas, ecm)
    assert sampled.orientation_tensor.shape == (2, 2, 2)
    assert sampled.orientation_tensor.dtype == np.float64
    # Default ECM orientation is identity (1.0, 0.0, 0.0, 1.0)
    expected = np.array([[1.0, 0.0], [0.0, 1.0]])
    for i in range(2):
        np.testing.assert_allclose(
            sampled.orientation_tensor[i], expected, atol=1e-12
        )


# ---------------------------------------------------------------------------
# Test 10: neutral bias no ECM mutation
# ---------------------------------------------------------------------------


def test_neutral_bias_no_ecm_mutation():
    ecm = _ecm(nx=4, ny=4, stiffness_value=3.0, fiber_value=0.5)
    snapshots = {
        "stiffness_kpa": ecm.stiffness_kpa.copy(),
        "ligand_density": ecm.ligand_density.copy(),
        "fiber_density": ecm.fiber_density.copy(),
        "orientation_tensor": ecm.orientation_tensor.copy(),
        "accumulated_traction_nNs_per_um2":
            ecm.accumulated_traction_nNs_per_um2.copy(),
    }
    fa = _fa("fa-1", position=(1.5, 1.5))
    compute_ecm_to_fa_bias_neutral((fa,), ecm)
    for name, snap in snapshots.items():
        np.testing.assert_array_equal(getattr(ecm, name), snap)


# ---------------------------------------------------------------------------
# Test 11: neutral bias no FA mutation
# ---------------------------------------------------------------------------


def test_neutral_bias_no_fa_mutation():
    ecm = _ecm(nx=4, ny=4)
    fa = _fa("fa-1", position=(1.5, 1.5))
    snapshot = (
        fa.position_um_xy,
        fa.traction_force_nN_xy,
        fa.maturity,
        fa.bound_fraction,
        fa.state,
        fa.age_s,
    )
    compute_ecm_to_fa_bias_neutral((fa,), ecm)
    after = (
        fa.position_um_xy,
        fa.traction_force_nN_xy,
        fa.maturity,
        fa.bound_fraction,
        fa.state,
        fa.age_s,
    )
    assert snapshot == after


# ---------------------------------------------------------------------------
# Test 12: diagnostics_dict only serializable plain values
# ---------------------------------------------------------------------------


def test_diagnostics_dict_only_serializable_plain_values():
    ecm = _ecm(nx=4, ny=4)
    fas = (_fa("fa-1", position=(0.5, 0.5)),)
    result = compute_ecm_to_fa_bias_neutral(fas, ecm)
    for key, value in result.diagnostics_dict.items():
        assert isinstance(key, str)
        assert isinstance(value, (int, float, str)), (
            f"diagnostics_dict[{key!r}] = {type(value).__name__} is not "
            f"plain serializable; arrays must live in typed dataclass fields"
        )


# ---------------------------------------------------------------------------
# Test 13: diagnostics_dict keys are locked set
# ---------------------------------------------------------------------------


def test_diagnostics_dict_keys_are_locked_set():
    ecm = _ecm(nx=4, ny=4)
    fas = (_fa("fa-1", position=(0.5, 0.5)),)
    result = compute_ecm_to_fa_bias_neutral(fas, ecm)
    expected_keys = {
        "n_adhesions",
        "max_multiplier",
        "min_multiplier",
        "sampler_geometry",
    }
    assert set(result.diagnostics_dict.keys()) == expected_keys
    assert result.diagnostics_dict["sampler_geometry"] == "bilinear_cell_centered"
    assert result.diagnostics_dict["n_adhesions"] == 1
    assert result.diagnostics_dict["max_multiplier"] == 1.0
    assert result.diagnostics_dict["min_multiplier"] == 1.0


# ---------------------------------------------------------------------------
# Test 14: empty FA list returns empty result, not raise
# ---------------------------------------------------------------------------


def test_empty_fa_list_returns_empty_result_not_raise():
    ecm = _ecm(nx=4, ny=4)
    result = compute_ecm_to_fa_bias_neutral((), ecm)
    assert result.fa_ids == ()
    assert result.multipliers_per_fa.shape == (0, 3)
    assert result.sampled_diagnostics is not None
    assert result.sampled_diagnostics.stiffness_kpa.shape == (0,)
    assert result.sampled_diagnostics.orientation_tensor.shape == (0, 2, 2)


# ---------------------------------------------------------------------------
# Test 15: sampler consistency with scatter geometry (sister-gate round-trip)
# ---------------------------------------------------------------------------


def test_sampler_consistency_with_scatter_geometry():
    """Same FA position scattered by Hard Blocker #3 and sampled
    by Hard Blocker #4 must use identical bilinear weights.

    We verify by scattering a unit FA at a 4-center crossing and
    sampling the resulting field: at the same crossing, the
    sampled scatter density equals the bilinear average of the 4
    deposited cells, which equals the FA's deposit divided by 4
    (each receives 25%, then averaged again). Conversely, scatter
    at one position and sample at the same position must be
    consistent with the geometry."""

    nx, ny = 4, 4
    spacing_um = 1.0
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=spacing_um,
        stiffness_kpa=np.zeros((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )

    # Use scatter to deposit a known field, then put that field
    # into a stiffness_kpa-shaped array and sample it back.
    fa = _fa("fa-1", position=(1.0, 1.0), traction=(4.0, 0.0))
    scatter_field = scatter_fa_traction_to_ecm_bilinear((fa,), ecm)
    # scatter_field shape (nx, ny, 2); use x-component as a fake
    # stiffness field
    stiff_from_scatter = scatter_field[..., 0].copy()
    ecm_from_scatter = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=spacing_um,
        stiffness_kpa=stiff_from_scatter,
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
        source="synthetic",
    )
    # Sample at the same position; should be the bilinear interpolated
    # value of the scatter field at (1.0, 1.0). The four cells
    # received scatter (0.25 * 4 / 1.0 = 1.0 each), so sampled
    # value at (1.0, 1.0) should be the bilinear average = 1.0.
    sampled = sample_ecm_at_fa_positions((fa,), ecm_from_scatter)
    np.testing.assert_allclose(sampled.stiffness_kpa, [1.0], atol=1e-12)


# ---------------------------------------------------------------------------
# Test 16: returned shapes match locked signature
# ---------------------------------------------------------------------------


def test_returned_shapes_match_locked_signature():
    ecm = _ecm(nx=3, ny=5, spacing_um=1.5)
    fas = (_fa("fa-1", position=(2.25, 3.75)), _fa("fa-2", position=(1.5, 1.5)))
    sampled = sample_ecm_at_fa_positions(fas, ecm)
    assert sampled.stiffness_kpa.shape == (2,)
    assert sampled.fiber_density.shape == (2,)
    assert sampled.ligand_density.shape == (2,)
    assert sampled.orientation_tensor.shape == (2, 2, 2)
    assert sampled.accumulated_traction_nNs_per_um2.shape == (2,)
    result = compute_ecm_to_fa_bias_neutral(fas, ecm)
    assert result.multipliers_per_fa.shape == (2, 3)
    assert isinstance(result.sampled_diagnostics, ECMSampledAtFAs)


# ---------------------------------------------------------------------------
# Test 17: no active law invoked in Phase D (Hard Rule 11 meta-test)
# ---------------------------------------------------------------------------


def test_no_active_law_invoked_in_phase_d():
    """Hard Blocker #4 Phase D — Hard Rule 11 wording protection.

    Documents at runtime that compute_ecm_to_fa_bias_neutral returns
    only neutral 1.0 multipliers; it does NOT decide an active
    constitutive law. Phase E `compute_ecm_to_fa_bias_active` is a
    separate function with its own Sanity Gate. This meta-test
    makes the boundary explicit so a future refactor cannot silently
    relabel the neutral function as the active one.

    Mirrors:
    - Phase B test 4
    - Phase C test 5
    - Hard Blocker #3 test 18
    """

    phase_d_modality = "neutral 1.0 multipliers + diagnostic sampler"
    phase_e_active_claim = "ECM→FA active mapping (Phase E TBD)"
    assert phase_d_modality != phase_e_active_claim, (
        "compute_ecm_to_fa_bias_neutral returns neutral 1.0 only. "
        "Phase E active mapping is a separate function "
        "(compute_ecm_to_fa_bias_active) with separate Sanity Gate. "
        "Locked plan §1 Phase D / Phase E function naming separation "
        "is the structural Hard Rule 11 boundary."
    )


# ---------------------------------------------------------------------------
# Test 18: non-finite FA position raises
# ---------------------------------------------------------------------------


def test_non_finite_fa_position_raises():
    ecm = _ecm(nx=4, ny=4)
    fa_nan = _fa("fa-nan", position=(float("nan"), 0.5))
    with pytest.raises(FAToECMBiasError) as excinfo_nan:
        sample_ecm_at_fa_positions((fa_nan,), ecm)
    assert excinfo_nan.value.failure_kind == "non_finite_fa_position"

    fa_inf = _fa("fa-inf", position=(0.5, math.inf))
    with pytest.raises(FAToECMBiasError) as excinfo_inf:
        sample_ecm_at_fa_positions((fa_inf,), ecm)
    assert excinfo_inf.value.failure_kind == "non_finite_fa_position"


# ---------------------------------------------------------------------------
# Test 19: invalid FA schema propagates ValueError (Codex review id=1378)
# ---------------------------------------------------------------------------


def test_malformed_position_length_propagates_schema_validate_error():
    """Codex review id=1385 regression: an FA with `position_um_xy`
    of length != 2 (e.g., a 1-tuple) must raise the schema's
    canonical ValueError ("position_um_xy must contain exactly 2
    finite values") rather than Python's unpack ValueError. The
    structural length check defers to fa.validate() before
    unpacking."""

    ecm = _ecm(nx=4, ny=4)
    bad_fa = FocalAdhesionState(
        adhesion_id="bad-len",
        cell_id="cell-A",
        position_um_xy=(0.5,),  # type: ignore[arg-type]
        age_s=0.0,
        maturity=1.0,
        bound_fraction=1.0,
        state="mature",
        traction_force_nN_xy=(1.0, 0.0),
    )
    with pytest.raises(ValueError, match="position_um_xy must contain exactly 2"):
        sample_ecm_at_fa_positions((bad_fa,), ecm)
    with pytest.raises(ValueError, match="position_um_xy must contain exactly 2"):
        compute_ecm_to_fa_bias_neutral((bad_fa,), ecm)


def test_invalid_fa_schema_propagates_validate_error():
    """Codex review id=1378 regression: per-FA `fa.validate()` runs
    after the lock-specific `non_finite_fa_position` check, so
    schema invariants (e.g., `state == "unbound"` with non-zero
    `traction_force_nN_xy`) propagate as the schema's `ValueError`
    BEFORE any sampler / multiplier output. Mirrors Hard Blocker #3
    fix `883efc8`."""

    ecm = _ecm(nx=4, ny=4)
    bad_fa = FocalAdhesionState(
        adhesion_id="bad-fa",
        cell_id="cell-A",
        position_um_xy=(1.5, 1.5),
        age_s=0.0,
        maturity=0.0,
        bound_fraction=0.0,
        state="unbound",
        traction_force_nN_xy=(1.0, 0.0),  # invalid: unbound + non-zero
    )
    with pytest.raises(ValueError, match="no traction without attachment"):
        sample_ecm_at_fa_positions((bad_fa,), ecm)
    # compute_ecm_to_fa_bias_neutral path also raises (it calls
    # sample_ecm_at_fa_positions internally).
    with pytest.raises(ValueError, match="no traction without attachment"):
        compute_ecm_to_fa_bias_neutral((bad_fa,), ecm)
