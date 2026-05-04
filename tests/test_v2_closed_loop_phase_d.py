"""Tests for Phase D no-op closed-loop scaffolding integrator.

Test catalog mirrors
``docs/v2_phase_d_no_op_scaffolding_locked.md`` §6 (commit
``7786b20``) + ``docs/v2_phase_d_no_op_scaffolding_sanity_gate.md``
§8 (commit ``42b8e34``). Each test maps to one Sanity Gate item or
one forbidden behavior.

Per Codex review id=1420 monkeypatch note, tests 5 + 9 patch the
wrapper-layer ``step_ecm_to_fa_bias`` symbol in
``acs.v2.dynamics.closed_loop_phase_d`` rather than the HB#4
primitive import — keeping the test robust against wrapper
internal-import-style changes.
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from acs.v2.dynamics import closed_loop_phase_d as cl_phase_d
from acs.v2.dynamics.closed_loop_phase_d import (
    FAToECMResponseResult,
    PhaseDNoOpStepResult,
    step_ecm_to_fa_bias,
    step_fa_to_ecm_response,
    step_phase_d_no_op,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMSampledAtFAs,
    ECMToFABiasResult,
    FAToECMBiasError,
    RATE_NAMES,
)
from acs.v2.dynamics.fa_to_ecm_scattering import (
    FAToECMScatteringError,
    scatter_fa_traction_to_ecm_bilinear,
)
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
# Test 1: returns typed dataclass with locked fields (frozen + slots)
# ---------------------------------------------------------------------------
def test_phase_d_step_returns_typed_dataclass():
    ecm = _ecm()
    fa = _fa("a1", position=(1.5, 1.5), traction=(2.0, -3.0))
    result = step_phase_d_no_op((fa,), ecm)
    assert isinstance(result, PhaseDNoOpStepResult)
    # Frozen dataclass: assignment must raise.
    with pytest.raises((AttributeError, TypeError)):
        result.fa_to_ecm = result.fa_to_ecm  # type: ignore[misc]
    # Locked fields exactly.
    assert set(PhaseDNoOpStepResult.__dataclass_fields__.keys()) == {
        "fa_to_ecm",
        "ecm_to_fa",
        "updated_ecm",
    }
    assert isinstance(result.fa_to_ecm, FAToECMResponseResult)
    assert isinstance(result.ecm_to_fa, ECMToFABiasResult)
    assert isinstance(result.updated_ecm, ECMSubstrateState)


# ---------------------------------------------------------------------------
# Test 2: triple-identity invariant
# ---------------------------------------------------------------------------
def test_phase_d_step_object_identity_updated_ecm():
    ecm = _ecm()
    fa = _fa("a1", position=(1.5, 1.5), traction=(2.0, -3.0))
    result = step_phase_d_no_op((fa,), ecm)
    assert result.updated_ecm is ecm
    assert result.fa_to_ecm.updated_ecm is ecm
    # Same for the leg called directly.
    leg_result = step_fa_to_ecm_response((fa,), ecm)
    assert leg_result.updated_ecm is ecm


# ---------------------------------------------------------------------------
# Test 3: no input mutation
# ---------------------------------------------------------------------------
def test_phase_d_step_no_input_mutation():
    ecm = _ecm(stiffness_value=2.5, fiber_value=0.7, ligand_value=0.9)
    fa = _fa("a1", position=(1.5, 1.5), traction=(2.0, -3.0))
    stiffness_before = ecm.stiffness_kpa.copy()
    fiber_before = ecm.fiber_density.copy()
    ligand_before = ecm.ligand_density.copy()
    orientation_before = ecm.orientation_tensor.copy()
    accumulated_before = ecm.accumulated_traction_nNs_per_um2.copy()
    fa_position_before = fa.position_um_xy
    fa_traction_before = fa.traction_force_nN_xy
    fa_state_before = fa.state

    _ = step_phase_d_no_op((fa,), ecm)

    assert np.array_equal(ecm.stiffness_kpa, stiffness_before)
    assert np.array_equal(ecm.fiber_density, fiber_before)
    assert np.array_equal(ecm.ligand_density, ligand_before)
    assert np.array_equal(ecm.orientation_tensor, orientation_before)
    assert np.array_equal(
        ecm.accumulated_traction_nNs_per_um2, accumulated_before
    )
    assert fa.position_um_xy == fa_position_before
    assert fa.traction_force_nN_xy == fa_traction_before
    assert fa.state == fa_state_before


# ---------------------------------------------------------------------------
# Test 4: empty FA list invariants (per locked §4)
# ---------------------------------------------------------------------------
def test_phase_d_step_empty_fa_list():
    ecm = _ecm(nx=5, ny=7)
    result = step_phase_d_no_op((), ecm)

    # FA→ECM leg invariants.
    assert result.fa_to_ecm.traction_density_xy.shape == (5, 7, 2)
    assert result.fa_to_ecm.traction_density_xy.dtype == np.float64
    assert np.all(result.fa_to_ecm.traction_density_xy == 0.0)
    assert result.fa_to_ecm.updated_ecm is ecm

    # ECM→FA leg invariants.
    assert result.ecm_to_fa.multipliers_per_fa.shape == (0, 3)
    assert result.ecm_to_fa.multipliers_per_fa.dtype == np.float64
    assert result.ecm_to_fa.rate_names == RATE_NAMES
    # Non-None empty ECMSampledAtFAs (per HB#4 sampler semantics —
    # never None in Phase D).
    assert result.ecm_to_fa.sampled_diagnostics is not None
    assert isinstance(result.ecm_to_fa.sampled_diagnostics, ECMSampledAtFAs)
    assert result.ecm_to_fa.sampled_diagnostics.fa_ids == ()
    assert result.ecm_to_fa.sampled_diagnostics.stiffness_kpa.shape == (0,)

    # Triple-identity.
    assert result.updated_ecm is ecm

    # Also accept list[] form.
    result_list = step_phase_d_no_op([], ecm)
    assert result_list.fa_to_ecm.traction_density_xy.shape == (5, 7, 2)
    assert np.all(result_list.fa_to_ecm.traction_density_xy == 0.0)


# ---------------------------------------------------------------------------
# Test 5: deterministic call order — out-of-grid surfaces as HB#3 first
# ---------------------------------------------------------------------------
def test_phase_d_step_call_order_deterministic(monkeypatch):
    ecm = _ecm()
    # FA outside the inclusive footprint [0, 4] x [0, 4].
    fa = _fa("a1", position=(10.0, 10.0), traction=(1.0, 0.0))

    bias_called = {"count": 0}

    def _spy_step_ecm_to_fa_bias(adhesions, ecm_arg):
        bias_called["count"] += 1
        return step_ecm_to_fa_bias(adhesions, ecm_arg)

    monkeypatch.setattr(
        cl_phase_d, "step_ecm_to_fa_bias", _spy_step_ecm_to_fa_bias
    )

    with pytest.raises(FAToECMScatteringError) as exc_info:
        step_phase_d_no_op((fa,), ecm)
    assert exc_info.value.failure_kind == "fa_position_outside_ecm_grid"
    # HB#4 must NOT have been called when HB#3 raises.
    assert bias_called["count"] == 0


# ---------------------------------------------------------------------------
# Test 6: structural meta-test — no Phase E kwargs
# ---------------------------------------------------------------------------
def test_phase_d_step_no_phase_e_kwargs():
    for fn in (
        step_phase_d_no_op,
        step_fa_to_ecm_response,
        step_ecm_to_fa_bias,
    ):
        sig = inspect.signature(fn)
        params = list(sig.parameters.values())
        assert (
            len(params) == 2
        ), f"{fn.__name__} must have exactly 2 parameters, got {len(params)}"
        names = [p.name for p in params]
        assert names == ["adhesions", "ecm"], (
            f"{fn.__name__} parameter names {names} must be exactly "
            f"['adhesions', 'ecm']"
        )
        for p in params:
            assert (
                p.default is inspect.Parameter.empty
            ), (
                f"{fn.__name__} parameter {p.name!r} must not have a "
                f"default value (silent-activation guard)"
            )
            assert p.kind in (
                inspect.Parameter.POSITIONAL_OR_KEYWORD,
                inspect.Parameter.POSITIONAL_ONLY,
            ), (
                f"{fn.__name__} parameter {p.name!r} kind {p.kind} not "
                f"allowed; must be positional-or-keyword or positional-only"
            )


# ---------------------------------------------------------------------------
# Test 7: legs match primitives exactly
# ---------------------------------------------------------------------------
def test_phase_d_step_legs_match_primitives():
    ecm = _ecm(stiffness_value=3.0, fiber_value=0.5)
    adhesions = (
        _fa("a1", position=(1.5, 1.5), traction=(2.0, -3.0)),
        _fa("a2", position=(2.5, 2.5), traction=(-1.0, 4.0)),
    )

    result = step_phase_d_no_op(adhesions, ecm)

    expected_traction = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    np.testing.assert_array_equal(
        result.fa_to_ecm.traction_density_xy, expected_traction
    )

    leg_fa_to_ecm = step_fa_to_ecm_response(adhesions, ecm)
    np.testing.assert_array_equal(
        leg_fa_to_ecm.traction_density_xy, expected_traction
    )

    leg_ecm_to_fa = step_ecm_to_fa_bias(adhesions, ecm)
    np.testing.assert_array_equal(
        result.ecm_to_fa.multipliers_per_fa, leg_ecm_to_fa.multipliers_per_fa
    )
    assert result.ecm_to_fa.rate_names == leg_ecm_to_fa.rate_names
    assert result.ecm_to_fa.fa_ids == leg_ecm_to_fa.fa_ids


# ---------------------------------------------------------------------------
# Test 8: HB#3 failure_kinds propagate (non_finite_fa_position)
# ---------------------------------------------------------------------------
def test_phase_d_step_propagates_hb3_failure_kinds():
    ecm = _ecm()
    fa = _fa("a1", position=(float("nan"), 1.5), traction=(1.0, 0.0))
    with pytest.raises(FAToECMScatteringError) as exc_info:
        step_phase_d_no_op((fa,), ecm)
    assert exc_info.value.failure_kind == "non_finite_fa_position"


# ---------------------------------------------------------------------------
# Test 9: HB#4 failure propagation when bias leg raises (mocked/stubbed)
# ---------------------------------------------------------------------------
def test_phase_d_step_propagates_hb4_failure_kinds_when_bias_leg_raises(
    monkeypatch,
):
    ecm = _ecm()
    fa = _fa("a1", position=(1.5, 1.5), traction=(1.0, 0.0))

    def _stubbed_bias(adhesions, ecm_arg):
        raise FAToECMBiasError(
            "fa_bias_position_outside_ecm_grid",
            "stubbed HB#4 failure for test 9 — natural geometry should "
            "match HB#3 sister-gate",
        )

    monkeypatch.setattr(cl_phase_d, "step_ecm_to_fa_bias", _stubbed_bias)

    with pytest.raises(FAToECMBiasError) as exc_info:
        step_phase_d_no_op((fa,), ecm)
    assert exc_info.value.failure_kind == "fa_bias_position_outside_ecm_grid"


# ---------------------------------------------------------------------------
# Test 10: HB#3 scatter conservation propagated through the wrapper
# ---------------------------------------------------------------------------
def test_phase_d_step_scatter_conservation_inherited():
    ecm = _ecm(nx=6, ny=6)
    cell_area = ecm.spacing_um * ecm.spacing_um
    adhesions = (
        _fa("a1", position=(2.5, 2.5), traction=(2.0, -3.0)),
        _fa("a2", position=(4.0, 1.5), traction=(-1.0, 4.0)),
        _fa("a3", position=(0.7, 5.3), traction=(0.5, 0.5)),
    )
    result = step_phase_d_no_op(adhesions, ecm)
    total_x = float(result.fa_to_ecm.traction_density_xy[..., 0].sum() * cell_area)
    total_y = float(result.fa_to_ecm.traction_density_xy[..., 1].sum() * cell_area)
    expected_x = sum(fa.traction_force_nN_xy[0] for fa in adhesions)
    expected_y = sum(fa.traction_force_nN_xy[1] for fa in adhesions)
    assert math.isclose(total_x, expected_x, abs_tol=1e-9)
    assert math.isclose(total_y, expected_y, abs_tol=1e-9)


# ---------------------------------------------------------------------------
# Test 11: Hard Rule 11 runtime meta-test — wrapper does NOT update ECM
# ---------------------------------------------------------------------------
def test_phase_d_step_does_not_satisfy_phase_e_response_law():
    """Phase D wrapper must NOT update ECM fields based on traction.

    A Phase E response law would mutate `stiffness_kpa`,
    `fiber_density`, `orientation_tensor`, or
    `accumulated_traction_nNs_per_um2` in response to the scatter
    output. Phase D must leave all fields byte-identical.
    """

    ecm = _ecm(stiffness_value=2.0, fiber_value=0.4)
    # Non-trivial scatter input — multiple FAs with different
    # tractions to ensure the wrapper has a real scatter array to
    # hypothetically respond to.
    adhesions = (
        _fa("a1", position=(1.5, 1.5), traction=(5.0, 0.0)),
        _fa("a2", position=(2.5, 2.5), traction=(0.0, -7.0)),
        _fa("a3", position=(0.7, 3.3), traction=(3.0, 4.0)),
    )
    stiffness_before = ecm.stiffness_kpa.copy()
    fiber_before = ecm.fiber_density.copy()
    ligand_before = ecm.ligand_density.copy()
    orientation_before = ecm.orientation_tensor.copy()
    accumulated_before = ecm.accumulated_traction_nNs_per_um2.copy()

    result = step_phase_d_no_op(adhesions, ecm)

    # Object identity invariant.
    assert result.updated_ecm is ecm
    assert result.fa_to_ecm.updated_ecm is ecm

    # Field-level byte-identity (no Phase E response law applied).
    assert np.array_equal(ecm.stiffness_kpa, stiffness_before)
    assert np.array_equal(ecm.fiber_density, fiber_before)
    assert np.array_equal(ecm.ligand_density, ligand_before)
    assert np.array_equal(ecm.orientation_tensor, orientation_before)
    assert np.array_equal(
        ecm.accumulated_traction_nNs_per_um2, accumulated_before
    )

    # Multipliers are the all-1.0 neutral default (no active bias
    # law applied either).
    assert np.all(result.ecm_to_fa.multipliers_per_fa == 1.0)


# ---------------------------------------------------------------------------
# Test 12: structural meta-test — no diagnostics_dict field
# ---------------------------------------------------------------------------
def test_phase_d_step_no_diagnostics_dict_field():
    assert set(FAToECMResponseResult.__dataclass_fields__.keys()) == {
        "traction_density_xy",
        "updated_ecm",
    }
