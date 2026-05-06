"""Tests for Phase E v1 closed-loop ECM gate composition.

Test catalog mirrors
``docs/v2/v2_phase_e_composition_locked.md`` §4 (commit ``c8b9550``)
and the impl-work Sanity Gate
``docs/v2/v2_phase_e_composition_sanity_gate.md`` §8 (commit
``b0b5baa`` + amendments ``a3bc8a6`` + ``e5a8186``). 17 tests
per locked §4.

Per Codex `id=1552` Sanity Gate PASS, this catalog covers:
- Composition invariants (5): tests 1, 2, 3, 4, 5
- Error propagation (4): tests 6, 7, 8, 9
- ECM-side evidence Items 1-4 (4): tests 10, 11, 12, 13
- Meta + guard (3): tests 14, 15, 16
- Exports (1): test 17
"""

from __future__ import annotations

import ast
import inspect
import math

import numpy as np
import pytest

from acs.v2 import dynamics as v2_dynamics
import acs.v2 as v2_pkg
from acs.v2.dynamics import closed_loop_phase_e as ce
from acs.v2.dynamics import ecm_constitutive_response as cr
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v1,
    step_closed_loop_phase_e_v2,
)
from acs.v2.dynamics.ecm_constitutive_response import (
    ECMConstitutiveResponseError,
    ECMOrientationResponseResult,
)
from acs.v2.dynamics.ecm_lyapunov_metric import (
    ECMOrientationLyapunovMetricResult,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    FAToECMBiasError,
)
from acs.v2.dynamics.fa_to_ecm_scattering import (
    FAToECMScatteringError,
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
    orientation_value: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 1.0),
):
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


def _fa(
    adhesion_id: str,
    *,
    position=(1.5, 1.5),
    cell_id: str = "cell-A",
    state="mature",
    bound_fraction: float = 1.0,
    traction=(2.0, 0.0),
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


def test_phase_e_v1_returns_phase_e_step_result_dataclass():
    """Test 1: return is PhaseEStepResult with all 5 expected fields
    populated."""

    ecm = _ecm()
    fa = _fa("a1")
    result = step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert isinstance(result, PhaseEStepResult)
    assert isinstance(result.traction_density_xy, np.ndarray)
    assert result.traction_density_xy.shape == (4, 4, 2)
    assert isinstance(result.orientation_response, ECMOrientationResponseResult)
    assert isinstance(result.ecm_to_fa_bias, ECMToFABiasResult)
    assert isinstance(result.lyapunov_metric, ECMOrientationLyapunovMetricResult)
    assert isinstance(result.updated_ecm, ECMSubstrateState)


def test_phase_e_v1_identity_invariant_updated_ecm_is_orientation_response_updated_ecm():
    """Test 2 (Y4): result.updated_ecm is
    result.orientation_response.updated_ecm — object identity, not
    just bytewise equality."""

    ecm = _ecm()
    fa = _fa("a1")
    result = step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert result.updated_ecm is result.orientation_response.updated_ecm


def test_phase_e_v1_module_does_not_define_wrapper_failure_kinds():
    """Test 3 (Y17): source/meta guard — no PhaseE*Error class
    definition, no failure_kind = assignment in module body."""

    src = inspect.getsource(ce)
    tree = ast.parse(src)

    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            assert not (
                node.name.startswith("PhaseE") and node.name.endswith("Error")
            ), (
                f"Phase E composition must not define wrapper failure kinds, "
                f"found class {node.name}"
            )

    src_no_spaces = src.replace(" ", "")
    assert "failure_kind=" not in src_no_spaces, (
        "Phase E composition must not introduce failure_kind labels; "
        "sub-call errors propagate as-is"
    )


def test_phase_e_v1_call_order_via_monkeypatch(monkeypatch):
    """Test 4 (Y5 + Y10 corrected pattern + Y12 monkeypatch where used):
    HB#3 → HB#1+#2 → HB#4 → HB#5 order; HB#4 receives (adhesions,
    post-update ECM); HB#5 receives (post-update ECM, original
    scatter)."""

    calls: list[tuple[str, tuple, dict]] = []

    real_scatter = ce.scatter_fa_traction_to_ecm_bilinear
    real_orient = ce.step_ecm_orientation_response
    real_bias = ce.compute_ecm_to_fa_bias_neutral
    real_lyap = ce.compute_ecm_orientation_lyapunov_metric

    def _spy_scatter(*args, **kwargs):
        out = real_scatter(*args, **kwargs)
        calls.append(("scatter", args, {"out": out}))
        return out

    def _spy_orient(*args, **kwargs):
        out = real_orient(*args, **kwargs)
        calls.append(("orient", args, {"out": out}))
        return out

    def _spy_bias(*args, **kwargs):
        out = real_bias(*args, **kwargs)
        calls.append(("bias", args, {"out": out}))
        return out

    def _spy_lyap(*args, **kwargs):
        out = real_lyap(*args, **kwargs)
        calls.append(("lyap", args, {"out": out}))
        return out

    monkeypatch.setattr(ce, "scatter_fa_traction_to_ecm_bilinear", _spy_scatter)
    monkeypatch.setattr(ce, "step_ecm_orientation_response", _spy_orient)
    monkeypatch.setattr(ce, "compute_ecm_to_fa_bias_neutral", _spy_bias)
    monkeypatch.setattr(ce, "compute_ecm_orientation_lyapunov_metric", _spy_lyap)

    ecm = _ecm()
    fa = _fa("a1")
    result = step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)

    assert [c[0] for c in calls] == ["scatter", "orient", "bias", "lyap"]

    scatter_out = calls[0][2]["out"]
    orient_out = calls[1][2]["out"]
    post_update_ecm = orient_out.updated_ecm

    bias_args = calls[2][1]
    assert bias_args[1] is post_update_ecm, (
        "HB#4 must receive post-HB#1+#2 ECM as second arg, not pre-update"
    )

    lyap_args = calls[3][1]
    assert lyap_args[0] is post_update_ecm, (
        "HB#5 must receive post-HB#1+#2 ECM as first arg, not pre-update"
    )
    assert lyap_args[1] is scatter_out, (
        "HB#5 must receive original HB#3 scatter output, not recomputed"
    )

    assert result.updated_ecm is post_update_ecm


def test_phase_e_v1_empty_adhesions_chain():
    """Test 5 (Y7 + Y11): empty FA list → all-zero scatter (cite
    fa_to_ecm_scattering.py:170); orientation bytewise unchanged but
    fresh array (no aliasing); HB#4 result shape (0, 3); V_active = 0;
    identity invariant holds; result.updated_ecm is not ecm."""

    ecm = _ecm()
    result = step_closed_loop_phase_e_v1((), ecm, dt_s=600.0)

    assert np.all(result.traction_density_xy == 0.0)
    assert result.traction_density_xy.shape == (4, 4, 2)

    np.testing.assert_array_equal(
        result.updated_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert not np.shares_memory(
        result.updated_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert result.updated_ecm is not ecm

    assert result.ecm_to_fa_bias.multipliers_per_fa.shape == (0, 3)

    assert result.lyapunov_metric.v_active_um2 == 0.0

    assert result.updated_ecm is result.orientation_response.updated_ecm


def test_phase_e_v1_propagates_hb3_error_stops_before_orient(monkeypatch):
    """Test 6 (Y6): monkeypatch HB#3 to raise; verify HB#1+#2 / HB#4 /
    HB#5 not called."""

    calls = []
    monkeypatch.setattr(
        ce,
        "scatter_fa_traction_to_ecm_bilinear",
        lambda *a, **k: (_ for _ in ()).throw(
            FAToECMScatteringError(
                "fa_position_outside_ecm_grid",
                "test stub HB#3 failure",
            )
        ),
    )
    monkeypatch.setattr(
        ce,
        "step_ecm_orientation_response",
        lambda *a, **k: calls.append("orient") or (_ for _ in ()).throw(
            AssertionError("HB#1+#2 must NOT be called when HB#3 raises")
        ),
    )
    monkeypatch.setattr(
        ce,
        "compute_ecm_to_fa_bias_neutral",
        lambda *a, **k: calls.append("bias") or (_ for _ in ()).throw(
            AssertionError("HB#4 must NOT be called when HB#3 raises")
        ),
    )
    monkeypatch.setattr(
        ce,
        "compute_ecm_orientation_lyapunov_metric",
        lambda *a, **k: calls.append("lyap") or (_ for _ in ()).throw(
            AssertionError("HB#5 must NOT be called when HB#3 raises")
        ),
    )

    ecm = _ecm()
    fa = _fa("a1")
    with pytest.raises(FAToECMScatteringError):
        step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert calls == []


def test_phase_e_v1_propagates_hb12_error_stops_before_bias(monkeypatch):
    """Test 7 (Y6): monkeypatch HB#1+#2 to raise; verify HB#4 / HB#5
    not called."""

    calls = []

    def _raise_orient(*a, **k):
        raise ECMConstitutiveResponseError(
            "dt_invalid", "test stub HB#1+#2 failure"
        )

    monkeypatch.setattr(ce, "step_ecm_orientation_response", _raise_orient)
    monkeypatch.setattr(
        ce,
        "compute_ecm_to_fa_bias_neutral",
        lambda *a, **k: calls.append("bias") or None,
    )
    monkeypatch.setattr(
        ce,
        "compute_ecm_orientation_lyapunov_metric",
        lambda *a, **k: calls.append("lyap") or None,
    )

    ecm = _ecm()
    fa = _fa("a1")
    with pytest.raises(ECMConstitutiveResponseError):
        step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert calls == []


def test_phase_e_v1_propagates_hb4_error_stops_before_lyap(monkeypatch):
    """Test 8 (Y6): monkeypatch HB#4 to raise; verify HB#5 not called."""

    calls = []

    def _raise_bias(*a, **k):
        raise FAToECMBiasError(
            "fa_bias_position_outside_ecm_grid", "test stub HB#4 failure"
        )

    monkeypatch.setattr(ce, "compute_ecm_to_fa_bias_neutral", _raise_bias)
    monkeypatch.setattr(
        ce,
        "compute_ecm_orientation_lyapunov_metric",
        lambda *a, **k: calls.append("lyap") or None,
    )

    ecm = _ecm()
    fa = _fa("a1")
    with pytest.raises(FAToECMBiasError):
        step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert calls == []


def test_phase_e_v1_propagates_hb5_error_at_end(monkeypatch):
    """Test 9 (Y6): monkeypatch HB#5 to raise; verify HB#3 / HB#1+#2 /
    HB#4 all called first."""

    calls = []

    real_scatter = ce.scatter_fa_traction_to_ecm_bilinear
    real_orient = ce.step_ecm_orientation_response
    real_bias = ce.compute_ecm_to_fa_bias_neutral

    def _spy_scatter(*a, **k):
        calls.append("scatter")
        return real_scatter(*a, **k)

    def _spy_orient(*a, **k):
        calls.append("orient")
        return real_orient(*a, **k)

    def _spy_bias(*a, **k):
        calls.append("bias")
        return real_bias(*a, **k)

    def _raise_lyap(*a, **k):
        calls.append("lyap_attempted")
        raise ValueError("test stub HB#5 failure")

    monkeypatch.setattr(ce, "scatter_fa_traction_to_ecm_bilinear", _spy_scatter)
    monkeypatch.setattr(ce, "step_ecm_orientation_response", _spy_orient)
    monkeypatch.setattr(ce, "compute_ecm_to_fa_bias_neutral", _spy_bias)
    monkeypatch.setattr(ce, "compute_ecm_orientation_lyapunov_metric", _raise_lyap)

    ecm = _ecm()
    fa = _fa("a1")
    with pytest.raises(ValueError):
        step_closed_loop_phase_e_v1((fa,), ecm, dt_s=600.0)
    assert calls == ["scatter", "orient", "bias", "lyap_attempted"]


def test_phase_e_v1_provides_item_1_response_monotonicity_evidence():
    """Test 10 (Y9): apply with two traction magnitudes; assert larger
    traction → larger convex_weight + faster orientation alignment.

    ECM-side evidence ONLY; not full closed-loop monotonicity."""

    ecm = _ecm()
    fa_small = _fa("small", traction=(0.5, 0.0))
    fa_large = _fa("large", traction=(5.0, 0.0))

    result_small = step_closed_loop_phase_e_v1((fa_small,), ecm, dt_s=600.0)
    result_large = step_closed_loop_phase_e_v1((fa_large,), ecm, dt_s=600.0)

    assert (
        result_large.orientation_response.diagnostics.max_convex_weight
        > result_small.orientation_response.diagnostics.max_convex_weight
    )

    initial_distance = math.sqrt(
        2.0
    )
    distance_small = float(
        np.linalg.norm(
            result_small.updated_ecm.orientation_tensor[1, 1]
            - np.array([[1.0, 0.0], [0.0, 0.0]])
        )
    )
    distance_large = float(
        np.linalg.norm(
            result_large.updated_ecm.orientation_tensor[1, 1]
            - np.array([[1.0, 0.0], [0.0, 0.0]])
        )
    )
    assert distance_large < distance_small <= initial_distance + 1e-9


def test_phase_e_v1_provides_item_2_saturation_evidence():
    """Test 11 (Y9): extreme traction → v_active ≤ v_active_max_bound
    (HB#1+#2 convex weight ∈ [0, 1] saturation built-in).

    ECM-side evidence ONLY."""

    ecm = _ecm()
    fa_extreme = _fa("ext", traction=(1000.0, 0.0))

    result = step_closed_loop_phase_e_v1((fa_extreme,), ecm, dt_s=3600.0 * 24)
    assert result.lyapunov_metric.v_active_um2 >= 0.0
    assert (
        result.lyapunov_metric.v_active_um2
        <= result.lyapunov_metric.diagnostics.v_active_max_bound_um2 + 1e-9
    )

    assert result.orientation_response.diagnostics.max_convex_weight <= 1.0


def test_phase_e_v1_provides_item_3_zero_traction_no_response_evidence():
    """Test 12 (Y9): zero traction → orientation unchanged + V_active
    = 0.

    ECM-side evidence ONLY."""

    ecm = _ecm()
    fa_zero = _fa("zero", traction=(0.0, 0.0), state="unbound", bound_fraction=0.0, position=(1.5, 1.5))

    result = step_closed_loop_phase_e_v1((fa_zero,), ecm, dt_s=600.0)

    np.testing.assert_array_equal(
        result.updated_ecm.orientation_tensor, ecm.orientation_tensor
    )
    assert result.lyapunov_metric.v_active_um2 == 0.0


def test_phase_e_v1_provides_item_4_ecm_side_boundedness_evidence():
    """Test 13 (Y9): apply 50 steps with rotating traction direction
    (θ uniform [0, 2π]); assert v_active_um2 ≤ v_active_max_bound_um2
    at every step AND multipliers_per_fa all-1.0 (neutral) at every
    step.

    ECM-side evidence ONLY, NOT closed-loop convergence (per locked
    Y9 wording)."""

    ecm = _ecm()
    state = ecm
    n_steps = 50

    for k in range(n_steps):
        theta = 2 * math.pi * k / n_steps
        fa = _fa(
            f"a_{k}",
            traction=(2.0 * math.cos(theta), 2.0 * math.sin(theta)),
        )
        result = step_closed_loop_phase_e_v1((fa,), state, dt_s=300.0)

        assert result.lyapunov_metric.v_active_um2 >= 0.0
        assert (
            result.lyapunov_metric.v_active_um2
            <= result.lyapunov_metric.diagnostics.v_active_max_bound_um2 + 1e-9
        )

        assert np.all(result.ecm_to_fa_bias.multipliers_per_fa == 1.0), (
            f"step {k}: HB#4 multipliers must remain neutral all-1.0 "
            f"in Phase E v1 (HB#4 hard-wired neutral)"
        )

        state = result.updated_ecm


def test_phase_e_v1_module_doc_does_not_overclaim_full_loop():
    """Test 14 (Y2 forward guard): module docstring + function
    docstring + return class docstring do NOT contain forbidden
    overclaim strings."""

    forbidden_phrases = [
        "Items 1-4 satisfied",
        "full closed-loop",
        "full FA→ECM→FA",
        "gate satisfied",
    ]

    module_doc = ce.__doc__ or ""
    function_doc = step_closed_loop_phase_e_v1.__doc__ or ""
    result_class_doc = PhaseEStepResult.__doc__ or ""

    for phrase in forbidden_phrases:
        assert phrase not in module_doc, (
            f"module docstring must not contain forbidden phrase {phrase!r}"
        )
        assert phrase not in function_doc, (
            f"function docstring must not contain forbidden phrase {phrase!r}"
        )
        assert phrase not in result_class_doc, (
            f"PhaseEStepResult docstring must not contain forbidden phrase "
            f"{phrase!r}"
        )


def test_phase_e_v1_module_does_not_reference_effective_stiffness():
    """Test 15 (Y13): both string-search and AST walk; the lock doc
    explains effective_stiffness in prose, but the .py file cannot."""

    src = inspect.getsource(ce)
    assert "effective_stiffness" not in src, (
        "closed_loop_phase_e.py source must not reference "
        "effective_stiffness anywhere (function body, docstring, "
        "comment, module docstring)"
    )

    tree = ast.parse(src)
    referenced = (
        {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
        | {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    )
    assert "effective_stiffness" not in referenced, (
        "closed_loop_phase_e.py AST must not have effective_stiffness "
        "as Attribute.attr or Name.id"
    )


def test_phase_e_v1_imports_hb12_constants_no_redefinition():
    """Test 16 (Y14): ce.TRACTION_REF_NN_PER_UM2 is
    cr.TRACTION_REF_NN_PER_UM2 (object identity); ce.K_ORIENT_PER_S
    is cr.K_ORIENT_PER_S; no surprise uppercase numeric constants."""

    assert ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2, (
        "TRACTION_REF_NN_PER_UM2 must be re-imported (object identity), "
        "not redefined"
    )
    assert ce.K_ORIENT_PER_S is cr.K_ORIENT_PER_S, (
        "K_ORIENT_PER_S must be re-imported (object identity), "
        "not redefined"
    )

    src = inspect.getsource(ce)
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    name = target.id
                    if (
                        name.isupper()
                        and "_" in name
                        and isinstance(node.value, ast.Constant)
                        and isinstance(node.value.value, (int, float))
                    ):
                        raise AssertionError(
                            f"closed_loop_phase_e.py defines uppercase "
                            f"numeric constant {name!r}; constants must be "
                            f"imported from sister modules per Y14"
                        )


def test_phase_e_v1_exports_through_both_init():
    """Test 17: PhaseEStepResult + step_closed_loop_phase_e_v1
    importable from BOTH acs.v2.dynamics AND acs.v2; identity match
    + __all__ inclusion."""

    expected_symbols = {"PhaseEStepResult", "step_closed_loop_phase_e_v1"}

    for sym in expected_symbols:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym!r}"
        assert hasattr(v2_dynamics, sym), f"acs.v2.dynamics missing {sym!r}"
        assert getattr(v2_pkg, sym) is getattr(ce, sym), (
            f"acs.v2.{sym} is not "
            f"acs.v2.dynamics.closed_loop_phase_e.{sym}"
        )
        assert getattr(v2_dynamics, sym) is getattr(ce, sym), (
            f"acs.v2.dynamics.{sym} is not "
            f"acs.v2.dynamics.closed_loop_phase_e.{sym}"
        )
        assert sym in v2_pkg.__all__, f"acs.v2.__all__ missing {sym!r}"
        assert sym in v2_dynamics.__all__, (
            f"acs.v2.dynamics.__all__ missing {sym!r}"
        )


def test_phase_e_v2_returns_phase_e_step_result_and_identity_invariant():
    """Phase E v2 reuses PhaseEStepResult and preserves the v1 Y4
    updated-ECM identity invariant."""

    ecm = _ecm(orientation_value=(1.0, 0.0, 0.0, -1.0))
    fa = _fa("a1", traction=(1.0, 0.0))

    result = step_closed_loop_phase_e_v2((fa,), ecm, dt_s=0.0, k_active=1.0)

    assert isinstance(result, PhaseEStepResult)
    assert result.updated_ecm is result.orientation_response.updated_ecm
    assert result.ecm_to_fa_bias.diagnostics_dict["mechanism"] == (
        "deviatoric_rayleigh_orientation"
    )
    assert result.ecm_to_fa_bias.diagnostics_dict["k_active"] == 1.0


def test_phase_e_v2_call_order_uses_active_bias_not_neutral(monkeypatch):
    """HB#3 -> HB#1+#2 -> HB#4-active -> HB#5, with HB#5 receiving
    the original scatter object."""

    calls: list[tuple[str, tuple, dict]] = []

    real_scatter = ce.scatter_fa_traction_to_ecm_bilinear
    real_orient = ce.step_ecm_orientation_response
    real_active = ce.compute_ecm_to_fa_bias_active
    real_lyap = ce.compute_ecm_orientation_lyapunov_metric

    def _spy_scatter(*args, **kwargs):
        out = real_scatter(*args, **kwargs)
        calls.append(("scatter", args, {"out": out}))
        return out

    def _spy_orient(*args, **kwargs):
        out = real_orient(*args, **kwargs)
        calls.append(("orient", args, {"out": out}))
        return out

    def _spy_active(*args, **kwargs):
        out = real_active(*args, **kwargs)
        calls.append(("active", args, {"kwargs": kwargs, "out": out}))
        return out

    def _neutral_should_not_run(*args, **kwargs):
        raise AssertionError("Phase E v2 must not call HB#4-neutral")

    def _spy_lyap(*args, **kwargs):
        out = real_lyap(*args, **kwargs)
        calls.append(("lyap", args, {"out": out}))
        return out

    monkeypatch.setattr(ce, "scatter_fa_traction_to_ecm_bilinear", _spy_scatter)
    monkeypatch.setattr(ce, "step_ecm_orientation_response", _spy_orient)
    monkeypatch.setattr(ce, "compute_ecm_to_fa_bias_active", _spy_active)
    monkeypatch.setattr(ce, "compute_ecm_to_fa_bias_neutral", _neutral_should_not_run)
    monkeypatch.setattr(ce, "compute_ecm_orientation_lyapunov_metric", _spy_lyap)

    ecm = _ecm(orientation_value=(1.0, 0.0, 0.0, -1.0))
    fa = _fa("a1", traction=(1.0, 0.0))
    result = step_closed_loop_phase_e_v2((fa,), ecm, dt_s=0.0, k_active=0.7)

    assert [c[0] for c in calls] == ["scatter", "orient", "active", "lyap"]
    scatter_out = calls[0][2]["out"]
    post_update_ecm = calls[1][2]["out"].updated_ecm

    active_args = calls[2][1]
    active_kwargs = calls[2][2]["kwargs"]
    assert active_args[1] is post_update_ecm
    assert active_kwargs["k_active"] == 0.7

    lyap_args = calls[3][1]
    assert lyap_args[0] is post_update_ecm
    assert lyap_args[1] is scatter_out
    assert result.updated_ecm is post_update_ecm


def test_phase_e_v2_invalid_k_active_raises_before_scatter(monkeypatch):
    """Valid ECM + invalid k_active must fail before HB#3 scatter."""

    ecm = _ecm()
    scatter_calls = {"count": 0}

    def _tracking_scatter(*args, **kwargs):
        scatter_calls["count"] += 1
        raise AssertionError("scatter should not run when k_active is invalid")

    monkeypatch.setattr(ce, "scatter_fa_traction_to_ecm_bilinear", _tracking_scatter)

    with pytest.raises(FAToECMBiasError) as excinfo:
        step_closed_loop_phase_e_v2((), ecm, dt_s=0.0, k_active=-1.0)

    assert excinfo.value.failure_kind == "k_active_invalid"
    assert scatter_calls["count"] == 0


def test_phase_e_v2_anti_collapse_under_anisotropic_ecm_with_two_perpendicular_fas():
    """Central v2 invariant: anisotropic ECM differentiates two FAs at
    the same position by traction direction."""

    ecm = _ecm(orientation_value=(1.0, 0.0, 0.0, -1.0))
    fas = (
        _fa("x", traction=(1.0, 0.0)),
        _fa("y", traction=(0.0, 1.0)),
    )

    result = step_closed_loop_phase_e_v2(fas, ecm, dt_s=0.0, k_active=1.0)
    multipliers = result.ecm_to_fa_bias.multipliers_per_fa

    np.testing.assert_allclose(multipliers[0, :], np.exp(1.0), rtol=1e-12)
    np.testing.assert_allclose(multipliers[1, :], np.exp(-1.0), rtol=1e-12)
    assert not np.allclose(multipliers, 1.0)
    assert not np.allclose(multipliers[0, :], multipliers[1, :])


def test_phase_e_v2_docstring_marks_composition_closure_without_item_overclaim():
    """Q4 hybrid wording: per-step composition closure only."""

    function_doc = step_closed_loop_phase_e_v2.__doc__ or ""

    assert "composition-structure" in function_doc
    assert "does **not** by itself establish Items 1-4" in function_doc
    for forbidden in (
        "Items 1-4 satisfied",
        "full closed-loop satisfied",
        "satisfies_item",
    ):
        assert forbidden not in function_doc


def test_phase_e_v2_exports_through_both_init():
    """Only the new v2 wrapper function is exported; result type is
    reused."""

    assert hasattr(v2_pkg, "step_closed_loop_phase_e_v2")
    assert hasattr(v2_dynamics, "step_closed_loop_phase_e_v2")
    assert v2_pkg.step_closed_loop_phase_e_v2 is ce.step_closed_loop_phase_e_v2
    assert v2_dynamics.step_closed_loop_phase_e_v2 is ce.step_closed_loop_phase_e_v2
    assert "step_closed_loop_phase_e_v2" in v2_pkg.__all__
    assert "step_closed_loop_phase_e_v2" in v2_dynamics.__all__

    assert not hasattr(v2_pkg, "PhaseEV2StepResult")
    assert not hasattr(v2_dynamics, "PhaseEV2StepResult")
    assert "PhaseEV2StepResult" not in v2_pkg.__all__
    assert "PhaseEV2StepResult" not in v2_dynamics.__all__
