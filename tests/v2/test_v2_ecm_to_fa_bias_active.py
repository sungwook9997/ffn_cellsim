"""Tests for ECM->FA active bias (HB#4-active step 1, Phase E v2 step 1).

Test catalog mirrors ``docs/v2/v2_hard_blocker_4_active_locked.md`` §4
(21 tests) and is cross-referenced in
``docs/v2/v2_hard_blocker_4_active_sanity_gate.md`` §5.

The catalog covers:

- Validation (5): k_active disciplines + ecm.validate() ordering
- Deviatoric form (3): isotropic neutral / schema-tight sqrt(2) / rank-1 half
- Multiplier behavior (3): aligned boost / orthogonal suppress / 3-rates-same
- Zero-traction + boundary (3): exact-neutral / count diagnostic / empty
- Diagnostics (3): exact 10-key set / no loose extras / sampler_geometry sister
- Failure-kind (1): active_multiplier_non_finite via NaN orientation_tensor
- Validation order (1): k_active raises before sampler invoked (monkeypatch)
- Effective-stiffness AST + string guard (1): function-scoped
- Exports + private constant (1): single public symbol; private hidden
"""

from __future__ import annotations

import ast
import inspect
import textwrap

import numpy as np
import pytest

from acs.v2.dynamics import ecm_to_fa_bias as ecm_to_fa_bias_mod
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMSampledAtFAs,
    ECMToFABiasResult,
    FAToECMBiasError,
    RATE_NAMES,
    compute_ecm_to_fa_bias_active,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


def _ecm_isotropic_half(
    nx: int = 4,
    ny: int = 4,
    *,
    spacing_um: float = 1.0,
    origin=(0.0, 0.0),
):
    """ECM with orientation_tensor = 0.5 * I (Item 5 IC).

    Under deviatoric form Q = T - 0.5*trace(T)*I:
        T = 0.5 * I  =>  trace(T) = 1  =>  Q = 0.5*I - 0.5*I = 0
    so alignment_score == 0 exactly for any traction direction.
    """
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 0.5
    orientation[..., 1, 1] = 0.5
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


def _ecm_uniform_orientation(
    t_xx: float,
    t_xy: float,
    t_yy: float,
    nx: int = 4,
    ny: int = 4,
    *,
    spacing_um: float = 1.0,
    origin=(0.0, 0.0),
):
    """ECM with uniform symmetric orientation tensor [[t_xx, t_xy], [t_xy, t_yy]]."""
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = t_xx
    orientation[..., 0, 1] = t_xy
    orientation[..., 1, 0] = t_xy
    orientation[..., 1, 1] = t_yy
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
    adhesion_id: str = "fa-0",
    *,
    position=(2.0, 2.0),
    cell_id: str = "cell-A",
    state: str = "mature",
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
# Validation (5)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_k_active_zero_or_negative_raises():
    ecm = _ecm_isotropic_half()
    fa = _fa()
    for bad in (0.0, -1.0, -1e-12):
        with pytest.raises(FAToECMBiasError) as excinfo:
            compute_ecm_to_fa_bias_active((fa,), ecm, k_active=bad)
        assert excinfo.value.failure_kind == "k_active_invalid"


def test_compute_ecm_to_fa_bias_active_k_active_nonfinite_raises():
    """Lock §4 #2 specified nonfinite cases (NaN/+inf/-inf); per Codex
    id=1625 BLOCKER fold-in directive (no test catalog increase),
    nonnumeric cases (str / None / list) are also folded into this
    bucket — all raise FAToECMBiasError(failure_kind="k_active_invalid")
    rather than leaking raw Python/NumPy TypeErrors."""
    ecm = _ecm_isotropic_half()
    fa = _fa()
    for bad in (float("nan"), float("inf"), -float("inf")):
        with pytest.raises(FAToECMBiasError) as excinfo:
            compute_ecm_to_fa_bias_active((fa,), ecm, k_active=bad)
        assert excinfo.value.failure_kind == "k_active_invalid"
    for bad_nonnumeric in ("1.0", None, [], (1.0,), {"x": 1.0}):
        with pytest.raises(FAToECMBiasError) as excinfo:
            compute_ecm_to_fa_bias_active((fa,), ecm, k_active=bad_nonnumeric)
        assert excinfo.value.failure_kind == "k_active_invalid"
        assert "numeric" in str(excinfo.value).lower() or "bool" in str(excinfo.value).lower()


def test_compute_ecm_to_fa_bias_active_k_active_bool_raises():
    ecm = _ecm_isotropic_half()
    fa = _fa()
    with pytest.raises(FAToECMBiasError) as excinfo:
        compute_ecm_to_fa_bias_active((fa,), ecm, k_active=True)
    assert excinfo.value.failure_kind == "k_active_invalid"
    assert "bool" in str(excinfo.value).lower()


def test_compute_ecm_to_fa_bias_active_k_active_too_large_raises_non_finite_bound():
    ecm = _ecm_isotropic_half()
    fa = _fa()
    # Threshold derivation: exp(k * sqrt(2)) > float64.max iff
    #   k > log(float64.max) / sqrt(2) ~= 709.78 / 1.4142 ~= 501.892.
    # Use 1000.0 to be safely above. (Codex id=1618 BLOCKER fix.)
    with pytest.raises(FAToECMBiasError) as excinfo:
        compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1000.0)
    assert excinfo.value.failure_kind == "k_active_invalid"


def test_compute_ecm_to_fa_bias_active_invalid_ecm_raises_at_validate():
    bad_orientation = np.zeros((4, 4, 2, 2), dtype=np.float64)
    bad_orientation[..., 0, 0] = 5.0  # violates |T_ij| <= 1
    bad_ecm = ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.zeros((4, 4), dtype=np.float64),
        ligand_density=np.zeros((4, 4), dtype=np.float64),
        fiber_density=np.zeros((4, 4), dtype=np.float64),
        orientation_tensor=bad_orientation,
        accumulated_traction_nNs_per_um2=np.zeros((4, 4), dtype=np.float64),
        source="synthetic",
    )
    fa = _fa()
    with pytest.raises(ValueError):
        compute_ecm_to_fa_bias_active((fa,), bad_ecm, k_active=1.0)


# ---------------------------------------------------------------------------
# Deviatoric form (3, Y1+Y10)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_isotropic_orientation_returns_neutral_multiplier():
    """Item 5 IC T = 0.5*I -> Q = 0 -> score = 0 -> multiplier == 1.0 exact
    for any traction direction. Y1 anti-pattern guard: raw Rayleigh
    n.T @ T @ n would give 0.5 here (false alignment claim)."""
    ecm = _ecm_isotropic_half()
    for traction in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (-1.0, 0.5), (3.0, -2.0)]:
        fa = _fa(traction=traction)
        result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.5)
        assert result.multipliers_per_fa.shape == (1, 3)
        assert result.multipliers_per_fa[0, 0] == 1.0
        assert result.multipliers_per_fa[0, 1] == 1.0
        assert result.multipliers_per_fa[0, 2] == 1.0
        assert result.diagnostics_dict["max_alignment_score"] == 0.0
        assert result.diagnostics_dict["min_alignment_score"] == 0.0


def test_compute_ecm_to_fa_bias_active_schema_tight_bound_score_approaches_sqrt2():
    """T = [[1, 1], [1, -1]] is already traceless (trace = 0) so Q = T.
    Eigenvalues from lambda^2 - 2 = 0: lambda = +/- sqrt(2). Eigenvector
    for +sqrt(2) at angle theta = pi/8 (since tan(2*theta) = 1/1).
    Y10 corrected math (id=1611 / id=1612 SEAL).
    """
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=1.0, t_yy=-1.0)
    theta = np.pi / 8.0
    nx_dir = float(np.cos(theta))
    ny_dir = float(np.sin(theta))
    fa = _fa(traction=(nx_dir, ny_dir))
    k_active = 0.7
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=k_active)

    expected_score = float(np.sqrt(2.0))
    expected_multiplier = float(np.exp(k_active * expected_score))
    assert result.multipliers_per_fa[0, 0] == pytest.approx(expected_multiplier, rel=1e-10)
    assert result.diagnostics_dict["max_alignment_score"] == pytest.approx(
        expected_score, rel=1e-10
    )
    assert result.diagnostics_dict["score_bound"] == pytest.approx(
        float(np.sqrt(2.0)), rel=1e-15
    )


def test_compute_ecm_to_fa_bias_active_rank1_aligned_score_half():
    """Rank-1 aligned T = n_outer @ n_outer.T with traction along same n.
    For n = (1, 0): T = [[1, 0], [0, 0]], trace = 1, Q = [[0.5, 0], [0, -0.5]],
    score = n.T @ Q @ n = 0.5. Y1 verification: rank-1 minus isotropic-half.
    """
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.0, t_yy=0.0)
    fa = _fa(traction=(1.0, 0.0))
    k_active = 1.0
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=k_active)
    expected_score = 0.5
    expected_multiplier = float(np.exp(k_active * expected_score))
    assert result.diagnostics_dict["max_alignment_score"] == pytest.approx(
        expected_score, rel=1e-12
    )
    assert result.multipliers_per_fa[0, 0] == pytest.approx(expected_multiplier, rel=1e-12)


# ---------------------------------------------------------------------------
# Multiplier behavior (3)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_aligned_target_boosts_above_neutral():
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.0, t_yy=0.0)
    fa = _fa(traction=(1.0, 0.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.0)
    assert result.multipliers_per_fa[0, 0] > 1.0


def test_compute_ecm_to_fa_bias_active_orthogonal_target_suppresses_below_neutral():
    # T = [[1, 0], [0, 0]] but FA traction along (0, 1) -> orthogonal to T's principal
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.0, t_yy=0.0)
    fa = _fa(traction=(0.0, 1.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.0)
    # n.T @ Q @ n with n=(0,1), Q=[[0.5,0],[0,-0.5]] -> -0.5
    assert result.multipliers_per_fa[0, 0] < 1.0
    assert result.diagnostics_dict["min_alignment_score"] == pytest.approx(-0.5, rel=1e-12)


def test_compute_ecm_to_fa_bias_active_all_three_rates_same_multiplier_per_fa():
    """Y4 scaffolding contract: all 3 rates per FA receive same multiplier."""
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.5, t_yy=-0.3)
    fas = (
        _fa("fa-0", position=(1.5, 1.5), traction=(1.0, 0.0)),
        _fa("fa-1", position=(2.5, 2.5), traction=(0.5, 0.5)),
        _fa("fa-2", position=(2.5, 1.5), traction=(0.0, 1.0)),
    )
    result = compute_ecm_to_fa_bias_active(fas, ecm, k_active=0.8)
    assert result.multipliers_per_fa.shape == (3, 3)
    for i in range(3):
        m0 = result.multipliers_per_fa[i, 0]
        m1 = result.multipliers_per_fa[i, 1]
        m2 = result.multipliers_per_fa[i, 2]
        assert m0 == m1 == m2


# ---------------------------------------------------------------------------
# Zero-traction + boundary (3, Y7)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_zero_traction_returns_exact_neutral():
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.5, t_yy=-0.5)
    fa = _fa(state="unbound", bound_fraction=0.0, traction=(0.0, 0.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=2.0)
    assert result.multipliers_per_fa[0, 0] == 1.0
    assert result.multipliers_per_fa[0, 1] == 1.0
    assert result.multipliers_per_fa[0, 2] == 1.0
    assert result.diagnostics_dict["max_alignment_score"] == 0.0
    assert result.diagnostics_dict["min_alignment_score"] == 0.0


def test_compute_ecm_to_fa_bias_active_zero_traction_count_in_diagnostics():
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.0, t_yy=0.0)
    fas = (
        _fa("fa-0", position=(1.5, 1.5), traction=(1.0, 0.0)),
        _fa("fa-1", position=(2.5, 1.5), state="unbound", bound_fraction=0.0,
            traction=(0.0, 0.0)),
        _fa("fa-2", position=(2.5, 2.5), state="unbound", bound_fraction=0.0,
            traction=(0.0, 0.0)),
    )
    result = compute_ecm_to_fa_bias_active(fas, ecm, k_active=1.0)
    assert result.diagnostics_dict["zero_traction_count"] == 2


def test_compute_ecm_to_fa_bias_active_empty_adhesions_returns_empty_consistent_shape():
    ecm = _ecm_isotropic_half()
    result = compute_ecm_to_fa_bias_active((), ecm, k_active=1.0)
    assert isinstance(result, ECMToFABiasResult)
    assert result.multipliers_per_fa.shape == (0, 3)
    assert result.fa_ids == ()
    assert result.sampled_diagnostics is not None
    assert isinstance(result.sampled_diagnostics, ECMSampledAtFAs)
    assert result.diagnostics_dict["n_adhesions"] == 0
    assert result.diagnostics_dict["max_multiplier"] == 1.0
    assert result.diagnostics_dict["min_multiplier"] == 1.0


# ---------------------------------------------------------------------------
# Diagnostics (3, Y6+Y13)
# ---------------------------------------------------------------------------


_EXPECTED_DIAGNOSTICS_KEYS = frozenset(
    {
        "n_adhesions",
        "max_multiplier",
        "min_multiplier",
        "sampler_geometry",
        "mechanism",
        "k_active",
        "max_alignment_score",
        "min_alignment_score",
        "zero_traction_count",
        "score_bound",
    }
)


def test_compute_ecm_to_fa_bias_active_diagnostics_dict_has_exact_10_keys():
    ecm = _ecm_uniform_orientation(t_xx=1.0, t_xy=0.0, t_yy=0.0)
    fa = _fa(traction=(1.0, 0.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.0)
    assert set(result.diagnostics_dict.keys()) == _EXPECTED_DIAGNOSTICS_KEYS
    assert len(result.diagnostics_dict) == 10


def test_compute_ecm_to_fa_bias_active_diagnostics_dict_no_loose_extra_keys():
    ecm = _ecm_isotropic_half()
    fa = _fa(traction=(1.0, 0.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.5)
    extra = set(result.diagnostics_dict.keys()) - _EXPECTED_DIAGNOSTICS_KEYS
    missing = _EXPECTED_DIAGNOSTICS_KEYS - set(result.diagnostics_dict.keys())
    assert extra == set(), f"unexpected diagnostic keys: {sorted(extra)}"
    assert missing == set(), f"missing diagnostic keys: {sorted(missing)}"
    assert result.diagnostics_dict["mechanism"] == "deviatoric_rayleigh_orientation"


def test_compute_ecm_to_fa_bias_active_diagnostics_sampler_geometry_matches_v1():
    """Y13: sampler_geometry string must match HB#4 v1 exactly."""
    ecm = _ecm_isotropic_half()
    fa = _fa(traction=(1.0, 0.0))
    result = compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.0)
    assert result.diagnostics_dict["sampler_geometry"] == "bilinear_cell_centered"


# ---------------------------------------------------------------------------
# Failure-kind: active_multiplier_non_finite (1)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_multiplier_non_finite_failure_kind(monkeypatch):
    """Defense in depth (Y11): even if ecm.validate() is bypassed (here via
    sampler monkeypatch returning corrupted orientation_tensor), the active
    function must raise active_multiplier_non_finite from the per-FA finite
    guard before returning."""
    ecm = _ecm_isotropic_half()
    fa = _fa(traction=(1.0, 0.0))

    def _corrupted_sampler(adhesions, _ecm):
        n = len(adhesions)
        bad_orientation = np.full((n, 2, 2), np.nan, dtype=np.float64)
        return ECMSampledAtFAs(
            fa_ids=tuple(a.adhesion_id for a in adhesions),
            stiffness_kpa=np.zeros(n, dtype=np.float64),
            fiber_density=np.zeros(n, dtype=np.float64),
            ligand_density=np.zeros(n, dtype=np.float64),
            orientation_tensor=bad_orientation,
            accumulated_traction_nNs_per_um2=np.zeros(n, dtype=np.float64),
        )

    monkeypatch.setattr(ecm_to_fa_bias_mod, "sample_ecm_at_fa_positions", _corrupted_sampler)
    with pytest.raises(FAToECMBiasError) as excinfo:
        compute_ecm_to_fa_bias_active((fa,), ecm, k_active=1.0)
    assert excinfo.value.failure_kind == "active_multiplier_non_finite"


# ---------------------------------------------------------------------------
# Validation order: k_active before sampler (1, Y12)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_invalid_k_active_raises_before_sampling(monkeypatch):
    """Y12: validation order is ecm.validate() -> _validate_k_active -> sampler.
    Bad k_active must raise before sampler is invoked, even with empty adhesions.
    """
    ecm = _ecm_isotropic_half()
    sampler_calls = {"count": 0}

    def _tracking_sampler(adhesions, _ecm):
        sampler_calls["count"] += 1
        raise AssertionError("sampler should not be called when k_active is invalid")

    monkeypatch.setattr(ecm_to_fa_bias_mod, "sample_ecm_at_fa_positions", _tracking_sampler)
    with pytest.raises(FAToECMBiasError) as excinfo:
        compute_ecm_to_fa_bias_active((), ecm, k_active=-1.0)
    assert excinfo.value.failure_kind == "k_active_invalid"
    assert sampler_calls["count"] == 0


# ---------------------------------------------------------------------------
# Effective-stiffness function-scoped AST + string guard (1, Y8)
# ---------------------------------------------------------------------------


_FORBIDDEN_FIELD_NAMES = frozenset(
    {
        "effective_stiffness",
        "stiffness_kpa",
        "fiber_density",
        "ligand_density",
        "accumulated_traction_nNs_per_um2",
    }
)


def test_compute_ecm_to_fa_bias_active_does_not_reference_effective_stiffness():
    """Y8: function-scoped guard. The active function body must not
    reference any forbidden mechanosensing field name (AST attribute /
    name access) and must not contain any forbidden literal substring.
    HB#5 Y11 sister-pattern.

    This is function-scoped (NOT module-scoped) because the surrounding
    module legitimately references some of these names in HB#4 v1
    docstrings, the ECMSampledAtFAs schema definition, and the shared
    sampler's sampled output.
    """
    src = textwrap.dedent(inspect.getsource(compute_ecm_to_fa_bias_active))
    tree = ast.parse(src)

    forbidden_attr_hits: list[str] = []
    forbidden_name_hits: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in _FORBIDDEN_FIELD_NAMES:
            forbidden_attr_hits.append(node.attr)
        if isinstance(node, ast.Name) and node.id in _FORBIDDEN_FIELD_NAMES:
            forbidden_name_hits.append(node.id)

    assert forbidden_attr_hits == [], (
        f"compute_ecm_to_fa_bias_active body references forbidden attributes: "
        f"{forbidden_attr_hits}"
    )
    assert forbidden_name_hits == [], (
        f"compute_ecm_to_fa_bias_active body references forbidden names: "
        f"{forbidden_name_hits}"
    )

    for forbidden_substr in _FORBIDDEN_FIELD_NAMES:
        assert forbidden_substr not in src, (
            f"compute_ecm_to_fa_bias_active source contains forbidden substring "
            f"{forbidden_substr!r}"
        )


# ---------------------------------------------------------------------------
# Exports + private constant (1, Y14)
# ---------------------------------------------------------------------------


def test_compute_ecm_to_fa_bias_active_only_function_exported():
    """Y14: single public symbol compute_ecm_to_fa_bias_active through both
    acs.v2.dynamics and acs.v2; private constant _DEVIATORIC_SCORE_BOUND
    and private helper _validate_k_active are NOT in either __all__."""
    import acs.v2 as v2_top
    import acs.v2.dynamics as v2_dynamics

    assert "compute_ecm_to_fa_bias_active" in v2_top.__all__
    assert "compute_ecm_to_fa_bias_active" in v2_dynamics.__all__
    assert hasattr(v2_top, "compute_ecm_to_fa_bias_active")
    assert hasattr(v2_dynamics, "compute_ecm_to_fa_bias_active")

    # Private constant + helper must NOT appear in either __all__
    assert "_DEVIATORIC_SCORE_BOUND" not in v2_top.__all__
    assert "_DEVIATORIC_SCORE_BOUND" not in v2_dynamics.__all__
    assert "_validate_k_active" not in v2_top.__all__
    assert "_validate_k_active" not in v2_dynamics.__all__

    # The constant lives at the module level of the implementation module
    # (private; importable via internal path but not exported through __all__).
    assert hasattr(ecm_to_fa_bias_mod, "_DEVIATORIC_SCORE_BOUND")
    assert ecm_to_fa_bias_mod._DEVIATORIC_SCORE_BOUND == float(np.sqrt(2.0))
