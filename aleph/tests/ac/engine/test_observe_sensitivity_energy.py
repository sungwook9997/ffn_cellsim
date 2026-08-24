"""Gates for the sensitivity / Fisher observables, the energy tests, and the run-record writer.

The reference systems are analytic: a power-law observable whose log-log Jacobian is known exactly, a
linear spring field whose loop work is exactly zero, and a deliberately non-conservative rotational
field whose loop work is exactly the enclosed circulation.  Testing against closed forms is what makes
these gates able to fail — an instrument checked only against itself cannot.
"""

from __future__ import annotations

import json
from typing import Mapping

import numpy as np
import pytest

from aleph.engine.contracts import (
    EvidenceLabel,
    EvidenceRung,
    GateVerdict,
    QuantitativeClaim,
    VoidCeiling,
)
from aleph.engine.observe.artifact import (
    ARTIFACT_SCHEMA,
    build_stamp,
    config_hash,
    observation_artifact,
    timing_block,
    write_artifact,
)
from aleph.engine.observe.energy import (
    EnergyBalance,
    balance_convergence,
    closed_loop_work,
    dissipated_work,
    loop_convergence,
    paired_path_work,
    path_work,
    step_energy_balance,
)
from aleph.engine.observe.sensitivity import (
    CENTRAL_DIFFERENCE_RELATIVE_STEP,
    fisher_report,
    log_log_jacobian,
)


# ── sensitivity ──────────────────────────────────────────────────────────────────────────────────────
def _power_law(exponents: np.ndarray):
    """Return an observable ``o_i = Π_j θ_j^{a_ij}``, whose log-log Jacobian is exactly ``a``."""

    def observable(parameters: Mapping[str, float]) -> np.ndarray:
        theta = np.array([parameters["a"], parameters["b"], parameters["c"]], np.float64)
        return np.exp(exponents @ np.log(theta))

    return observable


def test_log_log_jacobian_recovers_known_exponents() -> None:
    exponents = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, -1.0], [0.5, 0.5, 0.0]])
    sensitivity = log_log_jacobian(
        _power_law(exponents), {"a": 2.0, "b": 3.0, "c": 5.0}, names=("a", "b", "c"))
    assert sensitivity.jacobian == pytest.approx(exponents, abs=1e-7)
    assert sensitivity.parameter_names == ("a", "b", "c")
    assert sensitivity.relative_step == pytest.approx(CENTRAL_DIFFERENCE_RELATIVE_STEP)


def test_an_unresponsive_parameter_has_a_zero_column() -> None:
    exponents = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    sensitivity = log_log_jacobian(
        _power_law(exponents), {"a": 1.5, "b": 4.0, "c": 9.0}, names=("a", "b", "c"))
    norms = dict(zip(sensitivity.parameter_names, sensitivity.column_norms))
    assert norms["b"] == pytest.approx(0.0, abs=1e-7)
    assert norms["c"] == pytest.approx(0.0, abs=1e-7)
    assert norms["a"] > 1.0


def test_non_positive_observable_is_refused_rather_than_clamped() -> None:
    def observable(parameters: Mapping[str, float]) -> np.ndarray:
        return np.array([parameters["a"], -1.0])

    with pytest.raises(ValueError, match="non-positive"):
        log_log_jacobian(observable, {"a": 1.0})


def test_zero_parameter_is_refused() -> None:
    with pytest.raises(ValueError, match="relative perturbation is undefined"):
        log_log_jacobian(lambda p: np.array([1.0]), {"a": 0.0})


def test_missing_requested_parameter_is_named() -> None:
    with pytest.raises(ValueError, match="missing the requested names"):
        log_log_jacobian(lambda p: np.array([1.0]), {"a": 1.0}, names=("a", "z"))


def test_fisher_of_confounded_parameters_has_one_stiff_direction() -> None:
    # Both parameters enter only through their product, so only the combination is identifiable: the
    # information matrix must have exactly one resolvable direction.
    def observable(parameters: Mapping[str, float]) -> np.ndarray:
        product = parameters["a"] * parameters["b"]
        return np.array([product, product**2, product**3])

    sensitivity = log_log_jacobian(observable, {"a": 2.0, "b": 3.0}, names=("a", "b"))
    fisher = fisher_report(sensitivity)
    assert fisher.n_resolvable_directions == 1
    assert fisher.effective_dimension == pytest.approx(1.0, abs=1e-6)
    weights = fisher.eigenvectors[:, 0]
    assert abs(abs(weights[0]) - abs(weights[1])) < 1e-6  # the stiff direction is a+b, equally weighted


def test_fisher_of_independent_parameters_has_full_effective_dimension() -> None:
    exponents = np.eye(3)
    sensitivity = log_log_jacobian(
        _power_law(exponents), {"a": 2.0, "b": 2.0, "c": 2.0}, names=("a", "b", "c"))
    fisher = fisher_report(sensitivity)
    assert fisher.n_resolvable_directions == 3
    assert fisher.effective_dimension == pytest.approx(3.0, rel=1e-6)
    assert fisher.sloppiness_decades == pytest.approx(0.0, abs=1e-6)


def test_fisher_artifact_fields_name_the_leading_combination() -> None:
    exponents = np.array([[3.0, 0.0], [3.0, 0.0]])
    sensitivity = log_log_jacobian(
        lambda p: np.exp(exponents @ np.log(np.array([p["a"], p["b"]]))),
        {"a": 2.0, "b": 5.0}, names=("a", "b"))
    fields = fisher_report(sensitivity).as_artifact_fields(n_directions=1)
    leading = fields["leading_directions"][0]["parameter_weights"]
    assert list(leading)[0] == "a"


# ── energy ───────────────────────────────────────────────────────────────────────────────────────────
def _linear_spring_force(stiffness: float):
    """A conservative field ``F = −k x`` (potential ``½k|x|²``); its loop work is exactly zero."""

    def force(positions: np.ndarray) -> np.ndarray:
        return -stiffness * np.asarray(positions, np.float64)

    return force


def _rotational_force(omega: float):
    """A NON-conservative field: ``F_x = −ω y``, ``F_y = +ω x`` on node 0. Its curl is ``2ω``."""

    def force(positions: np.ndarray) -> np.ndarray:
        out = np.zeros_like(np.asarray(positions, np.float64))
        out[0, 0] = -omega * positions[0, 1]
        out[0, 1] = omega * positions[0, 0]
        return out

    return force


def test_conservative_field_closes_its_loop_exactly() -> None:
    origin = np.zeros((4, 3))
    origin[:, 0] = np.arange(4, dtype=np.float64)
    rng = np.random.default_rng(0)
    result = closed_loop_work(
        _linear_spring_force(12.0), origin, rng.standard_normal((4, 3)),
        rng.standard_normal((4, 3)), amplitude_um=1e-2, n_segments=8)
    assert abs(result.work) <= result.round_off_floor
    assert result.at_round_off_floor is True
    assert result.exchanged > 0.0


def test_non_conservative_field_is_caught_with_the_right_magnitude() -> None:
    # Circuit in the (x, y) plane of node 0 with edge ``a``: the enclosed circulation of a field with
    # curl ``2ω`` is exactly ``2ω·a²``.
    origin = np.zeros((1, 3))
    u = np.zeros((1, 3))
    u[0, 0] = 1.0
    v = np.zeros((1, 3))
    v[0, 1] = 1.0
    amplitude = 0.25
    omega = 3.0
    result = closed_loop_work(
        _rotational_force(omega), origin, u, v, amplitude_um=amplitude, n_segments=32)
    assert result.at_round_off_floor is False
    assert result.work == pytest.approx(2.0 * omega * amplitude**2, rel=1e-9)


def test_loop_work_of_a_conservative_field_is_amplitude_invariant() -> None:
    origin = np.zeros((3, 3))
    rng = np.random.default_rng(1)
    u, v = rng.standard_normal((3, 3)), rng.standard_normal((3, 3))
    for amplitude in (1e-1, 1e-2, 1e-3):
        result = closed_loop_work(
            _linear_spring_force(5.0), origin, u, v, amplitude_um=amplitude, n_segments=4)
        assert abs(result.work) <= result.round_off_floor


def test_path_work_of_a_spring_matches_its_potential_difference() -> None:
    start = np.zeros((1, 3))
    end = np.zeros((1, 3))
    end[0, 0] = 2.0
    signed, unsigned = path_work(_linear_spring_force(4.0), np.stack((start, end)))
    assert signed == pytest.approx(-0.5 * 4.0 * 2.0**2)
    assert unsigned == pytest.approx(abs(signed))


def test_loop_rejects_a_degenerate_circuit() -> None:
    origin = np.zeros((2, 3))
    with pytest.raises(ValueError, match="zero norm"):
        closed_loop_work(
            _linear_spring_force(1.0), origin, np.zeros((2, 3)), np.ones((2, 3)),
            amplitude_um=1.0)
    with pytest.raises(ValueError, match="positive-finite"):
        closed_loop_work(
            _linear_spring_force(1.0), origin, np.ones((2, 3)), np.ones((2, 3)),
            amplitude_um=0.0)


def _rest_length_spring_network(stiffness: float, rest_length: float, pairs):
    """A conservative NONLINEAR field: the rest-length pair spring the engine actually launches.

    Nonlinearity is the point — the midpoint rule is EXACT for a linear field, so only a nonlinear
    conservative one exercises the case the convergence study exists to separate: a real truncation
    residual that must NOT be read as a broken adjoint pair.  This particular nonlinearity is chosen
    because ``F = k(|d| − r0)·d̂`` is the pair force every tangent family in this engine differentiates,
    so the reference and the subject share their nonlinear structure.
    """

    def force(positions: np.ndarray) -> np.ndarray:
        array = np.asarray(positions, np.float64)
        out = np.zeros_like(array)
        for a, b in pairs:
            delta = array[b] - array[a]
            length = float(np.linalg.norm(delta))
            pull = stiffness * (length - rest_length) * (delta / length)
            out[a] += pull
            out[b] -= pull
        return out

    return force


def _loops(force, origin, u, v, amplitude, segments):
    """Walk the three refinements the convergence study compares."""
    return dict(
        base=closed_loop_work(force, origin, u, v, amplitude_um=amplitude, n_segments=segments),
        half_amplitude=closed_loop_work(
            force, origin, u, v, amplitude_um=amplitude * 0.5, n_segments=segments),
        refined_segments=closed_loop_work(
            force, origin, u, v, amplitude_um=amplitude, n_segments=segments * 2),
    )


def test_a_nonlinear_conservative_field_is_identified_as_quadrature() -> None:
    rng = np.random.default_rng(3)
    origin = rng.standard_normal((5, 3))
    u, v = rng.standard_normal((5, 3)), rng.standard_normal((5, 3))
    field = _rest_length_spring_network(7.0, 0.9, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0)))
    report = loop_convergence(**_loops(field, origin, u, v, 0.2, 8))
    assert report.verdict == "quadrature"
    # The segment exponent is the sharp discriminator: the midpoint rule is exactly second order in N,
    # while a circulation does not depend on N at all.
    assert report.segment_exponent == pytest.approx(2.0, abs=1e-3)
    # The amplitude exponent is at least the predicted 3 (higher-order terms can raise it, and here do:
    # ~3.7).  What matters for the verdict is that it is nearer 3 than 2.
    assert report.amplitude_exponent > 3.0
    assert report.amplitude_distance_to_quadrature < report.amplitude_distance_to_circulation
    assert report.consistent_with_circulation is False


def test_a_rotational_field_is_identified_as_circulation() -> None:
    origin = np.zeros((1, 3))
    u = np.zeros((1, 3))
    u[0, 0] = 1.0
    v = np.zeros((1, 3))
    v[0, 1] = 1.0
    report = loop_convergence(**_loops(_rotational_force(2.0), origin, u, v, 0.3, 8))
    assert report.verdict == "circulation"
    assert report.amplitude_exponent == pytest.approx(2.0, abs=1e-9)
    assert report.segment_exponent == pytest.approx(0.0, abs=1e-9)
    assert report.consistent_with_quadrature is False


def test_convergence_study_refuses_circuits_that_differ_in_more_than_one_respect() -> None:
    origin = np.zeros((1, 3))
    u = np.zeros((1, 3))
    u[0, 0] = 1.0
    v = np.zeros((1, 3))
    v[0, 1] = 1.0
    good = _loops(_rotational_force(1.0), origin, u, v, 0.2, 8)
    with pytest.raises(ValueError, match="half the base amplitude"):
        loop_convergence(
            base=good["base"],
            half_amplitude=closed_loop_work(
                _rotational_force(1.0), origin, u, v, amplitude_um=0.19, n_segments=8),
            refined_segments=good["refined_segments"])
    with pytest.raises(ValueError, match="MORE segments"):
        loop_convergence(
            base=good["base"], half_amplitude=good["half_amplitude"],
            refined_segments=good["base"])


def test_energy_balance_forms_a_dimensionless_closure_ratio() -> None:
    balance = EnergyBalance.from_terms(
        delta_potential=2.0, dissipated=3.0, event_jump=1.0, active_input=6.0)
    assert balance.residual == pytest.approx(0.0)
    assert balance.closure_ratio == pytest.approx(0.0)
    leaking = EnergyBalance.from_terms(
        delta_potential=2.0, dissipated=3.0, event_jump=1.0, active_input=4.0)
    assert leaking.residual == pytest.approx(2.0)
    assert leaking.closure_ratio == pytest.approx(0.5)


def test_energy_balance_rejects_a_non_finite_term() -> None:
    with pytest.raises(ValueError, match="finite"):
        EnergyBalance.from_terms(
            delta_potential=np.inf, dissipated=0.0, event_jump=0.0, active_input=0.0)


# ── artifact ─────────────────────────────────────────────────────────────────────────────────────────
def _record(**overrides) -> dict:
    base = dict(
        run_label="unit",
        evidence=EvidenceLabel(
            rung=EvidenceRung.CUDA_UNIT, quantitative=QuantitativeClaim.BLOCKED,
            basis="a probe was applied on a resolved accelerator"),
        config={"k": 1.0},
        census={"n_nodes": 3},
        t0={"n_bound": 0},
        timing=timing_block(wall_seconds=2.0, physical_time_s=0.5, n_steps=5),
        measurements={"x": 1.0},
        device="unit-device",
    )
    base.update(overrides)
    return observation_artifact(**base)


def test_artifact_stamps_its_own_build() -> None:
    record = _record()
    assert record["schema"] == ARTIFACT_SCHEMA
    assert set(record["build"]) >= {"commit", "dirty", "branch"}
    assert record["config_sha256"].startswith("sha256:")


def test_artifact_carries_both_evidence_axes() -> None:
    record = _record()
    assert record["evidence"] == EvidenceRung.CUDA_UNIT.value
    assert record["quantitative_claim_status"] == QuantitativeClaim.BLOCKED.value
    assert record["evidence_basis"]


def test_artifact_refuses_a_free_string_rung() -> None:
    with pytest.raises(TypeError, match="EvidenceLabel"):
        _record(evidence="CONNECTED")


def test_artifact_requires_a_t0_row() -> None:
    with pytest.raises(ValueError, match="t0 is required"):
        _record(t0={})


def test_a_partially_specified_verdict_is_refused() -> None:
    with pytest.raises(ValueError, match="all of void_ceiling"):
        _record(gate_passed=True, residual=1.0)


def test_a_run_whose_residual_swamps_its_signal_reports_void_not_pass() -> None:
    record = _record(
        void_ceiling=VoidCeiling(0.05, "declared before the run"),
        gate_passed=True, residual=50.0, signal=100.0)
    assert record["gate"]["verdict"] == GateVerdict.VOID.value
    assert record["gate"]["criterion_met"] is True
    assert record["gate"]["residual_over_signal"] == pytest.approx(0.5)


def test_a_clean_run_reports_pass_with_its_ceiling_recorded() -> None:
    record = _record(
        void_ceiling=VoidCeiling(0.05, "declared before the run"),
        gate_passed=True, residual=1.0, signal=100.0)
    assert record["gate"]["verdict"] == GateVerdict.PASS.value
    assert record["gate"]["void_ceiling_rationale"] == "declared before the run"


def test_config_hash_is_order_invariant_and_content_sensitive() -> None:
    assert config_hash({"a": 1, "b": 2}) == config_hash({"b": 2, "a": 1})
    assert config_hash({"a": 1}) != config_hash({"a": 2})


def test_build_stamp_never_fabricates_a_commit(tmp_path) -> None:
    stamp = build_stamp(tmp_path)
    assert stamp["commit"] == "unknown"
    assert stamp["source"] == "unavailable"
    assert "reason" in stamp


def test_a_declared_commit_is_recorded_as_declared_not_as_verified(tmp_path) -> None:
    # The native machine's tree is synced, not checked out, so a git-only stamp silently degrades to
    # "unknown" there — which is the missing-build defect the record exists to close.
    stamp = build_stamp(tmp_path, declared_commit="deadbeef")
    assert stamp["commit"] == "deadbeef"
    assert stamp["source"] == "declared"
    assert "not verified here" in stamp["reason"]


def test_artifact_carries_a_declared_commit_through(tmp_path) -> None:
    record = _record(declared_commit="cafebabe")
    assert record["build"]["commit"] in ("cafebabe",) or record["build"]["source"] == "git"


def test_write_artifact_round_trips(tmp_path) -> None:
    destination = write_artifact(tmp_path / "nested" / "record.json", _record())
    assert destination.exists()
    assert json.loads(destination.read_text())["schema"] == ARTIFACT_SCHEMA


# ── The step energy ledger ────────────────────────────────────────────────────────────────────────────
#
# The reference is a REST-LENGTH SPRING NETWORK, deliberately, and not a central field like `−k|x|²x`:
# a central field's symmetry drives the loop and path residuals to zero by cancellation, so it would
# exercise none of the quadrature these gates exist to bound.  A rest-length network is the nonlinearity
# the engine actually differentiates, and it makes the predicted exponents appear cleanly.

_SPRING_EDGES = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 3))
_SPRING_K = 1.7
_SPRING_R0 = 0.9


def _spring_force(positions: np.ndarray) -> np.ndarray:
    """Conservative rest-length spring network — the passive channel of the reference system."""
    force = np.zeros_like(positions)
    for a, b in _SPRING_EDGES:
        delta = positions[b] - positions[a]
        length = float(np.linalg.norm(delta))
        pull = _SPRING_K * (length - _SPRING_R0) * (delta / length)
        force[a] += pull
        force[b] -= pull
    return force


def _spring_potential(positions: np.ndarray) -> float:
    """The closed form the path integral must reproduce — the gate's independent reference."""
    return float(sum(
        0.5 * _SPRING_K * (float(np.linalg.norm(positions[b] - positions[a])) - _SPRING_R0) ** 2
        for a, b in _SPRING_EDGES
    ))


def _active_force(positions: np.ndarray) -> np.ndarray:
    """A non-conservative adjoint PAIR: rotational, so it has curl, and momentum-conserving."""
    force = np.zeros_like(positions)
    delta = positions[4] - positions[0]
    swirl = 0.3 * np.array([-delta[1], delta[0], 0.0])
    force[0] += swirl
    force[4] -= swirl
    return force


def _overdamped_trajectory(
    start: np.ndarray, *, mobility: float, iterations: int, total: object
) -> tuple[np.ndarray, float]:
    """Walk explicit-Euler overdamped dynamics, returning the full trajectory and the exact Σ|Δx|²."""
    position = np.array(start, np.float64)
    path = [position.copy()]
    squared = 0.0
    for _ in range(iterations):
        previous = position.copy()
        position = position + mobility * total(position)
        squared += float(np.sum((position - previous) ** 2))
        path.append(position.copy())
    return np.asarray(path, np.float64), squared


def test_path_work_recovers_the_closed_form_potential_difference() -> None:
    """``−∫F·dx`` along a wiggly path equals ``U(end) − U(start)`` for a conservative field."""
    rng = np.random.default_rng(0)
    start = rng.standard_normal((6, 3))
    end = start + 0.03 * rng.standard_normal((6, 3))
    steps = 33
    path = np.stack([
        start + (end - start) * (i / (steps - 1))
        + 0.004 * np.sin(np.pi * i / (steps - 1)) * rng.standard_normal((6, 3))
        for i in range(steps)
    ])
    path[0], path[-1] = start, end

    work, exchanged = path_work(_spring_force, path)
    exact = _spring_potential(end) - _spring_potential(start)
    assert exchanged > abs(work)                                   # traffic exceeds the net change
    assert abs(-work - exact) < 1e-4 * abs(exact)


def test_paired_path_work_evaluates_both_fields_at_the_same_midpoints() -> None:
    """A difference of integrals taken at different quadrature points would measure the quadrature."""
    seen: dict[str, list[np.ndarray]] = {"a": [], "b": []}

    def record(key: str):
        def force(positions: np.ndarray) -> np.ndarray:
            seen[key].append(np.array(positions, np.float64))
            return np.zeros_like(positions)
        return force

    rng = np.random.default_rng(1)
    path = rng.standard_normal((7, 4, 3))
    paired_path_work(record("a"), record("b"), path)
    assert len(seen["a"]) == len(seen["b"]) == path.shape[0] - 1
    for left, right in zip(seen["a"], seen["b"], strict=True):
        assert np.array_equal(left, right)


def test_the_step_energy_balance_closes_and_its_residual_is_the_integrator() -> None:
    """``ΔU + D − W_active`` is small, and halving the mobility halves it — first order, not a leak."""
    rng = np.random.default_rng(1)
    start = rng.standard_normal((6, 3)) * 0.5

    def total(positions: np.ndarray) -> np.ndarray:
        return _spring_force(positions) + _active_force(positions)

    ledgers = []
    for mobility, iterations in ((0.02, 400), (0.01, 800)):
        path, squared = _overdamped_trajectory(
            start, mobility=mobility, iterations=iterations, total=total)
        ledgers.append(step_energy_balance(
            path=path, force_total=total, force_passive=_spring_force,
            squared_displacement_um2=squared, mobility_um_per_pN=mobility))

    coarse, fine = ledgers
    # every term is nonzero, so the closure is not passing by everything being zero
    assert coarse.balance.delta_potential < 0.0 < coarse.balance.dissipated
    assert coarse.balance.active_input != 0.0
    assert coarse.balance.closure_ratio < 0.05
    # ΔU is path-independent, which is what licenses calling it a potential difference
    assert coarse.path_independence_ratio < 1e-3

    verdict = balance_convergence(base=coarse, refined=fine)
    assert verdict.verdict == "integrator_consistency"
    assert abs(verdict.residual_exponent - 1.0) < 0.1


def test_a_double_counted_force_channel_is_reported_as_a_LEAK() -> None:
    """The discriminator has to fire the other way too, or it only ever says one thing.

    The trajectory is driven by the real network, but the ledger's force evaluation counts three of its
    bonds TWICE — the "same physics twice" trap that all six of the PI framework's documented double
    counts are instances of.  The spurious work is a line integral in its own right, so refining the
    mobility does not move it, and the verdict must say so.
    """
    rng = np.random.default_rng(2)
    start = rng.standard_normal((6, 3)) * 0.5
    doubled_edges = _SPRING_EDGES + ((0, 3), (1, 2), (2, 3))

    def double_counted(positions: np.ndarray) -> np.ndarray:
        """The ledger's view: three bonds accumulated twice, as a duplicated launch would do."""
        force = np.zeros_like(positions)
        for a, b in doubled_edges:
            delta = positions[b] - positions[a]
            length = float(np.linalg.norm(delta))
            pull = _SPRING_K * (length - _SPRING_R0) * (delta / length)
            force[a] += pull
            force[b] -= pull
        return force

    blind = double_counted
    ledgers = []
    for mobility, iterations in ((0.02, 400), (0.01, 800)):
        path, squared = _overdamped_trajectory(
            start, mobility=mobility, iterations=iterations, total=_spring_force)
        ledgers.append(step_energy_balance(
            path=path, force_total=blind, force_passive=blind,
            squared_displacement_um2=squared, mobility_um_per_pN=mobility))

    coarse, fine = ledgers
    # Judged by the residual REFUSING to shrink when the mobility does — not by its size, which is what
    # makes this instrument different from a threshold on the residual.
    verdict = balance_convergence(base=coarse, refined=fine)
    assert verdict.verdict == "leak"
    assert abs(verdict.residual_exponent) < 0.3


def test_dissipated_work_is_dimensional_and_rejects_a_non_physical_mobility() -> None:
    assert dissipated_work(4.0, 2.0) == pytest.approx(2.0)
    with pytest.raises(ValueError, match="mobility"):
        dissipated_work(1.0, 0.0)
    with pytest.raises(ValueError, match="nonnegative"):
        dissipated_work(-1.0, 1.0)


def test_balance_convergence_refuses_a_comparison_that_measures_nothing() -> None:
    """Two ledgers at the same mobility differ by their relaxations, not by the integrator."""
    rng = np.random.default_rng(3)
    start = rng.standard_normal((6, 3)) * 0.5
    path, squared = _overdamped_trajectory(
        start, mobility=0.02, iterations=100, total=_spring_force)
    one = step_energy_balance(
        path=path, force_total=_spring_force, force_passive=_spring_force,
        squared_displacement_um2=squared, mobility_um_per_pN=0.02)
    with pytest.raises(ValueError, match="STRICTLY smaller mobility"):
        balance_convergence(base=one, refined=one)


# ── timing: a cost is a property of an engine AT A STATED CONVERGENCE ─────────────────────────────────
def test_timing_is_required_because_a_record_without_a_cost_answers_nothing() -> None:
    """The mirror of the t0 rule. This repo could not say whether anything got faster for a month."""
    with pytest.raises(ValueError, match="timing is required"):
        _record(timing={})
    with pytest.raises(ValueError, match="wall_seconds"):
        _record(timing={"n_steps": 5})


def test_timing_block_derives_the_rate_anyone_actually_asks_for() -> None:
    block = timing_block(wall_seconds=25.0, physical_time_s=0.2, n_steps=20,
                         n_inner_iterations=92_000)
    assert block["real_time_factor"] == pytest.approx(125.0)
    assert block["wall_seconds_per_step"] == pytest.approx(1.25)
    assert block["wall_seconds_per_inner_iteration"] == pytest.approx(25.0 / 92_000)


def test_timing_block_refuses_a_cost_that_cannot_be_one() -> None:
    with pytest.raises(ValueError, match="wall_seconds"):
        timing_block(wall_seconds=0.0)
    with pytest.raises(ValueError, match="physical_time_s"):
        timing_block(wall_seconds=1.0, physical_time_s=0.0)


def test_a_timing_is_marked_NON_comparable_unless_the_gate_passed() -> None:
    """The whole point: any solver is arbitrarily fast if allowed to stop early.

    The concrete failure this prevents already happened — `outer_step_wall_s = 1.13` sits in the same
    committed native artifact as `inner_converged: false`, and the two were read apart.
    """
    ceiling = VoidCeiling(0.01, "a residual above 1% of the signal makes the signal a transient")

    void = _record(void_ceiling=ceiling, gate_passed=False, residual=0.23, signal=1.0)
    assert void["gate"]["verdict"] == GateVerdict.VOID.value
    assert void["timing"]["comparable"] is False
    assert "VOID" in void["timing"]["not_comparable_reason"]

    no_gate = _record()                       # a driver that classified nothing
    assert no_gate["gate"] is None
    assert no_gate["timing"]["comparable"] is False
    assert "absent" in no_gate["timing"]["not_comparable_reason"]

    passed = _record(void_ceiling=ceiling, gate_passed=True, residual=1e-6, signal=1.0)
    assert passed["gate"]["verdict"] == GateVerdict.PASS.value
    assert passed["timing"]["comparable"] is True
    assert "same contract" in passed["timing"]["comparable_against"]


# ── the census may not claim more than it measures ───────────────────────────────────────────────────
def test_artifact_refuses_the_fraction_of_native_that_only_saw_the_cortex() -> None:
    """The banned key is banned by NAME, because its arithmetic was never the defect.

    Measured 2026-07-29: two runs whose `n_total_nodes` differ 1.73x and whose step cost differs 1.95x
    both stamped `fraction_of_native: 1.0`, since it was `cortex_filaments/70686` and the difference was
    the membrane. Every arithmetic check over 51 committed records passed; only the vocabulary was wrong.
    """
    with pytest.raises(ValueError, match="fraction_of_native"):
        _record(census={"n_nodes": 3, "fraction_of_native": 1.0})
    # the honest name is accepted, so the guard forces a rename rather than blocking the measurement
    assert _record(census={"n_nodes": 3, "fraction_of_native_cortex": 1.0})["census"]


def test_artifact_refuses_a_population_knob_whose_realized_size_is_unrecorded() -> None:
    """A knob that resizes a compartment must show up in the census, or a reduced run reads as full."""
    with pytest.raises(ValueError, match="membrane_vertices"):
        _record(config={"k": 1.0, "membrane_subdiv": 7}, census={"n_nodes": 3})
    ok = _record(config={"k": 1.0, "membrane_subdiv": 7},
                 census={"n_nodes": 3, "membrane_vertices": 163842})
    assert ok["census"]["membrane_vertices"] == 163842


def test_the_census_guard_is_not_vacuous() -> None:
    """A guard that cannot fail is not a guard — today's lesson, applied to today's guard.

    Both halves must be exercised: the banned key must actually raise, and a census with neither defect
    must actually pass. A test that only asserts the happy path would have accepted `fraction_of_native`.
    """
    assert _record(census={"n_nodes": 3})["census"] == {"n_nodes": 3}
    for bad_census, bad_config in (({"n_nodes": 3, "fraction_of_native": 0.5}, {"k": 1.0}),
                                   ({"n_nodes": 3}, {"k": 1.0, "membrane_subdiv": 8})):
        with pytest.raises(ValueError):
            _record(config=bad_config, census=bad_census)
