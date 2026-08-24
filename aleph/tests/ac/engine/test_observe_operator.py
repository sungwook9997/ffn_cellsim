"""Gates for the operator probe and the spectrum observables.

Every reference operator here is a pure-NumPy spring network built in this file — no Warp launch, no
device, and nothing imported from a physics module.  That is deliberate: the probe's contract is a
callable, so exercising it against an operator whose exact spectrum is known by hand tests the
INSTRUMENT rather than re-testing the physics the instrument will be pointed at.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.observe.operator_probe import assemble_dense_operator, symmetry_report
from aleph.engine.observe.spectrum import (
    lanczos_extremal_spectrum,
    lanczos_ritz,
    mode_localization,
    stiffness_spectrum,
)


def _chain_stiffness(n_nodes: int, k: float) -> np.ndarray:
    """Return the exact tangent of a 1-D chain of axial springs along x, free at both ends."""
    matrix = np.zeros((3 * n_nodes, 3 * n_nodes), np.float64)
    for node in range(n_nodes - 1):
        a, b = 3 * node, 3 * (node + 1)
        matrix[a, a] += k
        matrix[b, b] += k
        matrix[a, b] -= k
        matrix[b, a] -= k
    return matrix


def _matvec_of(matrix: np.ndarray, n_nodes: int):
    """Wrap a dense matrix as the ``(n, 3) -> (n, 3)`` callable the probe consumes."""

    def matvec(vector: np.ndarray) -> np.ndarray:
        return (matrix @ np.asarray(vector, np.float64).reshape(-1)).reshape(n_nodes, 3)

    return matvec


def test_column_probing_recovers_the_operator_exactly() -> None:
    reference = _chain_stiffness(6, 3.5)
    probed = assemble_dense_operator(_matvec_of(reference, 6), 6)
    assert probed.shape == (18, 18)
    assert np.array_equal(probed, reference)


def test_probe_refuses_a_population_it_cannot_afford() -> None:
    with pytest.raises(ValueError, match="SLICE instrument"):
        assemble_dense_operator(_matvec_of(np.zeros((9, 9)), 3), 3, max_dof=4)


def test_probe_rejects_a_matvec_that_lies_about_its_shape() -> None:
    def wrong(vector: np.ndarray) -> np.ndarray:
        return np.zeros((2, 3))

    with pytest.raises(ValueError, match="expected"):
        assemble_dense_operator(wrong, 4)


def test_probe_rejects_a_non_finite_response() -> None:
    def broken(vector: np.ndarray) -> np.ndarray:
        return np.full((4, 3), np.nan)

    with pytest.raises(ValueError, match="non-finite"):
        assemble_dense_operator(broken, 4)


def test_symmetric_operator_is_reported_conservative() -> None:
    report = symmetry_report(_chain_stiffness(5, 2.0), contribution_count=8)
    assert report.max_abs_asymmetry == 0.0
    assert report.relative_asymmetry == 0.0
    assert report.conservative is True
    assert report.round_off_floor > 0.0


def test_asymmetric_operator_is_reported_non_conservative() -> None:
    matrix = _chain_stiffness(5, 2.0)
    matrix[0, 3] += 1.0  # one broken adjoint pair: no potential exists for this field
    report = symmetry_report(matrix, contribution_count=8)
    assert report.max_abs_asymmetry == pytest.approx(1.0)
    assert report.conservative is False


def test_symmetry_floor_scales_with_the_declared_contribution_count() -> None:
    matrix = _chain_stiffness(5, 2.0)
    few = symmetry_report(matrix, contribution_count=2)
    many = symmetry_report(matrix, contribution_count=10_000)
    assert many.round_off_floor > few.round_off_floor
    assert few.round_off_floor > 0.0


def test_symmetry_report_rejects_a_non_square_operator() -> None:
    with pytest.raises(ValueError, match="square"):
        symmetry_report(np.zeros((3, 4)), contribution_count=1)


def test_free_chain_null_space_matches_its_rigid_body_count() -> None:
    # A free chain of axial-only springs along x has one translation constrained per axis pair: the
    # only stiff directions are the n-1 axial stretches, so the null space is 3n - (n-1).
    n_nodes = 6
    spectrum = stiffness_spectrum(_chain_stiffness(n_nodes, 4.0))
    assert spectrum.n_zero_modes == 3 * n_nodes - (n_nodes - 1)
    assert spectrum.n_negative_modes == 0
    assert spectrum.lambda_max > 0.0
    assert spectrum.lambda_min_nonzero is not None
    assert spectrum.condition_number == pytest.approx(
        spectrum.lambda_max / spectrum.lambda_min_nonzero)


def test_excess_null_space_over_the_rigid_body_bound_is_reported_not_absorbed() -> None:
    spectrum = stiffness_spectrum(_chain_stiffness(4, 1.0), rigid_body_zero_modes=6)
    assert spectrum.rigid_body_zero_modes == 6
    assert spectrum.null_space_at_least_rigid_body is True
    assert spectrum.n_excess_zero_modes == spectrum.n_zero_modes - 6
    assert spectrum.n_excess_zero_modes > 0  # the chain is floppy in every transverse direction


def test_fewer_zero_modes_than_free_rigid_bodies_is_a_hard_defect() -> None:
    # A fully constrained operator cannot be right for a system declared to hold two free rigid bodies.
    spectrum = stiffness_spectrum(np.eye(12), rigid_body_zero_modes=12)
    assert spectrum.n_zero_modes == 0
    assert spectrum.null_space_at_least_rigid_body is False
    assert spectrum.n_excess_zero_modes == -12


def test_stable_explicit_step_is_two_over_lambda_max() -> None:
    spectrum = stiffness_spectrum(_chain_stiffness(5, 7.0))
    assert spectrum.stable_explicit_mobility_step == pytest.approx(2.0 / spectrum.lambda_max)


def test_gershgorin_bound_is_computed_from_the_same_operator_and_never_undershoots() -> None:
    # Gershgorin bounds the spectrum, so the ratio can never fall below 1 — a value under 1 would mean
    # the bound and the eigenvalues were taken from different matrices, which is the drift this field
    # exists to make impossible.
    for stiffness in (1.0, 250.0, 4.6e5):
        spectrum = stiffness_spectrum(_chain_stiffness(7, stiffness))
        assert spectrum.gershgorin_bound >= spectrum.lambda_max
        assert spectrum.gershgorin_over_lambda_max >= 1.0
    # For a 1-D chain of equal springs the interior row sum is 4k while the true lambda_max is 4k too
    # only in the long-chain limit; the bound is therefore mildly conservative and never wrong.
    chain = stiffness_spectrum(_chain_stiffness(7, 3.0))
    assert chain.gershgorin_bound == pytest.approx(4.0 * 3.0)


def test_negative_eigenvalues_are_counted_not_clipped() -> None:
    matrix = np.diag([-5.0, 1.0, 2.0, 3.0, 4.0, 6.0])
    spectrum = stiffness_spectrum(matrix)
    assert spectrum.n_negative_modes == 1
    assert spectrum.eigenvalues[0] == pytest.approx(-5.0)


def test_zero_tolerance_is_derived_from_the_operator_norm() -> None:
    small = stiffness_spectrum(_chain_stiffness(4, 1.0))
    large = stiffness_spectrum(_chain_stiffness(4, 1.0e6))
    assert large.zero_tolerance == pytest.approx(small.zero_tolerance * 1.0e6, rel=1e-9)


def test_spectrum_artifact_fields_carry_the_low_end_verbatim() -> None:
    spectrum = stiffness_spectrum(_chain_stiffness(5, 2.0))
    fields = spectrum.as_artifact_fields(n_report=4)
    assert len(fields["lowest_eigenvalues_pN_per_um"]) == 4
    assert fields["lowest_eigenvalues_pN_per_um"] == [
        pytest.approx(float(v)) for v in spectrum.eigenvalues[:4]
    ]
    assert fields["n_dof"] == 15


def test_mode_localization_separates_a_local_mode_from_a_spread_one() -> None:
    n_nodes = 8
    positions = np.zeros((n_nodes, 3))
    positions[:, 0] = np.arange(n_nodes, dtype=np.float64)
    values = np.arange(3 * n_nodes, dtype=np.float64) + 1.0
    vectors = np.zeros((3 * n_nodes, 3 * n_nodes))
    vectors[0, 0] = 1.0  # a single-node mode
    spread = np.zeros(3 * n_nodes)
    spread[0::3] = 1.0 / np.sqrt(n_nodes)  # every node, x only
    vectors[:, 1] = spread
    local, delocal = mode_localization(values, vectors, positions, indices=(0, 1))
    assert local.participation_ratio == pytest.approx(1.0)
    assert local.gyration_radius_um == pytest.approx(np.std(positions[:, 0]) * 0.0 + 0.0, abs=1e-12)
    assert delocal.participation_ratio == pytest.approx(float(n_nodes))
    assert delocal.gyration_radius_um > local.gyration_radius_um


def test_mode_localization_rejects_an_index_outside_the_spectrum() -> None:
    positions = np.zeros((2, 3))
    values = np.ones(6)
    vectors = np.eye(6)
    with pytest.raises(ValueError, match="outside the spectrum"):
        mode_localization(values, vectors, positions, indices=(6,))


# ── Lanczos: the only path to the operator at a population column probing refuses ─────────────────────
def _dense_matvec(matrix: np.ndarray, n_nodes: int):
    """Return a node-shaped matrix-free application of a dense reference operator."""

    def matvec(vector: np.ndarray) -> np.ndarray:
        return (matrix @ np.asarray(vector, np.float64).reshape(-1)).reshape(n_nodes, 3)

    return matvec


def _random_spd(n_dof: int, *, seed: int, floor: float = 1e-3) -> np.ndarray:
    rng = np.random.default_rng(seed)
    root = rng.standard_normal((n_dof, n_dof))
    matrix = root @ root.T / n_dof + np.eye(n_dof) * floor
    return 0.5 * (matrix + matrix.T)


def test_lanczos_lambda_max_matches_the_dense_eigensolve_to_round_off() -> None:
    """The large end is what an explicit step is sized by, and it must be exact, not indicative."""
    n_nodes = 64
    matrix = _random_spd(3 * n_nodes, seed=3)
    truth = float(np.linalg.eigvalsh(matrix)[-1])
    result = lanczos_ritz(_dense_matvec(matrix, n_nodes), n_nodes, n_iterations=80)
    value, bound = result.extreme_high
    assert value == pytest.approx(truth, rel=1e-12)
    assert bound < 1e-6 * truth                    # and it says so itself, without a reference to check
    assert result.reorthogonalized is True
    assert result.breakdown is False


def test_lanczos_reports_an_unresolved_end_instead_of_a_confident_wrong_number() -> None:
    """Folding converges slowly at the small end of an ill-conditioned operator — say so, do not hide it.

    This is the honest limit of a solve-free small end, and the record has to carry it: a condition
    number built from an unresolved Ritz value would be a property of the iteration count.
    """
    n_nodes = 64
    matrix = _random_spd(3 * n_nodes, seed=3)
    gershgorin = float(np.max(np.sum(np.abs(matrix), axis=1)))
    record = lanczos_extremal_spectrum(
        _dense_matvec(matrix, n_nodes), n_nodes, upper_bound=gershgorin, n_iterations=80)

    truth = np.linalg.eigvalsh(matrix)
    assert record["lambda_max_pN_per_um"] == pytest.approx(float(truth[-1]), rel=1e-12)
    assert record["lambda_max_resolved"] is True
    assert record["upper_bound_over_lambda_max"] > 1.0     # the bound is an over-estimate, as expected
    if not record["lambda_min_resolved"]:
        assert record["condition_number"] is None
        assert "not resolved" in record["condition_number_note"]


def test_lanczos_recovers_the_small_end_when_the_operator_permits_it() -> None:
    """On a well-conditioned operator the folded pass resolves λ_min, so the method is not vacuous."""
    n_nodes = 24
    matrix = _random_spd(3 * n_nodes, seed=5, floor=2.0)
    gershgorin = float(np.max(np.sum(np.abs(matrix), axis=1)))
    record = lanczos_extremal_spectrum(
        _dense_matvec(matrix, n_nodes), n_nodes, upper_bound=gershgorin, n_iterations=3 * n_nodes)
    truth = np.linalg.eigvalsh(matrix)
    assert record["lambda_min_pN_per_um"] == pytest.approx(float(truth[0]), rel=1e-6)
    assert record["lambda_min_resolved"] is True
    assert record["condition_number"] == pytest.approx(float(truth[-1] / truth[0]), rel=1e-6)


def test_folding_against_a_bound_below_lambda_max_raises_instead_of_returning_the_wrong_end() -> None:
    n_nodes = 16
    matrix = _random_spd(3 * n_nodes, seed=7)
    with pytest.raises(ValueError, match="below the measured lambda_max"):
        lanczos_extremal_spectrum(
            _dense_matvec(matrix, n_nodes), n_nodes, upper_bound=1e-6, n_iterations=20)


def test_lanczos_breaks_down_exactly_on_an_invariant_subspace() -> None:
    """A rank-deficient operator exhausts its Krylov space; that is exact convergence, and is recorded."""
    n_nodes = 8
    n_dof = 3 * n_nodes
    basis = np.linalg.qr(np.random.default_rng(11).standard_normal((n_dof, n_dof)))[0]
    matrix = basis[:, :2] @ np.diag([5.0, 2.0]) @ basis[:, :2].T
    result = lanczos_ritz(_dense_matvec(matrix, n_nodes), n_nodes, n_iterations=n_dof)
    assert result.breakdown is True
    # the reachable space is rank 2, so the recurrence exhausts it and stops far short of the request
    assert result.n_iterations <= 4
    assert float(result.ritz_values[-1]) == pytest.approx(5.0, rel=1e-10)


def test_lanczos_refuses_a_krylov_space_larger_than_the_space_it_lives_in() -> None:
    matrix = _random_spd(12, seed=13)
    with pytest.raises(ValueError, match="exceeds the operator width"):
        lanczos_ritz(_dense_matvec(matrix, 4), 4, n_iterations=13)


def test_dense_probing_still_refuses_the_population_lanczos_exists_for() -> None:
    """The two instruments have to stay complementary: the refusal is what motivates the iterative path."""
    with pytest.raises(ValueError, match="Lanczos"):
        assemble_dense_operator(lambda v: v, 100_000)


# ── The projected operator: what the engine's CG actually inverts is P K P, not K ─────────────────────
def _projected_reference(n_nodes: int, *, n_constraints: int, seed: int):
    """Build ``(P, K, P K P)`` for a random orthogonal projector and a random SPD ``K``.

    Mirrors the engine's inner operator shape: the inextensibility projector removes one direction per
    segment, so ``P K P`` carries an exact null space of that dimension that has nothing to do with the
    mechanics.
    """
    n_dof = 3 * n_nodes
    rng = np.random.default_rng(seed)
    constraints = np.linalg.qr(rng.standard_normal((n_dof, n_constraints)))[0]
    projector = np.eye(n_dof) - constraints @ constraints.T
    stiffness = _random_spd(n_dof, seed=seed + 1, floor=1.5)
    return projector, stiffness, projector @ stiffness @ projector


def test_projected_lanczos_measures_the_restriction_not_the_constraint_null_space() -> None:
    """``P K P``'s small end must be the softest MECHANICAL mode, never the projector's own zeros.

    The engine solves ``aI + P K P``; at ``a = 0`` that operator is singular by construction, so the
    number that decides whether an implicit solve is warranted is the conditioning of its RESTRICTION to
    ``range(P)``.  Reading the raw small end instead reports zero for every operator ever built, which
    is why the projector has to enter the folding rather than be cleaned up afterwards.
    """
    n_nodes, n_constraints = 24, 17
    projector, _, projected = _projected_reference(n_nodes, n_constraints=n_constraints, seed=21)
    n_dof = 3 * n_nodes
    truth = np.linalg.eigvalsh(projected)
    # the exact restriction spectrum: drop the n_constraints structural zeros
    restricted = truth[n_constraints:]

    apply_p = _dense_matvec(projector, n_nodes)
    record = lanczos_extremal_spectrum(
        _dense_matvec(projected, n_nodes), n_nodes,
        upper_bound=float(np.max(np.sum(np.abs(projected), axis=1))),
        n_iterations=n_dof - n_constraints, projector=apply_p)

    assert record["lambda_max_pN_per_um"] == pytest.approx(float(restricted[-1]), rel=1e-9)
    assert record["lambda_min_pN_per_um"] == pytest.approx(float(restricted[0]), rel=1e-6)
    assert record["lambda_min_resolved"] is True
    assert record["condition_number"] == pytest.approx(
        float(restricted[-1] / restricted[0]), rel=1e-6)
    assert "range(P)" in record["subspace"]


def test_without_the_projector_the_same_operator_reports_its_constraint_zeros() -> None:
    """The control that makes the test above mean something: the defect it prevents is real and silent."""
    n_nodes, n_constraints = 24, 17
    _, _, projected = _projected_reference(n_nodes, n_constraints=n_constraints, seed=21)
    n_dof = 3 * n_nodes
    restricted = np.linalg.eigvalsh(projected)[n_constraints:]

    record = lanczos_extremal_spectrum(
        _dense_matvec(projected, n_nodes), n_nodes,
        upper_bound=float(np.max(np.sum(np.abs(projected), axis=1))),
        n_iterations=n_dof - n_constraints)

    # λ_max is unaffected — the null space is at the OTHER end — but the small end collapses to the
    # projector's zeros, so a condition number read from this record would be about the constraints.
    assert record["lambda_max_pN_per_um"] == pytest.approx(float(restricted[-1]), rel=1e-9)
    assert abs(record["lambda_min_pN_per_um"]) < 1e-8 * float(restricted[0])
    assert "the whole space" in record["subspace"]


def test_a_supplied_start_vector_is_used_and_a_degenerate_one_is_refused() -> None:
    n_nodes = 12
    n_dof = 3 * n_nodes
    matrix = _random_spd(n_dof, seed=31, floor=1.0)
    matvec = _dense_matvec(matrix, n_nodes)
    truth = np.linalg.eigvalsh(matrix)

    # An eigenvector start converges in ONE iteration, which a random start cannot do — so this shows
    # the argument reaches the recurrence rather than being accepted and ignored.
    top = np.linalg.eigh(matrix)[1][:, -1]
    one_step = lanczos_ritz(matvec, n_nodes, n_iterations=1, start_vector=top.reshape(n_nodes, 3))
    assert float(one_step.ritz_values[-1]) == pytest.approx(float(truth[-1]), rel=1e-10)

    with pytest.raises(ValueError, match="zero norm"):
        lanczos_ritz(matvec, n_nodes, n_iterations=4, start_vector=np.zeros((n_nodes, 3)))
    with pytest.raises(ValueError, match="the operator is"):
        lanczos_ritz(matvec, n_nodes, n_iterations=4, start_vector=np.ones(n_dof + 3))
