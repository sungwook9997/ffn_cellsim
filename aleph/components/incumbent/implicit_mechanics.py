"""All-device analytic projected-CG mechanics for the Active Cell Warp runtime.

This is the full-coupling alternative to the rejected local block and scalar-spectral experiments. Each
nonlinear mechanics iteration solves

``(a I + P K_analytic P) dx = P F``

with fixed-launch, device-predicated preconditioned conjugate gradient. ``K_analytic`` contains the exact
NF2007 bending Hessian plus current SPD tangents of crosslinks, WCA radial contact, explicit minifilament
bonds/crossbridges/angle harmonics, LINC, and unilateral ERM. A topology-aware additive coarse block restores
the per-fiber translation modes that node Jacobi suppresses by counting internal bending diagonals even though
those forces cancel under a collective fiber translation. The regularizer ``a`` is the maximum
already-declared stiffness
of force families not yet represented in ``K_analytic`` (surface area/bending, volume, and pressure), with a
machine-derived floor. It is therefore neither a biological parameter nor a value tuned to pass NG-3.

PCG scalars remain one-element CUDA arrays. The host launches a fixed maximum number of iterations; device
flags freeze recurrence updates after residual convergence or invalid curvature. Operator launches retain the
fixed schedule, and there is no per-CG or per-inner
device-to-host readback. A nonlinear candidate is never authorized from an unconverged PCG solve; the driver
then compares its exact projected residual with both the unchanged state and matched explicit-CFL candidate.

Sanity Gate:
    * Dimensions: ``A [pN/um] dx [um] = P F [pN]``.
    * Symmetry/positivity: ``a>0`` and every analytic tangent is PSD, so ``aI+PKP`` is SPD.
    * Constraints: operator input and stiffness output are projected, hence CG remains in the NF2007 tangent.
    * Boundary cases: zero RHS returns zero; compressed unilateral links contribute no force or stiffness.
    * Numerical: CG residual norm target is ``sqrt(eps64)`` relative, matching the engine's derived resolution;
      a fixed iteration cap is a reported compute budget, never interpreted as physical time.
    * Residency: vectors, dot products, recurrence scalars, convergence, and validity stay Warp-CUDA resident.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import warp as wp

from aleph.components.incumbent.implicit_mechanics_analytic import rigid_strain_coarse_basis

_LOG = logging.getLogger(__name__)

__all__ = [
    "ProjectedAnalyticCG",
    "add_angle_gauss_newton_stiffness_kernel",
    "add_bending_stiffness_kernel",
    "add_crossbridge_stiffness_kernel",
    "add_segment_crossbridge_stiffness_kernel",
    "add_pair_stiffness_kernel",
    "add_erm_stiffness_kernel",
    "add_erm_tension_side_preconditioner_kernel",
    "add_wca_stiffness_kernel",
    "add_fiber_translation_coarse_kernel",
    "build_fiber_block_cholesky_kernel",
    "apply_fiber_block_preconditioner_kernel",
    "build_erm_augmented_fiber_block_kernel",
    "apply_erm_augmented_fiber_block_kernel",
    "restrict_fiber_sum_kernel",
    "restrict_fiber_trace_majorizer_kernel",
    "prolong_fiber_translation_kernel",
    "assemble_coarse_matrix_kernel",
    "factor_coarse_cholesky_kernel",
    "coarse_restrict_kernel",
    "coarse_solve_kernel",
    "coarse_prolong_add_kernel",
    "compute_regularization_kernel",
    "conditional_displacement_kernel",
    "decide_better_line_search_trial_kernel",
    "conditional_copy_vec3_kernel",
    "commit_better_line_search_trial_kernel",
    "finalize_line_search_kernel",
    "copy_validity_kernel",
    "omitted_regularization_base",
]

_TWO_POW_1_6 = wp.constant(wp.float64(2.0 ** (1.0 / 6.0)))


@wp.func
def _central_action(delta: wp.vec3d, relative: wp.vec3d, stiffness: wp.float64,
                    rest_length: wp.float64) -> wp.vec3d:
    length = wp.length(delta)
    if length <= wp.float64(1.0e-12):
        return wp.vec3d(0.0, 0.0, 0.0)
    unit = delta / length
    axial = stiffness * wp.dot(unit, relative) * unit
    extension = length - rest_length
    transverse_k = wp.max(stiffness * extension / length, wp.float64(0.0))
    transverse = transverse_k * (relative - wp.dot(unit, relative) * unit)
    return axial + transverse


@wp.func
def _central_diagonal(delta: wp.vec3d, stiffness: wp.float64, rest_length: wp.float64) -> wp.vec3d:
    """Return the Cartesian diagonal of the same PSD central-spring tangent used by the operator."""
    length = wp.length(delta)
    if length <= wp.float64(1.0e-12):
        return wp.vec3d(0.0, 0.0, 0.0)
    unit = delta / length
    transverse_k = wp.max(stiffness * (length - rest_length) / length, wp.float64(0.0))
    difference = stiffness - transverse_k
    return wp.vec3d(
        transverse_k + difference * unit[0] * unit[0],
        transverse_k + difference * unit[1] * unit[1],
        transverse_k + difference * unit[2] * unit[2],
    )


@wp.kernel
def set_mass_action_kernel(
    vector: wp.array(dtype=wp.vec3d),
    regularization: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    out[i] = regularization[0] * vector[i]


@wp.kernel
def add_bending_stiffness_kernel(
    vector: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),
    alpha: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    a = triples[t, 0]
    b = triples[t, 1]
    c = triples[t, 2]
    curvature = vector[a] - wp.float64(2.0) * vector[b] + vector[c]
    action = alpha[t] * curvature
    wp.atomic_add(out, a, action)
    wp.atomic_add(out, b, -wp.float64(2.0) * action)
    wp.atomic_add(out, c, action)


@wp.kernel
def add_angle_gauss_newton_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.float64,
    theta0: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Add a PSD tangent of ``E = 0.5 k (theta-theta0)^2``.

    Away from a collinear triple this is the Gauss--Newton term ``k grad(theta) grad(theta)^T``. At the
    exactly straight backbone rest state, ``grad(theta)`` has no unique bending plane; the limiting quadratic
    energy has two transverse modes. That branch applies the isotropic straight-rod limit instead of silently
    returning zero stiffness. Both branches conserve translation by construction.
    """
    t = wp.tid()
    i = triples[t, 0]
    j = triples[t, 1]
    k = triples[t, 2]
    r1 = pos[i] - pos[j]
    r2 = pos[k] - pos[j]
    n1 = wp.length(r1)
    n2 = wp.length(r2)
    if n1 <= wp.float64(1.0e-12) or n2 <= wp.float64(1.0e-12):
        return
    u1 = r1 / n1
    u2 = r2 / n2
    c = wp.clamp(wp.dot(u1, u2), wp.float64(-1.0), wp.float64(1.0))
    sin_theta = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - c * c))
    straight_rest = wp.abs(theta0 - wp.float64(3.141592653589793)) <= wp.float64(1.0e-12)
    if straight_rest and c < wp.float64(0.0) and sin_theta < wp.float64(1.0e-8):
        # At theta=pi, d(u1+u2) is the two-dimensional transverse bend measure.
        axis = u2
        relative = (vector[i] - vector[j]) / n1 + (vector[k] - vector[j]) / n2
        bend = relative - wp.dot(axis, relative) * axis
        action_i = stiffness * bend / n1
        action_k = stiffness * bend / n2
    elif sin_theta >= wp.float64(1.0e-8):
        grad_i = -(u2 - c * u1) / (n1 * sin_theta)
        grad_k = -(u1 - c * u2) / (n2 * sin_theta)
        grad_j = -(grad_i + grad_k)
        directional_angle = (
            wp.dot(grad_i, vector[i])
            + wp.dot(grad_j, vector[j])
            + wp.dot(grad_k, vector[k])
        )
        action_i = stiffness * directional_angle * grad_i
        action_k = stiffness * directional_angle * grad_k
    else:
        return
    wp.atomic_add(out, i, action_i)
    wp.atomic_add(out, k, action_k)
    wp.atomic_add(out, j, -(action_i + action_k))


@wp.kernel
def add_pair_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.array(dtype=wp.float64),
    rest: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    action = _central_action(pos[j] - pos[i], vector[i] - vector[j], stiffness[t], rest[t])
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, j, -action)


@wp.kernel
def add_uniform_pair_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.float64,
    rest: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    action = _central_action(pos[j] - pos[i], vector[i] - vector[j], stiffness, rest)
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, j, -action)


@wp.kernel
def add_erm_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    stiffness: wp.float64,
    rupture_force: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    if bound[t] == 0:
        return
    i = membrane_idx[t]
    j = cortex_idx[t]
    delta = pos[j] - pos[i]
    length = wp.length(delta)
    extension = length - rest[t]
    if extension <= wp.float64(0.0) or stiffness * extension > rupture_force:
        return
    action = _central_action(delta, vector[i] - vector[j], stiffness, rest[t])
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, j, -action)


@wp.kernel
def add_linc_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    nucleus_idx: wp.array(dtype=wp.int32),
    anchor_idx: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    stiffness: wp.float64,
    stiffening: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    i = nucleus_idx[t]
    j = anchor_idx[t]
    delta = pos[j] - pos[i]
    length = wp.length(delta)
    extension = length - rest[t]
    if length <= wp.float64(1.0e-12) or extension <= wp.float64(0.0):
        return
    unit = delta / length
    relative = vector[i] - vector[j]
    tension = stiffness * extension * (wp.float64(1.0) + stiffening * extension * extension)
    axial_k = stiffness * (wp.float64(1.0) + wp.float64(3.0) * stiffening * extension * extension)
    action = axial_k * wp.dot(unit, relative) * unit
    transverse_k = tension / length
    action = action + transverse_k * (relative - wp.dot(unit, relative) * unit)
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, j, -action)


@wp.kernel
def add_crossbridge_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    stiffness: wp.float64,
    rest: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Add the fixed-KMC-state tangent of every bound, power-stroke-shifted crossbridge."""
    h = wp.tid()
    if bound[h] == 0 or anchor[h] < 0:
        return
    i = head_node[h]
    j = anchor[h]
    # Directional crossbridge tangent (rank-one K = k w wᵀ, matches segment_motor.py B; see the segment kernel).
    action = stiffness * wp.dot(walk_dir[h], vector[i] - vector[j]) * walk_dir[h]
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, j, -action)


@wp.kernel
def add_segment_crossbridge_stiffness_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    stiffness: wp.float64,
    rest: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Add the symmetric fixed-KMC tangent of a barycentric three-node crossbridge."""
    h = wp.tid()
    if bound[h] == 0 or seg_a[h] < 0 or seg_b[h] < 0:
        return
    i = head_node[h]
    a = seg_a[h]
    b = seg_b[h]
    t = bary_t[h]
    one_minus_t = wp.float64(1.0) - t
    relative = vector[i] - one_minus_t * vector[a] - t * vector[b]
    # Directional crossbridge tangent (2026-07-24, matches segment_motor.py crossbridge B): the force
    # ``f = k (d·w − r0) w`` with a frozen unit ``w = walk_dir`` is LINEAR in position, so its EXACT tangent is the
    # constant rank-one block ``K = k w wᵀ`` — force along ``w`` only. The old ``_central_action`` used ``d/|d|`` as
    # the axis and (for r0_xb=0) reduced to the isotropic ``k I``, injecting two spurious stiffness directions ⊥ w.
    action = stiffness * wp.dot(walk_dir[h], relative) * walk_dir[h]
    wp.atomic_add(out, i, action)
    wp.atomic_add(out, a, -one_minus_t * action)
    wp.atomic_add(out, b, -t * action)


@wp.kernel
def add_wca_stiffness_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    vector: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Add the symmetric PSD part of the current WCA tangent using the force kernel's live hash grid.

    In the repulsive branch the radial energy curvature is
    ``U'' = 24 eps/r^2 (26 (sigma/r)^12 - 7 (sigma/r)^6) > 0``. The transverse curvature is negative and
    is therefore omitted from this SPD solve. A force-capped pair has constant radial force and zero radial
    derivative, so it contributes zero rather than an invented core stiffness.
    """
    i = wp.tid()
    if active_nodes[i] < wp.int32(0):
        return
    zero = wp.float64(0.0)
    pi = pos[i]
    vi = vector[i]
    fi = fiber_id[i]
    cutoff = _TWO_POW_1_6 * sigma
    action = wp.vec3d(zero, zero, zero)
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fi:
            delta = pi - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < cutoff:
                ratio = sigma / length
                sr6 = ratio * ratio * ratio
                sr6 = sr6 * sr6
                force_magnitude = (wp.float64(24.0) * epsilon / length) * (
                    wp.float64(2.0) * sr6 * sr6 - sr6)
                capped = force_cap > zero and force_magnitude >= force_cap
                if not capped:
                    axial_k = (wp.float64(24.0) * epsilon / (length * length)) * (
                        wp.float64(26.0) * sr6 * sr6 - wp.float64(7.0) * sr6)
                    unit = delta / length
                    action = action + axial_k * wp.dot(unit, vi - vector[j]) * unit
    out[i] = out[i] + action


@wp.kernel
def compute_regularization_kernel(
    omitted_base: wp.float64,
    max_pressure_jump: wp.array(dtype=wp.float64),
    pressure_edge_um: wp.float64,
    regularization: wp.array(dtype=wp.float64),
) -> None:
    """Compose the live omitted-family stiffness scale ``a`` without host readback."""
    regularization[0] = wp.max(omitted_base, max_pressure_jump[0] * pressure_edge_um)


@wp.kernel
def conditional_displacement_kernel(
    pos: wp.array(dtype=wp.vec3d),
    displacement: wp.array(dtype=wp.vec3d),
    scale: wp.array(dtype=wp.float64),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
) -> None:
    """Apply one fixed line-search trial while both device validity latches remain set."""
    i = wp.tid()
    if active[0] != 0 and finite[0] != 0:
        pos[i] = pos[i] + scale[0] * displacement[i]


@wp.kernel
def decide_better_line_search_trial_kernel(
    trial_residual: wp.array(dtype=wp.float64),
    best_residual: wp.array(dtype=wp.float64),
    trial_finite: wp.array(dtype=wp.int32),
    best_finite: wp.array(dtype=wp.int32),
    linear_converged: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    take_trial: wp.array(dtype=wp.int32),
) -> None:
    """Decide one residual-monotone line-search comparison entirely on device."""
    take_trial[0] = wp.int32(0)
    if active[0] != 0 and linear_converged[0] != 0 and trial_finite[0] != 0:
        if best_finite[0] == 0 or trial_residual[0] < best_residual[0]:
            take_trial[0] = wp.int32(1)


@wp.kernel
def conditional_copy_vec3_kernel(
    source: wp.array(dtype=wp.vec3d),
    destination: wp.array(dtype=wp.vec3d),
    predicate: wp.array(dtype=wp.int32),
) -> None:
    """Copy a vector candidate under a device scalar predicate."""
    i = wp.tid()
    if predicate[0] != 0:
        destination[i] = source[i]


@wp.kernel
def commit_better_line_search_trial_kernel(
    trial_residual: wp.array(dtype=wp.float64),
    trial_finite: wp.array(dtype=wp.int32),
    trial_index: wp.int32,
    take_trial: wp.array(dtype=wp.int32),
    best_residual: wp.array(dtype=wp.float64),
    best_finite: wp.array(dtype=wp.int32),
    best_trial_index: wp.array(dtype=wp.int32),
) -> None:
    """Commit the scalar metadata after the predicated candidate copy has completed."""
    if take_trial[0] != 0:
        best_residual[0] = trial_residual[0]
        best_finite[0] = trial_finite[0]
        best_trial_index[0] = trial_index


@wp.kernel
def finalize_line_search_kernel(
    best_finite: wp.array(dtype=wp.int32),
    best_trial_index: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    implicit_accept_count: wp.array(dtype=wp.int32),
    explicit_accept_count: wp.array(dtype=wp.int32),
    stationary_accept_count: wp.array(dtype=wp.int32),
    scale_accept_counts: wp.array(dtype=wp.int32),
) -> None:
    """Install validity and classify the winning current, explicit, or accelerated candidate."""
    finite[0] = best_finite[0]
    index = best_trial_index[0]
    if index >= wp.int32(0):
        implicit_accept_count[0] += wp.int32(1)
        wp.atomic_add(scale_accept_counts, index, wp.int32(1))
    elif index == wp.int32(-1):
        explicit_accept_count[0] += wp.int32(1)
    else:
        stationary_accept_count[0] += wp.int32(1)


@wp.kernel
def copy_validity_kernel(source: wp.array(dtype=wp.int32), destination: wp.array(dtype=wp.int32)) -> None:
    """Copy one device validity latch without a host decision."""
    destination[0] = source[0]


@wp.kernel
def set_preconditioner_kernel(
    regularization: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    value = regularization[0]
    diagonal[i] = wp.vec3d(value, value, value)


@wp.kernel
def add_bending_preconditioner_kernel(
    triples: wp.array(dtype=wp.int32, ndim=2),
    alpha: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    value = alpha[t]
    isotropic = wp.vec3d(value, value, value)
    wp.atomic_add(diagonal, triples[t, 0], isotropic)
    wp.atomic_add(diagonal, triples[t, 1], wp.float64(4.0) * isotropic)
    wp.atomic_add(diagonal, triples[t, 2], isotropic)


@wp.kernel
def add_angle_gauss_newton_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.float64,
    theta0: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    """Add the Cartesian block diagonal of the angle Gauss--Newton/straight-limit tangent."""
    t = wp.tid()
    i = triples[t, 0]
    j = triples[t, 1]
    k = triples[t, 2]
    r1 = pos[i] - pos[j]
    r2 = pos[k] - pos[j]
    n1 = wp.length(r1)
    n2 = wp.length(r2)
    if n1 <= wp.float64(1.0e-12) or n2 <= wp.float64(1.0e-12):
        return
    u1 = r1 / n1
    u2 = r2 / n2
    c = wp.clamp(wp.dot(u1, u2), wp.float64(-1.0), wp.float64(1.0))
    sin_theta = wp.sqrt(wp.max(wp.float64(0.0), wp.float64(1.0) - c * c))
    straight_rest = wp.abs(theta0 - wp.float64(3.141592653589793)) <= wp.float64(1.0e-12)
    if straight_rest and c < wp.float64(0.0) and sin_theta < wp.float64(1.0e-8):
        axis = u2
        projection_diagonal = wp.vec3d(
            wp.float64(1.0) - axis[0] * axis[0],
            wp.float64(1.0) - axis[1] * axis[1],
            wp.float64(1.0) - axis[2] * axis[2],
        )
        inv_n1 = wp.float64(1.0) / n1
        inv_n2 = wp.float64(1.0) / n2
        value_i = stiffness * inv_n1 * inv_n1 * projection_diagonal
        value_k = stiffness * inv_n2 * inv_n2 * projection_diagonal
        value_j = stiffness * (inv_n1 + inv_n2) * (inv_n1 + inv_n2) * projection_diagonal
    elif sin_theta >= wp.float64(1.0e-8):
        grad_i = -(u2 - c * u1) / (n1 * sin_theta)
        grad_k = -(u1 - c * u2) / (n2 * sin_theta)
        grad_j = -(grad_i + grad_k)
        value_i = stiffness * wp.vec3d(
            grad_i[0] * grad_i[0], grad_i[1] * grad_i[1], grad_i[2] * grad_i[2])
        value_j = stiffness * wp.vec3d(
            grad_j[0] * grad_j[0], grad_j[1] * grad_j[1], grad_j[2] * grad_j[2])
        value_k = stiffness * wp.vec3d(
            grad_k[0] * grad_k[0], grad_k[1] * grad_k[1], grad_k[2] * grad_k[2])
    else:
        return
    wp.atomic_add(diagonal, i, value_i)
    wp.atomic_add(diagonal, j, value_j)
    wp.atomic_add(diagonal, k, value_k)


@wp.kernel
def add_pair_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.array(dtype=wp.float64),
    rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    value = _central_diagonal(pos[j] - pos[i], stiffness[t], rest[t])
    wp.atomic_add(diagonal, i, value)
    wp.atomic_add(diagonal, j, value)


@wp.kernel
def add_uniform_pair_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.float64,
    rest: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    i = links[t, 0]
    j = links[t, 1]
    value = _central_diagonal(pos[j] - pos[i], stiffness, rest)
    wp.atomic_add(diagonal, i, value)
    wp.atomic_add(diagonal, j, value)


@wp.kernel
def add_erm_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    stiffness: wp.float64,
    rupture_force: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    if bound[t] == 0:
        return
    extension = wp.length(pos[cortex_idx[t]] - pos[membrane_idx[t]]) - rest[t]
    if extension <= wp.float64(0.0) or stiffness * extension > rupture_force:
        return
    value = _central_diagonal(
        pos[cortex_idx[t]] - pos[membrane_idx[t]], stiffness, rest[t])
    wp.atomic_add(diagonal, membrane_idx[t], value)
    wp.atomic_add(diagonal, cortex_idx[t], value)


@wp.kernel
def add_erm_tension_side_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    stiffness: wp.float64,
    rupture_force: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Add the ERM TENSION-SIDE axial diagonal ``k (û⊗û)`` for every bound, sub-rupture tether.

    Unlike :func:`add_erm_preconditioner_kernel` (which skips a tether at ``extension<=0``, matching the
    unilateral FORCE), this preconditioner anticipates that a RESTING (force-free, ``extension==0``) ERM WILL
    engage under the outward turgor traction — the membrane must move outward to balance turgor, which stretches
    the tether. ``_central_diagonal`` already returns the pure axial ``k û_α²`` (its transverse part clamps to
    zero at/below rest), i.e. exactly the tension-side stiffness. Feeding it into the membrane and cortex node
    diagonals scales the membrane's balancing step to ``turgor/(a+γ_area+k_erm)`` instead of ``turgor/(a+γ_area)``
    (a ~k_erm/a ≈ 5× over-step that the line search then throttles to a crawl). Preconditioner-only: the ERM
    force stays exactly unilateral; this only conditions the resting membrane→ERM→cortex step and is used by the
    backbone-aware augmented block path.
    """
    t = wp.tid()
    if bound[t] == 0:
        return
    delta = pos[cortex_idx[t]] - pos[membrane_idx[t]]
    length = wp.length(delta)
    if length <= wp.float64(1.0e-12):
        return
    extension = length - rest[t]
    if extension > wp.float64(0.0) and stiffness * extension > rupture_force:
        return  # a genuinely over-stretched (ruptured) tether carries no stiffness
    value = _central_diagonal(delta, stiffness, rest[t])
    wp.atomic_add(diagonal, membrane_idx[t], value)
    wp.atomic_add(diagonal, cortex_idx[t], value)


@wp.kernel
def add_linc_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    nucleus_idx: wp.array(dtype=wp.int32),
    anchor_idx: wp.array(dtype=wp.int32),
    rest: wp.array(dtype=wp.float64),
    stiffness: wp.float64,
    stiffening: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    t = wp.tid()
    extension = wp.length(pos[anchor_idx[t]] - pos[nucleus_idx[t]]) - rest[t]
    if extension <= wp.float64(0.0):
        return
    delta = pos[anchor_idx[t]] - pos[nucleus_idx[t]]
    length = wp.length(delta)
    unit = delta / length
    tension = stiffness * extension * (wp.float64(1.0) + stiffening * extension * extension)
    axial = stiffness * (wp.float64(1.0) + wp.float64(3.0) * stiffening * extension * extension)
    transverse = tension / length
    difference = axial - transverse
    value = wp.vec3d(
        transverse + difference * unit[0] * unit[0],
        transverse + difference * unit[1] * unit[1],
        transverse + difference * unit[2] * unit[2],
    )
    wp.atomic_add(diagonal, nucleus_idx[t], value)
    wp.atomic_add(diagonal, anchor_idx[t], value)


@wp.kernel
def add_crossbridge_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    anchor: wp.array(dtype=wp.int32),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    stiffness: wp.float64,
    rest: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    h = wp.tid()
    if bound[h] == 0 or anchor[h] < 0:
        return
    i = head_node[h]
    j = anchor[h]
    # Directional crossbridge diagonal: diag(k w wᵀ) = k (w_x², w_y², w_z²) (matches the segment kernel).
    w = walk_dir[h]
    value = stiffness * wp.vec3d(w[0] * w[0], w[1] * w[1], w[2] * w[2])
    wp.atomic_add(diagonal, i, value)
    wp.atomic_add(diagonal, j, value)


@wp.kernel
def add_segment_crossbridge_preconditioner_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    stiffness: wp.float64,
    rest: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    """Block-Jacobi diagonal of the barycentric three-node crossbridge tangent."""
    h = wp.tid()
    if bound[h] == 0 or seg_a[h] < 0 or seg_b[h] < 0:
        return
    i = head_node[h]
    a = seg_a[h]
    b = seg_b[h]
    t = bary_t[h]
    one_minus_t = wp.float64(1.0) - t
    # Directional crossbridge diagonal (matches the rank-one K = k w wᵀ tangent): the Cartesian diagonal of
    # ``k w wᵀ`` is ``k (w_x², w_y², w_z²)`` — no delta-based axis, no clamped-transverse isotropy.
    w = walk_dir[h]
    value = stiffness * wp.vec3d(w[0] * w[0], w[1] * w[1], w[2] * w[2])
    wp.atomic_add(diagonal, i, value)
    wp.atomic_add(diagonal, a, one_minus_t * one_minus_t * value)
    wp.atomic_add(diagonal, b, t * t * value)


@wp.kernel
def add_wca_preconditioner_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    diagonal: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    if active_nodes[i] < wp.int32(0):
        return
    zero = wp.float64(0.0)
    pi = pos[i]
    fi = fiber_id[i]
    cutoff = _TWO_POW_1_6 * sigma
    value = wp.vec3d(zero, zero, zero)
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fi:
            delta = pi - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < cutoff:
                ratio = sigma / length
                sr6 = ratio * ratio * ratio
                sr6 = sr6 * sr6
                force_magnitude = (wp.float64(24.0) * epsilon / length) * (
                    wp.float64(2.0) * sr6 * sr6 - sr6)
                if not (force_cap > zero and force_magnitude >= force_cap):
                    axial = (wp.float64(24.0) * epsilon / (length * length)) * (
                        wp.float64(26.0) * sr6 * sr6 - wp.float64(7.0) * sr6)
                    unit = delta / length
                    value = value + wp.vec3d(
                        axial * unit[0] * unit[0],
                        axial * unit[1] * unit[1],
                        axial * unit[2] * unit[2],
                    )
    diagonal[i] = diagonal[i] + value


@wp.kernel
def apply_preconditioner_kernel(
    residual: wp.array(dtype=wp.vec3d),
    diagonal: wp.array(dtype=wp.vec3d),
    preconditioned: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    r = residual[i]
    d = diagonal[i]
    preconditioned[i] = wp.vec3d(r[0] / d[0], r[1] / d[1], r[2] / d[2])


@wp.kernel
def add_fiber_translation_coarse_kernel(
    residual: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    external_diagonal: wp.array(dtype=wp.vec3d),
    preconditioned: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Add an SPD block-Jacobi solve on each actin fiber's three translation modes.

    ``external_diagonal`` contains only the regularizer and couplings that do not cancel when every node of
    one fiber translates together (crosslinks, WCA, crossbridges, LINC, and ERM). For a translation basis
    vector ``z_f``, the diagonal coarse Rayleigh approximation is therefore the componentwise sum over the
    fiber. Adding ``z_f (z_f^T D_ext z_f)^-1 z_f^T`` to node Jacobi is symmetric positive semidefinite and
    introduces no relaxation parameter.
    """
    f = wp.tid()
    begin = fiber_off[f]
    end = fiber_off[f + 1]
    rhs = wp.vec3d(0.0, 0.0, 0.0)
    diagonal = wp.vec3d(0.0, 0.0, 0.0)
    for node in range(begin, end):
        rhs = rhs + residual[node]
        diagonal = diagonal + external_diagonal[node]
    correction = wp.vec3d(0.0, 0.0, 0.0)
    if diagonal[0] > wp.float64(1.0e-300):
        correction[0] = rhs[0] / diagonal[0]
    if diagonal[1] > wp.float64(1.0e-300):
        correction[1] = rhs[1] / diagonal[1]
    if diagonal[2] > wp.float64(1.0e-300):
        correction[2] = rhs[2] / diagonal[2]
    for node in range(begin, end):
        preconditioned[node] = preconditioned[node] + correction


@wp.kernel
def build_fiber_block_cholesky_kernel(
    external_diagonal: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    triple_off: wp.array(dtype=wp.int32),
    bending_alpha: wp.array(dtype=wp.float64),
    max_nodes: wp.int32,
    factor: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
) -> None:
    r"""Factor each actin fiber's ``D_external + D2^T alpha D2`` Cartesian block.

    One CUDA thread owns one short fiber. The three Cartesian blocks share the exact NF2007 bending stencil
    but retain their componentwise external diagonal. Dense Cholesky is used because the native discretization
    has only seven nodes per cortical filament; ``max_nodes`` is the allocated topology maximum, not a physics
    parameter.
    """
    f = wp.tid()
    begin = fiber_off[f]
    node_count = fiber_off[f + 1] - begin
    triple_begin = triple_off[f]
    block_stride = max_nodes * max_nodes
    for component in range(3):
        base = (f * wp.int32(3) + component) * block_stride
        for row in range(node_count):
            for column in range(node_count):
                factor[base + row * max_nodes + column] = wp.float64(0.0)
        for row in range(node_count):
            factor[base + row * max_nodes + row] = external_diagonal[begin + row][component]
        for local_triple in range(node_count - 2):
            alpha = bending_alpha[triple_begin + local_triple]
            for local_row in range(3):
                row_weight = wp.float64(1.0)
                if local_row == 1:
                    row_weight = wp.float64(-2.0)
                row = local_triple + local_row
                for local_column in range(3):
                    column_weight = wp.float64(1.0)
                    if local_column == 1:
                        column_weight = wp.float64(-2.0)
                    column = local_triple + local_column
                    index = base + row * max_nodes + column
                    factor[index] = factor[index] + alpha * row_weight * column_weight
        for row in range(node_count):
            for column in range(row + 1):
                value = factor[base + row * max_nodes + column]
                for inner in range(column):
                    value = value - (
                        factor[base + row * max_nodes + inner]
                        * factor[base + column * max_nodes + inner]
                    )
                if row == column:
                    if value <= wp.float64(1.0e-300) or not wp.isfinite(value):
                        finite[0] = wp.int32(0)
                        return
                    factor[base + row * max_nodes + column] = wp.sqrt(value)
                else:
                    factor[base + row * max_nodes + column] = (
                        value / factor[base + column * max_nodes + column])


@wp.kernel
def apply_fiber_block_preconditioner_kernel(
    residual: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    max_nodes: wp.int32,
    factor: wp.array(dtype=wp.float64),
    solve_work: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
    preconditioned: wp.array(dtype=wp.vec3d),
) -> None:
    """Overwrite actin Jacobi values with exact independent-fiber Cholesky block solves."""
    f = wp.tid()
    if finite[0] == wp.int32(0):
        return
    begin = fiber_off[f]
    node_count = fiber_off[f + 1] - begin
    block_stride = max_nodes * max_nodes
    for component in range(3):
        base = (f * wp.int32(3) + component) * block_stride
        for row in range(node_count):
            value = residual[begin + row][component]
            for column in range(row):
                value = value - (
                    factor[base + row * max_nodes + column]
                    * solve_work[wp.int32(3) * (begin + column) + component]
                )
            solve_work[wp.int32(3) * (begin + row) + component] = (
                value / factor[base + row * max_nodes + row])
        for reverse in range(node_count):
            row = node_count - wp.int32(1) - reverse
            value = solve_work[wp.int32(3) * (begin + row) + component]
            for column in range(row + 1, node_count):
                value = value - (
                    factor[base + column * max_nodes + row]
                    * solve_work[wp.int32(3) * (begin + column) + component]
                )
            solve_work[wp.int32(3) * (begin + row) + component] = (
                value / factor[base + row * max_nodes + row])
    for row in range(node_count):
        node = begin + row
        preconditioned[node] = wp.vec3d(
            solve_work[wp.int32(3) * node],
            solve_work[wp.int32(3) * node + wp.int32(1)],
            solve_work[wp.int32(3) * node + wp.int32(2)],
        )


@wp.func
def _aug_local_global(
    local: wp.int32,
    begin: wp.int32,
    node_count: wp.int32,
    erm_begin: wp.int32,
    erm_mem_global: wp.array(dtype=wp.int32),
) -> wp.int32:
    """Map an augmented-block local row (actin first, then ERM-tethered membrane) to its global node."""
    if local < node_count:
        return begin + local
    return erm_mem_global[erm_begin + (local - node_count)]


@wp.kernel
def build_erm_augmented_fiber_block_kernel(
    pos: wp.array(dtype=wp.vec3d),
    external_diagonal: wp.array(dtype=wp.vec3d),
    membrane_diagonal: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    triple_off: wp.array(dtype=wp.int32),
    bending_alpha: wp.array(dtype=wp.float64),
    erm_off: wp.array(dtype=wp.int32),
    erm_local_c: wp.array(dtype=wp.int32),
    erm_mem_global: wp.array(dtype=wp.int32),
    erm_index: wp.array(dtype=wp.int32),
    erm_bound: wp.array(dtype=wp.int32),
    erm_rest: wp.array(dtype=wp.float64),
    k_erm: wp.float64,
    rupture_force: wp.float64,
    max_aug: wp.int32,
    factor: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
) -> None:
    r"""Factor each cortical fiber's BACKBONE-AWARE block AUGMENTED with its ERM-tethered membrane DOFs.

    The plain :func:`build_fiber_block_cholesky_kernel` captures each fiber's backbone bending +
    external-diagonal but leaves the soft membrane↔ERM↔stiff-cortex coupling only diagonal, which is the
    conditioning gap diagnosed 2026-07-22e (a spring-pair Schwarz over-steps because it cannot see the
    backbone constraint). This block appends every membrane node tethered to this fiber (one ERM per membrane
    node ⇒ each membrane node is in exactly one fiber block, disjoint across fibers) and adds the ERM tangent
    as an explicit off-diagonal, so the block solves membrane→ERM→cortex→backbone consistently. The three
    Cartesian components stay decoupled (bending is isotropic; the external diagonals are Cartesian; the ERM
    coupling is taken per-component from the same PSD ``_central_diagonal`` the operator's diagonal uses, so
    the block remains a principal submatrix of the SPD ``aI+PKP`` preconditioner — diagonally dominant in the
    2-node ERM sub-blocks because the +H already lives in both node diagonals). ``max_aug`` is the allocated
    topology maximum (max fiber nodes + max ERM-per-fiber), not a physics parameter.
    """
    f = wp.tid()
    begin = fiber_off[f]
    node_count = fiber_off[f + 1] - begin
    triple_begin = triple_off[f]
    erm_begin = erm_off[f]
    n_erm = erm_off[f + 1] - erm_begin
    n_aug = node_count + n_erm
    block_stride = max_aug * max_aug
    for component in range(3):
        base = (f * wp.int32(3) + component) * block_stride
        for row in range(n_aug):
            for column in range(n_aug):
                factor[base + row * max_aug + column] = wp.float64(0.0)
        for row in range(node_count):
            factor[base + row * max_aug + row] = external_diagonal[begin + row][component]
        for e in range(n_erm):
            mrow = node_count + e
            m = erm_mem_global[erm_begin + e]
            factor[base + mrow * max_aug + mrow] = membrane_diagonal[m][component]
        for local_triple in range(node_count - 2):
            alpha = bending_alpha[triple_begin + local_triple]
            for local_row in range(3):
                row_weight = wp.float64(1.0)
                if local_row == 1:
                    row_weight = wp.float64(-2.0)
                row = local_triple + local_row
                for local_column in range(3):
                    column_weight = wp.float64(1.0)
                    if local_column == 1:
                        column_weight = wp.float64(-2.0)
                    column = local_triple + local_column
                    index = base + row * max_aug + column
                    factor[index] = factor[index] + alpha * row_weight * column_weight
        for e in range(n_erm):
            idx_e = erm_index[erm_begin + e]
            if erm_bound[idx_e] != wp.int32(0):
                lc = erm_local_c[erm_begin + e]
                mrow = node_count + e
                m = erm_mem_global[erm_begin + e]
                delta = pos[begin + lc] - pos[m]
                length = wp.length(delta)
                extension = length - erm_rest[idx_e]
                # TENSION-SIDE coupling: engage the ERM axial stiffness even at the force-free rest
                # (extension==0), because the resting membrane WILL stretch it under turgor. A ruptured
                # (over-stretched) tether is the only one dropped. Consistent with
                # add_erm_tension_side_preconditioner_kernel feeding the node diagonals, so each 2-node ERM
                # sub-block stays diagonally dominant / SPD.
                if (length > wp.float64(1.0e-12)
                        and not (extension > wp.float64(0.0) and k_erm * extension > rupture_force)):
                    coupling = _central_diagonal(delta, k_erm, erm_rest[idx_e])[component]
                    factor[base + lc * max_aug + mrow] = factor[base + lc * max_aug + mrow] - coupling
                    factor[base + mrow * max_aug + lc] = factor[base + mrow * max_aug + lc] - coupling
        for row in range(n_aug):
            for column in range(row + 1):
                value = factor[base + row * max_aug + column]
                for inner in range(column):
                    value = value - (
                        factor[base + row * max_aug + inner]
                        * factor[base + column * max_aug + inner]
                    )
                if row == column:
                    if value <= wp.float64(1.0e-300) or not wp.isfinite(value):
                        finite[0] = wp.int32(0)
                        return
                    factor[base + row * max_aug + column] = wp.sqrt(value)
                else:
                    factor[base + row * max_aug + column] = (
                        value / factor[base + column * max_aug + column])


@wp.kernel
def apply_erm_augmented_fiber_block_kernel(
    residual: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    erm_off: wp.array(dtype=wp.int32),
    erm_mem_global: wp.array(dtype=wp.int32),
    max_aug: wp.int32,
    factor: wp.array(dtype=wp.float64),
    solve_work: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
    preconditioned: wp.array(dtype=wp.vec3d),
) -> None:
    """Exact Cholesky solve of each backbone-aware ERM-augmented fiber block (actin + membrane rows)."""
    f = wp.tid()
    if finite[0] == wp.int32(0):
        return
    begin = fiber_off[f]
    node_count = fiber_off[f + 1] - begin
    erm_begin = erm_off[f]
    n_erm = erm_off[f + 1] - erm_begin
    n_aug = node_count + n_erm
    block_stride = max_aug * max_aug
    for component in range(3):
        base = (f * wp.int32(3) + component) * block_stride
        for row in range(n_aug):
            g_row = _aug_local_global(row, begin, node_count, erm_begin, erm_mem_global)
            value = residual[g_row][component]
            for column in range(row):
                g_col = _aug_local_global(column, begin, node_count, erm_begin, erm_mem_global)
                value = value - (
                    factor[base + row * max_aug + column]
                    * solve_work[wp.int32(3) * g_col + component]
                )
            solve_work[wp.int32(3) * g_row + component] = value / factor[base + row * max_aug + row]
        for reverse in range(n_aug):
            row = n_aug - wp.int32(1) - reverse
            g_row = _aug_local_global(row, begin, node_count, erm_begin, erm_mem_global)
            value = solve_work[wp.int32(3) * g_row + component]
            for column in range(row + 1, n_aug):
                g_col = _aug_local_global(column, begin, node_count, erm_begin, erm_mem_global)
                value = value - (
                    factor[base + column * max_aug + row]
                    * solve_work[wp.int32(3) * g_col + component]
                )
            solve_work[wp.int32(3) * g_row + component] = value / factor[base + row * max_aug + row]
    for row in range(n_aug):
        g_row = _aug_local_global(row, begin, node_count, erm_begin, erm_mem_global)
        preconditioned[g_row] = wp.vec3d(
            solve_work[wp.int32(3) * g_row],
            solve_work[wp.int32(3) * g_row + wp.int32(1)],
            solve_work[wp.int32(3) * g_row + wp.int32(2)],
        )


@wp.kernel
def restrict_fiber_sum_kernel(
    fine: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    coarse: wp.array(dtype=wp.vec3d),
) -> None:
    """Restrict a fine vector with the transpose of piecewise-constant fiber prolongation."""
    f = wp.tid()
    value = wp.vec3d(0.0, 0.0, 0.0)
    for node in range(fiber_off[f], fiber_off[f + 1]):
        value = value + fine[node]
    coarse[f] = value


@wp.kernel
def restrict_fiber_trace_majorizer_kernel(
    fine_diagonal: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    coarse_majorizer: wp.array(dtype=wp.vec3d),
) -> None:
    """Restrict a conservative isotropic majorizer for the vector-valued contact graph.

    Every represented pair tangent is PSD, so its trace bounds its largest Cartesian eigenvalue. Summing
    per-node traces therefore handles cross-component entries that a componentwise diagonal would miss and
    gives the graph-normalized spectral bound used by the fixed 1/2 Jacobi damping.
    """
    f = wp.tid()
    trace_sum = wp.float64(0.0)
    for node in range(fiber_off[f], fiber_off[f + 1]):
        value = fine_diagonal[node]
        trace_sum = trace_sum + value[0] + value[1] + value[2]
    coarse_majorizer[f] = wp.vec3d(trace_sum, trace_sum, trace_sum)


@wp.kernel
def prolong_fiber_translation_kernel(
    coarse: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    fine: wp.array(dtype=wp.vec3d),
) -> None:
    """Prolong one translation vector to every node of its disjoint actin fiber."""
    f = wp.tid()
    for node in range(fiber_off[f], fiber_off[f + 1]):
        fine[node] = coarse[f]


@wp.kernel
def coarse_jacobi_step_kernel(
    residual: wp.array(dtype=wp.vec3d),
    diagonal: wp.array(dtype=wp.vec3d),
    step: wp.array(dtype=wp.vec3d),
) -> None:
    """Apply the derived 1/2-damped Jacobi step for an SPD graph-Laplacian coarse operator."""
    f = wp.tid()
    r = residual[f]
    d = diagonal[f]
    half = wp.float64(0.5)
    step[f] = wp.vec3d(half * r[0] / d[0], half * r[1] / d[1], half * r[2] / d[2])


@wp.kernel
def add_coarse_step_kernel(
    step: wp.array(dtype=wp.vec3d),
    correction: wp.array(dtype=wp.vec3d),
) -> None:
    f = wp.tid()
    correction[f] = correction[f] + step[f]


@wp.kernel
def subtract_coarse_action_kernel(
    action: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
) -> None:
    f = wp.tid()
    residual[f] = residual[f] - action[f]


@wp.kernel
def set_initial_residual_kernel(
    rhs: wp.array(dtype=wp.vec3d),
    action: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    residual[i] = rhs[i] - action[i]


@wp.kernel
def assemble_coarse_matrix_kernel(
    basis: wp.array2d(dtype=wp.vec3d),
    action: wp.array2d(dtype=wp.vec3d),
    n_nodes: wp.int32,
    n_modes: wp.int32,
    matrix: wp.array(dtype=wp.float64),
) -> None:
    r"""Form the Galerkin coarse operator ``A_c[j, k] = b_j^T (A b_k)`` (row-major ``n_modes x n_modes``).

    ``action[k]`` is the already-computed operator apply ``A b_k``; one thread owns one ``(j, k)`` entry and
    reduces the full ``3N`` inner product. The result is symmetric to round-off because ``A`` is SPD.
    """
    entry = wp.tid()
    j = entry / n_modes
    k = entry - j * n_modes
    acc = wp.float64(0.0)
    for i in range(n_nodes):
        acc = acc + wp.dot(basis[j, i], action[k, i])
    matrix[entry] = acc


@wp.kernel
def factor_coarse_cholesky_kernel(
    matrix: wp.array(dtype=wp.float64),
    n_modes: wp.int32,
    factor: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
) -> None:
    r"""Dense lower-Cholesky factor of the tiny ``n_modes x n_modes`` SPD coarse operator (single thread).

    ``a I`` in the fine operator guarantees ``A_c = B^T A B`` is SPD for any full-rank ``B``; a non-positive
    pivot can therefore only be a non-finite input and latches the existing device validity predicate rather
    than emitting a spurious coarse correction.
    """
    for idx in range(n_modes * n_modes):
        factor[idx] = matrix[idx]
    for row in range(n_modes):
        for column in range(row + 1):
            value = factor[row * n_modes + column]
            for inner in range(column):
                value = value - factor[row * n_modes + inner] * factor[column * n_modes + inner]
            if row == column:
                if value <= wp.float64(1.0e-300) or not wp.isfinite(value):
                    finite[0] = wp.int32(0)
                    return
                factor[row * n_modes + column] = wp.sqrt(value)
            else:
                factor[row * n_modes + column] = value / factor[column * n_modes + column]


@wp.kernel
def coarse_restrict_kernel(
    basis: wp.array2d(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
    n_nodes: wp.int32,
    coarse_rhs: wp.array(dtype=wp.float64),
) -> None:
    """Restrict a fine residual onto the coarse space: ``coarse_rhs[k] = b_k^T r`` (one thread per mode)."""
    k = wp.tid()
    acc = wp.float64(0.0)
    for i in range(n_nodes):
        acc = acc + wp.dot(basis[k, i], residual[i])
    coarse_rhs[k] = acc


@wp.kernel
def coarse_solve_kernel(
    factor: wp.array(dtype=wp.float64),
    coarse_rhs: wp.array(dtype=wp.float64),
    n_modes: wp.int32,
    finite: wp.array(dtype=wp.int32),
    coarse_solution: wp.array(dtype=wp.float64),
) -> None:
    """Solve ``L L^T y = coarse_rhs`` by dense forward/back substitution (single thread)."""
    if finite[0] == wp.int32(0):
        for i in range(n_modes):
            coarse_solution[i] = wp.float64(0.0)
        return
    for row in range(n_modes):
        value = coarse_rhs[row]
        for column in range(row):
            value = value - factor[row * n_modes + column] * coarse_solution[column]
        coarse_solution[row] = value / factor[row * n_modes + row]
    for reverse in range(n_modes):
        row = n_modes - 1 - reverse
        value = coarse_solution[row]
        for column in range(row + 1, n_modes):
            value = value - factor[column * n_modes + row] * coarse_solution[column]
        coarse_solution[row] = value / factor[row * n_modes + row]


@wp.kernel
def coarse_prolong_add_kernel(
    basis: wp.array2d(dtype=wp.vec3d),
    coarse_solution: wp.array(dtype=wp.float64),
    n_modes: wp.int32,
    preconditioned: wp.array(dtype=wp.vec3d),
) -> None:
    """Prolong the coarse solution back and add it to the fine preconditioned residual."""
    i = wp.tid()
    correction = wp.vec3d(0.0, 0.0, 0.0)
    for k in range(n_modes):
        correction = correction + coarse_solution[k] * basis[k, i]
    preconditioned[i] = preconditioned[i] + correction


@wp.kernel
def _copy_zero_dx_kernel(
    rhs: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
    displacement: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    residual[i] = rhs[i]
    displacement[i] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _dot_kernel(a: wp.array(dtype=wp.vec3d), b: wp.array(dtype=wp.vec3d),
                out: wp.array(dtype=wp.float64)) -> None:
    i = wp.tid()
    wp.atomic_add(out, 0, wp.dot(a[i], b[i]))


@wp.kernel
def _cg_begin_kernel(rr: wp.array(dtype=wp.float64), rz: wp.array(dtype=wp.float64),
                     rr_initial: wp.array(dtype=wp.float64), rz_old: wp.array(dtype=wp.float64),
                     active: wp.array(dtype=wp.int32), converged: wp.array(dtype=wp.int32),
                     finite: wp.array(dtype=wp.int32), iterations: wp.array(dtype=wp.int32)) -> None:
    value = rr[0]
    rr_initial[0] = value
    rz_value = rz[0]
    rz_old[0] = rz_value
    iterations[0] = 0
    if value <= wp.float64(0.0):
        active[0] = 0
        converged[0] = 1
    elif rz_value <= wp.float64(1.0e-300) or not wp.isfinite(rz_value):
        active[0] = 0
        converged[0] = 0
        finite[0] = 0
    else:
        active[0] = 1
        converged[0] = 0


@wp.kernel
def _cg_alpha_kernel(rz_old: wp.array(dtype=wp.float64), p_ap: wp.array(dtype=wp.float64),
                     active: wp.array(dtype=wp.int32), finite: wp.array(dtype=wp.int32),
                     alpha: wp.array(dtype=wp.float64)) -> None:
    if active[0] == 0 or finite[0] == 0:
        active[0] = 0
        alpha[0] = wp.float64(0.0)
        return
    denominator = p_ap[0]
    if denominator <= wp.float64(1.0e-300) or not wp.isfinite(denominator):
        active[0] = 0
        finite[0] = 0
        alpha[0] = wp.float64(0.0)
    else:
        alpha[0] = rz_old[0] / denominator


@wp.kernel
def _cg_update_kernel(displacement: wp.array(dtype=wp.vec3d), residual: wp.array(dtype=wp.vec3d),
                      direction: wp.array(dtype=wp.vec3d), action: wp.array(dtype=wp.vec3d),
                      alpha: wp.array(dtype=wp.float64), active: wp.array(dtype=wp.int32)) -> None:
    i = wp.tid()
    if active[0] != 0:
        displacement[i] = displacement[i] + alpha[0] * direction[i]
        residual[i] = residual[i] - alpha[0] * action[i]


@wp.kernel
def _cg_finalize_kernel(rr_new: wp.array(dtype=wp.float64), rr_initial: wp.array(dtype=wp.float64),
                        eps64: wp.float64,
                        iteration: wp.int32, active: wp.array(dtype=wp.int32),
                        converged: wp.array(dtype=wp.int32), finite: wp.array(dtype=wp.int32),
                        iterations: wp.array(dtype=wp.int32)) -> None:
    if active[0] == 0:
        return
    value = rr_new[0]
    iterations[0] = iteration
    if not wp.isfinite(value):
        active[0] = 0
        converged[0] = 0
        finite[0] = 0
    elif value <= eps64 * rr_initial[0]:
        active[0] = 0
        converged[0] = 1


@wp.kernel
def _cg_beta_kernel(rz_new: wp.array(dtype=wp.float64), rz_old: wp.array(dtype=wp.float64),
                    active: wp.array(dtype=wp.int32), finite: wp.array(dtype=wp.int32),
                    beta: wp.array(dtype=wp.float64)) -> None:
    if active[0] == 0 or finite[0] == 0:
        beta[0] = wp.float64(0.0)
        return
    old = rz_old[0]
    new = rz_new[0]
    if old <= wp.float64(1.0e-300) or new <= wp.float64(0.0) or not wp.isfinite(new):
        active[0] = 0
        finite[0] = 0
        beta[0] = wp.float64(0.0)
    else:
        beta[0] = new / old
        rz_old[0] = new


@wp.kernel
def _cg_direction_kernel(preconditioned: wp.array(dtype=wp.vec3d), beta: wp.array(dtype=wp.float64),
                         active: wp.array(dtype=wp.int32), direction: wp.array(dtype=wp.vec3d)) -> None:
    i = wp.tid()
    if active[0] != 0:
        direction[i] = preconditioned[i] + beta[0] * direction[i]


def omitted_regularization_base(cell: object) -> float:
    """Return the largest declared stiffness not represented by the analytic matrix-free operator.

    This host calculation runs once when the solver workspace is configured. WCA is excluded because its
    current radial tangent is represented explicitly; the live pressure term is composed on device by
    :func:`compute_regularization_kernel`.
    """
    values = [float(np.sqrt(np.finfo(np.float64).eps) * cell.mechanical_kmax)]
    nucleus = cell.nucleus
    if nucleus is not None:
        edge = max(float(nucleus.mean_edge_um), np.sqrt(np.finfo(np.float64).eps))
        values.extend([
            float(nucleus.kappa_tilde) / edge**3,
            float(nucleus.k_soft),
            float(nucleus.k_ac),
            float(nucleus.k_vol),
        ])
    membrane = cell.membrane
    if membrane is not None:
        edge = max(float(membrane.mean_edge_um), np.sqrt(np.finfo(np.float64).eps))
        values.append(float(membrane.kappa_tilde) / edge**3)
        # γ_mem (area tension) is NO LONGER lumped here — it is now an explicit edge-spring term in the analytic
        # operator (ProjectedAnalyticCG.mem_edges_d), so counting it again would double-represent it.
    base = max(values)
    if not np.isfinite(base) or base <= 0.0:
        raise ValueError("omitted implicit regularization must be finite and positive")
    return base


class ProjectedAnalyticCG:
    """Fixed-launch, all-device CG workspace and analytic stiffness operator."""

    def __init__(
        self, cell: object, *, max_iterations: int = 32, coarse_iterations: int = 8,
        coarse_modes: int = 0, augment_fiber_block_erm: bool = False,
        erm_tension_side_precond: bool = False, disable_fiber_block: bool = False,
        fiber_quotient_coarse: bool = False, fq_coarse_iterations: int = 40,
        fq_coarse_mode: str = "pathA", multigrid: bool = False, mg_smooth: int = 2,
        mg_omega: float = 0.8, mg_fiber_quotient: bool = True, mg_min_fiber_nodes: int = 4,
        mg_max_levels: int = 6, mg_assembled_coarse: bool = False, cg_check_every: int = 0,
    ) -> None:
        if max_iterations <= 0:
            raise ValueError("max_iterations must be positive")
        if cg_check_every < 0:
            raise ValueError("cg_check_every must be nonnegative")
        if coarse_iterations < 0:
            raise ValueError("coarse_iterations must be nonnegative")
        if not 0 <= coarse_modes <= 12:
            raise ValueError("coarse_modes must be between 0 and 12")
        self.cell = cell
        self.device = cell.device
        self.n = int(cell.n_total)
        self.max_iterations = int(max_iterations)
        # SPEED: host early-exit cadence for the outer PCG. 0 == incumbent fixed-budget loop (byte-identical: the
        # Python loop always runs max_iterations, the device active/converged latch already no-ops the updates
        # after convergence). When > 0, read back the device `active` flag every `cg_check_every` iterations and
        # break once it clears (converged OR non-finite latch). The returned dx is BIT-IDENTICAL to running to
        # max_iterations (post-convergence launches only recompute no-op updates); this saves those launches —
        # each a full MG V-cycle at native scale. A control-flow readback for loop termination, NOT an
        # authoritative per-step physics roundtrip (the driver still globalizes with the exact residual).
        self.cg_check_every = int(cg_check_every)
        self.coarse_iterations = int(coarse_iterations)
        self.max_fiber_nodes = int(getattr(cell, "max_fiber_nodes", 0))
        # The per-fiber Cholesky block inverts a fiber's soft (bending/translation) modes without the
        # inextensibility projector or its inter-fiber crosslinks, so its direction OVER-STEPS the actin and a
        # single scalar line-search t (shared with the membrane) throttles the whole step to a crawl
        # (2026-07-23 diagnosis). ``disable_fiber_block`` drops it so the direction is pure PER-NODE scalar
        # Jacobi ``D^-1 r`` — each node gets its own stiffness scale (membrane turgor/(a+γ+k_erm), balanced
        # actin ~0), no over-step, so the coupled membrane→ERM→cortex Jacobi iteration converges.
        self.use_fiber_block = bool(
            cell.n_fibers and self.max_fiber_nodes > 0 and hasattr(cell, "toff_d")
            and not disable_fiber_block)
        # Backbone-aware ERM-augmented fiber block (fork-2 gate-closure design, 2026-07-23): extend each
        # per-fiber Cholesky block with its ERM-tethered membrane DOFs so the soft-membrane↔ERM↔stiff-cortex
        # coupling is solved consistently per fiber (see build_erm_augmented_fiber_block_kernel). Off by
        # default; the existing solvers are unchanged.
        self.use_erm_augmented_block = False
        self.erm_max_aug = 0
        # Tension-side ERM Jacobi diagonal WITHOUT block augmentation: give the membrane node its k_erm-scaled
        # scalar-Jacobi step (turgor/(a+γ_area+k_erm)) so the coupled membrane→ERM→cortex Jacobi iteration
        # converges (spectral radius ≈ k_erm/√(D_mem·D_cortex) ≪ 1) instead of over-stepping at
        # turgor/(a+γ_area). The membrane stays in scalar Jacobi (NOT coupled into the over-stepping fiber
        # block). Implied ON whenever the augmented block is used (its diagonals must match).
        self.erm_tension_side_precond = bool(erm_tension_side_precond)
        d = self.device
        with wp.ScopedDevice(d):
            self.r = wp.zeros(self.n, dtype=wp.vec3d)
            self.z = wp.zeros(self.n, dtype=wp.vec3d)
            self.projected_z = wp.zeros(self.n, dtype=wp.vec3d)
            self.p = wp.zeros(self.n, dtype=wp.vec3d)
            self.ap = wp.zeros(self.n, dtype=wp.vec3d)
            self.dx = wp.zeros(self.n, dtype=wp.vec3d)
            self.k_input = wp.zeros(self.n, dtype=wp.vec3d)
            self.k_output = wp.zeros(self.n, dtype=wp.vec3d)
            self.projected_k = wp.zeros(self.n, dtype=wp.vec3d)
            self.preconditioner = wp.zeros(self.n, dtype=wp.vec3d)
            self.coarse_external_diagonal = wp.zeros(self.n, dtype=wp.vec3d)
            factor_size = (
                int(cell.n_fibers) * 3 * self.max_fiber_nodes * self.max_fiber_nodes
                if self.use_fiber_block else 0
            )
            self.fiber_block_factor = wp.zeros(factor_size, dtype=wp.float64)
            self.fiber_block_work = wp.zeros(3 * self.n, dtype=wp.float64)
            self._build_erm_augmented_topology(cell, augment_fiber_block_erm)
            self.coarse_residual = wp.zeros(cell.n_fibers, dtype=wp.vec3d)
            self.coarse_diagonal = wp.zeros(cell.n_fibers, dtype=wp.vec3d)
            self.coarse_step = wp.zeros(cell.n_fibers, dtype=wp.vec3d)
            self.coarse_correction = wp.zeros(cell.n_fibers, dtype=wp.vec3d)
            self.coarse_action = wp.zeros(cell.n_fibers, dtype=wp.vec3d)
            self.coarse_fine = wp.zeros(self.n, dtype=wp.vec3d)
            self.coarse_fine_action = wp.zeros(self.n, dtype=wp.vec3d)
            self.diag_work = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64)
            self.rhs_work = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64)
            self.rr = wp.zeros(1, dtype=wp.float64)
            self.rr_initial = wp.zeros(1, dtype=wp.float64)
            self.rz = wp.zeros(1, dtype=wp.float64)
            self.rz_old = wp.zeros(1, dtype=wp.float64)
            self.p_ap = wp.zeros(1, dtype=wp.float64)
            self.alpha = wp.zeros(1, dtype=wp.float64)
            self.beta = wp.zeros(1, dtype=wp.float64)
            self.active = wp.zeros(1, dtype=wp.int32)
            self.converged = wp.zeros(1, dtype=wp.int32)
            self.iterations = wp.zeros(1, dtype=wp.int32)

            # Membrane in-plane area-tension tangent as EDGE SPRINGS (an exact analytic operator term replacing
            # the scalar-regularizer stand-in). Each triangle edge is a spring of stiffness γ_mem at its rest
            # length; at rest (rest == current edge) only the AXIAL stiffness survives (transverse tension = 0),
            # which is exactly the in-plane surface-tension tangent. So the analytic operator now sees the
            # membrane as a connected taut sheet — closing the MEMBRANE end of the coupled membrane→ERM→fiber
            # mode that per-fiber blocks miss. The force stays exact (``membrane_area_kernel``); this is
            # preconditioning only. Bending (κ̃) is still carried by the regularizer (a separate follow-up).
            self.mem_edges_d = None
            mem = getattr(cell, "membrane", None)
            if mem is not None and getattr(mem, "with_area_tension", False) and int(getattr(mem, "n_faces", 0)):
                faces = mem.faces_d.numpy()
                edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
                edges = np.unique(np.sort(edges, axis=1), axis=0)
                pos0 = cell.pos_d.numpy()
                rest = np.linalg.norm(pos0[edges[:, 0]] - pos0[edges[:, 1]], axis=1)
                self.mem_edges_d = wp.array(np.ascontiguousarray(edges, np.int32), dtype=wp.int32)
                self.mem_k_d = wp.array(np.full(edges.shape[0], float(mem.gamma_mem)), dtype=wp.float64)
                self.mem_rest_d = wp.array(np.ascontiguousarray(rest, np.float64), dtype=wp.float64)

            # Global rigid-body + constant-strain coarse space (l <= 2). The per-fiber block and node-local
            # Jacobi cannot see the whole-cell breathing / translation / ellipsoidal modes that couple the
            # membrane, ERM anchors, fibers and nucleus through one smooth global displacement — exactly the
            # residual component the preload stall carries as "global radial". A two-level additive Galerkin
            # correction ``B (B^T A B)^-1 B^T`` deflates that space each solve. ``A`` is rebuilt from the live
            # analytic operator (``coarse_modes`` operator applies), so this is preconditioning only: it cannot
            # move the fixed point or relax any residual gate, only accelerate the smooth modes. The basis is
            # geometry-derived (no tuned constant); orthonormalizing merely conditions ``A_c``.
            self.coarse_modes = 0
            self.coarse_basis_d = None
            requested_modes = min(int(coarse_modes), self.n)
            if requested_modes > 0:
                try:
                    basis = rigid_strain_coarse_basis(cell.pos_d.numpy(), n_modes=requested_modes)
                except ValueError as error:
                    _LOG.warning("global coarse space disabled (%s)", error)
                else:
                    self.coarse_modes = int(basis.shape[0])
                    self.coarse_basis_d = wp.array(
                        np.ascontiguousarray(basis, np.float64), dtype=wp.vec3d)
                    self.coarse_action_d = wp.zeros((self.coarse_modes, self.n), dtype=wp.vec3d)
                    self.coarse_matrix_d = wp.zeros(self.coarse_modes * self.coarse_modes,
                                                    dtype=wp.float64)
                    self.coarse_factor_d = wp.zeros(self.coarse_modes * self.coarse_modes,
                                                    dtype=wp.float64)
                    self.coarse_rhs_d = wp.zeros(self.coarse_modes, dtype=wp.float64)
                    self.coarse_solution_d = wp.zeros(self.coarse_modes, dtype=wp.float64)

            # Fiber-quotient inter-fiber coarse (Path A, matrix-free). The per-fiber block, translation
            # coarse, and global l<=2 space all miss the collective inter-fiber near-rigid modes that couple
            # fibers only through crosslinks — exactly the GATE-A resting plateau residual. This deflates that
            # space additively (M^-1 = S^-1 + P A_c^-1 P^T, SPD ⇒ PCG stays valid); it is preconditioning only
            # and defaults OFF. Enabled by the explicit flag OR the AC_FQ_COARSE=1 env var (the driver route,
            # which does not thread this flag). The prolongator is geometry-derived (no tuned constant).
            self.fq_coarse = None
            self._fq_reg = None
            # Mode: "pathA" = matrix-free P^T A P block-Jacobi inner CG (correctness-first, expensive SpMV);
            # "pathB" = explicit crosslink-weighted BSR A_c + DEEP block-Jacobi CG (cheap SpMV -> the strong
            # solve that closes the native plateau). Env overrides let the driver route (which does not thread
            # these) select without code changes: AC_FQ_MODE, AC_FQ_ITERS.
            self.fq_mode = os.environ.get("AC_FQ_MODE", str(fq_coarse_mode))
            fq_iters = int(os.environ.get("AC_FQ_ITERS", str(fq_coarse_iterations)))
            enable_fq = bool(fiber_quotient_coarse) or os.environ.get("AC_FQ_COARSE", "") == "1"
            if enable_fq and int(cell.n_fibers) > 0:
                from aleph.components.incumbent.fiber_quotient_coarse import (
                    FiberQuotientCoarse, FiberQuotientCoarsePathB)
                if self.fq_mode == "pathB":
                    self.fq_coarse = FiberQuotientCoarsePathB(cell, inner_iterations=fq_iters)
                else:
                    self.fq_coarse = FiberQuotientCoarse(cell, iterations=fq_iters)
                _LOG.info(
                    "fiber-quotient coarse ENABLED (mode=%s, N_c=%d, block_rank=%d, inner_iters=%d, ranks=%s)",
                    self.fq_mode, self.fq_coarse.n_coarse, getattr(self.fq_coarse, "block_rank",
                    self.fq_coarse.max_rank), fq_iters, self.fq_coarse.ranks_unique)

            # Fiber-arclength geometric-multigrid V-cycle preconditioner (design memo
            # FINE_MESH_MULTIGRID_DESIGN_2026-07-24 §2.2-2.4, §4.3). This is the ONE preconditioner that spans
            # the fine-mesh class-3 error (smooth-along-arclength AND coupled across fibers) that the per-fiber
            # block (no inter-fiber coupling) and the fiber-quotient rigid coarse (no arclength variation) both
            # miss. When ON it REPLACES the incumbent smoother+coarse stack in _precondition with a symmetric
            # V-cycle (SPD ⇒ the outer PCG is unchanged: preconditioning-only, the fixed point + residual gate
            # are IDENTICAL, CLAUDE.md no-gate-loosening). DEFAULT OFF; with it off every incumbent path below
            # is byte-untouched (self.mg stays None and the _precondition/build branches gate on self.multigrid).
            # Enabled by the explicit flag OR AC_MG=1 (the driver route, which does not thread this flag).
            self.multigrid = bool(multigrid) or os.environ.get("AC_MG", "") == "1"
            self.mg = None
            self._mg_reg = None
            if self.multigrid and int(cell.n_fibers) > 0:
                from aleph.components.incumbent.fiber_arclength_mg import FiberArclengthMultigrid
                mg_fq = bool(mg_fiber_quotient) and os.environ.get("AC_MG_FQ", "1") != "0"
                self.mg = FiberArclengthMultigrid(
                    cell, nu1=int(mg_smooth), nu2=int(mg_smooth), omega=float(mg_omega),
                    min_fiber_nodes=int(mg_min_fiber_nodes), max_levels=int(mg_max_levels),
                    attach_fiber_quotient=mg_fq, fq_inner_iterations=fq_iters,
                    assembled_coarse=(bool(mg_assembled_coarse) or os.environ.get("AC_MG_ASM", "") == "1"))
                _LOG.info("fiber-arclength MULTIGRID ENABLED (levels=%d, lengths=%s, coarsest=%s)",
                          self.mg.n_levels, self.mg.level_lengths,
                          "fiber-quotient" if self.mg.fq is not None else "line-smooth")
            elif self.multigrid:
                self.multigrid = False  # no fibers to precondition; leave the incumbent path untouched

    def _build_erm_augmented_topology(self, cell: object, requested: bool) -> None:
        """Group each bound ERM tether under the fiber owning its cortex node (static CSR, host, at setup).

        Each membrane node carries one ERM tether ⇒ every ERM-tethered membrane DOF is claimed by exactly one
        fiber block, disjoint across fibers, so the augmented block-diagonal preconditioner stays SPD. Live
        stiffness (position, bound flag, rupture) is evaluated in the build kernel; only the topology is fixed
        here. If augmentation is not requested, or there is no fiber block, or the membrane has no bound ERM,
        the plain per-fiber block is used unchanged.
        """
        self.erm_aug_off_d = None
        if not requested or not self.use_fiber_block:
            return
        mem = getattr(cell, "membrane", None)
        node_fiber = getattr(cell, "node_fiber_d", None)
        if mem is None or int(getattr(mem, "n_erm", 0)) == 0 or node_fiber is None:
            return
        n_fibers = int(cell.n_fibers)
        foff = np.asarray(cell.foff_d.numpy(), dtype=np.int64)
        node_fiber_h = np.asarray(node_fiber.numpy(), dtype=np.int64)  # (n_actin,) fiber id per actin node
        erm_m = np.asarray(mem.erm_m_d.numpy(), dtype=np.int64)        # global membrane node
        erm_c = np.asarray(mem.erm_c_d.numpy(), dtype=np.int64)        # global cortex (actin) node
        bound = np.asarray(mem.erm_bound_d.numpy(), dtype=np.int64).astype(bool)
        idx = np.nonzero(bound)[0]
        if idx.size == 0:
            return
        cortex = erm_c[idx]
        # cortex node must be an actin node with a valid fiber id
        valid = (cortex >= 0) & (cortex < node_fiber_h.shape[0])
        idx = idx[valid]
        cortex = cortex[valid]
        if idx.size == 0:
            return
        fiber_of = node_fiber_h[cortex]
        local_c = (cortex - foff[fiber_of]).astype(np.int32)
        order = np.argsort(fiber_of, kind="stable")
        fiber_of = fiber_of[order]
        erm_off = np.zeros(n_fibers + 1, dtype=np.int64)
        np.cumsum(np.bincount(fiber_of, minlength=n_fibers), out=erm_off[1:])
        per_fiber = np.diff(erm_off)
        max_erm_per_fiber = int(per_fiber.max()) if per_fiber.size else 0
        max_aug = int(self.max_fiber_nodes + max_erm_per_fiber)
        aug_local_c = local_c[order].astype(np.int32)
        aug_mem_global = erm_m[idx][order].astype(np.int32)
        aug_index = idx[order].astype(np.int32)
        d = self.device
        self.erm_aug_off_d = wp.array(np.ascontiguousarray(erm_off, np.int32), dtype=wp.int32, device=d)
        self.erm_aug_local_c_d = wp.array(np.ascontiguousarray(aug_local_c), dtype=wp.int32, device=d)
        self.erm_aug_mem_global_d = wp.array(np.ascontiguousarray(aug_mem_global), dtype=wp.int32, device=d)
        self.erm_aug_index_d = wp.array(np.ascontiguousarray(aug_index), dtype=wp.int32, device=d)
        self.erm_max_aug = max_aug
        self.erm_aug_factor = wp.zeros(n_fibers * 3 * max_aug * max_aug, dtype=wp.float64, device=d)
        self.use_erm_augmented_block = True
        self.erm_aug_n_tethers = int(idx.size)
        self.erm_aug_max_per_fiber = max_erm_per_fiber

    def _project(self, pos: wp.array, source: wp.array, destination: wp.array,
                 finite: wp.array) -> None:
        wp.copy(destination, source)
        if self.cell.n_fibers:
            from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
            wp.launch(project_constraint_forces_kernel, dim=self.cell.n_fibers,
                      inputs=[pos, source, self.cell.foff_d, self.cell.soff_d, destination,
                              self.diag_work, self.rhs_work, finite], device=self.device)

    def _stiffness(self, pos: wp.array, vector: wp.array, out: wp.array,
                   regularization: wp.array) -> None:
        cell = self.cell
        d = self.device
        wp.launch(set_mass_action_kernel, dim=self.n, inputs=[vector, regularization, out], device=d)
        if cell.n_tri:
            wp.launch(add_bending_stiffness_kernel, dim=cell.n_tri,
                      inputs=[vector, cell.tri_d, cell.alpha_d, out], device=d)
        if cell.n_xl:
            wp.launch(add_pair_stiffness_kernel, dim=cell.n_xl,
                      inputs=[pos, vector, cell.xl_d, cell.kxl_d, cell.r0xl_d, out], device=d)
        if self.mem_edges_d is not None:
            wp.launch(add_pair_stiffness_kernel, dim=self.mem_edges_d.shape[0],
                      inputs=[pos, vector, self.mem_edges_d, self.mem_k_d, self.mem_rest_d, out], device=d)
        steric = cell.steric
        if steric is not None:
            wp.launch(add_wca_stiffness_kernel, dim=steric.n,
                      inputs=[steric.grid.id, steric._qpts, pos, vector, steric.fiber_id, steric.active,
                              wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                              wp.float64(steric.f_cap), out], device=d)
        myosin = cell.myosin
        if myosin is not None:
            if int(myosin.backbone_bonds.shape[0]):
                wp.launch(add_uniform_pair_stiffness_kernel, dim=myosin.backbone_bonds.shape[0],
                          inputs=[pos, vector, myosin.backbone_bonds, myosin.k_backbone,
                                  myosin.r0_backbone, out], device=d)
            if int(myosin.head_bonds.shape[0]):
                wp.launch(add_uniform_pair_stiffness_kernel, dim=myosin.head_bonds.shape[0],
                          inputs=[pos, vector, myosin.head_bonds, myosin.k_head_spring,
                                  myosin.params.r0_head, out], device=d)
            if myosin.backbone_angles is not None and float(myosin.k_theta_backbone) > 0.0:
                wp.launch(add_angle_gauss_newton_stiffness_kernel, dim=myosin.backbone_angles.shape[0],
                          inputs=[pos, vector, myosin.backbone_angles, myosin.k_theta_backbone,
                                  wp.float64(np.pi), out], device=d)
            if myosin.head_arm_angles is not None and float(myosin.k_theta_arm) > 0.0:
                wp.launch(add_angle_gauss_newton_stiffness_kernel, dim=myosin.head_arm_angles.shape[0],
                          inputs=[pos, vector, myosin.head_arm_angles, myosin.k_theta_arm,
                                  wp.float64(0.5 * np.pi), out], device=d)
            state = myosin.state
            if int(state["bound"].shape[0]):
                if myosin.segment_runtime is None:
                    wp.launch(add_crossbridge_stiffness_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, vector, myosin.head_node, state["bound"], state["anchor"],
                                      state["abscissa"], state["walk_dir"], myosin.params.k_xb,
                                      myosin.params.r0_xb, out], device=d)
                else:
                    wp.launch(add_segment_crossbridge_stiffness_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, vector, myosin.head_node, state["bound"], state["seg_a"],
                                      state["seg_b"], state["bary_t"], state["abscissa"], state["walk_dir"],
                                      myosin.params.k_xb, myosin.params.r0_xb, out], device=d)
        nucleus = cell.nucleus
        if nucleus is not None and nucleus.n_linc:
            wp.launch(add_linc_stiffness_kernel, dim=nucleus.n_linc,
                      inputs=[pos, vector, nucleus.linc_n_d, nucleus.linc_a_d, nucleus.linc_rest_d,
                              wp.float64(nucleus.k_linc), wp.float64(nucleus.linc_stiffening), out], device=d)
        membrane = cell.membrane
        if membrane is not None and membrane.n_erm:
            wp.launch(add_erm_stiffness_kernel, dim=membrane.n_erm,
                      inputs=[pos, vector, membrane.erm_m_d, membrane.erm_c_d, membrane.erm_bound_d,
                              membrane.erm_rest_d, wp.float64(membrane.k_erm),
                              wp.float64(membrane.f_rupt), out], device=d)

    def _build_preconditioner(
        self, pos: wp.array, regularization: wp.array, finite: wp.array,
    ) -> None:
        """Assemble a scalar block-majorizer of the represented analytic stiffness families."""
        cell = self.cell
        d = self.device
        wp.launch(set_preconditioner_kernel, dim=self.n,
                  inputs=[regularization, self.preconditioner], device=d)
        wp.launch(set_preconditioner_kernel, dim=self.n,
                  inputs=[regularization, self.coarse_external_diagonal], device=d)
        if cell.n_tri:
            wp.launch(add_bending_preconditioner_kernel, dim=cell.n_tri,
                      inputs=[cell.tri_d, cell.alpha_d, self.preconditioner], device=d)
        if cell.n_xl:
            wp.launch(add_pair_preconditioner_kernel, dim=cell.n_xl,
                      inputs=[pos, cell.xl_d, cell.kxl_d, cell.r0xl_d, self.preconditioner], device=d)
            wp.launch(add_pair_preconditioner_kernel, dim=cell.n_xl,
                      inputs=[pos, cell.xl_d, cell.kxl_d, cell.r0xl_d,
                              self.coarse_external_diagonal], device=d)
        if self.mem_edges_d is not None:
            # The membrane in-plane area-tension edge springs are in the OPERATOR (see _stiffness) but were
            # previously absent from the preconditioner diagonal, so membrane nodes were under-stiffened and
            # over-stepped (contributing to the ~0.38 pN tangential membrane floor, 2026-07-22e). Adding their
            # PSD Cartesian diagonal makes the membrane node diagonal consistent with the operator; the
            # backbone-aware augmented block then reads it as each ERM-tethered membrane node's full diagonal.
            wp.launch(add_pair_preconditioner_kernel, dim=self.mem_edges_d.shape[0],
                      inputs=[pos, self.mem_edges_d, self.mem_k_d, self.mem_rest_d,
                              self.preconditioner], device=d)
        steric = cell.steric
        if steric is not None:
            wp.launch(add_wca_preconditioner_kernel, dim=steric.n,
                      inputs=[steric.grid.id, steric._qpts, pos, steric.fiber_id, steric.active,
                              wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                              wp.float64(steric.f_cap), self.preconditioner], device=d)
            wp.launch(add_wca_preconditioner_kernel, dim=steric.n,
                      inputs=[steric.grid.id, steric._qpts, pos, steric.fiber_id, steric.active,
                              wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                              wp.float64(steric.f_cap), self.coarse_external_diagonal], device=d)
        myosin = cell.myosin
        if myosin is not None:
            if int(myosin.backbone_bonds.shape[0]):
                wp.launch(add_uniform_pair_preconditioner_kernel, dim=myosin.backbone_bonds.shape[0],
                          inputs=[pos, myosin.backbone_bonds, myosin.k_backbone, myosin.r0_backbone,
                                  self.preconditioner], device=d)
            if int(myosin.head_bonds.shape[0]):
                wp.launch(add_uniform_pair_preconditioner_kernel, dim=myosin.head_bonds.shape[0],
                          inputs=[pos, myosin.head_bonds, myosin.k_head_spring, myosin.params.r0_head,
                                  self.preconditioner], device=d)
            if myosin.backbone_angles is not None and float(myosin.k_theta_backbone) > 0.0:
                wp.launch(add_angle_gauss_newton_preconditioner_kernel,
                          dim=myosin.backbone_angles.shape[0],
                          inputs=[pos, myosin.backbone_angles, myosin.k_theta_backbone,
                                  wp.float64(np.pi), self.preconditioner], device=d)
            if myosin.head_arm_angles is not None and float(myosin.k_theta_arm) > 0.0:
                wp.launch(add_angle_gauss_newton_preconditioner_kernel,
                          dim=myosin.head_arm_angles.shape[0],
                          inputs=[pos, myosin.head_arm_angles, myosin.k_theta_arm,
                                  wp.float64(0.5 * np.pi), self.preconditioner], device=d)
            state = myosin.state
            if int(state["bound"].shape[0]):
                if myosin.segment_runtime is None:
                    wp.launch(add_crossbridge_preconditioner_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, myosin.head_node, state["bound"], state["anchor"], state["abscissa"],
                                      state["walk_dir"], myosin.params.k_xb, myosin.params.r0_xb,
                                      self.preconditioner], device=d)
                    wp.launch(add_crossbridge_preconditioner_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, myosin.head_node, state["bound"], state["anchor"], state["abscissa"],
                                      state["walk_dir"], myosin.params.k_xb, myosin.params.r0_xb,
                                      self.coarse_external_diagonal], device=d)
                else:
                    wp.launch(add_segment_crossbridge_preconditioner_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, myosin.head_node, state["bound"], state["seg_a"], state["seg_b"],
                                      state["bary_t"], state["abscissa"], state["walk_dir"],
                                      myosin.params.k_xb, myosin.params.r0_xb, self.preconditioner], device=d)
                    wp.launch(add_segment_crossbridge_preconditioner_kernel, dim=state["bound"].shape[0],
                              inputs=[pos, myosin.head_node, state["bound"], state["seg_a"], state["seg_b"],
                                      state["bary_t"], state["abscissa"], state["walk_dir"],
                                      myosin.params.k_xb, myosin.params.r0_xb,
                                      self.coarse_external_diagonal], device=d)
        nucleus = cell.nucleus
        if nucleus is not None and nucleus.n_linc:
            wp.launch(add_linc_preconditioner_kernel, dim=nucleus.n_linc,
                      inputs=[pos, nucleus.linc_n_d, nucleus.linc_a_d, nucleus.linc_rest_d,
                              wp.float64(nucleus.k_linc), wp.float64(nucleus.linc_stiffening),
                              self.preconditioner], device=d)
            wp.launch(add_linc_preconditioner_kernel, dim=nucleus.n_linc,
                      inputs=[pos, nucleus.linc_n_d, nucleus.linc_a_d, nucleus.linc_rest_d,
                              wp.float64(nucleus.k_linc), wp.float64(nucleus.linc_stiffening),
                              self.coarse_external_diagonal], device=d)
        membrane = cell.membrane
        if membrane is not None and membrane.n_erm:
            # The augmented-block path uses the TENSION-SIDE ERM diagonal (engages the axial stiffness at the
            # force-free resting extension==0) so the membrane→ERM→cortex step is scaled by k_erm and the
            # 2-node ERM sub-blocks stay diagonally dominant; every other solver keeps the unilateral
            # (extension>0) diagonal that mirrors the ERM force exactly.
            erm_precond_kernel = (
                add_erm_tension_side_preconditioner_kernel
                if (self.use_erm_augmented_block or self.erm_tension_side_precond)
                else add_erm_preconditioner_kernel)
            wp.launch(erm_precond_kernel, dim=membrane.n_erm,
                      inputs=[pos, membrane.erm_m_d, membrane.erm_c_d, membrane.erm_bound_d,
                              membrane.erm_rest_d, wp.float64(membrane.k_erm), wp.float64(membrane.f_rupt),
                              self.preconditioner], device=d)
            wp.launch(erm_precond_kernel, dim=membrane.n_erm,
                      inputs=[pos, membrane.erm_m_d, membrane.erm_c_d, membrane.erm_bound_d,
                              membrane.erm_rest_d, wp.float64(membrane.k_erm), wp.float64(membrane.f_rupt),
                              self.coarse_external_diagonal], device=d)
        if self.use_erm_augmented_block:
            mem = cell.membrane
            wp.launch(build_erm_augmented_fiber_block_kernel, dim=cell.n_fibers,
                      inputs=[pos, self.coarse_external_diagonal, self.preconditioner,
                              cell.foff_d, cell.toff_d, cell.alpha_d,
                              self.erm_aug_off_d, self.erm_aug_local_c_d, self.erm_aug_mem_global_d,
                              self.erm_aug_index_d, mem.erm_bound_d, mem.erm_rest_d,
                              wp.float64(mem.k_erm), wp.float64(mem.f_rupt),
                              wp.int32(self.erm_max_aug), self.erm_aug_factor, finite], device=d)
        elif self.use_fiber_block and not self.multigrid:
            wp.launch(build_fiber_block_cholesky_kernel, dim=cell.n_fibers,
                      inputs=[self.coarse_external_diagonal, cell.foff_d, cell.toff_d,
                              cell.alpha_d, wp.int32(self.max_fiber_nodes),
                              self.fiber_block_factor, finite], device=d)
        if self.multigrid:
            # Stash the live regularizer (the V-cycle's matrix-free coarse operators need it) and refresh the
            # per-level smoother diagonals + coarsest A_c. The scalar node-Jacobi diagonal (self.preconditioner)
            # is still assembled above so the caller can precondition the non-actin DOFs.
            self._mg_reg = regularization
            self.mg.build(pos, self.coarse_external_diagonal, regularization, finite, self._operator)
        if self.coarse_modes:
            self._build_coarse_operator(pos, regularization, finite)
        if self.fq_coarse is not None:
            # Stash the live regularizer (Path A's matrix-free SpMV in _precondition has no reg argument).
            self._fq_reg = regularization
            if self.fq_mode == "pathB":
                # Assemble the crosslink-weighted BSR A_c at the live position + factor its block-Jacobi.
                self.fq_coarse.build(pos, regularization, finite)
            else:
                # Path A: factor the per-fiber block-Jacobi from the assembled coarse_external_diagonal.
                self.fq_coarse.build(self.coarse_external_diagonal, finite)

    def _build_coarse_operator(
        self, pos: wp.array, regularization: wp.array, finite: wp.array,
    ) -> None:
        """Form and factor the Galerkin coarse operator ``A_c = B^T A B`` at the live position."""
        d = self.device
        for mode in range(self.coarse_modes):
            self._operator(pos, self.coarse_basis_d[mode], self.coarse_action_d[mode],
                           regularization, finite)
        wp.launch(assemble_coarse_matrix_kernel, dim=self.coarse_modes * self.coarse_modes,
                  inputs=[self.coarse_basis_d, self.coarse_action_d, wp.int32(self.n),
                          wp.int32(self.coarse_modes), self.coarse_matrix_d], device=d)
        wp.launch(factor_coarse_cholesky_kernel, dim=1,
                  inputs=[self.coarse_matrix_d, wp.int32(self.coarse_modes),
                          self.coarse_factor_d, finite], device=d)

    def _precondition(self, pos: wp.array, finite: wp.array) -> None:
        wp.launch(apply_preconditioner_kernel, dim=self.n,
                  inputs=[self.r, self.preconditioner, self.z], device=self.device)
        if self.multigrid:
            # M^-1 = blockdiag(V-cycle on the actin fiber block, node-Jacobi on the disjoint remainder). The
            # node-Jacobi apply above already filled z everywhere (giving the non-actin DOFs their diagonal
            # preconditioner); the V-cycle now OVERWRITES the leading n_actin block with its SPD correction.
            self.mg.apply(pos, self.r, self.z, self._operator, self._mg_reg, finite)
            self._project(pos, self.z, self.projected_z, finite)
            return
        if self.use_erm_augmented_block:
            # Overwrite the scalar-Jacobi values on both the actin fiber nodes AND their ERM-tethered membrane
            # nodes with the exact augmented-block Cholesky solve (disjoint across fibers ⇒ no race, SPD
            # block-diagonal preconditioner). Membrane nodes with no ERM keep their scalar-Jacobi value.
            wp.launch(apply_erm_augmented_fiber_block_kernel, dim=self.cell.n_fibers,
                      inputs=[self.r, self.cell.foff_d, self.erm_aug_off_d, self.erm_aug_mem_global_d,
                              wp.int32(self.erm_max_aug), self.erm_aug_factor, self.fiber_block_work,
                              finite, self.z], device=self.device)
            wp.launch(add_fiber_translation_coarse_kernel, dim=self.cell.n_fibers,
                      inputs=[self.r, self.cell.foff_d, self.coarse_external_diagonal, self.z],
                      device=self.device)
        elif self.use_fiber_block:
            wp.launch(apply_fiber_block_preconditioner_kernel, dim=self.cell.n_fibers,
                      inputs=[self.r, self.cell.foff_d, wp.int32(self.max_fiber_nodes),
                              self.fiber_block_factor, self.fiber_block_work, finite, self.z],
                      device=self.device)
            wp.launch(add_fiber_translation_coarse_kernel, dim=self.cell.n_fibers,
                      inputs=[self.r, self.cell.foff_d, self.coarse_external_diagonal, self.z],
                      device=self.device)
        if self.coarse_modes:
            wp.launch(coarse_restrict_kernel, dim=self.coarse_modes,
                      inputs=[self.coarse_basis_d, self.r, wp.int32(self.n), self.coarse_rhs_d],
                      device=self.device)
            wp.launch(coarse_solve_kernel, dim=1,
                      inputs=[self.coarse_factor_d, self.coarse_rhs_d, wp.int32(self.coarse_modes),
                              finite, self.coarse_solution_d], device=self.device)
            wp.launch(coarse_prolong_add_kernel, dim=self.n,
                      inputs=[self.coarse_basis_d, self.coarse_solution_d,
                              wp.int32(self.coarse_modes), self.z], device=self.device)
        if self.fq_coarse is not None:
            # Additive fiber-quotient coarse (restrict self.r → inner block-Jacobi CG → prolong-add into
            # self.z), BEFORE the final projection. Additive with the smoother ⇒ composite stays SPD.
            if self.fq_mode == "pathB":
                # SpMV = the assembled BSR A_c (cheap ⇒ deep budget); no fine operator apply.
                self.fq_coarse.apply(pos, self.r, self.z, self._fq_reg, finite)
            else:
                # SpMV = P^T A P via the exact fine operator.
                self.fq_coarse.apply(pos, self.r, self.z, self._operator, self._fq_reg, finite)
        self._project(pos, self.z, self.projected_z, finite)

    def _operator(self, pos: wp.array, vector: wp.array, out: wp.array,
                  regularization: wp.array, finite: wp.array) -> None:
        self._project(pos, vector, self.k_input, finite)
        self._stiffness(pos, self.k_input, self.k_output, regularization)
        self._project(pos, self.k_output, self.projected_k, finite)
        # projected_k already contains a*P(v); add the tiny difference a*(v-Pv) so the exact operator is
        # aI + PKP even if a CG vector has accumulated round-off outside the tangent.
        wp.launch(_combine_operator_kernel, dim=self.n,
                  inputs=[vector, self.k_input, self.projected_k, regularization, out], device=self.device)

    def preconditioned_direction(
        self, pos: wp.array, rhs: wp.array, regularization: wp.array, finite: wp.array,
    ) -> wp.array:
        r"""Return ``P M^-1 rhs`` without interpreting the preconditioner as a converged linear solve.

        This is the source-independent nonlinear block-descent direction. ``M`` is rebuilt from the same live
        analytic tangents as PCG, including the exact short-fiber bending Cholesky block and additive fiber
        translation mode. The caller must globalize the direction with an exact nonlinear residual comparison;
        no linear-convergence latch is implied or required.
        """
        self._build_preconditioner(pos, regularization, finite)
        wp.copy(self.r, rhs)
        self._precondition(pos, finite)
        return self.projected_z

    def solve(self, pos: wp.array, rhs: wp.array, regularization: wp.array,
              finite: wp.array) -> wp.array:
        """Return the device displacement after the fixed CG launch budget."""
        d = self.device
        self._build_preconditioner(pos, regularization, finite)
        wp.launch(_copy_zero_dx_kernel, dim=self.n, inputs=[rhs, self.r, self.dx], device=d)
        if self.use_fiber_block and self.coarse_iterations:
            wp.launch(restrict_fiber_sum_kernel, dim=self.cell.n_fibers,
                      inputs=[rhs, self.cell.foff_d, self.coarse_residual], device=d)
            wp.launch(restrict_fiber_trace_majorizer_kernel, dim=self.cell.n_fibers,
                      inputs=[self.coarse_external_diagonal, self.cell.foff_d,
                              self.coarse_diagonal], device=d)
            self.coarse_correction.zero_()
            for _ in range(self.coarse_iterations):
                wp.launch(coarse_jacobi_step_kernel, dim=self.cell.n_fibers,
                          inputs=[self.coarse_residual, self.coarse_diagonal,
                                  self.coarse_step], device=d)
                wp.launch(add_coarse_step_kernel, dim=self.cell.n_fibers,
                          inputs=[self.coarse_step, self.coarse_correction], device=d)
                self.coarse_fine.zero_()
                wp.launch(prolong_fiber_translation_kernel, dim=self.cell.n_fibers,
                          inputs=[self.coarse_step, self.cell.foff_d, self.coarse_fine], device=d)
                self._operator(pos, self.coarse_fine, self.coarse_fine_action,
                               regularization, finite)
                wp.launch(restrict_fiber_sum_kernel, dim=self.cell.n_fibers,
                          inputs=[self.coarse_fine_action, self.cell.foff_d,
                                  self.coarse_action], device=d)
                wp.launch(subtract_coarse_action_kernel, dim=self.cell.n_fibers,
                          inputs=[self.coarse_action, self.coarse_residual], device=d)
            self.coarse_fine.zero_()
            wp.launch(prolong_fiber_translation_kernel, dim=self.cell.n_fibers,
                      inputs=[self.coarse_correction, self.cell.foff_d, self.coarse_fine], device=d)
            wp.copy(self.dx, self.coarse_fine)
            self._operator(pos, self.dx, self.coarse_fine_action, regularization, finite)
            wp.launch(set_initial_residual_kernel, dim=self.n,
                      inputs=[rhs, self.coarse_fine_action, self.r], device=d)
        self._precondition(pos, finite)
        wp.copy(self.p, self.projected_z)
        self.rr.zero_()
        wp.launch(_dot_kernel, dim=self.n, inputs=[self.r, self.r, self.rr], device=d)
        self.rz.zero_()
        wp.launch(_dot_kernel, dim=self.n, inputs=[self.r, self.projected_z, self.rz], device=d)
        wp.launch(_cg_begin_kernel, dim=1,
                  inputs=[self.rr, self.rz, self.rr_initial, self.rz_old,
                          self.active, self.converged, finite, self.iterations], device=d)
        for iteration in range(1, self.max_iterations + 1):
            self._operator(pos, self.p, self.ap, regularization, finite)
            self.p_ap.zero_()
            wp.launch(_dot_kernel, dim=self.n, inputs=[self.p, self.ap, self.p_ap], device=d)
            wp.launch(_cg_alpha_kernel, dim=1,
                      inputs=[self.rz_old, self.p_ap, self.active, finite, self.alpha], device=d)
            wp.launch(_cg_update_kernel, dim=self.n,
                      inputs=[self.dx, self.r, self.p, self.ap, self.alpha, self.active], device=d)
            self.rr.zero_()
            wp.launch(_dot_kernel, dim=self.n, inputs=[self.r, self.r, self.rr], device=d)
            wp.launch(_cg_finalize_kernel, dim=1,
                      inputs=[self.rr, self.rr_initial, wp.float64(2.220446049250313e-16),
                              wp.int32(iteration), self.active, self.converged, finite, self.iterations],
                      device=d)
            self._precondition(pos, finite)
            self.rz.zero_()
            wp.launch(_dot_kernel, dim=self.n, inputs=[self.r, self.projected_z, self.rz], device=d)
            wp.launch(_cg_beta_kernel, dim=1,
                      inputs=[self.rz, self.rz_old, self.active, finite, self.beta], device=d)
            wp.launch(_cg_direction_kernel, dim=self.n,
                      inputs=[self.projected_z, self.beta, self.active, self.p], device=d)
            # SPEED (opt-in, bit-identical): once the device latch clears (converged or non-finite), every further
            # loop body only recomputes no-op updates (all _cg_* kernels gate on `active`). Break to skip those
            # launches — each a full MG V-cycle. dx is unchanged; only the wasted post-convergence work is saved.
            if self.cg_check_every and iteration % self.cg_check_every == 0:
                if int(self.active.numpy()[0]) == 0:
                    break
        return self.dx


@wp.kernel
def _combine_operator_kernel(
    vector: wp.array(dtype=wp.vec3d),
    projected_vector: wp.array(dtype=wp.vec3d),
    projected_mass_and_stiffness: wp.array(dtype=wp.vec3d),
    regularization: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    out[i] = projected_mass_and_stiffness[i] + regularization[0] * (vector[i] - projected_vector[i])
