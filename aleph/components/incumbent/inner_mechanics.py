"""Device-side convergence machinery for the Active Cell inner mechanical solve.

The old integrator ran a fixed number of descent iterations and called the result "stable" without asking
whether the projected mechanics had reached a fixed point.  This module supplies whole-cell reductions and
convergence-guarded position/reshape kernels.  Host code may launch a fixed *maximum* number of iterations, but
after the device convergence flag is set every subsequent update is a no-op; no state-dependent D2H read is
needed to stop the mechanics.

Convergence requires the actual post-projection maximum nodal displacement, maximum NF2007 segment-length
constraint residual, and explicit-equivalent projected-force displacement ``dt_mu * max|P F|`` to fall below
tolerance, with every position/force component finite. An unconstrained force can coexist with a stationary
constrained state, so the force test uses ``P F`` rather than raw force. The tolerance is
``sqrt(float64 epsilon) * ell_ref`` where ``ell_ref`` is the actin segment discretization length. This is a
derived numerical-resolution criterion, grid-invariant and not tuned to make a biological gate pass. The
descent step also uses the live device maximum pressure jump, so the pressure-traction Jacobian cannot outrun
a CFL value frozen at resting Π₀.

Sanity Gate:
    * Dimensions: displacement/tolerance [um]; ``dt_mu [um/pN] * projected force [pN]`` is [um].
    * Boundary cases: zero force yields zero displacement and converges in one iteration; a non-zero update
      cannot set the flag; once set, both descent and reshape leave positions bit-identical.
    * Invariants: conditional reshape is the NF2007 COG-conserving projection copied exactly from the retained
      Warp component, with only the device convergence guard added.
    * Numerical: ``sqrt(eps64)`` is derived from machine precision; ``max_inner`` is a budget and a run that
      exhausts it reports ``converged=False`` rather than weakening the criterion.  Additional retry chunks
      extend only this numerical budget; they never advance physical time.
    * Residency: reductions, flag and first-converged iteration remain device arrays; no inner-loop readback.
"""

from __future__ import annotations

import warp as wp

_FLUID = wp.constant(1)

__all__ = [
    "begin_inner_attempt_kernel",
    "conditional_axpy_kernel",
    "conditional_rollback_vec3_kernel",
    "conditional_reshape_kernel",
    "compute_descent_step_kernel",
    "finite_vec3_kernel",
    "fire_adapt_kernel",
    "fire_mix_and_step_kernel",
    "fire_power_kernel",
    "fire_reset_state_kernel",
    "fire_sumsq_kernel",
    "invalidate_convergence_kernel",
    "inner_state_init_kernel",
    "max_abs_pressure_kernel",
    "max_constraint_error_kernel",
    "max_displacement_kernel",
    "max_force_kernel",
    "project_constraint_forces_kernel",
    "record_convergence_history_kernel",
    "reset_active_reduction_kernel",
    "update_convergence_kernel",
]

# ---------------------------------------------------------------------------
# FIRE control constants (Bitzek et al. 2006, PRL 97 170201 -- the published defaults).
# FIRE accelerates the SAME projected explicit descent (``pos += dt_mu * P F``) with a
# fictitious-inertia velocity and adaptive damping. It uses only force evaluations: no global
# tangent, no linear solve, no preconditioner -- so it never introduces the coupled-preload
# conditioning of a tangent accelerator, and it obeys the identical strict projected-force
# convergence gate (it cannot manufacture convergence). These five constants are universal,
# grid-invariant, and not tuned to pass any biological gate (Magic-Number Block: sourced
# closed-form algorithm constants, not empirical tuning knobs).
# ---------------------------------------------------------------------------
_FIRE_N_MIN = wp.constant(5)
_FIRE_F_INC = wp.constant(wp.float64(1.1))
_FIRE_F_DEC = wp.constant(wp.float64(0.5))
_FIRE_ALPHA0 = wp.constant(wp.float64(0.1))
_FIRE_F_ALPHA = wp.constant(wp.float64(0.99))
_FIRE_DT_MAX_MULT = wp.constant(wp.float64(10.0))


@wp.kernel
def project_constraint_forces_kernel(
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    seg_off: wp.array(dtype=wp.int32),
    projected: wp.array(dtype=wp.vec3d),
    diag_work: wp.array(dtype=wp.float64),
    rhs_lambda_work: wp.array(dtype=wp.float64),
    finite: wp.array(dtype=wp.int32),
) -> None:
    """Apply the exact NF2007 inextensibility projector to one variable-length fiber.

    For ``P = I - J^T (J J^T)^-1 J``, ``J J^T`` is symmetric tridiagonal. One CUDA thread owns one
    fiber and performs a Thomas solve in its disjoint slices of the global work arrays. ``projected`` must
    first be a D2D copy of ``force`` so non-fiber nodes pass through unchanged; this kernel overwrites every
    fiber node with ``P force``. A non-positive pivot is impossible for a nondegenerate fiber and latches the
    existing device finite/valid predicate instead of silently emitting an unconstrained update.
    """
    f = wp.tid()
    a = fiber_off[f]
    b = fiber_off[f + 1]
    p = b - a - 1
    s0 = seg_off[f]
    if p <= 0:
        return

    # Thomas forward elimination. ``rhs_lambda_work`` begins as J force and is overwritten by lambda below.
    for k in range(p):
        g = pos[a + k + 1] - pos[a + k]
        diagonal = wp.float64(8.0) * wp.dot(g, g)
        rhs = wp.float64(2.0) * wp.dot(g, force[a + k + 1] - force[a + k])
        if k > 0:
            g_prev = pos[a + k] - pos[a + k - 1]
            off = -wp.float64(4.0) * wp.dot(g_prev, g)
            previous_pivot = diag_work[s0 + k - 1]
            if previous_pivot <= wp.float64(1.0e-300):
                finite[0] = 0
                return
            multiplier = off / previous_pivot
            diagonal = diagonal - multiplier * off
            rhs = rhs - multiplier * rhs_lambda_work[s0 + k - 1]
        if diagonal <= wp.float64(1.0e-300):
            finite[0] = 0
            return
        diag_work[s0 + k] = diagonal
        rhs_lambda_work[s0 + k] = rhs

    # Back substitution: the work RHS becomes lambda = (J J^T)^-1 J force.
    last = p - 1
    rhs_lambda_work[s0 + last] = rhs_lambda_work[s0 + last] / diag_work[s0 + last]
    for reverse in range(p - 1):
        k = p - 2 - reverse
        g = pos[a + k + 1] - pos[a + k]
        g_next = pos[a + k + 2] - pos[a + k + 1]
        off = -wp.float64(4.0) * wp.dot(g, g_next)
        rhs_lambda_work[s0 + k] = (
            rhs_lambda_work[s0 + k] - off * rhs_lambda_work[s0 + k + 1]
        ) / diag_work[s0 + k]

    # J^T lambda at node j receives -2 g_j lambda_j + 2 g_{j-1} lambda_{j-1}.
    for j in range(p + 1):
        reaction = wp.vec3d(0.0, 0.0, 0.0)
        if j < p:
            g_right = pos[a + j + 1] - pos[a + j]
            reaction = reaction - wp.float64(2.0) * g_right * rhs_lambda_work[s0 + j]
        if j > 0:
            g_left = pos[a + j] - pos[a + j - 1]
            reaction = reaction + wp.float64(2.0) * g_left * rhs_lambda_work[s0 + j - 1]
        projected[a + j] = force[a + j] - reaction


@wp.kernel
def begin_inner_attempt_kernel(active: wp.array(dtype=wp.int32), attempts: wp.array(dtype=wp.int32)) -> None:
    """Count a numerical retry chunk only while the device solve remains active."""
    if active[0] != 0:
        attempts[0] += 1


@wp.kernel
def reset_active_reduction_kernel(
    value: wp.array(dtype=wp.float64), active: wp.array(dtype=wp.int32),
) -> None:
    """Reset a convergence diagnostic only while a retry chunk is still active."""
    if active[0] != 0:
        value[0] = wp.float64(0.0)


@wp.kernel
def inner_state_init_kernel(
    active: wp.array(dtype=wp.int32),
    converged: wp.array(dtype=wp.int32),
    iters: wp.array(dtype=wp.int32),
    attempts: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    max_total_inner: wp.int32,
    fluid_valid: wp.array(dtype=wp.int32),
) -> None:
    """Reset convergence state, disabling every position/projection update when fluid input is invalid."""
    finite[0] = fluid_valid[0]
    converged[0] = 0
    attempts[0] = 0
    if fluid_valid[0] == 1:
        active[0] = 1
        iters[0] = max_total_inner
    else:
        active[0] = 0
        iters[0] = 0


@wp.kernel
def max_force_kernel(force: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.float64)) -> None:
    """Whole-cell ``max_i |F_i|`` reduction [pN]."""
    i = wp.tid()
    wp.atomic_max(out, 0, wp.length(force[i]))


@wp.kernel
def max_displacement_kernel(
    before: wp.array(dtype=wp.vec3d),
    after: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.float64),
) -> None:
    """Whole-cell projected-step ``max_i |x_i^{k+1}-x_i^k|`` reduction [um]."""
    i = wp.tid()
    wp.atomic_max(out, 0, wp.length(after[i] - before[i]))


@wp.kernel
def max_constraint_error_kernel(
    pos: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    seg_off: wp.array(dtype=wp.int32),
    seg_rest: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.float64),
) -> None:
    """Maximum NF2007 segment-length constraint residual [um]."""
    f = wp.tid()
    a = fiber_off[f]
    b = fiber_off[f + 1]
    s0 = seg_off[f]
    for k in range(b - a - 1):
        err = wp.abs(wp.length(pos[a + k + 1] - pos[a + k]) - seg_rest[s0 + k])
        wp.atomic_max(out, 0, err)


@wp.kernel
def finite_vec3_kernel(values: wp.array(dtype=wp.vec3d), finite: wp.array(dtype=wp.int32)) -> None:
    """Latch ``finite=0`` if any component is NaN or infinite."""
    i = wp.tid()
    v = values[i]
    if not wp.isfinite(v[0]) or not wp.isfinite(v[1]) or not wp.isfinite(v[2]):
        finite[0] = 0


@wp.kernel
def invalidate_convergence_kernel(
    finite: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    converged: wp.array(dtype=wp.int32),
) -> None:
    """A non-finite final force/state invalidates an earlier displacement convergence latch."""
    if finite[0] == 0:
        active[0] = 0
        converged[0] = 0


@wp.kernel
def max_abs_pressure_kernel(
    pressure: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    p_ext: wp.float64,
    out: wp.array(dtype=wp.float64),
) -> None:
    """Maximum live pressure jump over FLUID cells [Pa == pN/um^2]."""
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID:
        wp.atomic_max(out, 0, wp.abs(pressure[i, j, k] - p_ext))


@wp.kernel
def compute_descent_step_kernel(
    mechanical_kmax: wp.float64,
    max_pressure_jump: wp.array(dtype=wp.float64),
    pressure_edge_um: wp.float64,
    dt_mu: wp.array(dtype=wp.float64),
) -> None:
    """Set ``dt_mu=0.1/max(k_mech, |Delta p| ell)`` entirely on device."""
    k_pressure = max_pressure_jump[0] * pressure_edge_um
    dt_mu[0] = wp.float64(0.1) / wp.max(mechanical_kmax, k_pressure)


@wp.kernel
def update_convergence_kernel(
    max_displacement: wp.array(dtype=wp.float64),
    max_constraint_error: wp.array(dtype=wp.float64),
    max_projected_force: wp.array(dtype=wp.float64),
    dt_mu: wp.array(dtype=wp.float64),
    tolerance_um: wp.float64,
    iteration: wp.int32,
    finite: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    converged: wp.array(dtype=wp.int32),
    iters: wp.array(dtype=wp.int32),
) -> None:
    """Latch a small update only when the explicit-equivalent projected force is also small.

    ``dt_mu [um/pN] * |P F| [pN] <= tolerance [um]`` is redundant for an unscaled explicit step at
    equilibrium. It prevents a backtracked accelerator from manufacturing convergence by taking an
    arbitrarily small displacement while the mechanical residual remains large.
    """
    if active[0] != 0:
        if finite[0] == 0:
            active[0] = 0
            converged[0] = 0
        elif (
            max_displacement[0] <= tolerance_um
            and max_constraint_error[0] <= tolerance_um
            and max_projected_force[0] * dt_mu[0] <= tolerance_um
        ):
            active[0] = 0
            converged[0] = 1
            iters[0] = iteration


@wp.kernel
def record_convergence_history_kernel(
    slot: wp.int32,
    iteration: wp.int32,
    max_displacement: wp.array(dtype=wp.float64),
    max_constraint_error: wp.array(dtype=wp.float64),
    max_projected_force: wp.array(dtype=wp.float64),
    iteration_history: wp.array(dtype=wp.int32),
    displacement_history: wp.array(dtype=wp.float64),
    constraint_history: wp.array(dtype=wp.float64),
    projected_force_history: wp.array(dtype=wp.float64),
) -> None:
    """Record one device-side diagnostic checkpoint without a hot-loop host read."""
    iteration_history[slot] = iteration
    displacement_history[slot] = max_displacement[0]
    constraint_history[slot] = max_constraint_error[0]
    projected_force_history[slot] = max_projected_force[0]


@wp.kernel
def conditional_axpy_kernel(
    pos: wp.array(dtype=wp.vec3d),
    dt_mu: wp.array(dtype=wp.float64),
    force: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=wp.int32),
) -> None:
    """Projected-descent position update, disabled after device convergence."""
    i = wp.tid()
    if active[0] != 0:
        pos[i] = pos[i] + dt_mu[0] * force[i]


@wp.kernel
def conditional_rollback_vec3_kernel(
    state: wp.array(dtype=wp.vec3d),
    step_start: wp.array(dtype=wp.vec3d),
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Restore a rejected mechanical candidate from the device step-start snapshot."""
    i = wp.tid()
    if accepted[0] == 0:
        state[i] = step_start[i]


@wp.kernel
def conditional_reshape_kernel(
    pos: wp.array(dtype=wp.vec3d),
    fiber_off: wp.array(dtype=wp.int32),
    seg_off: wp.array(dtype=wp.int32),
    seg_rest: wp.array(dtype=wp.float64),
    n_iter: wp.int32,
    active: wp.array(dtype=wp.int32),
) -> None:
    """NF2007 reshape with a convergence guard (same projection arithmetic as ``ff.network_warp``)."""
    if active[0] == 0:
        return
    f = wp.tid()
    a = fiber_off[f]
    b = fiber_off[f + 1]
    p = b - a - 1
    s0 = seg_off[f]
    for _it in range(n_iter):
        for k in range(p):
            g = pos[a + k + 1] - pos[a + k]
            d = wp.length(g)
            if d > wp.float64(1.0e-300):
                u = g / d
                e = d - seg_rest[s0 + k]
                n_a = wp.float64(k + 1)
                n_b = wp.float64(p - k)
                total = n_a + n_b
                d_a = e * n_b / total
                d_b = e * n_a / total
                for m in range(a, a + k + 1):
                    pos[m] = pos[m] + d_a * u
                for m in range(a + k + 1, b):
                    pos[m] = pos[m] - d_b * u


@wp.kernel
def fire_reset_state_kernel(
    dt_mult: wp.array(dtype=wp.float64),
    alpha: wp.array(dtype=wp.float64),
    n_positive: wp.array(dtype=wp.int32),
    reset_flag: wp.array(dtype=wp.int32),
) -> None:
    """Initialise the FIRE scalar controller at the start of one inner solve.

    The fictitious velocity array is zeroed separately by the caller. ``reset_flag=1`` makes the first
    step a pure projected descent (``v`` starts at zero), matching the ``explicit`` solver's first move.
    """
    dt_mult[0] = wp.float64(1.0)
    alpha[0] = _FIRE_ALPHA0
    n_positive[0] = 0
    reset_flag[0] = 1


@wp.kernel
def fire_sumsq_kernel(vec: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.float64)) -> None:
    """Whole-cell ``sum_i |vec_i|^2`` reduction for the FIRE L2 norms (caller zeroes ``out`` first)."""
    i = wp.tid()
    v = vec[i]
    wp.atomic_add(out, 0, wp.dot(v, v))


@wp.kernel
def fire_power_kernel(
    force: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.float64),
) -> None:
    """FIRE power ``P = sum_i F_i . v_i`` reduction (caller zeroes ``out`` first)."""
    i = wp.tid()
    wp.atomic_add(out, 0, wp.dot(force[i], velocity[i]))


@wp.kernel
def fire_adapt_kernel(
    power: wp.array(dtype=wp.float64),
    dt_mult: wp.array(dtype=wp.float64),
    alpha: wp.array(dtype=wp.float64),
    n_positive: wp.array(dtype=wp.int32),
    reset_flag: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
) -> None:
    """Advance the FIRE controller from the current power sign, entirely on device.

    Uphill motion (``P <= 0``) freezes the velocity and shrinks the step; sustained downhill motion
    (more than ``N_min`` positive steps) grows the step toward ``dt_max`` and reduces the mixing. The
    step multiplier is bounded above by ``_FIRE_DT_MAX_MULT`` so the FIRE step can never outrun the
    frozen resting CFL scale ``dt_mu`` by more than the published factor.
    """
    if active[0] == 0:
        return
    if power[0] > wp.float64(0.0):
        n_positive[0] = n_positive[0] + 1
        if n_positive[0] > _FIRE_N_MIN:
            dt_mult[0] = wp.min(dt_mult[0] * _FIRE_F_INC, _FIRE_DT_MAX_MULT)
            alpha[0] = alpha[0] * _FIRE_F_ALPHA
        reset_flag[0] = 0
    else:
        n_positive[0] = 0
        dt_mult[0] = dt_mult[0] * _FIRE_F_DEC
        alpha[0] = _FIRE_ALPHA0
        reset_flag[0] = 1


@wp.kernel
def fire_mix_and_step_kernel(
    pos: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    dt_mu: wp.array(dtype=wp.float64),
    dt_mult: wp.array(dtype=wp.float64),
    alpha: wp.array(dtype=wp.float64),
    vel_sumsq: wp.array(dtype=wp.float64),
    force_sumsq: wp.array(dtype=wp.float64),
    reset_flag: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
) -> None:
    """FIRE velocity mixing plus a semi-implicit Euler step on the projected force, convergence-guarded.

    ``v <- (1-alpha) v + alpha |v| Fhat`` then ``v <- v + dt F`` and ``x <- x + dt v`` with unit
    fictitious mass and ``dt = dt_mu * dt_mult``. The global norms ``|v|`` and ``|F|`` arrive as the
    pre-mix sums of squares. On a reset (``reset_flag != 0``) the velocity is zeroed first, so the update
    degrades exactly to the projected explicit descent ``x <- x + dt F`` -- never a larger move. The
    force here is the NF2007-projected ``P F``, so FIRE stays tangent to the inextensibility constraint.
    """
    i = wp.tid()
    if active[0] == 0:
        return
    f_i = force[i]
    v_i = velocity[i]
    if reset_flag[0] != 0:
        v_i = wp.vec3d(0.0, 0.0, 0.0)
    else:
        f_norm = wp.sqrt(force_sumsq[0])
        if f_norm > wp.float64(1.0e-300):
            scale = wp.sqrt(vel_sumsq[0]) / f_norm
            a = alpha[0]
            v_i = (wp.float64(1.0) - a) * v_i + a * scale * f_i
    dt = dt_mu[0] * dt_mult[0]
    v_i = v_i + dt * f_i
    velocity[i] = v_i
    pos[i] = pos[i] + dt * v_i
