"""Device-resident nonlinear acceleration for projected Active Cell mechanics.

The accelerators in this module change only the numerical path used to propose a displacement. Every proposal
is evaluated by the driver's exact nonlinear force assembly and competes with both the unchanged state and the
matched explicit-CFL update. The installed projected residual therefore cannot increase. Irreversible KMC
remains outside the inner solve.

Sanity Gate:
    * Dimensions: block directions and Anderson mappings are positions/displacements in um; every RKC force
      increment is ``dt_mu [um/pN] * projected_force [pN]``.
    * Boundary cases: the first Anderson call and a rank-deficient secant both reduce exactly to the supplied
      base direction; inactive/invalid solves cannot create history.
    * Conservation: the base direction is projected before entry; Anderson forms affine combinations of
      projected mappings and is followed by the same NF2007 reshape/check as every other candidate.
    * Numerical: depth-one Anderson has no history-size or damping knob. The secant coefficient is obtained by
      a device FP64 least-squares reduction; a zero denominator selects the base direction. RKC uses the exact
      first-order Chebyshev recurrence for its declared integer stage budget and no fitted damping parameter.
    * Residency: history, reductions, coefficient, validity, and proposal remain Warp-CUDA resident.
"""

from __future__ import annotations

import warp as wp

__all__ = [
    "AndersonDepthOne",
    "rkc1_recurrence_kernel",
    "displacement_between_kernel",
]


@wp.kernel
def _mapping_kernel(
    pos: wp.array(dtype=wp.vec3d),
    base_direction: wp.array(dtype=wp.vec3d),
    mapping: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    mapping[i] = pos[i] + base_direction[i]


@wp.kernel
def _difference_kernel(
    current_residual: wp.array(dtype=wp.vec3d),
    previous_residual: wp.array(dtype=wp.vec3d),
    current_mapping: wp.array(dtype=wp.vec3d),
    previous_mapping: wp.array(dtype=wp.vec3d),
    residual_difference: wp.array(dtype=wp.vec3d),
    mapping_difference: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    residual_difference[i] = current_residual[i] - previous_residual[i]
    mapping_difference[i] = current_mapping[i] - previous_mapping[i]


@wp.kernel
def _dot_kernel(
    left: wp.array(dtype=wp.vec3d),
    right: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.float64),
) -> None:
    i = wp.tid()
    wp.atomic_add(out, 0, wp.dot(left[i], right[i]))


@wp.kernel
def _coefficient_kernel(
    numerator: wp.array(dtype=wp.float64),
    denominator: wp.array(dtype=wp.float64),
    history_valid: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    coefficient: wp.array(dtype=wp.float64),
) -> None:
    """Solve the depth-one residual least-squares coefficient or select the base map."""
    coefficient[0] = wp.float64(0.0)
    if history_valid[0] == 0 or finite[0] == 0:
        return
    den = denominator[0]
    num = numerator[0]
    if den > wp.float64(1.0e-300) and wp.isfinite(den) and wp.isfinite(num):
        coefficient[0] = num / den


@wp.kernel
def _anderson_direction_kernel(
    pos: wp.array(dtype=wp.vec3d),
    base_direction: wp.array(dtype=wp.vec3d),
    current_mapping: wp.array(dtype=wp.vec3d),
    mapping_difference: wp.array(dtype=wp.vec3d),
    coefficient: wp.array(dtype=wp.float64),
    history_valid: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    direction: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    if active[0] != 0 and finite[0] != 0 and history_valid[0] != 0:
        direction[i] = current_mapping[i] - coefficient[0] * mapping_difference[i] - pos[i]
    else:
        direction[i] = base_direction[i]


@wp.kernel
def _record_history_kernel(
    current_residual: wp.array(dtype=wp.vec3d),
    current_mapping: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    previous_residual: wp.array(dtype=wp.vec3d),
    previous_mapping: wp.array(dtype=wp.vec3d),
) -> None:
    i = wp.tid()
    if active[0] != 0 and finite[0] != 0:
        previous_residual[i] = current_residual[i]
        previous_mapping[i] = current_mapping[i]


@wp.kernel
def _mark_history_kernel(
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    history_valid: wp.array(dtype=wp.int32),
) -> None:
    if active[0] != 0 and finite[0] != 0:
        history_valid[0] = wp.int32(1)


class AndersonDepthOne:
    r"""Depth-one Type-II Anderson acceleration with FP64 device reductions.

    For base fixed-point residual ``f_k = G(x_k)-x_k`` this computes

    ``gamma = <Delta f, f_k> / <Delta f, Delta f>`` and
    ``x_AA = G(x_k) - gamma * Delta G``.

    The returned displacement is ``x_AA-x_k``. A singular secant or first call returns ``f_k`` exactly.
    """

    def __init__(self, n: int, device: str) -> None:
        self.n = int(n)
        self.device = device
        with wp.ScopedDevice(device):
            self.previous_residual = wp.zeros(self.n, dtype=wp.vec3d)
            self.previous_mapping = wp.zeros(self.n, dtype=wp.vec3d)
            self.current_mapping = wp.zeros(self.n, dtype=wp.vec3d)
            self.residual_difference = wp.zeros(self.n, dtype=wp.vec3d)
            self.mapping_difference = wp.zeros(self.n, dtype=wp.vec3d)
            self.direction_d = wp.zeros(self.n, dtype=wp.vec3d)
            self.numerator = wp.zeros(1, dtype=wp.float64)
            self.denominator = wp.zeros(1, dtype=wp.float64)
            self.coefficient = wp.zeros(1, dtype=wp.float64)
            self.history_valid = wp.zeros(1, dtype=wp.int32)

    def reset(self) -> None:
        """Invalidate numerical history at the start of each outer transaction."""
        self.history_valid.zero_()
        self.coefficient.zero_()

    def direction(
        self,
        pos: wp.array,
        base_direction: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the accelerated displacement and record the current base mapping on device."""
        d = self.device
        wp.launch(_mapping_kernel, dim=self.n,
                  inputs=[pos, base_direction, self.current_mapping], device=d)
        wp.launch(_difference_kernel, dim=self.n,
                  inputs=[base_direction, self.previous_residual, self.current_mapping,
                          self.previous_mapping, self.residual_difference,
                          self.mapping_difference], device=d)
        self.numerator.zero_()
        self.denominator.zero_()
        wp.launch(_dot_kernel, dim=self.n,
                  inputs=[self.residual_difference, base_direction, self.numerator], device=d)
        wp.launch(_dot_kernel, dim=self.n,
                  inputs=[self.residual_difference, self.residual_difference, self.denominator], device=d)
        wp.launch(_coefficient_kernel, dim=1,
                  inputs=[self.numerator, self.denominator, self.history_valid, finite,
                          self.coefficient], device=d)
        wp.launch(_anderson_direction_kernel, dim=self.n,
                  inputs=[pos, base_direction, self.current_mapping, self.mapping_difference,
                          self.coefficient, self.history_valid, active, finite,
                          self.direction_d], device=d)
        wp.launch(_record_history_kernel, dim=self.n,
                  inputs=[base_direction, self.current_mapping, active, finite,
                          self.previous_residual, self.previous_mapping], device=d)
        wp.launch(_mark_history_kernel, dim=1,
                  inputs=[active, finite, self.history_valid], device=d)
        return self.direction_d


@wp.kernel
def rkc1_recurrence_kernel(
    previous: wp.array(dtype=wp.vec3d),
    current: wp.array(dtype=wp.vec3d),
    projected_force: wp.array(dtype=wp.vec3d),
    dt_mu: wp.array(dtype=wp.float64),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    next_state: wp.array(dtype=wp.vec3d),
) -> None:
    r"""One first-order Runge--Kutta--Chebyshev recurrence stage.

    With ``Y0=x`` and ``Y1=x+dt_mu F(x)``, later stages use
    ``Yj = 2 Y(j-1) - Y(j-2) + 2 dt_mu F(Y(j-1))``. For a frozen linear mode this evaluates
    ``T_s(1 + dt_mu lambda)`` and therefore extends one validated Euler step into an ``s^2``-slope
    super-step without a sparse solve. The driver evaluates the full nonlinear result and may reject it.
    """
    i = wp.tid()
    if active[0] != 0 and finite[0] != 0:
        next_state[i] = (
            wp.float64(2.0) * current[i]
            - previous[i]
            + wp.float64(2.0) * dt_mu[0] * projected_force[i]
        )
    else:
        next_state[i] = current[i]


@wp.kernel
def displacement_between_kernel(
    start: wp.array(dtype=wp.vec3d),
    end: wp.array(dtype=wp.vec3d),
    displacement: wp.array(dtype=wp.vec3d),
) -> None:
    """Return ``end-start`` for a candidate generated in a device position workspace."""
    i = wp.tid()
    displacement[i] = end[i] - start[i]
