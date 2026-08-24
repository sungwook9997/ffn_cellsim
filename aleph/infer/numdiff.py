"""Finite differences, kept in one place because they are the arbiter for every gradient here.

Three modules need a numerical derivative: the likelihoods need one to verify their analytic
gradients, the Fisher matrix needs a forward-model Jacobian, and the Laplace approximation needs a
Hessian.  Writing it three times would mean three step-size conventions and three sets of accuracy
claims, so it is written once.

Step size follows the standard optimum for the truncation/roundoff trade-off.  Central differences
have truncation error O(h^2 f''') and roundoff O(eps*|f|/h), which balance at h ~ eps^(1/3); forward
differences balance at h ~ eps^(1/2).  The step is scaled by ``max(|theta_i|, 1)`` so a parameter of
size 1e-9 does not get a step larger than itself and a parameter of size 1e9 does not get a step
that vanishes into its own rounding.

Central differences are the default and the only method used for accuracy claims.  Forward
differences exist because they cost half as much, and are marked as the lower-accuracy option.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Callable

__all__ = [
    "CENTRAL_STEP",
    "FORWARD_STEP",
    "gradient",
    "hessian",
    "jacobian",
    "max_relative_gradient_error",
]

#: Optimal relative step for central differences: eps^(1/3) ~ 6.06e-6 in float64.
CENTRAL_STEP: float = float(np.finfo(float).eps ** (1.0 / 3.0))
#: Optimal relative step for forward differences: eps^(1/2) ~ 1.49e-8 in float64.
FORWARD_STEP: float = float(np.finfo(float).eps ** 0.5)


def _steps(theta: NDArray, rel_step: float) -> NDArray:
    return rel_step * np.maximum(np.abs(theta), 1.0)


def jacobian(
    model: Callable[[NDArray], NDArray],
    theta: NDArray,
    *,
    rel_step: float | None = None,
    method: str = "central",
) -> NDArray:
    """Jacobian ``d model_i / d theta_j`` of a vector-valued forward model.

    Args:
        model: Parameters -> predicted observable.  May return a scalar or a 1-D array.
        theta: Point at which to differentiate.
        rel_step: Relative step.  Defaults to :data:`CENTRAL_STEP` or :data:`FORWARD_STEP`.
        method: ``"central"`` (default, second order) or ``"forward"`` (first order, half the cost).

    Returns:
        Array of shape ``(n_observations, n_parameters)``.
    """
    theta = np.asarray(theta, dtype=float)
    if theta.ndim != 1:
        raise ValueError(f"theta must be 1-D, got shape {theta.shape}")
    if method not in ("central", "forward"):
        raise ValueError(f"method must be 'central' or 'forward', got {method!r}")
    if rel_step is None:
        rel_step = CENTRAL_STEP if method == "central" else FORWARD_STEP

    base = np.atleast_1d(np.asarray(model(theta), dtype=float))
    if base.ndim != 1:
        raise ValueError(f"model must return a scalar or 1-D array, got shape {base.shape}")

    h = _steps(theta, rel_step)
    out = np.empty((base.size, theta.size), dtype=float)
    for j in range(theta.size):
        up = theta.copy()
        up[j] += h[j]
        if method == "central":
            down = theta.copy()
            down[j] -= h[j]
            # Recover the exactly representable step actually taken; theta[j] + h[j] is generally
            # not theta[j] + h[j] in floating point, and using the nominal h loses a digit or two.
            step = up[j] - down[j]
            out[:, j] = (
                np.atleast_1d(np.asarray(model(up), dtype=float))
                - np.atleast_1d(np.asarray(model(down), dtype=float))
            ) / step
        else:
            step = up[j] - theta[j]
            out[:, j] = (np.atleast_1d(np.asarray(model(up), dtype=float)) - base) / step
    return out


def gradient(
    scalar_fn: Callable[[NDArray], float],
    theta: NDArray,
    *,
    rel_step: float | None = None,
    method: str = "central",
) -> NDArray:
    """Gradient of a scalar function, as a 1-D array."""
    return jacobian(
        lambda t: np.array([float(scalar_fn(t))]), theta, rel_step=rel_step, method=method
    )[0]


def hessian(
    grad_fn: Callable[[NDArray], NDArray],
    theta: NDArray,
    *,
    rel_step: float | None = None,
) -> NDArray:
    """Hessian by central-differencing an analytic gradient, symmetrized.

    Differencing the gradient rather than the value costs one derivative order of accuracy less and
    keeps the result close to symmetric to begin with.  The explicit symmetrization removes the
    asymmetry that is left, which is pure numerical noise and would otherwise put a small imaginary
    part into eigenvalues downstream.
    """
    raw = jacobian(grad_fn, np.asarray(theta, dtype=float), rel_step=rel_step, method="central")
    return 0.5 * (raw + raw.T)


def max_relative_gradient_error(
    value_fn: Callable[[NDArray], float],
    grad_fn: Callable[[NDArray], NDArray],
    theta: NDArray,
    *,
    rel_step: float | None = None,
) -> float:
    """Largest disagreement between an analytic gradient and central differences.

    Scaled by the size of the gradient itself so the number means the same thing for a likelihood
    of magnitude 1 and one of magnitude 1e6.  This is the arbiter: an analytic gradient that
    disagrees with this is wrong, and the sign convention it implies is the one to keep.
    """
    theta = np.asarray(theta, dtype=float)
    analytic = np.asarray(grad_fn(theta), dtype=float)
    numeric = gradient(value_fn, theta, rel_step=rel_step)
    scale = max(float(np.max(np.abs(analytic))), float(np.max(np.abs(numeric))), 1e-300)
    return float(np.max(np.abs(analytic - numeric)) / scale)
