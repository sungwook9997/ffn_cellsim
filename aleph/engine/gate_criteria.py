"""The convergence gate's force floor, derived from a run's own numerics instead of frozen.

PI decision **D4**, approved 2026-07-28: replace the hardcoded ``0.21 pN`` literal with each run's own
``F_pred = inner_tolerance_um / inner_dt_mu`` and report the dimensionless ratio ``max|PF| / F_pred``.

**Why the literal had to go, stated precisely.** ``0.21`` was never a physics threshold. The executing
predicate is a DISPLACEMENT test — ``dt_mu * max|PF| <= tolerance_um``
(``ac/cell/inner_mechanics.py:296``) — with ``tolerance_um = sqrt(eps64) * ell`` (``driver.py:506``) and
``dt_mu = 0.1 / kmax`` (``assemble.py:964``). Rearranged, the force it admits is::

    F_pred = tolerance_um / dt_mu = 10 * sqrt(eps64) * ell * kmax

At the ``ell = 0.4999 um`` mesh that number is ``0.2107 pN``. Someone read it off a native run in July,
wrote "gate ~ 0.21", and it was hand-copied into dozens of sites — where it became a CONSTANT while the
quantity it came from is **proportional to the mesh spacing** ``ell``.

**So the frozen literal LOOSENS the gate on every refinement**, which is the opposite of what a reader
assumes a fixed threshold does. Scoring the 200 nm and 75 nm cortex rungs against ``0.21`` loosens them
by ``2.5x`` and ``6.6x`` respectively. Replacing it with this function is therefore a strict TIGHTENING
everywhere except the one mesh it was measured at, where it is a no-op to 0.34%.

The right fix was never to derive a better constant. It was to delete the constant and let each run
generate its own — the same principle as the D8 balance tolerance, which takes only a count.

Sanity Gate:
    * dimensional: ``tolerance_um [um] / dt_mu [um/pN]`` -> ``[pN]``. ``ell [um] * kmax [pN/um]`` -> ``[pN]``.
      The reported score is ``[pN]/[pN]`` = dimensionless, which is the point.
    * boundary: a non-finite or non-positive ``inner_dt_mu`` or ``inner_tolerance_um`` is rejected rather
      than producing an infinite or negative floor that would admit anything.
    * grid-invariance: ``F_pred`` scales linearly with ``ell`` BY CONSTRUCTION, so the dimensionless score
      is what may be compared across mesh rungs. Comparing raw ``max|PF|`` across rungs is the defect this
      module exists to remove.
    * measurement-protocol: ``max|PF|`` must be the PROJECTED residual (Dirichlet rows removed). Scoring a
      raw ``|F|`` against ``F_pred`` compares two different observables.
    * not-tuned: every input is read from the run that is being scored; nothing here can be chosen to make
      a gate pass.
"""

from __future__ import annotations

import math

import numpy as np

#: ``sqrt`` of float64 epsilon — the smallest position change resolvable at unit scale.
SQRT_EPS64: float = float(np.sqrt(np.finfo(np.float64).eps))

#: The inner predicate's step-size coefficient: ``dt_mu = DT_MU_COEFF / kmax`` (``ac/cell/assemble.py:964``).
DT_MU_COEFF: float = 0.1


#: The ``inner_tolerance_um`` recorded by the 2026-07 native resting run at ``ell = 0.4999 um``.
#:
#: Present ONLY so archived figures and probes over THAT run's data can draw a derived floor instead of
#: the hand-copied ``0.21``. It is a property of one recorded run, not a project constant: a new run must
#: use its OWN report via :func:`predicted_force_floor`. Source: ``native_resting_70686.json`` via
#: ``docs/v2_audit/GATE_A_THRESHOLD_DERIVATION_2026-07-25.md``.
NATIVE_L500_TOLERANCE_UM: float = 7.449163398500275e-09

#: The ``inner_dt_mu`` of the same recorded run. Same restriction as :data:`NATIVE_L500_TOLERANCE_UM`.
NATIVE_L500_DT_MU: float = 3.535943748517097e-08


def predicted_force_floor(*, inner_tolerance_um: float, inner_dt_mu: float) -> float:
    """Return the run's own force floor ``F_pred`` [pN] from its executing convergence predicate.

    Args:
        inner_tolerance_um: The run's displacement tolerance [um] — ``sqrt(eps64) * ell`` in the
            incumbent driver, but read from the run rather than recomputed so a driver that changes it
            is scored against what it actually used.
        inner_dt_mu: The run's inner mobility step [um/pN].

    Returns:
        ``inner_tolerance_um / inner_dt_mu`` [pN]: the largest residual force the predicate admits,
        because it is the force whose one step moves a node by exactly the tolerance.

    Raises:
        ValueError: If either input is non-finite or non-positive — both would produce a floor that
            admits an arbitrary residual, which is a silent gate removal rather than a gate.
    """
    if not math.isfinite(inner_tolerance_um) or inner_tolerance_um <= 0.0:
        raise ValueError(f"inner_tolerance_um must be finite and positive; got {inner_tolerance_um!r}")
    if not math.isfinite(inner_dt_mu) or inner_dt_mu <= 0.0:
        raise ValueError(f"inner_dt_mu must be finite and positive; got {inner_dt_mu!r}")
    return inner_tolerance_um / inner_dt_mu


def predicted_force_floor_from_mesh(*, seg_mean_um: float, kmax_pn_per_um: float) -> float:
    """Return ``F_pred`` [pN] from mesh spacing and stiffness, for callers with no run report.

    This is the same quantity as :func:`predicted_force_floor`, expanded through the incumbent driver's
    definitions ``tolerance_um = sqrt(eps64)*ell`` and ``dt_mu = 0.1/kmax``::

        F_pred = 10 * sqrt(eps64) * ell * kmax

    Prefer :func:`predicted_force_floor` whenever the run's own tolerance and step are recorded — this
    form re-derives them and so silently assumes the incumbent definitions still hold.

    Args:
        seg_mean_um: Realized mean segment length ``ell`` [um].
        kmax_pn_per_um: The stiffness scale ``kmax`` [pN/um] that set the step.

    Returns:
        The force floor [pN].

    Raises:
        ValueError: If either input is non-finite or non-positive.
    """
    if not math.isfinite(seg_mean_um) or seg_mean_um <= 0.0:
        raise ValueError(f"seg_mean_um must be finite and positive; got {seg_mean_um!r}")
    if not math.isfinite(kmax_pn_per_um) or kmax_pn_per_um <= 0.0:
        raise ValueError(f"kmax_pn_per_um must be finite and positive; got {kmax_pn_per_um!r}")
    return (SQRT_EPS64 * seg_mean_um) / (DT_MU_COEFF / kmax_pn_per_um)


def projected_residual_score(max_projected_force_pn: float, force_floor_pn: float) -> float:
    """Return the dimensionless gate score ``max|PF| / F_pred``.

    This is the number a gate should report and compare across mesh rungs. ``<= 1`` passes.

    Args:
        max_projected_force_pn: The run's ``max|PF|`` [pN] — the PROJECTED residual, Dirichlet rows
            removed. Passing a raw ``|F|`` scores a different observable and is a measurement-protocol
            error, not a stricter test.
        force_floor_pn: ``F_pred`` [pN] from :func:`predicted_force_floor`.

    Returns:
        The dimensionless ratio.

    Raises:
        ValueError: If the residual is negative or non-finite, or the floor is not positive.
    """
    if not math.isfinite(max_projected_force_pn) or max_projected_force_pn < 0.0:
        raise ValueError(
            f"max|PF| is a magnitude and must be finite; got {max_projected_force_pn!r}"
        )
    if not math.isfinite(force_floor_pn) or force_floor_pn <= 0.0:
        raise ValueError(f"force floor must be finite and positive; got {force_floor_pn!r}")
    return max_projected_force_pn / force_floor_pn


__all__ = [
    "DT_MU_COEFF",
    "SQRT_EPS64",
    "predicted_force_floor",
    "predicted_force_floor_from_mesh",
    "projected_residual_score",
]
