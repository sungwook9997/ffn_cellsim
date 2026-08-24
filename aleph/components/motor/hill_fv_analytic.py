r"""Hill-1938 force-velocity closed form — the analytic ground truth for the per-head NMII stepping.

Host-side acceptance oracle (pure NumPy, NO Warp, NO simulation runtime): the Warp-CUDA per-head Hill
stepping kernel is gated AGAINST this closed form, never the other way round (founding rule: analytic
ground truth is primary; peer engines are on-demand cross-checks only).

Hill 1938 muscle hyperbola ``(F + a)(v + b) = (F0 + a) b``. Normalizing by the stall force ``F0 = F_s`` and
writing the curvature ``kappa = a / F0`` and unloaded velocity ``v0 = b/kappa`` gives the one-line form used
throughout this track:

    v(F) = v0 * (1 - F/F_s) / (1 + (F/F_s) / kappa)                          (forward map)
    F(v) = F_s * (v0 - v) / (v0 + v / kappa)                                 (inverse map)

with the two anchors ``v(0) = v0`` (unloaded) and ``v(F_s) = 0`` (stall) holding for EVERY ``kappa > 0``.

FV-form reconciliation (⚠ params_i0b3.yaml::kappa_hill): the ``kappa -> inf`` limit is EXACTLY the linear
law ``v = v0 (1 - F/F_s)`` that PI ratified for non-muscle NMII on 2026-07-07 (Freedman 2017 AFINES;
Cytosim; Tam 2021), while ``kappa = 0.25`` is the Hill-1938 skeletal/cardiac MUSCLE hyperbola. So the linear
law is not a competing model — it is the flat-curvature limit of the SAME Hill machinery mandated by I0-A
(2026-07-16). The curvature is an I0-B3 GAP; this oracle SWEEPS it (it never chooses a value).

Sign convention: ``F`` is the resisting (opposing) load magnitude on the head, ``v`` is the shortening
(barbed-end-directed) stepping speed. Below stall the head does POSITIVE mechanical work (shortens against
load); past stall (F > F_s) NMII does not actively lengthen — ``v`` is clamped to 0 by default (matching the
retired ``minifilament_kernel``'s ``if Fmag < 0: Fmag = 0`` guard), so power stays >= 0.

Sanity Gate (of the oracle itself — self-tested in tests/ac/motor/test_hill_fv_oracle.py):
  * anchors: v(0) = v0 and v(F_s) = 0 for all kappa; F(0) = F_s, F(v0) = 0;
  * forward/inverse round-trip to machine precision on [0, F_s] x [0, v0];
  * monotone: v strictly decreasing in F on [0, F_s]; concave-down for finite kappa (Hill curvature);
  * linear limit: kappa -> inf reproduces v0 (1 - F/F_s) to tolerance; kappa = 0.25 matches Hill-1938 table;
  * force-free at v0: F(v0) = 0 (unloaded minifilament exerts ~0 contractile force);
  * power/work sign: P(F) = F v(F) >= 0 on [0, F_s], P(0) = P(F_s) = 0, single interior maximum.

References:
  Hill, A.V. 1938, Proc R Soc Lond B 126:136 (the heat/force-velocity hyperbola, a/F0 ~= 0.25 muscle).
  Freedman 2017 PLoS Comput Biol (AFINES linear stall); Tam 2021 (non-muscle FV). PI note 2026-07-07.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "hill_velocity",
    "hill_force",
    "linear_velocity",
    "power_per_head",
    "step_work",
]


def hill_velocity(
    force: npt.ArrayLike,
    v0: float,
    f_stall: float,
    kappa: float,
    *,
    clamp: bool = True,
) -> npt.NDArray[np.float64]:
    """Per-head Hill shortening velocity ``v(F)`` [µm/s] at resisting load ``force`` [pN].

    ``v = v0 (1 - F/F_s) / (1 + (F/F_s)/kappa)``, the Hill-1938 hyperbola with curvature ``kappa = a/F0``.

    Args:
        force: resisting load magnitude ``F`` [pN] (scalar or array).
        v0: unloaded stepping velocity ``v(0)`` [µm/s] (> 0).
        f_stall: per-head stall force ``F_s`` [pN] (> 0), where ``v(F_s) = 0``.
        kappa: Hill curvature ``a/F0`` (> 0); ``kappa -> inf`` is the linear law.
        clamp: if True (default) floor ``v`` at 0 past stall (no active lengthening); False keeps the
            raw hyperbola (v < 0 for F > F_s), useful for inspecting the mathematical form.

    Returns:
        Shortening velocity ``v`` [µm/s], same shape as ``force``.
    """
    if v0 <= 0.0:
        raise ValueError("v0 must be positive")
    if f_stall <= 0.0:
        raise ValueError("f_stall must be positive")
    if kappa <= 0.0:
        raise ValueError("kappa must be positive (kappa -> inf is the linear limit)")
    f = np.asarray(force, dtype=np.float64)
    fr = f / f_stall
    v = v0 * (1.0 - fr) / (1.0 + fr / kappa)
    if clamp:
        v = np.maximum(v, 0.0)
    return v


def hill_force(
    velocity: npt.ArrayLike,
    v0: float,
    f_stall: float,
    kappa: float,
) -> npt.NDArray[np.float64]:
    """Inverse Hill map: resisting load ``F(v)`` [pN] that produces shortening velocity ``velocity`` [µm/s].

    ``F = F_s (v0 - v) / (v0 + v/kappa)`` — the exact inverse of :func:`hill_velocity` on ``v in [0, v0]``.
    At ``v = 0`` (isometric) this returns ``F_s``; at ``v = v0`` (unloaded) it returns 0 (force-free).

    Args:
        velocity: shortening velocity ``v`` [µm/s] (scalar or array), physical range ``[0, v0]``.
        v0: unloaded stepping velocity [µm/s] (> 0).
        f_stall: per-head stall force [pN] (> 0).
        kappa: Hill curvature ``a/F0`` (> 0).

    Returns:
        Resisting load ``F`` [pN], same shape as ``velocity``.
    """
    if v0 <= 0.0:
        raise ValueError("v0 must be positive")
    if f_stall <= 0.0:
        raise ValueError("f_stall must be positive")
    if kappa <= 0.0:
        raise ValueError("kappa must be positive")
    v = np.asarray(velocity, dtype=np.float64)
    return f_stall * (v0 - v) / (v0 + v / kappa)


def linear_velocity(
    force: npt.ArrayLike,
    v0: float,
    f_stall: float,
    *,
    clamp: bool = True,
) -> npt.NDArray[np.float64]:
    """Linear (affine) force-velocity ``v = v0 (1 - F/F_s)`` [µm/s] — the ``kappa -> inf`` Hill limit.

    This is the PI-2026-07-07 non-muscle NMII law (Freedman 2017 / Cytosim / Tam 2021) and the isolated
    diagnostic form of the retired ``ff/myosin_linear.minifilament_kernel``. Provided so the ``kappa -> inf``
    limit of :func:`hill_velocity` can be gated against it directly.

    Args:
        force: resisting load ``F`` [pN].
        v0: unloaded velocity [µm/s].
        f_stall: stall force [pN].
        clamp: floor at 0 past stall (default True).

    Returns:
        Shortening velocity [µm/s], same shape as ``force``.
    """
    if v0 <= 0.0 or f_stall <= 0.0:
        raise ValueError("v0 and f_stall must be positive")
    f = np.asarray(force, dtype=np.float64)
    v = v0 * (1.0 - f / f_stall)
    if clamp:
        v = np.maximum(v, 0.0)
    return v


def power_per_head(
    force: npt.ArrayLike,
    v0: float,
    f_stall: float,
    kappa: float,
) -> npt.NDArray[np.float64]:
    """Instantaneous mechanical power ``P(F) = F * v(F)`` [pN·µm/s] delivered by one head at load ``force``.

    Non-negative on ``[0, F_s]`` with ``P(0) = P(F_s) = 0`` and a single interior maximum (the muscle
    max-power point) — the mechanical-work-sign gate. Uses the clamped Hill velocity.

    Args:
        force: resisting load ``F`` [pN].
        v0: unloaded velocity [µm/s].
        f_stall: stall force [pN].
        kappa: Hill curvature.

    Returns:
        Power ``P`` [pN·µm/s], same shape as ``force``.
    """
    f = np.asarray(force, dtype=np.float64)
    return f * hill_velocity(f, v0, f_stall, kappa, clamp=True)


def step_work(force: npt.ArrayLike, d_step: float) -> npt.NDArray[np.float64]:
    """Mechanical work per power stroke ``W = F * d_step`` [pN·µm] against resisting load ``force``.

    The ATP/step work-sign gate checks ``0 <= W < dG_ATP`` on ``[0, F_s]`` (positive work below stall;
    2nd-law efficiency ``eta = W/dG_ATP < 1``). ``d_step`` is the ~5-10 nm working stroke (I0-B3 GAP).

    Args:
        force: resisting load ``F`` [pN].
        d_step: per-ATP working-stroke displacement [µm].

    Returns:
        Work per step ``W`` [pN·µm], same shape as ``force``.
    """
    if d_step <= 0.0:
        raise ValueError("d_step must be positive")
    return np.asarray(force, dtype=np.float64) * d_step
