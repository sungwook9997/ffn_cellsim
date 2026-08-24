"""Closed-form acceptance oracle: linear viscoelastic relaxation & creep.

The S1 (viscoelastic sphere) / M2 gate. Provides the relaxation modulus ``G(t)``
(step-strain) and creep compliance ``J(t)`` (step-stress) for the three canonical
linear viscoelastic models — Maxwell, Kelvin-Voigt, and the Standard Linear Solid
(Zener) — plus a robust relaxation-time recovery from a measured ``G(t)`` curve.

VALIDATION oracle (runtime-import-forbidden). Pure functional forms, SI units.

Provenance / discipline
-----------------------
Standard linear viscoelasticity (Findley, Lai & Onaran, *Creep and Relaxation of
Nonlinear Viscoelastic Materials*; Christensen, *Theory of Viscoelasticity*). The FF
cortex is **emergently** viscoelastic: crosslink turnover relaxes stress with an
SLS-like ``G(t) = G_inf + G1 e^{-t/tau}`` (``docs/v2_audit/FF_RESULTS_LOG.md``;
``ff/ecm_mechanics.py:stress_relaxation`` MEASURES tau from a run). This module is
the analytic ground truth the emergent relaxation is gated against — the oracle is
a cross-check, never the runtime mechanism.

Model summary (spring modulus E [Pa], dashpot viscosity eta [Pa·s], tau = eta/E [s])
------------------------------------------------------------------------------------
- Maxwell (spring + dashpot in SERIES): a viscoelastic FLUID.
  G(t) = G0 e^{-t/tau};  J(t) = (1/G0)(1 + t/tau)  (unbounded creep).
- Kelvin-Voigt (spring + dashpot in PARALLEL): a viscoelastic SOLID, no instant strain.
  J(t) = (1/E)(1 - e^{-t/tau});  relaxation is a delta at t=0 (not a simple exponential).
- Standard Linear Solid / Zener (spring E1 ∥ [spring E2 + dashpot in series]): solid with
  both instantaneous and equilibrium moduli — the FF-relevant model.
  G(t) = E1 + E2 e^{-t/tau_sigma}, tau_sigma = eta/E2 (relaxation time);
  J(t) = 1/E1 - (E2/(E1(E1+E2))) e^{-t/tau_eps}, tau_eps = tau_sigma (E1+E2)/E1 (retardation time).
  G(0)=E1+E2 (glassy), G(inf)=E1 (equilibrium); J(0)=1/(E1+E2), J(inf)=1/E1; tau_eps > tau_sigma.

Sanity Gate
-----------
- Dimensional: G [Pa], J [1/Pa], t & tau [s]; G(t)·J(t) → 1 only at t=0 and t→inf (not for all t).
- Boundary: G monotone non-increasing; J monotone non-decreasing; both → their equilibrium value.
- Consistency: J(0) = 1/G(0) and J(inf) = 1/G(inf) for the SLS (checked in tests).
- Retardation > relaxation: tau_eps > tau_sigma for a genuine SLS (E2 > 0).
- Recovery: ``relaxation_time_from_curve`` recovers (G_inf, G1, tau) from a noise-free G(t).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "maxwell_relaxation",
    "maxwell_creep_compliance",
    "kelvin_voigt_creep",
    "sls_relaxation",
    "sls_creep_compliance",
    "sls_retardation_time",
    "relaxation_time_from_curve",
]

ArrayLike = npt.ArrayLike


def _check_time(t: ArrayLike) -> npt.NDArray[np.float64]:
    arr = np.asarray(t, dtype=np.float64)
    if np.any(arr < 0.0):
        raise ValueError("time t must be non-negative.")
    return arr


# --------------------------------------------------------------------------- #
# Maxwell (series) — viscoelastic fluid
# --------------------------------------------------------------------------- #
def maxwell_relaxation(t: ArrayLike, G0: float, tau: float) -> npt.NDArray[np.float64]:
    """Maxwell relaxation modulus ``G(t) = G0 e^{-t/tau}`` (step strain).

    Args:
        t: Time(s) (s, >= 0).
        G0: Instantaneous modulus (Pa, > 0).
        tau: Relaxation time ``eta/G0`` (s, > 0).

    Returns:
        Relaxation modulus ``G(t)`` (Pa), elementwise.

    Raises:
        ValueError: if ``G0 <= 0``, ``tau <= 0``, or any ``t < 0``.
    """
    if G0 <= 0.0 or tau <= 0.0:
        raise ValueError("G0 and tau must be strictly positive.")
    return G0 * np.exp(-_check_time(t) / tau)


def maxwell_creep_compliance(
    t: ArrayLike, G0: float, tau: float
) -> npt.NDArray[np.float64]:
    """Maxwell creep compliance ``J(t) = (1/G0)(1 + t/tau)`` (step stress).

    Unbounded (fluid) creep: the dashpot flows without limit.

    Args:
        t: Time(s) (s, >= 0).
        G0: Instantaneous modulus (Pa, > 0).
        tau: Relaxation time (s, > 0); ``eta = G0 * tau``.

    Returns:
        Creep compliance ``J(t)`` (1/Pa), elementwise.

    Raises:
        ValueError: if ``G0 <= 0``, ``tau <= 0``, or any ``t < 0``.
    """
    if G0 <= 0.0 or tau <= 0.0:
        raise ValueError("G0 and tau must be strictly positive.")
    return (1.0 / G0) * (1.0 + _check_time(t) / tau)


# --------------------------------------------------------------------------- #
# Kelvin-Voigt (parallel) — viscoelastic solid, delayed elasticity
# --------------------------------------------------------------------------- #
def kelvin_voigt_creep(
    t: ArrayLike, E: float, tau: float
) -> npt.NDArray[np.float64]:
    """Kelvin-Voigt creep compliance ``J(t) = (1/E)(1 - e^{-t/tau})`` (step stress).

    No instantaneous strain (``J(0) = 0``); bounded equilibrium ``J(inf) = 1/E``.

    Args:
        t: Time(s) (s, >= 0).
        E: Spring modulus (Pa, > 0).
        tau: Retardation time ``eta/E`` (s, > 0).

    Returns:
        Creep compliance ``J(t)`` (1/Pa), elementwise.

    Raises:
        ValueError: if ``E <= 0``, ``tau <= 0``, or any ``t < 0``.
    """
    if E <= 0.0 or tau <= 0.0:
        raise ValueError("E and tau must be strictly positive.")
    return (1.0 / E) * (1.0 - np.exp(-_check_time(t) / tau))


# --------------------------------------------------------------------------- #
# Standard Linear Solid (Zener) — the FF-relevant model
# --------------------------------------------------------------------------- #
def sls_relaxation(
    t: ArrayLike, G_inf: float, G1: float, tau: float
) -> npt.NDArray[np.float64]:
    """SLS relaxation modulus ``G(t) = G_inf + G1 e^{-t/tau}`` (step strain).

    ``G_inf`` is the equilibrium (relaxed) modulus ``E1``; ``G1`` is the relaxing
    part ``E2``; ``tau`` is the relaxation time ``tau_sigma = eta/E2``. Glassy modulus
    ``G(0) = G_inf + G1``.

    Args:
        t: Time(s) (s, >= 0).
        G_inf: Equilibrium modulus (Pa, > 0).
        G1: Relaxing modulus (Pa, >= 0).
        tau: Relaxation time (s, > 0).

    Returns:
        Relaxation modulus ``G(t)`` (Pa), elementwise.

    Raises:
        ValueError: if ``G_inf <= 0``, ``G1 < 0``, ``tau <= 0``, or any ``t < 0``.
    """
    if G_inf <= 0.0 or G1 < 0.0 or tau <= 0.0:
        raise ValueError("Require G_inf > 0, G1 >= 0, tau > 0.")
    return G_inf + G1 * np.exp(-_check_time(t) / tau)


def sls_retardation_time(G_inf: float, G1: float, tau: float) -> float:
    """SLS retardation (creep) time ``tau_eps = tau_sigma (E1+E2)/E1``.

    Always ``>= tau`` (the relaxation time); equal only if ``G1 = 0`` (purely elastic).

    Args:
        G_inf: Equilibrium modulus ``E1`` (Pa, > 0).
        G1: Relaxing modulus ``E2`` (Pa, >= 0).
        tau: Relaxation time ``tau_sigma`` (s, > 0).

    Returns:
        Retardation time ``tau_eps`` (s).

    Raises:
        ValueError: if ``G_inf <= 0``, ``G1 < 0``, or ``tau <= 0``.
    """
    if G_inf <= 0.0 or G1 < 0.0 or tau <= 0.0:
        raise ValueError("Require G_inf > 0, G1 >= 0, tau > 0.")
    return float(tau * (G_inf + G1) / G_inf)


def sls_creep_compliance(
    t: ArrayLike, G_inf: float, G1: float, tau: float
) -> npt.NDArray[np.float64]:
    """SLS creep compliance (step stress).

    ``J(t) = 1/E1 - (E2/(E1(E1+E2))) e^{-t/tau_eps}`` with ``E1 = G_inf``, ``E2 = G1``
    and retardation time ``tau_eps`` (see :func:`sls_retardation_time`). Satisfies
    ``J(0) = 1/(G_inf+G1)`` (glassy) and ``J(inf) = 1/G_inf`` (equilibrium).

    Args:
        t: Time(s) (s, >= 0).
        G_inf: Equilibrium modulus ``E1`` (Pa, > 0).
        G1: Relaxing modulus ``E2`` (Pa, >= 0).
        tau: Relaxation time ``tau_sigma`` (s, > 0).

    Returns:
        Creep compliance ``J(t)`` (1/Pa), elementwise.

    Raises:
        ValueError: if ``G_inf <= 0``, ``G1 < 0``, ``tau <= 0``, or any ``t < 0``.
    """
    if G_inf <= 0.0 or G1 < 0.0 or tau <= 0.0:
        raise ValueError("Require G_inf > 0, G1 >= 0, tau > 0.")
    tau_eps = sls_retardation_time(G_inf, G1, tau)
    amp = G1 / (G_inf * (G_inf + G1))
    return 1.0 / G_inf - amp * np.exp(-_check_time(t) / tau_eps)


# --------------------------------------------------------------------------- #
# Recovery: fit (G_inf, G1, tau) from a measured relaxation curve
# --------------------------------------------------------------------------- #
def _sls_separable_fit(
    tt: npt.NDArray[np.float64], g: npt.NDArray[np.float64], taus: npt.NDArray[np.float64]
) -> tuple[float, float, float, float]:
    """Best ``(ss_res, tau, G_inf, G1)`` over a grid of candidate relaxation times.

    For a fixed ``tau`` the model ``G = G_inf + G1 e^{-t/tau}`` is LINEAR in
    ``(G_inf, G1)``, so each candidate is a 2-parameter least squares (separable
    nonlinear LS); the grid picks the ``tau`` with the smallest residual.
    """
    best = (float("inf"), 0.0, 0.0, 0.0)
    for tau in taus:
        X = np.column_stack([np.ones_like(tt), np.exp(-tt / tau)])
        coef, _, _, _ = np.linalg.lstsq(X, g, rcond=None)
        ss = float(np.sum((g - X @ coef) ** 2))
        if ss < best[0]:
            best = (ss, float(tau), float(coef[0]), float(coef[1]))
    return best


def relaxation_time_from_curve(t: ArrayLike, G: ArrayLike) -> dict[str, object]:
    """Recover ``(G_inf, G1, tau)`` from a measured SLS relaxation curve ``G(t)``.

    Uses separable nonlinear least squares: since ``G(t) = G_inf + G1 e^{-t/tau}`` is
    linear in ``(G_inf, G1)`` at fixed ``tau``, a coarse log-spaced ``tau`` grid (each
    a 2-parameter linear solve) is minimised, then refined locally. This avoids the
    tail bias of a plateau-subtract + log-linear fit. Matches the SLS form the FF
    crosslink-turnover relaxation produces (``ff/ecm_mechanics.py:stress_relaxation``).

    Args:
        t: Sample times (s, >= 0), strictly increasing, length n >= 3.
        G: Measured relaxation modulus (Pa), length n, decreasing toward a plateau.

    Returns:
        Dict with ``G_inf`` (Pa), ``G1`` (Pa), ``tau`` (s), ``G0`` (= G_inf + G1),
        ``tau_eps`` (retardation time), and ``r_squared`` of the fit (in G-space).

    Raises:
        ValueError: on length mismatch, n < 3, or a non-decaying curve (fitted G1 <= 0).
    """
    tt = np.asarray(t, dtype=np.float64)
    g = np.asarray(G, dtype=np.float64)
    if tt.shape != g.shape:
        raise ValueError("t and G must have the same shape.")
    if tt.ndim != 1 or tt.size < 3:
        raise ValueError("Need at least 3 samples.")

    span = float(tt.max() - tt.min())
    if span <= 0.0:
        raise ValueError("t must span a positive interval.")

    coarse = np.geomspace(span / 500.0, span * 5.0, 400)
    _, tau_c, _, _ = _sls_separable_fit(tt, g, coarse)
    fine = np.geomspace(tau_c / 2.0, tau_c * 2.0, 400)
    ss_res, tau, g_inf, g1 = _sls_separable_fit(tt, g, fine)

    if g1 <= 0.0 or g_inf <= 0.0:
        raise ValueError("Curve is not a decaying relaxation (fitted G1 or G_inf <= 0).")

    ss_tot = float(np.sum((g - g.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0

    return {
        "G_inf": g_inf,
        "G1": g1,
        "tau": tau,
        "G0": g_inf + g1,
        "tau_eps": sls_retardation_time(g_inf, g1, tau),
        "r_squared": r_squared,
    }
