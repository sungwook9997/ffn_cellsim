r"""Arp2/3 branch junction — 3-body ANGLE-HARMONIC + thermal fluctuation (host NumPy oracle) — I4.

The mechanistic Arp2/3 branch (CLAUDE.md HARD worked-example: "Angle-harmonic with thermal fluctuation",
NOT a rigid 72° constraint). One branch junction is the mother-node → branch-node → daughter-node triple; the
angle ``theta`` at the branch vertex between the (branch→mother) and (branch→daughter) arms carries a harmonic
potential

    U(theta) = 1/2 * k_theta * (theta - theta0)^2 ,   theta0 ~ 70 deg (Faessler 2020 in-cell cryo-ET 68 +- 9),

so the angle FLUCTUATES and is not clamped. In thermal equilibrium the small-fluctuation (harmonic)
width is fixed by equipartition,

    <(theta - theta0)^2> = kT / k_theta   =>   sigma_theta = sqrt(kT / k_theta) ,

which is exactly how ``ff.architecture_spec.ARP23_BRANCH_K = 0.173 pN*um/rad^2`` was set from the observed
sigma_theta = 9 deg at T = 310 K (kT = 4.28e-3 pN*um): k_theta = kT / sigma_theta^2. This module is the
host-NumPy acceptance oracle for those two claims; the device realization is
``ac.weave.branch_angle_warp.branch_angle_kernel`` (a port of ``ff.network_warp.branch_angle_kernel``, whose
standard angle-force expression this oracle matches bit-for-formula so the native gate reproduces the CPU
gate).

engine units: length um, force pN, energy pN*um, angle rad; kT_310K = 4.28e-3 pN*um.

Sanity Gate (self-tested in tests/ac/weave/test_branch_angle_oracle.py):
  * dimensional: U in pN*um; the analytic torque -dU/dtheta = -k_theta*(theta-theta0); minimum at theta0.
  * FD-gradient sign arbiter: the analytic nodal force matches a finite-difference of U to tolerance, and it
    RESTORES toward theta0 (a bent junction is pulled back, an over-opened one pushed in).
  * equipartition round-trip: sigma_theta = sqrt(kT/k_theta) and k_theta = kT/sigma_theta^2 invert exactly;
    the ARP23 anchor (0.173 pN*um/rad^2, sigma = 9 deg, 310 K) round-trips.
  * thermal distribution: a Boltzmann ensemble exp(-U/kT) has SD ~ sigma_theta (measure-free harmonic limit
    exact; the physical solid-angle sin(theta) measure shifts the peak slightly, reported not tuned).
  * Newton's 3rd law: the three nodal forces sum to zero (an internal potential exerts no net force).

Reference: Faessler et al. 2020 EMBO J 39:e104254 (in-cell branch angle 68 +- 9 deg); Mueller 2017; CLAUDE.md
worked-example table; matches ff.network_warp.branch_angle_kernel (standard harmonic-angle force).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "KT_310K_PN_UM",
    "ARP23_THETA0_RAD",
    "ARP23_SIGMA_DEG",
    "ARP23_K_THETA",
    "branch_angle",
    "angle_energy",
    "angle_forces",
    "thermal_sigma",
    "k_theta_from_sigma",
    "sample_branch_angles",
]

# Thermal energy at the physiological setpoint (310 K). kT = k_B*T = 1.380649e-23 J/K * 310 K = 4.28e-21 J;
# 1 pN*um = 1e-18 J, so kT = 4.28e-3 pN*um. Sourced physical constant, not a tuning knob.
KT_310K_PN_UM: float = 1.380649e-23 * 310.0 / 1.0e-18

# Arp2/3 branch anchors (mirror ff.architecture_spec; lit-anchored, Faessler 2020). Recorded here so the
# oracle can round-trip the equipartition relation the ff constant was derived from — NOT re-chosen.
ARP23_THETA0_RAD: float = 1.2217304764        # 70.0 deg rest branch angle theta0
ARP23_SIGMA_DEG: float = 9.0                  # thermal angular spread sigma_theta (Faessler 2020)
ARP23_K_THETA: float = 0.173                  # pN*um/rad^2 = kT / sigma_theta^2 (equipartition; see round-trip)


def branch_angle(
    p_mother: npt.NDArray[np.float64],
    p_branch: npt.NDArray[np.float64],
    p_daughter: npt.NDArray[np.float64],
) -> float:
    """Angle ``theta`` [rad] at the branch vertex between the (branch->mother) and (branch->daughter) arms.

    Args:
        p_mother: (3,) mother-side node position [um] (the ``i`` of the ff triple).
        p_branch: (3,) branch-vertex node position [um] (the center ``j``).
        p_daughter: (3,) daughter node position [um] (the ``k``).

    Returns:
        The junction angle in [0, pi].
    """
    r1 = np.asarray(p_mother, float) - np.asarray(p_branch, float)
    r2 = np.asarray(p_daughter, float) - np.asarray(p_branch, float)
    d1 = float(np.linalg.norm(r1))
    d2 = float(np.linalg.norm(r2))
    if d1 < 1e-12 or d2 < 1e-12:
        raise ValueError("degenerate branch arm (coincident nodes)")
    c = float(np.clip(np.dot(r1, r2) / (d1 * d2), -1.0, 1.0))
    return float(np.arccos(c))


def angle_energy(theta: float | npt.NDArray[np.float64], k_theta: float, theta0: float) -> npt.NDArray[np.float64]:
    """Harmonic branch-angle potential ``U = 1/2 k_theta (theta - theta0)^2`` [pN*um]."""
    if k_theta < 0.0:
        raise ValueError("k_theta must be >= 0")
    dth = np.asarray(theta, float) - theta0
    return 0.5 * k_theta * dth * dth


def angle_forces(
    p_mother: npt.NDArray[np.float64],
    p_branch: npt.NDArray[np.float64],
    p_daughter: npt.NDArray[np.float64],
    k_theta: float,
    theta0: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Analytic nodal forces ``(F_mother, F_branch, F_daughter)`` [pN] of the harmonic angle potential.

    The standard harmonic-angle force (identical to ``ff.network_warp.branch_angle_kernel``):
    ``F_i = pref (r2/(d1 d2) - c r1/d1^2)``, ``F_k = pref (r1/(d1 d2) - c r2/d2^2)``,
    ``F_j = -(F_i + F_k)``, with ``pref = k_theta (theta - theta0) / sin(theta)``. Restores toward
    ``theta0`` and sums to zero (Newton's 3rd law for an internal potential).

    Returns:
        Tuple of three (3,) force vectors on (mother, branch, daughter).
    """
    r1 = np.asarray(p_mother, float) - np.asarray(p_branch, float)
    r2 = np.asarray(p_daughter, float) - np.asarray(p_branch, float)
    d1 = float(np.linalg.norm(r1))
    d2 = float(np.linalg.norm(r2))
    if d1 < 1e-12 or d2 < 1e-12:
        raise ValueError("degenerate branch arm (coincident nodes)")
    c = float(np.clip(np.dot(r1, r2) / (d1 * d2), -1.0, 1.0))
    s = float(np.sqrt(max(1.0 - c * c, 1e-12)))
    theta = float(np.arccos(c))
    pref = k_theta * (theta - theta0) / s
    f_i = pref * (r2 / (d1 * d2) - c * r1 / (d1 * d1))
    f_k = pref * (r1 / (d1 * d2) - c * r2 / (d2 * d2))
    f_j = -(f_i + f_k)
    return f_i, f_j, f_k


def thermal_sigma(k_theta: float, kT: float = KT_310K_PN_UM) -> float:
    """Equipartition angular width ``sigma_theta = sqrt(kT / k_theta)`` [rad] of the harmonic junction."""
    if k_theta <= 0.0:
        raise ValueError("k_theta must be > 0")
    return float(np.sqrt(kT / k_theta))


def k_theta_from_sigma(sigma_theta: float, kT: float = KT_310K_PN_UM) -> float:
    """Inverse equipartition: set ``k_theta = kT / sigma_theta^2`` [pN*um/rad^2] from an OBSERVED angular SD.

    This is the SETTING rule (how ``ARP23_K_THETA`` was fixed from the measured sigma_theta = 9 deg) — the
    angle stiffness is DERIVED from the measured thermal spread, never tuned to an outcome.
    """
    if sigma_theta <= 0.0:
        raise ValueError("sigma_theta must be > 0")
    return float(kT / (sigma_theta * sigma_theta))


def sample_branch_angles(
    theta0: float,
    k_theta: float,
    n: int,
    rng: np.random.Generator,
    kT: float = KT_310K_PN_UM,
    *,
    solid_angle_measure: bool = True,
) -> npt.NDArray[np.float64]:
    """Draw ``n`` branch angles from the Boltzmann equilibrium of the harmonic junction.

    With ``solid_angle_measure=True`` the physical distribution ``p(theta) ~ sin(theta) exp(-U/kT)`` is used
    (the free daughter arm's azimuth sweeps a solid angle ~ sin(theta)); with ``False`` the measure-free
    harmonic limit ``p(theta) ~ exp(-U/kT)`` is drawn, whose SD is EXACTLY ``sqrt(kT/k_theta)``. Because
    ``d(ln sin)/dtheta = cot(theta) > 0`` for ``theta0 = 70 deg < 90 deg``, the solid-angle measure shifts the
    peak slightly ABOVE ``theta0`` (toward pi/2) by ~ ``sigma^2 cot(theta0)`` — a physical effect, reported not
    tuned; the SD is essentially unchanged. Sampled by inverse-CDF on a fine grid (deterministic given ``rng``).

    Returns:
        (n,) sampled angles [rad] in (0, pi).
    """
    sigma = thermal_sigma(k_theta, kT)
    lo = max(1e-4, theta0 - 8.0 * sigma)
    hi = min(np.pi - 1e-4, theta0 + 8.0 * sigma)
    grid = np.linspace(lo, hi, 4096)
    logw = -angle_energy(grid, k_theta, theta0) / kT
    if solid_angle_measure:
        logw = logw + np.log(np.sin(grid))
    w = np.exp(logw - logw.max())
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    u = rng.random(n)
    return np.interp(u, cdf, grid)
