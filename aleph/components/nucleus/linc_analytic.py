"""LINC tether + pre-authored I7 load-path oracles — host references (pure NumPy, NO Warp).

Two roles:

  (I2, now) LINC tether force-extension — the nucleus↔cytoskeleton coupling (nesprin/SUN). Nesprin is
  a NONLINEAR spring (its spectrin/Ig domains stiffen/unfold under load), so the tether is a
  strain-stiffening finite-extensible link, force-free at its formation length. ⚠ k_linc is a GAP
  (I0-B7): "8 pN is a TENSION, not a stiffness" (§4C) — the flatten magnitude gate is INVALID until
  it is sourced; only the STRUCTURAL law (monotone, force-free at rest, dE/dx=f) is gated here.

  (I7, PRE-AUTHORED standalone) the load-path capstan + oblate-flatten oracles, so I7 (cap→LINC→nucleus,
  gated behind I3+I4) can be built against a green analytic suite the moment its live inputs land:
    * CAPSTAN ∮ T·κ ds — a contractile perinuclear-actin-cap FIBER of tension T draped over the
      curved nucleus presses INWARD by T·κ per unit length (tangential tension → normal pressure,
      NOT a radial strut). ∮ T·κ (n̂·d̂) ds over a wrap = the net directional inward force.
    * OBLATE VOLUME CONSERVATION — flattening at constant nucleoplasm volume: V=(4/3)π a²c=(4/3)π a³/A
      (A=aspect a/c) ⇒ at fixed V the equatorial radius grows a ∝ A^{1/3}. Adherent flatten is
      measured at the RESTING baseline with volume conserved (Khatau), NOT chased at AFM strain.

References: capstan/belt over a curve p=T/R=Tκ (classical); nesprin nonlinearity Arsenovic 2016 (nesprin
tension sensor ~8 pN), Bouzid 2019; oblate flatten Khatau 2009 PNAS 106:19017.
"""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "linc_tether_force",
    "linc_tether_energy",
    "capstan_normal_pressure",
    "capstan_line_integral",
    "oblate_equatorial_radius_at_constant_volume",
    "oblate_aspect_from_flatten",
]


# ---------------------------------------------------------------------------------------------------
# (I2) LINC tether — strain-stiffening nonlinear finite-extensible spring
# ---------------------------------------------------------------------------------------------------
def linc_tether_force(
    length: npt.ArrayLike, rest: float, k_linc: float, stiffening: float = 0.0
) -> npt.NDArray[np.float64]:
    """LINC tether tension along the link [pN] at end-to-end ``length`` (rest length ``rest``):

        f(Δ) = k_linc·Δ·(1 + stiffening·Δ²)          Δ = length − rest,   f(0) = 0

    A cubic-stiffening tether (``stiffening`` [1/µm²], 0 → Hookean). Only tension (Δ>0) transmits —
    a tether does not push. ``k_linc`` and ``stiffening`` are GAP (I0-B7): the magnitude is not gated,
    the STRUCTURE (monotone in Δ≥0, force-free at rest) is."""
    d = np.asarray(length, dtype=np.float64) - float(rest)
    d = np.maximum(d, 0.0)                                # tether: no compression
    return k_linc * d * (1.0 + stiffening * d * d)


def linc_tether_energy(
    length: npt.ArrayLike, rest: float, k_linc: float, stiffening: float = 0.0
) -> npt.NDArray[np.float64]:
    """LINC stored energy with ``dE/dΔ = f`` and E(rest)=0 (FD ground truth):
    ``E = k_linc(½Δ² + ¼·stiffening·Δ⁴)`` for Δ≥0, else 0."""
    d = np.asarray(length, dtype=np.float64) - float(rest)
    d = np.maximum(d, 0.0)
    return k_linc * (0.5 * d * d + 0.25 * stiffening * d ** 4)


# ---------------------------------------------------------------------------------------------------
# (I7, pre-authored) capstan — tangential cap tension → normal pressure on the nucleus
# ---------------------------------------------------------------------------------------------------
def capstan_normal_pressure(tension: npt.ArrayLike, curvature: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Inward normal force per unit fiber length ``p⊥ = T·κ`` [pN/µm] of a cap fiber of tension T [pN]
    following a curve of curvature κ=1/R [1/µm]. The capstan relation: tension over a curved surface
    is a NORMAL pressure, not a radial strut. (Frictionless/geometric belt-over-cylinder.)"""
    return np.asarray(tension, dtype=np.float64) * np.asarray(curvature, dtype=np.float64)


def capstan_line_integral(
    tension: npt.ArrayLike,
    curvature: npt.ArrayLike,
    ds: npt.ArrayLike,
    normal: npt.NDArray[np.float64],
    direction: npt.ArrayLike,
) -> float:
    """Net directional inward force ``∮ T·κ (n̂·d̂) ds`` [pN] from a cap draped along a discretized
    curve: per-segment normal pressure T·κ projected onto unit ``direction`` d̂, weighted by arc
    length ds. ``normal`` (Ns,3) are inward unit normals; ``tension``/``curvature``/``ds`` are (Ns,).
    For a full symmetric wrap the transverse components cancel and the axial component accumulates."""
    T = np.asarray(tension, dtype=np.float64)
    kappa = np.asarray(curvature, dtype=np.float64)
    dl = np.asarray(ds, dtype=np.float64)
    d_hat = np.asarray(direction, dtype=np.float64)
    d_hat = d_hat / (np.linalg.norm(d_hat) + 1e-300)
    proj = normal @ d_hat                                 # n̂·d̂ per segment
    return float(np.sum(T * kappa * proj * dl))


# ---------------------------------------------------------------------------------------------------
# (I7, pre-authored) oblate flatten at constant nucleoplasm volume
# ---------------------------------------------------------------------------------------------------
def oblate_equatorial_radius_at_constant_volume(v0: float, aspect: float) -> float:
    """Equatorial radius ``a`` [µm] of an oblate spheroid of aspect A=a/c enclosing a fixed volume V0:
    ``V0=(4/3)π a³/A`` ⇒ ``a = (3 V0 A / 4π)^{1/3}``. As the nucleus flattens (A↑) at conserved
    nucleoplasm volume, the equatorial radius grows as A^{1/3} (the I7 flatten kinematics)."""
    if v0 <= 0.0 or aspect < 1.0:
        raise ValueError("v0 must be > 0 and aspect ≥ 1")
    return float((3.0 * v0 * aspect / (4.0 * np.pi)) ** (1.0 / 3.0))


def oblate_aspect_from_flatten(height_ratio: float) -> float:
    """Aspect A=a/c from a measured polar-height reduction: if the resting sphere (A=1, height 2R)
    flattens to polar height ``height_ratio·2R`` at constant volume, then c=height_ratio·R·(a/R) and
    solving V-conservation gives A = height_ratio^{-3/2}. (Height ratio 1 → sphere; 0.5 → A≈2.83.)"""
    if not (0.0 < height_ratio <= 1.0):
        raise ValueError("height_ratio must be in (0, 1]")
    return float(height_ratio ** (-1.5))
