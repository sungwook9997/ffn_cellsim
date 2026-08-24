"""Cortical tension γ — method-of-planes estimator (FF engine, Stage 6d).

Faithful port of the validated HOOMD γ estimator
(``archive/hoomd_legacy/cortex/cortical_tension.py::_method_of_planes_gamma``) so the MD-free FF γ
is measured with the **same** instrument as the BAOAB-MD γ it will be compared against (the
decisive γ-floor experiment, ENGINE.md §4).

For each of ``n_planes`` near-isotropic diametral planes (Fibonacci normals) through the cell
centre, sum the tensile force ``T·|û·n̂|`` of every load-bearing element that crosses the plane,
divide by the cut circumference ``2πR``, and average ``|γ|`` over orientations:

    γ = mean_n̂ | (1/2πR) Σ_{elements crossing n̂}  T · |û·n̂| |

Elements are (endpoint A, endpoint B, signed scalar tension T): for the FF cortex these are the
inextensible actin segments (axial tension = constraint multiplier, ``constraints.segment_axial_
tension``), the crosslinker links, and the myosin links. ``|û·n̂|`` (absolute) is REQUIRED — it is
what makes an isotropic tension give the right √N→N scaling (see the source-module docstring).

Reference bands (literature; for OVERLAY only — γ is SWEPT as a controlled variable, never tuned to
hit a band, per the hard rule):
  * Salbreux, Charras & Paluch 2012:  γ ∈ [0.35, 0.65] mN/m  (= [3.5e-4, 6.5e-4] N/m).
  * MCF7 interphase (Hosseini 2020):  γ IQR ≈ [0.18, 0.40] mN/m.
In FF units γ has dimension force/length = pN/µm; **1 pN/µm = 1e-6 N/m = 1e-3 mN/m**, so the
Salbreux band [0.35, 0.65] mN/m = **[350, 650] pN/µm** (NOT 0.35–0.65 pN/µm — that 1000× slip is a
trap; γ values are large in pN/µm).
"""

from __future__ import annotations

import numpy as np

# Literature cortical-tension bands, in FF units (pN/µm). 1 pN/µm = 1e-3 mN/m, so 0.35 mN/m = 350.
# ⚠️ FRAMING (FF_STAGE6P, 2026-07-01): these bands are the TOTAL apparent cortical tension (AFM/Laplace),
# which is ~70% myosin-DEPENDENT (blebbistatin −68% Fischer-Friedrich 2016) + ~0.04 mN/m passive floor +
# a turgor Young-Laplace PARTNER (NOT additive — double-book) + minor membrane. The FF γ_active (method-
# of-planes, actomyosin only) must be compared against the ACTIVE FRACTION of the band (ACTIVE_FRACTION ·
# band ≈ 0.245–0.455 mN/m), NOT the full total band. Provenance NOTE: the (350,650) "Salbreux band" is a
# MISNOMER — its real empirical anchor is Chugh et al. 2017 Nat Cell Biol 19:689 (HeLa, interphase, AFM,
# T₀=230 pN/µm × peak ~1.6 ⇒ ~0.37 mN/m floor); Salbreux 2012 is a REVIEW carrying the broad 0.1–1 mN/m.
SALBREUX_BAND_PN_UM = (350.0, 650.0)   # TOTAL cortical tension band — Chugh 2017 (HeLa interphase AFM); not Salbreux
MCF7_IQR_PN_UM = (180.0, 400.0)        # MCF7 interphase SUSPENDED, Hosseini 2020 (0.18–0.40 mN/m); adherent = absent
ACTIVE_FRACTION = 0.70                  # myosin-dependent fraction of the total band (blebbistatin; FF_STAGE6P)


def active_band_pn_um(band=SALBREUX_BAND_PN_UM):
    """The myosin-ACTIVE target = ACTIVE_FRACTION · total band [pN/µm] — what FF γ_active is compared to
    (not the full total band; FF_STAGE6P). e.g. (350,650) → (245,455)."""
    return (ACTIVE_FRACTION * band[0], ACTIVE_FRACTION * band[1])


def fibonacci_plane_normals(n_planes: int) -> np.ndarray:
    """``n_planes`` near-isotropic unit normals on S² (Fibonacci lattice). Ported verbatim."""
    phi = (1.0 + 5.0**0.5) / 2.0
    i = np.arange(n_planes, dtype=np.float64)
    z = 1.0 - 2.0 * (i + 0.5) / n_planes
    rxy = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    theta = 2.0 * np.pi * i / phi
    return np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)


def method_of_planes_gamma(rA: np.ndarray, rB: np.ndarray, tension: np.ndarray,
                           R: float, *, n_planes: int = 50, centre: np.ndarray | None = None) -> float:
    """Method-of-planes cortical tension γ [force/length] from per-element scalar tensions.

    Args:
        rA: (M, 3) element endpoint A positions.
        rB: (M, 3) element endpoint B positions.
        tension: (M,) signed scalar element tensions (+ tension, − compression).
        R: cell radius (great-circle radius for the cut circumference 2πR).
        n_planes: number of diametral cut-plane orientations to average over.
        centre: cell centre (default = mean of all endpoints) — subtracted before crossing tests.

    Returns:
        mean over orientations of |γ|. 0.0 if there are no elements.
    """
    rA = np.asarray(rA, dtype=np.float64)
    rB = np.asarray(rB, dtype=np.float64)
    if rA.shape[0] == 0:
        return 0.0
    T = np.asarray(tension, dtype=np.float64)
    if centre is None:
        centre = 0.5 * (rA.mean(axis=0) + rB.mean(axis=0))
    a = rA - centre
    b = rB - centre
    d = b - a
    ln = np.linalg.norm(d, axis=1)
    safe = ln > 0
    u = np.zeros_like(d)
    u[safe] = d[safe] / ln[safe, None]

    normals = fibonacci_plane_normals(n_planes)
    circ = 2.0 * np.pi * R
    gammas = []
    for n_hat in normals:
        a_side = a @ n_hat
        b_side = b @ n_hat
        crossing = (a_side * b_side) < 0.0
        if not crossing.any():
            gammas.append(0.0)
            continue
        f_cut = T[crossing] * np.abs(u[crossing] @ n_hat)
        gammas.append(float(np.sum(f_cut)) / circ)
    return float(np.mean(np.abs(np.asarray(gammas))))


def gamma_to_mN_per_m(gamma_pn_um: float) -> float:
    """Convert γ from FF units (pN/µm) to mN/m. 1 pN/µm = 1e-3 mN/m."""
    return gamma_pn_um * 1e-3


def gamma_to_N_per_m(gamma_pn_um: float) -> float:
    """Convert γ from FF units (pN/µm) to SI N/m. 1 pN/µm = 1e-6 N/m."""
    return gamma_pn_um * 1e-6
