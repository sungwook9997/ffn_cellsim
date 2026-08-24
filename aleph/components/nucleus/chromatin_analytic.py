"""Closed-form chromatin-polymer references — the internal polymer net of the deformable nucleus.

Framework #6: the nucleus interior is a CHROMATIN polymer net (not a bulk modulus, not a radial
spring) governing the SMALL-strain response, with lamin-A/C taking over past the knee (``lamina_analytic``).
The fine-grained primitive is a worm-like chain (WLC); the internal net is a crosslinked population of
them tethered to the envelope. Host oracles (pure NumPy, NO Warp).

Marko–Siggia interpolation for a single WLC of contour length L and persistence length ℓ_p:

    f(x) = (k_BT/ℓ_p) [ x/L + 1/(4(1−x/L)²) − 1/4 ]      tension at end-to-end extension x

    small-x  → f ≈ (3 k_BT / (2 ℓ_p L)) x                entropic spring constant k_wlc
    x → L    → f → ∞                                     finite extensibility (never interpenetrates)

Energy E(x) = ∫₀ˣ f dx' (so dE/dx = f — the FD-gradient sign arbiter). k_BT in FF units: at 310 K,
k_BT ≈ 4.28e-3 pN·µm.

References: Marko & Siggia 1995 Macromolecules 28:8759; Bustamante 1994. Chromatin as a WLC network:
Stephens 2017 MBoC 28:1984 (chromatin governs small-strain nuclear rigidity; lamin-A/C large-strain).
"""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "KBT_PN_UM_310K",
    "wlc_tension",
    "wlc_energy",
    "wlc_small_strain_stiffness",
]

KBT_PN_UM_310K = 4.28e-3        # k_BT at 310 K [pN·µm] (1.380649e-23 · 310 / 1e-18·1e-12 units)


def wlc_tension(
    x: npt.ArrayLike, contour_length: float, persistence_length: float, kbt: float = KBT_PN_UM_310K
) -> npt.NDArray[np.float64]:
    """Marko–Siggia WLC tension f(x) [pN] at end-to-end extension ``x`` (0 ≤ x < L). Monotone
    increasing; diverges as x→L (finite extensibility). ``x`` may be an array."""
    xr = np.asarray(x, dtype=np.float64)
    L = float(contour_length)
    if L <= 0.0 or persistence_length <= 0.0:
        raise ValueError("contour_length and persistence_length must be positive")
    r = np.clip(xr / L, 0.0, 1.0 - 1e-12)                 # guard the x→L pole
    return (kbt / persistence_length) * (r + 1.0 / (4.0 * (1.0 - r) ** 2) - 0.25)


def wlc_energy(
    x: npt.ArrayLike, contour_length: float, persistence_length: float, kbt: float = KBT_PN_UM_310K
) -> npt.NDArray[np.float64]:
    """WLC stored energy E(x) [pN·µm] with ``dE/dx = f(x)`` and E(0)=0 (the FD ground truth):

        E(x) = (k_BT/ℓ_p) [ x²/(2L) + L/(4(1−x/L)) − L/4 − x/4 ].
    """
    xr = np.asarray(x, dtype=np.float64)
    L = float(contour_length)
    if L <= 0.0 or persistence_length <= 0.0:
        raise ValueError("contour_length and persistence_length must be positive")
    r = np.clip(xr / L, 0.0, 1.0 - 1e-12)
    return (kbt / persistence_length) * (
        L * r * r / 2.0 + L / (4.0 * (1.0 - r)) - L / 4.0 - L * r / 4.0
    )


def wlc_small_strain_stiffness(
    contour_length: float, persistence_length: float, kbt: float = KBT_PN_UM_310K
) -> float:
    """Entropic spring constant in the small-extension limit ``k_wlc = 3 k_BT / (2 ℓ_p L)`` [pN/µm].
    The chromatin net's soft small-strain stiffness feeds the lamina ``k_chrom`` (I0-B2)."""
    return float(3.0 * kbt / (2.0 * persistence_length * contour_length))
