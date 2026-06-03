"""Closed-form acceptance oracle: single-cell cortical tension -> tissue surface tension.

VALIDATION oracle (runtime-import-forbidden). Encodes the published *bridge* that links the
single-cell cortical-tension observable (KU-3.5, gamma) to the multicellular spheroid's
emergent aggregate surface tension and the ``A/A0 = a + b/R + c/R^2`` spreading law. The PI
thesis (2026-06-03) is the literature default: **single-cell cortical tension is the root of
spheroid/aggregate surface tension** — and three tissue-tension papers supply the exact
relations that make the chain quantitative + testable.

The chain (all literature-anchored; see ``docs/CORTICAL_TENSION_TRIAGE_2026-06-03.md``)::

    gamma  (single-cell cortical tension, KU-3.5; g_rigid native 0.57 mN/m IN band)
      - beta (E-cadherin adhesion energy density)        [Brodland DITH; Okuda 2026]
      = Gamma_cc  (interior cell-cell interfacial tension)
    sigma_tissue (aggregate FREE-surface tension) = gamma (surface cells' free cortex)
                                                    [Roffay 2021: outer ~1.6-2x interior]
    Young-Laplace:  dP = sigma (1/R + 1/R')           [Roffay 2021, 3D mean curvature]
    radius balance:  R = lambda / sigma                [Fastabend 2026, 2D analog]
    => the 1/R curvature scaling that the b/R term of A/A0 encodes.

Provenance / discipline
-----------------------
Every relation below is a PUBLISHED paper-model and therefore enters ONLY as a validation
oracle / literature anchor — NEVER a runtime mechanism (CLAUDE.md inversion rule). No physics
is computed in the runtime from this module; the Layer-2 CBM produces sigma *emergently* (via
the virial-pressure observable in ``ffn_sim.spheroid.observables``), and this oracle is the
acceptance target it is checked against.

Sources:
- Chugh et al. 2017, Nat Cell Biol 19:689 (DOI 10.1038/ncb3525) -- cortical-tension band /
  magnitude oracle (model peak ~0.37 mN/m = the KU-3.5 band floor; T0 = 230 pN/um).
- Roffay, Chan, Guirao, Hiiragi, Graner 2021, Development 148:dev192773 -- Young-Laplace
  ``dP = t*K``, 3D ``K = 1/R + 1/R'``; outer (cell-medium) tension ~1.6-2x interior cell-cell
  (mouse embryo, pipette-validated => PROXY ratio, absolute a.u. do not transfer).
- Fastabend, Bidan, Dunlop, Kollmannsberger 2026, Phys Rev E 113:024403 -- ``R = lambda/sigma``
  2D Young-Laplace analog; tissue growth-rate proportional to local interfacial curvature.
- Thiticharoentam, ..., Okuda 2026, bioRxiv 2026.03.17.712503 -- ``Gamma = cortical - adhesion``;
  monolayer -> 3D cap transition once free-surface tension exceeds ~0.2 * cell-cell tension.
- Tissue surface-tension magnitude PROXY: MCF10DCIS ~21 mN/m (Nagle 2022) -- NO MCF7 datum,
  flagged proxy. Single-cell MCF7 anchors are MCF7-specific (Wagner 2011 size, Iturri 2020
  de-adhesion); the aggregate-sigma band is proxy until an MCF7 tissue-tensiometry datum.

Sanity Gate
-----------
- Dimensional: tensions [N/m]; R, R' [m]; pressure [Pa = N/m^2]; energy density / adhesion
  tension [J/m^2 = N/m]; work [J]; area [m^2]. ``young_laplace_pressure`` returns Pa.
- Boundary: R -> inf  => dP -> 0 (flat interface, no Laplace overpressure). adhesion -> 0
  (free surface, no cohesive partner) => Gamma = cortical = aggregate surface tension.
- Sign-sense: aggregate surface tension > 0; cohesive interior tension Gamma_cc >= 0 (adhesion
  reduces but cannot invert it -- clamped at 0; adhesion > cortical would mean wetting/no
  surface, surfaced to the caller as Gamma_cc = 0). Young-Laplace dP > 0 (interior over-pressured).
- Inversion: ``surface_tension_from_pressure`` is the exact inverse of ``young_laplace_pressure``
  (round-trips to machine precision -- test).
- Convention (stated, not silent): the differential-interfacial-tension hypothesis (Brodland
  2002) writes a cell-cell interface tension as cortical tension minus adhesion. This module
  uses the *reduced single-cortex-per-face* form ``Gamma_cc = gamma - beta`` (``beta`` = adhesion
  energy density per face) so that the free surface (no adhesive partner) has ``sigma_tissue =
  gamma`` and the Roffay outer/interior ratio ``gamma / Gamma_cc = 1 / (1 - beta/gamma)``
  reproduces ~1.6-2x at ``beta/gamma ~ 0.375-0.5``. (Factor-of-2 conventions for ``beta`` exist
  in the literature; this module fixes one and documents it rather than leaving it ambiguous.)
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "CORTICAL_TENSION_BAND_MN_M",
    "G_RIGID_NATIVE_MN_M",
    "ROFFAY_SURFACE_INTERIOR_RATIO",
    "TISSUE_SIGMA_PROXY_MN_M",
    "OKUDA_CAP_THRESHOLD",
    "adhesion_tension",
    "interfacial_tension",
    "aggregate_surface_tension",
    "surface_interior_ratio",
    "young_laplace_pressure",
    "surface_tension_from_pressure",
    "laplace_radius",
    "is_three_d_cap",
]

# --- Literature anchors (cross-check only; provenance in module docstring) ----------------

#: KU-3.5 cortical-tension acceptance band [mN/m] (= [pN/um]/1000). Chugh model peak ~0.37.
CORTICAL_TENSION_BAND_MN_M: tuple[float, float] = (0.35, 0.65)

#: KU-3.5 structural/passive (g_rigid) native cortical tension [mN/m] -- IN band (STAGE-2).
G_RIGID_NATIVE_MN_M: float = 0.57

#: Roffay 2021 outer(cell-medium)/interior(cell-cell) tension ratio [dimensionless] -- PROXY
#: (mouse 8/16-cell embryo, pipette-validated; ratio transfers, absolute a.u. do not).
ROFFAY_SURFACE_INTERIOR_RATIO: tuple[float, float] = (1.6, 2.0)

#: Aggregate tissue surface-tension magnitude PROXY [mN/m] -- MCF10DCIS (Nagle 2022). NO MCF7
#: datum; do NOT treat as an MCF7-specific acceptance band (flagged proxy).
TISSUE_SIGMA_PROXY_MN_M: float = 21.0

#: Okuda 2026: monolayer -> 3D stratified cap once free-surface tension exceeds this fraction
#: of the cell-cell (lateral) tension [dimensionless].
OKUDA_CAP_THRESHOLD: float = 0.2


# --- Gamma = cortical - adhesion (Brodland DITH; Okuda 2026) ------------------------------

def adhesion_tension(work_of_deadhesion: float, contact_area: float) -> float:
    """Adhesion energy density (tension contribution) ``beta = W / A_contact`` [N/m].

    The energy released per unit contact area when two cells de-adhere. Equivalently the
    adhesion tension that opposes cortical tension at a cell-cell interface (Brodland 2002).

    Args:
        work_of_deadhesion: Work to separate one cell-cell contact [J] (> 0). E.g. MCF7-MCF7
            single-cell force-spectroscopy work (~200 fJ scale, Iturri 2020).
        contact_area: Cell-cell contact-zone area [m^2] (> 0).

    Returns:
        Adhesion energy density beta [J/m^2 = N/m].

    Raises:
        ValueError: if ``work_of_deadhesion < 0`` or ``contact_area <= 0``.
    """
    if work_of_deadhesion < 0.0:
        raise ValueError("work_of_deadhesion must be non-negative [J].")
    if contact_area <= 0.0:
        raise ValueError("contact_area must be strictly positive [m^2].")
    return work_of_deadhesion / contact_area


def interfacial_tension(cortical_tension: float, adhesion: float) -> float:
    """Interior cell-cell interfacial tension ``Gamma_cc = gamma - beta`` (DITH; Okuda 2026).

    The reduced single-cortex-per-face form of the differential-interfacial-tension hypothesis
    (Brodland 2002): cortical tension minus the adhesion energy density at the shared face. See
    the module Sanity Gate for the convention. Clamped at 0 -- adhesion exceeding cortical
    tension means wetting (no cohesive surface), surfaced as ``Gamma_cc = 0`` rather than a
    negative tension.

    Args:
        cortical_tension: Single-cell cortical tension gamma [N/m] (>= 0).
        adhesion: Adhesion energy density beta [N/m] (>= 0) -- see ``adhesion_tension``.

    Returns:
        Interior cell-cell interfacial tension Gamma_cc [N/m], >= 0.

    Raises:
        ValueError: on negative inputs.
    """
    if cortical_tension < 0.0 or adhesion < 0.0:
        raise ValueError("cortical_tension and adhesion must be non-negative [N/m].")
    return float(max(cortical_tension - adhesion, 0.0))


def aggregate_surface_tension(cortical_tension: float) -> float:
    """Aggregate FREE-surface tension ``sigma_tissue = gamma`` (Roffay 2021).

    The aggregate's outer boundary is the surface cells' free cortex (no adhesive partner on
    the medium side), so the tissue surface tension equals the single-cell cortical tension.
    Roffay's inference finds this outer tension ~1.6-2x the interior cell-cell tension, which
    this identity reproduces together with ``interfacial_tension`` (adhesion lowers the
    interior but not the free surface).

    Args:
        cortical_tension: Single-cell cortical tension gamma [N/m] (>= 0).

    Returns:
        Aggregate surface tension sigma_tissue [N/m].
    """
    if cortical_tension < 0.0:
        raise ValueError("cortical_tension must be non-negative [N/m].")
    return float(cortical_tension)


def surface_interior_ratio(cortical_tension: float, adhesion: float) -> float:
    """Predicted outer/interior tension ratio ``sigma_tissue / Gamma_cc`` (Roffay check).

    Returns ``inf`` if the interior tension clamps to 0 (full wetting). The Roffay-anchored
    target band is :data:`ROFFAY_SURFACE_INTERIOR_RATIO` (~1.6-2.0, PROXY).

    Args:
        cortical_tension: gamma [N/m] (>= 0).
        adhesion: beta [N/m] (>= 0).

    Returns:
        Dimensionless ratio (``inf`` if interior tension is 0).
    """
    gamma_cc = interfacial_tension(cortical_tension, adhesion)
    sigma = aggregate_surface_tension(cortical_tension)
    if gamma_cc <= 0.0:
        return float("inf")
    return float(sigma / gamma_cc)


# --- Young-Laplace bridge (Roffay 2021, 3D; Fastabend 2026) -------------------------------

def young_laplace_pressure(
    sigma: float, R: float, R2: float | None = None
) -> float:
    """Young-Laplace overpressure ``dP = sigma * (1/R + 1/R')`` [Pa] (Roffay 2021, 3D).

    Interior over-pressure of an aggregate of surface tension ``sigma`` and mean curvature
    ``K = 1/R + 1/R'``. For a sphere (``R2 = None`` => ``R' = R``) this is ``2*sigma/R``.

    Args:
        sigma: Aggregate surface tension [N/m] (>= 0).
        R: First principal radius of curvature [m] (> 0).
        R2: Second principal radius [m] (> 0); defaults to ``R`` (sphere).

    Returns:
        Laplace overpressure dP [Pa].

    Raises:
        ValueError: if ``sigma < 0`` or any radius ``<= 0``.
    """
    if sigma < 0.0:
        raise ValueError("sigma must be non-negative [N/m].")
    if R <= 0.0 or (R2 is not None and R2 <= 0.0):
        raise ValueError("radii must be strictly positive [m].")
    r2 = R if R2 is None else R2
    return float(sigma * (1.0 / R + 1.0 / r2))


def surface_tension_from_pressure(
    delta_P: float, R: float, R2: float | None = None
) -> float:
    """Invert Young-Laplace: ``sigma = dP / (1/R + 1/R')`` [N/m] (the measurement route).

    The exact inverse of :func:`young_laplace_pressure`. Given an aggregate's measured interior
    overpressure (e.g. from the virial-pressure observable) and its curvature, recover the
    emergent surface tension -- the validation observable for the bridge.

    Args:
        delta_P: Interior Laplace overpressure [Pa] (>= 0).
        R: First principal radius [m] (> 0).
        R2: Second principal radius [m] (> 0); defaults to ``R`` (sphere).

    Returns:
        Surface tension sigma [N/m].

    Raises:
        ValueError: if ``delta_P < 0`` or any radius ``<= 0``.
    """
    if delta_P < 0.0:
        raise ValueError("delta_P must be non-negative [Pa].")
    if R <= 0.0 or (R2 is not None and R2 <= 0.0):
        raise ValueError("radii must be strictly positive [m].")
    r2 = R if R2 is None else R2
    return float(delta_P / (1.0 / R + 1.0 / r2))


def laplace_radius(line_tension: float, sigma: float) -> float:
    """Fastabend 2026 radius-tension balance ``R = lambda / sigma`` [m] (2D Young-Laplace analog).

    The equilibrium radius at which a line/cortical tension ``lambda`` balances a surface
    tension / internal pressure ``sigma`` -- the structural origin of the 1/R curvature scaling
    in ``A/A0 = a + b/R + c/R^2``.

    Args:
        line_tension: lambda [N] (line tension) or [N/m] in the 2D analog (> 0).
        sigma: surface tension [N/m] or internal pressure proxy (> 0).

    Returns:
        Balance radius R [m] (units follow the inputs: [N]/[N/m] = [m]).

    Raises:
        ValueError: if ``line_tension <= 0`` or ``sigma <= 0``.
    """
    if line_tension <= 0.0 or sigma <= 0.0:
        raise ValueError("line_tension and sigma must be strictly positive.")
    return float(line_tension / sigma)


def is_three_d_cap(
    free_surface_tension: float, cell_cell_tension: float
) -> bool:
    """Okuda 2026 monolayer -> 3D-cap criterion: free-surface tension > 0.2 * cell-cell tension.

    Underwrites the L2.6 finding that cohesive MCF7 forms a 3D cap (not a monolayer). When the
    free-surface (cortical) tension exceeds :data:`OKUDA_CAP_THRESHOLD` of the cell-cell
    (lateral) tension, the aggregate stratifies into a 3D cap.

    Args:
        free_surface_tension: aggregate free-surface (cortical) tension [N/m] (>= 0).
        cell_cell_tension: interior cell-cell (lateral) tension [N/m] (> 0).

    Returns:
        True if the aggregate is predicted to be a 3D cap (stratified), False if a monolayer.

    Raises:
        ValueError: if ``free_surface_tension < 0`` or ``cell_cell_tension <= 0``.
    """
    if free_surface_tension < 0.0:
        raise ValueError("free_surface_tension must be non-negative [N/m].")
    if cell_cell_tension <= 0.0:
        raise ValueError("cell_cell_tension must be strictly positive [N/m].")
    return bool(free_surface_tension > OKUDA_CAP_THRESHOLD * cell_cell_tension)
