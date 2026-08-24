"""Sourced single-cell geometry as a measured DISTRIBUTION, not a point value (PI 2026-07-15).

Cell and nucleus size are STATISTICAL: real single-cell measurements give a distribution with
~15–30 % CV, not one number. Hard-coding ``R_cell``/``R_nuc = const`` (and, worse, the phantom
``Moore 2016 N:C 1.9`` citation that is NOT in the KB) is replaced here by the sourced MCF7 size
distributions + a deterministic sampler, so a quenched ensemble draws ``(R_cell, R_nuc)`` per
realization exactly as it draws network disorder (``ff/ENGINE.md`` §1, regime-A/B slow disorder →
parallel-realization Monte-Carlo). Bigger cells carry bigger nuclei (co-variation) at a ~conserved
nucleus:cell radius ratio (N:C).

Sourcing (MCF7, all PRIMARY — this replaces the un-registered "Moore 2016 N:C 1.9"):
  - Cell VOLUME 1760 µm³ = 1.76 pL → sphere-equivalent radius **7.49 µm**.  BioNumbers BNID 115154
    (Gamcsik, Millis & Colvin 1995, *Cancer Res* 55:2012, ³¹P/¹³C NMR; cited in Wagner, Venkataraman
    & Buettner 2011, *Free Radic Biol Med* 51:700, 10.1016/j.freeradbiomed.2011.05.024). The volume
    anchor — and volume is the conserved quantity under turgor — so it is the DEFAULT central radius.
  - Imaging flow cytometry, **n = 2164** suspended cells: cell dia 18.88 ± 2.86 µm (r 9.44 ± 1.43),
    nucleus dia 12.68 ± 1.94 µm (r 6.34 ± 0.97), **N:C radius ratio 0.68 ± 0.08**.  PMC7000884
    (10.1371/journal.pone... comparison study). Largest single distribution; the authors flag that
    IFC masking systematically OVER-estimates absolute size (→ the volume/photoacoustic anchor is
    ~1.25× smaller), so the SHAPE (CV, N:C) is taken from IFC and the central MEAN from volume.
  - Ultra-high-freq ultrasound / photoacoustic, n = 37, suspended-in-agarose: cell dia 15.2 ± 3.5 µm
    (r 7.6 — agrees with the volume anchor), nucleus dia 10.2 ± 3.5, **N:C 0.68 ± 0.19**.  Springer
    10.1007/s10765-016-2129-y. Independent confirmation of the 0.68 N:C ratio AND the ~7.5 µm radius.

Central defaults (MCF7, suspended): R_cell **7.5 µm** (volume anchor, matches photoacoustic), radius
CV **0.15** (IFC), N:C radius ratio **0.68 ± 0.08**. The IFC larger anchor (9.44 µm) is a selectable
method-upper (``MCF7_GEOMETRY_IFC``) so the ~1.25× method spread can be swept rather than hidden.

FF units: µm. Never sweep these to hit a mechanical target — they are measured geometry.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True, slots=True)
class CellGeometryDist:
    """A sourced single-cell geometry distribution (suspended sphere), FF units (µm)."""

    name: str
    R_cell_um: float        # central (mean) suspended cell radius [µm]
    R_cell_cv: float        # coefficient of variation of the cell radius (σ/µ)
    nc_ratio: float         # nucleus:cell RADIUS ratio (N:C), mean
    nc_ratio_sd: float      # SD of the N:C radius ratio
    R_cell_min_um: float    # truncation floor (reject non-physical tiny draws) [µm]
    R_cell_max_um: float    # truncation ceil  [µm]
    nc_min: float = 0.30    # N:C truncation (keep the nucleus inside the cell, above a floor)
    nc_max: float = 0.95

    @property
    def R_nuc_um(self) -> float:
        """Central nucleus radius = N:C · R_cell [µm]."""
        return self.nc_ratio * self.R_cell_um


# MCF7 suspended geometry — volume-anchored central radius, IFC-shaped spread, 0.68 N:C (all sourced above).
MCF7_GEOMETRY = CellGeometryDist(
    name="mcf7_volume_anchor", R_cell_um=7.5, R_cell_cv=0.15, nc_ratio=0.68, nc_ratio_sd=0.08,
    R_cell_min_um=5.5, R_cell_max_um=10.5)

# Method-upper alternative: the imaging-flow-cytometry central radius 9.44 µm (n=2164). Same N:C/CV.
# Use to sweep the ~1.25× volume-vs-IFC method spread (the authors attribute it to IFC masking).
MCF7_GEOMETRY_IFC = replace(MCF7_GEOMETRY, name="mcf7_ifc_anchor", R_cell_um=9.44,
                            R_cell_min_um=6.5, R_cell_max_um=13.5)

_REGISTRY = {g.name: g for g in (MCF7_GEOMETRY, MCF7_GEOMETRY_IFC)}
_ALIASES = {"mcf7": "mcf7_volume_anchor", "volume": "mcf7_volume_anchor", "ifc": "mcf7_ifc_anchor"}


def resolve_cell_geometry(name: str = "mcf7") -> CellGeometryDist:
    """Return the sourced :class:`CellGeometryDist` for ``name`` (registry key or alias)."""
    key = _ALIASES.get(name, name)
    if key not in _REGISTRY:
        raise ValueError(f"unknown geometry {name!r}; choose {sorted(_REGISTRY)} (+ aliases {sorted(_ALIASES)})")
    return _REGISTRY[key]


def central_geometry(dist: CellGeometryDist | str = "mcf7") -> tuple[float, float]:
    """The deterministic CENTRAL (mean) ``(R_cell, R_nuc)`` [µm] — single-realization / back-compat runs.

    This is the point that the old hard-coded ``R_cell=7.5, R_nuc=0.70·R`` approximated; here R_nuc uses
    the SOURCED N:C ratio 0.68 (was the phantom-cited 0.70), so R_nuc = 0.68·7.5 = 5.1 µm."""
    d = resolve_cell_geometry(dist) if isinstance(dist, str) else dist
    return d.R_cell_um, d.R_nuc_um


def sample_cell_geometry(rng: np.random.Generator, dist: CellGeometryDist | str = "mcf7") -> tuple[float, float]:
    """Draw one ``(R_cell, R_nuc)`` [µm] from the measured distribution (one ensemble realization).

    ``R_cell`` ~ truncated-normal(mean, cv) over ``[R_cell_min, R_cell_max]``; ``R_nuc = nc · R_cell`` with
    the N:C ratio ~ truncated-normal(nc_ratio, nc_ratio_sd) over ``[nc_min, nc_max]`` — so the nucleus
    CO-VARIES with the cell (bigger cell → bigger nucleus) at a ~conserved, measured N:C ratio. Deterministic
    in ``rng`` (pass a seeded ``np.random.default_rng`` per realization for a reproducible quenched ensemble)."""
    d = resolve_cell_geometry(dist) if isinstance(dist, str) else dist
    sigma = d.R_cell_um * d.R_cell_cv
    # rejection-free clamp is biased at the tails, so resample the (rare) out-of-range draws
    R_cell = float(rng.normal(d.R_cell_um, sigma))
    for _ in range(64):
        if d.R_cell_min_um <= R_cell <= d.R_cell_max_um:
            break
        R_cell = float(rng.normal(d.R_cell_um, sigma))
    R_cell = float(np.clip(R_cell, d.R_cell_min_um, d.R_cell_max_um))
    nc = float(rng.normal(d.nc_ratio, d.nc_ratio_sd))
    for _ in range(64):
        if d.nc_min <= nc <= d.nc_max:
            break
        nc = float(rng.normal(d.nc_ratio, d.nc_ratio_sd))
    nc = float(np.clip(nc, d.nc_min, d.nc_max))
    return R_cell, nc * R_cell


__all__ = ["CellGeometryDist", "MCF7_GEOMETRY", "MCF7_GEOMETRY_IFC", "resolve_cell_geometry",
           "central_geometry", "sample_cell_geometry"]
