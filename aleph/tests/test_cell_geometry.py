"""Sourced MCF7 cell/nucleus geometry distribution — sampler sanity gates.

The geometry is a MEASURED distribution (PI 2026-07-15), not a point value. These gates check the
sampler reproduces the sourced statistics (mean, CV, N:C ratio), stays inside the physical truncation,
co-varies nucleus with cell, and is deterministic in the rng (reproducible quenched ensemble).
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.laws.cell_geometry import (
    MCF7_GEOMETRY, MCF7_GEOMETRY_IFC, resolve_cell_geometry, central_geometry, sample_cell_geometry)


def test_central_is_volume_anchor():
    """Central geometry = the sourced volume anchor (7.5 µm) with the sourced N:C 0.68 (not 0.70)."""
    R_cell, R_nuc = central_geometry("mcf7")
    assert R_cell == pytest.approx(7.5)
    assert R_nuc == pytest.approx(0.68 * 7.5, abs=1e-9)      # 5.1 µm — sourced N:C, replaces phantom 0.70
    # sphere-equivalent volume matches BioNumbers 1760 µm³ (Gamcsik1995/Wagner2011) within rounding
    V = 4.0 / 3.0 * np.pi * R_cell ** 3
    assert V == pytest.approx(1767.0, rel=0.02)


def test_sample_recovers_distribution_statistics():
    """A large ensemble recovers the sourced mean radius, ~15 % CV, and 0.68 N:C ratio."""
    rng = np.random.default_rng(0)
    R = np.array([sample_cell_geometry(rng, MCF7_GEOMETRY) for _ in range(20000)])
    R_cell, R_nuc = R[:, 0], R[:, 1]
    assert R_cell.mean() == pytest.approx(7.5, abs=0.15)
    assert (R_cell.std() / R_cell.mean()) == pytest.approx(0.15, abs=0.03)   # IFC radius CV
    nc = R_nuc / R_cell
    assert nc.mean() == pytest.approx(0.68, abs=0.02)                        # sourced N:C radius ratio
    assert nc.std() == pytest.approx(0.08, abs=0.02)


def test_truncation_keeps_physical():
    """No draw leaves the physical band; the nucleus is always strictly inside the cell."""
    rng = np.random.default_rng(1)
    for _ in range(5000):
        R_cell, R_nuc = sample_cell_geometry(rng, MCF7_GEOMETRY)
        assert MCF7_GEOMETRY.R_cell_min_um <= R_cell <= MCF7_GEOMETRY.R_cell_max_um
        assert 0.0 < R_nuc < R_cell                                         # nucleus fits inside


def test_nucleus_covaries_with_cell():
    """Bigger cells carry bigger nuclei (positive R_cell↔R_nuc correlation), not an independent draw."""
    rng = np.random.default_rng(2)
    R = np.array([sample_cell_geometry(rng, MCF7_GEOMETRY) for _ in range(5000)])
    corr = np.corrcoef(R[:, 0], R[:, 1])[0, 1]
    assert corr > 0.7                                                       # co-variation, not independence


def test_deterministic_in_rng():
    """Same seed → identical draws (reproducible quenched ensemble)."""
    a = [sample_cell_geometry(np.random.default_rng(7), MCF7_GEOMETRY) for _ in range(3)]
    b = [sample_cell_geometry(np.random.default_rng(7), MCF7_GEOMETRY) for _ in range(3)]
    assert a == b


def test_ifc_anchor_is_larger():
    """The IFC method-upper anchor is ~1.25× the volume anchor (the real method spread, exposed not hidden)."""
    assert MCF7_GEOMETRY_IFC.R_cell_um / MCF7_GEOMETRY.R_cell_um == pytest.approx(9.44 / 7.5, rel=1e-6)
    assert resolve_cell_geometry("ifc").nc_ratio == pytest.approx(0.68)     # same N:C, only the mean shifts
