"""Pure-NumPy sanity gates for live membrane pressure traction (no Warp device execution)."""

from __future__ import annotations

import numpy as np

from aleph.components.incumbent.membrane_pressure import (
    signed_volume,
    uniform_pressure_force_reference,
)
from aleph.laws.membrane_surface import build_membrane_mesh


def _sphere():
    return build_membrane_mesh(7.5, subdivisions=2)


def test_zero_pressure_jump_is_bit_zero() -> None:
    mesh = _sphere()
    force = uniform_pressure_force_reference(mesh.verts, mesh.faces, 0.0)
    assert np.array_equal(force, np.zeros_like(force))


def test_positive_pressure_pushes_outward_and_reverses_sign() -> None:
    mesh = _sphere()
    force = uniform_pressure_force_reference(mesh.verts, mesh.faces, 40.0)
    radial_work = np.einsum("ij,ij->i", force, mesh.verts)
    assert radial_work.mean() > 0.0
    assert np.array_equal(
        uniform_pressure_force_reference(mesh.verts, mesh.faces, -40.0), -force)


def test_closed_surface_uniform_pressure_has_zero_net_force() -> None:
    mesh = _sphere()
    force = uniform_pressure_force_reference(mesh.verts, mesh.faces, 40.0)
    # Scale by total absolute nodal load so this is a topology/round-off statement, not a dimensional tolerance.
    closure = np.linalg.norm(force.sum(axis=0)) / np.linalg.norm(force, axis=1).sum()
    assert closure < 128.0 * np.finfo(np.float64).eps


def test_pressure_force_is_volume_work_gradient() -> None:
    mesh = _sphere()
    pressure = 40.0
    rng = np.random.default_rng(20260717)
    direction = rng.standard_normal(mesh.verts.shape)
    direction /= np.linalg.norm(direction)
    force = uniform_pressure_force_reference(mesh.verts, mesh.faces, pressure)

    # Central-difference optimum h ~ eps^(1/3)*length for a smooth scalar; derived, not fitted.
    h = np.cbrt(np.finfo(np.float64).eps) * mesh.R_mem
    vp = signed_volume(mesh.verts + h * direction, mesh.faces)
    vm = signed_volume(mesh.verts - h * direction, mesh.faces)
    fd_work = pressure * (vp - vm) / (2.0 * h)
    analytic_work = float(np.sum(force * direction))
    scale = max(abs(fd_work), abs(analytic_work), np.finfo(np.float64).tiny)
    assert abs(fd_work - analytic_work) / scale < 32.0 * np.cbrt(np.finfo(np.float64).eps)


def test_virial_and_young_laplace_identity() -> None:
    mesh = _sphere()
    pressure = 40.0
    force = uniform_pressure_force_reference(mesh.verts, mesh.faces, pressure)
    volume = signed_volume(mesh.verts, mesh.faces)
    virial = float(np.sum(force * mesh.verts))
    assert np.isclose(virial, 3.0 * pressure * volume, rtol=256.0 * np.finfo(np.float64).eps)

    gamma_equiv = pressure * mesh.R_mem / 2.0
    pressure_resultant = pressure * 4.0 * np.pi * mesh.R_mem**2
    curvature_tension_resultant = (2.0 * gamma_equiv / mesh.R_mem) * 4.0 * np.pi * mesh.R_mem**2
    assert pressure_resultant == curvature_tension_resultant
