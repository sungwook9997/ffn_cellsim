"""``laws.volume`` — the one enclosed-volume gradient, and the lemma that lets it be one.

``dV/dx`` was written three times in two different algebraic forms.  These tests pin the lemma that
makes those forms the same thing, AND the condition under which they stop being the same thing —
because a lemma whose precondition is never exercised is an assumption with better handwriting.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.laws import volume as volume_mod
from aleph.laws.surface_manifold import icosphere
from aleph.laws.volume import (consistent_load_reference, is_closed_manifold, mesh_volume,
                               pressure_force_reference, volume_gradient)

DELTA_P = 40.0          # pN/um^2; a probe value for the algebra, not a sourced physiological pressure
RADIUS = 5.0


@pytest.fixture()
def mesh():
    v, f = icosphere(2, RADIUS)
    return np.ascontiguousarray(v, np.float64), np.ascontiguousarray(f, np.int64)


def test_self_check_passes() -> None:
    volume_mod._demo()


def test_the_two_algebraic_forms_agree_on_a_closed_mesh(mesh) -> None:
    v, f = mesh
    own = pressure_force_reference(v, f, DELTA_P)
    fem = consistent_load_reference(v, f, DELTA_P)
    assert np.abs(own - fem).max() <= 64 * np.finfo(np.float64).eps * np.abs(own).max()


def test_they_do_not_agree_on_an_open_one(mesh) -> None:
    """The lemma's precondition, exercised.  Half an icosphere is not a closed manifold."""
    v, f = mesh
    patch = f[: f.shape[0] // 2]
    assert not is_closed_manifold(patch)
    d = np.abs(pressure_force_reference(v, patch, DELTA_P)
               - consistent_load_reference(v, patch, DELTA_P)).max()
    assert d > 1e-6 * np.abs(pressure_force_reference(v, f, DELTA_P)).max()


def test_assembled_gradient_is_translation_invariant(mesh) -> None:
    """``V`` is translation-invariant on a closed mesh, so its gradient must be — though no face's is.

    This is the check that fails first when a launch covers only part of a surface's face range.
    """
    v, f = mesh
    g = volume_gradient(v, f)
    far = volume_gradient(v + np.array([1.0e3, -7.0e2, 3.0e2]), f)
    assert np.abs(far - g).max() <= 1e-6 * np.abs(g).max()


def test_uniform_pressure_gives_zero_net_force_and_torque(mesh) -> None:
    v, f = mesh
    force = pressure_force_reference(v, f, DELTA_P)
    scale = np.abs(force).max()
    assert np.abs(force.sum(axis=0)).max() <= 1e-9 * scale
    assert np.abs(np.cross(v, force).sum(axis=0)).max() <= 1e-9 * scale * RADIUS


def test_virial_identity_fixes_the_sign_and_the_factor(mesh) -> None:
    """``sum_i f_i . x_i = 3 delta_p V`` arbitrates both at once; a picture arbitrates neither."""
    v, f = mesh
    force = pressure_force_reference(v, f, DELTA_P)
    target = 3.0 * DELTA_P * mesh_volume(v, f)
    assert abs(float(np.einsum("ij,ij->", force, v)) - target) <= 1e-9 * abs(target)


def test_discrete_volume_is_less_than_its_sphere(mesh) -> None:
    """Substituting the analytic sphere volume produced a phantom nucleus force at t0 once already."""
    v, f = mesh
    assert mesh_volume(v, f) < 4.0 / 3.0 * np.pi * RADIUS ** 3


def test_zero_pressure_is_bit_zero_and_the_law_is_odd(mesh) -> None:
    v, f = mesh
    assert np.array_equal(pressure_force_reference(v, f, 0.0), np.zeros_like(v))
    assert np.array_equal(pressure_force_reference(v, f, -DELTA_P),
                          -pressure_force_reference(v, f, DELTA_P))
