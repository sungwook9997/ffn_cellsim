"""Self-test of the nucleoplasm-volume oracle (pure NumPy — no Warp/CUDA).

Volume conservation is the incompressible-nucleoplasm (ν→½) invariant: the divergence-theorem signed
volume matches the continuum sphere/oblate, its exact vertex gradient matches a finite difference (the
volume-constraint force direction), and the closed-mesh volume is translation-invariant.
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.geometry import (
    build_oblate_mesh,
    mesh_volume,
    mesh_volume_gradient,
    oblate_volume,
)


@pytest.mark.parametrize("sub, tol", [(3, 0.02), (4, 0.005), (5, 0.0015)])
def test_sphere_volume_converges_to_continuum(sub: int, tol: float) -> None:
    """Signed-tet mesh volume → (4/3)πR³ as the icosphere refines (discretization deficit shrinks)."""
    v, f = build_oblate_mesh(3.0, 1.0, sub)
    v_cont = 4.0 / 3.0 * np.pi * 3.0 ** 3
    assert mesh_volume(v, f) == pytest.approx(v_cont, rel=tol)


def test_oblate_volume_matches_continuum() -> None:
    """Oblate mesh volume ≈ continuum (4/3)π a² c for a range of aspects (volume-conservation basis)."""
    for aspect in (1.5, 2.0, 3.0):
        v, f = build_oblate_mesh(3.0, aspect, 4)
        assert mesh_volume(v, f) == pytest.approx(oblate_volume(3.0, aspect), rel=0.01)


def test_volume_gradient_matches_fd() -> None:
    """∂V/∂x_i (the divergence-theorem gradient = the volume-constraint force direction) matches FD."""
    rng = np.random.default_rng(3)
    v, f = build_oblate_mesh(3.0, 1.6, 3)
    v = v + 0.05 * rng.standard_normal(v.shape)
    g_ana = mesh_volume_gradient(v, f)
    eps = 1e-7
    g_fd = np.zeros_like(v)
    for i in range(v.shape[0]):
        for d in range(3):
            vp = v.copy(); vp[i, d] += eps
            vm = v.copy(); vm[i, d] -= eps
            g_fd[i, d] = (mesh_volume(vp, f) - mesh_volume(vm, f)) / (2 * eps)
    assert np.max(np.abs(g_ana - g_fd)) < 1e-6 * np.max(np.abs(g_ana))


def test_volume_translation_invariant() -> None:
    """Closed-mesh volume is unchanged by a rigid translation (the shift terms telescope to 0)."""
    v, f = build_oblate_mesh(3.0, 2.0, 3)
    v0 = mesh_volume(v, f)
    v_shift = mesh_volume(v + np.array([12.3, -4.5, 7.7]), f)
    assert v_shift == pytest.approx(v0, rel=1e-10)
