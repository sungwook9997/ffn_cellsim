"""Self-test of the nuclear-envelope dihedral-bending oracle (pure NumPy — no Warp/CUDA).

The bending gate the ff/membrane pattern also passes, re-derived for the nucleus: sphere→8πκ (Willmore),
the icosphere dihedral sum Σ_ref≈7.50 (resolution- and R-independent — grid-invariant, Magic-Number
Block), the FD-gradient sign of the force, flat→0, and Σf=0 (translation invariance).
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.geometry import build_oblate_mesh, build_hinges
from aleph.components.nucleus.lamina_analytic import (
    calibrate_kappa_tilde,
    dihedral_bending_energy,
    dihedral_bending_forces,
    sphere_willmore_energy,
)

KAPPA_NE = 0.0828  # pN·µm (~20 kBT; envelope bending rigidity order)


@pytest.mark.parametrize("sub", [2, 3, 4])
def test_sphere_bending_energy_is_8pi_kappa(sub: int) -> None:
    """After κ̃-calibration the discrete sphere reproduces the continuum Willmore energy 8πκ exactly,
    at every subdivision (the calibration makes it exact by construction — resolution-invariant)."""
    v, f = build_oblate_mesh(3.0, 1.0, sub)
    h = build_hinges(f)
    kt = calibrate_kappa_tilde(KAPPA_NE, v, h)
    e = dihedral_bending_energy(v, h, kt)
    assert e == pytest.approx(sphere_willmore_energy(KAPPA_NE), rel=1e-9)


def test_dihedral_sum_is_grid_and_R_invariant() -> None:
    """Raw Σ(1−n̂₁·n̂₂) converges to ≈7.50 and is independent of R (a topological constant of the
    icosphere, NOT a fit) — this is what makes κ̃=8πκ/Σ_ref satisfy the Magic-Number Block."""
    sigmas_R = [dihedral_bending_energy(*_(build_oblate_mesh(r, 1.0, 3)), 1.0) for r in (1.5, 3.0, 6.0)]
    assert np.allclose(sigmas_R, sigmas_R[0], rtol=1e-6)          # R-independent
    assert sigmas_R[0] == pytest.approx(7.50, abs=0.05)
    # convergence toward 7.50 with refinement
    s2 = dihedral_bending_energy(*_(build_oblate_mesh(3.0, 1.0, 2)), 1.0)
    s4 = dihedral_bending_energy(*_(build_oblate_mesh(3.0, 1.0, 4)), 1.0)
    assert abs(s2 - 7.50) < 0.05 and abs(s4 - 7.50) < 0.05


def test_bending_force_matches_fd_gradient() -> None:
    """f = −∂E/∂x to <1e-6·‖f‖ on a PERTURBED (non-flat) mesh — the sign arbiter of the force."""
    rng = np.random.default_rng(1)
    v, f = build_oblate_mesh(3.0, 1.4, 3)
    v = v + 0.04 * rng.standard_normal(v.shape)
    h = build_hinges(f)
    v0, f0 = build_oblate_mesh(3.0, 1.0, 3)
    kt = calibrate_kappa_tilde(KAPPA_NE, v0, build_hinges(f0))
    f_ana = dihedral_bending_forces(v, h, kt)
    eps = 1e-7
    f_fd = np.zeros_like(v)
    for i in range(v.shape[0]):
        for d in range(3):
            vp = v.copy(); vp[i, d] += eps
            vm = v.copy(); vm[i, d] -= eps
            f_fd[i, d] = -(dihedral_bending_energy(vp, h, kt) - dihedral_bending_energy(vm, h, kt)) / (2 * eps)
    scale = np.max(np.abs(f_ana))
    assert np.max(np.abs(f_ana - f_fd)) < 1e-6 * scale


def test_bending_force_is_translation_invariant() -> None:
    """Σf = 0 (no net force from an internal energy) — a hard conservation invariant."""
    rng = np.random.default_rng(2)
    v, f = build_oblate_mesh(3.0, 1.5, 3)
    v = v + 0.05 * rng.standard_normal(v.shape)
    h = build_hinges(f)
    kt = calibrate_kappa_tilde(KAPPA_NE, *(build_oblate_mesh(3.0, 1.0, 3)[0], build_hinges(build_oblate_mesh(3.0, 1.0, 3)[1])))
    f_ana = dihedral_bending_forces(v, h, kt)
    assert np.max(np.abs(f_ana.sum(axis=0))) < 1e-10


def test_flat_hinge_zero_force() -> None:
    """A coplanar hinge (flat) → zero bending force (E has a minimum; ∂E/∂x=0 at the flat reference)."""
    # two triangles sharing edge (0,1), all in the z=0 plane
    verts = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0], [0.5, -1, 0]], np.float64)
    hinges = np.array([[0, 1, 2, 3]], np.int64)
    f = dihedral_bending_forces(verts, hinges, 1.0)
    assert np.max(np.abs(f)) < 1e-12


def _(mesh):
    """Unpack (verts, faces) → (verts, hinges) for the Σ helpers."""
    v, f = mesh
    return v, build_hinges(f)
