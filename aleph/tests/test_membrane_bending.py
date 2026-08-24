"""Helfrich membrane-bending kernel — MANDATORY acceptance gate (Sanity Gate Protocol).

The lipid membrane's DEFINING mechanic is bending; the production cell had only a lumped 2γ/R tension (no
bending anywhere). These tests validate the new Seung–Nelson dihedral-hinge kernel
(`ff.membrane_surface.helfrich_bending_kernel`) BEFORE it is wired into the whole-cell run:

  1. FINITE-DIFFERENCE gradient — the kernel force equals −∂E/∂x of E=Σκ̃(1−n̂₁·n̂₂). This is the ARBITER of
     the force sign/structure (the one thing the analytic formula can get wrong).
  2. Sphere energy → 8πκ_m, R-INDEPENDENT — catches winding/normal-orientation bugs (a flipped normal spikes
     the energy on the affected hinges).
  3. Flat patch → zero force — a planar sheet has no curvature ⇒ no bending force.
  4. Displaced vertex → a FINITE-WIDTH dimple, not a spike — the physical acceptance the lumped tension can't
     give (bending distributes the curvature over neighbours).

Run:  conda activate ffn_sim && python -m pytest aleph/tests/test_membrane_bending.py -q
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.laws.membrane_surface import (
    build_membrane_mesh, build_membrane_hinges, helfrich_bending_kernel,
    membrane_bending_energy, kappa_tilde, calibrate_kappa_tilde, SIGMA_ICOSPHERE,
)

KAPPA_M = 0.0828   # pN·µm (20 kBT, Rawicz 2000; KB-3.B1.2)


def _bending_force(verts: np.ndarray, hinges: np.ndarray, kt: float, device: str = "cpu") -> np.ndarray:
    """Launch the Warp kernel and return the (N,3) bending force at ``verts``."""
    N = verts.shape[0]
    pos_d = wp.array(np.ascontiguousarray(verts, np.float64), dtype=wp.vec3d, device=device)
    hin_d = wp.array(np.ascontiguousarray(hinges, np.int32), dtype=wp.int32, ndim=2, device=device)
    f_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if hinges.shape[0] > 0:
        wp.launch(helfrich_bending_kernel, dim=hinges.shape[0],
                  inputs=[pos_d, hin_d, wp.float64(kt), f_d], device=device)
        wp.synchronize_device(device)
    return f_d.numpy().astype(np.float64)


def test_finite_difference_gradient():
    """The kernel force = −∂E/∂x (central differences). Sign/structure arbiter."""
    m = build_membrane_mesh(3.0, subdivisions=2)
    hinges = build_membrane_hinges(m.faces)
    rng = np.random.default_rng(0)
    verts = m.verts + 0.05 * rng.standard_normal(m.verts.shape)   # break symmetry → robust nonzero forces
    kt = kappa_tilde(KAPPA_M)

    f_kernel = _bending_force(verts, hinges, kt)
    # central-difference −∂E/∂x on a sample of vertices (all of them; mesh is small)
    eps = 1e-6
    f_fd = np.zeros_like(verts)
    for i in range(verts.shape[0]):
        for a in range(3):
            vp = verts.copy(); vp[i, a] += eps
            vm = verts.copy(); vm[i, a] -= eps
            Ep = membrane_bending_energy(vp, hinges, kt)
            Em = membrane_bending_energy(vm, hinges, kt)
            f_fd[i, a] = -(Ep - Em) / (2.0 * eps)

    scale = np.abs(f_kernel).max()
    err = np.abs(f_kernel - f_fd).max()
    assert scale > 1e-6, f"forces vanished (scale={scale}); test would be vacuous"
    assert err < 1e-4 * scale, f"kernel force != −∂E/∂x: max |Δ|={err:.3e}, scale={scale:.3e}"


def test_translation_invariance():
    """Σf = 0 on every hinge ⇒ the total bending force is zero (no spurious net force / drift)."""
    m = build_membrane_mesh(4.0, subdivisions=2)
    hinges = build_membrane_hinges(m.faces)
    rng = np.random.default_rng(1)
    verts = m.verts + 0.05 * rng.standard_normal(m.verts.shape)
    f = _bending_force(verts, hinges, kappa_tilde(KAPPA_M))
    net = np.abs(f.sum(axis=0)).max()
    assert net < 1e-9, f"net bending force should vanish (translation invariance); got {net:.3e}"


def test_sphere_energy_8pi_kappa():
    """Closed sphere bending energy → 8πκ_m, R-independent (winding-bug catcher). Uses the CALIBRATED κ̃
    (=8πκ_m/Σ_ico); the SN 2κ/√3 factor is wrong for the icosphere (see the module header). The raw dihedral
    sum Σ(1−cosθ) is the geometric constant SIGMA_ICOSPHERE≈7.50 — checked here (that IS the calibration)."""
    target = 8.0 * np.pi * KAPPA_M
    for R in (3.0, 6.0):
        for sub in (3, 4):
            m = build_membrane_mesh(R, subdivisions=sub)
            hinges = build_membrane_hinges(m.faces)
            # Euler check: closed sphere has 3·Ntri/2 interior edges (all shared)
            assert hinges.shape[0] == 3 * m.faces.shape[0] // 2, "hinge count != 3·Ntri/2 (non-closed?)"
            # the raw dihedral sum is the icosphere geometric constant (resolution/R-independent)
            sigma = membrane_bending_energy(m.verts, hinges, 1.0)
            assert abs(sigma - SIGMA_ICOSPHERE) < 0.05 * SIGMA_ICOSPHERE, \
                f"icosphere Σ(1−cos)={sigma:.4f} vs constant {SIGMA_ICOSPHERE} (winding/topology bug?)"
            # with the DEFAULT (constant) calibration the sphere energy is 8πκ to the Σ spread (~0.3%)
            E_const = membrane_bending_energy(m.verts, hinges, kappa_tilde(KAPPA_M))
            assert abs(E_const - target) < 0.05 * target, f"sphere E(R={R},sub={sub})={E_const:.4f} vs 8πκ={target:.4f}"
            # per-mesh auto-calibration is EXACT
            E_exact = membrane_bending_energy(m.verts, hinges, calibrate_kappa_tilde(KAPPA_M, m.verts, hinges))
            assert abs(E_exact - target) < 1e-9 * target, "calibrate_kappa_tilde should give exactly 8πκ"


def test_flat_patch_zero_force():
    """A planar sheet has no curvature ⇒ zero bending force (to round-off)."""
    # a small flat triangulated square in the z=0 plane
    xs = np.linspace(-1.0, 1.0, 5)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    verts = np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], axis=1).astype(np.float64)
    faces = []
    n = 5
    for r in range(n - 1):
        for c in range(n - 1):
            a = r * n + c; b = r * n + c + 1; cc = (r + 1) * n + c; d = (r + 1) * n + c + 1
            faces.append((a, b, d)); faces.append((a, d, cc))    # consistent winding
    faces = np.asarray(faces, np.int64)
    hinges = build_membrane_hinges(faces)
    f = _bending_force(verts, hinges, kappa_tilde(KAPPA_M))
    assert np.abs(f).max() < 1e-9, f"flat patch bending force should be ~0; got {np.abs(f).max():.3e}"


def test_displaced_vertex_dimple_not_spike():
    """Pull one vertex outward on a sphere; bending resists AND spreads the curvature to a finite width
    (a dimple), not a single-node spike. Acceptance: the perturbed vertex feels a RESTORING (inward-tangent-to-
    the-radial) force, and its ring-1 neighbours feel a comparable force (non-local), unlike a spike."""
    m = build_membrane_mesh(3.0, subdivisions=3)
    hinges = build_membrane_hinges(m.faces)
    verts = m.verts.copy()
    centre = verts.mean(0)
    p = 0                                                    # perturb vertex 0 radially outward by 0.3 µm
    r_hat = (verts[p] - centre) / np.linalg.norm(verts[p] - centre)
    verts[p] = verts[p] + 0.3 * r_hat
    f = _bending_force(verts, hinges, kappa_tilde(KAPPA_M))

    # restoring: the bending force on the pulled vertex opposes the outward displacement (negative radial comp.)
    f_radial = float(np.dot(f[p], r_hat))
    assert f_radial < 0.0, f"bending should pull the bumped vertex back inward; f·r̂={f_radial:.3e}"
    # non-local (finite width): neighbours of p carry appreciable force too (not concentrated on p alone)
    neigh = set()
    for h in hinges:
        s = {int(h[0]), int(h[1]), int(h[2]), int(h[3])}
        if p in s:
            neigh |= s
    neigh.discard(p)
    neigh = np.array(sorted(neigh))
    fmag = np.linalg.norm(f, axis=1)
    assert fmag[neigh].max() > 0.1 * fmag[p], "curvature should spread to neighbours (dimple, not a spike)"


if __name__ == "__main__":
    test_finite_difference_gradient()
    test_translation_invariance()
    test_sphere_energy_8pi_kappa()
    test_flat_patch_zero_force()
    test_displaced_vertex_dimple_not_spike()
    print("membrane bending: ALL PASS")
