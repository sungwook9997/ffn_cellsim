"""CUDA-only kernel gates for the assembly milestone-2 nucleus + membrane compartments.

The old version launched these production physics kernels on Warp CPU, violating I0-A. They now select the
available CUDA device without an ordinal and skip as a group on the dev Mac. Pure-NumPy acceptance oracles live
in their own tests; every kernel execution happens on the gbook/A5000-or-better path.
"""
from __future__ import annotations

import numpy as np
import pytest

wp = pytest.importorskip("warp")

from aleph.components.incumbent.compartments import (  # noqa: E402
    build_membrane_compartment,
    build_nucleus_compartment,
    mesh_mean_edge,
)

wp.init()
_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)
pytestmark = pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
DEV = str(_CUDA_DEVICE) if _CUDA_DEVICE is not None else ""


def _cortex_shell(n: int = 1500, R: float = 7.5, seed: int = 0) -> np.ndarray:
    """A quasi-uniform actin shell at radius R (the LINC/ERM anchor candidates)."""
    rng = np.random.default_rng(seed)
    p = rng.standard_normal((n, 3))
    p *= R / np.linalg.norm(p, axis=1, keepdims=True)
    return np.ascontiguousarray(p, np.float64)


@pytest.fixture(scope="module")
def built():
    wp.init()
    pa = _cortex_shell()
    n_actin = pa.shape[0]
    nuc, nuc_v, prov = build_nucleus_compartment(pa, n_actin, DEV)
    mem, mem_v = build_membrane_compartment(pa, n_actin + nuc_v.shape[0], DEV, R_mem=7.5)
    pos_all = np.concatenate([pa, nuc_v, mem_v], axis=0)
    return dict(pa=pa, n_actin=n_actin, nuc=nuc, nuc_v=nuc_v, prov=prov, mem=mem, mem_v=mem_v,
               pos_all=pos_all)


# ── build invariants ────────────────────────────────────────────────────────────────────────────
def test_nucleus_build_invariants(built):
    nuc = built["nuc"]
    assert nuc.n_verts == 642 and nuc.n_faces == 1280 and nuc.n_hinges == 1920   # icosphere subdiv-3
    assert nuc.n_hinges == 3 * nuc.n_faces // 2                                   # closed-mesh Euler
    assert nuc.v0_ref > 0 and nuc.kappa_tilde > 0
    assert nuc.n_linc == nuc.n_verts                                             # all nodes paired at gap<reach


def test_membrane_build_invariants(built):
    mem = built["mem"]
    assert mem.n_verts == 642 and mem.n_faces == 1280 and mem.n_hinges == 1920
    assert mem.f_rupt > 0 and mem.kappa_tilde > 0
    # most membrane nodes pair to a cortex node within the ERM reach; on the sparse TEST cortex a few may miss
    # (the native cortex at 494k nodes is dense enough that all pair).
    assert 0.9 * mem.n_verts <= mem.n_erm <= mem.n_verts


def test_global_index_offsets(built):
    """The appended blocks address GLOBAL node indices (nucleus after actin, membrane after nucleus)."""
    n_actin = built["n_actin"]
    nuc, mem = built["nuc"], built["mem"]
    assert nuc.node_off == n_actin
    assert mem.node_off == n_actin + nuc.n_verts
    # face indices point into the block's own node range
    fmin = int(nuc.faces_d.numpy().min())
    fmax = int(nuc.faces_d.numpy().max())
    assert fmin >= n_actin and fmax < n_actin + nuc.n_verts


# ── resting force balance ───────────────────────────────────────────────────────────────────────
def _accumulate(built) -> np.ndarray:
    pos = wp.array(built["pos_all"], dtype=wp.vec3d, device=DEV)
    f = wp.zeros(built["pos_all"].shape[0], dtype=wp.vec3d, device=DEV)
    built["nuc"].accumulate(pos, f)
    built["mem"].accumulate(pos, f)
    wp.synchronize_device(DEV)
    return f.numpy()


def test_resting_forces_finite_and_newton_third_law(built):
    fh = _accumulate(built)
    assert np.all(np.isfinite(fh))
    # every family is an internal action-reaction pair ⇒ the global sum is machine-zero (no force leak)
    assert np.linalg.norm(fh.sum(axis=0)) < 1e-8


def test_resting_tethers_are_force_free_on_cortex(built):
    """At rest LINC/ERM rest = the formation length ⇒ the tethers exert ZERO force on the cortex nodes."""
    fh = _accumulate(built)
    f_cortex = np.linalg.norm(fh[: built["n_actin"]], axis=1)
    assert f_cortex.max() < 1e-9


def test_resting_volume_pressure_is_zero(built):
    """V = V0 at the resting mesh ⇒ the nucleoplasm pressure p_vol = 0 (no spurious inflation)."""
    _accumulate(built)   # runs the device volume→pressure chain
    assert abs(float(built["nuc"].p_vol_d.numpy()[0])) < 1e-6


def test_resting_no_rupture(built):
    nuc = built["nuc"]
    pos = wp.array(built["pos_all"], dtype=wp.vec3d, device=DEV)
    nuc.ruptured_d.zero_()
    nuc.rupture_step(pos, wp.ones(1, dtype=wp.int32, device=DEV))
    wp.synchronize_device(DEV)
    assert int(nuc.ruptured_d.numpy().sum()) == 0


# ── force-family sign / activation under a perturbation ──────────────────────────────────────────
def test_volume_penalty_resists_compression(built):
    """Shrink the nucleus mesh (V<V0) ⇒ p_vol>0 ⇒ vertices pushed OUTWARD (restoring) along +∂V/∂x."""
    nuc = built["nuc"]
    n_actin, nv = built["n_actin"], nuc.n_verts
    verts = built["nuc_v"]
    centre = verts.mean(axis=0)
    shrunk = built["pos_all"].copy()
    shrunk[n_actin:n_actin + nv] = centre + 0.9 * (verts - centre)   # 10% radial shrink
    pos = wp.array(shrunk, dtype=wp.vec3d, device=DEV)
    f = wp.zeros(shrunk.shape[0], dtype=wp.vec3d, device=DEV)
    nuc.accumulate(pos, f)
    wp.synchronize_device(DEV)
    assert float(nuc.p_vol_d.numpy()[0]) > 0.0                       # under-volume ⇒ inflating pressure
    fnuc = f.numpy()[n_actin:n_actin + nv]
    rhat = (shrunk[n_actin:n_actin + nv] - centre)
    rhat /= np.linalg.norm(rhat, axis=1, keepdims=True) + 1e-30
    radial = np.sum(fnuc * rhat, axis=1)
    assert radial.mean() > 0.0                                       # net outward (restoring)


def test_rupture_emerges_past_threshold(built):
    """Inflate a face past eps_rupture ⇒ it tears (ruptured→1); below threshold it does not."""
    nuc = built["nuc"]
    n_actin, nv = built["n_actin"], nuc.n_verts
    verts = built["nuc_v"]
    centre = verts.mean(axis=0)
    # inflate by 60% (areal strain ~1.56 ≫ eps_rupture 0.5) ⇒ every face should tear
    huge = built["pos_all"].copy()
    huge[n_actin:n_actin + nv] = centre + 1.6 * (verts - centre)
    pos = wp.array(huge, dtype=wp.vec3d, device=DEV)
    nuc.ruptured_d.zero_()
    nuc.rupture_step(pos, wp.ones(1, dtype=wp.int32, device=DEV))
    wp.synchronize_device(DEV)
    assert int(nuc.ruptured_d.numpy().sum()) == nuc.n_faces


# ── nucleus mask provider (fluid domain relative-no-flux inclusion) ───────────────────────────────
class _GridStub:
    def __init__(self, shape, dx, origin):
        self.shape = shape
        self.dx = dx
        self.origin = np.asarray(origin, np.float64)
        self.device = DEV


def test_mask_occupancy_matches_sphere_volume(built):
    prov = built["prov"]
    dx = 0.5
    half = 10.5
    nx = int(round(2 * half / dx)) + 1
    grid = _GridStub((nx, nx, nx), dx, (-half, -half, -half))
    inside = prov._inside_grid(grid)
    n_in = int(inside.sum())
    expected = (4.0 / 3.0 * np.pi * 5.0 ** 3) / dx ** 3               # nucleus r_eq=5
    assert 0.9 < n_in / expected < 1.02                              # inscribed icosphere ⇒ slight deficit


def test_mask_only_marks_inside(built):
    """No grid cell OUTSIDE the nucleus is marked inside (winding number is 0 there)."""
    prov = built["prov"]
    dx = 0.5
    half = 10.5
    nx = int(round(2 * half / dx)) + 1
    grid = _GridStub((nx, nx, nx), dx, (-half, -half, -half))
    inside = prov._inside_grid(grid)
    # a corner cell (|x|~half) must be OUTSIDE
    assert inside[0, 0, 0] == 0 and inside[-1, -1, -1] == 0


def test_mesh_mean_edge_positive(built):
    assert mesh_mean_edge(built["nuc_v"], built["nuc"].mesh.faces) > 0
