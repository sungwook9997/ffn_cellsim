"""Arp2/3 cortex sub-population overlap-free build — ``overlap_free=True`` WCA-relaxes the short leaf.

Certifies that ``_build_cortex_arp23_region(..., overlap_free=True)`` (the mixed formin+Arp2/3 cortex path)
resolves the dense sub-σ inter-filament interpenetrations of the short Arp2/3 sub-population — the artifact that
saturates the build-time all-fiber WCA StericForce and kept the mixed cortex from converging in GATE-A — while
preserving the θ₀=70° branch geometry and the rest-0 daughter-base↔branch-vertex anchors.

CUDA-free: the overlap-relax is host NumPy (the weave is pre-CUDA), so the geometric certifications run on the
Mac. A single optional check drives the real Warp ``StericForce`` on the CPU device to report the headline
max|F| the GATE-A gate reads; it is skipped if no Warp device is available.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial import cKDTree

from aleph.components.weave.regions import _build_cortex_arp23_region, cortex_arp23_region

SIGMA_EV_UM = 0.007                                   # F-actin steric diameter (7 nm; assemble.SIGMA_EV_UM)
R_C_UM = 2.0 ** (1.0 / 6.0) * SIGMA_EV_UM             # WCA contact distance
WELD_TOL_UM = 1.0e-6                                  # below this two nodes are the SAME welded junction point


def _build(overlap_free: bool, *, n_arp23: int = 3000, seed: int = 0):
    """Build a representative Arp2/3 cortex leaf (assemble defaults: L=0.15 µm, seg=0.05 µm, mother_frac=0.2)."""
    spec = cortex_arp23_region(
        n_arp23, length_um=0.15, seg_um=0.05, mother_fraction=0.2, R_um=7.5)
    return _build_cortex_arp23_region(spec, np.random.default_rng(seed), overlap_free=overlap_free)


def _cross_fiber_near_contacts(leaf) -> tuple[int, float]:
    """(# of NEAR different-fiber contacts with WELD_TOL<r<r_c, min NON-welded different-fiber distance [µm]).

    Exactly-coincident nodes (r < WELD_TOL) are the intentional welds (daughter-base↔branch-vertex anchor and
    the co-located daughters of one branch rosette); they carry ZERO WCA force (degenerate direction) and are
    excluded. Only NEAR partial overlaps (0<r<r_c) generate steric force — those must be gone after the relax.
    """
    pos = leaf.pos
    off = leaf.fiber_offsets
    fiber_of = np.repeat(np.arange(len(off) - 1), np.diff(off))
    pairs = cKDTree(pos).query_pairs(R_C_UM, output_type="ndarray")
    if pairs.shape[0]:
        pairs = pairs[fiber_of[pairs[:, 0]] != fiber_of[pairs[:, 1]]]
    if pairs.shape[0] == 0:
        return 0, np.inf
    r = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
    near = r > WELD_TOL_UM                              # drop exact welds
    n_near = int(near.sum())
    dmin = float(r[near].min()) if n_near else np.inf
    return n_near, dmin


def _branch_angles_deg(leaf) -> np.ndarray:
    """Signed branch angle at each junction from the (m_after, branch_vertex, daughter_arm) triple [deg]."""
    bt = leaf.branch_triples
    if not bt.size:
        return np.zeros(0)
    v1 = leaf.pos[bt[:, 0]] - leaf.pos[bt[:, 1]]
    v2 = leaf.pos[bt[:, 2]] - leaf.pos[bt[:, 1]]
    c = np.einsum("ij,ij->i", v1, v2) / (
        np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-12)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def test_off_seed_has_dense_near_interpenetration() -> None:
    """Without overlap_free the raw seed carries many NEAR sub-r_c inter-filament contacts (the artifact)."""
    n_near, dmin = _cross_fiber_near_contacts(_build(False))
    assert n_near > 100, f"the raw Arp2/3 seed must show the build interpenetration (got {n_near} near contacts)"
    assert dmin < SIGMA_EV_UM, f"raw seed min cross distance must be sub-σ (got {dmin*1e3:.3f} nm)"


def test_overlap_free_removes_sub_sigma_contacts() -> None:
    """overlap_free=True leaves NO near (0<r<r_c) different-fiber contact ⇒ min non-welded distance ≥ σ."""
    n_near, dmin = _cross_fiber_near_contacts(_build(True))
    assert n_near == 0, f"overlap-free Arp2/3 leaf must have zero near contacts (got {n_near})"
    assert dmin >= SIGMA_EV_UM, f"min non-welded different-fiber distance must be ≥ σ=7 nm (got {dmin*1e3:.3f} nm)"


def test_overlap_free_preserves_branch_geometry() -> None:
    """The θ₀=70° branch angle survives the sub-σ relax (all node moves ≪ 0.5 µm mesh)."""
    ang_off = _branch_angles_deg(_build(False))
    ang_on = _branch_angles_deg(_build(True))
    assert ang_off.size and ang_on.size
    med_off, med_on = float(np.median(ang_off)), float(np.median(ang_on))
    assert abs(med_on - 70.0) < 5.0, f"branch angle median must stay ~70° (got {med_on:.2f}°)"
    assert abs(med_on - med_off) < 2.0, f"relax must not move the branch median (off {med_off:.2f} on {med_on:.2f})"


def test_overlap_free_preserves_rest0_anchors() -> None:
    """Each daughter base stays coincident with its mother branch vertex (rest-0 α-actinin anchor)."""
    leaf = _build(True)
    a = leaf.branch_anchors
    assert a.size, "the Arp2/3 leaf must carry daughter-base↔branch-vertex anchors"
    d = np.linalg.norm(leaf.pos[a[:, 0]] - leaf.pos[a[:, 1]], axis=1)
    assert float(d.max()) < WELD_TOL_UM, f"anchor distance must stay ~0 (got max {d.max()*1e3:.3e} nm)"


def test_overlap_free_zeros_real_wca_steric_force() -> None:
    """The real Warp StericForce kernel: max|F| drops from thousands of pN (OFF) to ~0 (ON). CPU device is fine
    (host reference); skipped only if no Warp device can be created at all."""
    wp = pytest.importorskip("warp")
    from types import SimpleNamespace

    from aleph.components.solid.steric_warp import StericForce

    wp.init()
    dev = "cuda" if any(d.is_cuda for d in wp.get_devices()) else "cpu"

    def steric_max(leaf) -> float:
        fid = np.repeat(np.arange(leaf.fiber_offsets.size - 1), np.diff(leaf.fiber_offsets)).astype(np.int32)
        st = StericForce(fiber_id=fid, sigma=SIGMA_EV_UM, k_ev=1.0e5, device=dev, force_cap=1.0e3)
        with wp.ScopedDevice(dev):
            pos_d = wp.array(np.ascontiguousarray(leaf.pos, np.float64), dtype=wp.vec3d, device=dev)
            f_d = wp.zeros(leaf.pos.shape[0], dtype=wp.vec3d, device=dev)
            st.accumulate(SimpleNamespace(pos=pos_d), f_d)
            return float(np.linalg.norm(f_d.numpy(), axis=1).max())

    off_max = steric_max(_build(False))
    on_max = steric_max(_build(True))
    assert off_max > 100.0, f"the raw seed must saturate the WCA steric force (got {off_max:.1f} pN)"
    assert on_max < 1.0e-6, f"overlap-free Arp2/3 leaf must start with ~zero steric force (got {on_max} pN)"
