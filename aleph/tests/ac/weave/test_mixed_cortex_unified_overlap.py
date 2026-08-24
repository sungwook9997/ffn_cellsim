"""Unified whole-cortex cross-region overlap relax — mixed formin + Arp2/3 cortex converges (2026-07-24).

The mixed cortex (``cortex_arp23_fraction > 0``) is woven as TWO separately-relaxed leaves: the long-formin
cortex (relaxed by ``ff.weave.weave(..., overlap_free=True)``) and the short branched Arp2/3 leaf (relaxed
within-leaf by ``regions._relax_arp23_overlaps``). Co-placed on the SAME cortex shell, they were never relaxed
AGAINST each other, so ~arp23↔formin CROSS-REGION sub-σ interpenetrations survived and saturated the build-time
all-fiber WCA StericForce — native GATE-A started at max|PF|≈2233 (vs the formin-only ~0.14) and the weak pathA
descent diverged.

``weave_cell(..., overlap_free=True)`` now runs a UNIFIED WCA push-apart over ALL different-fiber node pairs of
the COMBINED cortex (formin↔arp23 + arp23↔arp23 together; formin↔formin already clean), holding every
daughter-base↔branch-vertex weld EXCLUDED from the push + re-welded each sweep so the θ₀=70° branch geometry and
the rest-0 anchors survive. It fires ONLY when an Arp2/3 (``sphere_dendritic``) population is present, so the
formin-only cell is untouched (Gate-1 byte-identical).

This test builds a representative mixed cortex (frac=0.33) and certifies:
  * BEFORE (the two per-leaf-relaxed leaves merely CONCATENATED, the pre-fix state) carries cross-region sub-σ
    contacts and a saturated (>100 pN) real Warp StericForce;
  * AFTER (``weave_cell(overlap_free=True)``, the unified relax) has NO different-fiber sub-σ contact across
    formin↔arp23 AND arp23↔arp23, the θ₀=70° branch median survives, the anchors stay rest-0, and the real Warp
    StericForce max|F| drops to ~0 (<1 pN).

CUDA-free: the relax is host NumPy (pre-CUDA weave), so the geometric certifications run on the Mac; the real
Warp ``StericForce`` check runs on whatever Warp device exists (CUDA if present, else the CPU host reference —
it is the same kernel), and is skipped only if Warp cannot be imported at all.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial import cKDTree

from aleph.components.incumbent.assemble import R_CORTEX_UM, _cortex_region
from aleph.components.weave.regions import build_region, cortex_arp23_region, derive_cortex_arp23_split
from aleph.components.weave.woven_cell import weave_cell

SIGMA_EV_UM = 0.007                                  # F-actin steric diameter (7 nm; assemble.SIGMA_EV_UM)
R_C_UM = 2.0 ** (1.0 / 6.0) * SIGMA_EV_UM            # WCA contact distance
WELD_TOL_UM = 1.0e-6                                 # below this two nodes are the SAME welded junction point

_N_FILAMENTS = 10_000                                # representative mixed-cortex budget (native is 70,686)
_ARP23_FRACTION = 0.33                               # Bovellan 2014 Arp2/3 mass fraction (assemble default)
_FORMIN_LEN, _ARP23_LEN = 3.0, 0.15                  # long-formin / short-Arp2/3 contour lengths [µm]
_SEED = 0


def _mixed_cortex_regions():
    """The two-leaf mixed cortex region list (long formin + short branched Arp2/3), assemble defaults."""
    n_formin, n_arp23 = derive_cortex_arp23_split(_N_FILAMENTS, _ARP23_FRACTION, _FORMIN_LEN, _ARP23_LEN)
    formin = _cortex_region(n_formin, seg_um=0.5, length_um=_FORMIN_LEN, density_per_fil=20.0)
    arp23 = cortex_arp23_region(
        n_arp23, length_um=_ARP23_LEN, seg_um=0.05, mother_fraction=0.2, R_um=R_CORTEX_UM)
    return [formin, arp23], n_formin, n_arp23


def _combined_before():
    """Pre-fix state: each leaf relaxed WITHIN itself (overlap_free=True) then merely CONCATENATED — no
    relax against each other. Returns ``(pos, node_fiber_global, region_of_node)``."""
    regions, _, _ = _mixed_cortex_regions()
    rng = np.random.default_rng(_SEED)
    leaves = [build_region(spec, rng, overlap_free=True) for spec in regions]
    pos, node_fiber, region = [], [], []
    fid = 0
    for rid, leaf in enumerate(leaves):
        pos.append(leaf.pos)
        node_fiber.append(np.repeat(np.arange(leaf.fiber_offsets.size - 1) + fid, np.diff(leaf.fiber_offsets)))
        region.append(np.full(leaf.pos.shape[0], rid, np.int64))
        fid += leaf.fiber_offsets.size - 1
    return np.concatenate(pos), np.concatenate(node_fiber).astype(np.int64), np.concatenate(region)


def _combined_after():
    """Post-fix state: ``weave_cell(overlap_free=True)`` runs the unified cross-region relax on the combined
    cortex. Returns the assembled :class:`WovenCell` (same seed/order as ``_combined_before``)."""
    regions, _, _ = _mixed_cortex_regions()
    return weave_cell(regions, rng=np.random.default_rng(_SEED), overlap_free=True)


def _cross_fiber_near(pos, node_fiber, region_of_node=None, require_cross_region=False):
    """(# NEAR different-fiber contacts WELD_TOL<r<r_c, min non-welded different-fiber distance [µm]).

    Exactly-coincident nodes (r < WELD_TOL) are the intentional welds (zero WCA force) and are excluded — only
    NEAR partial overlaps (0<r<r_c) generate steric force. With ``require_cross_region`` the pair must also span
    the two regions (formin↔arp23), isolating the cross-region interpenetration the unified relax targets."""
    pairs = cKDTree(pos).query_pairs(R_C_UM, output_type="ndarray")
    if pairs.shape[0]:
        pairs = pairs[node_fiber[pairs[:, 0]] != node_fiber[pairs[:, 1]]]       # different fibers only
    if pairs.shape[0] and require_cross_region and region_of_node is not None:
        pairs = pairs[region_of_node[pairs[:, 0]] != region_of_node[pairs[:, 1]]]
    if pairs.shape[0] == 0:
        return 0, np.inf
    r = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
    near = r > WELD_TOL_UM
    return int(near.sum()), (float(r[near].min()) if near.any() else np.inf)


def _branch_angles_deg(pos, branch_triples):
    """Branch angle at each junction from the (m_after, branch_vertex, daughter_arm) triple [deg]."""
    if not branch_triples.size:
        return np.zeros(0)
    v1 = pos[branch_triples[:, 0]] - pos[branch_triples[:, 1]]
    v2 = pos[branch_triples[:, 2]] - pos[branch_triples[:, 1]]
    c = np.einsum("ij,ij->i", v1, v2) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1) + 1e-12)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def _steric_max_force(pos, node_fiber, device) -> float:
    """max|F| of the real Warp ``StericForce`` kernel over the combined cortex (GATE-A build-time headline)."""
    from types import SimpleNamespace

    import warp as wp

    from aleph.components.solid.steric_warp import StericForce

    st = StericForce(fiber_id=node_fiber.astype(np.int32), sigma=SIGMA_EV_UM, k_ev=1.0e5,
                     device=device, force_cap=1.0e3)
    with wp.ScopedDevice(device):
        pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=device)
        f_d = wp.zeros(pos.shape[0], dtype=wp.vec3d, device=device)
        st.accumulate(SimpleNamespace(pos=pos_d), f_d)
        return float(np.linalg.norm(f_d.numpy(), axis=1).max())


def test_before_has_cross_region_interpenetration() -> None:
    """The pre-fix concatenation (leaves relaxed only WITHIN themselves) still carries cross-region contacts."""
    pos, node_fiber, region = _combined_before()
    n_cross, dmin_cross = _cross_fiber_near(pos, node_fiber, region, require_cross_region=True)
    assert n_cross >= 1, "the co-placed per-leaf-relaxed cortex must show ≥1 formin↔arp23 cross-region contact"
    assert dmin_cross < R_C_UM, f"cross-region contact must be sub-r_c (force-generating) (got {dmin_cross * 1e3:.3f} nm)"


def test_unified_relax_removes_all_sub_sigma_contacts() -> None:
    """AFTER the unified relax: NO different-fiber sub-r_c contact survives across formin↔arp23 AND arp23↔arp23
    ⇒ min non-welded different-fiber distance ≥ σ over the WHOLE combined cortex."""
    wc = _combined_after()
    region_of_node = wc.region_id[wc.node_fiber]
    n_all, dmin_all = _cross_fiber_near(wc.pos, wc.node_fiber)
    n_cross, _ = _cross_fiber_near(wc.pos, wc.node_fiber, region_of_node, require_cross_region=True)
    assert n_cross == 0, f"no formin↔arp23 cross-region sub-σ contact may survive (got {n_cross})"
    assert n_all == 0, f"no different-fiber sub-σ contact may survive anywhere in the cortex (got {n_all})"
    assert dmin_all >= SIGMA_EV_UM, f"min non-welded different-fiber distance must be ≥ σ (got {dmin_all * 1e3:.3f} nm)"


def test_unified_relax_preserves_branch_geometry_and_anchors() -> None:
    """The θ₀=70° branch median survives the sub-σ relax and every daughter base stays welded to its vertex."""
    wc = _combined_after()
    ang = _branch_angles_deg(wc.pos, wc.branch_triples)
    assert ang.size, "the mixed cortex must carry Arp2/3 branch triples"
    med = float(np.median(ang))
    assert abs(med - 70.0) < 5.0, f"branch-angle median must stay ~70° after the unified relax (got {med:.2f}°)"

    a = wc.branch_anchors
    assert a.size, "the mixed cortex must carry daughter-base↔branch-vertex anchors"
    d = np.linalg.norm(wc.pos[a[:, 0]] - wc.pos[a[:, 1]], axis=1)
    assert float(d.max()) < WELD_TOL_UM, f"anchors must stay rest-0 (got max {d.max() * 1e3:.3e} nm)"


def test_unified_relax_zeros_real_wca_steric_force() -> None:
    """The real Warp StericForce over the COMBINED cortex: max|F| drops from thousands of pN (before) to ~0 pN
    (after). Runs on the Warp CPU host reference if no CUDA is present (same kernel); skipped only if Warp is
    unimportable."""
    wp = pytest.importorskip("warp")
    wp.init()
    device = "cuda" if any(d.is_cuda for d in wp.get_devices()) else "cpu"

    pos_b, node_fiber_b, _ = _combined_before()
    before_max = _steric_max_force(pos_b, node_fiber_b, device)

    wc = _combined_after()
    after_max = _steric_max_force(wc.pos, wc.node_fiber, device)

    assert before_max > 50.0, f"the co-placed unrelaxed cortex must carry a force-generating interpenetration (got {before_max:.1f} pN)"
    assert after_max < 1.0, f"the unified-relaxed combined cortex must start with ~zero steric force (got {after_max} pN)"


def test_formin_only_cortex_is_byte_identical() -> None:
    """The guard is exact: with NO Arp2/3 population the unified relax never fires, so the formin-only cortex is
    byte-identical whether or not the mixed-cortex path exists (Gate-1 protection)."""
    formin = _cortex_region(2000, seg_um=0.5, length_um=_FORMIN_LEN, density_per_fil=20.0)
    a = weave_cell([formin], rng=np.random.default_rng(_SEED), overlap_free=True)
    b = weave_cell([formin], rng=np.random.default_rng(_SEED), overlap_free=True)
    assert np.array_equal(a.pos, b.pos)
    assert np.array_equal(a.xl_rest, b.xl_rest)
    # a formin-only cortex has no sphere_dendritic leaf ⇒ no branch anchors ⇒ the relax branch is skipped
    assert a.branch_anchors.size == 0
