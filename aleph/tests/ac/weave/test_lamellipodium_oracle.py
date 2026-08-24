"""Self-test of the dendritic Arp2/3 lamellipodium rebuild (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.emergence.nematic import nematic_order
from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD, branch_angle, thermal_sigma
from aleph.components.weave.lamellipodium import build_lamellipodium


def _measure_branch_angles(seed) -> np.ndarray:
    """Label-blind junction angles measured from geometry (branch vertex between mother and daughter arms)."""
    tr = seed.branch_triples
    return np.array([branch_angle(seed.pos[m], seed.pos[j], seed.pos[d]) for m, j, d in tr])


def test_dendritic_topology_is_a_tree() -> None:
    """Every daughter has exactly one parent branch; triples/anchors reference valid nodes; N stays fixed."""
    rng = np.random.default_rng(0)
    seed = build_lamellipodium(n_mothers=8, n_daughter_pool=40, length_um=1.0, seg_um=0.5,
                               patch_half_um=8.0, rng=rng, monomer_c=50.0, k_arp0=5.0, k_m=5.0, tau=1.0)
    n_fib = seed.fiber_offsets.shape[0] - 1
    assert n_fib == 8 + 40                              # N fixed (mothers + dormant pool), no new nodes
    assert seed.branch_triples.shape[0] == 40          # one branch per daughter (a tree)
    assert seed.branch_anchors.shape[0] == 40
    assert seed.branch_triples.max() < seed.pos.shape[0]
    # each daughter appears as the daughter-node of exactly one triple
    daughters = seed.branch_triples[:, 2]
    assert len(np.unique(daughters)) == 40


def test_branch_angle_distribution_matches_theta0_sigma() -> None:
    """The MEASURED junction angles match theta0 +- sigma_theta (the §3A.b branch-angle gate)."""
    rng = np.random.default_rng(1)
    seed = build_lamellipodium(n_mothers=20, n_daughter_pool=3000, length_um=1.0, seg_um=0.5,
                               patch_half_um=8.0, rng=rng, monomer_c=1e6, k_arp0=50.0, k_m=5.0, tau=5.0)
    ang = _measure_branch_angles(seed)
    sigma = thermal_sigma(ARP23_K_THETA)
    assert np.mean(ang) == pytest.approx(ARP23_THETA0_RAD, abs=0.03)   # ~70 deg (small sin-measure shift up)
    assert np.mean(ang) > ARP23_THETA0_RAD                            # sin-measure shifts toward pi/2
    assert np.std(ang) == pytest.approx(sigma, rel=0.15)              # ~9 deg thermal spread, not clamped


def test_flux_limit_activates_more_with_more_monomer() -> None:
    """A rich monomer field activates MORE dormant daughters than a depleted one (N-fixed; same geometry)."""
    lo = build_lamellipodium(n_mothers=10, n_daughter_pool=400, length_um=1.0, seg_um=0.5, patch_half_um=8.0,
                             rng=np.random.default_rng(7), monomer_c=0.5, k_arp0=5.0, k_m=5.0, tau=0.3)
    hi = build_lamellipodium(n_mothers=10, n_daughter_pool=400, length_um=1.0, seg_um=0.5, patch_half_um=8.0,
                             rng=np.random.default_rng(7), monomer_c=50.0, k_arp0=5.0, k_m=5.0, tau=0.3)
    n_active_lo = int(lo.active_mask[10:].sum())       # active daughters (mothers are [:10])
    n_active_hi = int(hi.active_mask[10:].sum())
    assert n_active_hi > n_active_lo
    assert lo.active_mask[:10].all() and hi.active_mask[:10].all()    # mothers always active


def test_isotropic_seed_is_not_pre_aligned() -> None:
    """With a UNIFORM (kappa=0) NPF field the mother orientation order S is low (the two-mode is NOT scripted);
    a strong NPF bias raises S. The branch topology is dendritic either way."""
    iso = build_lamellipodium(n_mothers=300, n_daughter_pool=0, length_um=1.0, seg_um=0.5, patch_half_um=8.0,
                              rng=np.random.default_rng(2), npf_kappa=0.0)
    aln = build_lamellipodium(n_mothers=300, n_daughter_pool=0, length_um=1.0, seg_um=0.5, patch_half_um=8.0,
                              rng=np.random.default_rng(2), npf_kappa=6.0)
    # mother axes as a fiber array: 2 nodes per fiber is enough for an end-to-end axis; reuse the nematic S
    s_iso = nematic_order(iso.mother_axes)
    s_aln = nematic_order(aln.mother_axes)
    assert s_iso < 0.6                                  # not a scripted aligned bundle (2D-isotropic baseline ~0.25)
    assert s_aln > s_iso + 0.15                         # the NPF bias raises order (a declared, ablatable seed)
