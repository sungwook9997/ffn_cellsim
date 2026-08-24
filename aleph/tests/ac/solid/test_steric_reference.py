"""Self-test of the excluded-volume ALGORITHM oracle (pure NumPy — no Warp/CUDA).

The Warp hash-grid ``steric_warp.StericForce`` must reproduce these host references. Here we certify the
references against each other and against the gate contract: hash-grid (cell-list) == brute-force
(grid-invariance of the neighbour accelerator = the Magic-Number-Block statement for k_EV's grid), OFF /
zero-overlap == zero-force, same-filament exclusion, Newton's third law, a compressed patch resisting
collapse, and the CFL k_EV bookkeeping being derived (not tuned to a crowding outcome).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.solid.steric_reference import (
    cfl_stiffness_term,
    steric_force_bruteforce,
    steric_force_celllist,
    total_steric_energy,
)
from aleph.components.solid.wca_analytic import (
    compressed_pair_separation,
    contact_stiffness,
    epsilon_from_contact_stiffness,
    wca_cutoff,
)

SIGMA = 0.05
K_EV = 200.0
EPS = epsilon_from_contact_stiffness(K_EV, SIGMA)


def _crowded_cluster(seed: int = 0, n_fiber: int = 6, per_fiber: int = 8):
    """A small crowded multi-filament node cloud in a ~4 r_c box so many cross-fiber pairs overlap."""
    rng = np.random.default_rng(seed)
    r_c = wca_cutoff(SIGMA)
    box = 4.0 * r_c
    pos = rng.uniform(0.0, box, size=(n_fiber * per_fiber, 3))
    fiber_id = np.repeat(np.arange(n_fiber), per_fiber).astype(np.int64)
    return pos, fiber_id


# --------------------------------------------------------------------------------------------------
# Gate 3a — hash-grid (cell-list) == brute-force : grid-invariance of the neighbour accelerator
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("cell_mult", [1.0, 1.5, 2.0, 3.0])
def test_celllist_matches_bruteforce(cell_mult: float) -> None:
    """The cell-list (CPU analogue of the device HashGrid) equals the O(N^2) sum for any cell >= r_c."""
    pos, fiber_id = _crowded_cluster()
    r_c = wca_cutoff(SIGMA)
    f_brute = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    f_cell = steric_force_celllist(pos, SIGMA, EPS, fiber_id, cell_size=cell_mult * r_c)
    assert np.allclose(f_brute, f_cell, rtol=1e-12, atol=1e-12)


def test_celllist_rejects_undersize_cell() -> None:
    """A cell smaller than r_c would miss within-cutoff pairs -> refused (not silently wrong)."""
    pos, fiber_id = _crowded_cluster()
    with pytest.raises(ValueError):
        steric_force_celllist(pos, SIGMA, EPS, fiber_id, cell_size=0.5 * wca_cutoff(SIGMA))


# --------------------------------------------------------------------------------------------------
# Gate 2 — OFF / zero-overlap == zero-force + same-filament exclusion
# --------------------------------------------------------------------------------------------------
def test_no_overlap_is_zero_force() -> None:
    """Nodes all farther apart than r_c feel exactly zero steric force (OFF == bit-identical)."""
    r_c = wca_cutoff(SIGMA)
    pos = np.array([[0.0, 0.0, 0.0], [3.0 * r_c, 0.0, 0.0], [0.0, 3.0 * r_c, 0.0]])
    fiber_id = np.array([0, 1, 2], dtype=np.int64)
    f = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    assert np.all(f == 0.0)
    assert total_steric_energy(pos, SIGMA, EPS, fiber_id) == 0.0


def test_same_fiber_pairs_excluded() -> None:
    """A single filament's own nodes never sterically self-repel (held by bending/inextensibility)."""
    r_c = wca_cutoff(SIGMA)
    # three closely-spaced nodes, ALL on fiber 0 (well within r_c) -> zero EV force
    pos = np.array([[0.0, 0.0, 0.0], [0.4 * r_c, 0.0, 0.0], [0.8 * r_c, 0.0, 0.0]])
    fiber_id = np.array([0, 0, 0], dtype=np.int64)
    f = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    assert np.allclose(f, 0.0, atol=1e-15)
    # relabel the middle node to a different fiber -> now it repels its neighbours
    fiber_id2 = np.array([0, 1, 0], dtype=np.int64)
    f2 = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id2)
    assert np.linalg.norm(f2) > 0.0


# --------------------------------------------------------------------------------------------------
# Newton's third law — momentum-clean internal force
# --------------------------------------------------------------------------------------------------
def test_newton_third_law_net_force_zero() -> None:
    """The internal steric force sums to ~0 (each pair pushes both nodes equal-and-opposite)."""
    pos, fiber_id = _crowded_cluster(seed=3)
    f = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    net = np.sum(f, axis=0)
    scale = np.sum(np.linalg.norm(f, axis=1)) + 1e-30
    assert np.linalg.norm(net) / scale < 1e-12


# --------------------------------------------------------------------------------------------------
# Gate 4 — a compressed patch resists volumetric collapse (finite steric volume)
# --------------------------------------------------------------------------------------------------
def test_compressed_patch_pushes_outward() -> None:
    """A crowded multi-fiber patch has positive steric energy and a net OUTWARD (expansive) virial."""
    pos, fiber_id = _crowded_cluster(seed=7)
    energy = total_steric_energy(pos, SIGMA, EPS, fiber_id)
    assert energy > 0.0                                   # overlaps store repulsive energy
    f = steric_force_bruteforce(pos, SIGMA, EPS, fiber_id)
    centroid = pos.mean(axis=0)
    virial = float(np.sum((pos - centroid) * f))          # sum r.F about the centroid
    assert virial > 0.0                                   # forces point outward -> resists collapse


def test_denser_patch_stores_more_energy() -> None:
    """Compressing the same patch (shrinking the box) raises the steric energy -> real steric volume."""
    pos, fiber_id = _crowded_cluster(seed=7)
    centroid = pos.mean(axis=0)
    loose = total_steric_energy(pos, SIGMA, EPS, fiber_id)
    squeezed_pos = centroid + 0.6 * (pos - centroid)      # pull every node 40% toward the centroid
    squeezed = total_steric_energy(squeezed_pos, SIGMA, EPS, fiber_id)
    assert squeezed > loose


# --------------------------------------------------------------------------------------------------
# Gate 3b — CFL k_EV bookkeeping : derived (not tuned to a crowding outcome), grid-invariant
# --------------------------------------------------------------------------------------------------
def test_cfl_stiffness_is_contact_stiffness_and_grid_free() -> None:
    """k_EV folded into kmax is the well curvature; it depends on (sigma, eps) only, not on any grid."""
    k = cfl_stiffness_term(SIGMA, EPS)
    assert k == pytest.approx(contact_stiffness(EPS, SIGMA), rel=1e-12)
    assert k == pytest.approx(K_EV, rel=1e-12)
    # capped-load variant gives the deeper-overlap stiffness, still finite and >= the contact value
    k_cap = cfl_stiffness_term(SIGMA, EPS, f_cap=100.0)
    assert np.isfinite(k_cap) and k_cap >= k


def _r_over_rc(k_ev: float, f_load: float) -> float:
    """Residual separation fraction r_eq/r_c of a crossing pair squeezed by ``f_load`` at stiffness ``k_ev``."""
    r_c = wca_cutoff(SIGMA)
    eps = epsilon_from_contact_stiffness(k_ev, SIGMA)
    return compressed_pair_separation(f_load, SIGMA, eps) / r_c


def test_no_interpenetration_insensitive_to_k_ev_above_threshold() -> None:
    """The 'held apart' outcome is INSENSITIVE to k_EV above the DERIVED threshold -> k_EV is not a fit.

    The Magic-Number-Block way to pick k_EV: choose it so the residual overlap under the operating load is
    within a numerical penetration tolerance (here <=10% of the contact distance, r_eq/r_c >= 0.90). That
    threshold is DERIVED from (load, tolerance), not tuned to a crowding density. Above it, a stiffer k_EV
    only pushes r_eq closer to r_c; the verdict (fibers do not interpenetrate) is a robust plateau. A value
    tuned to hit a crowding target would instead make the outcome swing with k_EV.
    """
    f_load = 20.0                 # pN compressive load on a crossing pair
    tol = 0.90                    # penetration tolerance: r_eq/r_c >= 0.90 (<=10% overlap)

    # DERIVE the threshold k_EV that just meets the tolerance (bisection on the monotone r_eq/r_c(k_ev)).
    lo, hi = 1.0, 1.0e6
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if _r_over_rc(mid, f_load) < tol:
            lo = mid
        else:
            hi = mid
    k_thr = hi
    assert _r_over_rc(k_thr, f_load) == pytest.approx(tol, abs=2e-3)

    # Sweep from the threshold up two decades; the outcome holds and saturates toward the contact distance.
    k_sweep = k_thr * np.array([1.0, 3.0, 10.0, 30.0, 100.0])
    frac = np.array([_r_over_rc(k, f_load) for k in k_sweep])
    assert np.all(frac >= tol - 1e-9)          # no interpenetration for EVERY k_EV above threshold
    assert np.all(np.diff(frac) > 0.0)         # stiffer -> nearer contact, monotonically
    assert frac[-1] - frac[-2] < 0.01          # high end saturates -> insensitive, not a tuned knife-edge
