"""Self-test of the stress-fiber prestress + traction-loop oracle (NumPy only; label-blind; GAP-guarded).

Covers the Sanity Gate declared in ``ac.weave.stress_fiber``:
  * dimensional/convention (unit axis; shape dimensionless; tension in pN when f_head given),
  * sign-sense (antiparallel FA<->FA is contractile; a balanced both-ends SF is a closed dipole),
  * measurement-protocol (ventral both-ends = uniform positive shaft prestress; dorsal one-end = open residual),
  * firewall (the wiring is invariant to a region_id permutation — it never reads the label),
  * ledger (bundles partition the anchored-motorized set; cortex fibers are excluded, not double-counted),
  * GAP guard (no f_head -> tension None + GAP status; never a tuned number).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.stress_fiber import (
    GAP_STATUS,
    stress_fiber_load_path,
    extract_sf_bundles,
    wire_stress_fibers,
)


# ── synthetic bundles ────────────────────────────────────────────────────────────────────────────
def _antiparallel_sarcomere():
    """Two anti-parallel filaments (barbed ends OUTWARD) along x, an FA at each outer end, one bipolar motor."""
    A = np.column_stack([np.linspace(0.0, 2.0, 5), np.zeros(5), np.zeros(5)])   # nodes 0..4
    B = np.column_stack([np.linspace(2.0, 4.0, 5), np.zeros(5), np.zeros(5)])   # nodes 5..9
    pos = np.vstack([A, B])
    foff = np.array([0, 5, 10], np.int64)
    polarity = np.array([-1, +1], np.int64)         # A barbed at x=0 (first node), B barbed at x=4 (last node)
    myo_i, myo_j = np.array([4]), np.array([5])      # interior overlap
    return pos, foff, polarity, myo_i, myo_j


def test_axis_and_shape_are_dimensionally_clean() -> None:
    pos, foff, pol, mi, mj = _antiparallel_sarcomere()
    lp = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0, 9]), np.zeros(0, np.int64), f_head_pN=0.5)
    assert abs(np.linalg.norm(lp.axis) - 1.0) < 1e-12          # unit axis
    assert lp.tension_per_fhead.shape == lp.s_bins.shape        # shape aligned to cross-sections
    assert lp.tension is not None and lp.tension.shape == lp.s_bins.shape
    np.testing.assert_allclose(lp.tension, lp.tension_per_fhead * 0.5)   # tension = f_head * shape (pN)
    assert lp.dipole_pN == pytest.approx(lp.dipole_per_fhead * 0.5)


def test_antiparallel_bundle_is_contractile_closed_dipole() -> None:
    """SIGN-SENSE: barbed ends outward + FA at both ends -> contractile, self-equilibrated dipole (residual 0)."""
    pos, foff, pol, mi, mj = _antiparallel_sarcomere()
    lp = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0, 9]), np.zeros(0, np.int64), f_head_pN=0.5)
    assert lp.contractile
    assert lp.dipole_per_fhead > 0.0
    assert lp.balance_residual_per_fhead < 1e-9                 # closed both-ends SF
    assert lp.closed_dipole
    thirds = lp.load_path_thirds()
    assert min(thirds.values()) > 0.0                          # positive prestress across the whole bundle (tension)


def test_dorsal_one_end_is_open_residual() -> None:
    """MEASUREMENT-PROTOCOL: a single anchored end hands a net residual to the network (not a closed dipole)."""
    pos, foff, pol, mi, mj = _antiparallel_sarcomere()
    lp = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0]), np.zeros(0, np.int64), f_head_pN=0.5)
    assert lp.dipole_per_fhead > 0.0
    assert lp.balance_residual_per_fhead == pytest.approx(lp.dipole_per_fhead)   # network carries the free end
    assert not lp.closed_dipole


def test_non_contractile_bundle_has_no_dipole() -> None:
    """A fascin-only (no NMII) bundle is anchored but generates no contractile prestress (filopodium taxonomy)."""
    pos, foff, pol, _, _ = _antiparallel_sarcomere()
    lp = stress_fiber_load_path(pos, foff, pol, np.zeros(0, np.int64), np.zeros(0, np.int64),
                                np.array([9]), np.zeros(0, np.int64))
    assert not lp.contractile
    assert lp.dipole_per_fhead == 0.0


def test_gap_guard_no_magnitude_without_f_head() -> None:
    """GAP guard: no provisional per-head force -> tension None + GAP status (never a tuned number)."""
    pos, foff, pol, mi, mj = _antiparallel_sarcomere()
    lp = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0, 9]), np.zeros(0, np.int64))
    assert lp.tension is None and lp.dipole_pN is None
    assert lp.magnitude_status == GAP_STATUS
    # the magnitude-independent SHAPE is still fully populated (the analytic gate is magnitude-free)
    assert np.any(lp.tension_per_fhead > 0.0)


def test_anchor_kind_taxonomy() -> None:
    """anchor_kind emerges from the anchor node set (fa / linc / fa+linc), never from a region label."""
    pos, foff, pol, mi, mj = _antiparallel_sarcomere()
    fa = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0, 9]), np.zeros(0, np.int64))
    linc = stress_fiber_load_path(pos, foff, pol, mi, mj, np.zeros(0, np.int64), np.array([0, 9]))
    both = stress_fiber_load_path(pos, foff, pol, mi, mj, np.array([0]), np.array([9]))
    assert (fa.anchor_kind, linc.anchor_kind, both.anchor_kind) == ("fa", "linc", "fa+linc")


# ── whole-cell wiring on a real woven network ─────────────────────────────────────────────────────
def _small_cell(with_cortex: bool = False):
    from aleph.components.weave.woven_cell import weave_cell
    from aleph.components.weave.regions import (
        VENTRAL_SF_REGION, DORSAL_SF_REGION, PERINUCLEAR_CAP_REGION, FILOPODIUM_REGION, CORTEX_REGION)
    regions = [VENTRAL_SF_REGION, DORSAL_SF_REGION, PERINUCLEAR_CAP_REGION, FILOPODIUM_REGION]
    if with_cortex:
        regions = [CORTEX_REGION] + regions
    return weave_cell(regions, rng=np.random.default_rng(0))


def test_bundles_partition_anchored_set_no_double_count() -> None:
    cell = _small_cell()
    w = wire_stress_fibers(cell)
    led = w.ledger()
    assert led["N_unique_sf_fibers"] == led["N_sf_fibers_summed"]     # a fiber lands in at most one bundle
    # every discovered bundle fiber carries an FA or LINC anchor (SF = anchored bundle)
    anchor_fibers = set()
    for sites in (cell.fa_sites, cell.linc_sites):
        if sites.size:
            anchor_fibers.update(int(f) for f in np.unique(cell.node_fiber[sites]))
    all_bundle_fibers = set(int(f) for arr in w.global_fibers for f in arr)
    assert all_bundle_fibers <= anchor_fibers


def test_cortex_fibers_are_excluded_not_double_counted() -> None:
    """The (un-anchored) cortex is NOT a stress fiber: its fibers must not enter any SF bundle."""
    cell = _small_cell(with_cortex=True)
    w = wire_stress_fibers(cell)
    cortex = cell.region_slice("cortex")                     # diagnostic-only, used by the TEST (not the module)
    cortex_fibers = set(range(cortex.start, cortex.stop))
    bundle_fibers = set(int(f) for arr in w.global_fibers for f in arr)
    assert cortex_fibers.isdisjoint(bundle_fibers)


def test_firewall_invariant_to_region_id_permutation() -> None:
    """FIREWALL: the wiring never reads region_id — permuting the label must not change the result."""
    import dataclasses

    cell = _small_cell()
    base = wire_stress_fibers(cell, f_head_pN=0.5)
    scrambled = dataclasses.replace(cell, region_id=cell.region_id[::-1].copy())
    perm = wire_stress_fibers(scrambled, f_head_pN=0.5)
    assert base.ledger()["N_unique_sf_fibers"] == perm.ledger()["N_unique_sf_fibers"]
    assert sorted(lp.dipole_pN for lp in base.load_paths) == \
        pytest.approx(sorted(lp.dipole_pN for lp in perm.load_paths))


def test_wiring_gap_default_and_finding_opt_in() -> None:
    cell = _small_cell()
    assert wire_stress_fibers(cell).magnitude_status == GAP_STATUS            # GAP by default
    w = wire_stress_fibers(cell, f_head_pN=0.5)
    assert all(lp.dipole_pN is not None for lp in w.load_paths)               # provisional FINDING opt-in
    # the seed is isotropic-conditional -> it has NOT condensed into a clean closed dipole yet (I5's proof)
    assert w.ledger()["N_sf_bundles"] >= 3
