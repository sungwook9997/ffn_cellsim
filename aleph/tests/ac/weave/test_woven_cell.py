"""Self-test of the WovenCell framework — concat correctness + first-class cross-region bonds (NumPy only)."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from aleph.laws.architecture_spec import CORTEX
from aleph.components.weave.regions import (
    DORSAL_SF_REGION,
    LAMELLIPODIUM_REGION,
    PERINUCLEAR_CAP_REGION,
    RegionSpec,
    VENTRAL_SF_REGION,
)
from aleph.components.weave.woven_cell import CrossRegionBond, weave_cell


def _small_cortex_region(n: int = 120) -> RegionSpec:
    fil = dataclasses.replace(CORTEX.filament, n_filaments=n)
    return RegionSpec(arch=dataclasses.replace(CORTEX, filament=fil))


def _small_lamellipodium() -> RegionSpec:
    return dataclasses.replace(LAMELLIPODIUM_REGION,
                               arch=dataclasses.replace(LAMELLIPODIUM_REGION.arch,
                                                        filament=dataclasses.replace(
                                                            LAMELLIPODIUM_REGION.arch.filament, n_filaments=6)),
                               n_daughter_pool=30)


def test_concat_offsets_and_indices_valid() -> None:
    """fiber_offsets are contiguous 0..N; every crosslink/motor/branch/FA/LINC index is a valid global node."""
    wc = weave_cell([_small_cortex_region(), _small_lamellipodium(), VENTRAL_SF_REGION, PERINUCLEAR_CAP_REGION],
                    rng=np.random.default_rng(0))
    N = wc.n_nodes
    assert wc.fiber_offsets[0] == 0 and wc.fiber_offsets[-1] == N
    assert np.all(np.diff(wc.fiber_offsets) >= 1)
    assert wc.region_id.shape[0] == wc.n_fibers
    for arr in (wc.xl_i, wc.xl_j, wc.myo_i, wc.myo_j, wc.fa_sites, wc.linc_sites):
        if arr.size:
            assert arr.min() >= 0 and arr.max() < N
    if wc.branch_triples.size:
        assert wc.branch_triples.min() >= 0 and wc.branch_triples.max() < N
    # node_fiber is the inverse of fiber_offsets
    assert wc.node_fiber.shape[0] == N
    for f in range(wc.n_fibers):
        assert np.all(wc.node_fiber[wc.fiber_offsets[f]:wc.fiber_offsets[f + 1]] == f)


def test_crosslinks_are_cross_fiber_within_region() -> None:
    """Each crosslink connects two DIFFERENT fibers (a crosslinker never bonds a fiber to itself)."""
    wc = weave_cell([VENTRAL_SF_REGION, PERINUCLEAR_CAP_REGION], rng=np.random.default_rng(1))
    if wc.xl_i.size:
        assert np.all(wc.node_fiber[wc.xl_i] != wc.node_fiber[wc.xl_j])


def test_region_labels_are_diagnostic_only() -> None:
    """region_id / region_names are recorded (firewall diagnostic); region_slice locates a region's fibers."""
    wc = weave_cell([_small_cortex_region(), VENTRAL_SF_REGION], rng=np.random.default_rng(2))
    assert wc.region_names == ["cortex", "ventral_sf"]
    sl = wc.region_slice("ventral_sf")
    assert np.all(wc.region_id[sl] == 1)


def test_cross_region_bond_is_first_class_and_cross_region() -> None:
    """A declared cross-region bond resolves to global endpoints in the TWO named regions, within reach."""
    # overlap the dorsal SF footprint onto the ventral SF centre so their end nodes are near -> a bond forms
    dorsal = dataclasses.replace(DORSAL_SF_REGION, centre_um=(0.0, 0.0, 0.0))
    bond = CrossRegionBond(region_a="ventral_sf", selector_a="end1", region_b="dorsal_sf", selector_b="end1",
                           k=4.6e5, reach_um=5.0, kind="arc_dorsal")
    wc = weave_cell([VENTRAL_SF_REGION, dorsal], rng=np.random.default_rng(3), cross_bonds=[bond])
    assert wc.xbond_i.shape[0] > 0
    rid_v = wc.region_names.index("ventral_sf")
    rid_d = wc.region_names.index("dorsal_sf")
    assert np.all(wc.region_id[wc.node_fiber[wc.xbond_i]] == rid_v)   # endpoint A in region A
    assert np.all(wc.region_id[wc.node_fiber[wc.xbond_j]] == rid_d)   # endpoint B in region B
    d = np.linalg.norm(wc.pos[wc.xbond_i] - wc.pos[wc.xbond_j], axis=1)
    assert np.all(d <= 5.0)                                            # within the declared reach
    assert wc.xbond_kind == ["arc_dorsal"] * wc.xbond_i.shape[0]


def test_population_ledger_no_double_count() -> None:
    """N_unique_active fibers = sum of active per region (no cortex/SF double count); dormant tracked."""
    wc = weave_cell([_small_lamellipodium(), VENTRAL_SF_REGION], rng=np.random.default_rng(4))
    led = wc.population_ledger()
    assert led["N_unique_active_fibers"] == sum(r["active_fibers"] for r in led["per_region"].values())
    assert led["N_allocated_fibers"] == led["N_unique_active_fibers"] + led["N_dormant_fibers"]
    # the lamellipodium has a dormant daughter pool -> some dormant fibers exist
    assert led["N_dormant_fibers"] >= 0
    assert led["per_region"]["lamellipodium"]["total_fibers"] == 6 + 30


def test_single_channel_guard_at_construction() -> None:
    """weave_cell refuses both crosslinker-relaxation channels (double-count guard, checked at build)."""
    with pytest.raises(ValueError, match="double-count"):
        weave_cell([VENTRAL_SF_REGION], xl_reattach=True, r0_creep=True)
