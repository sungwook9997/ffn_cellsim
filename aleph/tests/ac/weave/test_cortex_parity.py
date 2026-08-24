"""Gate-1: weave_cell([CORTEX]) all-OFF is BIT-IDENTICAL to ff.weave.weave(CORTEX) — the master's I4 gate.

The cortex region DELEGATES to ff.weave.weave(CORTEX) with untouched RNG order, so the assembled network is
array-for-array identical => identical gamma (a deterministic function of an identical network). Parity is
POPULATION-INDEPENDENT (it is a property of the delegation, not of a particular count), so this CPU gate runs a
SMALL cortex for speed; the FULL 70,686-filament parity is the lead's native gate (INTEGRATION.md §3).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from aleph.laws.architecture_spec import CORTEX
from aleph.laws.weave import weave as ff_weave
from aleph.components.weave.regions import RegionSpec
from aleph.components.weave.woven_cell import weave_cell


def _small_cortex(n: int = 200):
    """CORTEX with a reduced filament count (delegation parity is population-independent)."""
    fil = dataclasses.replace(CORTEX.filament, n_filaments=n)
    return dataclasses.replace(CORTEX, filament=fil)


def test_gate1_cortex_bit_identical() -> None:
    """Every assembled array (pos, crosslinks, motors) is bit-identical to ff.weave.weave(CORTEX)."""
    spec = _small_cortex(200)
    ref = ff_weave(spec, rng=np.random.default_rng(0))
    wc = weave_cell([RegionSpec(arch=spec)], rng=np.random.default_rng(0))

    assert np.array_equal(wc.pos, ref.net.pos)
    assert np.array_equal(wc.fiber_offsets, ref.net.fiber_offsets)
    assert np.array_equal(wc.xl_i, ref.xl_i)
    assert np.array_equal(wc.xl_j, ref.xl_j)
    assert np.array_equal(wc.xl_k, ref.xl_k)
    assert np.array_equal(wc.xl_rest, ref.xl_rest)
    assert np.array_equal(wc.myo_i, ref.myo_i)
    assert np.array_equal(wc.myo_j, ref.myo_j)


def test_gate1_crosslinked_cortex_roundtrip() -> None:
    """to_crosslinked_cortex() reproduces the ff CrosslinkedCortex network (the gamma harness input)."""
    spec = _small_cortex(150)
    ref = ff_weave(spec, rng=np.random.default_rng(0))
    wc = weave_cell([RegionSpec(arch=spec)], rng=np.random.default_rng(0))
    cc = wc.to_crosslinked_cortex()
    assert np.array_equal(cc.net.pos, ref.net.pos)
    assert np.array_equal(cc.xl_i, ref.xl_i)
    assert np.array_equal(cc.myo_i, ref.myo_i)
    # R0_mean (the turgor volume reference) matches to round-off (same node cloud)
    assert cc.R0_mean == pytest.approx(ref.R0_mean, rel=1e-12)


def test_gate1_all_other_regions_off() -> None:
    """The single-cortex weave has no branch/cross-region/anchor content (all other regions genuinely OFF)."""
    wc = weave_cell([RegionSpec(arch=_small_cortex(120))], rng=np.random.default_rng(0))
    assert wc.branch_triples.shape[0] == 0
    assert wc.xbond_i.shape[0] == 0
    assert wc.fa_sites.shape[0] == 0
    assert wc.linc_sites.shape[0] == 0
    assert wc.region_names == ["cortex"]
    assert wc.provenance["cortex"] == "delegated_ff_weave"
    assert wc.relaxation_channel == "kmc_reattach"     # the single selected crosslinker channel
