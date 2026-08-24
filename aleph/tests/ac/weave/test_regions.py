"""Self-test of the region builders — emergent seeds valid; demoted control; cross-region concat (NumPy only)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.regions import (
    DORSAL_SF_REGION,
    FILOPODIUM_REGION,
    PERINUCLEAR_CAP_REGION,
    TRANSVERSE_ARC_REGION,
    VENTRAL_SF_REGION,
    build_region,
)


@pytest.mark.parametrize("region", [
    VENTRAL_SF_REGION, DORSAL_SF_REGION, TRANSVERSE_ARC_REGION, PERINUCLEAR_CAP_REGION, FILOPODIUM_REGION])
def test_emergent_region_builds_valid_leaf(region) -> None:
    """Each bundle region seeds a valid RegionLeaf: contiguous offsets, per-fiber polarity, valid anchors."""
    leaf = build_region(region, np.random.default_rng(0))
    n = leaf.pos.shape[0]
    n_fib = leaf.fiber_offsets.shape[0] - 1
    assert leaf.fiber_offsets[0] == 0 and leaf.fiber_offsets[-1] == n
    assert leaf.polarity.shape[0] == n_fib
    assert leaf.provenance == "isotropic_conditional_seed"
    if leaf.xl_i.size:
        assert leaf.xl_i.max() < n and leaf.xl_j.max() < n
    for arr in (leaf.fa_sites, leaf.linc_sites):
        if arr.size:
            assert arr.max() < n


def test_anchor_roles() -> None:
    """Ventral SF -> FA at both ends; perinuclear cap -> LINC stubs; filopodium -> tip FA."""
    v = build_region(VENTRAL_SF_REGION, np.random.default_rng(0))
    assert v.fa_sites.size == 2 * (v.fiber_offsets.shape[0] - 1)      # both ends
    cap = build_region(PERINUCLEAR_CAP_REGION, np.random.default_rng(0))
    assert cap.linc_sites.size == (cap.fiber_offsets.shape[0] - 1)    # one LINC stub per cap fiber
    assert cap.fa_sites.size == 0
    fil = build_region(FILOPODIUM_REGION, np.random.default_rng(0))
    assert fil.fa_sites.size == (fil.fiber_offsets.shape[0] - 1)      # tip end only


def test_emergent_seed_is_not_a_scripted_bundle() -> None:
    """The emergent seed is isotropic-conditional (low nematic order), NOT the demoted parallel bundle."""
    from aleph.components.emergence.nematic import fiber_axes, nematic_order

    leaf = build_region(VENTRAL_SF_REGION, np.random.default_rng(0))
    s_seed = nematic_order(fiber_axes(leaf.pos, leaf.fiber_offsets))   # UNIT end-to-end axes
    assert s_seed < 0.9                                              # not a perfect parallel bundle (must condense)


def test_demoted_scripted_bundle_control() -> None:
    """emergent=False reproduces the DEMOTED ff scripted bundle (a labeled non-authoritative control)."""
    import dataclasses

    demoted = dataclasses.replace(VENTRAL_SF_REGION, emergent=False)
    leaf = build_region(demoted, np.random.default_rng(0))
    assert leaf.provenance == "demoted_scripted_bundle"
    n = leaf.pos.shape[0]
    assert leaf.fiber_offsets[0] == 0 and leaf.fiber_offsets[-1] == n
