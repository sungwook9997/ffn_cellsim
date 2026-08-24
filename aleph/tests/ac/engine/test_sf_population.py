"""Structural gates for the disjoint ``sf_arc`` F-actin population builder (NumPy only; native forbidden).

Covers the Sanity Gate declared in ``ac.engine.sf_population``:
  * ownership / no-double-count (unique sf_arc-block ids; disjoint from a cortex ledger; N_unique == N_summed),
  * connector conformance (feeds exactly the canonical fa_actin_anchor / sf_cortex_transient / actin_cap_linc /
    dorsal_arc_crosslink edges of ``reference_cell_architecture``; each bidirectional + adjoint),
  * emergent prestress (ventral both-ends FA = contractile CLOSED dipole; dorsal one-end = OPEN; no lumped k_SF),
  * anchor taxonomy (anchor_kind emerges from the anchor node set, not the class tag),
  * firewall (load paths invariant to a class-tag permutation),
  * GAP guard (no f_head -> tension None + GAP status; never a tuned magnitude),
  * no-resample (the fixed-capacity ledger raises rather than growing to a target).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from aleph.engine.contracts import ConnectorScope, reference_cell_architecture
from aleph.engine.population import PopulationLedger
from aleph.engine.sf_population import (
    ACTIN_CAP_LINC,
    DORSAL_ARC_CROSSLINK,
    FA_ACTIN_ANCHOR,
    SF_COMPONENT,
    SF_CORTEX_TRANSIENT,
    SFArcPopulation,
    SFBundle,
    SFClass,
    build_sf_arc_population,
    validate_sf_connectors,
)
from aleph.components.weave.stress_fiber import GAP_STATUS


def _pop(**kw) -> SFArcPopulation:
    return build_sf_arc_population(**kw)


def _by_class(pop: SFArcPopulation, sf_class: str) -> list:
    return [b for b in pop.bundles if b.sf_class == sf_class]


# ── ownership / no-double-count ───────────────────────────────────────────────────────────────────
def test_population_owns_unique_disjoint_ids() -> None:
    pop = _pop()
    pop.assert_partitioned()
    c = pop.census()
    assert c["N_unique_sf_fibers"] == c["N_sf_fibers_summed"]
    assert c["active_count"] == c["N_sf_fibers_summed"]
    lo, hi = pop.ledger.block
    for b in pop.bundles:
        assert all(lo <= int(fid) < hi for fid in b.global_fiber_ids)   # every id in the sf_arc block


def test_disjoint_from_cortex_population() -> None:
    """The PI 2026-07-22 no-double-count invariant: sf_arc shares no id block/active id with the cortex."""
    pop = _pop(id_base=1_000_000)
    cortex = PopulationLedger("cortex", id_base=0, capacity=70_686)
    cortex.seed_active(range(70_686))
    pop.assert_disjoint_from(cortex)                                    # non-overlapping blocks -> OK


def test_overlapping_block_is_rejected_as_double_count() -> None:
    pop = _pop(id_base=100)
    clash = PopulationLedger("cortex", id_base=0, capacity=200)         # block [0,200) overlaps sf_arc [100,..)
    clash.seed_active(range(1))
    with pytest.raises(ValueError, match="overlapping global"):
        pop.assert_disjoint_from(clash)


# ── connector conformance (the ONLY couplings; co-location is never a connection) ──────────────────
def test_feeds_exactly_the_four_canonical_sf_edges() -> None:
    pop = _pop()
    validate_sf_connectors(reference_cell_architecture(), pop)          # must not raise
    assert set(pop.connector_endpoints()) == {
        FA_ACTIN_ANCHOR, SF_CORTEX_TRANSIENT, ACTIN_CAP_LINC, DORSAL_ARC_CROSSLINK,
    }


def test_connector_endpoints_reference_real_sf_nodes() -> None:
    pop = _pop()
    ep = pop.connector_endpoints()
    n = pop.n_nodes
    for name in (FA_ACTIN_ANCHOR, SF_CORTEX_TRANSIENT, ACTIN_CAP_LINC):
        idx = ep[name]
        assert idx.size > 0 and idx.min() >= 0 and idx.max() < n        # valid global node indices
    joints = ep[DORSAL_ARC_CROSSLINK]
    assert joints.ndim == 2 and joints.shape[1] == 2                    # (dorsal-free, arc) pairs
    assert joints.max() < n


def test_validate_rejects_non_sf_incident_connector() -> None:
    arch = reference_cell_architecture()
    # break the FA edge so it no longer touches sf_arc
    broken = tuple(
        replace(c, component_a="lamellipodium") if c.name == FA_ACTIN_ANCHOR else c
        for c in arch.connectors
    )
    from aleph.engine.contracts import CellArchitecture

    with pytest.raises(ValueError, match="must join sf_arc"):
        validate_sf_connectors(CellArchitecture(arch.components, broken), _pop())


def test_fa_and_dorsal_arc_scopes() -> None:
    arch = reference_cell_architecture()
    by = {c.name: c for c in arch.connectors}
    assert by[FA_ACTIN_ANCHOR].scope is ConnectorScope.INTER_COMPONENT
    assert by[DORSAL_ARC_CROSSLINK].scope is ConnectorScope.INTERNAL
    assert by[DORSAL_ARC_CROSSLINK].component_a == by[DORSAL_ARC_CROSSLINK].component_b == SF_COMPONENT


# ── emergent prestress (NMII cumsum, no lumped k_SF) ───────────────────────────────────────────────
def test_ventral_is_contractile_closed_dipole() -> None:
    pop = _pop()
    ventral = _by_class(pop, SFClass.VENTRAL)[0]
    lp = ventral.load_path(f_head_pN=0.5)
    assert lp.contractile
    assert lp.anchor_kind == "fa"
    assert lp.closed_dipole
    assert lp.balance_residual_per_fhead < 1e-9                         # self-equilibrated both-ends SF
    assert min(lp.load_path_thirds().values()) > 0.0                    # positive prestress across the bundle


def test_dorsal_one_end_is_open_residual() -> None:
    pop = _pop()
    dorsal = _by_class(pop, SFClass.DORSAL)[0]
    lp = dorsal.load_path(f_head_pN=0.5)
    assert lp.dipole_per_fhead > 0.0
    assert not lp.closed_dipole
    assert lp.balance_residual_per_fhead == pytest.approx(lp.dipole_per_fhead)   # network carries the free end
    assert dorsal.free_local.size == 1                                  # exactly one dorsal free end


def test_cap_is_linc_anchored_and_arc_unanchored() -> None:
    pop = _pop()
    cap = _by_class(pop, SFClass.PERINUCLEAR_CAP)[0].load_path(f_head_pN=0.5)
    arc = _by_class(pop, SFClass.TRANSVERSE_ARC)[0].load_path(f_head_pN=0.5)
    assert cap.anchor_kind == "linc" and cap.contractile
    assert arc.anchor_kind == "none" and arc.contractile               # arc has motors but no FA/LINC


def test_dorsal_arc_joints_link_free_ends_to_arcs() -> None:
    pop = _pop(n_dorsal=2, n_arc=2)
    joints = pop.connector_endpoints()[DORSAL_ARC_CROSSLINK]
    assert joints.shape[0] == 2                                        # one joint per dorsal free end
    # each joint's first index is a dorsal free-end global node; the second is a transverse-arc node.
    dorsal_free_globals = {
        int(b.node_base + b.free_local[0])
        for b in _by_class(pop, SFClass.DORSAL) if b.free_local.size
    }
    arc_globals = {
        int(b.node_base + k)
        for b in _by_class(pop, SFClass.TRANSVERSE_ARC)
        for k in range(b.pos.shape[0])
    }
    for free_node, arc_node in joints:
        assert int(free_node) in dorsal_free_globals
        assert int(arc_node) in arc_globals


# ── firewall: physics never reads the class tag ────────────────────────────────────────────────────
def test_load_paths_invariant_to_class_tag_permutation() -> None:
    """FIREWALL: prestress emerges from motors + anchors; permuting the class tag must not change it."""
    pop = _pop()
    base = sorted(lp.dipole_per_fhead for lp in pop.load_paths(f_head_pN=0.5))
    scrambled_bundles = [replace(b, sf_class=SFClass.ALL[i % len(SFClass.ALL)])
                         for i, b in enumerate(pop.bundles)]
    scrambled = replace(pop, bundles=scrambled_bundles)
    perm = sorted(lp.dipole_per_fhead for lp in scrambled.load_paths(f_head_pN=0.5))
    assert base == pytest.approx(perm)


# ── GAP guard + no-resample ────────────────────────────────────────────────────────────────────────
def test_gap_guard_default_no_magnitude() -> None:
    pop = _pop()
    assert pop.magnitude_status == GAP_STATUS
    lps = pop.load_paths()                                             # no f_head -> GAP
    assert all(lp.tension is None and lp.dipole_pN is None for lp in lps)
    assert any(np.any(lp.tension_per_fhead > 0.0) for lp in lps)       # magnitude-free SHAPE still populated


def test_finding_opt_in_scales_all_bundles() -> None:
    pop = _pop(f_head_pN=0.5)
    assert pop.magnitude_status != GAP_STATUS
    assert all(lp.tension is not None for lp in pop.load_paths(f_head_pN=0.5))


def test_capacity_is_fixed_no_resample_to_target() -> None:
    pop = _pop()
    # the whole fixed block is seeded active -> a further allocate must raise, never grow to a target density
    with pytest.raises(RuntimeError, match="free-list exhausted"):
        pop.ledger.allocate()


def test_counts_scale_with_requested_instances() -> None:
    pop = _pop(n_ventral=3, n_dorsal=1, n_arc=1, n_cap=1)
    c = pop.census()
    assert c["N_bundles"] == 6
    assert c["N_sf_fibers_summed"] == 12                               # two filaments per SF instance
    assert c["N_by_class"][SFClass.VENTRAL] == 6
