"""CPU-green MEASURED gates for the SF-straddle NMII minifilament population.

Every physical claim below is a measurement against the landed pure-NumPy point-to-segment oracle
(:func:`aleph.components.motor.segment_query.point_to_segment_query_reference`), not an assertion about how the
builder was written:

* the ``±head_offset`` heads land ON the two anti-parallel filament LINES (perpendicular residual ≈ 0) —
  and the legacy co-located sarcomere provably does not (residual = ``head_offset``, i.e. binding sits at the
  very edge of the capture reach);
* the ``+`` and ``−`` head groups land on two DIFFERENT filaments whose barbed directions are anti-parallel,
  so the straddling minifilament is a contractile dipole (bipolar dot ≈ −1);
* the population owns its own particle array and its own disjoint ID block (no cortex/SF double count);
* a population whose sarcomere lacks the straddle overlap is REFUSED (no arbitrary-perpendicular fallback).

No Warp kernel is launched (CPU-Warp launches are forbidden in the test lane); the CUDA behaviour is the
Lead's native gate on the A5000.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.population import PopulationLedger
from aleph.engine.sf_mechanics import build_sf_mechanics_topology
from aleph.engine.sf_nmii_population import (
    barbed_ward_tangents,
    build_sf_straddle_nmii_population,
)
from aleph.engine.sf_population import build_sf_arc_population
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.segment_query import point_to_segment_query_reference

#: The NMII reference topology the sarcomere geometry must be derived FROM (I0-B3 layout; magnitudes are GAPs).
_TOPOLOGY = MinifilamentTopology(
    n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)
#: NMII head→actin capture reach [µm] (``ff.hand_kmc`` NMIIA, KU-3.5 reconciled) — used only to report margin.
_CAPTURE_UM = 0.210
_ID_BASE = 3_000_000


def _population(*, straddle: bool = True, n_per_fiber: int = 9):
    kwargs = {}
    if straddle:
        kwargs = {
            "sarcomere_lateral_um": 2.0 * _TOPOLOGY.head_offset_um,
            "sarcomere_overlap_um": _TOPOLOGY.backbone_length_um,
        }
    return build_sf_arc_population(
        n_ventral=4, n_dorsal=2, n_arc=2, n_cap=2, n_per_fiber=n_per_fiber, **kwargs)


def _segment_frame(pop):
    """Return ``(seg_a, seg_b, barbed, filament_id)`` for the SF axial segments, in link order."""
    topology = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0)
    seg_a = topology.links[:, 0].astype(np.int64)
    seg_b = topology.links[:, 1].astype(np.int64)
    delta = pop.pos[seg_b] - pop.pos[seg_a]
    barbed = delta * topology.seg_polarity[:, None].astype(np.float64)
    barbed /= np.linalg.norm(barbed, axis=1, keepdims=True)
    return seg_a, seg_b, barbed, topology.seg_filament_id


def _head_query(pop, nmii):
    seg_a, seg_b, barbed, filament_id = _segment_frame(pop)
    heads = nmii.pos[nmii.head_node]
    result = point_to_segment_query_reference(heads, pop.pos[seg_a], pop.pos[seg_b], barbed, _CAPTURE_UM)
    return result, barbed, filament_id


def test_straddle_puts_every_head_on_an_actin_line_and_legacy_does_not() -> None:
    """MEASURED: the straddle residual collapses to ~0 while the legacy sarcomere leaves it at head_offset."""
    pop = _population(straddle=True)
    nmii = build_sf_straddle_nmii_population(pop, _TOPOLOGY, id_base=_ID_BASE)
    result, _, _ = _head_query(pop, nmii)
    dist = result["dist"]

    assert dist.size == nmii.n_heads > 0
    # the physical claim: the heads sit ON the filament lines, not merely inside the capture ball.
    assert float(np.median(dist)) < 0.01 * _TOPOLOGY.head_offset_um
    assert float(dist.max()) < 0.5 * _TOPOLOGY.head_offset_um
    assert int((result["seg_id"] >= 0).sum()) == nmii.n_heads          # all within capture, with margin

    # the legacy co-located sarcomere cannot be straddled at all — the builder refuses it (next test); the
    # geometric reason is measurable here: its two filaments meet at one point, so any placement leaves the
    # heads a full head_offset away from the actin, i.e. at the edge of the capture reach.
    legacy = _population(straddle=False)
    bundle = legacy.bundles[0]
    node_a, node_b = int(bundle.myo_i[0]), int(bundle.myo_j[0])
    separation = float(np.linalg.norm(bundle.pos[node_b] - bundle.pos[node_a]))
    assert separation == pytest.approx(0.0, abs=1e-12)


def test_plus_and_minus_heads_land_on_the_two_antiparallel_filaments() -> None:
    """MEASURED bipolar straddle: the two head groups bind two DIFFERENT, anti-parallel filaments."""
    pop = _population(straddle=True)
    nmii = build_sf_straddle_nmii_population(pop, _TOPOLOGY, id_base=_ID_BASE)
    result, barbed, filament_id = _head_query(pop, nmii)
    seg_id, half = result["seg_id"], _TOPOLOGY.n_heads_per_side

    dots: list[float] = []
    for m in range(nmii.n_minifilaments):
        block = slice(m * 2 * half, (m + 1) * 2 * half)
        sid = seg_id[block]
        assert int((sid >= 0).sum()) == 2 * half                      # every head of this motor found actin
        plus = set(filament_id[sid[:half]].tolist())
        minus = set(filament_id[sid[half:]].tolist())
        assert len(plus) == 1 and len(minus) == 1                     # each side hugs ONE filament
        assert not (plus & minus)                                     # and they are DIFFERENT filaments
        w_plus = barbed[sid[:half]].mean(axis=0)
        w_minus = barbed[sid[half:]].mean(axis=0)
        dots.append(float(w_plus @ w_minus / (np.linalg.norm(w_plus) * np.linalg.norm(w_minus))))

    # Only STRAIGHT sarcomeres are motor stations, and a straight sarcomere's two filaments are EXACTLY
    # anti-parallel — so the bound is the geometric identity dot = −1, not a threshold placed around whatever
    # the run happened to produce.  (A curved bundle would sit at −0.600; it is excluded, see the next test.)
    assert max(dots) == pytest.approx(-1.0, abs=1e-9)


def test_curved_perinuclear_cap_is_excluded_from_motor_stations_and_counted() -> None:
    """A rigid bipolar minifilament cannot straddle an ARCHED sarcomere — it is excluded, never approximated.

    The cap's ``mid`` is an arch apex, so its two filaments meet at an angle instead of being anti-parallel.
    Placing a motor there would put the backbone on the angle bisector, so neither head group would walk along
    the actin it holds (bipolar dot −0.600, scale-invariant).  The exclusion must be STRUCTURAL and REPORTED.
    """
    pop = _population(straddle=True)
    by_class = {bundle.sf_class: bundle for bundle in pop.bundles}
    assert by_class["ventral"].motor_station_ready
    assert by_class["dorsal"].motor_station_ready
    assert by_class["transverse_arc"].motor_station_ready
    assert not by_class["perinuclear_cap"].motor_station_ready       # arched ⇒ not a station

    census = pop.motor_station_census()
    assert census["excluded_by_class"] == {"perinuclear_cap": 2}
    assert set(census["ready_by_class"]) == {"ventral", "dorsal", "transverse_arc"}
    assert census["N_stations"] == 4 + 2 + 2                          # the straight bundles only

    nmii = build_sf_straddle_nmii_population(pop, _TOPOLOGY, id_base=_ID_BASE)
    assert nmii.n_minifilaments == census["N_stations"]
    assert nmii.station_census["excluded_by_class"] == {"perinuclear_cap": 2}
    # and the exclusion is visible in the census a report/figure would print.
    assert nmii.census()["station_census"]["excluded_by_class"] == {"perinuclear_cap": 2}


def test_all_curved_population_is_refused_with_a_class_breakdown() -> None:
    """A population of ONLY curved bundles has no motor station at all — that must raise, not place nothing."""
    cap_only = build_sf_arc_population(
        n_ventral=0, n_dorsal=0, n_arc=0, n_cap=2, n_per_fiber=7, **{
            "sarcomere_lateral_um": 2.0 * _TOPOLOGY.head_offset_um,
            "sarcomere_overlap_um": _TOPOLOGY.backbone_length_um})
    with pytest.raises(ValueError, match="every bundle is a CURVED sarcomere"):
        build_sf_straddle_nmii_population(cap_only, _TOPOLOGY, id_base=_ID_BASE)


def test_overlap_shorter_than_the_backbone_is_refused() -> None:
    """REQUIRED-GEOMETRY: a shared span shorter than the backbone leaves outer heads facing only ONE filament."""
    short = build_sf_arc_population(
        n_ventral=2, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=9,
        sarcomere_lateral_um=2.0 * _TOPOLOGY.head_offset_um,
        sarcomere_overlap_um=0.2 * _TOPOLOGY.backbone_length_um)
    with pytest.raises(ValueError, match="SHORTER than the backbone contour"):
        build_sf_straddle_nmii_population(short, _TOPOLOGY, id_base=_ID_BASE)


def test_population_owns_its_own_array_and_a_disjoint_id_block() -> None:
    """Ownership / no-double-count: separate particle storage, disjoint ledger block from sf_arc and cortex."""
    pop = _population(straddle=True)
    nmii = build_sf_straddle_nmii_population(pop, _TOPOLOGY, id_base=_ID_BASE)

    assert nmii.pos.shape[0] == nmii.n_minifilaments * _TOPOLOGY.n_particles
    assert not np.shares_memory(nmii.pos, pop.pos)                    # NOT the sf_arc node array
    assert nmii.head_node.max() < nmii.n_particles                    # head indices address ITS OWN array

    cortex = PopulationLedger("cortex", id_base=0, capacity=70_686)
    cortex.seed_active(range(0, 128))
    nmii.assert_disjoint_from(pop.ledger, cortex)                     # raises if any block/id overlaps

    census = nmii.census()
    assert census["component"] == "nmii"
    assert census["N_minifilaments"] == nmii.n_minifilaments == census["active_count"]
    assert census["N_heads"] == nmii.n_minifilaments * _TOPOLOGY.n_heads
    assert nmii.station_nodes.shape == (nmii.n_minifilaments, 2)


def test_legacy_colocated_population_is_refused() -> None:
    """REQUIRED-GEOMETRY: no fallback placement onto an arbitrary perpendicular."""
    legacy = _population(straddle=False)
    assert not legacy.straddle_ready
    with pytest.raises(ValueError, match="REQUIRED-GEOMETRY"):
        build_sf_straddle_nmii_population(legacy, _TOPOLOGY, id_base=_ID_BASE)


def test_sarcomere_lateral_inconsistent_with_the_motor_is_refused() -> None:
    """A sarcomere not derived from THIS minifilament cannot host it (heads would miss the filament lines)."""
    mismatched = build_sf_arc_population(
        n_ventral=2, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=9,
        sarcomere_lateral_um=4.0 * _TOPOLOGY.head_offset_um,          # twice the separation the heads span
        sarcomere_overlap_um=_TOPOLOGY.backbone_length_um)
    with pytest.raises(ValueError, match="must EQUAL 2·head_offset"):
        build_sf_straddle_nmii_population(mismatched, _TOPOLOGY, id_base=_ID_BASE)

    # and a near-miss is refused too: the derived value is exact, so a "few percent" window would only serve to
    # admit a non-derived number (the no-magic-number rule).
    near = build_sf_arc_population(
        n_ventral=2, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=9,
        sarcomere_lateral_um=2.0 * _TOPOLOGY.head_offset_um * 1.02,
        sarcomere_overlap_um=_TOPOLOGY.backbone_length_um)
    with pytest.raises(ValueError, match="must EQUAL 2·head_offset"):
        build_sf_straddle_nmii_population(near, _TOPOLOGY, id_base=_ID_BASE)


def test_barbed_ward_tangents_follow_the_polarity_flag() -> None:
    """The tangent field points toward each filament's barbed end (the sarcomere's two run OPPOSITE ways)."""
    pop = _population(straddle=True)
    bundle = pop.bundles[0]
    tangent = barbed_ward_tangents(bundle.pos, bundle.fiber_offsets, bundle.polarity)

    assert np.allclose(np.linalg.norm(tangent, axis=1), 1.0)
    n_per_fiber = int(bundle.fiber_offsets[1])
    # fiber A: polarity −1 ⇒ barbed at its FIRST node ⇒ tangent runs node-index-DECREASING.
    step_a = bundle.pos[1] - bundle.pos[0]
    assert float(tangent[0] @ step_a) < 0.0
    # fiber B: polarity +1 ⇒ barbed at its LAST node ⇒ tangent runs node-index-INCREASING.
    step_b = bundle.pos[n_per_fiber + 1] - bundle.pos[n_per_fiber]
    assert float(tangent[n_per_fiber] @ step_b) > 0.0
    # barbed ends point OUTWARD (away from each other) — the precondition for inward contraction.
    assert float(tangent[0] @ tangent[n_per_fiber]) < -0.5


def test_sarcomere_geometry_is_byte_identical_without_the_straddle_arguments() -> None:
    """Additive: omitting the straddle geometry reproduces the landed end-to-end population exactly."""
    a = build_sf_arc_population(n_ventral=3, n_dorsal=2, n_arc=2, n_cap=1, n_per_fiber=7)
    b = build_sf_arc_population(n_ventral=3, n_dorsal=2, n_arc=2, n_cap=1, n_per_fiber=7,
                                sarcomere_lateral_um=None, sarcomere_overlap_um=None)
    assert np.array_equal(a.pos, b.pos)
    assert a.ledger.block == b.ledger.block
    for bundle_a, bundle_b in zip(a.bundles, b.bundles, strict=True):
        assert np.array_equal(bundle_a.myo_i, bundle_b.myo_i)
        assert np.array_equal(bundle_a.myo_j, bundle_b.myo_j)
        assert bundle_a.sarcomere_lateral_um is None and bundle_b.sarcomere_lateral_um is None


def test_straddle_requires_both_geometry_arguments() -> None:
    """Half-supplied geometry is a REQUIRED-PARAM error, never silently ignored."""
    with pytest.raises(ValueError, match="BOTH lateral_um and overlap_um"):
        build_sf_arc_population(n_ventral=1, n_dorsal=0, n_arc=0, n_cap=0,
                                sarcomere_lateral_um=0.4)
    with pytest.raises(ValueError, match="BOTH lateral_um and overlap_um"):
        build_sf_arc_population(n_ventral=1, n_dorsal=0, n_arc=0, n_cap=0,
                                sarcomere_overlap_um=0.3)


def test_overlap_that_does_not_fit_the_sarcomere_is_refused() -> None:
    """A shared span longer than the half-sarcomere is a build error, not a silently clipped geometry."""
    with pytest.raises(ValueError, match="does not fit"):
        build_sf_arc_population(
            n_ventral=1, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=5, cell_radius_um=0.05,
            sarcomere_lateral_um=2.0 * _TOPOLOGY.head_offset_um,
            sarcomere_overlap_um=_TOPOLOGY.backbone_length_um)
