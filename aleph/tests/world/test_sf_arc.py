"""The transverse-arc population — PI decision 8's *"different structure, not a rename"*, checked.

The module's own ``_demo()`` is the primary check and carries the refusals; it runs first here,
because a self-check reachable only through ``python -m`` is a self-check nobody runs.

The four cases added on top are the ones a reader has to be able to see FAIL:

1. the count is not typed anywhere — move the cell and it moves;
2. a derivation does not launder its inputs — a ``PI_GAP`` pitch may not come out ``DERIVED``;
3. the PI rule of 2026-08-20 holds across the axis ranges, not just at one lucky point — this is the
   defect that put 95.5% of the ventral population outside the cell while every structural test passed;
4. what is built is not the ventral bundle: lifted off the substrate and curved, both measured on the
   nodes rather than asserted in a docstring.

No device is touched.  The build refuses a device-less arena, which is one of ``_demo``'s assertions.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.world.bond import BondCount, SourceClass
from aleph.world.build import sf_arc as sf_arc_mod
from aleph.world.build.sf_arc import (
    ARC_PITCH_GAP,
    arc_count_from_footprint,
    host_arc_positions,
    plan_sf_arcs,
)
from aleph.world.build.stress_fiber import plan_stress_fibers
from aleph.world.geometry import assert_inside_membrane, footprint

_AXES = {"pitch_um": 0.6, "lift_um": 0.4, "seg_um": 0.05, "spacing_um": 0.012}
_BUNDLE = BondCount(basis="explicit", value=20.0, scope="transverse arc bundle, generic — band 7-30",
                    source_class=SourceClass.PI_GAP, provenance="no arc-specific bundle count exists")


def _count(fp, pitch_um: float = 0.6, cls: SourceClass = SourceClass.PI_GAP) -> BondCount:
    return arc_count_from_footprint(
        fp, pitch_um=pitch_um, spacing_um=_AXES["spacing_um"], lift_um=_AXES["lift_um"],
        n_filaments_per_arc=20, pitch_scope="MCF7 epithelial on glass",
        pitch_provenance=ARC_PITCH_GAP, pitch_source_class=cls)


def test_sf_arc_self_check_passes() -> None:
    """The module's own ``_demo`` is the primary check."""
    sf_arc_mod._demo()


def test_the_count_is_derived_from_the_cell_and_moves_when_the_cell_does() -> None:
    """A count typed at a call site cannot do this, and that is the 2026-08-20 defect."""
    small, large = footprint(7.5, -7.0), footprint(7.5, -6.0)
    assert large.r_footprint_um > small.r_footprint_um
    assert _count(large).value > _count(small).value

    # And the pitch is held while the count yields — the direction that keeps a sourced quantity sourced.
    assert _count(small, pitch_um=1.2).value < _count(small, pitch_um=0.6).value


def test_a_pi_gap_pitch_may_not_yield_a_derived_count() -> None:
    """A derivation is only as good as its inputs; laundering one is how a blank acquires a label."""
    fp = footprint(7.5, -7.0)
    assert _count(fp, cls=SourceClass.PI_GAP).source_class is SourceClass.PI_GAP
    assert _count(fp, cls=SourceClass.UNRATIFIED_PROXY).source_class is SourceClass.UNRATIFIED_PROXY
    # A sourced spacing is what turns this into a real derivation, and nothing else has to change.
    assert _count(fp, cls=SourceClass.SOURCED).source_class is SourceClass.DERIVED
    assert "ABSENT_FROM_CONTRACT_GRAPH" in _count(fp).provenance


@pytest.mark.parametrize("z_basal_um", [-7.2, -7.0, -6.5, -5.5])
@pytest.mark.parametrize("lift_um", [0.2, 0.4, 1.0])
def test_every_node_lies_inside_the_built_membrane(z_basal_um: float, lift_um: float) -> None:
    """PI 2026-08-20, on the nodes — a partition check says nothing about where a node is."""
    fp = footprint(7.5, z_basal_um)
    axes = _AXES | {"lift_um": lift_um}
    count = arc_count_from_footprint(
        fp, pitch_um=axes["pitch_um"], spacing_um=axes["spacing_um"], lift_um=lift_um,
        n_filaments_per_arc=20, pitch_scope="MCF7", pitch_provenance=ARC_PITCH_GAP,
        pitch_source_class=SourceClass.PI_GAP)
    inv = plan_sf_arcs(fp, count=count, bundle_count=_BUNDLE, **axes)
    pos = host_arc_positions(inv, (0.0, 0.0))

    class _P:
        def __init__(self, lo: int, hi: int) -> None:
            self.nodes = type("R", (), {"lo": lo, "hi": hi})()

    rep = assert_inside_membrane(pos, {"sf_arc": _P(0, pos.shape[0])}, fp.r_cell_um)
    assert rep["sf_arc"]["n_nodes"] == inv["n_strands"] * int(inv["nodes_per_strand"])
    assert rep["sf_arc"]["n_outside"] == 0


def test_what_is_built_is_not_the_ventral_bundle() -> None:
    """Lifted and curved, measured on the nodes. PI decision 8 is about the structure, not the name."""
    fp = footprint(7.5, -7.0)
    inv = plan_sf_arcs(fp, count=_count(fp), bundle_count=_BUNDLE, **_AXES)
    pos = host_arc_positions(inv, (0.0, 0.0))

    # A ventral fibre lies IN the basal plane; every arc node is above it.
    assert pos[:, 2].min() > fp.z_basal_um
    # A ventral fibre is straight; one arc filament has constant radius and a spread in y.
    one = pos[:int(inv["nodes_per_strand"])]
    assert np.ptp(np.hypot(one[:, 0], one[:, 1])) < 1e-9
    assert np.ptp(one[:, 1]) > 1e-3
    # It is a CHORD of the cell, so its end-to-end distance is shorter than its own arc length.
    chord = float(np.linalg.norm(one[-1] - one[0]))
    assert chord < float(inv["arc_length_um_outer"])
    # It spans the leading edge to within its own bundle width — a filament at radius r subtends
    # 2 r sin(half_angle), and r is the centreline plus its hex offset.
    assert abs(chord - float(inv["chord_um"])) <= float(inv["bundle_width_um"])

    # And the two builders answer different questions: neither census key set contains the other's.
    ventral = plan_stress_fibers(
        count=BondCount(basis="explicit", value=4.0, scope="x", source_class=SourceClass.PI_GAP,
                        provenance="fixture"),
        bundle_count=_BUNDLE, length_um=2.0, seg_um=0.05, pitch_um=0.5, spacing_um=0.012,
        sarcomere_um=1.0)
    assert ventral["kind"] == "stress_fiber" and inv["kind"] == "sf_arc"
    assert "half_angle_rad" in inv and "half_angle_rad" not in ventral
    assert "sarcomere_um_declared" in ventral and "sarcomere_um_declared" not in inv
    assert inv["dorsal_sf_built"] is False, "the dorsal fibre is named as unbuilt, never folded in"


def test_the_builder_writes_no_polarity_and_places_no_station() -> None:
    """Both are ratified open questions — PI decision 4 and STATE.md (e) 5 — not build details."""
    fp = footprint(7.5, -7.0)
    inv = plan_sf_arcs(fp, count=_count(fp), bundle_count=_BUNDLE, **_AXES)
    assert "strand_polarity" not in inv
    assert inv["n_stations_placed"] == 0
