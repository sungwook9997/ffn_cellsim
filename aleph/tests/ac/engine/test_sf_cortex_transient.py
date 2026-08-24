"""`sf_cortex_transient` topology — derived from a sourced radius, and honest when it derives nothing.

The value of this edge is that nothing in it is chosen. These tests pin that: the radius and stiffness
come from the sourced α-actinin record, the rest length is the measured separation, and a geometry that
supplies no pair returns no connector rather than an empty one.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.load_path import ActorRecord, ElementKind, EndpointRole
from aleph.engine.sf_cortex_transient import (
    SF_CORTEX_TRANSIENT,
    build_sf_cortex_transient_specs,
    pair_sf_to_cortex_segments,
)
from aleph.engine.sf_mechanics import ALPHA_ACTININ

_SF_POS = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
_SF_SEG = np.array([[0, 1]])


def _cortex(offset: float) -> tuple[np.ndarray, np.ndarray]:
    return np.array([[0.0, offset, 0.0], [1.0, offset, 0.0]]), np.array([[0, 1]])


def _pair(offset: float):
    pos, seg = _cortex(offset)
    return pair_sf_to_cortex_segments(
        sf_positions=_SF_POS, sf_segments=_SF_SEG, cortex_positions=pos, cortex_segments=seg,
    )


def test_the_radius_comes_from_the_sourced_record_not_from_here() -> None:
    assert _pair(0.03).capture_radius_um == pytest.approx(ALPHA_ACTININ.capture_radius_um)
    assert ALPHA_ACTININ.capture_radius_um == pytest.approx(0.06)


def test_a_segment_inside_the_radius_pairs_and_one_outside_does_not() -> None:
    inside = _pair(0.03)
    assert inside.n_pairs == 1
    assert inside.separation_um[0] == pytest.approx(0.03)

    outside = _pair(0.5)
    assert outside.n_pairs == 0, "beyond the sourced radius must not pair"


def test_just_outside_the_radius_is_still_outside() -> None:
    """The boundary is the sourced number, not a rounded one — 60 nm means 60 nm."""
    assert _pair(0.0599).n_pairs == 1
    assert _pair(0.0601).n_pairs == 0


def test_zero_pairs_yields_zero_specs_rather_than_a_placeholder_joint() -> None:
    """A connector with no joints carries no force. Emitting one would be the vacuity trap again."""
    record = ActorRecord("sf_arc", 0, 0, 0, 0, 1)
    cortex = ActorRecord("cortex", 1, 0, 0, 0, 1)
    specs = build_sf_cortex_transient_specs(_pair(0.5), sf_record=record, cortex_record=cortex)
    assert specs == ()


def test_rest_length_is_the_measured_separation_so_the_joint_starts_force_free() -> None:
    """The r0_bind convention: a placement-derived rest length manufactures an elastic artifact."""
    pairing = _pair(0.042)
    specs = build_sf_cortex_transient_specs(
        pairing, sf_record=ActorRecord("sf_arc", 0, 0, 0, 0, 1),
        cortex_record=ActorRecord("cortex", 1, 0, 0, 0, 1),
    )
    assert len(specs) == 1
    assert specs[0].rest_um == pytest.approx(0.042)
    assert specs[0].rest_um == pytest.approx(float(pairing.separation_um[0]))


def test_stiffness_is_the_sourced_alpha_actinin_value() -> None:
    specs = build_sf_cortex_transient_specs(
        _pair(0.03), sf_record=ActorRecord("sf_arc", 0, 0, 0, 0, 1),
        cortex_record=ActorRecord("cortex", 1, 0, 0, 0, 1),
    )
    assert specs[0].stiffness_pn_per_um == pytest.approx(ALPHA_ACTININ.link_k)
    assert ALPHA_ACTININ.link_k == pytest.approx(4.6e5)


def test_ports_address_segments_with_a_barycentric_coordinate() -> None:
    """`LoadPathJointRuntime` scatters through barycentric weights, so a NODE port would not bind."""
    specs = build_sf_cortex_transient_specs(
        _pair(0.03), sf_record=ActorRecord("sf_arc", 0, 0, 0, 0, 1),
        cortex_record=ActorRecord("cortex", 1, 0, 0, 0, 1),
    )
    spec = specs[0]
    assert spec.edge == SF_CORTEX_TRANSIENT
    assert spec.port_a.component == "sf_arc" and spec.port_b.component == "cortex"
    for port in (spec.port_a, spec.port_b):
        assert port.element_kind is ElementKind.SEGMENT
        assert 0.0 <= port.local_coordinates[0] <= 1.0
        assert isinstance(port.role, EndpointRole)


def test_multiplicity_above_one_refuses_rather_than_inventing_a_density() -> None:
    """More than one crosslink per SF segment needs a sourced density. There is none."""
    with pytest.raises(NotImplementedError, match="sourced"):
        pair_sf_to_cortex_segments(
            sf_positions=_SF_POS, sf_segments=_SF_SEG,
            cortex_positions=_cortex(0.03)[0], cortex_segments=_cortex(0.03)[1],
            max_pairs_per_sf_segment=2,
        )


def test_skew_segments_use_exact_closest_approach_not_midpoint_distance() -> None:
    """A midpoint prefilter that decided pairing would miss crossing segments; it only prefilters."""
    sf_pos = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    cx_pos = np.array([[0.0, -1.0, 0.02], [0.0, 1.0, 0.02]])
    pairing = pair_sf_to_cortex_segments(
        sf_positions=sf_pos, sf_segments=np.array([[0, 1]]),
        cortex_positions=cx_pos, cortex_segments=np.array([[0, 1]]),
    )
    assert pairing.n_pairs == 1
    assert pairing.separation_um[0] == pytest.approx(0.02)
    assert pairing.u_sf[0] == pytest.approx(0.5)
    assert pairing.u_cortex[0] == pytest.approx(0.5)


def test_the_summary_says_the_count_is_placement_derived() -> None:
    payload = _pair(0.03).as_dict()
    assert "NOT a sourced crosslink inventory" in payload["scope"]
    assert "not widened" in payload["scope"]
    assert payload["capture_radius_provenance"].startswith("SOURCED")


def test_binding_without_kinetics_is_refused_because_it_would_be_a_permanent_weld() -> None:
    """The charter forbids a permanent weld by name, and the connector only raises on propose_events.

    So a caller could bind a kinetics-less transient crosslink and simply never ask for events —
    carrying force through a bond that cannot break, with nothing in the record saying so. The builder
    refuses instead. Job 91 is what found this: the wiring ran and the runtime caught it.
    """
    from aleph.engine.load_path import ActorRegistry
    from aleph.engine.sf_cortex_transient import build_sf_cortex_transient_connector

    sf = ActorRecord("sf_arc", 0, 0, 0, 0, 1)
    cortex = ActorRecord("cortex", 1, 0, 0, 0, 1)
    with pytest.raises(ValueError, match="permanent weld"):
        build_sf_cortex_transient_connector(
            pairing=_pair(0.03), sf_actor=object(), cortex_actor=object(),
            registry=ActorRegistry((sf, cortex)), sf_record=sf, cortex_record=cortex,
            kinetics=None,
        )


def test_the_refusal_names_the_rates_as_sourced_so_it_reads_as_a_missing_delegate() -> None:
    """What is missing is a delegate, never a datum — the message must not send anyone hunting a value."""
    from aleph.engine.load_path import ActorRegistry
    from aleph.engine.sf_cortex_transient import build_sf_cortex_transient_connector

    sf = ActorRecord("sf_arc", 0, 0, 0, 0, 1)
    cortex = ActorRecord("cortex", 1, 0, 0, 0, 1)
    with pytest.raises(ValueError, match="sourced"):
        build_sf_cortex_transient_connector(
            pairing=_pair(0.03), sf_actor=object(), cortex_actor=object(),
            registry=ActorRegistry((sf, cortex)), sf_record=sf, cortex_record=cortex, kinetics=None,
        )
