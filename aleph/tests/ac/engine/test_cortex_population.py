"""Structural gates for the cortex component's DISJOINT population in cortex-local indices.

These run on the dev Mac at a small filament count.  They gate STRUCTURE only: a passing suite here says
the cortex owns a well-formed, disjoint, connector-conformant population — it says nothing about physics,
which is measured at the native 70,686-filament population on the A5000 (develop on a slice, conclude at
native).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.cortex_population import (
    CORTEX_COMPONENT,
    CORTEX_CONNECTORS,
    CortexEndpointDomain,
    CortexPopulation,
    build_cortex_population,
    cortex_endpoint_domains,
    validate_cortex_connectors,
)
from aleph.engine.population import PopulationLedger, assert_disjoint_populations

#: A development slice: small enough for the Mac, large enough to have crosslinks and multi-node fibers.
SLICE_FILAMENTS = 24
SLICE_DENSITY = 4.0


@pytest.fixture(scope="module")
def population() -> CortexPopulation:
    return build_cortex_population(
        SLICE_FILAMENTS, seg_um=0.5, length_um=3.0, density_per_fil=SLICE_DENSITY, seed=7)


def test_population_tables_are_internally_consistent(population: CortexPopulation) -> None:
    """Segment and bending row counts must follow from the fiber offsets, not from the table lengths."""
    t = population.topology
    assert t.n_fibers == SLICE_FILAMENTS
    assert t.n_segments == t.n_nodes - t.n_fibers
    assert t.n_triples == t.n_nodes - 2 * t.n_fibers
    assert int(t.fiber_offsets[0]) == 0 and int(t.fiber_offsets[-1]) == t.n_nodes
    assert t.n_links > 0, "a cortex with no crosslink is not a load-sharing network"
    census = population.census
    assert census["n_nodes"] == t.n_nodes
    assert census["n_active_filament_ids"] == t.n_fibers
    assert census["n_dormant_filament_ids"] == 0


def test_every_index_is_bounded_by_the_components_own_node_count(population: CortexPopulation) -> None:
    """The falsifiable form of ownership: no cortex table may address beyond the cortex's own nodes."""
    t = population.topology
    n = t.n_nodes
    for table in (t.seg_node_a, t.seg_node_b, t.links, t.bend_triples, t.barbed_node):
        array = np.asarray(table)
        if array.size:
            assert array.min() >= 0 and array.max() < n
    assert np.asarray(t.node_fiber).max() < t.n_fibers


def test_an_index_one_past_the_owned_nodes_is_rejected(population: CortexPopulation) -> None:
    """Negative control for the check above, and the reason it is a LENGTH check and not a 'looks local' one.

    The incumbent places actin at global offset 0, so under the aliased port a cortex index and a global
    index coincide and any "the indices look local" assertion passes vacuously.  What does NOT pass
    vacuously is an index reaching the first NON-cortex node: valid in the global ``n_total`` array,
    out of bounds here.
    """
    import dataclasses

    t = population.topology
    reaching_beyond = np.array(t.seg_node_a, np.int32, copy=True)
    reaching_beyond[0] = t.n_nodes  # the first myosin node in the incumbent's concatenated array
    broken = dataclasses.replace(t, seg_node_a=reaching_beyond)
    with pytest.raises(ValueError, match="reaching into another component"):
        broken.assert_local()


def test_a_truncated_table_is_rejected(population: CortexPopulation) -> None:
    """A silently short segment table must fail against the fiber offsets, not be accepted as smaller."""
    import dataclasses

    t = population.topology
    broken = dataclasses.replace(
        t,
        seg_node_a=np.asarray(t.seg_node_a)[:-1],
        seg_node_b=np.asarray(t.seg_node_b)[:-1],
        seg_polarity=np.asarray(t.seg_polarity)[:-1],
        seg_rest=np.asarray(t.seg_rest)[:-1],
    )
    with pytest.raises(ValueError, match="fiber offsets imply"):
        broken.assert_local()


def test_population_is_disjoint_from_every_other_component(population: CortexPopulation) -> None:
    """Each physical filament belongs to exactly one component — the no-double-count build invariant."""
    sf = PopulationLedger("sf_arc", id_base=1_000_000, capacity=8)
    sf.seed_active(range(1_000_000, 1_000_008))
    ecm = PopulationLedger("ecm", id_base=2_000_000, capacity=4)
    ecm.seed_active(range(2_000_000, 2_000_004))
    assert_disjoint_populations([population.ledger, sf, ecm])
    assert population.ledger.block == (0, SLICE_FILAMENTS)


def test_an_overlapping_id_block_is_rejected(population: CortexPopulation) -> None:
    """Negative control: a cortex block colliding with sf_arc must fail the build, not be reconciled."""
    colliding = build_cortex_population(
        4, seg_um=0.5, length_um=3.0, density_per_fil=SLICE_DENSITY, seed=7, id_base=1_000_000)
    sf = PopulationLedger("sf_arc", id_base=1_000_000, capacity=8)
    sf.seed_active(range(1_000_000, 1_000_008))
    with pytest.raises(ValueError, match="overlapping global"):
        assert_disjoint_populations([colliding.ledger, sf])


def test_declared_endpoint_domains_are_exactly_the_graphs_cortex_edges(
    population: CortexPopulation,
) -> None:
    """No invented connector, none missing — measured against ``reference_cell_architecture()``."""
    edges = validate_cortex_connectors(population)
    architecture = reference_cell_architecture()
    expected = {
        connector.name
        for connector in architecture.connectors
        if CORTEX_COMPONENT in (connector.component_a, connector.component_b)
    }
    assert set(edges) == expected == set(CORTEX_CONNECTORS)
    assert len(expected) == 9, "the T2 seven plus the two contact edges added PI 2026-07-25"


def test_an_invented_connector_is_rejected(population: CortexPopulation) -> None:
    """A component may not declare an edge the composition does not have."""
    with pytest.raises(ValueError, match="not a cortex-terminating connector"):
        CortexEndpointDomain("cortex_wishful_edge", "node", 4, "invented")


def test_a_missing_endpoint_domain_is_rejected(population: CortexPopulation) -> None:
    """Dropping a declared edge must fail loudly: an unaddressable connector is a dead connector."""
    import dataclasses

    reduced = dict(cortex_endpoint_domains(population.topology))
    reduced.pop("nmii_cortex_motor")
    broken = dataclasses.replace(population, endpoints=reduced)
    with pytest.raises(ValueError, match="no declared endpoint domain"):
        validate_cortex_connectors(broken)


def test_the_motor_edge_addresses_segments_and_the_rest_address_nodes(
    population: CortexPopulation,
) -> None:
    """A motor head binds a material coordinate on a polar filament — a segment, never a node."""
    domains = population.endpoints
    assert domains["nmii_cortex_motor"].domain == "segment"
    assert domains["nmii_cortex_motor"].size == population.topology.n_segments
    for name, domain in domains.items():
        if name == "nmii_cortex_motor":
            continue
        assert domain.domain == "node"
        assert domain.size == population.topology.n_nodes


def test_endpoint_domain_rejects_an_index_at_its_upper_bound(population: CortexPopulation) -> None:
    """``contains`` is half-open; the check a connector builder owes the component it binds into."""
    domain = population.endpoints["membrane_erm_cortex"]
    assert domain.contains(np.array([0, domain.size - 1], np.int64))
    assert not domain.contains(np.array([domain.size], np.int64))
    assert domain.contains(np.zeros(0, np.int64)), "an unbound connector offers no indices, not a failure"


def test_geometry_provenance_is_carried_and_names_itself_convenience(
    population: CortexPopulation,
) -> None:
    """The cortex mesh defaults are CONVENIENCE values; the population must say so, not imply sourcing."""
    provenance = population.provenance
    assert "CONVENIENCE" in provenance["geometry_provenance"]
    assert "weave_cell" in provenance["generator"]
    assert "double-count" in provenance["not_a_channel"]


def test_the_builder_defaults_are_the_incumbents_cortex_defaults() -> None:
    """If these drift apart, the "same cortex, private arrays" claim quietly becomes two cortices.

    The native parity gate passes the geometry explicitly on both sides, so a drift would not corrupt a
    measurement — it would make the DEFAULT engine cortex a different object from the default incumbent
    cortex, which is the harder mistake to notice.
    """
    import inspect

    from aleph.components.incumbent.assemble import CellConfig

    defaults = inspect.signature(build_cortex_population).parameters
    incumbent = CellConfig()
    assert defaults["seg_um"].default == incumbent.cortex_seg_um
    assert defaults["length_um"].default == incumbent.cortex_length_um
    assert defaults["density_per_fil"].default == incumbent.cortex_density_per_fil
    assert defaults["arp23_fraction"].default == incumbent.cortex_arp23_fraction
    assert defaults["overlap_free"].default == incumbent.overlap_free_cortex
    assert defaults["overlap_mode"].default == incumbent.cortex_overlap_mode
    assert defaults["overlap_span"].default == incumbent.cortex_overlap_span


def test_the_tables_are_the_generators_own_arrays_not_a_reimplementation() -> None:
    """The cortex is the SAME cortex: its tables must be the weave's own output, only re-typed.

    Host-side parity, runnable without CUDA.  The native gate checks the FORCE fields agree; this checks
    the tables they are built from agree, which is where a transposition or a dropped column would hide.
    """
    from aleph.components.incumbent.assemble import _cortex_region
    from aleph.components.weave.woven_cell import weave_cell

    seed, n = 11, 8
    population = build_cortex_population(
        n, seg_um=0.5, length_um=3.0, density_per_fil=SLICE_DENSITY, seed=seed)
    reference = weave_cell(
        [_cortex_region(n, seg_um=0.5, length_um=3.0, density_per_fil=SLICE_DENSITY)],
        rng=np.random.default_rng(seed), overlap_free=True, overlap_mode="transverse", overlap_span=2,
    )
    crosslinked = reference.to_crosslinked_cortex()
    t = population.topology
    assert np.array_equal(t.pos, crosslinked.net.pos)
    assert np.array_equal(t.fiber_offsets, crosslinked.net.fiber_offsets)
    assert np.array_equal(t.bend_triples, crosslinked.net.bend_triples)
    assert np.array_equal(t.seg_rest, crosslinked.net.seg_rest)
    assert np.array_equal(t.links[:, 0], crosslinked.xl_i)
    assert np.array_equal(t.links[:, 1], crosslinked.xl_j)
    assert np.array_equal(t.link_k, crosslinked.xl_k)
    assert np.array_equal(t.link_r0, crosslinked.xl_rest)
    assert np.array_equal(t.barbed_node, reference.barbed_node)


def test_a_zero_filament_population_is_refused() -> None:
    """Boundary: a component that owns nothing is a configuration error, not an empty success."""
    with pytest.raises(ValueError, match="positive int"):
        build_cortex_population(0)
