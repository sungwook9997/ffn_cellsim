r"""Build the ``cortex`` component's DISJOINT F-actin population in CORTEX-LOCAL indices (host).

WHAT THIS CLOSES.  ``cortex`` is the last actin component that does not own its own state.  Every native
cortex number so far was measured through a *bind-target port*: the incumbent
:func:`aleph.components.incumbent.assemble.build_cell` concatenates ``[actin | myosin | nucleus | membrane]`` into ONE
global ``pos_d``/``f_d`` (``assemble.py:809``), and the engine's cortex endpoint hands that whole array out
under the cortex's name (``scripts/ac_gate_b_cortex_motor_native.py:198,227``: "ALIAS the global array";
``cortex_motor_slice`` §Sanity Gate: "cortex is a bind-target PORT (never a participant)").  Co-location in
one array is exactly what the composition contract says is NEVER a connection (PI 2026-07-22), so while that
holds, **every connector terminating on the cortex is structurally unverifiable** — there is no cortex-owned
force array for an adjoint pair to close against.  Measured on ``reference_cell_architecture()``: **nine** of
the 36 connectors terminate on the cortex (the seven named in
``COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md`` T2 plus ``nucleus_cortex_contact`` and
``membrane_cortex_contact``, added PI 2026-07-25) — see :data:`CORTEX_CONNECTORS`.

This module is the host-side builder for the population a cortex state-owner owns.  It is the cortex
counterpart of :mod:`aleph.engine.sf_population` and deliberately mirrors its shape: generate the
population ONCE, give every filament a globally-unique id from a
:class:`~aleph.engine.population.PopulationLedger`, assert disjointness from every other component, and
declare — never invent — which cortex-owned index domain each connector is allowed to address.

SAME GENERATOR, PRIVATE ARRAYS.  The geometry comes from the SAME call chain the incumbent uses —
``weave_cell(<cortex regions>) -> to_crosslinked_cortex() -> ff.fiber_network`` — so this is not a second
cortex model, it is the same cortex with its own tables.  What changes is the index space: every table here
is indexed in ``[0, n_nodes)`` over cortex nodes ONLY, so a cortex index can no longer reach a myosin,
nucleus or membrane node by construction.  That the incumbent happens to place actin at global offset 0
(making local == global today) is precisely why the invariant this module asserts is the array LENGTH and
the index BOUND, never "the indices look local" — the latter is vacuously true and would gate nothing.

WHAT IS NOT DECIDED HERE.  Which cortex node a given connector binds is the connector's own geometry/kinetics
problem (ERM pairs radially, NMII binds by a device segment query, contact resolves against a live mesh).
Inventing a node subset per edge would be a made-up rule with no source, so this module declares only the
*domain* each edge may address (:class:`CortexEndpointDomain`) and checks the declared edge set against the
architecture.  Binding is the connector builder's job.

engine units: length µm, force pN, stiffness pN/µm.  Host NumPy only — CPU-importable, so the structural
gates drive it on the dev Mac at a small filament count; the native population (70,686 filaments / 494,802
nodes) is built on the gbook A5000 by the CUDA-lane owner in :mod:`aleph.engine.cortex_state`.

Sanity Gate (self-tested in ``tests/ac/engine/test_cortex_population.py``):
  * ownership / no-double-count: every cortex filament draws a unique id from the ``cortex`` block; the
    population is disjoint from an ``sf_arc``/``ecm``/``nmii`` ledger
    (:func:`~aleph.engine.population.assert_disjoint_populations`); ``active_count == n_fibers``.
  * boundary: EVERY index in EVERY table is in ``[0, n_nodes)`` — the falsifiable form of "the cortex cannot
    address another component's node".  Table row counts are cross-checked against the fiber offsets
    (segments ``= N − F``, interior bending triples ``= N − 2F``), so a silently-truncated table fails.
  * conservation: the node count is the sum of the per-fiber node counts and the fiber offsets are
    contiguous ``0..N`` — a gap would mean an orphan node owned by nobody.
  * connector conformance: the declared endpoint domains are exactly the cortex-terminating edges of
    ``reference_cell_architecture()`` — no invented edge, none missing — and each is bidirectional with
    adjoint transfer required.
  * measurement-protocol: the crosslink channel and the bending channel are the SAME two ``ff`` channels
    ``ac.cell.driver._accumulate_all`` launches on the cortex, in the same units, so a native parity check
    against the incumbent is meaningful rather than a comparison of two different models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np
import numpy.typing as npt

from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.population import PopulationLedger

__all__ = [
    "CORTEX_COMPONENT",
    "CORTEX_CONNECTORS",
    "CortexEndpointDomain",
    "CortexPopulation",
    "CortexTopology",
    "build_cortex_population",
    "cortex_endpoint_domains",
    "validate_cortex_connectors",
]

CORTEX_COMPONENT = "cortex"

#: Every connector in ``reference_cell_architecture()`` with the cortex as an endpoint, measured (not
#: recited): the seven of T2 plus the two contact edges added PI 2026-07-25.  Asserted against the
#: architecture by :func:`validate_cortex_connectors`, so this tuple cannot silently drift from the graph.
CORTEX_CONNECTORS: tuple[str, ...] = (
    "membrane_erm_cortex",
    "surface_porous_transfer",
    "sf_cortex_transient",
    "mt_cortex_capture",
    "lamellipodium_cortex_seam",
    "filopodium_cortex_root",
    "nmii_cortex_motor",
    "nucleus_cortex_contact",
    "membrane_cortex_contact",
)

#: Index domains a cortex-terminating connector may address.  ``node`` = a cortex node index (the
#: position/force arrays); ``segment`` = a cortex segment index (a material coordinate on a polar filament,
#: which is what a motor head binds).
NODE_DOMAIN = "node"
SEGMENT_DOMAIN = "segment"


@dataclass(frozen=True, slots=True)
class CortexEndpointDomain:
    """Which cortex-owned index space one connector is allowed to address, and how large it is.

    This is a DECLARATION, not a binding: it says a connector addresses cortex nodes (or cortex segments)
    and how many there are, so a connector builder can be checked for reaching outside the component.  The
    particular node/segment each bond lands on is decided by that connector from live geometry.

    Attributes:
        connector: the connector name; must be one of :data:`CORTEX_CONNECTORS`.
        domain: :data:`NODE_DOMAIN` or :data:`SEGMENT_DOMAIN`.
        size: half-open upper bound of the addressable index range ``[0, size)``.
        rationale: why this edge addresses that domain, in one line (kept so the choice is reviewable).
    """

    connector: str
    domain: str
    size: int
    rationale: str

    def __post_init__(self) -> None:
        if self.connector not in CORTEX_CONNECTORS:
            raise ValueError(f"{self.connector!r} is not a cortex-terminating connector")
        if self.domain not in (NODE_DOMAIN, SEGMENT_DOMAIN):
            raise ValueError(f"unknown cortex endpoint domain {self.domain!r}")
        if isinstance(self.size, bool) or not isinstance(self.size, int) or self.size <= 0:
            raise ValueError(f"{self.connector!r} endpoint domain must have a positive size")
        if not self.rationale.strip():
            raise ValueError(f"{self.connector!r} endpoint domain needs a stated rationale")

    def contains(self, indices: npt.NDArray[np.integer]) -> bool:
        """Whether every index is inside this domain (the check a connector builder owes the component)."""
        idx = np.asarray(indices)
        if idx.size == 0:
            return True
        return bool(idx.min() >= 0 and idx.max() < self.size)


@dataclass(frozen=True, slots=True)
class CortexTopology:
    """The cortex's own node/segment/crosslink/bending tables, ALL indexed in ``[0, n_nodes)``.

    The three force channels these tables feed are exactly the cortex channels the incumbent's
    ``_accumulate_all`` launches: Cytosim bending over ``bend_triples``/``bend_alpha``, the Hookean
    crosslink spring over ``links``/``link_k``/``link_r0``, and — only for a mixed formin+Arp2/3 cortex —
    the 70° branch-angle harmonic over ``branch_triples``/``branch_active``.  There is deliberately NO
    axial-spring channel: cortical inextensibility is an NF2007 constraint carried by the solver, not an
    accumulated force, and adding a spring here would double-count it.

    ``seg_rest`` is the NF2007 per-segment rest length that constraint uses; it is carried so the state
    owner can hand a complete component to the solver without reaching back into the incumbent.

    Attributes:
        pos: ``(N, 3)`` cortex node positions [µm].
        fiber_offsets: ``(F+1,)`` contiguous ``0..N`` node offsets per filament.
        polarity: ``(F,)`` per-filament barbed-end orientation flag (>=0 -> forward along the offsets).
        node_fiber: ``(N,)`` owning filament index of every node (steric same-filament exclusion).
        barbed_node: ``(F,)`` cortex-local barbed-end node of each filament (motor ``walk_dir`` hand-off).
        seg_node_a / seg_node_b / seg_polarity: ``(S,)`` adjacent-node segments + barbed direction sign.
        seg_rest: ``(S,)`` NF2007 segment rest lengths [µm].
        bend_triples: ``(T, 3)`` interior bending triples; ``bend_alpha``: ``(T,)`` α = κ/seg³ [pN/µm].
        links: ``(L, 2)`` crosslink node pairs; ``link_k`` [pN/µm]; ``link_r0`` [µm].
        branch_triples: ``(B, 3)`` Arp2/3 branch-angle triples; ``branch_active``: ``(B,)`` 0/1 mask.
        myosin_site_i / myosin_site_j: ``(M,)`` cortex node pairs the weave offers as NMII seeding sites.
            Offered, NOT owned: the minifilaments themselves belong to the ``nmii`` component and reach the
            cortex only through ``nmii_cortex_motor``.
    """

    pos: npt.NDArray[np.float64]
    fiber_offsets: npt.NDArray[np.int32]
    polarity: npt.NDArray[np.int32]
    node_fiber: npt.NDArray[np.int32]
    barbed_node: npt.NDArray[np.int32]
    seg_node_a: npt.NDArray[np.int32]
    seg_node_b: npt.NDArray[np.int32]
    seg_polarity: npt.NDArray[np.int32]
    seg_rest: npt.NDArray[np.float64]
    bend_triples: npt.NDArray[np.int32]
    bend_alpha: npt.NDArray[np.float64]
    links: npt.NDArray[np.int32]
    link_k: npt.NDArray[np.float64]
    link_r0: npt.NDArray[np.float64]
    branch_triples: npt.NDArray[np.int32]
    branch_active: npt.NDArray[np.int32]
    myosin_site_i: npt.NDArray[np.int32]
    myosin_site_j: npt.NDArray[np.int32]

    @property
    def n_nodes(self) -> int:
        """Number of cortex nodes this component owns."""
        return int(self.pos.shape[0])

    @property
    def n_fibers(self) -> int:
        """Number of cortex filaments this component owns."""
        return int(self.fiber_offsets.shape[0]) - 1

    @property
    def n_segments(self) -> int:
        """Number of cortex segments (the motor-bindable material coordinates)."""
        return int(self.seg_node_a.shape[0])

    @property
    def n_links(self) -> int:
        """Number of cortex crosslinks."""
        return int(self.links.shape[0])

    @property
    def n_triples(self) -> int:
        """Number of interior bending triples."""
        return int(self.bend_triples.shape[0])

    @property
    def n_branch(self) -> int:
        """Number of Arp2/3 branch-angle triples (zero for the formin-only default cortex)."""
        return int(self.branch_triples.shape[0])

    def assert_local(self) -> None:
        """Assert every table is bounded by this component's own node count.

        This is the falsifiable form of component ownership.  It is NOT "the indices look local" — the
        incumbent places actin at global offset 0, so local and global indices coincide today and any such
        check would pass vacuously.  What cannot pass vacuously is that no index reaches beyond ``n_nodes``:
        under the aliased port the addressable array is ``n_total`` long, and a myosin/membrane index is a
        perfectly valid index into it.
        """
        n = self.n_nodes
        f = self.n_fibers
        if n <= 0 or f <= 0:
            raise ValueError("cortex population must own at least one filament and one node")
        off = np.asarray(self.fiber_offsets, np.int64)
        if int(off[0]) != 0 or int(off[-1]) != n:
            raise ValueError("cortex fiber offsets must run contiguously from 0 to n_nodes")
        if np.any(np.diff(off) < 2):
            raise ValueError("every cortex filament must own at least two nodes (one segment)")
        expected_segments = n - f
        expected_triples = max(n - 2 * f, 0)
        if self.n_segments != expected_segments:
            raise ValueError(
                f"cortex segment table has {self.n_segments} rows but the fiber offsets imply "
                f"{expected_segments} — a truncated or concatenated table"
            )
        if self.n_triples != expected_triples:
            raise ValueError(
                f"cortex bending table has {self.n_triples} rows but the fiber offsets imply "
                f"{expected_triples}"
            )
        for label, table in (
            ("node_fiber", self.node_fiber),
            ("barbed_node", self.barbed_node),
            ("seg_node_a", self.seg_node_a),
            ("seg_node_b", self.seg_node_b),
            ("bend_triples", self.bend_triples),
            ("links", self.links),
            ("branch_triples", self.branch_triples),
            ("myosin_site_i", self.myosin_site_i),
            ("myosin_site_j", self.myosin_site_j),
        ):
            array = np.asarray(table)
            if array.size == 0:
                continue
            bound = f if label == "node_fiber" else n
            if int(array.min()) < 0 or int(array.max()) >= bound:
                raise ValueError(
                    f"cortex {label} addresses index {int(array.max())} outside its own "
                    f"[0, {bound}) — the component is reaching into another component's array"
                )
        for label, values, rows in (
            ("link_k", self.link_k, self.n_links),
            ("link_r0", self.link_r0, self.n_links),
            ("bend_alpha", self.bend_alpha, self.n_triples),
            ("seg_rest", self.seg_rest, self.n_segments),
            ("seg_polarity", self.seg_polarity, self.n_segments),
            ("branch_active", self.branch_active, self.n_branch),
            ("polarity", self.polarity, self.n_fibers),
        ):
            if int(np.asarray(values).shape[0]) != int(rows):
                raise ValueError(f"cortex {label} length does not match its table row count")


@dataclass(frozen=True, slots=True)
class CortexPopulation:
    """One built cortex component: its tables, its unique-ID ledger, and its connector domains.

    Attributes:
        topology: the cortex-local tables (:class:`CortexTopology`).
        ledger: the ``cortex`` :class:`~aleph.engine.population.PopulationLedger`, one active id per
            filament, drawn from the block ``[id_base, id_base + n_fibers)``.
        endpoints: connector name -> the cortex index domain that edge may address.
        provenance: build-time provenance strings (region geometry, generator, overlap policy).
    """

    topology: CortexTopology
    ledger: PopulationLedger
    endpoints: Mapping[str, CortexEndpointDomain]
    provenance: Mapping[str, str] = field(default_factory=dict)

    @property
    def census(self) -> dict[str, int]:
        """Standing inventory for the run artifact (never a rounded or 'representative' count)."""
        t = self.topology
        return {
            "n_filaments": t.n_fibers,
            "n_nodes": t.n_nodes,
            "n_segments": t.n_segments,
            "n_crosslinks": t.n_links,
            "n_bending_triples": t.n_triples,
            "n_branch_triples": t.n_branch,
            "n_myosin_sites_offered": int(self.topology.myosin_site_i.shape[0]),
            "n_active_filament_ids": self.ledger.active_count,
            "n_dormant_filament_ids": self.ledger.dormant_count,
        }


def cortex_endpoint_domains(topology: CortexTopology) -> dict[str, CortexEndpointDomain]:
    """Declare the cortex index domain each cortex-terminating connector may address.

    Only two of the nine address segments rather than nodes, and for the same physical reason: a motor
    binds a *material coordinate on a polar filament*, which is a point on a segment, not a node.  The
    other seven address nodes because their endpoint role in ``reference_cell_architecture()`` is a
    "cortex material point" or a surface quadrature over the cortical shell.

    No edge is given a node SUBSET here.  Which node an ERM tether pairs with, which segment a head binds,
    and which quadrature a contact resolves against are decided by those connectors from live geometry;
    narrowing the domain without a source would be an invented rule.
    """
    nodes = topology.n_nodes
    segments = topology.n_segments
    return {
        "membrane_erm_cortex": CortexEndpointDomain(
            "membrane_erm_cortex", NODE_DOMAIN, nodes,
            "ERM tethers a membrane vertex to a cortex material point (radial pairing is the connector's)",
        ),
        "surface_porous_transfer": CortexEndpointDomain(
            "surface_porous_transfer", NODE_DOMAIN, nodes,
            "porous transfer integrates over the cortical shell quadrature, i.e. the owned nodes",
        ),
        "sf_cortex_transient": CortexEndpointDomain(
            "sf_cortex_transient", NODE_DOMAIN, nodes,
            "a transient subcortical crosslink joins an SF node to a cortical actin material point",
        ),
        "mt_cortex_capture": CortexEndpointDomain(
            "mt_cortex_capture", NODE_DOMAIN, nodes,
            "cortical dynein capture acts at a cortex material point",
        ),
        "lamellipodium_cortex_seam": CortexEndpointDomain(
            "lamellipodium_cortex_seam", NODE_DOMAIN, nodes,
            "the rear seam crosslinks lamellipodial nodes to cortical actin material points",
        ),
        "filopodium_cortex_root": CortexEndpointDomain(
            "filopodium_cortex_root", NODE_DOMAIN, nodes,
            "a filopodium bundle roots on cortical actin material points",
        ),
        "nmii_cortex_motor": CortexEndpointDomain(
            "nmii_cortex_motor", SEGMENT_DOMAIN, segments,
            "an NMII head binds a live polar actin material coordinate = a point on a cortex segment",
        ),
        "nucleus_cortex_contact": CortexEndpointDomain(
            "nucleus_cortex_contact", NODE_DOMAIN, nodes,
            "envelope contact resolves against cortex material points",
        ),
        "membrane_cortex_contact": CortexEndpointDomain(
            "membrane_cortex_contact", NODE_DOMAIN, nodes,
            "membrane contact resolves against cortex material points",
        ),
    }


def validate_cortex_connectors(
    population: CortexPopulation,
    architecture: CellArchitecture | None = None,
) -> tuple[str, ...]:
    """Assert the declared endpoint domains are EXACTLY the architecture's cortex-terminating edges.

    Args:
        population: the built cortex population.
        architecture: the composition to check against; defaults to ``reference_cell_architecture()``.

    Returns:
        The cortex-terminating connector names, in architecture order.

    Raises:
        ValueError: if an edge is declared here but absent from the graph (an invented connector), if an
            edge terminates on the cortex but has no declared domain (a silently unaddressable connector),
            or if any such edge is not bidirectional with adjoint transfer required.
    """
    architecture = architecture if architecture is not None else reference_cell_architecture()
    architecture.component(CORTEX_COMPONENT)
    graph_edges = tuple(
        connector.name
        for connector in architecture.connectors
        if CORTEX_COMPONENT in (connector.component_a, connector.component_b)
    )
    declared = set(population.endpoints)
    invented = declared - set(graph_edges)
    if invented:
        raise ValueError(
            f"cortex declares endpoint domains for {sorted(invented)}, which are not connectors of the "
            "architecture — a component may not invent an edge"
        )
    missing = set(graph_edges) - declared
    if missing:
        raise ValueError(
            f"cortex-terminating connectors {sorted(missing)} have no declared endpoint domain, so they "
            "would be unaddressable on the component that owns their far side"
        )
    if set(graph_edges) != set(CORTEX_CONNECTORS):
        raise ValueError(
            f"CORTEX_CONNECTORS {sorted(CORTEX_CONNECTORS)} has drifted from the architecture "
            f"{sorted(graph_edges)}"
        )
    for connector in architecture.connectors:
        if connector.name not in declared:
            continue
        if not connector.bidirectional or not connector.adjoint_transfer_required:
            raise ValueError(
                f"cortex connector {connector.name!r} must be bidirectional with adjoint transfer — a "
                "one-sided scatter is the defect array ownership exists to expose"
            )
    return graph_edges


def build_cortex_population(
    n_filaments: int,
    *,
    seg_um: float = 0.5,
    length_um: float = 3.0,
    density_per_fil: float = 20.0,
    arp23_fraction: float = 0.0,
    arp23_length_um: float = 0.15,
    arp23_seg_um: float = 0.05,
    arp23_mother_fraction: float = 0.2,
    overlap_free: bool = True,
    overlap_mode: str = "transverse",
    overlap_span: int = 2,
    seed: int = 20260722,
    id_base: int = 0,
) -> CortexPopulation:
    """Build the cortex's DISJOINT population + its connector domains in cortex-LOCAL indices (host).

    The geometry generator is the one the incumbent uses — ``ac.weave.weave_cell`` over the cortex region
    specs, then ``to_crosslinked_cortex()`` — so this is the same cortex, not a second model.  The defaults
    reproduce ``ac.cell.assemble.CellConfig``'s cortex defaults, which are CONVENIENCE values, not sourced
    ones (``PARAM_PROVENANCE_AUDIT_2026-07-24.md``): ``seg_um`` 0.5 and ``length_um`` 3.0 are 5–10× and
    3–30× above the sourced cortical mesh and filament lengths, and ``density_per_fil`` 20.0 has no sourced
    value at all — it is the percolation floor that keeps the network single-spanning.  They are repeated
    here so a caller sees them; they are not endorsed by being defaults.

    Args:
        n_filaments: cortex filament count.  Native is 70,686 = 100 µm⁻² × 4π(7.5 µm)²; a smaller count is
            a development slice and no physics conclusion may be drawn from it.
        seg_um / length_um / density_per_fil: cortex mesh geometry (see above).
        arp23_fraction: Arp2/3 fraction of cortical actin BY MASS; 0.0 (default) = formin-only, which
            weaves zero branch triples.
        arp23_length_um / arp23_seg_um / arp23_mother_fraction: the branched sub-population's geometry.
        overlap_free / overlap_mode / overlap_span: build-time interpenetration policy.  The production
            default in the incumbent is ``overlap_free_cortex=True``; note the repo-level default flag is
            still False elsewhere (STATE.md (f)), which is why it is explicit here.
        seed: host RNG seed for the weave.
        id_base: this component's global filament-ID block base.  0 by convention; ``sf_arc`` starts at
            1,000,000, ECM at 2,000,000, NMII at 3,000,000 — so the blocks cannot overlap.

    Returns:
        The built :class:`CortexPopulation`, with :meth:`CortexTopology.assert_local` and the ledger
        invariants already checked.
    """
    # Imported here rather than at module scope: `ac.cell.assemble` pulls the whole incumbent build path
    # (and Warp) for what is, on this side, three region-spec constructors.
    from aleph.components.incumbent.assemble import R_CORTEX_UM, _cortex_region
    from aleph.components.weave.regions import cortex_arp23_region, derive_cortex_arp23_split
    from aleph.components.weave.woven_cell import weave_cell
    from aleph.laws.forces_warp import _per_triple_alpha

    if isinstance(n_filaments, bool) or not isinstance(n_filaments, int) or n_filaments <= 0:
        raise ValueError("cortex filament count must be a positive int")
    if not 0.0 <= arp23_fraction < 1.0:
        raise ValueError("Arp2/3 mass fraction must be in [0, 1)")

    if arp23_fraction <= 0.0:
        regions = [_cortex_region(
            n_filaments, seg_um=seg_um, length_um=length_um, density_per_fil=density_per_fil)]
    else:
        n_formin, n_arp23 = derive_cortex_arp23_split(
            n_filaments, arp23_fraction, length_um, arp23_length_um)
        regions = [
            _cortex_region(
                n_formin, seg_um=seg_um, length_um=length_um, density_per_fil=density_per_fil),
            cortex_arp23_region(
                n_arp23, length_um=arp23_length_um, seg_um=arp23_seg_um,
                mother_fraction=arp23_mother_fraction, R_um=float(R_CORTEX_UM)),
        ]

    wc = weave_cell(
        regions,
        rng=np.random.default_rng(seed),
        overlap_free=overlap_free,
        overlap_mode=overlap_mode,
        overlap_span=overlap_span,
    )
    crosslinked = wc.to_crosslinked_cortex()
    net = crosslinked.net

    # ── the cortex's own tables, derived exactly as `ac.cell.assemble.build_cell:727-752` derives them ──
    pos = np.ascontiguousarray(net.pos, np.float64)
    n_nodes = int(net.n_nodes)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    node_fiber = np.repeat(np.arange(net.n_fibers), np.diff(fiber_off)).astype(np.int32)
    seg_a, seg_b, seg_pol = wc.actin_segment_topology()
    links = (
        np.ascontiguousarray(np.stack([crosslinked.xl_i, crosslinked.xl_j], axis=1), np.int32)
        if crosslinked.xl_i.size
        else np.zeros((0, 2), np.int32)
    )

    topology = CortexTopology(
        pos=pos,
        fiber_offsets=fiber_off,
        polarity=np.ascontiguousarray(wc.polarity, np.int32),
        node_fiber=node_fiber,
        barbed_node=np.ascontiguousarray(wc.barbed_node, np.int32),
        seg_node_a=np.ascontiguousarray(seg_a, np.int32),
        seg_node_b=np.ascontiguousarray(seg_b, np.int32),
        seg_polarity=np.ascontiguousarray(seg_pol, np.int32),
        seg_rest=np.ascontiguousarray(net.seg_rest, np.float64),
        bend_triples=np.ascontiguousarray(net.bend_triples, np.int32),
        bend_alpha=np.ascontiguousarray(_per_triple_alpha(net), np.float64),
        links=links,
        link_k=np.ascontiguousarray(crosslinked.xl_k, np.float64),
        link_r0=np.ascontiguousarray(crosslinked.xl_rest, np.float64),
        branch_triples=np.ascontiguousarray(crosslinked.branch_triples, np.int32).reshape(-1, 3),
        branch_active=np.ascontiguousarray(wc.branch_active.astype(np.int32)).reshape(-1),
        myosin_site_i=np.ascontiguousarray(crosslinked.myo_i, np.int32),
        myosin_site_j=np.ascontiguousarray(crosslinked.myo_j, np.int32),
    )
    topology.assert_local()

    n_fibers = topology.n_fibers
    ledger = PopulationLedger(
        component=CORTEX_COMPONENT,
        id_base=int(id_base),
        capacity=n_fibers,
        n_nodes=n_nodes,
        n_explicit_states=topology.n_links,
    )
    ledger.seed_active(range(int(id_base), int(id_base) + n_fibers))
    ledger.assert_invariants()

    population = CortexPopulation(
        topology=topology,
        ledger=ledger,
        endpoints=cortex_endpoint_domains(topology),
        provenance={
            "generator": "ac.weave.weave_cell -> to_crosslinked_cortex (same chain as ac.cell.assemble)",
            "geometry": (
                f"n_filaments={n_filaments} seg_um={seg_um} length_um={length_um} "
                f"density_per_fil={density_per_fil} arp23_fraction={arp23_fraction}"
            ),
            "geometry_provenance": (
                "CONVENIENCE (coarse/percolation), not sourced — seg_um and length_um sit above the "
                "sourced cortical mesh and filament lengths and density_per_fil has no sourced value "
                "(PARAM_PROVENANCE_AUDIT_2026-07-24.md)"
            ),
            "overlap": f"overlap_free={overlap_free} mode={overlap_mode} span={overlap_span}",
            "seed": str(seed),
            "force_channels": "cytosim bending + Hookean crosslink (+ Arp2/3 branch angle when mixed)",
            "not_a_channel": (
                "axial inextensibility is an NF2007 solver constraint, not an accumulated force — a "
                "spring here would double-count it"
            ),
        },
    )
    validate_cortex_connectors(population)
    return population
