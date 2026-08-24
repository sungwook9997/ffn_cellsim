r"""Build the ``nmii`` minifilament population that STRADDLES a stress fiber's anti-parallel pairs (host).

WHAT THIS CLOSES.  :mod:`aleph.engine.sf_mechanics` made ``sf_arc`` ``KERNEL_BOUND`` for its PASSIVE
rod-cable laws, and said so explicitly: *"The ACTIVE contractile prestress is NOT a lumped ``k_SF`` baked in
here; it is the reaction of the discrete NMII heads, which enters through a separate MOTOR connector (the SF
counterpart of ``nmii_cortex_motor``) … that motor connector is left honestly SEAMED below because it needs the
NMII actuator + a live SF bind-target port."*  This module builds the first of those two missing pieces: the
**explicit head-resolved NMII minifilament population placed on the SF sarcomeres**, host-side.

WHY IT IS NOT A COPY OF THE CORTEX PATH.  ``ac/cell/assemble.py::_build_myosin`` places cortical minifilaments
by SEARCHING for an anti-parallel partner filament with a KD-tree (``_select_antiparallel_partner``), because a
woven cortex only knows its cross-fiber contact pairs.  A stress fiber does not need that search: a sarcomere
**is** the anti-parallel pair, built as one (:func:`aleph.engine.sf_population.build_sf_arc_population`),
so the partner is KNOWN by construction and :meth:`SFArcPopulation.motor_stations` hands it over exactly.  The
placement itself reuses the audited ``ac/motor`` primitives verbatim —
:func:`~aleph.components.motor.minifilament_topology.straddle_frame` for the frame and
:meth:`~aleph.components.motor.minifilament_topology.MinifilamentTopology.placed_positions` for the rigid geometry —
so this is a placement lane, never a second motor model.

ENGINE-LAYER PURE.  Host NumPy + ``ac/motor`` only: **no** ``ac/cell`` import (the frozen incumbent) and **no**
``warp`` import, so the module is CPU-importable and the structural gates drive it directly on the dev Mac.  The
device upload lives in :mod:`aleph.engine.sf_motor_slice`.

SPLIT OWNERSHIP IS REAL HERE.  The cortex GATE-B lane has the actuator and the cortex port ALIAS one global
``cell.pos_d`` (a logical split with global indices — honest, but one array).  These SF minifilaments own their
OWN particle array, disjoint from the SF nodes: the ``nmii`` component's ``position_d``/``force_d`` and
``sf_arc``'s are two never-merged arrays, and the only coupling is the connector's two-array adjoint
crossbridge scatter (Newton's 3rd law).  That is the ownership the composed graph contract asks for.

GEOMETRY IS DERIVED, NEVER CHOSEN (Magic-Number Block).  Nothing here picks a length.  The sarcomere the
minifilament sits in must have been built with the lateral separation and axial overlap DERIVED from this very
minifilament (``lateral = 2·head_offset_um``, ``overlap = backbone_length_um``; see
:func:`build_sf_arc_population`), and :func:`build_sf_straddle_nmii_population` re-asserts that the hosting
geometry is consistent with the topology it is handed — a mismatch raises instead of silently placing heads off
the actin.

MAGNITUDES ARE A PI-GAP (report-not-tune).  Head count per side, stall force, ``k_on``, and every mechanical
stiffness are unresolved KB/PI GAPs (``PI_GAP_EVIDENCE_CARDS`` N1–N9).  This module only places STRUCTURE; the
kinetic/mechanical constants are supplied by the caller through :class:`MinifilamentTopology` and the slice's
required-param config, and are never defaulted here.

engine units: length µm.

STRAIGHT SARCOMERES ONLY — the curved cap is EXCLUDED, not approximated.  A rigid bipolar minifilament can only
straddle two filaments that are ANTI-PARALLEL over the span it occupies.  The ventral / dorsal / transverse-arc
classes are straight (their ``mid`` is exactly the chord midpoint), so their two directions are exactly opposite
and the measured bipolar dot is **−1.0000**.  The perinuclear cap's ``mid`` is an arch APEX, so its two
filaments meet at an angle; a minifilament placed there gets its backbone along the ANGLE BISECTOR and neither
head group walks along the actin it holds (measured bipolar dot **−0.600**, scale-invariant — it does not
improve with a finer discretisation or a different nuclear radius).  Such a placement would be a scaffold
wearing a mechanism's name, so :meth:`SFArcPopulation.motor_stations` excludes curved bundles by construction
and :meth:`SFNMIIPopulation.station_census` reports the exclusion by class.  A curved-sarcomere motor needs its
own placement geometry — a separate design, surfaced to PI rather than silently approximated here.

Sanity Gate (self-tested in tests/ac/engine/test_sf_nmii_population.py):
  * ownership / no-double-count: the minifilament particles live in their OWN array (never an SF node index);
    the population carries its own :class:`PopulationLedger` block, asserted disjoint from ``sf_arc``/cortex.
  * placement (the physical claim): every head's perpendicular distance to the SF filament LINE it is meant to
    sit on is ≈0, and the ``+``/``−`` head groups land on the two DIFFERENT anti-parallel filaments — the
    bipolar straddle, measured, not asserted by construction.
  * sign/boundary: the two head groups' barbed directions are exactly anti-parallel (bipolar dot = −1 for every
    straight class), i.e. the pair contracts inward; curved bundles are excluded and counted, never placed.
  * REQUIRED-GEOMETRY: a population whose sarcomere was NOT built with the straddle overlap raises, and so does
    one whose lateral separation ≠ ``2·head_offset`` or whose overlap is shorter than the backbone contour (no
    fallback to a placement that puts heads off the actin).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.engine.sf_population import SFArcPopulation
from aleph.components.motor.minifilament_topology import MinifilamentTopology, straddle_frame

__all__ = [
    "NMII_COMPONENT",
    "SFNMIIPopulation",
    "barbed_ward_tangents",
    "build_sf_straddle_nmii_population",
]

NMII_COMPONENT = "nmii"

_EPS = 1.0e-12
#: Relative slack allowed between the hosting sarcomere's lateral separation and ``2·head_offset_um``.  It is
#: FLOATING-POINT slack, not a physics window: the separation must EQUAL twice the head offset, because the
#: head↔backbone arm rests at exactly ``head_offset`` — any mismatch δ leaves every head a perpendicular δ/2
#: off its filament line AND pre-strains the arm before a single crossbridge forms.  The correct value is
#: therefore not "within a few percent" but the DERIVED one
#: (:func:`aleph.engine.sf_motor_slice.sf_sarcomere_geometry_for`), which matches to ~1e-16.  A wider
#: window would exist only to admit a non-derived number, which is what the no-magic-number rule forbids.
_LATERAL_MATCH_REL_SLACK = 1.0e-9


def barbed_ward_tangents(
    pos: npt.NDArray[np.float64],
    fiber_offsets: npt.NDArray[np.int64],
    polarity: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Return the per-node unit tangent pointing toward each node's filament BARBED end.

    Central difference inside a filament, one-sided at its two ends, signed by the filament's barbed flag
    (``polarity >= 0`` ⇒ barbed at the LAST node ⇒ the tangent runs node-index-increasing; ``< 0`` ⇒ flipped).
    This is the SF-side analogue of the cortex tangent field, written here so the engine layer does not import
    the frozen ``ac/cell`` helper; it is pure geometry (no parameter, nothing to tune).

    Args:
        pos: ``(N, 3)`` node positions [µm] of one contiguous multi-filament block.
        fiber_offsets: ``(F+1,)`` contiguous per-filament node offsets into ``pos``.
        polarity: ``(F,)`` per-filament barbed-end flag.

    Returns:
        ``(N, 3)`` unit tangents; a degenerate (zero-length) filament segment yields a zero row.
    """
    pos = np.asarray(pos, np.float64)
    offsets = np.asarray(fiber_offsets, np.int64)
    tangent = np.zeros_like(pos)
    for f in range(int(offsets.size) - 1):
        lo, hi = int(offsets[f]), int(offsets[f + 1])
        if hi - lo < 2:
            continue
        block = pos[lo:hi]
        diff = np.zeros_like(block)
        diff[1:-1] = block[2:] - block[:-2]
        diff[0] = block[1] - block[0]
        diff[-1] = block[-1] - block[-2]
        norm = np.linalg.norm(diff, axis=1, keepdims=True)
        diff = np.divide(diff, norm, out=np.zeros_like(diff), where=norm > _EPS)
        tangent[lo:hi] = diff * (1.0 if int(polarity[f]) >= 0 else -1.0)
    return tangent


@dataclass(slots=True)
class SFNMIIPopulation:
    """The ``nmii`` minifilament population straddling a stress-fiber sarcomere set (host NumPy).

    All index arrays are LOCAL to :attr:`pos` — this population owns its own particle array, which is a
    different array from the ``sf_arc`` node array (split ownership; the crossbridge connector is the only
    coupling).  ``head_node`` indexes :attr:`pos`; ``station_nodes`` records, for provenance/diagnostics only,
    which ``sf_arc`` node pair each minifilament was placed between.
    """

    pos: npt.NDArray[np.float64]                  # (P,3) minifilament particle positions [µm]
    backbone_bonds: npt.NDArray[np.int32]         # (M·(n_bb-1), 2) backbone rod pairs
    head_bonds: npt.NDArray[np.int32]             # (M·2H, 2) head↔backbone arm pairs
    backbone_angles: npt.NDArray[np.int32]        # (M·(n_bb-2), 3) backbone bending triples (rest π)
    head_arm_angles: npt.NDArray[np.int32]        # (M·2H, 3) head-arm orientation triples (rest π/2)
    head_node: npt.NDArray[np.int32]              # (M·2H,) head particle index into pos
    head_side: npt.NDArray[np.int32]              # (M·2H,) 0 = '+' side, 1 = '−' side
    minifilament_offset: npt.NDArray[np.int32]    # (M+1,) first particle index of each minifilament
    station_nodes: npt.NDArray[np.int64]          # (M,2) the sf_arc node pair each minifilament straddles
    frames: npt.NDArray[np.float64]               # (M,3,3) [centre, e_x, e_y] per minifilament (diagnostics)
    ledger: PopulationLedger
    topology: MinifilamentTopology
    #: Which SF bundles could host a motor and which were EXCLUDED, by class (never a silent truncation).
    station_census: dict[str, object] = field(default_factory=dict)

    @property
    def n_minifilaments(self) -> int:
        return int(self.minifilament_offset.size) - 1

    @property
    def n_heads(self) -> int:
        return int(self.head_node.size)

    @property
    def n_particles(self) -> int:
        return int(self.pos.shape[0])

    def census(self) -> dict[str, object]:
        """Unique-active NMII census (each minifilament owned by exactly one component, no double count)."""
        return {
            "component": self.ledger.component,
            "N_minifilaments": self.n_minifilaments,
            "N_heads": self.n_heads,
            "N_particles": self.n_particles,
            "N_heads_per_side": int(self.topology.n_heads_per_side),
            "n_backbone_beads": int(self.topology.n_bb),
            "active_count": self.ledger.active_count,
            "dormant_count": self.ledger.dormant_count,
            "backbone_length_um": float(self.topology.backbone_length_um),
            "head_offset_um": float(self.topology.head_offset_um),
            "station_census": dict(self.station_census),
        }

    def assert_disjoint_from(self, *others: PopulationLedger) -> None:
        """Assert this population shares no global ID block or active id with ``others`` (sf_arc, cortex, …)."""
        assert_disjoint_populations([self.ledger, *others])


def build_sf_straddle_nmii_population(
    population: SFArcPopulation,
    topology: MinifilamentTopology,
    *,
    id_base: int,
) -> SFNMIIPopulation:
    """Place one bipolar NMII minifilament on every SF sarcomere station (host NumPy).

    For each station — the KNOWN anti-parallel node pair from
    :meth:`~aleph.engine.sf_population.SFArcPopulation.motor_stations` — the frame is taken from
    :func:`~aleph.components.motor.minifilament_topology.straddle_frame` (backbone along the shared filament line,
    ``+y`` toward filament A, centre at the pair midpoint) and the rigid particle layout from
    :meth:`MinifilamentTopology.placed_positions`.  Bond/angle topology comes from the same
    :class:`MinifilamentTopology`, offset into this population's own contiguous particle array.  Nothing is
    searched, sampled, or tuned.

    Args:
        population: a built ``sf_arc`` population whose sarcomeres carry the STRADDLE overlap geometry.
        topology: the minifilament reference topology (its ``head_offset_um`` / ``backbone_length_um`` are the
            values the hosting sarcomere geometry must have been derived from).
        id_base: inclusive lower bound of this population's global ``nmii`` ID block (must not overlap the
            ``sf_arc`` / cortex / ECM blocks — asserted by the caller via
            :func:`~aleph.engine.population.assert_disjoint_populations`).

    Returns:
        The built :class:`SFNMIIPopulation`.

    Raises:
        ValueError: if the population has no sarcomere station, if it was NOT built with the straddle overlap
            geometry (the co-located legacy sarcomere cannot be straddled — heads would land on an arbitrary
            perpendicular), or if its lateral separation is inconsistent with ``2·head_offset_um``.
    """
    if not population.bundles:
        raise ValueError("sf_arc population has no bundles to host an NMII minifilament")
    if not population.straddle_ready:
        raise ValueError(
            "REQUIRED-GEOMETRY: the sf_arc population was built with the LEGACY end-to-end sarcomere (its two "
            "anti-parallel filaments are co-located at mid), which a bipolar minifilament cannot straddle — "
            "rebuild it with build_sf_arc_population(sarcomere_lateral_um=2*head_offset_um, "
            "sarcomere_overlap_um=backbone_length_um).  No fallback placement is offered: placing the heads on "
            "an arbitrary perpendicular is the defect this lane exists to avoid."
        )
    expected_lateral = 2.0 * float(topology.head_offset_um)
    min_overlap = float(topology.backbone_length_um)
    for bundle in population.bundles:
        lateral = float(bundle.sarcomere_lateral_um)          # not None: straddle_ready checked above
        overlap = float(bundle.sarcomere_overlap_um)
        if abs(lateral - expected_lateral) > _LATERAL_MATCH_REL_SLACK * expected_lateral:
            raise ValueError(
                f"hosting sarcomere lateral separation {lateral!r} µm must EQUAL 2·head_offset = "
                f"{expected_lateral!r} µm: every head would otherwise sit a perpendicular "
                f"{abs(lateral - expected_lateral) / 2.0:.4g} µm off its filament line and its arm would be "
                f"pre-strained before any crossbridge forms.  Derive the population geometry from THIS "
                f"topology (sf_sarcomere_geometry_for)."
            )
        if overlap < min_overlap * (1.0 - _LATERAL_MATCH_REL_SLACK):
            raise ValueError(
                f"hosting sarcomere axial overlap {overlap!r} µm is SHORTER than the backbone contour "
                f"{min_overlap!r} µm it must contain: the heads are distributed along the whole backbone, so "
                f"those beyond the shared span would face only ONE filament (or none) and the bipolar dipole "
                f"would be unbalanced.  Derive it from THIS topology (sf_sarcomere_geometry_for)."
            )

    stations = population.motor_stations()
    if stations.shape[0] == 0:
        census = population.motor_station_census()
        raise ValueError(
            "sf_arc population exposes no motor station: every bundle is a CURVED sarcomere, which a rigid "
            f"bipolar minifilament cannot straddle ({census['excluded_by_class']}).  Include at least one "
            "straight class (ventral / dorsal / transverse_arc)."
        )

    # barbed-ward tangents over the flat population (per-bundle blocks, each with its own offsets/polarity).
    tangent = np.zeros_like(np.asarray(population.pos, np.float64))
    for bundle in population.bundles:
        base = int(bundle.node_base)
        n_local = int(bundle.pos.shape[0])
        tangent[base:base + n_local] = barbed_ward_tangents(
            bundle.pos, bundle.fiber_offsets, bundle.polarity)

    n_bb = int(topology.n_bb)
    n_part = int(topology.n_particles)
    n_side = int(topology.n_heads_per_side)
    n_heads_mf = int(topology.n_heads)
    bb_bonds_local = np.column_stack([np.arange(n_bb - 1), np.arange(1, n_bb)]).astype(np.int64)
    head_bonds_local = topology.head_backbone_bonds()
    bb_angles_local = topology.backbone_angle_triples()
    arm_angles_local = topology.head_arm_angle_triples()
    head_node_local = np.arange(n_bb, n_part, dtype=np.int64)
    side_local = np.concatenate([np.zeros(n_side, np.int32), np.ones(n_side, np.int32)])

    pos_parts: list[npt.NDArray[np.float64]] = []
    bb_parts: list[npt.NDArray[np.int64]] = []
    hb_parts: list[npt.NDArray[np.int64]] = []
    ba_parts: list[npt.NDArray[np.int64]] = []
    ha_parts: list[npt.NDArray[np.int64]] = []
    hn_parts: list[npt.NDArray[np.int64]] = []
    frames: list[npt.NDArray[np.float64]] = []
    kept_stations: list[tuple[int, int]] = []
    for (node_a, node_b) in stations:
        ia, ib = int(node_a), int(node_b)
        centre, e_x, e_y = straddle_frame(
            population.pos[ia], tangent[ia], population.pos[ib], tangent[ib])
        if float(np.linalg.norm(e_x)) <= _EPS or float(np.linalg.norm(e_y)) <= _EPS:
            continue                                          # degenerate station: skip, never place blindly
        base = len(pos_parts) * n_part
        pos_parts.append(topology.placed_positions(centre, e_x, e_y))
        bb_parts.append(bb_bonds_local + base)
        hb_parts.append(head_bonds_local + base)
        ba_parts.append(bb_angles_local + base)
        ha_parts.append(arm_angles_local + base)
        hn_parts.append(head_node_local + base)
        frames.append(np.vstack([centre, e_x, e_y]))
        kept_stations.append((ia, ib))

    if not pos_parts:
        raise ValueError("every sf_arc motor station was geometrically degenerate; nothing was placed")

    n_mf = len(pos_parts)
    ledger = PopulationLedger(component=NMII_COMPONENT, id_base=int(id_base), capacity=n_mf)
    ledger.seed_active(range(int(id_base), int(id_base) + n_mf))
    pos_all = np.ascontiguousarray(np.vstack(pos_parts), dtype=np.float64)
    ledger.n_nodes = int(pos_all.shape[0])
    ledger.n_heads = n_mf * n_heads_mf
    ledger.assert_invariants()

    return SFNMIIPopulation(
        pos=pos_all,
        backbone_bonds=np.ascontiguousarray(np.vstack(bb_parts), dtype=np.int32),
        head_bonds=np.ascontiguousarray(np.vstack(hb_parts), dtype=np.int32),
        backbone_angles=(np.ascontiguousarray(np.vstack(ba_parts), dtype=np.int32)
                         if bb_angles_local.shape[0] else np.zeros((0, 3), np.int32)),
        head_arm_angles=np.ascontiguousarray(np.vstack(ha_parts), dtype=np.int32),
        head_node=np.ascontiguousarray(np.concatenate(hn_parts), dtype=np.int32),
        head_side=np.ascontiguousarray(np.tile(side_local, n_mf), dtype=np.int32),
        minifilament_offset=np.arange(n_mf + 1, dtype=np.int32) * np.int32(n_part),
        station_nodes=np.asarray(kept_stations, np.int64).reshape(-1, 2),
        frames=np.asarray(frames, np.float64).reshape(-1, 3, 3),
        ledger=ledger,
        topology=topology,
        station_census=population.motor_station_census(),
    )
