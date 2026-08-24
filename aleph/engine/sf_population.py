r"""Build the ``sf_arc`` component's DISJOINT F-actin population + its canonical connector wiring (host).

WHAT THIS CLOSES.  The engine's :mod:`aleph.engine.stress_fiber` facade orchestrates SF/arc mechanics
and its four common edges, but it takes the SF state (device arrays, connectors) as GIVEN — nothing built the
population it orchestrates.  The whole-cell walkthrough's *biggest structural gap*
([[project-filament-subsystem-plan]] G1/G2: "SF not in the whole-cell loop") is exactly this: a stress fiber
is a **pre-stressed contractile load path that OWNS its filaments and couples to the rest of the cell only
through explicit connectors**.  This module is the host-side builder that generates that population.

SF-SEPARATE (PI 2026-07-22, [[project-ac-engine-canonical-component-connector]]).  This supersedes the
label-blind ``ac.weave.stress_fiber`` "SF = emergent partition of one WovenCell" approach FOR THE ENGINE.
``sf_arc`` is a **separate state-owning component owning a DISJOINT filament population**: every SF/arc/cap
filament belongs to ``sf_arc`` and to NO other component (no shared node with the cortex, no permanent weld).
Unique global filament IDs come from a :class:`~aleph.engine.population.PopulationLedger`; disjointness
from the cortex (and every other component) is a build-time assertion
(:func:`~aleph.engine.population.assert_disjoint_populations`), promoting the viz-time no-double-count
check to a build invariant.

EMERGENT, NOT LUMPED (CLAUDE.md hard rule; [[project-unified-actin-architecture]]).  The axial prestress is
NEVER a lumped ``k_SF`` bundle stiffness.  It is the load-path integral of the DISCRETE NMII head reactions
along each fiber (each head walks toward its actin barbed end -> Newton-3rd reaction on the actin; an
antiparallel bipolar minifilament -> net inward contraction), computed by the vetted
:func:`aleph.components.weave.stress_fiber.stress_fiber_load_path` (a cumulative sum, no bundle spring).  This
builder OWNS the population and the wiring; it reuses that oracle read-only for the per-SF load path.

MAGNITUDE IS A GAP — report-not-tune (CLAUDE.md; [[project-sf-nmii-forcescale-result]]).  Per-head stall
(I0-B3) and engaged-head density (I0-B6) are unresolved; ``f_head_pN`` defaults to ``None`` and the load path
returns its magnitude-independent SHAPE.  Filament COUNTS here are a small structural population for the Mac
gates, NOT the native SF inventory — the sourced count/geometry distribution arrives from
:class:`~aleph.engine.cell_state.CellState` in a biology phase, never chosen here to hit a band.

FIREWALL.  The emergent taxonomy (ventral / dorsal / transverse-arc / perinuclear-cap) is a build-time
GEOMETRY-SEEDING descriptor — WHERE a fiber is laid down and WHICH anchor kind it carries — never an input to
the force.  The prestress and traction dipole are pure functions of geometry + motors + anchors; a test
permutes the class tag and asserts the load paths are unchanged (the same "never read the label" firewall the
weave world enforces on ``region_id``).

engine units: length um, force pN, stiffness pN/um.  Host NumPy only (no Warp import) — CPU-importable so the
structural gates drive it directly on the dev Mac; the native force assembly (the I3 ``ac.motor`` Hill-NMII
kernels, on the gbook A5000) is the lead's and is NOT invoked here.

Sanity Gate (self-tested in tests/ac/engine/test_sf_population.py):
  * ownership / no-double-count: every fiber gets a unique global id from the ``sf_arc`` block; the population
    is disjoint from a cortex ledger (``assert_disjoint_populations``); ``N_unique == N_summed``.
  * connector conformance: the endpoint sets feed exactly the canonical ``fa_actin_anchor`` /
    ``sf_cortex_transient`` / ``actin_cap_linc`` / ``dorsal_arc_crosslink`` edges of
    ``reference_cell_architecture`` — each bidirectional + adjoint (co-location is never a connection).
  * emergent prestress: a ventral (both-ends-FA) SF is a contractile CLOSED dipole (residual -> 0); a dorsal
    (one-end-FA) SF is OPEN (network carries the free end); prestress comes from the NMII cumsum, no ``k_SF``.
  * anchor taxonomy: anchor_kind emerges from the anchor node set (fa / linc / none), not from the class tag.
  * firewall: load paths are invariant to a class-tag permutation.
  * GAP guard: ``f_head_pN=None`` -> ``tension is None`` and the GAP status (never a tuned magnitude).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.engine.contracts import CellArchitecture, ConnectorScope
from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.components.weave.stress_fiber import (
    GAP_STATUS,
    OK_STATUS,
    StressFiberLoadPath,
    stress_fiber_load_path,
)

__all__ = [
    "SF_COMPONENT",
    "SFClass",
    "SFBundle",
    "SFArcPopulation",
    "FA_ACTIN_ANCHOR",
    "SF_CORTEX_TRANSIENT",
    "ACTIN_CAP_LINC",
    "DORSAL_ARC_CROSSLINK",
    "build_sf_arc_population",
    "validate_sf_connectors",
]

SF_COMPONENT = "sf_arc"

#: canonical connector names this population feeds (all declared in ``reference_cell_architecture``).
FA_ACTIN_ANCHOR = "fa_actin_anchor"          # sf_arc <-> focal_adhesion (ventral/dorsal basal end)
SF_CORTEX_TRANSIENT = "sf_cortex_transient"  # sf_arc <-> cortex        (transient subcortical crosslink)
ACTIN_CAP_LINC = "actin_cap_linc"            # sf_arc <-> nucleus       (perinuclear cap LINC)
DORSAL_ARC_CROSSLINK = "dorsal_arc_crosslink"  # sf_arc <-> sf_arc (internal: dorsal free end <-> arc)


class SFClass:
    """Emergent SF taxonomy — a build-time GEOMETRY/ANCHOR descriptor, NEVER read to compute force.

    The class fixes WHERE a fiber is seeded and WHICH anchor kind it carries; the contractile prestress still
    emerges from motors + anchors (firewall-checked).  ``VENTRAL`` = basal, both ends FA (closed dipole);
    ``DORSAL`` = one basal FA end + one dorsal free end (opens into an arc joint); ``TRANSVERSE_ARC`` = elevated,
    no FA, couples to dorsal free ends + cortex; ``PERINUCLEAR_CAP`` = arches over the nucleus, LINC-anchored.
    """

    VENTRAL = "ventral"
    DORSAL = "dorsal"
    TRANSVERSE_ARC = "transverse_arc"
    PERINUCLEAR_CAP = "perinuclear_cap"

    ALL = (VENTRAL, DORSAL, TRANSVERSE_ARC, PERINUCLEAR_CAP)


_EPS = 1.0e-12


@dataclass(slots=True)
class SFBundle:
    """One built stress fiber: its local geometry + motors + anchors + the two global filament IDs it owns.

    A bundle is an antiparallel two-filament sarcomere unit (barbed ends outward, one bipolar NMII seed at the
    interior overlap) — the smallest unit that carries a self-consistent contractile load path.  ``fa_local`` /
    ``linc_local`` are the anchor node indices *within this bundle*; ``cortex_local`` are the mid-shaft nodes
    offered to the transient cortex connector.  ``global_fiber_ids`` are the unique ``sf_arc``-block ids.
    """

    sf_class: str
    pos: npt.NDArray[np.float64]                 # (n,3) local node positions [um]
    fiber_offsets: npt.NDArray[np.int64]         # (F+1,) local contiguous fiber offsets
    polarity: npt.NDArray[np.int64]              # (F,) barbed-end flag per fiber (+1 last / -1 first)
    myo_i: npt.NDArray[np.int64]                 # (M,) NMII seed endpoint local node index
    myo_j: npt.NDArray[np.int64]                 # (M,) partner local node index (antiparallel overlap)
    fa_local: npt.NDArray[np.int64]              # substrate-anchor local nodes (fa_actin_anchor endpoints)
    linc_local: npt.NDArray[np.int64]            # nuclear-anchor local nodes (actin_cap_linc endpoints)
    cortex_local: npt.NDArray[np.int64]          # mid-shaft local nodes (sf_cortex_transient endpoints)
    free_local: npt.NDArray[np.int64]            # dorsal free-end local nodes (feed the internal arc joint)
    global_fiber_ids: npt.NDArray[np.int64]      # (F,) unique sf_arc-block filament ids
    node_base: int = 0                           # this bundle's first node index in the flat population
    #: straddle geometry actually built (``None`` for the legacy end-to-end sarcomere): the anti-parallel
    #: lateral separation [µm] and the axial overlap length [µm] the two filaments share.  Both are DERIVED
    #: from the NMII minifilament the bundle must host (see :func:`build_sf_arc_population`), never chosen.
    sarcomere_lateral_um: float | None = None
    sarcomere_overlap_um: float | None = None
    #: True only for a STRAIGHT sarcomere, whose two filaments are exactly anti-parallel over the shared span
    #: (``t_a`` and ``t_b`` collinear) — the geometric precondition a bipolar minifilament's placement frame
    #: assumes.  A CURVED sarcomere (the perinuclear cap, whose ``mid`` is an arch apex rather than the chord
    #: midpoint) is deliberately NOT a motor station: see :attr:`motor_station_ready`.
    motor_station_ready: bool = False

    @property
    def straddle_ready(self) -> bool:
        """True iff this bundle was built with the anti-parallel lateral overlap an NMII straddle needs."""
        return self.sarcomere_lateral_um is not None and self.sarcomere_overlap_um is not None

    def load_path(self, *, f_head_pN: float | None = None, n_bins: int = 24) -> StressFiberLoadPath:
        """Emergent contractile load path of this bundle (NMII cumsum; magnitude GAP-guarded)."""
        return stress_fiber_load_path(
            self.pos, self.fiber_offsets, self.polarity, self.myo_i, self.myo_j,
            self.fa_local, self.linc_local, f_head_pN=f_head_pN, n_bins=n_bins)


@dataclass(slots=True)
class SFArcPopulation:
    """The built ``sf_arc`` component: its disjoint filament population, ledger, and connector endpoint sets.

    All node indices in :meth:`connector_endpoints` are GLOBAL (into :attr:`pos`).  The population is generated
    ONCE; there is no resample-to-density path (the ledger enforces fixed capacity).
    """

    bundles: list[SFBundle]
    ledger: PopulationLedger
    pos: npt.NDArray[np.float64]                       # (N,3) flat population node positions [um]
    fa_sites: npt.NDArray[np.int64]                    # global FA-anchor nodes (sf side of fa_actin_anchor)
    linc_sites: npt.NDArray[np.int64]                  # global LINC nodes (sf side of actin_cap_linc)
    cortex_sites: npt.NDArray[np.int64]                # global subcortical nodes (sf side of sf_cortex_transient)
    dorsal_arc_joints: npt.NDArray[np.int64]           # (J,2) internal (dorsal-free, arc) global node pairs
    f_head_pN: float | None = None
    magnitude_status: str = field(default=GAP_STATUS)

    @property
    def n_fibers(self) -> int:
        return int(sum(b.global_fiber_ids.size for b in self.bundles))

    @property
    def n_nodes(self) -> int:
        return int(self.pos.shape[0])

    def connector_endpoints(self) -> dict[str, npt.NDArray[np.int64]]:
        """Return the ``sf_arc``-side endpoint node sets for each canonical connector this population feeds.

        These are the ONLY mechanical couplings of ``sf_arc`` to the rest of the cell — co-location in a shared
        array is never a connection.  The FA/cortex/nucleus sides are owned by ``focal_adhesion`` / ``cortex`` /
        ``nucleus``; here we expose the SF material points that participate.
        """
        return {
            FA_ACTIN_ANCHOR: self.fa_sites,
            SF_CORTEX_TRANSIENT: self.cortex_sites,
            ACTIN_CAP_LINC: self.linc_sites,
            DORSAL_ARC_CROSSLINK: self.dorsal_arc_joints,
        }

    def load_paths(self, *, f_head_pN: float | None = None, n_bins: int = 24) -> list[StressFiberLoadPath]:
        """Per-bundle emergent load path (NMII cumsum, GAP-guarded) — the SF prestress entering the loop."""
        return [b.load_path(f_head_pN=f_head_pN, n_bins=n_bins) for b in self.bundles]

    @property
    def straddle_ready(self) -> bool:
        """True iff EVERY bundle carries the anti-parallel lateral overlap an NMII straddle placement needs."""
        return bool(self.bundles) and all(b.straddle_ready for b in self.bundles)

    def motor_stations(self) -> npt.NDArray[np.int64]:
        """Return the ``(M, 2)`` GLOBAL anti-parallel node pairs a bipolar NMII minifilament straddles.

        One row per motor station: ``(node on filament A, node on the anti-parallel filament B)`` in flat
        :attr:`pos` indices.  These are the geometric seeds the ``nmii`` actuator's straddle placement consumes —
        the SF analogue of the cortex cross-fiber motor pairs, except the anti-parallel partner is KNOWN by
        construction here (a sarcomere is built as the pair), so no KD-tree partner search is needed.

        Only bundles with :attr:`SFBundle.motor_station_ready` contribute.  A CURVED sarcomere (the perinuclear
        cap) is excluded by construction: its two filaments are not collinear over the shared span, so a rigid
        bipolar minifilament placed there would have its backbone along the ANGLE BISECTOR and neither head
        group would walk along the actin it holds (measured bipolar dot −0.600 instead of −1).  Placing a motor
        there anyway would be a scaffold wearing a mechanism's name; a curved-sarcomere motor needs its own
        placement geometry, which is a separate design (surfaced to PI, not silently approximated here).
        """
        rows = [(int(b.node_base + b.myo_i[m]), int(b.node_base + b.myo_j[m]))
                for b in self.bundles if b.motor_station_ready for m in range(int(b.myo_i.size))]
        return np.asarray(rows, np.int64).reshape(-1, 2)

    def motor_station_census(self) -> dict[str, object]:
        """Report which bundles can host a motor and which are excluded, BY CLASS (no silent truncation)."""
        ready: dict[str, int] = {}
        excluded: dict[str, int] = {}
        for bundle in self.bundles:
            target = ready if bundle.motor_station_ready else excluded
            target[bundle.sf_class] = target.get(bundle.sf_class, 0) + 1
        return {
            "N_stations": int(self.motor_stations().shape[0]),
            "ready_by_class": ready,
            "excluded_by_class": excluded,
            "exclusion_reason": "curved sarcomere (mid is an arch apex, not the chord midpoint): the two "
                                "filaments are not collinear over the shared span, so a rigid bipolar "
                                "minifilament cannot straddle them",
        }

    def census(self) -> dict[str, object]:
        """Unique-active SF census — no cortex/SF double count (each fiber owned by exactly one component)."""
        all_ids = np.concatenate([b.global_fiber_ids for b in self.bundles]) if self.bundles \
            else np.zeros(0, np.int64)
        by_class: dict[str, int] = {}
        for b in self.bundles:
            by_class[b.sf_class] = by_class.get(b.sf_class, 0) + int(b.global_fiber_ids.size)
        return {
            "component": self.ledger.component,
            "N_bundles": len(self.bundles),
            "N_unique_sf_fibers": int(np.unique(all_ids).size),
            "N_sf_fibers_summed": int(all_ids.size),          # == unique iff ids truly partition (asserted)
            "N_by_class": by_class,
            "active_count": self.ledger.active_count,
            "dormant_count": self.ledger.dormant_count,
            "N_fa_sites": int(self.fa_sites.size),
            "N_linc_sites": int(self.linc_sites.size),
            "N_cortex_sites": int(self.cortex_sites.size),
            "N_dorsal_arc_joints": int(self.dorsal_arc_joints.shape[0]),
            "magnitude_status": self.magnitude_status,
        }

    def assert_disjoint_from(self, *others: PopulationLedger) -> None:
        """Assert this population shares no global ID block or active id with ``others`` (e.g. the cortex)."""
        assert_disjoint_populations([self.ledger, *others])

    def assert_partitioned(self) -> None:
        """Assert the ledger is internally consistent and the fibers truly partition (no id counted twice)."""
        self.ledger.assert_invariants()
        c = self.census()
        if c["N_unique_sf_fibers"] != c["N_sf_fibers_summed"]:
            raise AssertionError(
                f"sf_arc fibers are double-counted: {c['N_sf_fibers_summed']} summed vs "
                f"{c['N_unique_sf_fibers']} unique"
            )


# ── geometry seeding (build-time descriptors; native counts/geometry arrive from CellState) ────────
def _unit(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Return ``v`` normalised, or ``v`` unchanged when it is degenerate (never raises, never renormalises 0)."""
    n = float(np.linalg.norm(v))
    return v / n if n > _EPS else np.asarray(v, np.float64)


def _lateral_axis(axis: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Deterministic unit vector perpendicular to ``axis`` — the sarcomere's lateral offset direction.

    WHICH transverse direction the two anti-parallel filaments are separated along is a geometric degree of
    freedom of the bundle, so it is fixed DETERMINISTICALLY rather than sampled: ``axis × ẑ``, i.e. the
    substrate-PARALLEL transverse direction, so a basal fiber's anti-parallel partner lies beside it IN the basal
    plane instead of being buried below the substrate.  When ``axis`` is itself ẑ-parallel that cross product
    degenerates and ``axis × x̂`` is used.

    HONEST SCOPE (measured 2026-07-25 — this is NOT a free gauge).  The choice does not change the STRAIGHT
    classes' load path, but it does move the connector ENDPOINTS by ``±lateral/2``: a ventral FA site shifts
    0.2 µm within the basal plane, and a cap LINC anchor leaves the nominal nuclear radius by ~7 nm.  For the
    CURVED cap it also changes the load path (``peak_frac`` 0.000 → 0.478) and can flip which arc apex a dorsal
    free end pairs with when two apices are equidistant.  That is one further reason the curved cap is not a
    motor station (:attr:`SFBundle.motor_station_ready`); the straddle geometry targets the straight classes, and
    the legacy default path is byte-identical either way.
    """
    axis = _unit(axis)
    lateral = np.cross(axis, np.array([0.0, 0.0, 1.0]))
    if float(np.linalg.norm(lateral)) <= 1.0e-8:            # axis ∥ ẑ → any perpendicular is equivalent
        lateral = np.cross(axis, np.array([1.0, 0.0, 0.0]))
    return _unit(lateral)


def _is_straight_sarcomere(
    a0: npt.NDArray[np.float64], mid: npt.NDArray[np.float64], b0: npt.NDArray[np.float64],
) -> bool:
    """True iff the two filament directions are COLLINEAR, i.e. the sarcomere is straight (not arched).

    A straight sarcomere has ``mid`` exactly on the ``a0``→``b0`` chord, so ``t_a = unit(mid−a0)`` equals
    ``t_b = unit(b0−mid)`` to floating-point precision and the two anti-parallel filaments are exactly
    anti-parallel over their shared span.  The perinuclear cap's ``mid`` is an arch APEX, so its two
    directions differ and the pair is not straddle-able by a rigid bipolar minifilament.

    The comparison is a floating-point exactness check on a geometric identity (collinearity), not a physical
    tolerance to tune: the straight classes are built from ``mid = ½(a0+b0)`` and satisfy it to ~1e-16.
    """
    t_a = _unit(np.asarray(mid, np.float64) - np.asarray(a0, np.float64))
    t_b = _unit(np.asarray(b0, np.float64) - np.asarray(mid, np.float64))
    return bool(np.allclose(t_a, t_b, rtol=0.0, atol=1.0e-12))


def _sarcomere(
    a0: npt.NDArray[np.float64], mid: npt.NDArray[np.float64], b0: npt.NDArray[np.float64],
    n_per_fiber: int,
    *,
    lateral_um: float | None = None,
    overlap_um: float | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64], npt.NDArray[np.int64],
           npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """One antiparallel two-filament sarcomere ``a0 -> mid -> b0`` (barbed ends OUTWARD, one bipolar NMII seed).

    Fiber A runs ``a0 -> mid`` with its barbed end at ``a0`` (polarity -1); fiber B runs ``mid -> b0`` with its
    barbed end at ``b0`` (polarity +1).  Both barbed ends point away from ``mid`` -> the interior NMII head pair
    pulls the two outer ends together (contraction).  Returns local pos/offsets/polarity/myo pairs.

    Two geometries are available:

    * **legacy end-to-end** (``lateral_um``/``overlap_um`` omitted): the two filaments MEET at ``mid`` — fiber
      A's last node and fiber B's first node carry the SAME coordinate.  A real bipolar minifilament cannot
      straddle that: with the two filaments co-located there is no anti-parallel pair to sit between, and
      :func:`aleph.components.motor.minifilament_topology.straddle_frame` falls back to an arbitrary perpendicular,
      throwing the heads off the actin (the SF instance of the cortex ``defect#3`` head-placement failure).
      Kept as the default so every landed load-path/connector gate stays byte-identical.
    * **straddle overlap** (both supplied): the two filaments run SIDE BY SIDE, separated by ``lateral_um``
      across the bundle axis and sharing an axial ``overlap_um`` around ``mid`` — each filament EXTENDS past
      ``mid`` along its OWN direction by ``overlap_um/2``, so a curved bundle (the perinuclear cap, whose
      ``mid`` is an arch apex, not the chord midpoint) keeps its arch instead of being straightened.  This is
      the geometry a real NMII minifilament needs: with ``lateral_um = 2·head_offset`` its ``±head_offset``
      heads land ON the two filament lines.  Both filaments keep their outward barbed ends, so the contractile
      SENSE is unchanged.  Note the lateral offset does move the anchor endpoints (see :func:`_lateral_axis` for
      the measured consequences, which include a changed load path for the curved cap) — only a STRAIGHT
      sarcomere is a motor station.

    Args:
        a0 / mid / b0: the sarcomere's outer end, interior station, and other outer end ``(3,)`` [µm].
        n_per_fiber: nodes per filament (>= 3).
        lateral_um: anti-parallel lateral separation [µm] (straddle mode; must be supplied with ``overlap_um``).
        overlap_um: axial overlap the two filaments share around ``mid`` [µm] (straddle mode).

    Returns:
        ``(pos, fiber_offsets, polarity, myo_i, myo_j)`` with local node indices; in straddle mode ``myo_i`` /
        ``myo_j`` are the A/B nodes AXIALLY CLOSEST to ``mid`` (the anti-parallel pair the motor straddles).

    Raises:
        ValueError: if only one of ``lateral_um``/``overlap_um`` is supplied, if either is non-positive, or if
            the requested overlap does not fit inside the shorter half of the sarcomere.
    """
    if (lateral_um is None) != (overlap_um is None):
        raise ValueError(
            "sarcomere straddle geometry needs BOTH lateral_um and overlap_um (they are two halves of one "
            "NMII-derived placement: lateral = 2·head_offset, overlap = backbone contour); supply both or "
            "neither"
        )
    if lateral_um is None:
        fiber_a = np.linspace(a0, mid, n_per_fiber)      # nodes 0..n-1, barbed at node 0 (=a0)
        fiber_b = np.linspace(mid, b0, n_per_fiber)      # nodes n..2n-1, barbed at node 2n-1 (=b0)
        pos = np.vstack([fiber_a, fiber_b])
        fiber_offsets = np.array([0, n_per_fiber, 2 * n_per_fiber], np.int64)
        polarity = np.array([-1, +1], np.int64)
        # NMII bipolar seed straddles the interior overlap: last node of A (n-1) <-> first node of B (n)
        myo_i = np.array([n_per_fiber - 1], np.int64)
        myo_j = np.array([n_per_fiber], np.int64)
        return pos, fiber_offsets, polarity, myo_i, myo_j

    if not (np.isfinite(lateral_um) and lateral_um > 0.0):
        raise ValueError(f"sarcomere lateral_um must be positive-finite, got {lateral_um!r}")
    if not (np.isfinite(overlap_um) and overlap_um > 0.0):
        raise ValueError(f"sarcomere overlap_um must be positive-finite, got {overlap_um!r}")

    t_a = _unit(np.asarray(mid, np.float64) - np.asarray(a0, np.float64))   # A's build direction (a0 -> mid)
    t_b = _unit(np.asarray(b0, np.float64) - np.asarray(mid, np.float64))   # B's build direction (mid -> b0)
    half_ov = 0.5 * float(overlap_um)
    reach_a = float(np.linalg.norm(np.asarray(mid, np.float64) - np.asarray(a0, np.float64)))
    reach_b = float(np.linalg.norm(np.asarray(b0, np.float64) - np.asarray(mid, np.float64)))
    if half_ov >= min(reach_a, reach_b):
        raise ValueError(
            f"sarcomere overlap_um={overlap_um} does not fit: half-overlap {half_ov:.4g} µm must be shorter "
            f"than each half-sarcomere ({reach_a:.4g} / {reach_b:.4g} µm)"
        )
    # The lateral direction is fixed on the OVERLAP axis (the mean of the two filament directions), so the two
    # filaments are separated across the segment they actually share.  For a straight sarcomere (ventral /
    # dorsal / arc, where mid is exactly the chord midpoint) t_a == t_b and this is exact.
    lateral = 0.5 * float(lateral_um) * _lateral_axis(_unit(t_a + t_b) if float(np.linalg.norm(t_a + t_b)) > _EPS
                                                      else t_a)
    fiber_a = np.linspace(a0 + lateral, mid + half_ov * t_a + lateral, n_per_fiber)
    fiber_b = np.linspace(mid - half_ov * t_b - lateral, b0 - lateral, n_per_fiber)
    pos = np.vstack([fiber_a, fiber_b])
    fiber_offsets = np.array([0, n_per_fiber, 2 * n_per_fiber], np.int64)
    polarity = np.array([-1, +1], np.int64)
    # the anti-parallel pair the minifilament straddles: the A/B nodes axially closest to the shared station.
    axis_ov = _unit(t_a + t_b) if float(np.linalg.norm(t_a + t_b)) > _EPS else t_a
    s_a = (fiber_a - np.asarray(mid, np.float64)) @ axis_ov
    s_b = (fiber_b - np.asarray(mid, np.float64)) @ axis_ov
    myo_i = np.array([int(np.argmin(np.abs(s_a)))], np.int64)
    myo_j = np.array([n_per_fiber + int(np.argmin(np.abs(s_b)))], np.int64)
    return pos, fiber_offsets, polarity, myo_i, myo_j


def _interior_nodes(n_per_fiber: int) -> npt.NDArray[np.int64]:
    """Local mid-shaft node indices of a two-fiber sarcomere (exclude the four outer/overlap ends)."""
    keep = []
    for f in range(2):
        base = f * n_per_fiber
        for k in range(1, n_per_fiber - 1):
            keep.append(base + k)
    return np.asarray(keep, np.int64)


def build_sf_arc_population(
    *,
    n_ventral: int = 4,
    n_dorsal: int = 2,
    n_arc: int = 2,
    n_cap: int = 2,
    n_per_fiber: int = 5,
    cell_radius_um: float = 7.5,
    nucleus_radius_um: float = 3.0,
    id_base: int = 1_000_000,
    f_head_pN: float | None = None,
    sarcomere_lateral_um: float | None = None,
    sarcomere_overlap_um: float | None = None,
) -> SFArcPopulation:
    """Build a structural ``sf_arc`` population: ventral/dorsal SF + transverse arcs + perinuclear cap.

    The counts/geometry are a small STRUCTURAL population for the Mac gates (native forbidden), NOT the sourced
    SF inventory — the physiological count/length distribution arrives from ``CellState`` in a biology phase and
    is drawn against this same fixed-capacity ledger (never resampled to a density here).  Every filament gets a
    unique id from the ``sf_arc`` block ``[id_base, id_base + 2*(n_ventral+n_dorsal+n_arc+n_cap))``.

    Args:
        n_ventral, n_dorsal, n_arc, n_cap: number of SF instances per taxonomy class (each owns two filaments).
        n_per_fiber: nodes per filament (>=3 so a mid-shaft cortex-transient node exists).
        cell_radius_um: basal-plane extent used to place the ventral/dorsal/arc geometry.
        nucleus_radius_um: nuclear radius used to arch the perinuclear cap over the nucleus.
        id_base: inclusive lower bound of this component's global filament-ID block.
        f_head_pN: provisional per-head NMII stall [pN] or ``None`` (GAP — the default; report-not-tune).
        sarcomere_lateral_um: anti-parallel lateral separation [µm] for the STRADDLE sarcomere geometry.  Supply
            it together with ``sarcomere_overlap_um`` when the population must HOST a bipolar NMII minifilament
            (the ``nmii_sf_motor`` edge): both are DERIVED from that minifilament's own geometry —
            ``lateral = 2·head_offset_um`` (so the ``±head_offset`` heads land ON the two filament lines) and
            ``overlap = backbone_length_um`` (so the whole backbone lies inside the shared span).  Omit both
            (the default) for the legacy end-to-end sarcomere, which is byte-identical to every landed gate but
            CANNOT be straddled (the two filaments are co-located at ``mid``).
        sarcomere_overlap_um: axial overlap the two anti-parallel filaments share around ``mid`` [µm].

    Returns:
        The built :class:`SFArcPopulation` (population + ledger + connector endpoint sets).
    """
    if n_per_fiber < 3:
        raise ValueError("n_per_fiber must be >= 3 so a bundle has a mid-shaft (cortex-transient) node")
    rng_free = float(cell_radius_um)
    z_dorsal = 0.6 * cell_radius_um                       # dorsal SF rise / arc elevation
    specs: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, str]] = []

    # ventral SF: basal plane (z=0), spanning x, FA at BOTH outer ends -> closed contractile dipole.
    for k in range(n_ventral):
        y = (-0.5 + (k + 0.5) / max(n_ventral, 1)) * rng_free
        a0 = np.array([-0.9 * rng_free, y, 0.0])
        b0 = np.array([+0.9 * rng_free, y, 0.0])
        specs.append((SFClass.VENTRAL, a0, 0.5 * (a0 + b0), b0, "fa_both"))

    # dorsal SF: basal FA end -> rises to a dorsal FREE end (which feeds a transverse-arc joint).
    dorsal_free_globals: list[int] = []
    for k in range(n_dorsal):
        x = (-0.4 + 0.8 * (k + 0.5) / max(n_dorsal, 1)) * rng_free
        a0 = np.array([x, -0.8 * rng_free, 0.0])         # basal FA end
        b0 = np.array([x, -0.2 * rng_free, z_dorsal])    # dorsal free end (top)
        specs.append((SFClass.DORSAL, a0, 0.5 * (a0 + b0), b0, "fa_low_free_high"))

    # transverse arc: elevated (z=z_dorsal), transverse (along x near the front), NO FA; couples dorsal + cortex.
    for k in range(n_arc):
        y = (-0.3 + 0.2 * k) * rng_free
        a0 = np.array([-0.7 * rng_free, y, z_dorsal])
        b0 = np.array([+0.7 * rng_free, y, z_dorsal])
        specs.append((SFClass.TRANSVERSE_ARC, a0, 0.5 * (a0 + b0), b0, "none"))

    # perinuclear cap: arches OVER the nucleus (dorsal), LINC-anchored at both ends near the envelope.
    for k in range(n_cap):
        ang = (k / max(n_cap, 1)) * np.pi
        ax, ay = np.cos(ang), np.sin(ang)
        a0 = np.array([-nucleus_radius_um * ax, -nucleus_radius_um * ay, nucleus_radius_um])
        b0 = np.array([+nucleus_radius_um * ax, +nucleus_radius_um * ay, nucleus_radius_um])
        mid = np.array([0.0, 0.0, nucleus_radius_um + 0.5 * nucleus_radius_um])   # arch apex above the nucleus
        specs.append((SFClass.PERINUCLEAR_CAP, a0, mid, b0, "linc_both"))

    n_total_fibers = 2 * len(specs)
    ledger = PopulationLedger(component=SF_COMPONENT, id_base=id_base, capacity=n_total_fibers)
    ledger.seed_active(range(id_base, id_base + n_total_fibers))

    bundles: list[SFBundle] = []
    flat_pos: list[np.ndarray] = []
    fa_sites: list[int] = []
    linc_sites: list[int] = []
    cortex_sites: list[int] = []
    dorsal_free_nodes: list[int] = []
    arc_apex_nodes: list[int] = []
    node_cursor = 0
    next_id = id_base
    n_nodes_total = 0
    n_heads_total = 0

    for (sf_class, a0, mid, b0, anchor_kind) in specs:
        lpos, loff, lpol, lmi, lmj = _sarcomere(
            a0, mid, b0, n_per_fiber,
            lateral_um=sarcomere_lateral_um, overlap_um=sarcomere_overlap_um)
        n_local = lpos.shape[0]
        gids = np.array([next_id, next_id + 1], np.int64)
        next_id += 2

        outer_a, outer_b = 0, n_local - 1                 # the two outer ends of the sarcomere
        mid_iface = n_per_fiber - 1                        # interior overlap (NMII station)
        interior = _interior_nodes(n_per_fiber)

        fa_local = np.zeros(0, np.int64)
        linc_local = np.zeros(0, np.int64)
        free_local = np.zeros(0, np.int64)
        if anchor_kind == "fa_both":
            fa_local = np.array([outer_a, outer_b], np.int64)
        elif anchor_kind == "fa_low_free_high":
            fa_local = np.array([outer_a], np.int64)       # basal FA end
            free_local = np.array([outer_b], np.int64)     # dorsal free end
        elif anchor_kind == "linc_both":
            linc_local = np.array([outer_a, outer_b], np.int64)
        # "none" (transverse arc): no FA/LINC; only cortex-transient + the internal arc joint

        # global-node bookkeeping
        fa_sites.extend(int(node_cursor + i) for i in fa_local)
        linc_sites.extend(int(node_cursor + i) for i in linc_local)
        cortex_sites.extend(int(node_cursor + i) for i in interior)
        if free_local.size:
            dorsal_free_nodes.append(int(node_cursor + free_local[0]))
        if sf_class == SFClass.TRANSVERSE_ARC:
            arc_apex_nodes.append(int(node_cursor + mid_iface))

        bundles.append(SFBundle(
            sf_class=sf_class, pos=lpos, fiber_offsets=loff, polarity=lpol, myo_i=lmi, myo_j=lmj,
            fa_local=fa_local, linc_local=linc_local, cortex_local=interior, free_local=free_local,
            global_fiber_ids=gids, node_base=node_cursor,
            sarcomere_lateral_um=sarcomere_lateral_um, sarcomere_overlap_um=sarcomere_overlap_um,
            motor_station_ready=(sarcomere_lateral_um is not None and _is_straight_sarcomere(a0, mid, b0))))
        flat_pos.append(lpos)
        node_cursor += n_local
        n_nodes_total += n_local
        n_heads_total += int(lmi.size)

    # internal dorsal-free <-> transverse-arc joints (dorsal_arc_crosslink): pair each dorsal free end to the
    # nearest arc apex node (a real geometric coupling, not a co-location weld).
    pos_flat = np.vstack(flat_pos) if flat_pos else np.zeros((0, 3), np.float64)
    joints: list[tuple[int, int]] = []
    if dorsal_free_nodes and arc_apex_nodes:
        arc_xyz = pos_flat[np.asarray(arc_apex_nodes, np.int64)]
        for fn in dorsal_free_nodes:
            d = np.linalg.norm(arc_xyz - pos_flat[fn], axis=1)
            joints.append((fn, arc_apex_nodes[int(np.argmin(d))]))
    dorsal_arc_joints = np.asarray(joints, np.int64).reshape(-1, 2)

    ledger.n_nodes = n_nodes_total
    ledger.n_heads = n_heads_total
    ledger.assert_invariants()

    status = OK_STATUS if f_head_pN is not None else GAP_STATUS
    return SFArcPopulation(
        bundles=bundles, ledger=ledger, pos=pos_flat,
        fa_sites=np.asarray(fa_sites, np.int64), linc_sites=np.asarray(linc_sites, np.int64),
        cortex_sites=np.asarray(cortex_sites, np.int64), dorsal_arc_joints=dorsal_arc_joints,
        f_head_pN=f_head_pN, magnitude_status=status)


def validate_sf_connectors(architecture: CellArchitecture, population: SFArcPopulation) -> None:
    """Assert the population feeds exactly the canonical ``sf_arc`` edges of ``architecture`` — nothing else.

    Every connector this population wires must be a registered, bidirectional + adjoint edge whose endpoints
    include ``sf_arc`` (co-location is never a connection).  The three inter-component edges join ``sf_arc`` to
    ``focal_adhesion`` / ``cortex`` / ``nucleus``; ``dorsal_arc_crosslink`` is the internal dorsal<->arc joint.
    This validates the wiring against the SoT graph without declaring any molecular joint here.

    Raises:
        ValueError: if a fed connector is missing, mis-scoped, one-way, or not incident to ``sf_arc``.
    """
    by_name = {c.name: c for c in architecture.connectors}
    expected = {
        FA_ACTIN_ANCHOR: ("focal_adhesion", ConnectorScope.INTER_COMPONENT),
        SF_CORTEX_TRANSIENT: ("cortex", ConnectorScope.INTER_COMPONENT),
        ACTIN_CAP_LINC: ("nucleus", ConnectorScope.INTER_COMPONENT),
        DORSAL_ARC_CROSSLINK: (SF_COMPONENT, ConnectorScope.INTERNAL),
    }
    fed = set(population.connector_endpoints())
    if fed != set(expected):
        raise ValueError(f"population feeds {sorted(fed)}, expected {sorted(expected)}")
    for name, (partner, scope) in expected.items():
        connector = by_name.get(name)
        if connector is None:
            raise ValueError(f"architecture must register {name!r}")
        if connector.scope is not scope:
            raise ValueError(f"{name} must be {scope.value} scope")
        endpoints = {connector.component_a, connector.component_b}
        if SF_COMPONENT not in endpoints or partner not in endpoints:
            raise ValueError(f"{name} must join sf_arc and {partner!r}, got {sorted(endpoints)}")
        if not connector.bidirectional or not connector.adjoint_transfer_required:
            raise ValueError(f"{name} must be bidirectional with adjoint transfer (Newton 3rd)")
        if connector.kinetics and not connector.commit_on_accept:
            raise ValueError(f"{name} is kinetic but does not commit only on an accepted step")
