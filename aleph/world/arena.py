"""The world arena: one fixed-capacity allocation, claimed as contiguous ID ranges.

WHAT THIS REPLACES, AND WHY.  Today six components each allocate their own ``position_d``/``force_d``
pair — ``surface_body``, ``cortex_state``, ``protrusion``, ``microtubule_rig``,
``intermediate_filament_rig`` and ``composed_native``'s SF owner — and six near-identical validators
re-derive the same shape.  The composed cell then addresses membrane, cortex and nucleus as *slice
views* of one array anyway (``STATE.md``: "Slice views, NOT private allocations"), so the private
allocation was already a fiction at the top of the stack while remaining real at the bottom.

This module makes the arrangement the declaration.  One allocation of fixed CAPACITY; a population is
a NAMED, CONTIGUOUS, HALF-OPEN ID RANGE claimed out of it.  Ownership becomes bookkeeping rather than
allocation, which is what ``PopulationLedger`` (``aleph/engine/population.py:19-22``) already says it
is: "a host-side accounting structure over integer ID ranges … it allocates no device memory".

THE THREE INVARIANTS, and each exists because something specific goes wrong without it.

1. **Capacity is fixed; LIVE is a contiguous prefix; every launch is over ``n_live``, never capacity.**
   A ten-million-node arena costs ~1.5 GiB, which is affordable — but a kernel launched at
   ``dim=capacity`` when a tenth of it is live wastes nine tenths of the memory traffic, and the
   dominant kernels in this engine are bandwidth-bound.  The prefix discipline is what makes the
   capacity free: launch ``dim=arena.n_live`` against the FULL array and no slicing is needed, because
   the live entries are exactly ``[0, n_live)``.

2. **Claims come from the tail and are never returned.**  There is deliberately no free list yet.  The
   build order is to weave populations IN, one at a time; nothing is removed, so nothing needs
   reclaiming, and a free list would be machinery for a case that does not exist.  When a builder first
   needs to remove a population, that is the moment to add one — and to decide then whether it
   compacts or fragments, which is a real question this module should not pre-answer.

3. **Every claim exposes its BYTE SPAN, and disjointness is asserted on spans, not pointers.**  This
   is not defensive; it is a bug the arena makes reachable.  ``coupled_solve.collect_force_arrays``
   deduplicates force arrays by device pointer, so a range whose ``id_base`` is 0 has the same pointer
   as the whole-world array and the two collapse into one entry — silently, in the ledger the balance
   gate depends on.  ``cortex_state.assert_component_state_disjoint`` already tests byte spans for
   exactly this reason and says so: "Pointer identity is not enough once any owner holds a SLICE VIEW".

WHAT THIS MODULE DELIBERATELY DOES NOT DO.  It computes no force, owns no physics, launches no kernel
and decides no acceptance.  It also allocates no solver workspace: whether the inner solve needs Krylov
vectors, per-node clock state, or velocities depends on a measurement that has not been taken yet, and
guessing would bake an answer into the layout.  Those arrays get claimed the same way everything else
does, once the answer exists.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``position``/``position_snap`` [µm], ``force`` [pN], ``mobility`` [µm/(pN·s)];
    ``strand_id`` and ``range_id`` are indices and carry no unit.
  * boundary — a claim of zero elements is rejected rather than returning an empty range that later
    reads as "present but silent"; a claim exceeding remaining capacity raises with both numbers.
  * conservation/invariant — claims partition ``[0, n_live)`` exactly: the ranges are contiguous, in
    claim order, with no gap and no overlap, and :meth:`WorldArena.assert_partitioned` checks it.
  * CFL/precision — no integration happens here; float64 throughout, matching the engine's arrays.
  * sign sense — not applicable; nothing is accumulated.
  * measurement protocol — ``census()`` is a host-side read of host-side bookkeeping. It touches no
    device array and is safe to call between accepted steps at any cadence, which is what the live
    floorplan view consumes.

engine units: length µm, force pN, drag pN·s/µm.  Runtime: NVIDIA Warp on CUDA for allocation; this
module is CPU-importable and allocates nothing until a device is given.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
import warp as wp

__all__ = ["Kind", "Claim", "WorldArena"]


class Kind(StrEnum):
    """The primitive kinds an arena allocates capacity for.

    Eight, and the set is closed by argument rather than by convenience: every cell structure this
    engine builds is a composition of these and nothing else.  Three candidates were considered and
    rejected — ``HEAD`` is a NODE plus a BOND whose far address is a point on a SEGMENT; ``SITE`` is an
    address resolved from live geometry, not an allocation; and a species pool is a named channel on a
    ``GRID_CELL`` field plus one conserved scalar per consuming range.

    ⚠ ``GRID_CELL`` was ADDED 2026-08-20, PI-approved, and the reason is written in the paragraph above:
    the rejection of a species pool cites ``GRID_CELL`` BY NAME, so the argument that closed a
    seven-member set already presupposed an eighth that had never been added.  Three sessions reached
    that sentence independently from different directions.  The cytosol is what forced it: the port
    source is Eulerian (``components/fluid/field_grid.py:9-11`` separates the PDE grid from the filament
    ``HashGrid`` as *"distinct allocations"*), ``FluidVolumeStateOwner`` validates ``ndim=3``, and the six
    ``*_cytosol_transfer`` edges couple through a Peskin-4 stencil around a POSITION — there is no
    cytosol node to bond to, and inventing one would replace the interpolation with a lumped surrogate.
    See ``docs/v2_audit/CYTOSOL_ARENA_REPRESENTATION_2026-08-20.md``.

    ⚠ Adding it costs one member and no allocation branch: this arena allocates arrays only for
    ``NODE`` (allocated in :meth:`WorldArena.__post_init__`).  Every other kind is capacity, claim and partition
    bookkeeping, and its arrays are owned by the builder — ``cortex.py``'s ``seg_node``/``seg_rest_um``,
    ``membrane.py``'s ``face_idx``/``hinge_idx``, and now ``cytosol.py``'s field arrays.  A ``GRID_CELL``
    claim counts ``nx·ny·nz`` cells against a 3-D array exactly as a ``SEGMENT`` claim counts ``S``
    against an ``(S, 2)`` array; the 1-D claim / n-D array shape difference is not new here.

    ``ANGLE3`` and ``ANGLE4`` stay separate because a dihedral is not a bending triple: they differ in
    arity, in stencil, and in which law reads them.  Merging them would promote scalar coefficients to
    arrays and add an in-kernel branch for no measured gain.
    """

    NODE = "node"
    STRAND = "strand"
    SEGMENT = "segment"
    ANGLE3 = "angle3"
    ANGLE4 = "angle4"
    FACE = "face"
    BOND = "bond"
    GRID_CELL = "grid_cell"


#: Bytes per element of the arena's per-NODE state, by array.  Used by :meth:`WorldArena.census` to
#: report a footprint that can be compared against a driver reading rather than trusted on its own —
#: the two disagreeing is allocator overhead, and is worth knowing rather than hiding.
_NODE_ARRAY_BYTES: dict[str, int] = {
    "position": 24, "force": 24, "position_snap": 24,
    "mobility": 8, "strand_id": 4, "range_id": 4,
}


@dataclass(frozen=True, slots=True)
class Claim:
    """One population's half-open range of one primitive kind.

    Attributes:
        population: the name the range belongs to, e.g. ``"cortex"``.  A label, never a type.
        kind: which primitive this range indexes.
        lo: inclusive lower bound.
        hi: exclusive upper bound.
    """

    population: str
    kind: Kind
    lo: int
    hi: int

    @property
    def count(self) -> int:
        """Number of elements in the range."""
        return self.hi - self.lo

    def owns(self, index: int) -> bool:
        """Whether ``index`` falls inside this range."""
        return self.lo <= index < self.hi

    def overlaps(self, other: "Claim") -> bool:
        """Whether two claims of the SAME kind share any index."""
        return self.kind == other.kind and self.lo < other.hi and other.lo < self.hi


@dataclass(slots=True)
class WorldArena:
    """One fixed-capacity world, claimed as contiguous ranges.

    Args:
        capacity: elements to reserve per :class:`Kind`.  Missing kinds reserve nothing.  The figure is
            a reservation, not a prediction — raise it when a compartment needs more, trim it when the
            inventory settles, exactly as a virtual address space is sized.
        device: the resolved CUDA device string, or ``None`` to build the bookkeeping alone.  ``None``
            is what makes this module CPU-importable and testable without a card; it allocates nothing
            and every range check still runs.

    Raises:
        RuntimeError: if ``device`` names a non-CUDA device.  There is no CPU simulation path and there
            must never be one, so resolving to ``cpu`` is an error rather than a fallback.
    """

    capacity: dict[Kind, int]
    device: str | None = None
    node_arrays: dict[str, wp.array] = field(default_factory=dict, repr=False)
    _claims: list[Claim] = field(default_factory=list, repr=False)
    _live: dict[Kind, int] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.capacity = {Kind(k): int(v) for k, v in self.capacity.items()}
        for kind, n in self.capacity.items():
            if n < 0:
                raise ValueError(f"capacity[{kind}] must be nonnegative; got {n}")
        self._live = {kind: 0 for kind in self.capacity}
        if self.device is None:
            return
        dev = wp.get_device(self.device)
        if not dev.is_cuda:
            raise RuntimeError(
                f"arena device resolved to {dev!r}, which is not CUDA. Warp CUDA is the only simulation "
                "runtime; a non-CUDA device raises rather than proceeding."
            )
        n = self.capacity.get(Kind.NODE, 0)
        if n:
            with wp.ScopedDevice(dev):
                self.node_arrays = {
                    "position": wp.zeros(n, dtype=wp.vec3d),
                    "force": wp.zeros(n, dtype=wp.vec3d),
                    "position_snap": wp.zeros(n, dtype=wp.vec3d),
                    "mobility": wp.zeros(n, dtype=wp.float64),
                    "strand_id": wp.full(n, -1, dtype=wp.int32),
                    "range_id": wp.full(n, -1, dtype=wp.int32),
                }

    # ── claiming ────────────────────────────────────────────────────────────────────────────────
    def claim(self, population: str, kind: Kind, count: int) -> Claim:
        """Claim ``count`` contiguous elements of ``kind`` from the tail, for ``population``.

        The claim extends the live prefix, so ``[0, n_live(kind))`` remains exactly the union of every
        claim of that kind, in claim order, with no gap.

        Raises:
            ValueError: on a non-positive count, or when the claim exceeds remaining capacity.  A
                zero-element claim is refused rather than returned empty: a population that is
                "present but has nothing" is indistinguishable downstream from one that failed to
                build, and the two are different facts.
        """
        kind = Kind(kind)
        if not population:
            raise ValueError("a claim needs a population name; an unnamed range cannot be attributed")
        if count <= 0:
            raise ValueError(
                f"claim({population}, {kind}, {count}) — a claim must be positive. An empty range reads "
                "downstream as 'present but silent', which is not distinguishable from a failed build."
            )
        cap = self.capacity.get(kind, 0)
        lo = self._live.get(kind, 0)
        hi = lo + count
        if hi > cap:
            raise ValueError(
                f"claim({population}, {kind}, {count}) exceeds capacity: {lo} live + {count} requested "
                f"> {cap} reserved. Raise the arena's capacity rather than lowering a population — the "
                "count comes from physiological density x geometry."
            )
        claim = Claim(population=population, kind=kind, lo=lo, hi=hi)
        self._claims.append(claim)
        self._live[kind] = hi
        return claim

    # ── reading ─────────────────────────────────────────────────────────────────────────────────
    def n_live(self, kind: Kind) -> int:
        """Live element count for ``kind`` — **the dimension every launch uses**, never the capacity."""
        return self._live.get(Kind(kind), 0)

    def claims(self, population: str | None = None, kind: Kind | None = None) -> tuple[Claim, ...]:
        """Claims, optionally filtered by population and/or kind, in claim order."""
        out = self._claims
        if population is not None:
            out = [c for c in out if c.population == population]
        if kind is not None:
            k = Kind(kind)
            out = [c for c in out if c.kind == k]
        return tuple(out)

    def population_of(self, kind: Kind, index: int) -> str | None:
        """The population owning ``index``, or ``None`` if it is unclaimed.

        This is the derivation that lets a bond record only ``(node_i, node_j)`` and have its component
        pair follow: the pair is read off the ranges rather than declared on the edge.
        """
        k = Kind(kind)
        for c in self._claims:
            if c.kind == k and c.owns(index):
                return c.population
        return None

    def byte_span(self, name: str, claim: Claim) -> tuple[int, int]:
        """Half-open ``(start, end)`` device byte span of ``claim`` within node array ``name``.

        Spans, not pointers, are what disjointness must be tested on once ranges share one allocation:
        a range at ``lo == 0`` has the same device pointer as the whole array, so pointer identity
        cannot separate them.

        Raises:
            KeyError: if ``name`` is not an allocated node array.
            ValueError: if ``claim`` is not a NODE claim, since only node arrays are spanned here.
        """
        if claim.kind is not Kind.NODE:
            raise ValueError(f"byte_span is for NODE claims; got {claim.kind}")
        arr = self.node_arrays[name]
        stride = _NODE_ARRAY_BYTES[name]
        base = int(arr.ptr) if arr.ptr else 0
        return base + claim.lo * stride, base + claim.hi * stride

    # ── device transfer ─────────────────────────────────────────────────────────────────────────
    def upload_nodes(self, claim: Claim, name: str, values) -> None:
        """Write ``values`` into node array ``name`` over ``claim``'s range.

        This is a HOST->DEVICE transfer and it belongs at build time or between accepted steps, never
        inside the inner loop. The engine's residency gate exists to catch exactly that, and it is
        non-vacuous only because a positive control puts one copy where it does not belong — so this
        method is deliberately not called by anything that runs per iteration.

        The write is a ``wp.copy`` into the claimed offset rather than an assignment to a slice, so the
        untouched remainder of the array keeps its contents and a short ``values`` cannot silently
        shorten the range.

        Args:
            claim: a NODE claim of this arena.
            name: an allocated node array.
            values: host data of exactly ``claim.count`` elements, in the array's dtype.

        Raises:
            RuntimeError: if the arena holds no device allocation.
            KeyError: if ``name`` is not an allocated node array.
            ValueError: if ``claim`` is not a NODE claim, exceeds the live prefix, or ``values`` has a
                length other than ``claim.count``.
        """
        arr = self._node_array(claim, name)
        host = wp.array(np.asarray(values), dtype=arr.dtype, device="cpu")
        if int(host.shape[0]) != claim.count:
            raise ValueError(
                f"upload_nodes({claim.population}, {name}): got {int(host.shape[0])} values for a range "
                f"of {claim.count}. A short write would leave the tail of the range at its previous "
                "contents while every launch still covers it."
            )
        wp.copy(arr, host, dest_offset=claim.lo, src_offset=0, count=claim.count)

    def download_nodes(self, claim: Claim, name: str):
        """Read ``claim``'s range of node array ``name`` back to the host, as a numpy array.

        For verification and rendering between accepted steps. Same rule as :meth:`upload_nodes`: never
        inside the inner loop.
        """
        arr = self._node_array(claim, name)
        return arr.numpy()[claim.lo:claim.hi]

    def _node_array(self, claim: Claim, name: str):
        """Resolve and validate a node array for a claim, or raise with the reason."""
        if not self.node_arrays:
            raise RuntimeError(
                "this arena holds no device allocation — it was built with device=None, which is the "
                "bookkeeping-only mode. Build it with a CUDA device to transfer."
            )
        if claim.kind is not Kind.NODE:
            raise ValueError(f"node transfer needs a NODE claim; got {claim.kind}")
        if claim.hi > self.n_live(Kind.NODE):
            raise ValueError(
                f"claim [{claim.lo}, {claim.hi}) reaches past the live prefix "
                f"[0, {self.n_live(Kind.NODE)}) — it is not a claim of this arena."
            )
        return self.node_arrays[name]

    # ── invariants ──────────────────────────────────────────────────────────────────────────────
    def assert_partitioned(self) -> None:
        """Assert the claims of every kind tile ``[0, n_live)`` exactly — no gap, no overlap.

        Contiguous tail-claiming makes this true by construction, which is the point: the assertion is
        cheap and it fails loudly the moment someone adds a code path that claims out of order.

        Raises:
            AssertionError: on a gap, an overlap, or a live count that disagrees with the claims.
        """
        for kind in self.capacity:
            spans = [c for c in self._claims if c.kind == kind]
            cursor = 0
            for c in spans:
                if c.lo != cursor:
                    raise AssertionError(
                        f"{kind} claims are not contiguous: {c.population} starts at {c.lo}, expected "
                        f"{cursor}. The live prefix is what makes launching over n_live correct."
                    )
                cursor = c.hi
            if cursor != self.n_live(kind):
                raise AssertionError(
                    f"{kind} live count {self.n_live(kind)} disagrees with the claims, which end at {cursor}"
                )

    def assert_disjoint_spans(self, name: str) -> None:
        """Assert every pair of NODE claims occupies disjoint byte spans of array ``name``.

        Redundant with :meth:`assert_partitioned` while claims are contiguous — deliberately so. It is
        the check that survives if compaction or a free list ever makes ranges non-contiguous, and it
        is the one that speaks the same language as the ledger's own guard.

        Raises:
            AssertionError: if two claims overlap in bytes.
        """
        node_claims = [c for c in self._claims if c.kind is Kind.NODE]
        spans = [(self.byte_span(name, c), c) for c in node_claims]
        spans.sort(key=lambda item: item[0][0])
        for (a_span, a), (b_span, b) in zip(spans, spans[1:]):
            if a_span[1] > b_span[0]:
                raise AssertionError(
                    f"{a.population} and {b.population} overlap in {name}: "
                    f"{a_span} vs {b_span}. Byte spans, not pointers — a range at index 0 shares the "
                    "whole array's pointer, which is how a ledger channel silently collapses."
                )

    # ── census ──────────────────────────────────────────────────────────────────────────────────
    def census(self) -> dict[str, object]:
        """Host-side census: what is claimed, by whom, and what it costs.

        This is the payload the live floorplan view reads — kilobytes, host-side, touching no device
        array, so it is safe to emit every accepted step.  ``node_bytes_arithmetic`` is a computed
        figure and is labelled one: it is meant to be compared against a driver reading, not trusted
        instead of it.
        """
        populations: dict[str, dict[str, int]] = {}
        for c in self._claims:
            populations.setdefault(c.population, {})[str(c.kind)] = (
                populations.setdefault(c.population, {}).get(str(c.kind), 0) + c.count
            )
        n_node_cap = self.capacity.get(Kind.NODE, 0)
        return {
            "capacity": {str(k): v for k, v in self.capacity.items()},
            "live": {str(k): self.n_live(k) for k in self.capacity},
            "utilisation": {
                str(k): (self.n_live(k) / v if v else 0.0) for k, v in self.capacity.items()
            },
            "populations": populations,
            "n_claims": len(self._claims),
            "device": self.device,
            "node_arrays": sorted(self.node_arrays),
            "node_bytes_arithmetic": n_node_cap * sum(_NODE_ARRAY_BYTES.values()),
        }


def _demo() -> None:
    """Self-check: the invariants this module exists to hold, asserted without a device."""
    arena = WorldArena(capacity={Kind.NODE: 1_000, Kind.SEGMENT: 900, Kind.BOND: 500})

    cortex = arena.claim("cortex", Kind.NODE, 600)
    membrane = arena.claim("membrane", Kind.NODE, 300)
    assert (cortex.lo, cortex.hi) == (0, 600)
    assert (membrane.lo, membrane.hi) == (600, 900)
    assert arena.n_live(Kind.NODE) == 900, "live is the claim total, not the capacity"
    assert arena.n_live(Kind.SEGMENT) == 0, "an unclaimed kind is live-zero, not live-capacity"

    # The derivation that lets a bond carry only its two node ids.
    assert arena.population_of(Kind.NODE, 0) == "cortex"
    assert arena.population_of(Kind.NODE, 599) == "cortex"
    assert arena.population_of(Kind.NODE, 600) == "membrane"
    assert arena.population_of(Kind.NODE, 950) is None, "beyond live is unclaimed, not the last owner"

    arena.assert_partitioned()
    assert not cortex.overlaps(membrane)

    # A claim past the reservation must name both numbers rather than silently truncating.
    try:
        arena.claim("nucleus", Kind.NODE, 200)
    except ValueError as exc:
        assert "exceeds capacity" in str(exc) and "900" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("over-capacity claim must raise")

    # An empty claim is refused: "present but silent" and "failed to build" are different facts.
    for bad in (0, -1):
        try:
            arena.claim("ecm", Kind.BOND, bad)
        except ValueError as exc:
            assert "must be positive" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"claim of {bad} must raise")

    # Out-of-order claiming is what assert_partitioned exists to catch.
    broken = WorldArena(capacity={Kind.NODE: 100})
    broken._claims.append(Claim("ghost", Kind.NODE, 10, 20))
    broken._live[Kind.NODE] = 20
    try:
        broken.assert_partitioned()
    except AssertionError as exc:
        assert "not contiguous" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a gap at the head must fail the partition check")

    c = arena.census()
    assert c["live"]["node"] == 900 and c["capacity"]["node"] == 1_000
    assert c["populations"]["cortex"]["node"] == 600
    assert c["node_bytes_arithmetic"] == 1_000 * 88
    print("arena self-check OK —", {k: c[k] for k in ("live", "n_claims", "node_bytes_arithmetic")})


if __name__ == "__main__":
    _demo()
