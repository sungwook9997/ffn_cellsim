"""Per-component population / ID / free-list ledger with the disjoint-ID no-double-count invariant.

Spec §3 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``) + PI 2026-07-22 (SF-separate): cortex,
SF/arc, lamellipodium, filopodium, … are separate state-owning components, **each owning a DISJOINT
filament population**.  Every physical filament belongs to exactly ONE component; components couple only
through explicit connectors, never a shared array.  :func:`~aleph.engine.dump_state` already assigns
a globally-unique filament-ID namespace and asserts disjoint blocks at *viz* time; this module promotes
that assert to a **build-time + accepted-step** invariant and adds the standing inventory + the free-list
discipline the ratified doc requires ("track unique active IDs, dormant allocated capacity, nodes,
explicit states, peak bytes").

The one hard mechanism: **count changes only via the free-list, never by resampling to a target.**
Nucleation pops a dormant slot (:meth:`PopulationLedger.allocate`); severing splits one filament into two
by popping one more slot (:meth:`PopulationLedger.sever`); depolymerisation/dissolution returns a slot
(:meth:`PopulationLedger.release`).  There is no "draw N filaments to hit a density" path — the persistent
population is generated ONCE (from :class:`~aleph.engine.cell_state.CellState`) and only ever
allocated/released against its fixed capacity.

This is a host-side accounting structure over integer ID ranges (the native lane mirrors it with device
bitmaps / free-lists); it allocates no device memory and is CPU-importable, so the structural gates drive
it directly.  Contour-length / mass conservation on a sever is a biology-phase concern layered on top; here
the invariant is purely the ID bookkeeping.

Sanity Gate:
    * ownership: ``active_ids`` and ``free_list`` partition the fixed block ``[id_base, id_base+capacity)``;
      their union is the whole block and their intersection is empty (checked by :meth:`assert_invariants`).
    * conservation: ``active_count + dormant_count == capacity`` every step — allocation is conserved, never
      resampled; :meth:`allocate` raises rather than growing capacity to hit a target.
    * boundary/sign: :func:`assert_disjoint_populations` rejects any overlap of blocks OR of active IDs across
      components — a filament counted under two components is a rejected build, not a silent double-count.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

__all__ = ["PopulationLedger", "assert_disjoint_populations"]


@dataclass(slots=True)
class PopulationLedger:
    """One component's unique-ID population, dormant free-list, and standing inventory census.

    Args:
        component: the owning component name (e.g. ``"cortex"``); must be non-empty.
        id_base: inclusive lower bound of this component's global filament-ID block.
        capacity: number of allocated slots; the block is ``[id_base, id_base + capacity)``.
        n_nodes: standing count of explicit nodes (host metadata, informational).
        n_heads: standing count of explicit motor heads.
        n_explicit_states: standing count of explicit kinetic states (bonds/clutches/…).
        n_field_cells: standing count of owned field cells.
        peak_gpu_bytes: exact peak device bytes for this component (informational; the ratified doc
            requires it be tracked, never used to justify lowering biological density).
    """

    component: str
    id_base: int
    capacity: int
    n_nodes: int = 0
    n_heads: int = 0
    n_explicit_states: int = 0
    n_field_cells: int = 0
    peak_gpu_bytes: int = 0
    _active: set[int] = field(default_factory=set, init=False, repr=False)
    _free: list[int] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        if not self.component.strip():
            raise ValueError("population ledger needs a component name")
        for label, value in (("id_base", self.id_base), ("capacity", self.capacity)):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"population ledger {label} must be an int")
        if self.id_base < 0:
            raise ValueError("id_base must be nonnegative")
        if self.capacity < 0:
            raise ValueError("capacity must be nonnegative")
        for label in ("n_nodes", "n_heads", "n_explicit_states", "n_field_cells", "peak_gpu_bytes"):
            value = getattr(self, label)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"population ledger {label} must be a nonnegative int")
        # The whole block starts dormant; the initial population is seeded explicitly via `seed_active`.
        self._free = list(range(self.id_base + self.capacity - 1, self.id_base - 1, -1))

    @property
    def block(self) -> tuple[int, int]:
        """Half-open global filament-ID block ``[id_base, id_base + capacity)`` owned by this component."""
        return (self.id_base, self.id_base + self.capacity)

    @property
    def active_ids(self) -> frozenset[int]:
        """Snapshot of the currently active (occupied) global filament IDs."""
        return frozenset(self._active)

    @property
    def active_count(self) -> int:
        """Number of active filaments."""
        return len(self._active)

    @property
    def dormant_count(self) -> int:
        """Number of dormant (allocated-but-free) slots."""
        return len(self._free)

    def owns(self, fid: int) -> bool:
        """Return whether ``fid`` falls in this component's global block."""
        lo, hi = self.block
        return lo <= fid < hi

    def seed_active(self, ids: Iterable[int]) -> None:
        """Mark an initial set of IDs active by drawing them from the dormant free-list (build-time only).

        This is the ONE-time generation of the persistent population (from ``CellState`` inventory).  Each
        seeded id must be a dormant slot inside this block; there is no path to activate an out-of-block id.
        """
        free_set = set(self._free)
        for fid in ids:
            if not self.owns(fid):
                raise ValueError(f"{self.component}: id {fid} is outside its block {self.block}")
            if fid not in free_set:
                raise ValueError(f"{self.component}: id {fid} is already active or unknown")
            free_set.remove(fid)
            self._active.add(fid)
        self._free = sorted(free_set, reverse=True)

    def allocate(self) -> int:
        """Nucleate one filament: pop a dormant slot into the active set.

        Raises ``RuntimeError`` when the free-list is empty — capacity is fixed and is NEVER grown by
        resampling to hit a target density (the ratified no-resample rule).
        """
        if not self._free:
            raise RuntimeError(
                f"{self.component}: population free-list exhausted; capacity is fixed and must not be "
                "resampled to a target (increase sourced capacity via CellState, not at runtime)"
            )
        fid = self._free.pop()
        self._active.add(fid)
        return fid

    def release(self, fid: int) -> None:
        """Depolymerise/dissolve one filament: return its slot to the dormant free-list."""
        if fid not in self._active:
            raise ValueError(f"{self.component}: id {fid} is not active and cannot be released")
        self._active.discard(fid)
        self._free.append(fid)

    def sever(self, fid: int) -> int:
        """Sever one filament into two: keep ``fid`` active and pop one more slot for the new fragment.

        Returns the new fragment's global id.  Contour-length/mass conservation across the split is a
        biology-phase concern layered on top of this ID bookkeeping.
        """
        if fid not in self._active:
            raise ValueError(f"{self.component}: id {fid} is not active and cannot be severed")
        return self.allocate()

    def assert_invariants(self) -> None:
        """Assert the active/dormant partition and the fixed-capacity conservation law."""
        if len(self._active) + len(self._free) != self.capacity:
            raise AssertionError(
                f"{self.component}: active+dormant ({len(self._active)}+{len(self._free)}) "
                f"!= capacity ({self.capacity}) — allocation was resampled, not conserved"
            )
        if self._active & set(self._free):
            raise AssertionError(f"{self.component}: an id is both active and dormant")
        lo, hi = self.block
        for fid in self._active:
            if not lo <= fid < hi:
                raise AssertionError(f"{self.component}: active id {fid} escaped its block {self.block}")


def assert_disjoint_populations(ledgers: Sequence[PopulationLedger]) -> None:
    """Assert every physical filament belongs to exactly one component (build-time + accepted-step).

    Rejects any overlap of global ID *blocks* (structural) OR of *active* IDs (runtime) across the given
    component ledgers — the PI 2026-07-22 no-double-count invariant, promoted here from a viz-time check.

    Raises:
        ValueError: if two components share a block range or a live active filament id, or if two ledgers
            declare the same component name.
    """
    names = [ledger.component for ledger in ledgers]
    if len(set(names)) != len(names):
        raise ValueError("population ledgers must each own a distinct component")
    ordered = sorted(ledgers, key=lambda ledger: ledger.block)
    for earlier, later in zip(ordered, ordered[1:], strict=False):
        if earlier.block[1] > later.block[0]:
            raise ValueError(
                f"components {earlier.component!r} and {later.component!r} share overlapping global "
                f"ID blocks {earlier.block} and {later.block} — a filament would be double-counted"
            )
    seen: dict[int, str] = {}
    for ledger in ledgers:
        for fid in ledger.active_ids:
            owner = seen.get(fid)
            if owner is not None:
                raise ValueError(
                    f"active filament id {fid} is claimed by both {owner!r} and "
                    f"{ledger.component!r} — components must own disjoint populations"
                )
            seen[fid] = ledger.component
