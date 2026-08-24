"""CellState — the slow-context state vector that conditions rates/inventory but owns NO physical state.

Spec §5 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``).  CellState is a small state vector, fixed
for a seconds-to-minutes run, over six axes: ``lineage · emt · cycle · adhesion · geometry ·
osmotic_state``.  The first ratified sourced instance is
``MCF7 × hybrid_E_M × G1 × collagen_spread × polarized × physiological_resting``
(:func:`mcf7_reference_state`).

CellState **sets ONLY** initial inventory distributions and reaction *rates*, and only conditionally
(one CellState changes many components at once).  It **MUST NOT** directly set cortex tension, total
contractile force, a bound-head count, a filament-count target, or crosslink topology — those *emerge*
from mechanics + events.  That boundary is enforced structurally: every value CellState emits is a
:class:`Rate` or a :class:`Distribution`, each of which is rejected if its target is a
force/pos/tension/count kind, and :func:`assert_cellstate_writes_no_force` proves the instance exposes no
device array and no force/count-writing setter.  ``CellState → force`` is a rejected build, mirroring the
engine's "co-location ≠ connection" enforcement.

Sourced-only discipline (thin → thick): every :class:`Rate`/:class:`Distribution` carries a non-empty KB
provenance (Magic-Number Block — a value with no source is rejected).  Add an axis value or a conditional
branch ONLY when it has a sourced value; :meth:`CellState.rate` / :meth:`CellState.initial_inventory`
raise for an unsourced ``(component, reaction)`` rather than returning a convenient default.  The first
ratified instance ships with EMPTY registries: the axes list is the Step-1 contract; the sourced rate and
inventory values arrive per-component in the biology phases (after GATE A banks).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

__all__ = [
    "Rate",
    "Distribution",
    "CellState",
    "FORBIDDEN_OUTPUT_KINDS",
    "mcf7_reference_state",
    "assert_cellstate_writes_no_force",
]

# Output kinds CellState may NEVER emit — these are mechanical/topological state that must EMERGE from
# the transaction, never be set from the slow-context state vector.
FORBIDDEN_OUTPUT_KINDS: frozenset[str] = frozenset(
    {
        "force",
        "pos",
        "position",
        "tension",
        "count",
        "filament_count",
        "bound_head_count",
        "crosslink_topology",
        "stress",
    }
)

# Method-name fragments that would let a caller treat CellState as a force/count writer.
_FORBIDDEN_WRITER_HINTS: frozenset[str] = frozenset(
    {"force", "pos", "position", "tension", "count", "stress", "topology"}
)


def _require_kb_source(kb_source: str, *, what: str) -> str:
    """Reject a value without KB provenance (Magic-Number Block)."""
    if not isinstance(kb_source, str) or not kb_source.strip():
        raise ValueError(f"{what} needs a non-empty KB provenance (no unsourced magic numbers)")
    return kb_source


def _reject_forbidden_kind(kind: str, *, what: str) -> None:
    """Reject a CellState output whose kind is emergent mechanical/topological state."""
    if kind in FORBIDDEN_OUTPUT_KINDS:
        raise ValueError(
            f"CellState may not emit {kind!r} as {what}: cortex tension / contractile force / "
            "bound-head or filament count / crosslink topology must EMERGE from mechanics + events"
        )


@dataclass(frozen=True, slots=True)
class Rate:
    """A reaction rate CellState conditions, with its mandatory KB provenance.

    Args:
        reaction: the reaction label (e.g. ``"nucleation"``, ``"nmii_activation"``) — never a forbidden kind.
        value: the rate value [per unit time], nonnegative.
        kb_source: non-empty KB provenance (``KB-x.y`` / DOI); an unsourced value is rejected.
    """

    reaction: str
    value: float
    kb_source: str

    def __post_init__(self) -> None:
        if not self.reaction.strip():
            raise ValueError("rate needs a reaction label")
        _reject_forbidden_kind(self.reaction, what="a rate")
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise TypeError("rate value must be a real number")
        if not (self.value >= 0.0) or self.value != self.value:  # reject negatives and NaN
            raise ValueError("rate value must be finite and nonnegative")
        _require_kb_source(self.kb_source, what="rate")


@dataclass(frozen=True, slots=True)
class Distribution:
    """An initial-inventory distribution CellState seeds, with its mandatory KB provenance.

    Args:
        kind: the inventory kind (e.g. ``"length"``, ``"spatial_density"``) — never a forbidden kind.
        params: distribution parameters (e.g. ``{"mean": ..., "scale": ...}``); provenance-carried, opaque here.
        kb_source: non-empty KB provenance; an unsourced distribution is rejected.
    """

    kind: str
    params: Mapping[str, float]
    kb_source: str

    def __post_init__(self) -> None:
        if not self.kind.strip():
            raise ValueError("distribution needs a kind")
        _reject_forbidden_kind(self.kind, what="an inventory distribution")
        if not isinstance(self.params, Mapping):
            raise TypeError("distribution params must be a mapping")
        object.__setattr__(self, "params", MappingProxyType(dict(self.params)))
        _require_kb_source(self.kb_source, what="distribution")


@dataclass(frozen=True, slots=True)
class CellState:
    """Slow-context state vector conditioning rates + inventory; owns no physical state.

    The six axes are fixed for a run.  ``rates`` maps ``(component, reaction) -> Rate``; ``inventory``
    maps ``component -> Distribution``.  Both are frozen at construction; CellState exposes only the
    read accessors below, so it can never write a mechanical or topological array.
    """

    lineage: str
    emt: str
    cycle: str
    adhesion: str
    geometry: str
    osmotic_state: str
    rates: Mapping[tuple[str, str], Rate] = field(default_factory=dict)
    inventory: Mapping[str, Distribution] = field(default_factory=dict)

    def __post_init__(self) -> None:
        axes = {
            "lineage": self.lineage,
            "emt": self.emt,
            "cycle": self.cycle,
            "adhesion": self.adhesion,
            "geometry": self.geometry,
            "osmotic_state": self.osmotic_state,
        }
        for name, value in axes.items():
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"CellState axis {name!r} must be a non-empty sourced label")
        for key, rate in self.rates.items():
            if not (isinstance(key, tuple) and len(key) == 2):
                raise TypeError("CellState.rates must be keyed by (component, reaction)")
            if not isinstance(rate, Rate):
                raise TypeError(f"CellState rate for {key!r} must be a Rate (carries KB provenance)")
            if rate.reaction != key[1]:
                raise ValueError(f"rate key {key!r} disagrees with its reaction label {rate.reaction!r}")
        for component, dist in self.inventory.items():
            if not isinstance(dist, Distribution):
                raise TypeError(f"CellState inventory for {component!r} must be a Distribution")
        object.__setattr__(self, "rates", MappingProxyType(dict(self.rates)))
        object.__setattr__(self, "inventory", MappingProxyType(dict(self.inventory)))

    @property
    def axes(self) -> tuple[str, str, str, str, str, str]:
        """The six ratified axis values in canonical order."""
        return (
            self.lineage,
            self.emt,
            self.cycle,
            self.adhesion,
            self.geometry,
            self.osmotic_state,
        )

    def rate(self, component: str, reaction: str) -> Rate:
        """Return the sourced :class:`Rate` for ``(component, reaction)`` or raise if unsourced.

        Raising (rather than defaulting) is the sourced-only discipline: an unsourced reaction must be
        added with a KB value and PI sign-off, not silently zero/one.
        """
        try:
            return self.rates[(component, reaction)]
        except KeyError as error:
            raise KeyError(
                f"no sourced rate for ({component!r}, {reaction!r}) under CellState {self.axes}; "
                "add a KB-sourced value (Magic-Number Block) + PI sign-off before use"
            ) from error

    def initial_inventory(self, component: str) -> Distribution:
        """Return the sourced initial-inventory :class:`Distribution` for ``component`` or raise."""
        try:
            return self.inventory[component]
        except KeyError as error:
            raise KeyError(
                f"no sourced initial inventory for {component!r} under CellState {self.axes}; "
                "add a KB-sourced distribution + PI sign-off before use"
            ) from error


def mcf7_reference_state() -> CellState:
    """Return the first ratified sourced CellState axes instance (empty registries; biology fills them).

    ``MCF7 × hybrid_E_M × G1 × collagen_spread × polarized × physiological_resting`` — the axes list is
    the Step-1 contract PI ratifies; the per-component sourced rates/inventory arrive after GATE A banks.
    """
    return CellState(
        lineage="MCF7",
        emt="hybrid_E_M",
        cycle="G1",
        adhesion="collagen_spread",
        geometry="polarized",
        osmotic_state="physiological_resting",
    )


def assert_cellstate_writes_no_force(state: object) -> None:
    """Assert a CellState instance exposes no device array and no force/count-writing setter.

    The structural proof of ``CellState → force`` being a rejected build: the object must carry no
    ``wp.array`` attribute and expose no callable whose name implies writing a force / position /
    tension / count / stress / topology array.  Its only value-emitting accessors are ``rate`` and
    ``initial_inventory``, which return :class:`Rate` / :class:`Distribution`.

    Raises:
        AssertionError: if the instance exposes a forbidden writer or holds a device array.
    """
    for attr in dir(state):
        if attr.startswith("__"):
            continue
        lowered = attr.lower()
        member = getattr(state, attr, None)
        if callable(member) and lowered.startswith(("set_", "write_", "assign_", "apply_")):
            if any(hint in lowered for hint in _FORBIDDEN_WRITER_HINTS):
                raise AssertionError(
                    f"CellState exposes forbidden mechanical/topological writer {attr!r}"
                )
        # A device array attribute would mean CellState owns physical state — it must not.
        device = getattr(member, "device", None)
        if getattr(member, "shape", None) is not None and device is not None:
            raise AssertionError(
                f"CellState holds a device array {attr!r}; it must own no physical state"
            )
