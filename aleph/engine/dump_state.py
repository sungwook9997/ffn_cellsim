"""Composed-world state dump for the component-first Active Cell engine.

This is the composed-engine-world analog of :mod:`aleph.components.incumbent.dump_state` (which dumps the legacy
monolithic ``AssembledCell``).  It walks a :class:`~aleph.engine.composition.ComposedCellWorld` and, for
each *registered* component state-owner and connector runtime, records the per-compartment geometry census the
per-compartment viz needs — node / face / connector-endpoint counts — together with a **global-unique
actor / filament ID namespace**.

The ID namespace is the anti-double-draw contract:

* Each distinct runtime object is one actor with one ``actor_id``; the composite α2β1–collagen FA joint (bound
  under both ``fa_actin_anchor`` and ``integrin_collagen_clutch``) is deduplicated by object identity, so it is
  a single actor rather than two.  The extra edge name(s) are recorded as ``aliases``.
* Stress fibres (``sf_arc``) are a separate actor from the cortex, the lamellipodium, etc., so each carries a
  disjoint ``[filament_id_base, filament_id_base + n_filaments)`` block of globally unique filament IDs — a
  per-compartment renderer keyed on these IDs never draws one physical filament under two labels.

This module allocates no device memory, launches no kernel, and performs **no device-to-host readback**: it
reads only host-side array *shape* metadata and plain integer count attributes, so it is safe to call inside a
zero-DtoH physical step and importable on a CUDA-free host.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from aleph.engine.composition import ComposedCellWorld

__all__ = [
    "ComposedWorldEntryDump",
    "ComposedWorldDump",
    "dump_composed_world_state",
]

# Documented public-attribute precedence used when a runtime does not implement ``dump_state_entry``.
# These read host-side metadata only (plain ints or ``wp.array.shape``); none triggers a device readback.
_FILAMENT_COUNT_ATTRS = ("n_filaments", "n_minifilaments", "n_fibers", "n_active_filaments")
_NODE_COUNT_ATTRS = ("n_nodes", "n_particles", "n_verts")
_NODE_ARRAY_ATTRS = ("position_d", "pos_d")
_FACE_COUNT_ATTRS = ("n_faces",)
_FACE_ARRAY_ATTRS = ("faces_d",)
_ENDPOINT_COUNT_ATTRS = ("n_endpoints", "n_heads", "n_segments")


@dataclass(frozen=True, slots=True)
class ComposedWorldEntryDump:
    """One composed-world actor's per-compartment census + its global-unique ID block."""

    kind: str
    name: str
    aliases: tuple[str, ...]
    runtime_type: str
    actor_id: int
    filament_id_base: int
    n_filaments: int
    n_nodes: int
    n_faces: int
    n_endpoints: int
    opaque: bool

    @property
    def filament_id_range(self) -> tuple[int, int]:
        """Return the half-open ``[base, base + n_filaments)`` global filament-ID block for this actor."""
        return (self.filament_id_base, self.filament_id_base + self.n_filaments)


@dataclass(frozen=True, slots=True)
class ComposedWorldDump:
    """The whole composed cell's per-compartment census with a disjoint global ID namespace."""

    entries: tuple[ComposedWorldEntryDump, ...]
    n_actors: int
    n_filaments_total: int
    n_nodes_total: int

    def by_name(self, name: str) -> ComposedWorldEntryDump:
        """Return the entry whose primary name or alias is ``name``."""
        for entry in self.entries:
            if entry.name == name or name in entry.aliases:
                return entry
        raise KeyError(name)

    def actor_ids_are_unique(self) -> bool:
        """Return whether every entry carries a distinct actor id (the no-double-draw invariant)."""
        ids = [entry.actor_id for entry in self.entries]
        return len(set(ids)) == len(ids)

    def filament_id_blocks_are_disjoint(self) -> bool:
        """Return whether the per-actor global filament-ID blocks never overlap."""
        blocks = sorted(
            entry.filament_id_range for entry in self.entries if entry.n_filaments > 0
        )
        return all(blocks[i][1] <= blocks[i + 1][0] for i in range(len(blocks) - 1))


def _shape0(array: object) -> int | None:
    """Return ``array.shape[0]`` (host-side metadata) or ``None`` if unavailable — never reads device data."""
    shape = getattr(array, "shape", None)
    if isinstance(shape, tuple) and shape:
        try:
            return int(shape[0])
        except (TypeError, ValueError):
            return None
    return None


def _first_int_attr(runtime: object, names: tuple[str, ...]) -> int | None:
    for name in names:
        value = getattr(runtime, name, None)
        if isinstance(value, int) and not isinstance(value, bool):
            return value
    return None


def _first_array_len(runtime: object, names: tuple[str, ...]) -> int | None:
    for name in names:
        length = _shape0(getattr(runtime, name, None))
        if length is not None:
            return length
    return None


def _census(runtime: object) -> tuple[int, int, int, int, bool]:
    """Return ``(n_filaments, n_nodes, n_faces, n_endpoints, opaque)`` for one runtime.

    A runtime may implement ``dump_state_entry() -> Mapping`` to declare its own census exactly; otherwise a
    documented public-attribute precedence is used, and ``opaque`` marks a runtime that exposed none of them.
    """
    hook = getattr(runtime, "dump_state_entry", None)
    if callable(hook):
        declared = hook() or {}
        return (
            int(declared.get("n_filaments", 0)),
            int(declared.get("n_nodes", 0)),
            int(declared.get("n_faces", 0)),
            int(declared.get("n_endpoints", 0)),
            False,
        )
    n_filaments = _first_int_attr(runtime, _FILAMENT_COUNT_ATTRS) or 0
    n_nodes = _first_int_attr(runtime, _NODE_COUNT_ATTRS)
    if n_nodes is None:
        n_nodes = _first_array_len(runtime, _NODE_ARRAY_ATTRS)
    n_faces = _first_int_attr(runtime, _FACE_COUNT_ATTRS)
    if n_faces is None:
        n_faces = _first_array_len(runtime, _FACE_ARRAY_ATTRS)
    n_endpoints = _first_int_attr(runtime, _ENDPOINT_COUNT_ATTRS)
    if n_endpoints is None:
        state = getattr(runtime, "state", None)
        if state is not None:
            n_endpoints = _first_int_attr(state, _ENDPOINT_COUNT_ATTRS)
    opaque = n_filaments == 0 and n_nodes is None and n_faces is None and n_endpoints is None
    return (n_filaments, n_nodes or 0, n_faces or 0, n_endpoints or 0, opaque)


def dump_composed_world_state(world: ComposedCellWorld) -> ComposedWorldDump:
    """Census every registered component/connector runtime with a disjoint global actor/filament ID namespace.

    Components are visited in declared architecture order, then connectors; each *distinct* runtime object is
    assigned one monotonic ``actor_id`` and a contiguous global filament-ID block, and any further edge name
    that shares the same object is recorded as an alias (the composite FA joint is one actor, not two).
    """
    architecture = world.architecture
    actor = world.actor

    ordered: list[tuple[str, str]] = [
        ("component", component.name)
        for component in architecture.components
        if component.name not in set(actor.missing_bindings()["components"])
    ]
    missing_connectors = set(actor.missing_bindings()["connectors"])
    ordered += [
        ("connector", connector.name)
        for connector in architecture.connectors
        if connector.name not in missing_connectors
    ]

    entries: list[ComposedWorldEntryDump] = []
    by_identity: dict[int, int] = {}
    filament_cursor = 0

    for kind, name in ordered:
        runtime = (
            actor.component_runtime(name)
            if kind == "component"
            else actor.connector_runtime(name)
        )
        identity = id(runtime)
        if identity in by_identity:
            index = by_identity[identity]
            existing = entries[index]
            entries[index] = ComposedWorldEntryDump(
                kind=existing.kind,
                name=existing.name,
                aliases=(*existing.aliases, name),
                runtime_type=existing.runtime_type,
                actor_id=existing.actor_id,
                filament_id_base=existing.filament_id_base,
                n_filaments=existing.n_filaments,
                n_nodes=existing.n_nodes,
                n_faces=existing.n_faces,
                n_endpoints=existing.n_endpoints,
                opaque=existing.opaque,
            )
            continue
        n_filaments, n_nodes, n_faces, n_endpoints, opaque = _census(runtime)
        entry = ComposedWorldEntryDump(
            kind=kind,
            name=name,
            aliases=(),
            runtime_type=type(runtime).__name__,
            actor_id=len(entries),
            filament_id_base=filament_cursor,
            n_filaments=n_filaments,
            n_nodes=n_nodes,
            n_faces=n_faces,
            n_endpoints=n_endpoints,
            opaque=opaque,
        )
        by_identity[identity] = len(entries)
        entries.append(entry)
        filament_cursor += n_filaments

    return ComposedWorldDump(
        entries=tuple(entries),
        n_actors=len(entries),
        n_filaments_total=filament_cursor,
        n_nodes_total=sum(entry.n_nodes for entry in entries),
    )
