"""The scene tree: components and connectors, kept apart.

Manuscript §11: *"Scene tree exposes components and connectors separately with ownership, endpoints,
force/work, events, and isolation."*

"Separately" is the requirement that does the work. The tempting design is one node type with an
optional second endpoint, because a connector is nearly a component with two ends. That design loses
the two things a viewer most needs to show: a component *owns* state and a connector *does not*, and a
connector's force is a pair that must sum to zero while a component's force is a resultant that need
not. Merge them and the tree can no longer say which of those two rules applies to a given node, so it
applies neither, and the render stops being able to show a broken adjoint pair at all.

So there are two node types with different fields, two accessors, and a constructor that refuses a
name used by both. The tree is built from :class:`~aleph.runtime.participant.EntityCensus` and
:class:`~aleph.runtime.participant.AdjointPair` — the runtime's own types — so ownership and force
come from what the step actually reported rather than from a parallel description of it.

Isolation is a first-class state rather than a boolean. A connector that is disengaged this step, a
connector that was never bound into the schedule, and a component deliberately excluded from the run
are three different situations, and a viewer that shows all three as "greyed out" hides the one that
matters: `NOT_SCHEDULED` is the defect this project was restarted over.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum, unique
from typing import Any, Iterable, Mapping, Sequence

from aleph.runtime.participant import AdjointPair, EntityCensus, norm

__all__ = [
    "ComponentNode",
    "ConnectorNode",
    "IsolationState",
    "SceneTree",
    "declare_unevaluated_connector",
    "scene_from_step",
]


@unique
class IsolationState(StrEnum):
    """Why a node is or is not participating in the load path this step."""

    COUPLED = "coupled"
    """Evaluated this step and carrying load."""

    ENGAGED_ZERO = "engaged_zero"
    """Evaluated this step and carrying exactly zero — a unilateral element on its inactive branch.

    Distinct from ``COUPLED`` with a small force and distinct from ``NOT_SCHEDULED``. A tensile-only
    tether under compression belongs here, and the fact that it reports an *exact* zero rather than a
    small number is the evidence that it is on the inactive branch rather than nearly balanced.
    """

    ISOLATED = "isolated"
    """Deliberately decoupled by the run's configuration — a break mode or a disabled connector."""

    NOT_SCHEDULED = "not_scheduled"
    """Registered and never evaluated. The failure this project exists to make visible.

    A viewer must render this differently from ``ISOLATED``: isolated is a decision, not scheduled is
    a hole in the schedule that nothing else in the system will report.
    """

    EXCLUDED = "excluded"
    """Declared out of scope, with a reason. Not a fault."""

    @property
    def carries_load(self) -> bool:
        """Whether this state means load is actually being transmitted."""
        return self is IsolationState.COUPLED

    @property
    def is_fault(self) -> bool:
        """Whether this state is a defect rather than a decision."""
        return self is IsolationState.NOT_SCHEDULED

    @property
    def requires_reason(self) -> bool:
        """Whether a node in this state must say why."""
        return self in (
            IsolationState.ISOLATED,
            IsolationState.NOT_SCHEDULED,
            IsolationState.EXCLUDED,
        )


@dataclass(frozen=True, slots=True)
class ComponentNode:
    """One state-owning component.

    Attributes:
        name: The owner's registration name.
        entities: The entities it owns.
        element_counts: Population of each entity, index-aligned with ``entities``.
        resultant_force_pn: Net force on the component [pN], when the step reported one. A resultant
            legitimately vanishes on a correct closed body, so this is displayed and never gated on.
        isolation: Participation state.
        isolation_reason: Required when :attr:`IsolationState.requires_reason`.
        events: Typed events at this step.
    """

    name: str
    entities: tuple[str, ...]
    element_counts: tuple[int, ...]
    resultant_force_pn: float | None = None
    isolation: IsolationState = IsolationState.COUPLED
    isolation_reason: str = ""
    events: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("a component node needs a name")
        if len(self.entities) != len(self.element_counts):
            raise ValueError(
                f"component {self.name!r} lists {len(self.entities)} entities and "
                f"{len(self.element_counts)} counts"
            )
        if self.isolation.requires_reason and not self.isolation_reason.strip():
            raise ValueError(
                f"component {self.name!r} is {self.isolation.value} and gives no reason; an "
                "unexplained isolation is indistinguishable from a scheduling defect"
            )

    @property
    def total_elements(self) -> int:
        """Sum over owned entities."""
        return sum(self.element_counts)

    @property
    def owns_state(self) -> bool:
        """Always ``True``. Present so the tree's two node types answer the same question."""
        return True

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able component node."""
        return {
            "kind": "component",
            "name": self.name,
            "owns_state": True,
            "entities": list(self.entities),
            "element_counts": list(self.element_counts),
            "total_elements": self.total_elements,
            "resultant_force_pn": self.resultant_force_pn,
            "isolation": self.isolation.value,
            "isolation_reason": self.isolation_reason,
            "isolation_is_fault": self.isolation.is_fault,
            "events": list(self.events),
        }


@dataclass(frozen=True, slots=True)
class ConnectorNode:
    """One coupling, with both endpoints and both sides of what it did.

    Attributes:
        name: The connector's registration name.
        endpoint_a: Owner on side a.
        endpoint_b: Owner on side b.
        internal_state_owner: Who owns the connector's own internal state (rest lengths, bound flags).
            A connector owns that and owns neither endpoint's state, and the tree says so explicitly
            because the reference failure was exactly a confusion about who owned what.
        force_a_pn: Magnitude of the force applied to endpoint a [pN].
        force_b_pn: Magnitude of the force applied to endpoint b [pN].
        force_residual_pn: ``|force_a + force_b|`` [pN]. Displayed, because a non-zero residual on a
            two-force member is the connector failing.
        endpoint_work_pn_um: Work done on the endpoints over the step [pN.um].
        stored_energy_delta_pn_um: Change in the connector's internal stored energy [pN.um].
        dissipation_pn_um: Energy converted to heat [pN.um].
        energy_residual_pn_um: ``W + dU + D - A`` [pN.um]. Zero when the books close.
        isolation: Participation state.
        isolation_reason: Required when :attr:`IsolationState.requires_reason`.
        events: Typed events at this step (bound, unbound, severed).
    """

    name: str
    endpoint_a: str
    endpoint_b: str
    internal_state_owner: str
    force_a_pn: float = 0.0
    force_b_pn: float = 0.0
    force_residual_pn: float = 0.0
    endpoint_work_pn_um: float = 0.0
    stored_energy_delta_pn_um: float = 0.0
    dissipation_pn_um: float = 0.0
    energy_residual_pn_um: float = 0.0
    isolation: IsolationState = IsolationState.COUPLED
    isolation_reason: str = ""
    events: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("name", "endpoint_a", "endpoint_b", "internal_state_owner"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"ConnectorNode.{name} is required")
        if self.endpoint_a == self.endpoint_b:
            raise ValueError(
                f"connector {self.name!r} names {self.endpoint_a!r} on both sides; a coupling with "
                "one endpoint is an internal force and belongs to the component that owns it"
            )
        for name in (
            "force_a_pn",
            "force_b_pn",
            "force_residual_pn",
            "endpoint_work_pn_um",
            "stored_energy_delta_pn_um",
            "dissipation_pn_um",
            "energy_residual_pn_um",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"connector {self.name!r} reports a non-finite {name}")
            object.__setattr__(self, name, value)
        if self.isolation.requires_reason and not self.isolation_reason.strip():
            raise ValueError(
                f"connector {self.name!r} is {self.isolation.value} and gives no reason"
            )
        if self.isolation is IsolationState.ENGAGED_ZERO and (
            self.force_a_pn != 0.0 or self.force_b_pn != 0.0
        ):
            raise ValueError(
                f"connector {self.name!r} claims ENGAGED_ZERO while reporting "
                f"({self.force_a_pn}, {self.force_b_pn}) pN. That state means the element is on its "
                "inactive branch and its force is identically zero, not merely small."
            )

    @property
    def owns_state(self) -> bool:
        """Always ``False`` for endpoint state.

        A connector owns internal state — that is what :attr:`internal_state_owner` names — and owns
        neither endpoint's arrays. Answering ``False`` here is the statement that it is not an owner
        in the census sense.
        """
        return False

    @property
    def endpoints(self) -> tuple[str, str]:
        """The two owners this connector couples."""
        return (self.endpoint_a, self.endpoint_b)

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able connector node."""
        return {
            "kind": "connector",
            "name": self.name,
            "owns_state": False,
            "endpoint_a": self.endpoint_a,
            "endpoint_b": self.endpoint_b,
            "internal_state_owner": self.internal_state_owner,
            "force_a_pn": self.force_a_pn,
            "force_b_pn": self.force_b_pn,
            "force_residual_pn": self.force_residual_pn,
            "endpoint_work_pn_um": self.endpoint_work_pn_um,
            "stored_energy_delta_pn_um": self.stored_energy_delta_pn_um,
            "dissipation_pn_um": self.dissipation_pn_um,
            "energy_residual_pn_um": self.energy_residual_pn_um,
            "isolation": self.isolation.value,
            "isolation_reason": self.isolation_reason,
            "isolation_is_fault": self.isolation.is_fault,
            "events": list(self.events),
        }


@dataclass(frozen=True, slots=True)
class SceneTree:
    """Components and connectors, side by side and never merged.

    Attributes:
        accepted_step: Which accepted step this tree describes.
        topology_epoch: Topology generation, so a tree cannot be drawn against the wrong connectivity.
        components: The state-owning components.
        connectors: The couplings.
    """

    accepted_step: int
    topology_epoch: int
    components: tuple[ComponentNode, ...] = ()
    connectors: tuple[ConnectorNode, ...] = ()

    def __post_init__(self) -> None:
        component_names = [node.name for node in self.components]
        connector_names = [node.name for node in self.connectors]
        for label, names in (("component", component_names), ("connector", connector_names)):
            if len(set(names)) != len(names):
                raise ValueError(f"a {label} name appears twice: {names}")
        overlap = sorted(set(component_names) & set(connector_names))
        if overlap:
            raise ValueError(
                f"{overlap} appear as both a component and a connector. The two carry different "
                "rules — a component's resultant may vanish, a connector's pair must — so one name "
                "meaning both makes the tree unable to say which rule applies."
            )
        for connector in self.connectors:
            unknown = [
                owner for owner in connector.endpoints if owner not in set(component_names)
            ]
            if unknown:
                raise ValueError(
                    f"connector {connector.name!r} names endpoint owner(s) {unknown} that are not "
                    "components of this tree; an endpoint with no owner is an edge to nowhere"
                )

    # -- separate accessors, by design --------------------------------------------------------
    def component(self, name: str) -> ComponentNode:
        """The component named ``name``.

        Raises:
            KeyError: If no component has that name. If a *connector* has it, the message says so,
                because reaching for a connector through the component accessor is the exact
                confusion this tree is shaped to prevent.
        """
        for node in self.components:
            if node.name == name:
                return node
        if any(node.name == name for node in self.connectors):
            raise KeyError(
                f"{name!r} is a connector, not a component; use SceneTree.connector(). A connector "
                "owns no endpoint state and its force is a pair, not a resultant."
            )
        raise KeyError(f"no component named {name!r}")

    def connector(self, name: str) -> ConnectorNode:
        """The connector named ``name``.

        Raises:
            KeyError: If no connector has that name.
        """
        for node in self.connectors:
            if node.name == name:
                return node
        if any(node.name == name for node in self.components):
            raise KeyError(f"{name!r} is a component, not a connector; use SceneTree.component()")
        raise KeyError(f"no connector named {name!r}")

    # -- the queries a viewer actually needs ---------------------------------------------------
    @property
    def faults(self) -> tuple[str, ...]:
        """Every node whose isolation state is a defect rather than a decision."""
        return tuple(
            f"{node.name}: {node.isolation.value} — {node.isolation_reason}"
            for node in (*self.components, *self.connectors)
            if node.isolation.is_fault
        )

    @property
    def load_path(self) -> tuple[str, ...]:
        """Connectors actually transmitting load this step."""
        return tuple(node.name for node in self.connectors if node.isolation.carries_load)

    def worst_force_residual(self) -> float:
        """Largest connector force-closure residual in the tree [pN]."""
        return max((node.force_residual_pn for node in self.connectors), default=0.0)

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able tree, with the two node kinds in two separate arrays."""
        return {
            "accepted_step": self.accepted_step,
            "topology_epoch": self.topology_epoch,
            "components": [node.as_json_obj() for node in self.components],
            "connectors": [node.as_json_obj() for node in self.connectors],
            "faults": list(self.faults),
            "load_path": list(self.load_path),
        }


def scene_from_step(
    *,
    accepted_step: int,
    topology_epoch: int,
    censuses: Iterable[EntityCensus],
    pairs: Iterable[AdjointPair] = (),
    isolation: Mapping[str, tuple[IsolationState, str]] | None = None,
    events: Mapping[str, Sequence[str]] | None = None,
    resultant_forces_pn: Mapping[str, float] | None = None,
    connector_internal_owners: Mapping[str, str] | None = None,
) -> SceneTree:
    """Build a scene tree from what one step reported.

    The inputs are the runtime's own types, which is what makes the tree an *observation* of the step
    rather than a second description of the model kept in step by hand.

    Args:
        accepted_step: Accepted-step index.
        topology_epoch: Topology generation.
        censuses: One :class:`EntityCensus` per state-owning component.
        pairs: The :class:`AdjointPair` each evaluated connector returned.
        isolation: Per-node ``(state, reason)`` overrides. A connector that appears in ``isolation``
            with ``NOT_SCHEDULED`` and does *not* appear in ``pairs`` is how the registered-but-never-
            evaluated case is rendered.
        events: Per-node event labels.
        resultant_forces_pn: Per-component resultant force magnitude, when the step measured one.
        connector_internal_owners: Which owner holds each connector's internal state. Defaults to the
            connector's own name, which is the usual case: a connector owns its own rest lengths.

    Returns:
        The tree.
    """
    isolation = dict(isolation or {})
    events = {key: tuple(value) for key, value in (events or {}).items()}
    resultants = dict(resultant_forces_pn or {})
    internal_owners = dict(connector_internal_owners or {})

    components: list[ComponentNode] = []
    for census in censuses:
        state, reason = isolation.get(census.owner, (IsolationState.COUPLED, ""))
        components.append(
            ComponentNode(
                name=census.owner,
                entities=tuple(census.entities),
                element_counts=tuple(census.element_counts),
                resultant_force_pn=resultants.get(census.owner),
                isolation=state,
                isolation_reason=reason,
                events=events.get(census.owner, ()),
            )
        )

    connectors: list[ConnectorNode] = []
    evaluated: set[str] = set()
    for pair in pairs:
        evaluated.add(pair.connector)
        force_a = norm(pair.force_a)
        force_b = norm(pair.force_b)
        default_state = (
            IsolationState.ENGAGED_ZERO
            if force_a == 0.0 and force_b == 0.0
            else IsolationState.COUPLED
        )
        default_reason = (
            "evaluated this step and reported an exact zero, i.e. a unilateral element on its "
            "inactive branch"
            if default_state is IsolationState.ENGAGED_ZERO
            else ""
        )
        state, reason = isolation.get(pair.connector, (default_state, default_reason))
        connectors.append(
            ConnectorNode(
                name=pair.connector,
                endpoint_a=pair.endpoint_a,
                endpoint_b=pair.endpoint_b,
                internal_state_owner=internal_owners.get(pair.connector, pair.connector),
                force_a_pn=force_a,
                force_b_pn=force_b,
                force_residual_pn=norm(pair.force_residual),
                endpoint_work_pn_um=pair.endpoint_work,
                stored_energy_delta_pn_um=pair.stored_energy_delta,
                dissipation_pn_um=pair.dissipation,
                energy_residual_pn_um=pair.energy_residual,
                isolation=state,
                isolation_reason=reason,
                events=events.get(pair.connector, ()),
            )
        )

    return SceneTree(
        accepted_step=int(accepted_step),
        topology_epoch=int(topology_epoch),
        components=tuple(components),
        connectors=tuple(connectors),
    )


def declare_unevaluated_connector(
    *,
    name: str,
    endpoint_a: str,
    endpoint_b: str,
    reason: str,
    internal_state_owner: str | None = None,
    events: Sequence[str] = (),
) -> ConnectorNode:
    """Build the node for a connector that was registered and never evaluated.

    Kept as a named constructor rather than left to callers, because every field except the isolation
    state has to be zero and it must be *visibly* zero rather than absent. A node built this way shows
    the viewer a connector with no force, no work and an explicit ``NOT_SCHEDULED``, which is the only
    honest rendering: a connector that did nothing did not thereby carry no load, it was never asked.
    """
    return ConnectorNode(
        name=name,
        endpoint_a=endpoint_a,
        endpoint_b=endpoint_b,
        internal_state_owner=internal_state_owner or name,
        isolation=IsolationState.NOT_SCHEDULED,
        isolation_reason=reason,
        events=tuple(events),
    )
