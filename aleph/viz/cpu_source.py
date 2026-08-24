"""A small CPU producer, so the viewer can be exercised before the trajectory closes.

Manuscript §11's dependency-boundary row says the packing, the viewer process and the UI *can be
prototyped before the whole-cell trajectory closes*. This module is the producer that makes that true
for Aleph: a real :class:`~aleph.runtime.transaction.Transaction` over a real
:class:`~aleph.runtime.pipeline.Pipeline`, with one owner doing arithmetic through the backend and one
connector that actually transmits load.

What it is **not**: it is not physics anyone should quote. The relaxation here is a linear pull toward
a target, chosen because it converges monotonically and deterministically, which makes it a good thing
to render and a useless thing to measure. Nothing in this module is an oracle, a model, or evidence,
and the viewer's job is to show whatever the producer holds — which means the producer's own honesty is
not what the viewer lane is testing.

It exists so that :mod:`aleph.viz` has a producer it owns. The runtime lane's test doubles are the
right thing to test the *runtime* against; a script that a person runs needs a producer inside the
package rather than inside a test tree.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from aleph.runtime.backend import NumpyBackend
from aleph.runtime.ledger import ForceWorkLedger
from aleph.runtime.participant import (
    AdjointPair,
    EndpointHandle,
    EntityCensus,
    Phase,
    StepContext,
    WorkReceipt,
)
from aleph.runtime.pipeline import Pipeline
from aleph.runtime.transaction import Clock, RngStream, Topology, Transaction

__all__ = [
    "AnchorSpring",
    "RelaxingField",
    "ViewerDemoWorld",
    "build_demo_world",
]


@dataclass(slots=True)
class RelaxingField:
    """One owner holding node positions and a per-node force accumulator.

    Each step it pulls every node a fixed fraction of the way toward its target, through the backend so
    the coverage gate's evaluation witness genuinely moves. Deterministic, monotone, and interesting
    enough to look at.
    """

    name: str
    entity: str
    positions: np.ndarray
    targets: np.ndarray
    stiffness: float = 40.0
    rate: float = 5.0e-3
    forces: np.ndarray = field(init=False)
    _snapshot: np.ndarray | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        self.positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        self.targets = np.ascontiguousarray(self.targets, dtype=np.float64)
        if self.positions.shape != self.targets.shape:
            raise ValueError("positions and targets must have the same shape")
        self.forces = np.zeros_like(self.positions)

    def owned_entities(self) -> EntityCensus:
        """The census of what this owner holds."""
        return EntityCensus(
            owner=self.name,
            entities=(self.entity,),
            element_counts=(int(self.positions.shape[0]),),
        )

    def snapshot(self) -> None:
        """Capture positions so a rejected candidate can be undone exactly."""
        self._snapshot = self.positions.copy()

    def commit(self) -> None:
        """Drop the snapshot; the candidate is now authoritative."""
        self._snapshot = None

    def rollback(self) -> None:
        """Restore positions exactly."""
        if self._snapshot is None:
            raise AssertionError(f"{self.name}: rollback without a snapshot")
        self.positions[:] = self._snapshot
        self._snapshot = None

    def accumulate(self, phase: Phase, ctx: StepContext) -> WorkReceipt:
        """Do this phase's work through the backend and report it."""
        backend = ctx.backend
        if phase is Phase.ACCUMULATE_INTERNAL:
            backend.scale(self.forces, 0.0, out=self.forces)
            pull = backend.scale(backend.subtract(self.targets, self.positions), self.stiffness)
            backend.add(self.forces, pull, out=self.forces)
            ops, detail = 4, "restoring force toward the target configuration"
        elif phase is Phase.SOLVE:
            backend.add(self.positions, backend.scale(self.forces, self.rate), out=self.positions)
            ops, detail = 2, "advanced candidate positions"
        else:
            ops = 1
            detail = f"summed the accumulator ({float(backend.sum(self.forces)):.3e})"
        return WorkReceipt(
            participant=self.name,
            phase=phase,
            touched_entities=(self.entity,),
            declared_ops=ops,
            detail=detail,
        )


@dataclass(slots=True)
class AnchorSpring:
    """A passive linear spring between one node of each owner.

    Each side's force is computed from that side's own geometry rather than by negating the other's, so
    the ledger's force-closure check remains a check.
    """

    name: str
    endpoints: tuple[str, str]
    stiffness: float = 120.0
    rest_length: float = 0.5
    last_pair: AdjointPair | None = field(default=None, init=False, repr=False)

    def owned_entities(self) -> EntityCensus:
        """A connector owns its internal state, never its endpoints'."""
        return EntityCensus(owner=self.name, entities=("spring",), element_counts=(1,))

    def snapshot(self) -> None:
        """No mutable internal state to capture."""

    def commit(self) -> None:
        """No mutable internal state to commit."""

    def rollback(self) -> None:
        """No mutable internal state to restore."""

    def accumulate(
        self,
        phase: Phase,
        ctx: StepContext,
        endpoint_a: EndpointHandle,
        endpoint_b: EndpointHandle,
    ) -> AdjointPair:
        """Apply the spring force to both endpoints and report both sides."""
        r_a = np.array(endpoint_a.position)
        r_b = np.array(endpoint_b.position)
        separation = r_a - r_b
        length = float(np.linalg.norm(separation))
        if length == 0.0:
            raise ValueError(
                f"connector {self.name!r} has coincident endpoints, so the line of action is "
                "undefined; refusing rather than returning a zero force that looks like equilibrium"
            )
        u_ab = separation / length
        extension = length - self.rest_length
        force_a = -self.stiffness * extension * u_ab
        force_b = -self.stiffness * extension * (-u_ab)

        if endpoint_a.force_sink is not None:
            ctx.backend.scatter_add(endpoint_a.force_sink, [endpoint_a.site_index], [force_a])
        if endpoint_b.force_sink is not None:
            ctx.backend.scatter_add(endpoint_b.force_sink, [endpoint_b.site_index], [force_b])

        pair = AdjointPair(
            connector=self.name,
            endpoint_a=endpoint_a.owner,
            endpoint_b=endpoint_b.owner,
            force_a=tuple(force_a),
            force_b=tuple(force_b),
            point_a=tuple(r_a),
            point_b=tuple(r_b),
            stored_energy_delta=0.0,
            dissipation=0.0,
            passive=True,
            detail=f"linear spring k={self.stiffness} L0={self.rest_length}",
        )
        self.last_pair = pair
        return pair


@dataclass(slots=True)
class ViewerDemoWorld:
    """The assembled producer, with every channel the parity fingerprint needs."""

    pipeline: Pipeline
    backend: NumpyBackend
    ledger: ForceWorkLedger
    rng: RngStream
    clock: Clock
    topology: Topology
    field_owner: RelaxingField
    anchor_owner: RelaxingField
    spring: AnchorSpring

    def transaction(self, **kwargs: Any) -> Transaction:
        """A transaction over this world's four channels."""
        return Transaction(
            self.pipeline,
            backend=self.backend,
            ledger=self.ledger,
            rng=self.rng,
            clock=self.clock,
            topology=self.topology,
            **kwargs,
        )

    def authoritative_arrays(self) -> dict[str, np.ndarray]:
        """Every authoritative array, labelled ``owner/array``, for parity fingerprinting."""
        return {
            f"{self.field_owner.name}/positions": self.field_owner.positions,
            f"{self.field_owner.name}/forces": self.field_owner.forces,
            f"{self.anchor_owner.name}/positions": self.anchor_owner.positions,
            f"{self.anchor_owner.name}/forces": self.anchor_owner.forces,
        }


def build_demo_world(*, n_nodes: int = 64, seed: int = 20260730) -> ViewerDemoWorld:
    """Assemble the demo producer.

    Args:
        n_nodes: Nodes on the relaxing field.
        seed: RNG seed. Nothing here draws from the stream; the seed is carried so the parity
            fingerprint has a stream position to compare, and an attachment that draws from it is
            caught.

    Returns:
        The world.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, int(n_nodes), endpoint=False)
    ring = np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    # Start perturbed and relax toward the clean ring, so the rendered field has visible structure
    # that decays monotonically. The perturbation is deterministic, not drawn.
    perturbation = 0.25 * np.stack(
        [np.sin(3.0 * theta), np.cos(5.0 * theta), 0.35 * np.sin(7.0 * theta)], axis=1
    )
    field_owner = RelaxingField(
        name="field",
        entity="field_nodes",
        positions=ring + perturbation,
        targets=ring,
    )
    anchor_owner = RelaxingField(
        name="anchor",
        entity="anchor_nodes",
        positions=np.array([[0.0, 0.0, -0.5], [0.0, 0.0, -1.5]]),
        targets=np.array([[0.0, 0.0, -0.5], [0.0, 0.0, -1.5]]),
    )
    spring = AnchorSpring(name="field_anchor_spring", endpoints=("field", "anchor"))

    pipeline = Pipeline()
    pipeline.register_participant(field_owner, phases=(Phase.ACCUMULATE_INTERNAL, Phase.SOLVE))
    pipeline.register_participant(anchor_owner, phases=(Phase.ACCUMULATE_INTERNAL, Phase.SOLVE))
    pipeline.register_connector(spring, phases=(Phase.ACCUMULATE_COUPLING,))
    pipeline.bind_participant("field", Phase.ACCUMULATE_INTERNAL)
    pipeline.bind_participant("anchor", Phase.ACCUMULATE_INTERNAL)
    pipeline.bind_participant("field", Phase.SOLVE)
    pipeline.bind_participant("anchor", Phase.SOLVE)

    def resolve(ctx: StepContext) -> tuple[EndpointHandle, EndpointHandle]:
        return (
            EndpointHandle(
                owner="field",
                entity="field_nodes",
                site_index=0,
                position=tuple(field_owner.positions[0]),
                force_sink=field_owner.forces,
            ),
            EndpointHandle(
                owner="anchor",
                entity="anchor_nodes",
                site_index=0,
                position=tuple(anchor_owner.positions[0]),
                force_sink=anchor_owner.forces,
            ),
        )

    pipeline.bind_connector("field_anchor_spring", resolve)

    return ViewerDemoWorld(
        pipeline=pipeline,
        backend=NumpyBackend(),
        ledger=ForceWorkLedger(),
        rng=RngStream(seed),
        clock=Clock(),
        topology=Topology(),
        field_owner=field_owner,
        anchor_owner=anchor_owner,
        spring=spring,
    )
