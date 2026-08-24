"""Replay: chunked topology, frames, events and operator outputs, kept re-measurable.

Manuscript §11, replay row: *"Chunked topology, state frames, events, annotations, and operator outputs
remain re-measurable."*

"Re-measurable" is a stronger claim than "stored", and it is the one this module is built around. An
archive that stores an operator's *output* lets a reader see what the number was. An archive that is
re-measurable lets a reader recompute it from the stored frames and find out whether it was right. So
:meth:`ReplayArchive.remeasure` takes the operator, runs it over the stored frames, and compares. A
disagreement is a finding — either the stored output was computed from something the archive does not
contain, or the operator has changed since.

**Why topology is chunked rather than per-frame.** Connectivity changes at typed transitions and not at
every step, so storing it per-frame would be almost all redundancy — but the redundancy is not the
reason to chunk it. The reason is that a frame *must not* be storable without the connectivity it was
drawn against. :meth:`append_frame` refuses a frame whose topology epoch has no chunk, so an archive
cannot end up holding values attached to elements it cannot enumerate. That is the situation where a
replay silently draws the new mesh with the old field.

**Why replay may not mint.** :class:`ReplayArchive` accepts only frames that were already minted by an
:class:`~aleph.viz.capture.EvidenceCapture`, and ``ViewerMode.REPLAY.publishes_evidence`` is ``False``.
A replay that could create an accepted frame would be a way to manufacture accepted state that no
transaction ever accepted — the frames would look exactly like real ones, because they would be real
ones structurally, and nothing downstream could tell.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from aleph.viz.frames import AcceptedFrame
from aleph.viz.modes import ViewerMode

__all__ = [
    "EventRecord",
    "MissingTopologyChunkError",
    "OperatorOutput",
    "RemeasureResult",
    "ReplayArchive",
    "TopologyChunk",
    "archive_from_capture",
]


class MissingTopologyChunkError(RuntimeError):
    """A frame or event refers to a topology epoch the archive does not hold."""


@dataclass(frozen=True, slots=True)
class TopologyChunk:
    """The connectivity in force over a contiguous span of accepted steps.

    Attributes:
        epoch: The topology generation this chunk describes.
        first_accepted_step: First accepted step at which it was in force.
        connectivity: Named integer arrays — faces, edges, endpoint index pairs. Frozen on the way in.
        component_names: State-owning components present in this epoch.
        connector_names: Couplings present in this epoch.
        transition_reason: What changed to open this epoch. Required for every epoch after the first:
            an unexplained connectivity change is the thing a replay most needs to show.
    """

    epoch: int
    first_accepted_step: int
    connectivity: Mapping[str, np.ndarray]
    component_names: tuple[str, ...]
    connector_names: tuple[str, ...]
    transition_reason: str = ""

    def __post_init__(self) -> None:
        if self.epoch < 0:
            raise ValueError(f"topology epoch cannot be negative; got {self.epoch}")
        if self.first_accepted_step < 0:
            raise ValueError("first_accepted_step cannot be negative")
        if self.epoch > 0 and not self.transition_reason.strip():
            raise ValueError(
                f"topology epoch {self.epoch} gives no transition reason. Every epoch after the "
                "first opened because something changed structurally, and a replay that cannot say "
                "what changed will draw the new connectivity as though it had always been there."
            )
        frozen: dict[str, np.ndarray] = {}
        for name, array in self.connectivity.items():
            copied = np.array(array, copy=True)
            copied.setflags(write=False)
            frozen[str(name)] = copied
        object.__setattr__(self, "connectivity", frozen)

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able chunk header; arrays are described, not inlined."""
        return {
            "epoch": self.epoch,
            "first_accepted_step": self.first_accepted_step,
            "connectivity": {
                name: {"shape": list(array.shape), "dtype": str(array.dtype)}
                for name, array in sorted(self.connectivity.items())
            },
            "component_names": list(self.component_names),
            "connector_names": list(self.connector_names),
            "transition_reason": self.transition_reason,
        }


@dataclass(frozen=True, slots=True)
class EventRecord:
    """One typed event at one accepted step.

    Attributes:
        accepted_step: When it happened.
        kind: Event class, for example ``"bond_formed"`` or ``"contact_engaged"``.
        subject: Which node it happened to.
        detail: Free text.
        topology_epoch: The epoch in force. Recorded because an event that changes connectivity sits
            exactly at an epoch boundary and a replay has to know which side it belongs to.
    """

    accepted_step: int
    kind: str
    subject: str
    topology_epoch: int
    detail: str = ""

    def __post_init__(self) -> None:
        for name in ("kind", "subject"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"EventRecord.{name} is required")

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able event."""
        return {
            "accepted_step": self.accepted_step,
            "kind": self.kind,
            "subject": self.subject,
            "topology_epoch": self.topology_epoch,
            "detail": self.detail,
        }


@dataclass(frozen=True, slots=True)
class OperatorOutput:
    """One number an observation operator produced, with what it was computed over.

    Attributes:
        name: Operator identity.
        accepted_step: Which accepted step it read.
        channel: Which derived channel it read.
        value: The result.
        units: Units of the result.
    """

    name: str
    accepted_step: int
    channel: str
    value: float
    units: str

    def __post_init__(self) -> None:
        for name in ("name", "channel", "units"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"OperatorOutput.{name} is required")

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able operator output."""
        return {
            "name": self.name,
            "accepted_step": self.accepted_step,
            "channel": self.channel,
            "value": float(self.value),
            "units": self.units,
        }


@dataclass(frozen=True, slots=True)
class RemeasureResult:
    """What re-running an operator over the stored frames produced.

    Attributes:
        name: Operator identity.
        agrees: Whether every recomputed value matched its stored value exactly.
        comparisons: ``(accepted_step, stored, recomputed)`` for each output.
        worst_absolute_difference: Largest ``|stored - recomputed|``.
    """

    name: str
    agrees: bool
    comparisons: tuple[tuple[int, float, float], ...]
    worst_absolute_difference: float

    def describe(self) -> str:
        """One line."""
        verdict = "agrees" if self.agrees else "DISAGREES"
        return (
            f"operator {self.name!r} {verdict} on {len(self.comparisons)} stored output(s); "
            f"worst |stored - recomputed| = {self.worst_absolute_difference:.6e}"
        )


@dataclass(slots=True)
class ReplayArchive:
    """Chunked topology, accepted frames, events and operator outputs, kept together.

    Attributes:
        run_label: What capture this archive replays.
        mode: Always ``REPLAY``. Non-publishing by contract.
    """

    run_label: str
    mode: ViewerMode = ViewerMode.REPLAY
    _chunks: dict[int, TopologyChunk] = field(default_factory=dict, init=False)
    _frames: list[AcceptedFrame] = field(default_factory=list, init=False)
    _events: list[EventRecord] = field(default_factory=list, init=False)
    _outputs: list[OperatorOutput] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        if self.mode is not ViewerMode.REPLAY:
            raise ValueError(f"a ReplayArchive carries the REPLAY contract; got {self.mode!r}")

    # -- topology ------------------------------------------------------------------------------
    def append_topology_chunk(self, chunk: TopologyChunk) -> TopologyChunk:
        """Store a connectivity chunk.

        Raises:
            ValueError: If the epoch is already stored. Two chunks for one epoch means one of them is
                describing connectivity that was not in force, and there is no way to tell which.
        """
        if not isinstance(chunk, TopologyChunk):
            raise TypeError("append_topology_chunk takes a TopologyChunk")
        if chunk.epoch in self._chunks:
            raise ValueError(
                f"topology epoch {chunk.epoch} is already stored; an epoch is one connectivity, so "
                "a second chunk for it describes a state that never existed"
            )
        self._chunks[chunk.epoch] = chunk
        return chunk

    @property
    def topology_epochs(self) -> tuple[int, ...]:
        """Stored epochs, ascending."""
        return tuple(sorted(self._chunks))

    def topology_at(self, accepted_step: int) -> TopologyChunk:
        """The chunk in force at ``accepted_step``.

        Raises:
            MissingTopologyChunkError: When no chunk had opened by that step.
        """
        candidates = [
            chunk
            for chunk in self._chunks.values()
            if chunk.first_accepted_step <= int(accepted_step)
        ]
        if not candidates:
            raise MissingTopologyChunkError(
                f"no topology chunk had opened by accepted step {accepted_step}; the archive holds "
                f"epochs {list(self.topology_epochs)}"
            )
        return max(candidates, key=lambda chunk: chunk.first_accepted_step)

    # -- frames --------------------------------------------------------------------------------
    def append_frame(self, frame: AcceptedFrame) -> AcceptedFrame:
        """Store an already-minted accepted frame.

        Args:
            frame: A frame an :class:`~aleph.viz.capture.EvidenceCapture` minted. Replay does not mint;
                ``ViewerMode.REPLAY.publishes_evidence`` is ``False`` and this method has no path to
                constructing an :class:`AcceptedFrame`.

        Raises:
            TypeError: If it is not an :class:`AcceptedFrame`.
            MissingTopologyChunkError: If its epoch has no chunk. A frame stored without the
                connectivity it was drawn against will be replayed against whatever connectivity
                happens to be nearest, which draws this step's field on another step's mesh.
        """
        if not isinstance(frame, AcceptedFrame):
            raise TypeError(
                f"a replay archive stores accepted frames only; got {type(frame).__name__}. Replay "
                "re-reads evidence and does not create it."
            )
        if frame.topology_epoch not in self._chunks:
            raise MissingTopologyChunkError(
                f"frame at accepted step {frame.accepted_step} carries topology epoch "
                f"{frame.topology_epoch}, which has no stored chunk (have "
                f"{list(self.topology_epochs)}). Store the chunk first."
            )
        self._frames.append(frame)
        return frame

    @property
    def frames(self) -> tuple[AcceptedFrame, ...]:
        """Stored frames, in insertion order."""
        return tuple(self._frames)

    def frames_for_epoch(self, epoch: int) -> tuple[AcceptedFrame, ...]:
        """Frames drawn against one connectivity epoch."""
        return tuple(frame for frame in self._frames if frame.topology_epoch == int(epoch))

    def frame_at(self, accepted_step: int) -> AcceptedFrame:
        """The frame at ``accepted_step``.

        Raises:
            KeyError: When no frame was captured there — which is normal off-cadence, and the message
                says so rather than implying the step did not happen.
        """
        for frame in self._frames:
            if frame.accepted_step == int(accepted_step):
                return frame
        raise KeyError(
            f"no frame at accepted step {accepted_step}; the capture cadence samples "
            f"{[f.accepted_step for f in self._frames][:8]}... — an absent frame here means the step "
            "was not sampled, not that it did not occur"
        )

    # -- events and operator outputs ------------------------------------------------------------
    def append_event(self, event: EventRecord) -> EventRecord:
        """Store a typed event.

        Raises:
            MissingTopologyChunkError: If its epoch has no chunk.
        """
        if event.topology_epoch not in self._chunks:
            raise MissingTopologyChunkError(
                f"event {event.kind!r} at accepted step {event.accepted_step} names topology epoch "
                f"{event.topology_epoch}, which has no stored chunk"
            )
        self._events.append(event)
        return event

    @property
    def events(self) -> tuple[EventRecord, ...]:
        """Stored events."""
        return tuple(self._events)

    def append_operator_output(self, output: OperatorOutput) -> OperatorOutput:
        """Store one operator result.

        Raises:
            KeyError: If no stored frame carries the step and channel the output claims to have read.
                An output over a frame the archive does not hold is not re-measurable, and storing it
                anyway is how an archive comes to contain numbers nobody can check.
        """
        frame = self.frame_at(output.accepted_step)
        if output.channel not in frame.buffers:
            raise KeyError(
                f"operator {output.name!r} claims to have read channel {output.channel!r} at "
                f"accepted step {output.accepted_step}, but that frame carries "
                f"{sorted(frame.buffers)}. An output over data the archive does not hold cannot be "
                "re-measured."
            )
        self._outputs.append(output)
        return output

    @property
    def operator_outputs(self) -> tuple[OperatorOutput, ...]:
        """Stored operator outputs."""
        return tuple(self._outputs)

    # -- the re-measurement --------------------------------------------------------------------
    def remeasure(
        self,
        name: str,
        operator: Callable[[np.ndarray], float],
        *,
        tolerance: float = 0.0,
    ) -> RemeasureResult:
        """Re-run ``operator`` over the stored frames and compare against the stored outputs.

        Args:
            name: Which operator's outputs to check.
            operator: ``f(active_values) -> float``, applied to the channel the output names.
            tolerance: Allowed absolute difference. Defaults to ``0.0`` — exact. A stored output and a
                recomputation over the *same stored bytes* by the *same operator* have no reason to
                differ at all, so any non-zero tolerance here is hiding something.

        Returns:
            The comparison.

        Raises:
            KeyError: If no output carries that name.
        """
        outputs = [output for output in self._outputs if output.name == name]
        if not outputs:
            raise KeyError(
                f"no stored operator outputs named {name!r}; stored: "
                f"{sorted({o.name for o in self._outputs})}"
            )
        comparisons: list[tuple[int, float, float]] = []
        worst = 0.0
        for output in outputs:
            frame = self.frame_at(output.accepted_step)
            recomputed = float(operator(frame.buffers[output.channel].active_values))
            comparisons.append((output.accepted_step, float(output.value), recomputed))
            worst = max(worst, abs(float(output.value) - recomputed))
        return RemeasureResult(
            name=name,
            agrees=worst <= float(tolerance),
            comparisons=tuple(comparisons),
            worst_absolute_difference=worst,
        )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able archive index."""
        return {
            "schema": "aleph-viz-replay@1",
            "run_label": self.run_label,
            "mode": self.mode.value,
            "publishes_evidence": self.mode.publishes_evidence,
            "topology_chunks": [
                self._chunks[epoch].as_json_obj() for epoch in self.topology_epochs
            ],
            "frames": [frame.as_json_obj() for frame in self._frames],
            "events": [event.as_json_obj() for event in self._events],
            "operator_outputs": [output.as_json_obj() for output in self._outputs],
        }


def archive_from_capture(
    capture: Any,
    *,
    chunks: Sequence[TopologyChunk],
    events: Sequence[EventRecord] = (),
) -> ReplayArchive:
    """Build a replay archive from a finished capture.

    Chunks go in first, which is not an ordering convenience: :meth:`ReplayArchive.append_frame`
    refuses a frame whose epoch has no chunk, so a caller who assembles them the other way round gets
    a refusal rather than an archive with a hole in it.
    """
    archive = ReplayArchive(run_label=capture.run_label)
    for chunk in chunks:
        archive.append_topology_chunk(chunk)
    for frame in capture.frames:
        archive.append_frame(frame)
    for event in events:
        archive.append_event(event)
    return archive
