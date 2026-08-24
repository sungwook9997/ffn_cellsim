"""The three frame types, which are three different kinds of claim.

Manuscript §11 requires that an accepted frame carry *accepted step, physical time, topology epoch,
active count, units, scalar channel, and dropped-frame count*, and that *rejected candidates can
appear only in a separately labeled diagnostic layer*.

Three classes, no shared base, on purpose:

* :class:`AcceptedFrame` is evidence. It cannot be constructed without every identifier §11 lists,
  including the dropped-frame count, and it cannot be constructed in a mode that is not
  ``EVIDENCE_CAPTURE``.
* :class:`LiveFrame` is not evidence and has no field that could be mistaken for one. It has no
  ``frame_id``, no census digest and no build stamp; it carries ``latest_accepted_step`` rather than
  ``accepted_step``, because a live frame is a picture of *whatever the producer had got to*, not a
  statement about a particular accepted state.
* :class:`DiagnosticFrame` carries a rejected candidate. It knows the candidate was rejected and why.

There is no common ancestor and no converter. That is the mechanism: :meth:`AcceptedFrameStream.append`
type-checks, and since ``LiveFrame`` is not a subclass of ``AcceptedFrame`` there is no
``isinstance``-shaped hole to fall through, no ``promote()`` to call by accident, and nothing a
copy-constructor could smuggle across. Filing a live frame as evidence requires re-deriving every
identifier it does not have, which is exactly the work that ought to be required.

A note on why the dropped count is mandatory rather than defaulted to zero. Defaulting it to zero
makes "no frames were dropped" and "nobody counted" the same value, and the second one is the state a
best-effort transport is in by default. A frame that does not know its drop count cannot support a
rate, and a rate is the most common thing anyone wants from a frame series.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any, Mapping

from aleph.artifacts.stamp import canonical_json
from aleph.viz.buffers import DerivedBuffer
from aleph.viz.channels import ColourBar
from aleph.viz.modes import ViewerMode

__all__ = [
    "AcceptedFrame",
    "AcceptedFrameStream",
    "DiagnosticFrame",
    "FrameIdentity",
    "LiveFrame",
    "NotEvidenceError",
    "frame_digest",
]


class NotEvidenceError(TypeError):
    """Something that is not an accepted frame was offered to the accepted stream."""


@dataclass(frozen=True, slots=True)
class FrameIdentity:
    """The build/config/census/topology identifiers §11 requires on an accepted frame.

    Attributes:
        build_commit: The commit the producing build was made on, or a declared commit. Not optional:
            a frame whose build is unknown cannot be reproduced, and an unreproducible frame is a
            picture.
        config_digest: Canonical hash of the run configuration.
        census_digest: Canonical hash of the population census. The census itself lives in the capture
            manifest; the frame carries its digest so a frame lifted out of its manifest still says
            which population it depicts.
        accepted_state_digest: The runtime-issued immutable digest of the accepted state this frame
            was derived from, when the producer computed one. ``None`` is permitted and is visible —
            it downgrades what the frame supports rather than being silently absent.
        schema_hash: Digest of the state schema.
    """

    build_commit: str
    config_digest: str
    census_digest: str
    schema_hash: str
    accepted_state_digest: str | None = None

    def __post_init__(self) -> None:
        for name in ("build_commit", "config_digest", "census_digest", "schema_hash"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"FrameIdentity.{name} is required; §11 lists it among the identifiers an "
                    "evidence frame carries, and a frame missing one cannot be traced to the state "
                    "it depicts"
                )

    @property
    def traces_to_accepted_state(self) -> bool:
        """Whether this frame names the accepted state it came from."""
        return bool(self.accepted_state_digest)

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able identity block."""
        return {
            "build_commit": self.build_commit,
            "config_digest": self.config_digest,
            "census_digest": self.census_digest,
            "schema_hash": self.schema_hash,
            "accepted_state_digest": self.accepted_state_digest,
            "traces_to_accepted_state": self.traces_to_accepted_state,
        }


def _buffer_payload_digest(buffers: Mapping[str, DerivedBuffer]) -> str:
    """Content digest over the derived values themselves, in name order.

    The values are what a reader would re-measure, so they are inside the frame's identity. Hashing
    the metadata alone would let two frames with identical labels and different pixels share an id.
    """
    hasher = hashlib.sha256()
    for name in sorted(buffers):
        buffer = buffers[name]
        hasher.update(name.encode("utf-8"))
        hasher.update(str(buffer.values.dtype).encode("utf-8"))
        hasher.update(repr(buffer.values.shape).encode("utf-8"))
        hasher.update(int(buffer.active_count).to_bytes(8, "big"))
        hasher.update(buffer.active_values.tobytes())
    return "sha256:" + hasher.hexdigest()


def frame_digest(payload: Mapping[str, Any]) -> str:
    """Return ``"sha256:<hex>"`` over the canonical JSON of ``payload``.

    Uses :func:`aleph.artifacts.stamp.canonical_json` rather than a local encoder so a frame id and a
    config digest are computed by the same rules. Two encoders is two conventions, and two conventions
    is a pair of digests that disagree about the same object.
    """
    return "sha256:" + hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class AcceptedFrame:
    """One published observation of one accepted step. This is the evidence object.

    Attributes:
        mode: Must be ``EVIDENCE_CAPTURE``. Present as a field rather than implied so a frame read
            back off disk says what contract produced it.
        identity: Build/config/census/schema identifiers.
        accepted_step: The accepted-step index, from the transaction's clock.
        attempted_step: The attempted-step index at the same instant. Kept separately because a state
            reached after rejections is not the state reached without them.
        physical_time_s: Simulated physical time [s].
        topology_epoch: Topology generation at this accepted step.
        active_count: Live population the frame depicts.
        colour_bar: The display mapping, carrying units, actual extrema and clamps.
        buffers: The derived render buffers, by name.
        dropped_frame_count: How many publishable frames were lost before this one. Required.
        cost_seconds: Wall time this frame's derivation and publication cost. §11 requires exact
            cost, and the cost of observing is part of the cost of the run.
        cadence_every_n_accepted_steps: The predeclared cadence this frame satisfies.
        events: Typed events at this step, for the scene tree's event column.
        notes: Anything else worth recording.
    """

    mode: ViewerMode
    identity: FrameIdentity
    accepted_step: int
    attempted_step: int
    physical_time_s: float
    topology_epoch: int
    active_count: int
    colour_bar: ColourBar
    buffers: Mapping[str, DerivedBuffer]
    dropped_frame_count: int | None
    cost_seconds: float
    cadence_every_n_accepted_steps: int
    events: tuple[str, ...] = ()
    notes: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.mode, ViewerMode):
            raise TypeError("AcceptedFrame.mode must be a ViewerMode")
        self.mode.require_publishing("building an AcceptedFrame")
        if not isinstance(self.identity, FrameIdentity):
            raise TypeError("AcceptedFrame.identity must be a FrameIdentity")
        if not isinstance(self.colour_bar, ColourBar):
            raise TypeError(
                "AcceptedFrame.colour_bar must be a ColourBar: a frame without units and clamps is a "
                "picture, and §11 forbids reading a quantity off a picture"
            )
        if self.dropped_frame_count is None:
            raise ValueError(
                "AcceptedFrame.dropped_frame_count is None. A frame that does not know how many "
                "publishable frames were lost before it cannot support a rate, and defaulting the "
                "count to zero would make 'none were dropped' and 'nobody counted' the same value."
            )
        if isinstance(self.dropped_frame_count, bool) or not isinstance(
            self.dropped_frame_count, int
        ):
            raise TypeError(
                f"dropped_frame_count must be an integer; got {self.dropped_frame_count!r}"
            )
        if self.dropped_frame_count < 0:
            raise ValueError(
                f"dropped_frame_count cannot be negative; got {self.dropped_frame_count}"
            )
        for name in ("accepted_step", "attempted_step", "topology_epoch", "active_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"AcceptedFrame.{name} must be an integer; got {value!r}")
            if value < 0:
                raise ValueError(f"AcceptedFrame.{name} cannot be negative; got {value}")
        if self.attempted_step < self.accepted_step:
            raise ValueError(
                f"attempted_step ({self.attempted_step}) is below accepted_step "
                f"({self.accepted_step}); a step is attempted before it can be accepted"
            )
        if not isinstance(self.physical_time_s, (int, float)) or isinstance(
            self.physical_time_s, bool
        ):
            raise TypeError(f"physical_time_s must be a real number; got {self.physical_time_s!r}")
        if float(self.physical_time_s) < 0.0:
            raise ValueError(f"physical_time_s cannot be negative; got {self.physical_time_s}")
        if float(self.cost_seconds) < 0.0:
            raise ValueError(
                f"cost_seconds cannot be negative; got {self.cost_seconds}. §11 requires exact "
                "cost, and a negative cost is a broken clock rather than a cheap frame"
            )
        if self.cadence_every_n_accepted_steps < 1:
            raise ValueError(
                "cadence_every_n_accepted_steps must be at least 1; a frame that satisfies no "
                "declared cadence was published opportunistically, which is the LIVE contract"
            )
        if not self.buffers:
            raise ValueError(
                "an accepted frame with no derived buffers depicts nothing; publish the reduction or "
                "do not publish the frame"
            )
        for name, buffer in self.buffers.items():
            if not isinstance(buffer, DerivedBuffer):
                raise TypeError(f"buffer {name!r} must be a DerivedBuffer")
            if buffer.accepted_step != self.accepted_step:
                raise ValueError(
                    f"buffer {name!r} was derived at accepted step {buffer.accepted_step} but the "
                    f"frame claims {self.accepted_step}. Mixing steps inside one frame renders two "
                    "different states as though they were simultaneous."
                )
            if buffer.topology_epoch != self.topology_epoch:
                raise ValueError(
                    f"buffer {name!r} carries topology epoch {buffer.topology_epoch}, frame claims "
                    f"{self.topology_epoch}; connectivity and values must be from the same epoch or "
                    "the render attaches values to elements that no longer exist"
                )
        if self.colour_bar.channel.name not in self.buffers:
            raise ValueError(
                f"the colour bar renders channel {self.colour_bar.channel.name!r}, which is not "
                f"among the frame's buffers {sorted(self.buffers)}"
            )

    # -- derived ------------------------------------------------------------------------------
    @property
    def units(self) -> str:
        """Units of the rendered channel."""
        return self.colour_bar.channel.units

    @property
    def scalar_channel(self) -> str:
        """Name of the rendered channel."""
        return self.colour_bar.channel.name

    @property
    def rate_is_computable(self) -> bool:
        """Whether a frame rate may be computed from this frame and its predecessor.

        True exactly when nothing was dropped before it. The count being *known* is what makes the
        question answerable at all; the count being zero is what makes the answer yes.
        """
        return self.dropped_frame_count == 0

    @property
    def supports_quantitative_claim(self) -> bool:
        """Whether a number may be read off this frame's channel.

        Requires an unclamped colour bar *and* a named accepted state. A clamped render cannot be
        inverted; an unnamed state cannot be retrieved.
        """
        return not self.colour_bar.is_clamped and self.identity.traces_to_accepted_state

    @property
    def frame_id(self) -> str:
        """Content-addressed frame identifier, derived rather than assigned.

        Derived so it cannot drift from the frame it names. An assigned id is a label somebody typed;
        this one changes if any depicted value changes.
        """
        return frame_digest(self._identity_payload())

    def _identity_payload(self) -> dict[str, Any]:
        return {
            "mode": self.mode.value,
            "identity": self.identity.as_json_obj(),
            "accepted_step": self.accepted_step,
            "attempted_step": self.attempted_step,
            "physical_time_s": float(self.physical_time_s),
            "topology_epoch": self.topology_epoch,
            "active_count": self.active_count,
            "colour_bar": self.colour_bar.as_json_obj(),
            "buffers": {name: self.buffers[name].as_json_obj() for name in sorted(self.buffers)},
            "payload_digest": _buffer_payload_digest(self.buffers),
            "dropped_frame_count": self.dropped_frame_count,
            "cadence_every_n_accepted_steps": self.cadence_every_n_accepted_steps,
            "events": list(self.events),
        }

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able frame header. Values are bulk data and are not inlined."""
        payload = self._identity_payload()
        payload["frame_id"] = frame_digest(payload)
        payload["cost_seconds"] = float(self.cost_seconds)
        payload["units"] = self.units
        payload["scalar_channel"] = self.scalar_channel
        payload["rate_is_computable"] = self.rate_is_computable
        payload["supports_quantitative_claim"] = self.supports_quantitative_claim
        if self.notes:
            payload["notes"] = dict(self.notes)
        return payload


@dataclass(frozen=True, slots=True)
class LiveFrame:
    """A best-effort observation. Not evidence, and structurally unable to pretend to be.

    It deliberately has no ``frame_id``, no :class:`FrameIdentity` and no cadence. What it has is
    ``latest_accepted_step``, whose name says the honest thing: this is a picture of whatever the
    producer had reached when the buffer was grabbed.

    Attributes:
        mode: Always ``LIVE``; refused otherwise, so an evidence-mode session cannot emit one of these
            and later be read as though it had published nothing.
        latest_accepted_step: The most recent accepted step the producer had reached.
        physical_time_s: Simulated time at that step.
        topology_epoch: Topology generation at that step.
        active_count: Live population.
        colour_bar: Display mapping. Required even here: a live viewer that shows an unlabelled field
            is where the habit of reading numbers off pictures starts.
        buffers: Derived buffers.
        dropped_since_previous: Best-effort drop count, or ``None`` when the transport does not know.
            Permitted to be ``None``, unlike on an accepted frame — that permission is precisely the
            difference between the two contracts.
    """

    mode: ViewerMode
    latest_accepted_step: int
    physical_time_s: float
    topology_epoch: int
    active_count: int
    colour_bar: ColourBar
    buffers: Mapping[str, DerivedBuffer]
    dropped_since_previous: int | None = None

    def __post_init__(self) -> None:
        if self.mode is not ViewerMode.LIVE:
            raise ValueError(
                f"a LiveFrame carries the LIVE contract; got mode {self.mode!r}. Publishing a "
                "best-effort frame under an evidence mode is how a dropped frame ends up inside an "
                "evidence series."
            )
        if not isinstance(self.colour_bar, ColourBar):
            raise TypeError("LiveFrame.colour_bar must be a ColourBar")
        if not self.buffers:
            raise ValueError("a live frame with no buffers depicts nothing")

    @property
    def is_evidence(self) -> bool:
        """Always ``False``. Present so the answer is queryable rather than remembered."""
        return False

    @property
    def not_evidence_reason(self) -> str:
        """Why this frame may not be cited."""
        return ViewerMode.LIVE.not_evidence_reason or ""

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able live-frame header, self-labelled as not evidence."""
        return {
            "mode": self.mode.value,
            "is_evidence": False,
            "not_evidence_reason": self.not_evidence_reason,
            "latest_accepted_step": self.latest_accepted_step,
            "physical_time_s": float(self.physical_time_s),
            "topology_epoch": self.topology_epoch,
            "active_count": self.active_count,
            "colour_bar": self.colour_bar.as_json_obj(),
            "buffers": {name: self.buffers[name].as_json_obj() for name in sorted(self.buffers)},
            "dropped_since_previous": self.dropped_since_previous,
        }


@dataclass(frozen=True, slots=True)
class DiagnosticFrame:
    """A rejected candidate, in the separately labelled layer §11 confines it to.

    Attributes:
        attempted_step: The attempted-step index of the rejected candidate.
        last_accepted_step: The most recent accepted step, i.e. what state is actually authoritative.
        physical_time_s: Physical time of the last accepted step, not of the candidate. A rejected
            candidate never had a time: the clock advances only on commit.
        topology_epoch: Topology generation of the candidate, which may differ from the accepted one
            and is exactly why a rejected candidate must not be rendered in the accepted layer.
        active_count: Population the candidate had.
        colour_bar: Display mapping.
        buffers: Derived buffers from the candidate.
        rejection_reasons: Which acceptance sub-conditions failed. Non-empty: a rejection with no
            stated failing condition is not a rejection anybody can act on.
        layer: Always ``"diagnostic"``.
    """

    attempted_step: int
    last_accepted_step: int
    physical_time_s: float
    topology_epoch: int
    active_count: int
    colour_bar: ColourBar
    buffers: Mapping[str, DerivedBuffer]
    rejection_reasons: tuple[str, ...]
    layer: str = "diagnostic"

    def __post_init__(self) -> None:
        if not self.rejection_reasons:
            raise ValueError(
                "a diagnostic frame must name the failing acceptance conditions; a rejected "
                "candidate with no stated reason is indistinguishable from an accepted one that was "
                "filed in the wrong layer"
            )
        if self.layer != "diagnostic":
            raise ValueError(
                f"a rejected candidate belongs in the 'diagnostic' layer, never {self.layer!r}; §11 "
                "confines rejected candidates to a separately labelled layer"
            )
        if not isinstance(self.colour_bar, ColourBar):
            raise TypeError("DiagnosticFrame.colour_bar must be a ColourBar")

    @property
    def is_evidence(self) -> bool:
        """Always ``False``: this depicts state that no transaction accepted."""
        return False

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able diagnostic-frame header."""
        return {
            "layer": self.layer,
            "is_evidence": False,
            "not_evidence_reason": (
                "this frame depicts a candidate the acceptance predicate rejected, so no "
                "authoritative state ever had these values"
            ),
            "attempted_step": self.attempted_step,
            "last_accepted_step": self.last_accepted_step,
            "physical_time_s": float(self.physical_time_s),
            "topology_epoch": self.topology_epoch,
            "active_count": self.active_count,
            "colour_bar": self.colour_bar.as_json_obj(),
            "buffers": {name: self.buffers[name].as_json_obj() for name in sorted(self.buffers)},
            "rejection_reasons": list(self.rejection_reasons),
        }


class AcceptedFrameStream:
    """The accepted stream. Nothing but an :class:`AcceptedFrame` gets in.

    Two invariants beyond the type check, both of which have caught real classes of mistake in other
    frame pipelines:

    * accepted steps are strictly increasing, so a replayed or duplicated frame cannot inflate a count;
    * a frame's cadence must match the stream's declared cadence, so a series cannot be assembled from
      frames published under two different sampling rules and then read as one rate.
    """

    def __init__(self, *, cadence_every_n_accepted_steps: int) -> None:
        if cadence_every_n_accepted_steps < 1:
            raise ValueError("a stream cadence must be at least 1 accepted step")
        self._cadence = int(cadence_every_n_accepted_steps)
        self._frames: list[AcceptedFrame] = []
        self._rejected_offers: list[str] = []

    @property
    def cadence_every_n_accepted_steps(self) -> int:
        """The declared cadence every member frame must satisfy."""
        return self._cadence

    @property
    def frames(self) -> tuple[AcceptedFrame, ...]:
        """The accepted frames, in publication order."""
        return tuple(self._frames)

    @property
    def rejected_offers(self) -> tuple[str, ...]:
        """Why each refused offer was refused. Kept so a silent drop is impossible."""
        return tuple(self._rejected_offers)

    @property
    def total_dropped(self) -> int:
        """Sum of the drop counts the member frames report."""
        return sum(int(frame.dropped_frame_count or 0) for frame in self._frames)

    def append(self, frame: object) -> AcceptedFrame:
        """Add a frame to the accepted stream.

        Args:
            frame: Must be an :class:`AcceptedFrame`.

        Returns:
            The frame.

        Raises:
            NotEvidenceError: If ``frame`` is a live frame, a diagnostic frame, or anything else. The
                message names what was offered, because the usual cause is a driver that publishes
                through one code path in both modes.
            ValueError: If the accepted step does not increase, or the cadence disagrees.
        """
        if not isinstance(frame, AcceptedFrame):
            offered = type(frame).__name__
            why = {
                "LiveFrame": (
                    "a LIVE frame is best-effort and does not know its drop count, its build, its "
                    "census or the accepted state it depicts; there is no conversion because there "
                    "is nothing to convert from"
                ),
                "DiagnosticFrame": (
                    "a diagnostic frame depicts a rejected candidate, which no accepted state ever "
                    "matched; it belongs in the diagnostic layer only"
                ),
            }.get(offered, "only an AcceptedFrame is evidence")
            self._rejected_offers.append(f"{offered}: {why}")
            raise NotEvidenceError(
                f"refusing to file a {offered} in the accepted stream — {why}"
            )
        if self._frames and frame.accepted_step <= self._frames[-1].accepted_step:
            why = (
                f"accepted step {frame.accepted_step} does not advance on "
                f"{self._frames[-1].accepted_step}; a repeated or reordered frame in an accepted "
                "series silently changes every rate computed from it"
            )
            self._rejected_offers.append(f"AcceptedFrame: {why}")
            raise ValueError(why)
        if frame.cadence_every_n_accepted_steps != self._cadence:
            why = (
                f"frame declares cadence {frame.cadence_every_n_accepted_steps}, stream declared "
                f"{self._cadence}; a series assembled from two sampling rules is not a series"
            )
            self._rejected_offers.append(f"AcceptedFrame: {why}")
            raise ValueError(why)
        self._frames.append(frame)
        return frame

    def offer(self, frame: object) -> bool:
        """Try to append, returning whether it was accepted instead of raising.

        For a driver that must keep stepping physics no matter what the viewer layer thinks.
        :meth:`append` has already recorded the reason in :attr:`rejected_offers` by the time this
        returns ``False``, so this is a non-raising path and not a silent one.
        """
        try:
            self.append(frame)
        except (NotEvidenceError, ValueError):
            return False
        return True
