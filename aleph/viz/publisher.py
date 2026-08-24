"""The accepted-step publisher, and the failure isolation that makes it safe to attach.

Manuscript §11: *"UI/viewer failure never blocks physics"* and *"Viewer off/on authoritative final
state is bit-identical."*

Those are two different guarantees and they need two different mechanisms.

**Never blocks** is this module. Every call out to anything that is not the producer — a derivation, a
transport, a consumer callback — goes through :meth:`ViewerAttachment._isolated`, which catches
``Exception``, records a typed :class:`ViewerFault`, and returns. The producer's step continues. What is
deliberately *not* caught is ``BaseException`` that is not an ``Exception``: a ``KeyboardInterrupt`` or
a ``SystemExit`` is the operator asking the process to stop, and swallowing it in the name of protecting
physics would make the run unkillable.

**Bit-identical** is :mod:`aleph.viz.parity`. Isolation alone does not give it: an attachment that
succeeds perfectly can still perturb the producer by consuming its RNG stream or moving its backend
operation counter, and the operation counter is read by the coverage predicate, so a viewer that ran one
extra backend op could flip an acceptance decision. Nothing about a try/except detects that. The
attachment therefore also refuses to be given the producer's backend or RNG at all: derivations receive
arrays and nothing else, which is why :class:`~aleph.viz.buffers.DerivationSpec.reduce` takes a mapping
of arrays rather than a step context.

The publisher does not drive the transaction. It is called *after* a decision has been finalized, by a
driver, and it has no reference to the transaction, no way to propose or commit, and nothing that could
make it look like part of the step. :func:`step_and_publish` is that driver, and it is thirty lines so
the ordering is checkable by eye: physics first, decision second, observation last.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

from aleph.viz.buffers import DerivedBuffer, ObservationRegistry
from aleph.viz.channels import ColourBar, ScaleKind, colour_bar_for
from aleph.viz.frames import AcceptedFrame, DiagnosticFrame, LiveFrame
from aleph.viz.modes import ViewerMode
from aleph.viz.transport import FrameTransport, PublishOutcome

__all__ = [
    "PublishReport",
    "PublishingDriver",
    "ViewerAttachment",
    "ViewerFault",
    "step_and_publish",
]


@dataclass(frozen=True, slots=True)
class ViewerFault:
    """One caught viewer-side failure.

    Attributes:
        where: Which stage failed — ``"derive"``, ``"colour_bar"``, ``"frame"``, ``"transport"``,
            ``"consumer"``.
        step: The step index being observed when it failed.
        exception_type: Class name of the exception.
        message: ``str(exc)``.
        detail: Which named thing failed, when the stage has more than one.
    """

    where: str
    step: int
    exception_type: str
    message: str
    detail: str = ""

    def render(self) -> str:
        """One line, for a log or a manifest."""
        tail = f" [{self.detail}]" if self.detail else ""
        return f"step {self.step}: {self.where}{tail} raised {self.exception_type}: {self.message}"


@dataclass(frozen=True, slots=True)
class PublishReport:
    """What one observation attempt produced. Never raises; always says what happened.

    Attributes:
        published: Whether a frame reached the transport.
        frame: The frame, when one was built.
        outcome: The transport's own outcome, when it was reached.
        faults: Faults caught during this attempt.
        skipped_reason: Why nothing was published, when that was a decision rather than a failure —
            off-cadence, for instance, which is not an error.
    """

    published: bool
    frame: Any = None
    outcome: PublishOutcome | None = None
    faults: tuple[ViewerFault, ...] = ()
    skipped_reason: str = ""


class ViewerAttachment:
    """Everything the producer holds on the viewer's behalf, and nothing more.

    It has no transaction, no backend, no RNG stream and no clock. It has a registry (which holds the
    authoritative arrays for derivation and hands out nothing), a transport, a mode, and a list of
    consumer callbacks that are treated as hostile code.
    """

    def __init__(
        self,
        *,
        registry: ObservationRegistry,
        transport: FrameTransport,
        mode: ViewerMode,
        channel: str,
        cadence_every_n_accepted_steps: int = 1,
        display_range: tuple[float, float] | None = None,
        scale: ScaleKind = ScaleKind.LINEAR,
        palette: str = "viridis",
        consumers: Sequence[Callable[[object], None]] = (),
        clock: Callable[[], float] = time.perf_counter,
        name: str = "viewer",
    ) -> None:
        """Attach a viewer.

        Args:
            registry: The observation registry.
            transport: Where frames go.
            mode: Which §11 contract this attachment operates under.
            channel: Which registered derivation is the rendered scalar channel.
            cadence_every_n_accepted_steps: Publish every N accepted steps.
            display_range: ``(display_min, display_max)``, or ``None`` for no clamping. Supplying a
                range narrower than the field is how a clamp gets into a frame, and the frame will
                carry it.
            scale: Colour-bar scale.
            palette: Colour map name.
            consumers: Viewer-side callbacks. Called with the frame. Assumed able to raise.
            clock: Monotonic seconds source, injectable so cost is deterministic in a test.
            name: Attachment identity, for faults and manifests.
        """
        if not isinstance(registry, ObservationRegistry):
            raise TypeError("ViewerAttachment needs an ObservationRegistry")
        if not isinstance(mode, ViewerMode):
            raise TypeError("ViewerAttachment.mode must be a ViewerMode")
        if cadence_every_n_accepted_steps < 1:
            raise ValueError("cadence must be at least 1 accepted step")
        if channel not in registry.derivation_names:
            raise ValueError(
                f"channel {channel!r} is not a registered derivation; registered: "
                f"{list(registry.derivation_names)}"
            )
        if mode.requires_predeclared_cadence and cadence_every_n_accepted_steps < 1:
            raise ValueError("evidence capture needs a predeclared cadence")
        self._registry = registry
        self._transport = transport
        self._mode = mode
        self._channel = channel
        self._cadence = int(cadence_every_n_accepted_steps)
        self._display_range = display_range
        self._scale = scale
        self._palette = palette
        self._consumers = list(consumers)
        self._clock = clock
        self.name = name

        self._faults: list[ViewerFault] = []
        self._published_count = 0
        self._skipped_off_cadence = 0
        self._cost_seconds = 0.0

    # -- introspection ------------------------------------------------------------------------
    @property
    def mode(self) -> ViewerMode:
        """The §11 contract this attachment operates under."""
        return self._mode

    @property
    def cadence_every_n_accepted_steps(self) -> int:
        """Declared cadence."""
        return self._cadence

    @property
    def faults(self) -> tuple[ViewerFault, ...]:
        """Every viewer-side failure caught so far, in order."""
        return tuple(self._faults)

    @property
    def fault_count(self) -> int:
        """How many viewer-side failures were caught and contained."""
        return len(self._faults)

    @property
    def published_count(self) -> int:
        """Frames that reached the transport."""
        return self._published_count

    @property
    def dropped_total(self) -> int:
        """Frames the transport discarded."""
        try:
            return int(self._transport.dropped_total)
        except Exception:  # noqa: BLE001 - a transport that cannot report is not a physics problem
            return 0

    @property
    def cost_seconds(self) -> float:
        """Total wall time spent observing. §11 requires exact cost, and this is Aleph's share."""
        return self._cost_seconds

    def wants(self, accepted_step: int) -> bool:
        """Whether the declared cadence publishes at this accepted step."""
        return int(accepted_step) % self._cadence == 0

    # -- isolation ---------------------------------------------------------------------------
    def _isolated(
        self, where: str, step: int, fn: Callable[..., Any], *args: Any, detail: str = "", **kwargs: Any
    ) -> tuple[bool, Any]:
        """Call ``fn``, containing any ``Exception`` it raises as a recorded fault.

        Returns:
            ``(ok, result)``. On failure ``ok`` is ``False`` and ``result`` is ``None``.

        Note what is *not* caught. ``KeyboardInterrupt`` and ``SystemExit`` derive from
        ``BaseException`` and not from ``Exception``, so they propagate: an operator stopping the run
        must not be defeated by the viewer's crash barrier.
        """
        try:
            return (True, fn(*args, **kwargs))
        except Exception as exc:  # noqa: BLE001 - containing viewer failure is this method's job
            self._faults.append(
                ViewerFault(
                    where=where,
                    step=int(step),
                    exception_type=type(exc).__name__,
                    message=str(exc),
                    detail=detail,
                )
            )
            return (False, None)

    def _notify(self, frame: object, step: int) -> tuple[ViewerFault, ...]:
        """Hand the frame to every consumer, containing each one's failure separately.

        Separately, so one broken consumer does not deprive the others of the frame. A single
        try/except around the loop would make consumer 1's exception silently skip consumer 2, which is
        a viewer failure affecting a viewer — still not physics, but still a bug that hides.
        """
        before = len(self._faults)
        for index, consumer in enumerate(self._consumers):
            self._isolated(
                "consumer", step, consumer, frame, detail=f"consumer[{index}]"
            )
        return tuple(self._faults[before:])

    # -- the observations ---------------------------------------------------------------------
    def _derive(
        self, *, step: int, topology_epoch: int
    ) -> tuple[Mapping[str, DerivedBuffer] | None, ColourBar | None, tuple[ViewerFault, ...]]:
        before = len(self._faults)
        ok, buffers = self._isolated(
            "derive",
            step,
            self._registry.derive_all,
            accepted_step=step,
            topology_epoch=topology_epoch,
        )
        if not ok or not buffers:
            return (None, None, tuple(self._faults[before:]))
        lo, hi = (None, None) if self._display_range is None else self._display_range
        ok, bar = self._isolated(
            "colour_bar",
            step,
            colour_bar_for,
            buffers[self._channel],
            display_min=lo,
            display_max=hi,
            scale=self._scale,
            palette=self._palette,
            detail=self._channel,
        )
        if not ok:
            return (buffers, None, tuple(self._faults[before:]))
        return (buffers, bar, tuple(self._faults[before:]))

    def on_accepted_step(
        self,
        *,
        accepted_step: int,
        attempted_step: int,
        physical_time_s: float,
        topology_epoch: int,
        capture: Any = None,
        events: Sequence[str] = (),
    ) -> PublishReport:
        """Observe one accepted step. Never raises.

        In ``EVIDENCE_CAPTURE`` mode the frame is built by ``capture`` — an
        :class:`~aleph.viz.capture.EvidenceCapture` — because only that object holds the predeclared
        cadence, the identifiers and the census. In ``LIVE`` mode a :class:`LiveFrame` is built here.

        Args:
            accepted_step: Accepted-step index after the commit.
            attempted_step: Attempted-step index at the same instant.
            physical_time_s: Simulated time after the commit.
            topology_epoch: Topology generation after the commit.
            capture: Required in ``EVIDENCE_CAPTURE`` mode.
            events: Typed events at this step.

        Returns:
            A :class:`PublishReport`. A viewer failure shows up as a fault in the report, never as an
            exception in the caller.
        """
        started = self._clock()
        try:
            if not self.wants(accepted_step):
                self._skipped_off_cadence += 1
                return PublishReport(
                    published=False,
                    skipped_reason=(
                        f"accepted step {accepted_step} is off the declared cadence of "
                        f"{self._cadence}"
                    ),
                )
            buffers, bar, faults = self._derive(
                step=accepted_step, topology_epoch=topology_epoch
            )
            if buffers is None or bar is None:
                return PublishReport(published=False, faults=faults)

            active = int(buffers[self._channel].active_count)
            if self._mode is ViewerMode.LIVE:
                ok, frame = self._isolated(
                    "frame",
                    accepted_step,
                    LiveFrame,
                    mode=ViewerMode.LIVE,
                    latest_accepted_step=int(accepted_step),
                    physical_time_s=float(physical_time_s),
                    topology_epoch=int(topology_epoch),
                    active_count=active,
                    colour_bar=bar,
                    buffers=buffers,
                    dropped_since_previous=self.dropped_total,
                )
            else:
                if capture is None:
                    self._faults.append(
                        ViewerFault(
                            where="frame",
                            step=int(accepted_step),
                            exception_type="ValueError",
                            message=(
                                f"mode {self._mode.value} needs an EvidenceCapture to mint a frame; "
                                "the identifiers, the census and the predeclared cadence live there "
                                "and cannot be defaulted here"
                            ),
                        )
                    )
                    return PublishReport(published=False, faults=tuple(self._faults[-1:]))
                ok, frame = self._isolated(
                    "frame",
                    accepted_step,
                    capture.mint,
                    accepted=True,
                    accepted_step=int(accepted_step),
                    attempted_step=int(attempted_step),
                    physical_time_s=float(physical_time_s),
                    topology_epoch=int(topology_epoch),
                    active_count=active,
                    colour_bar=bar,
                    buffers=buffers,
                    dropped_frame_count=self.dropped_total,
                    events=tuple(events),
                    cost_seconds=max(0.0, self._clock() - started),
                )
            if not ok:
                return PublishReport(published=False, faults=tuple(self._faults[-1:]))

            ok, outcome = self._isolated("transport", accepted_step, self._transport.publish, frame)
            if not ok:
                return PublishReport(published=False, frame=frame, faults=tuple(self._faults[-1:]))
            if outcome.published:
                self._published_count += 1
            consumer_faults = self._notify(frame, accepted_step)
            return PublishReport(
                published=bool(outcome.published),
                frame=frame,
                outcome=outcome,
                faults=faults + consumer_faults,
            )
        finally:
            self._cost_seconds += max(0.0, self._clock() - started)

    def on_rejected_candidate(
        self,
        *,
        attempted_step: int,
        last_accepted_step: int,
        physical_time_s: float,
        topology_epoch: int,
        rejection_reasons: Sequence[str],
        diagnostic_transport: FrameTransport | None = None,
    ) -> PublishReport:
        """Observe a rejected candidate into the diagnostic layer only. Never raises.

        A rejected candidate is *not* offered to the attachment's own transport unless that transport
        is explicitly passed as ``diagnostic_transport``. §11 confines rejected candidates to a
        separately labelled layer, and the cleanest reading of "separately" is a separate channel: a
        consumer subscribed to the accepted stream never sees one of these at all.
        """
        started = self._clock()
        try:
            reasons = tuple(str(reason) for reason in rejection_reasons if str(reason).strip())
            if not reasons:
                self._faults.append(
                    ViewerFault(
                        where="frame",
                        step=int(attempted_step),
                        exception_type="ValueError",
                        message=(
                            "a rejected candidate was offered with no failing condition named; "
                            "without one it cannot be distinguished from an accepted step filed in "
                            "the wrong layer"
                        ),
                    )
                )
                return PublishReport(published=False, faults=tuple(self._faults[-1:]))
            buffers, bar, faults = self._derive(
                step=last_accepted_step, topology_epoch=topology_epoch
            )
            if buffers is None or bar is None:
                return PublishReport(published=False, faults=faults)
            ok, frame = self._isolated(
                "frame",
                attempted_step,
                DiagnosticFrame,
                attempted_step=int(attempted_step),
                last_accepted_step=int(last_accepted_step),
                physical_time_s=float(physical_time_s),
                topology_epoch=int(topology_epoch),
                active_count=int(buffers[self._channel].active_count),
                colour_bar=bar,
                buffers=buffers,
                rejection_reasons=reasons,
            )
            if not ok:
                return PublishReport(published=False, faults=tuple(self._faults[-1:]))
            if diagnostic_transport is None:
                return PublishReport(
                    published=False,
                    frame=frame,
                    faults=faults,
                    skipped_reason=(
                        "no diagnostic transport is attached, so the rejected candidate was built "
                        "and not published; it is never offered to the accepted stream"
                    ),
                )
            ok, outcome = self._isolated(
                "transport", attempted_step, diagnostic_transport.publish, frame
            )
            return PublishReport(
                published=bool(ok and outcome and outcome.published),
                frame=frame,
                outcome=outcome if ok else None,
                faults=faults,
            )
        finally:
            self._cost_seconds += max(0.0, self._clock() - started)

    def report(self) -> dict[str, Any]:
        """Return the attachment's own manifest block."""
        return {
            "name": self.name,
            "mode": self._mode.value,
            "channel": self._channel,
            "cadence_every_n_accepted_steps": self._cadence,
            "published_count": self._published_count,
            "skipped_off_cadence": self._skipped_off_cadence,
            "dropped_total": self.dropped_total,
            "fault_count": self.fault_count,
            "faults": [fault.render() for fault in self._faults],
            "cost_seconds": self._cost_seconds,
            "derivations": list(self._registry.derivation_names),
        }


def step_and_publish(
    transaction: Any,
    dt: float,
    predicate: Any,
    attachment: ViewerAttachment | None,
    *,
    attempted_step: int,
    capture: Any = None,
    diagnostic_transport: FrameTransport | None = None,
    events: Sequence[str] = (),
) -> Any:
    """Drive one candidate step and then observe it. Physics first, observation last.

    Deliberately short, because the ordering is the safety property and a reader has to be able to
    check it without following a call graph:

    1. the transaction proposes, decides and finalizes — the viewer is not consulted and cannot
       intervene;
    2. only then is the attachment called, with values read *out of* the transaction;
    3. the attachment's return value is discarded here, because nothing about the step depends on it.

    ``attachment=None`` is the viewer-off path, and it is the same function. That matters for the
    off/on parity test: if the off path were a different function, the test would be comparing two
    code paths rather than the presence of a viewer.

    ``attempted_step`` is a parameter rather than something read off the transaction because the
    transaction's clock counts *accepted* steps only — ``Clock.advance`` runs on commit. Nothing in the
    runtime tracks attempts, so a driver has to, and inferring it here from the accepted count would
    silently report a run with rejections as a run without any. :class:`PublishingDriver` keeps the
    counter.

    Args:
        transaction: An :class:`aleph.runtime.transaction.Transaction`.
        dt: Step size.
        predicate: The acceptance predicate.
        attachment: The viewer, or ``None``.
        attempted_step: Attempted-step count *including* the step about to be taken.
        capture: The :class:`~aleph.viz.capture.EvidenceCapture`, in evidence mode.
        diagnostic_transport: Where rejected candidates go, if anywhere.
        events: Typed events to attach to this step's frame.

    Returns:
        The transaction's :class:`~aleph.runtime.transaction.Decision`, unchanged.
    """
    decision = transaction.step(dt, predicate)
    if attachment is None:
        return decision
    if decision.accepted:
        attachment.on_accepted_step(
            accepted_step=int(transaction.clock.step_index),
            attempted_step=int(attempted_step),
            physical_time_s=float(transaction.clock.time),
            topology_epoch=int(transaction.topology.generation),
            capture=capture,
            events=events,
        )
    else:
        attachment.on_rejected_candidate(
            attempted_step=int(attempted_step),
            last_accepted_step=int(transaction.clock.step_index),
            physical_time_s=float(transaction.clock.time),
            topology_epoch=int(transaction.topology.generation),
            rejection_reasons=tuple(
                str(condition) for condition in decision.verdict.failed_conditions
            )
            or (
                "the acceptance predicate evaluated no condition, so there was nothing to accept on",
            ),
            diagnostic_transport=diagnostic_transport,
        )
    return decision


class PublishingDriver:
    """Steps a transaction, counts attempts, and observes accepted steps through an attachment.

    The counter is here rather than in the runtime because attempts are a property of *this loop*, not
    of the transaction: two drivers sharing one transaction would each be attempting steps and neither
    would own the total. Keeping it local is also why the viewer-off path is exact — with
    ``attachment=None`` this class does nothing but call ``transaction.step``.
    """

    def __init__(
        self,
        transaction: Any,
        *,
        attachment: ViewerAttachment | None = None,
        capture: Any = None,
        diagnostic_transport: FrameTransport | None = None,
    ) -> None:
        self.transaction = transaction
        self.attachment = attachment
        self.capture = capture
        self.diagnostic_transport = diagnostic_transport
        self._attempted = 0
        self._accepted = 0
        self._rejected = 0

    @property
    def attempted_steps(self) -> int:
        """Candidate steps proposed."""
        return self._attempted

    @property
    def accepted_steps(self) -> int:
        """Candidate steps accepted."""
        return self._accepted

    @property
    def rejected_steps(self) -> int:
        """Candidate steps rejected."""
        return self._rejected

    def step(self, dt: float, predicate: Any, *, events: Sequence[str] = ()) -> Any:
        """Take one candidate step, then observe it if a viewer is attached."""
        self._attempted += 1
        decision = step_and_publish(
            self.transaction,
            dt,
            predicate,
            self.attachment,
            attempted_step=self._attempted,
            capture=self.capture,
            diagnostic_transport=self.diagnostic_transport,
            events=events,
        )
        if decision.accepted:
            self._accepted += 1
        else:
            self._rejected += 1
        return decision

    def run(self, dt: float, predicate: Any, n_steps: int) -> tuple[Any, ...]:
        """Take ``n_steps`` candidate steps."""
        return tuple(self.step(dt, predicate) for _ in range(int(n_steps)))
