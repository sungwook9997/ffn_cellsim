"""Frame transports: the CPU path built, the device path declared and refusing.

Manuscript §11's data path is

    authoritative state -> accepted-step publisher -> ObservationRegistry + derived buffers
    -> triple-buffered CUDA IPC -> independent native viewer -> recorder / replay / export

and §11 also says the packing and the viewer *can be prototyped before the whole-cell trajectory
closes*. That is what this module does: the transport is a Protocol, the triple-buffer semantics are
implemented and tested on the CPU, and the CUDA IPC implementation is a named seam that refuses with
its reason rather than a stub that returns something plausible.

The triple buffer is the interesting part and it is not an optimisation. With two buffers the producer
must wait for the consumer to release one, which means a slow viewer stalls physics — the exact thing
§11 forbids. With three, the producer always has a free slot: it writes into it, publishes it, and if
the consumer has not taken the previously published frame, that frame is *overwritten and counted as
dropped*. The producer never waits and never fails. Dropping is not a degradation of this design, it is
the design, which is why the count is mandatory on every accepted frame.

Two transports, two contracts:

* :class:`TripleBuffer` implements the LIVE contract. It drops silently-but-countably and never blocks.
* :class:`LosslessFrameQueue` implements the EVIDENCE_CAPTURE contract. It also never blocks, and when
  it cannot keep a frame it records a *capture failure* rather than dropping quietly. The distinction
  matters: a live viewer that misses a frame has missed a frame, and an evidence capture that misses a
  frame is no longer a capture at the cadence it declared. Physics continues in both cases.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Protocol, runtime_checkable

from aleph.viz.modes import ViewerMode

__all__ = [
    "CudaIpcTransport",
    "FrameTransport",
    "LosslessFrameQueue",
    "PublishOutcome",
    "TransportUnavailableError",
    "TripleBuffer",
    "transport_report",
]


class TransportUnavailableError(RuntimeError):
    """A transport was constructed that cannot run in this environment, and said so."""


@dataclass(frozen=True, slots=True)
class PublishOutcome:
    """What one publish attempt did.

    Attributes:
        published: Whether the frame is now visible to a consumer.
        dropped_now: How many frames this attempt discarded (0 or 1 for a triple buffer).
        dropped_total: Cumulative drops on this transport.
        detail: What happened, in words.
    """

    published: bool
    dropped_now: int
    dropped_total: int
    detail: str = ""


@runtime_checkable
class FrameTransport(Protocol):
    """A one-way channel from the producer to a viewer.

    Every method is non-blocking by contract. A transport that can block is a transport that can stall
    the authoritative loop, and §11's guarantee is that viewer behaviour never affects physics.
    """

    @property
    def name(self) -> str:
        """Stable transport identity, for the manifest."""
        ...

    @property
    def mode(self) -> ViewerMode:
        """Which contract this transport implements."""
        ...

    def publish(self, frame: object) -> PublishOutcome:
        """Offer a frame. Never blocks, never raises on a full channel."""
        ...

    def acquire(self) -> object | None:
        """Take the most recent available frame, or ``None``."""
        ...

    @property
    def dropped_total(self) -> int:
        """Cumulative frames this transport discarded."""
        ...


@dataclass(slots=True)
class TripleBuffer:
    """Three slots, a producer that never waits, and an honest drop count.

    The slot bookkeeping is deliberately explicit rather than a ring index: a reader can check by eye
    that the producer never writes the slot the consumer holds, which is the whole safety argument.

    Attributes:
        name: Transport identity.
        mode: Always ``LIVE``; refused otherwise. The lossy contract belongs to one mode only.
    """

    name: str = "triple_buffer"
    mode: ViewerMode = ViewerMode.LIVE
    _slots: list[object | None] = field(default_factory=lambda: [None, None, None], init=False)
    _published: int | None = field(default=None, init=False)
    _held_by_consumer: int | None = field(default=None, init=False)
    _dropped: int = field(default=0, init=False)
    _publish_count: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        if self.mode is not ViewerMode.LIVE:
            raise ValueError(
                f"a triple buffer drops frames by construction, so it implements the LIVE contract "
                f"only; got mode {self.mode!r}. Use LosslessFrameQueue for evidence capture."
            )

    @property
    def dropped_total(self) -> int:
        """Cumulative frames overwritten before a consumer took them."""
        return self._dropped

    @property
    def publish_count(self) -> int:
        """How many frames were offered."""
        return self._publish_count

    @property
    def pending(self) -> bool:
        """Whether a published frame is waiting to be taken."""
        return self._published is not None

    def _free_slot(self) -> int:
        for index in range(3):
            if index != self._published and index != self._held_by_consumer:
                return index
        raise AssertionError(
            "a triple buffer always has a free slot: at most one is published and at most one is "
            "held by the consumer. Reaching here means the slot bookkeeping is wrong."
        )

    def publish(self, frame: object) -> PublishOutcome:
        """Write ``frame`` into the free slot and make it the published one.

        The previously published frame, if the consumer never took it, is discarded and counted. The
        producer does not wait and cannot fail.
        """
        self._publish_count += 1
        target = self._free_slot()
        dropped_now = 0
        if self._published is not None:
            self._slots[self._published] = None
            self._dropped += 1
            dropped_now = 1
        self._slots[target] = frame
        self._published = target
        return PublishOutcome(
            published=True,
            dropped_now=dropped_now,
            dropped_total=self._dropped,
            detail=(
                "overwrote an unconsumed frame; the consumer is behind"
                if dropped_now
                else "published into a free slot"
            ),
        )

    def acquire(self) -> object | None:
        """Take the published frame, releasing whatever the consumer held before."""
        if self._held_by_consumer is not None:
            self._slots[self._held_by_consumer] = None
            self._held_by_consumer = None
        if self._published is None:
            return None
        self._held_by_consumer = self._published
        self._published = None
        return self._slots[self._held_by_consumer]

    def release(self) -> None:
        """Give back whatever the consumer holds. Idempotent."""
        if self._held_by_consumer is not None:
            self._slots[self._held_by_consumer] = None
            self._held_by_consumer = None


@dataclass(slots=True)
class LosslessFrameQueue:
    """A bounded queue that never drops silently and never blocks the producer.

    On overflow it does not discard the frame and does not raise into the producer's call stack.
    Instead it records a capture failure and reports ``published=False``. The capture that owns this
    queue then knows its declared cadence was not met and can say so in its manifest, which is the
    honest outcome: the physics ran, the evidence is incomplete, and the record says which.

    Attributes:
        name: Transport identity.
        capacity: Maximum frames held. ``0`` means unbounded, which is fine for a CPU capture writing
            to disk on the same thread.
        mode: ``EVIDENCE_CAPTURE`` or ``REPLAY``; ``LIVE`` is refused.
    """

    name: str = "lossless_frame_queue"
    capacity: int = 0
    mode: ViewerMode = ViewerMode.EVIDENCE_CAPTURE
    _frames: list[object] = field(default_factory=list, init=False)
    _taken: int = field(default=0, init=False)
    _failures: list[str] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        if self.mode is ViewerMode.LIVE:
            raise ValueError(
                "a lossless queue is for a contract that may not drop; LIVE may drop, so use "
                "TripleBuffer for it"
            )
        if self.capacity < 0:
            raise ValueError(f"capacity cannot be negative; got {self.capacity}")

    @property
    def dropped_total(self) -> int:
        """Always zero: this transport does not drop. Overflow becomes a failure instead."""
        return 0

    @property
    def capture_failures(self) -> tuple[str, ...]:
        """Every frame this queue could not accept, with the reason."""
        return tuple(self._failures)

    @property
    def depth(self) -> int:
        """Frames currently held."""
        return len(self._frames)

    def publish(self, frame: object) -> PublishOutcome:
        """Enqueue ``frame`` unless the queue is full.

        Returns an outcome with ``published=False`` on overflow. It does not raise, because raising
        here would propagate into the producer's step and §11 forbids the viewer affecting physics.
        """
        if self.capacity and len(self._frames) >= self.capacity:
            why = (
                f"queue is full at capacity {self.capacity}: the consumer is not draining, so the "
                "declared capture cadence was not met at this step"
            )
            self._failures.append(why)
            return PublishOutcome(
                published=False, dropped_now=0, dropped_total=0, detail=why
            )
        self._frames.append(frame)
        return PublishOutcome(
            published=True, dropped_now=0, dropped_total=0, detail="enqueued"
        )

    def acquire(self) -> object | None:
        """Take the oldest frame. Oldest, not newest: order is information in a capture."""
        if not self._frames:
            return None
        self._taken += 1
        return self._frames.pop(0)

    def drain(self) -> tuple[object, ...]:
        """Take every held frame, in order."""
        out = tuple(self._frames)
        self._taken += len(out)
        self._frames.clear()
        return out


@dataclass(slots=True)
class CudaIpcTransport:
    """The device-side transport §11 specifies. **Declared, not implemented.**

    This class exists so the seam is a named type with a signature rather than a sentence in a design
    document, and so a caller that reaches for it gets a refusal that says what is missing instead of a
    stub that returns something usable-looking.

    What is genuinely absent, as of this lane:

    * there are no device kernels to tap — no ``@wp.kernel`` exists in the package, so there is no
      authoritative device buffer for an IPC handle to point at;
    * ``warp`` is not installed in this environment;
    * no GPU authorization is on record, and none was created.

    Substituting the real implementation later is a transport swap and nothing more, because the
    :class:`FrameTransport` protocol is what the publisher and the capture are written against. The one
    thing the substitution must preserve is the direction of the IPC handle: the viewer receives a
    handle onto a *derived* buffer, read-only, never onto authoritative device memory. That is a
    property of the exporting side and cannot be checked here, so it is written into
    :attr:`SUBSTITUTION_CONTRACT` where whoever does the work will meet it.
    """

    SUBSTITUTION_CONTRACT: ClassVar[str] = (
        "the exported IPC handle must address a derived tap buffer, opened read-only on the viewer "
        "side; exporting a handle onto an authoritative array is a write path into physics even when "
        "the viewer intends only to read, and no test on the viewer side can detect it"
    )

    name: str = "cuda_ipc_triple_buffer"
    mode: ViewerMode = ViewerMode.LIVE

    def __post_init__(self) -> None:
        raise TransportUnavailableError(
            "CudaIpcTransport is a declared seam and is not implemented. There are no device "
            "kernels in aleph/ to tap, warp is not installed, and no GPU authorization is on "
            "record. The CPU path (TripleBuffer / LosslessFrameQueue) implements the same "
            "FrameTransport protocol, so this is a substitution and not a rewrite.\n"
            f"When it is written: {CudaIpcTransport.SUBSTITUTION_CONTRACT}."
        )

    def publish(self, frame: object) -> PublishOutcome:  # pragma: no cover - unreachable
        raise TransportUnavailableError("CudaIpcTransport is not implemented")

    def acquire(self) -> object | None:  # pragma: no cover - unreachable
        raise TransportUnavailableError("CudaIpcTransport is not implemented")

    @property
    def dropped_total(self) -> int:  # pragma: no cover - unreachable
        raise TransportUnavailableError("CudaIpcTransport is not implemented")


def transport_report(transport: Any) -> dict[str, Any]:
    """Return a manifest block describing ``transport``.

    Written as a function over the protocol rather than a method on each class, so adding a transport
    does not mean remembering to add a reporter.
    """
    return {
        "name": str(getattr(transport, "name", type(transport).__name__)),
        "type": type(transport).__name__,
        "mode": str(getattr(transport, "mode", ViewerMode.LIVE).value),
        "dropped_total": int(getattr(transport, "dropped_total", 0)),
        "capture_failures": list(getattr(transport, "capture_failures", ())),
    }
