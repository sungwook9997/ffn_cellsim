"""Viewer off/on parity: the headline guarantee, measured rather than asserted.

Manuscript §11: *"Viewer off/on authoritative final state is bit-identical; D2H on the authoritative
hot loop remains zero."*

The first half is checkable here and this module checks it. The second half is a device property and
there is no device in this lane, so it is stated as an open item rather than implied by a passing test.

**Why a digest and not ``np.allclose``.** Bit-identical means bit-identical. Two runs that agree to
1e-15 have *not* satisfied this contract, because the failure mode being excluded is not "the viewer
perturbed the physics a lot" — it is "the viewer perturbed the physics at all", and the way that
usually happens is through a shared counter rather than through arithmetic. So the fingerprint hashes
raw array bytes, and it also covers four things that are not arrays:

* the backend operation counter, because the coverage predicate reads it as an evaluation witness — a
  viewer that ran one extra backend op could flip an acceptance decision without touching a single
  physical value;
* the RNG draw count and the generator's exact state, because a viewer that drew one variate moves
  every subsequent stochastic proposal;
* the clock, both physical time and accepted-step index;
* the topology generation.

A fingerprint that covered only the arrays would pass on a run where the viewer had consumed the random
stream, and the divergence would appear later, in a different run, as irreproducibility.

**What a failure means.** :func:`assert_bit_identical` names the first component that differs. An array
difference is a viewer that wrote through a handle it should not have had. A counter difference is a
viewer that did work on the producer's substrate. They are different bugs and the message says which.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any, Callable, Mapping

import numpy as np

__all__ = [
    "ParityResult",
    "ProducerFingerprint",
    "ParityViolation",
    "assert_bit_identical",
    "fingerprint_producer",
    "run_off_on_parity",
]


class ParityViolation(AssertionError):
    """The producer's state differs between the viewer-off and viewer-on runs."""


def _array_digest(array: np.ndarray) -> str:
    """Digest of an array's exact bytes, dtype and shape.

    ``tobytes`` on a C-ordered copy, so two arrays that differ only in stride bookkeeping still hash
    alike while any difference in a stored value does not.
    """
    hasher = hashlib.sha256()
    hasher.update(str(array.dtype).encode("utf-8"))
    hasher.update(repr(array.shape).encode("utf-8"))
    hasher.update(np.ascontiguousarray(array).tobytes())
    return hasher.hexdigest()


@dataclass(frozen=True, slots=True)
class ProducerFingerprint:
    """Everything about a producer that a correct viewer must leave untouched.

    Attributes:
        arrays: ``owner/array -> byte digest``.
        backend_op_count: The backend's evaluation witness.
        rng_draw_count: Variates drawn from the stream.
        rng_state_digest: Digest of the generator's exact internal state.
        clock_time: Physical time.
        clock_step_index: Accepted-step index.
        topology_generation: Topology generation counter.
    """

    arrays: Mapping[str, str]
    backend_op_count: int
    rng_draw_count: int
    rng_state_digest: str
    clock_time: float
    clock_step_index: int
    topology_generation: int

    def differences(self, other: ProducerFingerprint) -> tuple[str, ...]:
        """Every component that differs, most structural first.

        Ordered so the first line of a failure is the most informative one: a counter difference
        explains an array difference, never the other way round.
        """
        out: list[str] = []
        for name in (
            "backend_op_count",
            "rng_draw_count",
            "rng_state_digest",
            "clock_step_index",
            "clock_time",
            "topology_generation",
        ):
            mine, theirs = getattr(self, name), getattr(other, name)
            if mine != theirs:
                out.append(f"{name}: {mine!r} != {theirs!r}")
        missing = sorted(set(self.arrays) ^ set(other.arrays))
        if missing:
            out.append(f"array set differs: {missing}")
        for key in sorted(set(self.arrays) & set(other.arrays)):
            if self.arrays[key] != other.arrays[key]:
                out.append(
                    f"array {key}: {self.arrays[key][:16]}... != {other.arrays[key][:16]}..."
                )
        return tuple(out)

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able fingerprint."""
        return {
            "arrays": dict(self.arrays),
            "backend_op_count": self.backend_op_count,
            "rng_draw_count": self.rng_draw_count,
            "rng_state_digest": self.rng_state_digest,
            "clock_time": self.clock_time,
            "clock_step_index": self.clock_step_index,
            "topology_generation": self.topology_generation,
        }


def fingerprint_producer(
    *,
    arrays: Mapping[str, np.ndarray],
    backend: Any,
    rng: Any,
    clock: Any,
    topology: Any,
) -> ProducerFingerprint:
    """Capture a producer's complete observable state.

    Args:
        arrays: ``"owner/array" -> ndarray``, the authoritative arrays.
        backend: Anything with ``op_count``.
        rng: Anything with ``draw_count`` and ``snapshot()``.
        clock: Anything with ``time`` and ``step_index``.
        topology: Anything with ``generation``.

    Returns:
        The fingerprint.
    """
    snapshot = rng.snapshot()
    state_repr = repr(sorted(dict(snapshot.bit_generator_state).items(), key=lambda kv: str(kv[0])))
    return ProducerFingerprint(
        arrays={name: _array_digest(array) for name, array in sorted(arrays.items())},
        backend_op_count=int(backend.op_count),
        rng_draw_count=int(snapshot.draw_count),
        rng_state_digest=hashlib.sha256(state_repr.encode("utf-8")).hexdigest(),
        clock_time=float(clock.time),
        clock_step_index=int(clock.step_index),
        topology_generation=int(topology.generation),
    )


def assert_bit_identical(
    off: ProducerFingerprint, on: ProducerFingerprint, *, context: str = ""
) -> None:
    """Raise unless the two fingerprints are identical in every component.

    Args:
        off: Fingerprint of the run with no viewer attached.
        on: Fingerprint of the run with a viewer attached.
        context: What was being compared, for the message.

    Raises:
        ParityViolation: On the first difference, with every difference listed.
    """
    diffs = off.differences(on)
    if diffs:
        head = f"viewer off/on parity failed{(' — ' + context) if context else ''}"
        raise ParityViolation(
            head
            + "\n  "
            + "\n  ".join(diffs)
            + "\n\nA counter difference means the viewer did work on the producer's substrate — the "
            "backend op count is read by the coverage predicate, so this can change an acceptance "
            "decision without changing a single physical value. An array difference means the viewer "
            "held a writeable handle on authoritative memory."
        )


@dataclass(frozen=True, slots=True)
class ParityResult:
    """The outcome of an off/on comparison.

    Attributes:
        off: Fingerprint with no viewer.
        on: Fingerprint with a viewer.
        identical: Whether every component matched.
        differences: The differences, empty when identical.
        frames_published: How many frames the attached run published — reported so a "parity holds"
            result cannot come from a viewer that was attached and did nothing.
    """

    off: ProducerFingerprint
    on: ProducerFingerprint
    identical: bool
    differences: tuple[str, ...]
    frames_published: int

    def describe(self) -> str:
        """One-paragraph summary."""
        head = "IDENTICAL" if self.identical else "DIFFERENT"
        return (
            f"viewer off/on: {head}; {self.frames_published} frame(s) published on the attached run"
            + ("" if self.identical else "\n  " + "\n  ".join(self.differences))
        )


def run_off_on_parity(
    build: Callable[[], Any],
    run: Callable[[Any, Any], int],
    fingerprint: Callable[[Any], ProducerFingerprint],
    make_attachment: Callable[[Any], Any],
) -> ParityResult:
    """Run the same trajectory twice — once with no viewer, once with one — and compare.

    The two runs must be built from scratch by ``build`` rather than reset, because a reset that misses
    a field would make the comparison meaningless in the direction that hides a bug.

    Args:
        build: ``() -> world``. Must produce an identically-seeded world every call.
        run: ``(world, attachment_or_None) -> frames_published``.
        fingerprint: ``(world) -> ProducerFingerprint``.
        make_attachment: ``(world) -> attachment``.

    Returns:
        The comparison.
    """
    off_world = build()
    run(off_world, None)
    off = fingerprint(off_world)

    on_world = build()
    attachment = make_attachment(on_world)
    published = int(run(on_world, attachment))
    on = fingerprint(on_world)

    diffs = off.differences(on)
    return ParityResult(
        off=off,
        on=on,
        identical=not diffs,
        differences=diffs,
        frames_published=published,
    )
