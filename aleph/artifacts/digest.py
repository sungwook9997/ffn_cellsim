"""A content digest for accepted state, built so that omitting a component is impossible by accident.

The reference project shipped a state digest that was the SHA-256 of whatever its serializer happened
to emit, and the serializer omitted two fields.  An external review found the consequence on
2026-07-29: two states with different provenance graphs and different attempted step indices produced
**identical digests**.  A content address that cannot see a difference is not a content address, and
the failure was structural rather than a typo — when the digest is "hash of the dict", every future
field is opt-in to the identity and the default is invisibility.

Aleph inverts the default.  :data:`DIGEST_COVERAGE` names the components the digest must cover, and
:func:`accepted_state_digest` takes every one of them as a **required keyword argument with no
default**.  A caller cannot omit one.  A future component cannot become quietly optional, because
adding it to the signature without adding it to :data:`DIGEST_COVERAGE` fails
``test_signature_matches_declared_coverage``.

The encoding is domain-separated and length-prefixed at every level.  Each component is hashed inside
a labelled frame carrying its own byte length, so no two differently-structured states can produce
the same byte stream by concatenation — the classic collision where owner names ``("ab", "c")`` and
``("a", "bc")`` flatten to ``abc``.
"""

from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "DIGEST_COVERAGE",
    "DIGEST_SCHEMA",
    "OwnerStateBlock",
    "RngState",
    "SimulationClock",
    "accepted_state_digest",
]

#: Identifier of this digest construction.  It is hashed into every digest, so changing the encoding
#: changes every digest — which is correct: a digest computed by a different rule is a different
#: identity, and letting the two share a namespace would be silent collision.
DIGEST_SCHEMA: str = "aleph-accepted-state-digest@1"

#: Every component of accepted state that this digest covers, in the order it is absorbed.
#:
#: This tuple is the contract, and it is checked against the real signature by a test.  Adding a
#: component to accepted state means adding it here *and* to the signature in the same change;
#: forgetting either half is what produced the reference project's colliding digests.
DIGEST_COVERAGE: tuple[str, ...] = (
    "schema_hash",
    "context_hash",
    "topology_generation",
    "owner_blocks",
    "rng_state",
    "clock",
)


def _frame(label: str, payload: bytes) -> bytes:
    """Wrap ``payload`` in a self-delimiting, labelled frame.

    The frame is ``len(label) || label || len(payload) || payload``, with both lengths as unsigned
    64-bit little-endian integers.  Length-prefixing is what makes concatenation unambiguous: without
    it, two adjacent fields can borrow bytes from each other and two different states can flatten to
    the same stream.

    Args:
        label: A short ASCII name for the field, so that the same bytes under a different name hash
            differently.
        payload: The field's bytes.

    Returns:
        The framed bytes.
    """
    encoded_label = label.encode("utf-8")
    return (
        struct.pack("<Q", len(encoded_label))
        + encoded_label
        + struct.pack("<Q", len(payload))
        + payload
    )


def _text(value: str) -> bytes:
    """Encode a string as UTF-8 for framing."""
    return value.encode("utf-8")


def _int(value: int) -> bytes:
    """Encode an integer of any width, sign included, without a fixed-size overflow.

    Fixed-width packing would silently wrap a step index past 2**64, so the value is written as its
    signed big-endian minimal byte form with an explicit length via the surrounding frame.

    Args:
        value: The integer.

    Returns:
        The encoded bytes.

    Raises:
        TypeError: If ``value`` is not an integer, or is a bool.  ``True`` and ``1`` must not hash
            alike, since a flag and a count are different state.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"expected an integer, got {type(value).__name__} ({value!r})")
    width = max(1, (value.bit_length() + 8) // 8)
    return value.to_bytes(width, "big", signed=True)


def _float(value: float) -> bytes:
    """Encode a float by its exact IEEE-754 double bits.

    Bit-exact rather than decimal, because the clock is state: two clocks one ULP apart are two
    different states and the digest must be able to say so.  ``-0.0`` and ``0.0`` therefore hash
    differently, which is intended.

    Args:
        value: The float.

    Returns:
        Eight bytes, little-endian.

    Raises:
        ValueError: If the value is NaN.  NaN is not equal to itself, so a state containing one has
            no stable identity and must be rejected before it acquires a digest.
    """
    as_float = float(value)
    if as_float != as_float:
        raise ValueError("NaN has no stable identity and cannot enter an accepted-state digest")
    return struct.pack("<d", as_float)


def _array_bytes(name: str, array: Any) -> bytes:
    """Frame one named state array, covering its dtype and shape as well as its contents.

    Two arrays with identical bytes but different shapes are different state, so the shape is hashed
    alongside the buffer rather than left implicit in the byte count.

    Args:
        name: The array's name within its owner block.
        array: A numpy array, a numpy scalar, or a nested sequence of numbers.

    Returns:
        The framed bytes.

    Raises:
        TypeError: If the object is neither array-like nor a sequence of numbers.
    """
    dtype = getattr(array, "dtype", None)
    shape = getattr(array, "shape", None)
    tobytes = getattr(array, "tobytes", None)
    if dtype is not None and shape is not None and callable(tobytes):
        # numpy path. `dtype.str` carries byte order, kind, and item size, so a float32 buffer and a
        # float64 buffer of the same length can never be confused.
        import numpy as np

        contiguous = np.ascontiguousarray(array)
        payload = (
            _frame("dtype", _text(str(contiguous.dtype.str)))
            + _frame("shape", b"".join(_int(int(n)) for n in contiguous.shape))
            + _frame("buffer", contiguous.tobytes(order="C"))
        )
        return _frame(f"array:{name}", payload)

    if isinstance(array, (int, float)) and not isinstance(array, bool):
        payload = _frame("scalar", _float(float(array)))
        return _frame(f"array:{name}", payload)

    if isinstance(array, Sequence) and not isinstance(array, (str, bytes, bytearray)):
        payload = _frame("shape", _int(len(array))) + _frame(
            "buffer", b"".join(_float(float(item)) for item in array)
        )
        return _frame(f"array:{name}", payload)

    raise TypeError(
        f"owner state array {name!r} is a {type(array).__name__}, which has no digestible form; "
        "supply a numpy array, a real scalar, or a sequence of reals"
    )


@dataclass(frozen=True, slots=True)
class RngState:
    """The exact position of a random stream, not merely the seed it started from.

    Recording only the seed is the shape of the reference project's bug: two states that started from
    the same seed but have drawn different numbers of variates are *different states*, and a digest
    over the seed alone calls them equal.  :attr:`position` is what distinguishes them.

    Attributes:
        algorithm: The generator's name, for example ``"PCG64"``.  Two generators at the same
            position produce different futures, so the name is part of the identity.
        seed: The seed the stream was created from.
        stream_id: Which independent stream of that seed this is.  Aleph gives each owner its own
            stream, so a state's identity depends on which stream it is reading.
        position: How far into the stream the state has advanced.  A counter of drawn variates or
            the generator's own counter word — either is fine as long as it is monotone and is the
            same convention across a comparison.
        state_words: Optional opaque generator words, when the backend can export them.  Included
            when present, so a state whose generator carries hidden internal state is still fully
            covered.
    """

    algorithm: str
    seed: int
    stream_id: int
    position: int
    state_words: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if not self.algorithm.strip():
            raise ValueError("an RNG state must name its algorithm")
        for name in ("seed", "stream_id", "position"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"RngState.{name} must be an integer; got {value!r}")
        if self.position < 0:
            raise ValueError(f"RNG stream position cannot be negative; got {self.position}")

    def digest_bytes(self) -> bytes:
        """Return the framed bytes this RNG state contributes to a digest."""
        payload = (
            _frame("algorithm", _text(self.algorithm))
            + _frame("seed", _int(self.seed))
            + _frame("stream_id", _int(self.stream_id))
            + _frame("position", _int(self.position))
            + _frame("state_words", b"".join(_frame("w", _int(w)) for w in self.state_words))
        )
        return _frame("rng_state", payload)


@dataclass(frozen=True, slots=True)
class SimulationClock:
    """Where a state sits in time, in both the physical and the bookkeeping sense.

    Attributes:
        step_index: How many candidate steps have been attempted.
        accepted_step_index: How many were accepted.  Kept separately because a state reached after
            three rejections is not the same state as one reached with none, even when the physical
            time matches — and conflating the two is exactly what let the reference project's
            differing attempted-step indices hash alike.
        time_seconds: Physical time [s], hashed by its exact double bits.
    """

    step_index: int
    accepted_step_index: int
    time_seconds: float

    def __post_init__(self) -> None:
        for name in ("step_index", "accepted_step_index"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"SimulationClock.{name} must be an integer; got {value!r}")
            if value < 0:
                raise ValueError(f"SimulationClock.{name} cannot be negative; got {value}")
        if self.accepted_step_index > self.step_index:
            raise ValueError(
                f"accepted steps ({self.accepted_step_index}) cannot exceed attempted steps "
                f"({self.step_index}); a step must be attempted before it can be accepted"
            )

    def digest_bytes(self) -> bytes:
        """Return the framed bytes this clock contributes to a digest."""
        payload = (
            _frame("step_index", _int(self.step_index))
            + _frame("accepted_step_index", _int(self.accepted_step_index))
            + _frame("time_seconds", _float(self.time_seconds))
        )
        return _frame("clock", payload)


@dataclass(frozen=True, slots=True)
class OwnerStateBlock:
    """One owner's local physical state.

    Attributes:
        owner: The owning entity's name, for example ``"membrane"`` or ``"cortex"``.
        arrays: Named state arrays.  Hashed in sorted key order so that two blocks assembled by
            different code paths agree, while each array's own element order is preserved because
            element order in a state array is information.
    """

    owner: str
    arrays: Mapping[str, Any]

    def __post_init__(self) -> None:
        if not self.owner.strip():
            raise ValueError("an owner state block must name its owner")
        for key in self.arrays:
            if not isinstance(key, str):
                raise TypeError(f"owner {self.owner!r} has a non-string array name: {key!r}")

    def digest_bytes(self) -> bytes:
        """Return the framed bytes this owner block contributes to a digest."""
        body = b"".join(
            _array_bytes(name, self.arrays[name]) for name in sorted(self.arrays)
        )
        payload = _frame("owner", _text(self.owner)) + _frame("arrays", body)
        return _frame("owner_block", payload)


def _owner_blocks_bytes(owner_blocks: Iterable[OwnerStateBlock]) -> bytes:
    """Frame the whole set of owner blocks, sorted by owner name.

    Sorting makes the digest independent of the order the caller happened to assemble the owners in,
    which is a presentation detail, while the length-prefixed frames keep two owners' names from
    running together.

    Args:
        owner_blocks: The blocks.

    Returns:
        The framed bytes.

    Raises:
        TypeError: If any element is not an :class:`OwnerStateBlock`.
        ValueError: If two blocks claim the same owner, or if there are none.
    """
    blocks = list(owner_blocks)
    if not blocks:
        raise ValueError(
            "accepted state with no owner blocks has nothing to be a digest of; if a run genuinely "
            "owns no state, it has no accepted state to address"
        )
    names: set[str] = set()
    for block in blocks:
        if not isinstance(block, OwnerStateBlock):
            raise TypeError(f"expected OwnerStateBlock, got {type(block).__name__}")
        if block.owner in names:
            raise ValueError(
                f"owner {block.owner!r} appears twice; each owner contributes exactly one block"
            )
        names.add(block.owner)
    blocks.sort(key=lambda block: block.owner)
    body = b"".join(block.digest_bytes() for block in blocks)
    return _frame("owner_blocks", _frame("count", _int(len(blocks))) + body)


def accepted_state_digest(
    *,
    schema_hash: str,
    context_hash: str,
    topology_generation: int,
    owner_blocks: Iterable[OwnerStateBlock],
    rng_state: RngState,
    clock: SimulationClock,
) -> str:
    """Return the immutable content digest of one accepted state.

    Every argument is keyword-only and has **no default**.  That is the load-bearing design choice:
    the reference project's digest silently omitted two fields because omission was expressible, and
    two genuinely different states came out with the same address.  Here, a caller who forgets a
    component gets a ``TypeError`` at the call, and a maintainer who adds a component without
    declaring it in :data:`DIGEST_COVERAGE` gets a test failure.

    Coverage, in absorption order — the same order as :data:`DIGEST_COVERAGE`:

    1. ``schema_hash`` — which state schema this is.  A state read under a different schema is a
       different state even if the bytes match.
    2. ``context_hash`` — the experimental protocol and external conditions.  Context is *not* part
       of physical state (master plan §9) but it is part of the state's identity: the same
       configuration under a different protocol is a different accepted state.
    3. ``topology_generation`` — increments on every typed topology transition, so a remesh cannot
       hide behind unchanged coordinates.
    4. ``owner_blocks`` — every owner's local arrays, with dtype and shape.
    5. ``rng_state`` — including stream *position*, not just the seed.
    6. ``clock`` — attempted steps, accepted steps, and physical time.

    Args:
        schema_hash: Digest of the state schema, in ``algorithm:hex`` form.
        context_hash: Digest of the context block, in ``algorithm:hex`` form.
        topology_generation: Monotone counter of typed topology transitions.
        owner_blocks: One :class:`OwnerStateBlock` per owning entity.
        rng_state: The random stream's exact position.
        clock: The simulation clock.

    Returns:
        ``"sha256:"`` followed by 64 lowercase hex characters.

    Raises:
        TypeError: If a component has the wrong type.
        ValueError: If a hash argument is empty, if ``topology_generation`` is negative, if
            ``owner_blocks`` is empty, or if two blocks claim the same owner.
    """
    for name, value in (("schema_hash", schema_hash), ("context_hash", context_hash)):
        if not isinstance(value, str):
            raise TypeError(f"{name} must be a string; got {type(value).__name__}")
        if not value.strip():
            raise ValueError(
                f"{name} is empty; an accepted state whose schema or context is unidentified has "
                "no identity to address"
            )
    if isinstance(topology_generation, bool) or not isinstance(topology_generation, int):
        raise TypeError(
            f"topology_generation must be an integer; got {topology_generation!r}"
        )
    if topology_generation < 0:
        raise ValueError(
            f"topology_generation counts transitions and cannot be negative; got "
            f"{topology_generation}"
        )
    if not isinstance(rng_state, RngState):
        raise TypeError("rng_state must be an RngState — a bare seed cannot distinguish positions")
    if not isinstance(clock, SimulationClock):
        raise TypeError("clock must be a SimulationClock")

    stream = (
        _frame("schema", _text(DIGEST_SCHEMA))
        + _frame("schema_hash", _text(schema_hash))
        + _frame("context_hash", _text(context_hash))
        + _frame("topology_generation", _int(topology_generation))
        + _owner_blocks_bytes(owner_blocks)
        + rng_state.digest_bytes()
        + clock.digest_bytes()
    )
    return "sha256:" + hashlib.sha256(stream).hexdigest()
