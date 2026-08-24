"""The ``aleph-run-record@1`` schema: one typed, self-stamping record per run.

Everything a later reader needs in order to decide whether to believe a number, in one object:
which state was accepted, which code produced it, which configuration it ran, what hardware it ran
on, how much of the population actually existed, what it cost, what the gate decided, where the
result sits on both evidence axes, and which files it wrote.

The record is a frozen dataclass rather than a ``dict`` builder.  The reference project assembled its
records as dictionaries, and the consequence showed up immediately: drivers hand-typed their evidence
rung as a free string, and one of them stamped a rung it had not earned.  A typed record cannot be
hand-typed wrong, because the fields will not accept a string where an enum belongs.

**The anti-overclaiming guard.**  :func:`reject_overclaiming_census` is the part of this module that
exists because of a specific incident rather than a principle.  In the reference project, two runs
1.73x apart in node count both stamped ``fraction_of_native: 1.0``, because the fraction was computed
over one compartment while a different compartment had been quietly reduced.  Every number in both
records was arithmetically correct.  Nothing short of a guard on the *vocabulary* and on the
*decomposition* catches that, so this module bans the names that promise more than they measure and
refuses any fraction that claims completeness without showing its per-owner breakdown.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from enum import Enum, unique
from pathlib import Path
from typing import Any, Mapping, Sequence

from aleph.artifacts.stamp import BuildStamp, config_hash
from aleph.evidence.ladder import EvidenceLabel, EvidenceRung, QuantitativeStatus

__all__ = [
    "ARTIFACT_SCHEMA",
    "ArtifactRef",
    "BANNED_CENSUS_KEYS",
    "DeviceDescription",
    "FRACTION_TOKENS",
    "OverclaimingCensusError",
    "PERCENT_TOKENS",
    "RecordIdentity",
    "RunRecord",
    "TimingBlock",
    "Verdict",
    "is_fraction_key",
    "reject_overclaiming_census",
]

#: Schema identifier, written into every record.  Bump it when the shape changes in a way a reader
#: must notice.  It starts at ``@1`` because Aleph inherits no record history.
ARTIFACT_SCHEMA: str = "aleph-run-record@1"


class OverclaimingCensusError(ValueError):
    """Raised when a census claims a completeness it has not demonstrated."""


@unique
class Verdict(Enum):
    """What one gate evaluation concluded.

    Members:
        PASS: The criterion was met and the measurement was interpretable.
        FAIL: The criterion was not met.  This is a real result about the physics.
        VOID: Nothing interpretable was measured — the numerical error was too large a fraction of
            the signal to separate physics from solver transient.  ``VOID`` is not a soft ``FAIL``.
            ``FAIL`` says the physics missed; ``VOID`` says no claim of either kind is available.
            Keeping it distinct is what stops a non-converged run from reporting a pass.
    """

    PASS = "PASS"
    FAIL = "FAIL"
    VOID = "VOID"


#: Census keys refused by name, each with the reason a record may not carry it.
#:
#: These are banned rather than corrected because the honest replacement is a *decomposition*, not a
#: better scalar.  A single number that claims whole-system completeness will always be computed over
#: whichever part the author was looking at.
BANNED_CENSUS_KEYS: Mapping[str, str] = {
    "fraction_of_native": (
        "a single scalar cannot say whether a multi-owner population is complete: it will be "
        "computed over whichever owner the author had in mind, and will read 1.0 while another "
        "owner sits at a quarter resolution. Record a per-owner breakdown and let the reader see "
        "the whole system"
    ),
    "fraction_complete": (
        "'complete' names a scope the census never declares. Record the realized count for every "
        "owner alongside the target count, so completeness is something the reader computes rather "
        "than something the record asserts"
    ),
    "percent_complete": (
        "same defect as fraction_complete, and additionally in percent, which invites a run at "
        "0.85 to be recorded as 85 and later read as a fraction"
    ),
    "is_native": (
        "a boolean cannot carry the population it is asserting about. Record the realized counts "
        "and the native targets separately"
    ),
    "all_owners_present": (
        "a boolean asserting completeness of an enumeration the census is already required to "
        "carry; if every owner is present, the enumeration shows it"
    ),
    "coverage_fraction": (
        "dispatch coverage is a structural fact and a census is a population fact. Mixing them is "
        "how a connector that is never evaluated passes a coverage gate"
    ),
}

#: Name tokens that make a census key a fraction, and therefore subject to the completeness rule.
FRACTION_TOKENS: frozenset[str] = frozenset(
    {"fraction", "frac", "ratio", "coverage", "share", "completeness"}
)

#: Name tokens Aleph refuses outright in a census.  A percentage in a machine-read record is a
#: scaling ambiguity waiting to happen: 0.85 and 85 are both plausible readings of the same field.
PERCENT_TOKENS: frozenset[str] = frozenset({"pct", "percent", "percentage"})

_TOKEN_SPLIT: re.Pattern[str] = re.compile(r"[^a-z0-9]+")

#: The suffix a fraction key's per-owner breakdown must carry.
_BREAKDOWN_SUFFIX: str = "_by_owner"


def _tokens(key: str) -> set[str]:
    """Split a census key into lowercase name tokens."""
    return {token for token in _TOKEN_SPLIT.split(key.lower()) if token}


def is_fraction_key(key: str) -> bool:
    """Return whether a census key names a fraction.

    Args:
        key: The census key.

    Returns:
        ``True`` when any of :data:`FRACTION_TOKENS` appears as a name token.  Token matching rather
        than substring matching, so ``"fractional_step"`` is caught but ``"refractory_count"`` — a
        legitimate biological name that happens to contain ``frac`` — is not.
    """
    return bool(_tokens(key) & FRACTION_TOKENS)


def _as_real(value: object) -> float | None:
    """Return ``value`` as a float when it is a real number, else ``None``."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    item = getattr(value, "item", None)
    if callable(item):
        try:
            unwrapped = item()
        except (TypeError, ValueError):
            return None
        if isinstance(unwrapped, (int, float)) and not isinstance(unwrapped, bool):
            return float(unwrapped)
    return None


def reject_overclaiming_census(census: Mapping[str, Any]) -> None:
    """Raise if a census claims more completeness than it demonstrates.

    Four checks, in order:

    1. **Banned names.** Any key in :data:`BANNED_CENSUS_KEYS` is refused, with the reason.
    2. **Percentages.** Any key carrying a token from :data:`PERCENT_TOKENS` is refused, because
       ``0.85`` and ``85`` are both plausible readings of the same field and nothing in the record
       says which was meant.
    3. **Fractions need a population.** A census containing any fraction key must carry
       ``owners``: a non-empty sequence of distinct owner names.  A fraction with no stated
       denominator population cannot be checked by anyone.
    4. **A fraction of 1.0 must show its work.** A fraction key at exactly ``1.0`` must be
       accompanied by ``<key>_by_owner``, a mapping covering every name in ``owners``, with every
       owner also at ``1.0``. An aggregate of 1.0 over a breakdown containing anything less is the
       exact defect this guard exists for: it is arithmetically consistent, and it is false.

    Fractions strictly below 1.0 are left alone. They are not claiming completeness, so they are not
    overclaiming, and forcing a breakdown on every partial fraction would make the guard so noisy
    that someone would turn it off.

    Args:
        census: The population census a run assembled.

    Raises:
        OverclaimingCensusError: On any of the four failures, naming the key and what to record
            instead.
        TypeError: If ``census`` is not a mapping, or ``owners`` is not a sequence of strings.
    """
    if not isinstance(census, Mapping):
        raise TypeError(f"census must be a mapping; got {type(census).__name__}")

    for key in census:
        if not isinstance(key, str):
            raise TypeError(f"census keys must be strings; got {key!r}")

    for banned, why in BANNED_CENSUS_KEYS.items():
        if banned in census:
            raise OverclaimingCensusError(
                f"census key {banned!r} may not be recorded: {why}"
            )

    for key in census:
        offending = _tokens(key) & PERCENT_TOKENS
        if offending:
            raise OverclaimingCensusError(
                f"census key {key!r} is a percentage ({sorted(offending)[0]!r}); record a fraction "
                "in [0, 1] instead, because a percentage in a machine-read record cannot say "
                "whether 0.85 or 85 was meant"
            )

    fraction_keys = sorted(
        key for key in census if is_fraction_key(key) and not key.endswith(_BREAKDOWN_SUFFIX)
    )
    if not fraction_keys:
        return

    owners = census.get("owners")
    if owners is None:
        raise OverclaimingCensusError(
            f"census carries fraction key(s) {fraction_keys} but no 'owners' enumeration; a "
            "fraction whose denominator population is unstated cannot be checked by a reader"
        )
    if isinstance(owners, (str, bytes)) or not isinstance(owners, Sequence):
        raise TypeError(f"census['owners'] must be a sequence of owner names; got {owners!r}")
    owner_names = [str(name) for name in owners]
    if not owner_names:
        raise OverclaimingCensusError(
            "census['owners'] is empty; a census of nothing cannot support a fraction"
        )
    if len(set(owner_names)) != len(owner_names):
        raise OverclaimingCensusError(
            f"census['owners'] repeats a name: {owner_names}; each owner appears exactly once"
        )

    for key in fraction_keys:
        value = _as_real(census[key])
        if value is None:
            raise OverclaimingCensusError(
                f"census fraction key {key!r} is not a real number; got {census[key]!r}"
            )
        if value != value or not (0.0 <= value <= 1.0):
            raise OverclaimingCensusError(
                f"census fraction key {key!r} must lie in [0, 1]; got {value!r}"
            )
        if value != 1.0:
            continue

        breakdown_key = f"{key}{_BREAKDOWN_SUFFIX}"
        breakdown = census.get(breakdown_key)
        if breakdown is None:
            raise OverclaimingCensusError(
                f"census key {key!r} claims 1.0 but there is no {breakdown_key!r} showing the "
                f"per-owner values for {owner_names}. A completeness claim computed over one owner "
                "reads 1.0 while another owner is reduced, and every number in the record stays "
                "arithmetically correct while the record is false"
            )
        if not isinstance(breakdown, Mapping):
            raise TypeError(
                f"census[{breakdown_key!r}] must be a mapping from owner name to fraction; got "
                f"{type(breakdown).__name__}"
            )
        missing = [name for name in owner_names if name not in breakdown]
        if missing:
            raise OverclaimingCensusError(
                f"census key {key!r} claims 1.0 but {breakdown_key!r} has no entry for {missing}; "
                "a breakdown that skips an owner is how the reduced owner stays invisible"
            )
        extra = [name for name in breakdown if name not in owner_names]
        if extra:
            raise OverclaimingCensusError(
                f"{breakdown_key!r} names {extra}, which are not in census['owners']; the "
                "breakdown and the owner enumeration must describe the same system"
            )
        for name in owner_names:
            share = _as_real(breakdown[name])
            if share is None or share != share or not (0.0 <= share <= 1.0):
                raise OverclaimingCensusError(
                    f"{breakdown_key!r}[{name!r}] must be a fraction in [0, 1]; got "
                    f"{breakdown[name]!r}"
                )
            if share != 1.0:
                raise OverclaimingCensusError(
                    f"census key {key!r} claims 1.0 while {breakdown_key!r}[{name!r}] is {share}; "
                    "an aggregate cannot be complete while an owner is not. Record the aggregate "
                    "as the minimum over owners, or drop the aggregate and keep the breakdown"
                )


@dataclass(frozen=True, slots=True)
class RecordIdentity:
    """What state this record is about.

    Attributes:
        accepted_state_digest: From :func:`aleph.artifacts.digest.accepted_state_digest`.  ``None``
            only when the run produced no accepted state, which a reader must be able to see.
        schema_hash: Digest of the state schema the run used.
        context_hash: Digest of the experimental context.  Carried alongside the state digest even
            though the digest already covers it, so a reader can group records by context without
            recomputing anything.
    """

    schema_hash: str
    context_hash: str
    accepted_state_digest: str | None = None

    def __post_init__(self) -> None:
        for name in ("schema_hash", "context_hash"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"RecordIdentity.{name} is required")

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able identity block."""
        return {
            "schema_hash": self.schema_hash,
            "context_hash": self.context_hash,
            "accepted_state_digest": self.accepted_state_digest,
            "has_accepted_state": self.accepted_state_digest is not None,
        }


@dataclass(frozen=True, slots=True)
class DeviceDescription:
    """What the run executed on, recorded verbatim rather than inferred.

    Attributes:
        kind: ``"cpu"`` or ``"gpu"``.
        name: The device name as the driver reported it.
        backend: The compute backend and version.
        precision: The compute/accumulation precision contract, for example
            ``"float32 compute / float64 accumulation"``.  Recorded because a magnitude at one
            precision is not the same measurement as the same magnitude at another.
        note: Anything a later reader would need — an instrumented run is not a benchmark and should
            say so here.
    """

    kind: str
    name: str
    backend: str
    precision: str
    note: str | None = None

    def __post_init__(self) -> None:
        if self.kind not in ("cpu", "gpu"):
            raise ValueError(f"device kind must be 'cpu' or 'gpu'; got {self.kind!r}")
        for name in ("name", "backend", "precision"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"DeviceDescription.{name} is required")

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able device block."""
        block: dict[str, Any] = {
            "kind": self.kind,
            "name": self.name,
            "backend": self.backend,
            "precision": self.precision,
        }
        if self.note:
            block["note"] = self.note
        return block


@dataclass(frozen=True, slots=True)
class TimingBlock:
    """What the run cost, with the derived rates computed here rather than by each reader.

    A cost is a property of a run *at a stated convergence*, never of an engine, because any solver
    is arbitrarily fast if allowed to stop early.  :meth:`RunRecord.as_json_obj` therefore stamps
    this block with the run's own verdict and marks it non-comparable unless that verdict was
    ``PASS``.

    Attributes:
        wall_seconds: Measured wall clock [s].  Positive and finite.
        physical_time_seconds: Physical time the run covered [s], if it stepped a clock.
        n_attempted_steps: Candidate steps proposed.
        n_accepted_steps: Candidate steps accepted.
    """

    wall_seconds: float
    physical_time_seconds: float | None = None
    n_attempted_steps: int | None = None
    n_accepted_steps: int | None = None

    def __post_init__(self) -> None:
        wall = self.wall_seconds
        if not isinstance(wall, (int, float)) or isinstance(wall, bool):
            raise TypeError(f"wall_seconds must be a real number; got {wall!r}")
        if not (wall > 0.0) or wall != wall or wall == float("inf"):
            raise ValueError(f"wall_seconds must be positive and finite; got {wall!r}")
        if self.physical_time_seconds is not None and not (self.physical_time_seconds > 0.0):
            raise ValueError(
                f"physical_time_seconds must be positive when given; got "
                f"{self.physical_time_seconds!r}"
            )
        for name in ("n_attempted_steps", "n_accepted_steps"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
                raise TypeError(f"TimingBlock.{name} must be an integer; got {value!r}")
            if value is not None and value < 0:
                raise ValueError(f"TimingBlock.{name} cannot be negative; got {value}")
        if (
            self.n_attempted_steps is not None
            and self.n_accepted_steps is not None
            and self.n_accepted_steps > self.n_attempted_steps
        ):
            raise ValueError(
                f"accepted steps ({self.n_accepted_steps}) cannot exceed attempted steps "
                f"({self.n_attempted_steps})"
            )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able timing block with derived rates."""
        block: dict[str, Any] = {"wall_seconds": float(self.wall_seconds)}
        if self.physical_time_seconds is not None:
            block["physical_time_seconds"] = float(self.physical_time_seconds)
            block["wall_seconds_per_physical_second"] = float(self.wall_seconds) / float(
                self.physical_time_seconds
            )
        if self.n_attempted_steps is not None:
            block["n_attempted_steps"] = int(self.n_attempted_steps)
        if self.n_accepted_steps is not None:
            block["n_accepted_steps"] = int(self.n_accepted_steps)
            if self.n_accepted_steps > 0:
                block["wall_seconds_per_accepted_step"] = float(self.wall_seconds) / float(
                    self.n_accepted_steps
                )
        if self.n_attempted_steps:
            block["acceptance_rate"] = (
                None
                if self.n_accepted_steps is None
                else float(self.n_accepted_steps) / float(self.n_attempted_steps)
            )
        return block


@dataclass(frozen=True, slots=True)
class ArtifactRef:
    """One file the run produced, addressed by content.

    Attributes:
        role: What the file is for, for example ``"spectrum"`` or ``"accepted_state"``.
        path: Path as written, relative to the data root where possible.
        digest: Content digest in ``algorithm:hex`` form.  Required: a path without a digest names a
            location, not a file, and locations get overwritten.
        size_bytes: File size, when known.
    """

    role: str
    path: str
    digest: str
    size_bytes: int | None = None

    def __post_init__(self) -> None:
        for name in ("role", "path", "digest"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"ArtifactRef.{name} is required")
        if ":" not in self.digest:
            raise ValueError(
                f"artifact digest must carry its algorithm, as 'algorithm:hex'; got {self.digest!r}"
            )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able artifact reference."""
        block: dict[str, Any] = {"role": self.role, "path": self.path, "digest": self.digest}
        if self.size_bytes is not None:
            block["size_bytes"] = int(self.size_bytes)
        return block


@dataclass(frozen=True, slots=True)
class RunRecord:
    """One run, fully stamped.  The ``aleph-run-record@1`` schema.

    The census is validated at construction, so an overclaiming record cannot be built at all — not
    built and then rejected at write time, which would leave the offending object alive in memory
    for something else to serialize.

    Attributes:
        run_label: Free-text identity of this run.
        identity: Which state this record is about.
        build: The build stamp from :func:`aleph.artifacts.stamp.build_stamp`.
        config: The configuration, recorded in full and hashed into ``config_digest``.
        device: What it ran on.
        census: Counts of every owned entity.  Validated by
            :func:`reject_overclaiming_census`.
        timing: What it cost.
        peak_memory_bytes: Peak device memory, when the backend reports it.
        evidence: The two-axis label carrying the rung and the quantitative status.
        verdict: The gate outcome, or ``None`` when the run evaluated no gate.  ``None`` rather than
            a defaulted ``PASS``: a run with no gate has not passed one.
        artifacts: Files the run produced.
        notes: Anything else worth recording.
    """

    run_label: str
    identity: RecordIdentity
    build: BuildStamp
    config: Mapping[str, Any]
    device: DeviceDescription
    census: Mapping[str, Any]
    timing: TimingBlock
    evidence: EvidenceLabel
    peak_memory_bytes: int | None = None
    verdict: Verdict | None = None
    artifacts: tuple[ArtifactRef, ...] = ()
    notes: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        if not self.run_label.strip():
            raise ValueError("a run record needs a run_label")
        if not isinstance(self.identity, RecordIdentity):
            raise TypeError("identity must be a RecordIdentity")
        if not isinstance(self.build, BuildStamp):
            raise TypeError("build must be a BuildStamp from build_stamp()")
        if not isinstance(self.device, DeviceDescription):
            raise TypeError("device must be a DeviceDescription")
        if not isinstance(self.timing, TimingBlock):
            raise TypeError("timing must be a TimingBlock")
        if not isinstance(self.evidence, EvidenceLabel):
            raise TypeError(
                "evidence must be an EvidenceLabel from classify_claim(), not a free string; a "
                "hand-typed rung records an intention rather than an observation"
            )
        if self.verdict is not None and not isinstance(self.verdict, Verdict):
            raise TypeError("verdict must be a Verdict or None")
        if self.peak_memory_bytes is not None and (
            isinstance(self.peak_memory_bytes, bool)
            or not isinstance(self.peak_memory_bytes, int)
            or self.peak_memory_bytes < 0
        ):
            raise ValueError(
                f"peak_memory_bytes must be a non-negative integer; got {self.peak_memory_bytes!r}"
            )
        for artifact in self.artifacts:
            if not isinstance(artifact, ArtifactRef):
                raise TypeError("every entry of artifacts must be an ArtifactRef")
        reject_overclaiming_census(self.census)

    @property
    def rung(self) -> EvidenceRung:
        """The structural rung this run earned."""
        return self.evidence.rung

    @property
    def quantitative_status(self) -> QuantitativeStatus:
        """Whether this run's magnitudes may be quoted."""
        return self.evidence.status

    @property
    def config_digest(self) -> str:
        """The canonical hash of the configuration, computed on demand from the config itself.

        Derived rather than stored so it cannot drift from the config it is supposed to address.
        """
        return config_hash(self.config)

    def _timing_block(self) -> dict[str, Any]:
        """Return the timing block stamped with this run's verdict.

        A cost lifted out of a run that did not converge is a measurement of not converging, so the
        block says explicitly whether it may be compared against another run.
        """
        block = self.timing.as_json_obj()
        if self.verdict is Verdict.PASS:
            block["comparable"] = True
            block["comparable_against"] = (
                "another run whose gate also PASSED under the same contract and observable"
            )
        else:
            block["comparable"] = False
            block["not_comparable_reason"] = (
                f"the gate verdict is {self.verdict.value if self.verdict else 'absent'}, so this "
                "run did not establish the convergence its cost would have to be quoted at"
            )
        return block

    def as_json_obj(self) -> dict[str, Any]:
        """Return the whole record as a JSON-able object.

        Returns:
            A mapping whose top-level keys are stable across ``@1``.  Both evidence axes are always
            present, and ``verdict`` is explicitly ``null`` when no gate was evaluated.
        """
        record: dict[str, Any] = {
            "schema": ARTIFACT_SCHEMA,
            "run_label": self.run_label,
            "identity": self.identity.as_json_obj(),
            "build": self.build.as_json_obj(),
            "device": self.device.as_json_obj(),
            "config": dict(self.config),
            "config_digest": self.config_digest,
            "census": dict(self.census),
            "timing": self._timing_block(),
            "peak_memory_bytes": self.peak_memory_bytes,
            "verdict": None if self.verdict is None else self.verdict.value,
            "artifacts": [artifact.as_json_obj() for artifact in self.artifacts],
        }
        record.update(self.evidence.as_artifact_fields())
        if self.notes:
            record["notes"] = dict(self.notes)
        return record

    def to_json(self, *, indent: int = 2) -> str:
        """Return the record as JSON text, newline-terminated.

        Args:
            indent: Indentation passed to :func:`json.dumps`.

        Returns:
            The JSON text.

        Raises:
            TypeError: If any value in ``config``, ``census``, or ``notes`` is not JSON-able.  There
                is no ``default=str`` fallback: a record that silently stringifies an object it does
                not understand is a record whose contents cannot be trusted to mean what they say.
        """
        return json.dumps(self.as_json_obj(), indent=indent, sort_keys=False) + "\n"

    def write_json(self, path: Path | str) -> Path:
        """Write the record to ``path`` as JSON, creating the parent directory.

        Args:
            path: Destination file.

        Returns:
            The written path.
        """
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination
