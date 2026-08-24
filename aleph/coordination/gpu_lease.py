r"""A single-holder, self-expiring lease on the one CUDA device, plus the run index beside it.

WHY THIS EXISTS.  There is exactly ONE A5000 and there are now several concurrent sessions.  Before
this module the protocol was "whoever launches first wins": two native runs would land on the same
device, contend for its 16 GB, and — worse — neither artifact would record that it had shared the
device with anything.  A run that swapped or was throttled by a neighbour is not the run its config
says it is, and nothing downstream could tell.

The second half of the same problem is duplication.  The engine's run records already stamp the build
commit and a canonical hash of the config (:func:`~aleph.engine.observe.artifact.config_hash`), so
"has this exact configuration already been measured on this exact source?" has had a mechanical answer
for some time — but nothing ever asked it, and the observed result was sessions re-running work that
existed and overwriting each other's output directories.  :class:`RunIndex` is the thing that asks.

DESIGN, and why each choice is the conservative one:

* **Self-expiring.**  A lease carries an absolute expiry.  A session that crashes, is interrupted, or
  simply forgets to release cannot wedge the device — the worst case is that the GPU idles until the
  expiry passes.  A lock without expiry would make a dead session an outage.
* **Atomic acquisition.**  The lease file is created with ``O_CREAT | O_EXCL`` and, if that fails,
  taken over ONLY after reading an expiry that has already passed.  Two simultaneous acquirers cannot
  both win, because exactly one ``O_EXCL`` create succeeds.
* **Holder-only release, but expiry always wins.**  Releasing someone else's live lease raises;
  reclaiming an EXPIRED lease is allowed and is the normal recovery path.
* **The clock is the holder's, and it is recorded.**  Both machines' clocks appear in the record, so a
  skew shows up as data rather than as a mysterious early expiry.
* **No physics, no device call.**  This module never imports Warp and never touches CUDA.  It is a
  coordination primitive over a filesystem path, which is why it is exercisable with no GPU and why
  the lease can live on the machine that HAS the GPU while the deciding code runs anywhere.

WHAT IT DELIBERATELY DOES NOT DO.  It does not schedule, queue, or prioritise.  A caller that cannot
acquire is told who holds it and for how long, and decides for itself whether to wait — the failure
mode of an automatic queue is a session blocking for hours on a lease whose holder died in a way the
expiry has not yet caught.

Sanity Gate:
    * boundary: a zero or negative duration is rejected; an unreadable or corrupt lease file is treated
      as ABSENT rather than as held, because a coordination primitive that fails closed on its own
      corruption turns one bad write into a permanent outage.
    * conservation/invariant: at most one holder at any instant — enforced by ``O_EXCL``, not by a
      read-then-write sequence, which would race.
    * numerical: times are UTC epoch seconds as float; no timezone arithmetic anywhere.
    * measurement-protocol: :meth:`RunIndex.find` matches on (build commit, config hash) TOGETHER.
      Either alone is insufficient — the same config on a different build is a different run, and the
      same build with a different config obviously is.
    * sign-sense: not applicable.
"""

from __future__ import annotations

import json
import os
import socket
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterator

__all__ = [
    "DEFAULT_LEASE_PATH",
    "DEFAULT_RUN_INDEX_PATH",
    "Lease",
    "LeaseHeld",
    "RunIndex",
    "RunRecord",
    "acquire",
    "read_lease",
    "release",
]

#: Where the lease lives on the machine that owns the GPU.  A dotfile under ``$HOME`` rather than
#: inside the repo, so a checkout, a branch switch, or a clean does not silently free a live lease.
DEFAULT_LEASE_PATH = Path.home() / ".ffn" / "gpu.lease"

#: Append-only index of runs that actually executed, beside the lease for the same reason.
DEFAULT_RUN_INDEX_PATH = Path.home() / ".ffn" / "runs.jsonl"


class LeaseHeld(RuntimeError):
    """Raised when the device is held by someone else and the lease has not expired.

    Carries the live :class:`Lease` so the caller can report who holds it and for how long rather
    than just failing.
    """

    def __init__(self, lease: "Lease") -> None:
        remaining = lease.seconds_remaining()
        super().__init__(
            f"GPU is leased by {lease.holder!r} on {lease.host} for another {remaining:.0f} s "
            f"(reason: {lease.reason or 'unstated'})"
        )
        self.lease = lease


@dataclass(frozen=True, slots=True)
class Lease:
    """One holder's claim on the device.

    Attributes:
        holder: Session identity — free text, but it is what a waiting session will read, so it
            should name something a human can go ask.
        host: Hostname the lease was written from.
        pid: Process id of the acquiring process, for the case where the holder is on this machine.
        acquired_at: UTC epoch seconds when the claim was made.
        expires_at: UTC epoch seconds after which anyone may take over.
        reason: What the device is being used for; recorded so a waiting session knows whether to
            wait or to do something else.
    """

    holder: str
    host: str
    pid: int
    acquired_at: float
    expires_at: float
    reason: str = ""

    def seconds_remaining(self, now: float | None = None) -> float:
        """Return seconds until expiry, floored at zero.

        Args:
            now: UTC epoch seconds to evaluate against; defaults to the current time.

        Returns:
            Remaining seconds, never negative.
        """
        return max(0.0, self.expires_at - (time.time() if now is None else now))

    def is_live(self, now: float | None = None) -> bool:
        """Whether this lease still holds the device."""
        return self.seconds_remaining(now) > 0.0


def _lease_path(path: Path | str | None) -> Path:
    return Path(path) if path is not None else DEFAULT_LEASE_PATH


def read_lease(path: Path | str | None = None) -> Lease | None:
    """Return the lease on disk, or ``None`` if absent or unreadable.

    An unreadable or malformed lease is reported as ABSENT, not as held: a coordination primitive that
    fails closed on its own corruption converts one bad write into an outage that only a human can
    clear. The cost of the other direction is bounded — a corrupt lease is takeable, and takeover is
    recorded.

    Args:
        path: Lease file; defaults to :data:`DEFAULT_LEASE_PATH`.

    Returns:
        The :class:`Lease`, or ``None``.
    """
    target = _lease_path(path)
    try:
        payload = json.loads(target.read_text())
        return Lease(**{k: payload[k] for k in Lease.__slots__ if k in payload})
    except (OSError, ValueError, TypeError, KeyError):
        return None


def acquire(
    holder: str,
    *,
    minutes: float,
    reason: str = "",
    path: Path | str | None = None,
) -> Lease:
    """Take the device lease, or raise :class:`LeaseHeld` if someone else holds it live.

    Acquisition is atomic: the file is created with ``O_CREAT | O_EXCL`` so exactly one of any number
    of simultaneous acquirers succeeds. If the create fails because a lease already exists, the
    existing one is read and taken over ONLY when it has expired.

    Args:
        holder: Session identity, non-empty.
        minutes: Lease duration. Must be positive. Pick the duration the work actually needs — a lease
            far longer than the job idles the device if the session dies, and one far shorter expires
            mid-run and lets a neighbour land on the same device.
        reason: What the device is for, recorded for whoever is waiting.
        path: Lease file; defaults to :data:`DEFAULT_LEASE_PATH`.

    Returns:
        The acquired :class:`Lease`.

    Raises:
        ValueError: If ``holder`` is empty or ``minutes`` is not positive-finite.
        LeaseHeld: If a live lease is held by anyone, including this same holder — re-acquiring your
            own live lease is a bug (two runs in one session), not an extension.
    """
    if not holder.strip():
        raise ValueError("a lease needs a holder a human can go ask")
    if not (minutes > 0.0) or minutes != minutes or minutes == float("inf"):
        raise ValueError(f"minutes must be positive-finite; got {minutes!r}")

    target = _lease_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    lease = Lease(
        holder=holder.strip(), host=socket.gethostname(), pid=os.getpid(),
        acquired_at=now, expires_at=now + minutes * 60.0, reason=reason.strip(),
    )
    payload = json.dumps(asdict(lease), indent=2) + "\n"

    try:
        fd = os.open(target, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        existing = read_lease(target)
        if existing is not None and existing.is_live(now):
            raise LeaseHeld(existing) from None
        # Expired or corrupt: take over. The write is still atomic via a temp file + rename, so a
        # reader never observes a half-written lease.
        temp = target.with_suffix(f".{os.getpid()}.tmp")
        temp.write_text(payload)
        os.replace(temp, target)
        return lease
    with os.fdopen(fd, "w") as handle:
        handle.write(payload)
    return lease


def release(holder: str, *, path: Path | str | None = None, force: bool = False) -> bool:
    """Release the lease if ``holder`` owns it, or if it has already expired.

    Args:
        holder: The releasing session's identity.
        path: Lease file; defaults to :data:`DEFAULT_LEASE_PATH`.
        force: Release even a live lease held by someone else. Reserved for a human clearing a lease
            whose holder is known dead before its expiry — it is exactly the operation the expiry
            exists to make unnecessary, so it is never the automatic path.

    Returns:
        ``True`` if a lease was removed, ``False`` if there was nothing to remove.

    Raises:
        LeaseHeld: If a live lease belongs to someone else and ``force`` is not set.
    """
    target = _lease_path(path)
    existing = read_lease(target)
    if existing is None:
        target.unlink(missing_ok=True)
        return False
    if existing.is_live() and existing.holder != holder.strip() and not force:
        raise LeaseHeld(existing)
    target.unlink(missing_ok=True)
    return True


@dataclass(frozen=True, slots=True)
class RunRecord:
    """One executed run, as the index remembers it.

    Attributes:
        run_label: What the run called itself.
        build_commit: The source commit it was measured on.
        config_sha256: Canonical hash of its configuration.
        artifact: Path the run wrote its record to.
        device: Device string the run resolved.
        started_at: UTC epoch seconds.
        holder: Which session ran it.
        population: Free-text population census summary — recorded because the 2026-07-28 retraction
            turned on nobody being able to see, at a glance, that a result came from 0.18% of native.
        closure_sha256: Digest of the CONTENT of every first-party file the driver imports. This is the
            code that actually ran; ``build_commit`` is only the repo HEAD at launch, so two runs of
            byte-identical code get different commits whenever an unrelated file was committed between
            them. A sweep comparing its points by commit therefore refuses runs that are in fact
            identical, and a sweep trusting the commit would miss a change confined to one file. Added
            2026-07-29 after a dose-response check refused two points on exactly that false difference.
    """

    run_label: str
    build_commit: str
    config_sha256: str
    artifact: str
    device: str
    started_at: float
    holder: str
    population: str = ""
    closure_sha256: str = ""


class RunIndex:
    """Append-only JSONL index answering "has this exact run already happened?".

    The question is (build commit, config hash) together — either alone answers a different question.
    """

    def __init__(self, path: Path | str | None = None) -> None:
        """Open (without creating) the index.

        Args:
            path: Index file; defaults to :data:`DEFAULT_RUN_INDEX_PATH`.
        """
        self.path = Path(path) if path is not None else DEFAULT_RUN_INDEX_PATH

    def __iter__(self) -> Iterator[RunRecord]:
        """Yield every well-formed record, skipping lines that are not.

        A partially written final line — the normal consequence of a process dying mid-append — must
        not make the whole index unreadable, so unparseable lines are skipped rather than raised on.
        """
        try:
            text = self.path.read_text()
        except OSError:
            return
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
                yield RunRecord(**{k: payload[k] for k in RunRecord.__slots__ if k in payload})
            except (ValueError, TypeError, KeyError):
                continue

    def find(self, *, build_commit: str, config_sha256: str) -> list[RunRecord]:
        """Return every prior run matching this build AND this config.

        Args:
            build_commit: Source commit.
            config_sha256: Canonical config hash, as produced by the engine's artifact writer.

        Returns:
            Matching records, oldest first. Empty means this run is new.
        """
        return [r for r in self
                if r.build_commit == build_commit and r.config_sha256 == config_sha256]

    def append(self, record: RunRecord) -> RunRecord:
        """Append one record. Creates the file and its parent if needed.

        Args:
            record: The run to remember.

        Returns:
            The record, unchanged, for chaining.
        """
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a") as handle:
            handle.write(json.dumps(asdict(record), default=str) + "\n")
        return record

    def summary(self) -> dict[str, Any]:
        """Return counts a dashboard can render without re-reading the file."""
        records = list(self)
        return {
            "n_runs": len(records),
            "n_distinct_configs": len({r.config_sha256 for r in records}),
            "n_builds": len({r.build_commit for r in records}),
            "holders": sorted({r.holder for r in records if r.holder}),
            "path": str(self.path),
        }
