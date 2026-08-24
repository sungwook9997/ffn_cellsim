"""Executing a declared campaign, resumably, without lying about what happened.

Three properties matter more than speed here.

**Resumable.** A sweep that dies at run 400 of 500 must continue, not restart. Results are appended
to a JSONL journal keyed by run address, and a restart skips addresses already present. Because
`aleph.campaign.rng` addresses streams by coordinates rather than by call order, a resumed run gets
exactly the bits it would have gotten in one pass.

**Failures are results.** A run that raises is recorded as a failed run with its traceback, not
dropped. A campaign that silently omits its failures reports a survivorship-biased summary, and the
omission is invisible precisely when it matters most — when the failures cluster in one corner of
parameter space.

**Refusals are not failures.** A model that correctly declines to answer (out of distribution,
unidentifiable, insufficient data) is behaving properly. Refusals are counted separately so that
"the sweep found nothing" can be distinguished from "the sweep was never able to look".
"""

from __future__ import annotations

import json
import pathlib
import time
import traceback
from collections.abc import Callable, Iterator
from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any

from aleph.campaign.rng import CampaignRNG, StreamAddress
from aleph.campaign.spec import SweepSpec


class RunStatus(StrEnum):
    OK = "ok"
    REFUSED = "refused"
    """The model declined, correctly. Not an error."""
    FAILED = "failed"
    """The run raised. Recorded, never dropped."""


@dataclass(frozen=True, slots=True)
class RunResult:
    point_index: int
    repeat_index: int
    parameters: dict[str, float]
    status: RunStatus
    value: dict[str, Any] = field(default_factory=dict)
    refusal_reason: str = ""
    error: str = ""
    seconds: float = 0.0
    rng_description: dict[str, Any] = field(default_factory=dict)

    @property
    def address(self) -> tuple[int, int]:
        return (self.point_index, self.repeat_index)


class Refusal(Exception):
    """Raise from a campaign task to record a principled decline rather than a failure."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


CampaignTask = Callable[[dict[str, float], StreamAddress, Any], dict[str, Any]]
"""(parameters, stream address, generator) -> result payload. Raise `Refusal` to decline."""


@dataclass(frozen=True, slots=True)
class CampaignSummary:
    spec_name: str
    spec_hash: str
    total: int
    ok: int
    refused: int
    failed: int
    seconds: float
    journal: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def headline(self) -> str:
        return (
            f"{self.spec_name} [{self.spec_hash[:12]}] "
            f"{self.ok}/{self.total} ok, {self.refused} refused, {self.failed} failed, "
            f"{self.seconds:.1f}s"
        )


class CampaignRunner:
    """Runs a `SweepSpec` against a task, journaling every run as it completes."""

    def __init__(self, spec: SweepSpec, journal_path: pathlib.Path) -> None:
        self.spec = spec
        self.journal_path = pathlib.Path(journal_path)
        self.journal_path.parent.mkdir(parents=True, exist_ok=True)
        self.rng = CampaignRNG(spec.seed, common_random_numbers=spec.common_random_numbers)

    def completed_addresses(self) -> set[tuple[int, int]]:
        """Which runs the journal already holds, so a resume skips them."""
        if not self.journal_path.exists():
            return set()
        done: set[tuple[int, int]] = set()
        with self.journal_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    # A partial final line means the process died mid-write. Everything before it
                    # is still good, so stop reading rather than discarding the whole journal.
                    break
                if record.get("spec_hash") != self.spec.spec_hash():
                    continue
                done.add((record["point_index"], record["repeat_index"]))
        return done

    def _append(self, result: RunResult) -> None:
        record = asdict(result)
        record["status"] = str(result.status)
        record["spec_hash"] = self.spec.spec_hash()
        record["spec_name"] = self.spec.name
        with self.journal_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, separators=(",", ":")) + "\n")
            handle.flush()

    def run(self, task: CampaignTask, *, resume: bool = True) -> CampaignSummary:
        started = time.monotonic()
        done = self.completed_addresses() if resume else set()
        counts = {RunStatus.OK: 0, RunStatus.REFUSED: 0, RunStatus.FAILED: 0}

        for point_index, repeat_index, parameters in self.spec.iter_runs():
            if (point_index, repeat_index) in done:
                continue
            address = StreamAddress(point_index=point_index, repeat_index=repeat_index)
            generator = self.rng.generator(address)
            run_started = time.monotonic()
            try:
                payload = task(parameters, address, generator)
                result = RunResult(
                    point_index=point_index,
                    repeat_index=repeat_index,
                    parameters=parameters,
                    status=RunStatus.OK,
                    value=payload,
                    seconds=time.monotonic() - run_started,
                    rng_description=self.rng.describe(address),
                )
            except Refusal as refusal:
                result = RunResult(
                    point_index=point_index,
                    repeat_index=repeat_index,
                    parameters=parameters,
                    status=RunStatus.REFUSED,
                    refusal_reason=refusal.reason,
                    seconds=time.monotonic() - run_started,
                    rng_description=self.rng.describe(address),
                )
            except Exception:  # noqa: BLE001 — a failed run is a result, not a reason to stop
                result = RunResult(
                    point_index=point_index,
                    repeat_index=repeat_index,
                    parameters=parameters,
                    status=RunStatus.FAILED,
                    error=traceback.format_exc(limit=12),
                    seconds=time.monotonic() - run_started,
                    rng_description=self.rng.describe(address),
                )
            counts[result.status] += 1
            self._append(result)

        # Count from the journal rather than from this invocation's tallies.
        #
        # The bug this replaces: `ok=counts[OK] + len(done)` treated every previously-journaled run
        # as a success, because `done` is keyed by address and carries no status. So a resumed
        # campaign over a journal holding 242 ok and 46 refused reported "288 ok, 0 refused" —
        # silently erasing the ok/refused/failed distinction that this module's own docstring says
        # it exists to preserve. It could only ever misreport in the flattering direction, and only
        # on resume, which is the run you are least likely to check.
        journal_counts = {RunStatus.OK: 0, RunStatus.REFUSED: 0, RunStatus.FAILED: 0}
        for result in self.read_results():
            journal_counts[result.status] += 1

        return CampaignSummary(
            spec_name=self.spec.name,
            spec_hash=self.spec.spec_hash(),
            total=self.spec.total_runs(),
            ok=journal_counts[RunStatus.OK],
            refused=journal_counts[RunStatus.REFUSED],
            failed=journal_counts[RunStatus.FAILED],
            seconds=time.monotonic() - started,
            journal=str(self.journal_path),
        )

    def read_results(self) -> Iterator[RunResult]:
        """Replay the journal. Results come back in completion order, which is not run order."""
        if not self.journal_path.exists():
            return
        with self.journal_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    break
                if record.get("spec_hash") != self.spec.spec_hash():
                    continue
                yield RunResult(
                    point_index=record["point_index"],
                    repeat_index=record["repeat_index"],
                    parameters=record["parameters"],
                    status=RunStatus(record["status"]),
                    value=record.get("value", {}),
                    refusal_reason=record.get("refusal_reason", ""),
                    error=record.get("error", ""),
                    seconds=record.get("seconds", 0.0),
                    rng_description=record.get("rng_description", {}),
                )
