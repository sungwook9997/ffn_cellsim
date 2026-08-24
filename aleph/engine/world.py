"""Global accepted-step transaction for the component-first cell engine.

The coordinator owns no physical state and never reads the device acceptance predicate on the
host.  Participants receive the same predicate and must predicate restore/commit in their Warp
kernels, preserving one physical clock across independently solved components.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import warp as wp

from aleph.engine.actor import CellActor

_TRANSACTION_HOOKS = ("snapshot_candidate", "rollback", "commit_irreversible")


@dataclass(slots=True)
class CellWorldTransaction:
    """Coordinate candidate snapshot, rollback, irreversible commit, and ledgers.

    Runtime objects are deduplicated by :class:`CellActor`, which matters for the composite
    FA series joint registered under two semantic connector edges.  ``rollback`` and
    ``commit_irreversible`` are both launched for every candidate; their device kernels select
    the rejected and accepted branches respectively without an authoritative GPU-to-host read.
    """

    actor: CellActor
    require_complete: bool = False
    _participants: tuple[Any, ...] = field(init=False, repr=False)
    _ledger_contributors: tuple[Any, ...] = field(init=False, repr=False)
    _registration_signature: tuple[int, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.require_complete:
            self.actor.assert_fully_bound()
        participants: list[Any] = []
        ledger_contributors: list[Any] = []
        runtimes = self.actor.unique_runtime_objects()
        for runtime in runtimes:
            implemented = tuple(callable(getattr(runtime, hook, None)) for hook in _TRANSACTION_HOOKS)
            if any(implemented) and not all(implemented):
                missing = tuple(
                    hook for hook, present in zip(_TRANSACTION_HOOKS, implemented, strict=True)
                    if not present
                )
                raise TypeError(
                    f"runtime {type(runtime).__name__} has an incomplete transaction API; "
                    f"missing {missing}"
                )
            if all(implemented):
                participants.append(runtime)
            if callable(getattr(runtime, "accumulate_ledger", None)):
                ledger_contributors.append(runtime)
        self._participants = tuple(participants)
        self._ledger_contributors = tuple(ledger_contributors)
        self._registration_signature = tuple(id(runtime) for runtime in runtimes)

    def _assert_registry_unchanged(self) -> None:
        signature = tuple(id(runtime) for runtime in self.actor.unique_runtime_objects())
        if signature != self._registration_signature:
            raise RuntimeError("cell actor bindings changed; rebuild the world transaction")

    @property
    def participants(self) -> tuple[Any, ...]:
        """Return unique state owners/joints in deterministic registration order."""
        return self._participants

    def begin_candidate(self) -> None:
        """Snapshot every authoritative state before candidate mechanics or kinetics."""
        self._assert_registry_unchanged()
        for participant in self._participants:
            participant.snapshot_candidate()

    def finalize_candidate(
        self,
        accepted_d: wp.array,
        *,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Launch rejected-state restore and accepted-state irreversible commit.

        ``accepted_d`` remains device-resident; this method must not call ``numpy()``, ``bool()``,
        or otherwise branch on its value from Python.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if isinstance(rng_seed, bool) or not isinstance(rng_seed, int) or rng_seed < 0:
            raise ValueError("rng_seed must be a nonnegative integer")
        self._assert_registry_unchanged()
        for participant in self._participants:
            participant.rollback(accepted_d)
        for participant in self._participants:
            participant.commit_irreversible(accepted_d, dt_phys, rng_seed)

    def accumulate_ledgers(self, ledger: object) -> None:
        """Collect component and connector invariants before the acceptance decision."""
        self._assert_registry_unchanged()
        for contributor in self._ledger_contributors:
            contributor.accumulate_ledger(ledger)
