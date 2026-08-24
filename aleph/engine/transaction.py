"""CellTransaction — the whole-cell accepted-step orchestrator (all 14 components / 36 connectors).

Spec §1 (``WHOLE_CELL_COMMON_CONTRACTS_SPEC_2026-07-23.md``).  ``ecm_world.py`` already runs a
multi-participant accepted-step transaction (snapshot ECM + crosslink + boundary + contact + clutch as ONE
candidate, one device ``accepted`` predicate to every ``rollback``, one ``commit_irreversible`` finalize).
This module lifts that pattern to the *whole* cell: one top-level :meth:`CellTransaction.step` that drives
every deduplicated state-owner and connector runtime under one device-resident acceptance predicate, in the
ratified order::

    snapshot_candidate (all)  ->  propose_events (all, into device candidate buffers)
                              ->  inner mechanical solve (caller-owned; accumulate -> converge)
                              ->  accumulate_ledger (all)  ->  GlobalCellLedger.assemble_balance(tol^2)
                              ->  finalize: rollback + commit_irreversible (all)   [device selects branch]
                              ->  advance event clock + RNG epoch (device-gated on the same predicate)
                              ->  assert disjoint populations (accepted-step no-double-count invariant)

No new transaction machinery is introduced: the snapshot/rollback/commit fan-out and participant
deduplication are :class:`~aleph.engine.world.CellWorldTransaction`; the balance gate is
:class:`~aleph.engine.ledger.GlobalCellLedger`; the clock is an
:class:`~aleph.engine.runtime.EventClock`.  This orchestrator owns no physical state, advances no
time itself, and never reads the acceptance predicate on the host — both ``rollback`` and
``commit_irreversible`` are launched every step and select their branch in-device, and the clock advance is
device-gated.  The inner mechanical solve is a caller-supplied callable (the engine's physics), so this
module is CPU-importable and the structural gates drive it with spies.

Sanity Gate:
    * ownership: participants + ledger contributors come from one :class:`CellWorldTransaction`; the composite
      FA joint is one participant, not two (deduplicated upstream).
    * boundary/sign: the acceptance predicate is device-resident; ``step`` never calls ``bool()``/``numpy()``
      on it — rollback and commit both run and select in-device, the clock advances only when it selects accept.
    * conservation: when population ledgers are supplied, :func:`assert_disjoint_populations` runs every
      accepted step, promoting the viz-time no-double-count check to a hot-loop invariant.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import warp as wp

from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.engine.world import CellWorldTransaction

__all__ = ["CellTransaction"]

_EVENT_HOOK = "propose_events"


@dataclass(slots=True)
class CellTransaction:
    """Top-level whole-cell physical-step transaction over a composed cell world.

    Args:
        world: the deduplicated snapshot/rollback/commit + ledger coordinator for every component
            state-owner and connector runtime.
        clock: the device-resident event clock / RNG epoch, advanced only on an accepted step.
        population_ledgers: optional per-component population ledgers; when given, the disjoint-ID
            no-double-count invariant is asserted at build time and after every accepted step.
    """

    world: CellWorldTransaction
    clock: Any
    population_ledgers: tuple[PopulationLedger, ...] = ()
    _proposers: tuple[Any, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not callable(getattr(self.clock, "advance", None)):
            raise TypeError("event clock must expose a callable advance(accepted, dt_phys)")
        base_seed = getattr(self.clock, "base_seed", None)
        if isinstance(base_seed, bool) or not isinstance(base_seed, int) or base_seed < 0:
            raise TypeError("event clock must expose a nonnegative int base_seed")
        self.population_ledgers = tuple(self.population_ledgers)
        if self.population_ledgers:
            assert_disjoint_populations(self.population_ledgers)  # build-time invariant
        self._proposers = self.world.actor.runtimes_with(_EVENT_HOOK)

    @property
    def proposers(self) -> tuple[Any, ...]:
        """Deduplicated runtimes exposing ``propose_events`` in deterministic registration order."""
        return self._proposers

    def step(
        self,
        *,
        dt_phys: float,
        solve: Callable[[], None],
        accepted_d: wp.array | None = None,
        ledger: Any | None = None,
        tol_sq_d: wp.array | None = None,
        balance_gamma_n_d: wp.array | None = None,
        rates: object | None = None,
        neighbors: object | None = None,
    ) -> None:
        """Run one whole-cell physical step under a single device acceptance predicate.

        Args:
            dt_phys: the outer physical timestep [s]; must be finite and positive.
            solve: the caller-owned inner mechanical solve (accumulate candidate forces -> converge).
                Runs after events are proposed and before ledgers are assembled; owns all physics.
            accepted_d: the device acceptance predicate.  If ``None``, it is taken from the ledger's
                ``balance_ok_d`` after :meth:`GlobalCellLedger.assemble_balance` — so a caller either
                supplies a composite predicate or lets the force-balance gate decide.
            ledger: optional :class:`GlobalCellLedger`; when given, every participant's ledger terms are
                accumulated and (if ``tol_sq_d`` is given) the force-balance gate is assembled.
            tol_sq_d: optional ``(1,)`` squared-tolerance device scalar for the balance gate.
            balance_gamma_n_d: optional ``(1,)`` float64 device scalar holding the float64 accumulation
                ratio ``γ_n`` from
                :func:`~aleph.engine.ledger.assemble_balance_tolerance_ratio` (PI **D8**).  When
                given, ``tol_sq_d`` is DERIVED on the device from the ledger's own just-accumulated
                resultants — after the accumulation and before the gate — instead of being a value the
                caller carried in.  Deriving it here rather than at the call site is what makes it this
                step's tolerance: the scale ``|reaction| + |traction|`` only exists once the step's own
                contributions have landed, so a host-side derivation could only ever use the previous
                step's magnitudes.
            rates: this step's local rate fields (set by ``CellState``), forwarded to ``propose_events``.
            neighbors: neighbour candidate lists, forwarded to ``propose_events``.

        The predicate stays device-resident throughout; this method never reads it on the host.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if not callable(solve):
            raise TypeError("solve must be a callable inner mechanical solve")
        if tol_sq_d is not None and ledger is None:
            raise ValueError("a balance tolerance requires a ledger to assemble against")
        if balance_gamma_n_d is not None and tol_sq_d is None:
            raise ValueError(
                "balance_gamma_n_d derives the squared tolerance IN PLACE, so tol_sq_d must be supplied "
                "as the device scalar to write it into"
            )

        rng_seed = int(self.clock.base_seed)

        # 1. snapshot every authoritative state as ONE candidate.
        self.world.begin_candidate()

        # 2. propose candidate events into device buffers (rate only; no force, no host state).
        for proposer in self._proposers:
            proposer.propose_events(rates, dt_phys, rng_seed, neighbors)

        # 3. inner mechanical solve to convergence (caller-owned physics).
        solve()

        # 4. accumulate ledgers and (optionally) assemble the force-balance acceptance gate.
        if ledger is not None:
            self.world.accumulate_ledgers(ledger)
            if tol_sq_d is not None:
                if balance_gamma_n_d is not None:
                    ledger.derive_balance_tolerance(balance_gamma_n_d, tol_sq_d)
                ledger.assemble_balance(tol_sq_d)

        predicate = accepted_d
        if predicate is None:
            if ledger is None or getattr(ledger, "balance_ok_d", None) is None:
                raise ValueError(
                    "no acceptance predicate: pass accepted_d, or a ledger with tol_sq_d so the "
                    "force-balance gate can decide"
                )
            predicate = ledger.balance_ok_d

        # 5. finalize: rollback the rejected branch and commit the accepted branch (device selects).
        self.world.finalize_candidate(predicate, dt_phys=dt_phys, rng_seed=rng_seed)

        # 6. advance biological time + RNG epoch, device-gated on the SAME predicate.
        self.clock.advance(predicate, dt_phys)

        # 7. accepted-step no-double-count invariant (host mirror of the device free-lists).
        if self.population_ledgers:
            assert_disjoint_populations(self.population_ledgers)
