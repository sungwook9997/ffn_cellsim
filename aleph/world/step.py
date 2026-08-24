r"""PHASE 4 — one physical step over the arena, and an HONEST answer about whether it was accepted.

**PHASE 4 is "the step, AND WHAT ACCEPTS IT". This module builds the first half and refuses to fake
the second.**

The port source already runs the step. ``engine/transaction.py``'s ``CellTransaction`` drives the
ratified order — snapshot every authoritative state as one candidate, propose events into device
buffers, run the inner mechanical solve, accumulate ledgers, assemble the balance gate, then launch
BOTH rollback and commit and let the device select the branch, advancing the clock only on accept. That
sequencing is not re-litigated here and this module does not reinvent it.

What this module adds is the thing the port source does not have: **a name for the fact that nobody
knows what "accepted" means yet.**

⚠ **BOTH acceptance predicates this engine has shipped are one-sided, and both were measured.**

* ``engine/ledger.py:239-259`` accepts when ``|reaction + traction|² ≤ tol²``, and ``tol²`` is
  ``assemble_balance_tolerance_ratio(n_terms)`` — the **Higham summation bound**, ``n_terms · eps64``.
  That is the floating-point error of the summation itself. ``sf_motor_slice.py:65`` says so in its own
  words: the predicate tests **adjoint closure**, not convergence. **It cannot fail for physics.**
* The ``descent_ratio`` gate proposed on 2026-08-20 read ``residual_end / residual_start``, and
  ``residual_end`` is a POST-ROLLBACK value pinned to the start on a rejected step. **It cannot
  succeed.** Job 69 was declared VOID for it.

One gate that cannot fail and one that cannot succeed are the same defect twice, from opposite sides:
a criterion defined by the process cannot certify the process.

And the amendment that would replace them, ``STATE.md`` (e) 1 — stationarity instead of convergence —
is still not decided. ⚠ **But the reason changed on 2026-08-21 and this paragraph used to give the old
one.** It said γ **is not emitted on the resting path at all**
(``PI_DECISION_e1_STATIONARITY_2026-08-16.md`` §2.2, *"a stationarity gate cannot be written against an
observable the run does not produce"*). That is no longer true: :mod:`aleph.world.observe_gamma` landed
``984c1292`` and ``world_phase4_native.py --emit-gamma`` produces a full-native γ trace — the τ records
carry 30,000 samples of it. ``GAMMA_EMITTED_e1_IS_NOW_DECIDABLE_2026-08-21.md`` measured what makes it
usable: frozen, γ's spread is **exactly zero** while the residual's is 2.111e-13 pN, and γ responds
monotonically under relaxation. **γ has no noise floor of its own; the residual does.**

**What remains is genuinely the PI's, and it is three things, not one**: WHICH statistic, over WHAT
window, judged against WHAT variance. So the amendment is now *undecided*, not *unwritable* — a
distinction that matters because the second cannot be put to the PI and the first can.

⚠ **And the window input is not ready, measured 2026-08-22 across three seeds.**
``world_phase4/tau3_cut_scan.py``: the equilibration search computes ``opens_on_transient`` **after**
choosing the window and never uses it in the choice, so widening its search bound alone lets seed 3
take the highest reliability of the three (124.8) on a window that opens on a rising transient. Today
every seed is refused and nobody is misled; widening alone would be **more confidently wrong**. Even
with the flag as an in-loop filter, only **1 seed of 3** yields a settled τ — a D-3 input, not D-2.

**So this module makes acceptance a declared, pluggable object with no default, and ships
:class:`UndefinedAcceptance` as the honest state of the art.** A step runs, every ledger is assembled,
the whole transaction sequence executes — and the verdict recorded is ``ACCEPTANCE_UNDEFINED``, which is
a different artifact from ``ACCEPTED`` and cannot be mistaken for one downstream.

⚠ **Why not just use the balance gate and call it accepted.** Because that is exactly how a
force-accepted trajectory gets recorded as an accepted physical trajectory, which the charter names as
a way results get promoted past what produced them. The gate would say yes to almost anything; writing
``ACCEPTED`` next to it would launder the Higham bound into a physics verdict. Every native run this
engine has recorded carries that laundering already, and this module is where it stops.

**The interface this reserves.** ``CellStatePosterior.accepted_step_evidence`` is typed
``AdmissibilityVerdict | None`` (``aleph/virtual_cell/cell_state_posterior.py``). It is ``None`` today
because nothing can fill it. :class:`StepVerdict` is the shape that will, and the day (e) 1 resolves,
the only change is which :class:`AcceptancePredicate` is passed in — not a rewrite of the loop.

Sanity Gate (recorded before first execution, per the charter):
    * **dimensions** — ``dt_phys`` is seconds and must be finite and positive; residuals are pN; the
      arena's own units are µm·pN·s throughout and no conversion happens here.
    * **boundary cases** — zero participants, zero events, and a solve that raises are each tested;
      a predicate that is undefined must not yield ``ACCEPTED`` under any of them.
    * **conservation invariants** — population ID ranges stay disjoint across a step
      (``arena.assert_partitioned``); the step asserts it AFTER the transaction, not only at build.
    * **CFL / precision** — not this module's: the inner solve owns its own step-size rule and this
      orchestrator neither picks nor checks one. Stated so the omission is deliberate.
    * **sign sense** — accept is ``1`` and reject is ``0`` in the device predicate, matching
      ``ledger.py``'s ``balance_ok``; ``UndefinedAcceptance`` writes NEITHER and reports that it did not.
    * **measurement protocol** — one ``StepVerdict`` per step, appended in order; the verdict carries
      the predicate's own name so an artifact can never be read without knowing what judged it.

engine units: length µm, force pN, time s. Runtime: NVIDIA Warp on CUDA. This module is CPU-importable
on purpose — it owns no physics and launches no kernel, so its gates can run on the dev machine, which
is the same reason ``engine/transaction.py`` is importable there.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Protocol, runtime_checkable

__all__ = [
    "ACCEPTANCE_UNDEFINED",
    "AcceptancePredicate",
    "AdjointClosurePredicate",
    "StepVerdict",
    "UndefinedAcceptance",
    "WorldStep",
]

#: The verdict a step carries when no acceptance criterion is defined for it. **Not a failure and not
#: an acceptance** — the run happened and nothing is entitled to say whether it counts.
ACCEPTANCE_UNDEFINED = "ACCEPTANCE_UNDEFINED"

#: What a predicate that HAS decided may report. Anything outside this set is refused at construction,
#: so a predicate cannot invent a third state that reads as approval.
_DECIDED = ("ACCEPTED", "REJECTED")


@runtime_checkable
class AcceptancePredicate(Protocol):
    """What decides whether a step counts.

    Deliberately tiny, and deliberately carrying its own ``name`` and ``tests``: an artifact that
    records a verdict without recording what produced it is how ``balance_ok`` came to be read as
    convergence for months. The predicate says what it examines, in its own words, and that string
    travels into every record.
    """

    #: Stable identifier, recorded in every verdict.
    name: str
    #: One sentence naming what this predicate ACTUALLY tests — not what a reader might assume.
    tests: str

    def decide(self, ctx: "StepContext") -> str:
        """Return one of ``_DECIDED``, or :data:`ACCEPTANCE_UNDEFINED`."""
        ...


@dataclass(frozen=True, slots=True)
class StepContext:
    """Everything a predicate is allowed to look at, and nothing else.

    Attributes:
        step_index: 0-based index of this step within the run.
        dt_phys_s: the outer physical timestep [s].
        residual_pn: the assembled force-balance residual [pN], or ``None`` when no ledger was given.
        tol_pn: the tolerance the residual was compared against [pN], or ``None``.
        n_events: candidate kinetic events proposed this step.
        observables: whatever the run chose to emit, by name. A predicate that needs an observable the
            run does not produce must say so rather than substitute one — that is the whole of (e) 1.
    """

    step_index: int
    dt_phys_s: float
    residual_pn: float | None
    tol_pn: float | None
    n_events: int
    observables: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class StepVerdict:
    """One step's outcome, with the thing that judged it attached.

    ``predicate_name`` and ``predicate_tests`` are not decoration. A verdict read without them is how a
    rounding bound gets quoted as convergence, so they are structural fields and
    :meth:`as_record` refuses to omit them.
    """

    step_index: int
    verdict: str
    predicate_name: str
    predicate_tests: str
    dt_phys_s: float
    residual_pn: float | None
    tol_pn: float | None
    n_events: int
    wall_s: float

    @property
    def is_accepted(self) -> bool:
        """True only for a literal ``ACCEPTED``.

        ⚠ ``ACCEPTANCE_UNDEFINED`` is falsy here ON PURPOSE. Downstream code that asks "was this
        accepted" must get ``False`` from a step nobody judged, because the alternative — treating an
        unjudged step as accepted — is the exact promotion the charter forbids.
        """
        return self.verdict == "ACCEPTED"

    @property
    def is_judged(self) -> bool:
        """Whether any criterion was applied at all. Distinct from :attr:`is_accepted`."""
        return self.verdict in _DECIDED

    def as_record(self) -> dict[str, object]:
        """The artifact row. Carries the judge, always."""
        return {
            "step_index": self.step_index,
            "verdict": self.verdict,
            "predicate": self.predicate_name,
            "predicate_tests": self.predicate_tests,
            "dt_phys_s": self.dt_phys_s,
            "residual_pn": self.residual_pn,
            "tol_pn": self.tol_pn,
            "n_events": self.n_events,
            "wall_s": round(self.wall_s, 6),
            "quantitative_claim_status": "BLOCKED" if not self.is_judged else "SEE_PREDICATE",
        }


@dataclass(frozen=True, slots=True)
class UndefinedAcceptance:
    """The honest default: run the step, judge nothing, and say that is what happened.

    ⚠ **This is the current state of the art and it is not a placeholder to be swapped for
    convenience.** `STATE.md` (e) 1 is still undecided — but ⚠ **for a different reason than this
    docstring gave until 2026-08-22.** It said γ "is not emitted on the resting path"; it is, since
    `984c1292`. The condition this class was written against was *"until γ is emitted AND the PI fixes
    the statistic, its window and the variance it is judged against"* — the **first conjunct is now
    satisfied and the other three are not**, so the class is unchanged and its reason is narrower.

    A predicate here would still only be asserting that the solver stopped, and one more thing is
    measured against writing it today: across three seeds only one yields a settled τ
    (`world_phase4/tau3_cut_scan.py`, 2026-08-22), so even the window has no D-2 input yet.
    """

    name: str = "undefined"
    tests: str = (
        "NOTHING. STATE.md (e) 1 is UNDECIDED, and as of 2026-08-21 that is no longer the same as "
        "UNWRITABLE: cortical tension is now emitted on the resting path (world/observe_gamma.py, "
        "984c1292), it has no noise floor of its own where the residual has 2.111e-13 pN, and it "
        "responds monotonically. What is still open is the PI's and it is three things: WHICH "
        "statistic, over WHAT window, judged against WHAT variance. Measured 2026-08-22 across three "
        "seeds: only 1 of 3 yields a settled tau, so the window has a D-3 input and not a D-2 one. "
        "This predicate exists so a step can run without a verdict being invented for it."
    )

    def decide(self, ctx: StepContext) -> str:
        """Always :data:`ACCEPTANCE_UNDEFINED`, regardless of what the step did."""
        return ACCEPTANCE_UNDEFINED


@dataclass(frozen=True, slots=True)
class AdjointClosurePredicate:
    """The engine's existing balance gate, **named for what it actually tests**.

    ``ledger.py:239-259`` accepts on ``|reaction + traction|² ≤ tol²`` where ``tol²`` is the Higham
    summation bound ``n_terms · eps64``. That is the rounding error of the summation, so passing it
    means the adjoint wiring closes — every force appears with both signs — and nothing more.

    ⚠ **Provided so it can be used HONESTLY, not so it can be used as convergence.** Its ``tests``
    string travels into every verdict it produces, and the sentence is deliberately blunt. If a future
    run wants a convergence gate, this is not it; ``sf_motor_slice.py:65`` said so first.
    """

    name: str = "adjoint_closure"
    tests: str = (
        "ADJOINT CLOSURE ONLY — that every force appears with both signs, to within the Higham "
        "summation bound n_terms*eps64. NOT convergence: this predicate cannot fail for physics, and "
        "a run that passes it has not been shown to have converged to anything."
    )

    def decide(self, ctx: StepContext) -> str:
        """``ACCEPTED`` when the residual is inside the summation bound."""
        if ctx.residual_pn is None or ctx.tol_pn is None:
            return ACCEPTANCE_UNDEFINED
        return "ACCEPTED" if ctx.residual_pn <= ctx.tol_pn else "REJECTED"


@dataclass(slots=True)
class WorldStep:
    """One physical step over the arena, under a declared acceptance predicate.

    Args:
        arena: the world whose populations take the step. Its partition invariant is re-asserted after
            every step, not only at build — ``assert_partitioned`` is cheap and a step that corrupts
            ranges should be caught by the step that did it.
        predicate: what decides whether a step counts. **No default.** A caller either states the
            criterion or passes :class:`UndefinedAcceptance` and gets an unjudged run; there is no
            path where one is chosen silently.
        on_verdict: optional callback per step, for a driver that wants to stream rather than collect.

    Raises:
        TypeError: if ``predicate`` does not satisfy :class:`AcceptancePredicate`, or its ``tests``
            string is empty — an unexplained judge is refused at construction.
    """

    arena: Any
    predicate: Any
    on_verdict: Callable[[StepVerdict], None] | None = None
    verdicts: list[StepVerdict] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.predicate, AcceptancePredicate):
            raise TypeError(
                f"{type(self.predicate).__name__} is not an AcceptancePredicate: it must carry `name`, "
                "`tests` and `decide(ctx)`. `tests` is required because a verdict recorded without "
                "what produced it is how a rounding bound came to be quoted as convergence."
            )
        if not str(getattr(self.predicate, "tests", "")).strip():
            raise ValueError(
                f"predicate {self.predicate.name!r} has an empty `tests` string. Say what it examines, "
                "in its own words; the string travels into every artifact it judges."
            )

    def step(
        self,
        *,
        dt_phys_s: float,
        solve: Callable[[], None],
        propose_events: Callable[[], int] | None = None,
        assemble_residual: Callable[[], tuple[float, float]] | None = None,
        observables: Callable[[], dict[str, float]] | None = None,
    ) -> StepVerdict:
        """Run one step in the ratified order and return its verdict.

        The order is ``engine/transaction.py``'s and is not re-derived: propose candidate events, run
        the inner mechanical solve, assemble the residual, judge, then re-assert the partition.

        Args:
            dt_phys_s: outer physical timestep [s]. Finite and positive.
            solve: the caller-owned inner mechanical solve. Owns all physics; this method owns none.
            propose_events: optional; returns how many candidate kinetic events were proposed.
            assemble_residual: optional; returns ``(residual_pn, tol_pn)``. Absent, the predicate is
                handed ``None`` for both and most predicates will answer
                :data:`ACCEPTANCE_UNDEFINED` — which is correct, not a degradation.
            observables: optional; whatever the run emits, for a predicate that reads state rather than
                residual. **This is the hook (e) 1 needs**: when γ is emitted on the resting path it
                arrives here, and a stationarity predicate becomes writable without touching this loop.

        Returns:
            The :class:`StepVerdict`, also appended to :attr:`verdicts`.

        Raises:
            ValueError: on a non-finite or non-positive ``dt_phys_s``.
            TypeError: if ``solve`` is not callable.
        """
        import time

        if not math.isfinite(dt_phys_s) or dt_phys_s <= 0.0:
            raise ValueError(f"dt_phys_s must be finite and positive; got {dt_phys_s!r}")
        if not callable(solve):
            raise TypeError("solve must be a callable inner mechanical solve")

        t0 = time.perf_counter()
        n_events = int(propose_events()) if propose_events is not None else 0
        solve()
        residual = tol = None
        if assemble_residual is not None:
            residual, tol = (float(x) for x in assemble_residual())

        ctx = StepContext(
            step_index=len(self.verdicts), dt_phys_s=float(dt_phys_s),
            residual_pn=residual, tol_pn=tol, n_events=n_events,
            observables=dict(observables()) if observables is not None else {},
        )
        decided = self.predicate.decide(ctx)
        if decided not in _DECIDED and decided != ACCEPTANCE_UNDEFINED:
            raise ValueError(
                f"predicate {self.predicate.name!r} returned {decided!r}, which is neither a decision "
                f"{_DECIDED} nor {ACCEPTANCE_UNDEFINED}. A third state would be read as approval by "
                "something downstream; there is no third state."
            )

        verdict = StepVerdict(
            step_index=ctx.step_index, verdict=decided,
            predicate_name=self.predicate.name, predicate_tests=self.predicate.tests,
            dt_phys_s=ctx.dt_phys_s, residual_pn=residual, tol_pn=tol, n_events=n_events,
            wall_s=time.perf_counter() - t0,
        )
        self.verdicts.append(verdict)
        if self.on_verdict is not None:
            self.on_verdict(verdict)

        # The no-double-count invariant, promoted from build time to every step. Cheap, and a step that
        # corrupts an ID range should be caught by the step that did it rather than at the next census.
        if hasattr(self.arena, "assert_partitioned"):
            self.arena.assert_partitioned()
        return verdict

    def as_record(self) -> dict[str, object]:
        """The run's step block, for a driver's artifact."""
        judged = [v for v in self.verdicts if v.is_judged]
        return {
            "schema": "ffn-world-step@1",
            "predicate": self.predicate.name,
            "predicate_tests": self.predicate.tests,
            "n_steps": len(self.verdicts),
            "n_judged": len(judged),
            "n_accepted": sum(1 for v in self.verdicts if v.is_accepted),
            "n_unjudged": len(self.verdicts) - len(judged),
            "steps": [v.as_record() for v in self.verdicts],
            # ⚠ Zero steps is NOT "every step judged". `len(judged) == len(verdicts)` is 0 == 0 for an
            # empty run, so the first version of this line reported evidence available for a run that
            # took no step at all — a vacuous pass, and the exact class this module exists to refuse.
            # Found on a fresh read 2026-08-21; the driver refuses a zero-step run, so it never
            # reached an artifact, but as_record() is the artifact interface and a future caller
            # would have.
            "accepted_step_evidence": (
                None if not self.verdicts or len(judged) != len(self.verdicts) else "SEE_PREDICATE"),
            "accepted_step_evidence_note": (
                "None means at least one step carries no criterion. CellStatePosterior."
                "accepted_step_evidence is typed AdmissibilityVerdict | None and stays None until "
                "every step in a run has been judged — an unjudged step may not be counted as "
                "admissible evidence, which is the whole reason this field is nullable. A run of ZERO "
                "steps is also None: nothing was judged because nothing happened, and 'no unjudged "
                "steps' is not the same claim as 'every step was judged'."),
        }


def _demo() -> None:
    """Sanity Gate execution — every boundary case in the module docstring, and the two that bite."""
    class _Arena:
        def __init__(self) -> None:
            self.checks = 0

        def assert_partitioned(self) -> None:
            self.checks += 1

    arena = _Arena()

    # 1. An unjudged step is NOT accepted, and says which is which.
    ws = WorldStep(arena=arena, predicate=UndefinedAcceptance())
    v = ws.step(dt_phys_s=0.05, solve=lambda: None)
    assert v.verdict == ACCEPTANCE_UNDEFINED
    assert not v.is_accepted and not v.is_judged
    assert arena.checks == 1, "the partition invariant runs every step, not only at build"

    # 2. Zero participants, zero events, and a residual that is absent are all fine and all unjudged.
    assert v.n_events == 0 and v.residual_pn is None

    # 3. The balance gate, used honestly: it decides, and its verdict carries what it tested.
    ws2 = WorldStep(arena=arena, predicate=AdjointClosurePredicate())
    a = ws2.step(dt_phys_s=0.05, solve=lambda: None, assemble_residual=lambda: (1e-9, 1e-8))
    r = ws2.step(dt_phys_s=0.05, solve=lambda: None, assemble_residual=lambda: (1.0, 1e-8))
    assert a.verdict == "ACCEPTED" and r.verdict == "REJECTED"
    assert "NOT convergence" in a.predicate_tests, "the caveat travels with the verdict"
    # ...and with no residual to look at, it declines rather than guessing.
    u = ws2.step(dt_phys_s=0.05, solve=lambda: None)
    assert u.verdict == ACCEPTANCE_UNDEFINED

    # 4. accepted_step_evidence stays None while ANY step is unjudged — this is the interface to
    #    CellStatePosterior and the reason that field is nullable.
    assert ws2.as_record()["accepted_step_evidence"] is None, "one unjudged step poisons the run"
    assert ws.as_record()["n_accepted"] == 0

    # 5. A predicate that invents a third state is refused. "MAYBE" reads as approval downstream.
    class _Sneaky:
        name = "sneaky"
        tests = "returns a third state"

        def decide(self, ctx: StepContext) -> str:
            return "MAYBE"

    try:
        WorldStep(arena=arena, predicate=_Sneaky()).step(dt_phys_s=0.05, solve=lambda: None)
    except ValueError as exc:
        assert "third state" in str(exc)
    else:
        raise AssertionError("a third verdict state must be refused")

    # 6. A judge with no explanation is refused at construction.
    class _Mute:
        name = "mute"
        tests = "   "

        def decide(self, ctx: StepContext) -> str:
            return "ACCEPTED"

    for bad, exc_type in ((_Mute(), ValueError), (object(), TypeError)):
        try:
            WorldStep(arena=arena, predicate=bad)
        except exc_type:
            pass
        else:
            raise AssertionError(f"{type(bad).__name__} must be refused")

    # 7. dt and solve are checked; a solve that raises is not swallowed into a verdict.
    for kwargs in ({"dt_phys_s": 0.0}, {"dt_phys_s": float("nan")}, {"dt_phys_s": -1.0}):
        try:
            ws.step(solve=lambda: None, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{kwargs} must be refused")

    def _boom() -> None:
        raise RuntimeError("the inner solve failed")

    n_before = len(ws.verdicts)
    try:
        ws.step(dt_phys_s=0.05, solve=_boom)
    except RuntimeError:
        pass
    else:
        raise AssertionError("a failing solve must propagate, not become a verdict")
    assert len(ws.verdicts) == n_before, "a failed step records no verdict at all"

    # 8. The observables hook — the seat (e) 1's stationarity gate will take without touching this loop.
    seen: dict[str, float] = {}

    class _NeedsGamma:
        name = "needs_gamma"
        tests = "would test stationarity of cortical tension; reports UNDEFINED when it is absent"

        def decide(self, ctx: StepContext) -> str:
            seen.update(ctx.observables)
            return ACCEPTANCE_UNDEFINED if "gamma_pn_per_um" not in ctx.observables else "ACCEPTED"

    ws3 = WorldStep(arena=arena, predicate=_NeedsGamma())
    assert ws3.step(dt_phys_s=0.05, solve=lambda: None).verdict == ACCEPTANCE_UNDEFINED
    got = ws3.step(dt_phys_s=0.05, solve=lambda: None,
                   observables=lambda: {"gamma_pn_per_um": 3.7})
    assert got.verdict == "ACCEPTED" and seen["gamma_pn_per_um"] == 3.7

    print(f"world.step self-check OK — {len(ws.verdicts) + len(ws2.verdicts) + len(ws3.verdicts)} "
          "verdicts, none of them invented")


if __name__ == "__main__":
    _demo()
