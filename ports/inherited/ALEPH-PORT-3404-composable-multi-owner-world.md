# ALEPH-PORT-3404 — a composable multi-owner world, so the engine can be stepped at all

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3404` |
| Lane | `f70de564` scenario lane |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | 30 connectors are wired and **none has ever been stepped inside a world** |

---

## 1. Aleph API

```python
from aleph.scenarios.world import (
    CellWorld,          # a multi-owner world satisfying the relax protocol
    OwnerSlot,          # one owner plus the phases it is scheduled in
    ConnectorSlot,      # one connector plus the resolver that finds its endpoints
    build_cell_world,   # compose a world from whichever owners a scenario needs
)
from aleph.scenarios.audit import render_report   # measured compartment/connector state
```

**`audit.py` reports and computes nothing physical.** It keeps three questions apart — declared,
implemented, exercised — because they had been collapsed into one and the third had never been asked
of the package as a whole. Run it with `python -m aleph.scenarios.audit`.

Its first run found something worth the entry: **`cortex`, `membrane` and `nucleus` declare no state
keys at all**, so nothing checks those three implementations against their registered contracts. The
other eleven do (`owned_state_keys()`, or `STATE_BLOCKS` in `nmii`'s case, which is the same
discipline under a different name and was a false positive until the audit was taught both spellings).

That gap matters more than its count suggests: **`membrane` and `cortex` are the verified vertical's
owners** — the two carrying this project's one verified physics result — and they are among the three
with no mechanical link to the registry. The convention arrived after the three oldest modules did.
Reported here rather than repaired: adding a declaration to a verified owner is that owner's lane.

**`aleph/vertical/assembly.py` is not modified.** That module holds the verified first vertical whose
Laplace and force-gradient results the project rests on, and it carries an **open PI decision** — the
turgor rewire (`HANDOFF.md` §C-④) moves the osmotic envelope from the cortex to the membrane, which
reverses the cortex↔membrane load-path direction and requires Laplace control to be *re-derived*
rather than re-run. Adding ten owners to that file at the same time would entangle two changes whose
failure modes are indistinguishable. This is a separate composition; the vertical keeps working.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` |
| Source path / symbol | **none named, none read** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened. |

`PLAN.md` §7 records that the reference project's `SOLVE_COUPLED` has a declaration and **zero
production tasks** — nothing there has ever converged as a whole cell, and its own documents say so.
There is nothing to port even if policy allowed it. That finding is also the warning: assembling many
owners and getting a number out is not the same as the assembly being right, and the reference
project is the worked example of the difference.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The composition is forced by `aleph/runtime/pipeline.py`,
which already defines `register_participant`, `bind_participant`, `register_connector` and
`bind_connector(name, resolver)`. This module supplies resolvers and a container.

## 4. What this is, and the one thing it must not become

`Pipeline` already schedules participants into phases and binds each connector to a **resolver** that
produces its two live `EndpointHandle`s from the step context. `build_vertical` does this for three
owners and two connectors. What is absent is a composition that does it for **any subset** of the
thirteen owners, so a scenario can ask for exactly the physics it needs.

`CellWorld` satisfies the same duck-typed protocol `relax_to_equilibrium` already consumes —
`transaction()`, `total_potential_energy_pn_um()`, `max_residual_force_pn()`, `evaluate_forces()` —
so the existing monotone-descent relaxation, its rejection logic and its report drive the wide world
**unchanged**. Reusing the verified driver rather than writing a second one is the whole design: a
second descent loop would be a second place for the step-acceptance rule to be right.

**The thing it must not become is a world that reports a plausible number because everything was
scheduled and nothing was evaluated.** `PLAN.md` fixes the standard as scheduled AND reached AND
*evaluated something*, and `Pipeline.coverage_verdict()` already computes it from `witness_delta`. A
`CellWorld` that assembles cleanly and produces an energy while half its connectors witnessed nothing
is precisely the reference project's `SOLVE_COUPLED` failure with better provenance. So:

**Every scenario must assert its coverage verdict, not merely that it ran.** A world is *complete*
only when every REQUIRED participant and connector was reached and witnessed non-zero work.

## 5. Units, domains, singular cases, invariants

Length µm, force pN, energy pN·µm, time s. This module introduces no physical law and therefore no
new unit.

Singular cases, each refused rather than absorbed:

- **An owner appearing twice**, or two owners claiming one `name` — refused. `Pipeline` refuses a
  duplicate registration, and a second owner under one name would silently shadow the first.
- **A connector whose declared endpoints are not both present** — refused with both names. Wiring a
  connector into a world missing one end is the configuration that produces forces attributed to a
  body that is not there.
- **A connector registered but never bound** — refused. `build_vertical` uses exactly this as its
  structural counterexample (`bind_contact=False`), and it must stay a refusal here.
- **An empty world** — refused. Zero owners assembles, steps, conserves energy perfectly and means
  nothing; `aleph/state/full_census.py` records the same trap for an empty manifest.

Invariants, each with the control that asserts it:

- **I1.** Every REQUIRED slot is reached and witnesses non-zero work (`coverage_verdict().ok`).
- **I2.** Force closure across every connector, against the **constituent** scale.
- **I3.** A rejected candidate restores every owner **bit-identically** — checked across the whole
  world, not per owner, because the failure this guards against is one owner rolling back and
  another not.
- **I4.** With one owner and no connectors, `CellWorld` reproduces that owner's own energy exactly.
- **I5.** Adding a connector whose two ends are already in equilibrium changes the total energy by
  exactly `0.0`, so a connector cannot inject energy merely by being present.

## 6. Source evidence class and known retractions

No claim about the reference implementation; none of it was read. Every magnitude a scenario supplies
is `UNSOURCED` and the scenarios say so. Evidence rung `ANALYTIC_ORACLE`, quantitative status
`BLOCKED`.

## 7. Independent oracle or derivation

1. **The single-owner reduction (I4).** A world containing one owner and no couplings must return
   that owner's own energy and forces bit-for-bit. This is an oracle for the composition itself,
   independent of every physical constant.
2. **The verified vertical.** A `CellWorld` built with membrane + cortex + turgor + the two original
   connectors must reproduce `build_vertical`'s relaxed tension to within the tolerance the vertical
   already asserts. **This is the strongest control available**, because the vertical's Laplace
   recovery is the project's one verified physics result — if the composition is wrong, this fails.
3. **Coverage as an oracle for the schedule**, not for the physics: `witness_delta > 0` per slot.

## 8. Positive control

**Measured this session on this tree.** Every number below came from a run in this session.

All controls in `tests/scenarios/test_cell_world.py`. **14 passed in 45.6 s.**

| Control | Asserts | Measured |
|---|---|---|
| `TestTheCompositionIsFaithful::test_it_reproduces_the_verified_vertical_before_any_step` | §7.2 before stepping | energy within **2 ULP**, `max_residual_force_pn` **exactly equal** |
| `TestTheCompositionIsFaithful::test_it_reproduces_the_verified_vertical_through_a_relaxation` | §7.2 through a full descent | see below |
| `TestTheCompositionIsFaithful::test_a_single_owner_world_reproduces_that_owner_s_own_energy` | I4 | energy `==`, forces `np.array_equal` |
| `TestEveryScheduledSlotWitnessedWork::test_the_coverage_verdict_is_complete_after_a_relaxation` | I1 | verdict `ok` |
| `TestTheWholeWorldRollsBackTogether::test_a_rejected_candidate_restores_every_owner_bit_identically` | I3 | `np.array_equal` on every owner |

**The load-bearing measurement**, at subdivision level 2, outside the suite — two independent
verticals, one driven as a `VerticalWorld` and one as a `CellWorld`:

| | converged | accepted | residual (pN) | energy (pN·µm) | σ (pN/µm) | wall |
|---|---|---|---|---|---|---|
| `VerticalWorld` | True | 3012 | 9.998478e-04 | 9301.989722999 | **29.999893750** | 25.4 s |
| `CellWorld` | True | 3012 | 9.998478e-04 | 9301.989722999 | **29.999893750** | 26.6 s |

**Identical to every digit printed**, including the accepted-step count — so the composed world is
not merely reaching the same answer, it is taking the same path there. The recovered tension is the
project's verified Laplace result against a declared 30.0 pN/µm.

**Why the energy is compared within a few ULP and not with `==`.** `CellWorld` sums with `math.fsum`,
which is exact; `VerticalWorld` sums with `+`, which is not. The composed world's energy is therefore
the *more* accurate of the two. At level 2 they agree to every digit; at level 1 they differ by 2
ULP. Asserting `==` would be asserting that fsum and repeated addition round identically, which is
false and says nothing about faithfulness.

## 9. Deliberately failing negative control

| Control | Asserts |
|---|---|
| `TestEveryScheduledSlotWitnessedWork::test_coverage_is_incomplete_before_anything_has_run` | Coverage `ok` is not a constant — it is `False` on a world that has never stepped, so the completeness test above cannot pass vacuously |
| `TestTheNegativeControlsFail::test_an_empty_world_is_refused` | Zero owners assembles, steps, conserves energy perfectly and means nothing |
| `TestTheNegativeControlsFail::test_a_connector_whose_endpoint_is_absent_is_refused` | Forces attributed to a body not in the world |
| `TestTheNegativeControlsFail::test_two_owners_under_one_name_are_refused` | The second shadows the first while every per-owner control still passes |
| `TestTheNegativeControlsFail::test_a_field_whose_load_lands_nowhere_is_refused` | The defect that produced a plausible `max|F|` from a cell with no pressure in it |
| `TestTheNegativeControlsFail::test_a_connector_without_accumulate_is_refused` | A coupling that cannot be scheduled — and one shipped that way and passed the wiring guard |

**Two defects in this module were caught by §7.2 before any scenario was written, and both produced
entirely plausible numbers rather than raising:**

1. Connector energy was read from a cache that is empty before the first step, so every coupling was
   invisible to `total_potential_energy_pn_um` — which is the function the relaxation *descends*.
   The world relaxed a membrane and a cortex that could not feel each other, reporting an energy
   46.7 pN·µm below the vertical's.
2. There was no hook for a field owner's applied load, so the turgor was scheduled, stepped, and
   contributed nothing to the force assembly. `max|F|` read 86.5 pN against the true 65.99 — a
   plausible residual from a cell with no pressure in it.

Neither raised, neither would have been caught by a smoke test, and both are the reason §7.2 compares
against the verified vertical rather than against a fixture.

**A rollback control this file does NOT have**, stated rather than implied: `I3` drives
`snapshot`/`rollback` directly rather than through a *rejected candidate step*, so it establishes
that the owners restore and not that the transaction rejects correctly across many owners. The
per-owner suites already cover the former. The latter needs a world whose candidate is rejected for a
physical reason, which is a scenario's job.

## 10. Numerical and precision envelope

float64. This module performs no arithmetic of its own beyond summing owner energies with
`math.fsum`, which is exact for the summation order. I3 and I4 are asserted with `==` and
`np.array_equal`, not tolerances: a composition that *nearly* reproduces one owner has a defect, and
a rollback that *nearly* restores accumulates over a rejected-step sequence.

`max_steps` is **per configuration and has no default here**. `PLAN.md` §2.5 measured an icosphere
converging in 3,426 steps where a Fibonacci-hull mesh at the same vertex count needed **189,419** —
a factor of 55 from tessellation alone. A shared default would be wrong for most worlds and
plausible for all of them.

## 11. Production-backend residency and transfer

Host, numpy, float64, `NumpyBackend`. **No GPU work was run by this lane and no authorization was
sought or held.** The composition is backend-agnostic by construction: it holds a `Backend` and hands
it to the transaction, exactly as `VerticalWorld` does.

## 12. Comments and docstrings to discard

Nothing to discard; no reference prose entered, because no reference file was read.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **To be filled after the controls run.** Status stays `PROPOSED` until then. |
| Reviewer | Agent-proposed. Unratified. No PI review. |
| Rollback | Delete `aleph/scenarios/` and `tests/scenarios/`. Nothing else imports either; `assembly.py` is untouched and the vertical is unaffected. |

## 14. Honest limits — what this entry does NOT establish

**A world that steps is not a cell that is right.** This entry establishes that many owners can be
scheduled, coupled, stepped and rolled back together. It establishes nothing about whether the
resulting configuration resembles a cell, and the reference project's `SOLVE_COUPLED` is the standing
warning that these are different claims.

**No acceptance predicate is decided.** `ALEPH-DQ-104` is open (`HANDOFF.md` §5.2), and the one
measurement that dominates it is recorded there: on a single fixed mesh, twenty stopping points gave
signed errors spanning **603×** with five sign changes, so a gate on `max|F|` admits states whose
answers differ by three orders of magnitude. Every scenario run under this module inherits that
uncertainty and must report where it stopped, not only what it found.

**The cortex fork is not resolved.** This world uses `ElasticCortexShell`, the status quo, because
whether the filament cortex replaces or coexists with it is an open PI decision. A world cannot hold
two `cortex` owners, so the choice is forced and it is recorded here rather than made.

**No magnitude is sourced.** Every card a scenario supplies is `UNSOURCED`.
