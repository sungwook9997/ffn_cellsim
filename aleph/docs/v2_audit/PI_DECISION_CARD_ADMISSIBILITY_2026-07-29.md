# PI Decision Card — accepted-step admissibility: conjoin convergence with adjoint closure (2026-07-29)

**Status: PROPOSED. This card changes nothing.** No engine file was edited. No gate threshold was
changed. No `ValidationGate` or `ModelContract` row was created (those are PI-authored). Lane A owns
exactly three files — `aleph/virtual_cell/admissibility.py`,
`aleph/tests/ac/virtual_cell/test_admissibility.py`, and this document — and has written to nothing
under `aleph/ac/`. Every number below was measured on **CPU**, in a sandbox oracle that imports no
engine runtime and no Warp. **A CPU number is not a physics result and none is claimed here**: what is
claimed is a statement about a *predicate's logic*, which is checkable without a device.

Repo state at authoring: branch `codex/ff-ac-codex`, HEAD `7daedb3d`.

---

## 1. The defect

**The canonical accepted-step predicate does not conjoin `inner_converged` with `balance_ok`.
`balance_ok` is an adjoint FORCE-CLOSURE check, not mechanical convergence — so the outer KMC/clock
can commit on top of an inner mechanical solve that did not converge.**

### 1a. The acceptance path, as written

`aleph/engine/transaction.py:154-161` — `CellTransaction.step` takes ONE predicate and, when the
caller omits `accepted_d`, falls back to the ledger's force-balance flag:

```python
predicate = accepted_d
if predicate is None:
    if ledger is None or getattr(ledger, "balance_ok_d", None) is None:
        raise ValueError(
            "no acceptance predicate: pass accepted_d, or a ledger with tol_sq_d so the "
            "force-balance gate can decide"
        )
    predicate = ledger.balance_ok_d
```

That same flag is then handed to *both* branches of `finalize_candidate` (`transaction.py:164`) and to
the clock (`transaction.py:167`). The inner mechanical solve is called at `transaction.py:143-145`:

```python
# 3. inner mechanical solve to convergence (caller-owned physics).
solve()
```

`solve` returns `None` and its convergence state is never consulted. There is no code path in
`CellTransaction` by which a non-converged inner solve can veto the commit.

### 1b. What `balance_ok` actually tests — the engine says it itself

`aleph/engine/ledger.py:242-258` (`_assemble_force_balance_kernel`) computes
`|Σreaction + Σtraction|²` and compares it to `tol_sq`. `ledger.py:412-421`
(`GlobalCellLedger.add_body_force`) states the invariant in its own docstring:

> "every force is either internal to one array — where the pairs cancel within that array's own
> resultant — or an adjoint pair scattered across the two, so the two resultants are equal and
> opposite **whatever the configuration and however unconverged the solve**. Which is the point: the
> gate then tests the **ADJOINT WIRING** (Newton's third law across the connector), **not the
> convergence**, and it fails exactly when a scatter is one-sided, double counted, or sign-flipped."

`aleph/engine/sf_motor_slice.py:65` repeats it in that slice's Sanity Gate:

> "The gate tests the adjoint wiring, not the convergence, and the driver's positive control shows it
> rejecting a deliberately unbalanced step."

`STATE.md` tier-(a) carries the same reading as a *quotable claim about the lane*:

> "**The predicate tests ADJOINT CLOSURE, not convergence** (`sf_motor_slice.py`'s own docstring says
> so), so a non-converged step is not what it rejects."

So this is not a discovered inconsistency between code and doc — code, module docstring and status
register all agree on what the predicate does. What is missing is the **conjunction**, and that
absence has not been surfaced as a decision.

### 1c. The two live paths use two DISJOINT halves of the conjunction

This is the sharpest form of the defect.

| Path | Acceptance predicate | Adjoint closure checked? | Convergence checked? |
|---|---|---|---|
| `ac/engine` `CellTransaction` (`transaction.py:154-161`) | `ledger.balance_ok_d` | **yes** | **no** |
| `aleph/components/incumbent/driver.py` (incumbent) — `outer_accepted = inner_converged` | `inner_d.converged_d` | **no** | **yes** |

`aleph/components/incumbent/driver.py` (`inner_solve.commit_irreversible(inner_d.converged_d)`, then
`outer_accepted = inner_converged` / `outer_rolled_back = not inner_converged` /
`committed_time = dt_phys if inner_converged else 0.0`) commits on convergence alone and never
assembles the balance gate. So a one-sided connector scatter — precisely what the SF lane's positive
control demonstrates the balance gate catching — is invisible on the driver path, and a stalled
mechanical solve is invisible on the engine path. **Neither path implements the conjunction, and the
two are not comparable to each other.**

Note that the driver *already* forms the conjunction, but only as a **post-hoc report** field, not as
the transaction predicate:

```python
stable = bool(
    pos_finite and force_finite and inner_finite and fluid_finite and inner_converged
    and outer_accepted and not outer_rolled_back ...
```

The clauses this card proposes therefore already exist as *reported* quantities; what is proposed is
promoting them from a report line to the gate.

### 1d. The cortex path in its current form

`aleph/scripts/ac_gate_b_cortex_motor_native.py:50-52` (module docstring):

> "FORCE-ACCEPT: the demo force-accepts every step (`accepted_d = ones`) so binding accumulates and
> the mechanism is visible; a resting-baseline PRODUCTION run would instead gate on the force-balance
> predicate (reject a non-converged candidate) once γ is high enough to hold turgor. Pass
> `--balance-gate` to switch."

`:696` — `slc.step(cortex_relax, dt_phys=float(args.dt), accepted_d=accepted_ones)  # force-accept`
`:794` — the run record self-stamps
`"acceptance": "force-accept (accepted_d = ones); --balance-gate is not implemented"`.

Two things follow. First, on that lane today acceptance is not a predicate at all. Second, the
docstring's own plan is inaccurate in one word: switching to the force-balance predicate would **not**
"reject a non-converged candidate" — per §1b it would reject a *miswired* candidate. The sentence
describes a gate the codebase does not have. That single word is what this card exists to resolve.

---

## 2. The counterexample — "balanced-but-stretched"

A step in which every connector's force pair cancels to floating-point round-off while the inner
mechanical residual sits six decades above the tolerance the solve was launched with.

**Measured** by `aleph/tests/ac/virtual_cell/test_admissibility.py`
(`test_control_1_balanced_but_stretched_is_accepted_by_adjoint_only_and_rejected_here`,
`test_control_1_numbers_are_the_ones_the_decision_card_quotes`), CPU, no GPU, no engine import:

| Quantity | Value |
|---|---|
| connector `nmii_sf_motor` — `\|Σf_a + Σf_b\|` | `2.273737e-13` pN |
| connector `nmii_sf_motor` — force traffic `Σ\|f_a,i\| + Σ\|f_b,i\|` | `2818.978155` pN |
| → relative closure | `8.066e-17` |
| connector `erm_cortex` — residual / traffic | `4.440892e-16` pN / `7.842194` pN → `5.663e-17` |
| adjoint-closure clause: worst-excess measured vs threshold | `4.440892e-16` pN vs `7.843194e-09` pN → **PASS** |
| inner mechanical residual | `2634.12` pN |
| inner tolerance the solve was judged against | `1.0e-3` pN |
| → residual / tolerance | `2.63412e+06` (**6.4 decades over**) |
| inner iterations / budget | `400 / 400` (budget exhausted) |

**Verdicts on the same record:**

- incumbent adjoint-only predicate (`adjoint_closure_only`, a faithful CPU model of
  `assemble_balance`): **ACCEPT**
- proposed conjunctive predicate (`evaluate_admissibility`): **REJECT**, `rejected_by =
  ('inner_converged',)`, expansion route `native cost failure -> solver / preconditioner / adaptive
  design`

**Provenance of the residual magnitude.** `2634.12` pN is not invented and is not re-measured here. It
is the value the archived full-native run recorded (`STATE.md` tier-(a): *"they also record
`inner_converged: false`, `residual_start = residual_end = 2634.12`"*), used as a **fixture magnitude**
so the counterexample sits at a scale the programme has actually seen rather than a scale chosen to
look bad. The `1.0e-3` pN tolerance is likewise a fixture, not a ratified gate value.

### 2b. Three further controls, all measured

| Control | Incumbent (adjoint-only) | Proposed (conjunctive) |
|---|---|---|
| **broken adjoint** — one connector scatters `+f` on A and `−f + [12, 0, 0]` pN on B (12 pN of the pair never lands) | REJECT (correct — this is the SF lane's positive control) | REJECT, `rejected_by = ('adjoint_closure',)`, measured residual `12.0` pN |
| **cancelling per-connector errors** — two connectors leak `+12` pN and `−12` pN | **ACCEPT** (the global resultant is exactly zero) | REJECT, `adjoint_closure` |
| **nonfinite** — `inner_residual = NaN`, forces clean | ACCEPT | REJECT, `state_finite` first; `adjoint_closure` / `inner_converged` / `iteration_budget` all `NOT_EVALUATED` with `measured = None` |
| **nonfinite state** — NaN or ±inf in `candidate_state.pos_um`, forces clean | **ACCEPT** (a NaN that never reaches the force resultant is invisible to the balance gate) | REJECT, `state_finite`, and the failing field is named |

The NaN trap, stated as the executable fact the test asserts: `float("nan") > 1.0e-3` is `False` **and**
`float("nan") <= 1.0e-3` is `False`. A predicate written as "reject if `residual > tol`" therefore
**accepts** NaN. This is why finiteness is clause zero and short-circuits: on a nonfinite record no
ordered comparison is performed at all, and the verdict proves it by reporting the downstream clauses
as `NOT_EVALUATED` with `measured = None` rather than as passes.

### 2c. Bit-exact rollback

On reject, all four channels — state, RNG, topology, clock — must be restored **bit-exactly**,
verified by SHA-256 over `ndarray.tobytes()`, never by float comparison. Measured negative controls,
one per channel (parametrized): a rollback that restores state/RNG/topology but leaves the **clock
advanced**; one that **drops a topology edge**; one that leaves the **RNG stream advanced**; one that
**misses a state array**. Each is detected, `bit_exact` is `False`, `transaction_sound` is `False`,
and the verdict names the leak — e.g. `"rollback leaked 1/4 channel(s): clock"`. **Restoring 3 of 4 is
a failure**, and the audit says which one.

Why bytes and not floats, as measured: `+0.0 == -0.0` is `True` (a float compare would call a
sign-flipped restore "restored") and `nan == nan` is `False` (a float compare would call a
byte-identical restore "leaked"). The digest gets both right. The audit additionally flags a
**vacuous** check — a rollback verified over a candidate that mutated nothing proves nothing.

### 2d. Scale safety of the closure tolerance — why an absolute pN threshold is wrong

Measured, two connectors whose force traffic differs by nine decades:

| Connector | residual | traffic `Σ\|f_i\|` | relative |
|---|---|---|---|
| heavy | `1.0e-4` pN | `2.0e+6` pN | `5.0e-11` |
| light | `1.0e-5` pN | `1.99e-3` pN | `5.025e-3` |

A single absolute threshold of `10^-4.5 = 3.162e-5` pN **passes the light pair and fails the heavy
one** — exactly backwards. The relative measure with an absolute floor,
`|Σf_a + Σf_b| ≤ rel_tol · (Σ|f_a,i| + Σ|f_b,i|) + abs_floor`, gets both right. The scale is the sum of
contribution *magnitudes*, not the magnitude of their sum, mirroring the engine's own D8 derivation
(`ledger.py:215-235`): a bipolar minifilament is a force dipole whose resultant is near zero while its
traffic is not, so bracketing by the resultant would put the tolerance orders of magnitude below the
arithmetic error it exists to admit. The absolute floor covers the genuinely idle connector, where a
relative test degenerates.

**Neither tolerance has a default in the oracle.** Both are required keyword arguments and the call
raises without them (asserted by test). Choosing their values is a **PI gate decision** and is *not*
part of this card's proposal.

---

## 3. What is PROPOSED

**P1 (primary).** Change the canonical accepted-step predicate from adjoint closure alone to the
**conjunction**, evaluated in this order with a hard short-circuit on the first clause:

```
admissible  ==  state_finite
            AND adjoint_closure   (per connector, scale-safe)
            AND inner_converged   (residual ≤ the tolerance THIS solve was launched with)
            AND iteration_budget  (iterations ≤ budget)
```

with the transaction additionally required to restore all four channels bit-exactly on reject
(`transaction_sound`).

**P2 (separable).** Judge adjoint closure **per connector** rather than on the global resultant. §2b
row 2 shows the global form accepting two connectors whose one-sided errors cancel. P2 can be adopted
or deferred independently of P1.

**P3 (companion, required if P1 is adopted).** An accepted-step predicate that can reject needs a
**dt controller**. Today a rejected step rolls back and the caller re-attempts at the same `dt_phys`
with the same state — which, for a deterministic solve, reproduces the same rejection forever. P1
without P3 converts a silent wrong commit into a live-lock. This is a design task, not a threshold.

**P4 (sequencing — recommended first move).** Adopt P1 in **shadow mode** before adopting it as the
gate: compute all four clauses every step, **keep committing on the incumbent predicate**, and record
both verdicts into the run record. This measures the reject fraction on real runs at zero behavioural
risk and turns §5's "cannot measure" into a number. Only after that number exists should P1 flip the
commit branch.

### What is explicitly NOT proposed

- **No threshold value.** `adjoint_rel_tol`, `adjoint_abs_floor_pn` and the inner tolerance are PI
  decisions and are not proposed here. Nothing in this card was chosen to make anything pass.
- **No change to the D8 tolerance derivation** — the scale-safe form above is the same principle.
- **No new gate row.** `ValidationGate` / `ModelContract` rows are PI-authored.
- **No claim that any existing result is wrong.** The claims `STATE.md` records for these lanes are
  already correctly hedged ("the predicate tests adjoint closure, not convergence"; "quote that it
  builds and steps — not that the step converged"). This card does not retract anything.

---

## 4. What would break if adopted — honest accounting

**Conjoining convergence will reject steps that currently commit.** That is the point, and it has a
cost that must be stated before, not after, the change.

1. **The GATE-B cortex-motor lane would stop committing entirely, by construction.** It force-accepts
   every step (`accepted_d = ones`), and its own record says `--balance-gate is not implemented`. Under
   P1 its reject fraction is whatever fraction of its steps fail to converge, and it currently has no
   code path that measures that. This is a fact about the script, cited from the script; it is not a
   measurement of a production run.
2. **The archived full-native run would commit zero steps.** `STATE.md` tier-(a) records
   `inner_converged: false` and `residual_start = residual_end = 2634.12` for the full-native
   build-and-step run — a flat residual, i.e. the solve was not progressing. Under P1 every step of
   that run is inadmissible and the outer clock advances by zero. **That is the correct verdict** — the
   run's own artifact already stamps `overall_status: INCOMPLETE_OR_FAIL` — but it means P1 does not
   make that configuration usable, it makes its unusability load-bearing.
3. **The honest worst case: the outer clock stalls and the cell makes no biological time.** If the
   native solve cannot reach tolerance, P1 turns a run that produced a trajectory into a run that
   produces nothing. That is a real loss of apparent output. It is not a reason to weaken the
   predicate: a trajectory assembled from non-converged steps was never a trajectory. It routes to
   **`native cost failure -> solver / preconditioner / adaptive design`** — the solver, the
   preconditioner, the timestep controller and the tolerance derivation are what must move.

### What I cannot measure, stated plainly

**I cannot quantify what fraction of a production run P1 would reject.** That number requires
instrumenting `CellTransaction` and running on the A5000, and this lane is CPU-only by mandate and
owns no engine file. Every number in §2 is a property of the *predicate*, measured on synthetic
records; **none of them is a reject rate.** The two figures in items 1–2 above are cited from other
lanes' recorded artifacts, not measured by me, and neither is a production reject rate:

- item 1 is 100% *by construction* (the script never evaluates a predicate), which is a statement about
  a demo script, not about the physics;
- item 2 is 100% for **one** archived run whose own status is `INCOMPLETE_OR_FAIL`.

Generalising either to "P1 would reject N% of production steps" would be exactly the kind of
unmeasured magnitude the programme has already had to retract. **P4 is the way to get the real
number.**

### Second-order effects to check before adoption

- Rejecting more steps changes the **KMC event statistics**: proposed events on rejected steps must be
  discarded and the RNG stream restored, or the accepted-step event distribution is biased. §2c's
  four-channel bit-exact requirement is the check for this, and the RNG channel is one of the four for
  exactly this reason.
- `iteration_budget` as a hard clause interacts with retry logic (`inner_max_attempts`): budget must be
  defined as *total across retries*, or the clause is trivially satisfiable.
- Per-connector closure (P2) requires each connector to expose **its own** two scattered force
  contributions, not a global reduction. Some connectors may not currently do so — that is a scoping
  question P2 must answer before it can be adopted.

---

## 5. Alternatives considered

| # | Option | Verdict |
|---|---|---|
| A | **Status quo** — adjoint closure alone | **Rejected.** §2 shows it accepting a step 6.4 decades from tolerance with the iteration budget exhausted. It also accepts a NaN that lives anywhere outside the force resultant. |
| B | **Convergence alone**, i.e. what `aleph/components/incumbent/driver.py` does | **Rejected.** Blind to a one-sided connector scatter — the exact defect the SF lane's positive control was built to catch. Adopting B would discard a working check. |
| C | **Conjunction (P1)** | **Proposed.** Each clause catches what the others cannot; no clause is redundant (each of A and B is a proper subset that demonstrably misses a control). |
| D | **Disjunction, or a weighted/soft score** | **Rejected without qualification.** Any predicate where a strong adjoint pass can offset a convergence failure is a licence to commit unconverged physics, and a tunable weight is a threshold chosen after seeing data. |
| E | **Keep adjoint-only as the gate, add convergence as a warning field** | **Rejected as an endpoint, adopted as a step** — this is P4. As a permanent arrangement it is the status quo with better bookkeeping: the run still commits. |
| F | **Per-connector vs global closure** | Separated out as **P2** so it can be decided independently. Global is strictly weaker (§2b row 2) but per-connector may need connector-side plumbing that does not exist yet. |
| G | **Loosen the inner tolerance until the current runs pass** | **Rejected — this is the forbidden move.** A threshold chosen to make a gate pass violates the hard rule. If the tolerance is wrong it must be re-derived from the run's own numerics (as D8 does), before the run, with PI sign-off. |

---

## 6. Evidence status of this card

| Item | Status |
|---|---|
| §1 defect statement (code citations) | **verified by reading the code at HEAD `7daedb3d`** |
| §2 control numbers | **measured, CPU, oracle-only** — properties of a predicate, not of the cell |
| §2 `2634.12` pN residual | **cited fixture** from `STATE.md` tier-(a); not re-measured |
| §4 items 1–2 | **cited from other lanes' artifacts**; not measured by Lane A |
| production reject fraction | **BLOCKED** — needs GPU instrumentation; route via **P4** |
| P1 / P2 / P3 adoption | **PI-DECISION** — nothing is adopted by this card |
| tolerance values | **PI-DECISION** — deliberately absent, and the oracle refuses to default them |

**Authority.** The oracle's verdict type is fixed to `authority = "unverified-observer"` and its
`evidence_source` is `ANALYTIC_ORACLE`; constructing it with any other authority raises. It cannot
emit `NATIVE_ACCEPTED`, by construction and by test.

## 7. Verification

- `aleph/tests/ac/virtual_cell/test_admissibility.py` — **37 passed**
- `aleph/tests/ac/virtual_cell/` (whole lane directory) — **125 passed**
- `ruff check` and `ruff format --check` on `aleph/virtual_cell/admissibility.py` and the test file —
  **clean**
- Warp is **not** imported by the module: asserted in a subprocess that imports
  `aleph.virtual_cell.admissibility` and checks `'warp' not in sys.modules`.
- No GPU was used. No GPU lease was taken. `ffn_gpu.py` was not run.

---

**This card CHANGES NOTHING until the PI signs.** Lane A has not touched the engine. The oracle is an
unverified observer; adopting its predicate inside `aleph/engine/transaction.py` is a
gate-contract change and requires PI sign-off.

`PI decision: ______  (P1 __ / P2 __ / P3 __ / P4 __)   date: ______   signature: ______`
