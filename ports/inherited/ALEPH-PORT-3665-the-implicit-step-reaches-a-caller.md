# ALEPH-PORT-3665 — the implicit step reaches a caller

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3665` |
| Lane | `b4fad06b` |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `SOURCE-DERIVED` |
| Completes | `ALEPH-PORT-3664`, whose §1 says *"Nothing is wired into `build_whole_cell` or `relax_to_equilibrium` by this entry"* — so the 708× it measured is reachable by no caller |
| Asked for by | the PI, this lane's transcript, priority table item 1: *"relax_implicit 배선 — 아키텍처와 무관, 12배 문제의 대부분. 제일 급함"* |

---

## 1. Aleph API

`aleph/vertical/relax_implicit.py` gains one driver:

- `relax_world_implicit(world, *, max_steps, dt, ...) -> ImplicitRelaxationReport`

`scripts/engine_clock_run.py` gains `--solver {explicit,implicit}`, defaulting to `explicit`.

**`aleph/vertical/relax.py` is still not edited.** Neither is any owner. The driver is a second
entry point beside `relax_to_equilibrium`, not a modification of it.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **`/Users/sw1/ffn_cellsim`** — READ ONLY, nothing written |
| Source commit | **`be0e58760baaddc04460bb6b34b0476b5d8aa6c5`** (HEAD, 2026-07-29) |
| Source path and symbol | `ffn_sim/ff/implicit_ff.py` — `implicit_step_current`'s `force_fn` contract; `ffn_sim/ff/relax.py` — `relax_implicit`'s driver shape |
| Working tree == commit? | Repository as a whole **NO** (five files modified); the two read here **YES, clean at HEAD**, re-checked with `git status --porcelain -- <paths>` |

## 3. Why source-derived porting beats clean-room

The **operator split** is the thing being taken, and it is a design decision rather than a formula.
`implicit_ff.py`'s header: *"Active/soft forces (turgor, membrane, nucleus, substrate, clutch,
protrusion, gravity, volume) enter `F_total` **EXPLICITLY** (they are soft; keeping them out of K
keeps M SPD and the factorization reusable)."* And `implicit_step_current`'s `force_fn` *"returns
the FULL force (elastic + active) at a given x"* — so the **implicit operator is elastic-only while
the right-hand side is everything**. A clean-room implementation would plausibly put the whole force
in both, lose SPD, and lose the stability argument with it.

## 4. Physical or mathematical law represented

No new law. This is NF2007 Eq (2) applied as an **IMEX operator split**: for the owner whose elastic
stiffness `K` is assembled, the step is

```text
( γ/dt · I  +  K ) · Δx  =  F_total(x)          (implicit in the elastic part)
```

and for every other owner it is the explicit overdamped rule `Δx = dt·μ·F` at the **same** `dt`, so
one global clock still means one thing.

## 5. Units, domains, singular cases, invariants

- **I1 — one `dt` for every owner.** Different owners advancing at different rates would make the
  world clock meaningless. Asserted by construction and checked.
- **I2 — the explicit driver is untouched.** `relax_to_equilibrium` on a default world is
  `np.array_equal` before and after this entry.
- **I3 — the implicit driver does NOT go through the acceptance predicate, and that is the point.**
  See §6.
- **I4 — the cortex sees the FULL force**, its own plus every connector's, taken from
  `world.evaluate_forces()`. An implicit step on the elastic force alone would be a different and
  wrong scheme.
- **I5 — a world with no cortex falls back to the explicit rule for everything**, rather than
  raising: the driver is then exactly the explicit one and says so in its report.

## 6. Source evidence class and known retractions

**This entry removes the accept/reject test from one code path, and that is `ALEPH-DQ-104`
territory. It is recorded as a finding, not decided.**

`MonotonePotentialDescent` exists because an **explicit** step can raise the energy: it is a
stability guard. NF2007 p17's whole content is that the semi-implicit step **cannot** be unstable in
the elastic modes for any `dt`, so the guard has nothing to guard there. `ALEPH-PORT-3637` reaches
the same conclusion from the residency side, in its own words: *"The adaptive step is why there is a
readback at all. **A fixed step needs no scalar and no sync**, but it is either unstable at the
start or uselessly slow at the end"* — and unconditional stability is exactly what removes the first
horn of that dilemma.

**What is NOT claimed:** that the guard is unnecessary in general. `CoveragePredicate` and
`BalancePredicate` check things a solver cannot make true by being stable — whether the connector
books close — and **this driver keeps neither**, which is a real loss and is stated in §14 rather
than argued away. The explicit path remains the ratified one and this is opt-in.

## 7. Independent oracle or derivation

1. **The explicit driver is the oracle.** At `dt` below the explicit stability limit the two must
   reach the same equilibrium, because both solve `F = 0`. Checked on a default world.
2. **`ALEPH-PORT-3664` measured 708× on the cortex in isolation. This entry's prediction is that
   the whole-cell figure is LOWER, and it says so before measuring.** The cortex is 93,679 of
   93,919 nodes but it is **not** the only stiff thing: `-3664` §7(1) measured that on the *default*
   construction the binding mode is the **microtubule** at `λ_max = 7.24e4`, and the microtubule
   stays explicit here. So the whole-cell `dt` is bounded by
   `min(accuracy limit, 2/(μ·λ_max of everything still explicit))`.
   **Predicted: a whole-cell speed-up well under 708×, and possibly under 10×.** If it comes out
   near 708× then the explicit remainder is softer than measured and that is a finding; if it comes
   out near 1× the split is buying nothing and that is a bigger one.
3. **The clock must still advance monotonically** and the residual must fall. A driver that is
   "fast" because it has stopped converging is the failure this control exists over.

## 8. Positive control

`tests/vertical/test_relax_world_implicit.py`:

- `test_both_drivers_reach_the_same_equilibrium` — §7(1).
- `test_one_dt_for_every_owner` — I1.
- `test_the_cortex_sees_the_connector_force` — I4, by checking the force handed to the solve differs
  from the cortex's own `internal_forces`.
- `test_a_world_without_a_cortex_is_all_explicit` — I5.

## 9. Deliberately failing negative control

- `test_the_explicit_driver_is_unchanged` — I2, `np.array_equal`.
- `test_the_implicit_driver_is_not_secretly_explicit` — the **vacuity control**: at a `dt` above the
  explicit stability limit the explicit driver's displacement blows up and this one's does not.
- `test_the_residual_falls` — §7(3).

## 10. Numerical and precision envelope

float64, host. CG tolerance is an argument and is reported. No bit-identity is claimed between the
two drivers — they are different schemes; §7(1) is an equilibrium statement.

## 11. Production-backend residency and transfer

**Host only, and deliberately so for now.** The PI's priority table item 4 freezes `resident_*`:
*"카드가 2배 느린 동안 디바이스 사본을 유지하지 않음."* This lane measured that the card is
**2.03× slower** than the laptop on this engine, so a device-resident implicit solve is deferred
until the card wins. `-3664` §11 already scoped that as a second port.

## 12. Comments and docstrings to discard

Nothing copied.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **`ACCEPTED` when the whole-cell speed-up is measured**
with `scripts/engine_clock_run.py --solver implicit --construction sourced` and written into
`docs/design/ENGINE_PERFORMANCE_TARGET.md` beside the explicit row — **including if it is far below
708×**, which §7(2) predicts it will be.

## 14. Honest limits

- **Coverage and balance are not checked on this path.** The explicit driver runs
  `AllOf(CoveragePredicate, BalancePredicate, MonotonePotentialDescent)` every step; this one runs
  none of them. Dropping the descent test is argued in §6; dropping the other two is **not**
  argued — it is a consequence of not going through the transaction, and it means a connector whose
  books stop closing would not be caught on this path. **Anyone using it for a result must run the
  explicit audit on the final state.**
- **No rollback.** The transaction rolls back state, clock, topology and RNG together on a rejected
  step. This driver has nothing to roll back to.
- **One owner is implicit.** Twelve are not.
- **The bending block inherits `INTERIOR_ONLY`'s deficit**, measured against Cytosim today at
  **0.497 of the continuum at `n = 3`** — the cortex's resolution. It is 0.33 % of `‖K‖` so it does
  not move `dt` today, and it would if the crosslinker ever stopped dominating.
