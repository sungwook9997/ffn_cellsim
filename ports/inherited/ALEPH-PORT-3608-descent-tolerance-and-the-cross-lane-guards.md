# ALEPH-PORT-3608 — a descent tolerance that is a statement about a precision, and two guards that make a remembered rule mechanical

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3608` |
| Lane | `46143f30` Lane G8 (CUDA track — decisions) |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | The PI answered `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md` questions 4, 5 and 6 at **2026-08-01 01:39 KST**. Three things follow, and none of them is a kernel: the descent criterion needs a tolerance that is a function of the precision its energy is read in; the oracle-contamination hazard of that proposal's §4 needs a guard rather than six lanes remembering; and the position mode has to be *in* a run artefact rather than in a note beside it. |

> **The instruction, quoted with coordinates so a session that did not receive it can check rather
> than trust** (`CLAUDE.md` §2 rule 2). Transcript
> `~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl`, user
> turn at **2026-08-01 01:39 KST**, verbatim:
>
> > *"ㅇㅇ gpu 권한 변경해줬고, 나머지는 다 너의 의견대로, 거부해야하는 이유는 없다고 생각함. 그리고
> > 수렴하는 것은 굉장히 어려운 것임. 오차를 주는 것 좋다고 생각함. 그리고 라플라스 다시 유도해봐도
> > 좋은데, 그거는 주말 동안 작업을 하지 못해도 gbook은 계속 돌릴 수 있으니깐 그렇게 먹여두는 것으로
> > 하자 토요일, 일요일 작업 중에 같이 보게 그리고 월요일날 같이 검토하면 되니깐"*
>
> **This entry carries no `decided_by` field and no agent may add one.** The record of what the
> sentence is taken to decide, and of the one clause in it that is genuinely ambiguous, is
> `docs/decisions/RATIFICATION-2026-08-01-position-precision-and-the-descent-tolerance.md`.

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered by it.

```python
# MODIFIED — aleph/vertical/relax.py, a VERIFIED module
from aleph.vertical.relax import (
    FLOAT32_UNIT_ROUNDOFF,              # NEW  5.960464477539063e-08 = eps32 / 2
    FLOAT64_UNIT_ROUNDOFF,              # NEW  1.1102230246251565e-16 = eps64 / 2
    DESCENT_SLACK_ROUNDOFFS,            # NEW  2 — derived in §4b, not chosen
    ENERGY_SCALE_FLOOR_PN_UM,           # NEW  1.0 pN.um — the old inline `max(|E|, 1.0)`, named
    MonotonePotentialDescent,           # CHANGED: `energy_unit_roundoff`, `roundoffs`,
                                        #          `relative_slack` becomes an explicit OVERRIDE
    RelaxationReport,                   # CHANGED: records the slack that was actually in force
    relax_to_equilibrium,               # CHANGED: `energy_unit_roundoff=`, `descent_relative_slack=`
)

# MODIFIED — aleph/runtime/parity.py
from aleph.runtime.parity import (
    LawCase,                            # CHANGED: `position_mode` field
    LawParity,                          # CHANGED: `position_mode` field
    LawParityReport,                    # CHANGED: `position_mode`, `position_modes` fields
    POSITION_MODE_COMPARISON_POLICY,    # NEW  the single declared flag, defaulting to refuse
    ReferenceInputError,                # NEW  the typed refusal §5 describes
    IncomparableArtefactsError,         # NEW  the typed refusal §6 describes
    refuse_incomparable_position_modes, # NEW
    compare_law_case,                   # CHANGED: refuses a narrowed reference input
    law_report_as_dict,                 # CHANGED: the mode and the policy are in the artefact
)

# MODIFIED — aleph/runtime/law_cases.py
from aleph.runtime.law_cases import (
    LAW_CASE_BUILDERS,                  # NEW  name -> (builder, its TRUE variant)
    every_law_case,                     # NEW  one built case per registered builder
)

# NEW — tests
tests/vertical/test_descent_tolerance.py
tests/runtime/test_reference_input_guard.py
tests/runtime/test_position_mode_in_artefacts.py
```

**Not authorised and not touched:** `aleph/runtime/law_kernels.py` (frozen by G2–G5 and driven by
G7), every module under `aleph/vertical/` **other than `relax.py`**, `aleph/scenarios/**`,
`aleph/viz/**`, `scripts/gpu_preflight.py`, `ports/ledger/INDEX.md`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

Everything read is Aleph's own: `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md`,
`ALEPH-PORT-3601`/`-3603`/`-3604`/`-3605`/`-3606`, `docs/results/2026-07-31-tension-sweep/`,
`PLAN.md` §2.5 and §6.1, and the five `⚠ G1…G5 LANDED` blocks of
`docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md`.

## 3. Why source-derived porting beats clean-room

**It does not, and this entry says so rather than manufacturing a reason.** A summation error bound
is textbook (Higham, *Accuracy and Stability of Numerical Algorithms*, §4.2), and the two guards are
answers to a hazard *this repository* discovered and wrote down. Importing a provider's stopping
rule would import their answer to the one question that matters here, which is what tolerance this
project's own oracle can afford.

The entry is `RE-DERIVED` and is kept because `PLAN.md` §0.2.5 wants the decision recorded, not only
the code — and because `aleph/vertical/relax.py` is a `VERIFIED` module whose published results have
to be **re-derived** rather than re-run when it changes.

## 4. Physical or mathematical law represented

### 4a. There is no physical law in the descent criterion, and that is stated first

`MonotonePotentialDescent` is scaffolding. Its own module docstring says it is **not** a candidate
answer to `ALEPH-DQ-104` (what makes a step physically acceptable), and nothing here changes that. A
tolerance on a scaffolding predicate is a numerical statement, not a physical one.

### 4b. What the tolerance is a statement about — derived, and the first derivation was wrong

The predicate compares `E_after` against `E_before`, each computed by the same code path in a
channel of unit round-off `u`. Two facts about *this* potential shape the derivation:

1. **Every term is non-negative.** `sigma*A`, `kappa*|K|^2/A`, every `k x^2 / 2` spring, the turgor's
   `K (V - V0)^2 / (2 V0)`, and every clamped connector spring. So the cancellation amplification
   `sum|E_i| / |sum E_i|` is **exactly 1** — the assembled total has no cancellation at all, and the
   only error is accumulation. (Asserted, not assumed: §8 control C2a.)
2. **The predicate reads a difference of two evaluations at different configurations**, whose
   round-off errors are independent once the displacement exceeds the rounding scale, so both enter.

#### The derivation in force

```
allowance   =  2 evaluations  x  1 unit round-off of the channel each

DESCENT_SLACK_ROUNDOFFS = 2                                (dimensionless, derived)

relative_slack(float64) = 2 * 1.1102230246251565e-16 = 2.220446049250313e-16
relative_slack(float32) = 2 * 5.9604644775390625e-08 = 1.1920928955078125e-07
```

**The same expression produces both**, which is the whole grounds for this option over widening the
energy array. A relative tolerance of this form is correct at *any* precision, so it introduces no
exception to `ALEPH-DQ-107`; widening one array would have introduced one.

#### The bound this entry first took, and why the measurement rejected it — this is the finding

For pairwise summation of `n` non-negative terms the classical bound (Higham §4.2) is
`|E_hat - E| <= gamma_k E` with `gamma_k = k u / (1 - k u)` and `k = ceil(log2 n)`. Declaring the
domain `n <= 4096` gives `k = 12`, and doubling for two evaluations gives a factor of **24**. This
entry was first written with that factor, and then with **4**.

**Neither survives contact with the thing it is a tolerance for.** Driven by `ALEPH-PORT-3603`'s
**real** float32-term membrane kernels on the vertical (§10.2), the window is **one binary order of
magnitude wide** and the derived value is the only point inside it:

| slack | `L1 sigma=10` | `L2 sigma=20` |
|---|---|---|
| `0` (exact) | stalls at 442 steps, step 6.4e-13 | stalls at 296 steps, step 6.1e-13 |
| `1 * u32` | **stalls** at 979 steps, step 5.7e-13 | converges, 2,806 steps, sigma 2.1e-05 |
| **`2 * u32` (derived)** | **converges, 2,742 steps, sigma 3.3e-06** | **converges, 4,610 steps, sigma 3.2e-05** |
| `4 * u32` | wanders 35,165 steps, sigma 6.2e-03 | wanders, sigma 1.1e-02 |
| `6 * u32` | wanders, sigma 2.2e-02 | wanders, sigma 6.3e-03 |
| `8 * u32` | wanders, sigma 2.0e-02 | wanders, sigma 2.6e-02 |
| `12 * u32` | wanders, sigma 2.4e-02 | wanders, sigma 2.8e-02 |
| `24 * u32` (the bound) | wanders, sigma **3.9e-02** | wanders, sigma **4.5e-02** |

The distinction this entry had to learn: **a tolerance is not only a bound on the evaluation error,
it is also the amount of uphill wander the descent is permitted**, and the two requirements pull in
opposite directions. A conservative error bound is the right thing to quote about accuracy and the
wrong thing to hand a stopping rule. "Tighter is safer" is the instinct, and it is wrong on the low
side too — `1 * u32` stalls.

#### What the assembled error actually is, measured

Rigid translation leaves the true energy exactly invariant (every term is a function of position
differences), so the computed drift under one is pure round-off with no model in between. Across
**all twelve** tension-sweep configurations the drift is **exactly one or two ULP of the result and
never more** (§10.1). The bound's `k = 12` is loose by an order of magnitude because the pairwise
errors do not align; what is left is the final representation.

**The float64 side of that has a consequence this entry states rather than buries.** Because the
host's error is a *quantisation* at whole ULP of the result, a step whose true decrease is below one
ULP yields a computed difference of **exactly `0.0`**, which a monotone test accepts at *any* slack
including zero. That is why the float64 scan is flat from `0` to `1e-10` (§10.3), and why the worst
drift — `2.64 u` at `L2_sigma10`, sitting outside the `2 u` allowance — is **not** a floor
violation. A float32-**term** channel perturbs the total instead of quantising it, and that
asymmetry is the entire reason this decision was needed.

`allowance = relative_slack * max(|E_before|, ENERGY_SCALE_FLOOR_PN_UM)`. The floor is the existing
inline `1.0`, now named: without it the allowance would collapse to zero for an assembly whose total
potential happens to pass through zero, and a relative tolerance on a quantity that can vanish is
not a tolerance.

### 4c. What this replaces, stated precisely, because "exact comparison" is not quite what was there

The predicate was **not** an exact `<=`. It carried `relative_slack = 1.0e-12`, described as "float64
dust". Three things are wrong with that and are what this entry fixes:

* **It is a round number.** `1e-12` is `9,007 u64`, and nothing derives 9,007.
* **It is a float64 constant with no float64 in its name.** Read the same energy in float32 and the
  constant is 5.4e8 times too tight; `ALEPH-PORT-3603` §13 measured the consequence — **0 of 12**
  configurations converge and recovered sigma degrades **331x**.
* **It says nothing about which channel produced the energy**, so a run artefact could not record
  the stopping rule that produced it. `RelaxationReport` now does.

### 4d. The one law-like statement in this entry, and it is about where a run stopped

`PLAN.md` records that the Laplace relative error **is exactly the dilational virial of the residual
force**, `(sum_i f_i . x_i) / (2 sigma A)`, matched to 2.2e-16. It is therefore a statement about
*where the relaxation stopped*, and loosening a stopping rule moves it directly. That is why §10
re-derives the twelve-configuration sweep instead of re-running it, and why a different number would
be a finding rather than an equivalent.

**The mechanism by which a slack could move it, derived so the measurement has something to falsify:**
an allowance `A` permits an accepted step that raises the energy by up to `A`. Near the minimum, an
energy excess `A` in a mode of stiffness `k` corresponds to a residual force `sqrt(2 k A)`. The
relaxation's convergence gate is a **force** tolerance (`force_tolerance_pn = 1e-3`), so the slack
cannot move the stopping point at all while `sqrt(2 k A) << 1e-3 pN`, and moves it abruptly once the
accepted wander exceeds that. §10 measures where that transition is rather than assuming it.

### 4e. The guards represent a hazard, not a law

`PROPOSAL-position-precision-for-the-cuda-port.md` §4: a law-parity case that hands its **reference**
float32-rounded input makes the oracle and the subject wrong together, and the gate prints a small
healthy number. That is not a law; it is a failure mode of the harness, and the correct form for it
is a **refusal in the harness** plus a **completeness check over the registry**, not a per-case test
six lanes have to remember.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `E_before`, `E_after` | pN·µm | 1e-18 J | finite; the predicate rejects non-finite |
| `relative_slack` | dimensionless | — | `>= 0`; `0.0` reproduces the exact comparison |
| `energy_unit_roundoff` | dimensionless | — | `(0, 1)`; `2^-53` or `2^-24` in practice |
| `accumulation_factor` | dimensionless | — | `4` = two evaluations x one ULP of the result each |
| `ENERGY_SCALE_FLOOR_PN_UM` | pN·µm | 1e-18 J | `> 0` |
| `allowance` | pN·µm | 1e-18 J | `>= 0` |
| `position_mode` | — | — | a `PositionPrecision` value, or `""` for a case that declares none |

Singular and boundary cases, each with the behaviour Aleph requires:

- **`E_before` absent from the candidate.** Unchanged: a named failing condition, never a crash.
- **The candidate's energy is not evaluable.** Unchanged: a rejection carrying the exception text.
- **`E_before = 0`.** The floor makes `allowance = relative_slack * 1.0 pN.um`, so the test does not
  degenerate to an exact comparison at the one configuration where it would be least meaningful.
- **`relative_slack = 0.0`.** Legal, and it is the *old* exact comparison. It is not the default and
  it is what the negative control N1 asserts still stalls a float32-term channel.
- **An assembly whose energy drifts by more than two ULP of the result.** Outside the domain the
  factor 2 was derived for. **Detected, not silently tolerated**: control C2 measures the drift in
  ULP of the result and asserts it is a *whole number* and `<= 2`, so an assembly whose accumulation
  really does grow with `n` goes red on a measurement rather than quietly widening the error the
  predicate cannot see past.
- **A case state carrying a float32 array.** `compare_law_case` raises `ReferenceInputError` and
  grades nothing. Refusal, not degradation.
- **Two artefacts at different position modes.** `refuse_incomparable_position_modes` raises
  `IncomparableArtefactsError` under the default policy.

Invariants that must hold, each with the test that asserts it:

- **I1.** The derived slack is exactly `accumulation_factor * energy_unit_roundoff`, in both channels.
  → `tests/vertical/test_descent_tolerance.py::test_the_slack_is_the_roundoffs_times_the_channel_roundoff`
- **I2.** The assembly's own energy-evaluation noise is inside the allowance in force, is a whole
  number of ULP of the result, is at most two, and the measurement is non-vacuous. →
  `::test_the_host_energy_error_is_a_whole_number_of_ulp_of_the_result`
- **I3.** A float32-term energy channel stalls under the float64 slack and converges under the
  float32 one. → `::test_a_float32_energy_channel_stalls_under_a_float64_slack`,
  `::test_a_float32_energy_channel_converges_under_the_derived_float32_slack`
- **I4.** The predicate still rejects a genuine increase larger than the allowance. →
  `::test_an_increase_larger_than_the_allowance_is_still_rejected`
- **I5.** The report records the slack that was in force. →
  `::test_the_report_records_the_stopping_rule_that_produced_it`
- **I6.** Every public `*_law_case` builder in `law_cases.py` is in `LAW_CASE_BUILDERS`. →
  `tests/runtime/test_reference_input_guard.py::test_every_law_case_builder_is_registered`
- **I7.** Every registered case's state is float64 throughout and carries at least one position
  array. → `::test_every_registered_case_feeds_its_reference_float64_positions`
- **I8.** Every registered case's reference is *measurably* sensitive to that precision. →
  `::test_every_registered_case_reference_is_sensitive_to_position_precision`
- **I9.** `compare_law_case` refuses a case whose state narrows the reference's input. →
  `::test_the_harness_refuses_a_case_that_feeds_its_reference_float32`
- **I10.** Every registered case declares a position mode, and it reaches the artefact. →
  `tests/runtime/test_position_mode_in_artefacts.py::test_every_case_declares_a_position_mode_and_it_reaches_the_artefact`
- **I11.** The comparison policy defaults to refuse, and one flag changes it. →
  `::test_the_comparison_policy_refuses_across_modes_and_one_flag_changes_it`

## 6. Source evidence class and known retractions

No source symbol is claimed, so there is nothing to retract from `ffn_cellsim`. What this entry
inherits and must state:

* `docs/results/2026-07-31-tension-sweep/README.md` is `ANALYTIC_ORACLE`. Laplace's law is an
  identity of the model, not a measurement of a cell. Every sigma below is `UNSOURCED` as physics.
* `ALEPH-PORT-3603` §13's 0-of-12 measurement was made on **warp's CPU device**, not on CUDA, and
  is `UNVERIFIED` as a statement about GPU behaviour. This entry does not change that and runs no
  device.
* `ALEPH-PORT-3605` §13.4 and `ALEPH-PORT-3606` measured **opposite signs** for whether
  `POSITIONS_F64` beats `LOCAL_F32`. Both stand; §14 records that this entry does not resolve them
  and must not be read as doing so.

## 7. Independent oracle or derivation

**The oracle is the twelve-configuration Laplace tension sweep**, and it is not agreement with
anything: `sigma_recovered = -(sum_i f_i . x_i) / (2 A)` is an exact identity of the discrete energy
(the bending term is degree zero in the vertex positions and the tension term degree two), so
recovering the *declared* sigma across a 6x tension range and two resolutions is a statement the
model can fail.

The derivation of the tolerance itself (§4b) is the second oracle, and it is checked against a
**model-free measurement**: every term of the assembly's potential is a function of position
*differences*, so a rigid translation changes the true energy by exactly zero. Whatever the computed
energy does under one is therefore pure round-off, measured at exactly the displacement scale the
predicate operates at, with no model in between. §10 reports it.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| C1 | `tests/vertical/test_descent_tolerance.py::test_the_slack_is_the_roundoffs_times_the_channel_roundoff` | The derived slack is exactly `4 * u` in both channels, the factor is exactly `2 * 2`, and the two constants are exactly `eps/2` of their dtypes |
| C1b | `::test_an_explicit_slack_overrides_the_derivation` | The override wins over the channel, so the two parameters cannot silently multiply; `0.0` is reachable |
| C2 | `::test_the_host_energy_error_is_a_whole_number_of_ulp_of_the_result` | Rigid-translation drift is `<= allowance` on four configurations **and is a whole number of ULP of the result, at most two** — the half that makes the derivation falsifiable rather than merely satisfied |
| C2a | `::test_the_assembly_potential_has_no_cancellation_so_the_amplification_is_one` | Every energy term is non-negative and the amplification is exactly 1, which is what licenses an accumulation-only bound |
| C2b | `::test_the_translation_probe_is_not_measuring_zero_everywhere` | At least one probed configuration moves at all (the vacuity half — 6 of 12 are bitwise invariant) |
| C3 | `::test_a_float32_energy_channel_stalls_under_a_float64_slack` and `::..._converges_under_the_derived_float32_slack` | On `ALEPH-PORT-3603`'s **real** float32-term kernels: the float64 slack collapses the step to `< 1e-11` in `< 1000` steps; the derived float32 slack converges with the step above `1e-6` |
| C3b | **UNTESTED** — measured once and never fixed by a control | `24 * u32` — the pairwise-summation bound — does **not** converge, which is the grounds for the factor being 4 |
| C4 | `::test_an_increase_larger_than_the_allowance_is_still_rejected` | A scripted world whose energy rises by `2 * allowance` is rejected; by `allowance / 2`, accepted; the reported `threshold` is the allowance |
| C5 | `::test_the_report_records_the_stopping_rule_that_produced_it` | `RelaxationReport.descent_relative_slack` and `.energy_unit_roundoff` equal what was passed, and the slack is in `describe()` |
| C5b | `::test_the_float64_default_still_relaxes_the_vertical_and_recovers_the_tension` | One sweep configuration in-suite, so a change to the default cannot pass the suite while costing convergence |
| C6 | `tests/runtime/test_reference_input_guard.py::test_every_registered_case_feeds_its_reference_float64_positions` | All 10 cases: every float array float64; `>= 1` array with a trailing axis of 3 |
| C7 | `::test_every_registered_case_reference_is_sensitive_to_position_precision` | For each case, float32-rounding its positions moves the reference by `> 1` ULP on at least one field — measured 2.07 to 519 — and **all ten take the numeric path**, so the control is not silently in its refusal branch |
| C7b | `::test_a_clean_case_is_not_refused` / `::test_the_refusal_ignores_integer_and_boolean_state` | The refusal fires on float32 and float16 and not on `int8`/`bool`, so it protects something and is not noise |
| C8 | `tests/runtime/test_position_mode_in_artefacts.py::test_every_case_declares_a_position_mode_and_it_reaches_the_artefact` | Every case's `position_mode` is non-empty, is a `PositionPrecision` value, and appears in `law_report_as_dict` and in `render()` |
| C8b | `::test_a_declared_mode_actually_follows_the_argument` | The field tracks the argument rather than being hard-wired, and `tether_law_case`'s single-dtype exception is asserted rather than skipped |
| C9 | `docs/results/2026-08-01-descent-tolerance/README.md` | The re-derived sweep, row for row, against the published table |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| N1 (must fail) | `tests/vertical/test_descent_tolerance.py::test_an_exact_comparison_stalls_the_same_way` | `relative_slack = 0.0` on a float32-term channel stalls at **the same accepted-step count and the same residual** as the float64 slack — the vacuity guard on C3, and the statement that `1e-12` was never doing any work in this channel |
| N2 (must fail) | `tests/runtime/test_reference_input_guard.py::test_the_harness_refuses_a_case_that_feeds_its_reference_float32` | `compare_law_case` raises `ReferenceInputError` on a case whose state was float32-rounded, and the message names the hazard |
| N3 (must fail) | `::test_the_registry_completeness_check_can_fail` | With one builder hidden from the registry, the completeness assertion is violated — so a new unregistered builder really would be caught |
| N4 (must fail) | `tests/runtime/test_position_mode_in_artefacts.py::test_the_comparison_policy_refuses_across_modes_and_one_flag_changes_it` | Two artefacts at different modes raise under the default policy; the same pair does not raise with the flag set to `"allow"`; and two artefacts at the *same* mode never raise |

## 10. Numerical and precision envelope

**Working precision.** The vertical is float64 throughout; the descent predicate's channel is
declared by `energy_unit_roundoff` and defaults to float64. Nothing in this entry moves
`ALEPH-DQ-107`'s float32 compute channel.

**Host for every figure below:** `gbook` for the twelve-configuration sweeps (16 threads,
`OMP_NUM_THREADS=1`, 10 workers), the PI's Mac for the single-configuration probes. §10.4 reports
what moving hosts did and did not change. **No GPU, on either host.**

### 10.1 What the assembled energy's error is — rigid-translation probe, all twelve configurations

The true energy is exactly invariant under a rigid translation, so this is pure round-off with no
model in between, at the displacement scale the predicate operates at (seven magnitudes, 1e-9 to
1e-3 um, eight directions each).

| configuration | `E` [pN·µm] | worst `\|dE\|` | relative | in `u64` | **in ULP of `E`** |
|---|---|---|---|---|---|
| L1_sigma10 | 2.9364e+03 | 4.5475e-13 | 1.5487e-16 | 1.39 | **1** |
| L1_sigma20 | 5.8763e+03 | 0 | 0 | 0 | **0** |
| L1_sigma30 | 8.8216e+03 | 1.8190e-12 | 2.0620e-16 | 1.86 | **1** |
| L1_sigma40 | 1.1772e+04 | 0 | 0 | 0 | **0** |
| L1_sigma50 | 1.4728e+04 | 0 | 0 | 0 | **0** |
| L1_sigma60 | 1.7689e+04 | 0 | 0 | 0 | **0** |
| L2_sigma10 | 3.1034e+03 | 9.0949e-13 | 2.9307e-16 | **2.64** | **2** |
| L2_sigma20 | 6.2109e+03 | 9.0949e-13 | 1.4643e-16 | 1.32 | **1** |
| L2_sigma30 | 9.3248e+03 | 1.8190e-12 | 1.9507e-16 | 1.76 | **1** |
| L2_sigma40 | 1.2445e+04 | 1.8190e-12 | 1.4617e-16 | 1.32 | **1** |
| L2_sigma50 | 1.5571e+04 | 1.8190e-12 | 1.1682e-16 | 1.05 | **1** |
| L2_sigma60 | 1.8703e+04 | 3.6380e-12 | 1.9452e-16 | 1.75 | **1** |

**Every entry is a whole number of ULP of the result, and never more than two.** The assembled error
is the final representation, not accumulation growth. Six of twelve are bitwise invariant, which is
why the vacuity control C2b exists.

### 10.2 The float32 window — `ALEPH-PORT-3603`'s real kernels, `max_steps = 40,000`

Membrane energy from the float32-term kernels, forces NumPy float64 — `ALEPH-PORT-3603` §13's
isolation. `ORDERED` scatter, warp **CPU** device.

| slack | L1 sigma=10 | | | L2 sigma=20 | | |
|---|---|---|---|---|---|---|
| | conv | acc | rel err | conv | acc | rel err |
| `0` (exact) | no | 442 | 3.36e-04 | no | 296 | 8.83e-04 |
| `1 * u32` | **no** | 979 | 1.64e-04 | yes | 2,806 | 2.13e-05 |
| **`2 * u32`** | **yes** | **2,742** | **3.30e-06** | **yes** | **4,610** | **3.24e-05** |
| `4 * u32` | no | 35,165 | 6.18e-03 | no | 35,166 | 1.08e-02 |
| `6 * u32` | no | 35,164 | 2.25e-02 | no | 35,166 | 6.33e-03 |
| `8 * u32` | no | 35,164 | 1.99e-02 | no | 35,166 | 2.60e-02 |
| `12 * u32` | no | 35,165 | 2.42e-02 | no | 35,165 | 2.82e-02 |
| `24 * u32` | no | 35,165 | 3.92e-02 | no | 35,166 | 4.50e-02 |

The `0` row **reproduces `ALEPH-PORT-3603` §13's stall signature exactly**: the step collapses to
`6.4e-13`, below `min_step`, after a few hundred accepted steps.

### 10.3 The float64 window — the same two configurations, host NumPy

| slack | /u64 | L1 sigma=10 | L2 sigma=20 |
|---|---|---|---|
| `0` (exact) | 0 | conv, 2683, 4.39e-06 | conv, 2732, 1.566e-05 |
| **`2 * u64`** | **2** | **conv, 2683, 4.39e-06** | **conv, 2732, 1.566e-05** |
| `1e-13` | 901 | conv, 2675, 2.84e-06 | conv, 2732, 1.566e-05 |
| `1e-12` (the old constant) | 9,007 | conv, 2689, 4.82e-06 | conv, 2732, 1.566e-05 |
| `1e-10` | 9.0e+05 | conv, 3307, 1.53e-05 | conv, 2742, 3.39e-05 |
| `1e-8` | 9.0e+07 | **REFUSED**, 8.78e-04 | conv, 3367, 3.47e-05 |
| `1e-6` | 9.0e+09 | **REFUSED**, 9.62e-03 | **REFUSED**, 3.73e-02 |
| `1e-4` | 9.0e+11 | **REFUSED**, 4.70e-01 | **REFUSED**, 4.11e-01 |

**Flat from `0` to `1e-10`** — seven orders — because the host's error is a quantisation and a
sub-ULP step gives a computed difference of exactly `0.0`. The float64 channel never needed a
tolerance, which is exactly why the float32 one was not noticed until `ALEPH-PORT-3603` §13.

### 10.4 The re-derived Laplace sweep, and what changing host did

Full table: `docs/results/2026-08-01-descent-tolerance/README.md`.

| | before (`1.0e-12`) | after (`2 * u64`) |
|---|---|---|
| **worst relative error** | **1.5663e-05** | **1.5663e-05** |
| **converged** | **12 / 12** | **12 / 12** |
| worst row | `L2_sigma20` | `L2_sigma20`, and it is **bit-identical** |
| total accepted steps | 38,100 | 37,888 (−0.6 %) |

**Changing host moved nothing that is reported and something that is not.** The same sweep at the
same slack, run on the Mac (10 cores, default BLAS threading) and on `gbook` (16 threads,
`OMP_NUM_THREADS=1`), gives **identical accepted-step counts and identical convergence on all twelve
rows**, and recovered sigma agreeing to `~3e-14` relative — the last two or three digits differ. That
is a BLAS/vectorisation difference on a float64 reduction; it is 9 orders below the 1.5663e-05 the
sweep is quoted by and does not touch any figure in this entry. **Stated rather than omitted**,
because "the sweep is deterministic" would have been the natural thing to assume and it is not quite
true across architectures.

**Range over which the stated tolerance holds:** an assembly whose potential terms are all
non-negative (§4b fact 1, asserted by C2a), whose assembled energy error is at most two ULP of the
result (§10.1, asserted by C2), and a convergence gate that is a **force** tolerance rather than an
energy one (§4d). Outside any of those the derivation does not apply and C2/C2a go red rather than
the run degrading silently.

## 11. Production-backend residency and transfer

`relax.py` runs on the **host** in float64 and this entry does not move it. What it does is make the
predicate able to *state* which channel produced the energy it reads, which is the prerequisite for
G6 keeping state resident on the device: a device-side energy is a float32-term channel, and until
now the predicate had a float64 constant with no way to say so.

The two guards run on the host and touch no device. `compare_law_case`'s refusal happens **before**
any launch, so a contaminated case never reaches a kernel.

**No CUDA device was touched by this lane.** A grant may exist for host `GBook` device 0; an agent
may never write that record and this one did not read it as authorization either.

## 12. Comments and docstrings to discard

Nothing is ported, so no provider prose crosses. What is **removed from Aleph's own text**, and what
replaces it:

| Removed | Why | Replaced by |
|---|---|---|
| `relative_slack: float = 1.0e-12` and the phrase *"so a step is not rejected on float64 dust"* | It is a round number, and "float64" is in the prose but not in the constant | `accumulation_factor * energy_unit_roundoff`, derived in §4b, with the channel a named parameter |
| the inline `max(abs(before), 1.0)` | An unexplained `1.0` in pN·µm | `ENERGY_SCALE_FLOOR_PN_UM`, with the reason it exists in its own docstring |

The module docstring's standing disclaimers — not dynamics, not thermal, not implicit, not a
minimiser with guarantees, and **not a proposal for `ALEPH-DQ-104`** — are kept verbatim and one
sentence is added saying that giving the predicate a tolerance does not change any of them.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **The twelve-configuration Laplace sweep re-derived on `gbook`, 2026-08-01: worst relative error `1.5663e-05` before and after, 12 of 12 converged both, and the worst row bit-identical.** `tests/runtime` + `tests/vertical` **1,808 passed / 3 skipped** (was 1,770/3). Mutation study **19 planted / 19 killed**, from a verified-empty bytecode tree |
| Reviewer | **agent-proposed; unratified.** No `decided_by` field exists in this entry and no agent may add one |
| Rollback | Revert `aleph/vertical/relax.py`, `aleph/runtime/parity.py`, `aleph/runtime/law_cases.py` and delete the three test modules. What breaks: nothing that runs today — the float64 path is behaviourally identical at the accepted slack (§13.1) — but the float32 energy channel returns to 0-of-12, and the §4e hazard returns to being six lanes' memory |

### 13.1 The re-derived tension sweep

`aleph/vertical/relax.py` is `VERIFIED`, so this is a **re-derivation** and not a re-run: the
baseline arm was measured first, on the same host, at the old `1.0e-12`, and it reproduces the
published table row for row including every accepted-step count. Without that check the comparison
would be uninterpretable.

| | before (`1.0e-12`) | after (`2 * u64 = 2.220446049250313e-16`) | published table |
|---|---|---|---|
| worst relative error | **1.5663e-05** | **1.5663e-05** | 1.566e-05 |
| converged | **12 / 12** | **12 / 12** | 12 / 12 |
| worst row | `L2_sigma20` | `L2_sigma20`, **bit-identical** (`19.99968673096192`, 2,732 steps) | `L2_sigma20` |
| total accepted steps | 38,100 | 37,888 | — |

**Individual rows moved, and that is a finding rather than a caveat.** `L1_sigma30` improved from
5.84e-06 to **8.99e-09**; `L2_sigma30` degraded from 3.54e-06 to **1.56e-05**. Seven improved, three
degraded, two are bit-identical. `PLAN.md` §2.5 is the reason this is expected and not alarming:
across twenty stopping points on one fixed mesh the signed errors spanned **603x** with five sign
changes, because a per-row error is a statement about *where the run stopped* and not about whether
the model is right. The sweep's claim was never a row — it was "sigma recovered = sigma declared
across a 6x tension range and two resolutions, to a worst 1.566e-05" — and that claim is unchanged.

**One number to be honest about:** the worst case is now attained twice, at `L2_sigma20`
(1.5663e-05, unchanged) and `L2_sigma30` (1.5640e-05, up from 3.54e-06). Under the old rule the
second-worst row was 9.87e-06. The distribution tightened at the top; the bound did not move.
Nothing here licenses a claim that the sweep got better.

**`ALEPH-PORT-3603` §13's kernel-driven figure is a different measurement and is not re-derived
here.** That lane recovered **1.5610e-05** with the membrane *forces* on kernels and the descent
energy in float64. This entry changed neither of those channels, so that number is untouched; what
it changed is the third channel that lane could not use at all.

## 14. Honest limits

1. **This entry does not answer questions 1, 2 or 3 of the proposal.** Which of A/B/C/D the port
   carries, whether it is global or per-law, and who owns a reference frame are all still open. G4
   and G5 measured **opposite signs** for `POSITIONS_F64` vs `LOCAL_F32`, and `ALEPH-PORT-3606`
   established that `LOCAL_F32` destroys family A's exact zero. Nothing here narrows that.
2. **The guard addresses the INPUT half of the hazard and not the LAW half.** `ALEPH-PORT-3606`'s M8
   set `bidirectional=False` on the frozen spring: both sides then computed a tether, in perfect
   agreement, with nine rows of physics gone — **a parity gate cannot see a defect in its own
   reference.** Nothing in §5 or §8 changes that, and this sentence exists so nobody reads the guard
   as more than it is. The series now stands at four: G1 — oracle and subject share a float32
   representation; G4 — the gate cannot see a defect in the state it grades on; G5 — nor in the
   reference it grades against; G8 — the input half of the first of those is now mechanical, and the
   other two are not.
3. **The float32 window is one binary order of magnitude wide, and nothing enforces staying in it.**
   `1 * u32` stalls and `4 * u32` wanders (§10.2), so the derived `2 * u32` has **no margin on
   either side** on the two configurations measured. That is a property of this descent controller
   (growth 1.10, shrink 0.5, `min_step` 1e-12) and of this energy scale, and it is `UNVERIFIED` for
   any other. A lane putting a different law's relaxation on a float32 energy should re-measure the
   window rather than inherit this number. **This is the single most fragile thing in the entry**,
   and the PI's own words for it — *"수렴하는 것은 굉장히 어려운 것임"* — are the accurate summary.
4. **The float32 controls run on warp's CPU device, not on CUDA.** They drive `ALEPH-PORT-3603`'s
   real float32-term kernels, so the channel is not modelled — but a CUDA launch has a different
   accumulation order and `ALEPH-PORT-3601` measured atomic scatter costing 10 ULP where the CPU
   permutation bound predicted 5.3. The window in §10.2 is `UNVERIFIED` on a device.

   **One modelling attempt is recorded here as rejected rather than deleted**, because it gets the
   answer backwards and the next lane will be tempted by it: rounding the assembled **total** to
   float32 is *not* a model of a float32-term channel. A quantised total makes a sub-ULP step's
   energy difference exactly `0.0`, which a monotone test accepts, so no stall appears at all and a
   *loose* slack is what breaks it — the opposite conclusion, reached cheaply and confidently.
5. **No GPU, and no device figure anywhere in this entry.**
6. **The tension sweep is re-derived on this host, not on the workstation the published table used.**
   The baseline run recorded in §13.1 reproduces the published table row for row before any change,
   which is the check that the two hosts agree; without that check the comparison would be
   uninterpretable.
7. **`refuse_incomparable_position_modes` is a function callers must call.** It is not a decorator on
   every read path, because a stored artefact is a dict that no Aleph code owns the loading of. A
   caller that never compares artefacts is unaffected; a caller that compares them without calling it
   is not protected, and that is a limit of the shape rather than of the implementation.
8. **`ports/ledger/INDEX.md` is not regenerated by this lane.** `test_index_is_not_stale` was already
   red from another lane's uncommitted work before this entry existed. Reported, not fixed.
