# ALEPH-PORT-3609 — the step made device-resident, and the host crossings counted rather than claimed

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3609` |
| Lane | `46143f30` Lane G6 (CUDA track — residency) |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md` Track G row **G6**: *"`CellWorld` step on device — keep state resident, stop round-tripping to host per step."* The 02:36–03:01 CUDA run's `j4` artefact wrote its own caveat — *"only the FORCE is on a kernel here"*, the descent energy being frozen NumPy float64 on the host — so the physics ran on a device but the *step* did not stay there. And `ALEPH-PORT-3608` §14.3–4 says in as many words that its float32 descent window was measured on **warp's CPU device** and that a device re-measurement is this lane's first job. |

> **This entry carries no `decided_by` field and no agent may add one.** It proposes work and
> records measurements; it settles none of the open questions in
> `docs/decisions/PROPOSAL-position-precision-for-the-cuda-port.md`.

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered by it.

```python
# NEW — aleph/runtime/residency_kernels.py
from aleph.runtime.residency_kernels import (
    DescentStepVariant,      # TRUE + three wrong on purpose
    load,                    # lazy: `import warp` happens here, never at module scope
    loaded_warp,
)

# NEW — aleph/runtime/residency.py
from aleph.runtime.residency import (
    TransferCounts,          # frozen record of what crossed the PCIe boundary
    CountingBackend,         # a Backend proxy that counts every crossing it can see
    ResidentMembraneState,   # device-resident positions, forces and every accumulator
    ResidentStepReport,      # what one device-resident step did
    ResidentDescent,         # propose / evaluate / decide / commit / rollback, on device buffers
    resident_membrane_forces_and_energy,   # the context manager the vertical sweep runs under
    main,                    # `python -m aleph.runtime.residency <job>` — the measurement driver
)

# NEW — tests
tests/runtime/test_residency.py

# NEW — results
docs/results/2026-08-01-device-descent-window/README.md
```

**Not authorised and not touched:** `aleph/vertical/**` (in particular `relax.py`, which
`ALEPH-PORT-3608` landed in hours ago, and `membrane.py`, the frozen parity reference),
`aleph/scenarios/**`, `aleph/viz/**`, `aleph/runtime/law_kernels.py`, `aleph/runtime/law_cases.py`,
`aleph/runtime/parity.py`, `aleph/runtime/backend.py`, `aleph/runtime/warp_kernels.py`,
`aleph/state/**`, `tests/firewall/**`, `aleph/represent/**`, `aleph/learn/**`,
`scripts/gpu_preflight.py`, `ports/ledger/INDEX.md`.

### 1a. Why the integrator kernel is not in `law_kernels.py`

Every kernel in `aleph/runtime/law_kernels.py` is a **constitutive law** with a parity row in a
`law_cases.py` builder, and `ALEPH-PORT-3608` §8 controls C6–C8 iterate that registry and require
each case to feed its reference float64 positions, to be measurably sensitive to position precision,
and to declare a position mode. A position **update** — `x <- x + dt m F` — is an integrator. It has
no constitutive content, and registering it as a law case to satisfy a completeness guard would make
the guard's own statement less true, not more. It therefore lives in a module of this lane's own and
is graded here, at a **budget of exactly `0.0` ULP** against its NumPy reference (§7), which is a
sharper claim than any law case in the tree can make.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

Everything read is Aleph's own: `ALEPH-PORT-3601`/`-3603`/`-3604`/`-3605`/`-3606`/`-3607`/`-3608`,
`aleph/runtime/{backend,law_cases,law_kernels,parity,participant,transaction,warp_kernels}.py`,
`aleph/vertical/{relax,membrane,assembly}.py` (read, never edited), `aleph/scenarios/world.py`
(read only — another session owns that package), `scripts/cuda_jobs.py`, and the eight
`⚠ … LANDED` blocks of `docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md`.

## 3. Why source-derived porting beats clean-room

**It does not, and this entry says so rather than manufacturing a reason.** Device residency is a
property of *this* project's `Backend` protocol, `Transaction` and `Participant` contracts. The
provider's residency arrangement would answer a question about the provider's ownership model. The
entry is `RE-DERIVED` and is kept because `PLAN.md` §0.2.5 wants the decision recorded.

## 4. Physical or mathematical law represented

### 4a. The integrator, and it is scaffolding

`x_out = x_base + h F`, with `h = dt · mobility` [µm/pN] and `F` [pN]. This is one explicit
overdamped descent step, the same expression `HelfrichMembrane.accumulate` applies in the `SOLVE`
phase. `aleph/vertical/relax.py`'s standing disclaimers apply unchanged and are not restated here:
**not dynamics, not thermal, not implicit, not a minimiser with guarantees, and not a candidate
answer to `ALEPH-DQ-104`.** Putting it on a device changes none of them.

**The one arithmetic decision.** `x` is float64 (`POSITIONS_F64`, the mode `j4`'s PASS was measured
under) and `F` arrives in the ratified float32 compute channel, so the kernel widens the force to
float64 and accumulates there:

```
out[i] = base[i] + step * float64(force[i])
```

Each of the widen, the multiply and the add is a single correctly-rounded IEEE-754 operation, and
NumPy's `base + step * force.astype(np.float64)` is the same three operations in the same order. So
the kernel is expected **bit-identical**, not merely close, and §8 control C1 asserts that with `==`.

**The alternative, rejected in writing:** accumulating in float32 (`out = float64(float32(base) +
float32(step)*force)`) is `ALEPH-PORT-3601` §14.7's float32 position storage — the defect that
measured ~100 ULP on the tether — reintroduced inside the integrator where no law-parity gate would
ever look at it. It ships as the wrong variant `WRONG_F32_ACCUMULATE` rather than as a comment.

### 4b. What "resident" means here, stated as an invariant rather than as an intention

> **Residency invariant.** Between `push_positions` and `pull_positions`, no array whose element
> count scales with the mesh crosses the host boundary. A step reads back **scalars only**: the
> assembled energy, the largest residual force, and the device-side refusal flags.

This is not a performance claim and this entry makes none. It is an *auditable* claim, and §8
control C5 audits it by counting `warp.array.numpy` invocations and their sizes during a real step
rather than by reading the source.

### 4c. The defect class residency creates, which did not exist before

A non-resident evaluation allocates its accumulators fresh on every call, so they are zero by
construction. A resident one does not. **Every scatter destination and every reduction accumulator
must be explicitly re-zeroed at the top of each evaluation**, and an omission is invisible in a
single-step test and grows linearly with step count. This is the reason five of this lane's thirteen
mutants are stale-accumulator mutants (§13a M4–M6, M8, M10) and the reason the headline control is a
**trajectory** comparison rather than a single-step one.

### 4d. The oracle for residency itself: it must be a bit-identical refactor

Residency changes **where the buffers live**, not the arithmetic. So:

> **The resident evaluation must agree with the existing non-resident evaluation
> (`law_cases.membrane_kernel_forces`) BITWISE, on the same positions, at every step of a
> trajectory.**

That is an exact identity rather than a tolerance, it needs no budget, and it is the sharpest
statement this lane can make. `ALEPH-PORT-3605` §13a M9 named the limit that makes it necessary — a
parity gate cannot see a defect in the state it grades on — and a resident buffer *is* state that
survives between gradings.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| `positions` | µm | 1e-6 m | finite; float64 on device under `POSITIONS_F64` |
| `forces` | pN | 1e-12 N | finite; float32 compute channel |
| `step` (`dt · mobility`) | µm/pN | 1e6 m/N | `> 0`, finite |
| `energy` | pN·µm | 1e-18 J | finite; float64 accumulator over float32 terms |
| `max_residual_force_pn` | pN | 1e-12 N | `>= 0`; the **row-norm** maximum, not the component maximum |
| `TransferCounts.*` | count / elements | — | `>= 0`, monotone non-decreasing within a session |

Singular and boundary cases, each with the behaviour Aleph requires:

- **`rollback` without a `snapshot`.** Raises, matching `HelfrichMembrane.rollback`'s
  `AssertionError` rather than silently doing nothing.
- **`snapshot` twice without a decision.** The second overwrites the first; the shadow is always the
  configuration the next `rollback` restores, and it is asserted to be a distinct buffer from the
  live positions (the aliasing mutant M9).
- **`push_positions` with a different vertex count.** Refuses. A resident buffer whose shape no
  longer matches the mesh is not resizable in place, and silently reallocating would discard the
  shadow with the rollback obligation still outstanding.
- **A zero-length mesh.** Refused at construction; every reduction below would return a vacuous
  `0.0` and every control over it would pass on an empty set (`ALEPH-PORT-3606` §13a M5's shape).
- **`NumpyBackend`.** Supported, and every control that does not require a device runs on it, so the
  suite is CPU-testable per `PLAN.md` §0.1. `TransferCounts` are then all zero **by definition, not
  by measurement**, and the control that would be vacuous there is skipped by name rather than
  passed.
- **The descent tolerance window.** If the device measurement (§10.1) shows `2·u32` outside the
  converging window, this entry **reports and stops**; it does not tune.

Invariants that must hold, each with the test that asserts it:

- **I1.** The step kernel is bit-identical to `base + step * force.astype(float64)`. →
  `::test_the_device_step_is_bit_identical_to_the_numpy_update`
- **I2.** The resident evaluation is bit-identical to `membrane_kernel_forces` on the same positions,
  for every step of a 25-step trajectory. →
  `::test_a_resident_trajectory_is_bit_identical_to_the_non_resident_one`
- **I3.** `rollback` restores the pre-step positions bitwise (`np.array_equal`). →
  `::test_rollback_restores_the_device_positions_bitwise`
- **I4.** `commit` advances the shadow, so the next proposal builds on the accepted configuration. →
  `::test_commit_advances_the_baseline_and_rollback_after_commit_refuses`
- **I5.** A resident step crosses **no** mesh-sized array; every readback is a scalar or a reduction
  partial, and the count is measured from `warp.array.numpy` rather than modelled. →
  `::test_a_resident_step_reads_back_scalars_only`
- **I6.** Every accumulator is re-zeroed per evaluation: evaluating twice at the same positions gives
  bit-identical energy and force. → `::test_evaluating_twice_at_one_configuration_is_idempotent`
- **I7.** `max_residual_force_pn` is the maximum **row norm** and equals the host's
  `np.linalg.norm(f, axis=1).max()` of the same force array. →
  `::test_the_residual_is_a_row_norm_maximum_and_not_a_component_maximum`
- **I8.** The three wrong step kernels each move something the true one does not. →
  `::test_each_wrong_step_kernel_is_caught`
- **I9.** `TransferCounts` is a measurement: the modelled count equals the observed
  `warp.array.numpy` count and observed `wp.array(...)` constructions. →
  `::test_the_transfer_count_matches_what_the_device_actually_did`

## 6. Source evidence class and known retractions

No source symbol is claimed, so there is nothing to retract from `ffn_cellsim`. What this entry
inherits and must state:

* `docs/results/2026-07-31-tension-sweep/README.md` is `ANALYTIC_ORACLE`. Laplace's law is an
  identity of the model, not a measurement of a cell. Every σ here is `UNSOURCED` as physics.
* `ALEPH-PORT-3608` §10.2's float32 descent window is **warp CPU device**, and that entry marks it
  `UNVERIFIED` on a device. §10.1 below is this lane's device measurement of the same thing and
  **supersedes nothing**: a different device is a different measurement, not a correction.
* `ALEPH-PORT-3605` §13.4 and `-3606` §13.4 measured **opposite signs** for `POSITIONS_F64` versus
  `LOCAL_F32`. Both stand. This entry runs `POSITIONS_F64` because that is the mode `j4`'s PASS was
  measured under, and it resolves nothing.
* No timing figure in this entry is a performance claim. Residency was built to remove a
  correctness-relevant round trip, and wall-clock is reported as context only.

## 7. Independent oracle or derivation

Three, and they are of different kinds on purpose.

1. **Exact identity — the integrator.** IEEE-754 makes `base + step*float64(force)` a fixed sequence
   of three correctly-rounded operations. The kernel and NumPy must agree **bitwise**, budget `0.0`.
2. **Exact identity — the refactor.** §4d: the resident chain must reproduce the non-resident chain
   bitwise on identical positions. No tolerance is involved, so no budget can be chosen to make it
   pass.
3. **`ANALYTIC_ORACLE` — the Laplace sweep.** σ recovered across a 6× tension range and two
   resolutions, worst relative error, twelve configurations. `PLAN.md` §2.5: the relative error *is*
   the dilational virial of the residual, so it is a statement about where the relaxation stopped —
   which is exactly why moving the descent's energy channel onto the device is a thing it can fail.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| C1 | `tests/runtime/test_residency.py::test_the_device_step_is_bit_identical_to_the_numpy_update` | The step kernel equals `base + step*force.astype(f64)` with `==`, over 4,096 randomised sites, at three step sizes and at a 50 µm offset |
| C1b | `::test_the_step_kernel_widens_the_force_before_it_multiplies` | An input whose float32 rounding differs from its float64 value gives the float64 answer — so the widen is real and not incidental |
| C2 | `::test_a_resident_trajectory_is_bit_identical_to_the_non_resident_one` | 25 steps, `L2` mesh: energy and both force fields bitwise equal to `membrane_kernel_forces` at every step |
| C3 | `::test_rollback_restores_the_device_positions_bitwise` | `np.array_equal` on the pulled positions after propose+rollback, over 10 rejected steps |
| C3b | `::test_the_snapshot_is_a_distinct_buffer` | Writing the live positions does not change the shadow — the aliasing half, which a passing rollback cannot distinguish on its own |
| C4 | `::test_commit_advances_the_baseline_and_rollback_after_commit_refuses` | After `commit`, `rollback` raises; the next proposal starts from the committed configuration |
| C5 | `::test_a_resident_step_reads_back_scalars_only` | With `warp.array.numpy` and the `wp.array` constructor spied on: **zero** crossings of `3·V` elements during a step, and every observed readback `<= ceil(n/CHUNK)` float64 partials |
| C5b | `::test_the_transfer_count_matches_what_the_device_actually_did` | `TransferCounts` equals the spy's tally — the vacuity half, so the count is a measurement and not a claim |
| C6 | `::test_evaluating_twice_at_one_configuration_is_idempotent` | Second evaluation bitwise equal to the first, which is what catches an un-zeroed accumulator in one call rather than in a hundred |
| C7 | `::test_the_residual_is_a_row_norm_maximum_and_not_a_component_maximum` | Equals `np.linalg.norm(f, axis=1).max()`, and on a state where the two differ by more than 10 % so the assertion discriminates |
| C8 | `::test_each_wrong_step_kernel_is_caught` | All three wrong variants differ from the true one, and the channel each is caught in is named in the test |
| C9 | `::test_the_resident_descent_decreases_the_energy_and_the_predicate_sees_it` | Over 50 resident steps the accepted energies are non-increasing and at least one step was rejected — so the accept/reject path is exercised, not merely present |
| C10 | `docs/results/2026-08-01-device-descent-window/README.md` | The device descent-tolerance window, slack by slack, beside `ALEPH-PORT-3608` §10.2's CPU-device table |
| C11 | The Laplace sweep on `cuda:0` with the descent energy read from the device | 12/12 converged, worst relative error in the 1.5–1.6e-05 band |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| N1 (must fail) | `::test_a_stale_accumulator_is_visible_in_the_trajectory_and_not_in_one_step` | An evaluation that skips the re-zero agrees with the non-resident chain on step 1 and disagrees by step 2 — the statement that a single-step control **cannot** see this defect class |
| N2 (must fail) | `::test_each_wrong_step_kernel_is_caught` | `WRONG_SIGN` ascends the energy, `WRONG_F32_ACCUMULATE` differs in the last bits, `WRONG_DROPS_BASE` moves the body to the origin |
| N3 (must fail) | `::test_an_aliased_snapshot_makes_rollback_a_no_op` | With the shadow aliased to the live buffer, `rollback` returns positions that are **not** the pre-step ones — so C3 is protecting something |
| N4 (must fail) | `::test_an_uncounted_readback_is_visible_to_the_spy` | A deliberately uncounted `to_host` makes C5b's equality fail — so the counter cannot drift from the device silently |

## 10. Numerical and precision envelope

**Working precision.** Positions float64 on device (`PositionPrecision.POSITIONS_F64`); every
membrane force term float32 (`ALEPH-DQ-107`'s ratified compute channel); every energy and reduction
accumulator float64 over float32 terms; `ScatterMode.ORDERED`, measured bit-identical on the A5000
by `j2` (0.0 spread, 16 repeats, 24,576 contributions into 4,096 rows).

**Host for every device figure below:** `gbook`, NVIDIA RTX A5000 Laptop GPU, `sm_86`, warp 1.15.0,
under the PI's authorization record (host `GBook`, device 0, expiring `2026-08-04T00:00+09:00`), with
`scripts/gpu_preflight.py` run immediately before every job. **An agent may never write that record
and this lane did not.**

### 10.1 The device descent-tolerance window — THIS LANE'S FIRST TASK, AND THE WINDOW MOVED

> **`2·u32` — the value `ALEPH-PORT-3608` derived and the PI ratified — is OUTSIDE the converging
> window on CUDA for `L2_sigma20`. The window moved down by one binary order of magnitude. The
> intersection over the two probed configurations is `1·u32`, not `2·u32`.**

Measured on `gbook`'s A5000 with the identical isolation `-3608` §10.2 used: the membrane **energy**
from the float32-term kernels, the **forces** NumPy float64, `max_steps = 40,000`,
`u32 = 2^-24 = 5.9604644775390625e-08`.

| slack | `L1 σ=10` **CUDA** | `L1 σ=10` warp CPU (`-3608`) | `L2 σ=20` **CUDA** | `L2 σ=20` warp CPU (`-3608`) |
|---|---|---|---|---|
| `0` (exact) | stall, 507 acc, step 6.1e-13 | stall, 442, 6.4e-13 | stall, 600 acc, step 5.3e-13 | stall, 296, 6.1e-13 |
| **`1·u32`** | **conv, 9,460 acc, 2.007e-06** | **STALL, 979 acc** | **conv, 3,113 acc, 2.838e-07** | conv, 2,806, 2.13e-05 |
| **`2·u32`** (derived) | conv, 3,938 acc, 1.664e-05 | conv, 2,742, 3.30e-06 | **WANDER, 35,166 acc, 5.323e-03** | **conv, 4,610, 3.24e-05** |
| `4·u32` | wander, 4.225e-03 | wander, 6.18e-03 | wander | wander, 1.08e-02 |
| `6·u32` | wander, 3.100e-03 | wander, 2.25e-02 | wander | wander, 6.33e-03 |
| `8·u32` | wander, 1.650e-02 | wander, 1.99e-02 | wander | wander, 2.60e-02 |
| `12·u32` | wander, 2.319e-02 | wander, 2.42e-02 | wander | wander, 2.82e-02 |
| `24·u32` | wander, 1.085e-02 | wander, 3.92e-02 | wander | wander, 4.50e-02 |

**Both edges moved down by one binary order**, and they moved together:

| | warp CPU device (`-3608`) | **CUDA `sm_86`** (this lane) |
|---|---|---|
| converging multiples, `L1 σ=10` | `{2}` | `{1, 2}` |
| converging multiples, `L2 σ=20` | `{1, 2}` | **`{1}`** |
| **intersection** | **`{2}`** | **`{1}`** |

`-3608` §14.3 called the one-binary-order width *"the single most fragile thing in the entry"* and
said in as many words that a lane putting a different law's relaxation on a float32 energy must
re-measure. It turns out that changing the **device** is enough.

#### The mechanism, and it is not the reduction

The two-stage reduction is **not** where this comes from. `warp_kernels.CHUNK` is a constant, so the
decomposition — and therefore the exact sequence of float64 additions — is identical on both devices
by construction (`warp_kernels.py`'s own docstring says so, and `j2` measured `ORDERED` scatter
bit-identical on this A5000). What differs is the **float32 terms**: CUDA contracts multiply-add
pairs into FMA and implements `sqrt`, `rsqrt` and division differently from the host's SSE/NEON
float32, so every `wp.length`, `wp.cross` and normalisation in the membrane chain lands on different
last bits. The assembled energy's *noise* is therefore a different number on the two devices, and
`-3608` §4b's insight is what turns that into a moved window: **a tolerance is simultaneously a
bound on the evaluation noise and the uphill wander the descent is permitted**, so a quieter channel
needs less slack *and* tolerates less before it climbs. Both edges move the same way, together, which
is what the table shows.

#### What this lane did about it: nothing, deliberately

`DESCENT_SLACK_ROUNDOFFS = 2` is in `aleph/vertical/relax.py`, which this entry's §1 does not
authorise and which `ALEPH-PORT-3608` derived and the PI ratified at 2026-08-01 01:39 KST. **It is
not changed here, and the Laplace sweep in §10.3 was run at the unchanged `2·u32`** so that the
oracle reports the consequence rather than being tuned past it. A derivation and a device
measurement now disagree, which `CLAUDE.md` §2 rule 5 makes a finding and not a licence. It is
raised in §14.6 for the PI.

### 10.2 What crosses the host boundary per step, before and after

Measured on the A5000 by `python -m aleph.runtime.residency transfers --host GBook --repeats 64`,
and **identically** on warp's CPU device — the counts are a property of the code path, not of the
hardware, which is why the control that asserts them is in the CPU-testable suite. Level-2 icosphere,
162 vertices / 320 faces, `POSITIONS_F64`, `ORDERED`, CSR cache warm.

| | calls | **elements** | what they are |
|---|---|---|---|
| **before** — one `membrane_kernel_forces` evaluation | 6 | **1,947** | 2 × 486 float64 positions up, 2 × 486 float32 forces down, 2 × ~1.5 float64 reduction partials down |
| … plus what the counter cannot see | 2 | **1,920** | `law_cases._device_triangles` uploads connectivity with `wp.array`, not `backend.array`, and says why. **So the before figure is a LOWER BOUND of 3,867 elements.** |
| **after** — one resident step (`propose` + `evaluate` + residual) | **3** | **4** | three float64 reduction partial reads; **zero** arrays |

**The headline is the zero, not the ratio.** `bulk_to_device_calls` and `bulk_to_host_calls` are
exactly `0` per resident step, asserted by
`::test_a_resident_step_reads_back_scalars_only` against a spy on `warp.array.numpy` and the
`wp.array` constructor rather than against the module's own bookkeeping. The four elements that do
cross are `ceil(320/256) + ceil(162/256) + ceil(162/256) = 2 + 1 + 1` float64 reduction partials —
the assembled areal energy, the assembled bending energy, and the largest residual force.

The relaxation makes **two to four** `membrane_kernel_forces` calls per step (the pre-step energy,
the SOLVE-phase force, the predicate's post-step energy, and `max_residual_force_pn`), so the
per-step *before* figure is 2–4× the row above.

**And the sweep's own figure, which is the honest one for the twelve configurations**, because there
the cortex, the osmotic envelope and the two connectors are still host NumPy and the positions
genuinely have to cross:

| | per accepted step |
|---|---|
| bulk crossings, `L1` | ~1.3 up + ~1.3 down, ~164 + ~161 elements |
| bulk crossings, `L2` | ~1.3 up + ~1.3 down, ~625 + ~614 elements |
| device evaluations per step | **~1.13**, because the context caches on exact position contents |

That is the number §14.2 is about: the membrane half of the sweep's step is resident and the rest is
not, and 1.3 crossings per step is what "the rest is not" costs. A fully resident step is the row
above it, and it is `0`.

### 10.3 The Laplace sweep with the descent energy on the device — IT DOES NOT REPRODUCE

`j4` reproduced the sweep on CUDA at **1.5748e-05, 12/12**, and wrote its own caveat: *"only the
FORCE is on a kernel here"*, the descent energy being frozen NumPy float64 on the host. G6 moved that
third channel onto the device. **The oracle does not survive the move**, at either slack measured:

| arm | descent energy | descent slack | converged | worst rel. err over converged |
|---|---|---|---|---|
| published NumPy | host float64 | `1e-12`, then `2·u64` | 12/12 | 1.566e-05 |
| `-3603` kernels, warp CPU | host float64 | `2·u64` | 12/12 | 1.5610e-05 |
| **`j4`, CUDA** | **host float64** | `2·u64` | **12/12** | **1.5748e-05** |
| **this lane, CUDA** | **DEVICE float32-term** | `2·u32` (derived, ratified) | **2 / 12** | 2.2429e-06 |
| **this lane, CUDA** | **DEVICE float32-term** | `1·u32` (the device window, §10.1) | **8 / 12** | **6.9551e-05** |

Per configuration, `2·u32` / `1·u32`:

| | `2·u32` | `1·u32` |
|---|---|---|
| `L1_sigma10` | wander, 3.4488e-03 | **conv, 1.1422e-06** |
| `L1_sigma20` | **conv, 2.2429e-06** | **conv, 1.2163e-06** |
| `L1_sigma30` | wander, 1.0400e-02 | **conv, 1.8771e-06** |
| `L1_sigma40` | wander, 3.8229e-03 | **conv, 2.2149e-06** |
| `L1_sigma50` | **conv, 1.9989e-06** | **conv, 1.7499e-06** |
| `L1_sigma60` | wander, 1.0820e-02 | **conv, 7.2397e-07** |
| `L2_sigma10` | wander, 9.6228e-03 | **conv, 6.9551e-05** |
| `L2_sigma20` | wander, 3.3573e-03 | **conv, 3.3424e-05** |
| `L2_sigma30` | wander, 4.3010e-03 | wander, 1.6872e-03 |
| `L2_sigma40` | wander, 4.1604e-03 | wander, 4.5592e-04 |
| `L2_sigma50` | wander, 8.6208e-03 | wander, 5.4970e-03 |
| `L2_sigma60` | wander, 1.6520e-03 | wander, 3.6488e-03 |

A non-converged row's recovered tension is **not a datum** and is printed only so the failure mode is
visible (`scripts/sweep_vertical_parameters.py`'s rule: a table cannot show you that one of its rows
stopped early).

#### What the two arms say, and what they do not

1. **No single slack converges all twelve.** `2·u32` converges the two `L1` rows whose window happens
   to contain it; `1·u32` converges every `L1` row and the two coarsest `L2` ones. The window is
   **per configuration**, and §10.1 measured that directly on two of these twelve.
2. **`L1` is comfortably inside the window at `1·u32`** — six of six converge at 7e-07 … 2e-06,
   *better* than the published 1.566e-05. That is not an improvement in the model; it is a different
   stopping point (§4d of `-3608`: the Laplace error **is** the dilational virial of the residual).
3. **`L2` is the half that does not reproduce.** Four of six wander at both slacks, and the two that
   converge land at 3.3e-05 and 7.0e-05 — a factor of 2–4 outside the 1.5–1.6e-05 band. **This lane
   did not go looking for a slack that would put them inside it**, which is what the G6 brief means
   by not tuning to reach the number.
4. **The `1·u32` arm is a measurement of §10.1's consequence, not a proposal.** It is reported beside
   the derived arm rather than instead of it, and `DESCENT_SLACK_ROUNDOFFS` is unchanged.

#### What this is NOT evidence of, stated because the natural reading is wrong

It is **not** evidence that the residency layer or the kernels regressed. The evidence against that
reading is in this entry's own §7 and §8: the resident chain reproduces the non-resident chain
**bitwise per force field over a 25-step trajectory**, and `j4` — the same kernels, the same
`POSITIONS_F64`, the same device, with the descent energy on the *host* — still gives 12/12 at
1.5748e-05. The one channel that changed is the one `ALEPH-PORT-3603` §13 already identified as the
one that breaks a descent, and `-3608` §14.3 already called its tolerance *"the single most fragile
thing in the entry"*. **G6's contribution here is the measurement that says how fragile: on this
device, one binary order of slack is the difference between 2 of 12 and 8 of 12, and neither is 12.**

## 11. Production-backend residency and transfer

This section is the entry's subject rather than a footnote to it.

**Allocated through the backend and staying there:** positions (float64), the position shadow, the
candidate positions, the assembled force, and every intermediate of the areal and bending chains —
face edges, per-face energy, corner forces, corner areas, cotangents, area gradients, cotangent
gradients, curvature contributions, the lumped vertex area, the vertex energy and the two vertex
coefficient fields. Connectivity (`triangles`, the corner and curvature scatter maps) is uploaded
**once** at construction; `ALEPH-PORT-3603` §11 requires scatter destinations to come from
`backend.zeros` and they do.

**What still crosses, and why:** the assembled energy (one float64), the largest residual force (one
float64), and — only when a host owner needs them — the positions and forces. In the twelve-
configuration Laplace sweep the cortex, the osmotic envelope and the two connectors are still host
NumPy, so that sweep's step is **not** fully resident and §14.1 says so rather than letting the
headline imply it.

**What an agent may never do here:** write `~/.aleph_data/gpu_authorization.json`. It is read by the
preflight and by `WarpBackend.__post_init__`, and never written.

## 12. Comments and docstrings to discard

Nothing is ported, so no provider prose crosses. Nothing is removed from Aleph's own text either:
this lane adds two modules and edits none.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **The residency layer landed and is graded by an exact identity; the sweep oracle did not survive moving the descent energy onto the device, and that is reported rather than tuned past.** `tests/runtime/test_residency.py` **32 passed**; `tests/runtime` **465 passed / 3 skipped** (was 416/3 after G5, 465/3 with `-3608`); `tests/runtime` + `tests/vertical` **1,841 passed / 3 skipped**; everything else **2,996 passed / 1 skipped / 2 failed** (both `tests/ports` INDEX guards, already red before this lane's first commit — §14.8). Mutation study **13 planted / 13 killed**, second pass from a verified-empty bytecode tree |
| Residency, measured | **0 bulk crossings per resident step**; 3 readbacks of 4 float64 elements, against 6 crossings / ≥3,867 elements per non-resident evaluation (§10.2), audited against a spy on `warp.array.numpy` |
| The refactor identity | resident chain **bitwise equal** to `law_cases`' non-resident chain, per force field, over a 25-step trajectory |
| The integrator | **bitwise equal** to `base + step*float64(force)` at a declared budget of **0.0 ULP**, at 4,096 sites, three step sizes, at the origin and at a 50 µm offset |
| The oracle | **2 / 12** at the derived `2·u32`; **8 / 12** at `1·u32`; worst over converged 2.2429e-06 and 6.9551e-05 respectively. `j4` — same kernels, same device, descent energy on the **host** — is 12/12 at 1.5748e-05 and is untouched |
| Reviewer | **agent-proposed; unratified.** No `decided_by` field exists in this entry and no agent may add one |
| Rollback | Delete `aleph/runtime/residency.py`, `aleph/runtime/residency_kernels.py` and `tests/runtime/test_residency.py`. What breaks: nothing that runs today — no existing module imports them, `aleph/vertical/` and `aleph/scenarios/` are byte-for-byte untouched, and `j4`'s path is unchanged. What is lost is §10.1's measurement, which is the part of this entry that matters most |

## 14. Honest limits

1. **This is not "`CellWorld` on device".** `aleph/scenarios/world.py` composes host NumPy owners
   through `OwnerSlot`/`ConnectorSlot`, and making its step resident means giving each owner device
   state — an edit inside `aleph/scenarios/**` and `aleph/vertical/**`, which another session owns
   and which this lane is forbidden to make. What lands here is the **residency layer** in
   `aleph/runtime/`, where a runtime facility belongs, plus the membrane demonstration and the
   sweep. The remaining owners are named in §11, not glossed.
2. **The Laplace sweep's step is partly resident.** The membrane half is; the cortex, turgor and the
   two connectors are host NumPy and the positions therefore still cross once per evaluation. §10.2
   reports both numbers separately for exactly this reason.
3. **A device measurement is one device.** Every figure is `sm_86` on one A5000 under warp 1.15.0.
   `ALEPH-PORT-3601` measured atomic scatter costing 10 ULP where the CPU permutation bound predicted
   5.3; nothing here should be read as a bound over hardware.
4. **The residency invariant is audited, not enforced.** §8 C5 counts what a step did. Nothing stops
   a later edit from adding a `to_host` inside an evaluation; the control would go red, which is the
   guard working, but it is a test and not a type.
5. **No performance claim.** Wall-clock is reported as context. A 162-vertex mesh is latency-bound on
   an A5000 and the interesting figure would be a mesh two orders larger, which this lane does not
   run.
6. **A DERIVATION AND A DEVICE MEASUREMENT NOW DISAGREE, AND THIS LANE DOES NOT RESOLVE IT.**
   `ALEPH-PORT-3608` §4b derives `DESCENT_SLACK_ROUNDOFFS = 2` — two evaluations, one unit round-off
   each — and measured it as the only converging point on warp's CPU device. §10.1 above measures
   that on CUDA `2·u32` **wanders** on `L2_sigma20` and `1·u32` converges on both probed
   configurations. Three things follow, and only the first two are this lane's to state:

   - The derivation is not wrong about what it derives. It bounds the *evaluation* error, and the
     evaluation error is smaller on CUDA because the float32 terms are contracted differently. What
     it does not bound is the *uphill wander*, which `-3608` §4b itself says pulls the other way.
   - **A per-device constant would be the wrong repair** and this lane proposes none. The window is
     a property of the controller, the energy scale, the law, *and* now the device; freezing four of
     those into one integer is how `1.0e-12` came to be a round number nobody could derive.
   - Which of "re-derive per device", "make the slack adaptive", or "leave `2` and accept that a
     float32 device descent needs a different controller" is right is **above a porting lane**, and
     it touches `aleph/vertical/relax.py`, which `-3608` owns and this entry does not. Raised, not
     taken.

   The Laplace sweep in §10.3 was therefore run at the **unchanged** `2·u32`, which is why its
   result is what it is. Reading §10.3 without §10.1 would make the sweep look like a regression in
   the port; it is a measurement of a tolerance that no longer sits inside its window on this device.
7. **Only two configurations were probed** for the window, the same two `-3608` used. `L1 σ=10` still
   converges at `2·u32` on CUDA; it is `L2 σ=20` that does not. So "the window moved" is established
   for the *intersection*, and the per-configuration picture is that the finer mesh is the one that
   moved. Whether that is about the mesh, the tension, or their product is not measured here — though
   §10.3's twelve rows say the same thing at a coarser grain: **every `L1` row converges at `1·u32`
   and four of six `L2` rows converge at neither slack.**
8. **`ports/ledger/INDEX.md` is red and it is foreign.**
   `test_the_index_is_not_missing_a_tracked_entry` and `test_index_is_not_stale` were **already red
   before this lane's first commit**, naming `-3606`, `-3607` and `-3608`. Committing `-3609` adds a
   fourth *name* to the same already-red assertions; it does not add a third red test.
   `test_named_controls_resolve_to_real_tests` passes — every control named here resolves to a real
   test. **Not regenerated**, on instruction: the generator reads the working tree and another lane
   still holds uncommitted `-3405`/`-3406`.
9. **The window and the sweep were measured with two different force channels, and the difference
   shows.** §10.1 holds the forces at NumPy float64 — `-3608` §10.2's isolation, so that only the
   energy channel is under test — while §10.3 puts *both* on the device, which is what a G6 lane is
   for. `L1_sigma10` converges at `2·u32` in the first and wanders in the second. So the window
   measured on the energy channel alone is an **upper bound** on the window of the full device path,
   and this entry does not measure the latter's edges directly. Naming it rather than leaving the two
   tables to be read as one.
