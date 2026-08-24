# ALEPH-PORT-3601 — a parity harness for force laws, not array operations

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3601` |
| Lane | `46143f30` Lane G1 (CUDA track) |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `ALEPH-DQ-107` ratified Warp-on-CUDA as the backend and **no line of physics has ever reached a kernel**. `aleph/runtime/parity.py` compares array operations; nothing compares a *law*. A CUDA kernel with no independent reference is unfalsifiable — a fast wrong answer and a fast right one are the same observation. |

---

## 1. Aleph API

```python
from aleph.runtime.parity import (
    ATOMIC_SCATTER_ULP_FLOOR,   # 10.00, measured on an A5000; see §10
    LawCase,                    # one law, one state, one declared ULP budget per output field
    LawFieldParity,             # the measurement for ONE output field of ONE case
    LawParity,                  # every field of one case, plus the verdict
    LawParityReport,            # every case, plus the backend context it was measured in
    compare_law_case,           # drive reference and candidate on identical state, measure
    law_parity_report,          # the gate
    law_report_as_dict,         # flatten for a run record
    ulp_at_scale,               # the ULP measure, defined in §4 and stated in the report header
)

from aleph.runtime.law_kernels import (
    TetherVariant,              # TRUE | WRONG_ENERGY_HALF | WRONG_UNILATERAL_GATE
    load,                       # compile the law kernels; lazy, mirrors warp_kernels.load()
    loaded_warp,                # the warp module or None; never imports anything
)

from aleph.runtime.law_cases import (
    tether_law_case,            # the ErmTether law as a LawCase, at a chosen kernel variant
    tether_state,               # reproducible membrane/cortex site geometry
)
```

Nothing outside this list is covered. In particular this entry does **not** authorise any change to
`aleph/runtime/backend.py`, to `aleph/runtime/warp_kernels.py`, or to any module under
`aleph/vertical/` — the NumPy laws are the frozen reference and are read, never edited.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

One structural observation about that project is already in Aleph's own record
(`docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md` §1.3, written by a prior read-only audit lane):
it keeps a pure-NumPy discrete oracle beside the Warp kernel it checks. **That is a shape, not code**,
and Aleph already has the shape thirty times over — every wired connector is a NumPy law with no
kernel beside it. This entry inverts the missing half: it writes the comparison, and the law that
already exists becomes the oracle. No line, no constant and no branch structure was taken.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** The harness is determined by Aleph's own surfaces:
`aleph.runtime.parity` fixes the vocabulary (`OpCase` → `LawCase`, `compare_case` →
`compare_law_case`, `parity_report` → `law_parity_report`), `aleph.runtime.backend.Backend` fixes
what a kernel driver may call, and `aleph.vertical.connectors.site_pair_forces` fixes the reference
law bit for bit. There is nothing here a second implementation would need help with.

## 4. Physical or mathematical law represented

Two things are represented and they must not be conflated.

**(a) The law under test — the ERM tether, re-derived.** `N` site pairs. Site `i` couples a membrane
point `a_i` to a cortex point `b_i`. Write `s_i = |a_i - b_i|` for the separation and `s_0` for the
rest length. The element is *unilateral in tension*: it stores energy only when stretched.

```
x_i  = max(s_i - s_0, 0)                     engaged stretch [um]
U    = sum_i (k/2) x_i^2                     stored energy [pN.um]
f_i^(a) = -(dU/ds_i) (a_i - b_i)/s_i = -k x_i (a_i - b_i)/s_i     [pN]
f_i^(b) = -f_i^(a)
```

`f^(a)` points from `a_i` toward `b_i` whenever `x_i > 0`: the tether pulls. On the inactive branch
`x_i = 0` and both the energy and the force are **exactly zero** — the number, not a small number.
That exact zero is a bit-level property, which is why it is worth a kernel-level gate at all: a
kernel that loses it still descends, still reports a stiffness, and turns a tether into a strut.

The assembled outputs are the scalar `U` and two per-vertex force arrays obtained by scattering the
per-site forces through index maps `p: site -> vertex` (`f^A_v = sum_{i: p_a(i)=v} f_i^(a)`, likewise
for B). Several sites may anchor on one vertex, so the scatter genuinely accumulates.

**(b) The measure — float32 ULP, normalised by the field's own scale.** For a reference field `r` and
a candidate field `c`:

```
ulp(c, r) = max|c - r| / (eps32 * max|r|),      eps32 = 1.1920929e-07
```

`ulp = 0` exactly when the two agree bitwise after the float64 widening; `ulp = 1` means they differ
by one float32 unit in the last place *of the field's largest entry*.

The denominator is the field's largest magnitude and **not** the element's own value. That choice is
forced, and it is the same choice `parity.py` already made and stated for `max_rel`: an elementwise
ULP distance is unusable on a unilateral force array, because a correct inactive site holds an exact
zero and any disagreement there — including a denormal — reports as an unbounded ULP count while
every physically meaningful digit is intact. A tether at rest would score infinity. Normalising by
the field's scale asks the question actually being asked: *how large is the disagreement compared
with the forces involved.*

It is also the definition the project's existing measurements are already in.
`ScatterDeterminism.ulp_spread` in `parity.py` is `max_rel_spread / EPS32`, which is this same
quantity, so the `0.00` / `10.00` figures in §10 are directly comparable to what this harness prints.
Adopting a second ULP convention would have made the two incomparable, silently.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| site positions `a`, `b` | µm | m | finite; `a_i != b_i` |
| separation `s`, rest length `s_0` | µm | m | `s > 0`, `s_0 > 0` |
| stiffness `k` | pN/µm | N/m | `k >= 0` finite |
| energy `U` | pN·µm | J | `>= 0` |
| force `f` | pN | N | finite |
| ULP figure | float32 ULP of the field's scale | dimensionless | `>= 0`; `inf` if the reference field is identically zero and the candidate is not |

Singular and boundary cases, each with the behaviour Aleph requires:

- **Coincident pair (`s_i = 0`).** The reference refuses with `ConnectorGeometryError`. This harness
  never constructs one and does not soften it; a state that would produce one is a state the case
  builder rejects before either side runs.
- **A reference field that is identically zero** (for example a fully slack tether). `ulp_at_scale`
  returns `0.0` if the candidate is also exactly zero and `inf` otherwise. `inf` exceeds every finite
  budget, so a kernel that invents force out of a dead law fails rather than dividing by zero.
- **Empty state (`N = 0`).** Not supported; the case builder refuses. A gate that passes vacuously on
  no sites is the failure mode this whole entry exists to prevent.
- **A candidate that returns a different set of output fields, or a field of a different shape.**
  Raised as `ValueError`, never graded. A contract break is not a precision number.

Invariants that must hold, each with the test that asserts it:

- **I1.** A law compared with itself is `0.00` ULP on every field, exactly —
  `test_a_law_compared_with_itself_is_exactly_zero_ulp`.
- **I2.** The verdict is the declared budget and nothing else: a budget of `0.0` fails the *true*
  kernel — `test_a_zero_ulp_budget_fails_the_true_kernel`.
- **I3.** Energy and forces are judged separately, so a kernel wrong in only one of them is caught in
  only that one — `test_the_wrong_energy_kernel_leaves_the_forces_bit_identical`.
- **I4.** The scatter destination is allocated through the backend, never as a host ndarray —
  `test_the_scatter_destination_is_allocated_through_the_backend`.
- **I5.** The report names the scatter mode, the device and the compute dtype it measured in —
  `test_the_report_names_the_scatter_mode_and_the_device`.
- **I6.** The atomic floor widens the budget of scatter-assembled fields only, never the energy —
  `test_the_atomic_floor_widens_only_the_scatter_assembled_fields`.

## 6. Source evidence class and known retractions

Nothing is inherited, so there is no inherited evidence class to carry. What *is* carried, and where
it was checked:

- **`docs/design/GPU_STATE_2026-07-31.md` and commit `636b0c8`** — the `ORDERED` = 0.00 ULP /
  `ATOMIC` = 10.00 ULP measurement on the A5000, against a CPU-derived permutation bound of ≤ 5.30.
  Read in this repository, not quoted from elsewhere.
- **Commit `636b0c8` is itself a retraction** and it is the reason for §11's warning:
  `278e6e5` reported a GPU-only defect that turned out to be a misuse of `WarpBackend.scatter_add`
  (a raw host ndarray passed as `dest`, which warp's CPU device tolerates and CUDA does not).
  The API was right; the caller was wrong. This entry's kernel driver allocates every scatter
  destination through `backend.zeros(...)`, and `test_the_scatter_destination_is_allocated_through_the_backend`
  is the control that keeps it that way.
- **No retraction was searched for in `/Users/sw1/ffn_cellsim`,** because nothing from that tree is
  used. `STATE.md` §(c)'s not-quotable list is therefore not engaged: no magnitude from that project
  appears in this entry or in the code it authorises.

## 7. Independent oracle or derivation

The oracle is **Aleph's own frozen NumPy law**: `aleph.vertical.connectors.site_pair_forces` driven by
`UnilateralSpring(engages_in_tension=True, law=ContactLaw.LINEAR)`, in float64, on the host, with no
backend involved. It is not the source project, it is not a second copy of the kernel, and it is not
modified by this entry. It already carries its own controls under `tests/vertical/`.

That is one oracle. It checks the *kernel*. It cannot check the *harness*, because a harness that
always returns "agree" would pass against any oracle. So there is a second, and it is the one this
task is actually about:

**The negative control is a pair of deliberately wrong kernels**, shipped in the package beside the
true one and driven through the same code path as production. See §9.

A third, weaker check: the analytic value of the law at a hand-computed configuration. `U` for a
single pair at `k = 200 pN/µm`, `s = 0.05 µm`, `s_0 = 0.02 µm` is `0.5 * 200 * 0.03^2 = 0.09 pN·µm`
and `|f| = 200 * 0.03 = 6 pN`. Used to confirm the reference is the law §4 derives, so that the
parity number is a statement about the kernel and not about two agreeing mistakes.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/runtime/test_law_parity.py::test_the_true_tether_kernel_agrees_with_the_frozen_numpy_law` | The float32 Warp tether kernel, driven on warp's CPU device with `ScatterMode.ORDERED`, agrees with the frozen NumPy law on **every** output field within that field's declared ULP budget, on a state where 40 of 96 sites are slack. Reports the measured figure for energy and for both force arrays. |
| Positive | `tests/runtime/test_law_parity.py::test_a_law_compared_with_itself_is_exactly_zero_ulp` | The harness's zero point: the NumPy law compared with itself is `0.00` ULP on every field, exactly. If this is not exact the harness has a defect of its own and every other number here is uninterpretable. |
| Positive | `tests/runtime/test_law_parity.py::test_the_report_names_the_scatter_mode_and_the_device` | The report carries `scatter_mode`, `device` and `compute_dtype`, so a stored measurement can never be read without knowing which accumulation path produced it. |

## 9. Deliberately failing negative control

Two wrong kernels ship in `aleph/runtime/law_kernels.py`. Both are *plausible* hand-translation
errors, not sabotage, and each is caught by a different output field — which is the argument for
grading energy and forces separately rather than reporting one number.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/runtime/test_law_parity.py::test_a_kernel_that_dropped_the_half_in_its_energy_is_caught` | `TetherVariant.WRONG_ENERGY_HALF` writes `k x^2` where the law says `(k/2) x^2`. Its **forces are bit-identical to the true kernel**. The harness must fail it, and must fail it on the `energy` field specifically. A gate that graded only the forces — which is what most parity gates grade, because forces are what the integrator consumes — would certify it. |
| Negative (must fail) | `tests/runtime/test_law_parity.py::test_a_kernel_that_dropped_the_unilateral_gate_is_caught` | `TetherVariant.WRONG_UNILATERAL_GATE` omits `x = max(x, 0)`, so a compressed tether pushes. This is the tether-becomes-a-strut defect the connector's own docstring names, one missing line, and it is invisible in a plot of a taut mesh. The harness must fail it on the force fields. |
| Negative (must fail) | `tests/runtime/test_law_parity.py::test_a_zero_ulp_budget_fails_the_true_kernel` | The vacuity guard on the verdict itself: with `ulp_budget` set to `0.0`, the *correct* kernel fails, because float32 storage is not bitwise float64. A verdict that survived a zero budget would not be reading the budget. |
| Negative (must fail) | `tests/runtime/test_law_parity.py::test_both_wrong_kernels_are_caught_in_atomic_mode_too` | Neither mutant escapes through the widened `ATOMIC` budget. A floor of 10 ULP is a floor, not an amnesty. |

And the finding that the negative controls produce, recorded as a control because it is the honest
limit of the method:

| Control | Location | Asserts |
|---|---|---|
| Blind spot (must hold) | `tests/runtime/test_law_parity.py::test_the_dropped_gate_is_invisible_on_an_all_taut_configuration` | On a state where **every** site is stretched, `WRONG_UNILATERAL_GATE` is bit-identical to the true kernel and the harness passes it at `0.00` ULP. The clamp never fires, so there is nothing to disagree about. **A law-parity case is only as strong as the branches its state visits**, and this is asserted rather than left for someone to discover on a kernel that mattered. |

## 10. Numerical and precision envelope

**Working precision.** Positions, per-site forces and the assembled force arrays are float32 on the
device (`ALEPH-DQ-107`'s compute channel). Per-site energies are widened to float64 *inside the
kernel* before the product is formed, and the total is reduced by `Backend.sum`, whose two-stage
float64 accumulation is already in `warp_kernels.py`. So the energy carries one rounding per input,
not one per addition. The reference is float64 end to end.

**Where the disagreement comes from, derived rather than assumed.** The dominant term is not the
arithmetic, it is that each position was rounded to float32 once before anything happened, and the
law then differences two nearby positions. Writing `u = eps32/2` for the unit round-off:

```
|d(a_i - b_i)| <~ u (|a_i| + |b_i|)            the difference inherits both roundings
|ds_i|         <~ u (|a_i| + |b_i|)
|dx_i|         =  |ds_i|                        x = s - s0 subtracts an exact constant
```

so, propagating into the two graded fields and normalising each by its own scale,

```
amp_energy = 2 * sum_i x_i (|a_i| + |b_i|) / sum_i x_i^2
amp_force  = max_i [ (|a_i| + |b_i|) (1 + x_i / s_i) ] / max_i x_i
bound_ulp  = (amp + 1) / 2
```

The `+1` covers the final rounding of the result back into float32, and the `/2` converts
`parity.py`'s relative bound `u (amp + 1)` into ULP, since `u = eps32/2`. `amp` is computed **from
the host state by the case**, so it is a property of the data and not a number chosen after seeing
the answer — the same discipline `parity.py` §"The round-off bound" already states for array ops.

**This is a large number for this law, and that is the finding.** A membrane vertex sits ~5 µm from
the origin and its tether is ~0.03 µm long, so the difference cancels about 300-fold before the
subtraction of `s_0` cancels it again. Float32 positions cannot express a 0.03 µm coordinate
difference at 5 µm to better than ~1e-5 relative. The measured figure and the declared budget for
each case are recorded in §13 after the run; the budget is set from the measurement with a stated
margin, and the report prints `bound_ulp` beside it so a reader can see whether the declared budget
is inside what round-off explains or wider than it. A budget wider than `bound_ulp` is flagged in the
rendered report, because a budget that exceeds its own round-off explanation can hide a defect.

**The consequence for the CUDA port, stated because it is the reason to measure at all:** for this
law at absolute cell coordinates, float32 position storage dominates the `ORDERED`-vs-`ATOMIC`
question by more than two orders of magnitude. The scatter mode is not what limits a tether's
accuracy; the coordinate representation is. That is a design question for G3 and is raised, not
resolved, here.

**The atomic floor.** `ATOMIC_SCATTER_ULP_FLOOR = 10.0`. Measured on the A5000 (`636b0c8`,
`docs/design/GPU_STATE_2026-07-31.md`): `ORDERED` scatter is bit-identical across 16 launches
(0.00 ULP) while `ATOMIC` measured 10.00 ULP against a CPU-derived permutation bound of ≤ 5.30 — the
CPU bound underestimates real CUDA atomic non-determinism by about 2×. When the candidate backend's
scatter mode is `ATOMIC`, the effective budget of a **scatter-assembled** field is
`max(declared, 10.0)`; the energy is not scatter-assembled and keeps its declared budget. The report
records which of the two supplied the effective number, so a widened budget can never be mistaken for
a declared one.

**Outside the envelope.** One device, one architecture, and `ATOMIC` has never been measured on this
Mac's CUDA-less build — `measure_scatter_determinism` already refuses to let a CPU zero be read as
evidence about CUDA, and this entry inherits that caveat verbatim rather than restating it weaker.

## 11. Production-backend residency and transfer

The kernel driver is production code and runs on the device. Per case:

| Array | Where | Precision | Transfer |
|---|---|---|---|
| site positions `a`, `b` | device | float32 | uploaded once per case by `backend.array` |
| per-site energy | device | float64 | never leaves the device before `backend.sum` reduces it |
| per-site forces | device | float32 | never leaves the device |
| assembled vertex forces | device, **allocated by `backend.zeros`** | float32 | read to host once, by the harness, for the comparison |

**The one rule this entry exists to keep.** `WarpBackend.scatter_add` takes a `dest` that was
allocated *through the backend*. Passing a raw host ndarray works on warp's CPU device and fails on
CUDA. Commit `278e6e5` reported that failure as a GPU-only defect in the backend and `636b0c8`
retracted it: the API was right and the caller was wrong. Every scatter destination here comes from
`backend.zeros(...)`, and I4's control asserts it by type rather than by comment.

The harness itself (`compare_law_case`, the ULP arithmetic, the report) is host-side and is not
production physics. It reads each field to the host exactly once per case. It is a gate, not a step.

## 12. Comments and docstrings to discard

Nothing was read, so there is no source prose to discard. The list is kept as a positive statement of
what must **not** appear in the modules this entry authorises, and
`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` is the
mechanical half of it:

- no provider repository name, module path, kernel name, gate name or branch name;
- no absolute path into another project;
- no provider datum as a literal;
- no claim that a number was measured on hardware this lane did not run on. Every ULP figure in the
  package's docstrings is either measured here on warp's CPU device and says so, or is cited to
  `docs/design/GPU_STATE_2026-07-31.md` and `636b0c8` and says that.

What replaces them: the derivation in §4 and §10, restated in the modules' own docstrings in Aleph's
vocabulary, and the measured numbers from this repository's own run records.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Measured 2026-07-31 on **warp's CPU device**, warp 1.15.0, float32 compute, `ScatterMode.ORDERED`, 96 site pairs (69 engaged / 27 slack) scattered into 32 + 32 vertices. See the table below. `tests/runtime/test_law_parity.py`: **21 passed**; `tests/runtime` as a whole: **269 passed, 3 skipped**. Mutation study: **7 of 7 killed** (see §13a). The status stays `PROPOSED` because a device run has not happened and no human has reviewed it. |
| Reviewer | **Agent-proposed. Unratified.** Written by a subagent of session `46143f30`; no PI review, no `decided_by` field, and none may be added by an agent. |
| Rollback | Delete `aleph/runtime/law_kernels.py`, `aleph/runtime/law_cases.py`, `tests/runtime/test_law_parity.py` and the law-level block appended to `aleph/runtime/parity.py`. Nothing else imports them, so nothing else breaks — and that is the point: **every downstream CUDA task (G2–G6) then has no way to tell a fast wrong kernel from a fast right one.** |

### The measurement

| kernel | `energy_pn_um` | `forces_a_pn` | `forces_b_pn` | verdict |
|---|---|---|---|---|
| `TRUE` | **20.78** ULP | **100.44** ULP | **72.45** ULP | pass |
| `WRONG_ENERGY_HALF` | **8,388,566.43** ULP | 100.44 (bit-identical) | 72.45 (bit-identical) | **fail, energy only** |
| `WRONG_UNILATERAL_GATE` | 643,808.20 ULP | **3,175,210.17** ULP | **3,226,080.09** ULP | **fail, all three** |
| declared budget | 64 | 256 | 256 | — |
| round-off explains | 520.23 | 376.45 | 376.45 | — |

On the all-taut state (96 engaged / 0 slack) the true kernel measures 10.75 / 63.96 / 51.58 ULP and
`WRONG_UNILATERAL_GATE` measures **exactly the same three numbers** — bit-identical, and the gate
passes it. That is §9's blind-spot control and §14.4's limit, observed rather than argued.

Three things worth reading off the table:

1. `WRONG_ENERGY_HALF`'s force figures are *identical to the true kernel's, to the last bit*. A
   parity gate that reported one aggregated number, or graded only the forces, certifies it.
2. Both mutants clear the budgets by four to five orders of magnitude, so the 2.5–3× margin between
   the measured figure and the declared budget is nowhere near the detection threshold.
3. The true kernel's ~100 ULP is **an order of magnitude larger than the 10.00 ULP that CUDA atomic
   scatter costs**. For this law at cell-scale coordinates the scatter mode is not what limits
   accuracy; float32 position storage is. §14.7.

### §13a. Mutation study — 7 planted defects, 7 killed

Run with `PYTHONDONTWRITEBYTECODE=1` (a same-length mutant with a sub-second edit/revert leaves a
`.pyc` CPython considers valid, which has already contaminated one study in this repository).

| # | Planted defect | Killed by |
|---|---|---|
| M1 | `within_budget=True` unconditionally | 5 tests |
| M2 | `ulp_at_scale` returns `0.0` | 8 tests |
| M3 | the atomic floor applied to every field, not only scatter-assembled ones | 1 test |
| M4 | the "every output field must carry a budget" check disabled | 1 test |
| M5 | the driver launches `WRONG_ENERGY_HALF` whatever variant was asked for | 6 tests |
| M6 | the scatter destination allocated as a host ndarray instead of `backend.zeros` | 13 tests |
| M7 | the true kernel's `max(x, 0)` clamp deleted | 4 tests |

M6 is worth a sentence: the raw-ndarray destination fails **here**, on warp's CPU device, because
`WarpBackend.scatter_add` branches on `dest.dtype is wp.float32` and a numpy dtype is not that
object. So this particular regression would be caught on this laptop. That does not weaken §11 — it
narrows it: the retracted report in `278e6e5` concerned a destination that the CPU device *did*
accept, and the type assertion in I4's control covers both shapes of the mistake.

## 14. Honest limits

What this entry does **not** establish:

1. **No CUDA device was touched.** `~/.aleph_data/gpu_authorization.json` is expired and names a
   different host; both preflights refuse with exit 2 and an agent may never write that file. Every
   number here is from **warp's CPU device**, where a kernel launch is a serial loop. It is evidence
   that the kernel's *algebra* matches the reference, and it is `UNVERIFIED` as evidence about CUDA
   — in particular the `ATOMIC` rows are reproducible here by accident of the schedule, exactly as
   `measure_scatter_determinism`'s caveat already says.
2. **One law, and the smallest one.** The tether is a central, pairwise, unilateral spring with no
   connectivity. It exercises energy, forces, an exact-zero branch and a scatter. It does **not**
   exercise a mesh operator, a curvature stencil, or anything with a face loop — so nothing here is
   evidence that the harness's shape survives the membrane's Helfrich term (G2).
3. **`ATOMIC_SCATTER_ULP_FLOOR = 10.0` is one measurement on one device**, not a bound. It is not
   safe to assume a different GPU stays under it, and this entry provides no way to find out without
   a device run.
4. **The negative controls are as strong as the state they run on**, and one of them is provably
   blind on a taut configuration. That is asserted by a control, but it generalises: this method
   cannot tell you that a *case's inputs* visit every branch of the law. Choosing the state remains a
   human judgement and there is no gate on it.
5. **The declared budgets are engineering numbers, not physics.** They are measurements plus a
   margin. No stiffness, no force and no energy in this entry may be reported as a property of a
   cell; the evidence class is `STRUCTURAL`, and `BLOCKED` for quoting.
6. **The energy field is compared as a scalar total.** A kernel that got two sites wrong in opposite
   directions would cancel and pass. A per-site energy comparison would catch that and is not
   implemented; it is named here so the gap is on the record rather than discovered later.
7. **The float32-position conditioning finding in §10 is measured but not acted on.** Whether the
   CUDA port should carry positions in float64, or in coordinates local to each connector, is a
   design decision above this lane and is raised for the PI rather than taken.
