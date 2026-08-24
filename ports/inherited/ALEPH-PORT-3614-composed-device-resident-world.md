# ALEPH-PORT-3614 — the composed multi-owner device-resident world

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3614` |
| Lane | `46143f30` Lane G8c (CUDA track — composed residency) |
| Status | `PROPOSED` |
| Written | `2026-08-01` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Exists because | `ALEPH-PORT-3609` made **one** owner resident and wrote its own limit into §14.1: *"This is not `CellWorld` on device"* — the residency layer landed, but the twelve-configuration sweep it drove still had the cortex, the osmotic envelope and both connectors on the host, so its step was only *partly* resident. `aleph/scenarios/world.py` is the host-side answer to the composition problem and its owners are NumPy. **The gap this entry closes is the device-side counterpart of that composition**: several owners' state resident at once, coupled through connectors that never leave the device, presenting the same duck-typed surface `aleph.vertical.relax.relax_to_equilibrium` already consumes. |

> **This entry carries no `decided_by` field and no agent may add one.** It proposes work and records
> measurements. It settles none of the open questions in `HANDOFF.md` §C (the descent tolerance
> disagreeing between derivation and device) or §C-0 (`ALEPH-PORT-3608`'s reference-input guard),
> and it changes no ratified constant.

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered by it.

```python
# NEW — aleph/runtime/world_kernels.py
from aleph.runtime.world_kernels import (
    WorldKernelVariant,      # TRUE + two wrong on purpose
    load,                    # lazy: `import warp` happens here, never at module scope
    loaded_warp,
    numpy_moment_about,      # the frozen host reference for the moment kernel
    numpy_arm_length,        # the frozen host reference for the arm-length kernel
)

# NEW — aleph/runtime/resident_world.py
from aleph.runtime.resident_world import (
    ResidentWorldError,          # the typed refusal; never a silent degradation
    ResidentSurfaceOwner,        # Helfrich areal + bending, device-resident, N instances per world
    ResidentCompartmentField,    # volume-penalty compartment loading ONE owner, device-resident
    ResidentContactCoupling,     # the unilateral contact, device-resident, two owners
    ResidentTetherCoupling,      # the tensile tether, device-resident, two owners
    ResidentOwnerSlot,           # one owner and the phases it is scheduled in
    ResidentFieldSlot,           # one field and the owner it loads
    ResidentConnectorSlot,       # one coupling and its two declared endpoints
    ResidentWorld,               # the composition; satisfies relax_to_equilibrium's surface
    build_resident_world,        # the builder, and every refusal lives in it
    build_shell_pair_world,      # the demonstration world: 3 owners, 2 couplings, all resident
    recovered_tension_pn_per_um, # the Laplace oracle, read off a resident world
    main,                        # `python -m aleph.runtime.resident_world <job>`
)
```

**Nothing in `aleph/scenarios/**` or `aleph/vertical/**` is edited by this entry.** Both are
imported read-only. `aleph/runtime/residency.py` is **composed with, not edited**:
`ResidentSurfaceOwner` *holds* a `ResidentMembraneState` rather than reimplementing its launch
chain, which is the whole of `ALEPH-PORT-3501` §3.1's argument against a second driver applied here.

### 1a. Why the two utility kernels are not in `law_kernels.py`

Identical to `ALEPH-PORT-3609` §1a and for the identical reason. Every kernel in
`aleph/runtime/law_kernels.py` is a **constitutive law** with a parity row in a `law_cases.py`
builder, and `ALEPH-PORT-3608`'s controls C6–C8 iterate that registry and require each case to feed
its reference float64 positions, to be measurably sensitive to position precision, and to declare a
position mode.

`k_moment_about` computes `(x_i − c) × f_i`. `k_arm_length` computes `|x_i − c|`. Neither has any
constitutive content: they are the **bookkeeping the `AdjointPair` ledger needs**, and registering
them as law cases in order to satisfy a completeness guard would make that guard's own statement
less true. So they are graded here instead, and graded harder — at a budget of exactly `0.0` ULP
against their NumPy references, which no law case in the tree can be graded at.

`k_zero_f32` / `k_zero_f64` are **imported from `aleph/runtime/residency_kernels.py`**, not
redefined. A second zeroing kernel would be a second place for the stale-accumulator defect of
§4c to be got wrong.

---

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) — cited for the record only |
| Source path / symbol | **none named, none read.** No file under that tree was opened by this lane. |
| Read from | **neither `git show` nor the working tree.** |
| Working tree == commit? | not applicable — nothing was read |

Everything read is Aleph's own: `ALEPH-PORT-3601`/`-3603`/`-3604`/`-3605`/`-3608`/`-3609`/`-3610`/
`-3611`/`-3613`; `aleph/runtime/{backend,law_cases,law_kernels,ledger,participant,pipeline,`
`residency,residency_kernels,transaction,warp_kernels}.py`;
`aleph/vertical/{relax,membrane,assembly,connectors,pressure,cortex,nucleus_interior}.py`
(read, never edited); **`aleph/scenarios/world.py` (read only — another session's package)**; and the
thirteen `⚠ … LANDED` blocks of `docs/design/BUILD_PLAN-2026-07-31-engine-gaps.md`.

| | |
|---|---|
| Aleph modules composed | `aleph/runtime/residency.py` (`ResidentMembraneState`, `TransferCounts`, `CountingBackend`), `aleph/runtime/law_kernels.py` (contact, tether-step, volume, membrane areal/bending), `aleph/runtime/law_cases.py` (`membrane_state`, `_corner_index`, `SWEEP_GRID`), `aleph/runtime/{pipeline,transaction,participant,ledger}.py` |
| Design read and deliberately mirrored | `aleph/scenarios/world.py` — **read only.** Its four refusals and its `coverage()` idea are reproduced here for a device world, with attribution on each. |

---

## 3. Why source-derived porting beats clean-room

The same argument `ALEPH-PORT-3609` makes, one level up. A clean-room device world would be a second
answer to a question `aleph/scenarios/world.py` has already answered *and paid for*: each of its four
refusals is a configuration that assembled, stepped, and produced a plausible number before somebody
noticed. The turgor that was scheduled and contributed nothing (`max|F|` 86.5 pN against the verified
65.99). The connector energy read from a cache that is `0.0` before the first step, making every
coupling invisible to the descent (46.7 pN·µm low, and it looked well). The adhesions dropped twice
over by a silent shape check, so a spreading world relaxed to a plausible equilibrium with no
adhesion in it.

**A device world reproduces every one of those failure modes and adds two of its own** (§4c, §4d).
Deriving the refusals from the module that already earned them is cheaper than re-earning them.

---

## 4. Physical or mathematical law represented

None that is new. This entry adds **no constitutive law**: every force and every energy in a
composed resident world comes from a kernel that already has a parity row and a measured ULP budget
in `law_cases.py`. What is new is the **composition** — which owners hold state on the device at
once, in what order their contributions are assembled, and what is refused.

### 4a. What "composed residency" means, stated as an invariant rather than as an intention

`ALEPH-PORT-3609` §4b states the single-owner invariant. The composed one is:

> Between `push_positions` on any owner and `pull_positions` on any owner, **no owner's positions,
> force field or intermediate accumulator crosses the host boundary** — not for the owners' internal
> mechanics, not for the couplings between them, not for the field loads, and not for the
> accept/reject decision.

Everything that crosses is a **fixed-length** scalar: per coupling two resultant forces, two
resultant moments and two centroids; per field one energy. The count does not grow with the mesh,
which §10.1 measures at two resolutions and
`test_the_transfer_count_per_step_does_not_grow_with_the_mesh` asserts.
`test_a_composed_step_reads_back_scalars_only` audits the bulk half by spying on
`warp.array.numpy` and the `wp.array` constructor during a real step rather than by bookkeeping.

**One term does scale, it is declared, and it is not this entry's.** `WarpBackend`'s reductions are
two-stage and end in `part.numpy()` over `ceil(n / CHUNK)` float64 partials.
`ALEPH-PORT-3609`'s `TransferCounts` already counts that in a separate channel, and this entry
reports it beside the bulk figure and never folds it in.

### 4b. The four refusals inherited from `aleph/scenarios/world.py`, each with what it would otherwise produce

Reproduced for a device world with attribution, because each is a configuration that **assembles,
steps and reports a number**:

| Refusal | What it would otherwise produce |
|---|---|
| **An empty world** | Zero owners assembles, steps, conserves energy perfectly and means nothing. `aleph/state/full_census.py` records the same trap one level up; `ALEPH-PORT-3606` §13a M5 records it as a control passing vacuously on an empty set. |
| **Two owners under one name** | The second shadows the first and every per-owner control still passes, because each owner tests itself. On a device this is worse than on the host: the two would also **share a force buffer under one dictionary key**, so the second owner's internal forces would be silently added into the first's residual. |
| **A connector whose two declared endpoints are not both present** | Forces attributed to a body that is not in the world — and they look reasonable. |
| **A connector registered and never bound** | `aleph/vertical/assembly.py` uses exactly this as its structural counterexample (`bind_contact=False`), which means `Pipeline` will hold it. Here it is refused at composition time. |

Two more, which are this entry's own and are not in `world.py` because a host world cannot have them:

| Refusal | Why a device world needs it |
|---|---|
| **A coupling between two owners on different backends** | Two `WarpBackend` instances are two devices. A scatter from one into the other's buffer is undefined behaviour on CUDA and *silently works* on warp's CPU device, which is the machine the suite runs on — the G6 hazard shape (invisible on the device, live on the test machine) in a new place. |
| **A field whose target owner is not in this world** | `world.py` refuses this on the host (`OwnerSlot.load`, the turgor that reported 86.5 pN). On a device the same defect additionally has nowhere to scatter, so it would either raise deep inside a launch or write into a buffer of the wrong length. |

### 4c. The defect class composition adds to the one residency already created

`ALEPH-PORT-3609` §4c: a resident evaluation does not get zeroed accumulators by construction, so
every scatter destination must be explicitly re-zeroed and an omission grows linearly with step
count. **Composition adds a second, sharper version of the same thing.**

In a composed world the owners' force buffers are written by *more than one author*: the owner's own
internal chain writes them, then every coupling scatters into them, then every field adds into them.
So the buffers have to be zeroed **exactly once per evaluation, before the first author runs**, and
each author after the first must **accumulate** rather than assign.

An author that assigns instead of accumulating **erases every contribution written before it** and
the world still steps, still descends, and still reports a plausible energy — because the *energy*
is summed from separate per-term reductions and is unaffected. That is the `WRONG_FROZEN_NORMAL`
shape from `ALEPH-PORT-3604` (invisible in the energy and in one of two force fields) applied to
composition. It is mutants **M2** (a coupling) and **M3** (a field) in §9.1 — and **M3 survived the
whole first pass**, for a reason that is worth more than the mutant: every force-channel control
this lane had was reading the *membrane*, and the field loads the *shell*.

### 4d. Scheduled is not evaluated — coverage on a device world

`aleph/scenarios/world.py` names this as *"the one thing a scenario must check that this module
cannot"*: a world can assemble, step, and report an entirely plausible energy while a connector
witnessed nothing at all.

`ResidentWorld` therefore composes a **real `Pipeline`**, not a private step loop.
`ALEPH-PORT-3609`'s `ResidentDescent` deliberately did not — it drove two device buffers directly
and said so — and the price was that `relax_to_equilibrium` could not drive it and `coverage()` did
not exist. Here every resident owner is a `Participant` and every resident coupling is a `Connector`
in the pipeline's sense, so:

* a coupling that is registered and never bound is caught by `assert_schedule_complete`;
* a coupling that is bound and launches nothing is caught by the **witness delta**, because the
  witness is `backend.op_count + rng.draw_count` and a kernel chain that launches nothing moves
  neither;
* the connector books are closed by `BalancePredicate` against the ledger, from **two independently
  accumulated sides** — `force_b := -force_a` appears nowhere.

**This is the property a device world is most likely to lose, and it is why the transaction is not
reimplemented.** A private device step loop would have no witness, and "the coupling ran" would be a
property of a function body again.

### 4e. The oracle for composition itself: it must be a bit-identical refactor

Composition changes **which buffers the same kernels write into**, not the arithmetic. So the
composed world's per-owner assembled force field must agree with the same terms evaluated
separately through `law_cases`' non-resident drivers, **bitwise, per field, at every step of a
trajectory**, in exactly the sense `ALEPH-PORT-3609` §4d establishes for one owner.

One arithmetic difference is declared rather than hidden, and it is the same one G6 declared: a
composed total is summed **on the device in float32** where a host assembly sums in float64 having
read each field back. The two differ by one rounding, so **each term is graded separately and
bitwise** and the total is not a channel. `ALEPH-PORT-3604`'s finding — `WRONG_FROZEN_NORMAL` was
invisible in the energy and in one of the two force fields — is why.

---

## 5. Units, domains, singular cases, invariants

Units throughout: length µm, force pN, energy pN·µm, time s, pressure pN/µm², stiffness pN/µm.

**Domains and refusals.**

| Condition | Behaviour |
|---|---|
| Zero owners | `ResidentWorldError` (§4b) |
| Duplicate owner or connector name | `ResidentWorldError` (§4b) |
| Connector endpoint absent from the world | `ResidentWorldError` (§4b) |
| Connector registered and not bound | `ScheduleError` from `Pipeline.assert_schedule_complete`, raised by the builder |
| Owners on different backend objects | `ResidentWorldError` (§4b) |
| Field target absent | `ResidentWorldError` (§4b) |
| A `NumpyBackend` anywhere | `ResidencyError` from `ResidentMembraneState` — the chain is warp kernels and there is nothing for state to be resident *in*. The warp **CPU** device is fully supported and is what the suite runs on. |
| Mesh with zero vertices or zero faces | `ResidencyError` — every reduction over it is vacuous (`ALEPH-PORT-3606` §13a M5) |
| A contact site at non-positive normal standoff | `ConnectorGeometryError`, **raised, never clamped**, per `ALEPH-PORT-3604` §5 — see §5a |
| Position precision other than `POSITIONS_F64` | `ResidencyError`, inherited from `ResidentMembraneState`; `j4`'s PASS was measured under `POSITIONS_F64` |

### 5a. The contact's device-side refusal is a reduction here, and that is a finding

`law_cases._contact_candidate` implements `ALEPH-PORT-3604` §5's second refusal by writing a
per-site `violation` flag on the device and reading **the whole `int32` array of length V** back to
the host to count the flagged sites.

**That readback is a bulk crossing of exactly the kind this entry exists to remove**, and it happens
inside a force evaluation, every step. So the resident contact allocates its violation channel as
`float32` and reduces it with `Backend.max_abs`, which crosses `ceil(V / CHUNK)` float64 partials
instead of `V` int32 elements. The refusal is unchanged in meaning: any non-zero flag raises, and
the kernel still writes exactly zero at the offending site rather than clamping.

**The general form is worth carrying:** *a refusal implemented as a per-site flag array is a bulk
crossing; residency forces it through a reduction.* A design that reports diagnostics per element
cannot be resident, and the fix is not to drop the diagnostic but to reduce it.

The float64 host pre-check (`_contact_host_standoff_um`) is **not** performed per resident step,
because it reads the positions on the host. It is performed once, at composition, when the positions
come from the host anyway. This is a declared narrowing and is in §14.

**Invariants asserted by controls.**

| # | Invariant |
|---|---|
| I1 | A composed step crosses **zero bulk elements**; every crossing is a fixed-length scalar (§4a). |
| I2 | The composed per-owner force field is **bitwise equal**, per term, to the same terms from `law_cases`' non-resident drivers, over a trajectory. |
| I3 | Rollback is **bitwise** for every owner in the world simultaneously — not one owner at a time. |
| I4 | `coverage().ok` is `True` for a correctly built world and `False`, with the naming fault, for each of the four structural counterexamples. |
| I5 | The two moment/arm kernels are **bit-identical** to their NumPy references at a declared budget of `0.0` ULP. |
| I6 | Evaluating twice at one configuration is idempotent — the composed force field is bitwise unchanged. This is what asserts that the zeroing in §4c is complete. |
| I7 | A coupling's two sides are accumulated independently: `force_b` is never computed from `force_a`. Asserted by mutating one side's kernel and requiring the ledger's force-closure channel to notice. |

---

## 6. Source evidence class and known retractions

**Evidence class of everything this entry composes:** each underlying kernel carries its own parity
row and measured ULP budget; nothing here loosens one. The **composition** carries only the evidence
this entry measures.

Retractions and findings this lane is built on, carried forward so they are not re-learned:

1. `ALEPH-PORT-3601` §14.7 — the binding constraint is **float32 position storage**, not the scatter
   mode. `POSITIONS_F64` only.
2. `ALEPH-PORT-3603` §13 — a float32 **descent energy** gives 0/12 while the forces are fine.
3. `ALEPH-PORT-3604` — grade **energy and each force field separately**; a total is not a channel.
4. `ALEPH-PORT-3605` §13a M9 — a gate cannot see a defect in **the state it grades on**, and a
   resident buffer is state that survives between gradings.
5. `ALEPH-PORT-3606` §13a — nor in **the reference it grades against**; and a control can pass
   **vacuously on an empty set**.
6. `ALEPH-PORT-3609` §14 — `warp.array.numpy()` on warp's **CPU** device returns a **zero-copy
   view**. Every pull in this module copies explicitly, and the reason is in the docstring at the
   call site rather than here.
7. `ALEPH-PORT-3610` §14 — an **analytic oracle can prefer the wrong kernel**; check an oracle at
   more than one resolution. §7's Laplace oracle is therefore driven at **two mesh levels**.
8. `ALEPH-PORT-3611` §14 — **a test that stops at its first failure does not tell you how many
   things are broken.** Every control in this lane that loops collects failures and asserts once.
9. `ALEPH-PORT-3613` §14 — **a count derived from a hand-written list cannot see what is not on the
   list.** The composed world's per-owner and per-coupling enumerations are derived from the
   `Pipeline`'s own registrations, never from a list written beside them; an unknown slot kind
   **refuses** rather than being partly graded.
10. `ALEPH-PORT-3611` §14 — **a ULP budget is a property of the state's assembly, not of the
    kernel**, and a float64 tolerance is not a float32 target.

---

## 7. Independent oracle or derivation

**The Laplace sweep, driven by a fully composed resident world.**

At membrane equilibrium the virtual work of every force on the membrane vanishes for any virtual
displacement, in particular a uniform dilation `dx_i = x_i`. The tension term contributes `−2σA`
because area is homogeneous of degree two in the vertex positions; the discrete bending term
contributes **exactly zero** because it is exactly scale-invariant (`ALEPH-PORT-1101` I5, and
conditional on `c₀ = 0`, which this entry sets and asserts). What is left is

    2 σ A = W_transmitted = Σ_i f_i^transmitted · (x_i − x̄)

so `σ_recovered = W_transmitted / (2 A)`, measured about the membrane's own centroid so a translated
world gives the same number. This is `aleph/vertical/assembly.py::recovered_tension_pn_per_um`'s
derivation, **re-derived here for a world whose transmitted load is the sum of the resident
couplings** rather than for `VerticalWorld` specifically. No continuum limit is taken anywhere, so
it is exact at any mesh resolution — which is what makes it usable as an oracle at **two** levels,
per §6 item 7.

**The world it is driven on is NOT `build_vertical`, and the difference is stated rather than
buried.** Three of the vertical's terms are not on kernels at all: `ElasticCortexShell`'s per-face
quadratic areal spring, its edge spring, and `aleph.vertical.pressure.OsmoticEnvelope`'s surface
load. Porting them is three new constitutive laws, three parity gates and three ULP budgets — a
different lane's work, and one that would also have to satisfy `ALEPH-PORT-3608`'s C6–C8. So the
demonstration world substitutes, **declared term by term**:

| `build_vertical` | `build_shell_pair_world` | Why |
|---|---|---|
| `HelfrichMembrane` (σ, κ) | the same, resident | on kernels since `ALEPH-PORT-3603` |
| `ElasticCortexShell` (areal + edge springs) | a **second Helfrich shell** (σ_shell, κ_shell) | the cortex's two springs have no kernel |
| `OsmoticEnvelope` (`VOLUME_PENALTY`, K in pN·µm) | `nucleus_interior`'s volume penalty (K in pN/µm²) | on kernels since `ALEPH-PORT-3605`. **The two moduli are not interchangeable** and differ by a factor of `V_ref` — `nucleus_interior.volume_energy_and_forces`'s own docstring says a reader who copies a number across is wrong by ~1e2 in this vertical, and this entry does not copy one. |
| `CortexMembraneContact` | the same, resident | on kernels since `ALEPH-PORT-3604` |
| `ErmTether` | the same, resident, evaluated at a **single** configuration | on kernels since `ALEPH-PORT-3604`; the step law with `start == end` takes its midpoint branch, which is the static force |

**Therefore the recovered σ from this world is not directly comparable with `j4`'s 1.5748e-05, and
§10.3 says so at the point of reporting.** What *is* comparable is the pair measured here: the same
composed world driven host-side through frozen `aleph.vertical` NumPy, and driven device-side
through the resident composition. That comparison is an identity, not a band.

---

## 8. Positive control

All in `tests/runtime/test_resident_world.py`. **30 controls, 30 passing.**

| Control | Asserts |
|---|---|
| `test_a_composed_world_assembles_and_reports_complete_coverage` | Five registrations compose, step through a real `Transaction`, and `coverage().ok` is `True` with every one of them covered and no zero witness delta. |
| `test_a_composed_trajectory_is_bit_identical_to_the_separate_law_drivers` | I2 — per owner, per coupling, per field, **bitwise**, over a six-step trajectory, against `law_cases`' own non-resident drivers. Failures collected, asserted once. |
| `test_a_composed_step_reads_back_scalars_only` | I1 — audited by a spy on `warp.array.numpy` and the `wp.array` constructor during a real step. |
| `test_the_transfer_count_per_step_does_not_grow_with_the_mesh` | I1's real content: the count at level 1 equals the count at level 2. |
| `test_rollback_restores_every_owner_bitwise_together` | I3. |
| `test_evaluating_twice_at_one_configuration_is_idempotent` | I6 — the §4c zeroing is complete in all three kinds of author. |
| `test_the_moment_kernel_is_bit_identical_to_its_numpy_reference` | I5, at `0.0` ULP, for both bookkeeping kernels. |
| `test_the_verified_relaxation_driver_can_drive_a_resident_world` | `aleph.vertical.relax.relax_to_equilibrium` — **unmodified** — drives the composed device world to a lower energy with coverage enforced and the books balanced every step, uploading zero elements. |
| `test_the_laplace_identity_holds_at_two_resolutions` | §7, at two mesh resolutions and two tensions, per §6 item 7. |
| `test_both_entry_points_reach_one_evaluation_implementation` | There is one assembly, not two. |
| `test_a_three_owner_world_composes_and_covers` | The layer is `n`-owner: three owners, three couplings, one field. |
| `test_the_world_energy_is_the_sum_of_every_term` | Every term reaches the descent. |
| `test_the_world_residual_is_the_maximum_over_every_owner` | Convergence is read off the world, driven on a world whose **first** owner is not the loudest. |
| `test_each_owner_total_is_the_sum_of_every_author` | I7's assembly half, for **every** owner rather than the one the oracle reads. |
| `test_the_two_sides_of_a_coupling_are_reduced_independently` | I7 — the two sides are near-opposite but **not bit-exact** opposites. |
| `test_the_float32_unit_roundoff_matches_the_verified_module` | The duplicated constant cannot drift from `relax.py`'s. |
| `test_a_slack_tether_contributes_exactly_zero` | Exactly zero, with the contact non-zero at the same configuration as the vacuity half. |
| `test_a_pulled_buffer_is_a_copy_and_not_a_view_of_the_device` | §6 item 6 — the vacuity half of every bitwise control above. |

## 9. Deliberately failing negative control

| Control (must fail) | Asserts |
|---|---|
| `test_an_empty_world_is_refused` | §4b, with the reason in the message. |
| `test_two_owners_under_one_name_are_refused` | §4b. |
| `test_a_coupling_to_an_absent_owner_is_refused` | §4b. |
| `test_a_registered_and_unbound_coupling_is_refused` | §4b — the counterexample `assembly.py` uses, refused at build time. |
| `test_a_coupling_across_two_backends_is_refused` | §4b, this entry's own. |
| `test_a_field_whose_target_is_absent_is_refused` | §4b, this entry's own. |
| `test_a_coupling_that_launches_nothing_fails_the_witness` | §4d — a stub returning a well-formed `AdjointPair` with closing books, caught by the witness delta and by nothing it returned. |
| `test_an_unknown_registration_kind_is_refused_rather_than_partly_graded` | §6 item 9, in the mechanism. |
| `test_a_wrong_moment_arm_is_invisible_to_the_ledger_residual_and_visible_in_its_scale` | **The finding of §14.3**: the ledger's moment channel is structurally unable to grade a moment reference point. Both halves asserted. |
| `test_an_assigning_coupling_erases_the_owner_internal_forces` | §4c — the counterexample that is **bit-identical in the energy**. |
| `test_a_violation_flag_does_not_survive_into_the_next_evaluation` | §5a — a resident flag that is not re-zeroed makes a refusal that fires once fire for ever. |
| `test_the_device_ledger_tolerance_is_not_inherited_from_float64` | §6 item 10, driven on a **jittered** world because the pristine one closes to exactly `0.0`. |

### 9.1 Mutation study

Planted one at a time, each reverted before the next, `PYTHONDONTWRITEBYTECODE=1`, second pass from
a **verified-empty** bytecode tree (`find aleph tests -name '*.pyc' | wc -l` → 0).

> **24 planted, 24 killed. Two survived the first pass and both were real gaps in the controls
> rather than equivalent mutants; the controls that close them are named below.**

| # | Mutant | Killed by — and **which channel** |
|---|---|---|
| M1 | both coupling sides land on `owner_a` | coverage + **force** (the Laplace identity, the relaxation) |
| M2 | a coupling assigns instead of accumulating | **force**, five controls; the **energy is bit-identical** |
| M3 | the field assigns instead of accumulating | **force**, `test_each_owner_total_is_the_sum_of_every_author` **only** — see below |
| M4 | the contact does not re-zero `normal_sum` | **bitwise parity against the non-resident driver** |
| M5 | the contact does not re-zero the violation flag | **the refusal channel** — it fires for ever |
| M6 | the compartment does not re-zero `triple_total` | energy + force + idempotence (5 controls) |
| M7 | the compartment does not re-zero its force buffer | **force**, and idempotence |
| M8 | the world residual is taken from one owner | `test_the_world_residual_is_the_maximum_over_every_owner` **only** — see below |
| M9 | the world energy drops every coupling and field | **energy** |
| M10 | the empty-world refusal removed | the refusal control |
| M11 | the duplicate-name refusal removed | the refusal control |
| M12 | the absent-endpoint refusal removed | the refusal control |
| M13 | the unbound-connector check removed | the refusal control |
| M14 | the one-device refusal removed from the builder | the refusal control |
| M15 | the absent-field-target refusal removed | the refusal control |
| M16 | every backend reports one device | the refusal control |
| M17 | the kind enumeration grades everything as an owner | the enumeration control **and** the three-owner control |
| M18 | the derived ledger band 8,000× too tight | the tolerance control |
| M19 | the coupling reports `force_b := -force_a` | the independence control **and** the tolerance control |
| M20 | the Laplace oracle divides by `A` instead of `2A` | **the oracle** |
| M21 | the moment kernel flips one component | **parity at 0.0 ULP** — and by nothing else |
| M22 | the arm kernel drops the `z` term | **parity at 0.0 ULP** |
| M23 | the flag reduction takes a minimum | the refusal control (it never fires) |
| M24 | the integer zeroing writes one | the refusal control, and five others through the raise |

**The two first-pass survivors, because they are the useful part of the study.**

* **M3 — a *field* that assigns instead of accumulating survived every control this lane had.** The
  field loads the **shell**, and every force-channel control was reading the **membrane**: the
  Laplace oracle is about the membrane, the bitwise parity control compares each author's own buffer
  rather than the owner's total, and the energy is unaffected. *An oracle that reads one owner is
  blind to a defect in another owner's assembly*, and in a multi-owner world that is a new shape of
  the same error the series has been collecting all night.
  `test_each_owner_total_is_the_sum_of_every_author` reconstructs **every** owner's total from its
  authors and is what closes it.
* **M8 — a residual taken from `owners[0]` survived, because on the demonstration world the
  first-registered owner is also the loudest one.** The control was correct, non-vacuous, and never
  exercised at the boundary it defends — **V2's finding exactly**. It now runs on a world whose
  first owner is the quiet one.

---

## 10. Numerical and precision envelope

Filled by measurement. `PositionPrecision.POSITIONS_F64` throughout, `ScatterMode.ORDERED`
throughout (`ALEPH-PORT-3601`: `ORDERED` is bit-identical across launches on the A5000 while
`ATOMIC` measured 10.00 ULP).

### 10.1 What crosses the host boundary per composed step — MEASURED

`python -m aleph.runtime.resident_world transfers --cpu --repeats 16`, warp CPU device, and
independently audited by the `warp.array.numpy` / `wp.array` spy in
`test_a_composed_step_reads_back_scalars_only`.

One composed step is one full `Transaction.step` — snapshot, `ACCUMULATE_INTERNAL` ×2,
`ACCUMULATE_COUPLING` ×2 with their `AdjointPair` books, `ACCUMULATE_FIELD`, `SOLVE` ×2,
`EVALUATE_GATES` ×2, commit or rollback — plus the two whole-world energy evaluations the
acceptance rule makes.

| | level 1 (84 vertices) | level 2 (324 vertices) |
|---|---|---|
| **owner positions / forces / accumulators crossing** | **0** | **0** |
| host→device elements | **0** | **0** |
| device→host **scalar** elements | **39** | **39** |
| device→host scalar calls | 15 | 15 |
| reduction partials (`WarpBackend`'s own two-stage `ceil(n/CHUNK)`) | 30 | 38 |

**The 39 is fixed and is the claim.** It itemises exactly: two couplings × (2 centroids + 2 resultant
forces + 2 resultant moments) × 3 = 36, plus the compartment's one-element energy buffer read once
per whole-world assembly, of which a step makes three.

**The partial count is the one term that grows**, as `ceil(n / 256)` per reduction, and it is
`WarpBackend`'s decomposition rather than this entry's — `ALEPH-PORT-3609` counts it in its own
channel for the same reason. 30 → 38 across a 3.9× mesh.

Against `ALEPH-PORT-3609` §10.2, and the comparison is **not** like for like, which is why both rows
carry what they are:

| | calls | elements |
|---|---|---|
| `-3609`, one non-resident `membrane_kernel_forces` evaluation, one owner, no coupling | 6 | **1,947** (≥3,867 counting the connectivity `CountingBackend` cannot see) |
| `-3609`, one resident step, one owner, **no transaction, no ledger** | 3 | **4** |
| **here**, one resident step, **2 owners + 2 couplings + 1 field, through a real `Transaction`** | 45 | **39 scalars + 30 partials, and 0 bulk** |

A single owner has no connector books; a world of five registrations has to produce them, and the
`AdjointPair` the ledger reads is 12 numbers per coupling by construction. The number that is
comparable between the two rows is the bulk one, and it is **zero** in both.

### 10.2 The composed world's own parity — MEASURED

`test_a_composed_trajectory_is_bit_identical_to_the_separate_law_drivers`, over six steps of a
trajectory, five force fields and four energies per step, against `law_cases._contact_candidate`,
`._tether_step_candidate`, `._volume_candidate` and `.membrane_kernel_forces` — each of which
uploads, allocates and reads back independently.

> **0.0 ULP on every field and every energy, at every step. It is an identity, not a budget.**

The one place a composed total is *not* bitwise is stated rather than hidden and is the same one
`ALEPH-PORT-3609` declares: an owner's **total** is summed on the device in float32, where a host
reconstruction sums the same contributions in float64. `test_each_owner_total_is_the_sum_of_every_author`
grades that at `64·u32` relative and measures well inside it; the per-author fields above are what
carry the exact claim.

### 10.3 The Laplace sweep on the composed resident world

To be measured by `python -m aleph.runtime.resident_world sweep`, reported beside — and explicitly
**not equated with** — the three published figures: `j4` 12/12 at worst 1.5748e-05 (forces on
kernels, descent on host); `ALEPH-PORT-3609` 2/12 at `2·u32` and 8/12 at `1·u32` (descent on
device). §7 states why the world is different. **No search is made for a slack that lands in any
band.**

---

## 11. Production-backend residency and transfer

Every buffer in a composed world is allocated through `Backend.zeros` / `Backend.array`, never as a
raw host `ndarray` handed to a launch — `636b0c8` is the retraction that made this a rule, and it is
the rule most easily broken by a composition, because a coupling allocates *between* two owners and
has no obvious home for its scratch.

The connectivity of each owner's mesh is uploaded **once**, at composition, with `wp.array` rather
than `backend.array` (`law_cases._device_triangles` says why: `WarpBackend` accepts float32 and
float64 only). `CountingBackend` cannot see it, so every non-resident figure this entry reports is a
**lower bound**, stated as such at the point of reporting.

### 10.4 The device run that was NOT made, and why

`gbook`'s A5000 was reachable and the grant covers it until 2026-08-04. It was **not used**, and the
reason is a rule rather than a preference: at 09:1x KST `nvidia-smi` on `GBook` reported one compute
application, pid `1553892`, 52 minutes in —

```
/home/sungwook/miniconda3/envs/aleph/bin/python .../scripts/cuda_jobs.py
    --job k4_descent_window_wide --out ~/.aleph_runner/artefacts/k4_descent_window_wide.json
```

which is **`ALEPH-PORT-3613`'s `k4`, still in flight**, with `systemctl --user is-active
aleph-cuda-runner` returning `active`. `CLAUDE.md` §3: *never start on a GPU with someone else's job
running.* Nothing was launched, nothing was queued, and no authorization record was written or read
by this lane. Every figure in §10 is therefore **warp CPU device**, and every device claim in this
entry is `UNVERIFIED` on CUDA.

That is not a small caveat here. `ALEPH-PORT-3609` §14.6 measured that the *descent window itself*
moves between the warp CPU device and CUDA, because CUDA contracts FMA and implements `sqrt`/`rsqrt`
differently and the float32 energy's noise floor is a different number. **A composed world has more
float32 terms than a single owner, not fewer**, so there is no reason to expect §10.3's numbers to
transfer, and this entry does not claim they will.

## 12. Comments and docstrings to discard

Nothing is copied from the reference project by this entry, so there is nothing to discard. Two
statements in Aleph's own tree are **read and deliberately not repeated**: `world.py`'s claim that
its refusals are about host NumPy owners (they are about composition, and transfer), and
`residency.py`'s §14.1 statement that composing several owners is out of scope — which this entry
supersedes rather than contradicts.

## 13. Acceptance

| | |
|---|---|
| `tests/runtime/test_resident_world.py` | **30 passed** |
| `tests/runtime` | **608 passed, 3 skipped, 1 failed** — the failure is `test_reference_input_guard.py::test_every_registered_case_reference_is_sensitive_to_position_precision`, which is `ALEPH-PORT-3608`'s guard, **red before this lane started**, foreign, and named in `HANDOFF.md` §C-0. Not fixed and not routed around. |
| `tests/ports` | 19 passed, 2 failed — both are `INDEX.md` staleness, already red for eight entries before this one and **not regenerated**, on instruction. `test_named_controls_resolve_to_real_tests` **passes**: every control this entry names is defined. |
| Mutation study | **24 planted, 24 killed**, second pass from a verified-empty bytecode tree |
| Bulk crossings per composed step | **0**, at both resolutions |
| Parity against the non-resident drivers | **0.0 ULP**, per field, per energy, over a trajectory |
| `aleph/scenarios/**`, `aleph/vertical/**`, `aleph/viz/**`, `scripts/**`, `aleph/state/**` | **byte-for-byte untouched** |

## 14. Honest limits

Every one of these is a limit of what was measured, not a thing to fix later without saying so.

### 14.1 The demonstration world is not the vertical, and no number here is a vertical number

`build_shell_pair_world` substitutes a second Helfrich shell for `ElasticCortexShell` and
`nucleus_interior`'s volume penalty for `OsmoticEnvelope`, because the cortex's areal spring, its
edge spring and the osmotic surface load **have no kernels at all**. Porting them is three
constitutive laws, three parity gates and three measured ULP budgets, plus `ALEPH-PORT-3608`'s C6–C8
for each. It is the single largest piece of remaining work between here and "the vertical runs on
the GPU", and this entry does not start it.

Consequence: §10.3's recovered tensions are **not comparable** with `j4`'s 1.5748e-05 or with
`ALEPH-PORT-3609`'s 2/12 and 8/12, and they are reported beside those numbers with that said at the
point of reporting rather than in a footnote.

### 14.2 The pipeline's coverage witness cannot see a raw `wp.launch` — and that is not this lane's to change

`StepContext.witness_count` is `backend.op_count + rng.draw_count`. A `wp.launch` moves neither. So a
device-resident participant whose phase work is *only* kernel launches is, to the coverage gate,
indistinguishable from a placeholder that launches nothing.

It fails **closed**, which is the safe direction — no stub can pass — but it is a false positive:
real work reported as none. Every phase handler here therefore reports a quantity it measured
*through the backend*, which is what a receipt is for and is what `HelfrichMembrane.accumulate`
already does. **That is a workaround at the call site, not a repair**, and the repair belongs to
whoever owns `aleph/runtime/pipeline.py`: a witness that counted kernel launches as well as backend
ops would make the gate true for device participants without weakening it for host ones. Raised,
not taken (`CLAUDE.md` §2 rule 5).

### 14.3 The ledger's moment channel is structurally unable to grade a moment reference point

Measured, and it was expected to go the other way. With force closure holding,

    sum_a (x - c) x f + sum_b (x - c) x f = [sum_a x x f + sum_b x x f] - c x (F_a + F_b)

and `F_a + F_b = 0` removes `c` from the residual **exactly**, at any displacement. The `WRONG_ARM_FROM_ORIGIN`
variant therefore changes the moment **scale** — inflating it on a displaced world, which makes the
band *wider* and hides real failures — while leaving the residual untouched.

**The only channel that kills it is the kernel's own parity gate at 0.0 ULP.** A world that shipped
`k_moment_about` without a bitwise reference would have had no grader for it at all. This is
`ALEPH-PORT-3606` §13a M8's shape — a gate cannot see a defect in the reference it grades against —
arriving at a *ledger* channel rather than at a parity one.

### 14.4 The device ledger band has ~160× of headroom over the measured residual

`DEVICE_LEDGER_ROUNDOFFS = 8` gives a relative band of `4.77e-07`. The worst closure ratio measured
on a jittered world is `3.01e-09`, i.e. **0.05 `u32`**. The derivation is a bound on the deepest
assembly path (seven roundings) and the measurement is far inside it, which means the band would
**not** catch a small systematic closure defect — anything under ~160× the round-off passes. What
catches those is the bitwise parity of §10.2, not the ledger. The band is a floor on noise, not a
test of correctness, and it should not be read as one.

**And on the *pristine* world the question does not arise at all**: two concentric icospheres give
every per-site force an exactly-opposite partner and the float64 resultant is **exactly `0.0`**. A
tolerance control written on the demonstration world would have reported a band that had never been
exercised — so `test_the_device_ledger_tolerance_is_not_inherited_from_float64` is driven on a
jittered world, and says so.

### 14.5 The contact's float64 host pre-check is not performed per resident step

`law_cases._contact_candidate` refuses twice: once on the host in float64 before the launch, once on
the device's own arithmetic after it. The host half reads the positions, which is a bulk crossing, so
a resident step performs it **only at composition**. The device half — §5a's reduced flag — runs
every step and is unchanged in meaning.

The narrowing is real: a configuration that becomes non-positive in float64 *during* a run is caught
by the device predicate rather than by the host one, and those two predicates are not the same
predicate (`ALEPH-PORT-3604` §5 says so explicitly, and whether they differ depends on
`PositionPrecision`). Declared here rather than discovered later.

### 14.6 One term in the transfer count still scales, and it is inherited

`WarpBackend`'s two-stage reductions read back `ceil(n / CHUNK)` float64 partials. 30 elements per
step at level 1, 38 at level 2. That is `ALEPH-PORT-3609`'s declared term, counted in its own channel
and never folded into the bulk figure.

### 14.7 Everything is `UNVERIFIED` on CUDA

§10.4. `gbook`'s A5000 was carrying `ALEPH-PORT-3613`'s `k4` when this lane finished, so no device
run was made. Given that `-3609` measured the descent window itself moving between the two devices,
and a composed world has **more** float32 terms than a single owner, §10.3 in particular should be
expected to move.

### 14.8 The two open PI decisions are untouched

`HANDOFF.md` §C (the descent tolerance disagreeing between derivation and device) and §C-0
(`-3608`'s reference-input guard, unclearable by at least three laws linear in the rounded input) are
neither settled nor worked around here. This entry runs at the ratified `2·u32` and reports what
happened.
