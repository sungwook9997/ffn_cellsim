# STATE LOG — the running record `STATE.md` (f) used to carry

Moved out of `STATE.md` on 2026-07-28. `STATE.md` (f) had become a session-by-session narrative and was
the second-largest thing in a file whose whole job is to be short. It keeps the one line a reader needs
— the date and what last changed — and appends the detail here, newest first.

**⚠ Everything here is status and most of it is superseded by the entry above it.** Nothing in this file
is quotable on its own; `STATE.md` (b) and (c) decide what may be cited.

---

## 2026-07-29 (overnight) 이전 (f) 전문


LAST VERIFIED: 2026-07-28 (night, 4) — **the stationarity detector is wired into the driver it was written for, and it has tests for the first time** (`4095882c`). `ac_gate_b_cortex_motor_native` now judges its own settling per observable and records the verdict, tau_int, discard point and contract; `--until-stationary` replaces
the step target with the series; `--gamma-every` buys the physical time that judgement needs — the per-step gamma readback is **53% of this driver's wall-clock** (3.85 s against a profiled 3.39 s subdiv-7 step). **The verdicts are RECORDED, not a gate**: PI held the amendment AND held `min_windows` for discussion today, so
nothing here rejects a run. Two findings bear on those held decisions. **The driving time was wrong three ways at once, and it is the denominator the contract is written in** — the driver stamped `1/0.4`; `NMII_KOFF0` is **0.35**; and the sweep ran `--catch-slip`, whose zero-load rate is the SUM `k_catch0+k_slip0` = 0.70/s ->
**1.43 s**. Re-judged against it, **all SIX blebbistatin points are `TOO_SHORT`, not five** ((c) 16): required 7.14 s, ran 0.60 s. **And the detector is optimistic exactly where it is used** — its "window never closed" branch is DEAD CODE (mean-subtraction forces `sum rho_k = -1/2` exactly, verified to 1e-9), so it cannot say a
series is too short for its own tau, while the measured bias at 60 samples of tau_int = 49.5 steps is **8% of the truth** (tau too small -> `n_eff` too large -> `sem` too tight). A contract in DRIVING lifetimes does not cover that, because tau's accuracy is set by the OBSERVABLE's own correlation time — **which nothing has yet
measured**, and a native run to measure it is on the device now. — **The framing error is still the largest item on this page**: the runtime was solving for a STATIC equilibrium that does not exist for an active cell. Measured control (§10 of `GATE_A_RESTING_CONVERGENCE_2026-07-23.md`): myosin OFF converges seed-robustly; 4,420
isometric heads inject a residual proportional to `f_head` that **no solver removes**. The force field is non-conservative ((b)), so there is no potential to minimise and the symmetric tangent assumed one. `CLAUDE.md` carries the PROPOSED gate amendment plus a new stage 0.5 (timescale separation); both are UNRATIFIED. — Budget
facts from the first native kernel profile (night, 3), which nothing else records: one inner iteration is **29.66 ms** (not the 192 ms wall/n_inner implied), filament physics is under 12% and the cost is fluid + moving boundary + membrane; **4x the filaments costs 1.5% more time**, so a slice does not buy its population ratio
and the sweep budget cannot be solved by shrinking the cell; membrane **subdiv 7 == 8 to 5 decimals** while subdiv 6 differs 4.17% (31 sigma), so 8->7 is free and worth 1.37x. `gamma_source` is the crossbridge STRESS TERM, not myosin's contribution — myosin acts by loading the network. — Standing rule from the PI pushback that
created (c) 14: **develop on a slice, conclude at native**. Earlier-session detail is not repeated here; it is in the commits, in (b)'s rows and in `STATE_NONQUOTABLE.md`. What remains live from it: **`.git` stays 7.1 GB** (only a history rewrite shrinks it, and that breaks the commit hashes this file rests on); **`ffn_gpu.py
run` refuses to launch when the remote tree does not hold the code the run would stamp** (`1192cc20`, AST closure vs ssh digests — it caught 2 of 79 files stale on first use, and a wrong build stamp cannot be corrected afterwards); **`CellState` is the 9th KB DB** (K=2 x M=2), so stage 5's denominator is representable, with
`ValidationGate` VOID and **`RunResult` population** still open — the latter is the absence behind (c) 14; and T10's exterior medium ran native at rung `CUDA_UNIT` but **its gate FAILED 3/10 and the FAIL STANDS**, `GC-0002` is filed **PENDING** with no re-run before signature, and **BLOCKS_CRAWL is NOT released** (connector
undispatched, `HALF_SPACE_BLAKE` raises, `mu_medium` a PI-GAP). — Still true, and the thing every roadmap item runs into: **no lane has a converging coupled solve** (`nmii_sf_motor` converges ONE of its two bodies; `nmii` plateaus at ~8e-02 pN) and **`ac/engine` has no coupled solver at all** — `SOLVE_COUPLED` is an empty phase,
so the 36/36 facade-owned connectors transfer force once and never reach balance ((c) 4). `overlap_free` is default-False. Branch `codex/ff-ac-codex`. Re-date at every milestone.


---

## (f) as of 2026-07-29 (overnight)

LAST VERIFIED: 2026-07-29 (overnight, autonomous) — **`ac/engine` can compute the cortex in a production run for the first time.** PI approved **`COMPARTMENT_VALIDATION_TRACKS` §6 D2** (`_accumulate_all(omit=)`; additive, bit-identical default), which four validation tracks were blocked behind; `make_inner_solve` forwards it, and behind `--engine-cortex` the
`CortexStateOwner` supplies the cortex channels while the incumbent omits them. **Honest scope: the arrays are still shared** — a private force array needs the inner solve to relax TWO arrays together, which is `SOLVE_COUPLED`, still an empty phase ((c) 4). So this closes "ac/engine fires no kernel" and gives the cortex an
owning object; it does NOT end co-location. The T2 device-side double-count control — which the gate script had to name as un-runnable — now runs: mask OFF must fail by an EXACT factor of two. Native A/B and the T2 gate are queued. — **The blebbistatin sweep is running to STATIONARITY** rather than to a step count; two points
have settled and `k_xb=1000` reproduces its own separate calibration run to four decimals. — **Cleanup, measured rather than asserted**: a native run's import closure went 129 files / 54,367 lines -> 119 / 50,635 and its actual import set 115 / 45,045 -> 101 / 37,201. The cuts were all ACCIDENTAL coupling: a seven-hop chain that
reached the parked `dcm/` engine to get two constants and a KD-tree query; five rig classes imported only to be dict KEYS in a lookup table; two eager subpackages. What remains is not core — of 101 imported modules **20 fire a kernel (22%)**, 16 define kernels and fire none, and 69 define none; the largest single item is a
2,093-line second solver imported unconditionally for a run that uses the explicit path. — **Two guards caught me, both correctly**: the I0-A static contract rejected a CPU-device literal in my own test (twice — the second time in a comment explaining why it is forbidden), and the KB coverage gate blocked a commit whose
REPORT.md headline no manifest declared. — **Two defects I caused**: five queued sweep points were silently discarded because I committed source without re-syncing the run host and the drainer recorded each refusal as "done" (fixed: exit 3/4 now requeue with a backoff and a bounded retry, since neither ever RAN), and I committed
once on a test run I had not read. — And the whole `aleph/tests/ac` suite had been collecting ZERO tests: a second `test_contracts.py` under an untracked `virtual_cell/` collided on module basename, and a basename collision INTERRUPTS collection rather than failing one file. — `CLAUDE.md` gained the prohibitions from the
PROPOSED virtual-cell master plan that need no signature (a force-accepted trajectory is not an accepted physical one; a surrogate is not reference physics; a gate scores a pre-declared projection; no model is scored against its own machinery's answers; nothing is "quantum" without quantum dynamics). Its IDENTITY claim is NOT
adopted — that plan is a proposal and says so. — Notion Dev Logs: four sections that each declared themselves stale were moved to a dedicated archive with per-section retirement markers, three removed from the live page. — arity detector is wired into the cortex-motor driver and has tests; the driving correlation time it is
judged against was wrong three ways, and re-judging with the corrected value retires **all six** blebbistatin points as `TOO_SHORT` ((c) 16). First settled native measurement in this lane landed the same night: four observables STATIONARY over 10.3 bound-head lifetimes. The boot context was cut from ~28,300 to ~14,000 tokens —
`CLAUDE.md` 524 -> 227 lines with the stage gates moved to `ROADMAP.md`, `MEMORY.md` stripped of the findings and magnitudes it was injecting ahead of their own retractions, and this section's running narrative moved to `aleph/docs/STATE_LOG.md`. **Session-by-session detail:
[`aleph/docs/STATE_LOG.md`](aleph/docs/STATE_LOG.md).** Re-date at every milestone; a stale `STATE.md` is a defect, not a formality.

## 2026-07-29 (overnight, later) — two "pending cleanup" items closed by measurement, not by editing

**The last import-cleanup item was a non-item, and the measurement says so.** `implicit_mechanics`
(2,093 lines, 55 kernels, none fired) and `contact_schwarz` (932 lines, 11 kernels, none fired) were left
as eager imports of the feature-frozen `ac/cell/driver.py`, deferred as needing a D2-style PI decision.
Measured before proposing one:

* **Cost:** 40.5 ms and 42.6 ms of a 605 ms driver import — 13% of import, and 0.003% of a 40-minute
  native run. Not an optimisation.
* **Closure:** removing both from the driver would delete **zero** files from the run's import closure
  (120 files). Both are imported independently by modules already in it — `erm_schwarz` and
  `erm_gauss_seidel` pull `contact_schwarz`; `medium_exterior` and `sf_implicit` pull
  `implicit_mechanics`. A lazy import in the driver would change the file list not at all.

So there is nothing to gain and a frozen file to touch. **Closed as measured, not deferred** — and the
general statement follows: the import cleanup is finished. The remaining 81 of 101 modules that fire no
kernel are imported for constants, dataclasses and helpers, which is ordinary Python dependency, not the
accidental coupling the earlier cuts removed. The 22% figure is not a backlog.

**The sweep landed three points and the ratio is still refused — correctly, on provenance.** `k_xb` =
1000 / 100 / 10 STATIONARY at 4.3256 / 3.4807 / 2.8839 pN/µm, bound fraction flat at 98.2–98.3%. Two
things the label hid are now printed and tested: `k_xb=10` passed by 2% (0.98 σ against 1.0, n_eff 21),
and its `γ_source` did not settle at all. Neither blocks R — it is a `γ_total` ratio — but a refusal the
reader never sees is the defect the analyser exists to prevent, one level down.

**Three sweep points were dropped by arithmetic and the arithmetic is recorded.** The transient grew
0.40 → 0.40 → 3.90 s across the measured decades, i.e. ≈ 1/`k_xb` once it leaves the detector's floor.
`k_xb` = 3, 1, 0.1 therefore need ≈ 2,300 / 4,600 / 46,000 steps against a 2,500-step cap and a 75-minute
wall: none can return a verdict at any affordable setting. `k_xb=3` was stopped three minutes in rather
than allowed to burn 75 minutes to a foregone `NOT_SETTLED`. The sweep spans 2 decades, not 4, for a
measured cost — and R will be an upper-bounded ratio at a stated softest `k_xb`, not a `k_xb`→0 limit.

**A boot-budget breach I caused and did not legislate away.** Adding the sweep result put the boot total
588 B over. Raising the budget to fit my own edit is the gate-loosening this repo forbids, so text came
out: the guards/defects paragraph (this file already owns it in full) and a closing line that restated
`CLAUDE.md`'s own carve-out verbatim. Worth knowing for next time: each file is individually UNDER budget
(+806 / +1,245 / +634 of slack) and the **total** binds at 46,000 against a 49,000 sum of parts. That is
deliberate — the parts were sized generously so the total is what forces a choice.


## 2026-07-29 (overnight, later) — the two boot questions, answered with measurements instead of adjectives

**"Is the boot budget resolved?" — NO, and it is at the wall: 45,989 / 46,000 B, ELEVEN bytes of slack.**
Per file there is room (CLAUDE.md +806, MEMORY.md +1,527, STATE.md +678 = 3,011 B), and the TOTAL still
binds because 46,000 was deliberately set below the 49,000 sum of parts. So the constraint is not "which
file is too big" but "which file may spend the next result". Six breaches tonight, six paid down in text,
every cut a real duplication — and the next result will breach it again. Raising the number to fit a
session's own edit is the gate-loosening this repo forbids, so it stays a PI decision.

**"Are the kernel-less modules that still get imported resolved?" — the question conflates three things,
and only one of them was ever a defect.** Measured on the current closure (120 files / 51,467 lines):

* **43 files (36%) call `wp.launch`.**
* **8 files define kernels and never launch them — `erm_tether` (7), `inner_mechanics` (21),
  `backbone_warp` (1), `hand` (3), `envelope` (8), `branch_angle_warp` (1), `presets` (1),
  `membrane_surface` (5), together 2,114 lines / 47 kernels.** These are LIBRARIES: `inner_mechanics`
  defines the 21 kernels `driver.py` launches. Separating a kernel's definition from its launch site is
  the intended architecture, not a defect, and calling them "kernel-0 modules" was my own bad framing.
* **69 files define no kernel at all** — constants, dataclasses, contracts, helpers. Ordinary dependency.
* The two genuinely deferred ones (`implicit_mechanics` 2,093 lines, `contact_schwarz` 932) were closed by
  measuring that removing them from the driver deletes **ZERO** closure files, because `erm_schwarz`,
  `erm_gauss_seidel`, `medium_exterior` and `sf_implicit` already import them.

**And the cost of all of it is 0.026%.** warp 0.176 s + driver 0.425 s = 0.6 s against a 2,322 s native
run. Deleting every importable line would save a quarter of a percent of a percent. **There is no speed
problem here to solve** — the reason to shrink the closure was never speed but PROVENANCE: every file in
it can silently change a run, must be rsynced, and is hashed by the build guard. 129 → 120 bought surface
area, not time, and saying otherwise would be inventing a benefit.

**A number in this file needed fixing.** The line above says "of 101 imported modules 20 fire a kernel
(22%)". That is the RUNTIME set (what reached `sys.modules` and actually launched); the 120 / 43 / 36%
here is the STATIC AST closure. Both are correct for their own definition and they are not comparable —
the earlier line never said which it was, which is how two true numbers become a contradiction for the
next reader.


## Moved out of `STATE.md` (f) on 2026-08-09, to make room for the merge entry

Both are 2026-07-29 status, superseded twice over. The **MEASUREMENT RULE** that sat between them
was deliberately left in `STATE.md`: it is a standing limit on every sweep, not a record of one.

**Blebbistatin sweep: NO R, and now no γ either.** Five points STATIONARY (γ 4.33/3.48/2.88/12.24/**49.80**
pN/µm at `k_xb` 1000→1); the γ ∝ `k_xb`⁻ⁿ exponent flips sign while every force FALLS. The cause was already
in the repo — `working_stroke_strain`; **only `k_xb`=100 is physiological** and
**nothing consulted it, including me**. A declared `dt` test refuted time-convergence → **(c) 17**. No R:
(e) 2. ⚠ **Five of my claims died here** — `outputs/ac/bleb_settled/REPORT.md`.

**Import cleanup is FINISHED and cost 0.026%** — 0.6 s of a 2,322 s run, so 129 → 120 bought PROVENANCE
surface, not time. Of 120 files: 43 launch kernels, 8 define them for others (libraries), 69 are constants.



## Moved out of `STATE.md` (f) on 2026-08-09 — the merge entry's detail

`STATE.md` hit its 20,000 B gate three times while this was being written, which is the gate
working. The pointer stays in (f); the prose is here, and the authoritative version of all of it
is `docs/decisions/PROPOSAL-the-merge-and-the-rename.md`.

**2026-08-09 — renamed Project Aleph, split into layers. No physics changed; no row above moved.**
`ffn_sim`→`aleph`, `ac/`→`engine/`+`components/`, `ff/`→`laws/`, `viz/` and `common/` resolved.
Decisions and every reversal: [`docs/decisions/`](docs/decisions/). Tree: [`STRUCTURE.md`](STRUCTURE.md).
⚠ **C-2 now precedes C-1** — PI, cited in the merge record §C.

⚠ **Six archive decisions were reversed by measurement** — `ac/cell`, `virtual_cell`, 23 of `ff/`,
`dcm`, `gamma_floor*`, `architecture_metrics`. Each rested on a document calling the code frozen or
unused; the code disagreed every time. Merge record §1a–§1i.

⚠ **Three findings, ratcheted not gated** (widening any blocks commits — PI calls): `gamma_floor_sweep`
defaults `device="cuda:0"`, which the charter forbids; `dcm/` imports its own oracles (the v1 pattern
the charter inverts); and **(c) reaches less far than it looks** — 0 unmarked repeats where it scans,
**46** in `aleph/docs/**`, and only **7 of 17** rows carry a fingerprint.



## Moved out of `STATE.md` (f) on 2026-08-09, second pass — the 07-29 PRIOR

Displaced by the resting-residual tier-(a) row. The engine-cortex A/B and the T2 arm are still true and
still cited by `outputs/ac/engine_cortex_ab/REPORT.md`; what moved is the paragraph, not the result. The
**MEASUREMENT RULE** that followed it stays in (f), because it is a standing limit on every sweep rather
than a record of one.

PRIOR: 2026-07-29 — **`ac/engine` computes the cortex in a production run for the first time**, behind
`--engine-cortex`, after PI approved **D2**. **Arrays are still shared** ((c) 4), so this does NOT end
co-location. **T2 and the A/B both PASS at full native** — T2's value is the arm D2 unblocked: mask OFF fails
by an EXACT factor of two, which makes the 2.80e-16 parity evidence of exclusivity rather than of a comparison
that cannot fail. Both are same-seed, so the scatter below does not touch them.
`outputs/ac/engine_cortex_ab/REPORT.md`.



## Moved out of `STATE.md` (f) on 2026-08-11 — two 08-09 entries, to make room without asking

The boot budget was raised to 22,000 B on 2026-08-10 with the PI's condition attached: *"Direction is still
down — move status out of the boot set before asking again."* `STATE.md` stood at 21,997 B, so the
connector device-run census row was landed by moving these two rather than by asking for a third raise.

Both are history and neither is load-bearing where it stood. The first is redundant with the tier-(a) row
it announces; the second is fully described by `STRUCTURE.md` and `docs/decisions/`.

PRIOR: 2026-08-09 (engine) — **the resting-residual result is now a tier-(a) row above**; the budget is
not the lever and the operator is what is left.

PRIOR: 2026-08-09 — **renamed Project Aleph, split into layers. No physics changed; no row moved.**
`ffn_sim`→`aleph`, `ac/`→`engine/`+`components/`, `ff/`→`laws/`. Tree: `STRUCTURE.md`; decisions, six
reversals and three ratcheted defects: `docs/decisions/`. ⚠ **C-2 now precedes C-1** (PI, merge record §C).
⚠ **Collection widened** 2,385→2,795, exposing **9 pre-existing failures**, xfail-strict.


## Moved out of `STATE.md` (f) on 2026-08-11, second pass — the 08-10 card5 entry

Displaced by the interface-residual tier-(a) row. The result itself is unchanged and still has its own
row above; what moved is the (f) paragraph, which restated it.

PRIOR: 2026-08-10 (card5) — **the interior column's three surfaces address DISJOINT blocks at full
native**. ⚠ **SLICE VIEWS, so (c) 4 is NARROWED, not closed.** Two guards were green for the wrong reason
and are fixed (pointer- vs byte-identity; a launch over the global array indexing per-component weights).

## 2026-08-11b — two blockers found by trying to bind, not by reading

**(e) 5 — the SF geometry refuses a motor station, in writing.** `SFArcPopulation.motor_station_census()`
returns `N_stations: 0` with every bundle excluded and states why:

> *"curved sarcomere (mid is an arch apex, not the chord midpoint): the two filaments are not collinear over
> the shared span, so a rigid bipolar minifilament cannot straddle them"*

`excluded_by_class` = ventral 8, dorsal 4, transverse_arc 4, perinuclear_cap 4 — all of them. So
`nmii_sf_motor` has been UNWIRED not because nobody wrote the connector (the class has been target-agnostic
since the cortex edge, and its docstring names this edge) but because the SF geometry as built cannot carry a
rigid bipolar minifilament. **The connector, the `sf_arc` port and head exclusivity all landed anyway**
(`b4095ace`) and are reusable the moment the geometry question is answered.

The decision is a model one: straighten the sarcomere over the shared span, or adopt a minifilament that
tolerates a curved one. **It must not be forced.** A proximity partition alone would have borrowed CORTEX
minifilaments — and the run says so plainly: with the NMII population placed for a 70,686-filament cortex,
every minifilament is nearest the cortex, the partition assigns 0 heads to `sf_arc`, and the "None stays
None" rule correctly registers nothing rather than a connector that cannot fail to balance.

**(e) 6 — four connectors are blocked by force-channel granularity, not by connector code.**
`membrane_erm_cortex`, `actin_cap_linc`, `membrane_cortex_contact` and `nucleus_cortex_contact` each have
both endpoints bound in the composed world AND a runtime class. What blocks them is that the incumbent
already computes them, inside `membrane_compartment_force` / `nucleus_compartment_force`, and
`OMITTABLE_CHANNELS` carries no finer entry — so binding the engine connector would double-count, and
omitting the whole compartment force would remove far more than the edge.

The kernels are visible in the step profile (`erm_tether_force_kernel`, `linc_tether_kernel`, both at 44
launches per step), so this is not speculative. Splitting those channels is the SAME mechanism D2 ratified
on 2026-07-29 — additive, bit-identical default, per-component — applied one level further down.

## Why `k_xb` is the wrong blebbistatin knob (moved from STATE.md (e) 2, 2026-08-11)

Only one dose sits inside the model's working-stroke window; `k_xb` also sets lever-arm stiffness, so the
sweep varies two things at once; and the bound fraction never moves across the whole range. Blebbistatin
lowers **duty ratio and `f_stall`**, neither of which `k_xb` reaches. Sweeping `f_stall` is the right knob
and needs PI sign-off, because `f_stall` is itself a PI-GAP (0.5 pN PROVISIONAL) — the sweep would vary an
unsourced constant. The instruction (**do not re-run the `k_xb` sweep**) stays in `STATE.md`; this is the
reasoning behind it.

## Moved out of `STATE.md` (f) on 2026-08-15 — the two 08-10/08-11 connector-census paragraphs

Displaced by the arena pre-study and the two T7 motor measurements. Both ⚠ lines were kept in `STATE.md`;
what moved is the narrative around them. Neither result is retracted.

PRIOR: 2026-08-10 (engine) — **connectors close at 38/38, MEASURED.** The last eleven edges were **two**
force laws, not eleven jobs; `connector_runtime_census()` maps every edge to its class and its gate
IMPORTS each symbol, catching **3 edges counted as served that could not have bound**.

PRIOR: 2026-08-11 (engine) — **"38 objects" and "1 object" are both true, of DIFFERENT worlds.**
Neither that line nor (c) 4 named a world; both stand, (c) 4 STRICTER never looser.

## 2026-08-15 — the arena pre-study, in detail

**Why these three ran before any arena code.** The PI directed a root redesign: one fixed-capacity node
arena claimed as contiguous ID ranges, seven primitives, interactions split four ways, and — the load-bearing
ruling — **every parameter a swept axis, because a value labelled physiological is physiological for some
particular cell and is not certain**. Each of the three measurements below could have changed the arena's
array list, so each ran first.

**W2 — capacity.** Ten million nodes plus segments, angles, bonds and faces allocate on the 4090 in five
stages with no OOM, ~1.53 GiB arithmetic against 23.99 GiB. ⚠ The driver memory read returned nothing, so
`exact_peak_gpu_bytes` — `null` everywhere in this repository — is STILL unmeasured. The allocation
succeeding answers "does it fit"; it does not answer "at what peak".

**W1 — conditioning.** `derive_gershgorin_lambda_max` (`overdamped_relax.py:148-161`) builds a per-node row
sum array and returns only its max, discarding the distribution. Kept, it reads 11.2 → 31,560,034 pN/µm:
max/min **2,827,026**, p99/p1 **19.8**. An explicit relaxation takes its step from the stiffest node while
the softest mode relaxes at the softest rate. ⚠ Gershgorin bounds λ_max but the minimum row is not a bound
on λ_min, so this is an indicator, not the condition number.

**A2 — the hypothesis that failed.** `hand_kmc.py:141` warns that relaxers over-stretch at the lit-anchored
crosslink stiffness "WITHOUT crosslink_turnover", and `gamma_floor.py:335-341` gives the mechanism: without
it, stiff crosslinks amplify the bending-settle residual into a large spurious passive γ_xl. Neither
`ac_resting_residual_curve.py` nor the incumbent driver references that flag. So: are the crosslinks
pre-stretched at t0? **No.** Extension is exactly 0.0 at every percentile across all 1,413,720, with zero
degenerate pairs — the `ac/` builder binds force-free by another route. The standing negative result is not
a pre-stress artefact.

⚠ **Incidental correction.** The median crosslink `k` is **820,000 pN/µm — FILAMIN** (Ferrer 2008 companion),
not α-actinin at 4.6e5. The cortex crosslink population is mixed and its median member is filamin. A claim
made earlier the same day that it was α-actinin is wrong.

**What the three together eliminate.** Two of the three candidate causes for the standing non-descent are
gone. It is not a pre-stress artefact (A2). It is not non-conservativity: at rest no head is bound, and the
inner solve freezes kinetics regardless, which is consistent with the tier-(a) row recording the assembled
tangent as **symmetric to round-off** while the force field ACROSS steps is not — the non-conservativity
lives in the irreversible abscissa advance, not in the instantaneous tangent. What remains is conditioning,
and the partitioned solid–fluid coupling, which separate with one run at zero pressure.

**The physical picture that leaves.** The cortex at t0 is force-free but nearly rigid. The residual lives in
the soft modes — bending, steric, pressure — and relaxing them requires moving nodes against crosslinks at
8.2e5 pN/µm. The softest mode that must relax is opposed by the stiffest element. That is physics, not a
defect, and it is what an IMEX split or a constraint treatment exists for.

## 2026-08-15 — the two T7 motor measurements, in detail

**M-A.** T7 derived that the accepted-step stall cap `da_max = (f_stall − F)/k_xb` overtakes the Hill advance
below a critical tick. Measured on the production kernel with no parameter changed: Δt* = f_stall/(k_xb·v₀) =
**4.1667 ms** against the drivers' documented **10 ms**, a factor 2.4 over. At that tick v_eff(F=0) = **0.050
µm/s** against v₀ = 0.120, and 0.050 = f_stall/(k_xb·Δt) to the digit — the tick sets the unloaded velocity,
not the motor. The crossover has a closed form: the (1−f) factor cancels on both branches, leaving

    f* = κ (Δt/Δt* − 1)

so f* ≤ 0 ⟺ Δt ≤ Δt* (Hill everywhere, admissible), and f* ≥ 1 ⟺ Δt ≥ Δt*(1+1/κ) = **12.5 ms**, above which
**Hill never runs at any load**. At the drivers' tick f* = 0.700 exactly, both branches equal 1.5e-4 µm there,
so the 70% is the true intersection and not an artefact of sampling the load at tenths. The peer `engine`
session verified the crossover in closed form independently and read `hand.py:198-200` directly.
⚠ f_stall is itself a PI-GAP at 0.5 pN PROVISIONAL, so Δt* rests on an unsourced constant — but the margin is
2.4×, and k_xb also sets lever-arm stiffness, which is why this cannot be fixed by moving k_xb.
**Framing that matters for whoever fixes it:** `hand.py:195-196` says the cap exists to stop the settled load
overshooting f_stall within a tick — a numerical guard. Nothing is wrong with the cap; the tick is too coarse
for it to remain one.

**M-H.** T7 called the per-head transverse residual "the strongest available test that the nmii↔sf_arc
connector is real rather than nominal", and noted it was absent. Measured: the transverse force on a bound
head is **0.0000 pN at every offset out to 1.0 µm** — 4.8× the 0.21 µm capture radius it could never have
bound outside of — and the **axial load does not change either**. Newton-3rd residual 0.0 pN at every point;
positive control (axial arm) gives exactly −k_xb·δ. `segment_motor.py:495-497` is why: only ⟨d, ŵ⟩ enters the
magnitude and ŵ alone sets the direction, so the transverse part of d enters neither. The bond is not merely
without a restoring force — it is blind to transverse displacement. Scale, as a counterfactual only: a true
3D spring would give 100 pN at 0.1 µm against f_stall = 0.5 pN.

**(e) 5 correction, surfaced but NOT applied — it is a PI item.** `motor_station_census()` returns ONE
hard-coded `exclusion_reason` naming curvature, for every exclusion regardless of cause.
`sf_population.py:566` shows `motor_station_ready` is a CONJUNCTION of two conditions and the string names
only the second. Measured both ways on the charter env: the default call returns `N_stations: 0`; supplying
`sarcomere_lateral_um=0.400` and `sarcomere_overlap_um=0.301` — both DERIVED from the minifilament's own
geometry — returns `N_stations: 1`. VENTRAL, DORSAL and TRANSVERSE_ARC all build `mid = 0.5*(a0+b0)` and are
straight; only PERINUCLEAR_CAP sets an arch apex. So on the (e) 5 census, **16 of 20 exclusions were
geometrically straight** and were excluded solely by the missing parameter, while the diagnostic told all 20
the same story. This changes whether `nmii_sf_motor` is blocked by a MODEL decision or by a build argument —
which is the PI's call, not a session's.

## 2026-08-15 — the solver diagnosis CLOSES on the operator, in detail

**What this session was for.** Two candidate causes for the standing resting non-descent survived the arena
pre-study: conditioning, and the partitioned solid–fluid coupling. Three measurements were directed to
separate them. Both are now resolved, and neither was what its name said.

**(2) The Gershgorin number was never a condition number, and the operator is not merely ill-scaled.**
`ac_observe_sf_operator_native.py` run at `n_ventral=1` with the straddle geometry — which is what turns
`motor_station_census()` from `N_stations: 0` to `1` — gives 18 SF nodes + 34 NMII particles = 52 nodes,
**156 DOF**, and the dense tangent assembled by 156 column probes:

| | t0 (nothing bound) | bound (20/20 heads) |
|---|---|---|
| λ_max | 450,642 pN/µm | 453,011 pN/µm |
| true κ | **4.507e7** | **2.284e12** |
| null space | 35 / 156 | **40 / 156** |
| bond-graph prediction | 6 | **0** |
| relative asymmetry | 1.936e-17 | 6.241e-17 |

The zero threshold is `n_dof·eps64·‖K‖` = 1.569e-8 pN/µm — a machine floor, not a chosen tolerance
(`spectrum.py:197`), and κ is taken over the NONZERO modes (`:208`). The symmetry column CONFIRMS the
committed tier-(a) row at a second population rather than contradicting it. So does the minifilament block
probe: 102 DOF, 19 internal zero modes, 95.6% on the heads — that row claimed the property transfers, and
it transferred.

**⚠ The 2,827,026 max/min row-sum span is an indicator, not a condition number.** Gershgorin bounds λ_max;
the minimum row bounds nothing. At the one scale where both can be computed the true κ is five orders larger.

**(2b) The obvious explanation is REFUTED, and the refutation is the useful part.** 40 = 20 bound heads × 2
transverse directions matches the committed M-H result arithmetically. Measured: the predicted per-head
transverse directions carry a Rayleigh quotient of **100.042 pN/µm**, which is `k_head_arm = 100` — an
identification, not a near miss. Containment residual against the measured null space 0.985 max / 0.493
median. **A bound head is not transversally free; the head-arm spring holds it.** M-H's "no transverse
channel" is a property of the CROSSBRIDGE BOND, not of the head's mobility — `9ce447d6` scopes that
sentence in the original record rather than rewriting it. The null space is collective: weight 0.524 head /
0.367 backbone / 0.109 sf, min 0.0014 on the whole NMII block, so at least one mode lives almost entirely
on SF nodes.

**What is subtractable is 3, not 6 and not 12.** The IMEX operator carries a Dirichlet mask, but
`_make_matvec` applies `solver._stiffness()` and holds regularization at zero — its docstring says "the
probe sees K and not aI + K" — so the assembled matrix is raw unmasked K and the mask contributes no zero
modes at all (Cartesian basis at pinned DOF: containment residual 0.990). Splitting rigid-body motion,
which `_rigid_body_basis` QRs together: **translation** |Rayleigh| 9.746e-13, containment 2.30e-5 — exact
zero modes; **rotation** 1.695e-2, containment 0.564 — NOT null, which is correct physics away from
equilibrium. **Anomaly = 37 of 156.** The rotation stiffness is also an independent measurement, from the
tangent alone, that the configuration is not relaxed.

**(2c) No single stiffness family owns them.** Re-assembling with each of the eight families scaled to zero
raises the zero count by at most +4 (`sf_axial`); `sf_bending` and `alpha_actinin_arc` by 0,
`nmii_backbone` +2, `nmii_angle_backbone` +1, `nmii_head_arm` / `nmii_angle_arm` / `crossbridge` +3. So the
37 are unconstrained by all eight at full strength — a property of the assembled MODEL, not of a disabled
term. `sf_bending` is genuinely wired and nonzero (14 triples, α 0.1296), so it is not a build error either.
λ_max is set by `nmii_angle_arm` (447,514 of 453,011). ⚠ `alpha_actinin_arc` contributes exactly 0 because
`n_arc_joints = 0` as built — the term is UNTESTED here, not harmless, and that zero is about the
configuration while reading like it is about the code. Same class as the (e) 5 exclusion.

**(3) IMEX descends where explicit does not, and still does not close.** Same slice, 10 steps × dt=0.01,
only the inner solver differs; the ceiling is `SF_MOTOR_VOID_CEILING = 0.01`, declared at module scope
before any run.

| explicit `--relax` | 50 | 200 | 1000 | 4000 |
|---|---|---|---|---|
| inner iterations | 750 | 3,000 | 15,000 | 60,000 |
| residual/signal | 2.021 | 1.290 | **3.375** | **3.282** |

| IMEX `--newton` | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|---|---|
| inner iterations | 12 | 24 | 48 | 96 | 192 | 384 | 768 | 1,536 |
| residual/signal | 2.147 | 0.853 | 0.459 | 0.285 | 0.214 | 0.133 | 0.0632 | **0.0596** |

**Explicit is NON-MONOTONE** — 20× the budget ends higher. That is the same signature as the standing
full-native negative result, reproduced at 156 DOF, so that failure is not about scale and a seconds-long
discriminator now stands in for a Slurm job. **IMEX halves the ratio per doubling and then KNEES at
~0.06**, six times above the ceiling; 64→128 buys 5.7%. Every arm is VOID, correctly, so no tension
magnitude from these runs is quotable. IMEX descends on a singular operator because the `aI` term
regularises the null space by construction — that is convergence, not force balance. And λ_min over the
nonzero modes is 1.98e-7 pN/µm, only **12.6×** the machine floor, so there is no clean kernel to project
out: any fix targeting exactly 37 or 40 directions leaves a tail numerically indistinguishable from what it
removed, which is the shape of the knee.

**(1) The solid–fluid coupling is EXCLUDED; the turgor LOAD was the whole thing.** Full native, 70,686
filaments / 511,114 nodes, one outer step, 4000 explicit inner iterations. The criterion was written down
before the data (binary: DESCENDS iff candidate < start, no threshold):

| arm | residual_start → candidate | | `n_biot_subcycles` |
|---|---|---|---|
| Π₀=40, Biot ON (physiological) | 47.405172 → 47.471519 | +0.1400% | 67 |
| no Biot, no load | 3.405036 → 3.403676 | −0.0399% | 0 |
| **Π₀=0, Biot ON** | **3.405036 → 3.403676** | **−0.0399%** | **67** |

The two-arm A/B pointed at the coupling and was confounded: `--no-pressure` removes the Biot substrate AND
Π₀ together, and `residual_start` falls 47.4 → 3.4, so **93% of the resting residual IS the turgor load**.
The third arm resolves it — Π₀=0 with the fluid fully composed reproduces the no-fluid arm to the digit.
**The vacuity check is what makes it readable**: if Π₀=0 had also stopped the Biot solver, "Biot ON" would
have been nominal. It did not — 67 subcycles, the same as the physiological arm.

⚠ **−0.04% is a sign flip, not convergence.** `inner_converged` is False and the step is rejected in every
arm, at Π₀=0 as much as at 40. Removing the load stops the residual going backwards; it does not close the
solve. And Π₀ = 40 Pa IS the physiological operating point, so what the solver cannot relax is the real
load, not an artefact.

**Π₀ is now a declared axis** (`CellConfig.turgor_pi0_pa`, `--turgor-pi0`), resolved once per build so the
grid seed, CFL term, preload and ledger cannot disagree; the driver's membrane osmotic balance reads the
ledger rather than the module constant, because a build seeded at one Π₀ balancing against another would
leak flux the arms then differed by. Default `None` = the gated `PI_0_PA`, verified byte-identical
(delta exactly 0.000e+00, not round-off).

**Four defects found, all four fixed.** `assemble.py` handed out an unbound `nucleus_mask_provider` whenever
pressure was off, making the exact configuration (1) needed unrunnable; `forces_manifest.py:289` still
pointed at the pre-rename `ac/cell/driver.py`, so the invariant its own docstring calls "what stops the
observer from rotting" had been DEAD since the rename and blocked every `observe` artifact write;
`ac_gate_b_sf_motor_native.py:508` cast `float()` on four PI-GAP optionals that default to `None`, so the
DEFAULT SLIP path could not start (`b0b2738b`); and `aleph/outputs/tag_kb/` was absent on the run host, so
gates refused with `GATE-CONTRACT MISSING` — **indistinguishable in a transcript from a gate that ran and
passed**, and still only partially populated by one session's copy.

**What is left.** The operator is the only cause still standing for the resting non-descent, and the 37 are
the sharpest thing known about it. Neither (1) nor (3) reaches them: they are measured on a slice with no
fluid, and IMEX regularises rather than resolves them. Whether they are the ordinary transverse floppiness
of a network whose only stiff terms are axial and angular — in which case the NF2007 projector this engine
already uses for backbone inextensibility is the first candidate that is not a hidden finding, a new
PI-GAP, or a double count — is NOT measured and is not asserted here.

## 2026-08-19 — carried in from `codex/cortical-tension-afm`, which was never merged

The AFM cortical-tension line wrote its own `STATE.md` (f) entry on 2026-08-10 and then sat unmerged for
nine days on a worktree. Its entry does NOT fit in `STATE.md`: that file stands at 21,998 B against a
22,000 B gate, and `CLAUDE.md` forbids moving a threshold to make something fit. So the entry lands here,
verbatim, and (f) keeps the one line it already had (*"2026-08-10 outer: NN/report v3 ready; quantitative
AFM BLOCKED"*) — which records the BLOCK but not the acceptance.

> **2026-08-10 (AFM) — C-2 ACCEPTED; quantitative sweep blocked.** Restored membrane subdivision 8.
> Full-native RTX 3090 seeds 0/1/2 each committed 0.01 s at maxPF
> 0.209269/0.205555/0.208651 pN < ~0.210669 pN. Indenter/transaction CUDA job 56 PASS;
> rejected candidates expose no partial curve; aggregation requires ≥3 seeds. Driver permits passive demo
> only. No F–δ/γ: physical inputs and post-induction dt-converged active-cortex handoff remain open;
> (c) 3/17 stand.

⚠ **NOT ratified by being merged here.** "C-2 ACCEPTED" is that session's own stamp on its own run; no
session has read it against the working branch's evidence, and merging a branch is not review. Whether it
earns a tier-(a) row — and therefore whether something else leaves `STATE.md` to make room — is a PI call
that has not been made. Until then it is quotable as what it is: an unreviewed entry from a merged branch,
with its artifacts under `aleph/outputs/ac/afm_cortical_tension_readiness/`.

## 2026-08-20 (restructure) — the acceptance predicate, and two stationarity modules

Would belong in `STATE.md` (f). It does not fit: that file stands at 21,998 B against a 22,000 B gate and
`CLAUDE.md` forbids moving a threshold to make something fit. Full argument:
[`v2_audit/ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md`](v2_audit/ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md).

**The step-acceptance predicate tests roundoff.** `ledger.py:250` accepts when
`|Σreaction + Σtraction|² ≤ n_terms·eps64` — the Higham summation bound. That is adjoint/wiring closure,
true at any configuration where the connectors are attached. `inner_converged` is computed, recorded, and
read by no gate. Census: 17 verdict strings, an 8-rung ladder, 3 evidence tiers, 7 ratchets — **one**
predicate gating a physical step, **zero** criteria stating what the resting cell is.

**And there is no target to converge to.** Two resting bound-myosin constants are PI-GAP with no default
(`driver.py:1592, 1597`), so no active hoop tension opposes turgor; the driver says so in its own docstring
(`driver.py:260`); explicit relaxation is non-monotone at 156 DOF and 1.65 M DOF alike. The five
2026-08-15 eliminations stand; the inference *"the operator is the only cause left"* does not, and was
already PROVISIONAL.

**Two `stationarity.py`, disagreeing on 30.7% of series.** `tau_int` is identical (ratio 1.00000), but the
verdict contracts differ in opposite directions — `min_windows` 5.0 vs 50.0, `drift_sigma` 1.0 vs 3.0.
Over 420 AR(1) series: 129 disagreements, dominant mode **engine STATIONARY / observe DRIFTING** — the
permissive direction, on the module the production driver imports. Pinned by
`tests/architecture/test_stationarity_modules_agree.py`. **Not quotable as a physics result** — it is a
property of two Python modules, measured on synthetic AR(1) series, not of any cell.

## 2026-08-20 — moved out of `STATE.md` (a) to make room for the arena promotion

`STATE.md` stood at 21,998 B against its 22,000 B gate, and `CLAUDE.md` forbids moving a threshold so
something fits. Two rows that had become history rather than status were moved here instead. Both are
closed; neither is load-bearing where it stood.

PRIOR: ~~`aleph/common/`~~ — **Dissolved 2026-08-09.** `surface_manifold`, `turgor_pi0`, `compartments`
→ `laws/`; `sim_realtime` → `dcm/`; four modules nothing imported → `archive/`.

PRIOR: `aleph/scripts/` — the bare-`ac.*` import defect (22 scripts double-loading the package) was
**fixed 2026-07-29; zero remain**.

---

## Moved out of `STATE.md` (e), 2026-08-21 — CLOSED items in an OPEN queue

⚠ Moved, not deleted. `STATE.md` hit its 22,000 B boot budget with 17 B of headroom, and the
charter's stated direction is *"move status out of the boot set"* rather than raise it again
(it was already raised 46,000 → 48,000 on 2026-08-10, forced by a 16th tier-(a) row). Closed
decisions in an open-decision queue are the clearest case: every session pays for them at boot
and none of them is a live question.

**CLOSED 2026-08-10 — the boot budget** (was item 3, proposed 2026-07-28): TOTAL **46,000 → 48,000**,
`STATE.md` **20,000 → 22,000**. Forced by a 16th row against 41 B of headroom. **Direction is still down**
— move status out of the boot set before asking again. Why (b)/(c) may not move behind a link: `STATE_LOG.md`.

**CLOSED.** The eight D1–D8 decisions were ratified and implemented 2026-07-28; `AC_EXECUTION_PLAN_2026-07-25.md`
§9 is the per-decision table. Separately, **`COMPARTMENT_VALIDATION_TRACKS` §6 D2 was approved 2026-07-29**
(`_accumulate_all(omit=)`, additive, bit-identical default), unblocking T2 and four validation tracks. ⚠ **Two
different decisions are both called D2** — §9's is membrane area-stiffness, already landed. Always name the
document. GC-0001 (the `observe` predicate scoped to the
passive t0 configuration) is PI-signed, and exposed that `observe` has no `GATE_*.yaml`, so that row pins no hash.
Per-parameter gaps: `PI_GAP_EVIDENCE_CARDS_2026-07-25.md` (33 cards + Card G-A); ~50 per-track questions in
`_raw_2026-07-25/compartment_plan_tracks.raw.json`.

---

## 2026-08-22 (world) — the archival ruling, and the two edges outside the guard's scope

**PI ruling, in session:** *"아레나가 엔진이고, 다른 거는 모두 이제 아카이브화된 상황이다"* — the arena is the
engine, and everything else is now archived. `STATE.md` (a) is the single place that says which tree is
what, so that is where it landed: `aleph/engine/` and `aleph/components/incumbent/` go **PORT SOURCE,
frozen → ARCHIVED**.

⚠ **`laws/` and `dcm/` are NOT covered by the ruling and were not touched.** Both exemptions are measured,
not preferred: the arena parity ran two `laws/` kernels over arena-addressed arrays unchanged, and
`laws/relax.py` calls `dcm`'s `implicit_overdamped_step`. `STATE.md` (a) states both.

⚠ **THE RULING'S REAL COST IS IN (b), NOT (a).** Most tier-(a) rows were measured by the two trees just
archived. They do not become false and they are not retracted — but they are now **history taken on an
archived tree**, and only `world/` re-earning them makes them status again. Today `world/` owns four:
builds-inside-membrane, bitwise determinism, 0.1620 s/step, and the PHASE 4 gate that accepts nothing.

### Moved out of `STATE.md` (a) to pay for the change

The boot total stood at 47,880 B of 48,000 and the ruling's own wording cost 217 B. `CLAUDE.md` forbids
moving a threshold so something fits, so detail that the ruling itself turned into history moved here:

PRIOR, `components/incumbent/`: *"`engine/` imports it **243 times**. Its plumbing is what 2026-08-20's L0
runs failed in: the ERM preload applies 39→79 pN along the axis carrying 92% of the residual and the
residual does not move. That defect is being ported PAST, not fixed."* Unchanged as a fact; it is simply
no longer a statement about a live tree.

### What the ruling exposed: the guard's scope is the PACKAGE

`world/**` was already clean and `test_layer_directions` had enforced it since `b203413d`. But that guard
reads `aleph/<layer>/**`, and **two files just outside the package were still importing the archived
tree**:

| | edge | how long |
|---|---|---|
| `aleph/tests/world/test_bond.py` | `from aleph.engine.forces_manifest import PROVENANCE` | since `b203413d` — that commit MOVED the symbol to `units/` and left a re-export, so the old spelling never stopped working |
| `aleph/scripts/world_phase1_native.py` | `from aleph.components.incumbent import gpu_memory_accounting` | the arena's only native PHASE 1 driver, since the driver was written |

Neither was hidden. Neither was noticed, because the thing that would have noticed was not looking there.
A layer's test and its driver are not commentary on the layer — they are how it is exercised and how it is
run, and **an arena whose only native driver imports the tree it replaces has not replaced it.**

**Fixed at the guard, not at the two files.** `LAYER_EXTRA_FILES` adds `tests/world/**` and
`scripts/world_*.py` to the `world` row, applied to the layer-direction test ONLY — the oracle firewall
keeps package scope deliberately, because a runtime that reads its own oracle can be tuned to agree with
it while a TEST reading one is what oracles are for. Widening both at once would have forbidden that as a
side effect. Verified in both directions: injected violations in each new scope FAIL, removing them PASSES.
A vacuity control asserts each glob still matches files — the `common` lesson, applied before it happens.

**And the widened guard immediately found that the rule was written imprecisely.** It flagged
`world_phase4_native.py` importing `scripts.run_provenance` — the module that stamps a run with the build
it ran. *"`world/` may not import `scripts/`"* is a statement about the arena PACKAGE not depending on
drivers; read against a driver it forbids a driver from importing a driver. Encoded as one condition —
**a file is never forbidden from the directory it lives in** — rather than a carve-out list, so a future
extra glob cannot need a new exception. It is a no-op for every package row.

### The NVML probe moved, and PHASE 1 gained a field it never had

`components/incumbent/gpu_memory_accounting.py` → **`aleph/world/gpu_memory.py`**. `world/build/__init__.py`
had reported `exact_peak_gpu_bytes: None` behind a status string that named the module and called the port
*"an open item"*; the port is done and `build/` now reads it. No re-export shim was left, unlike the
`PROVENANCE` move: the test monkeypatches `_load_nvml` on the module object, and a shim would take the
patch on itself while the function resolved the name in the real module — an inert patch against live NVML.
Four import sites in the whole tree, rewritten. The ledger KEY names (`gpu_memory_accounting_*`) are NOT
renamed: they name a record format, not a module.

⚠ **`exact_peak_gpu_bytes` is READ now, which is not the same as EXACT.** NVML accounting mode must be
enabled by the driver before the process starts and this repository never takes that authority, so on a
host where it is off the field stays null with an explicit non-exact status. A null there is a HOST
configuration fact, not a measurement that failed.

### ⚠ A PHASE 1 field that could never have been filled

`world_phase1_native.py` read `getattr(res, "peak_bytes", None) or getattr(res, "max_memory_usage_bytes", None)`.
**Neither attribute exists on `NvmlAccountingResult` — the field is `exact_peak_bytes`** — so both `getattr`
calls returned their default and **`exact_peak_gpu_bytes` was unconditionally `None` on every PHASE 1 run
ever recorded, including on a host where accounting was enabled.** It did not read as a bug: the note
printed beside it was the `repr` of a result object that carried the number. The charter requires exact peak
GPU bytes be tracked; it has been tracked as a null that could not have been anything else. Fixed, and the
read now respects `.exact`.

⚠ **No re-run is claimed and no artifact is corrected.** Every committed PHASE 1 record still carries the
null. What changes is that a future run can carry a number.

### Verification

`aleph/tests/{world,architecture,scripts,units,ac/cell}` — two failures, **both pre-existing at HEAD**
(`ba0160f9`) and confirmed by running them in a clean worktree at that commit without these changes:
`test_retired_claim_coverage` (12 of 20 (c) rows unfingerprinted against a ratchet of 10 — (c) 19 and 20
were added without fingerprints) and `test_argparse_help_is_formattable` (`world_turntable.py:93` carries a
bare `%`). Neither is touched by this work and neither is fixed by it.

⚠ **Committed against a moving tree.** `ba0160f9` (ruling 14's second half) landed from another session at
17:13, mid-edit, and it touches two of the same files. Checked rather than assumed: `cortex_shell()` is
present at all six of its call sites and the diff against `ba0160f9` shows only additions from this work.

### 2026-08-22b — the ruling reached the two files that had not heard it (PI approved both)

`STRUCTURE.md` and `ownership.yaml` were surfaced as PI items at 17:32 and approved in the same
breath. Both are now transcribed.

**`STRUCTURE.md`** — the two rows go `✅ PORT SOURCE, 동결` → `✅ ARCHIVED 2026-08-22
(archive-readonly)`, with the reason and the cost stated above the table so the row is not read as a
migration target. ⚠ **Line and file counts are deliberately unchanged**: archiving is not deleting
and that table describes the disk.

**`ownership.yaml` — this could NOT be done by adding to `archive-readonly`.** `owner_of` returns on
the FIRST matching glob and `engine` is declared first, so a duplicate `aleph/engine/**` under
`archive-readonly` would have loaded cleanly and enforced nothing — the precise failure the file's
own header warns about, and the shape of the seven findings in
`docs/decisions/FINDING-2026-08-09-checks-that-pass-for-the-wrong-reason.md`. The globs had to LEAVE
`engine`. Verified by resolution, not by reading:

```
archive-readonly   aleph/engine/transaction.py
archive-readonly   aleph/components/incumbent/driver.py
card5              aleph/engine/cortex_state.py        <-- see below
engine             aleph/laws/relax.py                 <-- not covered by the ruling
engine             aleph/dcm/dcm_ecm_clutch_host.py    <-- not covered by the ruling
```

**The narrowing broke a declared overlap, and the test that broke had asked for it.**
`test_a_declared_overlap_is_still_really_an_overlap` pins `(path, winner, loser)` triples and its
docstring says *"if `engine` is narrowed later so it no longer covers these, the exception should go
with it — otherwise this table becomes a list of things somebody once had to think about."* It
failed on the narrowing instead of passing stale. Losers updated `engine` → `archive-readonly`.

⚠ **AND THE NEW LOSER IS A SHARPER STATEMENT THAN THE OLD ONE.** `card5` wins two files inside a
tree nobody may edit: `engine/cortex_state.py` and `engine/interior_column_slice.py`. **Left open on
purpose.** `card5` is a specific PI approval of 2026-08-09 22:04, and retiring it is the PI's call,
not a side effect of transcribing a different ruling — even though its work has landed as the
interior-column tier-(a) row (`89eb2a60`). Surfaced beside the lane, not decided.

⚠ **The archived trees' TESTS did not move.** `aleph/tests/ac/{cell,motor,ecm,fluid,solid,nucleus,
weave,emergence,engine}/**` stay in `engine`, deliberately: some of them no longer test an archived
tree at all — `tests/ac/cell/test_gpu_memory_accounting.py` tests `aleph/world/gpu_memory.py` as of
`7e9855b6`. Sorting which tests followed their subject into the archive is its own pass, and doing
it by glob today would have archived a live test.

⚠ **`engine` is now a lane whose name is wrong.** What it still owns is `laws/` — the kernel library
the arena binds — plus the tests and drivers of two trees nobody may edit. Recorded in the lane's
own comment rather than renamed, because a lane rename changes what `$FFN_SESSION` values resolve.

### Moved out of `STATE.md` (f), 2026-08-22 — made history by the archival ruling, not by being fixed

⚠ **Moved, not closed. Neither defect is repaired and both are still true of the code.** The boot total
stood at 48,094 B against 48,000 and `CLAUDE.md` forbids moving a threshold so something fits, so the
question was which line had stopped being STATUS. This one had — and by today's ruling rather than by
anyone's work:

PRIOR, `STATE.md` (f): *"⚠ 08-10/11 (engine), both stand: `ac_composed_world_dump.py` raises and **no
lane owns it**; `NOT_BOUND` is weaker than a stub — no object at all."*

Both facts are about the trees archived this morning. `NOT_BOUND` is the `engine/` connector census —
36 of 38 edges with no runtime object — and that census belongs to a tree no new work lands in;
`STATE.md` (a) already says `engine/` reached 2 of 38 device-run, which is the same fact as status.
`ac_composed_world_dump.py` is an `ac_*` driver of that engine: it still raises, `aleph/scripts/` is
still not archived, and **no lane still owns it** — `90c8617d` did not change that, because the script
sits in the 44 of 72 `aleph/scripts/*.py` that no lane claims at all.

⚠ **So the ownership half of it is NOT resolved and must not be read as resolved.** It is now one
instance of a larger open item that `ownership.yaml` states beside the `engine` lane in its own words —
*"44 of the 72 `aleph/scripts/*.py` are claimed by no lane at all … Redistributing the other 43 is a PI
governance call"*. That is where it should be answered, not in a per-file line in the boot context.

### 2026-08-22c — the reason (e) 1 was blocked went stale, in nine places at once

**What was wrong.** `world/step.py` (3 sites), `scripts/world_phase4_native.py` (5) and
`tests/world/test_step.py` (1) all said cortical tension γ *"is not emitted on the resting path"* and
that (e) 1 is therefore **undecidable as written**. γ has been emitted since `984c1292` (2026-08-21).
`outputs/ac/world_phase4/tau3_seed{1,2,3}.json` each carry `gamma.trace_pn_per_um` with **30,000
samples — in the same record as the string saying the trace does not exist.**

**Same shape as this morning's two.** A description that outlived the thing it described:
`getattr(res, "peak_bytes")` beside a repr that carried the number; a build status string naming a port
as "an open item" after the port; and now a blocker whose stated cause was removed the day before. None
of the three is a wrong computation. All three are a claim that reads as current and is not.

**⚠ What did NOT change: the verdict.** `UndefinedAcceptance` stays, every step still records
`ACCEPTANCE_UNDEFINED`, and no predicate was written. The class's own condition was *"until γ is emitted
AND the PI fixes the statistic, its window and the variance it is judged against"* — **the first
conjunct is now satisfied and the other three are not.** (e) 1 moved from UNWRITABLE to UNDECIDED,
which matters only because the second can be put to the PI and the first cannot.

**What makes γ usable, and it is not that it exists** — `GAMMA_EMITTED_e1_IS_NOW_DECIDABLE_2026-08-21.md`:
frozen, γ's spread is **exactly 0.000e+00** across six steps while the residual's is **2.111e-13 pN with
nothing moving** (atomic-add ordering, irreducible); and under a declared relaxation γ falls
**monotonically**, 4.4e-05 against a zero floor. γ has no noise floor of its own; the residual does, and
the residual is also blind — 8.2 M added force terms left it unchanged while its own tolerance grew 86×.

### ⚠⚠ And the window input is NOT ready — measured by the peer session, 2026-08-22

Recorded here because it changes what (e) 1's option (a) means, and because the obvious reading of
"widen the search bound" is **worse than the status quo**. `outputs/ac/world_phase4/tau3_cut_scan.py`:

```
seed  variant                cut     tau   n_eff  reliab  opens_on_transient  bar 50
 1    current (//2)        14625  1280.2     6.0    12.0   True                REFUSED
 1    (a) widen only       20061    97.3    51.1   102.1   False               passes
 2    (a) widen only       24519   118.7    23.1    46.2   False               REFUSED
 3    (a) widen only       23776    49.9    62.4   124.8   True  ⚠             passes
 3    (a) widen + filter   21547    88.2    47.9    95.9   False               passes
```

**Seed 3 passes the bar with the best reliability of the three, on a window that opens on a rising
transient.** Today every seed is refused and nobody is misled; widening alone would be *more
confidently wrong.* The cause is structural — `aleph/observe/stationarity.py:605-615` picks the window
by `n_effective` alone and computes `opens_on_transient` **after** the loop, so the flag exists and is
not used in the choice. **So (e) 1 (a) is TWO changes, not one:** widen the bound AND make the flag an
in-loop filter.

⚠ **And even with the filter, seed 3 is not a measurement** — its τ keeps falling (89.2 → 70.2 → 64.7 →
50.6). The flag catches *"did the window open on a rise"*, not *"has τ settled"*. Seed 1 sits at
95.5–106.5 over 19,000–21,000, an 11% spread. **1 seed of 3 yields a settled τ: a D-3 input, not D-2.**

⚠ **And the cell those series came from has no NMII** (11 populations, `bond` live 0, build `f4572b39`).
The cell the gate must finally judge is one with motors running, and that cell exists only since
`28e45747`. Whether (e) 1 is decided on the current series or on a re-measure is the PI's.

**Not touched, deliberately: `aleph/observe/stationarity.py`.** Changing the equilibration search
changes what the gate accepts, and `CLAUDE.md` forbids editing a gate contract inline or after seeing
data. The measurement is here; the ruling is not.
