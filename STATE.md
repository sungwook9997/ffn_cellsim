# STATE.md — what is canonical, what is proven, what must not be cited

**Editable at every milestone** (`CLAUDE.md` §File ownership); the single place current status lives. Six sections,
nothing else, under 22,000 B (`check_state_md.py`, raised from 20,000 by PI 2026-08-10) — the moment it grows a
narrative it becomes audit doc #255.
Source of record if a row is disputed: `aleph/docs/v2_audit/AUDIT_WHOLE_REPO_2026-07-25.md`, `AUDIT_AC_ENGINE_2026-07-25.md`,
`TRAJECTORY_HOW_WE_GOT_TO_AC_2026-07-25.md`, `cell_engine/ROLLING_ROADMAP.md`.

## (a) Which engine is canonical, and what each tree is FOR

| Tree | What it is FOR, now |
|---|---|
| `aleph/world/` | **CANONICAL 2026-08-20.** The arena: one fixed-capacity node arena claimed as contiguous ID ranges; primitives STRAND / BOND / FACE / **GRID_CELL** (eighth, PI 2026-08-21 — `arena.py`'s own closure argument already cited it). **All new work lands here.** **25 `@wp.kernel`**, **19,356 lines** vs 68,970 of port source; **11 populations + a cytosol field at native**; **takes a physical step**, 0.1620 s/step over 4.56 M nodes / 8.49 M force terms. ⚠ **No step it takes is accepted** — (b). |
| `aleph/engine/` | **ARCHIVED 2026-08-22** by PI ruling — *the arena is the engine, the rest is archived*. Reached **2 of 38 connectors device-run**; read to port FROM. ⚠ It still owns most of (b), so **those rows are history on an archived tree** until `world/` re-earns them. |
| `aleph/components/incumbent/` | **ARCHIVED 2026-08-22**, same ruling. Owns every native physics number so far. Its ERM preload defect — ported PAST, never fixed — and the 243 `engine/` imports: `STATE_LOG.md`. |
| `aleph/laws/` | **NOT archived — the kernel library `world/` binds.** The arena parity ran two of these kernels over arena-addressed arrays UNCHANGED, so this is a proven import path, not a port target. **14 of 49 modules** imported (AST, 2026-08-09). Own corpus invalidated-pending-rerun, (c). |
| `aleph/dcm/` | Parked but **NOT removable** — `laws/relax.py` calls its `implicit_overdamped_step`, found by the layer guard 2026-08-09 when the archive move broke it. Also **imports its own oracles twice** (see (f)). HOOMD parity fixtures and `tests/dcm` were DELETED 2026-07-29. |
| `aleph/validation/oracles/` | Closed-form **acceptance oracles only**, never imported by a runtime path (still clean three engine generations later). Keep it that way. |
| `aleph/scripts/` | Production drivers and gate scripts, versioned as neither. |

## (b) Tier-(a) results — native full population on CUDA, with the build each was measured on

**Rule: no row without a build commit, and none without its POPULATION** — a conclusion far below native is a (c)
entry however clean. Only self-stamping records carry a real build; every other "build" is the implementing commit —
the audit's top finding, and why the γ-floor row sits in (c).

**Each row is the HEADLINE — what may be quoted and what may not.** The full claim with its evidence,
scoping and controls is the `claim` field of [`state_rows.yaml`](aleph/docs/state_rows.yaml), which renders
this table; read it when a row is challenged. Edit the YAML and run `make state`, never the table by hand.
**Artifacts are relative to `aleph/outputs/ac/`; self-stamped builds carry `source: declared`.**

<!-- STATE-TIER-A:BEGIN -->

| Result — only what the left column says is quotable | Population / device | Build commit | Committed artifact |
|---|---|---|---|
| **PHASE 4 runs at native — its gate is EXACT (0.000% over 5 runs) about a residual that swings 1037%**, and 626x looser than faces-only. **0 accepted.** | 4,361,496 nodes; 4090 — FULL NATIVE | 2026-08-21 `world/` | `world_phase4/*.json` + `BALANCE_GATE_MEASURED_2026-08-21.md` |
| **The engine is BITWISE DETERMINISTIC for a fixed seed** — same closure digest, same gamma. | 4,558,554 nodes; 4090 — NATIVE | `ccbda122` | `world_phase4/det_[ab].json` |
| **Every node lies inside the built membrane**, enforced on built positions. Geometry; no force evaluated. | 4,588,662 nodes / 11 populations + 63³ field; 4090 | 2026-08-21 `world/` | `world_phase1/phase1_inside.json` |
| **The `nmii_sf_motor` force field is non-conservative**, and the crossbridge is why. Quote the attribution and the scaling exponents; **no magnitude**. | 904 nodes / 2,712 DOF; 320 heads / 16 minifilaments; A5000 | `db841185` (self-stamped) | `observe/sf_operator_observables.json` + `.npz` + `REPORT.md` + 4 figs |
| The assembled tangent is **exactly multilinear** in its stiffness families and **symmetric to round-off**, which licenses the analytic Hellmann-Feynman gradient. | same run | `db841185` | same artifact (`multilinearity`, `symmetry`) |
| **A built bipolar minifilament has 19 INTERNAL zero-stiffness directions**, 95.6% of that subspace on the heads. A single-object property, so it transfers. | 102 DOF / 34 particles (14 backbone + 20 heads) per block, 16 blocks; A5000 | `49d2215a` (self-stamped) | `observe/sf_operator_observables.json` (`minifilament_internal_floppiness`) + `REPORT.md` §1b |
| **The `nmii_sf_motor` lane accepts on a DEVICE PREDICATE, and that predicate can say no** — a positive control injects one head's stall force one-sidedly and the same path rejects it. **The predicate tests ADJOINT CLOSURE, not convergence** (`sf_motor_slice.py`'s own docstring says so), so a non-converged step is not what it rejects. **Not** the tension, **not** rung `CONNECTED`. | 904 nodes / 2,712 DOF; 320 heads / 16 minifilaments — **0.18% of native**; A5000 | `49d2215a` (self-stamped) | `gate_b_sf_motor/sf_motor_native_gate.json` + `REPORT.md` §7 |
| GATE-B cortex-motor: NMII binding and cortical tension **EMERGE from `k_on` events**. Quote the **mechanism** — tension is event-generated, not a static prestress. Not the magnitude. | 70,686 cortical filaments / 494,802 actin nodes; A5000 `cuda:0` | `ec409751` + `bc3fef26` (`r0_bind` attach-unstrained fix) | **None in a clean checkout.** `cell_assembled/gate_b_dynamic_v2.npz` is on disk but gitignored (`.gitignore:158`) and self-stamps `"evidence": "CONNECTED"`, which is wrong |
| Native ECM constitutive force over genuinely private device SoA. **Existence and sign only.** | 200–220 collagen fibers, **zero crosslink bonds**; A5000 `cuda:0` | `8b2e7adb` + `f47c03aa` | None — 6 CUDA-gated tests, no run artifact (the merge's "22 CUDA tests" headline is not the count; see (c) 5) |
| **Negative result, do not re-queue:** the fiber-quotient coarse operator does not close the static resting residual at full native. | 70,686 filaments / 1,413,720 inter-fiber crosslinks; A5000 | `adea8103` + `d567eaf8` (fiber-arclength multigrid speed fix: `09c14eb8`) | None — numbers live in the commits + `GATE_A_RESTING_CONVERGENCE_2026-07-23.md` §10 |
| **Negative result:** the fine-75 nm residual plateau is a **construction** artifact, not a solver wall. | fine 75 nm cortex rung (~2.9 M nodes); A5000 | `21a56220` (measurement), `3a1f8866` (fix, default-off byte-identical), `831e51a7` → `6cbefe57` (native validation) | None |
| The mixed formin + Arp2/3 cortex **builds** natively with its branch-angle harmonic force-evaluated. Quote build and wiring; **not** convergence. | 302,286 actin nodes / 70,686 fibers (6,514 long formin + 64,172 short Arp2/3, ~⅓ by mass); A5000 | `e2ddda48` (branch force wired) + `9e694ec4` / `32afc076` (overlap + unified relax) + `0dae82c4` (recording) | None — `NATIVE_QUEUE_RESULTS_2026-07-24.md` §Arp2/3 |
| GPU residency, non-vacuously: zero hot-loop D2H copies **with** a positive control at one. **Quote the residency gate only.** | 70,686 filaments / 494,802 actin + 511,114 total nodes; `cuda` | added `34a3890d`, last touched `19103665` (2026-07-20) — i.e. **before** `bc5ff3b0` and `81805535` | `foundation-hardening/native_gates.json` (`ng6` PASS) + `native_resting_70686.json` |
| The full-native cell **builds and steps** with every compartment present — **not** that the step converged. | same run as the row above | `34a3890d` / `19103665` (same run as the row above) | same artifacts: they also record `inner_converged: false`, `residual_start = residual_end = 2634.12`, `N_crosslinks = 70,686` (1 per filament — the pre-`bc5ff3b0` cortex), `overall_status: INCOMPLETE_OR_FAIL`, `exact_peak_gpu_bytes: None` |
| **The SF motor DOES build tension** once the α-actinin crosslink is out of the operator and the clock passes the induction period — first non-`VOID` verdict and first `comparable: true` cost. **Not** the tension as a physical SF magnitude. | 200 ventral sarcomeres / 10,400 nodes / 4,000 heads / 0 crosslinks; A5000 — native SCALE, but a reduced architecture and an inventory that is still a PI-GAP | `2001c119` (the topology A/B + the PASS run) + `576ebf42` (the population scaling); `source: declared` — gbook is not a git checkout | `gate_b_sf_motor/inner_solve_convergence_2026-07-28/` — 16 `run-record@2` artifacts + 7 figs; `REPORT.md` §9 (superseded) + §10 |
| **The production explicit inner step is stable at full native, and `kmax` was never a bound** — the ratio to it IS the per-node bond accumulation, measured. **Not** `λ_min`, the condition number, or any eigenvalue as a material property. | 551,434 nodes / 1,654,302 DOF; 70,686 fibers / 1,413,720 crosslinks / 8,840 heads — **FULL NATIVE**; A5000 `cuda:0` | `a44a1de3` (self-stamped, `source: declared`; the driver's whole first-party import closure was verified 79/79 against the run host before launch) | `observe/native_cortex_spectrum/record.json` + `REPORT.md` + 3 figs |
| **The cortex owns its device arrays, and they are the incumbent's cortex** — node-by-node force parity with a control that demonstrably fails. **Not** rung `CONNECTED`, and not a physics magnitude. | 70,686 filaments / 494,802 actin nodes / 1,413,720 crosslinks / 353,430 bending triples; A5000 `cuda:0` | `576ebf42` (`source: declared` — gbook is not a git checkout, so the commit is the caller's assertion) | `cortex_ownership/native_record.json` (`run-record@2`, verdict PASS, rung CUDA_UNIT / QuantitativeClaim BLOCKED, `timing` 29.7 s `comparable: true`) |
| **Negative result, do not re-queue: 66× the inner budget ends HIGHER than the start** (60 iterations −0.52%, 4,000 **+0.14%**). With the same day's uniformity result, construction and budget are both excluded. Quote the **exclusion** — not a rate, not an operator property, **not a curve**. | 70,686 filaments / 494,802 actin / 511,114 total nodes; 4,000 inner iterations; RTX 4090-1 `cuda:0` | `4c165576` (`source: declared`) — **the only row whose closure was content-checked**: 69 files hash-compared against the run host pre-launch, 0 mismatched. | `resting_residual_curve/residual_curve.json` (`kind: diagnostic`, `timing.comparable: false`) @ `2e2a2c9d` |
| **Membrane, cortex and nucleus address DISJOINT blocks at full native**, control refuses the arrangement they replaced. **Slice views, NOT private allocations.** **Not** `CONNECTED`, no magnitude. | 70,686 filaments / 494,802 actin / 536,406 total nodes; 40,962 membrane + 642 nucleus verts — FULL NATIVE; RTX 4090-1 `cuda:0` | `89eb2a60` (`source: declared`; 84/84 closure files hash-matched). ⚠ Code landed in `54755105` under a mislabelled message — see the claim. | `interior_column_ownership/` — `native_record.json` (`run-record@2`, gate, PASS, 13.2 s `comparable: true`) + 1 fig @ `23474924` |
| **The composed native cell bound 2 of 38 connectors and 3 of 14 components**; the other 36 had no runtime object. Quote the census + both controls; **no magnitude**. | 70,686 filaments / 551,434 nodes / 8,840 heads — FULL NATIVE; RTX 4090-1 `cuda:0` | `73a93883` (re-run after the fix; first pass `d58bab03`; 122/122 hash-matched) | `connector_devicerun/` — 2 × `run-record@2` (`comparable: false`) + 1 fig |
| **Both wired power ports close their adjoints** — but 1.74e-15 at **95.0 pN** and 6.40e-17 at **9.0e-10 pN**, which is not a force. The second is near-vacuous until something moves (row below). | same native cell; 8,754 of 8,840 heads bound | `0673ad9e` (isolated re-measure; first pass `ae523fa5`; 126/126 hash-matched) | `interface_residual/` — isolated + steps30/steps0 (vacuity) + `gate_controls.json` |
| **A component now advances its own positions in the accepted step, and it is LIVE** (0.0 → 0.0516 µm, same perturbation). Its bound was NOT a bound: `dt` 459× over. **No** convergence. | composed native `sf_arc` + cortex PORT + nmii + ecm; 8,754/8,840 heads; 4090-1 | `b4a0e8f4` + `a9f38bf0` + `ed646e92`; Slurm jobs 104/106/107 | `relax_bound_rebind/gate.json`; `interface_residual/perturb_relax_{on,off}.json` |
| **6/11 components step in 243 ms; the real-time ceiling is MEASURED at ~0.64.** Quote the curve + the ceiling; **no** physics. | 70,686 cortical filaments + membrane/nucleus/Biot grid; 4090 | `d226b9ab`→`58b04dc8`; Slurm 110-117 | `fullcell_timing/` — 9 records (baseline, card5, profile, 4x ninner) |

<!-- STATE-TIER-A:END -->

One tier-(b) result (reduced population), the cleanest single number here: **DCM grid-invariant single-cell spread
A/A₀ = 3.38 / 3.62 / 3.44 at subdiv 2/3/4, V/V₀ → 1.034 / 1.024 / 1.009**, `cuda:0`, build `75f95a07`, three
committed jsons matching to the digit — carrying `accel_factor = 1000`.

**What checks the two lists, and what does not.** `make kb-check` fails on an undeclared KU-tagged constant or
`REPORT.md` headline, but its inventory is YAML-only, so constants in `ac/**` and `ff/` are out.
`forces_manifest` certifies intent, not realized force — no run yet supplies a composed `cell`. Two gates are
hash-stamped; every other is not. `run-record@2` self-stamps build + t0 + verdict and a REQUIRED `timing` marked
`comparable: false` unless the verdict is PASS. **Neither the SCALE nor the COMPARABILITY stamp is enforced
against prose** — (c) 14 and (c) 15 were both caught by hand. `make boot-check` (2026-07-29) does enforce the
boot budget, a no-magnitudes rule on the charter, `MEMORY.md`'s index shape, and (c) against the memory corpus.

## (c) NOT quotable — a blocklist. Reasons and evidence: [`STATE_NONQUOTABLE.md`](STATE_NONQUOTABLE.md)

**This list is PUSHED, not looked up.** Its whole purpose is that a session sees a retraction without
thinking to ask for one — a pull-only version is re-quoted by whoever does not think to query it, which
is how it came to exist. So it stays deliberately terse: the claim, and where its reason lives. It only
GROWS; never shorten it by dropping an entry.

| # | Do not quote | Reason |
|---|---|---|
| 1 | GATE A "CLOSED" / `max\|PF\|` 0.1398 | [c1](STATE_NONQUOTABLE.md#c1) — provisional; the residual is intrinsic to discrete myosin, not a solver target |
| 2 | γ-floor "~530× under band" + every full-native `ff/` headline | [c2](STATE_NONQUOTABLE.md#c2) — retired cortex build, unrecorded, not re-run |
| 3 | γ = 3.72 pN/µm and every SF/cortex motor tension | [c3](STATE_NONQUOTABLE.md#c3) — read inside an induction period; past it the tension is real and passes (b) |
| 4 | Any `CONNECTED` rung; "first composed native cell" | [c4](STATE_NONQUOTABLE.md#c4) — wiring is complete; nothing SOLVES, so no magnitude is physical |
| 5 | "22 CUDA tests"; seam-`CUDA_UNIT`; 40-step; the ensembles | [c5](STATE_NONQUOTABLE.md#c5) — no artifact exists |
| 6 | `engine_reval` "19/21 native PASS"; P3.1 percolation | [c6](STATE_NONQUOTABLE.md#c6) — ledger unreconstructable; a gate a broken build cannot fail |
| 7 | "6/6 ECM materials real-Pa" | [c7](STATE_NONQUOTABLE.md#c7) — recovers an INPUT; exponent and N₁ sign both wrong, both inherited |
| 8 | "Physiological crawl 28–45 nm/s" | [c8](STATE_NONQUOTABLE.md#c8) — coarse only; native is 0.12 nm/s and `ac/` has no crawl gate |
| 9 | The `dcm` July claim set | [c9](STATE_NONQUOTABLE.md#c9) — no artifact, disowned, or a superseded pre-fix FAIL |
| 10 | "266k particles"; "~530× short" | [c10](STATE_NONQUOTABLE.md#c10) — retired counts, still live on public gh-pages |
| 11 | Both KB self-reports (benchmark 100%; the scale figures) | [c11](STATE_NONQUOTABLE.md#c11) — scored with the SQL it was fed; scale figures wrong both ways |
| 12 | The `lane_d` D0/D1 gate results | [c12](STATE_NONQUOTABLE.md#c12) — coarse diagnostic rigs by construction |
| 13 | ~~"HOOMD is never imported"~~ — **CLOSED** | [c13](STATE_NONQUOTABLE.md#c13) — verified by AST; only caveat is that it is still installed |
| 14 | Every population-DEPENDENT number from the 2026-07-28 `observe` run | [c14](STATE_NONQUOTABLE.md#c14) — measured at 0.18% of native; re-measurable via Lanczos |
| 15 | The `nmii_sf_motor` cost multiplier "≈20×" | [c15](STATE_NONQUOTABLE.md#c15) — confounded by solver budget, not residual |
| 16 | The 2026-07-28 blebbistatin dose response — **all six** | [c16](STATE_NONQUOTABLE.md#c16) — 0.60 s vs a required 7.14 s; all `TOO_SHORT` |
| 18 | ECM row's **"zero crosslink bonds"** as a population | [c18](STATE_NONQUOTABLE.md#c18) — about the code, not the configuration |
| 17 | **Every cortex-motor γ at `dt`=0.01 — magnitudes AND ratios** | [c17](STATE_NONQUOTABLE.md#c17) — not time-converged (8.07%), though only **3.0 σ** once judged against the measured seed scatter, not the 137 σ first claimed from within-run sem |
| 19 | Any MT tensegrity or MT-nucleus result on this cell | [c19](STATE_NONQUOTABLE.md#c19) — MTOC at the nuclear centre, 68.6% of MT nodes inside it |
| 20 | ~~NMII binding / contractility on this cell~~ **CLOSED from `ba0160f9`** | [c20](STATE_NONQUOTABLE.md#c20) — 442/442 station now; results measured BEFORE it stay blocked |
| 21 | **Every γ and τ from the PHASE 4 cell** — and the reason is the CELL, not the statistic | [c21](STATE_NONQUOTABLE.md#c21) — 2 of 11 populations carried a law; cortex at 1/10 actin; no motor in it |

## (d) Minimal read set — one clause each; the reasoning is [`aleph/docs/READ_SET.md`](aleph/docs/READ_SET.md)

<!-- STATE-READSET:BEGIN -->
1. `ROADMAP.md` — the stage gates. What each stage must PROVE; never whether it has been proved.
2. `CLAUDE.md` — the charter. Holds only what does not change; every number in it would be a bug.
3. `aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md` — per-stage working detail for the reframe.
4. `aleph/docs/v2_audit/cell_engine/ROLLING_ROADMAP.md` — accurate per-component status; the only doc that says the GATE A/B cortex is coarse.
5. `aleph/docs/v2_audit/cell_engine/CELL_ENGINE_ARCHITECTURE.md` — the component/connector contract and the 8-state ladder.
6. `aleph/docs/v2_audit/PI_GAP_EVIDENCE_CARDS_2026-07-25.md` — what is blocked on PI.
7. `aleph/docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md` — decides which older cortex results are still citable.
8. `aleph/docs/v2_audit/PARAM_PROVENANCE_AUDIT_2026-07-24.md` — SOURCED / DERIVED / CONVENIENCE / PI-GAP per constant.
9. `aleph/outputs/tag_kb/param_audit_report.md` — the disk-grounded parameter verdicts.
10. `aleph/docs/v2_audit/SUCCESSOR_HANDOFF_2026-07-25b.md` — the later of two same-day handoffs; that there are two is itself the finding.
<!-- STATE-READSET:END -->

## (e) Open PI decision queue

**2026-08-21: PI answered 16 of 20** —
[`PI_DECISIONS_2026-08-21.md`](aleph/docs/v2_audit/PI_DECISIONS_2026-08-21.md). ⚠ Two did NOT close and
both read as closed if skimmed: **6** ratified the 30-heads test point and **not** a band (both endpoints
fail the citation audit); **12** stays open, so no pre-`world/` conclusion is promoted to the canonical tree.

**OPEN — six items.** New 2026-08-11 (both in [`STATE_LOG.md`](aleph/docs/STATE_LOG.md)): **5** SF
geometry REFUSES a motor station (`N_stations: 0`) so `nmii_sf_motor` is unbindable for a MODEL reason —
do not force it; **6** four bindable connectors blocked by force-channel granularity, not by code.
1. **The gate amendment** (`ROADMAP.md`) — stationarity instead of convergence. ⚠ **No longer
   UNWRITABLE, only UNDECIDED**: γ IS emitted on the resting path (`984c1292`), has NO noise floor where
   the residual has 2.111e-13 pN, and responds. Open is the PI's alone — WHICH statistic, WHAT window,
   against WHAT variance. ⚠ The window input is **D-3, not D-2**: 1 seed of 3 settles, and the search's
   `opens_on_transient` flag is computed AFTER the choice and unused in it, so widening its bound alone
   makes seed 3 pass on a contaminated window with the best reliability of three. ⚠ Those τ series ran on a cell
   with **no NMII**. Why it was held, and the scan: `STATE_LOG.md`.
2. **The blebbistatin analogue needs a different knob. Do not re-run the `k_xb` sweep.** Sweeping `f_stall`
   instead needs PI sign-off — it is itself a PI-GAP. Why `k_xb` fails: [`STATE_LOG.md`](aleph/docs/STATE_LOG.md).
3. **`STRUCTURE.md` names directories that no longer exist** — record or repair? [e4](STATE_NONQUOTABLE.md#e4).
4. **`FFN_PROBABILISTIC_VIRTUAL_CELL_MASTER_PLAN.md`** — a lane's PROPOSED reframe to a Probabilistic
   Mechanistic Virtual Cell emitting a `CellStatePosterior`; it changed no charter, engine or gate. Its
   PROHIBITIONS are already in the charter (no signature needed); its IDENTITY claim needs one.

## (f) Last verified

LAST VERIFIED: 2026-08-24 — advisory only.

**2026-08-24 (world) — THE CELL WAS NOT RUNNING BECAUSE NOTHING DROVE IT, and the parameter denominator was never counted.** Residual at step 0 is 2.3e-13 pN (built AT rest, nothing to relax); the only motion is a 20.7 pm/step thermostat, and 30,000 steps move the cell 3.6 nm. `dt` was never the limit — it sat 227× under the stability bound. **2 of 11 populations carried a law**; `build/nmii.py` and `build/cytosol.py` run in NO production cell. Counted rather than remembered: **378 magnitudes declared in Python** under `laws/`+`world/`+`observe/` against 79 assembled by hand, and of `laws/`'s 46 modules **11 BOUND / 5 WIREABLE / 4 blocked on acceptance / 26 with no kernel to bind**. `param_ledger.py`, `law_wiring.py`, `strand_coefficients.py`.

2026-08-22 (world): the archival ruling; two edges outside the layer guard's package scope; a PHASE 1 field that was **structurally always null**. `STATE_LOG.md`.

2026-08-10 outer: NN/report v3 ready; quantitative AFM BLOCKED.

**2026-08-15 (motor) — the drivers' dt runs the WRONG LAW. WIDENS (c) 17 from γ to the law.** Δt* = **4.167 ms** vs **10 ms**; a linear-stall proxy governs **70%** of the load. The
crossbridge BOND has **no transverse channel**; the head is not thereby free — scoped `9ce447d6`.

**2026-08-15 (arena)** — ⚠ median crosslink k = **8.2e5 = FILAMIN**; ⚠ **2,827,026 is NOT a condition
number** ; **37 anomalous nulls of 156**. **Fluid coupling EXCLUDED (turgor LOAD was all of it); the OPERATOR is the only cause left and SURVIVES a feasible shell** (subdiv 6, capacity 0.33→16.6). ⚠ Every resting number used subdiv 3, `preload_…_feasible` **False**, NOT the cause. `STATE_LOG.md`

⚠ **MEASUREMENT RULE — compare γ against the SEED SCATTER, not the within-run `sem`.** Three replicates give
**S = 0.0668 pN/µm = 1.95% of γ, 67× the sem**, so **one seed cannot resolve a γ difference below ~5%** — a
standing limit on every sweep, stage 4 included. Every σ quoted before this is inflated.
