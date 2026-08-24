# GATE-B `nmii_sf_motor` — SF active tension as an emergent myosin-event outcome (2026-07-25)

**Status: MECHANISM native-PASS · QUANTITATIVE solver-BLOCKED.** The declared-but-unimplemented
`nmii_sf_motor` MOTOR edge now has a real runtime; the tension MAGNITUDES it produces are not yet
quotable because the inner solve does not converge (see §4). Commit `fd482920` (+ this closeout).

> **2026-07-28 UPDATE — §7 below.** The lane no longer accepts steps on `accepted_d = ones`: the device
> force-balance flag is the acceptance predicate, its tolerance is derived on-device (D8), a positive
> control shows the gate *rejecting*, and the energy ledger closes with its residual identified as the
> integrator rather than a leak. Run on the same 904-node slice, build `49d2215a`. **Rung stays
> `CUDA_UNIT`**: this is a slice, and the standing rule since the 2026-07-28 retraction is *develop on a
> slice, conclude at native.* §§1–5 are unchanged.

## 1. What landed

`nmii_sf_motor` (MOTOR, `nmii` ⟷ `sf_arc`, bidirectional + adjoint, kinetic, commit-on-accept) has been a
declared edge of `reference_cell_architecture` since the architecture landed, with nothing behind it —
`sf_mechanics.SF_CONNECTOR_BINDING_STATUS` recorded it as *"left honestly SEAMED … needs the NMII actuator +
a live SF bind-target port"*. Both pieces now exist:

| Piece | File | What it is |
|---|---|---|
| straddle-placed NMII population | `ac/engine/sf_nmii_population.py` | host-NumPy builder placing explicit head-resolved bipolar minifilaments on the SF sarcomeres; reuses the audited `ac/motor` primitives (`straddle_frame`, `placed_positions`) verbatim |
| SF motor port + slice | `ac/engine/sf_motor_slice.py` | the `sf_arc` bind-target port over SF-owned arrays, the `nmii` owner over its OWN particle array, the connector, and ONE `CellTransaction` with THREE participants |
| target-agnostic connector | `ac/engine/cortex_motor_slice.py` | `CortexMotorConnector` → `FilamentMotorConnector` (`name`/`component_b` parameterised; cortex defaults byte-identical; alias kept) |
| sarcomere straddle geometry | `ac/engine/sf_population.py` | flag-gated, legacy byte-identical; both lengths DERIVED from the minifilament that must sit inside |
| per-segment metadata | `ac/engine/sf_mechanics.py` | identity / polarity / arclength emitted by the SAME loop that emits the links, so the port cannot drift out of segment order |

**Split ownership is physical here** — stricter than the cortex lane, where the actuator and the port alias one
global `cell.pos_d`. The minifilament particles and the SF nodes are two separate device arrays, so the
crossbridge is a genuine two-array adjoint transfer (`+f` on the head, `−(1−t)f` / `−t·f` on the two SF nodes).

**`sf_arc` is a real participant**, not a port-only target: it owns its arrays and launches its own passive
rod-cable mechanics, so this slice has three participants (SF owner + NMII owner + connector) each
snapshotting / rolling back / committing under one device predicate.

## 2. Geometry — why the heads reach the actin (measured, host)

A rigid bipolar minifilament can only straddle an anti-parallel PAIR. The sarcomere is built as that pair with
both lengths DERIVED from the motor: `lateral = 2·head_offset` (so the `±head_offset` heads land ON the two
filament lines and the arm rests force-free) and `overlap = backbone_contour` (so the whole backbone lies inside
the shared span). Measured against the landed point-to-segment oracle:

| | legacy end-to-end | straddle (this work) |
|---|---|---|
| head → filament-line residual | 0.200 µm (median = max) | **0.000 µm** median, 0.055 max |
| capture (0.210 µm) margin | ~5% — binding sits at the cliff edge | **100%** |
| filaments per head group | 1.20 (heads spill onto both) | **1.00** (clean 1:1) |
| bipolar dot, straight classes | −0.909 (population mean) | **−1.0000** (exact) |
| two filaments at `mid` | co-located (interpenetrating) | side by side |

**Correction to my first diagnosis:** I initially claimed the legacy placement put heads "off the actin". Measured
per SEGMENT rather than per NODE that is wrong — legacy heads are inside capture. The real gain is robustness and
physicality, as tabulated.

**The curved perinuclear cap is EXCLUDED from motor stations by construction.** Its `mid` is an arch apex, so its
two filaments are not collinear over the shared span; a minifilament placed there gets its backbone on the angle
bisector and neither head group walks along the actin it holds (bipolar dot **−0.600**, scale-invariant — it does
not improve with discretisation or nuclear radius). The exclusion is reported by class
(`motor_station_census`), never silently approximated. A curved-sarcomere motor needs its own placement geometry
— **PI decision, surfaced, not approximated.**

## 3. Native gate — MECHANISM PASS (gbook A5000)

`scripts/ac_gate_b_sf_motor_native.py`, record in `sf_motor_native_gate.json` / `.log`.
16 minifilaments / 320 heads on 16 straight stations (4 cap bundles excluded), 360 SF nodes, 320 links.

| Claim | Result |
|---|---|
| heads UNBOUND at t0 (no prestress seed) | ✅ 0/320 |
| rest is force-free (no lumped `k_SF`) | ✅ SF max\|F\| 5.7e-14 pN, T 5.6e-14 pN |
| heads bind by `k_on` Poisson EVENTS | ✅ 0 → 317/320 |
| tension is ZERO while no head is bound | ✅ (nothing else in the slice can make tension) |
| tension EXCEEDS any single head's load | ✅ (load is transmitted + accumulated along the fiber) |
| FA anchor reaction is INWARD | ✅ 100% of anchors (contractile, not sign-flipped) |
| clock advances once per accepted step | ✅ |
| REJECTED step changes nothing | ✅ binding SoA **and** both position arrays bit-restored |

## 4. QUANTITATIVE claim — BLOCKED by the inner solve (do not quote magnitudes)

The explicit overdamped relax does **not** converge: the free-node residual is **15.5% of the reported tension**.
Decisive control — a **10× longer relax moved `T_max` by 3.4×** (0.451 → 0.133 pN) at an *identical* bound
population, so every tension/traction magnitude in the record is a relaxation TRANSIENT.

Root cause is structural, not a tuning matter: the SOURCED α-actinin dorsal↔arc crosslink (4.6e5 pN/µm, several
meeting at one arc apex) dominates the Gershgorin bound, so λ_max ≈ 3.7e6 pN/µm and the CFL-stable explicit step
is ~1.4e-7 µm/pN, while the SF axial mode relaxes orders of magnitude slower. **Fix = an implicit/CG inner solve
for this slice, not a longer explicit relax** — the same wall the cortex hit at GATE A, where
`ProjectedAnalyticCG` + the fiber-arclength multigrid resolved it. That machinery is already in-tree.

An earlier run of this gate reported "T 0 → 1.92 pN, VERDICT PASS". **That was an unconverged transient and is
withdrawn**; the gate now reports the residual every step and refuses the quantitative claim below 1%.

## 5. Adversarial verification (independent agent, tasked to refute)

No sign / index / stale-array bug was found (polarity convention verified end-to-end; contraction sense verified;
split crossbridge momentum-conserving; head ordering consistent). Findings acted on: the cap arch (§2), the
missing convergence diagnostic (§4), verdict thresholds that were effectively no-ops plus a
structurally-guaranteed bound-vs-tension correlation (replaced with causal checks), an unvalidated sarcomere
overlap and a 5% lateral window that admitted non-derived values (both now check against the topology with
floating-point slack only), CLI defaults on GAP params while the docstring claimed none (now REQUIRED), and
several overstated docstrings — including my own claim that the lateral-direction choice is a free gauge, which
measurement contradicts (it moves connector endpoints and changes the curved cap's load path).

## 6. Figures

Regenerate with `PYTHONPATH=<repo> python ffn_sim/scripts/ac_sf_arc_vis.py`; browser-verified with
`scripts/browser_check.py` (a grep of the HTML is not verification). Viewer:
`outputs/ac/sf_arc/sf_arc.html` (6 scenes, full native resolution, no downsampling).

| Figure | Caption |
|---|---|
| `../sf_arc/figs/shot_nmii_sf_motor_one_sarcomere_x15.png` | ONE sarcomere at ×15 isotropic magnification (a microscope objective — uniform on all three axes, no axis truncation; quoted lengths are true µm). The two anti-parallel filaments run side by side separated by 2·head_offset, the NMII backbone sits BETWEEN them, the head↔backbone arms reach from the backbone to BOTH filaments, and the barbed-direction stubs point OUTWARD at each end — the geometry that makes the pull contractile. |
| `../sf_arc/figs/shot_nmii_sf_motor_bipolar_straddle.png` | All 13 motors in the population: `+` heads (yellow) and `−` heads (violet) resolve into two groups at each station, and the perinuclear cap is drawn desaturated with its exclusion in the label (curved sarcomere — cannot be straddled). |

The lane has a **second** figure generator, because the accounting shares no data, axes or scene with the
geometry: `python ffn_sim/scripts/ac_sf_motor_energy_vis.py` (matplotlib only, runs on the dev Mac, reads only
the committed run record). Both are listed here so nothing is left un-regenerated.

| Figure | Caption |
|---|---|
| `figs/balance_gate_margin.png` | The acceptance gate's own margin: the measured `\|reaction + traction\|` and the tolerance derived on-device from `γ_n·Σ\|f_i\|`, on ONE log force axis so the margin is the thing you see — 3–4 decades, every step. The positive control (0.5 pN injected one-sided, rejected) is drawn on the same axis, which is what shows the gate has a working range rather than a floor. |
| `figs/energy_ledger.png` | Left: the balance term by term on a linear axis, so `ΔU + dissipated` visibly tracks `W_active` step to step rather than netting out once; nothing is truncated (headroom is added above the bars, not taken out of them). Right: the residual at `μ` and `μ/2` against **both** predicted power laws — it lies on `μ¹` (integrator consistency) and not on the flat leak line, which is the verdict stated as a picture. |

## 7. Physical acceptance + a closed energy ledger (2026-07-28, build `49d2215a`)

Plan §10 step 2 asks one lane for *a converging solve and a physical acceptance predicate*. The predicate
half landed; the convergence half did not, and the record says which is which.

**Slice: 360 `sf_arc` nodes + 544 `nmii` particles = 904 nodes / 2,712 DOF, 320 heads on 16 minifilaments,
20 pinned FA anchors; 20 accepted steps at `dt_phys` 0.01, explicit relax 4,000. A5000.** That is 0.18% of
the native cell, so nothing here is a magnitude claim — see the rung note in the banner.

### 7.1 Acceptance is decided by the device, on physics

`accepted_d = ones` is gone. Each body reduces its OWN force array into one half of
`GlobalCellLedger`'s balance gate — `sf_arc` into the reaction channel, `nmii` into the traction channel —
so the two channels are the two never-merged arrays whose only coupling is the split crossbridge. Their
vector sum must vanish by Newton's third law, and `balance_ok_d` is what `CellTransaction` accepts on.

| Quantity | Measured over 20 steps |
|---|---|
| accumulated contributions (both channels) | 904 |
| float64 summation ratio `γ_n` (D8, from the COUNT alone) | 1.0025e-13 |
| derived tolerance `γ_n·Σ|f_i|` [pN] | 3.5e-12 → 1.03e-11 |
| measured `|reaction + traction|` [pN] | 3.0e-16 → 3.6e-15 |
| steps accepted by the gate | **20 / 20** |

The residual sits **~3–4 orders below** the round-off floor it is judged against, so the adjoint scatter is
exact to far better than the arithmetic could require. **What the gate tests is the ADJOINT WIRING, not the
convergence** — it holds at any configuration and fails when a scatter is one-sided, sign-flipped or double
counted. Stating that is the point: a reader must not mistake 20/20 for a converged solve.

Two things had to be fixed for this to mean anything:

* **D8's device side had no caller.** `derive_balance_tolerance` now runs inside the transaction, between
  the accumulation and the gate — the only point at which this step's own scale exists.
* **The Higham bound's scale was wrong for a motor lane.** `|Σ̂−Σ| ≤ γ_n·Σ|x_i|` takes the sum of
  MAGNITUDES; bracketing it by `|reaction| + |traction|` is fine for channels that do not cancel, but a
  bipolar minifilament is a force *dipole* — its resultant is near zero while its arithmetic traffic is
  not. That bracket would have set the tolerance orders of magnitude below the round-off it exists to
  admit and rejected every physically correct step. `balance_scale_d` accumulates the real scale; the old
  bracket stays as the fallback, so no existing lane changes.

### 7.2 The positive control — the gate can say no

A gate that accepted every step is indistinguishable, in an artifact, from a gate that cannot fail. So one
head's stall force is injected one-sidedly onto a single `nmii` particle (the run's own `f_stall`, not a
chosen number), through the same predicate path:

| Check | Result |
|---|---|
| `balance_ok` | **0 — rejected** |
| `|reaction + traction|` | **0.49999999999999883 pN** against an injected 0.5 pN |
| against tolerance | 1.06e-11 pN |
| SF + NMII positions bit-restored | ✓ ✓ |
| binding SoA restored, clock did not advance | ✓ ✓ |

The residual reproduces the injection to 14 digits, so the gate is not merely firing — it is measuring the
right quantity.

### 7.3 The energy ledger closes, with NO energy function

No potential is implemented anywhere. Writing one per force family would be a second implementation of the
same physics, free to drift from the forces it claims to differentiate — the "gate a broken build cannot
fail" pattern. Every term comes from launches the step already performs:

| Term | How it is formed |
|---|---|
| `ΔU` | `−∫F_passive·dx` along the step's own trajectory — legitimate because that field was *measured* conservative, and cross-checked against the straight chord |
| `W_active` | the same walk with the force evaluated twice per midpoint, heads bound and heads detached, differenced — the control that attributed the non-conservative content in the first place |
| `dissipated` | `Σ|Δx|²/μ`, accumulated ON DEVICE every inner iteration (a line integral survives subsampling; a sum of squares does not) |
| `event_jump` | **0 by derivation** — the only binding-dependent force here is the crossbridge, which *is* the active channel, so an accepted commit at fixed configuration moves no passive potential |

| step | `ΔU` | `dissipated` | `W_active` | residual | closure ratio | ΔU path-independence |
|---|---|---|---|---|---|---|
| 5 | +3.742e-03 | 2.855e-02 | +3.228e-02 | +1.406e-05 | 4.36e-04 | 2.9e-10 |
| 10 | +5.364e-03 | 3.158e-02 | +3.693e-02 | +1.557e-05 | 4.22e-04 | 9.0e-10 |
| 15 | +5.860e-03 | 3.108e-02 | +3.693e-02 | +1.427e-05 | 3.87e-04 | 6.3e-10 |
| 20 | +5.929e-03 | 2.970e-02 | +3.562e-02 | +1.254e-05 | 3.52e-04 | 2.1e-09 |

All 16-segment paths; every state restored bit-exactly (verified, not assumed). `ΔU` is path-independent to
~1e-9, which is what licenses calling it a potential difference at all.

### 7.4 The residual is the integrator, and it says so itself

The balance does not close to zero, and the reason is structural rather than an accounting failure: an
explicit step takes `Δx = μ·F(x_old)`, so `Σ|Δx|²/μ` is a LEFT-endpoint quadrature of the same integral the
work terms take by the midpoint rule. The two differ by `O(μ)`. A genuine leak — an omitted channel, a
double count, a one-sided adjoint — is a work integral in its own right and does not move with `μ`. So the
residual identifies itself, and no tolerance is chosen anywhere:

| | residual [pN·µm] | closure ratio |
|---|---|---|
| at `μ` | +1.2508e-05 | 3.50e-04 |
| at `μ/2`, 2× iterations (same relaxation) | +6.2429e-06 | 1.75e-04 |
| **measured exponent** | **1.00251** | (1 = integrator consistency, 0 = leak) |

**Verdict `integrator_consistency`** — first order to 0.25%. Both passes are pure observations from the same
configuration with everything restored bit-exactly between them; no event fires and nothing commits.

### 7.5 What this run does NOT establish

* **Not convergence.** The free-node residual is **22.9% of the reported tension**, so the D7 verdict is
  **`VOID`** against the 1% ceiling declared before the run: the tension is a relaxation transient and no
  magnitude may be quoted. Unchanged from §4 — the balance gate does not address it, and does not pretend to.
* **Not rung `CONNECTED`.** That needs this same machinery at the native inventory. `--native-population` is
  an explicit assertion rather than a node-count threshold, because "native" is a claim about where the
  numbers came from, which no count can check.
* One instrumentation defect was found and fixed by this run, and is worth recording because it is a class:
  the energy observer perturbs the position AND force arrays, and restoring only the positions left the
  step's residual being read off the observer's last force evaluation. Recomputing instead of restoring is
  also wrong — an `accumulate` after the step sees the POST-commit binding, while every unmeasured step
  reports the PRE-commit force its solve left. Measured, either mistake inflated the residual on exactly the
  measured steps (23% → 99% → 214%) and flipped this gate's verdict. **An observer must restore the
  observable, bit-exactly, not just the state.**

## 8. Open / next

1. **PI-GATED — implicit/CG inner solve for the SF motor slice** (unblocks every quantitative SF claim). This is
   the critical path; the mechanism is done.
2. **PI decision — curved-sarcomere motor placement** (the perinuclear cap currently carries no motor).
3. **PI-GAP magnitudes** — `f_stall`, `k_on`, catch-slip constants, `k_axial`, NMII stiffnesses, minifilament
   layout (cards N1–N9, S1). Every one is passed explicitly and reported as PROVISIONAL/PROXY; none is defaulted.
4. **Composition** — folding `nmii_sf_motor` into `build_native_composed_cell_world` alongside
   `nmii_cortex_motor` additionally needs head-exclusivity between the two MOTOR edges (one head may bind only
   one target), a separate composition step.

## 9. The inner solve was run, and it retracts this lane's tension (2026-07-28 evening, builds `80473c65` / `a9861bed`)

§8 item 1 named the implicit/CG inner solve "the critical path; the mechanism is done". It has now been
**executed for the first time** — `sf_implicit.py` landed structurally verified but had never run on a GPU —
and the result inverts what the item expected. The solve works. What it converges *to* is the finding.

Nine runs, all on the same 904-node slice (0.18% of native), all `D7 = VOID`, all constants PI-GAP. Records +
figures: `inner_solve_convergence_2026-07-28/`. The lane's own files are byte-identical across the two build
stamps (the intervening commits belong to other lanes), so the nine are one comparable series.

### 9.1 `--implicit-step-scale 1.0` is not the implicit solve

The first run was *worse* than the explicit control (residual/signal 106.7 vs 22.9). The reason is not the
solver: at scale 1.0 the IMEX mobility equals the explicit CFL step, so 4 Newton iterations advance as far as
4 explicit iterations — against the explicit path's 4,000. The IMEX operator `(γ/dt)I + K` is unconditionally
stable precisely so the step may be raised, and the module's own docstring says the equilibrium must then be
invariant to it. **The scale is the whole point of the path, and its default makes the path a no-op.**

### 9.2 Raised, the solve converges six decades — and the tension follows it down

| inner solve | `res` sf_arc [pN] | `res` nmii [pN] | `T_max` [pN] | max head load [pN] | traction [pN] | bound | s/step |
|---|---|---|---|---|---|---|---|
| EXPLICIT, relax 4000 (the prior record) | 3.42e-01 | 1.63e-01 | 1.50e+00 | 2.03e-01 | 3.34e-01 | 316 | 1.24 |
| IMEX 1e2×, 8 Newton | 1.04e+00 | 3.64e-01 | 1.97e+00 | 3.81e-01 | 1.25e+00 | 316 | 1.17 |
| IMEX 1e4×, 8 Newton | 2.83e-02 | 2.20e-01 | 4.55e-01 | 8.22e-02 | 9.96e-01 | 316 | 1.21 |
| IMEX 1e5×, 8 Newton | 2.57e-03 | 5.45e-01 | 5.09e-02 | 1.22e-02 | 1.51e-01 | 316 | 1.17 |
| IMEX 1e6×, 8 Newton | 1.31e-04 | 5.10e-01 | 4.74e-04 | 1.17e-03 | 8.73e-04 | 316 | 1.27 |
| IMEX 1e6×, 32 Newton, 2000 CG | 5.89e-05 | 2.47e-01 | 1.21e-04 | 2.57e-04 | 1.00e-04 | 316 | 24.42 |
| **IMEX 1e8×, 32 Newton, 2000 CG** | **1.89e-07** | 7.84e-02 | **3.75e-06** | 4.57e-04 | 1.48e-06 | 316 | 24.93 |

**The bound-head count never moves.** 316 of 320 heads are bound in every row, so nothing here is a binding
or a kinetics artifact. What moves is the mechanics: as the inner solve is driven six decades closer to
equilibrium, the SF axial tension and the FA traction fall *with* it, essentially 1:1 across seven decades
(`figs/tension_is_the_residual.png`, which draws the 1:1 line as an exact reference, not a fit).

So the previously reported SF tension was **not** a relaxation transient in the innocent sense — a transient
is a real force on its way somewhere. It was the **unconverged residual itself**, read through the tension
diagnostic. `T = k_axial·(L − r0)` is a perfectly good estimator; it was estimating an un-relaxed geometry.

### 9.3 Three objections, measured rather than argued

* **"It would accumulate if you ran longer."** The same converged configuration for **10× the physical
  time** (200 steps, 2.0 s): `T_max` reaches 4.12e-03 pN with 314 heads bound — still **two decades below
  ONE head's `f_stall`** (0.5 pN), and the residual rises alongside it. `figs/duration_control.png`.
* **"The minifilament's 19 internal zero-stiffness modes let the crossbridge relax for free."** The obvious
  mechanism, given the tier-(a) floppiness row, and it is **REFUTED**: stiffening the head arm ×100 and the
  backbone ×100 at the most converged configuration leaves `T_max` at 3.80e-06 pN against the baseline's
  3.75e-06 — no change. Recorded as a diagnostic, not a tuning: no gate was scored against it.
* **"The fiber is anchored on one side."** It is not. `VENTRAL` bundles are built `fa_both` — both outer
  sarcomere ends carry an FA anchor node — and the implicit path puts the Dirichlet mask in the OPERATOR, so
  `dx` at those 20 nodes is exactly zero. A both-ends-pinned dipole with 316 bound heads carries 3.75e-06 pN.

### 9.4 The D7 ratio ceiling cannot certify this lane, and that is a contract question

Residual and signal vanish **together**, so `residual/signal` never falls under the 1% ceiling — the best it
reaches is 1.21%, while the residual underneath it falls six decades. A ratio ceiling presumes a signal that
survives convergence; here the signal *is* the residual, up to a factor of a few. **Per the no-gate-loosening
rule this is surfaced, not edited**: the fix is a denominator that does not vanish (the active channel's own
Σ|f|, or `f_stall × n_bound`), and that is a gate-contract change requiring PI sign-off.

### 9.5 What this DOES deliver for stage 2 — CORRECTED the same day

**The first version of this section was wrong and is retracted.** It read the series as "1.24 s/step at
residual 3.4e-01 → 24.9 s/step at 1.9e-07, ≈20× for six decades" and offered that as the convergence-cost
multiplier stage 2 needs. It is **confounded**: the cheap and expensive runs differ in the solver BUDGET the
operator was given, not only in the residual they reached. Read off the nine `config.argv`, wall clock tracks
`cg_iterations × newton` and is almost independent of residual — and every one of the nine records had already
stamped its own timing `comparable: false`, because every gate is `VOID`. The row quoted a ratio the artifacts
themselves refuse to compare. Caught by Lead review on the day it landed; retired to `STATE.md` (c) 15.

| budget = `cg × newton` | s/step | residual reached [pN] |
|---|---|---|
| 800 (`cg` 200 × `newton` 4) | 0.313 | 4.57e+00 |
| 3,200 (400 × 8) | 1.170 – 1.267 | 1.04e+00 … 1.31e-04 (five decades, at ONE budget) |
| 64,000 (2000 × 32) | 24.4 – 25.1 | 5.89e-05 … 1.89e-07 |

**Read correctly the same data says the opposite, and that is the better result.** At a FIXED budget of 3,200,
raising `--implicit-step-scale` 1e2 → 1e6 moved the residual **1.035e+00 → 5.569e-04, a factor 1,859, for
+0.2% wall clock** (1.171 → 1.174 s/step). Spending 20× more iterations at a fixed scale bought only 9.5×
(5.569e-04 → 5.891e-05) for 20.8× the time. **The cheap lever is the IMEX step scale; the iteration budget is
the expensive one.** That is a statement about which knob to turn, not a cost, so it survives the confound —
but it is still a slice, and conditioning worsens with density (2026-07-22e: 8k 2× → native 8×), so no cost
law is quotable from 904 nodes with 4 α-actinin joints. The native cost measurement stands undone.

### 9.6 What is now blocked, and on whom

1. **PI — the physics.** A both-ends-anchored sarcomere with 316 bound heads holds no load at equilibrium.
   The two candidate mechanisms left standing are (i) the walked `abscissa` being absorbed by minifilament
   **rigid-body translation** (its only tie to the world is its own crossbridges), and (ii) the contractility
   asymmetry Belmonte 2017 names — a minifilament on antiparallel filaments is predicted *neutral* without
   passive crosslinkers, and this slice has 4 α-actinin dorsal↔arc joints against native's 1,413,720.
   **(ii) is a CONNECTIVITY property, i.e. exactly the thing a slice cannot answer.**
2. **PI — the gate contract** (§9.4).
3. **Task 1 of this lane (cost of a converged native run) is not blocked by the solver any more.** It is
   blocked by there being no converged observable worth costing until 1 is answered.

### 9.7 Figures

Regenerate with `python ffn_sim/scripts/ac_gate_b_sf_motor_convergence_viz.py` (matplotlib only, runs on the
dev Mac, reads only the committed run records — an absent record is reported and skipped, never plotted as
zero). This is the lane's **third** figure generator; it shares no data or axes with the other two.

| Figure | Caption |
|---|---|
| `inner_solve_convergence_2026-07-28/figs/tension_is_the_residual.png` | SF tension and FA traction against the free-node residual they were measured with, log-log over seven decades, both panels on identical square limits with an **exact** 1:1 reference (not a regression) and one bound head at `f_stall` drawn as the model's own force scale. The explicit control and the NMII-stiffening diagnostic are on the same axes as distinct markers. |
| `inner_solve_convergence_2026-07-28/figs/per_step_traces.png` | Every run's per-step trace, drawn per realisation rather than as a summary scalar: each plateaus within a few steps at a level set by its OWN residual, so the near-zero tension is a steady state and not a ramp caught early. |
| `inner_solve_convergence_2026-07-28/figs/two_body_residuals.png` | Both owners' residuals against the IMEX step scale. The honest limit of the convergence claim: `sf_arc` falls six decades, `nmii` does not — split ownership means the gate must see both, and the nmii body is the open end. |
| `inner_solve_convergence_2026-07-28/figs/cost_vs_convergence.png` | Wall clock per accepted step against the residual that step reached, x-axis inverted so "more converged" is to the right. A cost is a property of an engine at a stated convergence; the two axes therefore belong on one plot, and the caption says every point comes from a VOID run. |
| `inner_solve_convergence_2026-07-28/figs/duration_control.png` | The duration control: the same converged configuration for 10× the physical time, with tension, residual and one head's `f_stall` on ONE log force axis so the remaining distance is the thing you see. |

## 10. NATIVE, and §9's conclusion is RETRACTED: the tension emerges — the runs were too short (2026-07-28 night, builds `576ebf42` / `2001c119`)

§9 concluded from nine runs that the converged equilibrium of this lane carries no tension. **Every one of
those runs was 20 steps = 0.20 s of physical time, and 0.20 s is inside an induction period.** Run natively at
five seconds on a topology that actually converges, the tension emerges properly and the gate PASSES. Seven
further native runs, records in the same directory.

### 10.1 The α-actinin crosslink is what blocked convergence — ten decades, measured A/B

Identical 10,400 nodes, 4,000 heads, 200 motor stations, identical solver budget; the ONLY difference is the
dorsal↔arc crosslink (4.6e5 pN/µm), at the same 0.20 s:

| topology | α-actinin joints | free-node residual [pN] | D7 |
|---|---|---|---|
| `topo_ventral_only` — pure closed ventral dipoles | **0** | **1.3334e-12** | FAIL (interpretable) |
| `topo_joint_dense` — dense joints | 160 | 1.4568e-02 | VOID |

**Ten decades.** §9 hunted for the absorber of the myosin work in the wrong place: the problem was never that
the model is neutral, it is that one stiffness family makes the inner solve unable to reach equilibrium, and
an unconverged residual then reads out through the tension diagnostic. `sf_implicit.py` says as much in its
own docstring — the α-actinin joint dominates λ_max — but it was written as a reason the EXPLICIT path fails,
and the same term turns out to defeat the implicit path at every budget tried here too.

### 10.2 The run that PASSES — the first in this lane

`topo_v_long500`: 200 ventral sarcomeres, 10,400 nodes, 4,000 heads, 500 accepted steps = **5.0 s** physical.

| t [s] | bound heads | `T_max` [pN] | max head load [pN] | free-node residual [pN] |
|---|---|---|---|---|
| 0.20 | 3,928 | **3.74e-06** | 1.19e-04 | 1.33e-12 |
| 0.50 | 3,921 | 0.947 | 0.627 | 6.21e-10 |
| 1.00 | 3,925 | 1.829 | 0.863 | 1.63e-11 |
| 2.00 | 3,933 | 2.991 | 1.072 | 8.96e-09 |
| 5.00 | 3,933 | **5.001** | 1.773 | **1.18e-07** |

**D7 verdict = PASS**, `residual/signal = 2.36e-08` against the 1% ceiling — eight decades of margin, and the
first non-`VOID` verdict this lane has produced. All twelve mechanism predicates true, including
`tension_exceeds_single_head_load`. The timing block is stamped **`comparable: true`**.

The 0.20 s row is the whole lesson: it is the converged value, it is ~0, and it is where §9 stopped.

### 10.3 What §9 got right, and what is retracted

* **RIGHT, with its scope now stated**: at t = 0.20 s the converged tension really is ~1e-6 pN, so the 1.50 pN
  the explicit relax reported at that instant WAS the residual. The 1:1 tracking across seven decades stands —
  **as a statement about runs sampled at 0.20 s**, which is why none of those magnitudes meant anything.
* **RETRACTED**: "the converged equilibrium has no tension", "a both-ends-anchored sarcomere with 316 bound
  heads holds no load", and the whole §9.6 item 1 framing that made this a Belmonte 2017 connectivity question.
  It is not a connectivity question. Removing crosslinks **entirely** is what let the tension appear.
* **Also retracted**: §9.3's minifilament and anchoring "refutations" were argued from 0.20 s runs. Re-run on
  the converged topology they are unchanged at 0.20 s (`topo_v_stiffnmii` 3.74e-06 → 5.85e-06 pN; ×100 axial →
  2.22e-05), but that now says nothing at all, because the correct answer at 0.20 s is ~0 either way.

### 10.4 Stage-2 cost — the first `comparable: true` measurement in this lane

`wall_seconds_per_step = 1.122 s` · `n_inner_iterations = 4,008` · `real_time_factor = 112` · residual
1.18e-07 pN · D7 PASS. A cost is a property of an engine **at a stated convergence**, and this is the first
time this lane has had both halves in one record. It is **not** the full-native cell (that is 511,114 nodes;
this is 10,400) and the native SF inventory is still a KB gap — see §10.5.

Population scaling, measured on the way (same budget, same 0.20 s): 904 → 9,040 nodes changed `T_max` by
**+0.03%** and the residual by +0.8% while the FA traction resultant scaled ×10.2 — intensive and extensive
quantities separating exactly as they must. At 90,400 nodes the same budget no longer converges (residual
4.5e-02), which is the conditioning-with-density effect, not a physics change.

### 10.5 Open, and what changes for the PI items in §9.6

1. **§9.6 item 1 is WITHDRAWN as posed.** The replacement question is narrower and sharper: *why does the
   α-actinin dorsal↔arc joint cost ten decades of convergence, and what preconditioner fixes it?* The cortex
   lane solved its own version of this with a fiber-arclength multigrid; `sf_implicit.py` deliberately does not
   use it, and its stated reason (the V-cycle presumes no intra-fiber axial operator) should now be re-examined.
2. **§9.6 item 2 stands, and is now sharper.** D7's ratio ceiling worked exactly as designed here — it called
   `VOID` on every run whose signal was residual and `PASS` on the one that was not. The defect is only that a
   ratio cannot certify a run whose signal is genuinely near zero, which is the 0.20 s case.
3. **NEW — no run in this lane may be 20 steps again.** The induction period is ~0.3 s at `dt = 0.01`; a
   duration below it cannot distinguish "no tension" from "not yet". This is a protocol item for stage 0.
4. The `topo_*` runs are a deliberately **reduced architecture** (ventral only). They establish the mechanism
   and the cost at native scale; they are not the cell, and the SF inventory remains a PI-GAP.

### 10.6 Figures

Same generator (`ac_gate_b_sf_motor_convergence_viz.py`), two more:

| Figure | Caption |
|---|---|
| `inner_solve_convergence_2026-07-28/figs/tension_emerges_when_it_converges.png` | The PASS run: tension and residual over 5.0 s on ONE log force axis, with `f_stall` for scale and a vertical marker at **0.20 s** showing where every earlier run stopped. The title states the late-time slope in pN/s because a log axis flatters a linear rise. |
| `inner_solve_convergence_2026-07-28/figs/crosslink_blocks_convergence.png` | The A/B at identical size, head count and budget: residual and reported tension for 0 vs 160 α-actinin joints, log axis, with `f_stall` overlaid. |
