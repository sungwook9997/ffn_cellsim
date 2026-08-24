# C-1 — the one experiment, declared in advance

**Status: PI-SET 2026-08-10** (session `47e953f4-93d6-4c67-a901-08ca181222f6`, 00:35 KST).
**Branch: `codex/c1-probe-scale-relaxation`**, worktree `/Users/sw1/ffn_cellsim_c1_relaxation`, cut from
`22e404f6`. Held separately from `codex/cortical-tension-afm`, which a Codex session is running.

> C-1 is *"the **one experiment** the sweep targets is declared in advance — name, protocol, observable,
> uncertainty"* (`docs/decisions/PROPOSAL-the-merge-and-the-rename.md` §C-1). This file is that declaration.
> It is a **gate contract**: it is written BEFORE the run and may not be edited inline or re-thresholded
> after seeing data. Its final home is `docs/decisions/`, which belongs to the `contracts` lane; it is
> written here (`v2_audit/`, `free`) so that landing it does not require a lane fight, and it should be
> moved, not rewritten, when `contracts` ratifies.

---

## 1. Name

**PSR — Probe-Scale Relaxation.**

*Which clock owns the cell's measured stress relaxation, and can any single-probe method tell?*

## 2. Why this experiment and not cortical tension

Cortical tension was considered first and **rejected as circular**. In this model γ = ΔP·R/2 with ΔP set
by the osmotic turgor Π₀, and `aleph/laws/turgor_pi0.py` records that the tree holds three competing Π₀
values, **none of them a measurement of the target cell line**, one of which (133 Pa) is the cortical-tension
band inverted through Young-Laplace — *"circular by construction"*, in the module's own words. Measuring γ
with Π₀ set that way scores the model against an answer its own machinery produced, which the charter
forbids outright.

PSR does not have that structure. Its inputs are the **solvent** viscosity (`ETA_SOLVENT = 1.0e-3`
pN·s·µm⁻², water at 37 °C) and explicit mesh geometry; the composite cytoplasm viscosity is never an input
on this path. `aleph/laws/units.py:60` states the design intent as a prediction: *"the measured effective
cytoplasm η=65.9 Pa·s IS the poroelastic composite (fiber friction + Darcy resistance of solvent
percolating the ~ξ mesh); with the Biot fluid explicit, the bare fibers feel ONLY the solvent and the 65.9
EMERGES from the coupling."* That is falsifiable and non-circular.

## 3. What the engine contains that makes this the right target

Four relaxation clocks are already implemented, spanning three orders of magnitude, and an experimentalist
sees all four lumped into one word — "viscoelasticity". These are **code parameters, not results**; the
`laws/` result corpus is invalidated-pending-rerun (`STATE.md` (a)).

| Clock | Where | Scale | Scales with probe radius `a`? |
|---|---|---|---|
| poroelastic `τ_p ~ L²/D` | `network_warp.py:708`, `D_um2_s` def 50 (Moeendarbary) | ~1 s | **YES — `∝ a²`** |
| crosslink turnover `1/k_off` | `xl_koff_per_s`; α-actinin `k_off0` 0.066 s⁻¹ | ~15 s | no |
| membrane drainage `τ_osm = (V₀−v_min)/(L_p·A·Π_in)` | `network_warp.py:840` | ~30–200 s | no |
| viscous press transient | `v_press_um_s`, `eta_bulk_Pa_s` | set by rate | no |

**Exactly one of the four carries a length scale.** That is the discriminator, and it is why the experiment
is a *sweep in `a`* rather than a measurement at one `a`.

The engine can also do what no experiment can: switch clocks off individually. Each switch already has a
regression gate asserting that OFF is bit-identical (`tests/ff/test_poroelastic_cytoplasm.py`:
`test_turnover_off_is_backward_compatible`, `test_press_ramp_off_is_backward_compatible`,
`test_backward_compat_all_off_is_pure_turgor`).

## 4. Protocol

**Cell state.** Rounded / non-adherent. Held by the exterior Stokes medium
(`aleph/engine/medium_exterior.py`), whose six rigid modes are a *physical* holder — not the numerical `aI`
regulariser it replaced. **Not** the adherent, spread state: there SF and FA dominate the response and the
measurement is no longer about the cytoplasm.

**Load.** Step indentation through the unilateral contact law
(`aleph/engine/connector_joints.py::_unilateral_contact_force_kernel` — force-free at and beyond `rest`,
repulsive inside it, sign gated on the gap so a contact can never pull), then **hold**. Record the reaction
on the indenter, the quantity an AFM cantilever reads.

**Sweep.** Contact radius `a` over at least four values spanning ≥ one decade, at **fixed δ/a** so the
strain field is self-similar and `a` is the only thing that changes. Dwell ≫ the longest clock under test.

**Rate.** A physical press speed (`v_press_um_s`) with η-calibrated drag, never an instant strain at a
numerical `dt`. Reading force before the cell has finished deforming is what produced the 857× artifact on
2026-07-08; that failure is the reason this line exists in the protocol.

**Controls, all in silico, all pre-declared:** (i) `Lp` off, (ii) `xl_koff_per_s` off, (iii) `biot_fsi` off.
Each removes exactly one clock. A control that changes nothing is as informative as one that changes τ.

## 5. Observable — the scalar, and the projection, declared here

For each `a`: fit `F(t)/F(0)` to the SLS form `G(t) = G_inf + G_1·exp(−t/τ)` using the **existing oracle**
`aleph/validation/oracles/mechanics/viscoelastic_relaxation.relaxation_time_from_curve`. The oracle is a
cross-check and scorer; it is never imported by a runtime path.

**The scored scalar is the exponent `p` in `τ(a) ∝ a^p`**, from a log-log fit across the sweep.

    poroelastic-dominated   ->  p = 2
    turnover / drainage / viscous dominated  ->  p = 0

An exponent, not a magnitude, and deliberately so: this tree's tier-(a) rows that survived scrutiny are the
ones that quote scaling and refuse magnitude (`STATE.md` (b), row 1: *"Quote the attribution and the scaling
exponents; no magnitude"*), while every magnitude claim about tension is currently in (c).

**Secondary, and NOT the gate:** the effective composite viscosity emerging from `D` + `M_biot`, compared
against Dessard 2024 MCF7 65.9 Pa·s. A magnitude, so it is reported as a cross-check with its scatter and
may not close C-1 on its own.

## 6. Uncertainty — declared before the run

1. **Seed scatter is measured first, and it is the yardstick.** The standing MEASUREMENT RULE (`STATE.md`
   (f)) is that within-run `sem` understates by ~67× and one seed cannot resolve a 5% difference. The
   analogous scatter for τ is unknown, so the first run of this line is **three replicates at one `a`** to
   measure it. No `p` is quoted before that number exists.
2. `p` is scored against a bootstrap CI over seeds, and the verdict is the *separation*: whether the CI
   excludes 0, excludes 2, or excludes neither. **"Excludes neither" is a real and publishable outcome** —
   it is the statement that the observation cannot identify the mechanism, which is the product this
   project says it is building.
3. **A fit quality floor is declared here**: an `a` whose SLS fit residual exceeds the pre-declared bound is
   reported as unfittable and excluded from the log-log fit *with the exclusion recorded*, never dropped.

## 7. What C-1 may not claim

- Not a rung. Not `CONNECTED`.
- Not a cytoplasm viscosity magnitude, unless §6.1 has been done and the scatter supports it.
- Not a statement about adherent cells — the protocol is explicitly the rounded state.
- Not a transfer claim: one archetype is not a universal law until shown to transfer.

## 8. Known blockers, in the order they bind

1. **C-2.** The resting step does not converge — 66× the inner budget ends *higher* than it began
   (`STATE.md` (b)), and a coupled 3-step run at full native on 2026-08-10 returned `maxPF = 0.7902`
   identical to four figures on every step with zero membrane flux. A relaxation curve read off a
   configuration that never reaches equilibrium is not a material response. **C-1 is declared now, but it
   cannot be measured before C-2.** Declaring it now is the point: C-2 needs a tolerance, and a tolerance
   needs the experiment named.
2. **Wiring.** `biot_fsi` lives in `aleph/laws/network_warp.py`, the legacy whole-cell compression driver.
   It has to become an `engine/` component + connector before it can run inside an accepted-step
   transaction.
3. **PI-GAP parameters.** `D_um2_s = 50` and `M_biot_Pa = 300` — provenance unverified for the target cell
   line. The contact law's `rest_um` and stiffness are declared GAP by `connector_joints.py` itself, which
   states that no rate constant, stiffness or contact distance is defined in that module.
4. **Indenter provenance.** Tip radius and cantilever stiffness come from whichever protocol is being
   reproduced; they are not ours to choose.

## 8a. STEP 1, FULL DECADE — the definitive run. RTX 4090, build `2753580a`

`outputs/ac/c1_biot_timescale_decade/record.json`, verdict **PASS**. Eleven probe sizes from 0.40 to
4.00 µm — exactly a factor of 10 — across four boxes (L = 24 / 48 / 96 / 192 µm) each with its own dx,
holding `a/dx` between 3.67 and 7.33 so an `a`-dependence cannot be a discretisation-ratio dependence
in disguise.

| sweep | declared | measured | evidence |
|---|---|---|---|
| `tau` vs `a` | +2 | **1.9997781** | r² = 0.99999893, 11 points over a **decade**, four boxes |
| `tau` vs `D` | −1 | **−1.0000000** | r² = 1.00000000, 6 points over a 32× range |
| `tau` vs `dx` | invariant | spread **0.18 %** | 4 spacings over 3× |

Oracle ratio median **0.999896** (min 0.997361, max 1.001227). Peak pore-mass drift **3.2e-11**. Three
box overlaps, each agreeing to **0.387 %** while changing box *and* mesh at once.

**And the residual is a pure function of `a/dx` — to MACHINE PRECISION:**

| `a/dx` | oracle ratio | probe sizes sharing it | spread within the group |
|---|---|---|---|
| 3.67 | 1.001227 | 0.55, 1.10, 2.20 | **8.1e-13** |
| 5.00 | 0.999896 | 0.75, 1.50, 3.00 | **2.7e-12** |
| 7.33 | 0.997361 | 0.55, 1.10, 2.20 | **9.7e-13** |

Three probe sizes spanning a factor of four give *the same* deviation from the closed form to **twelve
decimal places** whenever they share a mesh ratio. The dependence on `a` is not small — **it is absent
at machine precision.** Everything remaining in the residual is discretisation.

That is a much stronger statement than "the exponent is 2 with 0.1 % of something unexplained", and it
is only visible because the design put several probe sizes at each mesh ratio. A decade of `a` at a
single `dx` spans a decade of `a/dx` too, and no fit could then have separated them.

Per-observer over the decade:

| observer | exponent | prefactor [a²/D] | closed form | r² |
|---|---|---|---|---|
| whole-field L2 | 1.999778 | 1.396039 | 1.396834 | 0.99999893 |
| single point | 1.999684 | 0.475535 | 0.473867 | 0.99999816 |
| finite contact | 1.999880 | 0.579423 | — | 0.99999792 |

**Exponent spread 0.000196. Prefactor ratio 2.93572.**

## 8b. STEP 1, first pass — measured 2026-08-10, RTX 4090, build `ee4a26bb`

`outputs/ac/c1_biot_timescale/record.json`, verdict **PASS** on all three sweeps.

| sweep | declared | measured | evidence |
|---|---|---|---|
| `tau` vs `a` | +2 | **1.997519** | r² = 0.999999, 6 points over a 4× range, two boxes |
| `tau` vs `D` | −1 | **−1.000000** | r² = 1.000000, 5 points over a 16× range |
| `tau` vs `dx` | invariant | spread **0.18 %** | 4 spacings over 3×; scored on spread, not r² (§ why in the driver) |

**The poroelastic clock carries a length scale, and it is the one the theory says.** Against the
closed form `tau = a²(e^(4/3)−1)/2D`, the oracle ratio has median **0.99983** (min 0.99782, max
1.00174) — so the prefactor is right, not only the exponent.

**Controls that could have failed and did not.** Box overlap: `a` = 1.0 µm measured in a 48 µm and a
96 µm box agrees to **2.0e-3**, so the domain is not in the answer. Peak pore-mass drift across every
point is **1.1e-10**, four orders under the 1e-6 guard the setup render forced into existence.

**The multi-method result, at the field level.** Three observation operators on the same states:

| observer | exponent | prefactor [a²/D] | closed form |
|---|---|---|---|
| whole-field L2 | 1.99752 | 1.396271 | 1.396834 |
| single point | 1.99470 | 0.475203 | 0.473867 |
| finite contact patch | 1.98663 | 0.581867 | none |

**Exponents agree to 0.011. Prefactors span 2.938×.** The scaling is a property of the medium; the
magnitude is a property of the instrument. Note which observer is furthest from 2 and has the lowest
r²: the finite-contact one — the only one that resembles an instrument, and the only one whose
footprint is discretised. The reading that a laboratory could actually take is the least robust of
the three.

⚠ **THIS IS NOT A TIER-(a) ROW AND MUST NOT BECOME ONE.** (b) admits results at the native full
population; this is the Biot field alone — no mechanics, no contact law, no compartments, no cell. It
is a component-level gate and its record says so in `may_not_be_quoted_for`. It licenses the C-1
*design* — that an exponent is the right scored scalar — and nothing about a cell.

## 8g. All three methods over the decade, and how mesh-sensitive each one is

Once the residual is known to be a function of `a/dx` alone, two things follow that no single-mesh
sweep can give: **how steep that function is, is a property of the method** — a number for "how fine a
mesh does this protocol need", which is the practical question and is essentially never reported — and
the prefactor can be extrapolated toward an infinitely fine mesh.

| method | C spread across `a/dx` 3.67→7.33 | extrapolated continuum C | fit r² | oracle |
|---|---|---|---|---|
| step relaxation | **0.39 %** | 1.38855 | 0.938 | 1.396834 |
| oscillatory | **0.36 %** | 0.86229 | **0.018** | — |
| creep | **9.87 %** | 0.31742 | 0.934 | 1/π = 0.318310 (**−0.28 %**) |

**Creep is 25× more mesh-sensitive than step relaxation for the same measurement.** That is the honest
verdict on it — not "wrong" but "slowly convergent", and the difference matters because a study running
creep at a mesh that is perfectly adequate for step relaxation would be 10 % off and have no signal
saying so.

**And creep does reach its closed form.** Extrapolating C against 1/(`a/dx`) gives 0.31742 against
1/π = 0.318310 — **0.28 %**. The method that failed every grid check converges to the analytic answer;
it just needs a mesh nobody would have thought to use.

⚠ **Three different extrapolation statuses, and they are not interchangeable.** The oscillatory fit has
r² = 0.018: there is no trend, the residual is scatter at the converged level, and its "continuum"
value is a line through noise. The step method has a real trend but the extrapolation lands *outside*
the measured range and *away* from the oracle, which says the assumed order is wrong — a second-order
scheme errs as (dx/a)², and a fit linear in (dx/a) overshoots. Only creep has a trend both strong
enough and clean enough for the extrapolated number to mean anything. The synthesis flags all three
rather than reporting one column of continuum values as if they were equivalent.

## 8f. The synthesis — five readings of one medium

`outputs/ac/c1_synthesis/` (`synthesis.json`, one figure, one HTML page verified in a browser).

| method | reading | exponent p | prefactor C [a²/D] | r² | pts | verdict |
|---|---|---|---|---|---|---|
| step relaxation | whole-field L2 | 1.9975 | 1.3963 | 1.00000 | 6 | PASS |
| step relaxation | single point | 1.9947 | 0.4752 | 1.00000 | 6 | PASS |
| step relaxation | finite contact | 1.9866 | 0.58187 | 0.99998 | 6 | PASS |
| oscillatory | −3 dB bandwidth | 2.0030 | 0.86246 | 1.00000 | 3 | PASS |
| creep | flux decay | 0.7641 | 0.70786 | 0.98263 | 3 | **FAIL** |

**Exponent spread 0.0164** (bound 0.10) across the four valid readings; **prefactor ratio 2.9383**
(minimum 1.5). Two of the five prefactors have closed forms and both are reproduced: 1.3963 against
1.396834, 0.4752 against 0.473867.

**The rule the first synthesis was missing.** A reading whose *own* method reported a failing internal
validity check is evidence about the method, not about the medium. Creep's r² is **0.98263 — above its
own declared floor** — and its grid-invariance check failed by a factor of 21; letting it into "do the
exponents agree?" answers a question about the medium with a measurement that is not of the medium.
It is a partition, not a filter: creep keeps its row, keeps its curve on the log-log panel, and the
headline states that one reading was excluded and why.

## 8c. A limitation of this gate's own design, found by running it — THREE times

**An exponent can be exactly right while the measurement under it is corrupt.** This is now the most
robustly established thing in the whole C-1 line, because it turned up three independent times:

1. Run 34's frequency-domain `D` sweep sat entirely on `a` = 1.1 µm — the one probe size the
   self-similarity guard rejects — and returned exponent **−1.0000000** with r² at the floor's ceiling.
2. Run 38's creep `D` sweep returned exponent **−1.0000** on a method whose own grid-invariance check
   FAILED by a factor of 2.36 (below).
3. Run 38's oscillatory `a` sweep gives exponent **2.00195** with the corrupted point excluded and
   **2.00023** with it *included*. The corrupted point makes the exponent look **better**.

That third one is the sharpest. A sweep contaminated by a bad point did not merely survive the
contamination — it was flattered by it. No amount of staring at exponents and r² values would have
found the defect; only the physics-guaranteed self-similarity invariant did.

Run 34's frequency-domain `D` sweep was taken entirely at `a` = 1.1 µm — the one probe size the
self-similarity guard rejects (§8d). It nonetheless returned `tau` = 0.08464 s at `D` = 12.5 and
0.04232 s at `D` = 25: a ratio of exactly 2, an exponent of −1, and an r² that would have passed the
floor comfortably.

The reason is structural and worth stating plainly: **the corruption is scale-covariant.** Whatever
went wrong at `a` = 1.1 went wrong the same way at every `D`, because `D` enters only through the tone
grid, so it divided out of the ratio. A sweep whose systematic error scales with the sweep variable
returns the declared exponent and tells you nothing.

This does not retract §8b. That result carries something the `D` sweep here does not: agreement with a
**closed-form prefactor** to 0.017 %, on six points, in two boxes, with an independent mass invariant
at 1e-10. A prefactor cannot be right by covariance — it has to be the number.

It does qualify the C-1 headline, and the qualification belongs in the gate: *the exponent is what
transfers across methods, but agreement on an exponent is not by itself evidence that a measurement is
sound.* Where a closed form exists, the gate must check the prefactor too, and where none exists — the
finite-contact observer, the oscillatory method, creep — the exponent must be defended by internal
controls instead. That is why §8d's self-similarity check exists and why it is scored per point.

## 8d. Method 2 — two wrong diagnoses before the right one, and one defect still open

The 45-degree inversion fired at exactly one of four probe sizes in run 34. I gave two explanations
and both were wrong; both are kept here because a retraction that erases what was retracted teaches
nothing.

1. *"45° is the asymptote of a diffusive phase lag."* Refuted by data already in the same run —
   `a` = 1.1 reached 65.5°.
2. *"The inversion is ill-posed on this medium."* Also wrong. Widening the tone band from
   (0.25 … 4.0) to (0.125 … 16.0) × `D/a²` takes the phase to **82.7°** and the crossing exists at
   every probe size.

**The band was mis-centred, and that is the finding.** The dimensional guess `ω_c = D/a²` is ~4× low,
so the original top tone stopped at 43°, just short of the level being looked for. This is precisely
what a real instrument does when it selects its frequency range from a dimensional estimate, and what
it returns is either no measurement or — at one probe size — a corrupted one. Neither looks wrong.

**The corruption is real, reproducible, and still undiagnosed** — but its *consequence* is now the
most important result in this section. It affects one or two points per run, moves between runs, and
always carries the same signature: the phase saturates near 64° and the tail turns over
(64.1 → 60.1°) instead of climbing to ~82°.

Run 52 (six probe sizes, two corrupted), prefactors in `a²/D`:

| a | C(−3 dB) | C(phase-45) | self-similarity dev |
|---|---|---|---|
| 0.60 | 0.8631 | **0.7675** | 21.82° ← corrupted |
| 0.80 | 0.8619 | 0.8735 | 0.49° |
| 1.00 | 0.8625 | 0.8745 | 0.31° |
| 1.30 | 0.8631 | 0.8760 | 0.19° |
| 1.70 | 0.8644 | **0.7668** | 22.65° ← corrupted |
| 2.20 | 0.8645 | 0.8651 | 0.19° |

**The −3 dB inversion cannot see the corruption at all.** Its two corrupted readings, 0.8631 and
0.8644, fall *inside* the clean spread of 0.8619–0.8645. Phase-45 puts the same two points 12 % out.
A study measuring only the −3 dB point would have six values agreeing to 0.3 %, two of them corrupt,
and nothing anywhere in its own data to say so.

So the trade between the two inversions is not precision against transferability. **It is precision
against the ability to detect that you are wrong, and the more precise instrument is the blind one.**

### DIAGNOSED — it was one floating-point comparison

The lock-in projected against `t − t_window_start` while the drive is defined on **absolute** `t`. That
adds a phase offset `ω·t_start` to every tone, which is a multiple of 2π only when `t_start` is an
exact whole number of periods of every tone. The discard boundary is *designed* to land there — sample
4096 of 8192, every tone completing an integer count — but the comparison that finds it,
`t >= N_DISCARD * period_slow`, sits on a knife edge, and one ulp decides whether that sample is in or
out. Out by one shifts the phase reference by exactly one sample.

Predicted error: `−ω·sample_dt = −(2π/256)·m = −1.4062·m` degrees.

| tone multiple m | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| measured (corrupt − clean) | −0.9 | −1.6 | −3.0 | −5.7 | −11.3 | −22.3 |
| predicted | −0.70 | −1.41 | −2.81 | −5.62 | −11.25 | −22.50 |
| ratio | 1.28 | 1.14 | 1.07 | **1.01** | **1.00** | **0.99** |

**The signature was diagnostic for three runs and I did not read it.** Amplitudes were always clean —
a time shift does not touch them, and that table was in front of me. The phase error grew *linearly*
with frequency, which is what a time shift always looks like. And which probe size was hit moved
between runs, because a knife-edge comparison is decided by rounding, not by physics. I chased the box
transient, the patch footprint, `a/dx` and `steps_per_sample` instead — all correctly ruled out, none
of them the cause.

**Fix:** project on absolute `t`. The window start then cannot enter the answer at all. That is the
transferable lesson: *a phase reference tied to where you happened to start looking is a reference
that depends on rounding.*

Previously ruled out, and still correct: lock-in windowing arithmetic (identical at every point —
`T_slow/sample_dt` = 1024, `n_samples` = 8192), CFL, wall margin, `a/dx`, `steps_per_sample`.

**Ruled out — lock-in windowing.** Reproduced on the host with the field replaced by an analytic
first-order response, isolating the analysis from the physics: identical at all four probe sizes —
2048 samples, `T_slow / sample_dt` = 256.000 exactly, exactly 4.0000 slow periods in the window, 0.13°
phase error at every tone. (`aleph/scripts/c1_diagnose_lockin_windowing.py`.) **Also ruled out** — CFL
(`dt/dt_max` 0.2002–0.2005 everywhere) and wall margin (monotone in `a`, nothing special at 1.1).

## 8e. Method 3 (creep) FAILED its own grid check — which is why the check was written

`tau_creep` across `dx` = 0.30 / 0.20 / 0.15 µm at fixed `a` = 1.0: **0.03405 / 0.01830 / 0.01444 s**.
A 2× refinement moves the answer by **2.36×**. The `a` exponent comes out **0.764** against a declared
2, and the method's verdict is FAIL.

This was predicted in the driver's own docstring before the run: the patch is loaded instantaneously
into an empty medium, so the initial flux is set by the gradient at the patch edge, which is `~1/dx`.
The `1/e` crossing sits inside that grid-dominated region. **The creep number is mostly mesh.**

It is reported, not repaired. A method that fails a grid-invariance check is a real answer to "what do
the other ways of measuring this give" — some of them give you the mesh — and a poroelastic medium
with no intrinsic compliance is arguably the wrong material for a creep protocol in the first place.
What the run establishes is that you cannot tell from the creep data alone: its `D` exponent is
−1.0000, and its `a`-sweep r² is high enough to publish.

**Ruled out — lock-in windowing.** The windowing arithmetic was reproduced on the host with the field
replaced by an analytic first-order response, so any deviation would be the windowing in isolation. It
is *identical* for all four probe sizes: 2048 samples, `T_slow / sample_dt` = 256.000 exactly, exactly
4.0000 slow periods inside the measurement window, and a phase error of 0.13° at every tone for every
`a`. Tone leakage from non-integer periods is not the cause. (`scratchpad/diagnose_a11.py`.)

**Also ruled out — the CFL and the wall.** `dt/dt_max` is 0.2002–0.2005 at every point, and the wall
margins (10.6 / 7.7 / 5.7 / 4.2) are monotone in `a` with nothing special at 1.1.

**Still open.** The cause is on the physics side of the run, not the analysis side. The remaining
candidates are the box transient — the measurement window is `[2.011a², 4.021a²]` s and the box's own
diffusive time is 4.67 s, so `a` = 1.1 is the only point whose window straddles it — and the
discretised patch footprint. Neither is established, and the ordering of the transient overlap does
not by itself single out 1.1, so the leading candidate does not yet explain the observation either.
Isolating it needs a run varying box size at fixed `a`, which is not on the critical path: the
self-similarity guard rejects the point on an invariant that holds regardless of cause.

## 9. The method inventory — what each one assumes to turn a curve into a number

Added 2026-08-10 after the PI asked that the *other* ways of measuring this be built and swept too. The
point of the table is not that there are several methods; it is that each carries an **inversion**, the
inversion carries **assumptions**, and the assumptions are invisible in the number that gets published.

| # | Method | Driver | What the instrument records | The inversion, and what it assumes |
|---|---|---|---|---|
| 1a | step relaxation, whole-field | `c1_biot_timescale.py` | `‖p−p̄‖₂` after a step | `1/e` of the norm. Assumes you can see the whole field — **no instrument can**; this is the reference no laboratory owns |
| 1b | step relaxation, point probe | same | `p` at one point | `1/e` of a local reading. Assumes the point is representative — a tracer bead's assumption |
| 1c | step relaxation, finite contact | same | mean `p` over a contact patch | `1/e` of a patch mean. Assumes the contact is small compared with the gradient — the AFM assumption, and the one with **no closed form** |
| 2 | oscillatory microrheology | `c1_biot_microrheology.py` | amplitude + phase vs frequency | `τ = 1/ω₄₅`. Assumes the response is **first order**; a medium with two clocks has no single crossing and the method reports one anyway |
| 3 | creep under held load | `c1_biot_creep.py` | flux needed to hold the load | `1/e` of `dQ/dt`. Assumes the load was applied instantaneously — **its initial reading is set by the mesh**, which is why it alone needs a grid-invariance verdict |

Read the last column as a list of ways the same medium can produce different published numbers without
anyone making an error. That is the C-1 claim, and it is why the scored scalar is an exponent.

**What none of them can do, and the engine can.** Turn off one clock at a time. The three switches
(`Lp`, `xl_koff_per_s`, `biot_fsi`) each carry a regression gate asserting OFF is bit-identical, so the
attribution is a controlled experiment rather than an argument about which term dominates.

**The honest negative outcome is declared with the positive one.** If the exponents disagree across
methods, the scaling is *not* a property of the medium and C-1's choice to score an exponent is wrong.
`c1_method_synthesis.py` reports that case explicitly rather than filtering it.

## 10. Reference geometry note

The pressurised-shell oracle (`oracles/mechanics/thin_shell.py`) applies to the *tension*-dominated limit
and carries the Reissner zero-pressure stiffness `k₀ = 4Eh²/(R√(3(1−ν²)))` (Vella 2012). PSR must confirm
which regime each `a` sits in before reading τ, because bending-dominated and tension-dominated indentation
of the same shell are different measurements. Hertz (`oracles/mechanics/hertz.py`) is the solid-sphere
form and is **not** the right oracle for a pressurised shell; it is retained for the ECM-gel line.
