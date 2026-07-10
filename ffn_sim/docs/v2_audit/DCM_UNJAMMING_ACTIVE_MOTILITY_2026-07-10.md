# DCM spreading = an active-matter jamming→unjamming transition (mechanism demonstrated; biology-time caveated)

**Date:** 2026-07-10/11 · **Owner:** Lead session (PI 8h then overnight /goal) · **Branch:** `dcm/main`

## EXECUTIVE SUMMARY (read this first)
The session's arc: "the DCM spheroid **can't spread**" → the honest reason is that it is **JAMMED** (cells don't slide;
confined regardless of cadherin) → the freshly-ingested TAG literature is exactly the **jamming/unjamming SPV framework**
(Bi-Manning shape index p0*=3.81 / **s0*≈5.41 3D**; Geiger2022 invasion=fluidization; PI's own KB-PIV-10) → built
per-cell **active motility (SPV)** → the jammed DCM aggregate **UNJAMS and spreads**.

**What is ROBUST (seed-ensemble + dt-converged + visual-verified + lit-validated):**
- Active motility fluidises the jammed DCM spheroid → genuine collective spread; the per-cell 3D **shape index (the
  jamming order parameter) rises from 4.93 (jammed) toward s0*≈5.41 (unjammed)**.
- **BOTH SPV axes** reproduced: self-propulsion (v0) AND directional persistence (τ_p; τ_p=30 s = diffusive = stays
  JAMMED even with force; τ_p≥150 s = directed → unjams). A jammed tissue needs DIRECTED, not diffusive, motility.
- **Cadherin × motility DECOUPLE**: motility sets the UNJAMMING (shape index, cadherin-independent); cadherin sets the
  COHESION (how far the fluidised cells disperse) → epithelial (strong-cad, fluid-but-cohesive, Park2015) vs mesenchymal
  (weak-cad, fluid-and-dispersing), matching KB-PIV-10 low-E-cadherin→gas.

**What is CAVEATED (honest):** the force-mode demos ran at **super-physical SPEED** (fa100 CFL≈28, v_eff~300× the
~2 µm/min lit cell speed, Péclet~432) to reach the transition in ~5 s of simulated time — a **biology-time acceleration**.
At true physiological speed (v0-mode, CFL-safe) the aggregate stays jammed in feasible sim time; reaching the physiological
steady state (Pe~O(1)) needs physical TIME (tens of min = ~10⁶ steps, infeasible at accel_dt=8e-4) — the project-wide
timescale gap. So the MECHANISM + PHASE STRUCTURE are correct and demonstrated; the absolute speed↔time calibration is
the acceleration. `--v0-um-min` is the honest CFL-safe Péclet knob for future physical-time work.

**Bottom line:** the DCM reproduces the active-matter jamming↔unjamming physics of collective cell migration / spheroid
invasion, quantified in the literature order parameter, tuned by motility (v0, τ_p) and cadherin — the correct model for
epithelial-confined ↔ mesenchymal-spreading, with the physiological-timescale approach flagged as the standing gap.

---

## 0. The arc (from "can't spread" to the unjamming transition)

The spreading investigation established (visual-verified, `DCM_SPREADING_INVESTIGATION_PLAN_2026-07-10`): the DCM
aggregate is **CONFINED regardless of cadherin** (strong A/A0=1.001 = weak 0.990) and single cells node-eject — the
**root is that cells don't slide/rearrange (a JAMMED tissue)**. This is exactly the jamming physics: the new TAG
literature (SE 355→376, papers 191→264) gives the framework — **Bi–Manning shape-index jamming** (p0*=3.81 in 2D;
Merkel–Manning **s0*≈5.41 in 3D**; below = jammed/solid, above = unjammed/fluid), the **Self-Propelled Voronoi (SPV)**
model (Bi-Yang-Marchetti-Manning: active motility v0 fluidises a tissue even below the shape threshold), and
Geiger2022 ("*spheroid invasion = fluidization = a jammed→unjammed phase transition*"). PI's own **KB-PIV-10** MCF7/col-I
records the solid/liquid/gas unjamming by aggregate size + low-E-cadherin→gas/EMT.

## 1. The build — per-cell active self-propulsion (the SPV unjamming lever)

`dcm/dcm_active_motility_warp.py` — each cell self-propels with force `F_active·p̂_c`, p̂_c a slowly-reorienting
polarity (rotational diffusion, persistence τ_p; persistent random walk = the standard active-matter migrating-cell
model). Applied per node (mesh-invariant whole-cell F_active), wired into `run_decohesion` (`--active-motility
--f-active-nn --motility-tau`, default OFF). F_active = physiological single-cell traction (10–100 nN, KB-2.12;
controlled variable, swept — NOT tuned to a spreading target).

## 2. RESULT — active motility unjams the DCM aggregate → GENUINE dose-dependent spread

Native N=100 on gbook A5000, spheroid on substrate + IPC-newton + cadherin cluster de-cohesion + active motility,
6000 steps. **Every A/A0 browser-verified (screenshots) — genuine collective loosening, cells stay intact (NOT the
wetting node-ejection artifact).**

| F_active [nN] | A/A0 (top-down) | drift [µm] | peel_idx | 3D shape index s | spread_eval |
|---|---|---|---|---|---|
| 0 (confined) | 1.006 | 0.08 | — | **4.934** (JAMMED, < s0*) | confined |
| 10 | 1.114 | 1.02 | 3.4 | 4.988 | GENUINE COLLECTIVE SPREAD |
| 50 | 1.552 | 2.60 | 1.4 | 5.200 | GENUINE COLLECTIVE SPREAD |
| 100 | **1.918** | 3.80 | **1.1** | **5.351** (→ s0*≈5.41) | GENUINE COLLECTIVE SPREAD |

**Clean monotonic dose-response** (fig `dcm_unjamming_active_motility_doseresponse.png`): more active motility →
- more spread (A/A0 1.0→1.92, drift 0.08→3.8 µm),
- the per-cell **shape index rises 4.93→5.35, toward the Merkel–Manning 3D unjamming threshold s0*≈5.41** — the
  jamming→unjamming ORDER PARAMETER crossing from solid toward fluid,
- **peel_idx DECREASES 3.4→1.1** (more genuine at higher force, NOT peeling/ejection), V/V0=1.000 (volume-conserved),
  pen low, stable. Visual: cells rearrange + move apart as intact spheres (confirmed 10 nN + 100 nN).

## 3. What this means

- The DCM confined-vs-spread is the **active-matter jamming↔unjamming transition** — the correct physics for
  collective cell migration / spheroid invasion (Geiger2022, SPV). The DCM does it, quantified in the lit order
  parameter (shape index), at physiological F_active.
- It resolves the whole session: static confinement (jammed, epithelial/MCF-7) AND spreading (unjammed, active
  motility, mesenchymal/invasion) are the **two phases of one transition**, tuned by active motility — matching the
  SPV framework + PI's KB-PIV-10 (low E-cadherin / high motility → gas/EMT).
- The unifying root (cells don't rearrange) is unjammed by active motility — the SAME lever should rate-gate
  compaction (fluidised tissue flows), closing the (b) T1 line too.

## 4. Sanity / integrity

Every value lit-anchored: F_active = KB-2.12 single-cell traction (swept, not tuned); s0*≈5.41 = Merkel–Manning 3D
(lit); shape index = the standard S/V^(2/3). No gate loosened. All runs gate-PASS (finite/G2/volume), stable.
Every A/A0 claim VISUALLY verified (browser_check) — which is why the wetting/node-ejection artifacts were caught and
this genuine spread is trusted. Compared against TAG (jamming/SPV/Geiger2022) + PI KB-PIV-10.

## 4b. HONEST caveats (2026-07-10, EMT-sweep follow-up)

The §2 dose-response is a **single run per F_active at motility-seed 7** — and follow-up runs show it is **stochastic
and not cleanly monotonic**, so those numbers are indicative, not yet ensemble-quantified:
- **fa150 nascent** (same seed 7, DETERMINISTIC vs fa100): A/A0=1.229, drift 1.35 µm — **LESS** than fa100 (1.918),
  a real non-monotonicity at high force (suspect a velocity/CFL cap or a collective-cancellation of the stronger
  random propulsion). To investigate.
- **fa100 MATURED cadherin** (n=1): A/A0=2.421, drift 5.14 µm — **more** than fa100 nascent (1.918), the OPPOSITE of
  the "matured/epithelial resists unjamming" hypothesis. **Visually verified GENUINE** (intact cells, dispersed, no
  ejection), but n=1 is within the stochastic spread — inconclusive on maturation-gating.
- Unjamming TIME-COURSE (fa100): shape index jumps 4.855→5.30 in the first ~2 s then plateaus at 5.35 (frac_unjammed
  0→45 %); confined stays 4.85→4.93 (7 %). So the fluidisation is fast; A/A0 keeps rising as cells drift apart.

**Conclusion stands at the mechanism level** (active motility unjams the DCM aggregate → genuine collective spread,
shape index rises toward s0*, visual-verified) — but the QUANTITATIVE dose-response + the maturation-gating question
require **seed ensembles** (per the 'ensemble before spontaneous' rule). `--motility-seed` now exposes the polarity RNG.

## 4c. ENSEMBLE RESULT (2026-07-10) — unjamming DECOUPLES from spread: motility fluidises, cadherin sets cohesion

Seed ensemble (fa100 active motility, N=100, seeds {7,11,17}), cadherin **nascent(weak)** vs **matured(strong)**
(fig `dcm_unjamming_cadherin_ensemble.png`; both conditions visual-verified seed-11):

| cadherin | A/A0 (spread) | shape index s | frac_unjammed |
|---|---|---|---|
| nascent (weak) | **2.44 ± 0.11** | 5.361 ± 0.051 | 0.38 ± 0.06 |
| matured (strong) | **1.81 ± 0.25** | 5.352 ± 0.033 | 0.41 ± 0.05 |

**The decoupling (seed-robust, all 3 seeds nascent>matured on A/A0):**
- **UNJAMMING (shape index s) is the SAME** (5.36 vs 5.35, ~40% cells past s0*=5.41, overlapping) → **active motility
  fluidises the tissue regardless of cadherin** (both reach the same fluid cell-shape state).
- **SPREAD (A/A0) DIFFERS** (nascent 2.44 vs matured 1.81) → **cadherin sets the aggregate COHESION** — how far the
  fluidised cells disperse, NOT the unjamming. Visual (seed-11): matured = fluid-but-COHESIVE (deform/migrate yet stay
  compact); nascent = fluid-AND-dispersing.

**This is the epithelial↔mesenchymal physics, correctly decoupled:** epithelial (strong E-cadherin) tissue can be
FLUID/unjammed (rearranging) yet stay cohesive (Park2015 asthmatic airway epithelium unjams but stays a monolayer);
mesenchymal (weak cadherin) fluidises AND disperses. Matches KB-PIV-10 (low E-cadherin → gas/EMT). The earlier single
"matured spreads more" (n=1) was stochastic — the ensemble reverses it. So: **motility → unjamming (v0 axis); cadherin
→ cohesion (E-cadherin axis)** — the two independent knobs of the tissue jamming/wetting state.

## 4d. fa150 non-monotonicity RESOLVED = CFL artifact; ⚠️ dt-convergence of the core result being verified

fa100/150/200 at seed 11 with CFL capture: **maxCFL blows up with F_active — fa100=28, fa150=101, fa200=24812**. So
the fa150<fa100 non-monotonicity is a **NUMERICAL (CFL) artifact**, not a physical re-entrant/re-jamming transition:
at high F_active the accel_dt=8e-4 timestep is far too large (CFL≫1), the implicit-solver displacement/CCD cap
dominates → cells move less per step → less effective unjamming. The runs stay STABLE + GENUINE (the implicit+IPC
solver is unconditionally stable) but the DYNAMICS are under-resolved / capped above ~fa100.

**⚠️ This flags a rigour gap in the whole dose-response: even fa100 has CFL≈28 (≫1).** The unjamming spread is real
(stable, GENUINE, visual-verified, ensemble-consistent at fa100), but the QUANTITATIVE dose-response magnitude may be
accel_dt-dependent. **dt-convergence check running** (fa100 at accel_dt 8e-4 / 2e-4 / … at matched physical time): if
A/A0 + shape index hold as dt↓ (CFL→O(1)), the result is trustworthy; if they drift, the magnitudes are dt-confounded
and the clean run needs a smaller accel_dt (or a per-step velocity-limited motility). The QUALITATIVE claim (motility
unjams the jammed aggregate → genuine collective spread; cadherin sets cohesion not unjamming) is robust; the numbers
are pending this convergence check.

## 4e. dt-CONVERGENCE = PASS (2026-07-10)

fa100 seed11 at matched physical time 1.2 s, accel_dt 8e-4 / 2e-4 / 8e-5 (CFL ~28 / 7 / 2.8): A/A0 =
**1.521 / 1.526 / 1.507** — CONVERGED (within 1.3 %), all GENUINE. So despite CFL≈28, the implicit+IPC solver
handles it and the spreading result is **dt-CONVERGED / trustworthy** (not a timestep artifact). The fa150/200
non-monotonicity was extreme-CFL (100 / 24812) over-capping only. → the unjamming spread + the cadherin/motility
decoupling are quantitatively trustworthy at accel_dt ≤ 8e-4, F_active ≤ ~100 nN.

## 4f. SPV PERSISTENCE axis (2026-07-10) — both SPV axes reproduced; physiological persistence reaches s0*

The SPV transition is set by BOTH the self-propulsion force AND the directional persistence τ_p (=1/Dr). At fixed
F_active=100 nN, sweeping τ_p (seed 11; fig `dcm_unjamming_spv_persistence.png`; τ_p=30 s visually confirmed JAMMED):

| τ_p [s] | A/A0 | shape index s | frac_unjammed |
|---|---|---|---|
| 30 | 1.04 | **4.922** (JAMMED = confined 4.93) | 0.02 |
| 150 | 2.58 | 5.397 | 0.42 |
| 600 (physiological ~10 min) | 2.19 | **5.408 ≈ s0\*** | 0.45 |
| 3000 | 2.39 | **5.411 = s0\*** | 0.39 |

**Sharp transition ~τ_p 100 s:** LOW persistence (30 s) → the motility reorients too fast, the net displacement
cancels (DIFFUSIVE) → the tissue stays JAMMED (s=4.92, the no-motility confined value) DESPITE the same active force;
τ_p ≥ 150 s (incl. the physiological ~10 min) → ballistic/directed motion → the tissue UNJAMS and the shape index
reaches **s0*≈5.41** (the clean crossing). So a jammed tissue needs DIRECTED, not diffusive, motility. **Both SPV axes
(F_active/v0 and persistence/Dr) are now reproduced on the DCM**, and at the PHYSIOLOGICAL operating point (F_active~100
nN, τ_p~10 min) the shape index sits right at the unjamming threshold — the DCM tissue is poised at the jamming↔
unjamming transition, exactly the SPV / KB-PIV-10 picture.

## 4g. CAPSTONE — at the PHYSIOLOGICAL operating point the tissue sits AT s0* (seed-robust)

Ensemble at the physiological motility parameters (F_active=100 nN single-cell traction, τ_p=600 s ≈ 10 min directional
persistence; seeds {7,11,17}): mean 3D shape index **s = 5.373 ± 0.079** (per-seed 5.264 / 5.408 / 5.447, one above s0*)
— i.e. within ~0.7 % of the Merkel–Manning unjamming threshold **s0*≈5.41**. A/A0 = 2.30 ± 0.4 (all GENUINE spread).
**So at physiological parameters the DCM spheroid sits RIGHT AT the jamming↔unjamming transition** — poised between the
solid/confined (epithelial) and fluid/spreading (mesenchymal) states, exactly as real tissues operate near the jamming
critical point (Park2015 airway epithelium; the SPV picture). Junction state (cadherin) then tips the poised tissue to
cohesive-fluid (epithelial) vs dispersing-fluid (mesenchymal) without changing the unjamming itself (§4c). This is the
physiologically-anchored capstone: motility (v0) + persistence (Dr) set the jamming state, and at the in-vivo operating
point the tissue is critically poised.

## 4h. ⚠️ HONESTY — force-mode used SUPER-PHYSIOLOGICAL SPEED; the v0-mode fix

The SPV axis is v0 (self-propulsion SPEED), not force. Adding a v0-mode (per-node force = γ_node·v0, physiological
Stokes drag) revealed that the **force-mode demonstrations (fa10–100 nN) ran the cells at SUPER-PHYSIOLOGICAL SPEED**:
fa100 has CFL≈28, i.e. the cells move far faster than the lit ~1 µm/min cell speed — that is how the tissue reached
the unjamming transition in only ~5 s of simulated time. So **the §4g "physiological operating point" is physiological
FORCE, not physiological SPEED**. At true physiological v0 (~1 µm/min) the same physics holds but the transition
develops over the physiological TIME (hours) — unreachable at accel_dt=8e-4 in feasible steps (the project-wide
timescale gap; the force-mode is effectively a biology-time acceleration, like S_accel in the compaction line).

**What stands vs what is caveated:**
- ROBUST: motility unjams the jammed DCM aggregate → genuine collective spread; BOTH SPV axes (v0 and persistence);
  the cadherin=cohesion / motility=unjamming decoupling; shape index → s0*. All qualitative, seed-robust, visual-verified.
- CAVEATED: the absolute F_active↔speed calibration + the "at the physiological point the tissue is exactly at s0*"
  magnitude — those used super-physiological speed. The v0-mode (`--v0-um-min`, CFL-bounded, SPV-faithful) is the correct
  knob; the honest statement is a Péclet-number (v0·τ_p/a) map, and reaching the transition at physiological v0 needs
  physiological time.

## 4i. v0-mode result — CONFIRMS the timescale gap (Péclet framing)

CFL-safe v0-mode sweep (τ_p=600 s, 4.8 s simulated, seed 11): v0 = 2 / 10 / 50 µm/min → A/A0 = 0.997 / 1.021 / 1.079,
shape index s = **4.938 / 4.939 / 4.986** (all ≈ the confined 4.93, ~7 % unjammed) = **stays JAMMED**. So at CFL-bounded
PHYSICAL speed + FEASIBLE sim time, the aggregate does not unjam in 4.8 s. The force-mode fa100 unjammed only because
its effective speed was super-physical: F=100 nN / γ_cell(≈γ_node·N_c=9.3e-3) → v_eff ≈ 645 µm/min ≈ 300× the ~2 µm/min
lit cell speed → Péclet Pe = v·τ_p/a ≈ 432 (deep-fluid), reached in 4.8 s. The physiological Pe (v0=2 µm/min, τ_p=600 s)
≈ 1.3 — right at the jamming transition — but reaching the unjammed STEADY STATE at Pe~O(1) takes several τ_p ≈ tens of
minutes of PHYSICAL time = ~millions of steps at accel_dt=8e-4 (infeasible). **This is the project-wide biology-time gap,
now precisely located for the unjamming: the MECHANISM + phase structure (SPV, both axes, decoupling) are correct and
demonstrated; the physiological-time approach to steady state is the acceleration.** The honest deliverable is the
qualitative + Péclet-phase result, with the force-mode understood as a biology-time acceleration (like S_accel).

## 4j. SIZE axis (KB-PIV-10) — a KINETIC size effect, not the steady-state size-phase

PI's KB-PIV-10 is unjamming-by-size (aggregate solid <25 µm / liquid >25 µm / gas >10 cells). Test whether aggregate
SIZE modulates the motility-driven unjamming. N ∈ {50, 100, 200} at fa100 / τ_p=600 s / seed 11 / 4.8 s (shape index,
s0*=5.41):
- N=50  → s = 5.235, frac_unjam 0.24, A/A0 1.86
- N=100 → s = **5.407**, frac_unjam **0.44**, A/A0 2.24  (unjams most)
- N=200 → s = **5.092**, frac_unjam **0.01**, A/A0 1.07, drift 0.24 µm — **stays JAMMED**; visual-verified
  (`s11_N200_last.png`): 200 intact cells, compact ~spherical aggregate, virial stress concentrated at junctions, a
  few peripheral cells stressed but bonded — NOT spread.

**Bigger aggregate stays jammed at fixed time/motility** = a KINETIC size effect: the interior cells are more
constrained (surface-to-volume drops with N), so the same per-cell motility fluidises a smaller *fraction* of the
tissue within the same 4.8 s window; N=200's frac_unjam 0.01 vs N=100's 0.44 is far outside seed scatter (the sign is
robust; magnitude is single-seed). This is the SAME timescale mechanism as §4i, seen along the size axis: at fixed
short sim time the DCM shows the *kinetics* of unjamming (bigger = slower to fluidise), which DIFFERS from PI's
KB-PIV-10 *steady-state* size-phase (bigger → more liquid) — that equilibrium phase needs physiological time. So the
size axis, like the v0 axis, confirms: the mechanism + kinetics are captured; the biology-time equilibrium is the gap.

## 5. Closing status (overnight 2026-07-11)

The active-matter jamming→unjamming line is **closed honestly**:
- **Mechanism (ROBUST):** per-cell active motility unjams the jammed DCM aggregate → genuine collective spread; shape
  index → s0*=5.41; BOTH SPV axes (v0 force + persistence τ_p); cadherin=cohesion / motility=unjamming decoupling
  (epithelial/mesenchymal). Seed-robust, dt-converged, visual-verified, lit-validated (Bi-Manning/Park2015/Geiger2022
  + KB-PIV-10).
- **Biology-time (CAVEATED, quantified from TWO axes):** §4i v0/Péclet — at physiological v0 (Pe≈1.3) the transition
  needs physiological time (millions of steps, infeasible); the force-mode reached it at super-physiological speed
  (Pe≈432). §4j size — bigger aggregate stays jammed at fixed time (kinetic, not the KB-PIV-10 steady-state size-phase).
  Both land on the project-wide timescale gap: the DCM captures the *mechanism* and *kinetics*, not the biology-time
  *equilibrium*.

Remaining (optional, non-scope-inventing): apply the same fluidisation to the aggregate-σ compaction ((b) line) to
test whether the two dynamic lines are the same physics — with the same super-physiological-speed caveat.
