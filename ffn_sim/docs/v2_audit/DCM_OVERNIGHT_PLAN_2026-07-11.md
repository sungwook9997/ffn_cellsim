# DCM overnight autonomous plan (PI asleep, self-deciding) — 2026-07-11

**Directive:** PI set an overnight /goal — work autonomously, /loop ~1 h health checks, self-decide when the session
would wait for PI, and when a thread finishes write a new plan and keep going. Backups every commit (nothing lost).

## State at start of the night
The DCM active-matter **jamming→unjamming** line is comprehensively delivered + rigorously validated + honestly
caveated (doc `DCM_UNJAMMING_ACTIVE_MOTILITY_2026-07-10`, ~8 figures, 2 Notion milestones, memory, backup
`dcm-backup-2026-07-10` on origin+github):
- Active per-cell motility unjams the jammed DCM spheroid → genuine collective spread; shape index → s0*=5.41.
- BOTH SPV axes (v0 force + persistence τ_p) reproduced; cadherin×motility DECOUPLE (motility=unjamming,
  cadherin=cohesion → epithelial/mesenchymal). Seed-robust, dt-converged, visual-verified, lit-validated
  (Bi-Manning/Park2015/Geiger2022) + PI KB-PIV-10.
- ⚠️ HONESTY: the force-mode demos ran at SUPER-PHYSIOLOGICAL SPEED (CFL 28); a v0-mode (`--v0-um-min`, SPV-faithful,
  CFL-bounded) was added — the honest Péclet knob. A CFL-safe v0-axis sweep is running.

## Overnight phases (prioritised; self-decide, re-plan as results land)

**Phase 1 — the HONEST v0-mode / Péclet phase diagram (primary).** Finish the CFL-safe v0-axis sweep; then map the
SPV transition in the physiologically-faithful (v0, τ_p) parametrisation:
- v0-mode is CFL-bounded, so sweep v0 (2…50 µm/min accelerated) × τ_p → the shape index / A/A0 phase surface, all
  numerically clean (no CFL confound). The order parameter is the Péclet number Pe = v0·τ_p / a (a = cell size).
- Anchor: physiological v0 ~1–2 µm/min, τ_p ~10 min → the physiological Pe; show where the transition sits.
- Gate: shape index crosses s0*=5.41 as Pe rises; jammed (s≈4.93) below. Ensemble the transition point (3 seeds).
- Deliver: a clean Péclet phase-diagram figure + a doc section that REPLACES the super-physiological-speed capstone
  with the honest v0/Pe statement. This closes the honesty caveat properly.

**Phase 2 — consolidate into a thesis-grade summary.** One clean SUMMARY doc + a flagship interactive HTML/animation
(confined → unjamming → spread, shape-index annotated, browser-verified) — the complete jamming/unjamming/EMT story
with the honest caveats. Every claim cross-referenced to TAG + a figure.

**Phase 3 — KB-PIV-10 SIZE axis (PI's own experiment).** PI's KB-PIV-10: aggregate unjamming by size (solid <25 µm /
liquid >25 µm / gas >10 cells). Test whether aggregate SIZE (N cells) affects the motility-driven unjamming in the
DCM (it may not — motility is per-cell; a null is itself an honest finding vs PI's size-dependence). Run N∈{50,100,200}.

**Phase 4 — (stretch) fluidisation ⊗ compaction ((b) line).** Does motility-fluidisation let the aggregate-σ
compaction become rearrangement-mediated (vs the S4 gap-closing bypass)? Uncertain; measure, don't force.

## Process (every step — HARD)
1. **Visual-verify** every A/A0/unjamming claim (browser_check screenshot, READ it) — has caught artifacts all night.
2. **TAG-compare** every quantity; **stability gate** every run (finite/G2/volume, watch CFL); **seed-ensemble** any
   motility magnitude claim (single runs are stochastic).
3. **No magic numbers** (v0/τ_p/s0* all lit-anchored, swept not tuned); no gate-loosening.
4. **Backup** every commit: `git commit` + refresh `dcm-backup-2026-07-10` (origin + github).
5. **Re-plan**: when a phase finishes, append the result + the next plan here; keep the loop productive.
6. If gbook drops (it did once tonight — network), the work is safe (all backed up); the monitors retry; resume on recovery.

## One-line
Close the unjamming line honestly (Phase 1 Péclet map replacing the speed-confounded capstone), consolidate it
thesis-grade (Phase 2), then extend to PI's size axis (Phase 3) — visually + TAG-verified, seed-robust, backed up, all night.

---

## NIGHT RESULTS (self-logged as phases landed)

### Phase 1 — v0-mode / Péclet (CFL-safe) — DONE, the honesty result
CFL-bounded v0-mode sweep (per-node force = γ_node·v0, physiological Stokes drag; τ_p=600s, N=100, 4.8s sim):
- v0 = 2 µm/min (physiological, 33 nm/s) → A/A0 = **0.997** (no spread; slight peel)
- v0 = 10 µm/min → A/A0 = **1.021**
- v0 = 50 µm/min (super-physiological) → A/A0 = **1.079**

**At physiologically-faithful speed the aggregate does NOT reach the unjamming transition in the accessible sim
time.** Clean quantification of *why*: at v0 = 2 µm/min a cell moves v0·t = 0.033 µm/s × 4.8 s ≈ **0.16 µm** over the
whole run — a small fraction of a cell diameter, far too little for a T1 neighbour-exchange. To rearrange (~1 diameter,
~15 µm) at physiological speed takes ~7–8 min of *simulated* time = ~5×10⁵ steps at accel_dt=8e-4 → ~15 h wall/run
(infeasible). The force-mode "unjamming" demos were REAL physics but reached the transition only by running cells at
**super-physiological speed** (compressing physiological minutes–hours of migration into ~5 s). This is the project
timescale gap, stated precisely; the v0-mode is the honest Péclet knob (Pe = v0·τ_p/a). **The mechanism is robust; the
absolute biology-time is not accessible.** Fig source: `bmbtvli46` v0-sweep log.

### Phase 3 — KB-PIV-10 SIZE axis — DONE, a kinetic (not steady-state) size effect
N∈{50,100,200} at fa100 / τ_p=600 / seed 11 / 4.8 s, shape index (s0*=5.41):
- N=50  → s = 5.235, frac_unjam 0.24, A/A0 1.86
- N=100 → s = **5.407**, frac_unjam **0.44**, A/A0 2.24  (unjams most)
- N=200 → s = **5.092**, frac_unjam **0.01**, A/A0 1.07, drift 0.24 µm  (**stays JAMMED**) — visual-verified
  (`s11_N200_last.png`): 200 intact cells, compact ~spherical, stress at junctions, NOT spread.

**Bigger aggregate stays jammed at fixed time/motility** — a KINETIC size effect (interior cells more constrained →
slower to rearrange; the periphery-to-volume ratio drops with N so the same per-cell motility fluidises less of the
tissue in the same window). This DIFFERS from PI's KB-PIV-10 *steady-state* size-dependence (bigger → more liquid),
and for the same reason as Phase 1: the steady-state phase needs physiological TIME; at fixed short sim time the DCM
shows the *kinetics* (bigger = slower to unjam), not the equilibrium size-phase. Single-seed (trend robust in sign;
magnitude would need an ensemble — but N=200's frac 0.01 vs N=100's 0.44 is far outside seed scatter).

### Convergent conclusion (Phases 1 + 3)
Both independent axes land on the SAME honest statement already in the main doc's caveat: **the unjamming MECHANISM
is robust and lit-validated (SPV both axes, shape index → s0*, cadherin/motility decoupling), but the physiological
STEADY STATE — whether at physiological v0 or the size-dependent equilibrium phase — requires physiological time that
is not accessible at the feasible integration step.** The DCM does the *mechanism* and the *kinetics*; it cannot reach
the *biology-time equilibrium*. This is the project timescale gap, now quantified from two directions. Nothing here
loosens or overstates; it tightens the caveat.

### Re-plan (next, autonomous)
Phases 1–3 are the honest close of the unjamming line. Remaining productive threads, in priority:
- **(a) Consolidate (Phase 2 finish):** fold Phase-1 v0/Péclet + Phase-3 size into the main doc as the two closing
  honesty sections (replacing/augmenting the speed-confounded capstone framing). Update memory + Notion milestone.
  → do now (non-GPU, safe).
- **(b) Phase 4 stretch (fluidisation ⊗ compaction):** does super-physiological-speed motility make the aggregate-σ
  compaction rearrangement-mediated (vs S4 gap-closing)? Same super-physio-speed caveat, but tests whether the (b)
  compaction line and the unjamming line are the same physics. GPU run; launch after (a) is committed.
- **(c) If both done:** the unjamming line is closed; hand back to PI a clean summary at wake. Do NOT invent new
  scope overnight beyond (a)/(b) — the honest close is the deliverable, not more experiments.
