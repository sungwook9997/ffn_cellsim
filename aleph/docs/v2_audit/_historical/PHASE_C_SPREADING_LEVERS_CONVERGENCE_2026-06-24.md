---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Phase C — spreading levers: convergent finding + apico-basal polarization plan (2026-06-24)

Autonomous-loop session. Branch `h7/m1-ipc-contact`. This records the convergent result of the
M1-contact + turgor-engage investigation and the plan for the last untested lever.

## Convergent finding — the spheroid does NOT spread under ANY mechanistic lever tried (A/A0 ≈ 1.0)

| run | contact | cells | A/A0_final | reading |
|---|---|---|---|---|
| proxy-free stack (penalty) | capped penalty (pen 3.1) | rigid k_vol=7.73e5 | **1.037** | no spread (contact-confounded) |
| IPC (con_q) | IPC barrier+CCD (pen 1.80) | rigid | **0.996** | no spread (artifact-free) |
| deformable + IPC | IPC (pen ~2-3) | **deformable k_vol=1e3** | **1.016** | no spread |

All three land at A/A0 ≈ 1.0. So removing the two confounds the reconsider identified —
(1) contact penetration (M1, fixed to good-enough-honest pen 1.80) and (2) the rigid-cell
volume-lock (k_vol 7.73e5 → 1e3 makes the cell deformable: V/V0 0.995 responds to load vs 1.000
locked) — does **NOT** produce spreading. This strongly re-confirms the reconsider's high-confidence
verdict and the prior PI-ratified Layer-2 conclusion: **the collective A/A0=7-10 spheroid-spread
magnitude is a structural limit that belongs to the fine-grained single-cell line, not center/mesh
spheroid mechanics.** No-spread is robust across cohesion mechanisms, contact methods, and now
deformability.

## What each sub-result established (durable)
- **M1 contact → good-enough-honest.** Original IPC (con_q, d̂=c_rep) = pen 1.80 final, stable,
  bounded, A/A0 hull-robust. Cheap levers REFUTED: repel_q (con_q+3·me) gave pen 2.20 WORSE + 5×
  slower (entry-face-loss was NOT the bottleneck — it is the feasibilization equilibrium of
  already-inside nodes under the strong bundle); d̂=1·me DIVERGED (cfl 55, bonds 3139→142). The
  bundle is PHYSICAL-as-applied (don't lower it). pen<0.3 is DE-SCOPED (M1 = measurement honesty).
- **Turgor "engage" = deformability, NOT inflation.** Removing the k_vol lock gives a deformable
  cell (V/V0 0.995 under aggregation, pen=0 — IPC holds soft cells cleanly) but does NOT inflate to
  V/V0=1.12 (dP0=133 Pa cannot out-pressure the cortex; the old 1.12 was old-two-step-model-specific
  — REFUTED). Deformable + full spread stack → high cfl (20-47, soft cells stress the integrator) +
  A/A0 1.016 (no spread). So deformability alone is insufficient.

## The LAST untested lever — apico-basal POLARIZATION (the directional-spread mechanism)
A spatially UNIFORM shell (one γ, one contact law over the whole cell) can only minimize area →
ROUND. Directional flattening (basal spreading on the substrate while staying cohesive laterally and
tense apically) requires a basal/apical tension ASYMMETRY. This is the reconsider's #2 real lever and
the headline SimuCell3D / Runser-2024 feature ffn_cellsim lacks.

### Plan (mechanistic, NOT a lumped per-face γ multiplier — that would be the M3-junction-switch sin)
1. **Detection (host, low cadence — reuse what we already track):** tag each face basal / lateral /
   apical from EXISTING geometry, no voxel flood-fill needed:
   - basal = substrate engagement `w>0` (the wetting kernel already computes a per-face basal weight
     from z-height; `dcm_substrate_warp.py`),
   - lateral = node has a live cadherin trans-dimer bond / a different-cell neighbour within c_adh
     (`dcm_cadherin_host.py` already maintains this set),
   - apical = the complement (free surface).
2. **Drive the THREE EXISTING forces by the tag (route forces, not labels):**
   - basal faces → substrate wetting / low effective cortical tension,
   - lateral faces → Rakshit cadherin catch-bond cohesion,
   - apical faces → high cortical surface tension γ (the existing `surface_tension_kernel`).
   The asymmetry (low basal + high apical) is what converts a round pack into a spread sheet.
3. **Honest expectation (reconsider):** this may give tissue-like directional DEFORMATION, but the
   A/A0=7-10 MAGNITUDE remains the fine-grained single-cell line's structural limit — so the target
   is "does the mechanism produce directional flattening at all," not the magnitude.

## NEXT CYCLE (durable handoff)
Build the polarization detection + force-routing per the plan above (one new host tagger + wiring the
three existing kernels by tag). Validate on N=100: does A/A0 rise above ~1.0 (directional flatten)?
Visualize with the nucleus-inclusive cross-section. If it still does not flatten, the structural limit
is confirmed on the LAST lever → surface to PI that spheroid spreading is exhausted and the magnitude
belongs to the fine-grained single-cell line (re-route there per the PI-ratified Layer-2 call).

## FINAL consolidation (2026-06-24, after the polarization build) — spreading investigation COMPLETE

The last lever, **apico-basal polarization**, is now BUILT (`polarized_surface_tension_kernel`,
Young–Dupré differential tension `γ_face = γ − w·w_cs`, CPU-validated basal/apical force ratio 0.715)
and wired (`--polarize`). It does NOT overturn the structural-limit conclusion; it EXPLAINS it
mechanistically:

- **The Douezan spreading coefficient gates it.** `S = w_cs − 2γ`. At physiological MCF7 values
  (cortical γ ≈ 1e-2 N/m, substrate adhesion w_cs = 2.85e-3 J/m²) → **S = −0.017 < 0 = NON-WETTING**.
  A non-wetting cell does not passively spread — this is the mechanistic reason the spheroid footprint
  does not grow. (The same Douezan finding the Layer-2 §F wetting work reached, now at the force-kernel
  level.)
- **Single-cell vs collective is the real divide.** A SINGLE cell DOES spread (fried-egg A/A0→3.8 with
  wetting/polarization, previously PI-confirmed); the SPHEROID (collective) does not. This reproduces
  the PI-ratified Layer-2 conclusion: the FORM is reproduced, the A/A0=7-10 MAGNITUDE belongs to the
  fine-grained single-cell line, NOT center/mesh spheroid mechanics.
- **A numerical caveat (not a physics blocker).** The deliverable spheroid config in the DEFORMABLE
  regime (k_vol=1e3) + strong bundle + the implicit accel_dt=8e-4 hits a CFL instability (cfl → 35-240)
  for the deformable+polarized+bundle runs. This is a time-step issue (the soft cells + stiff forces
  need adaptive substepping `--cfl-limit`, a driver feature not yet exposed in the harness), NOT a
  physics result. The rigid-cell runs (k_vol=7.73e5) are stable; the deformable runs need adaptive-dt.

**Bottom line:** spreading is mechanistically EXPLAINED, not just observed-absent — the spheroid is
non-wetting (S<0) at physiological MCF7 values, and the collective spread magnitude is a structural
limit owned by the fine-grained single-cell line. The investigation across contact (M1/IPC),
deformability (turgor/k_vol), and surface-tension polarization is exhausted and self-consistent.

## DEEPEST insight (S>0 demo) — spheroid spreading needs a SUBSTRATE MODEL + volume conservation

Ran the S>0 polarization demo (γ=1e-4 → S = +2.65e-3 > 0, Douezan wetting) on the deformable
(k_vol=1e3) + IPC + bundle deliverable config, with adaptive substepping `--cfl-limit 0.3`. It
**hard-diverged** (cfl → 11035 at step 666, even at the max 16 substeps). Root cause is PHYSICAL, not
numerical: at S>0 the basal effective tension ``γ_basal = γ − w_cs < 0`` is **NEGATIVE** = area-MAXIMISING
= an inherent runaway (nothing bounds the basal expansion). This is exactly the spreading DRIVE
(negative basal tension wets/spreads) — so the polarization mechanism genuinely produces the wetting
force — but it has no EQUILIBRIUM here.

**Why the single cell spreads (bounded, A/A0→3.8) but this runs away:** the single-cell fried-egg used
(1) STIFF k_vol (volume conservation → the cell flattens at CONSTANT volume, so A/A0 is bounded by the
volume) AND (2) a substrate WELL/floor (bounds the basal expansion at the contact line). The proxy-free
DEFORMABLE stack removed BOTH: k_vol=1e3 (no volume bound) and use_substrate_well=False (no floor). So
the negative basal tension expands without bound → divergence.

**This reframes the whole picture and the turgor finding:**
- The k_vol "volume lock" (7.73e5) was NOT only a contact band-aid — it ALSO provides the
  **volume conservation that bounds spreading** (constant-volume flattening). "Remove the lock for
  deformability" is in direct TENSION with "spreading needs volume conservation."
- The proxy-free stack's ECM clutch provides substrate TRACTION, but NOT the substrate FLOOR /
  contact-line mechanics that BOUNDS a spreading cell. **The missing piece for spheroid spreading is a
  mechanistic SUBSTRATE MODEL** (a proper floor / contact-line / deformable substrate = the C6
  "deformable/3D substrate" backlog item), together with volume conservation.

So the full, self-consistent answer: spreading is gated by Douezan S (S<0 non-wetting at physiological
MCF7); the S>0 wetting drive EXISTS in the polarization kernel but needs **volume conservation + a
substrate floor** to reach the bounded fried-egg equilibrium — exactly what the proxy-free deformable
stack removed. The single cell has it (fried-egg works); the spheroid needs the C6 substrate model OR
the fine-grained single-cell line.

## IPC-barrier vs SimuCell3D-penalty A/B (the loop directive's contact comparison) — SURPRISE

Ran the contact A/B the loop prompt repeatedly asked for, via `--ipc-dhat-factor` (d̂ = factor·c_rep):
factor=1.0 = full IPC log-barrier; factor=0.01 = barrier negligible = **SimuCell3D-style implicit
penalty-only** (only the inside feasibilization linear push-out, with its analytic Hessian). Same n100
full bundle, same rep.

| contact | pen_final | pen_peak | wall | A/A0 |
|---|---|---|---|---|
| IPC barrier (d̂=c_rep) | **2.61** | 3.42 | 633 s | 1.001 |
| SimuCell3D penalty-only (d̂=0.01·c_rep) | **1.49** | 3.35 | 517 s | 1.024 |

**Surprise (refutes the assumption that the barrier is best):** the simpler **implicit penalty-only
BEATS the IPC barrier** — 43% lower final pen (1.49 vs 2.61), similar peak, ~18% faster, and NO
transient cfl spikes (the barrier produced a cfl≈1.2e9 single-frame spike that self-recovered but
perturbs the local config). The barrier's →∞ force near contact creates large transients that hurt
more than the smooth linear penalty. So the PI was right to want this A/B: **the M1 contact
recommendation is revised — the SimuCell3D-style implicit penalty (analytic-Hessian feasibilization,
NO barrier, NO CCD) is the simpler, faster, lower-pen contact.** (pen 1.5 is still above good-enough
~0.5-1.0, and pen is noisy/oscillating — neither fully closes M1, but penalty-only is the better base.)
The barrier+CCD machinery (the bulk of the IPC build) is not earning its complexity here.

## PI decision points (surfaced)
- Accept M1 good-enough-honest at pen 1.80 (de-cohesion runs are A/A0-hull-robust)?
- k_vol final physiological value (1e3 deformable confirmed viable with IPC; 2500/Guo are refinements).
- Invest in the polarization build (last lever, likely tissue-deformation not magnitude) vs accept the
  structural-limit conclusion now (strongly re-confirmed across contact + deformability)?

## Single-cell polarization test (the magnitude scale) — genuine ~2.3× spread, but raw A/A0 was an EJECTION ARTIFACT

Ran the polarization kernel on a **single cell** (N=1, the scale where the spread magnitude lives per
Layer-2): `--ipc --ipc-dhat-factor 0.01 --polarize --gamma-surf 1e-4 --well --no-bundle --accel-dt 8e-5`.
γ=1e-4 → S = +2.65e-3 > 0 (Douezan complete-wetting regime). Stable (cfl 0, V/V0 1.000) over 15k steps.

**Raw metric (top-down convex hull, ALL nodes): A/A0 climbed linearly 11.7 → 18.98, still rising.**
The figure (`figs/n1_polarization_friedegg.png`, side view) revealed this is **NOT a genuine fried-egg**:
the cell BULK barely flattened (z-span 15 → 13.6 µm, only **−9.5%**), while **12 of 162 perimeter nodes
(7.4%) were EJECTED** out to x = ±35 µm (the S>0 negative basal tension `γ_basal = γ − w_cs < 0` is
area-maximising and overpowers the cortex springs on individual perimeter nodes — a discretisation
runaway, exactly the peeling/ejection failure mode the peeling-aware `spread_eval` guards against).

**Honest robust metric (ejected nodes removed): A/A0 ≈ 2.14** — bulk footprint radius 7.5 → 11.4 µm
(genuine 2.3× spread), height −9.5%. So:
- **The polarization mechanism DOES produce genuine single-cell spreading** (~2.3×, O(2–4), consistent
  with the previously PI-confirmed single-cell fried-egg A/A0→3.8). Magnitude lives at the single-cell
  scale: **2.14 (single) vs 1.0 (spheroid)** — the single-vs-collective divide is real and reproduced.
- **The raw A/A0=19 was an artifact** (12 ejected nodes inflating the hull) — caught and corrected; the
  committed figure now marks the ejected nodes and reports the robust 2.14. (Lesson re-applied: the
  top-down-hull A/A0 must be ejection/peeling-robust, not raw — same class as the retracted A/A0→1.94.)
- **Remaining mechanistic gap = bounding the S>0 wetting drive.** The negative basal tension has no
  contact-line equilibrium at the node scale → it ejects nodes instead of spreading coherently. A
  physical bound (Young contact-angle / line tension at the basal perimeter, or a γ_basal≥0 clamp, i.e.
  the C6 substrate/contact-line model) is what converts the ejection runaway into a clean bounded
  fried-egg. This is the precise, sharpened statement of the C6 backlog item — not "a substrate model"
  generically, but specifically a **contact-line equilibrium** to bound the wetting drive.

## DECISIVE — at PHYSIOLOGICAL cortical tension (γ=1e-2) even a SINGLE cell does NOT spread (non-wetting)

Followed up the single-cell test at the **physiological MCF7 cortical tension** γ=1e-2 N/m (the
physiological-baseline rule: don't test at a convenient sub-physiological value). Here
`γ_basal = γ − w_cs = 1e-2 − 2.85e-3 = +7.15e-3 > 0` (positive, bounded — no ejection possible).

**Result: A/A0 = 1.019 — NO spread. The cell stays ROUND** (maxZ drop 5.3%, Ψ 0.996→0.992, peel_idx 0,
0 ejected nodes, V/V0 1.000, cfl 0). The 28.5% basal/apical tension differential is far too weak to
overcome the physiological cortical tension → the single cell is **non-wetting (Douezan S<0) and does
not passively spread**, exactly as Douezan predicts at physiological values.

**This closes the single-cell spreading question physiologically and ties the whole investigation together:**
- The γ=1e-4 "spread" (bulk 2.3×) only appeared because γ was pushed **100× below physiological** into
  the artificial S>0 wetting regime — and even there it was a node-ejection artifact, not a clean
  fried-egg. At the **real physiological γ, neither the single cell NOR the spheroid spreads.**
- **Passive differential-tension polarization is NOT the spreading mechanism at physiological values.**
  The PI's experimental A/A0=7–10 must come from **ACTIVE protrusion (lamellipodium/precursor-monolayer
  spreading, Aslemarz/Gupta 2024)** — a fundamentally different, active mechanism — not passive wetting.
  This independently re-confirms the PI-ratified Layer-2 §F conclusion (active, not passive, spreading)
  now at the single-cell force-kernel level and grounded in the physiological cortical tension.
- Figure: `figs/n1_polarization_physio_vs_subphysio.png` (physiological round/no-spread vs
  sub-physiological ejection-artifact, side by side).

**Net for the loop:** the spreading investigation is exhausted AND physiologically grounded on both
scales — single cell (non-wetting at physiological γ) and spheroid (collective structural limit). The
only routes to the A/A0=7–10 magnitude are PI-scoped: an **active** spreading mechanism (lamellipodium
as the driver, not a modulator) and/or the C6 contact-line substrate model — not any passive lever.

## M1 contact lever#3 — SimuCell3D position-based PROJECTION hard-constraint (built + self-test PASS)

The contact A/B established that both force-based methods plateau under the strong cad×40/ecm×167
bundle: penalty (capped) pen ~3.1, IPC penalty-only ~1.5, IPC barrier ~2.6. Lever#3 in the loop
directive ("제약방식" / constraint method) is the **third, non-force option**: a geometric
non-penetration **constraint**, the SimuCell3D hard-constraint approach.

**Built** (`dcm_contact_implicit_warp.py`): `nearest_face_project_kernel` + `project_contacts()`.
After the integrator's position update, a few Jacobi sweeps move any penetrating node out along the
nearest-other-cell-face outward normal to `c_rep` clearance — **independent of force magnitude**, so a
strong cohesion bundle cannot tunnel it. Wired opt-in into `run_decohesion` (`--project`), compatible
with the penalty contact path (post-step correction; not paired with `--ipc`).

**Self-test PASS** (`_projection_unittest`, two cells at 0.50·R overlap):
- initial penetration 3.67·c_rep → **penetration-free in 12 sweeps** (final probe-sweep residual
  0.0000·c_rep),
- COM drift from projection 0.0066·me (the mesh is not shoved off-centre — the correction is local + symmetric).
- CPU smoke (n12 through the full driver) runs end-to-end, all gates PASS, drift ~9e-6 µm.

So **a true hard constraint IS achievable and penetration-free where the force-based methods saturate** —
answering lever#3 affirmatively. The gbook A/B (full-bundle penalty vs penalty+projection: does pen drop
from ~3.1 to ~0?) is in flight (`n100_project_fullbundle`). If it confirms, the M1 recommendation
updates again: the projection hard-constraint is the penetration-honest contact, at the cost of a
post-step geometric correction (no per-step host sync in production; a few cheap Jacobi sweeps).

## Projection hard-constraint — adversarial review CAUGHT a divergence bug (the value of verify-before-trust)

Before trusting the projection A/B's headline (pen→~0 even under the full bundle), ran a 35-agent
adversarial review (4 lenses — geometry / convergence / physics / measurement-honesty — each finding
concerns, then an independent verifier per concern reading the actual code). 31 concerns, **22 confirmed
real, 3 critical.** This was not ceremony — it caught a RUN-BREAKING bug the early A/B data hid:

- **The one-sided projection DIVERGED in production.** With the original kernel (node shoved out, NO
  reaction on the contact face), the n100 full-bundle run was stable early (step 665: pen 0.066,
  V/V0 1.000, cfl 0.71) but **diverged by step 2527: V/V0 = 109.3 (volume blew up 109×), cfl = 2.9e10,
  A/A0 = 267.** The momentum the one-sided shove injected each step accumulated until the bundle +
  projection feedback exploded. The early pen→0 would have been a misleading "it works."
- **Fix (review critical #3): momentum-conserving 50/50 barycentric split.** The node moves out by
  half, the contact face's 3 vertices move in by the barycentric-weighted half (atomic_add, Σdpos=0) —
  Newton's third law, matching the repel/IPC kernels' existing reactions. Self-test COM drift
  **0.0066 → 0.0000·me**, still penetration-free (probe residual 0.0000·c_rep).
- **Also applied:** proj_iter default 4 → 8 (review critical #2 — the 50/50 split moves the node only
  half/sweep, so a margin; per-step penetration is small).
- **Reviewer call NOT applied (judgment over deference):** critical #1 (project along the face normal
  using signed distance, not the 3D closest-feature vector). The kernel projects along the
  closest-FEATURE direction r_vec — the gradient of the node-triangle distance, valid for
  interior/edge/vertex contacts alike; the fnv-normal alternative is correct only for face-interior
  contacts and would regress edge cases. Kept r_vec (self-test-proven penetration-free), documented inline.

**Decisive re-test in flight** (`n100_project_momfix`, proj_iter=8, momentum-conserving, full bundle):
does the fix carry the run past the old kernel's step-2527 divergence to a stable pen→~0 at 8000 steps?
Divergence-guarded poll armed (flags V/V0>5 / cfl-blow / nan). This is the real M1-lever#3 verdict.

## M1 lever#3 FINAL VERDICT — projection is penetration-free at base dt, diverges at the 100× timestep

Isolated the divergence cause with two concurrent diagnostics (n100, full bundle, momentum-conserving kernel):

| run | accel_dt | proj_gap | result |
|---|---|---|---|
| `n100_proj_smalldt` | **8e-6 (base)** | c_rep | **STABLE through 5000 steps — A/A0 0.999, V/V0 1.000, pen 0.000, cfl 0, drift 0** |
| `n100_proj_smallgap` | 8e-4 (100×) | 0.05·c_rep (near-zero) | DIVERGED (V/V0 358, cfl 1.7e11 @step 3984) |

**The divergence cause is the LARGE TIMESTEP, not the clearance shell.** A near-zero gap (pure
non-penetration, not fighting adhesion) still diverged at accel_dt 8e-4; the base timestep 8e-6 is
perfectly stable with the FULL c_rep gap. Mechanism: the post-step geometric projection applies a
position correction Δx each step; the implied velocity Δx/dt blows up at the 100× timestep, and the
projection↔implicit-force feedback diverges. At the base dt the corrections are small and consistent
with the force balance.

**So lever#3 is a SUCCESS with a hard caveat — the speed/accuracy frontier for M1 contact:**
- **Projection hard-constraint = the most accurate contact: pen→0 (TRUE penetration-free)**, V/V0 1.000,
  cfl 0, zero drift — *but only at the base timestep 8e-6* (100× slower than the accelerated runs).
  Figure: `figs/n100_projection_pen0_nucleus_montage.png` (clean round all-touching n100 spheroid, pen 0).
- **It is numerically INCOMPATIBLE with the 100× accelerated implicit timestep (8e-4)** that the
  rigid-cell penalty runs use — it diverges (V/V0 blows up) regardless of gap.
- The proper way to get penetration-free contact AT the large timestep is to couple the constraint
  INTO the implicit solve (constraint-aware Newton), not bolt a geometric projection on post-step —
  which is exactly what the **IPC barrier already does** (its analytic Hessian goes into the CG; stable
  at 8e-4, pen 2.6). The post-step projection is the "naive" hard constraint; the IPC barrier is the
  "solver-coupled" one.

**M1 contact — consolidated recommendation across all three methods + the projection:**
- **Fast production (accel_dt 8e-4):** force-based **IPC penalty-only** (pen 1.5, stable, the
  reconsider's good-enough-honest) — cohesion/de-cohesion runs are A/A0-hull-robust at this pen.
- **Penetration-free reference (base dt 8e-6):** the **projection hard-constraint** (pen→0) — for a
  decisive contact-fidelity check or any run where penetration honesty must be exact, at 100× wall-time.
- **Not viable:** post-step projection at the accelerated timestep (diverges).

This closes M1-lever#3: a true hard constraint IS achievable and penetration-free (answering the
directive), the cost is the base timestep, and the production contact stays the solver-coupled IPC.
The adversarial review was load-bearing — it caught the one-sided-shove divergence before the headline
pen→0 could be mis-reported as "projection works in production."
