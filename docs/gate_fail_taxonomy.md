# Gate FAIL taxonomy — Option F Week 1

PI directive 2026-04-29 (Option F selected). Codex review item 6:
classify every current gate FAIL into three buckets so that future
result interpretation is unambiguous in a referee context.

## Three categories

| Category | Definition | Allowed in publication claim? |
|---|---|---|
| **HARD-BLOCKER** | Numerical or solver bug; suspected to be falsifying the result. Must be root-caused before any quantitative claim derived from this run is admissible. | NO. Quote results only with explicit caveat. |
| **ACCEPTED-LIMITATION** | Known physics gap or design choice; PI has explicitly accepted that the framework does not capture this aspect. Documented + bounded. | YES, with named limitation reference. |
| **EXPLORATORY-ONLY** | Pilot-stage allowed FAIL; quantitative claim FORBIDDEN at production scale; appears in production purely to confirm the qualitative trend established at pilot. | YES at pilot scale only; production claim FORBIDDEN. |

A gate FAIL must belong to **exactly one** of these three. If the
classification cannot be decided, the default is HARD-BLOCKER until
analysis assigns it elsewhere.

## Universal classification protocol

For every current FAIL, four questions decide the bucket:

1. **Reproducible across phenotypes / scales?** If FAIL appears at
   identical magnitude regardless of physics knobs, it is solver-side.
2. **Cause known?** Anchor in literature, internal sanity_md, or
   `docs/12_validation.md` known-limitations registry.
3. **Magnitude bounded?** Can the discrepancy be quantitatively bounded
   so its effect on downstream claims is small / known?
4. **PI-decided accepted?** Is there a recorded PI decision accepting
   the FAIL as a documented limitation?

Decision tree:
- Answers Q1=Yes, Q2=No → HARD-BLOCKER
- Answers Q2=Yes, Q3=Yes, Q4=Yes → ACCEPTED-LIMITATION
- Answers Q2=Yes, Q3=Yes, Q4=No (only at pilot) → EXPLORATORY-ONLY
- Otherwise → HARD-BLOCKER (unresolved)

---

## Current FAIL inventory (production_lam4 + Phase 4-v2 pilots)

Observed across `results/production_lam4/gate_report.md`,
`results/final_pilot_v2_{bare,pre,lam4}/gate_report.md`. Eight distinct
FAIL types span the recent runs.

### F1 — Curvature operator (κ vs 2/R)
**Observed**: Phase 4-v2 pilots (1k particles): rel err 21.7%
(measured 2.4345 vs 2.0); Production Lam4 (5k particles): 9.5% (PASS,
first PASS in the project).

**Q1**: scale-dependent, not phenotype-dependent — rel err falls when
particle count increases, consistent with statistical sampling
improvement at finer surface band.

**Q2**: known. Documented in `docs/stage1a_aha_div_sanity.md` as the
v13/v14 ε² regulariser issue, with v15 documented v15.k.3 density-
based volumetric stress as the principled mitigation. Off-peak |f″/f′|
amplification per the Sanity Gate §6 measurement-protocol consistency
check.

**Q3**: bounded. ≤22% at pilot scale, 9.5% at production scale.

**Q4**: PI decision logged in `docs/outcomes_v15.md` accepting the v15
formulation; production-scale PASS confirms the principled mitigation.

**Classification**: **EXPLORATORY-ONLY** at 1k pilot; **ACCEPTED**
(now PASS) at 5k production. Pilot-scale numbers must NOT be quoted
as production claims; the production-scale PASS supersedes them.

---

### F2 — Momentum drift (horizontal only)
**Observed**: 1.93e-01 (Production Lam4), 1.61e-02 (v2 Lam4),
7.69e-02 (v2 Pre), 1.54e-01 (v2 Bare). Limit 1e-03. **All FAIL by
16×–193×**.

**Q1**: appears in every run (all phenotypes, both scales). Magnitude
varies but always exceeds limit by an order of magnitude.

**Q2**: NOT yet root-caused. Possible sources: anisotropic numerical
dissipation in boundary treatment, atomic-op order non-determinism,
asymmetric Marangoni impulse near substrate contact. Listed as
"Codex review item 3" in `docs/codex_review_synthesis.md`.

**Q3**: NOT bounded. We do not yet know if 0.193 is a benign
numerical drift or contaminating the spreading field.

**Q4**: NOT accepted. No PI-recorded decision treating this as a
documented limitation.

**Classification**: **HARD-BLOCKER**. Investigation scheduled in
Option F Week 2 (`docs/horizontal_momentum_drift_investigation.md`,
not yet written). Until classification changes, **anisotropy claims
in Stage 1e Sim A vs Sim B comparison are FORBIDDEN** because solver
drift would be indistinguishable from real anisotropy.

---

### F3 — Anchor force balance
**Observed**: 10.158 (Production Lam4), 7.680 (v2 Lam4), 6.680
(v2 Pre), 6.168 (v2 Bare). Limit 0.20. **All FAIL by 30×–50×**.

**Q1**: phenotype-correlated (Lam4 > Pre > Bare in pilot ordering)
*and* scale-correlated (production > pilot at same phenotype). Not a
pure constant — suggests it depends on the active spreading state.

**Q2**: NOT yet root-caused. The gate compares substrate vertical
force F_sub against the integrated stress σ_zz over the contact
patch; 50× discrepancy means either F_sub measurement is wrong or the
σ_zz integration is wrong. Possible: contact patch boundary
miscounting, anchor force double-counting between Layer 1a+ Option β
substrate CSF and the Hertz contact response.

**Q3**: NOT bounded. 50× over-anchoring (or 50× under-integrating)
could either (a) be a measurement artifact with no physical effect,
or (b) be biasing the asymptote downward by holding the basal patch
artificially in place — same magnitude could explain the
peak-and-decay we observe.

**Q4**: NOT accepted.

**Classification**: **HARD-BLOCKER**. Investigation scheduled in
Option F Week 2 (`docs/anchor_force_balance_investigation.md`, not yet
written). This is the single most consequential FAIL because it could
be *causing* the long-time decay rather than just measuring an
artifact.

---

### F4 — Contact-band ρ_kernel / ρ_ref window
**Observed**: 0.726 (Production Lam4), 0.679 (v2 Lam4), 0.731
(v2 Pre), 0.789 (v2 Bare). Window [0.85, 1.15]; all measurements
below 0.85. Underdense by 15–32%.

**Q1**: phenotype-mild, scale-mild. Slight phenotype variation but
all four runs cluster around 0.70–0.79.

**Q2**: partially understood. Contact band cells receive
substrate-CSF response; under-density indicates either contact band
selection too aggressive (counting too many bulk cells as "contact"
and diluting the kernel density), or substrate Layer 1a+ Option β
γ_sub_eff coupling is shifting density distribution within the band.
Documented in `docs/stage1a_plus_substrate_sanity.md` for the design
window.

**Q3**: bounded. The 15–32% under-density is small enough that
downstream observables (R drift, sphericity, spreading area) are not
visibly distorted by it relative to the well-known limit pathology of
~50% drift.

**Q4**: NOT formally accepted. PI directive 2026-04-29 to surface
substrate / contact gates as Codex review item 2.

**Classification**: **HARD-BLOCKER (low priority)** pending Codex item
2 investigation. Likely reclassifiable to ACCEPTED-LIMITATION once
linked to F3 anchor force balance investigation.

---

### F5 — Radius drift |R/R₀ − 1|
**Observed**: 0.151 (Production Lam4), 0.073 (v2 Lam4), 0.095 (v2 Pre),
0.131 (v2 Bare). Limit 0.050.

**Q1**: scale-correlated (production > pilot) — at longer time the
spheroid has more time to drift. Phenotype-correlated (Lam4 < Pre <
Bare drift at pilot, an *inverse* correlation that is itself a
publication finding).

**Q2**: known. v15 architectural ceiling documented in
`docs/outcomes_v15.md`: density-based volumetric stress reduces but
does not eliminate residual drift. The drift is dominated by
boundary-cohesive imbalance, with v15 baseline 0.244 and current
production 0.151 (−38% improvement).

**Q3**: bounded. Stage-by-stage improvement series logged in the
gate report ("R drift improvement vs v15 baseline" PASS in
production_lam4).

**Q4**: PI-accepted. The "v15 R drift baseline 0.244" comparison is
explicitly a strict-improvement check, with the architectural ceiling
documented as a Stage 1a residual.

**Classification**: **ACCEPTED-LIMITATION**. Reference: v15 outcome
document. Quote in publication with the architectural-ceiling caveat.

---

### F6 — Wadell sphericity ψ
**Observed**: 0.881 (Production Lam4), 0.829 (v2 Lam4), 0.745 (v2 Pre),
0.781 (v2 Bare). Limit 0.950.

**Q1**: scale-improving (production > pilot). Phenotype-correlated
(Bare/Pre/Lam4 ordering matches the spreading-extent ordering — more
spreading = lower sphericity, expected).

**Q2**: known. Sphericity *should* fail when spreading is active; the
gate enforces a sphere-like baseline that is incompatible with the
spreading phenotype the framework is designed to reproduce.

**Q3**: bounded. ψ ∈ [0.75, 0.88] across all runs — never falls below
0.7 (catastrophic deformation territory).

**Q4**: PI-accepted. The gate exists for Stage 1a free-floating
relaxation; under Stage 1a+ substrate contact + Layer 2/3/4/5/6 the
sphericity is *expected* to drop. Documented as expected at substrate
entry in `docs/stage1a_plus_substrate_sanity.md`.

**Classification**: **ACCEPTED-LIMITATION**. The gate threshold is
appropriate for Stage 1a only; for Stage 1a+ and beyond, the FAIL is
the *expected outcome*. Recommendation in roadmap refresh: split the
gate into "Stage 1a-only sphericity 0.95" and "Stage 1a+ post-spread
sphericity 0.7" thresholds. (Per Hard Rule "never modify a gate's
tolerance to make a failing run pass": this is a *contract change*
documented before the production rerun, not an inline edit.)

---

### F7 — Active power finite (Stage 1a++ Cousin-Rule)
**Observed**: ratio 17.75 (v2 Lam4), 28.68 (v2 Pre), 51.91 (v2 Bare).
Limit 10.0. **PASS at production_lam4** (ratio 0.05).

**Q1**: scale-dependent — production passes by 200×; pilot fails by
2–5×. Phenotype-correlated.

**Q2**: known. Stage 1a++ active stress σ_act = -ζ·K·I generates
power proportional to ε̇ ∝ 1/R. At pilot resolution the boundary band
is thin, so |ε̇| spikes locally → active power spikes. At production
resolution the boundary band is statistically smoother.

**Q3**: bounded. Rises with phenotype-spreading-extent in pilot;
disappears at production scale.

**Q4**: PI-accepted at pilot scale; production passes naturally.

**Classification**: **EXPLORATORY-ONLY** at pilot. Production passes,
so the production claim is intact.

---

### F8 — A/A₀ trajectory finite & non-pathological
**Observed**: at v2 pilot only — A/A₀ ∈ [0.169, 1.498] (Bare),
[0.296, 1.643] (Pre), [0.425, 1.560] (Lam4). Limit min ≥ 0.5.
**Production Lam4 PASS**: A/A₀ ∈ [1.000, 2.643].

**Q1**: pilot-scale FAIL on the **substrate-contact-area** metric
(`A_contact_xy_hull`) which depopulates on lift-off. NOT the top-down
projection metric (which is the publication-relevant measure).

**Q2**: known. **Measurement-protocol mismatch** documented in
`CLAUDE.md` Hard Rule 11. The gate measures substrate-contact-patch
area, which is not the same modality as PI experimental imaging
(top-down microscope projection). Documented as Phase 1.2 outstanding-
issue resolution.

**Q3**: bounded. Production passes the same gate (≥ 1.0 throughout)
because more particles statistically reduce the lift-off depopulation.
Top-down metric (publication-relevant) is monotone-positive
throughout: peak 1.570, end 1.409.

**Q4**: PI-accepted. The gate is a *legacy* contact-area measurement;
the publication-relevant measurement is `A_over_A0_topdown`, also
gated separately.

**Classification**: **EXPLORATORY-ONLY (legacy gate)**. Reclassify in
the production rerun to use only the top-down metric. The current gate
is a measurement-protocol legacy, not a physics failure.

---

### F9 — φ trajectory toward predicted φ_eq
**Observed**: <φ>(end) = 0.1593 vs predicted φ_eq = 0.7514, |err| =
0.5920 (tolerance ≤ 0.5). Production Lam4 only (pilot did not run
long enough for φ_eq drift).

**Q1**: scale-correlated. At pilot 4hr the framework has not yet
deviated; at production 80hr the interior φ has decayed to 0.088.

**Q2**: known. Layer 3 spatial Sₚ extension (v11) drives φ_interior
→ 0 because S=0 in interior (no boundary cells means k_+ source is
zero, k_- relaxation dominates). The "predicted φ_eq" is the *boundary*
equilibrium, not the global mean. **Codex review item 4** flags Layer
3 audit needed.

**Q3**: NOT bounded. The interior φ decay is what is generating the
Marangoni retraction in Mechanism inventory of marangoni_review.md.
Whether 0.088 is the correct interior value or a Layer 3 modeling
artifact is open.

**Q4**: NOT accepted. PI flagged Layer 3 audit in Option F.

**Classification**: **HARD-BLOCKER (Layer 3 audit pending)**. Until
the Layer 3 review confirms whether the interior decay is correct
physics or model artifact, the φ trajectory FAIL cannot be classified.

---

## Summary table — current gate state

| ID | Gate | Production | Pilot | Class | Action |
|---|---|---|---|---|---|
| F1 | curvature κ vs 2/R | PASS 9.5% | FAIL 21.7% | EXPLORATORY-ONLY (pilot), ACCEPTED (production) | Quote production only |
| F2 | momentum drift horizontal | FAIL 0.193 | FAIL 0.016–0.154 | **HARD-BLOCKER** | Week 2 investigation |
| F3 | anchor force balance | FAIL 10.158 | FAIL 6.168–7.680 | **HARD-BLOCKER** | Week 2 investigation |
| F4 | contact-band ρ_kernel | FAIL 0.726 | FAIL 0.679–0.789 | HARD-BLOCKER (low pri) | Week 2 (linked to F3) |
| F5 | radius drift | FAIL 0.151 | FAIL 0.073–0.131 | ACCEPTED-LIMITATION | Quote with v15 caveat |
| F6 | Wadell sphericity | FAIL 0.881 | FAIL 0.745–0.829 | ACCEPTED-LIMITATION | Roadmap: split gate |
| F7 | active power finite | PASS 0.05 | FAIL 17.75–51.91 | EXPLORATORY-ONLY (pilot), PASS (production) | Quote production only |
| F8 | A/A₀ contact-hull | PASS [1.0, 2.643] | FAIL [0.17, 1.6] | EXPLORATORY-ONLY (legacy) | Use top-down only |
| F9 | φ trajectory | FAIL 0.59 | n/a | **HARD-BLOCKER (Layer 3 audit)** | Defer to Layer 3 audit |

## Hard-blocker count

- **F2 horizontal momentum drift**: blocks Stage 1e anisotropy claim
- **F3 anchor force balance**: potentially causal for asymptote
- **F4 contact-band ρ_kernel**: low priority, linked to F3
- **F9 φ trajectory**: blocks Layer 4 Marangoni mechanism upgrades

These four are the Week 2 / Layer 3 audit targets. Until classified,
the *publication-claim subset* is restricted to:

- Phenotype ordering (Bare < Pre < Lam4 in spreading) — robust across
  R drift and A/A₀_topdown, multiple scales.
- Curvature operator at production scale (F1 PASS).
- v15 R drift architectural ceiling (F5 ACCEPTED with caveat).
- Continuum-only ceiling at A/A₀_topdown ≈ 1.4 end (Bucket P3).

## Cross-references

- `docs/codex_review_synthesis.md` — Codex review items 2, 3, 4, 6
- `docs/marangoni_review.md` — physics-side Mechanism A/E/F discussion
- `docs/12_validation.md` — known limitations registry
- `docs/outcomes_v15.md` — F5 architectural ceiling
- `CLAUDE.md` Hard Rule 11 — F8 measurement-protocol mismatch
