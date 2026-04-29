# Codex external review — synthesis with Marangoni review (2026-04-29)

External reviewer (Codex) provided a 10-item gap analysis. Cross-
referenced with `docs/marangoni_review.md` (committed d0a7d99) and
current Production Lam4 gate state. This document is the consolidated
PI decision input. Review-only; no code changes.

---

## Reviewer's four-bucket framing

1. **Physics mechanisms** (items 1, 4, 5)
2. **Validation contracts** (items 2, 3, 6, 8, 10)
3. **Code structure** (item 9)
4. **Documentation coherence** (item 7)

Reviewer's bottom line:
> 지금 제일 부족한 건 "세포가 진짜로 가장자리에서 기어 나가는
> 사건성"이고, 그다음은 검증 실패들을 연구 결과와 코드 결함으로
> 구분하는 정리 작업.

This **agrees with the marangoni_review.md** Mechanism A/E/F finding
on the physics axis (continuum cannot capture discrete edge events
above the Pearson 1958 Bénard threshold in our overdamped regime), and
**adds** four engineering-discipline items the Marangoni review did
not cover: gate fail taxonomy, momentum drift root-cause, doc
coherence, monolithic file split, parameter literature registry.

---

## Item-by-item synthesis

### Item 1 — Discrete boundary biology missing (HIGHEST priority)

Reviewer: spreading to A/A₀ = 8–33 needs lamellipodia / filopodia /
leader cells / discrete focal adhesion / edge-cell traction asymmetry /
neighbor exchange — the continuum framework saturates around 2–3.

Marangoni review §4: this is **Option β (Stage 1a++.b)** — orthogonal
to the Marangoni A/E/F path.

**Synthesis**: reviewer's item 1 + Marangoni review's Option β +
Mechanism D (Bénard-Marangoni cannot occur in overdamped continuum)
form a single coherent theoretical case. The continuum framework's
peak-and-decay is *not a bug* — it is the linear consequence of an
incomplete theory. The Pearson 1958 instability that would generate
discrete edge events spontaneously in inviscid thin films is suppressed
in our overdamped Stokes regime; therefore we need either (a)
explicit stochastic edge events (Stage 1a++.b) or (b) a much more
elaborate γ-dynamics model (Mechanisms A+E+F) — possibly both, sequenced.

**Action implication**: Stage 1a++.b is justified independently of
Marangoni mechanism upgrades. They are not redundant.

### Item 2 — Substrate / contact mechanics still failing gates

Reviewer flags four production_lam4 FAILs:
- anchor force balance (median 10.158 vs limit 0.20)
- contact-band ρ_kernel/ρ_ref (median 0.726 vs window [0.85, 1.15])
- radius drift (max 0.151 vs limit 0.050)
- sphericity (min 0.881 vs limit 0.950)

These are **not addressed in the Marangoni review**. The Marangoni
review tacitly accepts them as background; reviewer is right to flag
that this acceptance has not been formally documented. The substrate
contract has not closed.

**Synthesis**: at least one of these (anchor force balance, 50× over
limit) is severe enough that it could be a *cause* of long-time decay
rather than just a symptom — if the substrate is not balancing
vertical force properly, lateral spreading driven by the same anchor
would also be biased.

**Action implication**: substrate force balance investigation has to
precede Stage 1a++.b. If anchor force balance is broken at the
solver level, adding stochastic events on top will not fix the
asymptote — and worse, the events would be calibrated against a
biased baseline.

### Item 3 — Horizontal momentum drift

Reviewer: |Δp_xy|/(m·v_rms) = 0.193 (vs limit 0.001) means the
horizontal momentum drift is 193× the gate. Vertical drift is
substrate-absorbed and ungated; horizontal should be conserved by
symmetry.

This is a **possible solver bug**, not a physics omission. None of the
Marangoni review's six mechanisms (A–F) generates horizontal
asymmetry; if anything they preserve isotropy. The 193× drift suggests
either (i) anisotropic numerical dissipation in the boundary
treatment, (ii) accumulation in atomic ops without reduction-order
control, or (iii) asymmetric Marangoni impulse application near the
substrate contact.

**Risk**: in Stage 1e (Sim A vs Sim B comparison), this drift would be
*indistinguishable* from a real anisotropy finding. If not fixed before
the comparison, the radial-vs-full conclusion is contaminated.

**Action implication**: this is a **Hard Blocker** for Stage 1e and
for any anisotropy claim in the paper. It is an *Exploratory-only* fail
at the current stage where we are not making anisotropy claims, so
production-stage progress is allowed to continue. But it must be
classified explicitly.

### Item 4 — Layer 3 φ dynamics too simple / spatially awkward

Reviewer: φ_predicted vs measured |err| = 0.59 (gate FAIL). Interior
S=0 driving φ → 0 in interior, generating ∇γ that retracts. The
reviewer questions whether the single-scalar φ E-cad/integrin
representation is adequate.

The Marangoni review's Mechanism A/E/F **inherits this Layer 3 issue**
— γ(φ) is only as good as φ. If φ is structurally wrong, fixing γ
dynamics on top of it is fitting on sand.

Reviewer's specific gaps:
- boundary/interior φ source definition unclear
- phenotype-specific initial φ mapping needs literature anchoring
- no spatial diffusion / cell-cell communication
- E-cad/Integrin-β1 ratio compressed into single scalar

**Synthesis**: φ is a deeper structural issue than Marangoni. The
Cho 2020 ODE was anchored at the cell-population level, not for
spatially-resolved spheroid interiors. The S_p spatial extension was a
v11 patch, not a first-principles derivation.

**Action implication**: before adding Mechanism A or E to Layer 4,
audit Layer 3. If the Layer 3 φ field has structural issues, the
γ(φ_p) = γ_max(1−φ) + γ_min·φ map is *the* mechanism propagating
those issues into Layer 4.

### Item 5 — Osmotic / turgor coupling insufficient at long time

Reviewer: ρ_osm desaturates from pilot 1.600 → production 1.224. No
sustained turgor amplification.

Marangoni review Mechanism F couples Layer 5 → Layer 4 (γ depends on
ρ_osm). Reviewer goes deeper:
- water efflux/influx timescale needs revisit
- local strain-rate coupling at edge insufficient
- cortex stiffening ↔ boundary traction coupling missing

**Synthesis**: Mechanism F is the right *direction* but reviewer
correctly notes the Layer 5 model itself may be incomplete. Guo 2017
turgor model is for *single cells*; spheroid-level turgor with cell-
cell water transport (gap junctions, paracellular) is not present.

**Action implication**: Mechanism F should be paired with a Layer 5
audit similar to the Layer 3 audit above.

### Item 6 — Gate fail taxonomy missing (HIGHEST engineering priority)

Reviewer's three categories:
- **Hard blocker**: must fix before next-stage interpretation
- **Accepted limitation**: PI-decided documented failure
- **Exploratory-only fail**: pilot allowed, production-claim FORBIDDEN

Currently the production_lam4 gate report has **5 FAILs** (anchor
balance, contact-band ρ ratio, radius drift, sphericity, φ trajectory)
and the run is treated as a complete result anyway. This is not
defensible at publication time.

**Synthesis**: this is the single most important *non-physics* item in
the reviewer's list. Without taxonomy, every gate FAIL is ambiguous,
and the cumulative claim "framework reproduces phenotype ordering" is
unfalsifiable in a referee context.

**Action implication**: write `docs/gate_fail_taxonomy.md` *before*
any further code. Classify every current FAIL into the three
categories. This is a 1–2 day documentation task with no code risk.

### Item 7 — Documentation / actual state mismatch

Reviewer: `docs/10_dev_roadmap.md` Stage 0 says "not yet started"
(verified line 7 says Stage 0 — but lines 23–28 do mark it complete;
the rest of the roadmap is stale though). `SESSION_HANDOFF.md` and
stage outcomes documents are scattered.

**Synthesis**: this is correct. After Stage 1a → 1a+ → 1a++ → 1b →
1c → 1d → 1e → 2, the roadmap document is months of work behind. New
collaborators would be confused immediately.

**Action implication**: roadmap consolidation is a 1-day documentation
task. Should happen alongside Item 6.

### Item 8 — Tests are shallow (Taichi-skipped)

Reviewer: 11/12 tests pass but Taichi-dependent ones skip — meaning
the core solver is not under CI.

Reviewer's missing-test inventory:
- curvature static sphere regression
- top-down projection area regression
- φ ODE bounds / trajectory unit test
- Layer 5 osmotic bounds
- Layer 6 ECM/MMP monotonicity
- substrate force balance synthetic test
- horizontal momentum conservation synthetic symmetry test

**Synthesis**: this is engineering-grade discipline. The horizontal
momentum drift (Item 3) is the canonical case where a synthetic
symmetry test would have flagged the bug at Stage 1a+.

**Action implication**: add the symmetry-conservation tests *first*
(they are the cheapest and most disciplinary); the Layer-specific
bound tests are nice-to-have.

### Item 9 — `mlsmpm.py` monolith

Reviewer: 100k+ bytes in one file containing Layer 1–6 + diagnostics +
substrate + measurement.

Marangoni review Mechanism A/E adds two new state variables (γ_p, Γ_p)
and at least one new ODE / surface-diffusion kernel. Adding these into
the current monolith makes it strictly worse.

**Synthesis**: file split is a *prerequisite* for any Mechanism A/E
implementation, not a parallel task. The Layer 4 logic in mlsmpm.py
(lines 307–947) is the natural target file to extract first.

Reviewer's split proposal:
- `constitutive.py` (Maxwell + Neo-Hookean)
- `csf.py` (curvature, surface tension grid)
- `substrate.py` (Layer 1a+ Option β + anchor)
- `boundary_events.py` (Layer 2; Stage 1a++.b ready)
- `adhesion_phi.py` (Layer 3)
- `osmotic.py` (Layer 5)
- `chemistry.py` (Layer 6)
- `diagnostics.py` (invariants, force-balance, mass)

**Action implication**: file split is the *first* code change that
follows Items 6+7 documentation work, before any new physics.

### Item 10 — Literature anchoring partial

Reviewer's 4-tier registry:
- direct literature value
- derived from dimensional analysis
- PI-approved phenomenological placeholder
- accepted limitation, not used for fitting

Currently Magic-Number-Block PARTIAL pattern (v15 ρ_floor, etc.) is
the closest thing we have, but it is per-decision not per-parameter.

**Synthesis**: writing `docs/parameter_registry.md` consolidates Magic-
Number-Block PARTIAL entries into a single auditable table. 1 day of
work, blocks nothing.

---

## Reviewer's priority ordering vs Marangoni review's PI options

Reviewer's order (1 → 5):
1. Gate fail taxonomy (Item 6)
2. Horizontal momentum + substrate force balance hard-blocker decision (Items 2, 3)
3. Stage 1a++.b discrete boundary events design (Item 1)
4. Lam4 4hr pilot — does the plateau break under Stage 1a++.b? (validation step)
5. Then production / paper writing

Marangoni review's PI options:
- α: Stage 1d.b (Mechanisms A/E/F)
- β: Stage 1a++.b (discrete events)
- γ: sequenced α → β
- δ: paper-as-is

**The reviewer's order maps to a refined Option β-first sequencing**,
with two new prerequisites the Marangoni review did not surface:
documentation taxonomy + substrate / momentum bug investigation.

---

## Consolidated PI decision options (revised)

### Option E — Reviewer-aligned engineering-first (RECOMMENDED)

Sequence:
1. **Documentation pass** (1 week): gate fail taxonomy, roadmap
   refresh, parameter registry, Marangoni review label fix
   (production_lam4_finding.md). No code.
2. **Solver hardening** (1–2 weeks): root-cause horizontal momentum
   drift; root-cause anchor force balance; if either is a solver bug,
   fix; if a physics gap, document and reclassify into "accepted
   limitation".
3. **Test discipline** (1 week): add the 7 missing tests from Item 8.
   Symmetry-conservation tests first.
4. **File split** (1 week): extract Layer 2 / 3 / 4 / 5 / 6 modules
   from `mlsmpm.py`.
5. **Stage 1a++.b** (2–3 weeks): design + implement discrete boundary
   events on the cleaned-up baseline.
6. **Pilot 4hr Lam4** (1 day): validate plateau-breaking at pilot scale.
7. **Decision point**: production rerun OR continue to Marangoni A/E/F.

**Pros**: every step is a low-risk, well-scoped commit. The solver-bug
investigation in step 2 *could* explain part of the asymptote gap by
itself. The cleaned-up baseline makes Marangoni A/E/F much easier
later.

**Cons**: 6–8 weeks before any new physics. Slow on the publication
clock.

### Option F — Reviewer-Marangoni hybrid (FAST)

Sequence:
1. Gate fail taxonomy + roadmap refresh + parameter registry (3–5 days
   documentation only).
2. Single-issue solver investigation: horizontal momentum drift (1 week).
3. Stage 1a++.b *and* Marangoni Mechanism A in parallel (3 weeks each
   with file-split as joint prereq).
4. Production rerun.

**Pros**: trims 3–4 weeks of pure infrastructure off Option E without
abandoning the discipline items.

**Cons**: Layer 3/5 audits deferred; substrate / contact gate FAILs
remain unresolved.

### Option G — Paper-as-is + parallel infra cleanup (LOWEST RISK)

Sequence:
1. Documentation + taxonomy + registry (1 week).
2. Manuscript draft on existing 27-simulation narrative + Production
   Lam4 P3 finding + Mechanism D (Pearson 1958 overdamped suppression
   analysis, no code).
3. *In parallel*: solver hardening + tests + file split (background
   over the writing weeks).
4. Stage 1a++.b only after submission, as the second paper.

**Pros**: maximum publication speed on the existing finding. The P3
"continuum-only ceiling" result is genuinely paper-worthy.

**Cons**: leaves the most interesting upgrade (Mechanism A Yadav
reaccumulation) for paper 2.

### Option H — Status-quo Stage 1d.b (DEPRECATED by reviewer items)

This is the original Marangoni review Option α (Mechanisms A/E/F
without the engineering items). **Reviewer's items 2, 3, 6, 7
effectively veto this option** — building new physics on top of the
unresolved gate FAILs and broken roadmap is too high-risk.

---

## Reviewer's points that the Marangoni review missed entirely

1. **Substrate force balance 50× over limit** — `anchor_force_balance`
   median 10.158 vs limit 0.20. Marangoni review did not flag this
   as potentially causal.
2. **Horizontal momentum drift 193× over limit** — possible solver
   bug, threatens Stage 1e anisotropy interpretation.
3. **Layer 3 audit needed before Layer 4 upgrade** — Mechanisms A/E/F
   inherit any Layer 3 issues via γ(φ_p).
4. **Layer 5 audit needed before Mechanism F** — Guo 2017 is single-
   cell; spheroid-level turgor needs cell-cell water transport.
5. **Gate fail taxonomy** — most important non-physics item.
6. **mlsmpm.py monolith is a Mechanism A/E prerequisite** — file split
   needed before adding γ_p / Γ_p state.

---

## My recommendation to PI

**Option F (reviewer-Marangoni hybrid)** with one variation: do
substrate force balance investigation *before* deciding whether to
proceed with Stage 1a++.b or with Mechanism A. Reasoning: if the
anchor balance FAIL is a solver bug, fixing it could shift the
Production Lam4 asymptote upward enough that Stage 1a++.b's
incremental contribution becomes measurable; if it is a physics
limitation, the bug investigation closes the question and Stage 1a++.b
can proceed cleanly.

Concretely:
- **Week 1**: documentation pass (taxonomy, roadmap, parameter
  registry, Marangoni labeling fix)
- **Week 2**: anchor force balance + horizontal momentum drift
  investigation. Outcome: solver bug *or* accepted limitation.
- **Week 3 onward**: PI re-decides between Mechanism A vs Stage 1a++.b
  with cleaner baseline.

PI input required. No code changes pending decision. Auto-STOP.

---

## Action items if Option F selected (with no further input)

These would proceed without PI action since they are documentation
and investigation only:

- [ ] `docs/gate_fail_taxonomy.md` (new)
- [ ] `docs/parameter_registry.md` (new)
- [ ] `docs/10_dev_roadmap.md` refresh through Stage 2
- [ ] `docs/SESSION_HANDOFF.md` consolidation
- [ ] `docs/production_lam4_finding.md` labeling fix (top-down vs
  contact_xy_hull)
- [ ] Anchor force balance root-cause analysis (no code change yet)
- [ ] Horizontal momentum drift root-cause analysis (no code change yet)

The investigation outputs are sanity_md analogs; if they identify
solver bugs, those become separate fix commits with sanity gate.
