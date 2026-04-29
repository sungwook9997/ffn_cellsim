# Session Handoff — 2026-04-29 (Option F Week 2, end-of-day)

This document records where the current session stopped. Keep this
file ≤ 2 pages; for archive material see `docs/SESSION_HANDOFF_archive_v10.md`
(the original handoff written at v10 stage, kept for historical context).

---

## Where we are now

- **Active option**: **Option F** (engineering-Marangoni hybrid) per
  `docs/codex_review_synthesis.md`. PI selected this 2026-04-29.
- **Current week**: Option F Week 2 — solver investigations DONE.
- **Last code commit**: d93b4ee (Week 1 docs); Week 2 commit pending.
- **Last simulation**: Production Lam4 (commit b506b57), result
  reclassified to **Bucket P3 strict** (peak-and-decay, not plateau).
- **Hard blockers**: 4 → 1 after Week 2. Only F9 (φ trajectory, Layer
  3 audit pending) remains.

## Documentation-pass deliverables (Week 1)

| File | Status | Notes |
|---|---|---|
| `docs/marangoni_review.md` | DONE (commit d0a7d99) | 6 candidate mechanisms A–F, paper framings |
| `docs/codex_review_synthesis.md` | DONE (commit 6aa3b27) | 10-item Codex review, Options E/F/G/H |
| `docs/gate_fail_taxonomy.md` | DONE (this commit) | F1–F9 classification; 4 HARD-BLOCKERs identified |
| `docs/parameter_registry.md` | DONE (this commit) | 4-tier (L1/L2/L3/L4) registry; 11 L3 placeholders |
| `docs/10_dev_roadmap.md` Current Status | REFRESHED (this commit) | now reflects Stage 1a → 2 → Production Lam4 → Option F |
| `docs/production_lam4_finding.md` labeling | FIXED (this commit) | P2 plateau → P3 peak-and-decay |
| `docs/SESSION_HANDOFF.md` | CONSOLIDATED (this file) | original moved to archive |

## Hard blockers tracked (post-Week 2)

| ID | Gate | Status |
|---|---|---|
| ~~F2~~ | horizontal momentum drift | ACCEPTED-LIMITATION (gate too tight in overdamped equilibrium; Stage 1e seed-averaging required) |
| ~~F3~~ | anchor force balance | ACCEPTED-LIMITATION (kernel-truncation diagnostic formula bug; NOT causal) |
| ~~F4~~ | contact-band ρ_kernel/ρ_ref | ACCEPTED-LIMITATION (linked to F3) |
| F9 | φ trajectory vs predicted φ_eq | **REMAINING** — Layer 3 audit pending |

## Week 2 outcomes (this commit)

Both investigations completed read-only:
- `docs/horizontal_momentum_drift_investigation.md` — F2 root-caused
  to gate normalization too tight in overdamped equilibrium + initial-
  pack asymmetry. Absolute drift bounded ~1e-3 across 4 runs.
- `docs/anchor_force_balance_investigation.md` — F3 root-caused to
  Adami-Hu-Adams §3 kernel truncation making F_pressure_down formula
  report negative (tensile) values when contact band is under-densified.
  F_substrate_up is small and stable (~0.038); substrate is
  *under-deflected*, not over-anchoring. Not causal for asymptote.

Three of four hard blockers retired. Marangoni / asymptote
interpretation unaffected (independent root cause per
`docs/marangoni_review.md`).

## Where to find what

- Project mission + framing: `CLAUDE.md`, `docs/00_project_vision.md`
- Stage roadmap + current status: `docs/10_dev_roadmap.md`
- Stage outcome ledgers: `docs/outcomes_*.md` (one per stage)
- Sanity gate analyses: `docs/*_sanity.md` (one per stage / decision)
- Active review docs (Option F Week 1):
  - `docs/marangoni_review.md`
  - `docs/codex_review_synthesis.md`
  - `docs/gate_fail_taxonomy.md`
  - `docs/parameter_registry.md`
- Production Lam4 result + finding:
  - `results/production_lam4/` (artefacts)
  - `docs/production_lam4_finding.md` (revised P3 reading)

## What the next session should do

1. **PI re-decision**: with 3 of 4 hard blockers retired and the
   Marangoni / asymptote diagnosis intact, the choice between Option α
   (Stage 1d.b: Mechanism A/E/F continuum upgrade), Option β (Stage
   1a++.b: discrete boundary events), Option γ (sequenced), or Option
   G (paper-as-is) is now well-posed. Surface to PI.
2. Layer 3 audit (Codex item 4) — only remaining hard blocker (F9). A
   Layer 3 audit unblocks both Mechanism A/E/F (which use γ(φ) and
   thus inherit Layer 3 issues) AND any φ-related claim in the paper.
3. Once a path is selected, prerequisites:
   - **mlsmpm.py file split** (Codex item 9) before adding any new
     state variable (γ_p for Mechanism A, Γ_p for Mechanism E,
     stochastic-event state for Stage 1a++.b)
   - **Symmetry-conservation tests** (Codex item 8 subset) at
     minimum, ideally before file split, to lock in the current
     behavior as a regression baseline
4. Continue Option F discipline: sanity gate protocol on any code
   change, Magic-Number Block on any new constant.

## What is NOT current state (historical context only)

The original `SESSION_HANDOFF.md` (now `SESSION_HANDOFF_archive_v10.md`)
described the v10 era when the curvature operator was failing 192%
over the limit and reference calibration drift was an open root cause.
**Both are resolved**:
- v15 density-based volumetric stress (commit 8e06d8e) closed the
  reference-calibration drift via principled ρ_ref harmonic-mean over
  the well-resolved population matched to F-scale application.
- Production Lam4 (5k particles) achieved curvature first PASS at 9.5%
  rel err (commit b506b57). The Sanity Gate Protocol §6 (measurement-
  protocol consistency) caught the v13 anti-pattern on the way.

The reader of the v10 handoff document should treat it strictly as
project history, not state.
