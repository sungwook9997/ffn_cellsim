# Session Handoff — 2026-04-29 (Option F Week 1, end-of-day)

This document records where the current session stopped. Keep this
file ≤ 2 pages; for archive material see `docs/SESSION_HANDOFF_archive_v10.md`
(the original handoff written at v10 stage, kept for historical context).

---

## Where we are now

- **Active option**: **Option F** (engineering-Marangoni hybrid) per
  `docs/codex_review_synthesis.md`. PI selected this 2026-04-29.
- **Current week**: Option F Week 1 — documentation pass (in progress).
- **Last code commit**: 6aa3b27 (Codex review synthesis, no code).
- **Last simulation**: Production Lam4 (commit b506b57), result
  reclassified to **Bucket P3 strict** (peak-and-decay, not plateau)
  in this session.

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

## Hard blockers tracked (from gate_fail_taxonomy.md)

| ID | Gate | Severity | Week 2 owner |
|---|---|---|---|
| F2 | horizontal momentum drift (193× over) | blocks Stage 1e anisotropy claim | `horizontal_momentum_drift_investigation.md` |
| F3 | anchor force balance (50× over) | potentially causal for asymptote | `anchor_force_balance_investigation.md` |
| F4 | contact-band ρ_kernel/ρ_ref (15–32% under) | linked to F3 | (rolled into F3 investigation) |
| F9 | φ trajectory vs predicted φ_eq | Layer 3 audit pending | deferred to Week 3+ |

## Week 2 investigations (no code changes)

Both are read-only diagnostic analyses producing sanity_md-style
documents. Outputs determine whether the FAIL is a solver bug
(→ fix commit with sanity gate) or accepted physics
(→ documentation-only update to `docs/12_validation.md` + reclassify
in `docs/gate_fail_taxonomy.md`).

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

1. Resume Option F Week 2: anchor force balance + horizontal momentum
   drift investigations. These are read-only — open the relevant code
   paths in `acs/physics/mlsmpm.py` and `acs/runner.py`, examine the
   FAIL's metric definitions, write the two `_investigation.md` docs.
2. Do NOT start Stage 1a++.b or Mechanism A/E/F implementation. Those
   wait for Week 2 outputs to inform the PI re-decision.
3. Do NOT modify `mlsmpm.py` until Week 2 outputs surface a solver bug
   to fix; the file split (Codex item 9) is Week 3+.
4. Continue Option F discipline: review-only commits, sanity gate
   protocol on any code change, Magic-Number Block on any new constant.

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
