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
