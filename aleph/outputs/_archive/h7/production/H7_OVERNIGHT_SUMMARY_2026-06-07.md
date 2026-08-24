# H.7 overnight autonomous run — morning summary (2026-06-07)

Per PI /goal: ran production + worked autonomously past KU-3.5, parallel subagents, viz +
theoretical/novelty analysis. gbook RTX A5000, full ×40 scale. All committed on
`h7/full-cell-integration`; gbook code backup `~/ffn_cellsim_2026-06-07_preH7_codebackup.tgz`.

## Headline (honest, PROVISIONAL — KU-3.5 conclusion is a PI call; NO gate-loosening)
The full physiological MCF7 cell assembles + runs at full ×40 scale, and its EMERGENT cortical
tension **floors ~6-10× below the [0.35,0.65] mN/m band**. No clean emergent channel reaches it:
- **γ_active (myosin) = 0.0001 ± 0 mN/m** — n=8 seeds, constrained AND unconstrained, turnover ON
  AND OFF. Robust, scale-independent (×40 converged).
- **γ_rigid (structural, free pressurised cortex) = 0.06-0.09 mN/m** — n=3 + ×40-converged; turgor
  IMPLIES 0.50 but transmits only ~0.06 to the backbone (a transmission gap). Below band.
- **γ_passive = 0.499 mN/m** — circular identity with the B3 turgor setpoint, not a measurement.
- The one ever-in-band number (0.391) was an FA integrin-overload ARTIFACT (r0=0 bug), discarded.

**Root cause = a TIMESCALE GAP**: tension-generating processes are seconds-scale (myosin
contraction; turnover τ½=10 s) while feasible fine-grained constrained-MD reaches only ~ms — the
active-tension steady state is unreachable. Directly confirmed (turnover-ON didn't help) +
cross-checked (active~0 in both integrator paths → real generation limit, not an M-SHAKE artifact).
This is the active-gel continuum's seconds-scale advantage and the novelty-analysis's sharpest
risk ("emergent γ is a promise, not a result"), now quantified in the complete cell. Details:
`H7_GATEB_FINDINGS_2026-06-07.md`.

## What was delivered (all committed)
- **Full-cell GATE-B pipeline** (`scripts/h7_gate_b.py`): FA-adhered build → equilibration →
  3-channel γ (`cortex/cortical_tension.py`, turgor never folded in). Levers: --softstart,
  --dt-safety, --no-constrained, --turnover, --no-fa.
- **6 GPU production rounds** (rounds 1-6): γ ensembles, turnover test, clean γ_rigid, ×40
  convergence. Summaries in `outputs/h7/production/H7_*.md`.
- **Both single-cell lamellipodium geometries** (basal_ring isotropic, polarized_patch migrating)
  build in the full cell + A/A0 comparison (`h7_spreading_compare.py`) + **cell-movement animation**
  `outputs/h7/figs/h7_spreading_movement.gif` (ring +24% isotropic, patch directional).
- **3 research reviews** (CrossRef/Consensus-verified citations): direction-review,
  novelty-analysis (novelty = integration-level + the 3-channel γ discipline; narrow),
  paper-scout (16 SE candidates, NOT auto-registered).
- **Stability map + fa.py fix spec** (`H7_FA_INTEGRIN_OVERLOAD_FIX`): unconstrained-soft robust;
  constrained+FA integrin-overload-fragile (force-free-integrin fix spec'd, NOT implemented — see below).

## PI DECISIONS NEEDED (I did not force any of these)
1. **Path forward for active-γ** (the core question): (a) accelerated dynamics — raise myosin/
   turnover rates to reach steady state in feasible steps, BUT breaks Hill-validity (read γ at
   F/F_stall ≤ 1) → needs an explicit validity gate; I deliberately did NOT run this overnight
   because it would produce a Hill-invalid in-band number that could be misread as success.
   (b) multiscale active-gel γ-seam (hand coarse stress to a seconds-scale continuum). (c) re-target:
   a spread adherent MCF7 is traction/stress-fibre dominated → maybe traction, not cortical γ.
2. **Lamellipodium patch-vs-rim** — both built; pick by matching your A/A0 data (still open).
3. **Band target** — [0.35,0.65] is a rounded/de-adhered non-MCF7 proxy; confirm the right MCF7 target.
4. **De-circularize γ_passive** — independent osmotic datum, or keep labeling it a setpoint.
5. **fa.py fix** — implement the force-free-integrin fix (delicate multi-file core-FA change; spec'd,
   test-gated) to get a clean FA-adhered γ_rigid? (Not blocking — clean γ_rigid already via --no-fa.)
6. **Citation fix** — `validation/pereverzev.py:85-87` KU-2.18 has the wrong title for Bangasser 2013
   (metadata drift, like the 2026-06-02 audit) — recommend an SE title fix.

## Commits (this session, on h7/full-cell-integration)
Phase A (a50ee24→eb6c23a) → viz/Codex/B3 → B4/KU-5.x → triple-head geometries (2b89144) → Track-1
A/A0 (fb3bc1e) → gbook unblock + 6 production rounds → consolidated findings (bdf53d6) + this summary.
