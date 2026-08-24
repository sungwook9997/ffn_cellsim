---
kb_record:
  topic: KU-3.5-cortical-tension
  claim: KB-3.5
  gate: VG-H3-KU35-cortex-tension
  status: authoritative
  authoritative_as_of: 2026-06-04
  aliases: [KU-3.5, KU3.5, cortical tension, g_soft, g_rigid]
  supersedes:
    - KU35_FLOOR_ROOT_CAUSE_2026-05-31
    - CORTICAL_TENSION_TRIAGE_2026-06-03
    - STAGE2-v0accel-RESOLVED@18ab034
  conclusion: >-
    KU-3.5 is NOT a single resolved-bug gate. It splits into KU-3.5-ACTIVE
    (g_soft; blocker = force generation/aggregation; still OPEN), KU-3.5-PASSIVE
    (g_rigid in-band but myosin-INDEPENDENT structural tension, plus the turgor
    Pi0 physiological baseline), and KU-3.5-BRIDGE (g_rigid/gamma as the Layer-2
    spheroid surface-tension anchor). The 5/31-6/3 "KU-3.5 resolved" and
    --v0-accel-based conclusions are SUPERSEDED; do not let them override the 6/4
    record.
---

# KU-3.5 cortical tension — gate structure (2026-06-04 AUTHORITATIVE)

**PI-ratified 2026-06-04.** KU-3.5 is NOT a single "resolved bug" gate. The
2026-06-04 record SUPERSEDES the 5/31–6/3 "KU-3.5 resolved" conclusion and splits
KU-3.5 into an **active** gate, a **passive/composite** gate, and a **Layer-2
bridge** anchor. Treating any single tension number as "KU-3.5 ✅" is invalid.

## Why split (the 2026-06-04 findings)

1. **The `--v0-accel` dead flag was a real bug** (verified: `run_arm` threads
   `v0_accel` but `_build` never applies it to `p_myo` — myosin runs at literal
   NMIIA v0). BUT the old "KU-3.5 resolved" conclusion built on it is
   **SUPERSEDED**, not the answer.
2. **`g_rigid` (constrained-backbone Lagrange tension), native ≈ 0.57 mN/m,
   is IN-BAND but myosin-INDEPENDENT structural/passive tension** — it must NOT
   be credited as active KU-3.5.
3. **`g_soft` is the active myosin readout**, and in the Hill-valid /
   F_stall-valid regime (F/F_stall ≤ 1) it is still BELOW band.
4. **The current blocker is force GENERATION / AGGREGATION, not
   binding/percolation/turnover/coherence**: per-head myosin force does not
   aggregate into a sustained shell tension.
   - Cross-check (2026-06-04 connectivity A/B, `figs/connected_gamma_ab.png`):
     connectivity rebuild (z 1.3→3.3) raised `g_soft` only ~1.5× (3.0e-4 vs
     2.3e-4) and did NOT move `g_total` (g_rigid-dominated) → connectivity is
     NOT the active blocker; aggregation is.

## The three gates + anchor

### KU-3.5-ACTIVE  (the real, still-OPEN gate)
- **Readout**: `g_soft` (compliant crosslink-channel tension = active myosin).
- **Validity window**: motors ON−OFF *delta* (subtract the passive baseline),
  F/F_stall ≤ 1 ONLY (Hill-valid / F_stall-valid; reject super-stall artefacts).
- **Band**: KU-3.5 [0.35, 0.65] mN/m as the *active myosin* contribution.
- **Status**: BELOW band. **Blocker = force generation/aggregation** (per-head →
  sustained shell tension). Real targets: (a) `--v0-accel` dead flag (myosin
  under-drive), (b) per-head force → shell-stress aggregation (why the sum
  doesn't sustain), (c) Hill/F_stall validity of the operating point.

### KU-3.5-PASSIVE / COMPOSITE  (structural, NOT active credit)
- **Readout**: `g_rigid` (constrained-backbone Lagrange tension) + membrane /
  enclosed-volume Young-Laplace (γ = ΔP·R/2).
- **Physiological baseline (CLAUDE.md HARD rule)**: the resting cell carries a
  turgor Π₀ that PRE-TENSIONS the cortex. The enclosed-volume law
  `ΔP = Π₀ − K_vol·(V−V0)/V0` has `turgor_dP0 = Π₀` defaulting to **0** — an
  unpressurised floppy shell (the structural floor). Set Π₀ to the resting
  turgor (interphase ≈ 40 Pa Stewart 2011 → γ=0.15; band [0.35,0.65] ⇔ Π₀ ≈
  93–173 Pa). Driver: `mcf7_fullcell_stage1 --turgor-pa <Π₀>`
  (sets `turgor_dP0`, 2026-06-04). This gate validates the PASSIVE/structural
  tension only — it is explicitly NOT active myosin credit.

### KU-3.5-BRIDGE  (Layer-2 anchor)
- Use `g_rigid` / γ as the **spheroid surface-tension anchor** for the Layer-2
  spheroid/cohesion line (the composite cell-surface tension a multicellular
  aggregate presents), NOT as single-cell active KU-3.5.

## Superseded findings (TAG/RAG must tag, do NOT let 5/31–6/3 override 6/4)
- "KU-3.5 resolved" (single-bug closure) — **SUPERSEDED**.
- 5/31–6/3 hypotheses that the blocker is binding / percolation / turnover /
  coherence — **SUPERSEDED** by the 6/4 conclusion (blocker = generation/
  aggregation). The connectivity rebuild is correct + necessary structure but is
  NOT the active-γ lever.
- The `--v0-accel`-based "resolved" conclusion — **SUPERSEDED** (flag was dead;
  fixing it alone does not close KU-3.5-active).

## Aggregation-estimator result (2026-06-04) — generation, NOT aggregation

The force-AGGREGATION estimator is built + fixture-validated
(`scripts/h3_ku35_aggregation.py`, `tests/test_ku35_aggregation.py`,
`--aggregation` flag on `mcf7_fullcell_stage1`). It decomposes the chain
*per-head generation → coherent ceiling → realized g_soft* into orthogonal
efficiencies (η_agg = radial-projection loss; η_medium = xlink/spring
cancellation), built on the project's IK-virial 8πR² calibration (corrected MOP
≈ IK; the estimator audit's convention). It emits one of three verdicts:
GENERATION-LIMITED / AGGREGATION-LIMITED(RADIAL) / AGGREGATION-LIMITED(MEDIUM).

**Verdict: GENERATION-LIMITED — robust at dev (n_fil=400, CPU) AND production
scale (n_fil=1000, n_motors=200, gbook GPU, all compartments ON;
`outputs/h3/production/ku35_active/AGG_n{400_cpu,1000_gpu}.log`,
`figs/ku35_aggregation_motorsON.png`).** n_fil=1000 plateau (sample 6):
n_engaged 1242, Σ|F_head| 2.43 nN, η_agg 0.84, η_medium 1.05, η_total 0.88,
ceiling/band 0.0010, gate g_soft 2.90e-4 ≈ mop_attach 2.49e-4 ≈ ik_attach
3.01e-4 ≈ ceiling 3.56e-4 mN/m (all within ~1.4×). Identical verdict to
n_fil=400 (ceiling/band 0.0012) → scale-invariant. **Motor count verified
physiological + self-consistent**: native_n_motors = 2120.6 (Salbreux 3/µm² ×
706.9 µm²); n_motors=200 × factor 10.60 = native (invariant); total stall
capacity 21.21 nN is n_motors-invariant (≈0.45 mN/m over 2πR ⇒ the full native
complement at stall IS band-capable — the floor is the sub-stall operating point
+ short lever, NOT motor undercount); the estimator reads scaled k/F_stall from
the live run so the historical param-map mis-pin is impossible. **Connected mesh
verified**: 0 same-filament staples, giant component 1.000, z 3.65 at n_fil=1000.

- **Aggregation is NOT the wall.** η_agg ≈ 0.91, η_medium ≈ 1.02, η_total ≈ 0.94
  — the per-head forces aggregate ~94% efficiently. Cross-check: gate g_soft
  ≈ internal MOP ≈ IK ≈ coherent ceiling, all within ~10% (sample 4: gate
  3.99e-4, mop 3.92e-4, ik_attach 3.81e-4, ceiling 4.16e-4 mN/m). The connected
  mesh + Arp2/3 structure does its job; the realized g_soft ≈ what the attach
  forces *can* produce. **No hidden leak/cancellation** — what is generated is
  what is measured. (Refines the earlier "doesn't aggregate" wording: it DOES
  aggregate, faithfully, into a correspondingly small tension.)
- **The wall is GENERATION.** ceiling/band ≈ 0.0012 → the coherent ceiling of
  the attach channel is itself ~840× below band. Two compounding deficits:
  (i) raw force Σ|F_head|/band ≈ 0.15 (per-head 2.37 pN, deeply sub-stall —
  series F/F_stall ≈ 0.005, 100% Hill-valid; ~1038/4000 ≈ 26% heads recruited);
  (ii) short lever — the soft channel only sees the ~240 nm head→bead attach
  bond on a 7.5 µm shell (the rest of the motor reaction is absorbed by the
  rigid M-SHAKE backbone → the motor-insensitive g_rigid). Recruitment grew
  586→1038 but is slowing; ~4× more heads ≈ ~4× ceiling = still ~200× under.
- **Levers (active gate)**: per-head delivered force/lever + recruitment +
  whether the active stress needs a soft long-range transmission path (the
  rigid backbone shunts it locally). NOT aggregation geometry, NOT medium
  cancellation (both ~1).

## Action items
- [ ] Tag the above as `superseded` in RAG/TAG (RunResult/DecisionLedger/
      KnowledgeClaim) so 5/31–6/3 does not override the 6/4 record.
- [x] KU-3.5-active diagnosis: force generation/aggregation (per-head → shell
      tension) — **aggregation estimator built + validated; verdict =
      GENERATION-LIMITED (aggregation ~94% efficient, generation ~840× under)**.
- [x] KU-3.5-active follow-up: production-scale (n_fil=1000) `--aggregation` on
      gbook GPU — **DONE, verdict GENERATION-LIMITED confirmed scale-invariant**
      (ceiling/band 0.0010 @ n1000 vs 0.0012 @ n400). NEXT: quantify the
      per-head-force / recruitment / soft-transmission levers (PI-gated core
      myosin per the 6/4 CORTICAL_TENSION_RECORD).
- [ ] KU-3.5-passive: wire `turgor_dP0` = resting Π₀ as the production baseline.
- [ ] KU-3.5-bridge: expose g_rigid/γ as the Layer-2 spheroid surface-tension
      anchor.
