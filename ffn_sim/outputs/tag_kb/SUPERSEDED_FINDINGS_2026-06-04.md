# Superseded KU-3.5 / γ-floor findings — CANDIDATE list for KB tagging (2026-06-04)

> **STATUS: READ-ONLY CANDIDATE LIST — NOT auto-applied.** This file proposes
> which 2026-05-29 → 2026-06-03 KU-3.5 cortical-tension / γ-floor findings the
> **2026-06-04 authoritative record** supersedes, so the KB (Notion Contract-Graph
> → TAG duckdb / Obsidian) can be tagged `superseded` in a controlled, PI-gated
> pass. **Nothing here edits Notion, the TAG duckdb, or any KB file.** A human /
> PI must review each row before it is marked `superseded` in the SoT and the read
> layers refreshed.
>
> **Authoritative 6/4 record (the overriding source):**
> [`ffn_sim/docs/KU_3_5_GATE_STRUCTURE_2026-06-04.md`](../../docs/KU_3_5_GATE_STRUCTURE_2026-06-04.md)
> + [`ffn_sim/docs/CORTICAL_TENSION_RECORD_2026-06-04.md`](../../docs/CORTICAL_TENSION_RECORD_2026-06-04.md)
> (commits `68995a8` "KU-3.5 gate STRUCTURE split", `e1f66ab` "γ A/B connectivity-insensitive",
> `97f6b3e` "Track-1 force-aggregation").

## The 6/4 overriding points (the "WHY superseded" references)

- **P1** — KU-3.5 is NOT a single "resolved" gate; it splits into **ACTIVE**
  (`g_soft`, still OPEN, below band), **PASSIVE/COMPOSITE** (`g_rigid` + Young-Laplace,
  structural — NOT active credit), and **BRIDGE** (Layer-2 anchor).
- **P2** — `g_rigid` (native ≈ 0.57 mN/m, in-band) is **myosin-independent
  structural/passive** tension and must **NOT** be credited as active KU-3.5.
- **P3** — The active blocker is **force GENERATION / AGGREGATION** (per-head myosin
  force not aggregating into sustained shell tension), **NOT**
  binding / percolation / turnover / coherence / transmission. Connectivity rebuild
  is correct + necessary structure but is NOT the active-γ lever (A/B: z 1.3→3.3
  raised `g_soft` only ~1.5×, did not move `g_total`).
- **P4** — The `--v0-accel`-flag-based "RESOLVED" closure is superseded: the flag
  *was* a real dead-flag bug, but wiring it does **not** reach band — best case
  0.030 mN/m (~11.6× under floor) and only by over-driving per-head force to 2.7×
  stall (non-Hill / super-stall artefact). v0_accel is an un-derived tuning knob;
  chasing band with it violates no-magic-number / no-gate-loosening.
- **P5** — KU-3.5 tension is measurement-protocol / timescale / adhesion-context
  dependent (~3 orders of magnitude across protocols); the band [0.35,0.65] is the
  *slow-micropipette-rounded-cell* oracle, not a single universal number. (Context
  refinement, retained band — listed for completeness; gates the same items above.)

---

## Candidate superseded items

| # | Identifier (dated heading / KB row / commit) | Where it lives (path:line · branch / KB table:row) | Now-superseded claim (one line) | Why superseded (6/4 point) |
|---|---|---|---|---|
| S1 | **STAGE-2 "RESOLVED — dead --v0-accel flag was the floor"** (commit `18ab034`, 2026-06-03) | `ffn_sim/docs/STAGE2_TRANSMISSION_ELIMINATION_2026-06-03.md` (TL;DR) · branch `test/mcf7-fullcell` (NOT on current `layer2/spheroid-cbm`) | "The floor was a dead flag, not a physics wall… wiring `--v0-accel` makes `g_soft` rise 17–150× off the floor → the active channel works → RESOLVED." | **P4** — flag was real but does NOT reach band (0.030 mN/m, 11.6× under, only via 2.7× super-stall over-drive). "RESOLVED" formally superseded by `97f6b3e`; was self-contradicted (same doc body: "extrapolates to ~40 M steps — a real floor"). |
| S2 | **STAGE-2 figure "dead v0-accel flag → g_soft 17-150x climb"** (commit `1b00e3e`) + `fig_stage2_v0accel_breakthrough.png` | `ffn_sim/outputs/h3/figs/fig_stage2_v0accel_breakthrough.png`, `fig_stage2_v0accel_long.png` | The v0-accel "breakthrough" climb closes the active floor. | **P4** — "breakthrough" framing superseded; the climb is bought with super-stall (non-physical) force and never reaches band. The companion `0f62d47` already walked it back ("CH2 over-driving revealed"). |
| S3 | **STAGE-2 elimination chain conclusion** (commit `7e416c1` "transmission/turnover/coherence eliminated → upstream", 2026-06-03) | `ffn_sim/docs/STAGE2_TRANSMISSION_ELIMINATION_2026-06-03.md` (§eliminations) · `test/mcf7-fullcell`; summarized in `outputs/h3/REPORT.md` | Floor blocker isolated by eliminating **percolation, turnover, channel, coherence/meridional, buckling, soft-coupling-ceiling** as levers. | **PARTIALLY superseded by P3**: the *eliminations remain valid* (those are NOT the lever), but the framing that this localizes the cause to "transmission-upstream / soft-coupling ceiling / STAGE-1" is superseded — 6/4 names the cause **force generation/aggregation** specifically (binding/percolation/turnover/coherence explicitly listed as superseded hypotheses). |
| S4 | **CORTICAL_TENSION_TRIAGE §0 + §0.1** "the floor is UPSTREAM of transmission… force GENERATION / soft-coupling CEILING (STAGE-1) / model fidelity" (2026-06-03) | `ffn_sim/docs/CORTICAL_TENSION_TRIAGE_2026-06-03.md:23-77` | Wall narrowed to "force generation / soft-coupling ceiling (STAGE-1) / model fidelity"; binding-fraction hypothesis refuted (Track-1). | **P3** — the binding/percolation/turnover/coherence refutations stand, but the residual attribution to "soft-coupling ceiling / STAGE-1" is superseded by the sharper 6/4 "force generation/AGGREGATION" conclusion. (Doc itself is marked TRIAGE-not-contract; flag, don't delete.) |
| S5 | **`c3fafed` "native γ diagnosis: force-isotropy = the wall"** (2026-06-02) | commit `c3fafed`; logged in `ffn_sim/docs/AUTONOMOUS_LOG_2026-06-02.md:337-341` | "Lead diagnosed **force-isotropy / √N random-dipole cancellation** as the fundamental γ-floor blocker ('the wall')." | **P3** — isotropy/cancellation/coherence is explicitly the *coherence* hypothesis, now superseded; STAGE-2's meridional test (`|axis·meridian|=1.000`) already retired it empirically, and 6/4 confirms the blocker is generation/aggregation, not arrangement. |
| S6 | **AUTONOMOUS_LOG "soft-flag RESOLVED"** (iter 4, 2026-06-02) | `ffn_sim/docs/AUTONOMOUS_LOG_2026-06-02.md:169` | "soft-flag RESOLVED" — implies the γ-soft floor was closed by the production run advancing. | **P4 / P1** — this refers to the same v0-accel/soft-flag line; the run only confirmed warm-up, not a resolved floor (g_tot stayed at the known γ-floor ~4.7e-3 throughout). No "resolved" closure stands. |
| S7 | **KU35_FLOOR_ROOT_CAUSE banner** "the γ floor… is now FULLY DIAGNOSED… It is the myosin contractile mechanism (binned-r0 ratchet lumped proxy)" (2026-05-31) | `ffn_sim/docs/v2_audit/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md:1-11` | The floor was *fully* root-caused to the binned-r0 ratchet (relabel-not-transport); grip-walk fix would lift γ into band. | **P3 / P4** — the binned-r0→grip-walk fix landed (STAGE-1 `5bfb146`) and the floor **persisted**, so the "fully diagnosed / fix closes it" claim is superseded; the remaining blocker is force generation/aggregation, not just the proxy mechanism. (The architectural-principle critique stays valid.) |
| S8 | **REPORT.md 2026-05-29 "Post-R1 re-run… should produce γ_total in/near band"** | `ffn_sim/outputs/h3/REPORT.md:1043-1052` | Exposing `g_rigid` (R1 Lagrange tension) "should produce γ_total in or near [0.35,0.65]" → implies KU-3.5 ✅ candidate. | **P1 / P2** — `g_rigid`/`γ_total` being in-band is **structural/passive**, NOT active KU-3.5 credit; "γ_total in band ⇒ KU-3.5 ✅ candidate" (also REPORT.md:1131 "γ_total ∈ band → KU-3.5 ✅ candidate") is superseded — the active `g_soft` is the gate and it is below band. |
| S9 | **REPORT.md 2026-05-29 sweep "If v3 γ_total ∈ band → KU-3.5 ✅ candidate, H.3 ✅ DONE unblocked"** | `ffn_sim/outputs/h3/REPORT.md:1119-1133` | A single γ_total-in-band number from the v3 sweep would close KU-3.5 and unblock H.3 DONE. | **P1 / P2** — treating any single tension number (esp. g_rigid-dominated γ_total) as "KU-3.5 ✅" is explicitly invalid per 6/4; the active gate stays open. |
| S10 | **DecisionLedger `DEC-2026-05-31-grip-walk`** "Replace binned_r0 proxy with AFINES grip-walk … to address the KU-3.5 tension floor" (PI-ratified) | TAG duckdb `DecisionLedger` row `DEC-2026-05-31-grip-walk` (Notion SoT) | Adopting grip-walk stepping addresses / resolves the KU-3.5 tension floor. | **P3 (qualify, do NOT retract)** — decision was correct + ratified and grip-walk is the right mechanism, BUT the implied "this addresses the floor" is superseded: grip-walk landed and the floor persisted; the floor is force generation/aggregation. Tag as *superseded-as-floor-resolution*, keep as a valid mechanism decision. |

### Items reviewed and judged NOT superseded (kept for auditor context — do NOT tag)

- **`g_rigid` native ≈ 0.57 mN/m in-band** (TRIAGE §0; REPORT) — **retained**, but its
  *role* is reclassified by P2 to passive/composite + Layer-2 bridge (not active credit).
  This is a re-attribution, not a falsification.
- **STAGE-2 individual lever falsifications** (percolation z-sweep, turnover ON/OFF,
  motors ON/OFF channel ×340, meridional |axis·meridian|=1.000, buckling θ≈179°, k-sweep) —
  **the negatives remain valid** (those are genuinely NOT the active lever). Only the
  *overall causal attribution* built on top of them (S3/S4) is superseded.
- **TAG `KnowledgeClaim KB-3.5`** "γ_cortex = active(myosin)+passive(elastic), 0.1–1 mN/m" —
  literature claim, **not** superseded (it is consistent with the 6/4 active/passive split).
- **Chugh 2017 magnitude oracle** (band peak ~0.37 mN/m; T₀=230 pN/µm) and the Layer-2
  bridge papers (Fastabend/Okuda/Roffay) — **literature anchors, not superseded**.
- **L_p / KU-1.1 "4.7σ resolved"** (REPORT.md:233,793) — unrelated to KU-3.5; not in scope.

---

## Provenance notes (for the controlled tagging pass)

- The TAG `run_result` table currently returns **0 rows** for the KU-3.5 γ-floor
  run-results (`tag_query.py` confirmed) — the dated claims are **not yet ingested**
  into the structured KB; they live in commits + docs + REPORT.md + AUTONOMOUS_LOG +
  the figures. So most tagging targets are **doc/commit/figure artefacts**, plus the
  single `DecisionLedger:DEC-2026-05-31-grip-walk` row and `KnowledgeClaim:KB-3.5`
  (KB-3.5 = keep). When the γ-floor RunResults are eventually ingested, S1–S9 should be
  tagged `superseded` at ingest.
- **S1/S2/S3** live on branch `test/mcf7-fullcell`, not the current working branch —
  the STAGE-2 doc is absent from `layer2/spheroid-cbm`. Tag against the commit/branch,
  not a working-tree path.
- Recommended tag value: `superseded` with `superseded_by = KU_3_5_GATE_STRUCTURE_2026-06-04
  (68995a8 / e1f66ab / 97f6b3e)`. PI-gate each row (esp. S10, which is a valid decision
  whose *floor-resolution implication* — not the decision — is what's superseded).

*Generated read-only 2026-06-04. No KB / Notion / duckdb writes performed.*
