# DCM spreading investigation — plan (PI pivot: prioritise spreading over active assembly)

**Date:** 2026-07-10 · **Owner:** Lead session (8 h autonomous /goal) · **Branch:** `dcm/main`
**Directive:** PI — "스페로이드 조립을 너무 액티브하게 하는 게 잘 안 되면 퍼지는 거를 더욱이 우선적으로." Every step MUST
(1) compare against TAG/lit data and (2) be VISUALLY verified (HTML viewer / figure, checked in a real browser).

---

## 0. Why the pivot (two confirmed negatives on active assembly)

- **S4** (`DCM_CADHERIN_CLUSTER_REARRANGEMENT_DESIGN_2026-07-09` §6c): the aggregate-σ compaction RATE is NOT
  τ_mature-gated (τ=30 s ≈ τ=600 s). The radial liquid-drop densifies by drag-limited gap-closing, bypassing the junction.
- **(b)S1** (`DCM_T1_REARRANGEMENT_COMPACTION_DESIGN_2026-07-10` §6c): the tangential DAH (differential-adhesion) driver
  does NOT compact a confluent aggregate — stable run (CFL 0.24, gates PASS) porosity **0.302 → 0.301 flat**. Differential
  tension alone cannot densify (confirms the earlier §2e finding at native scale).

Active spheroid *compaction* is not junction-rate-gatable with the current drivers (aggregate-σ or DAH). Per PI, pivot to
**spreading**. The (b)S2 active-T1-contraction option is documented + deferred; the S1–S3 cadherin cluster+maturation
mechanism is CORRECT and is REUSED below (it is the fine-grained de-cohesion lever spreading needs).

## 1. What the DATA (TAG) says — the physiological target

TAG corpus (queried 2026-07-10; Douezan/Beaune/González wetting papers are NOT in the KB — flagged for ingest):
- **Warmt2021_NewJPhys** (the one quantitative breast datum): aggregates on 1.5 g/L collagen, 60 h — **MDA-MB-231 covered
  area ×8, MDA-MB-436 ×3 (mesenchymal = strong wetting/spread); MCF-10A stays in the aggregate = NON-wetting (epithelial).**
- **Mangani2025**: MCF-7 spheroids stay confined (limited dissemination) vs MDA-MB-231 rapid disseminate.
- **Chaudhuri2015** (U2OS): single-cell spread area ≈ **750–1500 µm²**, rises with stiffness (1.4→9 kPa) + RGD ligand.
- ⚠️ **Moazzeni2021 δ=a/a0−1 is ELECTRODEFORMATION shape, NOT substrate spreading** — do not conflate.

**Physiological verdict**: epithelial (MCF-7 / MCF-10A) = **non-wetting / confined** (a spheroid does NOT spread — matches
the prior DCM structural conclusion + Douezan S<0 at MCF7 γ=1e-2); **mesenchymal (MDA-MB-231) = spreading via single-cell
dispersal / crawl-out** (×8 area). So "spreading" in this system is **NOT bulk flattening** — it is **de-cohesion + single
cells crawling out** (an EMT/invasion phenotype), and the cell-type difference is set by **junction adhesion** (which the
S1–S3 cadherin cluster+maturation models directly).

## 2. What the CODE says — machinery + prior limits (spreading recon 2026-07-10)

- **Aggregate bulk flattening = confirmed structural DEAD-END** (`DCM_SPREADING_LIMIT_CONCLUSION_2026-07-01`): all 5
  mechanisms give A/A0 = 1.000–1.003. Reasons: FA clutch is an ANCHOR not a motor (no outward force); protrusion reaches
  only ~3 % rim cells (11/400 basal); turgor (V/V0=1) + cohesion resist flattening.
- **Single-cell spreading WORKS**: A/A0 → 3.8 (fried-egg) at Douezan S>0; but at physiological γ=1e-2, S<0 → even a single
  cell needs the ACTIVE driver (lamellipodium/`spreading_push`), stiff volume conservation (k_vol), and a substrate FLOOR
  (z-well) to bound it (else area-maximising node-ejection).
- **Faithful aggregate route = basal DE-COHESION → single-cell crawl-out** — needs (i) a working IPC contact (per-face
  penalty tunnels under lit cadherin bundles; **the finished #1 IPC-newton now provides this**), (ii) per-cell active
  protrusion (`all_cell_protrusion` / cryptic followers), (iii) de-cohesion (cadherin catch-bonds rupturing under traction
  — **the S1–S3 cluster gives this, fine-grained**).
- Machinery available: `dcm_substrate_warp` (z-well anchor, u-bottom ULA, wetting drive, compliant/mechano);
  `dcm_lamellipodium_ratchet_warp` (Brownian-ratchet rim protrusion + tether clutch, `--lamellipodium --lamel-clutch`);
  `dcm_ecm_clutch_host` (Pereverzev catch-slip FA anchor); A/A0 = `_topdown_area_um2` (top-down silhouette, PI rule);
  `spread_eval.py` (GENUINE spread vs PEELING vs no-spread, peeling_index guard); production driver `run_spread_overnight.py`.

## 3. Direction — staged, in the DCM lane (FF is the parallel session's)

**Stage 1 — DCM single-cell active spreading (fast, decisive, bounded).** Establish the working single-cell spread with the
ACTIVE lamellipodium driver at physiological conditions (γ=1e-2, stiff k_vol, z-well floor, FA catch-slip clutch). Reproduce
A/A0 > 1 (top-down, peeling-aware) and compare the spread AREA to Chaudhuri (750–1500 µm²) + the fried-egg A/A0→3.8. This
is the honest single-cell phenomenon and iterates fast (1 cell).

**Stage 2 — DCM aggregate spreading via de-cohesion + crawl-out (the PI aggregate target).** Substrate-adhered spheroid +
finished IPC-newton + S1–S3 cadherin cluster (de-cohesion) + per-cell active protrusion (`all_cell_protrusion`). Question:
do peripheral cells DE-COHERE and CRAWL OUT (A/A0 up via dispersal, peeling-aware — GENUINE, not rim-peel), and does the
extent scale with cadherin adhesion/maturation (strong/matured = confined = epithelial/MCF-7; weak/nascent = spreading =
mesenchymal/MDA-MB-231, matching Warmt2021 ×8 vs non-wetting)? This reuses ALL of S1–S3 and turns the "correct but
compaction-irrelevant" junction mechanism into the driver of a data-matching spreading phenotype.

## 4. Gates (write before running; do not loosen; TAG-anchored)

- **G-S1a — single cell spreads**: A/A0 (top-down) > 1.10 with maxZ drop > 3 % (GENUINE per `spread_eval`, peeling_index < 8),
  stable (finite, V/V0≈1, pen low, CFL < 1). Spread area in the Chaudhuri 750–1500 µm² order for a single MCF7-sized cell.
- **G-S1b — active-driven, bounded**: spreading vanishes with the lamellipodium OFF (it is the driver, not passive), and does
  NOT node-eject (z-well floor + k_vol bound it).
- **G-S2a — aggregate de-cohesion spreads by crawl-out**: A/A0 rises via peripheral single-cell dispersal (GENUINE, not
  bulk flatten, not rim-peel), with actual cadherin bond rupture at the periphery.
- **G-S2b — junction sets the phenotype (the data match)**: strong/matured cadherin → confined (A/A0≈1, epithelial/MCF-7);
  weak/nascent → dispersing (A/A0 up, mesenchymal) — reproducing Warmt2021's epithelial-non-wetting / mesenchymal-×8 split.
  Cadherin params at LIT values (no tuning to an A/A0 target).
- **G-stability — every run**: finite / G2 interpenetration / volume-conservation gates PASS; CFL < 1; sim does not blow up.
- **G-visual — every milestone**: interactive HTML cell-shape viewer (peel/cut/top-down), browser-verified (real Chrome).
- **G-TAG — every result**: compared against the §1 lit values before any claim.

## 5. Process (PI mandate — bake in)

1. **TAG compare** every result (the §1 data; query tag_query.py for any new quantity; confirm SourceEvidence verdict before citing).
2. **Visual verify** every milestone myself — generate the HTML cell-shape viewer + browser_check.py (screenshot ground truth),
   and READ the figure/screenshot. No claim from grep/log alone.
3. **Stability gate** every run (no blow-up). **Backups**: commit frequently + push the dated backup branch.
4. **Continuously replan**: update this doc + write new plans as directions emerge. After 4 h, re-query TAG for newly-ingested
   spreading data → refine the plan on it.

## 6. Risks / open questions

1. Single-cell spread without a contact-line bound node-ejects (`PHASE_C_SPREADING_LEVERS`) — the z-well floor + stiff k_vol
   is the bound; watch node-ejection (maxZ, V/V0, CFL).
2. Aggregate de-cohesion needs the IPC contact to hold under lit bundles — the finished IPC-newton should; verify pen stays low.
3. Douezan/Warmt values are sparse in TAG — ingest the wetting papers (Douezan 2011, Beaune 2014, Warmt2021 is present) as
   SourceEvidence so the comparison is KB-grounded (a TAG-ingest task, PI-gated for new SE rows).
4. Do NOT tune cadherin/γ to hit an A/A0 target (magic-number rule) — the phenotype must EMERGE from lit adhesion values.

## 7. One-line summary
Active spheroid compaction is a confirmed dead-end (S4 + (b)S1); pivot to spreading — establish the working single-cell active
spread (Stage 1, TAG-compared to Chaudhuri/fried-egg), then turn the S1–S3 cadherin cluster into the fine-grained de-cohesion
driver of **aggregate spreading by single-cell crawl-out** (Stage 2), where junction adhesion/maturation sets the
epithelial-confined (MCF-7) vs mesenchymal-spreading (MDA-MB-231, Warmt2021 ×8) phenotype — every step TAG-compared + visually
verified + stability-gated.
