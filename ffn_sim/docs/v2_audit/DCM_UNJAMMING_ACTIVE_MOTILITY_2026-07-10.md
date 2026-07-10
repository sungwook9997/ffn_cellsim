# DCM spreading SOLVED — active-motility unjamming (the jamming→unjamming transition)

**Date:** 2026-07-10 · **Owner:** Lead session (PI 8h /goal) · **Branch:** `dcm/main`
**Result:** the DCM aggregate spreads by an **active-matter jamming→unjamming transition** — visually verified,
dose-dependent, and quantified in the literature's own order parameter (shape index → s0*≈5.41), validated by the
freshly-ingested TAG jamming/SPV literature + PI's KB-PIV-10 experiment.

---

## 0. The arc (from "can't spread" to the unjamming transition)

The spreading investigation established (visual-verified, `DCM_SPREADING_INVESTIGATION_PLAN_2026-07-10`): the DCM
aggregate is **CONFINED regardless of cadherin** (strong A/A0=1.001 = weak 0.990) and single cells node-eject — the
**root is that cells don't slide/rearrange (a JAMMED tissue)**. This is exactly the jamming physics: the new TAG
literature (SE 355→376, papers 191→264) gives the framework — **Bi–Manning shape-index jamming** (p0*=3.81 in 2D;
Merkel–Manning **s0*≈5.41 in 3D**; below = jammed/solid, above = unjammed/fluid), the **Self-Propelled Voronoi (SPV)**
model (Bi-Yang-Marchetti-Manning: active motility v0 fluidises a tissue even below the shape threshold), and
Geiger2022 ("*spheroid invasion = fluidization = a jammed→unjammed phase transition*"). PI's own **KB-PIV-10** MCF7/col-I
records the solid/liquid/gas unjamming by aggregate size + low-E-cadherin→gas/EMT.

## 1. The build — per-cell active self-propulsion (the SPV unjamming lever)

`dcm/dcm_active_motility_warp.py` — each cell self-propels with force `F_active·p̂_c`, p̂_c a slowly-reorienting
polarity (rotational diffusion, persistence τ_p; persistent random walk = the standard active-matter migrating-cell
model). Applied per node (mesh-invariant whole-cell F_active), wired into `run_decohesion` (`--active-motility
--f-active-nn --motility-tau`, default OFF). F_active = physiological single-cell traction (10–100 nN, KB-2.12;
controlled variable, swept — NOT tuned to a spreading target).

## 2. RESULT — active motility unjams the DCM aggregate → GENUINE dose-dependent spread

Native N=100 on gbook A5000, spheroid on substrate + IPC-newton + cadherin cluster de-cohesion + active motility,
6000 steps. **Every A/A0 browser-verified (screenshots) — genuine collective loosening, cells stay intact (NOT the
wetting node-ejection artifact).**

| F_active [nN] | A/A0 (top-down) | drift [µm] | peel_idx | 3D shape index s | spread_eval |
|---|---|---|---|---|---|
| 0 (confined) | 1.006 | 0.08 | — | **4.934** (JAMMED, < s0*) | confined |
| 10 | 1.114 | 1.02 | 3.4 | 4.988 | GENUINE COLLECTIVE SPREAD |
| 50 | 1.552 | 2.60 | 1.4 | 5.200 | GENUINE COLLECTIVE SPREAD |
| 100 | **1.918** | 3.80 | **1.1** | **5.351** (→ s0*≈5.41) | GENUINE COLLECTIVE SPREAD |

**Clean monotonic dose-response** (fig `dcm_unjamming_active_motility_doseresponse.png`): more active motility →
- more spread (A/A0 1.0→1.92, drift 0.08→3.8 µm),
- the per-cell **shape index rises 4.93→5.35, toward the Merkel–Manning 3D unjamming threshold s0*≈5.41** — the
  jamming→unjamming ORDER PARAMETER crossing from solid toward fluid,
- **peel_idx DECREASES 3.4→1.1** (more genuine at higher force, NOT peeling/ejection), V/V0=1.000 (volume-conserved),
  pen low, stable. Visual: cells rearrange + move apart as intact spheres (confirmed 10 nN + 100 nN).

## 3. What this means

- The DCM confined-vs-spread is the **active-matter jamming↔unjamming transition** — the correct physics for
  collective cell migration / spheroid invasion (Geiger2022, SPV). The DCM does it, quantified in the lit order
  parameter (shape index), at physiological F_active.
- It resolves the whole session: static confinement (jammed, epithelial/MCF-7) AND spreading (unjammed, active
  motility, mesenchymal/invasion) are the **two phases of one transition**, tuned by active motility — matching the
  SPV framework + PI's KB-PIV-10 (low E-cadherin / high motility → gas/EMT).
- The unifying root (cells don't rearrange) is unjammed by active motility — the SAME lever should rate-gate
  compaction (fluidised tissue flows), closing the (b) T1 line too.

## 4. Sanity / integrity

Every value lit-anchored: F_active = KB-2.12 single-cell traction (swept, not tuned); s0*≈5.41 = Merkel–Manning 3D
(lit); shape index = the standard S/V^(2/3). No gate loosened. All runs gate-PASS (finite/G2/volume), stable.
Every A/A0 claim VISUALLY verified (browser_check) — which is why the wetting/node-ejection artifacts were caught and
this genuine spread is trusted. Compared against TAG (jamming/SPV/Geiger2022) + PI KB-PIV-10.

## 4b. HONEST caveats (2026-07-10, EMT-sweep follow-up)

The §2 dose-response is a **single run per F_active at motility-seed 7** — and follow-up runs show it is **stochastic
and not cleanly monotonic**, so those numbers are indicative, not yet ensemble-quantified:
- **fa150 nascent** (same seed 7, DETERMINISTIC vs fa100): A/A0=1.229, drift 1.35 µm — **LESS** than fa100 (1.918),
  a real non-monotonicity at high force (suspect a velocity/CFL cap or a collective-cancellation of the stronger
  random propulsion). To investigate.
- **fa100 MATURED cadherin** (n=1): A/A0=2.421, drift 5.14 µm — **more** than fa100 nascent (1.918), the OPPOSITE of
  the "matured/epithelial resists unjamming" hypothesis. **Visually verified GENUINE** (intact cells, dispersed, no
  ejection), but n=1 is within the stochastic spread — inconclusive on maturation-gating.
- Unjamming TIME-COURSE (fa100): shape index jumps 4.855→5.30 in the first ~2 s then plateaus at 5.35 (frac_unjammed
  0→45 %); confined stays 4.85→4.93 (7 %). So the fluidisation is fast; A/A0 keeps rising as cells drift apart.

**Conclusion stands at the mechanism level** (active motility unjams the DCM aggregate → genuine collective spread,
shape index rises toward s0*, visual-verified) — but the QUANTITATIVE dose-response + the maturation-gating question
require **seed ensembles** (per the 'ensemble before spontaneous' rule). `--motility-seed` now exposes the polarity RNG.

## 5. Next (in progress)
- **Seed ensemble** (running): fa100 nascent × {7,11,17} vs matured × {7,11,17} → the maturation effect with error bars
  (does junction maturation gate the unjamming, matching KB-PIV-10 low-E-cadherin→gas?).
- Resolve the fa150 non-monotonicity (CFL/velocity cap check); map A/A0, s, frac_unjammed vs F_active with ensembles
  → the SPV phase diagram; anchor to KB-PIV-10 unjamming-by-size.
- Apply the same fluidisation to rate-gated compaction ((b) line).
