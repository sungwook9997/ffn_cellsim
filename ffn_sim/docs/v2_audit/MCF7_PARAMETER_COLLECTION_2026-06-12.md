# MCF7 physiological parameter collection — DCM re-anchoring (2026-06-12)

Branch `h7/compartment-platform`. Collected via an **Ultracode multi-agent workflow** (17 agents:
8 parameter dimensions × search → adversarial verify → synthesis; 283 tool-uses over web +
PubMed + bioRxiv + Consensus). Every value is from a verified source with a checked DOI; PROXY /
UNVERIFIED provenance is flagged inline. Where no MCF7-specific datum exists, that is stated — no
value is invented (PI param rule, 2026-06-12 HARD). This doc anchors the node-face contact re-tune
(brief §1) and Phase 3 (γ + area). **PI decisions in §4 gate the re-tune.**

## 0. Headlines (read first)

- ⭐ **The platform γ = 2.7×10⁻⁴ N/m is NOT a measured MCF7 value.** It is a downstream SimuCell3D
  pipeline config number *below that paper's own generic epithelial floor* (γ_default = 1×10⁻³,
  range 0.5–2.5×10⁻³). The **only direct MCF7 single-cell cortical tension is ~1×10⁻² N/m
  (~10 mN/m)** (Moazzeni 2021, electrodeformation, n=47, SUSPENDED). ⇒ a **~40× upward re-anchor**
  (PI decision, §2/§4).
- ⭐ **The Nagle 21–45 mN/m band is aggregate/tissue surface tension on MCF10A/MCF10DCIS — NOT MCF7,
  NOT single-cell.** It must not be equated with cortical tension. The 0.27 / 0.35–0.65 mN/m
  proxies all sit in the weak micropipette-aspiration window (MPA under-reports — suction remodels
  the cortex), ~20–40× below the real MCF7 value.
- ⭐ **Contact re-tune confirmed quantitatively**: repulsion ξ̄ = 2.42 vs target 0.48 ⇒ **ξ ≈ 4×10⁷
  Pa/m (5× softer than current 2e8)**; smeared adhesion ω̄ = 9.7 vs the mcf7_p0.xml design ω = 0
  (**cohesion belongs in a discrete cadherin module, not a smeared ω**).
- ⚠ **Codebase citation fixes**: (a) cytoplasm viscosity 65.9 Pa·s is **Dessard, Manneville &
  Berret 2024** (Nanoscale Adv 6:1727, 10.1039/d4na00003j) — the "Hu 2024" attribution is WRONG
  (PMC ID was right, value 65.9±11.4 confirmed). (b) SimuCell3D = **Runser, Vetter, Iber 2024**.

## 1. MCF7 parameter table

SI units. "MCF7?" = the *measurement* was on MCF-7 (vs a proxy line / model default).

| Parameter | Best value (SI) | Band | Condition | Source (DOI) | MCF7? | Verdict |
|---|---|---|---|---|---|---|
| **Cortical/surface tension γ** | **~1×10⁻² N/m** (10 mN/m) @0.01 s naive cortex | 10⁻²–10⁻¹ N/m | SUSPENDED rounded cell, electrodeformation | Moazzeni 2021, PRE 103:032409, 10.1103/PhysRevE.103.032409 | **Y** (only direct) | CONFIRMED |
| — γ adherent operating point | **no direct MCF7 datum** | proxies 0.27–45 mN/m (§2) | adherent/spread | — | N | OPEN |
| **Cell–cell adhesion** (whole-pair force) | ~1–2 nN @5 s → ~4–5 nN @120 s | rises w/ contact time | MCF7–MCF7 pair, E-cad⁺ | Hyler 2020, Cells, PMC7227807 | **Y** | PLAUSIBLE (graph) |
| — ranking | MCF-7 highest of MCF-7>T47D>MDA-MB-231 | — | homotypic SCFS | Omidvar 2014, 10.1016/j.jbiomech.2014.08.002 | **Y** | CONFIRMED |
| — absolute benchmark (**PROXY**, rescale by density²) | ~20 nN@30 s → ~200 nN@60 min; **SF ∝ (cadherin density)²** | Ca²⁺-dep (nil<100 µM, max ~2 mM) | S180-Ecad, dual pipette | Chu 2004, JCB 167:1183, 10.1083/jcb.200406044 | N | CONFIRMED |
| — per-bond clutch | catch peak ~30 pN; slip τ₀~0.63 s; rupture ~30–70 pN | loading-rate dep | recombinant ectodomain | Rakshit 2012, 10.1073/pnas.1208349109; Perret 2004, 10.1073/pnas.0402085101 | N | CONFIRMED |
| **Elastic (Young) modulus E** | **357 ± 31 Pa** (G = 115 ± 20 Pa, E/G=3.09) | 0.3–2 kPa; pooled 0.4–87 kPa | ADHERENT colloidal-probe AFM | Eldridge 2019, BpJ 117:696, 10.1016/j.bpj.2019.07.008 | **Y** | CONFIRMED |
| — corroboration | 300–400 Pa | — | adherent, spherical model | Mirzaluo 2023, 10.1007/s00249-023-01642-3 | **Y** | CONFIRMED |
| — relative | 1.4–1.8× softer than MCF-10A @37°C | abs ~0.1–0.4 kPa | adherent 37°C | Li 2008, 10.1016/j.bbrc.2008.07.078 | **Y** | CONFIRMED (fold) |
| **Cytoplasm viscosity η** | **65.9 ± 11.4 Pa·s** (all-wire) / **56.4 ± 16.6** (L=3 µm); G=32.9±6.0 Pa | 10–70 Pa·s | ADHERENT 37°C, active MRS | **Dessard, Manneville & Berret 2024**, 10.1039/d4na00003j | **Y** | CONFIRMED |
| — probe-scale caveat | 1–6× water (~0.001–0.015 Pa·s) | — | passive submicron OT, 25°C | Roy 2023, 10.3389/fphy.2022.1099958 | **Y** | CONFIRMED (diff scale) |
| **Osmotic turgor ΔP** | **~40 Pa** interphase; ~320–400 Pa mitotic | 10–100 Pa interphase | PROXY HeLa, AFM-confine | Fischer-Friedrich 2014, 10.1038/srep06213 | N | PROXY |
| — current code | turgor_dP0 = 133 Pa (Laplace 2γ/R) | 80–160 Pa | rounded (NOT spread) | synthesis | N | PLAUSIBLE |
| **Osmotic bulk modulus K** | **NOT constant**: B=N·k_BT·V/(V−Vmin)²; Pa→kPa→MPa | 100 Pa–1 MPa | PROXY osmotic P–V | Guo 2017, 10.1073/pnas.1705179114 | N | PROXY |
| — model default | K = 2500 Pa (log law p=−K·ln(V/V₀)) | — | generic epithelium | Runser/Vetter/Iber 2024, 10.1038/s43588-024-00620-9 | N | CONFIRMED (default) |
| **Radius / volume** | **V₀ = 1760 µm³; R_sph ≈ 7.5 µm** | susp diam 15–21 µm | spheroid NMR per-cell est | Gamcsik 1995 (PMID 7743493); Wagner 2011, 10.1016/j.freeradbiomed.2011.05.024 | **Y** | CONFIRMED |
| **Traction + spread area** | **102 ± 59 nN**; **1822 ± 886 µm²**; energy 0.0052 pJ | — | ADHERENT 9 kPa FN-PDMS | Gil-Redondo 2023, 10.1002/jemt.24368 | **Y** | CONFIRMED |
| — traction stress | **~60–65 Pa** (Liew prints "kPa" = ~1000× unit err) | — | adherent 14 kPa PA | Liew 2024, 10.1007/s12195-024-00811-4 | **Y** | PLAUSIBLE (corrected) |
| — morphology | round on ~0.1 kPa → spread on 4–17 kPa | ~hours | stiffness series | Gil-Redondo 2021, 10.1016/j.jmbbm.2021.104979 | **Y** | CONFIRMED |
| **Membrane area modulus k_a** (effective) | **~0.13–0.14 N/m** | reservoir-buffered | PROXY MDCK II resting | Brückner 2015, 10.1038/srep14700 | N | CONFIRMED |
| — stiff bilayer (reservoir gone) | **0.243 N/m** ±10% | 0.1–1.0 N/m | PC bilayer GUV | Rawicz 2000, 10.1016/S0006-3495(00)76295-3 | N | CONFIRMED |
| — reservoir / lysis | excess area >20–40%; **~2–3% strain then lysis** | A₀×1.2–1.4 | mixed mammalian | Brückner 2015; Figard & Sokac 2014 | N | CONFIRMED |

## 2. Surface-tension resolution (the ~50–100× discrepancy)

Three numbers are being conflated — they are **different physical quantities under different
conditions**, not three estimates of one thing:

| Value | What it actually is | Verdict |
|---|---|---|
| **2.7×10⁻⁴ N/m** | SimuCell3D downstream config, **below the paper's own 1×10⁻³ floor**, not measured, not MCF7 | **NOT literature-anchored — reject** |
| **0.35–0.65 mN/m** (HeLa/L929) | de-adhered proxy cortical tension in the **weak MPA window** (suction remodels cortex) | PROXY (not MCF7) |
| **21–45 mN/m** (Nagle 2022) | **tissue AGGREGATE** surface tension (cortex + adhesion), on **MCF10A/MCF10DCIS** | PROXY (line + scale) |
| **~10 mN/m = 1×10⁻² N/m** | **the only direct MCF7 single-cell cortical tension** (Moazzeni 2021, n=47, suspended) | **CONFIRMED, MCF7** |

**Recommendation:** do NOT use 2.7×10⁻⁴. Anchor the cortical-generation setpoint at the MCF7
suspended value **γ ≈ 1×10⁻² N/m** (band down to ~0.5 mN/m for the spread/remodeled state). This
is a **~40× upward re-anchor** and respects the physiological-baseline rule. **PI decision** (§4):
route (A) use Moazzeni suspended value as the adherent setpoint anchor, or (B) treat γ as EMERGENT
from the full compartment build (consistent with the project's active-γ-floor reframe). No adherent
MCF7 γ exists — cannot be resolved from data alone.

## 3. Dimensionless ratios for the contact re-tune

ℓ = V₀^(1/3) = (1.760×10⁻¹⁵)^(1/3) = **1.208×10⁻⁵ m** (≈12.1 µm, between R=7.5 µm and diam 15 µm ✓).
γ̄ = γ/(Kℓ); ω̄ = ω·ℓ/K; ξ̄ = ξ·ℓ/K (ω, ξ in Pa/m). **All ratios ∝ 1/K — fix K (§4.3) before
freezing stiffnesses.**

| Quantity | K = 1×10³ Pa (current) | K = 2500 Pa (SimuCell3D) | target (mcf7_p0.xml) |
|---|---|---|---|
| γ̄ (= γ=1e-2 N/m, Moazzeni) | **0.83** | 0.33 | — (tension absent from node-face form, MUST add) |
| ξ̄ (current ξ=2e8 Pa/m) | **2.42** | 0.97 | **0.48** → ξ ≈ **4×10⁷ Pa/m** (5× softer) |
| ω̄ (current ω=8e8 Pa/m) | **9.68** | 3.87 | **0** (cohesion → cadherin module) |
| turgor ΔP₀/K (133 Pa) | 0.13 | 0.053 | — |

At the physiological point tension dominates turgor ~6–20×. **Re-tune (turgor:tension:adhesion:
repulsion):** keep ΔP₀/K≈0.13; **add** γ̄≈0.33–0.83 (γ=1e-2 N/m); cohesion ω=0 (cadherin) or
ω̄≈0.5–1 if smeared (anchored to Chu density², NOT per-bond pN); **ξ≈4×10⁷ Pa/m** (ξ̄≈0.48).

## 4. PI decision points (gate the re-tune)

1. **γ condition & value (§2):** Moazzeni suspended 1×10⁻² N/m as adherent anchor (route A, ~40× up)
   OR γ emergent (route B)? Reject 2.7×10⁻⁴ either way.
2. **Where does cohesion live?** Smeared node-face ω OR discrete cadherin module? Findings favor
   **discrete cadherin (Rakshit kinetics, already in repo) + Chu density²-rescaled whole-cell SF**
   for any smeared term — never per-bond pN inside ω. (mcf7_p0.xml uses ω=0.)
3. **Bulk modulus K: 1×10³ → ?** keep 1e3 (dilute), 2500 (SimuCell3D), or crowding law (Guo 2017,
   B volume-dependent kPa–MPa)? **Rescales all §3 ratios (∝1/K) — freeze first.**
4. **Turgor operating point:** dP0=133 Pa implies a ROUNDED cell; adherent-spread baseline should
   sit ≤~40 Pa (Fischer-Friedrich interphase). Which is the production baseline?
5. **η: 65.9 (all-wire) vs 56.4 (L=3 µm)?** + fix attribution to Dessard 2024.
6. **Radius: volume-anchored 7.5 µm vs suspended-diameter central ~9 µm** (~1.7× volume diff; sets ℓ).
7. **Membrane k_a:** soft reservoir 0.13 N/m + 1.2–1.4× excess-area headroom → stiff 0.243 N/m +
   2–3% lysis cap; correct MDCK proxy for MCF7?

## 5. Open data gaps (no MCF7-specific measured source)

Adherent γ · osmotic turgor ΔP · osmotic bulk modulus K · membrane k_a · cell–cell adhesion ENERGY
(J/m²; all MCF7 data are *forces*, force→energy needs a model assumption — do not invent) — every
one is PROXY/geometry-inherited. The one MCF7 membrane datum (Pradhan 2021, BBRC 587:126) reports
apparent *tension*, value paywalled (retrieve via gbook/KAIST).

## 7. Empirical re-tune validation (gbook A5000, N=7, node-face + remesh)

Tested the §3 lit-anchored softer contact via `--rep-strength/--adh-strength`:

| config | dt | outcome |
|---|---|---|
| un-tuned ξ=2e8, ω=8e8 | 1e-4 | STABLE but **COMPRESSES**: V/V0 0.87→**0.54** (adhesion≫turgor) |
| **re-tuned ξ=4e7, ω=5e7** | 1e-4 | STABLE, **V/V0 HELD ~0.87–1.0** (compression FIXED), A/A0 ~1.0, 120 swaps |
| **re-tuned ξ=4e7, ω=5e6** (ω̄~0.06) | 1e-4 | STABLE, **V/V0 ~1.0 + A/A0 creeps 1.0→1.10** (cell holds volume + begins to spread), 31 splits |
| re-tuned ξ=4e7, ω=5e7 | 1e-3 / 3e-3 | **DIVERGES** (A/A0→379×, maxZ→323µm, V/V0→±1e4) — *identical* to un-tuned |

**Two separable results:**
1. ✅ **The lit-anchored re-tune WORKS at the stable dt** — it fixes the adhesion-driven
   compression/collapse (V/V0 0.54→~1.0) and the cell begins to spread (A/A0→1.10). The §3
   dimensionless analysis is empirically validated.
2. ⛔ **The dt ceiling is a SEPARATE, contact-INDEPENDENT instability.** dt=1e-3 diverges with
   the softened contact *exactly* as with the stiff one (V/V0 goes NEGATIVE = mesh inversion /
   large-dt overshoot-tangling). So contact stiffness is NOT what gates dt — a distinct
   integration instability (likely large-dt node overshoot → mesh tangle → turgor sign-flip) is.
   This must be resolved (soft-start ramp, velocity cap, or sub-stepping) to reach the dt≥3e-3
   "5× faster" runtime; until then the spread runs correctly at dt=1e-4.

⇒ **NEXT (post-PI-decisions): (a) add γ as a node-face surface-tension term (anchored ~1e-2 N/m,
PI route A/B) and freeze ξ/ω/K per §4; (b) diagnose + fix the contact-independent dt-ceiling
instability** (instrument per-force contributions at dt=1e-3, check mesh-inversion onset, add
soft-start). The aggregation re-derivation (brief §6) remains the other gate to a multi-cell A/A0.

## 6. Citation-integrity fixes (carry into the KB + code)

| # | Fix | Where |
|---|---|---|
| a | η 65.9 Pa·s: **Dessard, Manneville & Berret 2024** (10.1039/d4na00003j), NOT "Hu 2024" | code comment + memory |
| b | SimuCell3D = **Runser, Vetter, Iber 2024**, NOT "Runeberg-Roos/Villoutreix" | spec/KB |
| c | "Moon et al. 19.74 µm" MCF7 diameter — **unresolved provenance** (not in Alshareef 2013); do not cite | KB |
| d | Liew 2024 traction stress printed "kPa" is physically **~Pa** (~1000× unit error) | note |
| e | Eldridge 2019 second line is **BT474**, not MCF-10A | note |
| f | Guo 2017 uses **MCF10A**, not A549 | note |
| g | 2–3% lysis cap from **Brückner 2015**, not Figard & Sokac | note |
