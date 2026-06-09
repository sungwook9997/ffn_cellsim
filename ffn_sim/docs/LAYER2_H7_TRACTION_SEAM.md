# Layer-2 ↔ H.7 traction scale-bridge (seam) — CONSUMER contract

**Date** 2026-06-10 · **Branch** `layer2/spheroid-cbm` · **Owner** Layer-2 session (consumer side)
**Pairs with** H.7 `docs/v2_audit/H7_SF_ARRAY_TRACTION_SCALEUP_2026-06-10.md` §5 (supplier side)

## 0. One line

H.7 derived the single-cell active traction mechanistically (sarcomeric ventral-SF array,
**+131 ± 8 pN per SF**, +2.67 pN/head ≈ native F_stall); this seam consumes that value as the
Layer-2 CBM per-cell `f_active`, replacing the heuristic anchor that the (closed) Layer-2 line
flagged as its one remaining unanchored scale — **without fitting to the PI A/A0** (overlay-only).

## 1. The gap this fills

The Layer-2 center-based CBM line is CLOSED (2026-06-05): it reproduced the PI law
`A/A0 = a + b/R + c/R²` in FORM (r² = 0.998, zero-calibration, full PI R₀ range) but its MAGNITUDE
is bounded ~5–9× under the PI medians (7–10). The `b/R` term is the **single-cell active traction**.
Layer-2 only ever set it from HEURISTICS:

| heuristic anchor | value | where | nature |
|---|---|---|---|
| whole-cell migration speed × γ | ~1.6 nN | `motility_bridge.resolve_active_traction` | v₀(MCF7 0.32 µm/min) × γ_cell |
| protrusion velocity × γ | ~9.4 nN | `motility_bridge.resolve_active_traction` | v₀(Bieling barbed-end) × γ_cell |
| clutch-scale guess `T_ref` | 2.5 nN | `ligand_traction.py` | *"no direct MCF7 single-cell traction in lit"* |

H.7 now supplies the missing **mechanistic** quantity.

## 2. Supplier → consumer interface

### 2.1 What H.7 supplies

- **Now (DECISIVE single fiber):** `outputs/h7/production/h7_sf_2c_sarcomeric_gpu.json` →
  `coherent_differential_pN = 131.08`, `coherent_differential_sem_pN = 8.20`, `engaged_heads ≈ 49`.
  Vendored read-only into `outputs/layer2/seam_inputs/` (see that dir's README).
- **Pending (SF-array aggregate + seam record):** H.7 §5 designs
  `outputs/h7/production/h7_traction_seam.json` with the schema below, written by the SF-array GPU
  run (`h7_sf_array_n4_gpu.json`, IN-PROGRESS on gbook). When it lands, drop it into
  `outputs/layer2/seam_inputs/` (or leave it in the H.7 production dir — the consumer checks both).

```jsonc
// h7_traction_seam.json (H.7 supplies; Layer-2 reads — DO NOT edit H.7 code)
{
  "single_cell_total_traction_N":   <array aggregate at N_SF, physiological density>,
  "single_cell_traction_stress_Pa": <aggregate / contact area>,
  "contact_area_m2":                1822e-12,
  "per_SF_traction_N":              <per-fiber, e.g. 131.08e-12>,
  "per_SF_traction_sem_N":          8.20e-12,
  "n_SF":                           20,
  "model_density_gap_factor":       <model aggregate ÷ platform 102 nN>,
  "provenance":                     "H.7 sarcomeric ventral-SF array (continuous_stroke myosin)"
}
```

### 2.2 What Layer-2 consumes (this side)

`ffn_sim/spheroid/h7_traction_seam.resolve_h7_traction()` → `ResolvedH7Traction`:

| field | value (provisional) | meaning |
|---|---|---|
| `per_sf_traction` | 131.08 pN | decisive single sarcomeric-SF coherent contractile traction |
| `n_sf` | 21.5 | ventral SF/cell = **FA_count/2** (43/2; each SF spans 2 FAs; Hotulainen-Lappalainen 2006) |
| `f_cell_test` | **2.82 nN** | per_sf × n_sf — raw H.7 mechanism at the validated TEST-fiber density |
| `f_cell_physio` | **102 nN** | Gil-Redondo 2023 MCF-7 total contractility (physiological-density anchor) |
| `stress_physio_Pa` | 56 Pa | f_cell_physio / 1822 µm² (≈ Gil-Redondo 63 Pa ✓) |
| `density_gap_factor` | ×36 | physiological ÷ test-density (≈ sanctioned ×40 mesoscale; the honest residual) |

**Reading rule:** if `h7_traction_seam.json` is present → use its measured array aggregate; else
PROVISIONAL from the single-SF × N_SF (status `provisional=True`). Until the array lands, the
per-cell value is provisional but the *order* (few-nN test, ~10²-nN physiological) is fixed.

## 3. The dimensional bridge (no magic number)

```
per-SF [N]  ──×  N_SF (= FA_count/2, anchored)  ──►  f_cell_test [N]
                                                 ──►  f_cell_physio [N] = Gil-Redondo 102 nN
                                                      (= mechanism × physiological minifilament density;
                                                       density_gap = f_physio / f_test ≈ ×36)
```

Every factor is derived or cited: per-SF from the H.7 artifact; N_SF from the FA-pairing anchor
(KU-2.4 FA count + Hotulainen-Lappalainen SF/FA topology); physiological scale from Gil-Redondo
2023 TFM. **Nothing is fitted to the PI A/A0.**

## 4. CBM reduction (the consumer's coarse-graining, NOT H.7's)

H.7 supplies a *fully-spread isolated single-cell* traction. The CBM `f_active` is a *per-edge-cell*
outward traction with a B1 ≤ 3 nN stability ceiling (above it the overdamped step `v=F/γ·dt` flings
a marginally-detaching edge cell). The two differ by spread-state and rim-vs-whole-cell projection.
The decisive experiment (§5) feeds the seam value directly as `f_active` and maps `A/A0(f_active)`
to read the fork empirically rather than assuming a rim-reduction factor.

## 5. Decisive fork — resolved (see `scripts/layer2_h7_decisive.py`, §6 below)

`A/A0(f_active)` curve at matched R₀ vs the PI band [7–10]:
- if the curve ENTERS [7–10] at the H.7 mechanistic value → gap was MISSING MAGNITUDE (H.7 supplies);
- if it SATURATES ~1–2 then EJECTS before [7–10] → gap is the 1-particle CBM ABSTRACTION LIMIT,
  and H.7's mechanistic traction confirms route a (fine-grained single-cell) owns the magnitude.

*(Result table filled in §6 after the run.)*

## 6. RESULT (2026-06-10, `layer2_h7_decisive.py` N0=300, 2 seeds, CPU, matched R₀≈60 µm)

A/A0(f_active) connected-core, ensemble mean of 2 seeds (figure
`outputs/layer2/figs/fig_layer2_h7_decisive.png`):

| f_active | provenance | A/A0_core | ejected |
|---|---|---|---|
| 0 nN | baseline (substrate only) | 1.99 | no |
| 1.6 nN | heuristic whole-cell | 2.07 | no |
| **2.82 nN** | **H.7 mechanistic, TEST density** | **2.17 (peak)** | no |
| 3.0 nN | B1 ceiling | 2.14 | no |
| 6.0 nN | curve-fill | 1.71 | **EJECT** |
| 9.4 nN | heuristic protrusion | 1.39 | **EJECT** |
| **102 nN** | **H.7 mechanistic, PHYSIOLOGICAL (Gil-Redondo)** | **0.31 (fractured)** | **EJECT** |

**⭐ FORK RESOLVED → CBM 1-particle ABSTRACTION LIMIT (now with a MECHANISTIC anchor).**
A/A0 PEAKS at ~2.2 right at the H.7 mechanistic test-density value (2.82 nN), then DECLINES and
EJECTS as f_active rises; the physiological 102 nN value catastrophically FRACTURES the spheroid
(core → 0.3, scattered debris inflates the raw hull). **No f_active — including either mechanistic
H.7 value — brings A/A0 into the PI band [7–10].** The earlier decisive verdict (structural limit)
held a HEURISTIC bracket (1.6 / 9.4 nN); this run replaces it with H.7's mechanistically-DERIVED
single-cell traction and reaches the SAME verdict, closing the open caveat the closed-line flagged.

**Interpretation.** The center-based 1-particle CBM has a structural ceiling A/A0 ≈ 2.2: a point
cell driven by an edge-outward traction either stays cohesion-locked (low f) or detaches (high f) —
it cannot crawl-while-attached or distribute contact-line traction to spread 7–10×. The H.7
mechanism SUPPLIES a correct, mechanistic single-cell traction; the magnitude that traction should
produce lives in the FINE-GRAINED single-cell representation (route a), exactly as the closed-line
predicted — not in a parameter the CBM was missing. H.7's value at validated test density lands the
CBM at its own ceiling (~2.2); the physiological scale just fractures the abstraction.

This is the mechanistic twin of the cortical-γ density story: the per-SF minifilament density gap
(×36, §3) is the documented fine-grained-prediction residual, and the CBM cannot host the
physiological traction regardless. The seam itself is VALIDATED — it cleanly turns +131 pN/SF into
a dimensionally-grounded per-cell traction (2.82 nN test / 102 nN physiological) that the CBM
consumes — it is the CBM representation, not the seam or the traction value, that bounds A/A0.

## 7. Sanity gates (CLAUDE.md — recorded before the run)

1. **Dimensional:** per_sf [N] · n_sf [–] = f_cell [N]; stress = f_cell / area [Pa]. ✓ (§3)
2. **Boundary:** per_sf = 0 ⇒ f_cell = 0 ⇒ A/A0 reduces to the cohesion-locked G1 baseline.
3. **Sign:** only the SARCOMERIC (+131 pN, contractile) value bridges; the mixed-polarity SMOKE
   (−, slackening) is NOT consumed.
4. **No fit:** N_SF, per-SF, physiological scale are all anchored/cited; f_active is swept, never
   tuned to land A/A0 in the PI band. Band [7–10] LOCKED.
5. **Integrator freeze:** dt-substep for f_active > ceiling is the PI-authorized Layer-2-only
   params copy (dt/factor + steps×factor, same biological time) — the frozen integrator is untouched.
6. **Provisional flagged:** until the H.7 array aggregate / seam record lands, the per-cell value is
   provisional (single-SF × N_SF) and labelled as such on every result.
