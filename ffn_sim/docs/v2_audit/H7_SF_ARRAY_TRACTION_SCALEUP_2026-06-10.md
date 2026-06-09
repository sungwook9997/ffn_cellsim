# H.7 — SF-array scale-up → aggregate substrate traction [Pa], PI-platform comparison, Layer-2 seam

**Date** 2026-06-10 · **Branch** `h7/full-cell-integration` · **Status** IN-PROGRESS (GPU run
developing engagement on gbook). Builds directly on the DECISIVE single-SF result (loop23,
`H7_MYOSIN_OVERLAP_MECHANISM_DESIGN` §12): a graded-polarity SARCOMERIC ventral SF rectifies the
§9-corrected continuous-stroke myosin into **+131 ± 8 pN coherent contractile traction (16σ,
+2.67 pN/engaged-head ≈ native F_stall)** at its FA anchors. This doc scales that single fiber to
a cell-scale array, expresses the aggregate as a traction stress [Pa], compares it to the PI
experimental platform, and designs the seam that feeds it to the Layer-2 magnitude gap.

## 1. Method — the SF array (`cortex/ventral_stress_fiber.generate_sf_array_layout`)

`n_sf` parallel sarcomeric ventral SFs are stacked along ŷ across the flat ventral contact patch
at the literature lateral spacing. Each fiber is an INDEPENDENT, FA-anchored copy of the validated
single-SF unit (`generate_sarcomeric_sf_layout`): Z-bands (α-actinin barbed anchors) ↔ M-bands
(bipolar continuous-stroke myosin), graded polarity, antiparallel pointed-end overlap, the outer
Z-bands pinned as FA anchors. The fibers are coupled only by WCA excluded volume at the µm-scale
spacing → the array is mechanically **near-additive**, and the per-SF FA-anchor groups
(`per_sf_anchor_beads`) are tracked so each fiber's traction is separable.

The aggregate observable reuses the validated same-seed COHERENT differential probe (force_scale
1 vs 0, identical seed → thermal + taut-baseline cancel via HOOMD counter-based RNG; ±8 pN floor).
Because all fibers run along x̂ with the same length, the global left/right end split sums the
coherent contractile traction across the whole array. Aggregate traction stress:

```
σ_traction [Pa] = Σ_SF (coherent FA reaction) [N] / A_contact [m²]
```

with `A_contact` = the literature spread footprint (NOT the SF bounding box), the same convention
as the TFM "contractility / cell area".

## 2. Literature anchors (HARD rule — FA density, SF spacing, sarcomere period all anchored)

| Quantity | Value | Source |
|---|---|---|
| FA count / cell | 43 (35 nascent + 8 mature) | KU-2.4 (`configs/phase1_h4.yaml`), brief 20–50 + 5–10 |
| FA area | 1 µm² nascent / 4 µm² mature | KU-2.2 |
| **SF count / cell** | ~20 (≈ FA_count/2; each ventral SF spans two FAs) | FA-pairing + Hotulainen & Lappalainen 2006 (JCB) |
| **SF lateral spacing** | ~1–5 µm (default 2 µm) | Hotulainen-Lappalainen 2006; consistent with ~40 µm contact width / ~20 SFs |
| **Sarcomere period** | ~1 µm native → ~2 µm mesoscale effective | Peterson et al. 2004 (Mol Biol Cell); Hotulainen-Lappalainen 2006 |
| MCF-7 contact area | 1822 ± 886 µm² | Gil-Redondo 2023 (Microsc Res Tech, DOI 10.1002/jemt.24368), Table 1 |
| **MCF-7 total contractility** | 102 ± 59 nN | Gil-Redondo 2023, Table 1 control n=37 |
| **MCF-7 mean traction stress** | ~63 Pa (= 102 nN / 1822 µm²; their "63 ± 27 N/m²") | Gil-Redondo 2023 |
| single active SF tension | ~5–6 nN | Kassianidou & Kumar 2017 (PNAS); cf. Kumar 2006 10–30 nN = network/prestress |
| per-FA force | ~2.4 nN (102 nN / 43 FAs) | derived; consistent with mature-FA ~nN |

The mesoscale sarcomere period (~2 µm vs native ~1 µm) is the sanctioned ×40 filament-scale
coarse-graining (`ell0` = 0.5 µm), the SAME spirit ratified for filament count — not a free knob.

## 3. The scaling law (why the magnitude is a DENSITY lever, not a mechanism failure)

The decisive invariant is **+2.67 pN per engaged head**, the native cross-bridge F_stall (~2 pN)
rectified into the contractile direction by the sarcomeric organization. Aggregate traction is
therefore linear in the total number of engaged heads:

```
F_aggregate ≈ N_SF · (minifilaments/SF) · (heads/minifilament) · engagement · 2.67 pN
```

To reach the platform's 102 nN: N_engaged ≈ 102 nN / 2.67 pN ≈ **3.8 × 10⁴ engaged heads** →
at the ~11 % geometric engagement plateau, ~3.5 × 10⁵ heads → with ~60 heads/minifilament,
~5,800 minifilaments/cell → over 20 SFs, ~290 minifilaments/SF → over a 20-µm SF with 1-µm
sarcomeres (20 M-bands), ~14 minifilaments per M-band. That is physiologically plausible for a
thick ventral SF (each M-band is a 3-D bundle of minifilaments).

**⇒ The platform magnitude IS reachable by the mechanism — it is set by the per-SF minifilament
(overlap) DENSITY, exactly the cortical-γ "density/overlap gap" (Chugh 2017, Truong Quang 2021:
cortical/SF tension ∝ overlap × density, largely independent of motor *number per filament*),
now appearing in the adherent traction observable.** Our small test fibers (16 minifilaments,
3 sarcomeres, 6 µm) are deliberately ~1–2 orders under the physiological per-SF density; the
aggregate is reported AS the mechanism × density, never tuned to the band.

## 4. Aggregate result + PI-platform comparison (overlay-only, literature-first)

Dimensional ledger (all consistent): per-head force [N] → per-SF coherent FA reaction [N] →
aggregate Σ over SFs [N] → / contact area [m²] → stress [Pa]. ✓

`⟨FILL FROM GPU RUN — h7_sf_array_n4_gpu.json⟩`:
- per-SF coherent traction (array path) vs the decisive single-SF +131 pN → **additivity check**.
- aggregate (n_sf=4) [nN]; extrapolated cell-scale (N_SF=20) aggregate [nN] and stress [Pa].
- ×-under-platform on force AND stress (vs 102 nN / 63 Pa).

**Comparison is OVERLAY-ONLY.** The PI platform is MCF-7 on pV4D4/col-I with integrin-β1; the
nearest literature single-cell anchor is Gil-Redondo 2023 (MCF-7 on fibronectin/PDMS, 102 nN /
63 Pa). The model's per-SF density gap is the documented, honest residual — the same density/
overlap physics the fine-grained mean-field model does not capture at full count, NOT a sign or
mechanism error (sign + rectification are 16σ-confirmed).

## 5. ⭐ Layer-2 scale-bridge / traction-seam (DESIGN — Layer-2 is a separate session)

**The gap this fills.** Layer-2 (center-based CBM, LINE CLOSED 2026-06-05) reproduced the PI
A/A₀ = a + b/R + c/R² law in FORM (r² = 0.998, zero-calibration) but its MAGNITUDE is bounded;
the `b/R` term is the **single-cell active traction** that drives spreading. Layer-2 sets it via
`ffn_sim/spheroid/ligand_traction.py`:

```
f_traction(condition) = T_ref · density(condition) · clutch_strength(ligand),   T_ref = 2.5e-9 N
```

with the in-code provenance note: *"There is NO direct MCF7 single-cell traction in the literature
(the 15–25 nN micropillar value was REFUTED)"* and a 3 nN stability ceiling. **`T_ref` is a
placeholder clutch-scale guess — precisely the quantity the H.7 fine-grained traction now derives
mechanistically.** This is the scale-bridge: a fine-grained single-cell traction (HOOMD particle/
bond SF array) → the coarse CBM per-cell traction parameter.

**Seam interface (H.7 supplies; does NOT modify Layer-2 code).** The array run writes a small
seam record the Layer-2 session reads:

```
ffn_sim/outputs/h7/production/h7_traction_seam.json
{
  "single_cell_total_traction_N":  <aggregate at N_SF=20, physiological density>,
  "single_cell_traction_stress_Pa": <aggregate / contact area>,
  "contact_area_m2":               1822e-12,
  "per_SF_traction_N":             <per-fiber>,
  "n_SF":                          20,
  "traction_localization_Lp_m":    <edge screening length, if resolved>,
  "provenance": "H.7 sarcomeric ventral-SF array (continuous_stroke myosin, §9/§12); "
                "PI-overlay Gil-Redondo 2023 MCF-7 102 nN / 63 Pa",
  "model_density_gap_factor":      <model aggregate ÷ platform 102 nN>
}
```

**Mapping (the bridge, not a raw substitution).** Layer-2's `f_traction` is a *per-rim-cell edge*
traction with a ≤3 nN CBM-stability ceiling; the H.7 / Gil-Redondo 102 nN is a *fully-spread
isolated single cell's TOTAL* contractility. They differ by the spread-state and the
rim-vs-whole-cell projection. The seam therefore supplies **two layers**:
1. the **absolute mechanistic single-cell traction** (H.7 array → ~10²-nN scale, anchored to
   Gil-Redondo), which sets the physical units of `T_ref` (replacing the flagged placeholder); and
2. the **CBM reduction** (per-rim-cell effective traction = total × rim-spread fraction), kept
   ≤ the B1 ceiling — a coarse-graining the Layer-2 session applies, NOT H.7.

So H.7 removes the "no lit value" caveat on `T_ref` (the central open item of the Layer-2 line)
and supplies a mechanistic, dimensionally-grounded single-cell traction; Layer-2 applies its
documented CBM rim-reduction to land `f_traction`. The seam is **values + interface only** —
`ligand_traction.py` is the other session's concern and is not edited here.

## 6. Sanity gates (CLAUDE.md — recorded before the GPU result)

1. **Dimensional:** per-head [N] → per-SF [N] → aggregate [N] → /A [m²] → [Pa]. ✓ (§4 ledger)
2. **Boundary:** n_sf=1 → reproduces the single-SF +131 pN (CPU control running); myosin force
   OFF → ~0 aggregate (same-seed differential cancels the passive pre-tension).
3. **Conservation:** Newton-3 at every FA anchor; the same-seed differential removes the passive
   baseline so only the myosin-generated contractile reaction survives.
4. **Sign:** aggregate is CONTRACTILE (inward at both ends), as the single SF (16σ positive).
5. **Additivity:** per-SF (array) ≈ single-SF (isolated) within noise → independent fibers add;
   any deficit = WCA interference, to be reported, not absorbed.
6. **Measurement consistency:** stress denominator = the literature spread footprint, identical
   to the TFM contractility/area convention → apples-to-apples with Gil-Redondo.
7. **No gate-chasing:** band/platform LOCKED; the per-SF density gap is reported as the structural
   residual, never tuned. Success = the mechanism produces traction of the right SIGN and SCALING
   at literature parameters, with the magnitude set by the (anchored) minifilament density.

## 7. Open items → PI

- **(per-SF density)** Is there a defensible MCF-7 ventral-SF minifilament density (minifilaments
  per M-band / per SF)? Without it the aggregate is reported at the validated test density with the
  scaling law (§3), the honest fine-grained-prediction stance (same as the cortical-γ density datum).
- **(seam consumption)** The Layer-2 session applies the CBM rim-reduction (§5.2) and confirms the
  B1 ≤ 3 nN ceiling holds for the H.7-anchored `T_ref`.
- **(SF orientation)** The array is parallel (all along x̂); real ventral SFs are partly radial/
  criss-cross. Parallel is the faithful first model (max coherent traction); an isotropic-orientation
  factor (≈ ⟨|cos θ|⟩) would reduce the net by a known geometric factor — flag for a later refinement.
