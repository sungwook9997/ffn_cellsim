# FF Stage 6X — plasma membrane (H.8) wired into the whole-cell AFM; the full cell gives the same cortical tension

**Date:** 2026-07-01  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**PI directives:** "진행" + "진짜 전부 다 필요할텐데" — complete the whole-cell compartment set (nucleus ✓ 6W,
membrane here, cytoplasm already shared). **Follows** FF_STAGE6W (nucleus wired).

## One line

The plasma membrane is now an additive, lit-anchored channel in the whole-cell AFM (a distinct surface
tension the FF cortex lacks), and a native run with **all compartments ON** (cortex + turgor + nucleus +
membrane) gives the **same confinement-independent cortical tension** (~0.15 mN/m) as cortex + turgor alone —
i.e. the cortical-tension result (γ-floor, 6V) is robust to the full compartment stack.

## Part 1 — lit-anchored params (parallel-researched + unit-corrected)

A workflow researched the three membrane constants (corpus BM25 + web) with citations, plus a code-level
design check. Values in FF units (**verified against the cortex ground truth** γ_cortex = ½·40 Pa·7.5 µm =
150 pN/µm ≡ 0.15 mN/m ⇒ **1 N/m = 1e6 pN/µm, 1 µN/m = 1 pN/µm; energy 1 J = 1e18 pN·µm**):

| param | FF value | lit range | primary source (DOI) | in corpus |
|---|---|---|---|---|
| γ_mem (bilayer in-plane tension) | **10 pN/µm** | 3–40 µN/m | Diz-Muñoz 2013 (10.1016/j.tcb.2012.09.006); KB-3.B1.1 | yes |
| K_A (area-expansion modulus) | **2.35e5 pN/µm** | 0.2–0.3 N/m | Rawicz 2000 (10.1016/S0006-3495(00)76295-3); KB-3.B1.3 | yes |
| κ (bending rigidity) | **0.0828 pN·µm** | 10–30 kBT | Rawicz 2000; KB-3.B1.2 | yes |

**Unit-error caught:** a research draft used a wrong 1 N/m = 1e3 pN/µm (γ_mem→0.01, K_A→235) — off by 1000×.
The cortex ground-truth cross-check fixed it (correct γ_mem = 10 pN/µm, K_A = 2.35e5 pN/µm). γ_mem is the
**bilayer-only** tension (NOT the apparent 30–300 µN/m, which includes membrane-cortex adhesion γ_MCA already
carried by the cortex — using it would double-count the cortex).

## Part 2 — design: no double-count; buffered-plateau, not bare-K_A

A code-level design agent confirmed (file:line cited): **the FF cortex has NO area-elastic term** —
`reshape_kernel` constrains segment length only (a "fishnet": inextensible edges, free area), bending
penalizes curvature, crosslinks are force-free at construction, and only turgor resists *volume* (not area).
So the membrane area/tension channel is **genuinely distinct** — additive (`γ_total = γ_membrane + γ_cortex`,
Sens & Plastino 2015), not a duplicate.

**Buffered-plateau, not bare-K_A (empirically decided).** A bare fixed-A0 K_A law (γ=γ_mem+K_A·(A−A0)/A0) is
**knife-edge**: because K_A is huge, any area shrink below A0 drives γ negative → clamps to 0 (slack), and any
growth spikes it explosively. Under our regulated turgor the cell *sheds* area under compression (A/A0 →
0.845 at strain 0.8, measured), so bare-K_A is **always slack** and cannot hold the physiological baseline
tension (verified: γ_mem_channel → 0). The membrane physically holds a **reservoir-buffered ~constant
tension γ_mem** over normal deformations (the folds/microvilli/caveolae reservoir unfolds to keep in-plane
tension near-constant; Raucher & Sheetz 1999). So we wire the **constant baseline γ_mem** (buffered plateau);
the K_A elastic upturn past reservoir capacity (areal strain > f_excess) is **deferred** — the reservoir DOF
is PI-blocked (f_excess unknown for MCF7).

## Part 3 — wiring

`simulate_whole_cell_compression_on_device(..., membrane=ResolvedMembrane)`: an inward **isotropic** Young-
Laplace pressure ΔP_mem = 2γ_mem/R applied on the cortex nodes (implemented as the turgor kernel with a
negative, inward dP·area/Nc; cortex-only `dim=Nc`) — a shell-baseline representation, not a true tangential
in-plane area-tension force. γ_mem is a **separate additive force channel**; the turgor Laplace partner
ΔP·R/2 is NOT re-reported as a passive tension (the one double-count risk, against turgor — avoided). The
`τ_lysis` clamp is present but inert (γ_mem = 10 ≪ 5000). **Metric note:** `gamma_apparent` (turgor Laplace,
½·ΔP·R) and `gamma_mem_channel` (bilayer 10 pN/µm) are reported as SEPARATE fields and must NOT be summed —
γ_mem is the bilayer-only channel; the cortex-adhesion tension is already inside the turgor-borne γ_apparent.

## Part 4 — native result (N=70686, regulated ΔP=40 Pa, all compartments ON)

| strain | F cortex+turgor | F full cell | γ cortex+turgor | γ full cell | γ_mem channel |
|---|---|---|---|---|---|
| 0.10 | 8.459e6 | 8.459e6 | 0.1502 | 0.1502 | 10.00 |
| 0.30 | 7.720e7 | 7.720e7 | 0.1502 | 0.1502 | 10.00 |
| 0.50 | 2.146e8 | 2.146e8 | 0.1508 | 0.1508 | 10.00 |

- **γ_apparent is compartment-independent** (nucleus + membrane change it by −0.000%) and confinement-
  independent (~0.15 mN/m ≈ Fischer-Friedrich interphase 0.17). The full whole cell gives the SAME cortical
  tension as cortex + turgor — the 6V γ-floor resolution is robust to the complete compartment stack.
- The membrane holds its **baseline γ_mem = 10 pN/µm** (ΔP_mem ≈ 2.7 Pa inward) — the physiological membrane
  tension is ON.
- **The membrane also leaves F_plate essentially unchanged** — Δ = −0.0008% / −0.0002% / −0.0001% at strain
  0.1/0.3/0.5 (adversarial-verified). ΔP_mem is ~7% of the turgor pressure (2.68–2.89 Pa vs 40 Pa), but that
  pressure share does **NOT** reduce the plate force: under regulated turgor ΔP is held at the 40 Pa setpoint,
  so F_plate ≈ ΔP·contact_area is invariant (R_eq and contact radius bit-identical between arms); the membrane
  only slightly lowers V/V0 (0.98985→0.98946). (Earlier framing that "the membrane bears a ~7% share of the
  plate force" was a misattribution of dP_mem/dP_turgor and is retracted.) *Confound note:* the full arm
  bundles the nucleus (n_nuc=642), so the ~0.001% F delta is not a membrane-isolated number — a membrane-only
  arm would isolate it if a cleaner value is ever needed.

## Part 5 — honest caveats (PI-roadmap)

1. **K_A elastic upturn + reservoir deferred.** The membrane's area-incompressibility (K_A) engages only past
   the reservoir capacity (areal strain > f_excess) — relevant for **spreading / blebbing** (area growth), not
   the resting/compression regime here. The reservoir DOF is PI-blocked (MCF7 f_excess unknown) → surface to
   PI before any large-area-strain production run.
2. **Regulated-ΔP masks tension partition on F_plate** (see Part 4) — a force-derived γ (real AFM extraction)
   would resolve the membrane share; our turgor-Laplace γ (a proxy, 6W caveat) does not.
3. **γ_mem is bilayer-only**; the apparent membrane tension (incl. cortex adhesion) is represented by the
   cortex side, cell-type specificity (MCF7 ~2× lower apparent tension, Tsujita 2021) belongs there.

## Sanity gate
- No magic number: γ_mem/K_A/κ are lit values (Rawicz/Diz-Muñoz; KB-3.B1), band-guarded in `resolve_membrane`,
  in the corpus; units verified against the cortex ground truth. ✓
- Physiological baseline: membrane ON at γ_mem = 10 pN/µm from the start (buffered plateau). ✓
- No double-count: cortex has no area term (code-verified); γ_mem is a separate channel; ΔP·R/2 not re-reported. ✓
- No regression: additive `membrane=None` default; nucleus/cortex paths unchanged; tests + kb-check green. ✓
- Honest: buffered-plateau chosen over bare-K_A (the native run has area < A0 at every strain → bare-K_A
  permanently slack); K_A/reservoir deferral stated; the F_plate claim was corrected by verification (below). ✓

## Adversarial verification

A 4-lens skeptic panel (units/physics, code, model/physiology, literature) tried to refute claims M-A..M-E;
a synthesizer aggregated. **M-A (units), M-C (no double-count), M-E (code/sign) confirmed (4/4); M-B (model)
confirmed** (the native area<A0-at-every-strain data backs the knife-edge argument). **M-D was REFUTED and the
doc corrected**: an earlier draft claimed the membrane lowers F_plate ~7%; the true native ΔF is ~0.001% —
the 7% was a misattribution of dP_mem/dP_turgor (the membrane's *pressure* share). Under regulated turgor
F_plate ≈ ΔP·contact_area is invariant; the membrane only nudges V/V0. Part 4 reflects the corrected finding.
The units cross-check (150.177 pN/µm ground truth) and the sign/dim/clamp wiring were independently reproduced.

## Figures
- `outputs/ff/figs/fullcell_membrane_native.png` — Panel A: F_plate(strain) cortex+turgor vs full cell
  (overlap). Panel B: γ_apparent(strain) flat ~0.15 mN/m for both = compartment-independent, interphase band
  overlaid, membrane baseline γ_mem annotated.

## Files
- `common/compartments.py` — `resolve_membrane` + `ResolvedMembrane` (lit-anchored, band-guarded).
- `ff/network_warp.py` — membrane channel in `simulate_whole_cell_compression_on_device` (commit f034952).
- `tests/ff/test_compartments_shared.py` (resolve_membrane bands), `tests/ff/test_compression.py`
  (whole-cell membrane additive/small/stable).

Related: [[project-cell-mechanics-extend]], FF_STAGE6W (nucleus), FF_STAGE6V (pressure-borne tension). Sources:
Rawicz 2000 (10.1016/S0006-3495(00)76295-3); Diz-Muñoz 2013 (10.1016/j.tcb.2012.09.006); Raucher & Sheetz 1999.
