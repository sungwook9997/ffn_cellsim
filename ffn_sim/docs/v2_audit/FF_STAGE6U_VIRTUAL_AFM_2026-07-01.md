# FF Stage 6U — virtual parallel-plate (AFM): measurement-protocol-consistent cortical tension

**Date:** 2026-07-01  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**Prompted by PI:** "cortical tension is MEASURED with the cell DEFORMED (indenter/AFM). We've been
measuring a RESTING cortex — shouldn't we press it the same way the experiments do?" — and the follow-up
that Nie is a DENSITY (imaging) measurement, not a tension (pressing) measurement.

This is the CLAUDE.md **measurement-protocol-consistency sanity gate**, which the γ-floor comparison had
been violating: we compared a RESTING method-of-planes γ_myo to a band measured by DEFORMING the cell.

## Two axes, not to be conflated (PI clarification)

| axis | quantity | method (does it press?) | sources |
|---|---|---|---|
| **density** | minifil/µm² | imaging / mass-spec — **NO pressing** | Nie 2015 (0.6), proteomics (10–90) |
| **tension** | γ [mN/m] | **AFM / parallel-plate / MPA — PRESSING/deforming** | Chugh, Fischer-Friedrich, Hosseini = the band |

The cortical-tension BAND comes ONLY from pressing experiments. Nie (our density input) never pressed.
So to compare our model's tension to the band we must **press the model the same way** — this stage.

## Build

`network_warp.plate_kernel` (two rigid parallel plates at z = centroid ± R0(1−strain), stiff penalty +
reaction accumulation) + `simulate_compressed_shell_on_device` (turgor uses the true **convex-hull volume**
for the non-spherical compressed shape; extracts γ_apparent = ΔP·R_eq/2 the liquid-drop way — exactly how
the experiment infers γ). Tests: `test_compression.py`. Self-consistency: the AFM force-derived pressure
F/(πa²) ≈ the turgor ΔP (≈535 vs 493 Pa at 10%).

## A5000 result — pressing reaches the band; the ~1000× floor is largely a PROTOCOL artifact

Moderate cortex, myosin ON vs OFF (blebbistatin analog), resting method-of-planes γ_myo = **2e-4 mN/m**:

| strain | γ_apparent (myosin ON) | (OFF) | myosin contrib | F [nN] | ΔP [Pa] |
|---|---|---|---|---|---|
| 0% | 0.01 | 0.01 | ~0 | 0 | 0 |
| 5% | **0.39–0.70** | same | 0.0003 | 2–5 | 180 |
| 10% | 1.7–3.0 | same | 0.0003 | 19–36 | 790 |
| 15% | 4.4 | 4.4 | 0.0004 | 70 | 1900 |
| 20% | 8.7 | 8.7 | 0.0007 | 175 | 3700 |

(MCF7 band 0.18–0.40; Chugh/Salbreux 0.35–0.65 mN/m.)

**Pressed γ_apparent enters the band at ~3–5% strain (~0.4–0.7 mN/m) — vs resting γ_myo 2e-4 mN/m,
i.e. ~10³–10⁵× higher.** So the ~1000× "floor" was substantially a **measurement-protocol artifact**:
resting prestress (method-of-planes) is NOT the quantity the band measures (deformation response).

## ⚠️ The honest reversal — pressed tension is TURGOR-borne, not myosin-borne

Myosin ON ≈ OFF (contribution ~0.0003 mN/m). In our model the pressed apparent tension is entirely the
**turgor/osmotic** response, NOT myosin. This is physically correct: **tangential cortical myosin does not
resist NORMAL (plate) compression — the turgor pressure does.** Myosin can only feed the pressed tension via
Laplace (P = 2γ/R), i.e. by raising the resting tension/pressure — but our myosin cannot generate band-scale
prestress (the generation floor), so it does not move the pressure and contributes ~nothing under compression.
The REAL band is ~70% myosin (blebbistatin) → our model reproduces the band **MAGNITUDE (via turgor)** but
NOT the **MECHANISM (myosin)**.

## Net reframe of the γ-floor (this stage relocates it)

- NOT "the model can't reach band-scale cortical tension" — it CAN, when measured the experimental way
  (pressed), via turgor: γ_apparent 0.4–8.7 mN/m at 5–20% strain.
- The floor is precisely: **myosin cannot GENERATE band-scale prestress, so it neither (a) sets a
  band-scale resting tension (method-of-planes γ_myo) nor (b) feeds the pressed AFM tension (Laplace).**
  The band magnitude in our model is carried by turgor; in the real cell it is carried ~70% by myosin.
- Combined with FF_STAGE6T (proteomics: cortex is NOT myosin-protein-limited; Nie undercounts ~20–150×)
  and FF_STAGE6H (network screening ~55–74×): the myosin-generation gap is **density-undercount ×
  network-screening**, and the "band" was additionally a **wrong-protocol / wrong-reference** comparison
  (resting vs pressed; isolated-cell vs tissue, per the DCM assembly evidence).

## Visualization (PI request: "I want to see how the filaments bend")

`outputs/ff/figs/afm_compression_montage.png` (+ `afm_compression.mp4`): cortex filaments colored by
curvature, strain 0→20%. The cell FLATTENS top/bottom and BULGES equatorially; **individual filaments bend
little** (the inextensible cortex accommodates compression by shape change, not local bending), with
curvature rising only locally near the plate contacts. γ_apparent annotated per frame (0.01 → 8.75 mN/m).

## Sanity gate
- **Protocol consistency:** γ now extracted the SAME way the band is (press → liquid-drop). ✓ (the whole point)
- **Self-consistency:** AFM force pressure F/(πa²) ≈ turgor ΔP. ✓
- **No magic number:** strain is the swept controlled variable; plates are a geometric constraint; γ_apparent
  is derived from ΔP + R_eq (measured), not tuned. Myosin ON/OFF is the blebbistatin control. ✓
- **Honest:** the result REVERSES a hoped-for closure — the pressed band-tension is turgor-borne, not the
  myosin mechanism; reported as such. ✓

## Files
- `ff/network_warp.py` — plate_kernel + simulate_compressed_shell_on_device (committed f327dbe).
- `tests/ff/test_compression.py`; figs afm_compression_montage.png + afm_compression.mp4.

Related: [[project-gamma-floor-likely-deficit]], FF_STAGE6T (proteomics density), FF_STAGE6P (band definition),
DCM_GAMMA_CONTROLLED_SWEEP (assembly fails at band-γ). Sources: Fischer-Friedrich 2014; Chugh 2017; Nie 2015 (density, NOT tension).
