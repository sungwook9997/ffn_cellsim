# MCF7-ification of the main single-cell line — relay spec (2026-06-02)

> Read-only audit by the Layer-2 session (no main-line files edited). For the **main-line
> session** to apply. Goal: build the single-cell line on an **MCF7** basis so it shares ONE
> MCF7 anchor set with the Layer-2 multicellular line (no generic-vs-MCF7 drift). Principle:
> **no global cancer modulus** — MCF7 identity is a per-compartment vector of literature-
> anchored presets, orthogonal to the substrate/ligand axis. Anchors cross-ref:
> `PI_EXP_VALIDATION_MAP.md`, `CANCER_CELLTYPE_PARAM_MAP.md`, `LAYER2_ANCHORS_2026-06-02.md`.

## Change table

| Parameter | Current (generic) + file:line | → MCF7 value + source | Change type |
|---|---|---|---|
| Cell radius R_cell | `10.0e-6` m — `configs/phase1_h3.yaml:40` | **7.5e-6 m** (diam 14.8 µm) — Wagner 2011, PMC3147247 (Coulter, measured) | **config-edit** ⚠️ re-validation (§B) |
| Cytoplasm η_eff | `cell_type: null` (water) — `configs/phase1_h10.yaml:31` | **65.9 Pa·s** via `cell_type: MCF7` — Hu 2024, PMC10929591 (already in `cytoplasm.py` table) | **config-edit** (set `enabled: true`, `cell_type: MCF7`) |
| Col-I clutch (ligand) | FN-like α5β1 catch-slip default — `configs/phase1_h4.yaml:78` | **collagen_I → α2β1 slip**: k_off0 **1.3/s**, x_β **2.3 Å** — Taubenberger 2007; already in `ligand_species.LIGAND_REGISTRY["collagen_I"]` | **needs-code** (wire `ligand_species` → `bridge/fa.py`; registry currently unconnected) |
| Intracellular shear G | Tier-2 GLE not implemented — `phase1_h10.yaml:70` | 32.9 Pa (Hu 2024) | **needs-code, PI-gated** (frozen integrator; defer) |
| Nucleus E_nuc / lamin | generic `5.0e3` Pa — `cell/nucleus.py:233` | **MCF7 value MISSING** | **needs-research** (don't invent; keep generic + tag) |
| Cortical tension γ | **emergent** in H.3 (no knob); myosin `n_motors_per_cell:100` | **no MCF7 value** (only pooled/refuted) | **needs-research** (documented absence) |
| Membrane tension γ_mem | `3.0e-5` N/m generic — `phase1_h8.yaml:51` | no MCF7-specific value | **needs-research** (keep generic) |
| Whole-cell modulus E₀ | n/a (emergent) | 0.26 kPa (Yubero 2020) | **overlay-only — NEVER a knob** |
| E-cadherin status | n/a (no cell-cell layer) | MCF7 = **E-cad POSITIVE** | no single-cell edit (Layer-2 only); don't add N-cad/EMT |

## A. Do-now (config-only, high-confidence, already MCF7-anchored + shared with Layer-2)
1. `phase1_h3.yaml:40` R_cell `10.0e-6` → `7.5e-6` (**but triggers §B re-validation**).
2. `phase1_h10.yaml` cytoplasm: `enabled: true`, `cell_type: MCF7` (→ η 65.9 Pa·s). Validation target is the MCF7/MDA **ratio** (~5×), not absolute Pa·s (η∝L²; state L≈3 µm).

## B. ⚠️ R_cell 10→7.5 µm is a RE-VALIDATION EVENT, not a free scalar
1. **Box**: `L_box = L_box_over_R_cell·R_cell` (factor 3.0). Absolute filament-fluctuation cushion (σ_perp≈1.56 µm, ℓ_0=0.5 µm) does NOT shrink → re-derive that factor 3.0 still clears `half_box ≥ R_cell+σ_perp+ℓ_0` (it does: 11.25 ≥ 9.6 µm) and fix the now-stale worked numbers in the header comment.
2. **Osmotic ΔP / enclosed volume** (`phase1_h3.yaml:233-257`): the 200/400 Pa Laplace cross-check was at 10 µm; smaller R → higher ΔP=2γ/R → re-evaluate.
3. **Cortex budgets**: surface 4πR² drops 44% (1257→707 µm²); confirm myosin density (3/µm²) + filament areal density + ×40 meso count + bead caps stay consistent.
4. **All landed H.3 results (`outputs/h3/`) are at 10 µm → invalidated for MCF7.** Re-run H.3 gates + refresh `outputs/h3/figs/` (visualize-at-closeout) before claiming the MCF7 line validated.
5. **FA capture geometry**: south-pole-to-substrate gap shrinks ~25%; re-confirm `fa_clutch_capture_radius` override still seeds the intended bond count.

## C. Implement as a `cell_type: MCF7` preset (not scattered edits)
Per-compartment vector, keyed off `cell_type` (the project's existing design,
`CANCER_CELLTYPE_PARAM_MAP.md:92`). Cytoplasm already supports `cell_type: MCF7`. Have the
preset set R_cell from the SAME Wagner-2011 anchor Layer-2 uses (one diameter source).
Carry nucleus/cortex/membrane as `null`/generic with an explicit "MCF7 value absent" tag so
generic is never silently passed off as MCF7. **Keep ligand/ECM OUT of the preset** (it's a
per-run substrate axis — col-I for the PI exp). Suggested home: `configs/celltype_mcf7.yaml`
overlay or a `cell_type:` key threaded through the H.3/8/9/10 resolvers.

## D. Over-claim flags (must accompany any MCF7 run)
- Whole-cell stiffness **ordering is method-disputed** (AFM: MCF7 softer 0.26 kPa; electrodef: MCF7 ~10× stiffer). Carry dual method-tagged bands; never encode "cancer=softer". 0.26 kPa is single-method.
- Nucleus E_nuc, cortical tension, membrane tension: MCF7 values are documented gaps/absences — generic-but-tagged, do not invent.
- col-I α2β1 rupture force is loading-rate dependent (38→90 pN over 180–8800 pN/s) — state the rate. Slip bond → KU-2.5 catch-peak gate must skip slip ligands.

## Net
2 config-only (R_cell, η — both already MCF7-anchored, shared with Layer-2) · 1 substantial
code task (wire col-I α2β1 slip into `fa.py`) · 1 PI-gated deferral (G/Tier-2) · 3
documented-absence flags (nucleus, cortical, membrane). Whole-cell modulus stays overlay.
Bundle R_cell + η + col-I as ONE coherent MCF7 preset + re-validation, not piecemeal.
