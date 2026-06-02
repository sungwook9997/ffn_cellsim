# Cancer cell-type → simulation-parameter map (MCF10A / MCF7 / MDA-MB-231)

> **STATUS: REFERENCE map — 2026-06-02. Prep, NOT a contract, NO code.** Parallel
> session (`AUTONOMOUS_LOG_2026-06-02.md`). This is the *how-to-parameterize* layer —
> "verified band → which sim knob → value, per cell-type" — **distinct** from the
> Notion KU v2 KnowledgeClaim rows (which hold the claims themselves; cross-referenced,
> not duplicated). Cell-type *contracts* are deliberately deferred until cell-type
> instantiation (Notion 2026-06-01); this map is the literature consolidation that
> instantiation will draw on. Companion to `H10_CYTOPLASM_DESIGN.md`,
> `LIGAND_IDENTITY_FA_SPEC.md`, the EXTEND memo, and the PI-exp validation map.

## Governing principle — there is NO global "cancer modulus"

Cell-type identity is **per-compartment**, not a single stiffness scalar. The platform
is compartmentalized (H.3 cortex, H.8 membrane, H.9 nucleus, H.10 cytoplasm), so a
cell-type is a *vector* of per-compartment parameters, each from its own measurement.
The **robust discriminator is cytoplasm viscosity** (MCF7 ≈5× MDA — clean, 3-0
verified); whole-cell stiffness *ordering* is **method-dependent and disputed** and
must NOT become a "cancer = softer" knob (see Contested below). Parameterize from the
robust quantities; carry the disputed ones as dual-method overlays only.

## The map (verified bands → sim knob → per-cell-type value)

All values are **acceptance/overlay targets, never fitting targets** (CLAUDE.md). SI.

### Cytoplasm viscosity η — H.10 Tier-1, `gamma_map` per cytoplasm-immersed type

`γ_b = 6π η_eff R` (the H.10 Tier-1 design). This is the **primary cell-type
discriminator** and the highest-confidence axis.

| Cell type | η (MRS, L≈3 µm) | η (Hu 2024 refined) | sim knob |
|---|---|---|---|
| MCF10A (normal) | 41.6 Pa·s | — | `η_eff` for interior types |
| **MCF7** (low-invasion) | 56.4 Pa·s | **65.9 ± 11.4 Pa·s** | `η_eff` |
| **MDA-MB-231** (high-invasion) | 10.7 Pa·s | **12.0 ± 5.7 Pa·s** | `η_eff` |
| **contrast** | MCF7 ≈ **5.3×** MDA | MCF7 ≈ **5.5×** MDA | the validation target |

Sources: H.10 brief / KU-3.B3.1; Hu 2024 *Nanoscale Adv* PMC10929591 (MRS). State the
probe length L (η ∝ L²). **Validation = reproduce the ~5× MCF7/MDA viscosity ratio**,
not the absolute Pa·s (which is L-dependent).

### Cytoplasm shear modulus G — H.10 Tier-2 (GLE storage; PI-gated)

| Cell type | intracellular G (Pa) | note |
|---|---|---|
| MCF10A | 79.3 | the elastic storage the GLE Prony series targets |
| MCF7 | 32.9 | ≈10–20× below whole-cell AFM E₀ (distinct quantity) |
| MDA-MB-231 | 38.6 | band 33–79 Pa |

Source: Hu 2024 PMC10929591. **Note the ordering differs from η** (MCF7 has the
*highest* η but the *lowest* G) — they are independent compartmental properties, a good
guard against collapsing "cancer mechanics" to one number.

### Nucleus — H.9 `NucleusConfinement` k_nuc (lamin-A/C set)

| Quantity | band | sim knob | note |
|---|---|---|---|
| Nuclear modulus E_nuc | 1–10 kPa (in-situ ~5) | `k_nuc` shell stiffness | KU-3.B2.1 |
| Lamin-A/C trend (cancer) | lamin-A/C ↓ → softer / more deformable → invasion ↑ | scale `k_nuc` down for high-invasion | context-dependent |

**GAP (cheap lit follow-up):** MCF7- vs MDA-specific lamin-A/C levels / E_nuc are not
yet anchored here — fill from literature before instantiating the nuclear axis (do NOT
invent). Trend only for now: MDA (high-invasion) ⇒ lower `k_nuc` than MCF7.

### Cortex (H.3) + membrane (H.8) — existing knobs

| Quantity | band | sim knob | source |
|---|---|---|---|
| Cortical tension | ~10⁻² N/m (10⁻²–10⁻¹), short-pulse | H.3 γ_cortex | PMC10625872 |
| Membrane surface tension γ_mem | KU-3.B1 (T 0.03–0.3 mN/m) | H.8 `MembraneSurfaceTension` | existing module |
| Cancer membrane overlay | metastatic ~2× lower apparent T | γ_mem / γ_MCA | H.8 brief KU-3.B1.5 |

### Whole-cell AFM modulus E₀ — OVERLAY ONLY (composite, not a single knob)

| Cell type | E₀(1 s) | fluidity exp β | note |
|---|---|---|---|
| MCF10A | 1.14 kPa | 0.186 | composite cortex+membrane+nucleus+cytoplasm |
| MCF7 | 0.26 kPa | 0.234 | NOT an input — an *emergent* whole-cell readout |
| MDA-MB-231 | 0.46 kPa | 0.147 | compare to the simulated AFM-indentation response |

Source: Yubero 2020 *Commun Biol* s42003-020-01330-4. **E₀ is an output to validate
against, never a knob** — it emerges from the compartment parameters above.

## Contested — do NOT encode as a single ordering

⚠️ **MCF7-vs-MDA whole-cell stiffness *ordering* is method-dependent:** AFM power-law →
MCF7 *softer* (0.26 vs 0.46 kPa); electrodeformation → MCF7 *~10× stiffer* (7.1 vs
0.7 kPa). Carry **both** as method-tagged dual bands; never collapse to "cancer is
softer/stiffer." The robust separators are **cytoplasm viscosity** and contractility,
not whole-cell modulus. (PI-exp map §Contested; matches memory "no global cancer modulus".)

## How a cell-type preset assembles (illustrative — NOT a config schema)

A cell-type = a vector over the compartment knobs, default-off-compatible (omit a
compartment ⇒ that module stays off / at its current value):

```
cell_type: MCF7
  cytoplasm:  { eta_eff: 65.9 Pa·s }     # H.10 Tier-1 (PI-gated magnitude)
  nucleus:    { k_nuc: <from E_nuc, lamin> }  # H.9 (GAP: MCF7 lamin)
  cortex:     { gamma_cortex: <H.3 band> }
  membrane:   { gamma_mem: <H.8 KU-3.B1> }
  # ligand identity is a SUBSTRATE-condition axis (LIGAND_IDENTITY_FA_SPEC), not a
  # cell-type property — keep the two axes orthogonal.
```

## Validation framing

Cell-type validation = **reproduce the per-compartment contrasts** (above all the ~5×
MCF7/MDA cytoplasm-viscosity ratio), with the PI poster + whole-cell E₀ as **overlays
at comparison time only**. A preset that needs a PI value to pass = a violation
(literature-first). Keep cell-type (modulus vector) and substrate-condition (ligand,
ECM) as orthogonal axes.

## Open items for PI

- [ ] Ratify the **per-compartment (no global modulus)** parameterization.
- [ ] Fill the **MCF7/MDA nuclear (lamin-A/C, E_nuc)** gap from literature (cheap).
- [ ] Confirm η probe-length convention (L≈3 µm) for the viscosity axis (η ∝ L²).
- [ ] Method-tag the contested whole-cell stiffness ordering as a dual band.
- [ ] Sequence cell-type instantiation after the H.9/H.10 modules are enabled (the
      knobs they expose) — this map is the input, not a trigger.
