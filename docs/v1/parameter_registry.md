# Parameter registry — Option F Week 1

PI directive 2026-04-29 (Option F selected). Codex review item 10:
classify every model parameter into a literature-anchoring tier so
publication can quote each with a defensible provenance.

## Four tiers

| Tier | Definition | Publication policy |
|---|---|---|
| **L1 — Direct literature value** | Number quoted from a peer-reviewed source (preferably IF ≥ 15) within an order of magnitude. | Quote with citation. |
| **L2 — Derived from dimensional analysis** | Number set by combining other L1 parameters via a physics relation (e.g., CFL bound, capillary number, Péclet number). | Quote with derivation chain. |
| **L3 — PI-approved phenomenological placeholder** | Number set in a sanity_md decision with a documented rationale; not directly anchored but bounded by sensitivity analysis. | Quote with PI-decision reference + sensitivity caveat. |
| **L4 — Accepted limitation, not used for fitting** | Number set to a bound that is known not to match physical reality but used to keep the simulation tractable; never used in claims. | Mention only as "framework limitation"; do not quote in results. |

A parameter must belong to **exactly one** tier. Re-tiering across
revisions is allowed but must be logged in this file.

## Universal classification protocol

For every parameter, three questions:

1. **Source**: literature value, derivation, or sanity_md decision?
2. **Sensitivity**: bounded by an explicit sweep or argument?
3. **Use in publication**: quoted as a result-driver, or only as a
   numerical input?

---

## Layer 1 — Bulk hydrodynamics

| Param | Value | Tier | Anchor | Notes |
|---|---|---|---|---|
| K (cortical bulk) | 1 kPa | L1 | Fischer-Friedrich 2014 *Sci Rep* IF 4.4 (cortical tension series); Krieg 2008 *Nat Cell Biol* IF 17 (HEK/cortical mechanics) | mid-range value; sensitivity not gated |
| μ (shear) | 300 Pa | L1 | Same chain as K; μ/K ≈ 0.3 dimensional consistency (incompressible-soft tissue) | derived from K via L2 actually |
| τ_relax (Maxwell) | 60 s | L1 | **Moeendarbary 2013 *Nat Mater* IF 47** ("cytoplasm of living cells behaves as poroelastic material") | canonical anchor for the dimensionless time unit |
| ρ_density | 1050 kg/m³ | L1 | water + soft tissue; standard biophysics value | not result-driver |
| γ_cohesive | 1 mJ/m² | L1 | Foty & Steinberg 2005 *Dev Biol* IF 4.5 (tissue surface tension); dimensional with cell-cell adhesion energy | mid-range value |

**Layer 1 sub-decision**: v15 density-based volumetric stress
σ_vol = K·(ρ_ref/ρ_kernel − 1)·I

| Param | Value | Tier | Anchor |
|---|---|---|---|
| ρ_ref construction | harmonic mean over well-resolved population | L2 | derived from Hu et al. 2018 §4.3 reference state + Jensen's inequality fix (v3b sanity) |
| ρ_floor | 0.1 · ρ_ref | **L3** | `docs/v1/outcomes_v15.md` PI decision; conservative bound, no IF≥15 anchor; sensitivity tested by ρ_floor sweep documented in v15 sanity |

---

## Layer 1a+ — Substrate contact (Option β)

| Param | Value | Tier | Anchor |
|---|---|---|---|
| γ_sub_alpha | 1.0 | **L3** | `docs/v1/outcomes_stage1a_plus.md` PI decision; α=1.0 chosen after α∈{0.1, 0.3, 1.0} sweep gave best R drift. Original derivation hinted γ_sub ≈ γ_cc/2 from Col1 wettability literature but final value PI-tuned. |
| n_contact_band | 3 cells | L2 | dimensional; matches CSF colour-smoothing radius (Brackbill 1992 §V) |
| substrate enabled | true | n/a | model toggle |

**Note**: γ_sub_alpha = 1.0 is the closest current parameter to a
fitted value in the project. The Magic-Number-Block PARTIAL pattern
applies — sweep documented, PI authorized, but it would benefit from
an independent literature anchor (e.g., Beningo 2002 for Col1-coated
substrate cell adhesion strength) before publication.

---

## Layer 2 — Active boundary stress

| Param | Value | Tier | Anchor |
|---|---|---|---|
| ζ_star (zeta_star, scalar) | 0.0 (active stress is overridden by Layer 3 ζ(φ)) | L4 | placeholder; superseded by Layer 3 |
| ζ_min (Layer 3 controlled) | 0.1 | **L3** | sanity_md `stage1a_plus_plus_layer2_sanity.md`; chosen so ζ(φ_int=0) > 0 to keep interior cells active |
| ζ_max (Layer 3 controlled) | 0.6 | **L3** | sanity_md + Phase 2.2 sweep (`docs/v1/stage2_sanity.md` ζ ∈ [0.1, 1.0] 13-point response curve); ζ=0.6 is the first PASS of R drift gate |

---

## Layer 3 — Adhesion φ-ODE (Cho 2020)

| Param | Value | Tier | Anchor |
|---|---|---|---|
| φ_initial | 0.80 | L1 | Cho 2020 Western blot timecourse, MCF7 baseline E-cad fraction |
| k_+ | 1.39e-3 s⁻¹ | L1 | Cho 2020 Fig. 4 Slope (E-cad → Int-β1 transition) |
| k_- | 4.6e-4 s⁻¹ | L1 | Cho 2020 reverse rate; ratio k_+/k_- = 3 at literature value gives φ_eq = 0.75 |
| spatial S extension | boundary-only S=1, interior S=0 | **L3** → **HARD-BLOCKER** flagged | v11 sanity decision; **Codex review item 4** flags Layer 3 audit pending |

The "spatial S extension" is currently a HARD-BLOCKER in
`docs/v1/gate_fail_taxonomy.md` (F9). This is the single most fragile
Layer 3 design choice; pending Layer 3 audit.

---

## Layer 4 — Marangoni γ(φ)

| Param | Value | Tier | Anchor |
|---|---|---|---|
| γ_max_star | 0.020 | L1 → L3 in current implementation | Pajic-Lijakovic & Milivojevic 2022 *Eur Biophys J* IF 2.6 (Marangoni for cell spreading); dimensionalised against γ_cc |
| γ_min_star | 0.003 | L1 → L3 in current implementation | same source; γ_max/γ_min ≈ 6.7× (consistent with Cazabat-style binary mixture spread) |
| γ(φ) form | γ_max(1-φ) + γ_min·φ | L1 | linear interpolation, Pajic-Lijakovic & Milivojevic 2022 sign convention |
| layer3_spatial_S | true | L3 | inherits Layer 3 issue |

**Layer 4 has Mechanism A/E/F upgrades pending** per
`docs/v1/marangoni_review.md`. Their parameters (τ_γ, D_s, α_osm-γ
coupling) are not yet in the registry — will enter as L1 (Yadav
2022, Stone 1990) / L2 (CFL diffusion bound) when the upgrades are
implemented.

---

## Layer 5 — Mechano-osmotic Tier 2

| Param | Value | Tier | Anchor |
|---|---|---|---|
| ρ_osm_initial | 1.0 | L1 | dimensionless reference, Guo 2017 *PNAS* IF 11 |
| α_osm_star | 8.33e-3 | L1 | Guo 2017 + dimensional analysis; spreading-induced volume loss timescale |
| β_osm_star | 0.1 | L1 | Guo 2017 osmotic relaxation timescale; β/α ≈ 12 |
| ρ_osm_min | 0.5 | L4 | numerical clamp; not physical |
| ρ_osm_max | 1.6 | **L3** | sanity_md decision based on osmotic saturation literature; Codex review item 5 flags Layer 5 audit |

**Codex review item 5** flags Layer 5 audit needed: Guo 2017 is
single-cell; spheroid-level turgor with cell-cell water transport
(gap junctions, paracellular flow) is not present. Until audit,
Mechanism F upgrades (Layer 5 ↔ Layer 4 coupling) deferred.

---

## Layer 6 — Chemistry / ECM degradation

| Param | Value | Tier | Anchor |
|---|---|---|---|
| α_mmp_star | 1.0e-4 | L1 | Egeblad & Werb 2002 *Nat Rev Cancer* IF 78 (MMP secretion rates); dimensionalised |
| β_deg_star | 2.0e-3 | L1 | Egeblad & Werb 2002 (Col1 degradation kinetics) |
| ecm_strength_initial | 1.0 | L1 | dimensionless reference |
| ecm_strength_min | 0.1 | **L3** | sanity_md decision; conservative floor preventing γ_sub_eff = 0 (substrate detachment singularity) |

ecm_strength floor is hit at production scale; this is a documented
limitation (`docs/v1/production_lam4_finding.md`).

---

## Path C — Effective gravity

| Param | Value | Tier | Anchor |
|---|---|---|---|
| gravity_star | 0.01 | **L3** | `docs/v1/path_c_sanity.md` PI decision; chosen via dimensional argument g_star · H ≈ γ · κ at characteristic spread thickness |

This is the canonical Hard Rule 10 (dimensional comparison
verification) anti-pattern caught at the 2026-04-29 Path C g_star
investigation: original sanity-md proposed g_star=0.1 ("same order as
γ_star=0.01") but failed the per-volume-vs-per-area unit check; the
corrected comparison g_star · H vs γ · κ gave 7× over-anchoring,
confirmed by the resulting pancake regime. Final value 0.01.

---

## Numerics

| Param | Value | Tier | Anchor |
|---|---|---|---|
| dt_star | 0.01 | L2 | overdamped MLS-MPM stability bound; dimensional with τ_relax |
| dx_um | 4.0 | L2 | grid_n = 64 → dx = domain/grid_n = 600/64 ≈ 9.4 μm; dx_um=4.0 used for spheroid radius scale (R₀=100μm → R₀/dx = 25 cells) |
| frame_interval_star | 15.0 | L2 | imaging-cadence-matched (PI experimental 15 min imaging) |
| domain_star | 6.0 | L2 | spheroid R₀=1 (star), domain 6× wider per side, prevents grid boundary touching |
| free_surface_density_threshold | 0.6 | L1 | Brackbill 1992 / Hu 2018 (kernel density boundary criterion) |
| n_material_points | 1000 (pilot) / 5000 (production) | L2 | VRAM-tier sizing per CLAUDE.md performance budget |

---

## Summary statistics

- **L1 (direct literature)**: 22 parameters
- **L2 (dimensional)**: 9 parameters
- **L3 (PI-approved phenomenological)**: 11 parameters
- **L4 (accepted limitation)**: 3 parameters

**L3 count is the audit target**. 11 phenomenological placeholders
across the framework. Each has either a sweep or a sanity_md decision,
but only 4 of them have IF ≥ 15 indirect anchors:

| L3 param | IF ≥ 15 indirect anchor exists? | Audit priority |
|---|---|---|
| ρ_floor (v15) | partial (Hu 2018 §4.3 + Jensen) | LOW |
| γ_sub_alpha | NO (Beningo 2002 IF≈10) | HIGH (Codex item 10) |
| ζ_min, ζ_max | partial (sweep + sanity) | LOW |
| spatial S extension | NO (audit pending) | **HIGH (Codex item 4)** |
| γ_max, γ_min (Layer 4) | partial (Pajic-Lijakovic IF 2.6) | MODERATE |
| ρ_osm_max | NO (Guo 2017 single-cell only) | HIGH (Codex item 5) |
| ecm_strength_min | NO (numerical clamp) | LOW |
| gravity_star | dimensional (Path C) | LOW |

Three HIGH-priority audits: γ_sub_alpha (substrate adhesion strength),
spatial S extension (Layer 3 audit), ρ_osm_max (Layer 5 audit). These
align with Codex review items 2/4/5 already in Option F Week 2 scope.

---

## Cross-references

- `docs/codex_review_synthesis.md` — Codex review item 10
- `docs/v1/gate_fail_taxonomy.md` — gate FAIL classification
- `docs/v1/marangoni_review.md` — Mechanism A/E/F upgrade parameters
- `docs/12_validation.md` — known limitations
- `docs/v1/outcomes_v15.md` — ρ_floor sweep
- `docs/v1/stage1a_plus_substrate_sanity.md` — γ_sub_alpha sweep
- `docs/v1/path_c_sanity.md` — gravity_star derivation
