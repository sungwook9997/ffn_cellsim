# SourceEvidence / KnowledgeClaim candidate — cytoplasm Biot coupling modulus M (FSI)

**Date:** 2026-07-16 · **Origin:** PI decision (B) — use the physically-motivated undrained-derived Biot
coupling modulus for the FSI fluid→compartment coupling, with a KB citation (not the conflated K_drained).
**Status:** DRAFT — not yet Notion rows. PI sign-off required before `harvest_ops`/registration.

## Why this exists

The spatial Biot FSI (`network_warp.py` `biot_fsi`) needs a **coupling modulus M** = the pore pressure
generated per unit volumetric strain in the undrained limit (`p ≈ M·α·ε`, α≈1). The first landing defaulted
`M_fsi = K_drained = 300 Pa`, which a 17-agent verification workflow (wf_9e422e1c) flagged as a
**conflation**: K_drained is the DRAINED skeleton stiffness (already carried by the Terzaghi `dP_solid`
channel), a physically distinct quantity from the Biot modulus M. This registers the sourced value for M.

## Physical basis (Terzaghi/Biot)

Under undrained loading the total stress splits into effective (skeleton) + pore pressure:
`σ = K_u·ε`, `σ' = K_drained·ε`, so the pore pressure `p = σ − σ' = (K_u − K_drained)·ε`. Hence the
pore-pressure coefficient **M ≈ K_undrained − K_drained** (α≈1 for a saturated soft cell).

## Sourced anchors

- **Drained modulus** (KB-3.B3.2, `Moeendarbary2013_NatMater`, DOI 10.1038/nmat3517, PMC3925878): HeLa
  ≈ 0.9 kPa, HT1080 ≈ 0.4 kPa, MDCK ≈ 0.4 kPa (equilibrium/long-time, cross-cell — flagged NOT breast-specific).
  Poroelastic diffusivity Dp ≈ 40–60 µm²/s, mesh ξ ≈ 14 nm, poroelasticity dominates for events < ~0.5 s.
- **Undrained (instantaneous) modulus** — NOT in the KB (confirmed via TAG); would be pulled from the paper.
  Cell poroelastic **undrained/drained ratio ~2–3** (standard for the instantaneous-vs-equilibrium AFM
  response). → K_u ~ 1–2.7 kPa → **M ~ 0.6–1.8 kPa**.

## Value (draft)

**M_biot = 1.0 kPa (1000 Pa)** for the MCF7 cytoplasm FSI coupling — the mid of the undrained-derived range,
consistent with K_drained~300 Pa (code) and K_u~1–1.3 kPa. Used as the `biot_fsi` `M_biot_Pa` default.
**Never tuned to a nucleus-response band** (no-param-tuning HARD rule).

## Honest caveat (report with the value)

Even at this physical M, the **fluid-mediated interior (nucleus) deformation is small (~sub-nm)**: the
*angular* (quadrupole) pore-pressure deviation the nucleus surface feels is ~Pa (fast diffusion equilibrates
it; the nucleus sits deep where the field is smooth), and the nucleus shell + incompressibility (K_vol~66 kPa)
is stiff. The fluid genuinely COUPLES to every compartment now (the CFD requirement), but it is NOT the
dominant nucleus-deforming force — that is DIRECT cytoskeletal coupling (IF cage / LINC), the next build.

## Registration target

New SourceEvidence pointer (Moeendarbary 2013 already SE'd via KB-3.B3.2) + a new KnowledgeClaim
"cytoplasm Biot coupling modulus M ≈ 1 kPa (undrained-derived)" → ModelContract for the FSI coupling.
PI sign-off required.
