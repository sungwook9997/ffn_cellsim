# SMOKE-ONLY Constants Provenance (2026-06-09)

Every constant the ENABLED-PATH smoke harnesses (`scripts/compartment_smoke/`)
used for a None-gated / PI-pending parameter, with its
`COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09` confidence tag, citation, and
ratification status. This is a single **PI ratification checklist**: a value here
makes the *smoke plumbing* run; it does NOT make the *production* path runnable
(production keeps the constant `None` → the enabled build raises, guarded by
`tests/test_compartment_unratified_raises.py`).

Confidence tags: **SOLID** = literature-anchored, in-band, not pending ·
**DERIVED** = computed from anchored inputs/geometry (grid-invariant) ·
**ORDER** = order-estimate, PI to ratify · **LAYOUT** = discretisation/layout knob
(grid-invariant, free) · **NONE/blocked** = deliberately not supplied (path raises).

## microtubules
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| EI | 2.2e-23 N·m² | SOLID | Gittes 1993 JCB 120:923 (module default) | no |
| L_p | 5.2e-3 m | SOLID | Gittes 1993 | no |
| Y_stretch | 2.3e-7 N | ORDER | Kis 2002 PRL / Pampaloni 2006 (E·A; the dt lever) | **yes** |
| L_mt | 5.0e-6 m | ORDER | interphase aster few-µm (Howard 2001); host geometry | **yes** |
| n_mt | 20 | ORDER | interphase centrosomal count (no MCF7 KU) | **yes** |
| beads_per_mt | 25 | LAYOUT | Dmitrieff 125 nm segmentation; grid-invariant | no |

## intermediate_filaments
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| E_if | 6.0e6 Pa | ORDER | Kreplak 2005 / Guo 2013 few-MPa (module default) | **yes** (keratin pick) |
| Lp | 0.5e-6 m | SOLID | Lichtenstern 2012 K8/K18 (module default) | no |
| d_if | 10e-9 m | SOLID | Herrmann/Aebi canonical assembled IF | no |
| max_stretch_ratio | 3.0 | SOLID | Kreplak 2005 (2-3.5×) | no |
| nonlinear Table curve | — | NONE | NOT supplied → `register_if_bond_params(nonlinear=True)` raises | **yes** |
| n_filaments / beads_per_fil | 60 / 8 | LAYOUT | layout knobs | no |

## linc
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| k_linc | 1.0e-2 N/m | ORDER | Rief 1999 folded-rod pre-unfolding secant (route A) | **yes** |
| r0 | 50e-9 m | SOLID | Crisp 2006 perinuclear span | no |
| f_rest (oracle) | 8.0e-12 N | SOLID | Déjardin 2020 (oracle; not a smoke *input*) | provenance |

## osmotic_regulation
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| Lp | 1.0e-12 m/(s·Pa) | ORDER | Jung 2011 MCF7/AQP5 (default 1e-13 likely wrong) | **yes** |
| Δc (shock) | -100 mol/m³ | ORDER | Jung 2011 100 mM sorbitol protocol | run-knob |
| turgor_dP0 | 133 Pa | SOLID | physiological baseline (registry; LIVE) | no |
| c_phys | 300 mol/m³ | SOLID | Lodish/Alberts mammalian osmolarity | no |

## cadherin_junction
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| k_trans | 2.92e-3 N/m | DERIVED | f0/contact_zone (Rakshit 2012 f0=29.2 pN / 10 nm) | no (overridable) |
| r_bind / r0_trans | 1.0e-8 m | ORDER | EC1 strand-swap reach ≈ contact zone | no (default) |
| n_cad_per_cell | 40 (smoke) | LAYOUT | smoke count — NOT the Iturri-2020 N_cad=223 bridge | n/a |
| Iturri SourceEvidence | — | NONE | unregistered in KB; register before any deliverable cite | **yes** |

## ventral_stress_fibers
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| EA_single | 4.3e-8 N | SOLID | Gittes 1993 / Kojima 1994 (audit-OK) | no |
| N_filaments | 20 | ORDER | Cramer 1997 (10-30) → μ_SF = N·EA_single | **yes** |
| k_actin | μ_SF/ell0 | DERIVED | grid-invariant backbone (for measure consistency) | follows N_filaments |
| sf_myosin_* prefix | — | NONE | NMII off — BLOCKER #2 (shares cortex_myosin_*) | **yes** |
| n_SF / n_beads_per_SF | 6 / 24 | LAYOUT | layout knobs | no |

## membrane_reservoir
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| W_MCA | 1.0e-5 J/m² | SOLID | KU-3.B1.4 (Hochmuth 1996 / Derényi 2002) band-mid | no |
| k_tether | 0.1 N/m | DERIVED | ERM radial pinning order (Charras 2008) | no |
| max_tether_dist | 200e-9 m | SOLID | KU-3.17 cortex thickness | no |
| σ_crit_bleb | — | NONE | NOT supplied → bleb updater raises | **yes** |
| f_excess | — | NONE | NOT supplied → reservoir release blocked | **yes** |

## junctional_actin (scalar only; HOOMD build is a STUB)
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| x_catch | 4.0e-9 m | SOLID | Buckley 2014 / KB-4.17 catch arm | no |
| x_slip | 0.4e-9 m | SOLID | Buckley 2014 / KB-4.17 slip arm | no |
| k_catch0 | 1.0 /s | ORDER | matched to Buckley τ(0)~1 s | **yes** |
| k_slip0 | 0.02 /s | ORDER | Pereverzev structure (k_slip0 ≪ k_catch0) | **yes** |
| k_couple / k_anchor / k_on / max_couple_dist / anchor_r0 | 1e-6 / 1e-5 / 10 / 6e-8 / 5e-9 | DERIVED | sanctioned H.3 transfers | confirm |
| HOOMD build path | — | NONE | reserved STUB → `extend_*` raises even anchored | **yes** |

## surface_manifold (geometry-only)
| const | value | tag | provenance | PI-pending? |
|---|---|---|---|---|
| R_cell | 7.5e-6 m | SOLID | MCF7 (Wagner 2011) | no |
| reach | 1.0e-6 m | LAYOUT | broad-phase candidate reach (geometry knob) | no |
| k_conf (soft normal confinement) | — | NONE | NOT used (no force) — Magic-Number Block [3e-6,3e-5] N/m PI-gated | **yes** (if a force is ever added) |

---
*Companion: `PLATFORM_PI_QUEUE.md` (blockers), `COMPARTMENT_SMOKE_REPORT_2026-06-09.md`
(results), `COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09.md` (full candidate rationale).*
