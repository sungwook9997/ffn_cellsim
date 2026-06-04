# SourceEvidence registration candidates — 2026-06-04 (collective-migration / SPP plithotaxis)

Anchors for the **Layer-2 §E coherent-collective-traction (SPP plithotaxis)** mechanism
(`outputs/layer2/DESIGN_plithotaxis.md`). Sourced by the deep-research pull `wf_e0df20b3-470`
(106 agents, 23 sources → 104 claims → 25 adversarially verified, 20 confirmed / 5 refuted).
All values below are **3-0 confirmed** (verbatim-quote-checked). **None is in the references
corpus yet** (no PDF, no BM25, no Notion SE row) → register the PDFs + SE rows before these
become citeable. **Overlay-only HARD rule:** every number anchors a model PARAMETER or is an
emergent-validation TARGET; the PI A/A₀ is never fitted.

> ⚠️ **Citation-integrity (2026-06-02 audit: 3 hallucinations found).** Verify each DOI on
> CrossRef before creating the SE row. Two integrity flags carried from the verification pass
> are marked ⚠️ below. SourceEvidence is edited in Notion (SoT), never here.

> ⚠️ **No MCF7-specific data exists** for any of these observables (deep-research confirmed).
> The closest breast-epithelial proxy is **MCF10A** (Smeets 2016); everything else is MDCK /
> HBEC / NRK / endothelial / generic SPP-SPV theory. Each row states its system. The PI's own
> MCF7-spheroid data remains the overlay target; these proxies anchor the *mechanism*.

---

## A. SPP / CIL parameters — the DIRECTLY-ANCHORED model inputs (register first)

| paper | DOI | system | value → model param | verdict |
|---|---|---|---|---|
| **Smeets, Alert, Pešek, Pagonabarraga, Ramon, Vincent 2016**, *Emergent structures and dynamics of cell colonies by contact inhibition of locomotion*, PNAS 113(51):14621-14626 | 10.1073/pnas.1521151113 | **MCF10A** (breast epithelial, MSD fit) | **D_r ≈ 0.05 min⁻¹** (→ τ=1/D_r=20 min, polarity persistence); **v_m ≈ 1 µm/min** (self-propulsion v0); **f_cil ≈ 0.1 min⁻¹** (CIL repolarization rate); **ψ=f_cil/(2D_r) ≈ 1** | 3-0 ✓ |
| ↑ same — CIL formalism | — | generic SPP | polarity rule **θ̇ᵢ = −f_cil(θᵢ − θᵢ^free) + √(2D_r)·η**, θ^free = away from weighted mean position of contacting cells; NP_N>1 ⇒ global outward polarity | 3-0 ✓ |
| **Bi, Yang, Marchetti, Manning 2016**, *Motility-driven glass and jamming transitions in biological tissues*, Phys. Rev. X 6:021011 | 10.1103/PhysRevX.6.021011 | SPV theory | **τ = 1/D_r** (white angular noise, definitional); triplet {v0, p0, D_r}; Pe~v0/D_r; p0*=3.81 | 3-0 ✓ |

**Provenance caveats (record on the Smeets row):** only `v_m` and `D_r` are MCF10A-direct MSD
fits; `f_cil≈0.1/min` is borrowed from hemocytes/fibroblasts/MDCK (no MCF10A CIL-contact
measurement); `ψ≈1` assumes MCF10A≈MDCK similarity. Still the closest breast-epithelial anchor.

## B. Velocity correlation length ξ — EMERGENT-VALIDATION target (not a fitted input)

| paper | DOI | system | value | verdict |
|---|---|---|---|---|
| **Petitjean, Reffay, … Silberzan 2010**, *Velocity fields in a collectively migrating epithelium*, Biophys J 98(9):1790-1800 | 10.1016/j.bpj.2010.01.030 | **MDCK** (PIV) | **ξ ≈ 200 µm** collective epithelial; **≈40 µm** NRK fibroblast contrast; fingers influence field over ~200 µm | 3-0 ✓ |
| **Garcia, Hannezo, … Gov 2015**, *Physics of active jamming during collective cellular motion in a monolayer*, PNAS 112(50):15314-15319 | 10.1073/pnas.1510973112 | **HBEC** (+MDCK,3T3) | ξvv **NOT** density-controlled (\|Pearson\|<0.3); controlled by velocity/junction maturation; epithelial polarity τ ~ **30 min**; **velocity correlations (~200µm) EMERGE from per-cell traction persistence WITHOUT explicit Vicsek alignment** | 3-0 ✓ |
| **Angelini, Hannezo, … Weitz 2011**, *Glass-like dynamics of collective cell migration*, PNAS 108(12):4714-4719 | 10.1073/pnas.1010059108 | **MDCK** | dynamic-heterogeneity ξh grows ~10→30 cell-bodies with density, declines past glass σ_g=2800 cells/mm² | 3-0 ✓ |
| **Henkes, Sknepnek, Marchetti 2020**, Nat Commun (corroborating) | (verify) | SPP theory | persistent uncoordinated motility + collective elastic modes ⇒ swirl/velocity correlations, **no Vicsek alignment needed** | corroborating |

**Decisive for the mechanism:** Garcia 2015 + Henkes 2020 ⇒ ξ~200µm is an **emergent** property
of persistence + cohesion, **not** an imposed alignment length. So ξ is a **diagnostic to
measure and overlay-validate**, NOT a parameter to fit a Vicsek J to. ⚠️ ξ has 3 operationally
different definitions (Petitjean PIV / Angelini ξh / Garcia ξvv) — do not conflate; the model's
emergent-ξ estimator must match whichever it is compared to (use Petitjean PIV-style).

## C. Plithotaxis / monolayer stress — OVERLAY-VALIDATION targets (emergent stress field)

| paper | DOI | system | value | verdict |
|---|---|---|---|---|
| **Tambe, Hardin, … Trepat, Fredberg 2011**, *Collective cell guidance by cooperative intercellular forces*, Nat Mater 10:469-475 | 10.1038/nmat3025 | RPME/MDCK/**MCF10A** | intercellular stress predominantly **tensile, >300 Pa over tens of cells**; **tension builds cumulatively edge→inward (tug-of-war)**; plithotaxis = migrate along **max-principal = min-shear** stress | 3-0 ✓ |
| **Trepat, Wasserman, … Fredberg 2009**, *Physical forces during collective cell migration*, Nat Phys 5:426-430 | 10.1038/nphys1269 | MDCK | ORIGINATOR of tug-of-war: traction forces pile up many rows behind the leading 5-10 rows | 3-0 ✓ |
| **Tambe et al. 2013**, *Monolayer Stress Microscopy: Limitations…*, PLoS One | 10.1371/journal.pone.0153471 (verify) | validation | hundreds-of-Pa stresses genuine (R≥0.66); plithotaxis robust to artifacts | 3-0 ✓ |

⚠️ **Citation-integrity flag (carried from verification):** the verbatim *"average tensile stress
increased systematically with increasing distance from the advancing front… global tug-of-war"*
sentence is **Tambe et al. 2011 Nat Mater**, NOT Trepat 2009 Nat Phys (Trepat originated the
*concept* and co-authors Tambe 2011, but the quoted sentence is the 2011 paper). 300 Pa is a
regional/peak tensile value over a heterogeneous field, not a uniform monolayer tension.

## D. Spreading-front kinematics — OVERLAY-VALIDATION targets

| paper | DOI | system | value | verdict |
|---|---|---|---|---|
| **Poujade, … Silberzan 2007**, *Collective migration of an epithelial monolayer in response to a model wound*, PNAS 104(41):15988-15993 | 10.1073/pnas.0705062104 | **MDCK** | front velocity accelerates **0→10±5 µm/h over ~15 h**, position ⟨s⟩~t^**1.8±0.4** (super-linear, NOT constant-velocity); **leader cells 18±2 µm/h constant**; **~70% fingers normal/outward**; ≤5 fingers/mm | 3-0 ✓ |

⚠️ **Do NOT use** (refuted 0-3 / 1-2 on adversarial review): ξ~100µm inside-MDCK; ξ~14µm
one-cell-diameter; B16F10 ~0.5-0.7 min autocorrelation; t^1.29/t^1.02 colony-radius growth.

---

## How these map into the model (DESIGN_plithotaxis.md §3)

| model param | anchored value | source | role |
|---|---|---|---|
| `D_r` (rotational diffusion) | **0.05 min⁻¹** = 8.33e-4 s⁻¹ | Smeets 2016 (MCF10A) | persistence τ=20 min input |
| `f_cil` (CIL repolarization rate) | **0.1 min⁻¹** = 1.67e-3 s⁻¹ | Smeets 2016 | free-edge outward drive |
| `f_active` self-propulsion arms | **1.6 / 5.0 / 9.4 nN** | motility_bridge (1.6 net MCF7, 9.4 protrusion) + **5.0 = v_m 1µm/min·γ_cell (MCF10A Smeets)** | the 3 bracket arms; 5.0 nN = the principled breast-epithelial self-propulsion anchor |
| Vicsek alignment J | **OMITTED (default-off)** | Garcia 2015 / Henkes 2020 | correlations emerge from persistence+cohesion; imposing J would be an unanchored tuned term (anti-magic-number) |
| ξ emergent | **~200 µm target** (overlay) | Petitjean 2010 | self-consistency diagnostic, NOT fitted |
| intercellular stress | **>300 Pa tensile, edge→in build-up** (overlay) | Tambe 2011 | emergent stress-field validation |
| ψ = f_cil/(2 D_r) | **= 0.1/(2·0.05) = 1.0** ✓ reproduces Smeets MCF10A | Smeets 2016 | derived consistency check (no free param) |

**Net:** every model INPUT (D_r, f_cil, f_active) is a measured value; ξ and stress are emergent
TARGETS; the Vicsek term is dropped on literature grounds (reduces free parameters, anti-magic-
number). v_m·γ_cell = 5.0 nN sits below cohesion 6.5 nN ⇒ admissible even with the brittle bond.
