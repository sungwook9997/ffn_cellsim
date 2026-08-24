# NMII production parameter GAP — PI decision packet (2026-07-21)

Status: **PI DECISION RECEIVED (2026-07-21): first-baseline paralog = pure NM2B (fewest-gap direction).**
Production magnitudes remain PARTIALLY blocked: the NM2B-sourced fields (`L_bb`, `N_side`, `duty`, `k_off0`,
`d_step`) now have a selected direct source, but value insertion into `params_i0b3.yaml` still waits on the
SourceEvidence → KnowledgeClaim → Parameter audit chain (step 5 below). The paralog-independent GAPs
(`F_stall_head`, `v0`, `kappa_hill`, `k_xb`, `k_on`, `r0_head`, backbone `L_p`, cortical minifilament density)
are unchanged. No value was selected to improve NG-1, cortical tension, or implicit convergence.

## PI DECISION 2026-07-21 — pure NM2B baseline

PI ratified NM2B because it leaves the fewest open GAPs: Nagy 2013 directly measures its duty, detachment/ADP
rate, power stroke, and dimers-per-side at the single-molecule level, and Melli 2018 gives its backbone length.
NM2A lacks single-molecule `k_off0`/`d_step` assays, and a heterotypic mix would additionally require sourced
population weights — both strictly more open than pure NM2B. MCF7 does express multiple isoforms (Dey 2017), so
this is a defensible **first**-baseline simplification to be revisited, not a claim that the MCF7 cortex is
NM2B-only.

| Field | NM2B ratified value | Direct source | Runtime status |
|---|---|---|---|
| `L_bb` | 0.323 µm (323±24 nm) | Melli 2018 (Billington meas.) | select pending SE-param audit |
| `N_side` | ~14 dimers/side (≈28 heads) | Nagy 2013 | select pending SE-param audit; retire AFINES-10 |
| `duty` | 0.20–0.25 | Nagy 2013 | cross-check only (engaged fraction stays emergent) |
| `k_off0` | 0.4 s⁻¹ detachment (NOT the 0.35 ADP-release) | Nagy 2013 | select pending SE-param audit; fixes the YAML mislabel |
| `d_step` | ~6 nm | Nagy 2013 | select pending SE-param audit; retire 11 nm muscle value |

Still GAP under NM2B (needs new sourcing/experiment, unaffected by the paralog pick): `F_stall_head`, `v0`,
`kappa_hill`, `k_xb`, head `k_on`, `r0_head`, backbone `L_p` (Adamovic/Kaufmann are proxies only), NM2B central
bare-zone geometry (Liu 2017 measured NM2A only), and MCF7 cortical minifilament density.

## Contract-Graph/TAG result

Three TAG queries were run against the refreshed DuckDB read layer.

- No Contract-Graph record quantifies MCF7 NMIIA/NMIIB absolute abundance or cortical minifilament density.
- No registered claim/parameter quantifies mature NMII minifilament molecule count, head count per bipolar
  half, contour length, bare-zone geometry, head reach, or mature-backbone flexural rigidity.
- Single-head mechanochemical values occur only in `KB-DRAFT-*` claims whose defining Stam/Kovács/Finer/Hill
  sources are absent from `SourceEvidence`; aggregate motor-clutch values are verified but are not admissible
  substitutes for a resolved NMII head.
- `KB-3.18` must not be promoted as a clean single-head source without edge cleanup: the project-wide citation
  audit records a hallucinated `Yao2011_NatCommun` edge into KB-3.18 even though other attached sources may be
  real. Claim-level provenance is therefore mixed.

## Primary-source candidates

| Runtime field | Direct evidence found | Applicability | PI action |
|---|---|---|---|
| `N_side` | Human NM2A/NM2B bipolar filaments contain about 30 molecules total; Melli et al. also state about 30 motors are available at one filament end. Nagy et al. describe NMIIB as about 14 dimers per side. | Direct human NM2 topology, isoform-specific but not MCF7 composition-specific. Supports a resolved head population near 28–30 heads per side, not the AFINES 10-head surrogate. | Ratify NM2A, NM2B, or heterotypic baseline before selecting the exact count. |
| `L_bb` | NM2A `301 ± 24 nm`; NM2B `323 ± 24 nm` in the Billington measurements summarized by Melli and Liu. | Direct purified human paralog filament geometry. | Choose paralog/composition; 0.301 µm is defensible only for an NM2A baseline. |
| bare-zone topology | Liu et al.: unphosphorylated NM2A bare zone `167.9 ± 29.7 nm`; phosphorylated NM2A `149.4 ± 18.8 nm`; ATP changes unphosphorylated bare-zone length. | Direct but biochemical-state dependent. The current no-bare-zone topology is incomplete. | Specify phosphorylation/ATP state, then add an explicit central bare zone; do not infer it from total length. |
| `duty` | Human NMIIB duty ratio about `0.20–0.25`; Melli reports NM2A single-molecule duty ratio `0.05`. | Direct paralog-specific measurements; not interchangeable. | Select paralog mixture. Reject the generic 0.1 placeholder. |
| `k_off0` | Human NMIIB detachment about `0.4 s^-1` at physiological 1 mM ATP; ADP release about `0.35 s^-1`. | Direct NMIIB assay. Does not establish NM2A or a single MCF7 mixture value. | May ratify for a pure-NM2B contract only; otherwise keep GAP. |
| `d_step` | Human NMIIB single- and double-headed HMM power stroke about `6 nm`. | Direct NMIIB optical-trap result. | May ratify for pure NM2B; do not use the generic 11 nm muscle value. |
| MCF7 isoform role | Dey et al. detect NMIIA/B/C in MCF7 and show longer cortical high-FRET dwell for NMIIB plus differential filament formation after cortex ablation. | Direct MCF7, but relative/phenotypic rather than absolute abundance or density. | Use to require a multi-paralog or explicitly chosen paralog contract; it cannot set copy number. |
| `v0` | Aggregate motor-clutch 100–120 nm/s is KB-verified; NM2A/B filament motility is paralog- and viscosity-dependent. | Aggregate clutch value is the wrong primitive; current single-head candidates are not closed in the KB. | Keep GAP until assay temperature, ATP, viscosity, and paralog are fixed. |
| `F_stall_head` | Aggregate clutch 2 pN is verified; primary NM2 filament studies explicitly note that individual nonmuscle filament force output remains unresolved. | Aggregate value cannot be copied into an explicit head. | Keep GAP; do not use the 2 pN clutch value as a head measurement. |
| `kappa_hill` | No numeric human NM2 Hill curvature found. | Missing for the mandated per-head Hill law. | PI must choose a documented proxy/form contract; keep magnitude GAP. |
| `k_xb` | No NMII crossbridge stiffness record in the Contract-Graph or direct human NM2 source set. | Missing master force/compliance parameter. | Keep GAP; 100–1000 pN/µm remains only a broad crossbridge-class plausibility band. |
| `k_on` | No dedicated NMII head attachment rate; verified on-rates are clutch rates. | Wrong mechanism if reused. | Keep GAP. |
| `r0_head` | No source for the current 0.2 µm head-arm reach. | Geometric proxy, not a measured lever-arm length. | Keep GAP and do not conflate head reach with bare-zone length. |
| backbone `L_p` | Adamovic et al. model an approximately 10 nm scallop-muscle S2 fragment and compare its stiffness with a 130–170 nm range reported across coiled coils; Kaufmann–Schwarz assume 130 nm for individual NM2 rods. | Neither is a measured mature multi-tail human NMIIA/B minifilament rigidity. | Keep production bending off; these values remain labeled diagnostic proxies only. |
| cortical minifilament density | No absolute MCF7 density/copy-number record found. | Missing population contract. | PI/experiment/quantitative proteomics required; do not back-calculate density from a tension target. |

## Proposed PI decision order

1. Choose the first-baseline paralog contract: NM2A, NM2B, or an explicitly sourced heterotypic mixture.
2. Ratify the matching topology (`N_side`, total length, biochemical-state-specific bare zone).
3. Ratify the matching measured kinetics (`duty`, `k_off0`, `d_step`) only where the assay applies.
4. Commission/source the still-blocking mechanochemical values (`F_stall_head`, `v0`, Hill curvature, `k_xb`,
   head `k_on`, head reach, mature-backbone rigidity) and the MCF7 cortical minifilament density.
5. Create/link atomic SourceEvidence → KnowledgeClaim → Parameter rows in Notion, run source/parameter audits,
   then refresh TAG/Obsidian. Only after that should non-null production values enter the runtime config.

## Citation integrity — all 7 candidate DOIs verified (2026-07-21)

Every source this packet leans on was resolved against the Crossref REST API before being proposed for
registration; none is a hallucinated edge. **KB status** = whether the exact paper is already an audited
`SourceEvidence` row (queried against `outputs/tag_kb/kb.duckdb`).

| Candidate | DOI | Crossref title match | KB status |
|---|---|---|---|
| Melli 2018, *eLife* | `10.7554/eLife.32871` | OK — "Bipolar filaments of human nonmuscle myosin 2-A and 2-B have distinct motility…" | **absent** (Melli hits are homonyms) |
| Liu 2017, *PNAS* | `10.1073/pnas.1702375114` | OK — resolves, PNAS 2017 | **absent** |
| Nagy 2013, *JBC* | `10.1074/jbc.M112.424671` | OK — "Kinetic Characterization of Nonmuscle Myosin IIB at the Single Molecule Level" | **absent** (0 hits) |
| Kovács 2003, *JBC* | `10.1074/jbc.M305453200` | OK — "Functional Divergence of Human Cytoplasmic Myosin II" | **absent** (Kovacs hit is a homonym) |
| Dey 2017, *MBoC* | `10.1091/mbc.E16-07-0524` | OK — "Differential role of nonmuscle myosin II isoforms during blebbing of MCF-7…" | **absent** |
| Adamovic 2008, *Biophys J* | `10.1529/biophysj.107.122028` | OK — "The Elastic Properties of the Structurally Characterized Myosin II S2 Subdomain" | **absent** (0 hits) |
| Kaufmann & Schwarz 2020, *PLoS CB* | `10.1371/journal.pcbi.1007801` | OK — "Electrostatic and bending energies predict staggering and splaying in nonmuscle myosin II…" | **absent** (0 hits) |

Related sources already ingested but **unclassified** (Source Type + Claims blank, `PI to classify`), so they
cannot yet anchor a Parameter: `Billington2013` (`10.1074/jbc.m113.499848`), `Veigel2003`,
`StamAlbertsGardelMunro2015`, `Freedman2017_BiophysicalJournal`. The whole evidentiary basis of the NMII
magnitudes is therefore **either unregistered or unclassified** — no runtime value can close until these are
staged into Notion and audited.

## The single unblocking decision — pick the first-baseline paralog

Nearly every magnitude cascades from ONE PI choice: the first-baseline NMII paralog contract. Dey 2017 shows
MCF7 expresses NMIIA, IIB and IIC with **NMIIB** holding the longer cortical high-FRET dwell — so a NM2B or
NM2B-weighted heterotypic baseline is the MCF7-defensible starting point, but this is a PI call, not set here.

| Runtime field | NM2A baseline | NM2B baseline | Direct source | If heterotypic |
|---|---|---|---|---|
| `L_bb` | 0.301 µm (301±24 nm) | 0.323 µm (323±24 nm) | Melli/Liu (Billington meas.) | per-population mix |
| `N_side` | ~15/side (≈30 total) | ~14 dimers/side (≈28 heads) | Melli / Nagy | ≈30 total |
| bare zone | 167.9±29.7 nm unphos / 149.4±18.8 nm phos | GAP | Liu (NM2A only) | choose per paralog |
| `duty` | 0.05 | 0.20–0.25 | Melli (NM2A) / Nagy (NM2B) | population-weighted |
| `k_off0` | GAP (no NM2A assay) | 0.4 s⁻¹ (ADP release 0.35) | Nagy (NM2B) | NM2B rows only |
| `d_step` | GAP (reject 11 nm muscle) | ~6 nm | Nagy (NM2B) | NM2B rows only |

**Still GAP regardless of paralog** (needs new sourcing or experiment, cannot be closed by the paralog pick):
`F_stall_head`, `v0`, `kappa_hill` (Hill curvature), `k_xb` (master force knob), head `k_on`, `r0_head` (lever
reach), mature-minifilament backbone rigidity `L_p`, and **MCF7 cortical minifilament density**. Backbone
bending stays production-off until `L_p` is a measured mature-minifilament value (Adamovic/Kaufmann are S2 /
single-rod proxies only).

## Reconciliation against the runtime GAP registry (`ac/motor/params_i0b3.yaml`)

The YAML registry predates this source set and carries older claims; where they diverge, the newer direct human
NM2 assays above supersede, but selection still awaits PI:

- `duty`: YAML holds `0.1` (Kovács, provisional cross-check). Packet: paralog-specific 0.05 (NM2A) / 0.20–0.25
  (NM2B). The YAML `0.1` is neither paralog's measured value.
- `k_off0`: YAML holds provisional `0.35 s⁻¹` labeled "Stam-Hocky/Tam". Nagy resolves this as NM2B **ADP-release**
  0.35 vs **detachment** 0.4 s⁻¹ — the YAML value is the ADP-release rate mislabeled as the Bell detachment
  prefactor. PI must pick which rate feeds `k_off0`.
- `N_side`: YAML `claim_b = 28` (Billington/Niederman-Pollard) is consistent with Nagy's ~14 dimers/side; the
  AFINES `claim_a = 10` surrogate should be retired for the production count.
- `d_step`: YAML band `[5,10] nm` is consistent with Nagy NM2B ~6 nm; bind the exact value with the paralog.

## Candidate SourceEvidence staging (PI-gated — NOT auto-created)

Per CLAUDE.md, new `SourceEvidence` rows are PI-authored in Notion (SoT); this table is the staging manifest to
approve, then create in Notion and refresh TAG/Obsidian. All DOIs are Crossref-verified above.

| Proposed `citation_key` | DOI | Claims it would anchor |
|---|---|---|
| `Melli2018_eLife` | `10.7554/eLife.32871` | NM2A/B filament composition (~30 molecules), paralog motility, NM2A duty 0.05 |
| `Liu2017_PNAS` | `10.1073/pnas.1702375114` | NM2 filament + bare-zone geometry vs ATP/RLC-phos state |
| `Nagy2013_JBC` | `10.1074/jbc.M112.424671` | NMIIB duty 0.20–0.25, detachment 0.4 / ADP 0.35 s⁻¹, ~6 nm power stroke, ~14 dimers/side |
| `Kovacs2003_JBC` | `10.1074/jbc.M305453200` | NMIIA transient kinetics (currently only an unregistered claim source) |
| `Dey2017_MBoC` | `10.1091/mbc.E16-07-0524` | MCF7 NMIIA/B/C paralog-specific cortical behavior (relative, not absolute abundance) |
| `Adamovic2008_BiophysJ` | `10.1529/biophysj.107.122028` | S2-fragment elasticity — **proxy only** for backbone `L_p`, flag as diagnostic |
| `KaufmannSchwarz2020_PLoSCB` | `10.1371/journal.pcbi.1007801` | 130 nm single-rod persistence length assumption — **proxy only**, flag as diagnostic |

Also finish classifying the four unclassified auto-ingests (`Billington2013`, `Veigel2003`,
`StamAlbertsGardelMunro2015`, `Freedman2017_BiophysicalJournal`) so they can anchor Parameters.

## Sources checked

- Melli et al. 2018, *eLife* 7:e32871, DOI `10.7554/eLife.32871` — human NM2A/B filament composition and
  paralog-dependent motility.
- Liu et al. 2017, *PNAS* 114:E6516–E6525, DOI `10.1073/pnas.1702375114`, PMCID `PMC5559010` — human NM2
  filament and bare-zone geometry versus ATP/RLC phosphorylation.
- Nagy et al. 2013, *JBC* 288:709–722, DOI `10.1074/jbc.M112.424671` — human NMIIB duty ratio,
  detachment/ADP-release rates, and approximately 6 nm power stroke.
- Kovács et al. 2003, *JBC* 278:38132–38140, DOI `10.1074/jbc.M305453200` — human NMIIA transient kinetics;
  candidate source, not yet registered/audited in the Contract-Graph.
- Dey et al. 2017, *Molecular Biology of the Cell* 28:1034–1042, DOI `10.1091/mbc.E16-07-0524` — MCF7
  NMII paralog-specific cortical behavior, without an absolute abundance/density measurement.
- Adamovic et al. 2008, *Biophysical Journal* 94:3779–3789, DOI `10.1529/biophysj.107.122028` — atomistic
  elasticity of a scallop-muscle myosin-II S2 fragment; proxy only for the production primitive.
- Kaufmann & Schwarz 2020, *PLoS Computational Biology* 16:e1007801, DOI
  `10.1371/journal.pcbi.1007801` — NM2 rod-assembly model that assumes a 130 nm individual-rod persistence
  length; not a mature-minifilament rigidity measurement.
