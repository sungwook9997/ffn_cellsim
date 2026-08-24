---
id: J1_bulk_pressure_unjamming
topic: J1 — bulk pressure / solid stress in a compacting spheroid + jamming↔unjamming + pressure-driven escape/spreading
kind: research-synthesis (junction track)
year: 2026
ffn_relevance: High (bridge + parameter-source + validation-oracle)
ffn_themes: [spheroid-scale-context, junction, cortex, bulk-pressure, jamming-unjamming, validation-oracle, parameter-source, physiological-baseline]
entities: [solid-stress, osmotic-pressure, turgor-pressure, interstitial-fluid-pressure, cell-shape-index, SPV-vertex-model, gap-junction, dextran, p27Kip1, cadherin, cortical-tension]
methods: [osmotic-compression-dextran, stress-clamp, residual-stress-cutting, SPV-vertex-simulation, AFM, optical-diffraction-tomography]
measurables: [bulk-pressure-Pa-kPa, solid-stress-kPa, critical-shape-index, cell-volume-change, migration-velocity, proliferation-rate, doubling-time]
keywords: [osmotic compression, Delarue, Montel, Dolega, solid stress, Stylianopoulos, jamming, unjamming, shape index 3.81, Bi Manning, Han cell swelling, gap junction fluid flow, pressure release, reversible G1 arrest, cell escape, invasion]
tags: ["#bulk-pressure", "#jamming-unjamming", "#solid-stress", "#osmotic-compression", "#shape-index", "#validation-oracle", "#parameter-source", "#junction", "#physiological-baseline"]
has_transferable_params: true
membrane_vs_cortical: SUPRACELLULAR/bulk (aggregate-scale isotropic stress + tissue-fluidity), not single-cortex tension
---

# [J1] Bulk pressure / solid stress in a compacting spheroid, the jamming↔unjamming transition, and pressure-driven cell escape

**Tags:** #bulk-pressure #jamming-unjamming #solid-stress #osmotic-compression #shape-index #validation-oracle #parameter-source #junction #physiological-baseline

> ⭐ **The pressure ledger for a compacting multicellular aggregate.** Three numeric regimes
> matter to ffn_cellsim: (1) resting baseline (~tens of Pa hydrostatic / cortical turgor),
> (2) growth-induced bulk solid stress that *builds* as a spheroid compacts (~0.1–20 kPa,
> with proliferation-gating onset already at ~0.5 kPa), and (3) the jamming↔unjamming
> structural transition controlled by the dimensionless cell **shape index** (q* ≈ 3.81 in 2D
> / ~5.4 in 3D). Compression jams + arrests; **release re-fluidizes and re-licenses spreading**.
> These are aggregate-scale **acceptance oracles + physiological-baseline anchors**, NOT a
> runtime mechanism — but they set the operating point and the multicell validation gate.

---

## 0. Why this matters for ffn_cellsim

The PI "physiological operating point" rule (CLAUDE.md, 2026-06-04) demands every compartment
start at its real in-vivo value: a cell is **not** an unpressurised floppy bag. This note
collects the *bulk/aggregate* pressure numbers (the spheroid-scale analogue of the single-cell
~40 Pa resting turgor) and the structural order parameter (shape index) that decides whether a
confluent cohort behaves like a jammed solid (no escape) or an unjammed fluid (collective
invasion / single-cell escape). The junction track needs these because **cadherin junction
tension and bulk pressure together set the shape index**, and the shape index is the candidate
multicell validation observable — directly in the same family as the L2 hull/shape metrics and
the Roffay [52] Young-Laplace `ΔP = σ(1/R+1/R')` oracle.

---

## 1. Pressure ledger — magnitudes that build in a compacting spheroid

| Regime | Magnitude | What it is | Source |
|---|---|---|---|
| **Resting single-cell hydrostatic / cortical "turgor"** | ~**40–100 Pa** (project baseline ~40 Pa) | Cortex-contraction + ECM-confinement back-pressure of an animal cell at rest; the physiological setpoint a cell starts from | Mammalian-cell mechanics reviews (Current Biology 2024 turgor review, doi:10.1016/j.cub.2024.07.084; JCS 2020 doi:10.1242/jcs.240341) |
| **Total trans-membrane osmotic pressure** | ~**0.1–1 MPa** | Balanced across the membrane by turgor — NOT a net force; do not confuse with the net solid stress | osmotic-shock estimates (reviews above) |
| **Onset of proliferation gating (applied osmotic)** | from ~**500 Pa**, saturating ~**5 kPa** | Compressive stress at which spheroid growth/motility starts to slow | Dolega 2021 eLife (doi:10.7554/eLife.63258); Delarue 2014 (below) |
| **Stress clamp that drastically halts growth** | ~**5–10 kPa** | Applied isotropic osmotic stress that arrests proliferation, mainly in the core | Montel 2011 PRL (doi:10.1103/PhysRevLett.107.188102); Delarue 2014 BpJ (doi:10.1016/j.bpj.2014.08.031) |
| **Growth-induced *endogenous* solid stress, murine tumors** | ~**0.37–8.0 kPa** (≈2.8–60 mmHg) | Residual mechanical stress that a tumor builds in itself as it compacts | Stylianopoulos 2012 PNAS (doi:10.1073/pnas.1213353109) |
| **Growth-induced solid stress, human tumors** | ~**2.2–19.0 kPa** (≈16.5–142.5 mmHg) | Same, human resected tumors (cut-and-relax) | Stylianopoulos 2012 PNAS |

**Take:** a spheroid that compacts as it grows builds *internal* solid stress of order **0.1–20 kPa**,
and the cell-cycle machinery already feels it at **~0.5 kPa** — i.e. the bulk pressure scale that
matters for proliferation/escape is **3 orders of magnitude above** the resting single-cell turgor.
The resting baseline (~40 Pa) is the *starting* operating point; the compaction physics happens
in the **sub-kPa → tens-of-kPa** band.

---

## 2. How applied compression changes growth & spreading (the osmotic-compression assays)

### Delarue, Montel, Vignjevic, Prost, Joanny, Cappello — *Biophys. J.* 107(8):1821–1828 (2014), doi:10.1016/j.bpj.2014.08.031
"Compressive stress inhibits proliferation in tumor spheroids through a volume limitation."
- **Method:** 100 kDa Dextran in medium (excluded from the aggregate) ⇒ osmotic stress on the outer
  shell, mechanically transmitted inward as a quasi-isotropic compressive stress. Applied **5 kPa**
  (55 g/L) and **10 kPa** (80 g/L).
- **Minutes (volume):** cell-to-cell distance fell ~**20%** in the MCS center within 5 min at 10 kPa
  (p<0.002) — a fast, near-elastic volume reduction.
- **Hours:** reversible induction of cell-cycle inhibitor **p27Kip1** propagating center→periphery.
- **Days:** cells arrest at the **late-G1 restriction point** (pRb-Thr373); bulk growth rate *k*
  drops by **≥ ×2**. Conserved across 5 lines (HT29, CT26, BC52, FHI, AB6).
- **★ Reversibility:** "the compression-induced proliferation arrest … is **completely reversible**" —
  silencing p27Kip1 antagonizes the pressure effect. Spheroid: ~200 µm initial diam, ~4×10⁶ µm³,
  10-day assay.

### Montel et al. — *Phys. Rev. Lett.* 107:188102 (2011), doi:10.1103/PhysRevLett.107.188102
"Stress clamp experiments on multicellular tumor spheroids." A clamped **10 kPa** isotropic stress
**drastically reduces growth by inhibiting proliferation mainly in the spheroid core** (apoptosis +
proliferation suppression) — the foundational stress-clamp result Delarue 2014 dissects.

### Dolega, Monnier, Brunel et al. — *eLife* 10:e63258 (2021), doi:10.7554/eLife.63258
"ECM in multicellular aggregates acts as a pressure sensor controlling cell proliferation and motility."
- Global compression with **big Dextran at 5 kPa** ⇒ mean cell **migration velocity −50%**, and
  spheroid **doubling time 36±1 h → 68±4 h** (~×1.9). Effect visible from **Πd = 500 Pa**, saturates
  ~**5 kPa**.
- **Mechanism subtlety (important for fine-grained modeling):** the ECM is **100–1000× more
  compressible than cells**, so it acts as the deformable pressure sensor — *selective* cell-only
  compression (small Dextran) produced **almost no effect**; the long-timescale compressive stress is
  transmitted through the dehydrated matrix. ⇒ in a particle model the **ECM cross-link compliance**,
  not the cell bulk modulus, is the dominant pressure transducer.

**Synthesis (compression direction):** isotropic bulk pressure of order 0.5–10 kPa → cell-volume
↓ → late-G1 arrest (p27Kip1) + motility ↓ → **jamming / spreading suppression**, reversibly.

---

## 3. The jamming ↔ unjamming transition — the structural order parameter

The transition is set by the **dimensionless cell shape index**, not by density alone (this is the
"density-independent rigidity" insight).

### Bi, Lopez, Schwarz, Manning — *Nat. Phys.* 11:1074–1079 (2015), doi:10.1038/nphys3471
- Vertex-model energy `E = Σ [K_A (A−A0)² + K_P (P−P0)²]`; effective **target shape index**
  `p0 = P0/√A0`. **Rigidity transition at p0* ≈ 3.81** (2D, shape index `p = P/√A`).
- **p0 < 3.81 ⇒ JAMMED / solid** (finite shear modulus, finite energy barriers to T1 rearrangement,
  cells caged → no escape). **p0 > 3.81 ⇒ UNJAMMED / fluid** (barriers vanish, cells flow → escape).
- The balance is **cell–cell adhesion vs cortical tension**: more adhesion / less cortical tension
  → larger p0 → unjammed. This is the junction lever (cadherin) directly setting the bulk phase.

### Bi, Yang, Marchetti, Manning — *Phys. Rev. X* 6:021011 (2016), doi:10.1103/PhysRevX.6.021011
- Self-Propelled-Voronoi (SPV): three control parameters — **self-propulsion speed v0, shape index
  p0, rotational noise / persistence Dr (τ = 1/Dr)**.
- Measured (snapshot) shape index `q = ⟨p/√a⟩`: **q < 3.81 solid (jammed), q > 3.81 fluid (unjammed)**;
  precise dynamical glass line `q ≈ 3.813`. Self-diffusivity D_eff turns on across q ≈ 3.8 (Fig. 2).
- **★ q is a STATIC, snapshot-measurable order parameter** — you can read jamming state straight off a
  geometry without tracking dynamics. This is the candidate ffn_cellsim multicell observable.

### Park, Kim, Bi, … Fredberg — *Nat. Mater.* 14:1040–1048 (2015), doi:10.1038/nmat4357
- Experimental confirmation in **primary human bronchial epithelium**: measured median shape index
  crosses ~**3.81** at the jamming transition; asthmatic cells stay **unjammed** (q above threshold)
  longer → more migratory. Direct in-vitro validation that q* ≈ 3.81 is real, not just a model artifact.

### 3D caveat (don't mix up the numbers)
The ~3.81 threshold is the **2D** perimeter/√area convention. In **3D** the analogous index is
`SI = surface area / volume^{2/3}`, with a critical value ~**5.4** (Merkel & Manning 2018,
New J. Phys. 20:022002, doi:10.1088/1367-2630/aaaa13; and Han iScience 2021 below uses 5.4). A
3D spheroid model must compare to ~5.4, a 2D monolayer slice to ~3.81.

---

## 4. Pressure RELEASE → sudden spreading / escape (the unjamming mechanism)

This is the load-bearing answer to the J1 question. Two complementary mechanisms:

**(a) Direct re-fluidization (compression is reversible).** Because the compression-induced arrest is
*purely reversible* (Delarue 2014: "completely reversible"; p27Kip1 decays, G1 block lifts), removing
the bulk stress lets cell volume re-expand, lowers density, raises the shape index back above q*, and
the cohort **re-enters the unjammed/fluid regime** — proliferation and motility resume. Release is the
mirror image of §2 compression.

**(b) Stress-gradient-driven swelling that pre-positions the periphery to unjam (Han mechanism).**
- **Han, Pegoraro, Li et al. — *Nat. Phys.* 16:101–108 (2020), doi:10.1038/s41567-019-0680-8**
  "Cell swelling, softening and invasion in a 3D breast cancer model" (MCF-10A organoids):
  - A **gradient of intra-tumor compressive stress** (high in core, low at periphery) drives
    **supracellular fluid flow through gap junctions** from compressed core → relaxed periphery.
  - Peripheral / invasive-branch cells become **larger (swollen), softer, and more dynamic**: branch
    cells show **~5× larger cytoplasmic force fluctuations at 1 Hz** vs core; peripheral cells migrate
    faster. Core cells = stiffest, branch cells = softest.
  - **★ Stress-release experiment:** releasing the external stress produced a **volume INCREASE in the
    core and DECREASE in the periphery** — direct proof the volume/stiffness gradient is driven by the
    *pressure gradient*, i.e. the periphery is locally "pressure-released" relative to the core.
  - Blocking gap junctions (carbenoxolone 500 µM / connexin-mimetic peptides) abolishes the gradient
    and **delays the invasive transition**. Artificially re-stiffening peripheral cells (osmotic
    compression, daunorubicin, jasplakinolide) **reduces invasion** — confirming softening/swelling is
    causal, not a passive marker.
- **Han follow-up — *iScience* 24(11):103252 (2021), doi:10.1016/j.isci.2021.103252**
  "A novel jamming phase diagram links tumor invasion to non-equilibrium phase separation."
  - Measured **3D shape index**: late-stage spheroid **core SI ≈ 5.84 ± 0.32 (jammed)**, **periphery
    SI ≈ 6.6 ± 0.79 (unjammed)**, approaching/exceeding the critical **SI ≈ 5.4**.
  - Phase diagram axes: cell motility (effective temperature T_eff) × confinement pressure P_conf
    (collagen density), with solid / fluid / gas-like phases. **Collective invasion = the unjammed
    fluid-like periphery splitting off the jammed solid core** (non-equilibrium phase separation).

**Mechanistic chain for "release → sudden spreading":**
relieved bulk pressure → cell volume re-expansion (water influx via gap junctions / channels) →
density drop + cortical softening → shape index climbs above q* (3.81 / 5.4) → energy barriers to T1
rearrangement vanish → unjamming → collective fluid-like spreading + single-cell escape. The
suddenness is the **sharp (near-critical) shape-index threshold**: a small drop in pressure / rise in
adhesion-to-cortical-tension ratio flips the cohort across the rigidity transition discontinuously.

**Developmental corroboration (non-cancer):** Mongera et al. *Nature* 561:401–405 (2018),
doi:10.1038/s41586-018-0479-2 — vertebrate body-axis elongation runs on the *same* physics: an
N-cadherin-dependent **yield-stress gradient** keeps the posterior progenitor zone **unjammed (fluid,
spreading/extending)** and the anterior presomitic mesoderm **jammed (solid)**. Confirms cadherin
junctions ↔ yield stress ↔ jamming as a general tissue-fluidity control knob, not a cancer artifact.

---

## 5. Relevance to ffn_cellsim — how to use these numbers

- **Physiological-baseline anchor (PI rule).** A multicell run must initialize the aggregate at its
  resting operating point (single-cell ~40 Pa turgor; aggregate near-zero *net* solid stress at small
  size) and then let solid stress **build emergently** to the 0.1–20 kPa band as it compacts — never
  start from an unphysical relaxed bag. The proliferation/motility gate should engage near **~0.5 kPa**
  (Dolega/Delarue), saturating ~5 kPa.
- **Validation oracle (multicell tier).** The **shape index q** (2D q*≈3.81 / 3D SI≈5.4) is a static,
  snapshot-measurable order parameter readable from a GSD trajectory — a clean acceptance gate for
  "did the emergent cohort jam or unjam?" Pairs with Roffay [52] `ΔP = σ(1/R+1/R')` and the [11]
  DAH/DITH surface-tension oracle. Han SI core/periphery (5.84 / 6.6) gives a *gradient* target.
- **Junction lever is the mechanism.** Adhesion-vs-cortical-tension ratio (cadherin catch-bond KU-4.2
  + actomyosin cortex H.1) is exactly what sets p0 in the vertex/SPV picture — so the fine-grained
  runtime *should* emergently reproduce the jamming line if the junction + cortex units are right.
  Unjamming = raise adhesion / lower cortical tension / lower bulk pressure.
- **ECM compliance, not cell bulk modulus, transduces bulk pressure (Dolega).** When the ECM track is
  on, cross-link compliance is the dominant pressure sensor — model the matrix as 100–1000× more
  compressible than the cell, or the pressure-gating physics won't appear.
- **Release → escape is reversible re-fluidization.** Reproduce by lowering applied confinement and
  checking q crosses back above threshold + motility resumes (Delarue reversibility; Han stress-release
  volume flip).

## 6. Classification & caveats
- **Class:** aggregate/supracellular-scale **physiological-baseline anchor + acceptance oracle +
  parameter source**. Paper-models (vertex/SPV) are **oracles only — never the ffn_cellsim runtime
  mechanism** (CLAUDE.md hard rule). q is an observable to *measure*, not a force law to *impose*.
- **Caveats:** (1) 2D (3.81) vs 3D (5.4) shape-index conventions must not be mixed. (2) Applied
  osmotic Dextran stress ≈ but ≠ endogenous growth-induced solid stress (Dolega: route matters — ECM
  vs cell). (3) Absolute solid-stress magnitudes are line- and host-dependent (Stylianopoulos spread
  0.37–19 kPa). (4) q is a *necessary structural* signature; active motility (v0) and persistence (Dr)
  also move the transition — a snapshot q near threshold is not sufficient without dynamics. (5) Han
  gap-junction fluid-flow mechanism is breast-organoid-specific; transfers as a *mechanism class*, not
  exact connexin numbers.

## 7. Primary citations (verify DOI verdict before deliverable use)
- Delarue M, Montel F, Vignjevic D, Prost J, Joanny JF, Cappello G. *Biophys J* 107(8):1821-1828 (2014). doi:10.1016/j.bpj.2014.08.031
- Montel F, Delarue M, Elgeti J, et al. *Phys Rev Lett* 107:188102 (2011). doi:10.1103/PhysRevLett.107.188102
- Dolega ME, Monnier S, Brunel B, et al. *eLife* 10:e63258 (2021). doi:10.7554/eLife.63258
- Stylianopoulos T, Martin JD, Chauhan VP, et al. *PNAS* 109(38):15101-15108 (2012). doi:10.1073/pnas.1213353109
- Bi D, Lopez JH, Schwarz JM, Manning ML. *Nat Phys* 11:1074-1079 (2015). doi:10.1038/nphys3471
- Bi D, Yang X, Marchetti MC, Manning ML. *Phys Rev X* 6:021011 (2016). doi:10.1103/PhysRevX.6.021011
- Park JA, Kim JH, Bi D, et al. *Nat Mater* 14:1040-1048 (2015). doi:10.1038/nmat4357
- Merkel M, Manning ML. *New J Phys* 20:022002 (2018). doi:10.1088/1367-2630/aaaa13
- Han YL, Pegoraro AF, Li H, et al. *Nat Phys* 16:101-108 (2020). doi:10.1038/s41567-019-0680-8
- Han YL, Ronceray P, Xu G, et al. *iScience* 24(11):103252 (2021). doi:10.1016/j.isci.2021.103252
- Mongera A, Rowghanian P, Gustafson HJ, et al. *Nature* 561:401-405 (2018). doi:10.1038/s41586-018-0479-2
