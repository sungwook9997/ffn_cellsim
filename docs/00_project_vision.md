# 00 — Project Vision & Framing

## The Core Idea (One Sentence)
Build an **independent first-principles 3D mechanobiological simulation** of MCF7 spheroid spreading on Col1-coated substrate, leveraging fluid-dynamics intuition from thin-film drying research, to generate new physical insight rather than to fit existing experimental data.

## Why This Matters
The PI's current research uses a phenomenological radial-symmetric model:

`A/A₀ = a + b/R + c/R²`

This fits experimental spreading curves well, but suffers from a fundamental epistemic limitation: **the model assumes radial-only forces**, while actual cellular spreading involves tangential traction, internal vortices, anisotropic active stress, and stochastic protrusions. Fitting the radial model to data cannot tell us whether the radial assumption is **valid** or just **incidentally correlated**.

This simulation addresses that gap by:
1. Building a full 3D anisotropic active hydrodynamic simulation from first principles (no fitting to PI's data)
2. Deriving (separately) a properly-reduced radial approximation through continuum mechanics, not by ad-hoc force projection
3. Comparing the two to quantify when/where the radial form is valid
4. Probing additional phenomena (Marangoni-like internal flow, coffee-ring analog density gradients, mechano-osmotic spreading-induced volume loss) inspired by the PI's own thin-film drying research

If results align with experiments → strong independent validation.
If results disagree → that **is** a finding (missing physics identified).

Either way, this is academic contribution.

## Scope
**Stage 1 (current)**: MCF7 cell line + Col1-coated rigid substrate. Single cell line, single substrate, parameter sweep over adhesion network composition (5 points spanning ULA-like ↔ pV4D4-like).

**Future stages**: documented in `10_dev_roadmap.md`, but explicitly out of scope for Stage 1.

## Why MCF7 + Col1 (Choice Rationale)
- **MCF7**: epithelial, strong E-cadherin, low invasion → spheroid behaves most like a coherent droplet → continuum/active-fluid framework is most valid here
- **Col1 coating**: most-studied baseline ECM, well-characterized mechanical properties, baseline state (no Lam supplementation) for clear interpretation
- **Together**: the regime where the simulation framework should work BEST — successful demonstration here is necessary before extending to invasive lines/conditions

## Framing Hierarchy
1. **Primary goal**: build a valid, publishable computational mechanobiology framework
2. **Secondary goal**: probe radial-vs-full discrepancies (methodological contribution)
3. **Tertiary goal**: test thin-film-drying-inspired hypotheses (Marangoni, coffee-ring, drying-concentration analogs) in cellular context
4. **Bonus**: alignment with PI's experimental data → validation

The simulation must be defensible **independently of the PI's experimental results**.

---

## Explicit Assumptions & Limitations (Stage 1)

This section is critical for academic honesty. Every assumption below is intentional; each has explicit justification and a defined "tier escalation" path if relaxation becomes necessary.

### Tier-0 (always assumed, won't be relaxed in this project)
- Newtonian-Cauchy continuum mechanics (relativistic and quantum effects ignored — trivially valid)
- Constant temperature 37 °C, constant pH, constant medium osmolality (matches culture conditions)
- No external mechanical perturbations (no shear flow, no pressure pulses)

### Tier-1 (Stage 1 simplifications, justified)
| Assumption | Justification | Escalation if needed |
|---|---|---|
| Cell-equivalent material points (~1 point ≈ 1 cell) | Bulk continuum + cell-level granularity sufficient for spheroid scale | Reduce cell count, refine to subcellular elements (Tier A model) |
| Constant cell mass (no proliferation/death) | 80 hr ≈ 1.5 MCF7 doubling cycles, secondary effect | Add proliferation as stochastic event (low cost) |
| All cells viable (no necrotic core) | Spreading thins spheroid → diffusion path shortens → necrosis suppressed | Activate Layer 6 chemistry (Stage 2) |
| Uniform cell properties (no heterogeneity) | First-principles baseline; heterogeneity is emergent perturbation | Add stiffness/size distributions in Stage 2 |
| Rigid substrate (no compliance) | Glass/TCPS effectively rigid (~GPa) compared to cells (~kPa) | Add substrate elasticity layer for PA gel comparison |
| Static ECM (no remodeling) | Col1 coating, no secreted matrix dynamics in Stage 1 | Add ECM remodeling for 3D collagen invasion studies |
| No chemical signaling (TGF-β, Wnt, etc.) | Mechanics-only baseline; adhesion remodeling captured via φ ODE | Couple φ ODE to chemical fields in Stage 3 |
| No nucleus mechanics | Too sub-cellular; nucleus stiffness absorbed into effective bulk modulus | Optional Stage 4 if mechano-genomics studied |
| Mechano-osmotic coupling: phenomenological (Tier 2) | Full poroelastic + ion channel kinetics has parameter explosion | Tier 3 (Biot poroelasticity + Pump-Leak) for water dynamics studies |
| External medium: Stokes drag only | Far-field hydrodynamic interactions weak at spheroid scale | Full Navier-Stokes external solver if microfluidic confinement |

### Tier-2 (would invalidate Stage 1 results if violated)
- Spheroid initial radius < 250 μm (otherwise necrotic core forms within 80 hr)
- Spreading is "slow" compared to cortical actin turnover (~30 s) → quasi-static cytoskeletal remodeling assumption
- No external active perturbation (drug, optical tweezers, etc.)
- Substrate fully wettable (effective spreading coefficient S > 0)

### What this simulation does NOT claim
- Does not predict **specific cell line** behavior beyond MCF7 in Stage 1
- Does not predict **drug response** (no pharmacokinetic layer)
- Does not predict **gene expression dynamics** (φ is mechanical proxy, not genetic state)
- Does not capture **single-cell migration** with high accuracy (cell-equivalent point granularity)

---

## Data Sources & Their Roles

### Primary literature (model construction)
- Active matter / continuum: Marchetti et al. Rev Mod Phys 2013; Banerjee & Marchetti
- Active wetting: Pérez-González et al. Nat Phys 2019
- Spreading dynamics: Douezan & Brochard-Wyart PNAS 2011; Beaune et al. PNAS 2014
- Mechano-osmotic coupling: Venkova et al. eLife 2022; Moeendarbary et al. Nat Mater 2013; Guo et al. PNAS 2017
- Cellular Marangoni: Pajic-Lijakovic & Milivojevic Eur Biophys J 2022; Phys Rev Fluids 2022
- Cell-cell adhesion energy: Maître et al. Science 2012
- Focal adhesion mechanics: Plotnikov et al. Cell 2012
- Active nematic order: Saw et al. Nature 2017; Duclos et al. Nat Phys 2017
- E-cad / Int-β1 crosstalk: Friedl & Alexander Cell 2011; Canel et al. J Cell Sci 2013
- Numerical method: Hu et al. SIGGRAPH 2018 (MLS-MPM); Taichi documentation

Full BibTeX in `references.bib`.

### PI's own experimental data
**Location**: `data/experimental/260313_{Bare,Pre,Lam4}.csv`
**Role**: Comparison overlay only (visualizing simulation results alongside experimental curves).
**Forbidden uses**: Parameter fitting, model selection, error minimization against this data.

Data summary (verified):
- Bare: 8 spheroids, 82 hr, 60 min interval, A/A₀ final 4.6–13.4
- Pre: 25 spheroids, 82 hr, 60 min interval, A/A₀ final 4.8–26.3
- Lam4: 26 spheroids, 83 hr, 60 min interval, A/A₀ final 8.0–33.1
- Initial effective radii: 87–419 μm

### Cho et al. 2020 (PI's lab paper, IF~5)
**Role**: Cross-check for adhesion remodeling timescales (φ ODE rates), NOT primary parameter source.
**Workflow**: Derive `k+`, `k-` from independent literature physical timescales (actin turnover, junction lifetime, etc.) → compare predicted Western blot trajectory to Cho's 12/24/48/72/96 hr timepoints → if matches, validation; if not, document as limitation.

---

## Authorship Note
This simulation framework, when published, should credit the thin-film drying intuition lineage from the PI's prior work even though no thin-film paper is cited directly in the model — the conceptual transfer (Marangoni, coffee-ring, drying concentration → cellular analogs) is the unique angle.
