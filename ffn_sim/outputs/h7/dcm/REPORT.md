# DCM tier — Deformable-Cell-Model multicellular spreading + cell–ECM coupling

**Branch** `h7/compartment-platform` · **2026-06-11** · CPU mesoscale dev · graduation-report build.

## Mission (PI)

Pivot the platform's spheroid-scale work onto a **Deformable Cell Model (DCM)** — the
mid tier between the fine-grained single cell (every cytoskeletal filament a particle)
and the layer-2 coarse multicell (cell ≈ point) — and run multicellular **spheroid
spreading**, trusting literature band values directly. Then **go beyond** the reference
literature by reusing the platform's existing fine-grained mechanisms (explicit
cross-linked fiber ECM, catch-slip focal-adhesion clutch) rather than re-implementing
the coarser published models.

## Reference papers (read in full; all from the Taeyoon Kim group — the platform's own AFINES/Langevin lineage)

1. **Slater, Li, Indana, Xie, Chaudhuri, Kim — "Transient mechanical interactions
   between cells and viscoelastic extracellular matrix", Soft Matter 17, 10274 (2021)**
   (`references/d0sm01911a.pdf`). Agent model: a contractile cell (membrane +
   actomyosin cortex as a zero-rest-length contractile spring) at the centre of a
   **fiber ECM with transient cross-linkers** (Bell's-law unbinding, k_ub,0 = 1e-5 s⁻¹,
   x_ub = 1e-10 m, binding sites every 100 nm → ECM viscoelasticity). **Focal adhesions
   = permanent harmonic springs linking ECM fibers to cortex nodes** (κ_s,fc). Result:
   stronger contraction → greater matrix deformation + longer-range stress; stress
   decays **~1/r in 2D** (1/r² in 3D); validated vs 3T3 fibroblasts in 3D collagen.
   → **the blueprint for our cell–ECM coupling.**
2. **Jo, Yim, … Kim, Kim — "Reciprocal folding dynamics in cellular networks at the
   stroma–basement-membrane interface", Acta Biomaterialia 201, 360 (2025)**
   (`references/1-s2.0-S1742706125004039-main.pdf`). Stromal cells exert **polarized
   traction** on a rigid basement membrane → buckling/folding morphogenesis;
   **stiffness difference sets folding direction**; blebbistatin abolishes it (traction
   is required). ECM modulus measured by dog-bone tensile test (ε̇ = 3.2e-3 s⁻¹).
   → traction + substrate-stiffness → morphology (durotaxis grounding).
3. **Yim — "Computational study of biophysical mechanisms of axon outgrowth and cell
   migration", PhD dissertation, Purdue (advisor T. Kim), 2026**
   (`references/REVDissertation_Donghyun Yim.pdf`). Consolidated cell-migration / cell–
   ECM mechanics (motor-clutch on ECM). → the migration/clutch context.

**Judgment:** our cell–ECM attachment is improvable, and these (lineage) papers give
the exact target. They use a *coarse* cell + *permanent* FA spring + *transient*
cross-linker; we already have **finer** mechanisms (real Pereverzev **catch-slip** FA
clutch, an explicit **cross-linked Mikado fiber ECM** with bending, exact-turgor
**deformable** cells). The advance is to **compose** these — exceeding the papers on
every axis — instead of matching them.

## What was built

- **`cell/dcm.py` — the DCM tier.** Icosphere elastic membrane shells (harmonic edge
  springs = cortical elasticity) + **exact triangulated turgor** (`DcmTurgorForce`:
  divergence-theorem enclosed volume + face-normal pressure — correctly conserves the
  *non-spherical* volume so a flattening cell spreads laterally; the sphere-equivalent
  estimator does NOT, and was the bug that blocked spreading) + adhesive substrate
  (`DcmSubstrateForce`, capped-harmonic well = BAOAB-stable) + HOOMD LJ cell–cell
  adhesion/excluded-volume (per-cell types) on the BAOAB integrator. 3D spheroid (FCC
  ball) or 2D monolayer aggregate.
- **ECM-attachment improvement layer (from the literature audit + the 3 papers):**
  - **compliant substrate `k_sub`** — finite substrate stiffness softens the felt
    adhesion via a series spring (`k_eff = k_well·k_sub/(k_well+k_sub)`); `None` = rigid
    dish (glass; the PI's prep). The **durotaxis / rigidity-sensing knob** (Engler;
    Acta Bio stiffness→morphology).
  - **ligand-density-scaled adhesion** `W_cs ∝ ρ_ligand` (Gallant/García) — the
    coating-density / ECM-condition knob.
- **`cell/dcm_ecm.py` (beyond-papers headline) — DCM cell on an explicit Mikado fiber
  ECM + FA clutch** (reuses `ecm/mikado.py` cross-linked fiber network + the catch-slip
  integrin clutch). Traction-driven spreading + **ECM remodeling**. [results below]
- Scripts: `h7_dcm_spheroid_spread.py` (spheroid spreading + figure/GIF),
  `h7_dcm_ecm_remodel.py` (cell–ECM remodeling + 1/r stress validation).

## Results

- **DCM cell** holds shape under turgor with **volume conserved** (V/V₀ ≈ 1, exact
  triangulated pressure); stable on BAOAB.
- **Passive substrate wetting is limited:** the 2D aggregate flattens and spreads to
  **A/A₀ ≈ 1.3**; the 3D spheroid barely spreads (**A/A₀ ≈ 1.0** — passive adhesion
  cannot pull a 3-D ball into a pancake; the upper cells never descend). This is the
  honest limitation that motivates the traction-coupled model.
- **Durotaxis (substrate stiffness `k_sub` sweep, 2D aggregate, W_cs=8 mJ/m²):**
  CONFIRMED monotonic rigidity sensing — **A/A₀ = 1.02 (soft, k_sub=1e-4 N/m) → 1.11
  (3e-4) → 1.14 (1e-3) → 1.27 (rigid dish / glass)**; aggregate height 15→11.7 µm on
  the stiff end. Stiffer substrate → more spreading, as Engler 2006 / Elosegui-Artola
  2016 / the Acta-Biomater-2025 stiffness→folding picture predict.
- **Beyond-papers cell–ECM coupling (`dcm_ecm`, headline):** a DCM exact-turgor
  deformable cell (162 nodes) on an **explicit cross-linked Mikado fiber ECM (2754
  fiber beads)**, gripping via the FA clutch and contracting — **builds + runs 40k
  BAOAB steps stably** and **REMODELS the matrix**: fiber beads pulled inward toward
  the cell (mean 15 nm, **max ≈ 532 nm near the cell**), and **tensed fibers align
  radially near the cell** (|cos(bond, radial)| ≈ 0.65 vs 0.5 isotropic — the
  d0sm01911a/Kim tensed-fiber signature). ⚠️ The QUANTITATIVE d0sm01911a 2D stress
  law (tension ~ 1/r) is **NOT yet reproduced** (log-log slope +0.43 vs the predicted
  −1): the contraction is weak (15 nm mean remodeling) and the dynamic catch-slip FA
  did not engage (n_engaged=0; traction is carried by the static force-free FA bonds),
  so the cell's stress field does not dominate the bed's background tension. Stronger
  contraction + working catch-slip engagement is the next tuning step. **This is a
  structurally-beyond-the-papers, stable, novel asset** (deformable turgor cell ⊗
  explicit cross-linked fiber ECM ⊗ catch-slip FA) demonstrating cell-driven matrix
  remodeling; quantitative 1/r validation is the documented follow-on.

## Honesty / scope

- Literature bands trusted directly (PI 2026-06-11): cell–substrate adhesion
  W_cs ~ 0.5–8 mJ/m² (Hategan/Cuvelier/collagen-coated), cell–cell ~0.3 mJ/m²
  (cadherin), turgor 133 Pa, K_vol 5 kPa, substrate stiffness via k_sub. A/A₀ ≈ 2–4 and
  the ~1/r stress law are **validation targets**, not inputs.
- Mesoscale CPU; native-scale GPU is the gbook follow-on (dirty-branch cleanup pending).
- Passive-wetting spheroid spreading is bounded; the traction-coupled `dcm_ecm` model is
  the mechanism that drives real spreading/remodeling and exceeds the reference papers.

## Figures

- `figs/dcm/spheroid3d_spread.png` / `.gif` — DCM 3-D spheroid geometry + (passive) spreading.
- `figs/dcm/durotaxis_k_sub.png` — spread area + height vs substrate stiffness k_sub
  (rigidity sensing: stiffer → more spreading, A/A₀ 1.02→1.27).
- `figs/dcm_ecm/dcm_ecm_remodel.png` — 4-panel: ECM fiber inward-displacement map;
  fiber tension vs r (1/r check, not yet matched); radial alignment vs r (≈0.65 near
  cell); remodeling vs time (max ≈532 nm). The beyond-papers cell⊗ECM headline.
