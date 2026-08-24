---
id: 01_a-computational-model-for-early-cell-spreading-migration-and
paper_n: 1
title: "A computational model for early cell spreading, migration, and competing taxis"
authors: "J. Betorz et al. (Betorz, Bokil, Deshpande, Kulkarni, Araya, Venturini, Saez)"
year: "2023"
venue: "Journal of the Mechanics and Physics of Solids (JMPS), vol. 179, 105390"
doi: "10.1016/j.jmps.2023.105390"
paper_type: continuum-model
ffn_relevance: Medium
ffn_themes: [motor/myosin, cortex, FA/clutch, membrane, cell-spreading, ECM/collagen, numerics/methods, parameter-source]
entities: [actin, f-actin, g-actin, myosin-ii, gtpase, rac1, rhoa, cdc42, arp2-3, integrin, talin, vinculin, focal-adhesion, lamellipodium, cell-membrane, nucleus, keratocyte]
methods: [continuum-model, finite-element, reaction-diffusion, advection-diffusion, active-gel, crank-nicolson, supg-stabilization]
measurables: [retrograde-flow, migration-speed, membrane-tension, traction-stress, cell-spreading-radius, actin-density, polymerization-velocity, elastic-modulus]
keywords: [cell-motility, retrograde-flow, actin-polymerization, durotaxis, chemotaxis, wave-pinning, active-gel, membrane-tension, mesenchymal-migration, amoeboid-migration, molecular-clutch, stall-tension]
tags: ["#cell-motility", "#continuum-model", "#active-gel", "#durotaxis", "#chemotaxis", "#retrograde-flow", "#parameter-source", "#1d-finite-element", "#validation-context"]
has_transferable_params: true
---

# [1] A computational model for early cell spreading, migration, and competing taxis

**Tags:** #cell-motility #continuum-model #active-gel #durotaxis #chemotaxis #retrograde-flow #parameter-source #1d-finite-element #validation-context

| Field | Value |
|---|---|
| Authors | J. Betorz et al. (P. Sáez corresponding, UPC-BarcelonaTech) |
| Year / Venue | 2023 / J. Mech. Phys. Solids 179, 105390 |
| DOI / ID | 10.1016/j.jmps.2023.105390 |
| Type | continuum-model (1D finite-element active-gel + reaction–diffusion) |
| Pages | 18 |
| ffn_cellsim relevance | Medium — a LUMPED continuum competitor model, not a particle simulator, but a rich source of phenomenology, target numbers (retrograde flow, membrane tension, migration speed), and a validation/cross-check oracle for the emergent single-cell motility of ffn_cellsim |

## 1. Summary
The authors build a single, unified 1D continuum (finite-element) model of single-cell motility that couples (i) two-species advection–diffusion–reaction transport of F-actin / G-actin, (ii) transport of bound/unbound myosin motors, (iii) a wave-pinning Rho-GTPase reaction–diffusion polarization module, (iv) an active viscous "active-gel" momentum balance for the actomyosin retrograde flow with myosin contractility and ECM friction, and (v) a membrane-tension-limited actin protrusion velocity at the moving cell boundaries. With **one parameter set**, the model reproduces the three phases of symmetric cell spreading (P0/P1/P2), mesenchymal migration, and amoeboid (confined/bleb-like) migration, then is used to study chemotaxis vs durotaxis and their competition ("tug-of-war"). Main biological finding: chemotaxis dominates durotaxis in most regimes, but a strong stiffness gradient on a soft substrate can diminish (not reverse) chemotaxis; inhibiting GTPase control of front actin polymerization lets durotaxis prevail on soft matrices. A free platform is promised.

## 2. Problem & motivation
Cell motility underlies development, wound healing, and cancer invasion/metastasis. Prior models each captured only one phase (spreading), one mode (mesenchymal OR amoeboid), or one taxis (chemo- OR duro-) in isolation. No single mechanistic model reproduced ALL phases/modes with one parameter set, and the direct *competition* of chemical vs mechanical guidance cues — which coexist in vivo — had not been studied theoretically or experimentally. This paper unifies the mechanisms and resolves the competition question.

## 3. Methods / model
- **Model class:** 1D continuum on a moving domain Ω = [l_r(t), l_f(t)]; cell length L(t)=l_f−l_r; migration velocity V_cell = (l̇_r + l̇_f)/2. Boundary protrusion velocity l̇_{r,f} = v_p − v_F (polymerization minus retrograde flow). Lumped/coarse-grained — densities are continuous fields, NOT explicit filaments/motors.
- **F-actin / G-actin transport (Eqs. 1–2):** two advection–diffusion–reaction PDEs with polymerization rate k_p and depolymerization k_d; G-actin diffuses only; cell-frame velocity w = v_F − V_cell; zero-flux BCs. A reduced single-species F-actin form is mainly used (G-actin assumed spatially constant, ρ_GL(t)=1−∫ρ_F dx).
- **Myosin transport (Eqs. 3–5):** bound ρ_M / unbound ρ_m two-species; with fast-binding assumption (k_uM/k_bM ≪ 1) reduces to a single advection–diffusion equation with effective diffusion D = (k_uM/k_bM)D_M + D_m.
- **GTPase polarization (Eqs. 6–8):** Mori et al. (2008) wave-pinning: bistable reaction–diffusion of active ρ_R / inactive ρ_S forms, with D_R ≪ D_S, positive feedback k_on = k0 + γ ρ_R²/(K²+ρ_R²), k_off = 1; external stimulus f_S = k_S(x,t) ρ_S.
- **Actomyosin mechanics (Eqs. 9–11):** active viscous gel — overdamped momentum balance ∂_x σ_F = η_F v_F (friction with ECM); constitutive σ_F = μ_F ∂_x v_F + ζ ρ_F ρ_M (shear viscosity + myosin active contraction ζ). Neumann BC σ = τ(L) where membrane tension is a Hookean spring τ(L)=k(L−L_b), L_b=L_0+L_r (resting + membrane-reservoir buffer).
- **Protrusion velocity (Eq. 12):** v_p = v_p0 [1 − τ(L)/τ_stall]^γ, free velocity v_p0 = k_on δ ρ_R (δ = monomer size); decay exponent γ = 8 (keratocytes, Keren 2008).
- **Numerics:** FEM in space (element size h=L/N), implicit 2nd-order Crank–Nicolson in time, staggered solve of Eqs. (1),(5),(6),(9); SUPG stabilization for convection-dominated (Péclet Pe = h v_F / 2D > 1) regimes; no remeshing.
- **Taxis modeling:** chemotaxis = imposed external GTPase-activating signal k_S; durotaxis = spatially varying effective friction η(E) fit from clutch-model traction data (Elosegui-Artola 2016; Saez & Venturini 2023) over substrate stiffness E ∈ 0.1–100 kPa, sample lengths 100/200/2000 μm.

## 4. Key results (quantitative)
- **Spreading kinetics:** fast phase P1 lasts ~3 min; cell radius increases ~3-fold by end of P1; steady state (P2) reached at ~10 min (Fig. 2c, p.7). Radius follows a power law early — matches Cuvelier/Dubin-Thaler/Gauthier data.
- **Membrane tension (spreading):** rises to ~0.05 nN/μm, which suppresses front polymerization (Fig. 3c, p.7); in agreement with Shi et al. 2018.
- **Retrograde flow (spreading):** max actin flow v_F = 0.045 μm/s at the cell sides, vanishes at center (Fig. 3f, p.8); matches Wilson 2010 / Raz-Ben Aroush 2017 / Barnhart 2015.
- **Actomyosin stress (spreading):** symmetric, peak ~30 Pa at cell center (Fig. 3g, p.8) — imposes traction on the nuclear region.
- **Steady-state density accumulation (spreading):** actin ~1.08, myosin ~1.5 at center (normalized) (Fig. 3e, p.8).
- **Mesenchymal migration:** retrograde flow v_F ≈ 0.1 μm/s at front and rear (p.9).
- **Inhibition test:** inhibiting Rac1→actin polymerization gives V_cell ≈ 35 nm/s; inhibiting RhoA→myosin gives V_cell ≈ 65 nm/s (Fig. 5, p.10) → front actin polymerization has the STRONGER effect on migration speed than myosin contractility.
- **Amoeboid migration:** steady migration velocity 0.13 μm/s (matches Bergert 2015, Stroka 2014); actin/myosin densities ~2-fold/6-fold higher at rear vs front; self-polarizes via imposed-stress BC, no GTPase signal needed (Fig. 6, p.10–11).
- **Durotaxis:** cells migrate up positive stiffness gradients with V_cell ≈ 25–45 nm/s; maximum 45 nm/s for cells starting at the softest location of the steepest-gradient (L=100 μm) sample (Fig. 8, p.11–12). Cells on stiff regions do NOT durotax.
- **Competition:** aligned chemo+duro cues are roughly ADDITIVE (total ≈ sum of single velocities); opposed cues → chemotaxis dominates in most cases but is diminished (not reversed) by strong durotaxis on soft substrate; uncoupling GTPase control of front polymerization lets durotaxis cancel chemotaxis on soft matrices (Figs. 9–10, p.12–15).

## 5. Parameters & constants of interest

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Protrusion stall-decay exponent γ | 8 (dimensionless), keratocytes | Eq. 12, p.7 (Keren 2008) |
| Free polymerization velocity form | v_p0 = k_on δ ρ_R | §2.5, p.6 |
| GTPase basal conversion rate k0 | 1 s⁻¹ | Eq. 8, p.6 |
| GTPase max feedback rate γ | 1 | Eq. 8, p.6 |
| GTPase saturation K | 1 | Eq. 8, p.6 |
| GTPase inactivation rate k_off | 1 | Eq. 8, p.6 |
| Max actin retrograde flow v_F | 0.045 μm/s (spreading); ~0.1 μm/s (mesenchymal) | Fig. 3f p.8; p.9 |
| Membrane tension (spreading, P1 end) | ~0.05 nN/μm | Fig. 3c, p.7 |
| Actomyosin stress peak (spreading) | ~30 Pa at center | Fig. 3g, p.8 |
| Mesenchymal migration speed V_cell | ~35 nm/s (Rac1-inhib) … ~65 nm/s (RhoA-inhib) | Fig. 5, p.10 |
| Amoeboid migration speed | 0.13 μm/s | p.10 |
| Durotaxis migration speed | ~25–45 nm/s | Fig. 8, p.12 |
| Spreading P1 duration / radius gain | ~3 min / ~3-fold radius | Fig. 2c, p.7 |
| Steady-state reached | ~10 min | Fig. 2c, p.7 |
| Substrate stiffness range (durotaxis) | E ∈ 0.1–100 kPa | §3.3.2, p.11 |
| Sample lengths | 100 / 200 / 2000 μm | §3.3.2, p.11 |
| Diffusivity ordering (GTPase) | D_R ≪ D_S (active slow, inactive fast) | §2.3, p.6 |
| Myosin fast-binding ratio | k_uM/k_bM ≪ 1 | §2.2, p.5 |

Note: the full numerical parameter table (Appendix B, Table B.1: D_F, D_G, D_M, k_p, k_d, μ_F, ζ, η_F, k, L_0, L_r, δ, etc.) lives in the Supplementary Material, not in this main-text PDF; the values above are those stated inline in the body.

## 6. Relevance to ffn_cellsim  (MOST IMPORTANT)
This is a **continuum (lumped) competitor model**, the exact opposite of ffn_cellsim's architectural principle: here actin, myosin, GTPases, and adhesion are continuous density fields with effective diffusion/friction/contractility coefficients, solved by 1D FEM — none of them are explicit particles/bonds. So it is **NOT a mechanism source for the runtime** and **NOT a numerics source** (FEM/Crank–Nicolson/SUPG are irrelevant to a HOOMD BAOAB particle sim). Its value to ffn_cellsim is threefold:

- **(a) Validation oracle / emergent-behavior cross-check (primary).** ffn_cellsim's whole bet is that motility *emerges* from the fine-grained AFINES filaments + Stam-Hocky myosin + clutch adhesions. This paper gives the **macroscopic targets** that emergence must reproduce: retrograde flow ~0.045–0.1 μm/s, peak actomyosin stress ~30 Pa, membrane tension ~0.05 nN/μm, single-cell migration speeds ~25–130 nm/s, ~3-fold spreading-radius gain in ~3 min, steady state ~10 min, and the qualitative ordering "front polymerization > myosin contractility for migration speed." These are clean numbers to band a ffn_cellsim single-cell run against. This directly supports the H.5/H.7 single-cell motility milestones and the PI-experimental traction overlay track.
- **(b) Parameter source / sanity anchor (secondary).** The stall-velocity exponent γ=8, the v_p = v_p0[1−τ/τ_stall]^γ force–velocity protrusion law (note: this is a *membrane-tension* stall law, complementary to the Hill 1938 motor force–velocity used at runtime), stiffness range 0.1–100 kPa, and the η(E) friction-vs-rigidity curve from Elosegui-Artola 2016 are usable as oracle anchors / Magic-Number-Block literature citations.
- **(c) Phenomenology / spheroid-scale context (tertiary).** Clear, well-cited articulation of the spreading P0/P1/P2 phases, mesenchymal vs amoeboid modes, GTPase front/back polarization, and the chemo/duro competition. Useful narrative scaffolding for the EXTEND membrane track (membrane-tension-limited protrusion) and for framing multicellular/taxis behaviors the simulator might later reach.

Maps to themes: **motor/myosin** (contractility ζ ρ_F ρ_M), **cortex/actomyosin** (retrograde flow), **FA/clutch** (η(E) friction abstracts the molecular clutch — explicitly derived from clutch models in their own Saez & Venturini 2023), **membrane** (tension-limited protrusion, reservoir buffer L_r — relevant to EXTEND H.8), **cell-spreading**, **ECM/collagen** (substrate stiffness/durotaxis). It is partly off-target as a *mechanism* reference (it lumps everything ffn_cellsim resolves), but on-target as a *validation oracle and parameter anchor*. Hence Medium, not High.

## 7. Limitations & caveats
- **1D and lumped.** Real geometry, 2D/3D shape, lamellipodial branching network, individual filament orientation, and discrete clutch stochasticity are all collapsed into scalar fields and effective coefficients — precisely the fidelity ffn_cellsim is built to recover, so this model cannot validate *mechanism*, only *emergent macro-observables*.
- **Two-species actin / single-myosin reductions** explicitly acknowledged as simplifications of true actin turnover (no capping, no explicit Arp2/3 branching, no ATP/ADP cycle beyond effective k_p/k_d).
- **Adhesion is a friction coefficient η(E)**, not Bell-Evans catch/slip bonds; clutch behavior is imported as a fitted curve, not simulated.
- **Amoeboid mode is a boundary-condition trick** (imposed stress BC), not a bleb/water-flux mechanism.
- **No explicit data fit in this paper** ("No data was used"); validation is against literature numbers — so transferable values are themselves literature-derived, not new measurements.
- Full quantitative parameter table is offloaded to Supplementary Material (not in this PDF).

## 8. Key figures / tables
- **Fig. 3 (p.8):** spreading steady-state — actin/myosin densities, retrograde flow (v_F=0.045 μm/s), membrane tension (0.05 nN/μm), actomyosin stress (30 Pa). The richest single source of ffn_cellsim validation targets.
- **Fig. 5 (p.10):** Rac1-vs-RhoA inhibition → 35 vs 65 nm/s; establishes front-polymerization > contractility for speed.
- **Fig. 6 (p.10–11):** amoeboid migration at 0.13 μm/s with rear-accumulated actomyosin.
- **Fig. 8 / Fig. 9 (p.12–13):** durotaxis (25–45 nm/s) and chemo-vs-duro competition maps across stiffness gradients and sample lengths.
- **Appendix B, Table B.1 (Supplementary):** complete parameter set (not in main PDF).

## 9. Notable quotes / citable claims
- "The maximum actin flow velocity is v_F = 0.045 μm/s at the cell sides … in agreement with previous data." (p.8)
- "the stress of the actomyosin network is also symmetric with a maximum at the cell center of 30 Pa … imposes traction on the nuclear region." (p.8)
- "actin polymerization at the cell front has a stronger effect on cell migration than myosin activity." (p.10)
- "the cell reaches a steady state migration velocity of 0.13 μm/s, in agreement with previous data" (amoeboid). (p.10)
- "chemotaxis directs cell migration in most situations … inhibiting the effect of GTPases … makes durotaxis prevalent in soft matrices." (Conclusions, p.14)
- "v_p = v_p0 [1 − τ(L)/τ_stall]^γ … γ … equal to 8 in keratocytes." (Eq. 12, p.7)
