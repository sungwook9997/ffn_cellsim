# TAG-validation digest (all 77 papers, sorted by score)


## p74 [4.0] Taeyoon Kim · research — 1Living Matter Department, AMOLF, Amsterdam, The Netherlands. 2Institu
- headline: TAG-validation lens — p74  (Mulla et al. 2022, "Weak catch bonds make strong networks", Nat Mater 21:1019)

## p61 [4] Taeyoon Kim · research — Disordered actomyosin networks are sufﬁcient to
- headline: Kim-lab disordered-actomyosin paper adds cooperativity (percolation ρc≈0.56 µm⁻²) + telescopic-contraction scaling laws as new FF cross-check gates; no existing gate covered, none tightened.
  - Telescopic law: boundary contraction velocity rises with activation radius (~0.025→0.07 µm/s, S/M/L radii 9.0/23.1/37.1 µm) while εmax and dε/dt stay ~constant (Fig 4c–e)
  - Cooperative contraction: εmax vs myosin density is a sharp percolation-like transition, εmax<0.5 below vs >1 (up to ~2.5) above critical density ρc≈0.56 µm⁻²; strain rate LINEAR in ρ (Fig 3)
  - Hill coefficient nH≈11 for εmax vs myosin density (FIT param, overlay only — not a portable band; Fig 3c)
  - Bending-invariance (SIMULATION): telescopic contraction retained as F-actin persistence length rises 20 µm→2.6 mm — filament buckling not required (Fig 6)
  - Filamin crosslinking suppresses contraction: strain ~2.4→0.3 and rate ~0.015→0.002 s⁻¹ at RC=[FLNa]/[actin]=0.035 (Fig 5h,i)
  - Active-stress buildup law σₐ(r,t)=σ₀(1−e^(−t/τₐ)); ε(t) three-phase kinetics (lag 0–100s / linear 100–300s / plateau >300s); ∇·v min ~−7×10⁻³ s⁻¹ (Methods, Fig 2c/3a)
  - FF network geometry: F-actin 280nm segments/7nm dia, myosin 42nm backbone segments, 7nm binding-site spacing, lp=20µm, ⟨l⟩≈7.1–7.5µm (corroborates FF KB-3.18 length dist)
  · n_new_gates=4
  · n_tighten=0
  · n_params=9
  · paper_kind=experimental

## p57 [4] Makito Miyazaki · experimental — Biphasic Effect of Proﬁlin Impacts the Formin mDia1
- headline: Single-molecule mDia1 formin force-velocity: an FF-only gap — closed-form Bell-gated elongation law (Eq.1, d=5.5 nm) as candidate analytic oracle plus in-vitro kon/koff/p_c⁰ constants; 3 new FF gate candidates, no DCM relevance, nothing superseded.
  - Force-dependent formin elongation law v = kon·C_A·[1−p_c⁰·e^(−fd/kT)] − koff·p_c⁰·e^(−fd/kT), d=5.5 nm (Eq.1, oracle candidate)
  - mDia1 kon = 21 µM⁻¹s⁻¹ (ATP) / 3.6 (ADP); koff = 1.3 (ATP) / 3.9 s⁻¹ (ADP), 26 °C
  - p_c⁰ zero-force closed-state prob. = 0.29 (ATP) / 0.61 (ADP)
  - Tension accelerates ADP elongation 3.4 → 10.9 subs/s (~3.2×), 0.13 → ~4 pN, 3 µM ADP-actin
  - Open-state probability p_open → ~1 above ~3 pN (elongation saturates by 3 pN)
  - Profilin force-switch: net −1.8 (depoly) → +40.2 subs/s, 3 µM ADP-actin + 3 µM profilin, low→high tension
  - Working distance d = 5.5 nm (FH2 dimer-half, single Bell force-sensing length); dumbbells break > 6 pN
  · n_new_gates=3
  · n_tighten=0
  · n_params=8
  · paper_kind=experimental

## p53 [4] Taeyoon Kim · research — Rapid assembly of a polar network architecture
- headline: Kim-lab agent-based actomyosin paper: rich FF-native formin/actin/crosslinker kinetics (v_elong 1.1-1.3 um/s, k_off 0.11/s, ~6um filaments, ACP Bell constants) + a filament-length analytic-oracle candidate; instantiates existing Bell-Evans oracle and flags a Ferrer-2008 ACP off-rate conflict (0.115 vs KB 0.066 s^-1).
  - Formin (CYK-1) barbed-end elongation speed 1.1-1.3 um/s (~450 monomers/s), C. elegans cortex
  - Formin turnover k_off=0.11 s^-1 (smPReSS); actin monomer turnover 0.05-0.15 s^-1
  - Mean formin-elongated F-actin length ~6 um via L=V/(1/tau_actin+1/tau_formin)
  - ACP (alpha-actinin) Bell off-rate k0_u=0.115 s^-1, x_u=1.04e-10 m (Ferrer 2008) — CONFLICTS with KB 0.066 s^-1
  - Actin bending stiffness kappa_b=2.64e-19 N.m (ell_p~15-17 um, Isambert 1995); extensional 1.69e-2 N/m
  - Myosin minifilament stoichiometry N_a=8 arms x N_h=4 heads = 32 heads (parallel-cluster-model)
  - Barbed-end half-life t1/2=ln2/(k_on[F])~3 s at k_on=0.029 nM^-1s^-1, [F]=10 nM
  · n_new_gates=3
  · n_tighten=0
  · n_params=9
  · paper_kind=experimental

## p46 [4] Taeyoon Kim · research — The nature of cell division forces in epithelial
- headline: First division/cytokinesis data for the KB (currently no such gate): a Taeyoon-Kim same-lineage DCM shows ~5% volume loss drives ~60% axial elongation — a strong DCM-to-DCM reproduction cross-check plus experimental volume/curvature overlays; no new oracle.
  - Cell volume drops only ~5% (−4.5 to −6%) across cytokinesis, MDCK n=22 confocal (measured, Fig 7a)
  - Equatorial ring constriction at ~5% volume loss yields ~60% axial elongation; ~40% loss needed to abolish (BD DCM-lineage sim, Fig 7d,e)
  - Membrane area reservoir expands +22% during division (sim, Fig S4d); κ_a swept 10⁻⁴/10⁻⁵/10⁻⁶ N/m
  - FEM-inferred division-elongation force 300–500 nN (best SSIM, Fig 5h,i)
  - Cell-cell boundary curvature 0.035→0.092 µm⁻¹ metaphase→cytokinesis (measured, Fig 4g,h)
  - Perpendicular stress <25% of axial during elongation (~75% in rounding/spreading), measured p3
  - Cell elongation 13.4% (meta→late ana) + 19.5% (late ana→cyto), n=94 (Fig 6b)
  · n_new_gates=4
  · n_tighten=0
  · n_params=2
  · paper_kind=experimental

## p22 [4] Taeyoon Kim · research — Cellular Pushing Forces during Mitosis Drive Mitotic
- headline: Kim/Chaudhuri mitosis-in-collagen paper: strongest new content is a measured G'(concentration) curve (10/100/400 Pa) extending VG-U1-G0, and a DCM division-force protocol (ring+spindle 100 pN/side → 40% elongation, volume-conserved) — no new mechanism (Bell's law is a duplicate); E_fibril/ell_p conflict-flags for PI; full sim block blocked on absent Table S1.
  - Collagen G' = ~10/100/400 Pa at 1/3/5 mg/mL (measured, AR-G2 rheometer)
  - Mitotic cell elongation ~40% (norm 1.4) at 1 mg/mL, MDA-MB-231 (exp + sim)
  - Ring-only elongation ~20% vs ring+spindle ~40% (additive-mechanism decomposition)
  - F_Cyto (ring) ~100 pN/side and F_Is (spindle) ~100 pN/side (sim-fitted inputs, never-fit)
  - Cell volume conserved through mitosis (n.s. Meta/Ana/T-C) — volume-constraint invariant
  - Collagen fibril E ~30 MPa (Graham 2004) + fibril radius 6.5 nm + specific volume 0.73 mL/g
  - Max matrix deformation ~0.4-0.8 um, control vs blebbistatin 0.75 vs 0.4 (ring pushes)
  - Fiber persistence length ell_p ~100 um (model bending input)
  · n_new_gates=5
  · n_tighten=1
  · n_params=8
  · paper_kind=experimental

## p14 [4] ? · ? — Reconstitution of actomyosin networks in cell-sized
- headline: Taeyoon Kim's canonical fine-grained actomyosin+deformable-membrane model (reconstituted liposome blebbing) confirms our Young-Laplace + Bell-Evans oracles, adds a new Bell+Laplace cortex-detachment oracle, and supplies a fully-specified native-scale FF/DCM parameter set — all cross-check/overlay grade (minimal reconstituted system + simulation, not MCF-7).
  - Analytic cortex-detachment oracle: γ_C^D* ≃ 0.35·(k_BT·ρ_H·R/r₀)·ln(ρ_A/K_m^2D), γ_C^D* ∝ R_C (Bell + Young-Laplace, Eqs 6–8)
  - F-actin severing threshold F_sev = 300 pN (range 200–600, Tsuda 1996) — FF severing gate candidate
  - Actin persistence length l_p ~9 µm (κ_b,A = 2.64e-19 N·m, Isambert 1995) — FF actin-stiffness anchor, distinct from collagen ell_p 17 µm
  - Reconstituted cortex network tension plateau ~0.2 mN/m (200 pN/µm; sweep 0–400 ∝R_C) — corroborates interphase γ~0.15 mN/m, overlay-only
  - NM-II minifilament length 0.64 ± 0.30 µm (TIRF, n=550; 294 nm modeled) + per-arm stall 5.7 pN, v₀ 140 nm/s (Erdmann-Schwarz PCM)
  - K_d His-α-actinin/F-actin = 1.4 µM (co-sedimentation); α-actinin rupture force 1.4–80 pN (Ferrer/Miyata) — consistent with our Ferrer-2008 k_off⁰ anchor
  - Young-Laplace ΔP ≃ 2γ_c/R confirmed in silico (γ_c≫γ_m near bleb threshold); extended with membrane term 2(γ_m+γ_c)/R (Eq 2)
  · n_new_gates=2
  · n_tighten=1
  · n_params=6
  · paper_kind=experimental

## p07 [4] Taeyoon Kim · research — Covalent cross-linking of basement
- headline: Taeyoon-Kim ECM paper opens a real FF validation gap — matrix VISCOPLASTICITY (creep-recovery permanent strain, 80%->10%) decoupled from stiffness — uncovered by any current gate; 2 FF gate candidates staged, crosslinker off-rates are simulation-origin (PI-flagged, not registered).
  - rBM Young's modulus ~100 Pa, ns with tTg (8 mg/mL Matrigel) — CONFIRMS VG-U1-G0 upper band
  - rBM creep-recovery permanent strain ~80% -> ~10% with 500 ug/mL tTg (plasticity drops, stiffness invariant)
  - IPN permanent strain ~30% (HP) -> <10% (LP-CC covalent) at ~2 kPa Young's, ns stiffness
  - Plasticity vs tTg dose = exponential decay, R2=0.95
  - Kim-model crosslinker Bell slip-bond reference: k*0,u=3e-6 s^-1, lambda*u=100 pm (SIMULATION-origin, not registered)
  - Model creep permanent strain ~0->~70% as k0,u rises; initial modulus ~500-750 Pa flat (stiffness-plasticity decoupling)
  - Model invasion max/perm strain: HP ~40%/~30% vs LP ~10%/~0%; protrusion ~1 nN
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p67 [3.5] Taeyoon Kim · research — Cytoskeletal Deformation at High Strains and the Role of Cross-link
- headline: Filamin-A single-molecule rupture (70±23 pN unbind / 57±19 pN unfold) + Bell slip-bond params fill the currently-empty FF cross-linker validation slot; no existing gate covers it.
  - Filamin-A unbinding critical force 70 ± 23 pN (single cross-link, LR 400–2000 pN/s)
  - Ig-domain unfolding critical force 57 ± 19 pN (same OTFS assay)
  - Unfolding sawtooth period 28 ± 5 nm (theory ~34 nm, 96 aa × 3.5 Å)
  - Rupture force drops sharply above ~45° inter-filament pulling angle (α-actinin ~70°)
  - Bell slip-bond sim inputs k_off0 = 0.115 s⁻¹, x = 0.416 nm
  - HS-fit k_off0: 1.70 s⁻¹ unbind / 0.52 s⁻¹ unfold (high LR); 0.087 s⁻¹ low LR
  - Reconstituted network discrete transition forces 20–56 pN
  - Barrier heights ΔG 5.6 (unbind) / 6.8 (unfold) k_BT via ½κ_m x‡²
  · n_new_gates=2
  · n_tighten=0
  · n_params=5
  · paper_kind=experimental

## p59 [3.5] Makito Miyazaki · research — Cell-sized spherical conﬁnement induces the
- headline: Experimental actomyosin-ring reconstitution: net-new validation lane (no existing gate overlaps) with 4 FF-engine emergent-benchmark gate candidates, 3 analytic-oracle seeds (duty^2 effective-conc law, constant-volume invariant, elastic-rod confinement), and a distinct actin Lp param; all overlay/cross-check, not MCF7 fit targets.
  - Fast-phase contraction rate 0.39 um/min per um perimeter (matches cytokinetic rings in cells), constant ring volume
  - Critical contractable perimeter P* ~4.5 um ~ mean actin filament length 5.4 um; velocity proportional to (perimeter - P*)
  - Critical confinement diameter R* ~15 um ~ actin persistence length Lp (equatorial ring vanishes above)
  - Ring/droplet diameter ratio R_ring/R_droplet = 0.96 +/- 0.05 (ring hugs inner equator)
  - Effective myosin concentration law [myosin]_eff = [bulk] x (duty ratio)^2, explains ~100x myoV vs myoII offset
  - Actin persistence length Lp = 5-18 um (distinct from registered collagen ell_p = 17 um)
  - Myosin-II head duty ratio 0.04-0.05; myosin-V 0.5-0.7
  - Constant-volume contraction invariant V = pi*R_ring*(pi*d^2/4) = const
  · n_new_gates=4
  · n_tighten=0
  · n_params=4
  · paper_kind=experimental

## p26 [3.5] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: One genuinely new FF mechanism — tension-induced actin severing at F_frag=500 pN (Tsuda 1996) with a ≥2× fragment-retraction diagnostic — plus a per-head-vs-per-motor stall clarification (5 vs 2 pN); all contraction/severing curves are sibling agent-based-model output (comparison-only) and Table S1 provenance is missing.
  - F_frag = 500 pN tensile actin-severing threshold (Tsuda 1996) — new FF mechanism, no current gate/param
  - Myosin single-head stall f_stall ≈ 5 pN vs our PARAM-F_stall_motor default 2 pN 'per motor' (per-head vs per-motor conflict)
  - Fragment retraction speed jump ≥2× (S_after/S_before>2) in 8/10 events — tension-release diagnostic (Fig 5c)
  - ~300 pN avg tension (subset >500 pN) on ⟨Lf⟩=2.6 µm filament, ~60 heads half-bound — worked estimate p.12
  - Contraction biphasic vs ATP: optimum ~100 µM, ~0 at 1000 µM (in-vitro, Fig 3c) — overlay only
  - Bipolar minifilament geometry: 8 arms × 8 heads = 64 heads (model, Table S1 absent)
  - Gelsolin-triggered contraction 2→120 µm²/s on stable R_ANLN=0.24 network (Fig 4c) — overlay
  · n_new_gates=2
  · n_tighten=1
  · n_params=2
  · paper_kind=experimental

## p13 [3.5] Taeyoon Kim · research — Swirling motion of breast cancer cells radially aligns collagen fibers
- headline: Independent Taeyoon-Kim BD fiber network confirms + anisotropizes our ~1/r force-propagation gate (radial n=1.2 vs tangential n=2.2) and adds a negative-normal-stress sign contract for the FF collagen ECM; experimental/tissue numbers are overlay-only, model params are reference-not-empirical.
  - Radial deformation decay exponent n=1.2 vs tangential n=2.2 (annular dipole-shear Kim BD fiber sim, Fig.4p) — confirms + anisotropizes VG-U1-force-propagation ~1/r
  - Negative normal stress |sigma_N| ~= sigma_shear (contractile, negative) in HP col1 acIPN + IBC tissue at 20% strain; positive in healthy/LP (Fig.4c-f)
  - Alginate-only control: radial deformation decays as fast as tangential when fibrous col1 removed (Fig.5e-f) — long-range ~1/r is a fibrous-network property
  - Human breast tissue Young's modulus 1.4 kPa healthy / 4.7 kPa IBC; plasticity 15% / 32% (Fig.1d,f) — overlay context only
  - Kim fiber-model constants: kappa_s,f=4e-3 N/m, kappa_b,f=8.27e-20 N.m, Bell k0_off=1e-6 s^-1, lambda=0.4 nm (Suppl. Table 2) — model reference, not empirical
  - Interface cell velocity ~0.05-0.2 um/min, tangential >> radial (swirling); 4T1 intravital anisotropy p=0.0274, abolished by Ecad KD (Fig.3)
  · n_new_gates=2
  · n_tighten=1
  · n_params=5
  · paper_kind=experimental

## p75 [3.0] Taeyoon Kim · research — Nature  |  Vol 626  |  15 February 2024  |  635
- headline: TAG-validation lens — p75  (Fan et al. 2024, "Matrix viscoelasticity promotes liver cancer", Nature 626:635)

## p69 [3] Taeyoon Kim · research — Comp. Part. Mech. (2015) 2:317–327
- headline: Kim-lab agent-based actomyosin simulation: high FF-cortex mechanistic relevance but low validation authority (peer sim, not experiment) — yields filamin-A/myosin param provenance and qualitative sigma_max~R_M^0.6 + crosslinker-clutch trend cross-checks, zero hard gate bands.
  - Filamin-A crosslinker off-rate k0_u,ACP = 0.115 s^-1, lambda = 1.04e-10 m (Ferrer 2008, single-molecule) — feeds Bell-Evans oracle
  - NMII thick-filament stall force f_M^stall ~5.7 pN and unloaded walking velocity ~140 nm/s (vs PARAM defaults 2 pN / 100 nm/s)
  - Peak network contractile stress sigma_max ~600 Pa plateau (R_M=0.014, R_ACP=0.1; motor/stall-limited) — sim comparison only
  - sigma_max scales as R_M^0.6 for R_M<0.06, then rolls off (Fig 5a) — qualitative FF cross-check
  - Sustainability S -> 1 as crosslinker density R_ACP rises (crosslinkers as molecular clutch stabilizing stress)
  - Storage modulus E' ~ f^0.09; critical frequencies 1 / 3.16 / 10-31.6 Hz across conditions; E' span 1e2-1e5 Pa
  - 64 heads per non-muscle myosin II thick filament (8 arms x 8 heads) — FF minifilament structural param
  · n_new_gates=0
  · n_tighten=0
  · n_params=5
  · paper_kind=simulation

## p64 [3] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:1143–1155
- headline: Kim-lab BD crosslinked-actin simulation — FF engine cross-check + ACP/actin parameter donor; 3 candidate FF gates (strain-stiffening γc~0.5, crosslinker-stiffness gating law, actin lₚ), but simulation-tier (comparison-only, no new oracle, no experimental overlay except imported Kang2012 lₚ)
  - F-actin lₚ ~3.5–4 µm baseline (sim; LOW vs canonical 9–17 µm — flag, do not adopt)
  - Actin segment extensional kₛ ~0.055–0.075 N/m (RTB-derived, FF param donor)
  - ACP crosslinker stiffness: compliant 0.002 N/m (filamin/α-actinin) vs rigid 0.2 N/m (scruin)
  - Network strain-stiffening critical strain ~0.5 (corroborates VG-U1 γ window)
  - Network shear stress @60% strain: ~200–230 Pa soft ACP, ~800–2500 Pa rigid ACP
  - Crosslinker-stiffness gating law: soft ACP insensitive (p=0.88), rigid ACP sensitive (p=0.03), ~3× stress spread over lₚ 2.1→12.7 µm
  - Kang2012 EXPERIMENTAL anchor: lₚ 2.1→12.7 µm over [MgCl₂] 0.5→5 mM (only exp. datum, overlay-tier)
  · n_new_gates=3
  · n_tighten=0
  · n_params=9
  · paper_kind=simulation

## p63 [3] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:345–355
- headline: Same-lineage (Taeyoon Kim) FF-class actomyosin-bundle simulation: rich comparison-only cross-checks for an FF bundle-contractility test (χ, F_max~R_M^0.6, sustainability) plus one analytic Euler-buckling oracle; Ferrer-2008 ACP off-rate 0.115 s⁻¹ conflicts with our 0.066 re-anchor — surface to PI.
  - Bundle F_max ~2 nN sustained ~100 s at R_M=0.025/R_ACP=0.032 (S≈1); motors-only collapses <5 s
  - Contractile efficiency χ ≈ 0.005 (no ACP) → ≈0.15 (with ACP), ~30-40× crosslinker amplification
  - F_max ∝ R_M^1 without ACPs → ∝ R_M^0.6 with ACPs (exponent change = clutch signature)
  - Motor per-head stall f_M^stall ≈ 4.5 pN, unloaded velocity ≈ 550 nm/s (PCM output)
  - ACP off-rate k0*_u,ACP = 0.115 s⁻¹, λ* = 1.04e-10 m (Ferrer 2008) — CONFLICT with KB 0.066 s⁻¹ re-anchor
  - Euler buckling oracle F = π²·κ_b,A·r0,A/L_seg² (Eq.7, symmetry-breaking analytic law)
  - Substrate biphasic-tension transition ~300 Pa (saturating, not peaked — confirms KB-2.8)
  - Parallel vs anti-parallel ACP efficiency η ≈ 0.61 vs 0.15
  · n_new_gates=3
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p62 [3] Makito Miyazaki · research — Processive Nanostepping of Formin mDia1 Loosely Coupled with
- headline: FF-only single-molecule formin mDia1 elongation data: 2.7-2.8 nm ≈ one-actin-subunit step invariant (strong sanity gate) + concentration-dependent stepping rates 0.92-1.4 s⁻¹; overlay/cross-check only, no new oracle, rates are near-critical-conc calibration not physiological setpoints.
  - Formin mDia1 forward unitary step +2.8 nm / backward -2.4 nm (≈ one 2.7 nm actin subunit), 4-6 pN, 26°C (Fig 3b)
  - Forward stepping rate k = 0.92 s⁻¹ (40 nM) → 1.4 s⁻¹ (75 nM), concentration-dependent (Fig 3d)
  - Backward stepping rate k = 0.56/0.62 s⁻¹, concentration-INdependent (Fig 3d)
  - Multi-subunit steps: fwd 2×=4.9 nm, 3×=8.7 nm; area fractions 60/38/2.2% (40nM) → 45/46/9.3% (75nM) = loose coupling
  - Actin critical concentration 64 nM (both ends, Fujiwara 2002); assay [G-actin] 40-75 nM ≪ physiological
  - Applied tension 4-6 pN, elongation ~force-independent (not a force-velocity law)
  · n_new_gates=3
  · n_tighten=0
  · n_params=6
  · paper_kind=experimental

## p60 [3] Taeyoon Kim · research — Interplay of active processes modulates tension
- headline: Taeyoon-Kim BD paper (FF engine's own physics class): actin turnover stabilizes the cortex against motor-driven aggregation and sets sustained prestress — a strong FF emergent-behavior comparison target plus experimental foci time constants (tau_C=28s, tau_D=155s), but nearly all magnitudes are simulation output (comparison-only, no new bands or oracles).
  - Peak internal prestress ~300-400 Pa at turnover 0 s^-1, decaying monotonically to ~0 Pa at 300 s^-1 (3D BD, 25 uM actin, 5% ACP, 1% motor) — SIMULATION, comparison-only
  - Critical turnover rate gates homogeneous<->aggregated phase transition; ->0 at ACP:motor ~=100, rises with motor% (sim)
  - Cluster-formation tau_C = 28 +/- 5 s (R2 0.77), MDA-MB-231, 2.5 uM Cytochalasin D, n=15 — EXPERIMENTAL overlay
  - Declustering tau_D = 155 +/- 42 s (R2 0.83), CytoD washout, n=12 — EXPERIMENTAL overlay
  - Pulsed actomyosin-foci / constriction timescale ~100 s (sim + Drosophila lit)
  - Rate-ratio master curve: peak-stress-vs-turnover collapses when turnover scaled by relative motor-walking rate (dimensionless, sim)
  - Physiological traction-stress reference ~100-1000 Pa (TFM, lit refs 6,48) — context, not fit target
  · n_new_gates=4
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p50 [3] ? · ? — Molecular Biology of the Cell • 34:ar67, 1–10, June 1, 2023
- headline: Taeyoon Kim parallel-focal-adhesion motor-clutch paper: no new numeric bands (mechanism already in MC-U2, constants SI-locked/NO_SI), but 3 FF emergent cross-checks — ~37 FAs/lamellipodium overlay, width→force↑/stress↓ monotone law, and PFAC parallel-adhesion lifetime rescue that mirrors the DCM cadherin load-sharing cluster.
  - FAs per lamellipodium: mean ~37, max 89, min 6 (SIM-TIRF, H4 glioma) — FF adhesion-density overlay
  - Mean force ~30→110 pN rises while substrate stress ~1700→200 pN/um2 falls monotonically as adhesion width A_FA 0.05→1 (agent-based)
  - PFAC: single filopodial adhesion load-and-fails; ≥10 parallel adhesions prolong lifetime, NO failure at 30 FA (Fig 4H)
  - FA lifetime ~11 s (A_FA 0.1) → ~2–3 s (wide lamellipodium) — seconds-order nascent FA (Fig 6F)
  - TFM traction control ~400 Pa vs Arp2/3-inhibited ~500–600 Pa on 90 kPa (overlay only, H4 glioma)
  - Migration speed deficit Cohen's d 0.9 (500 Pa) → 0.5 (4 kPa) → 0.25 (90 kPa) — overlay only
  - Arp2/3 branch angle ~70° confirms existing angle-harmonic convention
  · n_new_gates=3
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p44 [3] Taeyoon Kim · research — Morphological Transformation and Force
- headline: Kim-lab agent-based actomyosin network→bundle simulation: comparison-only for FF (its own tensions/times are model output), but supplies clean sourced params — filamin-A Bell slip-bond (Ferrer), NM-II PCM f_stall 5.7 pN / v 140 nm/s conflicting current defaults, actin L_p 9 µm.
  - Actin persistence length L_p = 9 µm (Isambert 1995) — actin bending anchor, distinct from collagen 17 µm
  - Filamin-A ACP slip-bond: k0_off = 0.115 s⁻¹, λ = 1.04e-10 m (Ferrer 2008) — Bell-Evans oracle datum for a cortex crosslinker
  - NM-II myosin f_stall ≈ 5.7 pN (PCM) — CONFLICTS PARAM-F_stall_motor default 2 pN (~3×)
  - NM-II unloaded velocity ≈ 140 nm/s — EXCEEDS PARAM-v_unloaded range 10–100 nm/s
  - NM-II reference unbinding k_u = 0.049 s⁻¹ (PCM catch-bond)
  - Network→bundle compaction ~10–20 s (sim), anchored by experimental transverse-arc ~20 s (Heath 1983)
  - Emergent bundle max tension ~0.8 nN (fail) vs ~4 nN (success) — SIM OUTPUT, FF comparison only
  - Buckling criterion end-to-end/contour < 0.6; suppressing buckling (100×κ_b,A) abolishes isotropic-network tension
  · n_new_gates=3
  · n_tighten=0
  · n_params=6
  · paper_kind=simulation

## p43 [3] Taeyoon Kim · research — Computational Analysis of Viscoelastic Properties of
- headline: Kim-2009 crosslinked-actin BD opens an uncovered FF network-rheology surface: the analytic ω^3/4 semiflexible exponent is a clean oracle candidate and the in-vitro G',G'' 0.1–10 Pa curve an overlay target, but all prestress/magnitude numbers are simulation-grade and rest on ÷40 stiffness + permanent (unbinding-free) crosslinks — comparison-only.
  - G'∝f^0.75 high-freq exponent, crosslink-free F-actin (analytic semiflexible ω^3/4; EV on)
  - G'∝f^0.3 low-freq exponent, well-crosslinked (R=0.021), cell-like 0.15–0.3
  - In-vitro G',G'' ≈ 0.1–10 Pa across 1–10 Hz (filamin-A/gelsolin, ⟨L_f⟩=1.5µm, C_A=12.1µM) — overlay-only
  - Prestress stiffening: G'∝τ₀^0.85 above threshold τ₀≈0.1 Pa (sim; softens to ~0.7 for soft k_s,A; in-vitro ref ~1.0)
  - G' ~100× enhancement at prestrain γ=0.55, near frequency-independent
  - Supportive framework: top-25% ACP (~28% filaments) carries ~70% of stress
  - Crosslink spacing l_c=0.393 µm at R=0.021; thermal fluctuations negligible unless l_p ≲ 3 µm
  - k_b,A=1.056e-18 N·m → actin l_p≈10–20 µm (physiological, non-scaled); k_s,A used is ÷40 unphysical
  · n_new_gates=3
  · n_tighten=0
  · n_params=3
  · paper_kind=simulation

## p42 [3] Taeyoon Kim · research — F‑Actin Fragmentation Induces Distinct Mechanisms of Stress
- headline: Kim-model actin-severing simulation: contributes a canonical Euler-buckling analytic oracle for FF and a directionality-of-cyclic-relaxation behavior gate, but its stress magnitudes are simulation-only (comparison, not fit); severing is a mechanism gap absent from FF (SI needed for a port).
  - Euler buckling F_buckle = π²·κ_b,A·r_0,A/L_c² = 4.1 pN at L_c=300 nm (analytic, eq 2) — candidate FF oracle
  - Peak shear stress ~100–120 Pa at 50% cyclic strain, R_ACP=0.032; ~150 Pa at R_ACP=0.1 (sim, 3µm cube)
  - Cyclic stress drop rises 0→~75% with amplitude 0.1→0.5 bidirectional; ~0% unidirectional (sim)
  - Max cyclic stress drop up to ~70% for long ⟨L_f⟩ + high R_ACP (sim)
  - Severing orientation ~135° vs high-load (≥8 pN) ACP unbinding ~45° — orthogonal, independent channels
  - Bell off-rate rises only 1.1× k_u⁰ at 4.1 pN — confirms existing Bell-Evans oracle
  - First-max stress spans >3 orders (0.1→>100 Pa) across ⟨L_f⟩ 0.7–1.3µm × R_ACP 0.005–0.05 sweep
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p41 [3] Taeyoon Kim · research — Ding, Chou et al. eLife 2025;14:RP105236. DOI: https://doi.org/10.7554
- headline: Kim-lab agent-based (simulation) actomyosin paper: FF's near-direct sibling — stages 2 FF cross-check gate candidates (cooperative-overlap force bound; bundle∝N_M vs network∝√N_M scaling) and myosin-architecture param anchors, but as a simulation it is comparison-only, not a truth oracle.
  - Two-motor cooperative-overlap force bound 0.5·F_M^max ≤ Ftot ≤ F_M^max (Eq 1): no-ACP→2× single motor, ACP-between→1×
  - Bundle force scaling Ftot ∝ N_M (slope≈1) vs 2D network Ftot ∝ N_M^0.65 (≈√N_M)
  - Efficiency η ∝ 1/N_M with high-N_M plateau ≈0.08; η=0.25→0.12 as bundle thickens (N_F 2→7)
  - Single myosin-head stall F_st≈5.7 pN, N_h=8 heads/arm → arm stall ≈45.6 pN; unloaded v≈140 nm/s (NMII, PCM)
  - Critical overlap length L_c = 2·L_sp·(N_a/4−1) (Eq 6); analytic force estimator F_est = F_st·N_h·N_a·Ξ/2 (Eq 7)
  - Bundle Ftot = 0.6/1.3/2.9/3.9 nN for N_F=2/4/6/7 (R_ACP=0.04, R_M=0.005, N_a=24)
  - Actin persistence length L_p=9 µm (distinct from KB collagen ell_p=17 µm); NMII minifilament 0.3–1.5 µm, bare zone ~160 nm
  · n_new_gates=2
  · n_tighten=0
  · n_params=7
  · paper_kind=simulation

## p39 [3] Taeyoon Kim · research — 10274 |  Soft Matter, 2021, 17, 10274–10285
- headline: Kim-lab agent-based cell-in-viscoelastic-ECM sim: CONFIRMS the passing r^-1 stress-decay gate (VG-U1) and supplies its 2-D/3-D geometric rationale + 2 ECM-crosslinker Bell params; relaxation/magnitude curves are sim-only comparison, ESI numbers absent.
  - Far-field transmitted radial stress decays as sigma proportional to r^-1 (2-D cylindrical slab) / r^-2 (3-D) — analytic, geometry = force/lateral-area
  - Peak radial stress ~55-60 Pa at r~20 um (strongest contraction); ~5-10 Pa weakest; peak linear in k_s,c (Kim-engine sim)
  - Peak stress independent of crosslinker off-rate k_ub,0 (stress peaks before unbinding); >50% relaxation even at zero unbinding (poroelastic split)
  - Stress relaxation ~80% (low fiber conc 0.5x) vs ~20% (high 2x) over 10 min/1 h; relaxation rate proportional to R_xl^-1
  - ECM crosslinker Bell slip-bond: k_ub,0 = 1e-5 s^-1, x_ub = 1.0e-10 m (0.1 nm), benchmarked Nam 2021 (ESI Table S1)
  - Tensed fibers align radial (~0 deg), buckled fibers align circumferential (~90 deg) under central contraction
  - In-vitro 3T3/collagen (overlay-only): remodeling decreases 1->5 mg/mL; reduced by tTG 500 ug/mL (p<0.05) and by ML-7 25 uM / blebbistatin 10 uM (p<0.0001)
  · n_new_gates=2
  · n_tighten=1
  · n_params=2
  · paper_kind=simulation

## p24 [3] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: FF-relevant Kim-lab actomyosin simulation: registrable literature anchors (3.8 pN/head stall, filamin-A k_off 0.115 s^-1, NMII <140 nm/s) but headline curves are simulation-output/overlay-only; stall value challenges the 2 pN param and informs the F*=7 pN debt.
  - Myosin-II stall force 3.8 pN/head (Erdmann/Schwarz PCM) — CHALLENGES PARAM-F_stall_motor default 2 pN
  - Filamin-A cross-linker zero-force off-rate k0_u,ACP = 0.115 s^-1 (Ferrer 2008); distinct from KB's alpha-actinin 0.066 s^-1
  - NMII unloaded walking velocity < 140 nm/s (Cuda 1997) — CONFIRMS PARAM-v_unloaded 100 nm/s
  - SKMM walking velocity < 5 um/s and ~4% duty ratio (Harris & Warshaw 1993)
  - Motor MSD anomalous exponent: 1<alpha<=2 superdiffusive with F-actin turnover ON, alpha<1 stalled in high-connectivity/no-turnover limit
  - F-actin turnover thresholds erasing confinement (sim): k_t,A=20 s^-1 local stalling, up to 70 s^-1 aggregation (comparison-only)
  · n_new_gates=3
  · n_tighten=1
  · n_params=6
  · paper_kind=simulation

## p21 [3] Taeyoon Kim · research — Dr. V. Yadav, Dr. A. P. Tabatabai, Prof. M. P. Murrell
- headline: Reconstituted formin/actin "model cortex" gives FF-scale cross-check data (150/10 sub/s elongation, l_p≈15 µm, FRAP τ 1.7–21.8 min) and an analytic filament-bending-energy oracle; confirms but does not tighten the G-actin-pool contract — all overlay-only, not MCF7.
  - formin(mDia1)-bound elongation ~150 subunits/s (≈0.37 µm/s) vs formin-free ~10 subunits/s — two-mode barbed-end kinetics (Fig 1B,C)
  - Filament length set by nucleation: ~100 µm @10 nM formin → ~5–10 µm @1 µM formin (Fig 1F,G)
  - F-actin FRAP turnover τ₁ = 1.7 / 21.8 / 12.6 min at 10 / 100 / 1000 nM formin — non-monotonic (Fig 2D)
  - Actin persistence length l_p ≈ 15 µm (Gittes 1993), distinct from collagen ell_p 17 µm
  - Filament-arc bending energy E=κL/2R², κ=k_BT·l_p → 10⁻²⁰–10⁻¹⁹ J (~100× k_BT) — analytic oracle candidate (Eq S2, Fig 3C)
  - Biphasic network growth crossover R₁→R₂ at ~20 min from finite shared G-actin pool depletion (Fig 1E) — confirms VG-U3-gactin-pool
  - Frank bend constant K₃₃ = 1.38–4.6 pN for the actin bundle director field (Fig 3B,D)
  · n_new_gates=3
  · n_tighten=0
  · n_params=3
  · paper_kind=experimental

## p19 [3] Taeyoon Kim · research — Mechanical Checkpoint for Cell Division in
- headline: Taeyoon-Kim ECM division-checkpoint sim: numbers are model output (comparison-only), but yields one solid FF analytic oracle (fiber Euler-buckling ~30 pN cap), a DCM mitotic-division load protocol, and a flag that its ℓ_p~100µm conflicts with KB PARAM-ell_p=17µm.
  - Fiber compressive-buckling cap ~30 pN (analytic Euler load, ∝κ_b,f/L²) — polar proximal fibers, mitotic elongation (Fig.4E/4I)
  - Cell aspect ratio 1.41 (2D) / ~1.75-2.0 (3D, extend ~60%), stop at 58.74% elongation (sim, Fig.2G/5B)
  - Anchoring ratio T_mat/T_pole ~0.37 reference (checkpoint engaged); <0.5 → failure (sim, Fig.3F)
  - Shell tension T_pole~6.5 / T_mat~2.5 nN (2D), T_pole~10 nN (3D); shell area +10%(2D)/+20%(3D) (sim)
  - Fiber persistence length ℓ_p ~100 µm (fiber/bundle) — CONFLICTS with PARAM-ell_p=17 µm (fibril, KU-1.1)
  - Bell slip-only cross-linker unbinding k₋=k₋⁰·exp(λ|F|/kBT); ×100 turnover → T_mat/T_pole→0 (Eq.7) — DUPLICATE of KB Bell-Evans oracle
  - Borrowed constitutive anchors: cell E=324 Pa (Eldridge2019), collagen E_mat=832 Pa/ν=0.313 (Lane2018) — overlay-only
  - Experimental MDA-MB-231/1mg-mL collagen PIV radial deformation ±4 µm (Fig.1B, the only measured datum)
  · n_new_gates=2
  · n_tighten=0
  · n_params=3
  · paper_kind=simulation

## p17 [3] ? · ? — Kinetic Control of Out-Of-Equilibrium
- headline: Cortical actomyosin turnover kinetics (myosin k_off 0.025-0.10 s^-1, k_app, ~4.5 s cascade delay) plus two closed-form analytic oracles (first-order RC leaky-integrator k_app/k_off and high-pass sqrt(1+w^2/k_off^2)) that match FF myosin bind/unbind physics; all experimental numbers are cross-species (C. elegans) overlay-only, no existing gate overlaps.
  - Myosin II cortical off-rate k_off ~0.025-0.10 s^-1 (C. elegans embryo, per pulse, Fig 2E)
  - Myosin k_app ~0.5-1.5 molecule/s; F-actin k_app 5x swing 0.5->2.5 (Fig 2D/4C)
  - Analytic high-pass turnover attenuation sqrt(1+omega^2/k_off^2) (Eq.4) + first-order step relaxation N(t)=N_eq2+(N_eq1-N_eq2)e^(-k_off t) (Eq.7)
  - RhoA/ROCK/formin->Myosin delay ~4.5 s; F-actin precedes myosin ~0.5 s; ~4.25 s is turnover-kinetic (Fig 1/2/4)
  - Pulse period 31 s WT / 42 s cofilin-KD (unc-60 RNAi) (Fig S4F)
  - Agent-based net contractile force non-monotonic vs delay, peak ~38 nN at d_M~5 s over 20x20x0.1 um slab (Fig 5D, sim)
  - Formin barbed-end elongation 1.2 um/s -> ~12 um in ~10 s (Methods, Costache 2022 anchor)
  · n_new_gates=3
  · n_tighten=0
  · n_params=4
  · paper_kind=experimental

## p15 [3] Taeyoon Kim · research — Dissecting Molecular Origins of the Mechano-Adaptive
- headline: Canonical Taeyoon-Kim agent-based actomyosin-bundle model + isolated single-SF micromanipulation: exact FF-engine physics (explicit actin beams, ACP Bell crosslinkers, multi-head myosin PCM catch-bond). Rich FF parameter set + 3 candidate FF bundle-buckling/recovery cross-check gates; data is sibling-simulation (comparison-only, not truth), experiment overlay-only.
  - Myosin minifilament: unloaded v0 ~56 nm/s, stall F_s ~24 pN (=6*N_h, 32-head, N_a=8/N_h=4), 30-pN slip cutoff (PCM, Fig.S6/p40)
  - Rate-dependent buckling law: bundle buckles iff compression rate > unloaded motor walking speed; no buckling at |ec|=0.01/s (Fig.4, exp Fig.1)
  - Recovery time set by motor WALK SPEED not density: k20 4->16/s gives t_r ~32->4.5 s; R_M 0.002->0.01 gives t_r ~14->12.5 s flat (Fig.5B/E)
  - Steady-state contractile force scales with motor density: F_B,ss ~0.3/0.82/1.38/1.62 nN for R_M 0.002/0.004/0.008/0.01 (Fig.5F)
  - Actin bending kb,A=2.64e-19 N*m (L_p=9 um); ref sims deliberately used 10x (L_p~90 um) (Table S1/p38)
  - ACP Bell off-rate k0_u=0.05/s, x_u=0.9 nm; generic crosslinker (Table S1)
  - Peak stretch force F_B ~5.5 nN, relaxes to <half before stretch end via motor/ACP unbinding (Fig.2C)
  - Cylinder-segment drag closed form: zeta=3*pi*mu*r_c*(3+2r0/r_c)/5, mu=8.6e-2 Pa*s (Eq.S3/p37)
  · n_new_gates=3
  · n_tighten=0
  · n_params=9
  · paper_kind=simulation

## p10 [3] ? · ? — Fine-tuning of material properties by catch bonds
- headline: Sibling-simulation (Kim/Cytosim-lineage) of α-actinin catch-slip actin crosslinkers: strong FF comparison target + a genuine catch-slip lifetime oracle gap (Eq.1), but all numbers are simulation-output comparison-only, not fit targets; actual force constants need SI Table S1.
  - Catch-slip reference lifetime: τ0=9.5 s, τmax=53 s, Fmax=14.7 pN (α-actinin-like crosslinker, area-matched to slip; Fig.1B) — peak force distinct from integrin ~30 pN
  - Single-motor contractile stress: catch-slip ~90 Pa vs slip ~25 Pa → ~3-4× (Fig.6B, 512-head myosin-II)
  - Passive shear rheology: slip peaks ~100 Pa, catch-slip ~200 Pa still rising at γ~1.5 (Fig.1C; γ̇=0.001 s⁻¹, CF=10µM, RCL=0.01)
  - Catch-slip beats slip only at high connectivity: σ advantage >1 (up to ~1.3×) at high RCL, <1 at RCL≈0.002; strain advantage up to ~2.5× (Fig.2)
  - Eq.1 catch-slip lifetime law k_u=k^s0·e^{+λ^s F/kT}+k^c0·e^{-λ^c F/kT} — closed-form oracle candidate (constants in SI Table S1, not in main PDF)
  - CL turnover essential + connectivity: 5µm hop default; no-turnover flips catch-slip below slip (Fig.S3)
  - Global CL turnover rate ~0.156-0.164 s⁻¹; critical-CL avg force ~1.7-1.8 pN vs ~1.2-1.3 global near yield (Fig.1G/H)
  · n_new_gates=3
  · n_tighten=0
  · n_params=8
  · paper_kind=simulation

## p03 [3] Taeyoon Kim · research — Determinants of Fluidlike Behavior and Effective Viscosity in Cross-Li
- headline: Canonical Taeyoon-Kim agent-model rheology paper (simulation, not experiment): opens an uncovered observable class — cortex network fluidization/effective viscosity — via a biphasic creep law (η_eff⁰≈35 Pa·s, σ_th≈1 Pa, ε̇∝turnover, 5%/99% load-bearing) that FF should reproduce as a model-vs-model cross-check.
  - Master creep law η_eff·n = η_eff⁰·max(σ,σ_th)/σ_th (Eq.6)
  - η_eff⁰ ≈ 35 Pa·s (effective-viscosity prefactor, sim)
  - σ_th ≈ 1 Pa critical/yield stress (biphasic onset)
  - Turnover scaling: ε̇ ∝ n, η_eff ∝ n⁻¹; ε̇→0 at n=0
  - Biphasic σ: ε̇∝σ below σ_th, ε̇∝σ⁰ plateau above (~2 decades)
  - ~5% of filaments carry ~99% of network tension (50 pN norm)
  - Treadmilling turnover WITH filaments → η_eff ~4× lower than pure ACP unbinding
  - Ferrer-2008 Bell ACP off-rate k⁰_ub=0.115 s⁻¹, λ=1.04e-10 m (filamin-A)
  · n_new_gates=3
  · n_tighten=0
  · n_params=4
  · paper_kind=simulation

## p02 [3] Taeyoon Kim · research — Dynamic Role of Cross-Linking Proteins in Actin Rheology
- headline: Simulation ancestor of the FF crosslinker layer — supplies lit-anchored filamin/ACP Bell-Evans slip constants (k0_ub=0.115 s^-1, lambda_ub=0.104 nm; Ferrer 2008) and the mechanistic result that unbinding, not unfolding, governs network rheology at physiological strain rate; constitutive curves are comparison-only.
  - ACP/filamin zero-force unbinding rate k0_ub = 0.115 s^-1 (Ferrer 2008, single-molecule) with Bell slip compliance lambda_ub = 1.04e-10 m (0.104 nm)
  - ACP/filamin unfolding rate k0_uf = 3.0e-5 s^-1, lambda_uf = 6e-10 m (0.6 nm); ACP arm WLC p=0.33 nm, r0=105 nm, l0,i=140+30i nm
  - Unbinding-vs-unfolding dominance crossover at gamma_dot_eff ~= 1 s^-1; unbinding governs rheology at physiological ~0.1 s^-1
  - Two-regime actin-network strain-stiffening tau = 2.5*gamma^0.4 + 500*gamma^4 Pa, differential modulus K ~ gamma^3, K_m ~ tau_m^(3/4) (simulation, comparison-only)
  - Plastic-flow exponent x ~= 0.17 at gamma_dot_eff=0.1 s^-1, matching Gardel 2006 in-vitro (gamma 0.08->0.02); high-load ACPs (>20 pN) rebind to different filament
  - Filamin-A saw-tooth unfolding peak force 100-200 pN at ~30 nm intervals (from AFM ref 7); in-vitro stress cap >100 Pa rarely seen, reorg endpoint ~10 Pa, relaxation plateau ~1 Pa
  · n_new_gates=2
  · n_tighten=0
  · n_params=6
  · paper_kind=simulation

## p58 [2.5] Makito Miyazaki · research — Myosin-Driven Advection and Actin Reorganization Control the
- headline: Confined Xenopus actomyosin gives FF cross-check numbers (myosin-driven advection ~1.2 µm/s, wave period ~101 s, boundary-shape→gel-shape transfer) but at bulk-reconstitution scale — overlay-only, 2 PI-gated candidate gates, no existing band challenged.
  - Inward actin advection speed v = 1.20 ± 0.41 µm/s (Xenopus extract, confined, PIV, n=6)
  - Contractile-wave period T = 101.0 ± 7.9 s (circular well D=350 µm, n=6)
  - Cluster roundness R=m/M: ~1.0 (W/D=0) → ~0.25–0.30 (W/D=1), size-independent
  - Cluster area set by well VOLUME not shape (n.s. across shapes at A_μ=4.8×10⁴ µm²)
  - Square-cluster quartic parameter 0.069 → 0.088 under CalA (star-like)
  - Confinement efficiency: well/gap F-actin intensity ratio 9.5×
  - Perturbation signs: CytoD↓k_p→rounder, VCA↑k_p→larger, CalA↑Pe→smaller
  · n_new_gates=2
  · n_tighten=0
  · n_params=2
  · paper_kind=experimental

## p48 [2.5] Taeyoon Kim · research — PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1013
- headline: Agent-based lamellipodium BD paper (FF's exact paradigm): outputs are comparison-only, but it independently CONFIRMS KB-2.8 traction-saturation and anchors 2 FF cross-check gate candidates (retrograde flow, flow-vs-friction) plus a v_unloaded 100→120 nm/s range flag.
  - Retrograde flow (reference steady) ~20-30 nm/s, t>100s — MODEL OUTPUT, physiological range (Fig 2E,F)
  - Total traction force PLATEAUS ~0.2-0.25 nN at motor density R_M>=0.005 — CONFIRMS KB-2.8 'saturating not peaked' (Fig 5G)
  - Flow biphasic vs depolymerization (peak ~20-25 nm/s at k-,A=4-6 s^-1) and vs severing (window 10^-5.5–10^-4.5) — qualitative laws (Figs 3F,4F)
  - NMII unloaded walking velocity 120 nm/s (Johnson 2019, Cuda 1997 — experimental) vs PARAM-v_unloaded range 10-100: range-ceiling conflict
  - Arp2/3 branch angle 70° (Amann-Pollard 2001 — experimental); likely duplicate of FF baked 70°±35°
  - Flow falls monotonically with FA-region size, frustrated at A_FA>0.4 — friction-vs-contraction competition (Fig 6F)
  · n_new_gates=2
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p37 [2.5] Taeyoon Kim · research — Molecular Biology of the Cell • 35:ar47, 1–12, April 1, 2024
- headline: Taeyoon-Kim-corpus modified motor-clutch simulation (Aplysia growth cone) — mostly confirms existing FF Bell-Evans/linear-F-V/v_u contracts; net-new value is a 6 pN NMII-stall conflict vs KB default 2 pN plus a candidate force-threshold reinforcement oracle (Mekhdjian Eq 6); all headline curves are simulation output on a non-MCF7 cell → overlay/cross-check only, no hard bands.
  - F_s single NMII stall = 6 pN (Lohner 2019) — CONFLICTS KB PARAM-F_stall_motor default 2 pN
  - v_u unloaded actin retrograde flow = 100 nm/s (Lin-Forscher 1995) — CONFIRMS PARAM-v_unloaded
  - Strong-coupling actin-flow threshold = 20 nm/s (80% reduction from v_u)
  - Latency time biphasic: 14 min @2.5 → ~2 min @4 → ~8-9 min @≥14 pN/nm (Aplysia, model vs Athamneh2015)
  - Substrate deformation Δx_sub: 1.6 µm @2.5 → ~0.8 µm plateau @≥14 pN/nm
  - Optimum substrate stiffness 4 pN/nm (min latency + min parameter sensitivity)
  - Reinforcement force threshold F_t=10 pN; K_clutch≈0.3 pN/nm & F_b≈4 pN are model fits (do NOT register)
  · n_new_gates=1
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p36 [2.5] Taeyoon Kim · research — 1548 | Soft Matter, 2020, 16, 1548--1559
- headline: FF-core actomyosin simulation (Kim-lab): net-new network-scale gliding/contractility phenomenology (biphasic speed-vs-k₂₀, ~3% crosslinker percolation arrest, intermediate-crosslink contractility optimum) — 4 candidate FF cross-comparison gates + 2 citable rate constants (k⁰_u,ACP=0.115 s⁻¹, k₂₀=20 s⁻¹); simulation output so comparison-only, never fit; exact force-law constants are in ESI (harvest before quoting).
  - Unloaded myosin gliding speed ~140 nm/s (0.14 µm/s) at k₂₀=20 s⁻¹, tuned to non-muscle myosin II (Fig 1/2b)
  - Gliding speed biphasic vs ATP-unbinding rate k₂₀ (5–640 s⁻¹); bound-motor fraction falls ~65%→~14% (Fig 2b,c)
  - ~3% of actins in crosslinks (R_ACP≈0.03) arrests network gliding at C_A=60 µM, permanent crosslinks (Fig 5a)
  - Contractility/heterogeneity maximal at INTERMEDIATE R_ACP AND intermediate k⁰_u,ACP; tensile force buckling-limited (Fig 5)
  - ACP zero-force unbinding k⁰_u,ACP = 0.115 s⁻¹ (Bell's-law reference); PCM k₂₀ = 20 s⁻¹ reference
  - Peak gliding speed ~2 µm/s at k₂₀=640 s⁻¹; motor stall force ∝ 1/k₂₀
  - Fast gliding requires ⟨L_f⟩ ≳ mean motor spacing (length-threshold scaling, Fig 1c,d)
  · n_new_gates=4
  · n_tighten=0
  · n_params=4
  · paper_kind=simulation

## p35 [2.5] Taeyoon Kim · research — Soft Matter, 2017, 13, 3213--3220 | 3213
- headline: Simulation (Taeyoon-Kim BD, same lineage as FF): adds a new severing-limited-contraction dimension the KB lacks — best value is Eq. 1 curvature-severing law as an FF mechanism/oracle plus 4 candidate FF-parity gates; all comparison-only, dimensional constants sit in the un-ingested ESI.
  - ξ_ss vs bending rigidity κ_b,A is NON-MONOTONIC — max at intermediate κ_b,A when severing is on (Fig. 3D)
  - Buckling (hence contraction) suppressed above critical R_ACP ≈ 0.32 molar ratio to actin (Fig. 2C)
  - Large contraction ξ_ss>0.50 requires R_M ≥ critical AND intermediate R_ACP window (Fig. 2C)
  - Critical buckling-cutoff κ_b,A ∝ R_ACP/R_M and independent of k⁰_s,A (Fig. 6B)
  - Eq. 1 curvature-severing law k_s,A = k⁰_s,A·exp(θ_s,A/λ_s,A) — importable FF mechanism/oracle
  - Steady-state focus count 1 → ~10–11 grows with severing; t_ss ∝ κ_b,A, ~independent of k⁰_s,A
  - 64× κ_b,A ↔ coarse-grained segment L_p 0.58 µm (model, not physiological actin L_p)
  · n_new_gates=4
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p27 [2.5] Taeyoon Kim · research — Cytoskeleton, 2026; 83:407–426
- headline: Plant (Arabidopsis, 1µM) agent-based actin paper from the Taeyoon-Kim FF lineage: rich FF filament-turnover verification targets + 3 analytic oracle laws (severing, Clift drag, FDT), but zero DCM/MCF7 contract value — all candidates plant-context, PI-gate review.
  - Barbed-end polymerization velocity: exp 1.6-1.8 µm/s (free ~1.25, formin 2-2.25) / model 2.2±0.4 µm/s — Arabidopsis cortex, C_A=1µM
  - Pointed-end depolymerization: exp 0.17-0.20 / model 0.2±0.01 µm/s
  - Severing frequency: exp 0.008-0.014 / model 0.01±0.002 events/µm/s via angle-dependent law k_sev=k0,sev·exp(λ_sev·θ), k0,sev=8e-5 s⁻¹
  - Nucleation frequency (total): exp 2e-4 / model 1.2e-4 events/µm²/s (branched 0.8e-4, de novo 0.4e-4)
  - Filament length: exp 10-17 / model 9.9±0.85 µm; formin enhancement 2-3× (α_form=3, Zhang 2016)
  - Actin bending stiffness κ_b,A=2.64e-19 N·m ⇒ l_p≈8.9 µm (Isambert 1995); extensional κ_s,A=1.69e-2 N/m; segment 140nm/7nm
  - Cylindrical-segment drag ζ=3πμ·r_c(3+2r0/r_c)/5 (Clift 1978); FDT thermal ⟨FᵀFᵀ⟩=2kBTζ/Δt (Eq 2/3)
  - 99% profilin-bound monomers; ~50% branched/20% de novo/30% formin filament classes; steady state ~500s
  · n_new_gates=4
  · n_tighten=0
  · n_params=6
  · paper_kind=simulation

## p76 [2.0] Makito Miyazaki · research — BIOPHYSICS AND COMPUTATIONAL BIOLOGY
- headline: TAG-validation lens — p76  (Sakamoto et al. 2022, actomyosin droplet motility, PNAS 119:e2121147119)

## p73 [2.0] Makito Miyazaki · research — Optogenetic actin network assembly on lipid
- headline: TAG-validation lens — p73  (Yamamoto & Miyazaki 2025, OptoVCA, Nat Commun 16:7583)

## p70 [2] Makito Miyazaki · research — Tug-of-war between actomyosin-driven
- headline: Xenopus egg-extract tug-of-war reconstitution: no single-cell gates, but a 10× α-actinin k_off conflict to surface and an FF force-percolation cross-check
  - α-actinin k_off=0.66 s⁻¹ (Wachsstock 1994) — CONFLICTS 10× with our KB anchor 0.066 s⁻¹ (Ferrer 2008)
  - α-actinin k_on=1.2×10⁶ M⁻¹s⁻¹; crosslinker turnover τ=1/(k_off+C0·k_on)≈0.54 s (C0≈1 µM)
  - F-actin contour length L≈5.7 µm control (4.3 gelsolin / 10 mDia2 fit) — CONFIRMS KB-3.18 1–10 µm range
  - Actomyosin wave period T≈46 s, contraction velocity 1–5 µm/s (droplet 150–350 µm)
  - Transition diameter D_c=85 µm measured (≈73 µm analytic, Eq.3 Rc/L=log₂(T/τ))
  - Force-percolation: longer filaments (larger L) + more crosslinkers (smaller τ) raise contractile-transmission threshold (Fig 7 sweeps L=4/8/12 µm, C0=0.1/1/10 µM)
  - Egg-extract context: [actin]~20 µM, [α-actinin]~1 µM (Wühr 2014) — reconstitution/overlay only, not MCF7
  · n_new_gates=2
  · n_tighten=0
  · n_params=6
  · paper_kind=experimental

## p66 [2] Taeyoon Kim · research — Bulletin of Mathematical Biology (2019) 81:3301–3321
- headline: Lumped agent-based PRW/CIL/nematic simulation (Taeyoon Kim) — comparison-only behavioral targets (single-cell speed ~35-60 µm/hr, MSD ballistic→diffusive crossover, jamming density ~2300-2400 cells/mm², h_p>h_v); no registrable parameters (constants in absent ESM), no direct oracle, but motivates an analytic PRW-MSD (Fürth/OU) oracle sourced elsewhere.
  - Single-cell mean crawl speed ~50-60 µm/hr (θ_F=180°), ~35-40 µm/hr (θ_F=240°) [Fig 3f]
  - MSD log-slope crossover: ~2 (ballistic) short-lag → ~1 (diffusive) long-lag (τ=400 min) [Fig 3b,d]
  - Critical PRW→CIL/jamming density ~2300-2400 cells/mm² (fibroblast, rigid substrate, no cadherin) [Fig 5c]
  - Nematic order: h_p asymptote ~0.8 vs h_v ~0.4-0.6, always h_p>h_v, saturating after ~20 h [Fig 7]
  - Final MSD slope biphasic in lamellipodium duration T_F (min at intermediate T_F) [Fig 3d]
  - Multi-cell anomalous long-tailed speed distribution up to ~150 µm/hr, mean drops with density [Fig 5e]
  · n_new_gates=4
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p65 [2] Taeyoon Kim · research — Computational Analysis of a Cross-linked
- headline: Ancestral Taeyoon-Kim BD crosslinked-actin sim — FF-engine blueprint yielding 3 sim-vs-sim cross-check morphology observables (crosslink angle, pore ∝ C_A^{−1/3}, L_c periodicity) but no truth contracts; params are lit-restatements or deliberately non-physiological.
  - Crosslink angle: bundler (α-actinin/fascin-type) −2.39±3.98° vs orthogonal (filamin-type) 87.5±10.5° — sim output, archetype discriminator
  - Mean pore size ∝ C_A^{−1/3} (∝ inter-filament spacing r₀), slope ~2.25, range ~4–10 σ_A — analytic geometric scaling
  - Inter-crosslink spacing L_c bimodal at ~5 σ_A & ~10 σ_A (half/full 74-nm helical turn) for bundlers vs monotonic for orthogonal
  - F-actin persistence length L_p≈20 µm calibration target (k̃_b,f=3000; matches Gittes/Isambert) — cite primary, not Kim
  - Bell force-dependent ACP off-rate k_ACP,−=k⁰·exp(γF/kT) confirms existing Bell-Evans oracle in a new (crosslinker) domain
  - NON-physiological: k_s,f softened ~25× (strain 0.05 vs 0.002@60pN) and k_A,e=1.41e9 M⁻¹s⁻¹ accelerated (real ~1.2e7) — do NOT register
  - Bundle effective diameter σ_b=√(N_f,b)·σ_A, measured 1.5–2.3 σ_A — analytic identity
  - Polymerization/length scaling collapses t̃_p~Φ^{−5/6}Da_n^{1/2}, ⟨L_f⟩~(Da_n Φ^{1/3})^{1/2} — sim kinetics, not FF-actionable
  · n_new_gates=3
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p55 [2] Makito Miyazaki · experimental — Directional Bleb Formation in Spherical Cells under Temperature Gradie
- headline: Experimental blebbing paper: no new params/oracles, but stages 3 candidate overlay-only gates for a directional-bleb / cortex-membrane-cohesion contract (DCM + FF) currently uncovered by the KB; confirms Young-Laplace γ·2/R-vs-ΔP balance already an oracle.
  - Directional bleb extension length 8.1 ± 3.0 µm (trypsin-rounded interphase, ΔT=12 °C, R≈10 µm, n=43; metaphase 7.5 ± 5.3 µm)
  - Bleb angle relative to asymmetry axis 8.4 ± 5.2° (n=26) — strongly co-axial/directional
  - Sphericity required: flat adherent interphase blebs only 29% (4/14) at 0.5 ± 0.8 µm vs 88–98% for rounded cells
  - Cortex F-actin & myosin(MRLC) density ↓ at high-activity pole, ↑ at opposite pole, in-phase with membrane, preceding detachment/rupture
  - Actomyosin is the effector: LatB/CytoD/blebbistatin/Y-27632 suppress; ML-7(MLCK) does not; blebbistatin blue-light restores length 0.8→3.4→6.2 µm
  - Convection-flow control gives 0.0 ± 0.0 µm bleb, then 6.3 ± 2.7 µm on heating — rules out flow artefact
  - Bleb onset latency ~3 s (timescale-fragile, not proposed as a gate)
  · n_new_gates=3
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p47 [2] Taeyoon Kim · research — Role of actin filaments and cis binding in
- headline: Simulation (2D force-free cadherin-clustering BD model) — all headline numbers are model outputs (comparison-only, never a gate/fit target); orthogonal to our junction-force gates, adds no analytic oracle, but is a strong mechanism/param reference for the ongoing DCM cadherin cluster+maturation redesign.
  - K_D^cis ~500 µm⁻² (range ~10¹–>10³) — MODEL OUTPUT, comparison-only (Fig 2B)
  - cluster lifetime τ ∝ [cad]² quadratic scaling — emergent, NOT an analytic oracle (Fig 2F)
  - max cadherin cluster ~600–700 monomers at [cad]=1600–2000, weak cis — sim output (Fig 2C)
  - cortical-actin corralling → sub-diffusive MSD (slope<1) — best native-FF target (Fig 3C)
  - cad-actin off-rates 0.1 s⁻¹ (α-catenin) / 10 s⁻¹ (vinculin) — CITED secondary [4,8]
  - cadherin-actin Mander's colocalization 0.076 — CITED expt [31], overlay only
  - cluster nearest-neighbour spacing ~31–45 nm across φ (Fig 5H); cadherin radius 2.5 nm
  - k_ass^cad/actin 10¹–10⁵ s⁻¹ flagged ARBITRARY/no-ref — do NOT register
  · n_new_gates=3
  · n_tighten=0
  · n_params=6
  · paper_kind=simulation

## p45 [2] Taeyoon Kim · research — Characterization, Enrichment, and Computational
- headline: Taeyoon-Kim-family agent-based CLAN paper: gives portable FF/DCM test-rigs (virial-stress network rheology; LINC nucleus compression) plus two qualitative simulation cross-checks (tension-dominant anisotropy 2.4:1:0.5; LINC nucleus stiffening x2.7), but all key numbers are simulation outputs or off-target glaucoma-TM measurements and the force constants live in the unavailable supplement — no new hard gates, no new params, no oracle.
  - Actin-network stiffness anisotropy tension 2.4e4 : shear 1.0e4 : compression 0.5e4 Pa (|epsilon|=0.05, sim, Fig6B)
  - Nucleus resistant stress ladder 2300 (no CLAN) / 3200 (CLAN no-link) / 6200 Pa (CLAN+LINC) at epsilon=-0.5 (sim, Fig6D) -> LINC coupling x2.7
  - AFM cell modulus CLAN+ ~3.5-4 kPa vs CLAN- ~0.3 kPa, GTM3L, Hertz nu=0.5 (off-target cell)
  - ACP crosslinker molar ratio R_ACP=0.1; actin C_A=100 uM (physiological, ref37)
  - LINC-linked F-actin fraction 20%; model nucleus diameter 5 um
  - CLAN incidence ~47% on glass vs ~5% on 2.36 kPa GelMA hydrogel (mechanosensing)
  - CLAN persistence 58.9+/-61.7 min (PEG withdrawn) vs 717.1+/-453.9 min (maintained)
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p40 [2] Taeyoon Kim · research — Durability of Aligned Microtubules Dependent
- headline: MT–kinesin collective-gliding paper adds microtubule/kinesin parameters (Lp 0.65–3.54 mm, v 601 nm/s, stall 5–7 pN, depletion 0.11 pN) for a future FF MT compartment plus a candidate cantilever-mode bending oracle, but yields no single-cell validation gate (off the MCF7 roadmap).
  - MT persistence length Lp = 0.65/1.93/3.54 mm (softer-GTP/stiffer-GTP/GMPCPP MTs, taxol, 25C) — distinct from collagen ell_p=17um
  - Kinesin-1 gliding velocity v = 601 ± 0.04 nm/s (all groups equal)
  - Single-kinesin-1 stall force = 5–7 pN (vs myosin F_stall=2 pN in KB)
  - Collision dispersion η_d = 16.3 ± 5.4°, Lp- & MC-independent (N=1200, p>0.05)
  - MC 0.3 wt% depletion force = 0.11 pN
  - Bundle length L_B = 48.7/164.3/141.5 µm (soft/stiff/GMPCPP, 50 min, ρ=5.5 µm⁻²); nematic S(50min)=0.54/0.62/0.80
  - Sim (comparison-only): stiffer N_B>130 & >200s vs softer N_B<45 & <20s; formation equal, disintegration differs
  · n_new_gates=0
  · n_tighten=0
  · n_params=6
  · paper_kind=experimental+simulation

## p38 [2] Makito Miyazaki · research — Controlling Physical and Biochemical Parameters of Actin
- headline: Experimental in-vitro nucleation-patterning paper: FF-relevant only — supplies one registerable cross-check param (NPF areal density 10-30k/um^2) and one qualitative architecture sanity gate (Arp2/3-branched-perpendicular vs formin+fascin-bundled-parallel); no scalar bands at our single-cell scale, no oracle, no fit targets.
  - NPF (His-VVCA/N-WASP) surface density 10,000-30,000 molecules/um^2 (0.66%/2% Ni-lipid), matching in-cell 6,500-25,000/um^2
  - Architecture dichotomy: Arp2/3+VVCA -> branched pillar perpendicular to membrane; mDia1DeltaN+fascin -> bundled vortex parallel to contour
  - Barbed-end polarity: new actin incorporates at pillar base (membrane) only (two-color pulse-chase)
  - Pillar length ~25-40 um (10 um micropattern); length ns low->med, high>med (p<0.001)
  - Reconstitution mix: actin 5 uM, Arp2/3 50 nM, profilin 15 uM, CP 25 nM, fascin 50 nM
  - NPF density vs Ni-lipid calibration FI=157.7x-2.32 (R^2=0.88) - instrument-specific, not a physical law
  - Area per POPC lipid 0.64 nm^2; bilayer height 4 nm
  · n_new_gates=1
  · n_tighten=0
  · n_params=1
  · paper_kind=experimental

## p32 [2] Taeyoon Kim · review — Weldon School of Biomedical Engineering,
- headline: Taeyoon-Kim-lab review (FF's ancestor program) — all numbers secondary/overlay; corroborates existing H9/H10/strain-stiffening gates and stages 5 mechanistic cross-check candidates (best: cortex G″-minimum at crosslinker off-rate), but supplies no MCF7 datum and tightens no band.
  - F-actin persistence length L_p ~9 µm, d 5-9 nm (semiflexible, ≈contour length), →[27]/[26]
  - Epithelial-dome superelastic surface-tension plateau σ ~1-3 mN/m at areal strain up to 1000% (Latorre MDCK, →[8])
  - Fascin-actin G″ has a local minimum AT the crosslinker unbinding rate; low/high-ω slopes ~0.06/0.5 (→[96]) — direct FF Bell-Evans off-rate cross-check
  - Poroelastic→viscoelastic crossover at t~0.05 s; water-flow relaxation window 0.1-10 s (→[101]/[98])
  - Nuclear chromatin↔lamin-A load-bearing crossover at extension ~3 µm; lamin network ~2× area reserve (→[104]/[105])
  - Crosslinked-actin differential modulus K′ ~1→100 Pa rising above critical prestress, myosin:actin 0.02/0.005/0.001 (→[115])
  - Poroelastic diffusion law D_p = E_solid·ξ²/η_fluid — candidate analytic oracle for H10 cytoplasm (Dp 40-60 µm²/s band)
  - Isolated-nucleus (HEK-293) AFM peak force rises with loading rate 1→16 µm/s — poroelastic rate-dependence (→[159])
  · n_new_gates=5
  · n_tighten=0
  · n_params=5
  · paper_kind=review

## p30 [2] Makito Miyazaki · research — PHYSICAL REVIEW RESEARCH 5, 013208 (2023)
- headline: Confined Xenopus-extract actomyosin flow↔wave↔rotation pattern-selection paper; genuinely new observable class (actin-wave period, cortical flow, Pe control) but at 200–300 µm droplet scale far from MCF-7 single cell — no current gate confirmed/tightened/challenged, no registrable param, no oracle; stage 2 pattern gates as overlay-only PI-gate candidates for a future DCM active-cortex module.
  - Periodic actin-wave period T ≈ 1.5 min (1/T ~0.8–1.0 min⁻¹), set by polymerization not contractility, high-CalA Xenopus extract droplet (Fig.1h)
  - Steady-flow→wave transition threshold: waves emerge above ~0.3 pM Calyculin A (contractility axis) (Fig.1h)
  - Steady inward actomyosin flow speed ~0.7 µm/s (PIV, control, droplet midplane) — droplet-scale, not single-cell (Fig.1c)
  - Two control knobs = myosin contractility σ_act & F-actin polymerization rate kp; Flow=lower-right / Wave=upper-left of (σ_act,kp) plane (Fig.2)
  - Simulation-only: steady flow at Pe=38, periodic waves at Pe=115, Pe=(ζµ)₀/(Dγ) (comparison, never truth) (Fig.3)
  - Simulation-only: waves require αp=kbulk_p/ksurf_p < 0.5; rotational wave stable only on the flow↔wave phase boundary, τd-independent (Fig.3h,4d)
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p29 [2] Taeyoon Kim · research — Mechanics and Morphology of Proliferating Cell Collectives with Self-I
- headline: Analytic proliferation-mechanics theory (bacterial/dimensionless): clean closed-form Darcy-growth oracles (p(r)=(R²−r²)/4, R∼(2/λ)^{1/2}t, Λ∝λ^{−1/2}) but only usable if DCM gains growth+division; nothing registerable now, no MCF7 params, no existing gate confirmed/challenged.
  - Colony radius crossover: R(t)=R₀e^{t/2} early → linear R∼(2/λ)^{1/2}t asymptote (Eq.8, nondim, discrete+continuum agree)
  - Darcy-growth pressure profile (λ=0): p(r)=(R²−r²)/4, velocity u(r)=r/2 (Eq.5, analytic closed form)
  - Concentric size-ring wavelength ∝ λ^{−1/2}, decreasing with radial distance, oscillation period ≈ln2 (Fig.4c)
  - Stress-inhibited growth law: ℓ̇/ℓ=(1/τ)e^{−λσ}, g=(1/τ)e^{−λp} (Eq.7/2b) — mechanistic pressure→growth coupling
  - Max radially-averaged stress at colony center, ~10× drop λ=0→10⁻² (peak ~6000–7000 nondim at λ=0, Fig.3a,c)
  - E. coli estimate (from ref [58], not measured): τ≈54 min, ℓ₀≈1 µm, dR/dt≈1.3×10⁻² µm/s → λ≈1.1×10⁻³ (dimensionless)
  · n_new_gates=4
  · n_tighten=0
  · n_params=0
  · paper_kind=analytic

## p28 [2] Makito Miyazaki · review — Molecular Crystals and Liquid Crystals
- headline: Review (all values re-cited, overlay-only); ~80% out of scope (muscle SPOC + mitotic spindle). In-scope residue is narrow: actin L_p 5-15µm (FF bending cross-check, likely already known) + confined-actin-ring self-organization <20µm as a DCM cytokinesis candidate. No new oracle, no gate/param conflicts.
  - Actin persistence length L_p = 5-15 µm (fig 5-18), F-actin (Gittes 1993 [65]) — FF bending anchor, distinct from collagen ell_p=17µm
  - Contractile actin ring self-organizes when confinement diameter < 20 µm (≈ L_p) — DCM cytokinesis cross-check (Miyazaki 2015 [5])
  - Ring constriction rate ∝ initial ring diameter, HMM 25 µM (reconstituted, [5])
  - Mitotic-spindle rigidity 10-100 Pa, anisotropy ~2×, Poisson ratio ~2.0 — OUT OF SCOPE (no spindle model)
  - Spindle viscoelastic crossover: elastic <1 s / viscous ~10 s / active-elastic ~min; |G*| 0.1→1 nN/µm — out of scope
  - Cardiac SL 2.00→1.72 µm (~14%), ΔLVP=287.3·ΔSL (R=0.60) in-vivo mouse — out of scope (sarcomere)
  - Giant uni-lamellar liposome yield >90% (inverted emulsion) — DCM membrane program method note
  · n_new_gates=2
  · n_tighten=0
  · n_params=1
  · paper_kind=review

## p25 [2] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: Kim-lab agent-based motility-assay simulation — comparison-only for FF's filament+motor collective regime (network/flock/band/ring, C<60µM band-collapse); no analytic oracle, no sourced params (all in unavailable SI Table S1), zero gate overlap with existing KB.
  - Emergent structure sequence network→flock→band→ring as κr↑ / κb↓ (4 classes), motility assay, 100 s (Figs 1–6)
  - Thick coherent bands vanish for actin C < 60 µM (collectivity threshold), reference κr/κb (Fig 2)
  - Velocity correlation ⟨|cos α|⟩ ~0.5–0.95, endpoint pairs <0.5 µm, avg 50–100 s (Figs 1c–3c)
  - Velocity-field curl ~0.018–0.033 s⁻¹, peaks where many rings (Figs 1d–6d)
  - Filament curvature ordering ring>band>flock>network (Fig S7f)
  - Ring size set by rigidity κb & repulsion κr, NOT by C or ⟨L⟩; ring count ∝ ⟨L⟩ (Figs S8–S12)
  - Implicit motors form NO rings + broader low-speed fraction vs explicit motors (Fig S6)
  - l_p = κb·r0/(k_BT); Δt=1.15e-5 s; effective myosin 96 µM (8 heads/arm, CM=12 µM); ⟨L⟩ref 2.5 µm
  · n_new_gates=3
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p18 [2] ? · ? — Durotactic Migration Driven by Anisotropic Matrix
- headline: Peer agent-based durotaxis model (Taeyoon-Kim lineage, same Mikado-ECM + SimuCell3D physics as FF/DCM): all numbers are simulation output → comparison-only, no oracle, no quantitative bands; contributes qualitative sign/trend/null cross-checks (durotactic-index positivity + frozen-matrix null; adhesion-driven F_r/F_θ>1) for the FF⊗ECM durotaxis and DCM aggregate-escape tracks.
  - Durotactic index max ~0.6-0.7 (↑ with strain ε, ↓ with porosity p); ≈0 in frozen-matrix vᵢ=0 control (sim, Fig.3C/E)
  - Radial/circumferential force ratio F_r/F_θ ≈ 2 (2-3 with cell-cell adhesion) in contracting cluster (sim, Fig.4D)
  - Cells with radial-bias>0.8: 28% single-cell / 73% multi-no-adh / 82% +adhesion / 83% +1.5×contractility — non-additive (sim, Fig.4F)
  - Stiffness anisotropy E∥/E⊥ peak ~1.05-1.07 at ε≈1-1.5, p=0, κ_b,M=0.001 nN·µm (sim, Fig.2B)
  - Matrix contraction strain ~0.5 sustained to ~550 h with adhesion vs decay after 100-200 h without (sim, Fig.4C)
  - Cell-cell adhesion κ_a,int = 1e-5 N/m (SimuCell3D tent form, Eq.8) — coarser than our KU-4.2 catch-bond
  - Fiber buckling asymmetry: κ_s,M ÷10 in compression (tension-stiff/compression-soft, Eq.13)
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p16 [2] Taeyoon Kim · research — Computational Modeling of Motor-Driven
- headline: Simulation sibling of FF (Kim cylindrical-element lineage) on axonal MT bundles — no existing gate touched; net-new value is one analytic Euler buckling oracle (F_cr∝κ_b/L_uns²) plus preprint-flagged dynein/MT params, all gated on a future FF MT module.
  - Euler critical buckling load F_cr ∝ κ_b/L_uns² (Eq.10, analytic) — candidate FF oracle FF currently lacks
  - Dynein unloaded walking speed ~1.1 µm/s → ~0 at stall; k_u,M rising to ~1 s⁻¹ under load (Monzon 2018, Fig.1F)
  - Crosslinker Bell slip-bond k_u ~1 s⁻¹ (0 pN) → ~4.5 s⁻¹ (60 pN), Tau-like (Bell 1978, Fig.1E)
  - Buckling curvature threshold 0.2 rad/µm per-node (Memet 2018, Eq.9)
  - MT persistence length l_p ≈ 1 mm (Isambert), MT diameter 25 nm — distinct from collagen ℓ_p=17 µm, do not cross-anchor
  - Sim outputs (comparison-only): ε_max ~0.08 at ~100 s no-XL; biphasic crosslink optimum R_X≈0.25 µm⁻¹; motor force ~5–7 pN
  · n_new_gates=3
  · n_tighten=0
  · n_params=5
  · paper_kind=simulation

## p04 [2] Taeyoon Kim · research — Balance between Force Generation and Relaxation
- headline: Sibling FF-class actomyosin simulation: mechanism-confirmatory (PCM/Hill motor, Bell ACP, treadmilling) with 2 candidate FF gates (pulsed-cluster 2–4 µm/40–100 s overlay; turnover-homogenizes sign gate) — no new oracle, no sourced params.
  - Ensemble motor speed ⟨v_M⟩ ~0.6–1.8 µm/s (R_M=0.04, network-emergent; simulation → soft FF cross-check vs Hill/PCM)
  - In-vivo pulsed-cluster size ~2–4 µm dia (Nishikawa ref 36; OVERLAY-only target)
  - In-vivo pulsed-cluster duration ~40–100 s (Nishikawa ref 36; OVERLAY-only)
  - Emergent weak-pulse cluster lifetime <50 s, size 0.5–3% of actin (balance regime, no severing; simulation)
  - Turnover monotonically homogenizes network (k_t,A 15→240 s⁻¹ ⇒ heterogeneity Q_A ↓; qualitative FF cross-check)
  - Contraction maximal at intermediate connectivity (R_ACP 0.0025–0.16; ⟨L_f⟩ 0.7–5 µm; qualitative)
  - Angle-dependent severing law Eq.1 k_s,A=k⁰·exp(θ/λ), >50% severing inside clusters (candidate FF mechanism, not oracle)
  · n_new_gates=2
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p77 [1.5] Makito Miyazaki · research — Spatial confinement of active microtubule networks
- headline: TAG-validation lens — p77  (Suzuki, Miyazaki et al. 2017, confined MT rotational flow, PNAS 114:2922)

## p52 [1.5] Makito Miyazaki · research — Microscopic Temperature Control Reveals Cooperative Regulation
- headline: Experimental HMM-actin gliding-assay paper: supplies motor sliding-speed + T-dependence (~7-15 um/s) and reusable actin constants (36 nm half-pitch, 1:5 stoichiometry), but its drebrin-E cooperativity (Hill 1->4, Kd 6/32 nM) is a neuronal regulator outside MCF7 engine scope; no gate overlap, no band challenged, 1 dormant FF cross-check candidate + 4 staged params.
  - HMM-actin sliding velocity ~7-8 um/s at 28-30 C, rising to ~15 um/s at 36.5-37.5 C (in-vitro gliding assay, Fig 4B)
  - Drebrin-E inhibition Hill coefficient n ~1 below 34 C jumping to 4.0 (3.3 excl. Qdot) at 36.5-37.5 C (Fig 4D)
  - IC50 of drebrin-E inhibition ~2-3 nM (low T) rising to >30 nM (37 C)
  - Kd drebrin-E-actin = 6 nM at 25 C / 32 nM at 37 C (luminescence densitometry; supersedes old 120 nM)
  - Actin helical half-pitch 36 nm bare -> 40 nm with drebrin (AFM, ref 22)
  - Drebrin:actin stoichiometry 1:5 subunits
  - Inhibition-recovery time constant 5.4 s (30 nM drebrin, heat-off, Fig S3)
  - Drebrin stays bound during sliding: Qdot F.I. unchanged, Friedman p=0.344 NS (n=15)
  · n_new_gates=1
  · n_tighten=0
  · n_params=4
  · paper_kind=experimental

## p51 [1.5] Taeyoon Kim · research — Mechanical Model for Durotactic Cell Migration
- headline: Kim-lab durotaxis simulation — durotaxis emerges from pure mechanics (DI biphasic in κ^L, falls with viscosity, transition at κ^H≈0.15 nN/µm); a behavioral FF cross-check, not a truth source; 3 candidate FF gates, 0 params.
  - Durotaxis emerges purely mechanically (no decision rule), DI ~1 vs 0.5 unbiased on a 100:1 stiffness step
  - DI is biphasic in soft-region stiffness κ^L and decreases monotonically with substrate viscosity ξ
  - DI high-level transition at κ^H_s,M ≈ 0.15 nN/µm (elasticity overtakes viscosity)
  - Migration speed 5–35 µm/hr, set by LOCAL stiffness independent of orientation (speed is not the durotaxis driver)
  - Distance-before-reorientation ∝ 1/gradient (~200 µm at low → ~50 µm at 200 pN/µm²)
  - Substrate chain forces up to ~0.5 nN in asymmetric front>rear dipole; network is linear-then-strain-stiffening
  · n_new_gates=3
  · n_tighten=0
  · n_params=0
  · paper_kind=simulation

## p09 [1.5] Taeyoon Kim · research — Reciprocal folding dynamics in cellular networks at the stroma-basemen
- headline: Tissue-scale continuum-FEM + in-vitro folding paper: 0 new single-cell gates; only a soft E_Str→G0 cross-check plus SourceEvidence-grade overlay moduli — weak for our single-cell KB.
  - E_Str = 354 Pa (fibroblast + 2 mg/ml collagen-I stroma, tensile) → G≈118 Pa, CONFIRMS VG-U1-G0-linear [15,200] Pa band
  - E_BM = 756 kPa (engineered nm-thick collagen-I basement membrane, dog-bone tensile)
  - E_Cell = 8.6 kPa (3T3 fibroblast AFM, cited) — non-MCF7 overlay only
  - P_Cell = 200 kPa imposed FEM boundary stress — model input, NOT a single-cell traction (do not register)
  - +1.91 µm evagination at measured stiffnesses — FEM simulation output, comparison-only (never truth)
  - Cell–cell coupling threshold ~90 µm (72.3/91.8/226.7 µm); tissue-scale, no single-cell analog
  - Fold angles: joining 16.15±7.94°, branching 100.20±10.26° — multicell network, not registerable
  - Cell–fold / collagen–pseudopod–fold co-alignment 93.5–97.8% — mechanistic narrative for FF only
  · n_new_gates=0
  · n_tighten=0
  · n_params=3
  · paper_kind=experimental

## p05 [1.5] Taeyoon Kim · research — Mobility of Molecular Motors Regulates Contractile
- headline: Kim-lab agent-based (simulation) actomyosin paper: strong FF cross-check partner but low new-contract yield — supplies motor-mobility contractility trends (1-arm biphasic vs 2-arm monotonic) + one wet-lab calibration velocity (140 nm/s NMII); its Q_A/force outputs are comparison-only, never truth.
  - Unloaded NMII walking velocity 140 nm/s @ reference k20 (experimental anchor, refs 39/40); ~2 um/s at highest k20
  - Avg tensile force per motor arm: one-arm up to ~13 pN, two-arm up to ~8 pN (=> ~1-3 pN/head over N_h=4-8) — consistent with F_stall_motor=2 pN
  - F-actin average speed ~0.02-0.08 um/s (one-arm), ~0.02-0.06 um/s (two-arm)
  - Heterogeneity Q_A vs motor mobility: BIPHASIC (peak at intermediate mobility) for one-arm motors; MONOTONIC-increasing for two-arm motors
  - Motor arm binding rate 40*N_h /s; N_h=8 (one-arm)/4 (two-arm); stall force ∝ 1/k20; F_max=N_h*stall*#arms
  - Reference network: <L_f>=1.5 um, actin L_p=9 um, R_M=0.8, R_ACP=0.1, k20=20/s, zeta*_M=8.10e-8 kg/s, 10x10x0.1 um quasi-2D slab
  · n_new_gates=3
  · n_tighten=1
  · n_params=6
  · paper_kind=simulation

## p68 [1] Taeyoon Kim · research — sharing, adaptation, distribution and reproduction in any medium or fo
- headline: Plant (Arabidopsis) agent-based cortical-MT stress-ordering sim — orthogonal to our MCF7 actin/ECM contracts; 0 new gates, only donor MT-DI params staged for a future (mammalian) microtubule module.
  - Plant cortical MT plus-end DI rates (Shaw2003 via cite): v_grow=3.69, v_shrink+=5.80, v_shrink-=0.53 µm/min
  - State-transition freqs 0.47/0.97/0.51/0.24/0.87/0.23 min⁻¹; nucleation 10 µm⁻²·min⁻¹
  - MT-MT collision: 40° zipper/catastrophe threshold; P(catastrophe)=0.2-0.8; spacing 25-50 nm
  - Single cortical MT length 2-4 µm (lit), model 2.0-2.9 µm; lifetime 2.2-3.0 min
  - Nematic order S_p ~0.4-0.5 isotropic vs ~0.8-0.85 anisotropic; alignment threshold 1.5:1
  - Array reorientation ~50 min (10 min-2 h lit) under 90° stress flip; wall stress 3-15 MPa
  - MT dynamic-instability tension scale ~10 pN (Akiyoshi2010 via cite)
  · n_new_gates=0
  · n_tighten=0
  · n_params=6
  · paper_kind=simulation

## p49 [1] Makito Miyazaki · research — Kinetic Scheme of Myosin Phosphorylation by ZIP Kinase
- headline: Upstream ZIPK→MRLC phosphorylation kinetics (kcat 0.40–0.60 s⁻¹, Km 4.2 µM, KI 0.30 µM); clean parameter/oracle source for a future myosin-activation layer but touches no current force/tension gate.
  - ZIPK-FL Km(MRLC) = 4.2 µM, both Ser19/Thr18 sites equal (5 nM ZIPK, 25°C)
  - kcat Ser19 = 0.60 s⁻¹ vs Thr18 = 0.40 s⁻¹ (ratio 1.5×), ZIPK-FL
  - Sequential Ser19-first order: monophospho Ser19:Thr18 ≈ 35:1 (LC-MS/MS)
  - ZIPK–MHC dissociation KI = 0.30 µM (~10× tighter than Km) → SMM sequestration
  - Native SMM activates 3.5× slower than isolated MRLC: Total-P 31.5% vs 111%
  - ZIPK-ΔC (1–291): Km 2.5 µM, kcat 0.95/0.63 s⁻¹, weaker KI 0.65 µM
  - Activation timescale ~1–10 s from kcat 0.40–0.60 s⁻¹ at physiological setpoint
  · n_new_gates=3
  · n_tighten=0
  · n_params=7
  · paper_kind=experimental

## p33 [1] Makito Miyazaki · experimental — Accurate polarity control and parallel alignment of actin
- headline: Experimental actin-alignment nanofab paper; validation value minimal — actin l_p 15-20 µm and 6 nm diameter are secondary cited constants, myosin-V 0.62 µm/s is wrong-isoform; no new gates, no oracle.
  - Actin persistence length l_p 15-20 µm (cited refs 16-19, single filament)
  - Actin filament diameter 6 nm (cited ref 15)
  - Myosin-V bead velocity 0.62 ± 0.10 µm/s (n=29, in vitro isopolar array)
  - Isopolar-array alignment σ = 10.8° Gaussian (n=134) — device metric, no cell gate
  - Polarity orientation accuracy 99.6% (n=264) — device metric
  - Microtubule l_p 5000-6000 µm / diameter 25 nm (comparison track)
  - Biotin-streptavidin K_d ~1e-14 M (cited ref 32)
  · n_new_gates=0
  · n_tighten=0
  · n_params=2
  · paper_kind=experimental

## p31 [1] Taeyoon Kim · research — PHYSICAL REVIEW RESEARCH 7, 013312 (2025)
- headline: Kim-lab agent-based motility-assay simulation (mobile-motor F-actin clustering + PCM motors); off our single-cell scale — only an analytic slender-body drag closed-form and PCM motor-kinetics anchors are FF cross-check-worthy; nothing confirms/tightens/challenges any existing gate band.
  - Reference motor diffusion D*_M = 5.11e-14 m2/s; mobility swept D̄_M = {0, 0.0625, 0.25, 1, 4} (model)
  - Cylindrical-segment drag closed-form ζ = 3πμ·r_c·(3 + 2 r0/r_c)/5 (Eq. A3, analytic — candidate FF drag oracle)
  - Motor binding rate k+,M = 40·N_H = 160 s^-1; binding-site spacing 7 nm; N_H=4 heads lumped/segment (model)
  - Actin-density onset for collective behavior ≈903 monomers/µm² (negligible) → strong clusters at 6020–12040 (model)
  - Filament-length crossover: clusters at ⟨L_f⟩≤2.58 µm vs closed rings at 7.43 µm (model)
  - Polar order P ~0.5–0.6 (mobile D̄_M=1,4) vs ~0.2–0.3 (fixed D̄_M=0); phase-sep M/M0 ~3 vs ~1 (model output)
  - Euler timestep Δt = 1.15e-5 s; motor rest length / actin-plane gap 13.5 nm (numerical/geometry)
  - Exp TIRF (qualitative, no fit target): short 0.5–1 µm + high conc → clusters; long 1–2 µm + low conc → bundles
  · n_new_gates=1
  · n_tighten=0
  · n_params=2
  · paper_kind=simulation

## p56 [0.5] Makito Miyazaki · method — Quantitative Analysis of the Lamellarity of Giant Liposomes Prepared b
- headline: Method/imaging paper: certifies inverted-emulsion GUVs are ~97% unilamellar (justifies single-bilayer membrane assumption) but supplies zero mechanistic force/rate/modulus data — 0 new gates, 0 params, 0 oracles; reference-only.
  - Inverted-emulsion GUVs 98.0% unilamellar (1 mM egg PC) vs hydration 49.8%; ~97.2% mean across all conditions
  - Fit focal-leakage depth sigma=0.913 um and pixel size delta=0.143 um (microscope-optics constants, Eq.1 fit; NOT physiological)
  - Unilamellar fraction stays >90% (92.1% @10 mM) across lipid conc 0.1-10 mM, +/-20% DOPE/DOPG, +/-50% cholesterol, and F-actin/BSA/egg-extract cargo
  - alpha-hemolysin stem ~5.2 nm ~= single bilayer thickness, pore diameter ~2.8 nm (both cited from Song 1996, not measured here)
  - alpha-HL actin-bundling proof: 94.3% bundle+ in lowest (unilamellar) band vs 8.7% in higher bands
  - Eq.1 = geometric-optics fluorescence-vs-radius curve (Gaussian out-of-focus kernel) — an imaging forward model, not a mechanobiology oracle
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=method

## p54 [0.5] Makito Miyazaki · experimental — on the number and size of liposomes formed
- headline: In-vitro inverted-emulsion liposome-formation methods study; all data are protocol-yield/count statistics with zero single-cell mechanical observables — no gate, param, or oracle candidates, methods-provenance only.
  - 10% DOPA raises liposome yield 13.7× vs egg-PC control (natural sedimentation, Fig. 2D)
  - 10% DOPG / DOPS yield 8.83× / 8.84× (natural sed.); centrifugal DOPG peaks ~10× sharply at 10% then collapses (Fig. 5B)
  - Yield vs centrifugation T: y=68.6·e^(−0.0855·T °C), ~2× rise 16→0 °C (empirical fit, Fig. 4B)
  - Yield vs incubation time: y=102−25.6·e^(−0.249·t min), y(0)≈75% of plateau (empirical fit, Fig. 4D)
  - Prob(liposome >10 µm) correlates positively with additive head MW; W/O droplet size distribution unchanged by composition (Fig. 3)
  - egg PC mean MW 770.123 Da; cell-sized liposomes counted only if spherical and >5 µm diameter (metadata)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p23 [0.5] Makito Miyazaki · method — Bio-protocol 16(8): e5656. DOI: 10.21769/BioProtoc.5656
- headline: Method/protocol paper (optogenetic Arp2/3 actin on a lipid bilayer) — no mechanistic parameters, no closed-form laws, no observables at our single-cell scale; zero gates/params/oracles, at most a PI-gated method SourceEvidence pointer.
  - iLID light-intensity control range 0.04–0.18 nW/µm² (apparatus calibration, 6% DGS-NTA(Ni) SLB)
  - SspB-VCA recruitment saturates ~0.18 nW/µm² (above → density uncontrollable)
  - Network density ∝ light intensity (monotonic, n=6) with transient depletion overshoot at ≥0.08 nW/µm² — qualitative trend from companion Nat Commun
  - Network thickness ∝ illumination duration (10/20/30 min), plateaus ~30–40 µm (n=12)
  - SspB membrane-recruitment delay ~20 s after light on; 500 ms/20 s stim duty cycle
  - Reconstitution recipe: actin 5 µM, Arp2/3 100 nM, profilin 15 µM, CP 25 nM, SspB-VCA 150 nM (setup values, not params)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=method

## p20 [0.5] ? · ? — PAPERS |  SEPTEMBER 01 2022
- headline: Physics-education outreach paper (elementary-student Brownian-motion workshop); only quantitative content is textbook MSD=4Dt + hand-clicked bead diffusion — zero mechanobiology validation data, no gates/params/oracles to register.
  - 2D Einstein diffusion law MSD = 4Dt (Eq.1, p.480) — canonical textbook analytic relation, no D evaluated
  - √MSD ∝ t^(1/2) diffusive scaling vs ballistic ∝t (Eq.2, Fig.4a) — qualitative demo
  - √MSD of 1-µm beads in water ≈6–7 µm at t≈20 s (Fig.4a) — hand-clicked, defocus-contaminated, not quantitative
  - Diffusion coefficient D: NOT numerically reported (symbol only)
  - Dice 2D 6-direction random walk √MSD ∝ N^(1/2) (Fig.8) — Monte-Carlo teaching demo
  - Survey: interesting 4.69±0.64, difficulty 3.27±1.09 (n=26, Fig.10) — pedagogy, out of scope
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=review

## p11 [0.5] Makito Miyazaki · research — RESEARCH ARTICLE |  APRIL 07 2011
- headline: Single-molecule hidden-trajectory inference METHOD (Onsager-Machlup DRP / Go-and-Back) on synthetic toy models — zero cell-mechanics validation targets; no gates, params, or oracles for DCM/FF. File-only in the Miyazaki methods corpus.
  - Compute-cost reduction >10^3-fold vs gradient descent (method benchmark, dimensionless)
  - Converged Onsager-Machlup action S = -3.04e4 (Go-and-Back) vs -2.91e4 (grad. descent), Model A, model energy units
  - BC-diminishing timescale tau_dim ~ 0.005 (fit) / ~0.001 analytic = gamma/(a1+k), model time units, <0.1% of trajectory affected
  - Initial-condition sensitivity <|x_hat_noise - x_hat|> = 1e-4 model length units (Model A, 1e4 pts)
  - Thermal SD of x = sqrt(2 k_BT dt/gamma) = 0.2 model length units (standard FDT identity, already covered)
  - Cited context only: probe temporal res up to 9.1 us, spatial res 0.1 nm, F1-ATPase 120deg step
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=method

## p08 [0.5] Makito Miyazaki · research — Chimeras of kinesin-6 and kinesin-14 reveal head-
- headline: Experimental fission-yeast spindle/kinesin-6-14 domain-swap paper; MT-motor observables sit entirely outside the actomyosin-cortex/adhesion/ECM validation surface — 0 new gates, 0 params, 0 oracles.
  - WT Klp9 anaphase-B spindle velocity 0.68 ± 0.07 µm/min (≈11 nm/s), S. pombe 25°C, n=20 (Fig 2C)
  - Chimera P1h:K9nt velocity 1.19 ± 0.24 µm/min, ~1.5× WT (Fig 2C); no-neck P1h:K9t 1.51 ± 0.32 µm/min (Fig 2H)
  - Final anaphase spindle length: WT 9.75 ± 1.11 µm vs P1h:K9nt 13.80 ± 1.61 µm (over-elongated → buckle/break) (Fig 2D)
  - Midzone FRAP recovery 47 ± 18 s (WT) vs 29 ± 6 s (P1h:K9nt); motor front speed 1.62 vs 2.39 µm/min (Fig 3B–C)
  - Aneuploidy: ~2% via spindle breakage (P1h:K9nt) and ~8–11% via MT-protrusion/misposition (Pkl1 chimeras) vs WT <0.1% (Fig 3E/4G)
  - ~25 microtubules per spindle pole, cited context (p423)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=experimental

## p06 [0.5] Taeyoon Kim · review — Contents lists available at ScienceDirect
- headline: Plant cell-wall FE morphogenesis review (Kim co-author) — all numbers are MPa-scale turgor/wall mechanics, 10^3–10^6× off our animal single-cell regime; zero gates/params/oracles, LOW relevance.
  - Turgor pressure 0.5–1 MPa (plant epidermal cells) — 10^4× above our ~40 Pa animal baseline; out-of-regime
  - Cell-wall stress ~10 MPa, scaling σ~Π·R/t (≈10× turgor); thin-wall pressure-vessel cousin of Young-Laplace
  - FE trichome-branch wall stress 5.0–12.6 MPa (Yanagisawa 2015 simulation, Fig 4c)
  - Max principal strain 0.023–0.039 per FE growth cycle (simulation)
  - New wall material 3.0–5.3 nm/cycle, base-biased (FE simulation)
  - Microtubule persistence length ~mm (generic cited value, not measured here)
  - Trichome branch angular spacing ~120° over ~14 h growth phase (plant morphology/timescale)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=review

## p01 [0.5] ? · ? — RESEARCH ARTICLE |  FEBRUARY 28 2011
- headline: Statistical-inference method paper (Bayesian/Onsager-Machlup/WKB) on a nondimensional 2-bead Langevin toy; zero dimensional cell-mechanics validation targets, no gate/param overlap, only a method-scoped closed form (Eq.20) as a possible future FF inference-tool self-check.
  - Estimator accuracy gain vs power-spectrum fitting: 1-2 orders of magnitude (up to x100), precise phase (dimensionless)
  - Sample-averaged variance jump at loss-of-precision critical curve: >4 orders of magnitude, abrupt (Gamma* large, h* small)
  - Systematic error scales as tau^-1 (log-log slope -1) in precise phase; ~independent of Delta_t
  - Analytic MAP hidden-trajectory error Eq.20: <(x-x_hat)^2> = (1/(beta*h))[2/kappa - 1/sqrt(kappa^2+g)], kappa=1+k/h, g=gamma/Gamma
  - H_infinity minimum lies exactly at true Pi* (ergodic, noiseless, tau->inf)
  - Noise-robustness divergence cutoffs sigma_noise = 5.62e-1, 1.78e-1, 3.16e-1 (nondim, SD of x=1)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=method

## p34 [0] ? · ? — Integr. Biol., 2015, 7, 1093--1108 | 1093
- headline: Review paper — no original data, equations, or laws; every number is a second-hand citation, so zero new gates/params/oracles. Confirmatory/roadmap value only (points to FF network-rheology primaries: Kim refs 4/47, Bidone ref 30). Priority 0.
  - Talin H1–H12 rod SMD: 300 pN applied, ΔL 3.2 nm ref → 6.2 nm at rupture (cited ref 2, sim, below our CG floor)
  - Filamin A domain-20 unfolding force >~4 pN, single-molecule optical tweezer (cited ref 19)
  - Actin+ACP Brownian-dynamics network stress ~10–40 Pa, oscillatory strain 1–6 ms (cited ref 4, read off Fig.1b axes)
  - In-vivo talin extension up to ~several hundred nm, myosin-II-driven live cell (cited ref 23)
  - Motor-driven cytoplasmic fluctuation processivity ~10 s, intracellular tracking (cited ref 132)
  - Chromatin decondensation force ~nN order at membrane, seconds (cited ref 135)
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=review

## p12 [0] Makito Miyazaki · review — Honmachi, Sakyo-ku, Kyoto 606-8501, Japan. ORCID iD: https://orcid.org
- headline: p12 is a numbers-free symposium editorial (review); no validation data, no gate/param/oracle candidates — register only as a low-priority provenance/context node.
  - No quantitative observables — 3-page symposium editorial with zero equations, parameters, or measurements
  - Only numerals are administrative (journal Vol.19/e190031, meeting dates, PRESTO grant IDs JPMJPR20ED/JPMJPR20E6)
  - No force laws, potentials, or closed-form relations anywhere in the text
  - No DCM/FF-relevant mechanistic data despite Miyazaki cytoskeleton-lab byline
  - Qualitative thesis only: 'creating and manipulating' supramolecular assemblies to reveal design principles
  · n_new_gates=0
  · n_tighten=0
  · n_params=0
  · paper_kind=review