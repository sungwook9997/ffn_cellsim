# DCM digest (all 77 papers, sorted by score)


## p46 [4.5] Taeyoon Kim · research — The nature of cell division forces in epithelial
- headline: A published T. Kim triangulated-surface-mesh DCM shows cytokinetic-ring constriction plus enclosed-volume conservation yields ~60% axial elongation from only ~5% volume loss — a drop-in division module, equations, and validation oracle for our DCM.
  - Cytokinetic-ring division module: prescribed inward constriction of a 1-µm equatorial node band (constant v_c, freeze at 0.5 µm from axis) — no ring line-tension magic number
  - Enclosed-volume constraint (U_v) as the elongation engine: validate/tune κ_vol against the 5%-vol-loss→~60%-axial-elongation coupling (40% to abolish)
  - Soft/asymmetric surface-area law with membrane reservoir (κ_a 1e-4…1e-6 N/m) permitting +22% area growth — enables division, blebs, membrane-cortex coupling
  - Node-area-proportional drag ζ_i=λA_i (grid-invariant, derivable from η) supersedes uniform per-node γ
  - Dead-band extensional edge spring (slack window 0.05–1.25 µm) for remesh-tolerant neck stretch
  - Division is neighbor-transmitted: build the scenario in the N≈400 aggregate (single isolated cell gives no elongation traction)
  - Spindle NOT required for elongation: ring pinch + volume incompressibility alone suffice — big scope reduction for the division module
  · n_parameters=14
  · n_equations=8
  · n_validation_targets=16

## p22 [4.5] Taeyoon Kim · research — Cellular Pushing Forces during Mitosis Drive Mitotic
- headline: A validated fine-grained division protocol (ring + spindle forces under volume conservation) that fills DCM's biggest open gap — cell division — with concrete forces (100 pN/side) and elongation targets (20% ring-only, 40% ring+spindle).
  - Division module: cytokinetic-ring inward force + spindle outward polar force under volume conservation (fills DCM's open division/cytokinetic-ring gap)
  - Ring as a mechanistic line-tension law (T_ring·curvature) replacing the paper's constant 100pN — self-intensifying furrow
  - Volume-conservation-driven elongation: equatorial pinch at constant V => ~20% axial elongation (direct test of DCM turgor/enclosed-volume)
  - ECM confinement via node-face excluded-volume penalty (Eq 8 == Ericson) — divide against real collagen G'=10/100/400 Pa
  - Spindle force routed through astral-MT tips on DCM microtubule compartment (mechanistic upgrade over raw membrane force)
  - Spindle-ablation additivity test: F_Is off => partial matrix-deformation relaxation
  · n_parameters=13
  · n_equations=6
  · n_validation_targets=11

## p59 [4] Makito Miyazaki · research — Cell-sized spherical conﬁnement induces the
- headline: Complete mechanistic blueprint for DCM's open cytokinetic-ring/division module: emergent equatorial furrow, constriction-rate proportional to perimeter law, constant-volume constriction, and duty-ratio-squared effective-tension.
  - Cytokinetic-ring/division module: emergent equatorial furrow placement (energy-min great circle, R_ring/R_droplet=0.96) + line-tension constriction
  - Constriction law dP/dt = -k_c(P-P*), k_c=0.39 um/min per um, P*=4.5um set by filament length
  - Constant-volume constriction (V=piR*pid^2/4): remesh/mass-conservation on the furrow band
  - Duty-ratio-squared effective-tension law [myosin]_eff=[bulk]*(duty)^2 for cortex gamma generation
  - Cell-sized spherical confinement BC + criterion R_cell<=L_p (persistence length) for confinement response
  - Two-plate compression -> furrow reorients plate-parallel (perpendicular to compression)
  - Motor-escalation life cycle (weak->strong->complete R_ring->0, 44%) to drive complete division
  · n_parameters=13
  · n_equations=6
  · n_validation_targets=11

## p14 [4] ? · ? — Reconstitution of actomyosin networks in cell-sized
- headline: Cell-sized deformable-membrane bleb model hands DCM a portable surface-mesh potential set, a Bell-kinetics membrane-cortex detachment law with closed-form tension oracle, and a bottom-up ~0.2 mN/m cortex-tension anchor.
  - Membrane-cortex coupling as explicit Bell-kinetics breakable linkers (Eq.5) → mechanistic bleb / decoupling module
  - Two-mechanism bleb selection (Detachment vs Rupture) set by orthogonal R_C (coupling) and R_X (connectivity) dimensionless ratios — derivable, non-magic config axes
  - Closed-form detachment tension oracle γ_C^{D*} ≃ 0.35·(k_BT·ρ_H·R/r₀)·ln(ρ_A/K_m) ∝ R_C (Eqs 7-8) as a bleb validation gate
  - Directly-portable membrane-mesh potential set: per-triangle area (S9), enclosed-volume (S8), dihedral bending (S6), steepened node-face self-contact (S10)
  - Young-Laplace split-tension closure ΔP=2(γ_m+γ_c)/R with γ_c≫γ_m at onset (Eqs 2-3) — same relation as DCM γ=ΔP·R/2, usable as a gate
  - Cortex Rupture via local tension-threshold severing cascade (F_sev=300 pN) as mechanistic bleb-nucleation / cortex-failure mode
  - 2D-shell vs 3D-fill cortex geometry (27.1% shell) governing single vs multiple protrusions
  - Local R_C/R_X patterning → deterministic bleb positioning = emergent mechanical polarity scenario
  · n_parameters=18
  · n_equations=11
  · n_validation_targets=14

## p55 [3.5] Makito Miyazaki · experimental — Directional Bleb Formation in Spherical Cells under Temperature Gradie
- headline: Turnkey bleb + membrane-cortex-coupling target: rounded cell held by 2gamma/R=DeltaP blebs directionally when the cortex-tension field is locally broken; port the induced gamma asymmetry + breakable linker + cortex-damage/advection, not the thermal driver.
  - Bleb module: nucleate from local cortex-tension/membrane-adhesion imbalance (staged detach->rupture->inflate, Fig.4E)
  - Membrane-cortex coupling as explicit breakable linker layer (Bell-Evans detachment), replacing DCM's co-located membrane==cortex
  - Per-face heterogeneous gamma field + cortex-damage scalar D (local tension collapse = rupture)
  - Cortex-density Marangoni advection to low-tension side (conserved rho_cortex, gamma tied to actomyosin density) -> near-down/far-up asymmetry
  - Imposed axial gamma asymmetry gamma(theta)=gamma_base[1+A cos theta] as the DCM perturbation, A anchored to measured dF/F0 (no thermal model)
  - Directional-vs-multiple bleb dichotomy (asymmetry -> single polar bleb; uniform -> multiple) as discriminating gate
  - Cytokinesis coupling: localized high-gamma furrow band + ectopic contractility -> mis-partition (secondary)
  · n_parameters=6
  · n_equations=5
  · n_validation_targets=11

## p30 [3.5] Makito Miyazaki · research — PHYSICAL REVIEW RESEARCH 5, 013208 (2023)
- headline: Confined active-fluid model gives DCM a mechanistic recipe to make cortical tension γ an emergent actomyosin-density field with in-plane flow, plus flow→contractile-ring wave selection.
  - Cortex tension: replace lumped γ with emergent density field σ_act=(ζμ)₀·ρ/(ρ₀+ρ) (Eq.3)
  - Cortical actomyosin FLOW: surface-Stokes force balance ∇²v+λ∇(∇·v)+∇f=v (Eq.6) driving in-plane flow
  - Cortical density transport: reaction-advection-diffusion ∂ρ/∂t+∇·(ρv)=D∇²ρ+k_p−k_d·ρ (Eq.4)
  - Cytokinetic ring as emergent cortical wave: flow→periodic-ring transition vs contractility (T≈1.5min)
  - Pe=(ζμ)₀/(Dγ) as grid-invariant contractility regime diagnostic (flow 38 / wave 115)
  - Membrane-shell polymerization contrast α_p=k_p^bulk/k_p^surf (<0.5→waves)
  - Volume-penalty U′=φ(φ−1)[φ−½−α₀(V/V_tar−1)] as cross-check of DCM enclosed-volume constraint
  · n_parameters=9
  · n_equations=6
  · n_validation_targets=7

## p19 [3.5] Taeyoon Kim · research — Mechanical Checkpoint for Cell Division in
- headline: Hands DCM a mechanistic cell-division protocol (V0-ramp growth + CRC/ISE mitotic force fields, stop at 58.74%) with AR 1.41(2D)/~2(3D) shape targets; the checkpoint science itself is ECM/FF.
  - Division module: CRC inward-equatorial + ISE outward-polar nodal force fields (the open cytokinetic-ring problem, mechanistic BC)
  - Growth-via-V0: ramp enclosed-volume target x2 over G1 to grow the cell before division
  - Elongation stop at 58.74% + furrow-close cutoff (1um from axis) as division protocol geometry
  - Soft-core differentiable node-segment repulsion (Eq.8) as membrane self-contact / soft cell-cell kernel
  - Bell tensile-only off-rate (Eq.7) as force-gated junction turnover to escalate aggregate T1 compaction
  - Harmonic enclosed-volume penalty (Eq.9) confirms/backs DCM turgor formulation, Pi = (kv/V0)(1-V/V0)
  · n_parameters=14
  · n_equations=6
  · n_validation_targets=9

## p76 [3.0] Makito Miyazaki · research — BIOPHYSICS AND COMPUTATIONAL BIOLOGY
- headline: DCM lens — p76  (Sakamoto, Izri, Shimamoto, Miyazaki, Maeda 2022, "Geometric trade-off ... actomyosin droplet motility", PNAS 119:e2121147119)

## p47 [3] Taeyoon Kim · research — Role of actin filaments and cis binding in
- headline: Sub-µm cadherin-AJ kinetic model: parameterizes DCM's lumped junction (cluster size/lifetime laws, weak-cis regime) and independently confirms cadherin kinetics don't rate-limit 24-48h compaction — a calibration/validation reference, not a force or scenario source.
  - Junction-as-load-sharing-cluster with emergent lifetime τ∝[cad]² (calibrates DCM's Rakshit catch-bond cluster)
  - External confirmation that cadherin junction kinetics (~6 s max) do NOT gate 24-48h compaction — corroborates S4 '(b) T1 needed'
  - Counterintuitive weak-cis/high-turnover → larger, more persistent clusters (escalate by density/recruitment, not by lowering off-rate)
  - Actin corralling → sub-diffusive cadherin confinement sets junction density ceiling (cortex-density → junction recruitment coupling)
  - Cad/actin tethering non-monotonic: high-affinity+small-footprint cooperates, ≥50% footprint depletes trans (junction↔cortex anchoring optimum)
  - p_ass=1−exp(−kΔt) baseline rate law + linker off-rates 0.1/10 s⁻¹ (α-catenin/vinculin) as zero-force junction parameters
  · n_parameters=17
  · n_equations=4
  · n_validation_targets=11

## p39 [3] Taeyoon Kim · research — 10274 |  Soft Matter, 2021, 17, 10274–10285
- headline: Contractile-cortex→ECM traction template + σ∝r⁻¹/r⁻² validation and emergent-viscoelasticity crosslink-turnover mechanism; corroborates DCM's Kim 1/r work but touches only the ECM/traction side, not cell interior.
  - Permanent FA harmonic link (Eq.1) as the cortex-node→fiber traction interface — mechanistic cell→ECM force path
  - σ(r) radial-stress diagnostic with derivable σ∝r⁻¹ (2-D) / r⁻² (3-D) dimensionality self-check
  - Emergent ECM viscoelasticity from Bell slip-bond crosslinker turnover (Eq.2), not a lumped modulus
  - Plastic/irreversible matrix remodeling as a DCM ECM validation behavior (little retraction on cell removal)
  - Counter-intuitive stronger-contraction→faster-relaxation caution for active-cortex + turnover-limited coupling
  - Tensional-homeostasis: finite steady traction sustained even in a relaxing/plastic matrix
  · n_parameters=9
  · n_equations=4
  · n_validation_targets=16

## p32 [3] Taeyoon Kim · review — Weldon School of Biomedical Engineering,
- headline: Kim-lab review is a DCM design brief: names the mechanistic origin of four missing inelastic behaviors (superelastic cortex plateau, two-regime nucleus, viscoplastic compaction, poroelasticity) plus a validation-target table, but ships no portable equation.
  - Superelastic cortex: γ tension-plateau via fluidization + IF backstop (γ-generation ceiling)
  - Nucleus two-regime stiffening: chromatin soft <3µm, lamin-A load-bearing above + 2× lamina area reserve
  - Viscoplastic irreversible offset from force-dependent crosslinker unbinding / turnover = physiological compaction-rate (T1) clock
  - Poroelastic cytoplasm: biphasic drainage, exp(<0.05s)->power-law relaxation, D_p=E·ξ²/η
  - Power-law (soft-glassy) rheology observable G'(ω),G''(ω) instead of single-time Maxwell
  · n_parameters=12
  · n_equations=3
  · n_validation_targets=10

## p28 [3] Makito Miyazaki · review — Molecular Crystals and Liquid Crystals
- headline: Review supplies a cytokinetic-ring closure law (rate ∝ diameter), a cell-sized confinement threshold, and spindle active-viscoelastic targets — hitting DCM's division and confinement open problems, but validation-grade not core-moving.
  - Cytokinetic furrow: constant-strain-rate ring closure giving dD/dt ∝ D₀ (division module, currently absent)
  - Furrow force-balance vs turgor: λ_ring/r ≳ Π₀ constriction condition on the physiological baseline
  - Cell-sized spherical confinement BC (R_conf ≈ L_p ≈ 20 µm threshold) as a new DCM scenario
  - Three-timescale active viscoelasticity target (elastic<1s/viscous~10s/active~min, |G*| 0.1→1 nN/µm) for stiff/mitotic inclusion
  - Anisotropic + anomalous-Poisson (ν≈2.0, ~2×) constitutive pattern for active inclusions
  - Uni-lamellar giant-liposome shell as membrane single-bilayer provenance reference
  · n_parameters=13
  · n_equations=4
  · n_validation_targets=10

## p07 [3] Taeyoon Kim · research — Covalent cross-linking of basement
- headline: Bell off-rate+rebind bond turnover is what makes deformation permanent vs elastic — the missing mechanistic driver for DCM's aggregate-compaction-RATE / cadherin T1 problem, plus a confinement/protrusion scenario and clean shape observables.
  - Cadherin junction PLASTICITY = Bell off-rate + rebind → permanent (not elastic) aggregate deformation, supplying the bond-turnover driver DCM's open compaction-RATE / T1-rearrangement problem needs
  - Confinement scenario: single cell in a tunable-plasticity surround; protrusion extension gated by environment YIELD not stiffness
  - Stiffness ⟂ plasticity decoupling as a hard validation invariant (turnover moves permanent-strain, leaves initial modulus flat)
  - Membrane edge length-window constraint [4 nm, 800 nm] as an extensional cap/floor enabling force-concentrating protrusions
  - Localized self-concentrating radial protrusive-force patch (25 pN × 40 nodes ≈ 1 nN) as a bleb/invadopodium driver
  - 3D shape index/sphericity Ψ (from 2D circularity 4π·area/perimeter²) + boundary strain metric ε=(d−r0)/2r0 as DCM shape observables
  · n_parameters=15
  · n_equations=4
  · n_validation_targets=11

## p03 [3] Taeyoon Kim · research — Determinants of Fluidlike Behavior and Effective Viscosity in Cross-Li
- headline: Turnover-set effective viscosity + biphasic yield-stress creep gives DCM a derivable cortex-rheology closure and a turnover-limited plastic-flow rate law for the compaction-rate open problem.
  - Cortex surface viscosity SET BY TURNOVER (η_eff ∝ 1/turnover, Eq.6) — removes a magic constant
  - Biphasic YIELD-STRESS creep closure for the cortex (linear-viscous below σ_th≈1 Pa, plastic above)
  - Turnover-limited plastic-flow RATE law (ε̇=k⁻·Δε_avg) for the loose→compact compaction-rate open problem
  - Turnover-gated T1/junction rearrangement clock instead of stress-escalated pull
  - Nonaffine 5%/99% load-pathway warning → reject naive linear-Maxwell cortex
  - Uniaxial cortex creep-test protocol as a direct DCM rheology validation gate
  · n_parameters=8
  · n_equations=4
  · n_validation_targets=12

## p70 [2.5] Makito Miyazaki · research — Tug-of-war between actomyosin-driven
- headline: In-vitro extract paper: no DCM constitutive numbers, but donates an exponential-in-size percolation-maturation RATE law (τp=2^N·τ) that fits the DCM compaction-timescale gap, plus confinement + emergent nucleus-positioning scenarios.
  - Percolation-gated maturation kernel τp=2^N·τ (exponential-in-size slow timescale) as the mechanistic RATE law for the DCM loose→compact / cadherin-maturation open problem
  - Nucleus positioning as EMERGENT cortex-vs-bulk force balance (tug-of-war) — gives the passive nucleus inclusion a real positioning role + a size-dependent centered↔edge transition to reproduce
  - Contractile centripetal-wave / ring line-tension (T≈46s, v≈1–5 µm/s) as a cytokinetic-ring / cortical-flow seed
  - Cell-sized rigid confinement BC scenario (node-face wall, quasi-2D h/D 0.3–0.6) — the missing DCM confinement-response test
  - Boundary–cortex coupling strength as the first-class control knob that flips internal-body polarity and converts contraction into net migration (PEG-center vs Arp2/3-edge)
  - Percolation occupancy probability p=1−[1−p_site^N]^{T/τ} as a fine-grained replacement for any binary junction latch (forbidden cad_mult switch)
  - Timescale-crossing criterion Rc/L=log₂(T/τ) as a derivable size-dependent regime boundary (no tuning)
  · n_parameters=10
  · n_equations=6
  · n_validation_targets=6

## p61 [2.5] Taeyoon Kim · research — Disordered actomyosin networks are sufﬁcient to
- headline: Continuum active-gel half gives DCM an exponential cortical-tension buildup law, a per-face localized-activation upgrade (cap/ring→cytokinesis), and cortical-flow validation targets — but cooperativity constants (nH=11, ρc=0.56) are in-vitro and stay out of the runtime.
  - Active-stress buildup law σₐ=σ₀(1−e^(−t/τₐ)): exponential γ ramp on activation instead of step-jump
  - Promote cortex γ to a per-face field with a localized activated cap (heterogeneous cortex)
  - Ring/equatorial-band activation = cytokinetic furrow protocol for the open division problem
  - E/σₐ ~ O(1) balanced-regime dimensionless check coupling γ to membrane K_A/κ_m
  - Bending-invariance result citably justifies DCM's non-filament-resolving continuum cortex
  - Crosslink stiffening (filamin) attenuates contraction — effective-E modulator of cortical strain
  - Surface-divergence observables ε=∇·u, dε/dt=∇·v for cortical-flow convergence
  · n_parameters=7
  · n_equations=6
  · n_validation_targets=9

## p60 [2.5] Taeyoon Kim · research — Interplay of active processes modulates tension
- headline: Kim-lab BD actomyosin paper gives DCM a turnover-gated dynamic cortex-tension law plus stress-scale (100–1000 Pa) and pulsed-constriction (~100 s, τC=28s/τD=155s) validation anchors — no direct kernel port.
  - Cortex γ as turnover-gated dynamic variable via active-kinetic-spring law (γ ∝ n_bound·f_motor), not a static set-point
  - Actin turnover as the anti-aggregation stabilizer: nondimensional control Π=k_turnover/v_motor collapses stress profiles (master curve)
  - Sustainability S diagnostic: cortex-collapse rate ∝ (1−S), classifies sustained vs collapsing regime
  - Turnover-pulse protocol for reversible foci → cytokinetic-ring / apical-constriction pulsation at ~100 s
  - Force-sensitive crosslinker feedback (Bell k_off) coupling face tension back to n_bound softening
  - Sigmoid kinetics (Eqs 1,2) to fit DCM foci/compaction τ against τC=28s, τD=155s
  - 100–1000 Pa TFM stress-scale sanity anchor for DCM cortical tension
  - Critical crosslinker:motor≈100 boundary where critical turnover→0
  · n_parameters=13
  · n_equations=4
  · n_validation_targets=7

## p58 [2.5] Makito Miyazaki · research — Myosin-Driven Advection and Actin Reorganization Control the
- headline: Turns DCM's uniform cortex γ into a mechanistic advection+turnover density field (γ_f=γ_0·ρ_f/ρ_0), adds a confinement-response scenario, cross-checked by emergent flow ≈1.2 µm/s and pulse period ≈101 s — but phase-field/bulk-collapse/dimensionless knobs do not transfer.
  - Advection-driven cortical-tension redistribution: make γ a per-face field γ_f=γ_0·(ρ_f/ρ_0) advected by tangential cortical flow (Eq 1 closure)
  - Cortex-density source/sink turnover k_p − k_d·ρ ported to Laplace–Beltrami surface mass balance (Eq 2 on the mesh)
  - Curvature/corner tension-concentration as an emergent M1 signature under confinement
  - Pulsatory cortical-flow instability (Kumar–Bois–Jülicher–Grill) → emergent cortex oscillation
  - Cell-sized rigid confinement scenario: boundary asymmetry transferred to cortex/shape (P1)
  - Contractility↑ vs polymerization↑ mechanistic perturbation battery (CalA/CytoD/VCA analogues)
  · n_parameters=6
  · n_equations=2
  · n_validation_targets=8

## p53 [2.5] Taeyoon Kim · research — Rapid assembly of a polar network architecture
- headline: Hands DCM a mechanistic, excitable pulsatile recipe for cortical-tension γ generation (RhoA→delayed-myosin→transient γ, 0.3 µm/s flow) plus polarity-driven force efficiency — but its filament-scale Table S4 stiffnesses must NOT be ported into the surface mesh.
  - Cortex γ generation: replace static uniform γ with an excitable RhoA-pulse surface field γ(x,t) (period 30s, myosin delay d_M=5s, duration τ_M=15s)
  - Polar cortical architecture → anisotropic/oriented surface tension tensor; barbed-out polarity transmits force farther and avoids aggregate collapse
  - Cytokinetic ring as an equatorial RhoA band reusing the same γ(x,t) pulse machinery (DCM division open problem)
  - Contraction onset lag ~3s → γ ramp rise-time law after each pulse
  - Cell-cycle γ_0 modulation (interphase 0.15 ↔ mitotic 1.5 mN/m) with temporal profile, contraction drop >10% at cytokinesis
  · n_parameters=8
  · n_equations=5
  · n_validation_targets=8

## p51 [2.5] Taeyoon Kim · research — Mechanical Model for Durotactic Cell Migration
- headline: Coarse 2-point cell model — no DCM cell mechanism, but its viscoelastic Kelvin-Voigt triangulated ECM (stretch+bend+per-edge drag) and decision-free emergent-durotaxis scenario/validation suite are portable to DCM's substrate/ECM layer.
  - Viscoelastic Kelvin-Voigt triangulated ECM/substrate: stretch (Eq9) + triplet bending (Eq10) + per-edge relative-velocity drag ξ (Eq11) — the per-edge dashpot is the novel ingredient DCM's ECM/contact layer lacks
  - Constant-strain-energy emergent traction: traction magnitude solved from substrate resistance under a fixed active energy budget (torque, Eq7), not prescribed — stiffer substrate yields larger traction for free
  - Decision-free emergent durotaxis as a DCM validation/design constraint (no stiffness-conditional parameter, mechanism must EMERGE)
  - Spatially-heterogeneous ECM stiffness field (sharp step + gradient) as a durotaxis test-bed scenario for an adherent DCM cell
  - Kelvin-Voigt short-vs-long-timescale force sustainment as the elasticity-vs-viscosity control of whether a substrate holds cell traction
  - Cell-free substrate calibration protocol (5% ±x stretch, center fixed) as an acceptance gate for the viscoelastic-ECM module
  · n_parameters=11
  · n_equations=4
  · n_validation_targets=12

## p45 [2.5] Taeyoon Kim · research — Characterization, Enrichment, and Computational
- headline: Triangulated-mesh nucleus with V/A conservation plus a NEW LINC cortex-to-nucleus tether and a 2300/3200/6200 Pa compression-stress ladder to validate DCM nuclear-coupling.
  - LINC-complex tether: attractive harmonic springs coupling 20% of cortex/perinuclear nodes to the nucleus surface (NEW cortex→nucleus mechanical coupling DCM lacks; ~2x nuclear compression resistance)
  - Nucleus compartment constraint potentials in SimuCell3D form: enclosed-volume U_V=1/2 kappa_V(V-V0)^2/V0 + surface-area U_A conservation on the triangulated nucleus sub-mesh, anchored to a measured stress-strain curve
  - Whole-cell AFM indentation scenario (10um sphere, <=400nm, 800nm/s, Hertz nu=0.5) to map cytoskeletal/LINC state -> effective Young's modulus
  - Osmotic deswelling (PEG300) as a physiological perturbation of the enclosed-volume/turgor constraint (lower V0 -> compaction/stiffening)
  - Nucleus-compression validation scenario S1: three arms (nucleus only / +perinuclear network / +LINC) reproducing the 2300->3200->6200 Pa stress ladder
  - Cell-sized microfluidic constriction (7um @30deg) as a DCM confinement-response scenario
  · n_parameters=20
  · n_equations=6
  · n_validation_targets=7

## p29 [2.5] Taeyoon Kim · research — Mechanics and Morphology of Proliferating Cell Collectives with Self-I
- headline: Rigid-rod colony paper offers DCM no new surface force, but a portable stress-inhibited growth law g=(1/τ)e^{−λp}, a large-dt complementarity contact solver validating the IPC direction, and a closed-form Darcy oracle for growing aggregates.
  - Stress-inhibited growth law dV0/dt=(V0/τ)exp(−λ·p_local) for aggregate proliferation (Eq.7/2b)
  - Complementarity/Lagrange-multiplier contact solver (min frictional dissipation, no penalty) — independent validation of DCM's IPC large-dt path
  - Darcy tissue-flow closure ξu+∇p=0, ∇·u=g as an emergent-physics audit of coarse-grained node velocities/pressures
  - Closed-form radial pressure/velocity p(r)=(R²−r²)/4, u(r)=r/2 (Eq.5) as a no-fit analytic oracle
  - exp→linear colony-radius crossover R~(2/λ)^{1/2}t (Eq.8) as a growing-aggregate validation target
  - Lagrangian inheritance division rule V̄0=2^{frac[∫g dt/ln2]} — stress-history-driven division timing (Eq.4)
  · n_parameters=7
  · n_equations=6
  · n_validation_targets=9

## p18 [2.5] ? · ? — Durotactic Migration Driven by Anisotropic Matrix
- headline: Same-lineage (SimuCell3D Eq.8) collective-contraction/escape story maps onto DCM's open aggregate-compaction-RATE problem — adhesion sustains contraction by delaying escape — but the coarse 2-point-agent, 2D, model-output cell offers oracles/protocols, not portable cell-shape physics.
  - Eq.8 SimuCell3D cell–cell repulsion+tent-adhesion as coarse-grained ORACLE for our fine-grained Rakshit cadherin junction (κ_a,int=1e-5 N/m cross-check, NOT a runtime setter)
  - Collective-remodeling gate: single cell cannot compact, cohesive cluster can — adhesion sustains contraction by delaying escape → maps to open aggregate-compaction-RATE / cadherin T1-rearrangement problem
  - Compaction rate is escape/turnover-kinetics-limited (junction lifetime), not contractile-force-magnitude-limited — supports making DCM compaction rate emerge from cadherin maturation/detachment
  - Contractility × adhesion are NON-ADDITIVE, saturating levers (1.5Fc≈adh alone) — guard for the DCM γ×κ_cadherin sweep
  - F_r/F_θ radial/circumferential tension decomposition (Eqs.19–20) as a new DCM aggregate compaction-anisotropy observable (target ≈2)
  - Circular-cluster 2×2 (±adhesion ×±1.5 contractility) + single-cell-vs-cluster control as ready-made N=400 aggregate-compaction protocols
  · n_parameters=4
  · n_equations=3
  · n_validation_targets=5

## p15 [2.5] Taeyoon Kim · research — Dissecting Molecular Origins of the Mechano-Adaptive
- headline: Sub-cellular SF paper (FF-native), but yields a mechanistic rate-dependent cortex γ(t) law — plateau∝myosin density, recovery-rate∝motor speed — that supersedes DCM's static γ.
  - Rate-dependent cortex γ(t): plateau ∝ myosin density, recovery rate ∝ motor walking speed (decoupled) — replaces static/lumped-viscoelastic γ
  - Buckling/wrinkle criterion: cortex wrinkles when compressive strain rate exceeds v_walk/ℓ_ref (grid-invariant, no fitted constant)
  - Stress relaxation sourced from force-dependent bond turnover (Bell/PCM), not a surface-viscosity dashpot
  - New DCM parameter v_walk (~200 nm/s in-cell) sets cortex recovery timescale τ_rec=ℓ_ref/v_walk
  - Sharp falsification tests: recovery rate must be density-INDEPENDENT; recovery time must be bending-stiffness-INDEPENDENT
  · n_parameters=2
  · n_equations=2
  · n_validation_targets=7

## p10 [2.5] ? · ? — Fine-tuning of material properties by catch bonds
- headline: Molecular catch-slip cross-linker paper (FF-native); for DCM it supplies conceptual guidance for the cadherin junction turnover/rearrangement and a ~3.6x cortex-tension amplifier, but no direct surface-mesh constitutive law.
  - Cadherin junction: force-biased turnover/redistribution of catch-bond clusters onto high-tension nodes (T1-like), not static maturation
  - Cortical-tension amplifier: catch-slip cross-linker redistribution multiplies motor-transmitted stress ~3.6x (candidate FF-derived multiplier on myosin->gamma for the gamma-floor)
  - Connectivity switch: catch-slip benefit is OFF at low connectivity, ON above RCL~0.01 -> emergent loose->compact compaction switch tied to contact density, not a hand-set escalation time
  - Turnover is essential: remove redistribution and catch-slip falls below slip (M1 falsification test)
  - Three-knob escalation map: raise Fmax (not lifetime) for compliant-yet-strong high-yield-strain contacts needed for large-deformation compaction
  - Eq.1 catch-slip off-rate ported per junction node-pair via (tau0, tau_max, Fmax) reparameterization
  · n_parameters=8
  · n_equations=2
  · n_validation_targets=7

## p04 [2.5] Taeyoon Kim · research — Balance between Force Generation and Relaxation
- headline: Cortical tension is a generation-vs-relaxation balance: add a turnover/curvature-gated relaxation channel to DCM's γ field to fix spurious aggregation and irreversible compaction.
  - per-face γ field with turnover/severing relaxation channel
  - curvature-gated local tension relaxation (continuum image of Eq. 1 severing)
  - turnover as homogeneity-preserving cortex stabilizer
  - reversible-vs-irreversible diagnosis for DCM aggregate-σ compaction
  - per-face tension-evolution ODE (generation minus relaxation)
  - cortical-patch pulsed-contraction validation scenario
  · n_parameters=4
  · n_equations=3
  · n_validation_targets=5

## p02 [2.5] Taeyoon Kim · research — Dynamic Role of Cross-Linking Proteins in Actin Rheology
- headline: Fine-grained network paper, indirect for DCM — but supplies the exact missing mechanism (turnover→partner-swap plastic reorganization) for the cadherin-junction T1 rearrangement that must drive loose→compact aggregate compaction at the real rate, plus a turnover-clocked stress-relaxation closure and a cortex strain-stiffening constitutive law.
  - Cadherin-junction T1 rearrangement = unbind-above-~20pN then rebind to a DIFFERENT partner (turnover=1), the fine-grained plastic-reorg engine for aggregate compaction RATE
  - Aggregate-σ / compaction relaxation is cross-link-TURNOVER-limited not drag-limited: clock = 1/(2 k0_ub), finite residual tension plateau b/a (Eqs 5-7)
  - Load-sharing supportive framework: 25% of junction bonds bear ~85% of load and preferentially rupture — cluster load concentration
  - Cortex nonlinear surface constitutive closure: strain-stiffening sigma=c1*g^0.4+c2*g^4, K~g^3, plastic gamma(t) (Eqs 4,9,10)
  - Negative result: SKIP cadherin/crosslinker domain-unfolding — at physiological gamma_eff~0.1/s unbinding dominates, unfolding cannot relax stress or change topology
  - Bell slip off-rate (k0_ub=0.115/s Ferrer2008, lambda_ub=1.04e-10 m) as CROSS-CHECK for the ratified catch-bond junction, never a replacement
  · n_parameters=13
  · n_equations=7
  · n_validation_targets=12

## p75 [2.0] Taeyoon Kim · research — Nature  |  Vol 626  |  15 February 2024  |  635
- headline: DCM lens — p75  (Fan, ..., Taeyoon Kim, ..., Török 2024, "Matrix viscoelasticity promotes liver cancer", Nature 626:635)

## p74 [2.0] Taeyoon Kim · research — 1Living Matter Department, AMOLF, Amsterdam, The Netherlands. 2Institu
- headline: DCM lens — p74  (Mulla, Jung, Kim, Tans, Koenderink 2022, "Weak catch bonds make strong networks", Nat Mater 21:1019)

## p69 [2] Taeyoon Kim · research — Comp. Part. Mech. (2015) 2:317–327
- headline: FF-native reconstituted-cortex paper; for DCM a coarse-graining anchor that derives interphase γ≈0.3 mN/m from emergent σ_max=600 Pa and independently corroborates the γ-ceiling=motor-generation-limit.
  - Cortical γ anchor: γ=σ_max·h_cortex → 600 Pa × 0.5µm ≈ 0.30 mN/m corroborates DCM interphase 0.15
  - γ ceiling is motor-count/stall-limited (plateau ~600 Pa) — external corroboration of project's γ-floor=generation-limit finding
  - γ as a relaxing state (σ_max generation vs S sustainment), not a frozen constant
  - Crosslinker-as-molecular-clutch → load-dependent (Bell) cortex surface viscosity
  - Myosin-density scaling law γ ∝ R_M^0.6 for interphase→mitotic switch (flags 10× gap needs multi-factor)
  - Cortex-patch clamp-and-read + oscillatory surface-rheology validation protocols
  · n_parameters=13
  · n_equations=5
  · n_validation_targets=9

## p68 [2] Taeyoon Kim · research — sharing, adaptation, distribution and reproduction in any medium or fo
- headline: Plant cortical-MT stress-feedback paper donates a computed-stress projection law and a per-face nematic order-tensor template for anisotropic cortex remodeling, plus MT kinetics for DCM's microtubule compartment — but wall-regime (3-15 MPa, 2D) and cytoskeleton mismatch make it a secondary conceptual donor.
  - Eq-1 stress-projection onto face principal-stress direction (Mohr transform, computed not imported)
  - Stress-feedback -> per-face nematic order tensor Q for anisotropic cortex/MT reinforcement
  - Linear->FF-derived saturating constitutive law biasing filament kinetics under tension
  - Transient homogenization (S_p->0, density halves) signature during ~50min reorientation
  - Nematic order parameter S_p as DCM anisotropy observable
  - Three-state MT dynamic-instability rate constants for --microtubules compartment
  - 1.5:1 stress-anisotropy alignment threshold as derived remodeling trigger
  · n_parameters=13
  · n_equations=3
  · n_validation_targets=12

## p67 [2] Taeyoon Kim · research — Cytoskeletal Deformation at High Strains and the Role of Cross-link
- headline: Filamin single-molecule crosslink paper is mostly FF; for DCM the one real lever is a shear/torsion angle-knockdown + two-barrier form on the cadherin junction off-rate, touching the open T1/compaction-rate problem.
  - Junction shear/torsion knockdown: rupture force flat below ~45deg then drops sharply (Fig 5a) -> makes shear-loaded sliding junctions rupture first, seeding emergent T1 rearrangement in the open compaction-rate problem
  - Two-barrier loading-rate-dependent off-rate (Fig 5b, ln(rate) two-slope) so the cadherin junction has correct lifetime across physiological-slow vs perturbation-fast regimes
  - Bell force-dependent off-rate k_off=k_off0*exp(F*x/kT) as the junction unbinding form (functional form only; filamin magnitudes stay FF)
  - Barrier-sanity guard dG=0.5*kappa*x^2 to keep any junction parameters physically sane (few kT), enforcing the magic-number block
  - Two-cell junction shear-vs-tension pull + loading-rate sweep as DCM micro-validation protocols for the above
  · n_parameters=12
  · n_equations=4
  · n_validation_targets=2

## p63 [2] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:345–355
- headline: FF-class paper; for DCM it supplies a χ≈0.005–0.15 cortex-tension efficiency oracle, a biphasic durotaxis law (knee ~300–1000 Pa), and filament-side confirmation that real aggregate compaction is junction-kinetics-rate-limited.
  - Cortex γ efficiency oracle: sustained tension is χ≈0.15 of Σ stall (0.005 without crosslinkers) — reframes γ-floor as an efficiency ceiling, not a missing datum
  - Biphasic rigidity-sensing law γ_active∝E_sub then plateaus above E*≈300–1000 Pa → durotaxis scenario
  - Loose→compact compaction is crosslinker/junction-kinetics-rate-limited (not surface-tension-drag) — corroborates cadherin T1 escalation for compaction RATE
  - Crosslinker connectivity, not motor density, gates SUSTAINED cortical tension (S≈1 vs <5 s collapse)
  - Buckling (Euler, Eq.7) rectifies net-contractile sign — justifies DCM always-contractile γ from first principles
  · n_parameters=11
  · n_equations=3
  · n_validation_targets=6

## p50 [2] ? · ? — Molecular Biology of the Cell • 34:ar67, 1–10, June 1, 2023
- headline: Weak for DCM overall (FF/substrate-migration paper), but the PFAC load-sharing principle directly validates and constrains DCM's cadherin junction-cluster robustness and reframes the compaction-rate open problem.
  - PFAC: parallel catch-bond load-sharing gives a junction cluster a robustness threshold (~10 bonds crossover, ~30 = never fails/self-healing) — validates & constrains DCM cadherin cluster
  - Reframe compaction-RATE open problem as cluster-level PFAC re-engagement dynamics, not a lumped single-bond off-rate
  - Contact stress = Σ(bond force)/patch area: total force ↑ to plateau while stress ↓ as contact area widens (add force-vs-stress split to node-face/junction reporter)
  - Shared-compliant-load requirement (E2): parallel bonds must pull on one deformable patch or PFAC vanishes
  - Catch-slip vs slip is qualitatively insensitive (Bangasser 2013) — keep Rakshit catch-bond, robustness conclusions hold
  - Parallel-bond load-sharing unit test (P1) sweeping N=1,3,6,10,15,20,30 with constant bond budget as a DCM junction-cluster validation
  · n_parameters=6
  · n_equations=6
  · n_validation_targets=6

## p48 [2] Taeyoon Kim · research — PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1013
- headline: Filament-scale lamellipodium BD paper — peripheral for DCM; usable only as a Bell slip-bond clutch form, a contraction-vs-adhesion flow/traction-plateau validation lens, and a compliant-substrate durotaxis scaffold.
  - Bell slip-bond off-rate (Eq S10) as the force-dependent form for DCM's substrate FA clutch (distinct from cadherin catch-bond)
  - Contraction-vs-adhesion-friction competition sets retrograde/cortical flow speed — validation lens for a future DCM cortex-flow term (20-50 nm/s)
  - Turnover keeps surface homogeneous: adopt CoV-of-density (0.2-0.3 healthy, >1.0 stalled) as a DCM remesh/compaction diagnostic
  - Force-dependent traction PLATEAU (~0.2-0.3 nN saturating) for adherent single-cell runs
  - Compliant elastic-substrate abstraction (kappa_eff 2.4e-4 N/m) as a durotaxis/adherent scenario scaffold
  - Curvature-triggered remesh weight as the continuum analog of buckling-induced severing (Eq S9)
  · n_parameters=8
  · n_equations=4
  · n_validation_targets=7

## p44 [2] Taeyoon Kim · research — Morphological Transformation and Force
- headline: Sub-grid FF paper whose one DCM lever is a rigorous load-sharing→per-bond-force→off-rate→lifetime law that validates the cadherin cluster and argues compaction rate is bond-turnover-gated (supports T1 escalation, not aggregate-σ).
  - Junction cluster stability law: load-sharing → per-bond force → force-dependent off-rate → cluster lifetime (validates DCM's explicit N-bond cadherin cluster over a lumped single bond)
  - Compaction RATE gated by bond turnover/remodeling, not driving force — fine-grained confirmation that DCM's loose→compact aggregate needs cadherin T1 rearrangement escalation, not more aggregate-σ pull
  - Catastrophic vs graceful junction failure: below a critical bond count the cluster disintegrates non-linearly (a lumped junction cannot reproduce this — it's the acceptance test)
  - Cortex γ generation has a mechanistic density-limited ceiling (~4 nN, ∝ motor·crosslinker density); buckling-derived tension is sub-grid to DCM and must come from FF, not tuning
  - Reusable σ-based, rate-of-change-thresholded compaction-time protocol (definition transfers to DCM aggregate/cortex timing; the ~20 s number does NOT)
  - Non-monotonic turnover-protection optimum (intermediate ξ_d≈0.6) if DCM adds a cortex/junction dissolution-protection term
  · n_parameters=8
  · n_equations=3
  · n_validation_targets=5

## p43 [2] Taeyoon Kim · research — Computational Analysis of Viscoelastic Properties of
- headline: A fine-grained FF-ancestor actin-BD paper; for DCM it supplies a cortical-material-law cross-check — prestress-stiffening exponent 0.85, threshold 0.1 Pa, cell G'~10³ Pa — plus a cited justification for deterministic overdamped mesh dynamics and an oscillatory-shear rheology validation protocol.
  - Cortex prestress-stiffening law G'∝τ₀^0.85 above threshold τ₀≈0.1 Pa (~100× stiffening relaxed→pretensioned) as a falsifiable acceptance test on the cortical-tension γ field
  - First-principles justification for DCM's deterministic noise-free mesh dynamics: thermal fluctuation negligible when l_c≪l_p (0.393µm vs 10-20µm)
  - Oscillatory-shear whole-cell rheology scenario + Eq-10 modulus extraction (G*=|τ|/|γ|, G'=|G*|cosφ) as a new DCM validation protocol
  - Loss modulus G''/phase-delay targets to calibrate DCM membrane surface viscosity
  - Supportive-framework heterogeneity: ~25% linkers carry ~70% stress (affine) — argues against uniform-γ/uniform aggregate-σ, aligns with load-sharing cadherin clusters
  - Yielding sanity: permanent-crosslink stress catastrophe (>100 Pa artifact vs 1-30 Pa real) = what DCM must NOT reproduce, motivates remesh/T1 release valve
  · n_parameters=13
  · n_equations=4
  · n_validation_targets=6

## p37 [2] Taeyoon Kim · research — Molecular Biology of the Cell • 35:ar47, 1–12, April 1, 2024
- headline: 1-D growth-cone motor-clutch: narrow but real DCM payoff — a force-gated, pool-saturating reinforcement law to replace age-based cadherin junction maturation, plus a compliant-substrate durotaxis BC and validation shapes; keep Rakshit catch-bond, do not import the fitted bimodal k_add0.
  - Force-gated + pool-saturating reinforcement law (Eq 6) to replace age-gated cadherin junction maturation
  - f_available=(n_max-n)/n_max derivable saturation → junction equilibrium size with no magic-number cap
  - Compliant elastic-substrate BC (K_sub spring, Eq 9) → DCM durotaxis scenario
  - Reinforcement converts load-and-fail oscillation → stable steady state (junction stability regression)
  - Kinematic maturation criterion: interface slip drops to 20% of free (grid-invariant ratio, not a clock)
  - Equal-load-sharing Hookean trans-bond (Eq 4) per contact face-pair
  · n_parameters=8
  · n_equations=4
  · n_validation_targets=10

## p35 [2] Taeyoon Kim · research — Soft Matter, 2017, 13, 3213--3220 | 3213
- headline: FF-lineage actomyosin-network paper; for DCM it supplies a reusable foci-area compaction metric (Eq.2), a single-vs-multiple-foci over-drive classifier, and mechanistic corroboration of the non-monotone σ sweet spot and ~100 s mechanical compaction limit — no new force ports (severing stays in FF).
  - Foci-area contraction extent ξ=1−ΣπRi(t)²/ΣπRi(0)² (Eq.2) reused as DCM aggregate compaction observable on cell centroids
  - Single-focus vs multiple-foci end-state classifier as a DCM aggregate-σ over-drive failure gate
  - Non-monotone compaction vs drive/stiffness corroborates DCM σ-sweep sweet spot (σ≈5 mN/m), warns against monotone-σ
  - γ/σ generation set by motor–crosslink RATIO (critical κ∝R_ACP/R_M), not motor drive alone — guardrail for any γ-generation sub-model
  - Assemble-then-drive two-stage protocol maps to DCM --from-resting baseline then activate tension
  - Literature-standard steady-state rule (ξ̇<0.1ξ̇₀, rate=ξ_ss/t_ss) for the DCM compaction-RATE open problem
  - ~100–500 s mechanical contraction timescale confirms DCM ~100 s is the mechanical limit, not the 24–48 h cadherin-rearrangement time
  · n_parameters=3
  · n_equations=2
  · n_validation_targets=6

## p24 [2] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: F-actin turnover is the cortex's plastic stress-relaxation channel — give DCM's cortical tension a turnover-set Maxwell relaxation time, fixing elastic-only γ and explaining why aggregate-σ compaction jams without a rearrangement channel.
  - Cortex Maxwell viscoelasticity: relaxation time τ_relax≈1/k_t,A from F-actin turnover, replacing elastic-only γ
  - Turnover = the missing plastic-relaxation channel that explains why aggregate-σ compaction is drag-limited/jammed (T1-rearrangement analog)
  - Cytokinetic ring = localized cortex-aggregation instability (equatorial band of raised γ / suppressed turnover)
  - Cortical-flow speed bounded by unloaded NMII walking velocity <140 nm/s
  - Connectivity→cortex-stiffness ordering (⟨L_f⟩, R_ACP) as a moduli-derivation sanity check
  - Preassemble-then-activate protocol reinforces DCM --from-resting checkpoint discipline
  · n_parameters=2
  · n_equations=1
  · n_validation_targets=4

## p17 [2] ? · ? — Kinetic Control of Out-Of-Equilibrium
- headline: Cortex tension γ becomes a first-order turnover field γ(t)=γ_ref·k_app/(k_off·m_ref) — a leaky-integrator/high-pass upgrade to DCM's static γ set-point, giving out-of-equilibrium and pulsed contractility for free.
  - Make cortex tension a per-face turnover field: dm/dt=k_app−k_off·m, γ=γ_ref·m/m_ref (Eq.5) instead of a constant set-point
  - High-pass rate ceiling on γ changes: realized amplitude = commanded/√(1+ω²/k_off²), phase-lag arctan(ω/k_off) (Eq.4) — physical limit on how fast the cortex can ramp tension
  - Instantaneous set-point γ_target=k_app/k_off (Eq.3) reuses DCM's pressure-borne γ=ΔP·R/2 as γ_ref, kinetically gated
  - Step-response τ_cortex=1/k_off≈10–40s (Eq.7) as the cortex-tension charging time and a unit test
  - Assembly-before-contraction sequencing: ramp γ over d_M≈5s on fresh/remeshed cortex, never spike it (soft stability rule from Fig 5 non-monotonicity)
  - Out-of-equilibrium diagnostic Δ_γ (Eq.8): amplitude ratio of commanded vs realized tension, ~5 at 30s forcing
  · n_parameters=7
  · n_equations=5
  · n_validation_targets=7

## p13 [2] Taeyoon Kim · research — Swirling motion of breast cancer cells radially aligns collagen fibers
- headline: Not new DCM physics (engine is ECM/FF), but a strong emergent-collective-motion validation: swirling gives DCM's cadherin+cortex aggregate a tangential-shear observable and reframes the T1 compaction-rate problem as shear-driven neighbor exchange.
  - Swirling: emergent tangential>>radial interfacial collective motion as a validation mode for the cadherin-junction + cortex aggregate
  - Reframe compaction-RATE / T1 open problem as shear-driven tangential neighbor exchange (not radial pulling, not aggregate-sigma)
  - Junction+cortex+adhesion requirement triad as a multi-knockout mechanistic gate matrix (Ecad/p120, myosin/ROCK, beta1/FAK)
  - Basement-membrane-as-mechanical-insulator: stiff area-incompressible outer shell that screens cortical shear (membrane-cortex coupling scenario)
  - Aggregate-at-ECM-interface equatorial polar velocity decomposition as a runnable DCM scenario
  - Radial alignment index |cos(theta_f-theta_n)| as a DCM boundary-stress diagnostic
  · n_parameters=6
  · n_equations=2
  · n_validation_targets=10

## p09 [2] Taeyoon Kim · research — Reciprocal folding dynamics in cellular networks at the stroma-basemen
- headline: Tissue-scale continuum FEM of stiffness-contrast folding: no force law to port, but supplies a compliant-substrate/BM scenario, volume-gated buckling sanity, and stiffness anchors (E_Cell 8.6/E_Str 354/E_BM 756 kPa) as bounds, not fits.
  - Compliant stiff-sheet/BM-on-soft-substrate scenario: let DCM cortex+FA traction deform an elastic substrate (emergent, not prescribed) to unlock durotaxis/folding
  - Stiffness-contrast buckling intuition: cortical-shell wrinkling should be NON-MONOTONIC in cortex stiffness and volume-gated (peaks ~E_Cell 14 kPa)
  - Volume-gating sanity: surface folding vanishes as enclosed volume->0, corroborating physiological turgor/Pi0 must be ON
  - Blebbistatin contraction-OFF logic -> DCM myosin-OFF mechanism-attribution control run
  - Corroborates compaction-rate problem: >24h topographic reorganization confirms aggregate-sigma (~100s drag-limited) is too fast, needs slow rearrangement
  - REJECT (violates fine-grained rule): contractility-as-boundary-stress P_Cell=200kPa and single E_Cell modulus
  · n_parameters=9
  · n_equations=2
  · n_validation_targets=13

## p05 [2] Taeyoon Kim · research — Mobility of Molecular Motors Regulates Contractile
- headline: Fine-grained FF-native motor paper: no DCM equations or parameters, but supplies two roadmap guardrails — connectivity/T1 (not drag-σ) sets compaction rate, and γ is membrane-mobility-gated via two-arm myosin-II.
  - Compaction RATE is connectivity/rearrangement-limited, not drag-pull — reinforces T1 rearrangement (not more aggregate-σ) as the real DCM compaction lever
  - Cortical tension γ is mobility-gated: dominated by two-arm/antiparallel myosin-II whose force is drag-INSENSITIVE — guardrail against tuning cytosol η to fix γ
  - Cytokinetic-ring assembly is gated by node MOBILITY: physiological baseline is anchored/low-mobility membrane-coupled nodes (too mobile → ring fails)
  - Membrane surface viscosity is a legitimate physiological modulator of the γ source term (myosin-I membrane anchoring) — couple membrane shell to cortex, don't leave it inert
  - Bell's-law crosslinker off-rate confirms the catch/slip form DCM's cadherin junction already uses (cross-check only)
  - Q_A grid-density-SD reusable as a surface cortex-uniformity diagnostic on the mesh (uniformity-not-V/V0)
  · n_parameters=0
  · n_equations=1
  · n_validation_targets=0

## p73 [1.5] Makito Miyazaki · research — Optogenetic actin network assembly on lipid
- headline: DCM lens — p73  (Yamamoto & Miyazaki 2025, OptoVCA branched-actin density gates ABP function, Nat Commun 16:7583)

## p66 [1.5] Taeyoon Kim · research — Bulletin of Mathematical Biology (2019) 81:3301–3321
- headline: 2D rigid-substrate 2-point migration machine — no surface-mesh physics to port; yields only behavioral validation targets (PRW speed 50-60 µm/hr, jamming ~2300 cells/mm², h_p>h_v) and weak corroboration of cadherin-as-friction.
  - Cadherin junction ≈ effective inter-cell viscous friction (Eq.2 β·Δv) — external corroboration for our load-sharing cadherin cluster's low-frequency limit, NOT a lumped replacement
  - CIL as shared-adhesion-site competition/partitioning (substrate-mediated, no explicit repulsion) — conceptual contrast to our node-face + junction contact; do NOT add as a rule
  - MSD ballistic→diffusive PRW analysis (Eq.6-7 s=d lnMSD/d lnτ) as reusable DCM post-processing observable
  - Nematic Q-tensor + Saupe order parameter (Eq.8-9, h_p>h_v) as monolayer-ordering diagnostic
  - REJECT: constant-torque non-centripetal propulsion + instantaneous FA (violates fine-grained hard rule)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=8

## p65 [1.5] Taeyoon Kim · research — Computational Analysis of a Cross-linked
- headline: FF-native BD actin-network paper; for DCM it only supplies mechanistic FORMS — a rupture-limited power-law reformation law that reframes the open aggregate-compaction-RATE problem, plus a thermal-fluctuation caveat on junction off-rates — with no portable surface force, parameter, or validation number.
  - compaction-rate reframe: rupture-limited power-law reformation
  - Bell off-rate on fluctuating force (thermal amplification)
  - cortex gamma as percolating orthogonal-network shadow
  - persistence-bending calibration form
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=0

## p42 [1.5] Taeyoon Kim · research — F‑Actin Fragmentation Induces Distinct Mechanisms of Stress
- headline: Mostly an FF paper; for DCM it justifies one constitutive upgrade — an irreversible, reversal-gated plastic relaxation term for cortex γ, mirroring the T1-rearrangement rate-limiter the compaction problem needs.
  - Cortex constitutive upgrade: add irreversible, reversal-gated plastic relaxation to γ (distinct from reversible viscous η_s)
  - Two-channel relaxation split: within-cycle recoverable (ACP-unbinding analogue) vs across-cycle irreversible (severing analogue)
  - Mechanistic justification that aggregate-compaction rate-limiter must be an irreversible reversal/T1-gated event, not more σ-pull
  - Rate-dependent stochastic relaxation (must be a derived rate, not a per-step probability)
  - Eq (1) cyclic stress-drop as a DCM cortex-fatigue diagnostic + no-remodeling null gate
  - Bidirectional-vs-unidirectional cyclic-loading DCM scenario (relaxation only under strain reversal)
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=4

## p41 [1.5] Taeyoon Kim · research — Ding, Chou et al. eLife 2025;14:RP105236. DOI: https://doi.org/10.7554
- headline: Fine-grained Kim-lab myosin paper — DCM home is FF; only lever is a mechanistic efficiency ceiling (η≈0.08) and √-density scaling to guardrail future cortical-γ generation.
  - Cortical-γ generation: impose myosin-density efficiency ceiling η≈0.08 (Eq 3) as upper bound, reframes γ-floor as architecture/efficiency effect
  - Network vs bundle scaling: cortex active tension ∝ ρ_M^0.65 (√density, Eq 8/Fig 7B), NOT linear in myosin
  - Single-motor force quantum F_st≈5.7 pN/head × N_h=8 (≈45.6 pN/arm) as γ-estimator cross-check
  - Mitotic γ ramp sanity: 10× tension needs ~32× myosin recruitment under super-linear (0.65) scaling
  - Isometric section-averaged tension read-out at steady-state plateau confirms DCM γ measurement protocol
  · n_parameters=6
  · n_equations=3
  · n_validation_targets=3

## p40 [1.5] Taeyoon Kim · research — Durability of Aligned Microtubules Dependent
- headline: Wrong-engine paper (2D MT-kinesin gliding = FF, not DCM); the one durable DCM takeaway is the durability-not-formation principle that reinforces cadherin-T1 as the compaction-rate limiter.
  - Durability-not-probability principle: junction lifetime (not formation) sets aggregate compaction rate — reframes the open cadherin-T1 problem
  - Two-observable gate: compaction rate tracks junction disintegration time, invariant to engagement/formation knobs (Fig 7e port)
  - Coupling-time t0->t1 protocol as cadherin-junction lifetime measurement (Fig 6a)
  - MT persistence length 0.65-3.54 mm -> kappa_bend = Lp*kBT, only if explicit microtubule-inclusion beam law is added
  - REJECT depletion force (0.11 pN) and 2D nematic/skewness metrics — no MCF7 analog, non-transferable
  · n_parameters=4
  · n_equations=1
  · n_validation_targets=2

## p38 [1.5] Makito Miyazaki · research — Controlling Physical and Biochemical Parameters of Actin
- headline: In-vitro 2D actin-nucleation paper with near-zero direct DCM content; its one usable lever is making cortex/protrusion generation a membrane-surface-localized, face-normal, NPF-density-weighted source (sigma_NPF ~6500-25000/um2), everything else is FF-scale.
  - Cortex/protrusion generation as a membrane-surface-localized, face-normal, NPF-density-weighted source field (per-face sigma_NPF) instead of a global scalar gamma
  - Physiological membrane-nucleator areal density anchor 6500-25000 molecules/um2 for a future generation field (baseline-rule setpoint)
  - Barbed-end-at-membrane polarity: generation/protrusion term is added at the membrane and points along the face normal (Fig 1f)
  - Nucleator-type selects architecture: Arp2/3 normal-isotropic vs formin+fascin tangential-anisotropic cortex (future gamma-tensor two-mode switch)
  - Patterned heterogeneous sigma_NPF(face) as a DCM polarized-cortex/protrusion/durotaxis scenario, footprint copies the patch
  - Geometry-size sensing: >=4um patch gives directed growth, 2um does not (qualitative sanity for patch-protrusion)
  · n_parameters=5
  · n_equations=2
  · n_validation_targets=4

## p36 [1.5] Taeyoon Kim · research — 1548 | Soft Matter, 2020, 16, 1548--1559
- headline: FF-family motility-assay paper with low direct DCM yield; usable only as an acceptance-oracle/calibration for the lumped cortical-tension abstraction (intermediate-crosslink optimum, percolation arrest, turnover-limited compaction) — no force equation ports to the mesh.
  - Cortical-tension γ as emergent, intermediate-optimum network quantity (oracle for a future sub-grid γ generator, not a free monotone knob)
  - Turnover-gated 'frozen network' arrest corroborates DCM finding that compaction RATE is bond-turnover/T1-rearrangement-limited, not tension-magnitude-limited
  - Percolation-gated arrest (~3% cross-links) as acceptance oracle: cortex closure must not rise monotonically with cross-linker density
  - Peak-then-nonzero-asymptote tensile relaxation shape as validation target for aggregate-σ + junction-turnover
  - Eq.5 heterogeneity metric Q_A portable as a DCM aggregate-uniformity diagnostic (nodes/cells per grid), not a force
  - Keep DCM catch-bond junction — do NOT downgrade to the paper's slip-only Bell's-law cross-linker
  · n_parameters=2
  · n_equations=1
  · n_validation_targets=5

## p26 [1.5] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: FF-regime filament paper; for DCM only one real lever — irreversible topology reduction drives catastrophic contraction, backing the T1-rearrangement (not aggregate-σ) hypothesis for the compaction-RATE problem.
  - Aggregate compaction RATE: irreversible topology-reduction (T1 rearrangement) unlocks catastrophic contraction of a stuck high-connectivity aggregate — independent support that more force (aggregate-σ) does NOT set the rate
  - Diagnostic signature for the T1 escalation: reversible unbinding = gradual/stalls, irreversible severing = sharp self-accelerating collapse; trigger on max-loaded junction not the mean
  - Cortex tension as connectivity-gated, stress-relaxing material (biphasic in myosin turnover, non-monotone in crosslink density, capped by a rupture ceiling) — forward-looking note for γ-generation redesign only
  - Cytokinetic ring closure may stall on line-tension alone without a connectivity-reduction/severing term (anillin = contractile-ring crosslinker)
  - Stable-until-perturbed → catastrophic-compaction trigger protocol (gelsolin analog) as a DCM T1-escalation test scenario
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=4

## p21 [1.5] Taeyoon Kim · research — Dr. V. Yadav, Dr. A. P. Tabatabai, Prof. M. P. Murrell
- headline: Mostly FF-turf actin-nucleation paper; for DCM it yields one real lever — a turnover-set finite cortical stress-relaxation time (Maxwell surface rheology) — plus a minutes-scale turnover anchor and a mechanical-memory caution for the compaction-rate problem.
  - Cortex surface tension gains a turnover-set finite relaxation time (Maxwell surface rheology, τ_c) — DCM's cortex currently has no stress-relaxation clock
  - FRAP turnover time τ₁ → cortex remodeling/relaxation rate k_turn = 1/τ₁ (physiological order: seconds-to-minutes)
  - 'Mechanical memory': relaxation-vs-turnover decoupling (τ₂/τ₁≫1) as literature support that slow non-dissipative cortical stress storage is real — a conceptual handle (not a fit) for the compaction-RATE open problem
  - Per-face Maxwell surface-stress law dσ_c/dt = −σ_c/τ_c + K_A·dε_A/dt ported onto the existing K_A area-elastic surface
  - Bending-energy functional form E=(κ/2)(L/R²) as an order-of-magnitude sanity check on DCM κ_m (validation-only)
  · n_parameters=1
  · n_equations=3
  · n_validation_targets=4

## p77 [1.0] Makito Miyazaki · research — Spatial confinement of active microtubule networks
- headline: DCM lens — p77  (Suzuki, Miyazaki, Takagi, Itabashi, Ishiwata 2017, "Spatial confinement of active MT networks induces rotational cytoplasmic flow", PNAS 114:2922)

## p64 [1] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:1143–1155
- headline: FF-engine paper (crosslinked-filament rheology) with only a tangential DCM lever: a strain-stiffening σ(γ) shape (γ_c≈0.5) for an optional effective cortex law and a retro-justification of DCM's lumped-γ abstraction.
  - Cortex effective law: nonlinear strain-stiffening σ(γ) with critical strain γ_c≈0.5 (optional K_A_eff(γ), FF-imported not fitted)
  - Softer-element-dominates principle retro-justifies DCM's lumped-cortex γ/K_A abstraction for compliant crosslinkers
  - σ=ΣF/A_z shear-rheometry protocol as an FF→DCM upscaling template, not a native DCM scenario
  - High-strain cortex stress order-of-magnitude band ~200–2500 Pa as a sanity reference
  · n_parameters=2
  · n_equations=2
  · n_validation_targets=4

## p62 [1] Makito Miyazaki · research — Processive Nanostepping of Formin mDia1 Loosely Coupled with
- headline: Single-molecule formin nanostepping: an FF barbed-end input with zero surface-mesh/continuum content; DCM relevance is provenance-only.
  - Formin barbed-end step kinetics — belongs to FF, not DCM
  - Loose coupling multi-sized steps — FF ratchet refinement
  - Cortical-actin renewal timescale — provenance-only future hook
  - No DCM-portable equation/parameter/geometry/validation target
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p57 [1] Makito Miyazaki · experimental — Biphasic Effect of Proﬁlin Impacts the Formin mDia1
- headline: Single-filament formin/actin optical-tweezers kinetics — an FF paper, not DCM; only thin thread is that cortical γ generation is tension-accelerated (feed via FF→DCM, never lump into DCM).
  - Cortical γ generation is mechanosensitive: tension ACCELERATES formin actin assembly (p_open→1 above ~3 pN) — motivates FF-computed, tension-fed γ rather than static γ
  - Bell force-factor exp(−f·d/k_BT), d=5.5nm — the only transferable idea, but only as an FF per-filament gate, NOT a DCM node-force term
  - Route formin/barbed-end force-sensing to FF engine (filament-resolved), explicitly NOT DCM
  - Physiological-baseline guard: in-vitro 26°C purified rate constants must never enter DCM cortex config
  · n_parameters=0
  · n_equations=1
  · n_validation_targets=0

## p56 [1] Makito Miyazaki · method — Quantitative Analysis of the Lamellarity of Giant Liposomes Prepared b
- headline: Liposome-lamellarity optics/method paper: no DCM mechanics — only licenses the single-bilayer membrane assumption, ~5.2 nm thickness, and osmotic-baseline protocol; register as reference, not a validation gate.
  - Membrane: keep single-bilayer (unilamellar) idealization — this paper is its quantitative license (~97% unilamellar GUVs)
  - Membrane thickness anchor ~5.2 nm (alpha-hemolysin stem = bilayer thickness) confirming DCM's lipid-sheet constant
  - Osmotic-baseline protocol: matched-osmolality relaxed vesicle as a two-step init template separating passive bilayer from active Pi_0 turgor
  - Cell-sized spherical confinement + encapsulated-cortex reconstitution as a reduced dev-only DCM scenario (reconfirm on full stack)
  - Gateway citation to Sykes/Murrell cortex-in-GUV papers (refs 12/16/17/18) that carry the cortical-tension gamma DCM's open problem needs
  · n_parameters=2
  · n_equations=1
  · n_validation_targets=0

## p54 [1] Makito Miyazaki · experimental — on the number and size of liposomes formed
- headline: In-vitro liposome-yield fabrication study with no mechanistic membrane model, no moduli/tension, and no DCM observable — essentially no transferable DCM physics.
  - Latent-only: headgroup-dependent spontaneous curvature c0 as a passive Helfrich term (not needed, not parameterizable now)
  - Cell-sized (5-40 um) passive-membrane vesicle as a degenerate cortex-off DCM baseline geometry
  - Reconstitution-protocol metadata for a future cytoskeleton-in-liposome scenario (KB provenance only)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p52 [1] Makito Miyazaki · research — Microscopic Temperature Control Reveals Cooperative Regulation
- headline: In-vitro drebrin-E actomyosin paper offers one portable DCM idea — cortical γ activation is a cooperative Hill switch, not a linear duty — but no continuum, geometry, or MCF7-usable parameters.
  - Cortex γ generation: replace linear active-tension activation with a cooperative Hill sigmoid gate (form only, not n=4)
  - Regulator acts by bound-state allosteric duty modulation (no bind/unbind latch), mirroring anti-lumping cortex design
  - Occupancy p=[R]/(Kd+[R]) as the physical argument of the activation gate
  - Reversible step-response activation with ~5.4 s recovery as a dynamic cortex-gate test template
  - Temperature confirms 37 °C physiological setpoint sits at maximal actomyosin cooperativity
  · n_parameters=10
  · n_equations=2
  · n_validation_targets=4

## p49 [1] Makito Miyazaki · research — Kinetic Scheme of Myosin Phosphorylation by ZIP Kinase
- headline: Enzyme-kinetics paper (no mechanics/geometry) — supplies only a derived seconds-scale myosin-activation time constant to gate a γ cortex-tension ramp; low DCM impact, no fix for γ magnitude or compaction.
  - Cortex γ activation ramp: gate γ by an active-myosin fraction φ_act, γ(t)=γ_passive+(γ_max−γ_passive)·φ_act
  - Sequential 2-state MRLC phospho ODE (Eqs 1–2) reduced to a scalar cortex activation-fraction field
  - Derived seconds-scale activation time constant τ_activate≈1/kcat≈1.7–2.5 s replacing an ad-hoc/instant γ ramp
  - Saturated (pseudo-first-order) regime: activation rate ~kcat, not substrate-concentration-limited
  · n_parameters=7
  · n_equations=3
  · n_validation_targets=6

## p34 [1] ? · ? — Integr. Biol., 2015, 7, 1093--1108 | 1093
- headline: Review paper — no portable equations/parameters/targets; a citation-graph pointer that only corroborates DCM's active-compaction and single-junction directions.
  - Confirms aggregate compaction is active-contractility+cadherin driven (not passive-σ) — supports DCM's OPEN T1-rearrangement direction (via ref 128, DCIS-CPM ref 12)
  - Frames cell-cell junction as ONE mechanosensitive catch/slip adhesion interface — corroborates DCM's cadherin catch-bond junction design
  - Cortex-γ frontier: stiffen-then-soften (bending→stretching + transient-crosslink + motor percolation) picture, but only conceptual — chase primary refs 4/29/47
  - Nucleus mechanotransduction frontier (LINC stiffening, nN chromatin decondensation) as future E_nuc direction (refs 135/136)
  - Points to confinement, durotaxis-on-deformable-substrate, cell-pair-through-ECM scenarios — qualitative only
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p31 [1] Taeyoon Kim · research — PHYSICAL REVIEW RESEARCH 7, 013312 (2025)
- headline: Flat 2D myosin-motility assay with no membrane/turgor/cortex-tension/junction — zero direct DCM transfer; route to FF engine, DCM impact is a single advisory caution on engaged-motor-fraction dynamics.
  - Advisory only: engaged-motor fraction is emergent + self-limiting (Fig.1j) — if DCM ever GENERATES gamma (vs imposed gamma=DP*R/2), engaged density must be dynamical, not a setpoint
  - Pointer (not adoptable): active-nematic cortex direction — but flat 2D, no curvature/surface coupling; use the cited cortex refs [2] Grill, [3] Salbreux instead
  - No DCM module, equation, parameter, geometry, or validation gate transfers — all primitives are FF-engine (filaments/PCM motors/Langevin/drag)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p27 [1] Taeyoon Kim · research — Cytoskeleton, 2026; 83:407–426
- headline: Plant fine-grained actin BD model — FF lineage, not DCM; only a curvature-gated turnover-law analogy is portable, all numbers are plant/non-MCF7.
  - Curvature-gated cortex remodeling: adapt Eq.9 severing k=k0·exp(λθ) to local mesh curvature to give the open T1/remesh escalation a mechanistic (not lumped) trigger
  - Cortex-as-turnover-steady-state as a future sub-grid substrate for spatially-varying γ generation (density→tension), not the uniform γ=ΔP·R/2
  - Homeostatic-steady-state + single-mechanism-sweep discipline (measure at turnover steady state, seed-ensembled) for any DCM remodeling model
  - FDT node-noise amplitude (Eq.3, 2kBTζδ/Δt) as a cross-check only — inert at 7.5µm cell scale
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=0

## p25 [1] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: FF-engine paper; for DCM only a conceptual lever — the implicit-motor failure independently validates keeping cadherin-T1 over aggregate-σ as the compaction driver, plus two shelved surface-flow metrics.
  - Implicit-motor caution → validates keeping cadherin-T1 (not aggregate-σ) as the compaction DRIVER; σ matches end-state, loses rate/mechanism
  - Eq 4 morphology-homogeneity Q as a DCM cortical surface-density uniformity observable (grid-invariant form only)
  - Eq 1 velocity-correlation, remapped node→node, as a cortical surface-flow alignment metric (shelved until membrane–cortex coupling exists)
  - Curl-of-velocity-field as a cortical rotational-flow / swirl metric for polarization & division (needs surface flow field first)
  - Steady-state onset criterion (density-field self-correlation ≥0.8) as a DCM measurement-window gate
  - Noise-suppressed velocity sampling (~2 s cadence) for any DCM surface-flow observable
  · n_parameters=0
  · n_equations=4
  · n_validation_targets=0

## p23 [1] Makito Miyazaki · method — Bio-protocol 16(8): e5656. DOI: 10.21769/BioProtoc.5656
- headline: Wet-lab optogenetic branched-actin-on-flat-bilayer protocol — no equations, forces, or continuum params; effectively zero DCM impact, route to FF.
  - Conceptual only: cortex mechanical output tracks upstream areal regulator density (γ-generation intuition, UNQUANTIFIED here)
  - Conceptual only: patterned reversible on/off cortex actuation is physically realizable (scenario idea for bleb/ring/confinement)
  - No portable mechanism — molecular single-filament Arp2/3 assembly sits below DCM coarse-graining (FF territory)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p16 [1] Taeyoon Kim · research — Computational Modeling of Motor-Driven
- headline: Filament-bundle (MT/dynein) paper — essentially non-applicable to DCM's surface-mesh continuum; yields only one optional confinement scenario and two cross-check equations, no new mechanism or parameter.
  - Force-gated moving-wall confinement scenario (fixed + advancing end-cap, advance only when internal push > F_R) adapted to a pressurised MCF7 cell
  - Cylinder drag closed form (Eq. 3) as a cross-check oracle, not a runtime replacement for DCM's eta-derived node drag
  - Bell slip-bond off-rate (Eq. 8) as a high-force-limb cross-check on the DCM cadherin catch-bond (KU-4.2), never a replacement
  · n_parameters=0
  · n_equations=3
  · n_validation_targets=0

## p08 [1] Makito Miyazaki · research — Chimeras of kinesin-6 and kinesin-14 reveal head-
- headline: In-vivo fission-yeast spindle/kinesin genetics paper — no equations, force laws, or continuum params; sole DCM carry-over is a magnitude-free hook toward a load-bearing nuclear-envelope shell.
  - Nuclear-envelope resistive tension + buckling (qualitative only): motivates a future load-bearing nucleus sub-mesh with its own K_A/kappa_m/enclosed-volume instead of a lumped elastic inclusion
  - No adoptable spindle/MT/kinesin mechanism — that physics lives below DCM's continuum surface (FF motor layer, not DCM)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p06 [1] Taeyoon Kim · review — Contents lists available at ScienceDirect
- headline: Plant cell-wall FE review (Kim co-author) — no equations/params/targets transferable to MCF7 DCM; only a turgor-shell relax-remesh-grow loop concept is worth filing for the division open problem.
  - Relax→remesh→grow→reset-material loop (E2) as a magic-number-free surface-growth operator for the division/volume-doubling open problem
  - Confirms DCM's γ = ΔP·R/2 as the correct animal-cell counterpart to thin-wall σ ≈ ΠR/t (nothing new)
  - REJECT: cellulose fiber-reinforced-matrix anisotropy on the MCF7 cortex (non-physiological wall)
  - REJECT: plant cortical-MT self-organization (FF-scale at best, and plant-specific)
  - Remesh energy-conservation sanity gate: reset per-face K_A rest-area + per-node γ setpoint, no work injected across remesh
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=0

## p72 [0.5] Makito Miyazaki · research — Protein design of two-component tubular
- headline: Protein-design nanotube paper with zero cell-surface physics — nothing to port into the DCM engine at our scale.
  - No DCM-applicable mechanism — de novo protein-tube reconstitution, no cell-scale surface physics
  - No membrane/cortex/turgor/junction/division/confinement content
  - WLC persistence-length (Eq.1) is single-filament FF cross-check, not DCM continuum
  - Coiled-coil binding (Kd 268nM, Tm 35C) is not a cadherin catch-bond junction analogue
  - Recommend no port and no DCM KB mechanism registration
  · n_parameters=0
  · n_equations=1
  · n_validation_targets=0

## p33 [0.5] Makito Miyazaki · experimental — Accurate polarity control and parallel alignment of actin
- headline: In-vitro isopolar actin-array nanofabrication paper with zero DCM overlap — no cell, surface, volume, equation, or cortical-tension datum to port; import nothing.
  - NONE — in-vitro actin-array nanofab, no cell/surface/volume mechanics
  - gelsolin capping / streptavidin crosslink / shear-flow alignment are sub-grid to DCM's lumped cortex γ
  - myosin-V is a CARGO motor (0.62 µm/s), NOT tension-generating myosin-II — do not feed γ/contractility
  - no governing equations, potentials, or constitutive laws to port
  - no DCM scenario (no confinement/aspiration/AFM/division/aggregate/durotaxis)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p11 [0.5] Makito Miyazaki · research — RESEARCH ARTICLE |  APRIL 07 2011
- headline: Single-molecule path-inference method (Go-and-Back / Onsager-Machlup); zero surface-mesh cell mechanics — irrelevant to DCM.
  - (none directly applicable) speculative-only: Onsager-Machlup dominant-reaction-path optimizer as an out-of-scope post-hoc analysis of rare DCM shape transitions (T1/bleb/division) — NOT an engine change
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p01 [0.5] ? · ? — RESEARCH ARTICLE |  FEBRUARY 28 2011
- headline: Single-molecule Bayesian inference toy (2-bead Langevin); no cortex/membrane/junction physics - zero mechanism/equation/parameter for DCM, only a methodological calibration guardrail.
  - Differentiable Bayesian calibration of DCM's unmeasured constitutive params (E_nuc, eta, gamma, kappa_m) via Warp autodiff instead of WKB - tooling only, literature cross-check not PI-data fitting
  - Loss-of-precision / non-identifiability guardrail: stiff weakly-coupled inner compartments (nucleus) are unidentifiable from outer-membrane motion; flag, don't report spurious precision
  - FDT noise convention 2*gamma*kBT + Stratonovich OM discretization - reference for FF thermal Langevin, NOT a DCM force
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p20 [0] ? · ? — PAPERS |  SEPTEMBER 01 2022
- headline: Elementary-school Brownian-motion education workshop; zero mechanistic content for DCM — only the textbook MSD=4Dt diffusion baseline, no cell-mechanics parameters or targets.
  - None — outreach/education paper, no DCM-applicable mechanism
  - Only content: textbook 2D Einstein diffusion MSD=4Dt (Eq.1), a thermostat FDT baseline DCM already satisfies
  - No cortex γ / turgor / membrane / nucleus / junction / division / confinement physics
  - No D, η, or motor parameters reported — nothing to calibrate
  - Corpus provenance value only (Miyazaki-lab outreach item)
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=0

## p12 [0] Makito Miyazaki · review — Honmachi, Sakyo-ku, Kyoto 606-8501, Japan. ORCID iD: https://orcid.org
- headline: 3-page symposium editorial (Miyazaki & Kosugi 2022); no equations, parameters, geometry, or observables — zero DCM-applicable content.
  - none — editorial with no mechanism, equation, or parameter
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0