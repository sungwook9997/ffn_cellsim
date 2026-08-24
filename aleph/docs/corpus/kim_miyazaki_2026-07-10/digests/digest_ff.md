# FF digest (all 77 papers, sorted by score)


## p63 [5] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:345–355
- headline: Sibling FF-class BD engine: crosslinker-clutch friction (not motor force) is the missing lever for sustained actomyosin tension — direct mechanistic route into FF's γ-floor plus a full oracle-grade contraction validation suite (F_max∝R_M^0.6, χ 0.005→0.15, S≈1 at R_M=0.025/R_ACP=0.032, buckling symmetry-break).
  - Crosslinker-as-molecular-clutch: ACP density/off-rate ratio (not motor force) governs SUSTAINED tension — direct γ-floor lever
  - F-actin buckling between crosslink points = emergent tensile/compressive symmetry-breaker (κ_b sensitivity test, Eq.7)
  - No artificial barbed-end stall: motors slide off → turnover-limited processivity, f_M<1pN unless crosslinker friction holds filament
  - Bell slip crosslinker off-rate with compression cutoff (Eq.6): k_u=k0·exp(λ|F|/kT), r<r0→k0
  - Catch vs ideal vs slip bond sets stability/density tradeoff; catch stabilizes at ~150x lower crosslinker density via reverse-sorting
  - Biphasic rigidity sensing: tension ∝ E below ~300Pa knee, flat above (stiff substrate lets motors reach stall fast)
  - PCM motor outputs as Hill-law targets: v0~550nm/s, stall~4.5pN, catch-bond softening of k_w/k_u under load
  · n_parameters=21
  · n_equations=7
  · n_validation_targets=17

## p44 [5] Taeyoon Kim · research — Morphological Transformation and Force
- headline: Kim-lab agent-based actomyosin model is FF's exact physics twin: ports Bell filamin off-rate + PCM catch-bond motor + buckling-required contraction, and hands FF a crosslinker-density sustainability law as a mechanistic lever for its open cortical-tension/contraction problem.
  - Crosslinker load-sharing → tension sustainability law (Bell slip-bond feedback): max-tension ∝ R_M·R_ACP, sustainability ∝ R_ACP/R_M — the mechanistic lever for FF's open γ/contraction problem
  - Filament buckling REQUIRED for contractile tension: 100×κ_b suppresses buckling → tension collapses to ~0 (validates FF bending-beam as necessary, not gold-plating)
  - PCM catch-bond NM-II motor: k_w and k_u both decrease with load; v_unl≈140 nm/s, f_stall≈5.7 pN, k_u,M=0.049 s⁻¹, 8 arms×8 heads — load-dependent detachment, not lumped off-rate
  - Bell ACP slip-bond off-rate k_u=k⁰·exp(λ|F|/kT), k⁰=0.115 s⁻¹, λ=1.04e-10 m (filamin-A, Ferrer 2008) — port as-is with tension-only gate
  - Depolymerization-protection turnover: k_d=k⁰_d(1−ξ_d), sweet spot ξ_d≈0.6 extends motor runway → sustains tension (state-coupled, not lumped decay)
  - Polarity-sorting as a second bundling channel distinct from buckling: bundles form but tension ≈0/compressive when buckling suppressed
  - Cylinder per-segment drag ζ=3πμr_c(3+2r₀/r_c)/5 (use physiological μ, not in-vitro)
  - Buckling diagnostic: filament 'buckled' when end-to-end/contour < 0.6
  · n_parameters=22
  · n_equations=11
  · n_validation_targets=16

## p74 [4.5] Taeyoon Kim · research — 1Living Matter Department, AMOLF, Amsterdam, The Netherlands. 2Institu
- headline: FF lens — p74  (Mulla, Jung, Kim, Tans, Koenderink 2022, "Weak catch bonds make strong networks", Nat Mater 21:1019)

## p69 [4.5] Taeyoon Kim · research — Comp. Part. Mech. (2015) 2:317–327
- headline: Kim-lab FF-class actomyosin model hands FF a directly-portable crosslinker/myosin force-law + parameter set and the key mechanism for our γ program: transient cross-linker friction (off-rate), not motor count, decides whether motor tension is sustained or dissipated.
  - Crosslinker as transient molecular clutch: force-accelerated Bell-slip off-rate sets tension SUSTAINABILITY (not motor count) — direct lever for our γ-floor/stress-relaxation problem
  - Myosin minifilament acts as TEMPORARY crosslinker during fast σ_max loading; passive ACPs govern slow holding — two-timescale split
  - Catch-bond myosin arm kinetics: k_w and k_u both DECREASE with load, letting the minifilament dwell and reach stall (add alongside Hill)
  - ACP Bell unbinding with compression gate: k_u=k0·exp(λ|F|/kT) only when arm stretched (r≥r0); F_b=kT/λ≈39.8pN, k0=0.115/s filamin-A
  - Instrument emergent remodeling (motor aggregation + mesh-coarsening decay-constant) as quantitative readouts distinguishing mechanistic relaxation from structural-instability artifact
  - Cortex-slab stress-sensor scenario: clamp filaments to elastic boundary E=3e4 Pa, read σ(t)=Eε, compute σ_max and sustainability S
  - Frequency-sweep rheology harness (5% sinusoidal strain, 0.1-1000 Hz) → E'(ω),E''(ω), critical frequency
  · n_parameters=24
  · n_equations=6
  · n_validation_targets=18

## p67 [4.5] Taeyoon Kim · research — Cytoskeletal Deformation at High Strains and the Role of Cross-link
- headline: Filamin-A single-molecule spectroscopy: adds an unfolding compliance channel, an angle/torsion rupture law, a two-barrier off-rate, and directly-usable Bell/HS crosslinker constants plus the parent BD bead-drag rheology scenario for FF.
  - Crosslinker Ig-domain UNFOLDING channel (28nm/domain, 57pN)
  - Angle-dependent shear/torsion rupture (~45deg filamin / ~70deg alpha-actinin)
  - Two-barrier ln(loading-rate) off-rate landscape
  - Filamin-A Bell/HS constants (k_off0=0.115s-1, x=0.416nm)
  - Bead-drag network microrheology scenario (Kim BD protocol)
  - Filamin link stiffness kappa_m 266-820 pN/nm
  · n_parameters=22
  · n_equations=5
  · n_validation_targets=18

## p61 [4.5] Taeyoon Kim · research — Disordered actomyosin networks are sufﬁcient to
- headline: Kim-lab agent-based disordered-actomyosin model (FF's own lineage) supplies bipolar-myosin geometry, spatially-gated processivity, and a strong cortex-contraction validation battery (n_H=11, ρc=0.56 µm⁻², telescopic law, bending-invariance).
  - Myosin bipolar minifilament geometry: hinged 3-seg×42nm backbone, 4 arms, 7nm binding lattice (Craig-Megerman/Kim) — pins FF's Stam-Hocky myosin discretization
  - Spatially-gated myosin off-rate (low/processive inside ROI, high/non-processive outside) → mechanistic local cortical-flow & telescopic-patch FF scenario
  - Cooperativity as connectivity percolation: emergent Hill ε_max(ρ) with n_H=11, ρc=0.56 µm⁻², but LINEAR dε/dt(ρ) — dual FF acceptance oracle
  - Telescopic law: boundary velocity ↑ with activation size ξ at ~constant mean strain/strain-rate — cortex-contraction validation target
  - Crosslinker (filamin) over-stiffening SUPPRESSES contraction (ε 2.4→0.3): contraction ∝ σ_a/E, non-monotonic connectivity lesson for α-actinin/filamin
  - Bending/buckling is a RATE modulator, not the mechanism: ℓ_p 20µm→2.6mm slows but preserves contraction — falsification guardrail for FF's buckling attribution
  - Strain=∇·u, strain-rate=∇·v with P2-slope/P3-plateau protocol — port verbatim as FF cortex-flow analysis pipeline
  · n_parameters=14
  · n_equations=5
  · n_validation_targets=16

## p60 [4.5] Taeyoon Kim · research — Interplay of active processes modulates tension
- headline: Same-lineage (Taeyoon-Kim) Brownian-dynamics actomyosin network: adds the actin-turnover stabilizer that keeps FF's motor-driven cortex homogeneous with sustained prestress instead of collapsing into foci — with clean emergent stress/turnover phase-map and live-cell τ_C=28s/τ_D=155s validation targets.
  - Actin treadmilling turnover (barbed-polym = pointed-depolym) with spontaneous ACP/motor release on segment removal — the decisive homogeneous-vs-foci stabilizer FF's turnover module currently lacks
  - Governing control = ratio (turnover rate)/(motor walking rate); peak-stress-vs-turnover curves collapse to a master curve — single dimensionless cortex-stability knob
  - Turnover bounds sustained prestress from above (motors relax before walking far as filaments disintegrate) — mechanistic connection to FF's cortical-tension γ-floor
  - Head count as motor lever: 64 heads/TF default vs 256 for fast aggregation; adopt PCM load-sharing-across-bound-heads unbinding while keeping Hill f-v
  - Pulsed turnover across phase boundary → reversible actomyosin foci with ~100 s transient stress spikes — new FF pulsed-constriction scenario
  - Permanent-bond (fascin/scruin) crosslinker variant + cofilin/buckling severing channel — dual role of turnover
  - Bell off-rate k_off=k_off0 exp(F x_b/kBT) anchored to Ferrer 2008 — confirms FF's existing crosslinker law and 0.066 s-1 anchor
  · n_parameters=14
  · n_equations=5
  · n_validation_targets=10

## p57 [4.5] Makito Miyazaki · experimental — Biphasic Effect of Proﬁlin Impacts the Formin mDia1
- headline: mDia1 barbed-end elongation is a Bell-type force sensor (exp(−fd/kBT), d=5.5nm) with ATP/ADP-resolved k_on/k_off/p_c⁰ — the exact force-dependent formin law FF's constant-velocity growth kernel lacks, plus a tension-switchable profilin biphasic gate.
  - Formin mDia1: replace force-blind barbed-end ratchet with force-dependent elongation law v_elong = k_on·C_A·(1−p_c) − k_off·p_c, p_c = p_c⁰·exp(−fd/kBT), d=5.5nm — tension ACCELERATES (mDia1≠Bni1p)
  - Force-dependent critical-concentration drop: C_crit(f) falls with tension (strongly for ADP) — emergent, not hard-coded
  - Nucleotide-resolved kinetics: ATP k_on=21 k_off=1.3 p_c⁰=0.29 vs ADP 3.6/3.9/0.61 — couples filament aging to production rate
  - Profilin as explicit force-coupled biphasic gate: excess profilin depolymerizes, ~4pN switches −1.8→+40.2 subs/s (all-or-none, optional channel off by default)
  - Processive formin cap ≈ permanent below 6pN (~1/20000 detach) — no fast Bell off-rate on the cap; detachment channel only above 6pN
  - mDia1 gating factor p_o=1 (supersedes assumed 0.5 for ADP) — stepping-second collapses to stair-stepping Eq.1
  · n_parameters=8
  · n_equations=4
  · n_validation_targets=13

## p53 [4.5] Taeyoon Kim · research — Rapid assembly of a polar network architecture
- headline: A Cytosim-sibling agent-based actomyosin engine with a complete parameter table: adopt formin barbed-end-pool depletion + emergent length law + Bell ACP kinetics, and reproduce the decisive result that barbed-out network polarity (not just density) sets contractile force transmission — a new lever on FF's gamma-floor.
  - Polarity as force lever: barbed-ends-out aster transmits ~2x more summed tension (0.76 vs 0.37 nN) than pointed-out despite less local contraction — new handle on gamma-floor/stress transmission
  - Finite depletable barbed-end pool (t1/2~3s, formins in excess) making formin elongation self-limiting and competitive instead of imposed
  - Emergent filament-length law L=v_elong/(1/tau_actin+1/tau_formin) ~6um — derivable exponential distribution, kills a magic-number
  - ACP Bell's-law off-rate k_u=0.115*exp(F*1.04e-10/kBT), F_b~40pN (Ferrer 2008) — reconcile vs KB 0.066
  - Pulsed RhoA activation protocol: local barbed-end boost rho_f=10 tau_f=10s, then myosin d_M=5s tau_M=15s, period 30s — mechanistic contraction driver
  - Formin barbed-end speed 1.1-1.3 um/s (~450 mon/s), k_off 0.11/s — physiological v0 candidate (vs FF 0.623)
  - Myosin minifilament stoichiometry Na=8 arms x Nh=4 heads (32 heads) with parallel-cluster-model kinetics as Hill cross-check
  - Formin-tip persistent-random-walk seeding (P_L=15um) to build a polar barbed-end orientation field on the cortex manifold
  · n_parameters=32
  · n_equations=6
  · n_validation_targets=15

## p48 [4.5] Taeyoon Kim · research — PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1013
- headline: Taeyoon-Kim lamellipodia BD model (FF's exact paradigm) delivers a drop-in angle-dependent buckling-severing turnover channel that unjams myosin-compacted networks, plus a self-consistent parameter set and emergent flow/traction validation targets.
  - Angle-dependent buckling-induced severing k_sev=k0_sev·exp(λ_sev·θ) — NEW turnover channel that unjams myosin-compacted bundles so contraction/flow can proceed (Eq S9)
  - Front-assembly / rear-disassembly balance + closed monomer pool (≥10% floor) → sustained retrograde flow instead of one-shot protrusion
  - Contractile-force-vs-FA-friction competition with force-limited (Bell slip) traction plateau — reframes FF's grid-drag crawl-speed problem
  - PCM/Hill motor cross-check anchors: NMII 120 nm/s unloaded, 5.33 pN stall, catch-like load dependence at k20=17 s⁻¹
  - ACP cross-link density as transverse percolation control on whole-network contraction (threshold R_ACP≥0.02)
  - Nascent-FA transient elastic clutch, Bell slip-off under tension only (Eq S10), 200 nm capture
  - Motor unbinding as active disassembly promoter (buckling→severing coupling), not just force generator
  · n_parameters=30
  · n_equations=6
  · n_validation_targets=13

## p43 [4.5] Taeyoon Kim · research — Computational Analysis of Viscoelastic Properties of
- headline: Direct methodological ancestor of FF: hands the FF crosslinker its missing angular-bending mechanics, an η/geometry-derived drag + excluded-volume law, two turnkey G'/G'' rheology protocols, and a rich table of network-rheology validation targets — with the caveat that k_s is 1/40-scaled and permanent crosslinks make it a passive-backbone unit test, not a physiological-cortex target.
  - Crosslinker as two-armed element with THREE elastic modes: arm stretch k_s,ACP + arm-arm bending k_b,ACP,1 + arm-filament-axis bending k_b,ACP,2 (θ2=π/2) — FF currently lumps this into one scalar link stiffness
  - Bundler vs orthogonal crosslinker as equilibrium arm-arm angle θ1,eq: π (bundle, filopodium) vs 0.417π (orthogonal mesh, cortex) — orthogonal percolates/gels far better, sets cortex crosslinker baseline
  - Oscillatory bulk-shear rheometry harness (Eq 9-10): prestrain + 5% superposed sinusoid → G'/G'' — gives FF its first quantitative viscoelastic-modulus test
  - Segment-tracking MSD microrheology (modified Mason, Eq 11) as an independent second modulus estimator to cross-check bulk shear
  - Anisotropic cylinder drag closed-form (Eq 6, ζ_perp/ζ_par=1.64) and rod-min-distance one-sided-harmonic excluded volume (Eq 4-5) with lever-arm force split for mid-segment binding
  - Prestress-stiffening three-regime emergence (low-γ filament bending, mid-γ crosslinker bending, high-γ stretch) as emergent readout not fitted law
  - Supportive-framework diagnostic: rank crosslinkers by bending force, top 25% carry 70% stress — FF cortex load-path check
  - Bell-Evans crosslinker unbinding is the paper's OWN missing physics (permanent crosslinks → >100Pa stress catastrophe vs experimental 1-30Pa) — FF's differentiator to deliver stress relaxation
  · n_parameters=22
  · n_equations=11
  · n_validation_targets=20

## p42 [4.5] Taeyoon Kim · research — F‑Actin Fragmentation Induces Distinct Mechanisms of Stress
- headline: Canonical Taeyoon-Kim agent model (FF's own family) supplies the missing FF mechanism: buckling-curvature-gated stochastic F-actin severing as an orthogonal, dominant source of cyclic stress relaxation — portable as a fine-grained turnover channel with a like-for-like rheology validation scenario.
  - filament severing (buckling-curvature-gated Poisson)
  - Euler buckling gate F_buckle∝κ/L_c²
  - severed ends feed formin/Arp2/3 growth (turnover balance)
  - Bell ACP unbind confirm (1.1× at F_buckle)
  - orthogonal severing~135° vs unbinding~45°
  - cyclic stress-drop rheology validation
  - reconstituted 3µm PBC shear scenario
  · n_parameters=14
  · n_equations=4
  · n_validation_targets=11

## p41 [4.5] Taeyoon Kim · research — Ding, Chou et al. eLife 2025;14:RP105236. DOI: https://doi.org/10.7554
- headline: Kim-lab agent-based NMII paper is a near-direct FF sibling: supplies explicit bipolar-minifilament geometry (bare zone, arm count/spacing), PCM catch-bond kinetics finer than Hill, an emergent ACP-partitioned contractile-unit mechanism FF should reproduce as a bug-check, and the bundle(N_M)-vs-network(sqrt N_M) force-scaling law that directly targets FF's cortical-tension gamma-floor open problem.
  - Myosin minifilament: make bare-zone L_bz + arm-count N_a + arm-spacing L_sp explicit DOF (L_c=2L_sp(N_a/4-1) controls force 2-3x)
  - PCM catch-bond walking+unbinding (Erdmann-Schwarz) as finer-grained runtime than Hill closed-form; both rates decrease under load
  - Emergent ACP-partitioned contractile units: no-ACP-between-motors adds forces (2x), ACP-between buckles actin and cancels (1x) - a free acceptance test FF should already pass
  - Two-spring crossbridge arm attachment: transverse r0=13.5nm + longitudinal r0=0 (right-angle)
  - Bundle-vs-network myosin-density scaling: F_tot proportional N_M (bundle) vs N_M^0.65 = sqrt(N_M) (isotropic slab) - direct lever on cortical-tension gamma-floor open problem
  - ACP slip-vs-catch sets contractile-unit length -> validates FF Bell-Evans crosslinker force response
  - Actin turnover as a tension regulator (competes with motors, unstalls them), not just a length channel
  · n_parameters=22
  · n_equations=15
  · n_validation_targets=9

## p36 [4.5] Taeyoon Kim · research — 1548 | Soft Matter, 2020, 16, 1548--1559
- headline: Kim-lab sibling of the FF engine: hands FF a load-coupled PCM motor kinetic scheme (unloaded 140 nm/s, biphasic-k20, stall∝1/k20), a Bell-law ACP off-rate cross-check, a soft-harmonic excluded-volume morphogen, and a percolation-gated / buckling-limited contractility law directly targeting the cortical-tension γ problem — with exact stiffness/timestep/PCM-rate numbers gated behind an ESI harvest.
  - Myosin: adopt Parallel-Cluster-Model (Erdmann-Schwarz) load-coupled walk+unbind with stall∝1/k20 as the mechanistic core / oracle for Stam-Hocky+Hill
  - Crosslinker: percolation-gated contractility — max tension at INTERMEDIATE density (~3% percolation) AND intermediate Bell off-rate, non-monotonic (informs γ generation)
  - Contraction driver = filament buckling rectification; buckling SUPPRESSED at very high crosslink density → candidate γ ceiling mechanism
  - Excluded volume as first-class morphogen: soft-harmonic Eq.6 repulsion self-organizes gliding filaments into bundles↔rings by strength κ_r and C_A
  - Bell's-law ACP off-rate k_u=k0·exp(Fλ/kBT) — cross-check FF's existing k_off0/F_b (do NOT overwrite Ferrer α-actinin 0.066)
  - One-arm motor = N_h=8 cooperative heads: sanctioned processivity/duty/stall coarse-graining keeping the cross-bridge cycle
  - Diagnostics to port: velocity autocorr ⟨cosφ⟩, collective ⟨cosα(r)⟩, heterogeneity Q_A grid, spring-energy-density map, tangent-tangent effective ℓ_p
  · n_parameters=17
  · n_equations=6
  · n_validation_targets=15

## p35 [4.5] Taeyoon Kim · research — Soft Matter, 2017, 13, 3213--3220 | 3213
- headline: Supplies the one mechanism FF lacks — buckling-induced curvature-activated F-actin severing (Kim lineage) — turning FF filaments from unbreakable beams into a real fragmentation/turnover + stress-relaxation channel, with a full contraction phase diagram as parity targets.
  - severing
  - fragmentation
  - contraction
  - buckling
  - crosslinker-rebind
  - polarity-sort
  · n_parameters=14
  · n_equations=2
  · n_validation_targets=12

## p27 [4.5] Taeyoon Kim · research — Cytoskeleton, 2026; 83:407–426
- headline: Taeyoon-Kim-lineage agent-based actin model donates a complete, mutually-calibrated stochastic turnover kit (angle-severing law, profilin monomer-pool gating, ATP-aging, capping, 4-angle Arp2/3, formin dwell) with rate constants for FF's turnover/nucleation/severing layer — plant-context, mechanism-donor not MCF7 gate.
  - Angle-dependent severing k_sev=k0,sev·exp(λ_sev·θ) — curvature-gated fragmentation, the natural motor/crosslinker-buckling coupling FF lacks
  - Profilin-gated finite monomer pool (bound/free interconversion) making cortical density self-limiting via pool depletion
  - ATP→ADP aging (t_ATP=3s) gating pointed-end depoly + severing to reproduce non-treadmilling barbed/pointed asymmetry
  - Barbed-end capping/uncapping (k+,cap/k−,cap, ~1000× asymmetric) with severing→immediate-cap coupling
  - Full 4-angle Arp2/3 branch geometry (mother 90°, daughter 20°, branch 70°, coplanar torsion) with asymmetric mother/daughter off-rates — supersedes single-angle stub
  - Processive formin as stochastic dwell state (t_form=5s, ×3, profilin-bound monomers) — loose-coupled elongation
  - One-sided harmonic excluded-volume U_r=½κ_r(r₁₂−r_c)² as Kim-lineage alternative to LJ
  · n_parameters=34
  · n_equations=9
  · n_validation_targets=12

## p26 [4.5] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: Same-family Kim-lineage actomyosin BD model + in-vitro validation delivers a directly portable tension-severing channel (F>500 pN) that attacks FF's stress-relaxation and cortical-tension-floor gaps at high connectivity, plus a rich biphasic contraction/severing-frequency validation suite.
  - tension-severing
  - fragmentation-vs-unbinding
  - connectivity-regime-map
  - k20-ATP-proxy
  - catch-bond-PCM-myosin
  - gelsolin-trigger
  · n_parameters=12
  · n_equations=6
  · n_validation_targets=15

## p24 [4.5] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: Kim-lab fine-grained actomyosin BD model (FF's direct sibling): F-actin treadmilling turnover is the mechanistic stress-relaxation channel that keeps the cortex ergodic — the missing lever for FF's gamma-floor/contraction/stall problems — with portable Bell-law crosslinker off-rate, PCM 3.8 pN/head motor kinetics, and quantitative turnover thresholds.
  - F-actin treadmilling turnover (k_-,A locked to k_+,A) as explicit stress-relaxation channel — releases ALL bound ACP/motor on a removed segment; more efficient than crosslinker unbinding
  - Bell's-law ACP slip-bond off-rate (Eq 1) gated on r>=r_0 — only stretched crosslinker arms feel force-accelerated unbinding
  - PCM force-dependent myosin walk+unbind with 3.8 pN/head stall + explicit slide-off, hopping, free-diffusion states
  - One-arm vs two-arm motor: emergent tension set by antiparallel geometry each half engages (near-zero on parallel filaments)
  - Turnover-rescue threshold k_t,A~=20 s-1 restores superdiffusive/ergodic motor motion — quantitative FF validation
  - Motor-driven aster formation + polarity sorting + buckling as confinement/condensation symmetry-breaking target
  - MSD instrumentation (TE/E/tau-MSD, alpha exponent, g(r) aggregation metric) as FF post-processing observables
  · n_parameters=13
  · n_equations=5
  · n_validation_targets=16

## p15 [4.5] Taeyoon Kim · research — Dissecting Molecular Origins of the Mechano-Adaptive
- headline: Canonical Taeyoon-Kim agent-based actomyosin model at the single-stress-fiber scale — a near drop-in, provenance-tagged parameter set + buckling/recovery validation suite for FF's myosin (PCM catch-bond), crosslinker (Bell ACP), and filament-bending modules.
  - Myosin PCM catch-bond force-velocity kernel (per-head Markov, v0=56nm/s, F_s=24pN, 30pN slip cutoff) as mechanistic alternative to Hill
  - Load-free motor walking as a turnover-independent contraction/stress-relaxation mode (recovery rate set by k20, not density)
  - ACP Bell slip-bond confirmed: k0_u=0.05/s, x_u=0.9nm -> F_b~4.6pN; stretch-induced unbinding = permanent lengthening channel
  - Buckling threshold = compression rate exceeds unloaded motor walking speed (emergent symmetry-breaking rule)
  - Density-vs-speed separation invariant: F_B,ss proportional to R_M; recovery time/rate proportional to k20 only
  - Single anti-parallel contractile-unit FF scenario (N_F=50, stretch/hold/compress/recover protocol) as native-mechanistic dev regression
  - Filament bending sets buckling shape but not recovery time (decoupling invariant); do NOT copy their non-physical 10x kb inflation
  · n_parameters=30
  · n_equations=9
  · n_validation_targets=24

## p14 [4.5] ? · ? — Reconstitution of actomyosin networks in cell-sized
- headline: Taeyoon Kim's own Brownian-dynamics actomyosin model is FF's exact paradigm — supplies a NEW force-triggered severing channel, a parallel-cluster motor kernel, cylindrical drag law, and a fully-specified native-scale parameter set with a ~200 pN/µm cortical-tension validation target.
  - Force-triggered F-actin severing channel (F_sev=300 pN, Tsuda 1996) — NEW turnover/stress-relaxation capability FF lacks
  - Erdmann-Schwarz parallel-cluster NM-II motor: emergent f-v (v0≈140 nm/s, stall 5.7 pN/arm) from per-head 3-state kinetics vs FF's imposed Hill
  - Underhill-Doyle cylindrical per-segment drag ζ=3πμd(3+2l/d)/5 (S3) — principled form for FF grid-drag/crawl-speed problem
  - Bell-kinetics cortex-membrane detachment oracle γ_C^D* ∝ R_C (Eqs 4-8) as coupled-run validation gate
  - α-actinin rupture-force spectrum (1.4-80 pN) to anchor Bell-Evans crosslinker F_b + severing-dominates-unbinding sanity gate
  - Confined-shell assemble-then-gate-myosin bleb scenario with independent R_C/R_X control and cue-free bleb positioning (symmetry-breaking)
  - Cortical network tension plateau ~200 pN/µm + ∝R_C scaling as γ-floor datum
  · n_parameters=28
  · n_equations=15
  · n_validation_targets=14

## p10 [4.5] ? · ? — Fine-tuning of material properties by catch bonds
- headline: Directly portable catch-slip α-actinin off-rate + a novel crosslinker-turnover/redistribution mechanism that raises both yield stress and strain and amplifies myosin contraction ~3–4× — FF's exact native regime from a Kim-lab (Cytosim-sibling) model.
  - Catch-slip crosslinker off-rate (Eq.1: slip Bell term + decaying catch term) supersedes FF's slip-only α-actinin — breaks stiffness↔extensibility trade-off
  - Crosslinker TURNOVER: unbind→hop (5µm)→re-bind — NEW FF capability (absent in AFINES/MEDYAN/ALENS/Cytosim); necessary for the catch-slip advantage (no-turnover control flips it negative)
  - Force-dependent CL redistribution amplifies myosin contraction ~3–4× (~90 vs ~25 Pa) — CLs arrive at motor before filament aggregation
  - Store {k_s0,k_c0,λ_s,λ_c} primitives, derive {τ0=9.5s, τmax=53s, Fmax=14.7pN}; slip-only recovered by k_c0=0 (back-compat, K255E control)
  - Connectivity gate: catch-slip beats slip only at high CF/RCL (CF=10µM, RCL=0.01) — validate at dense cortex regime, not loose
  - CL size/stiffness (α-actinin 35nm stiff vs filamin 160–250nm flexible) sets catch→slip transition rate — derive from FF per-CL geometry not magic number
  - γ̇-scaling: τ0,τmax ∝ 1/(strain rate ≈ motor walking rate) — rescale in-vitro lifetimes to FF physiological motor velocity
  · n_parameters=22
  · n_equations=3
  · n_validation_targets=15

## p05 [4.5] Taeyoon Kim · research — Mobility of Molecular Motors Regulates Contractile
- headline: This is FF's own Kim-lab/Cytosim physics lineage: it hands FF an explicit motor-mobility DOF (drag ζ_M) plus one-vs-two-arm connectivity that mechanistically explain cortical contraction/γ, with directly portable Bell/PCM kinetics and 140 nm/s + 8-13 pN validation targets — deep numeric tables are in the SI (follow-up needed).
  - Motor mobility as explicit DOF: per-minifilament drag ζ_M (ref 8.10e-8 kg/s) separating membrane-anchored vs cytosol-free myosin
  - One-arm vs two-arm connectivity: FF bipolar minifilament must sit on the MONOTONIC (mobility↑→contraction↑) branch via antiparallel-pair pulling
  - Single ATP knob k₂₀ (ref 20/s) driving stall∝1/k₂₀, walk & off-rate ∝k₂₀ self-consistently (PCM/catch-bond)
  - PCM cross-bridge scheme (3 states/5 rates, load-suppressed catch-bond) as mechanistic alternative/augment to FF Hill heads
  - Bell's-law ACP off-rate k_off0·exp(Fλ/kT) — confirms FF crosslinker parity (Ferrer 2008 λ)
  - Contraction is drag/elasticity-rate-limited to ~100s (scalars plateau ~20s) — cross-engine timescale corroboration
  - Q_A areal-density heterogeneity metric (N_G=20) as FF cortex clustering/condensation diagnostic
  - v_max=140 nm/s NMII calibration + ~8-13 pN per-arm force validation targets
  · n_parameters=20
  · n_equations=8
  · n_validation_targets=12

## p04 [4.5] Taeyoon Kim · research — Balance between Force Generation and Relaxation
- headline: Sibling Kim-lab actomyosin model hands FF two ready-to-port relaxation kernels — treadmilling turnover and buckling-selective severing — that close its network-contraction/stress-relaxation and fragmentation open problems, plus in-vivo pulse targets (2-4 um, 40-100 s).
  - Segment-resolved treadmilling turnover (barbed-polymerize = pointed-depolymerize) as explicit stress-relaxation channel — the missing dissipation behind FF's contraction/relaxation open problem
  - Angle/buckling-dependent F-actin severing k_s=k0·exp(theta/lambda) using FF's existing bending-junction angles — lets irreversible clusters disassemble
  - Both channels REQUIRED together: severing recycles fragments via treadmilling; neither alone sustains physiological pulses
  - Non-monotonic contractility vs connectivity (crosslinker density / filament length), max at INTERMEDIATE — a falsifiable gamma-generation prediction
  - PCM 4-arm x 8-head motor force-velocity + load-dependent unbinding as cross-check oracle for FF Stam-Hocky/Hill minifilament (do not overwrite)
  - Bell's-law ACP off-rate corroborates FF's existing Bell-Evans crosslinker law
  - Rheology diagnostics to port: g(r), heterogeneity Q_A, motor speed (Eq.2), velocity-persistency autocorrelation (Eq.5)
  · n_parameters=14
  · n_equations=5
  · n_validation_targets=12

## p03 [4.5] Taeyoon Kim · research — Determinants of Fluidlike Behavior and Effective Viscosity in Cross-Li
- headline: Canonical Taeyoon-Kim BD cortex-rheology paper = FF's own physics family: hands FF a turnover-driven (not force-driven) cross-link relaxation mechanism plus a clean emergent constitutive law (η_eff∝1/turnover, biphasic yield σ_th≈1Pa, η_eff0≈35 Pa·s, 5%/99% percolation) to reproduce natively as cortex-creep validation.
  - Turnover-driven cross-link detachment as PRIMARY relaxation channel (link dies when host segment depolymerizes + rebind lockout), superseding Bell-only unbinding — 4x lower η_eff, force-dependence sub-dominant below 10 Pa
  - Balanced steady-state treadmilling with single global scale factor n (Eq.5: k+_B=6n, k+_P=0.6n, k-_B=0.6n, k-_P=6n; rate=5.4n/s) tuning turnover independent of length/density
  - Per-segment cylinder Stokes drag ζ=3πη_m·r_c·(3+2r0/r_c)/5 (Eq.4) — port exactly, at physiological η_m
  - Emergent biphasic yield-creep from nonlinear subsegment force-extension: ε̇∝σ below σ_th≈1Pa, ε̇=const above; master law η_eff·n=η_eff0·max(σ,σ_th)/σ_th, η_eff0≈35 Pa·s
  - Nonaffine load-bearing percolation observable: ~5% of filaments carry ~99% of tension; creep = cyclic pathway formation→buildup→turnover→release
  - Force-force autocorrelation R(τ) with decay slope ∝ n as the turnover-timescale probe
  - Bell slip-bond k_ub=k0·exp(λF/kT) as SECONDARY channel (k0=0.115/s, λ=1.04e-10m → F_b≈39.5pN, Ferrer 2008 — already in KB)
  · n_parameters=16
  · n_equations=8
  · n_validation_targets=11

## p02 [4.5] Taeyoon Kim · research — Dynamic Role of Cross-Linking Proteins in Actin Rheology
- headline: Kim-Kamm 2011 (FF crosslinker-layer ancestor) hands FF a complete literature-anchored Bell-Evans unbinding + WLC-arm crosslinker model and proves UNBINDING (not domain unfolding) governs cortex rheology at physiological strain rates -- validating FF's off-rate design and giving 4 no-fit native-scale acceptance targets (K~gamma^3, ~1 Pa relaxation plateau, plastic exponent 0.17, 25%/85% load framework).
  - Bell-Evans crosslinker UNBINDING is the validated rheology driver: k_ub=k0_ub*exp(lambda_ub*|F|/kT), k0_ub=0.115/s (filamin, Ferrer2008), lambda_ub=1.04e-10 m -> F_b~=39.5 pN
  - Domain UNFOLDING negligible at physiological gamma_dot_eff~=0.1/s -- FF can defer the multi-step unfolding module (unfolding dominance only above ~1/s)
  - Rebinding must be ON at physiological rate: sets finite ~1 Pa stress-relaxation plateau; without it stress->0 (unphysical over-relaxation)
  - Emergent supportive-framework reorganization: ~25% of network bears ~85% of >20 pN ruptures, rebinding to a DIFFERENT filament above 20 pN
  - WLC crosslinker-arm force law (Eq.1) upgrades harmonic bonds: p=0.33 nm arm, r0=105 nm filamin, linear compression below r0
  - Bulk-shear-box validation scenarios (2.8 um, C_A=12.1 uM, R=0.02): strain-stiffening K~gamma^3, plastic exponent ~0.17, stress cap ~100 Pa
  - Register filamin as a DISTINCT crosslinker species -- do NOT overwrite FF alpha-actinin k_off0=0.066/s (same Ferrer2008 paper, different isoform)
  · n_parameters=17
  · n_equations=10
  · n_validation_targets=19

## p65 [4] Taeyoon Kim · research — Computational Analysis of a Cross-linked
- headline: Ancestral blueprint of FF: supplies the missing filament-torsion+helical-registry mechanism and two ACP archetype geometries that make bundle-vs-mesh crosslink angle (−2.4° vs 87.5°) emerge, plus rich morphology acceptance targets — but its kinetics are accelerated and stretch deliberately soft, so port forms/geometry, reject the rate constants.
  - Filament torsion + helical binding-site registry (κ_t=4.14e-18 N·m, 36°/bead, 10 beads/turn) — the missing potential that makes cross-link angle emergent
  - Two-gate crosslinker binding (attach→cross-link) with ±10° angular acceptance windows + dihedral φ gate
  - Two ACP archetypes as pure geometry parameter sets: bundler ACP^B (θ₁=π,φ₁=0,stiff κ_t→0° parallel) vs orthogonal ACP^C (φ₁=π/2,soft→90° mesh)
  - Bell's-law force-dependent unbinding k_off=k0·exp(γF/kT) — corroborates FF slip-bond form (keep FF's sourced k0, reject Kim's accelerated 0.4967 s⁻¹)
  - Emergent late-time crosslink remodeling (N_active ~ t^-0.05) as a mechanistic stress-relaxation/fluidization channel
  - L_p↔κ_b calibration via tangent correlation e^{-s/Lp}, L_p≈17-20 µm cross-check
  - Damköhler/Φ scaling laws: ⟨L_f⟩~k_n^{-1/2}, pore~C_A^{-1/3} as predictive checks on FF polymerization
  · n_parameters=18
  · n_equations=9
  · n_validation_targets=13

## p59 [4] Makito Miyazaki · research — Cell-sized spherical conﬁnement induces the
- headline: Cell-sized confinement reconstitution gives FF a parameter-free self-organization benchmark: equatorial contractile-ring emergence, duty² motor scaling, and P*≈filament-length contraction — all at our R<ℓ_p cortex scale.
  - Confinement-driven equatorial ring self-organization as a parameter-free FF emergence test (ℓ_p≫R elastic-rod localization)
  - Two-heads-bound + turnover requirement validates Stam–Hocky bipolar minifilament over single-head/rigor proxies
  - [myosin]_eff = [bulk]×(duty ratio)² — motor-model check via myoV/myoII 100× duty² offset
  - Contraction velocity ∝ (P−P*) with P*≈mean filament length 5.4 µm as emergent contractile unit
  - Constant-volume, depolymerization-free contraction invariant (V=πR·πd²/4) for FF motor-only sliding
  - α-actinin↔myosin bundling-vs-remodeling balance controls ring assembly
  - Cell-sized spherical/quasi-2D confinement scenario at our exact scale (R<ℓ_p) for symmetry-breaking/self-organization
  - Contraction rate anchor 0.39 µm/min per µm perimeter (matches living cytokinetic rings)
  · n_parameters=14
  · n_equations=6
  · n_validation_targets=9

## p21 [4] Taeyoon Kim · research — Dr. V. Yadav, Dr. A. P. Tabatabai, Prof. M. P. Murrell
- headline: Formin nucleation on a finite shared G-actin pool is a single mechanistic knob that sets emergent filament length, turnover, and polymerization-force buckling — replacing FF's imposed length distribution and adding a motor-free symmetry-breaking channel, with a full stack of in-vitro validation numbers.
  - Finite shared G-actin pool conservation (Σ L_i·ρ_sub + G_free = N_p) → emergent length control, replacing FF's imposed exponential length dist
  - Formin two-mode elongation (150→10 subunits/s) + length-triggered release + nucleator recycle
  - Polymerization growth-pressure buckling: stalled barbed-end axial force f_stall~1pN vs Euler f_crit=π²κ/L² → motor-free buckling/templating (new FF channel)
  - Assembly-limited turnover: τ varies 7× at constant disassembly rate purely from pool depletion — τ as derived observable not input
  - Bending-energy readout E_bend=(κ/2)(L/R²) + Frank K33=1.38-4.6pN as defect diagnostics
  - Templating/mechanical-memory: τ_relax/τ_turnover peaks at intermediate nucleation as a native-scale integration test
  · n_parameters=13
  · n_equations=5
  · n_validation_targets=17

## p17 [4] ? · ? — Kinetic Control of Out-Of-Equilibrium
- headline: Kim agent-based actomyosin (same physics class as FF) gives FF a NEW emergent lever — the crosslink-before-contract timing delay that non-monotonically maximizes contractile force (peak ~5 s, ~38 nN) — plus scarce in-vivo cortical myosin turnover rates (k_off 0.03-0.10 s-1) and a pulsed-cortex activation protocol.
  - Actin-crosslink-before-motor DELAY (d_M) as a first-class contractility lever — force non-monotonic, peaks ~5 s (NEW to FF)
  - Disassembly-dominated myosin turnover: modulate minifilament k_off (0.025-0.10 s-1), not just recruitment
  - RhoA-gated local/transient/stochastic subdomain activation protocol -> emergent pulsed cortex (~31 s period)
  - Formin boosted barbed-end elongation: 1.2 um/s, rho_f=10, tau_f=10 s -> ~10-12 um filaments
  - Leaky-integrator turnover dN/dt=k_app-k_off*N with high-pass 1/sqrt(1+w^2/k_off^2) as FF frequency-response check
  - F-actin ~0.5 s substrate lead before myosin loading (sequencing constraint)
  - Kim Bell + parallel-cluster-model network as cross-check oracle for FF Bell-Evans + Stam-Hocky (no downgrade)
  · n_parameters=18
  · n_equations=7
  · n_validation_targets=10

## p13 [4] Taeyoon Kim · research — Swirling motion of breast cancer cells radially aligns collagen fibers
- headline: Kim-lab BD semiflexible-fiber network is a Cytosim-class twin of FF's collagen side: hands FF a full oracle parameter set plus a first-principles negative-normal-stress + annular swirl-shear→radial-alignment (exponents 1.2/2.2) validation scenario at our fine-grained scale.
  - Negative normal stress: shear→contractile radial stress as EMERGENT falsification test of FF semiflexible-network buckling/tension-asymmetry (|σ_normal|≈σ_shear, ~50/50 stretch/compress)
  - Annular dipole-shear → radial fiber alignment (TACS-3) as an FF collagen-ECM self-organization scenario (radial exponent 1.2 vs tangential 2.2)
  - Bell's-law slip-bond crosslinker with tension-only off-rate acceleration (compression→base rate), Eq.9 — mechanistic refinement over |F|-based Bell-Evans
  - Clift cylinder drag ζ=3πμ·r_c(3+2r0/r_c)/5, length-dependent per-segment drag (form portable, μ not)
  - Full Kim-lab BD collagen-fiber parameter set (κ_s,f=4e-3 N/m, κ_b,f=8.27e-20 N·m→ℓ_p~20µm, k0_off=1e-6/s, Bell λ=0.4nm) as FF Mikado-ECM parity oracle
  - Overdamped Langevin + fluctuation-dissipation noise normalization as BAOAB high-friction cross-check
  · n_parameters=24
  · n_equations=9
  · n_validation_targets=10

## p75 [3.5] Taeyoon Kim · research — Nature  |  Vol 626  |  15 February 2024  |  635
- headline: FF lens — p75  (Fan, ..., Taeyoon Kim, ..., Török 2024, "Matrix viscoelasticity promotes liver cancer", Nature 626:635)

## p64 [3.5] Taeyoon Kim · research — Biomech Model Mechanobiol (2015) 14:1143–1155
- headline: Taeyoon-Kim BD crosslinked-actin engine (FF's direct ancestor) - a ready-made strain-stiffening rheology validation scenario + two-arm ACP geometry, not a new active mechanism.
  - Two-arm ACP crosslinker with explicit arm-alignment bending (k_f,ACP1=1.45e-25) + filament-attachment perpendicular spring (k_f,ACP2=5.8e-25 N.m2) - adds geometric bending FF's stretch-only crosslinker lacks
  - Crosslinked-network shear rheology scenario as an FF sub-engine validation gate (3um box, 20uM actin, R_ACP=0.01, dgamma/dt=0.1/s, sever+clamp shear)
  - Soft-element-dominates principle: with compliant alpha-actinin/filamin (k_s,ACP=0.002 N/m << filament), network is crosslinker-limited, insensitive to filament kappa - a prior for the gamma/stress-relaxation problem
  - Per-filament tunable persistence length (->kappa) from nucleotide/cation state, physiological default Mg2+-ADP; sweep 3.5-12.7um
  - lp<->kappa conversion kappa=kBT*lp to set bending rigidity from measured persistence length
  - Harmonic stretch/bend potentials V_s=0.5*k_s(l-l0)^2, V_b=0.5*k_f(theta-theta0)^2 - parity check with FF force inventory
  · n_parameters=17
  · n_equations=5
  · n_validation_targets=7

## p50 [3.5] ? · ? — Molecular Biology of the Cell • 34:ar67, 1–10, June 1, 2023
- headline: Kim-lab Arp2/3 lamellipodium-clutch agent model: adopt tip-gated 300 nm clutch seeding + reproduce emergent PFAC load-sharing crossover; highest-value lever is the cited buckling-fragmentation channel for contraction. Constants are NO_SI.
  - Tip-gated clutch seeding: nucleate clutch when substrate/ECM node is within 300 nm of an F-actin barbed end; one tip seeds multiple parallel Bell/catch springs
  - PFAC — emergent parallel-adhesion load-sharing: surviving clutch ramps force when neighbors fail; time-to-failure crossover (fail<10FA, no-fail@30FA) should EMERGE, not be parameterized
  - Buckling-induced F-actin fragmentation as contraction facilitator (cited Li/Kim 2017) — curvature-triggered severing to unblock FF network contraction & gamma-floor
  - Confirms FF runtime laws: Bell-Evans slip off-rate k_off=k_off0*exp(F/F_b), Arp2/3 70deg branch; keep Hill+catch over paper's linear Chan-Odde+slip (cross-check only)
  - Substrate-stress observable sigma=Sigma f_clutch/area reproducing force-up-while-stress-down as adhesion footprint widens
  - Keep F-actin polymerization ON to avoid the paper's 20 s transient-window limitation (FF is already ahead)
  · n_parameters=12
  · n_equations=8
  · n_validation_targets=12

## p39 [3.5] Taeyoon Kim · research — 10274 |  Soft Matter, 2021, 17, 10274–10285
- headline: Ports a Bell-slip-bond transient cross-link kinetics onto FF's collagen/ECM substrate so matrix viscoelasticity, plasticity and stress-relaxation emerge mechanistically (k_ub,0=1e-5 s⁻¹, x_ub=1e-10 m); ECM-side only — the paper's lumped zero-rest-length cortex is discarded in favor of FF's native actomyosin.
  - ECM collagen cross-links as explicit transient Bell slip-bonds: k_ub=k_ub,0·exp(x_ub·|F|/kBT), one-sided (tension-only)
  - Emergent viscoelastic/viscoplastic remodeling & irreversibility from load-biased cross-link turnover (no lumped modulus)
  - Connectivity-gated relaxation: load-sharing (C_f, R_xl) sets effective off-rate → sustained vs fast-relaxing stress
  - ECM cross-link constants k_ub,0=1e-5 s⁻¹, x_ub=1e-10 m (F_b≈41 pN), 100-nm binding-site spacing, ~10 µm fibers (Nam 2021 ref45)
  - Poroelastic guardrail: >50% relaxation persists at zero turnover — do not attribute all relaxation to Bell kinetics
  - Buckling-orientation readout: tensed fibers radial (0°), buckled fibers circumferential (90°)
  - σ∝r⁻¹ (2-D slab) / r⁻² (3-D) stress-transmission diagnostic + relaxation validation ladder
  · n_parameters=9
  · n_equations=3
  · n_validation_targets=9

## p32 [3.5] Taeyoon Kim · review — Weldon School of Biomedical Engineering,
- headline: Kim-lab review is the design brief for FF: its stated gap (no discrete model with turnover + force-dependent crosslinker unbinding + myosin together) IS the FF thesis; delivers network-level validation gates (K' 1->100 Pa vs prestress; G'' dip at crosslinker off-rate) and flags fragmentation as the missing FF channel — but prints no rate constants (chase primaries).
  - Crosslinker off-rate validated by network G'' minimum at k_off0 (Fig 5f) — turns Bell-Evans constant into a falsifiable observable
  - Buckling-induced F-actin fragmentation/severing channel as the root of stress relaxation + viscoplasticity (missing in FF)
  - Myosin-II prestress as the resting operating point (filaments pre-tensioned from t=0) — mechanistic basis for cortex tension gamma
  - Strain-stiffening as emergent bending->extension crossover, not a fitted differential-modulus law
  - Microtubule-facilitated stiffening: rigid struts suppress actin bending, force extension
  - Full catch-slip (not slip-only) crosslinker off-rate for load-bearing alpha-actinin
  - Power-law relaxation must EMERGE from turnover+unbinding timescale spread, never a lumped fractional element
  · n_parameters=7
  · n_equations=0
  · n_validation_targets=8

## p22 [3.5] Taeyoon Kim · research — Cellular Pushing Forces during Mitosis Drive Mitotic
- headline: Kim-lab collagen model is a Cytosim-native FF sub-engine: port the tension-gated Bell off-rate and geometry-derived cylindrical drag (a candidate grid-drag fix), plus a ready G'(density) stress-relaxation validation harness and mg/mL calibration for the FF collagen ECM.
  - Tension-gated Bell's-law crosslinker off-rate (Eq 9): apply exp(F/F_b) only when the bond is in tension, baseline k_0u under compression — mechanistic correctness fix over magnitude-based off-rate
  - Geometry-derived per-element cylindrical drag ζ=3πμr_c(3+2r_0/r_c)/5 (Eq 4) — candidate physical reframing of FF's count-scaling grid-drag / crawl-speed problem, exact for ECM fibers
  - FDT noise rebuilt from the same geometric ζ_i (Eq 3) — consistency constraint if Eq-4 drag adopted
  - Two-arm, two-angle crosslinker with perpendicular (pi/2) fiber-binding geometry — upgrade FF filamin/alpha-actinin from axial-only to orientationally constrained
  - Fiber-stiffness-from-modulus (Eq 10: k_s=piE r_c^2/r_0) and length-to-mg/mL calibration (Eq 11, nu=0.73) for the FF collagen ECM
  - Stress-relaxation / shear-rheology validation harness (50x50x1um, strain-rate 0.1/s, hold) reproducing G'(density) + strain-stiffening + strain-enhanced relaxation
  - Cell-pushes-collagen division scenario (F_Cyto/F_Is=100 pN/side applied loads + constant-volume) as an FF ECM scenario with emergent fiber-buckling observable
  · n_parameters=14
  · n_equations=9
  · n_validation_targets=12

## p19 [3.5] Taeyoon Kim · research — Mechanical Checkpoint for Cell Division in
- headline: Taeyoon-Kim-lab collagen fiber-network + tension-gated Bell crosslinker model: a complete ECM force-law set (Eqs.1-9), a growth-prestress→pericellular-shell scenario, and a quantitative buckling/anchoring validation battery for FF's ECM substrate half.
  - Tension-gated (one-sided) Bell unbinding: k_off=k_off0·exp(λF/kT) only when r≥r0, else k_off0 — compressive branch clamp for FF crosslinkers/clutches
  - Two-arm crosslinker with two bending stiffnesses + proximity-weighted mid-segment force distribution (also for myosin heads/clutches)
  - Cylinder per-element drag ζ=(3πµr_c/5)(3+2r0/r_c) — principled alternative to bulk-drag for FF's grid-drag/crawl-speed problem
  - Soft-core bounded excluded-volume repulsion (Eq.8) for fiber-fiber and membrane-fiber contact vs LJ singularity
  - Growth-prestress protocol: ramp membrane V0 ×2 to emergently build a densified, load-bearing pericellular ECM shell (physiological baseline)
  - Collagen ECM as self-assembled (nucleation+bidirectional elongation) semiflexible fiber network, ℓ_p=100µm, emergent ~30pN buckling
  - Discrete-vs-continuum discriminator: three coexisting resistance modes reproducible only at fiber resolution, not Neo-Hookean FE
  · n_parameters=20
  · n_equations=9
  · n_validation_targets=14

## p07 [3.5] Taeyoon Kim · research — Covalent cross-linking of basement
- headline: Taeyoon-Kim Bell slip-bond crosslinked-ECM model gives FF a tension-gated crosslinker off-rate (k0,u=3e-6 s⁻¹, F_b=41 pN), a covalent=parameter-limit trick, and ready-made matrix-plasticity + protrusion-into-ECM validation scenarios.
  - Bell slip-bond crosslinker off-rate k_u=k0,u·exp(λu·F/kBT), tension-ONLY gating (base rate under compression) — fidelity fix for FF crosslinker unbinding
  - F_b = kBT/λu = 41.1 pN and k_off0 = 3e-6 s⁻¹ as ECM/BM-crosslink Bell anchors (NOT cortex α-actinin)
  - Covalent crosslink = parameter limit (k0,u→0, λu→0) of existing Bell bond — no new module
  - Discrete 140-nm crosslinker binding-site lattice → plastic set emerges from unbind+rebind-in-deformed-state (no lumped viscoplastic term)
  - Creep→recovery matrix-plasticity assay (100 Pa 1h→0 Pa): flat modulus + collapsing permanent strain as k0,u↓ = sharp FF ECM validation gate
  - Protrusion-into-ECM invasion scenario (1 nN, 40 nodes) mapped to FF filopodium-into-collagen crawl line
  - Short-fiber (~4 µm) crosslinked Mikado builder as basement-membrane environment preset
  · n_parameters=13
  · n_equations=3
  · n_validation_targets=13

## p77 [3.0] Makito Miyazaki · research — Spatial confinement of active microtubule networks
- headline: FF lens — p77  (Suzuki, Miyazaki et al. 2017, "Spatial confinement of active MT networks induces rotational cytoplasmic flow", PNAS 114:2922)

## p68 [3] Taeyoon Kim · research — sharing, adaptation, distribution and reproduction in any medium or fo
- headline: Plant cortical-MT kinetic model donates a three-state dynamic-instability rate set, an angle-dependent zippering/catastrophe collision scheme, and a stress->plus-end mechanotransduction feedback (best cast as Bell-Evans catch-bond) for FF's underdeveloped --microtubules compartment.
  - MT three-state plus-end dynamic instability (6 transition freqs + constant minus-end shrink) for FF --microtubules compartment
  - Angle-dependent MT-MT collision: <40 zippering/bundle, >40 catastrophe-or-crossover (P_cat 0.2-0.8) - emergent nematic ordering FF lacks
  - Stress->plus-end-kinetics mechanotransduction feedback: project cortex stress onto MT (Eq1), bias polymerization/catastrophe - loop closes natively in FF
  - Implement the coupling as Bell-Evans exponential (10 pN / 300 nm2 catch-bond sensor), NOT the paper's admitted linear approximation
  - Plant cortical MT DI rate set (Shaw2003): v_grow 3.69, v_shrink+ 5.80, v_shrink- 0.53 um/min; nucleation 10 um^-2 min^-1; 100 nm segment
  - Catch-bond MT-membrane anchor as the directional stress-sensor (tension lowers dissociation) - fits FF's Rakshit catch-bond vocabulary
  - Nematic order parameter S_p (length-weighted cos2theta) as FF MT-order observable
  · n_parameters=18
  · n_equations=3
  · n_validation_targets=15

## p62 [3] Makito Miyazaki · research — Processive Nanostepping of Formin mDia1 Loosely Coupled with
- headline: mDia1 single-molecule nanostepping gives FF a mechanistic two-channel (association/dissociation) discrete 2.7 nm barbed-end growth kernel with direct step-size and dwell-time rate anchors, replacing the current lumped formin drift velocity.
  - formin: discretize barbed-end growth into explicit ~2.7 nm Poisson stepping (replace lumped drift v0)
  - formin: two-channel kinetics — conc-dependent association k_on·[G-actin] vs conc-independent dissociation k_off → emergent critical concentration
  - formin: couple on-rate to FF cytoplasm G-actin reservoir so monomer depletion mechanistically gates growth
  - formin: loose coupling — optional 2x/3x multi-subunit step channel (intermediate states enable CP cooperation)
  - measurement: dwell-time Gamma-mixture (eq1) + merge-artifact probability model (eqs 2-5) as sampling-protocol sanity guard
  - validation: single-formin single-filament ratchet benchmark reproducing step 2.7 nm + k_fwd 0.92/1.4 s-1
  · n_parameters=9
  · n_equations=5
  · n_validation_targets=7

## p58 [3] Makito Miyazaki · research — Myosin-Driven Advection and Actin Reorganization Control the
- headline: Purified confined-actomyosin gives FF emergent validation targets (inward flow v=1.2 um/s, pulsatile T=100s, contractility=size vs polymerization=shape) but its continuum active-fluid model is a lumped cross-check oracle only, not portable FF mechanism.
  - Myosin (Stam-Hocky) emergent inward advective flow benchmarked to v=1.2 um/s
  - Contraction-wave pulsatility at T=100s from coupled myosin+turnover
  - Turnover (k_p source / k_d sink) as controller of contracted steady-state shape/size
  - Surface-biased Arp2/3 nucleation (k_p^surf>k_p^bulk, membrane-templated)
  - Confined actomyosin condensation as native-scale FF self-organization scenario
  - Double dissociation: contractility sets SIZE, polymerization sets SHAPE (perturbation ordering)
  - Corner/flat-side stress concentration in shaped confinement
  · n_parameters=13
  · n_equations=4
  · n_validation_targets=7

## p45 [3] Taeyoon Kim · research — Characterization, Enrichment, and Computational
- headline: Taeyoon Kim's sibling actin+ACP engine hands FF a portable boundary-clamped bulk-rheology harness, an R_ACP=0.1 cross-linker connectivity floor, and a tension≫shear>compression stiffness-anisotropy target — but it is passive/motor-free, coarser than FF, and the real contraction lever lives in its ref 33.
  - Bulk-rheology validation harness: assemble network, sever+clamp filament ends at y-faces, strain in shear/tension/compression at 0.001/s to |ε|=0.05, measure σ from clamped-end reaction forces (virial boundary stress)
  - Cross-linker connectivity floor R_ACP = C_ACP/C_A = 0.1 as a build-time percolation gate for FF cross-linker seeding (guards against under-linked floppy networks)
  - Network stiffness anisotropy target: tension ~2.4e4 > shear ~1.0e4 > compression ~0.5e4 Pa — nonaffine signature testing FF's cross-linker+EV+bending stack
  - CLAN triangular-lattice / hub-spoke architecture as a new emergent FF network class (must self-assemble, not hard-code the lattice)
  - POINTER ref 33 (Li/Kim Soft Matter 2017): buckling-induced F-actin fragmentation channel that unjams motor-driven contraction — the genuine fine-grained lever for FF's contraction/stress-relaxation open problem (pull separately)
  - LINC-mimicking harmonic tether coupling 20% of filament ends to nucleus (~2x nuclear compression stiffening) — mostly DCM but FF supplies the tethered filament ends
  - Keep FF's mechanisms over theirs: BAOAB over forward-Euler, Bell-Evans ACP off-rate over fixed harmonic link, Hill/Stam-Hocky motors over passive network
  · n_parameters=12
  · n_equations=8
  · n_validation_targets=9

## p40 [3] Taeyoon Kim · research — Durability of Aligned Microtubules Dependent
- headline: Calibrates FF's microtubule compartment (Lp→κ, length, kinesin glide/stall) and hands a purely-mechanistic emergent-order regression scenario whose key law — alignment durability ∝ bending rigidity, not binding probability — directly tests FF's bending term + BAOAB thermostat.
  - MT compartment as explicit bending beam: κ=Lp·k_BT from measured Lp 0.65–3.54 mm (2.7e-24…1.4e-23 N·m²)
  - Directed motor-propelling force per segment (substrate-motor lawn law) calibrated to glide velocity 601 nm/s, stall 5–7 pN
  - Depletion/crowding lateral attraction ~0.11 pN that damps tip fluctuation (raises effective Lp) — soft potential, no off-rate
  - Cantilever mode-decomposition self-consistency gate (Eqs 1–3, q_n roots): recover input κ from tip fluctuation ∝ 1/Lp
  - Emergent principle: alignment DURABILITY (∝Lp), not PROBABILITY (θ_in-set, Lp-invariant), sets order — bending/thermostat is the lever
  - 2D self-propelled-filament regression scenario (S, Sk, N_B, η_d=16.3°) exercising bending+EV+motor+noise together
  · n_parameters=14
  · n_equations=6
  · n_validation_targets=11

## p37 [3] Taeyoon Kim · research — Molecular Biology of the Cell • 35:ar47, 1–12, April 1, 2024
- headline: Lumped Aplysia motor-clutch model — an acceptance-oracle for FF, not a runtime port: delivers single-molecule anchors (F_s=6, F_b=4 pN, K_clutch=0.3 pN/nm, v_u=100 nm/s), a mechanistic force-threshold adhesion-reinforcement lever, and a durotaxis validation scenario (biphasic latency 14->2->8 min, deformation 1.6->0.8 um, optimum 4 pN/nm).
  - Force-triggered adhesion reinforcement (F_t=10 pN) — add mechanistically via talin-unfold/integrin recruitment (5-10 pN), flips traction from oscillation to stable plateau
  - Bell slip-bond clutch off-rate F_b~4 pN anchor (apCAM regime; keep FF integrin catch-bond)
  - Single-NMII stall F_s=6 pN anchor (Lohner 2019) for Hill per-head force-velocity
  - Clutch spring stiffness K_clutch=0.3 pN/nm anchor
  - Substrate-stiffness durotaxis scenario + latency-time & substrate-deformation observables for ff_crawl_on_substrate
  - Retrograde-flow strong-coupling threshold: 20 nm/s = 80% drop from 100 nm/s unloaded (Aplysia)
  - Equal-force-sharing -> longer bond lifetime as cluster grows (already emergent in FF; sanity check)
  - Load-dependent myosin duty-ratio / Stam-Hocky adaptable NMIIA endorsed over lumped linear force-velocity
  · n_parameters=7
  · n_equations=5
  · n_validation_targets=11

## p25 [3] Taeyoon Kim · research — R E S E A R C H A R T I C L E
- headline: Kim-lab motility-assay model gives FF the Parallel Cluster Model motor-kinetics oracle plus an explicit-motors-only ring/band emergent-pattern benchmark, but its load-bearing constants are trapped in an unavailable SI.
  - Parallel Cluster Model (Erdmann–Schwarz) 8-head myosin kinetics as acceptance oracle for FF's Stam–Hocky/Hill minifilament
  - Explicit-motors-beat-implicit: rings only form with kinetically-resolved motors — external validation of FF's no-lumped-propulsion axiom
  - Collision-induced alignment via prompt LJ repulsion (κ_r) drives network→flock→band→ring; Δt small enough to engage repulsion at overlap onset
  - Geometric detachment: myosin head must release on reaching barbed end, not only force-dependently
  - l_p = κ_b·r0/k_BT (Eq2) unit-consistency gate + tangent-correlation ℓ_p estimator (Eq3) as FF buckling diagnostic
  - Buckling-induced F-actin fragmentation channel (Li 2017 pointer) for facilitated contraction / stress relaxation
  - FF gliding-assay benchmark scenario (10×10×0.1µm, fixed motors, sweep κ_r/κ_b/C/⟨L⟩) as mechanism-isolation regression
  · n_parameters=13
  · n_equations=4
  · n_validation_targets=12

## p18 [3] ? · ? — Durotactic Migration Driven by Anisotropic Matrix
- headline: Not a cytoskeleton paper — its FF value is a fully-specified tension/compression-asymmetric semiflexible collagen SUBSTRATE (buckling ÷10, bending, dilution p) plus anisotropic-stiffening / F_r/F_θ≈2 / durotaxis validation targets for FF-on-Kim-ECM scenarios.
  - ECM bond: tension-stiff/compression-soft ÷10 buckling asymmetry (Eq.13) — port to FF collagen kernel
  - ECM plasticity via random chain dilution p (0.3) — sub-isostatic Mikado connectivity
  - ECM semiflexible bending κ_b,M=0.001 nN·µm angular-harmonic on collagen triads
  - F_r/F_θ radial-tension decomposition (Eqs.19-20) as FF remodeling observable, target ≈2
  - E∥/E⊥ stiffness-probe protocol (F_ext=0.01 nN, Eq.17) as FF substrate characterization
  - Cell-matrix coupling as emergent viscous FA friction λ (Eq.9) — cross-check FF clutch durotactic sign
  · n_parameters=8
  · n_equations=8
  · n_validation_targets=11

## p16 [3] Taeyoon Kim · research — Computational Modeling of Motor-Driven
- headline: Kim-lineage MT-bundle BD paper: not actin, but a direct sibling of FF that hands over portable buckling↔crosslink-spacing physics, a biphasic cross-linker-density protocol, a cylindrical-segment drag law relevant to FF's grid-drag problem, and a curvature buckling observable — as cross-checks/validation, with no actin numbers transferable.
  - Buckling threshold tied to cross-link spacing: F_cr ∝ κ_b/L_uns² — emergent, magic-number-free; instrument in native cortex
  - Biphasic network response vs cross-linker density (buckling-resistance vs mobility-friction) → interior optimum; reusable contraction-sweep protocol/signature
  - Cylindrical-segment drag law ξ=3πμr_c(3+2r0/rc)/5 (Eq.3) — audit FF's per-segment γ, direct lever on the Σγ∝Nc crawl grid-drag problem
  - Per-node discrete curvature observable + 0.2 rad/µm buckling threshold — add to common/filament_math as a shared curved-fraction diagnostic
  - Multi-motor cooperativity required to exceed single-motor stall and buckle a cross-linked filament — load-sharing sanity check
  - Conditional: full dynein force-velocity + load-dependent detachment (Monzon 2018) recipe for activating the native --microtubules compartment (MT-MT sliding)
  · n_parameters=11
  · n_equations=5
  · n_validation_targets=8

## p76 [2.5] Makito Miyazaki · research — BIOPHYSICS AND COMPUTATIONAL BIOLOGY
- headline: FF lens — p76  (Sakamoto et al. 2022, "Geometric trade-off ... actomyosin droplet motility", PNAS 119:e2121147119)

## p73 [2.5] Makito Miyazaki · research — Optogenetic actin network assembly on lipid
- headline: FF lens — p73  (Yamamoto & Miyazaki 2025, OptoVCA branched-actin density gates ABP function, Nat Commun 16:7583)

## p55 [2.5] Makito Miyazaki · experimental — Directional Bleb Formation in Spherical Cells under Temperature Gradie
- headline: Experimental live-cell paper: gives FF the mechanistic bleb-nucleation sequence (cortex-membrane detachment + myosin-driven cortex rupture/fragmentation → directional bleb) and clean validation targets, but ZERO force laws or rate constants — trigger is a temperature gradient FF must NOT model, only its downstream cortical asymmetry.
  - cortex-membrane detachment bleb pathway
  - myosin-driven F-actin fragmentation/rupture
  - explicit ERM cortex-membrane linker Bell-Evans bond
  - duty-ratio/engaged-head local contractility knob
  - sharp-vs-uniform bleb bifurcation
  - local Laplace γ/R vs ΔP nucleation
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=9

## p51 [2.5] Taeyoon Kim · research — Mechanical Model for Durotactic Cell Migration
- headline: Cell side is lumped (reject); the real FF value is a viscoelastic ECM upgrade (per-bond bending + Kelvin–Voigt dashpot + sub-isostatic connectivity) plus a marquee no-decision durotaxis validation target for FF's clutch traction.
  - ECM bending springs per fiber junction (eq10, κ_b,M=1 nN·µm) → strain-stiffening — reuse Arp2/3 angle-harmonic kernel on collagen junctions
  - Per-bond Kelvin–Voigt dashpot on ECM fibers (eq11, F=-ξΔv, ξ 0.025–0.25 nN·hr/µm) → tunable network relaxation time τ=ξ/κ, the key durotaxis lever
  - Sub-isostatic ECM connectivity knob (random bond-removal p=0.4, Broedersz–MacKintosh) — coordination 0–6, currently missing in FF Mikado
  - Emergent durotaxis as a validation scenario for FF's existing motor-clutch traction — directed migration up a stiffness gradient with NO decision rule
  - Non-centripetal traction direction as an emergent cross-check (not imposed) against FF's explicit retrograde-flow field
  - Constant substrate strain energy vs stiffness (Gardel 2010) as an emergent-behavior falsification target for FF traction
  · n_parameters=7
  · n_equations=5
  · n_validation_targets=10

## p38 [2.5] Makito Miyazaki · research — Controlling Physical and Biochemical Parameters of Actin
- headline: Membrane-templated reconstitution gives FF a cell-matched nucleation-site density anchor (6,500-30,000/µm²) and two clean emergent-architecture validation scenarios (Arp2/3 perpendicular pillar vs formin+fascin horizontal vortex), but zero force laws or rate constants.
  - Nucleation-site areal density seeding (6,500-30,000 molecules/µm²)
  - Nucleator-type→architecture emergent switch (Arp2/3 pillar vs formin+fascin vortex)
  - Barbed-end-at-membrane protrusion polarity
  - Fascin parallel-bundling crosslinker (presence only)
  - Patch-size effect on growth + ≤2µm symmetry breaking
  - Area-not-shape architecture null test
  · n_parameters=10
  · n_equations=1
  · n_validation_targets=11

## p31 [2.5] Taeyoon Kim · research — PHYSICAL REVIEW RESEARCH 7, 013312 (2025)
- headline: Kim-lab explicit-motor BD assay: a motor-kinetics (PCM) and numerics ORACLE for FF, not a cortex mechanism — port the PCM rate cross-check, 7nm/40s^-1 binding, cylindrical-drag law, and collective-order metrics; the 2D clustering biology stays in-vitro-only.
  - myosin: adopt Parallel Cluster Model (Erdmann-Schwarz, 3-state/5-rate) as an ensemble unbinding+walking-rate ORACLE for the Stam-Hocky minifilament under load
  - myosin: quantize head binding to 7 nm F-actin lattice at k_+=40*N_H s^-1 (40/head)
  - filament drag: adopt closed-form cylindrical-segment drag zeta=3*pi*mu*r_c*(3+2*r0/r_c)/5 as a dimensional sanity check (swap mu to 65 Pa.s)
  - analysis: add nematic S, polar P, velocity-correlation, and grid M/M0 phase-separation metrics to FF's collective-behavior/cortex-flow toolkit
  - scenario: build a decoupled quasi-2D gliding-assay bench to unit-test FF explicit-motor+EV+bending against published polar/nematic ordering BEFORE trusting it in the 3D cortex
  - corroboration: explicit-motor > implicit-motor result directly endorses FF's mechanistic-over-lumped hard rule
  - length bookkeeping: adopt 2.6 nm/monomer helical-pitch conversion
  · n_parameters=12
  · n_equations=12
  · n_validation_targets=11

## p30 [2.5] Makito Miyazaki · research — PHYSICAL REVIEW RESEARCH 5, 013208 (2023)
- headline: Confined actomyosin flow→wave→rotation is a strong emergent-pattern BENCHMARK and two-knob (myosin σ_act × polymerization k_p) scenario for FF's fine-grained cortex — but the paper is lumped continuum, so it yields targets and a Péclet regime map, NOT micro-parameters.
  - Confined-cortex pattern-selection scenario (flow↔wave↔rotation) as a NEW FF emergent-observable test class
  - Two-knob mechanistic control: σ_act←myosin duty/density, k_p←Arp2/3+barbed-end kinetics — reproduce every arrow of the CalA/VVCA/CytoD/Y27632/CK666 perturbation map
  - Membrane-localized Arp2/3 nucleation: implement surface WASP/NPF enrichment (0.9R<r<R shell) mechanistically so the α_p<0.5 wave criterion EMERGES, not the lumped k_surf/k_bulk=10 ratio
  - Turnover sets the clock: explicit depolymerization/severing (k_d) makes wave period T≈1.5 min emerge, invariant to motor strength
  - Péclet regime map Pe=(ζμ)₀/(Dγ) — compute FF-derived Pe from micro-state to predict flow(≈38) vs wave(≈115) territory before native runs
  - Rotational symmetry breaking via system-wide staggered myosin onset (τd), stable only on the flow↔wave boundary — structural robustness test
  · n_parameters=9
  · n_equations=8
  · n_validation_targets=8

## p28 [2.5] Makito Miyazaki · review — Molecular Crystals and Liquid Crystals
- headline: Mostly muscle/spindle (DCM lane), but yields FF a confined actomyosin-ring self-organization scenario, a lattice-spacing-dependent myosin-binding extension, and an L_p=5-15 um bending anchor to reconcile with FF's 17 um.
  - mechanosensitive myosin binding gated by transverse head-to-actin spacing
  - confinement-driven emergent contractile-ring scenario (R < L_p)
  - rigid-rod bending-energy ring-vs-bulk criterion
  - L_p 5-15 um bending anchor vs FF 17 um reconciliation
  - depletion bundling interchangeable with alpha-actinin
  - constriction rate proportional to ring diameter
  · n_parameters=7
  · n_equations=3
  · n_validation_targets=7

## p52 [2] Makito Miyazaki · research — Microscopic Temperature Control Reveals Cooperative Regulation
- headline: In-vitro HMM motility assay (no sim) gives FF a myosin unloaded v_max sanity anchor (~15 µm/s at 37°C, skeletal not NMII) plus a mechanistic template for an allosteric actin-conformation ABP regulator with emergent Hill n≈4 cooperativity — no force law, moderate-low impact.
  - Myosin v_max sanity anchor: unloaded HMM sliding ~7-8 µm/s (28-30°C) → ~15 µm/s (37°C), ~2× Arrhenius scaling — cross-check Hill F=0 intercept (skeletal, NOT NMII)
  - Optional new module: allosteric actin-conformation ABP regulator (drebrin/tropomyosin-class) that gates myosin duty by changing actin pitch, NOT by steric competition — regulator stays bound (Qdot p=0.344)
  - Cooperativity via nearest-neighbor pitch coupling: actin half-pitch 36→40 nm propagating ≤2 helical repeats → 1-D Ising-like backbone state; emergent Hill n must go 1→4, coupling-set NOT tuned
  - Regulator binding kinetics: langmuir occupancy p=[D]/(K_d+[D]), K_d 6nM(25°C)→32nM(37°C), 1:5 regulator:actin stoichiometry, per-site stochastic on/off
  - Hill dose-response v=1−C^n(1−b)/(k^n+C^n) as validation ORACLE only (never runtime): n≈4, IC50 ~3→>30 nM at 37°C
  - FF-motility-assay scenario: 2-D sparse gliding assay of 1-5µm filaments over anchored heads, displacement/0.5s window — dev-scale motor calibration bench, non-authoritative
  · n_parameters=16
  · n_equations=5
  · n_validation_targets=12

## p49 [2] Makito Miyazaki · research — Kinetic Scheme of Myosin Phosphorylation by ZIP Kinase
- headline: Supplies myosin's upstream activation-switch kinetics (sequential Ser19→Thr18 MRLC phosphorylation, kcat 0.4–0.95 s⁻¹, Km 2.5–4.2 µM) — a candidate activation-fraction gate for FF's Hill minifilaments and the cortical-γ open problem, but no force law, off-rate, or mechanical datum.
  - Myosin: add explicit MRLC phospho activation state (0P/1P/2P) gating Hill-head engagement — activation fraction becomes derived, not implicitly 1.0
  - Sequential Michaelis-Menten activation kernel (Ser19 then Thr18) as O(N_minifil) stochastic per-minifilament updater alongside binding kinetics
  - Activation timescales for cortical-tension γ modulation: kcat 0.60/0.40 s⁻¹ (τ_act≈1.7/2.5 s), Km 4.2 µM — sets seconds-scale contractile ramp
  - Sequential ordering constraint: monophospho ≈35:1 Ser19:Thr18 (1P is Ser19, not random 50/50)
  - 2P enhances filament assembly but NOT unloaded gliding velocity — assembly-vs-force decoupling in activation logic
  - Pair on-limb phosphorylation with a PP1/MYPT1 phosphatase off-rate (NOT in paper) for a genuine steady-state active fraction
  - EXCLUDE Model-2 ZIPK-MHC competition (KI=0.30 µM) as in-vitro artifact; Arrhenius-correct 25°C rates to 37°C
  · n_parameters=9
  · n_equations=5
  · n_validation_targets=7

## p47 [2] Taeyoon Kim · research — Role of actin filaments and cis binding in
- headline: Rigid-actin cadherin-clustering paper: no force law for FF, but supplies cortex actin length/concentration/turnover bands (tip 20 s + stalk 2 min), a 3.5 nm F-actin excluded-volume radius, and a corralling->sub-diffusion target FF should reproduce emergently.
  - Bimodal cortex F-actin turnover: fast bundle-tip pool ~20 s + slow stalk half-life ~2 min (vs single-rate FF turnover)
  - F-actin excluded-volume radius 3.5 nm (7 nm diameter) as LJ sigma / contact-cutoff anchor
  - Corralling -> sub-diffusion of membrane tracers: reproduce EMERGENTLY on native cortex (MSD slope <1), never hard-code rejection
  - Cortex actin-concentration band (non-muscle 46-70 uM) as density cross-check vs ratified KB-3.18
  - Filament tilt (phi) controls effective mesh: near-normal->fine/punctate, near-tangential->open/linear
  - REJECT: constant force-independent k_dis (keep FF Bell-Evans k_off=k_off0 exp(F/F_b))
  · n_parameters=11
  · n_equations=3
  · n_validation_targets=4

## p46 [2] Taeyoon Kim · research — The nature of cell division forces in epithelial
- headline: DCM paper, low FF value: no filament-scale force laws to port, but it specifies a new native-scale emergent-cytokinetic-ring FF scenario and hard oracle numbers (300-500 nN ring force, 5%-vol->60%-elongation, myosin -60%, r=0.6 anillin) plus a drag~area idea for the grid-drag problem.
  - Emergent cytokinetic actomyosin ring as a NEW native-scale FF scenario
  - Anillin as a ring-specific crosslinker channel (graded contraction, r=0.6)
  - Ring force 300-500 nN as whole-structure contraction magnitude test
  - Near-volume-conservation coupling (5% vol loss -> 60% elongation) validates turgor ON
  - Drag ~ local area (zeta=lambda*A) template for the FF grid-drag problem
  - Myosin-necessity test: myosin OFF abolishes ~60% elongation
  · n_parameters=7
  · n_equations=1
  · n_validation_targets=10

## p09 [2] Taeyoon Kim · research — Reciprocal folding dynamics in cellular networks at the stroma-basemen
- headline: Tissue-scale FEM/in-vitro folding paper — no molecular force laws or rate constants, but a clean single-cell traction phenotype (polarized F-actin-tip/myosin-body/paxillin-FA zipper, blebbistatin-abolished, >90% fiber-protrusion-fold alignment) that FF's existing crawl/clutch/myosin chain should emergently reproduce as validation, not a parameter source.
  - zipper traction layout emergent target
  - blebbistatin myosin-off null gate
  - FA-polarity analysis metric (Eq.2)
  - alignment-% metric (Eq.1)
  - durotaxis fiber-align-before-protrusion
  - plastic proteolytic ECM future channel
  - FF-on-stiff-substrate scenario
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=7

## p66 [1.5] Taeyoon Kim · research — Bulletin of Mathematical Biology (2019) 81:3301–3321
- headline: Lumped 2-point over-damped agent CIL/PRW model — no fine-grained physics to port; useful to FF only as behavioral crawl oracles (MSD ballistic→diffusive crossover, ~10-17 nm/s speed) and analysis conventions, zero parameters registrable.
  - MSD log-slope crawl-persistence metric (s=d ln MSD/d ln τ): ballistic→diffusive analysis harness for FF's emergent crawl
  - Behavioral speed oracle: single-cell crawl ~10-17 nm/s (50-60 µm/hr) cross-checks FF's 28-45 nm/s
  - Finite-lifetime lamellipodium turnover (T_F~6 min) confirms FF clutch-turnover crawl-emergence story — as derived observable, not imposed timer
  - Redirection-angle distribution peaked at 0° + biphasic persistence vs protrusion lifetime as FF crawl targets
  - CIL as substrate-anchor partitioning — parked idea for future multi-cell FF integrin-occupancy competition (not single-cell)
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=6

## p34 [1.5] ? · ? — Integr. Biol., 2015, 7, 1093--1108 | 1093
- headline: Narrative review, no original physics — its FF value is a roadmap to co-author Kim's primary BD-actin oracles plus two clean native-scale FF scenarios (subnetwork-ablation, motor:ACP percolation); no parameter or force law is FF-usable or supersedes any current value.
  - Transient-ACP stress-relaxation + emergent load-bearing subnetwork (Kim BD ref 4/47) — validates FF crosslinker turnover
  - Bending→stretching strain-stiffening transition as a native FF emergent target (Fig. 3 three-regime K(σ))
  - Motor:crosslinker percolation phase transition (connected active-stiffening ↔ disconnected) as an FF density-sweep scenario
  - Fig. 1b subnetwork-ablation protocol: 25-75% ACP removal, most stress retained — FF load-path robustness scenario
  - ~10 s non-thermal motor-driven fluctuation correlation time as an FF active-noise validation target
  - Cation/nucleotide-tuned single-filament κ → network stiffening (Bidone ref 30) — optional state-dependent-κ refinement
  · n_parameters=7
  · n_equations=2
  · n_validation_targets=5

## p33 [1.5] Makito Miyazaki · experimental — Accurate polarity control and parallel alignment of actin
- headline: In-vitro actin nanofab/motility method paper — no dynamics model; useful to FF only as ℓ_p/diameter cross-checks (actin 15–20 µm, MT 5000–6000 µm) plus two soft mechanism pointers (gelsolin capping, α-actinin flow-bundling); zero force laws, rate schemes, or cell-mechanics validation targets.
  - Actin ℓ_p = 15–20 µm + d = 6 nm — cross-check for FF filament bending κ = k_BT·ℓ_p (confirms FF ℓ_p≈17µm, no change)
  - MT ℓ_p = 5000–6000 µm / d = 25 nm — bending-stiffness cross-check for the native microtubule compartment
  - Gelsolin plus-end capping + two-monomer nucleation — optional future turnover hook to make filament length distribution emergent (no rate constants given)
  - α-actinin bundling under an aligning field → parallel isopolar bundle — corroborates FF filopodium/stress-fiber crosslinker motif
  - Myosin-V unloaded speed 0.62 µm/s — out-of-model motor-speed reference ONLY (wrong isoform; not a myosin-II Hill/stall input)
  · n_parameters=5
  · n_equations=0
  · n_validation_targets=2

## p23 [1.5] Makito Miyazaki · method — Bio-protocol 16(8): e5656. DOI: 10.21769/BioProtoc.5656
- headline: Wet-lab optogenetic Arp2/3-network protocol: no importable physics, but yields a membrane-NPF-density nucleation-gating idea, a clean in-vitro branched-network unit-test scenario, and depletion/saturation validation shapes for FF's Arp2/3 + monomer-reservoir modules.
  - Arp2/3 nucleation rate gated by membrane NPF surface-density field σ_NPF (spatially patterned, density-graded cortex/lamellipodium seeding)
  - Monomer-depletion-limited transient overshoot: finite G-actin reservoir must reproduce rise-then-decline + apical-high z-gradient under strong sustained nucleation
  - Capping-protein rate (not length caps) as the mesh-density/branch-density regulator; profilin biases to Arp2/3-templated barbed-end growth
  - In-vitro branched-network unit-test scenario: flat membrane patch, patterned NPF, actin:Arp2/3:profilin:CP ≈ 200:4:600:1
  - Validation shapes: density saturates with NPF; thickness grows-then-plateaus (~30-40 µm in-vitro); polymerization reversible on trigger removal
  · n_parameters=6
  · n_equations=0
  · n_validation_targets=5

## p72 [1] Makito Miyazaki · research — Protein design of two-component tubular
- headline: Protein-design paper with only one FF-transferable item: an independent actin persistence-length datum (12.5 µm; 9–20 µm lit) plus the 2D WLC relation — a sanity cross-check for FF's bending stiffness, not a new mechanism.
  - actin persistence-length cross-check (L_p 12.5 µm / 9–20 µm)
  - 2D worm-like-chain relation for L_p fitting
  - κ=L_p·k_BT bending-stiffness anchor
  - single-filament thermal-bending unit test
  - 2D-vs-3D WLC dimensionality caveat
  · n_parameters=2
  · n_equations=1
  · n_validation_targets=2

## p56 [1] Makito Miyazaki · method — Quantitative Analysis of the Lamellarity of Giant Liposomes Prepared b
- headline: Membrane-fabrication/imaging methods paper: licenses FF's single-bilayer containment and supplies a turgor-free encapsulated-actin confinement scenario, but contributes zero filament/motor/crosslinker mechanism, force law, or rate constant.
  - Encapsulated-actin-in-GUV confinement scenario (50 µM F-actin, single bilayer, iso-osmotic turgor-free control)
  - Single-bilayer membrane assumption empirically licensed (97% unilamellar inverted-emulsion GUVs)
  - PEG-depletion actin bundling as in-vitro comparator only — keep FF bundling explicit-crosslinker (do NOT adopt AO proxy)
  - α-hemolysin membrane ruler (bilayer ~5.2 nm) for membrane-geometry sanity
  · n_parameters=2
  · n_equations=1
  · n_validation_targets=0

## p29 [1] Taeyoon Kim · research — Mechanics and Morphology of Proliferating Cell Collectives with Self-I
- headline: Proliferating-rigid-rod colony paper with zero cytoskeletal content; FF gains only two logged cross-engine reminders (complementarity contact, length-resolved rod drag) and no mechanism, force law, parameter, scenario, or validation target — the real payoff is DCM's.
  - Constraint-based no-overlap contact vs FF LJ excluded-volume stiffness (cross-engine note)
  - Load-inhibited elongation exp(-lambda*sigma) — FF's ratchet already finer-grained
  - Length-resolved rod drag 1/(xi*l), 12/(xi*l^3) for FF grid-drag audit
  - No portable cytoskeletal mechanism/force law/parameter present
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p11 [1] Makito Miyazaki · research — RESEARCH ARTICLE |  APRIL 07 2011
- headline: Single-molecule hidden-trajectory INFERENCE method (Go-and-Back / Onsager–Machlup), not a filament/motor force law — near-zero FF runtime residue; only a thermostat unit test and a shelved offline estimator for future bead-assay validation.
  - Overdamped Langevin with FDT noise Var=2γk_BT — restates FF's existing BAOAB substrate
  - Free-bead thermostat sanity check SD=√(2k_BTΔt/γ)
  - Go-and-Back offline hidden-trajectory estimator for future bead-assay validation (downstream only)
  - No forward mechanism/force-law/rate-constant to adopt
  · n_parameters=0
  · n_equations=1
  · n_validation_targets=1

## p08 [1] Makito Miyazaki · research — Chimeras of kinesin-6 and kinesin-14 reveal head-
- headline: In-vivo yeast mitotic-kinesin chimera paper with NO force law, bond kinetics, or rheology; only conceptual motor head/neck/tail factorization and coarse spindle-sliding velocities (0.35-2.39 um/min) that don't touch FF's actin-cortex modules or gates.
  - Factorize motor state into orthogonal knobs: directionality/localization vs speed vs force-generating stepping
  - Bipolar/tetrameric plus-end motor sliding antiparallel filaments at an overlap zone (future MT reuse of Stam-Hocky path)
  - Emergent excess-motor-drive -> buckling -> breakage instability as qualitative corroboration for FF buckling agenda
  - Resisting-boundary tensile force opposing motor-driven sliding as a future boundary motif
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p01 [1] ? · ? — RESEARCH ARTICLE |  FEBRUARY 28 2011
- headline: Bayesian single-molecule inference method (2-bead F1-ATPase toy) — no cytoskeletal mechanism, force law, or validation number for FF; value is only an optional Langevin-parameter calibration/inference tool and an integrator cross-check.
  - Overdamped Langevin bead-spring + FDT noise (2γk_BT) as integrator correctness cross-check — FF already exceeds via BAOAB
  - Onsager–Machlup midpoint/Stratonovich action Eq.(12) as a unit-test oracle for the discretization
  - Optional out-of-band Bayesian hidden-parameter calibration (integrate out hidden DOF, minimize Hτ=−τ⁻¹lnP([y]|Π)) to infer effective coarse-grained bead stiffness/drag
  - Loss-of-precision guardrail: only trust inferred params when k*≪h* and Δt<τr≈Γ*/h*, report sample-averaged variance
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p54 [0.5] Makito Miyazaki · experimental — on the number and size of liposomes formed
- headline: Inverted-emulsion liposome-formation methods paper — zero cytoskeletal/force-law content; nothing to port, parity-check, or validate in the Cytosim-physics FF engine (null result recorded).
  - NONE — no cytoskeletal filament/motor/crosslinker mechanism in paper
  - Liposome-yield levers (charge repulsion, headgroup steric neck-narrowing) are membrane-manufacture only, not FF force laws
  - Marginal out-of-scope tie only: cell-sized 5–40 µm vesicle container geometry for a hypothetical future encapsulated-cortex experiment
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p06 [0.5] Taeyoon Kim · review — Contents lists available at ScienceDirect
- headline: Plant cell-branching review (turgor cell-wall FE + microtubule cortical-array self-organization) with no transcribed equations and no actomyosin content — essentially not applicable to the FF actin/myosin/crosslinker engine.
  - MT dynamic instability 3-state channel
  - angle-dependent MT-MT zipper/catastrophe rule
  - katanin crossover-selective severing
  - turgor-shell load-relax-remesh growth loop (DCM analogy)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0

## p20 [0] ? · ? — PAPERS |  SEPTEMBER 01 2022
- headline: Elementary-school Brownian-motion outreach paper — zero FF-adoptable mechanism, parameter, or validation target; only the canonical MSD=4Dt diffusion law appears.
  - none — K-12 Brownian-motion outreach paper, no cytoskeletal content
  - only equation is textbook Einstein 2D diffusion MSD=4Dt (Eq.1)
  - BAOAB thermostat already satisfies MSD=2nDt — invariant, not a new mechanism
  - molecular motors mentioned qualitatively only (no duty ratio / stall / F-v)
  · n_parameters=0
  · n_equations=2
  · n_validation_targets=0

## p12 [0] Makito Miyazaki · review — Honmachi, Sakyo-ku, Kyoto 606-8501, Japan. ORCID iD: https://orcid.org
- headline: p12 is a 3-page BSJ symposium editorial (Miyazaki & Kosugi 2022) with zero equations, parameters, or cytoskeletal mechanisms — nothing applicable to the FF engine.
  - NONE — symposium editorial, no cytoskeletal mechanism
  - No filament/crosslinker/myosin/formin/Arp2/3/clutch physics present
  - Only lead worth following is a citation pointer (ActuAtor ref [7], primary paper not this one)
  · n_parameters=0
  · n_equations=0
  · n_validation_targets=0