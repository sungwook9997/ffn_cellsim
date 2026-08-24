# Literature Dossier — Intracellular Transport, Stress-Fiber Assembly & Field Actuation of the Cytoskeleton

> Compiled 2026-07-16 for the ffn_cellsim CFD-based engine research program. **All DOIs machine-verified** (crossref/PubMed) — citation-integrity HARD rule. 66 papers across 15 topics; 52 newly registered to the Notion KB (SourceEvidence), the rest already present.

Legend: 🆕 = newly registered SourceEvidence this session · 📚 = already in KB.

## Actin retrograde flow + molecular clutch

- 🆕 **Case & Waterman 2015** — *Integration of actin dynamics and cell adhesion by a three-dimensional, mechanosensitive molecular clutch*. Nat Cell Biol 17(8):955-963 (Review). [10.1038/ncb3191](https://doi.org/10.1038/ncb3191)
  - **Finding:** Review: the clutch is a conserved 3D nano-architecture of mechanosensitive protein-protein interactions transmitting cytoskeletal force to ECM
  - **Relevance:** Authoritative clutch-architecture review; frames the multi-protein force-transmission stack the FF adhesion module abstracts
- 📚 **Chan & Odde 2008** — *Traction dynamics of filopodia on compliant substrates*. Science 322(5908):1687-1691. [10.1126/science.1163595](https://doi.org/10.1126/science.1163595)
  - **Finding:** Stochastic motor-clutch model predicts frictional-slippage (stiff, fast flow, low force) vs load-and-fail oscillation (soft, slow flow, high force); stiffness switch near 1 kPa
  - **Relevance:** Canonical motor-clutch mechanism the project cites as the mechanistic (not lumped) clutch model; KB oracle for stiffness sensing
- 🆕 **Gardel, Sabass, Ji, Danuser, Schwarz & Waterman 2008** — *Traction stress in focal adhesions correlates biphasically with actin retrograde flow speed*. J Cell Biol 183(6):999-1005. [10.1083/jcb.200810060](https://doi.org/10.1083/jcb.200810060)
  - **Finding:** Traction vs F-actin flow is biphasic with a threshold flow speed of 8-10 nm/s independent of FA density/age; inverse near edge, direct in mature FAs
  - **Relevance:** Primary quantitative acceptance target (biphasic traction-flow law) for the FF molecular-clutch/traction actuation module
- 🆕 **Mitchison & Kirschner 1988** — *Cytoskeletal dynamics and nerve growth*. Neuron 1(9):761-772. [10.1016/0896-6273(88)90124-9](https://doi.org/10.1016/0896-6273(88)90124-9)
  - **Finding:** Founding conceptual framework: actin assembles at the leading membrane edge and undergoes rearward (retrograde) flow / treadmilling to drive growth-cone advance
  - **Relevance:** Origin of the retrograde-flow + clutch paradigm that the FF molecular-clutch / traction module implements; conceptual root of the whole topic
- 🆕 **Ponti, Machacek, Gupton, Waterman-Storer & Danuser 2004** — *Two distinct actin networks drive the protrusion of migrating cells*. Science 305(5691):1782-1786. [10.1126/science.1100533](https://doi.org/10.1126/science.1100533)
  - **Finding:** Fluorescent speckle microscopy resolves two colocalized networks: a lamellipodium that disassembles within 1-3 um and a lamella where actomyosin contraction couples to adhesion
  - **Relevance:** Defines the two-zone kinematics (lamellipodium vs lamella retrograde flow) the FF cortex/protrusion transport model must reproduce
- 🆕 **Thievessen, Thompson, Berlemont, ... Campbell & Waterman 2013** — *Vinculin-actin interaction couples actin retrograde flow to focal adhesions, but is dispensable for focal adhesion growth*. J Cell Biol 202(1):163-177. [10.1083/jcb.201303129](https://doi.org/10.1083/jcb.201303129)
  - **Finding:** Vinculin acts as the molecular clutch: it slows F-actin flow in maturing FAs, generates high ECM traction, but FA growth rate tracks flow speed independently of vinculin
  - **Relevance:** Identifies the specific clutch protein (vinculin) coupling flow->traction; informs which linker the FF clutch bond should represent
- 🆕 **Vallotton, Gupton, Waterman-Storer & Danuser 2004** — *Simultaneous mapping of filamentous actin flow and turnover in migrating cells by quantitative fluorescent speckle microscopy*. Proc Natl Acad Sci USA 101(26):9660-9665. [10.1073/pnas.0300552101](https://doi.org/10.1073/pnas.0300552101)
  - **Finding:** qFSM maps coupled F-actin flow + net assembly/disassembly; ~90% of polymer disassembles at the lamellipodium-lamellum junction, contraction-coupled depolymerization in convergence zone
  - **Relevance:** Quantitative flow+turnover field the transport-limited actin and retrograde-flow modules are validated against

## Stress-fiber assembly

- 🆕 **Burnette, Manley, Sengupta, ... Kachar & Lippincott-Schwartz 2011** — *A role for actin arcs in the leading-edge advance of migrating cells*. Nat Cell Biol 13(4):371-381. [10.1038/ncb2205](https://doi.org/10.1038/ncb2205)
  - **Finding:** Lamellipodial actin is condensed by myosin II into an actin arc that moves rearward and slows at FAs; arcs mechanically bridge lamellipodium and lamella to advance the edge
  - **Relevance:** Links protrusion actin -> arc -> adhesion; couples the retrograde-flow and SF modules structurally
- 🆕 **Colombelli, Besser, Kress, ... Schwarz & Stelzer 2009** — *Mechanosensing in actin stress fibers revealed by a close correlation between force and protein localization*. J Cell Sci 122(10):1665-1679. [10.1242/jcs.042986](https://doi.org/10.1242/jcs.042986)
  - **Finding:** SF laser nanosurgery + computational SF model: zyxin relocalizes to high-tension intermediate anchor points; measured localization matches computed force distribution -> direct force sensing along SFs
  - **Relevance:** Force-along-SF computational model + mechanosensing datum; template for FF SF force-distribution readout
- 📚 **Hotulainen & Lappalainen 2006** — *Stress fibers are generated by two distinct actin assembly mechanisms in motile cells*. J Cell Biol 173(3):383-394. [10.1083/jcb.200511093](https://doi.org/10.1083/jcb.200511093)
  - **Finding:** Canonical taxonomy: dorsal SFs assembled by formin (mDia1) at FAs; transverse arcs by annealing of myosin bundles + Arp2/3 bundles; both convert to ventral SFs
  - **Relevance:** Defines the dorsal/arc/ventral SF architecture the stress-fiber target build must instantiate
- 🆕 **Kumar, Maxwell, Heisterkamp, ... Mazur & Ingber 2006** — *Viscoelastic retraction of single living stress fibers and its impact on cell shape, cytoskeletal organization, and extracellular matrix mechanics*. Biophys J 90(10):3762-3773. [10.1529/biophysj.105.071506](https://doi.org/10.1529/biophysj.105.071506)
  - **Finding:** Laser-severed single SFs recoil as actomyosin-tensed viscoelastic cables (prestress); recoil abolished by MLCK inhibition; ECM behaves as a physical extension of the cytoskeleton
  - **Relevance:** Quantitative viscoelastic-prestress SF datum; acceptance oracle for single-SF mechanics in the SF target
- 🆕 **Oakes, Beckham, Stricker & Gardel 2012** — *Tension is required but not sufficient for focal adhesion maturation without a stress fiber template*. J Cell Biol 196(3):363-374. [10.1083/jcb.201107042](https://doi.org/10.1083/jcb.201107042)
  - **Finding:** SF template (formin/alpha-actinin-1 dependent) is needed for FA compositional maturation; maturation still occurs with cellular tension reduced ~80% -> tension necessary but not sufficient
  - **Relevance:** Constrains FA maturation logic (structural SF template vs pure tension) for the adhesion/SF coupling
- 🆕 **Tanner, Boudreau, Bissell & Kumar 2010** — *Dissecting regional variations in stress fiber mechanics in living cells with laser nanosurgery*. Biophys J 99(9):2775-2783. [10.1016/j.bpj.2010.08.071](https://doi.org/10.1016/j.bpj.2010.08.071)
  - **Finding:** Peripheral MLCK-dependent SFs show lower effective elasticity and higher plateau retraction than central ROCK-dependent SFs; SF populations share load and can absorb one another's roles
  - **Relevance:** Regional SF mechanical heterogeneity + load-sharing for spatially resolved SF modeling
- 🆕 **Tojkander, Gateva & Lappalainen 2012** — *Actin stress fibers - assembly, dynamics and biological roles*. J Cell Sci 125(8):1855-1864 (Review). [10.1242/jcs.098087](https://doi.org/10.1242/jcs.098087)
  - **Finding:** Review consolidating SF subtypes (dorsal SF, transverse arc, ventral SF), their assembly, turnover, and contractile/mechanosensing roles
  - **Relevance:** The J Cell Sci 2012 review named in the seed; reference map for SF-subtype definitions and dynamics
- 🆕 **Tojkander, Gateva, Schevzov, ... Gunning & Lappalainen 2011** — *A molecular pathway for myosin II recruitment to stress fibers*. Curr Biol 21(7):539-550. [10.1016/j.cub.2011.03.007](https://doi.org/10.1016/j.cub.2011.03.007)
  - **Finding:** Four tropomyosins + Dia2 formin nucleate Tm4-decorated cortical filaments; myosin II is recruited Tm4-dependently, then anneals with Arp2/3 filaments to form contractile transverse arcs
  - **Relevance:** Mechanistic transverse-arc assembly pathway (myosin recruitment) for the SF module; the actual mechanism paper behind the '2012 arc' seed
- 🆕 **Vallenius 2013** — *Actin stress fibre subtypes in mesenchymal-migrating cells*. Open Biol 3(6):130001 (Review). [10.1098/rsob.130001](https://doi.org/10.1098/rsob.130001)
  - **Finding:** Review: SF subtypes (dorsal, transverse arc, ventral) differ in molecular composition and spatially/temporally partition myosin-II contractility during mesenchymal migration
  - **Relevance:** Subtype composition reference for parameterizing distinct SF populations in the model

## Formin + polymerization kinetics

- 🆕 **Cao, Kerleau, Suzuki, ... Romet-Lemonne & Jegou 2018** — *Modulation of formin processivity by profilin and mechanical tension*. eLife 7:e34176. [10.7554/eLife.34176](https://doi.org/10.7554/eLife.34176)
  - **Finding:** Piconewton tensile force on the filament accelerates formin (mDia1/mDia2) dissociation by orders of magnitude, overriding profilin stabilization; two distinct force-selected dissociation pathways
  - **Relevance:** Force-dependent formin detachment - the mechanistic, force-coupled formin kinetics for the field-actuation program
- 🆕 **Drenckhahn & Pollard 1986** — *Elongation of actin filaments is a diffusion-limited reaction at the barbed end and is accelerated by inert macromolecules*. J Biol Chem 261(27):12754-12758. [10.1016/S0021-9258(18)67157-1](https://doi.org/10.1016/S0021-9258(18)67157-1)
  - **Finding:** Barbed-end elongation is diffusion-limited (k+ ~ 1e7 M^-1 s^-1 per centipoise, extrapolating through origin vs viscosity); inert crowders raise k+ via excluded-volume activity increase
  - **Relevance:** Justifies diffusion/viscosity dependence of elongation and macromolecular-crowding acceleration - directly ties polymerization to the cytoplasm-rheology/transport program
- 🆕 **Footer, Kerssemakers, Theriot & Dogterom 2007** — *Direct measurement of force generation by actin filament polymerization using an optical trap*. Proc Natl Acad Sci USA 104(7):2181-2186. [10.1073/pnas.0607052104](https://doi.org/10.1073/pnas.0607052104)
  - **Finding:** Optical trap: growth of ~8 parallel filaments stalls at ~1 pN (single-filament stall load), with large length fluctuations indicating only the longest filament touches the barrier -> dynamic instability limits bundle force
  - **Relevance:** Direct experimental stall-force calibration for polymerization actuation in bundle geometry
- 🆕 **Kovar & Pollard 2004** — *Insertional assembly of actin filament barbed ends in association with formins produces piconewton forces*. Proc Natl Acad Sci USA 101(41):14725-14730. [10.1073/pnas.0405902101](https://doi.org/10.1073/pnas.0405902101)
  - **Finding:** Filament growth between an immobilized formin and an anchor buckles segments as short as 0.7 um, demonstrating single-filament polymerization forces >1 pN (near theoretical max) via insertional assembly
  - **Relevance:** Direct polymerization-force datum for the actin-polymerization actuation (protrusion) force scale
- 📚 **Mogilner & Oster 1996** — *Cell motility driven by actin polymerization*. Biophys J 71(6):3030-3045. [10.1016/S0006-3495(96)79496-1](https://doi.org/10.1016/S0006-3495(96)79496-1)
  - **Finding:** Elastic Brownian ratchet: thermal bending fluctuations of polymerizing filaments (not just the load) rectify to produce directed protrusive force; quantitatively fits Listeria and lamellipodia
  - **Relevance:** Foundational polymerization force-generation theory; oracle for protrusion force-velocity
- 📚 **Mogilner & Oster 2003** — *Force generation by actin polymerization II: the elastic ratchet and tethered filaments*. Biophys J 84(3):1591-1605. [10.1016/S0006-3495(03)74969-8](https://doi.org/10.1016/S0006-3495(03)74969-8)
  - **Finding:** Tethered-ratchet extension: attached (tethered) plus working filaments yield the force-velocity relation for Listeria/lamellipodial protrusion, resolving symmetry breaking and low-load motion
  - **Relevance:** Refined force-velocity law for polymerization-driven protrusion actuation
- 🆕 **Paul & Pollard 2008** — *The role of the FH1 domain and profilin in formin-mediated actin-filament elongation and nucleation*. Curr Biol 18(1):9-19. [10.1016/j.cub.2007.11.062](https://doi.org/10.1016/j.cub.2007.11.062)
  - **Finding:** Elongation rate rises with the number of FH1 polyproline/profilin tracks; FH1 binding of profilin-actin is rate-limiting for subunit delivery up to rates >=88 s^-1; profilin inhibits FH2 nucleation
  - **Relevance:** Quantitative FH1-profilin delivery model for formin processive elongation parameters
- 📚 **Pollard & Borisy 2003** — *Cellular motility driven by assembly and disassembly of actin filaments*. Cell 112(4):453-465 (Review). [10.1016/s0092-8674(03)00120-x](https://doi.org/10.1016/s0092-8674(03)00120-x)
  - **Finding:** Dendritic-nucleation / treadmilling model: Arp2/3 branching, capping termination, ATP hydrolysis + ADF/cofilin debranching, profilin nucleotide exchange reconstitute motility
  - **Relevance:** Master reference for the actin turnover cycle the FF fine-grained filament dynamics reproduces
- 📚 **Pollard 1986** — *Rate constants for the reactions of ATP- and ADP-actin with the ends of actin filaments*. J Cell Biol 103(6 Pt 2):2747-2754. [10.1083/jcb.103.6.2747](https://doi.org/10.1083/jcb.103.6.2747)
  - **Finding:** Measures association/dissociation rate constants for ATP- and ADP-actin at barbed and pointed ends (linear elongation up to 20 uM); the canonical actin polymerization kinetics
  - **Relevance:** Ground-truth on/off rate constants for the FF actin polymerization/turnover kinetics
- 🆕 **Romero, Le Clainche, Didry, ... Pantaloni & Carlier 2004** — *Formin is a processive motor that requires profilin to accelerate actin assembly and associated ATP hydrolysis*. Cell 119(3):419-429. [10.1016/j.cell.2004.09.039](https://doi.org/10.1016/j.cell.2004.09.039)
  - **Finding:** FH1-FH2 formin is a processive barbed-end motor; profilin required; increases profilin-actin association rate constant ~15-fold and couples ATP hydrolysis to processive polymerization
  - **Relevance:** Mechanistic basis for formin-accelerated filament elongation in the FF formin/protrusion kinetics
- 🆕 **Vavylonis, Kovar, O'Shaughnessy & Pollard 2006** — *Model of formin-associated actin filament elongation*. Mol Cell 21(4):455-466. [10.1016/j.molcel.2006.01.016](https://doi.org/10.1016/j.molcel.2006.01.016)
  - **Finding:** Kinetic model with FH2 gating (open/closed) + FH1 multi-site profilin-actin transfer; predicts an optimal profilin concentration maximizing elongation, with high profilin suppressing growth
  - **Relevance:** Closed-form formin-elongation model usable as an acceptance oracle for the FF formin kinetics

## Myosin-II force

- 🆕 **Finer, Simmons & Spudich 1994** — *Single myosin molecule mechanics: piconewton forces and nanometre steps*. Nature 368(6467):113-119. [10.1038/368113a0](https://doi.org/10.1038/368113a0)
  - **Finding:** Feedback optical-trap measurement of a single myosin-actin interaction: ~11 nm steps at low load and 3-4 pN single-molecule isometric force transients
  - **Relevance:** Single-motor step size and force calibration for the mechanistic myosin head model
- 📚 **Hill 1938** — *The heat of shortening and the dynamic constants of muscle*. Proc R Soc Lond B 126(843):136-195. [10.1098/rspb.1938.0050](https://doi.org/10.1098/rspb.1938.0050)
  - **Finding:** Hyperbolic force-velocity relation of muscle (the Hill equation) and heat-of-shortening constants; foundational contractile force-velocity law
  - **Relevance:** Project's sanctioned mechanistic myosin force-velocity oracle (Hill 1938) per the architectural-principle table
- 📚 **Kovacs, Wang, Hu, Zhang & Sellers 2003** — *Functional divergence of human cytoplasmic myosin II: kinetic characterization of the non-muscle IIA isoform*. J Biol Chem 278(40):38132-38140. [10.1074/jbc.M305453200](https://doi.org/10.1074/jbc.M305453200)
  - **Finding:** Non-muscle myosin IIA transient kinetics: high actomyosin-ADP affinity, ADP release >> steady-state ATPase, so NM2A spends only a small fraction of its cycle strongly actin-bound (low duty ratio)
  - **Relevance:** Isoform-specific NM2A duty-ratio/kinetics for the mechanistic myosin-II minifilament parameterization
- 📚 **Stam, Alberts, Gardel & Munro 2015** — *Isoforms Confer Characteristic Force Generation and Mechanosensation by Myosin II Filaments*. Biophys J 108(8):1997-2006. [10.1016/j.bpj.2015.03.030](https://doi.org/10.1016/j.bpj.2015.03.030)
  - **Finding:** Cross-bridge minifilament model: motor force output on an elastic load is set by two timescales - F-actin attachment (varying sharply with ensemble size, duty ratio, load) and force build-up - conferring isoform-specific mechanosensation
  - **Relevance:** The 'Stam-Hocky' multi-head bipolar minifilament model the project prefers over AFINES's 2-head spring for mechanistic myosin II

## Cytoplasm rheology + transport-limited actin

- 🆕 **Appalabhotla, Butler, Bear & Haugh 2023** — *G-actin diffusion is insufficient to achieve F-actin assembly in fast-treadmilling protrusions*. Biophys J 122(18):3816-3829. [10.1016/j.bpj.2023.08.022](https://doi.org/10.1016/j.bpj.2023.08.022)
  - **Finding:** FRAP + constrained modeling of steady treadmilling protrusions shows molecular diffusion alone cannot resupply G-actin to the leading edge; a vectorial anterograde transport (velocity <1 um/s, with diffusivity ~2 um^2/s) is required
  - **Relevance:** The 2023 Biophysical Journal transport/advection-limited monomer-resupply analysis; core evidence that leading-edge actin assembly is transport-limited (CFD-transport program)
- 📚 **Dessard, Manneville & Berret 2024** — *Cytoplasmic viscosity is a potential biomarker for metastatic breast cancer cells*. Nanoscale Adv 6(6):1727-1738. [10.1039/d4na00003j](https://doi.org/10.1039/d4na00003j)
  - **Finding:** Magnetic rotational spectroscopy: cytoplasmic viscosity 10-70 Pa*s and elastic modulus 30-80 Pa across MCF-10A/MCF-7/MDA-MB-231; MCF-7 viscosity ~5x higher than MDA-MB-231 (viscosity discriminates metastatic potential)
  - **Relevance:** Source of the physiological MCF7 cytoplasm viscosity setpoint (~65 Pa*s) the baseline HARD rule mandates
- 📚 **Moeendarbary, Valon, Fritzsche, ... Mahadevan & Charras 2013** — *The cytoplasm of living cells behaves as a poroelastic material*. Nat Mater 12(3):253-261. [10.1038/nmat3517](https://doi.org/10.1038/nmat3517)
  - **Finding:** Microindentation validates a poroelastic (biphasic) cytoplasm: a porous elastic solid meshwork bathed in cytosol, where short-timescale deformation rate is set by intracellular water redistribution through the mesh
  - **Relevance:** Physical basis for the Biot poroelastic CFD (fluid-solid) coupling in the FF cytoplasm/FSI module
- 🆕 **Novak, Slepchenko & Mogilner 2008** — *Quantitative analysis of G-actin transport in motile cells*. Biophys J 95(4):1627-1638. [10.1529/biophysj.108.130096](https://doi.org/10.1529/biophysj.108.130096)
  - **Finding:** Diffusion-drift-reaction modeling in realistic 3D cell geometry: rear-disassembly/front-assembly + diffusion set up a G-actin gradient transporting monomer forward; convective cytoplasmic flow contributes little; relaxation time ~10-100 s
  - **Relevance:** Quantitative diffusion-vs-advection G-actin transport model underpinning the transport-limited actin resupply analysis
- 🆕 **Vitriol, McMillen, Kapustina, Gomez, Vavylonis & Zheng 2015** — *Two functionally distinct sources of actin monomers supply the leading edge of lamellipodia*. Cell Rep 11(3):433-445. [10.1016/j.celrep.2015.03.033](https://doi.org/10.1016/j.celrep.2015.03.033)
  - **Finding:** Two molecularly distinct G-actin pools feed lamellipodia: cytosolic thymosin-beta4-bound monomer targeted to formins/leading edge (sets elevated G/F ratio for protrusion) and locally recycled lamellipodial F-actin
  - **Relevance:** Establishes non-homogeneous monomer pools + spatial targeting - constrains the transport-limited actin resupply model

## 1. Low-Reynolds physics

- 🆕 **Purcell 1977** — *Life at low Reynolds number*. American Journal of Physics 45(1):3-11. [10.1119/1.10903](https://doi.org/10.1119/1.10903)
  - **Finding:** Scallop theorem: at Re<<1 (bacteria ~1e-4-1e-5) reciprocal motion yields zero net displacement; inertia irrelevant, only shape-sequence non-reciprocity propels.
  - **Relevance:** Sets the inertialess Stokes regime that governs all intracellular cytosol/organelle transport and the force-balance (no inertia) assumption underlying the FF/DCM solvers.

## 2. Cytoplasmic streaming / advective transport

- 🆕 **Goldstein, Tuval & van de Meent 2008** — *Microfluidics of cytoplasmic streaming and its implications for intracellular transport*. PNAS 105(10):3663-3667. [10.1073/pnas.0707223105](https://doi.org/10.1073/pnas.0707223105)
  - **Finding:** First quantitative advection-diffusion theory of rotational streaming (Chara/Nitella); mixing enhancement depends strongly on helical flow pitch and peaks at Nitella's exponential growth phase.
  - **Relevance:** Canonical model coupling motor-driven fluid entrainment to intracellular scalar transport — the physics the CFD/advective-transport channel must reproduce.
- 🆕 **Haraguchi, Tamanaha, ... Ito 2022** — *Discovery of ultrafast myosin, its amino acid sequence, and structural features*. PNAS 119(8):e2120962119. [10.1073/pnas.2120962119](https://doi.org/10.1073/pnas.2120962119)
  - **Finding:** Four Chara myosin XI cloned; XI-1/XI-2 motor domains 3.2x/2.8x faster than prior CcXI; chimeric XI-1 reaches 60 um/s (~10x fast skeletal-muscle myosin), explaining the ~70 um/s streaming; 2.8-A crystal structure implicates actin-binding loops.
  - **Relevance:** Upper-bound motor velocity that sets the streaming forcing amplitude in the advective-transport model; corrects the seed's 'Kashiyama 2021' to a verified DOI.
- 🆕 **Tominaga et al. 2013** — *Cytoplasmic streaming velocity as a plant size determinant*. Developmental Cell 27(3):345-352. [10.1016/j.devcel.2013.10.005](https://doi.org/10.1016/j.devcel.2013.10.005)
  - **Finding:** Chimeric high-/low-speed myosin XI-2 (swapping Chara vs human motor domains) accelerated/decelerated streaming and correspondingly enlarged/shrank the plant — causal link between streaming velocity and size.
  - **Relevance:** Direct myosin-XI actuation-of-streaming experiment; the seed 'Tominaga myosin XI' paper tying motor speed to advective transport output.
- 🆕 **Verchot-Lubicz & Goldstein 2010** — *Cytoplasmic streaming enables the distribution of molecules and vesicles in large plant cells*. Protoplasma 240(1-4):99-107. [10.1007/s00709-009-0088-x](https://doi.org/10.1007/s00709-009-0088-x)
  - **Finding:** Review: myosin-XI-driven organelle movement entrains cytosol to drive cyclosis; streaming velocity varies with cell type, developmental stage and species and shapes intracellular molecular gradients.
  - **Relevance:** Biological grounding for the streaming/advective-transport module — connects motor mechanics to emergent cytosol flow and gradient formation.
- 🆕 **van de Meent, Tuval & Goldstein 2008** — *Nature's microfluidic transporter: rotational cytoplasmic streaming at high Peclet numbers*. Physical Review Letters 101(17):178102. [10.1103/PhysRevLett.101.178102](https://doi.org/10.1103/PhysRevLett.101.178102)
  - **Finding:** Chiral helical wall forcing produces transverse Dean-like vortices; strongly enhanced lateral transport and longitudinal homogenization emerge above a critical transverse Peclet number (boundary-layer scaling).
  - **Relevance:** The advective-mixing theory paper: shows how a chiral flow field boosts transport — directly informs how streaming should be modelled as an actuated CFD transporter.

## 3. Galvanotaxis / electrotaxis

- 🆕 **Allen, Mogilner & Theriot 2013** — *Electrophoresis of cellular membrane components creates the directional cue guiding keratocyte galvanotaxis*. Current Biology 23(7):560-568. [10.1016/j.cub.2013.02.047](https://doi.org/10.1016/j.cub.2013.02.047)
  - **Finding:** Rules out ion-transport hypotheses; galvanotaxis is insensitive to membrane potential/ion fluxes but fails below pH 6 and slows with aqueous viscosity - i.e. lateral electrophoresis of charged membrane components is the primary sensing mechanism.
  - **Relevance:** Mechanistic basis for E-field actuation acting via surface-charge electrophoresis rather than currents - directly relevant to modelling field effects on the membrane.
- 🆕 **Cortese, Palama, D'Amone & Gigli 2014** — *Influence of electrotaxis on cell behaviour*. Integrative Biology (Camb) 6(9):817-830. [10.1039/c4ib00142g](https://doi.org/10.1039/c4ib00142g)
  - **Finding:** Review comparing electrotaxis studies across cell types; implicates electro-migration of surface receptors/ion channels and highlights EFs as a precisely controllable guidance cue vs chemical gradients.
  - **Relevance:** Comparative synthesis of electrotaxis for choosing field-actuation protocols; complements the mechanism papers.
- 🆕 **McCaig, Rajnicek, Song & Zhao 2005** — *Controlling cell behavior electrically: current views and future potential*. Physiological Reviews 85(3):943-978. [10.1152/physrev.00020.2004](https://doi.org/10.1152/physrev.00020.2004)
  - **Finding:** Authoritative review: endogenous DC fields (~40-200 mV/mm) exist in all developing/regenerating tissues; lays out physics and cellular mechanisms of field-directed cell behaviour.
  - **Relevance:** The canonical review framing physiological EF magnitudes and mechanisms for the field-actuation program's baseline conditions.
- 🆕 **Mycielska & Djamgoz 2004** — *Cellular mechanisms of direct-current electric field effects: galvanotaxis and metastatic disease*. Journal of Cell Science 117(Pt 9):1631-1639. [10.1242/jcs.01125](https://doi.org/10.1242/jcs.01125)
  - **Finding:** Review of galvanotactic mechanisms (Ca2+ push-pull, voltage-gated Na+ channels, surface charge, protein electrophoresis); strongly metastatic prostate/breast cells are more galvanotactic and can move oppositely to weakly metastatic ones.
  - **Relevance:** Cancer-cell electrotaxis mechanisms and the Nav-channel angle relevant to MCF7/MDA field-response modelling.
- 🆕 **Nuccitelli 2003** — *A role for endogenous electric fields in wound healing*. Current Topics in Developmental Biology 58:1-26. [10.1016/s0070-2153(03)58001-2](https://doi.org/10.1016/s0070-2153(03)58001-2)
  - **Finding:** Wounds instantly generate injury currents producing lateral fields of 40-200 mV/mm; removing the field slows healing ~25%, and clinical EF stimulation increases healing rate 13-50%.
  - **Relevance:** Quantitative endogenous-field magnitudes and the causal wound-field evidence underpinning the electrotaxis program.
- 🆕 **Zhao et al. 2006** — *Electrical signals control wound healing through phosphatidylinositol-3-OH kinase-gamma and PTEN*. Nature 442(7101):457-460. [10.1038/nature04925](https://doi.org/10.1038/nature04925)
  - **Finding:** Physiological-strength DC fields are the dominant directional cue in wound healing; genetic loss of PI(3)Kgamma abolishes electrotaxis while PTEN deletion enhances it, placing electrotaxis on the chemotaxis PI3K/PTEN axis.
  - **Relevance:** Foundational E-field-actuation-of-migration paper; anchors the electrotaxis signalling that field-actuation experiments perturb.

## 4. E/B field effects on filaments

- 🆕 **Gartzke & Lange 2002** — *Cellular target of weak magnetic fields: ionic conduction along actin filaments of microvilli*. American Journal of Physiology - Cell Physiology 283(5):C1333-C1346. [10.1152/ajpcell.00167.2002](https://doi.org/10.1152/ajpcell.00167.2002)
  - **Finding:** Proposes the polyelectrolyte microvillar F-actin bundle as the cellular interaction site for weak EM/magnetic fields via nonlinear cable-like cation conduction through condensed ion clouds (stochastic-resonance/Brownian-motor framing).
  - **Relevance:** Bridges topics 4 and 6 - a mechanistic hypothesis for how weak external fields could couple to charged actin, informing field-actuation modelling.
- 🆕 **Stracke, Bohm, Wollweber, Tuszynski & Unger 2002** — *Analysis of the migration behaviour of single microtubules in electric fields*. Biochemical and Biophysical Research Communications 293(1):602-609. [10.1016/S0006-291X(02)00251-6](https://doi.org/10.1016/S0006-291X(02)00251-6)
  - **Finding:** Single microtubules migrate toward the anode with electrophoretic mobility ~2.6e-4 cm2/Vs (~0.19 e- net charge per tubulin dimer, IEP ~pH 4.2); DC fields also reorient kinesin-gliding MTs.
  - **Relevance:** Quantitative MT charge/mobility for modelling electric-field actuation and orientation of the microtubule compartment.
- 🆕 **Tang & Janmey 1996** — *The polyelectrolyte nature of F-actin and the mechanism of actin bundle formation*. Journal of Biological Chemistry 271(15):8556-8563. [10.1074/jbc.271.15.8556](https://doi.org/10.1074/jbc.271.15.8556)
  - **Finding:** F-actin behaves as a polyelectrolyte: multivalent cations above a valence-dependent threshold bundle filaments (ionic-strength dependent), analogous to DNA condensation - a generic charge-driven mechanism.
  - **Relevance:** Establishes the fixed charge density of actin that any E-field-on-filament actuation model requires.
- 🆕 **Tang, Ito, Tao, Traub & Janmey 1997** — *Opposite effects of electrostatics and steric exclusion on bundle formation by F-actin and other filamentous polyelectrolytes*. Biochemistry 36(41):12600-12607. [10.1021/bi9711386](https://doi.org/10.1021/bi9711386)
  - **Finding:** Electrostatic (polycation) vs steric (PEG) bundling of F-actin, DNA, microtubules and vimentin have opposite dependence on ionic strength and filament concentration - the shared polyelectrolyte nature dictates behaviour.
  - **Relevance:** Generalizes the actin polyelectrolyte result across all cytoskeletal filaments; useful for a unified charge-based field-response treatment.

## 4. Filament rotation

- 🆕 **Nishizaka, Yagi, Tanaka & Ishiwata 1993** — *Right-handed rotation of an actin filament in an in vitro motile system*. Nature 361(6409):269-271. [10.1038/361269a0](https://doi.org/10.1038/361269a0)
  - **Finding:** Sliding force carries a right-handed torque: a fixed-front actin filament forms a left-handed superhelix before supercoiling, revealing axial rotation of the filament during actomyosin sliding.
  - **Relevance:** Seed control paper - filament rotation from motor torque (NOT external field), needed to disambiguate field-induced vs motor-induced filament reorientation.
- 🆕 **Sase, Miyata, Ishiwata & Kinosita 1997** — *Axial rotation of sliding actin filaments revealed by single-fluorophore imaging*. PNAS 94(11):5646-5650. [10.1073/pnas.94.11.5646](https://doi.org/10.1073/pnas.94.11.5646)
  - **Finding:** Single-fluorophore polarization imaging shows sliding actin rotates once per ~1 um, far exceeding the 72 nm actin helical pitch - myosin 'runs' skipping protomers, with each kick primarily axial.
  - **Relevance:** Quantifies motor-driven filament axial rotation (baseline vs any electro/magneto-orientation), and demonstrates single-molecule orientation readout.

## 5. Membrane electrical shielding / dielectrics

- 🆕 **Schwan 1957** — *Electrical properties of tissue and cell suspensions*. Advances in Biological and Medical Physics 5:147-209. [10.1016/b978-1-4832-3111-2.50008-0](https://doi.org/10.1016/b978-1-4832-3111-2.50008-0)
  - **Finding:** Foundational treatment of tissue/cell-suspension dielectrics and the alpha/beta/gamma dispersions - basis for the membrane-capacitance shielding and the Schwan equation for field-induced transmembrane potential.
  - **Relevance:** Origin of the membrane dielectric/beta-dispersion model that governs how much of an applied external field is shielded vs reaches the cytoplasm.

## 5. Membrane electrical shielding / electroporation

- 🆕 **Marszalek, Liu & Tsong 1990** — *Schwan equation and transmembrane potential induced by alternating electric field*. Biophysical Journal 58(4):1053-1058. [10.1016/S0006-3495(90)82447-4](https://doi.org/10.1016/S0006-3495(90)82447-4)
  - **Finding:** Experimentally validates the Schwan equation for AC-field-induced transmembrane potential; measured critical breakdown potentials 0.33-0.53 V and inferred cytoplasm resistivity 910-1100 Ohm-cm.
  - **Relevance:** Quantitative induced-transmembrane-potential and cytoplasm-conductivity values for an E-field-shielding/electroporation actuation model.
- 🆕 **Weaver & Chizmadzhev 1996** — *Theory of electroporation: A review*. Bioelectrochemistry and Bioenergetics 41(2):135-160. [10.1016/S0302-4598(96)05062-3](https://doi.org/10.1016/S0302-4598(96)05062-3)
  - **Finding:** Comprehensive transient-aqueous-pore theory of electroporation: transmembrane potential ~0.5-1 V drives pore formation; sets reversible vs irreversible thresholds and pore energetics.
  - **Relevance:** The reference theory for the field-strength threshold above which membrane shielding fails and the field couples directly into the cytoplasm.

## 6. Magnetic actuation

- 🆕 **Etoc, Lisse, Bellaiche, Piehler, Coppey & Dahan 2013** — *Subcellular control of Rac-GTPase signalling by magnetogenetic manipulation inside living cells*. Nature Nanotechnology 8(3):193-198. [10.1038/nnano.2013.23](https://doi.org/10.1038/nnano.2013.23)
  - **Finding:** Functionalized magnetic nanoparticles inside living cells act as displaceable nanoscopic hot spots that, under magnetic force, locally trigger Rac-GTPase signalling and remodel the actin cytoskeleton / cell morphology.
  - **Relevance:** The intracellular SPION actuation paper - demonstrates direct in-cytoplasm magnetic force application and downstream cytoskeletal response, the core of the field-actuation concept.
- 🆕 **Wang, Butler & Ingber 1993** — *Mechanotransduction across the cell surface and through the cytoskeleton*. Science 260(5111):1124-1127. [10.1126/science.7684161](https://doi.org/10.1126/science.7684161)
  - **Finding:** Magnetic twisting cytometry through integrin-bound beads: cytoskeletal stiffness rises in direct proportion to applied stress and requires intact microtubules + intermediate filaments + microfilaments (tensegrity).
  - **Relevance:** Founding magnetic-bead actuation + tensegrity result - directly supports the MT-compression/actin-tension tensegrity work in the codebase.

## 6. Magnetic actuation / bead microrheology

- 🆕 **Bausch, Ziemann, Boulbitch, Jacobson & Sackmann 1998** — *Local measurements of viscoelastic parameters of adherent cell surfaces by magnetic bead microrheometry*. Biophysical Journal 75(4):2038-2049. [10.1016/S0006-3495(98)77646-5](https://doi.org/10.1016/S0006-3495(98)77646-5)
  - **Finding:** Magnetic-bead microrheometer (forces up to 1e4 pN on 4.5 um beads) yields local surface shear modulus ~2-4e-3 Pa*m and apparent cytoplasmic viscosity ~2e3 Pa*s, with a ~7 um screened elastic-field cutoff.
  - **Relevance:** Provides bead-force ranges and local viscoelastic/cytoplasm-viscosity numbers to calibrate magnetic actuation and the cytoplasm drag/FSI parameters.
- 🆕 **Fabry, Maksym, Butler, Glogauer, Navajas & Fredberg 2001** — *Scaling the microrheology of living cells*. Physical Review Letters 87(14):148102. [10.1103/PhysRevLett.87.148102](https://doi.org/10.1103/PhysRevLett.87.148102)
  - **Finding:** Magnetic-bead microrheology collapses onto a single weak-power-law master curve g*(f) ~ f^x (soft glassy rheology), unifying cell mechanical response across drugs and conditions.
  - **Relevance:** The scaling/soft-glass framework for interpreting bead-actuation microrheology data against the FF/DCM viscoelastic response.

## 7. Amoeboid / bleb cytoplasmic flow

- 🆕 **Charras, Yarrow, Horton, Mahadevan & Mitchison 2005** — *Non-equilibration of hydrostatic pressure in blebbing cells*. Nature 435(7040):365-369. [10.1038/nature03550](https://doi.org/10.1038/nature03550)
  - **Finding:** Local cortex relaxation does not stop blebbing on the far side, so hydrostatic pressure does NOT equilibrate over ~10 um / ~10 s; cytoplasm is a poroelastic contractile network infiltrated by cytosol.
  - **Relevance:** Foundational poroelastic-cytoplasm evidence - directly underpins the Biot/FSI two-way coupling and spatial pressure-field work in the FF engine.
- 📚 **Liu et al. 2015** — *Confinement and low adhesion induce fast amoeboid migration of slow mesenchymal cells*. Cell 160(4):659-672. [10.1016/j.cell.2015.01.007](https://doi.org/10.1016/j.cell.2015.01.007)
  - **Finding:** Without focal adhesions and under confinement, mesenchymal cells spontaneously switch to fast amoeboid migration via either a local protrusion or a myosin-II cortical-instability global cortical flow; a phase diagram is set by confinement, adhesion and contractility.
  - **Relevance:** Confinement/adhesion/contractility phase diagram for the mesenchymal-amoeboid switch - directly relevant to modelling confined migration and cortical flow.
- 🆕 **Petrie, Koo & Yamada 2014** — *Generation of compartmentalized pressure by a nuclear piston governs cell motility in a 3D matrix*. Science 345(6200):1062-1065. [10.1126/science.1256965](https://doi.org/10.1126/science.1256965)
  - **Finding:** Actomyosin pulls the nucleus forward via nesprin-3, physically partitioning cytoplasm into fore/aft compartments and pressurizing the front to drive lamellipodia-independent lobopodial 3D migration; nesprin-3 knockdown equalizes intracellular pressure.
  - **Relevance:** Nucleus-as-piston compartmentalized-pressure mechanism - a load path (nucleus <-> cortex) and pressure-compartmentalization the internal-coupling/CFD model must represent.
- 📚 **Ruprecht et al. 2015** — *Cortical contractility triggers a stochastic switch to fast amoeboid cell motility*. Cell 160(4):673-685. [10.1016/j.cell.2015.01.008](https://doi.org/10.1016/j.cell.2015.01.008)
  - **Finding:** Increasing myosin-II activity reversibly converts cells to a fast 'stable-bleb' amoeboid mode; cortical contractility fluctuations trigger a stochastic switch and rearward cortical flow drives migration in adhesive and non-adhesive environments.
  - **Relevance:** Stable-bleb cortical-flow motility - the actuation-by-contractility regime the cortex/cytosol model should capture.
- 📚 **Stroka et al. 2014** — *Water permeation drives tumor cell migration in confined microenvironments*. Cell 157(3):611-623. [10.1016/j.cell.2014.02.052](https://doi.org/10.1016/j.cell.2014.02.052)
  - **Finding:** Osmotic Engine Model: polarized Na+/H+ pumps + aquaporins create leading-edge water inflow and trailing-edge outflow, propelling confined tumor-cell migration even when actin polymerization and myosin-II contractility are inhibited.
  - **Relevance:** Water-flux/osmotic actuation of migration - the advective/osmotic transport channel independent of the cytoskeleton, complementary to the CFD-transport program.
