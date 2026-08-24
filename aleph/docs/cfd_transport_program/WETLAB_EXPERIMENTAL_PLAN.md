# Wet-Lab Experimental Plan — Cytoskeletal Transport, Stress-Fiber Assembly, and Field Actuation

**Program:** ffn_cellsim single-cell mechanobiology — bench validation campaign
**Companion:** fine-grained HOOMD/Warp simulation (ffn_sim `ff/` Filament-FEM engine); in-silico predictions are pre-registered per aim and used as the *leading model* below. In-silico is a companion, **not** a substitute — every hypothesis is decided at the bench.
**Document status:** proposal draft (v1). Numbers are targets for power/feasibility; final SOPs to be locked before first run.

---

## 0. Hypotheses under test

| ID | Hypothesis (leading / mechanistic model) | Competing null |
|----|------------------------------------------|----------------|
| **H1** | Cytoskeletal responses to a DC electric field are **signaling-gated** (PI3K/PTEN → polarity → actomyosin), not a direct electrostatic torque on filaments at physiological ionic strength. | EF acts by **direct electro-mechanical force** on charged filaments/membrane; blocking signaling should not abolish the response. |
| **H2** | Local intracellular torque/force (magnetic) produces **local remodeling + local stirring only**; it does **not** generate directed long-range net transport of cytoplasm/monomer without a rectifying cytoskeletal template. | Imposed torque/gradient drives **directed net transport** (a "molecular pump" from stirring alone). |
| **H3** | G-actin delivery to fast protrusions is **advection/reaction co-limited** near the leading edge (measurable depletion zone), not pure Fickian diffusion from a well-mixed pool. | G-actin is **diffusion-limited from a well-mixed pool**; no depletion zone; formin/monomer titration only rescales rate, not spatial profile. |
| **H4** | Dorsal stress fibers assemble **de novo by formin (mDia1) at focal adhesions**, condense from arcs, and bear **prestress with semi-autonomous mechanics** (recoil ∝ internal tension, weakly coupled to neighbors). | SFs are passively templated; ablation recoil reflects only global cell tension, not fiber-autonomous prestress. |
| **H5** | The adhesion **clutch is biphasic in substrate stiffness** (retrograde flow ↓ and traction peaks at intermediate stiffness), set by motor–clutch load-and-fail dynamics (Chan–Odde / Bangasser). | Flow and traction are **monotonic** in stiffness (stiffer = always more traction, always less flow). |
| **H1-d1** | A DC field **can** misalign reconstituted actin *only* below a Debye-screening threshold (low ionic strength); at physiological ionic strength the effect **collapses**, so in-cell EF steering must be signaling-mediated (feeds H1). | Actin electro-orientation persists at physiological ionic strength and directly explains galvanotaxis. |

**Cross-aim logic.** Aim 6 (in-vitro, cell-free) isolates the *physics* of direct field action; Aim 1 (in-cell) isolates the *biology*. If Aim 6 shows Debye collapse at 150 mM AND Aim 1's PI3K/PTEN block abolishes steering, H1 is supported and H1-null rejected. Aims 3–5 constitute the transport/assembly/clutch backbone the field aims perturb. **Aim 7 (TFM) is the cross-cutting mechanical readout** that binds them all — its field on/off traction *kinetics* (gradual+drug-sensitive vs instantaneous) is the single cleanest signaling-vs-direct-force discriminator, and its traction field is the direct in-silico validation target (H8).

---

## 1. Shared methods, materials, and statistics

- **Cell models (rationale per aim):** MCF7 (program reference line; epithelial, weakly motile — conservative galvanotaxis/transport baseline), primary/ hTERT **human dermal fibroblast** (robust dorsal SFs + FAs — SF/clutch aims), **fish epidermal keratocyte** (*Cyprinus carpio* scale explants; canonical fast, steady galvanotaxis and retrograde flow — mechanism dissection), U2OS (SF/TIRF workhorse for FA-formin imaging). Each aim states its primary line + one orthogonal line for generality.
- **Labels:** mEmerald/EGFP-β-actin (low-expression, FACS-gated <1.5× endogenous), **PA-GFP-actin** (photoactivation), mDia1-GFP / F-tractin-mScarlet, paxillin-mApple (FA marker), Lifeact (verify no SF artifact vs phalloidin controls).
- **Substrates:** polyacrylamide (PAA) gels, tunable **0.3–150 kPa** (Young's modulus), fibronectin-conjugated (sulfo-SANPAH), fluorescent-bead traction layer (0.2 µm, 2-color for high-density TFM).
- **Microscopy:** spinning-disk confocal + TIRF (60–100× 1.4–1.49 NA); QFSM-grade EMCCD/sCMOS; environmental control 37 °C / 5% CO₂ / humidified; anti-drift (hardware focus lock; fiducials).
- **Analysis:** QFSM (Danuser lab pipeline) for speckle flow; PIV/optical flow cross-check; TFM by Fourier-transform traction cytometry (FTTC) + Bayesian regularization; blinded, automated ROI extraction; **pre-registered** analysis scripts.
- **Statistics:** effect sizes with 95% CI, mixed-effects models (cell = random effect, nested in biological replicate). Power target **80% at α=0.05**; N stated per aim from pilot effect sizes; **≥3 biological replicates** (independent passages/days) per condition unless noted. Report both per-cell and per-replicate means (avoid pseudoreplication).
- **Blinding/randomization:** condition-blinded image analysis; randomized acquisition order; vehicle-matched (DMSO) drug controls.
- **Biosafety/ethics:** BSL-1/2 as line requires; keratocyte scales from euthanized fish under approved animal protocol; nanoparticle (SPION) handling per institutional nano-safety SOP.

---

## Aim 1 — Galvanotaxis rig: signaling-mediated vs direct-force EF response

**Objective (tests H1, feeds H1-d1).** Determine whether a physiological-to-supraphysiological DC field steers migration, retrograde flow, and SF orientation via a **signaling cascade (PI3K/PTEN)** or via **direct electrostatic force** on the cytoskeleton.

**Cell model.** Primary: **fish keratocyte** (large, fast, steady, textbook galvanotactic — clean readouts). Orthogonal: **MCF7** (program line; expect weak but measurable directedness) and hTERT fibroblast (for SF reorientation, which keratocytes lack).

**Technique + apparatus.** Zhao/McCaig-style **galvanotaxis chamber**: glass-bottom channel (cross-section ~ 1 mm × 10 mm × 22 mm) sealed with coverslip + vacuum grease; **agarose–Steinberg's/PBS salt bridges** (2–4% agarose, ≥5 cm long) connecting the medium reservoirs to Ag/AgCl electrodes in beakers of Steinberg's solution — electrodes kept **remote** so electrolysis products, pH shifts, and metal ions never reach cells. Field delivered by constant-current source; **field measured in-chamber** by probe electrodes (mV across a known gap), not assumed from voltage. Perfusion + Peltier hold 37 °C / pH. Live imaging: phase for tracks; TIRF + QFSM for actin speckle flow; confocal z for SF.

**Perturbation.** DC field **0.1, 0.4, 1, 2, 4, 6 V/cm** (0.4–2 V/cm = physiological wound range; 4–6 V/cm supraphysiological to probe threshold/damage). **Field reversal** mid-experiment (repolarization latency is a signaling signature). Pharmacology, orthogonal arms:
- **PI3K:** LY294002 (10–50 µM) or wortmannin (100 nM); genetic: PTEN-null vs PTEN-reconstituted (Zhao 2006 predicts PI3Kγ⁻/⁻ and PTEN-null flip/abolish directedness).
- **Myosin II:** **blebbistatin 10–50 µM** (para-nitro-blebbistatin to avoid phototoxicity/photoinactivation).
- **Actin polymerization control:** low-dose cytochalasin D / CK-666 (Arp2/3) as specificity checks.

**Quantitative readout.**
- **Directedness (cosine of trajectory angle to field), cos θ**, and **directional velocity** (µm/min) over 60–120 min; tracked N ≥ **60 cells/condition** across ≥3 replicates.
- **Retrograde-flow speed** (QFSM/fluorescent-speckle) as a function of angular sector relative to field: expect asymmetry (cathodal-vs-anodal) under leading model. Keratocyte lamellipodial flow baseline **~2–10 µm/min**.
- **SF reorientation** (fibroblasts): order parameter S = ⟨2cos²φ−1⟩ of SF long axis vs field, over 1–3 h; expected reorientation **perpendicular** to field at ≥2 V/cm.
- **Repolarization latency** after field reversal (min).

**Controls.** Zero-field (sham chamber, same salt bridges/heat), vehicle (DMSO) match, **field-only heat/pH control** (measure medium pH & temp at cells; must be within ±0.05 pH, ±0.3 °C of sham), photobleaching control for QFSM, drug wash-out reversibility, current-density audit (confirm no bubble/heat artifact at 6 V/cm).

**Expected result — leading (H1) vs null.**
- *Leading (signaling):* Directedness rises with field, saturating ~2 V/cm (cos θ ≈ 0.8–0.9 keratocyte; ≈0.2–0.4 MCF7). **PI3K inhibition / PTEN loss abolishes or reverses steering** while cells stay motile; blebbistatin decouples steering from speed; **repolarization on reversal takes minutes** (cascade timescale). SF reorientation is Rho/myosin-dependent (blebbistatin-sensitive).
- *Null (direct force):* Steering scales linearly with field with **no drug-sensitivity of direction**, and reversal is near-instantaneous (electrostatic). Retrograde-flow asymmetry would track field polarity independent of PI3K.

**Pitfalls.** Electrode products/pH/temperature masquerading as "field response" (mitigated by remote salt bridges + in-chamber pH/T logging); Joule heating at 4–6 V/cm; blebbistatin blue-light photoinactivation + phototoxicity (use nitro-blebbistatin, minimize 488 dose); keratocyte-to-MCF7 generalization gap (report per-line); QFSM speckle density/flow-vector validity at cell edges; drug off-target (confirm with genetic PTEN arm).

---

## Aim 2 — Magnetic actuation: directed transport vs local stirring

**Objective (tests H2).** Impose controlled intracellular torque/gradient and determine whether it produces **directed net transport** of cytoplasm/tracers or **only local stirring + local mechanotransductive remodeling** (SF growth at the loaded adhesion).

**Cell model.** hTERT fibroblast / U2OS (robust SFs + FAs, well-characterized MTC response). Orthogonal: MCF7.

**Technique + apparatus.** Three complementary loaders:
1. **Magnetic Twisting Cytometry (MTC):** RGD-coated **4.5 µm ferromagnetic beads** bound to integrins; magnetize (strong transient pulse), then apply twisting field (**up to ~50 G**, sinusoidal 0.1–10 Hz). Specific torque calibrated (~**17.5 Pa/G** bead stress constant); measure bead rotation/displacement (apparent stiffness) and downstream SF response (Wang/Ingber/Riveline paradigm).
2. **SPION + rotating field:** internalized **100–250 nm superparamagnetic iron-oxide** particles; rotating magnetic field (10–100 mT, 0.5–20 Hz) to impose local micro-torque/stirring.
3. **Magnetic tweezers:** single **~1–5 µm** paramagnetic bead, force **~10 pN–10 nN** with imposed gradient for directed pulling.

Imaging: confocal/TIRF with **tracer beads** (40–200 nm) + PA-GFP-actin/dextran for transport mapping; SF/paxillin live channels.

**Perturbation.** Torque amplitude/frequency sweep (MTC 5–50 G, 0.3–3 Hz); rotating-field on/off; tweezer force ramp 10 pN → 5 nN. Biochemical arms: **blebbistatin** (does local remodeling need myosin?), **Rho-kinase inhibitor Y-27632** and **formin inhibitor SMIFH2** (is FA-local SF growth formin/Rho-mediated?), cytochalasin D (does stirring transport change when the actin template is disrupted — isolates "template-rectified" vs "pure stirring" transport).

**Quantitative readout.**
- **Local tracer flow field** around the bead (PIV of 40–200 nm tracers): decompose into **rotational (curl) vs net-translational (divergence-free directed)** components; report net displacement of tracer centroid over 5–30 min at radii 1, 5, 10, 20 µm from bead.
- **SF growth at the loaded FA:** paxillin area growth + F-actin/mDia1 intensity vs distance from bead; force-induced FA reinforcement (Riveline: local growth in direction of applied stress).
- **Bead creep/stiffness** (MTC) as mechano-readout.
- N ≥ **30 beads/cells per condition**, ≥3 replicates.

**Controls.** Non-magnetized bead (same bead, no field), poly-L-lysine (non-integrin) bead (adhesion-independent control), field-only no-bead (heating/vibration artifact), tracer-only diffusion baseline (no torque → pure Brownian, sets divergence-free null), photobleach/drift correction, iron-dose cytotoxicity (viability + ROS assay for SPION).

**Expected result — leading (H2) vs null.**
- *Leading:* Tracer field is **dominated by local rotational (curl) flow that decays steeply (~1/r² or faster)** with **near-zero net centroid transport** beyond ~10 µm; the durable output is **local, formin/Rho-dependent SF growth at the loaded FA** (blebbistatin/Y-27632/SMIFH2-sensitive). Disrupting actin (cytochalasin) *increases* passive stirring reach but still yields no directed transport — confirming the cytoskeleton is a *template*, not a pump.
- *Null:* Measurable **directed net tracer displacement** persists at 10–20 µm, scaling with torque, independent of myosin — a stirring-driven pump.

**Pitfalls.** Bead heterogeneity (integrin number, internalization depth) → wide stiffness scatter (bin by binding); SPION endosomal confinement mimics "no transport" trivially (verify particle is cytoplasmic, not vesicle-trapped); rotating-field heating; distinguishing true directed transport from tracking bias/drift (rigorous divergence/curl decomposition + drift fiducials); phototoxicity over 30 min; SMIFH2 known off-targets (confirm with mDia1 knockdown arm).

---

## Aim 3 — Monomer transport: diffusion vs advection, and transport- vs reaction-limitation

**Objective (tests H3).** Measure G-actin transport mode (Fickian diffusion vs advective delivery) and whether protrusion growth is **transport-limited** (monomer supply) or **reaction-limited** (nucleation/elongation), by mapping depletion zones and titrating the monomer pool and formin activity.

**Cell model.** Fibroblast + keratocyte (fast protrusions with clear leading edge); MCF7 for reference.

**Technique + apparatus.** (i) **FRAP** of EGFP-β-actin: bleach a leading-edge ROI (2–5 µm), fit recovery to diffusion vs diffusion+flow models; apparent G-actin **D ~ 3–6 µm²/s**. (ii) **PA-GFP-actin photoactivation** pulse-chase: activate a defined spot (in cytoplasm and at edge), track centroid drift (advection velocity) vs spread (diffusion) over seconds–tens of seconds. (iii) High-frame-rate confocal/TIRF (**5–20 fps** for FRAP kinetics; 1–2 fps for chase). Ratiometric G/F-actin sensor optional.

**Perturbation.**
- **Monomer pool titration:** low-dose **latrunculin A (10–100 nM)** to graded-sequester G-actin (below the dose that stops protrusion); **profilin** and **thymosin-β4** manipulation (overexpression/knockdown) to shift the sequestered vs polymerization-competent pool.
- **Formin activity:** **SMIFH2 (10–30 µM)** to inhibit; **constitutively active mDia1 (mDiaΔN3 / ΔGBD-ΔDAD)** to boost barbed-end elongation demand at the edge.
- Arp2/3 arm (CK-666) to separate branched vs formin-fed pools.

**Quantitative readout.**
- **Depletion-zone depth & width** at fast protrusions: G-actin intensity profile from edge inward (µm), expect a supply-limited dip under leading model; correlate depletion with **protrusion velocity** (µm/min).
- **FRAP-derived D and advection v** (µm²/s, µm/s); fraction mobile.
- **Protrusion rate vs monomer availability** dose–response: transport-limited → protrusion rate saturates as pool is depleted with a *spatial* signature; reaction-limited → uniform rescaling, no depletion zone.
- N ≥ **40 protrusions / 20 cells per condition**, ≥3 replicates.

**Controls.** Free-GFP diffusion (upper-bound D, cytoplasmic viscosity reference), fixed-cell photobleach controls, expression-level gating (FACS), vehicle match, immobile-fraction correction, monomer-sensor calibration in vitro.

**Expected result — leading (H3) vs null.**
- *Leading (co-limited):* A **measurable G-actin depletion zone** appears at the fastest protrusions; boosting formin demand (CA-mDia1) **deepens** the depletion and reveals transport limitation (protrusion fails to scale with added elongation capacity); mild latrunculin **narrows protrusions from the tip** (spatial signature), not uniform slowdown; PA-GFP shows **net advective drift toward the edge** superimposed on diffusion.
- *Null (diffusion-limited, well-mixed):* **No depletion zone**; FRAP fits pure diffusion (no flow term); monomer titration and formin boost rescale protrusion rate **uniformly** with no spatial gradient.

**Pitfalls.** GFP-actin can perturb dynamics (validate against low-expression + phalloidin morphology); FRAP conflates diffusion with treadmilling turnover (use PA-GFP chase to separate); photodamage at high frame rate; depletion-zone contrast at the thin edge (TIRF geometry, background subtraction); latrunculin dose that silently stops protrusion (titrate with live viability of the edge); SMIFH2 off-targets (myosin) — cross-check with formin knockdown.

---

## Aim 4 — Stress-fiber assembly at focal adhesions: formin-driven growth, prestress, semi-autonomy

**Objective (tests H4).** Establish that dorsal SFs grow **de novo by mDia1 at FAs**, condense from transverse arcs, and carry **fiber-autonomous prestress** (recoil reflects internal tension, weakly coupled to neighbors), correlated with traction.

**Cell model.** **U2OS** and hTERT fibroblast (canonical dorsal-SF / arc / ventral-SF architecture; Hotulainen–Lappalainen paradigm).

**Technique + apparatus.** (i) **TIRF live imaging** of mDia1-GFP + F-tractin/Lifeact + paxillin-mApple to watch **de novo dorsal-SF elongation from FAs** and **arc condensation** (frame rate 0.1–1 fps over 30–120 min). (ii) **Femtosecond laser nano-ablation** (Kumar-style) to sever a single SF and record **elastic recoil** (retraction of cut ends). (iii) **Traction-force microscopy** (fibronectin PAA gel, 0.2 µm beads, FTTC) simultaneous/serial with ablation to correlate fiber tension with substrate traction.

**Perturbation.** **SMIFH2 / mDia1 knockdown** (blocks de novo dorsal SF growth — test formin-specificity); **blebbistatin** (removes myosin prestress — recoil should collapse if prestress is myosin-borne); **Y-27632** (Rho-kinase); calyculin A (hyper-contract, increases prestress). Ablate single vs paired adjacent fibers to test coupling (semi-autonomy).

**Quantitative readout.**
- **De novo growth rate** of dorsal SFs from FAs (µm/min) and mDia1 tip-tracking velocity; **arc condensation** kinetics.
- **Ablation recoil:** initial retraction velocity (**~0.1–0.5 µm/s**), plateau recoil distance (**~1–3 µm**), viscoelastic **recoil half-time (~5–20 s)** (Kumar 2006). Recoil amplitude ∝ prestress.
- **Semi-autonomy index:** displacement of the *adjacent, uncut* fiber upon ablation (near-zero → autonomous; large → coupled).
- **Traction correlation:** FA traction (Pa or nN) vs recoil-inferred fiber tension; expect positive correlation.
- N ≥ **30 fibers / 15 cells per condition**, ≥3 replicates.

**Controls.** Sham illumination at ablation power minus cut (photobleach/heat control), pre/post ablation morphology, blebbistatin recoil-collapse as positive prestress control, fixed-cell zero-recoil control, TFM bead-detachment/regularization checks, mDia1-GFP expression gating.

**Expected result — leading (H4) vs null.**
- *Leading:* Dorsal SFs **nucleate and elongate from FAs with mDia1 at the growing tip** (SMIFH2/mDia1-KD abolishes; Arp2/3 loss does not); arcs condense and fuse into ventral SFs; **ablation recoil is substantial and blebbistatin-sensitive**, adjacent uncut fibers barely move (**semi-autonomous**), and recoil amplitude **correlates with FA traction** (r > 0.5).
- *Null:* Recoil is uniform regardless of local fiber identity, adjacent fibers displace strongly (globally coupled tension), and formin inhibition does not selectively block dorsal-SF birth.

**Pitfalls.** Laser cut collateral damage / cavitation vs clean sever (titrate energy, verify single-fiber cut by morphology); recoil masked by adhesion drag; distinguishing dorsal vs ventral vs arc SF by geometry (blind classification); TFM spatial resolution vs single-fiber tension attribution; mDia1 overexpression artifacts; blebbistatin photoinactivation during 488 imaging (nitro-blebbistatin).

---

## Aim 5 — Retrograde flow ↔ clutch across substrate stiffness (biphasic motor-clutch)

**Objective (tests H5).** Map retrograde-flow speed and traction vs substrate stiffness to test the **biphasic clutch** (Chan–Odde / Gardel / Bangasser): flow decreasing and traction peaking at intermediate stiffness, set by motor–clutch load-and-fail.

**Cell model.** Fibroblast + U2OS on tunable gels; keratocyte for clean steady flow.

**Technique + apparatus.** **Simultaneous TFM + QFSM speckle microscopy** on fibronectin-PAA gels spanning **0.3, 1, 3, 8, 20, 50, 150 kPa**. TIRF speckle of low-level EGFP-actin for **retrograde-flow vector fields** (QFSM); 2-color bead TFM beneath for traction. Frame every **5–10 s over 15–30 min**. FA marker (paxillin) to register flow/traction to adhesion zones.

**Perturbation.**
- **Myosin:** blebbistatin (10–50 µM) dose–response and washout (motor arm of the clutch).
- **Adhesion:** integrin function-blocking antibody / RGD-peptide competition / graded fibronectin density (clutch-number arm); talin or vinculin knockdown as reinforcement arm.
- Actin polymerization rate (mild CK-666) to shift flow.

**Quantitative readout.**
- **Retrograde-flow speed** (µm/min) vs stiffness — expect monotonic decrease from **~5–10 µm/min (soft)** toward **~1–2 µm/min (stiff)** in lamellipodia, OR biphasic "load-and-fail" oscillation on soft gels.
- **Traction stress** (Pa) and total force (nN) vs stiffness — expect **peak at intermediate stiffness** (predicted optimum ~ few–tens kPa; Bangasser optimal-stiffness).
- **Flow–traction coupling** (local): inverse relation in the engaged-clutch regime (Gardel).
- N ≥ **25 cells per stiffness × condition**, ≥3 replicates; ≥3 gel batches (Young's modulus QC by AFM/rheometry).

**Controls.** Rigid glass reference (upper stiffness bound), bead-only no-cell TFM (noise floor), gel-modulus QC per batch, blebbistatin washout reversibility, speckle-density validity, flat-field/drift correction, matched fibronectin coating across stiffness (measure ligand density — a classic confound).

**Expected result — leading (H5) vs null.**
- *Leading (biphasic):* **Traction peaks at intermediate stiffness** and falls on both very soft and very stiff gels; **retrograde flow decreases** with stiffness (or shows load-and-fail oscillations on soft); **blebbistatin flattens both curves** (flow rises toward free-flow, traction collapses) — motor-dependence confirmed; reducing clutch number (integrin block / low FN) **shifts the traction optimum to higher stiffness**.
- *Null (monotonic):* Traction rises monotonically with stiffness; flow falls monotonically with no optimum; perturbations only rescale, never shift the optimum.

**Pitfalls.** **Ligand-density confound** (softer PAA can present different effective FN density — must measure and match); gel modulus drift/nonlinearity at extremes; QFSM flow accuracy on very soft substrates (large deformation); coupling TFM and speckle registration; cell-to-cell heterogeneity swamping the optimum (need N and stiffness resolution near the peak); blebbistatin phototoxicity; distinguishing lamellipodial vs lamellar flow regimes.

---

## Aim 6 — Field-induced filament misalignment (H1-d1): in-vitro electro-orientation threshold + Debye collapse, then in-cell test

**Objective (tests H1-d1, feeds H1).** In a **cell-free reconstituted** system, measure whether a DC field can electro-orient actin filaments, find the **field threshold**, and show it **collapses under Debye screening** at physiological ionic strength — then test in-cell whether comparable fields disrupt migration by direct action.

**Cell model / system.** **Reconstituted actin** (purified rabbit-muscle or recombinant human actin), **TIRF single-filament imaging** (Alexa/Atto-phalloidin or SNAP-actin), in a flow-cell electro-orientation chamber. In-cell arm: keratocyte/fibroblast in the Aim-1 rig.

**Technique + apparatus.** Microfabricated **electro-orientation flow cell** with in-plane electrodes (defined gap, e.g. 100–500 µm), DC/low-AC field; TIRF to visualize individual filaments and quantify angular distribution. Ionic strength tuned by buffer: **salt-free / low-I (1–10 mM)** vs **physiological (150 mM KCl/Mg-ATP)**. Debye length λ_D = 0.304/√I nm (I in M): **~9.6 nm at 1 mM → ~0.78 nm at 150 mM** — screening kills the electrostatic torque at physiological I.

**Perturbation.** Field sweep (**~1 to ~100 V/cm** in low-I; note this is *in-buffer* field, higher than in-cell because cytoplasm is conductive/screened); ionic-strength ladder (1, 10, 50, 150 mM); filament length classes (short vs long — torque ∝ length). In-cell arm: apply Aim-1 fields (0.4–6 V/cm) and look for *direct* (drug-insensitive, instantaneous) reorientation.

**Quantitative readout.**
- **Orientational order parameter** S = ⟨2cos²θ−1⟩ of filaments vs field, and **threshold field** E_th (V/cm) at which S rises above noise, as a function of ionic strength and filament length.
- **Debye-collapse curve:** E_th vs I — expect E_th to diverge (effect vanishes) approaching 150 mM.
- In-cell: SF/actin order parameter change per field, and its drug-sensitivity (cross-referenced to Aim 1).
- N ≥ **200 filaments per condition**; ≥3 preparations.

**Controls.** Zero-field angular distribution (isotropic baseline S≈0), flow-only alignment control (buffer flow can align filaments — decouple from field), electrolysis/heating/pH control at electrodes, photobleach control, filament-length calibration, no-electrode sham.

**Expected result — leading (H1-d1) vs null.**
- *Leading:* Electro-orientation is **real but only at low ionic strength** with high threshold fields; increasing I to physiological **collapses the effect (E_th → impractically high)**. Therefore in-cell galvanotaxis at 0.4–6 V/cm **cannot** be direct filament torque — the in-cell arm shows **no drug-insensitive, instantaneous reorientation**, consistent with H1 (signaling).
- *Null:* Actin orients at accessible fields even at 150 mM, and the in-cell arm shows fast, drug-insensitive reorientation — supporting the direct-force null and undercutting H1.

**Pitfalls.** Buffer-flow-induced alignment mistaken for electro-orientation (rigorous flow-off controls); electrode fouling/electrolysis in the tiny gap at high field; joule heating altering filament dynamics; surface tethering biasing angles (passivate; use crowding agent for free filaments); translating in-buffer field to in-cell field (cytoplasm screens — cannot naively compare V/cm); phalloidin stabilization altering charge/flexibility; single-filament tracking density.

---

## Aim 7 — Traction-Force Microscopy (TFM) as the cross-cutting mechanical output

**Objective (tests H1/H3/H6/H7/H8, and validates Aims 1–2 & 4–5).** Use TFM as the **single unifying quantitative mechanical readout** — the force the cell actually delivers to its substrate — to (a) close the mechanistic loop on the field hypotheses with a *dynamical* signature, and (b) supply the direct validation target for the CFD engine's clutch + stress-fiber modules. The physical logic: every hypothesis in this plan ultimately makes a claim about **where, when, and how much traction** the cell generates. TFM measures exactly that, and its **temporal kinetics discriminate signaling-mediated from direct-force actuation** in a way that migration end-points cannot.

**Technique + apparatus.**
- **Substrate:** fibronectin-conjugated polyacrylamide gel (matched to the Aim-5 stiffness set; default **8 kPa** for the field arms), 0.2 µm 2-color fluorescent beads as the deformation sensor; or **PDMS micropost arrays** (direct per-post force) as an orthogonal cross-check.
- **Reference/relaxed state:** post-hoc **cell detachment** (trypsin / SDS) or micropost null to recover the unstressed bead lattice.
- **Inversion:** regularized / **Bayesian Fourier-transform traction cytometry (FTTC)**; report the regularization-parameter sweep. Co-register **paxillin-mApple** for per-FA traction.
- **Temporal resolution 5–30 s/frame** — fast enough to separate an *instantaneous* (direct-force) from a *minutes-scale* (signaling) traction response.

**Perturbations & what TFM measures in each.**
1. **Galvanotaxis (H3/H10 discriminator).** Apply the Aim-1 DC field (0.4–6 V/cm) on the TFM gel and track the **traction-polarity vector** (first moment / net contractile dipole) as the field turns on, reverses, and turns off. *Signaling-mediated* → traction repolarizes **gradually over minutes**, is **abolished by PI3K/PTEN inhibition (LY294002)**, and relaxes over minutes when the field is removed. *Direct-force* → would repolarize **within seconds, drug-insensitively** (but H3 predicts this is impossible after membrane shielding). **This on/off kinetic is the cleanest single discriminator in the whole plan.**
2. **Magnetic actuation (H4/H13).** During SPION torque/gradient (Aim 2), measure whether traction changes are **local** (stress-stiffening near beads — the tensegrity signature, Wang–Ingber) or **global**, and confirm that stirring produces **no net directed traction** (no locomotor force), consistent with "local mixing, not propulsion."
3. **Field-induced SF misalignment (H6).** After imposing filament/SF disorder (Aim 6 in-cell arm), quantify the drop in **organized** traction: strain-energy change and **traction anisotropy** (dipole eigenvalue ratio) should fall as directed force generation degrades, even if total scalar traction is partly preserved.
4. **Baseline validation (H7, ties Aims 4–5).** Anchor SF tension (Aim 4 ablation recoil) and clutch traction (Aim 5 biphasic curve) to the same TFM pipeline so the mechanical readouts are cross-consistent.

**Readouts.**
- **Total strain energy** U = ½∫ **t**·**u** dA (pJ) — the most robust scalar; expect ~**0.01–1 pJ/cell**.
- **RMS/peak traction stress** (Pa) and **total contractile force** (nN; MCF7/fibroblast ~**10–100 nN**).
- **Contractile-moment tensor** M_ij → **net force direction + anisotropy** (the polarity metric for field arms).
- **Per-FA traction** (Pa or nN) co-registered to paxillin.
- **Response kinetics:** traction-vector reorientation **half-time** after field on/off (the signaling-vs-direct discriminator).

**Controls.** Detached-cell zero-traction reference; **bead-only, no-cell** noise floor; regularization-parameter sensitivity; gel-modulus QC per batch; **electrode-sham** (leads on, field off) for the galvanotaxis arm; **LY294002 / blebbistatin** pharmacology; drift/flat-field correction; matched FN density across conditions; SPION cytoplasmic-localization confirmation (Aim 2 shared).

**Expected result — leading vs null.**
- *Leading:* Field-driven traction repolarization is **gradual (minutes), drug-sensitive, cathode-biased**, and reversible on field-off — a signaling signature (supports H3/H10). Magnetic stirring gives **local** traction perturbation with **zero net directed force** (H4/H13). SF misalignment **lowers traction anisotropy / strain energy** (H6). Baseline TFM reproduces the **biphasic** stiffness curve (H5) and **SF-tension ↔ FA-traction** correlation (H7).
- *Null:* Traction changes are **instantaneous and drug-insensitive** under the field (reopens direct-force actuation, contra H3); or magnetic stirring yields a **net directed traction** (would imply internal-flow propulsion, contra the scallop/momentum argument H5).

**In-silico mirror (H8).** The CFD engine's clutch + SF modules predict a **traction field**; TFM **strain energy, polarity, and per-FA magnitude** are the direct acceptance targets. The engine's prediction — that external fields alter traction **only** by reorganizing the SF/clutch machinery (slow, indirect), never by a direct cytoskeletal body force — is falsified precisely by the *instantaneous-traction* null above.

**Pitfalls.** Traction-inversion is **ill-posed** — over-regularization erases polarity, under-regularization amplifies bead noise (fix the regularizer by cross-validation, report sensitivity); **electrode products/pH** on the gel corrupting beads or cell state (perfusion, salt bridges); **gel nonlinearity** at high traction; separating a genuine traction-polarity change from **cell translocation** during the field (register to the cell frame, or use non-migrating confined cells); FTTC spatial resolution vs single-FA attribution; SPION-loaded beads perturbing local mechanics independent of the field.

---

## 7. Cross-aim decision matrix (go / no-go)

| Outcome pattern | Interpretation |
|---|---|
| Aim 6 shows Debye collapse at 150 mM **AND** Aim 1 PI3K/PTEN block abolishes steering | **H1 supported**, direct-force null rejected. EF steering is signaling-mediated. |
| Aim 6 electro-orientation persists at 150 mM **AND** Aim 1 steering is drug-insensitive/instant | **H1 rejected** — reopen direct-force model; re-parameterize ffn_cellsim membrane/filament electrostatics. |
| Aim 2 tracer field is pure curl, no net transport, SF growth formin-dependent | **H2 supported** — cytoskeleton is a template, not a stirring pump. |
| Aim 3 depletion zone + formin-boost fails to scale protrusion | **H3 supported** — transport co-limitation; feed advection term to in-silico monomer model. |
| Aim 4 recoil blebbistatin-sensitive + semi-autonomous + traction-correlated | **H4 supported** — formin-FA de novo SF with fiber-autonomous prestress. |
| Aim 5 traction peaks at intermediate stiffness, blebbistatin flattens | **H5 supported** — biphasic motor-clutch; lock optimal-stiffness parameter. |
| Aim 7 field-driven traction repolarizes **gradually (min), drug-sensitive, reversible** | **H3/H10 supported** — EF actuation is signaling-mediated, not a direct cytoskeletal force. |
| Aim 7 traction change is **instantaneous & drug-insensitive** under the field | **H3 rejected** — direct-force actuation is live; re-parameterize engine electrostatics + revisit shielding estimate. |
| Aim 7 magnetic stirring → **local** traction perturbation, **zero net directed force**; SF misalignment **drops traction anisotropy/strain energy** | **H4/H6/H13 supported** — stirring ≠ propulsion; disorder degrades directed force. Feed the traction field to the in-silico clutch/SF validation (H8). |

Each **wet-lab outcome is pre-registered against the ffn_cellsim prediction**; disagreements trigger a model-revision loop (surface to PI), never a post-hoc gate loosening.

---

## 8. Timeline, throughput, and resources (rough)

- **Phase A (months 0–4):** rig builds + calibration (galvanotaxis chamber, magnetic loaders, electro-orientation cell), PAA gel + TFM/QFSM pipeline validation, cell-line/label QC.
- **Phase B (months 4–12):** Aims 3, 4, 5 (transport / SF assembly / clutch) — the mechanistic backbone.
- **Phase C (months 10–18):** Aims 1, 2, 6 (field actuation) once backbone readouts are trusted.
- **Replication:** ≥3 biological replicates/condition throughout; blinded analysis; pre-registered scripts.
- **Key equipment:** TIRF + spinning-disk confocal with fs-laser ablation, QFSM/TFM analysis stack, constant-current galvanotaxis source + Ag/AgCl electrodes + agarose bridges, MTC/magnetic-tweezers/rotating-field rig, SPION + purified-actin reagents, tunable PAA gel fabrication + AFM/rheometer for modulus QC.
- **Personnel:** 1 imaging lead, 1 biochem/reconstitution lead, 1 device/field-engineering lead, shared analysis.

---

## 9. Risk register (top items)

1. **Field artifacts** (pH/heat/electrolysis) confounding Aims 1/6 → remote salt bridges, in-chamber logging, sham controls.
2. **Label perturbation** (GFP-actin, phalloidin) → low-expression gating, orthogonal labels, fixed-cell cross-checks.
3. **Ligand-density confound** in stiffness sweep (Aim 5) → measure and match FN density per gel.
4. **Off-target inhibitors** (SMIFH2, blebbistatin) → genetic (KD/CA) confirmatory arms.
5. **Single-fiber attribution** in ablation/TFM (Aim 4) → energy titration, blinded classification.
6. **SPION vesicle-trapping** trivially mimicking "no transport" (Aim 2) → confirm cytoplasmic localization.
7. **In-buffer vs in-cell field translation** (Aim 6) → do not directly compare V/cm across conductive media; interpret via screening physics.

*All quantitative targets (fields, bead sizes, frame rates, N, effect sizes) are pilot-derived estimates for power/feasibility and are to be finalized against Phase-A calibration before first production run.*