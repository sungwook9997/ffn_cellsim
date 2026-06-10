---
id: J2_cadherin_integrin_switch
topic: junction / adhesion-switch
title: "J2 — The cadherin (cell–cell) ↔ integrin (cell–ECM) adhesion switch under mechanical tension"
scope: "Mechanistic + quantitative synthesis: how rising junctional/bulk tension shifts adhesion from E-cadherin (adherens junction) to integrin (focal adhesion), in epithelial & breast-cancer (MCF7) context. Cadherin↔integrin mechanical crosstalk, force/tension thresholds, timescales, per-bond and per-cell force numbers."
ffn_relevance: High
ffn_themes: [junction, FA/clutch, cadherin, integrin, mechanotransduction, EMT, tension-homeostasis, parameter-source, validation-oracle]
entities: [e-cadherin, n-cadherin, ve-cadherin, integrin, alpha5beta1, integrin-beta1, talin, vinculin, alpha-catenin, beta-catenin, p120-catenin, focal-adhesion, adherens-junction, actomyosin, rhoa, fak, src, mcf7]
methods: [molecular-tension-fret-sensor, magnetic-tweezers, afm-force-spectroscopy, traction-force-microscopy, molecular-clutch-model, dna-tension-probe, single-cell-force-spectroscopy]
measurables: [single-bond-rupture-force-pN, molecular-tension-pN, traction-stress-Pa, intercellular-force-nN, force-threshold-pN, adhesion-energy, vinculin-binding-threshold, tension-ratio]
keywords: [cadherin-integrin-crosstalk, adhesion-switch, tensional-homeostasis, mechanotransduction, EMT, cadherin-switch, catch-bond, molecular-clutch, cell-cell-vs-cell-matrix, breast-cancer-invasion]
tags: ["#junction", "#cadherin-integrin-switch", "#mechanotransduction", "#tension-homeostasis", "#molecular-clutch", "#EMT", "#parameter-source", "#validation-oracle"]
has_transferable_params: true
status: synthesis (literature-derived, no project PDF in flat library covers this directly)
---

# J2 — The cadherin ↔ integrin adhesion switch under mechanical tension

**Tags:** #junction #cadherin-integrin-switch #mechanotransduction #tension-homeostasis #molecular-clutch #EMT #parameter-source #validation-oracle

> Provenance note. The flat reference library (papers #1–49) contains FA/clutch + integrin parameter sources (#13 curved-ECM integrin clutch, #29 catch-slip FA cluster oracle, #1 spreading active-gel) and MCF7/breast-cancer spheroid biology (#7, #15, #17), but **no paper directly on the cadherin↔integrin switch**. This file is therefore a literature synthesis built from PubMed + open-access full text, cross-linked to the in-library FA/clutch params. All numbers are cited with author/year/journal/DOI.

---

## 1. The core question and the one-sentence answer

Epithelial cells run **two competing adhesion systems** wired to the *same* actomyosin cytoskeleton: E-cadherin–based **adherens junctions (AJ)** (cell–cell) and integrin–based **focal adhesions (FA)** (cell–ECM). Both are force-activated catch-/clutch-type linkages whose maturation is gated by piconewton-scale forces. Because they **compete for a shared, finite pool of F-actin and the same force-sensitive adaptors (vinculin, α-catenin/talin)**, raising contractile/bulk tension does not simply strengthen one — it **reallocates the cytoskeletal "clutch budget."** Under EMT / spreading / stiff-ECM conditions the partition tips toward integrin–FA (and toward weaker N-cadherin), driving scattering and invasion; under soft-ECM / confluent conditions it tips toward cadherin–AJ. The switch is therefore best modeled as a **shared-resource (clutch-competition) + force-thresholded reinforcement** process, not a single hard threshold.

Authoritative conceptual anchors:
- **Mui, Chen & Assoian 2016**, *J Cell Sci* 129:1093–1100 — "The mechanical regulation of integrin–cadherin crosstalk organizes cells, signaling and forces." [DOI 10.1242/jcs.183699](https://doi.org/10.1242/jcs.183699) (PMC4813297). The crosstalk is **reciprocal and homeostatic**, not a clean antagonism.
- **Canel, Serrels, Frame & Brunton 2013**, *J Cell Sci* 126:393–401 — "E-cadherin–integrin crosstalk in cancer invasion and metastasis." [DOI 10.1242/jcs.100115](https://doi.org/10.1242/jcs.100115). Src/FAK/ILK + Rho-GTPase wiring of the switch in tumor cells.
- **Zuidema, Wang & Sonnenberg 2020**, *BioEssays* 42:e2000119 — "Crosstalk between cell adhesion complexes in regulation of mechanotransduction." [DOI 10.1002/bies.202000119](https://doi.org/10.1002/bies.202000119).

*(According to PubMed; DOIs linked above.)*

---

## 2. The mechanism, layer by layer

### 2.1 Shared cytoskeleton = a competition node
Both adhesions terminate on F-actin through **analogous force-sensing adaptor chains**:
- FA: integrin → **talin → vinculin** → actin.
- AJ: E-cadherin → β-catenin → **α-catenin → vinculin** → actin.

Vinculin and actin are the shared, rate-limiting currency. **Barcelona-Estaje et al. 2024**, *Nat Commun* 15:8824 [DOI 10.1038/s41467-024-53107-6](https://doi.org/10.1038/s41467-024-53107-6) (PMC11479646) make this explicit and quantitative with a **modified molecular-clutch model**: "talin–vinculin and α-catenin–vinculin would compete with each other to bind to actin, given the limited availability of actin filaments." Engaging cadherin (HAVDI N-cadherin peptide) on a substrate also bearing integrin ligand (RGD) **weakens integrin adhesion**:
- single-cell force spectroscopy: RGD-only adhesion force **up to ~6 nN**; adding HAVDI **lowers** the detachment force;
- focal-adhesion length **shrinks** when HAVDI is present, but **only above a substrate-viscosity/stiffness threshold** where the clutch is actually engaged (gel-phase DPPC ~1×10⁻⁴ Pa·s·m and glass, not fluid DOPC ~1×10⁻⁶ Pa·s·m);
- YAP nuclear translocation (mechanotransduction read-out) **drops** with cadherin ligation.
Interpretation: cadherin homophilic bonds are individually **weaker** than integrin–ECM bonds, so adding them creates "more, weaker bonds" that **siphon actin/vinculin away from FA**, capping integrin force transmission. **This is the cleanest mechanistic template for the ffn_cellsim switch: two clutch populations drawing on one actin pool.**

### 2.2 Force-thresholded reinforcement (the pN gates)
Each adhesion *matures* only once its adaptor crosses a force threshold that exposes a cryptic vinculin-binding site:
- **Talin (FA side):** the talin R3 rod domain unfolds at **~5 pN**, which is the **force threshold for vinculin binding and adhesion progression** — **Yao et al. 2014**, *Sci Rep* 4:4610 [DOI 10.1038/srep04610](https://doi.org/10.1038/srep04610) (PMC3980218). Vinculin binding then **locks talin open**, latching the FA.
- **α-catenin (AJ side):** a force-dependent conformational switch; the vinculin-binding M-domain unfurls reversibly at **~5 pN** (and irreversibly at 10–15 pN), and **5 pN triggers vinculin binding with nanomolar affinity**, converting transient force into a sustained junction-reinforcing signal — **Yao et al. 2014**, *Nat Commun* 5:4525 [DOI 10.1038/ncomms5525](https://doi.org/10.1038/ncomms5525).

So **both** adhesions share a ~5 pN per-molecule "arming" threshold. Whichever adhesion sustains supra-threshold load (set by where the actomyosin contractility is anchored and how stiff/ligand-dense the opposing partner is) **wins reinforcement**; the other relaxes and disassembles. Tension is the *selector*, not just an amplifier.

### 2.3 Tensional homeostasis / force balance (the cell-scale law)
At the whole-cell scale the two systems are **load-coupled and conserved**:
- **Maruthamuthu et al. 2011**, *PNAS* 108:4708–4713 [DOI 10.1073/pnas.1011123108](https://doi.org/10.1073/pnas.1011123108) — for an epithelial (MDCK) cell pair, the **cell–cell (cadherin) force ≈ 100 ± 40 nN**, directed **perpendicular** to the contact (88 ± 18°), and crucially the **intercellular force is a constant fraction of total cell–ECM traction, ratio ≈ 0.47 ± 0.07.** Increasing traction (e.g., stiffer ECM) **proportionally raises** junctional tension. → the systems do **not** trade off linearly bond-for-bond at the cell scale; they co-scale toward a homeostatic set point, while the *molecular* clutch competition (§2.1) sets *where* the actin goes.
- **Mertz et al. 2013**, *PNAS* 110:842–847 [DOI 10.1073/pnas.1217279110](https://doi.org/10.1073/pnas.1217279110) — in keratinocyte colonies, cadherin junctions **mechanically integrate the colony**: total traction grows with colony size and **localizes to the colony periphery** (interior cells offload to neighbors via cadherin). Loss of cadherin coupling redistributes traction back under every cell. → cadherin engagement **redirects** integrin traction spatially (rim-ward), it does not merely add to it.
- **E-cadherin molecular tension** itself is real and tunable: **Borghi et al. 2012**, *PNAS* 109:12568–12573 [DOI 10.1073/pnas.1117011109](https://doi.org/10.1073/pnas.1117011109) — a FRET tension sensor (sensitivity window **~1–6 pN**) shows E-cadherin is under **constitutive actomyosin tension (~1–2 pN per molecule)** even off-junction, **increasing** at contacts upon external stretch; requires αE-catenin + the catenin-binding domain.

---

## 3. The switch under specific perturbations

### 3.1 Rising bulk/substrate tension (stiff ECM, spreading) → toward integrin
- Cells **tune intracellular tension proportionally to ECM stiffness** (Mui 2016). On stiff polyacrylamide, FA-protein phosphorylation rises; **the FAK–p130Cas–Rac pathway is activated and *increases N-cadherin expression*** (a softer, more mesenchymal cadherin) (Mui 2016, citing Assoian lab work).
- Stiff ECM + spreading → larger FAs, more integrin clutch engagement, **partial AJ disassembly** ("focal adhesions disappear underneath cell–cell contacts; the effect depends on substrate stiffness and spreading extent" — Mui 2016). Switch direction: **cadherin → integrin.**

### 3.2 Soft ECM / confluence → toward cadherin
- On soft substrates (e.g., ~1.2 kPa) epithelial cells show **stronger E-cadherin / β-catenin membrane localization** and weakened ECM adhesion → cells aggregate (Mui 2016 and references therein). Switch direction: **integrin → cadherin.**
- E-cadherin rigidity sensing has its own stiffness transition during initial cell–cell contact (Bazellières-adjacent literature; **PNAS 2017** "Changes in E-cadherin rigidity sensing regulate cell adhesion," PMC5530647): % attachment on E-cadherin-Fc gels rose from **57% → 83%** as modulus increased.

### 3.3 EMT / breast cancer (MCF7 context) → integrin up, cadherin down & "softened"
- The classic **cadherin switch**: EMT lowers/loses **E-cadherin** and upregulates **N-cadherin** — a *qualitatively different, weaker, pro-migratory* cadherin — concomitant with **β1-integrin activation** and higher traction. **Canel 2013** [DOI 10.1242/jcs.100115](https://doi.org/10.1242/jcs.100115) reviews the Src/FAK/ILK/Rho wiring; loss of E-cadherin function de-represses integrin/FAK signaling and motility.
- **Figueiredo et al. 2021**, *Gastric Cancer* 25:124–137 [DOI 10.1007/s10120-021-01239-9](https://doi.org/10.1007/s10120-021-01239-9) (PMC8732838): E-cadherin dysfunction **raises traction forces and activates integrin β1**; **integrin β1 synergizes with E-cadherin loss** to drive scattering/invasion; clinically, **low E-cadherin + high integrin β1 = worse grade/survival.** Direct evidence the switch is *causal* for invasion.
- MCF7 specificity: MCF7 is an epithelial, E-cadherin-**positive**, relatively non-invasive luminal line (contrast MDA-MB-231, mesenchymal/E-cad-low). In MCF7 the balance sits on the **cadherin** side; pushing it toward integrin (stiff ECM, EMT induction, E-cad knockdown) is what unlocks dispersal. (Library refs #7, #15 use exactly the MCF7 vs MDA-MB-231 contrast; #15 profiles EMT/MMP/syndecan in these lines.)

---

## 4. Quantitative table — per-bond / per-molecule forces (the "adhesion energy/force" numbers requested)

| Quantity | Value | Source |
|---|---|---|
| **E-cadherin** single-bond unbinding (initial, weak state) | **< 12 pN** | Wang et al. 2016, *Sci Rep* 6:21584, [DOI 10.1038/srep21584](https://doi.org/10.1038/srep21584) |
| **E-cadherin** single-bond, matured/reinforced state | **> 43 pN** | Wang et al. 2016 (same) — shows force *escalation* during adhesion maturation |
| **E-cadherin** constitutive molecular tension (in cell) | **~1–2 pN/molecule** (sensor window 1–6 pN) | Borghi et al. 2012, *PNAS* 109:12568, [DOI 10.1073/pnas.1117011109](https://doi.org/10.1073/pnas.1117011109) |
| **N-cadherin** unbinding force, non-malignant cell | **26.1 ± 7.1 pN** | Lekka et al. 2011, *J Mol Recognit* 24:833, [DOI 10.1002/jmr.1123](https://doi.org/10.1002/jmr.1123) |
| **N-cadherin** unbinding force, malignant cell | **61.7 ± 14.6 pN** | Lekka et al. 2011 (same) — cancer *stabilizes* N-cadherin bonds |
| **Integrin α5β1–fibronectin** catch-bond optimum (lifetime ↑ with force) | **10–30 pN** window | Kong et al. 2009, *J Cell Biol* 185:1275, [DOI 10.1083/jcb.200810002](https://doi.org/10.1083/jcb.200810002) |
| **Talin** R3 unfolding = FA vinculin-binding / maturation threshold | **~5 pN** | Yao et al. 2014, *Sci Rep* 4:4610, [DOI 10.1038/srep04610](https://doi.org/10.1038/srep04610) |
| **α-catenin** unfolding = AJ vinculin-binding threshold | **~5 pN** (rev.); 10–15 pN irrev. | Yao et al. 2014, *Nat Commun* 5:4525, [DOI 10.1038/ncomms5525](https://doi.org/10.1038/ncomms5525) |
| **Cell–cell (cadherin) force**, epithelial pair | **~100 ± 40 nN** | Maruthamuthu et al. 2011, *PNAS* 108:4708, [DOI 10.1073/pnas.1011123108](https://doi.org/10.1073/pnas.1011123108) |
| **Intercellular force ÷ total traction** (constant fraction) | **0.47 ± 0.07** | Maruthamuthu et al. 2011 (same) |
| **In-library FA-clutch anchors:** single integrin–ligand bond spring k_LR | **~1.0 pN/nm**; eq. length **30 nm**; contact cutoff **300 nm** | library #13 (Table 1) — `13_rsc-ib-c2ib20159c-1-12.md` |

**Take-home force ordering:** matured integrin–ECM bonds (catch-bond ≥30 pN, FA clusters sum to nN/cell) generally **outcompete** individual cadherin bonds (E-cad ~12 pN initial, N-cad 26 pN normal). Both *individual* receptor classes "arm" at a **shared ~5 pN adaptor threshold**, but FAs build to far larger aggregate force because integrin catch-bonds *strengthen* over 10–30 pN whereas cadherin reinforcement is more modest unless cancer-stabilized (N-cad → 62 pN). This asymmetry is *why* rising tension biases the budget toward FA.

---

## 5. Tension thresholds & timescales (best current estimates)

**Thresholds**
- **Per-molecule arming:** ~5 pN (talin R3 / α-catenin) — universal gate for *both* adhesions (Yao 2014 ×2).
- **Integrin catch-bond reinforcement band:** 10–30 pN (Kong 2009) — the regime where ECM adhesion actively *wins* by lifetime extension.
- **Substrate-stiffness selector:** clutch only engages (and cadherin can then poach actin) **above a stiffness/viscosity threshold** (Barcelona-Estaje 2024: effect present on glass + gel-phase DPPC 10⁻⁴ Pa·s·m, absent on fluid DOPC 10⁻⁶ Pa·s·m). On the cell scale, soft (~1 kPa) favors cadherin; stiff (≥ tens of kPa) favors integrin (Mui 2016).
- There is **no single sharp "switch tension" in nN** in the literature — the transition is graded and set by the *ratio* of opposing ligand stiffness/density, consistent with the 0.47 constant-fraction homeostasis (Maruthamuthu 2011).

**Timescales**
- **Per-bond / conformational:** μs–s. Talin/α-catenin unfolding and vinculin capture happen within the lifetime of a loaded bond (sub-second to seconds at ~5 pN; Yao 2014). Integrin catch-bond lifetimes peak at **~seconds** in the 10–30 pN band (Kong 2009).
- **Adhesion (FA/AJ) remodeling:** **seconds–minutes** — FA assembly/disassembly and AJ vinculin recruitment respond to force on the **minutes** scale; FAs "disappear" under new cell–cell contacts over minutes (Mui 2016).
- **Spreading-driven shift:** spreading-phase mechanics resolve over **~3–10 min** (library #1 Betorz 2023: fast spreading P1 ~3 min, steady state ~10 min) — the window over which a spreading cell's integrin engagement ramps and can outcompete nascent junctions.
- **EMT / transcriptional cadherin switch (E→N):** **hours–days** — this is the slow, gene-expression-level arm (Snail/Zeb/Twist programs), distinct from the fast mechanical reallocation above. The mechanical switch can *precede and trigger* the transcriptional one (stiff-ECM → FAK–Rac → N-cadherin upregulation; Mui 2016).

---

## 6. How to model this in ffn_cellsim (fine-grained, no abstractions)

Per the project's "no-abstractions" rule (PI 2026-05-19), implement the switch as **two explicit particle-resolved clutch populations sharing one actin pool**, not a phenomenological switch variable:

1. **Two clutch types as bonds:** integrin–ECM clutches (already the H.4/H.5 FA track) **and** cadherin–cadherin clutches (new junction unit). Each is a Bell-Evans/catch-slip bond particle with its own force-dependent off-rate. Use the §4 numbers as per-bond anchors: E-cad k_off slip-dominated, low rupture (~12 pN); integrin α5β1 **catch-bond** (Pereverzev two-pathway, library #29 already adopts this family) peaking 10–30 pN; N-cad as a tunable stronger-when-malignant variant (26→62 pN).
2. **Shared actin / shared vinculin as the competition:** both clutch ensembles bind a **finite, explicit F-actin / vinculin pool** (Barcelona-Estaje 2024 mechanism). The "switch" then **emerges** from competition for binding sites, exactly the fine-grained behavior wanted — no hand-coded threshold.
3. **~5 pN maturation gate** as the adaptor unfolding event (talin/α-catenin) that *unlocks* extra actin-binding capacity for whichever clutch is loaded (Yao 2014 ×2) — implement as a force-gated state transition on the adaptor particle, not a global flag.
4. **Validation oracles (never runtime):**
   - cell-pair force balance: intercellular force ≈ 0.47 × total traction (Maruthamuthu 2011) — a clean emergent acceptance test.
   - colony-periphery traction localization (Mertz 2013).
   - stiffness-dependence: cadherin-favored < ~1 kPa, integrin-favored at high stiffness (Mui 2016); clutch-engagement threshold (Barcelona-Estaje 2024).
   - MCF7 vs MDA-MB-231 partition contrast (library #7, #15) for the EMT/breast-cancer arm.

**Caveats / open numbers.** (i) No literature gives a single "switch tension" in nN — it is ratio-/stiffness-set; encode it as competition, not a constant. (ii) Per-molecule *adhesion energies* (k_BT / ΔG) are rarely reported directly; the field uses rupture *forces* (§4) — convert via your bond potential, do not invent ΔG. (iii) The fast mechanical switch (s–min) and the slow transcriptional cadherin switch (h–d) are different processes; ffn_cellsim's single-cell mechanical scope captures the former, and should treat E→N cadherin identity change as a parameter swap, not a dynamic the mechanics produce on their own.

---

## 7. Source list (all PubMed-attributed)

According to PubMed and open-access full text:
1. Mui, Chen, Assoian 2016, *J Cell Sci* 129:1093 — [10.1242/jcs.183699](https://doi.org/10.1242/jcs.183699) (PMC4813297). Crosstalk review / homeostasis.
2. Canel, Serrels, Frame, Brunton 2013, *J Cell Sci* 126:393 — [10.1242/jcs.100115](https://doi.org/10.1242/jcs.100115). E-cad–integrin in cancer invasion.
3. Zuidema, Wang, Sonnenberg 2020, *BioEssays* 42:e2000119 — [10.1002/bies.202000119](https://doi.org/10.1002/bies.202000119). Adhesion-complex crosstalk.
4. Barcelona-Estaje et al. 2024, *Nat Commun* 15:8824 — [10.1038/s41467-024-53107-6](https://doi.org/10.1038/s41467-024-53107-6) (PMC11479646). N-cad/integrin clutch competition (the model template).
5. Barcelona-Estaje et al. 2021, *Adv Healthc Mater* 10:e2002048 — [10.1002/adhm.202002048](https://doi.org/10.1002/adhm.202002048). Adhesive-crosstalk biomaterials review.
6. Yao et al. 2014, *Sci Rep* 4:4610 — [10.1038/srep04610](https://doi.org/10.1038/srep04610) (PMC3980218). Talin R3 ~5 pN vinculin threshold.
7. Yao et al. 2014, *Nat Commun* 5:4525 — [10.1038/ncomms5525](https://doi.org/10.1038/ncomms5525). α-catenin ~5 pN vinculin switch.
8. Kong et al. 2009, *J Cell Biol* 185:1275 — [10.1083/jcb.200810002](https://doi.org/10.1083/jcb.200810002) (PMC2712956). Integrin α5β1 catch bond 10–30 pN.
9. Wang et al. 2016, *Sci Rep* 6:21584 — [10.1038/srep21584](https://doi.org/10.1038/srep21584) (PMC4753514). E-cad tension <12 pN → >43 pN.
10. Lekka et al. 2011, *J Mol Recognit* 24:833 — [10.1002/jmr.1123](https://doi.org/10.1002/jmr.1123). N-cad 26 pN (normal) vs 62 pN (malignant).
11. Maruthamuthu et al. 2011, *PNAS* 108:4708 — [10.1073/pnas.1011123108](https://doi.org/10.1073/pnas.1011123108). Cell–cell ~100 nN; intercellular = 0.47 × traction.
12. Mertz et al. 2013, *PNAS* 110:842 — [10.1073/pnas.1217279110](https://doi.org/10.1073/pnas.1217279110). Cadherin organizes colony traction (periphery).
13. Borghi et al. 2012, *PNAS* 109:12568 — [10.1073/pnas.1117011109](https://doi.org/10.1073/pnas.1117011109) (PMC3411997). E-cad constitutive ~1–2 pN molecular tension (1–6 pN sensor).
14. Figueiredo et al. 2021, *Gastric Cancer* 25:124 — [10.1007/s10120-021-01239-9](https://doi.org/10.1007/s10120-021-01239-9) (PMC8732838). Integrin β1 synergizes with E-cad loss → invasion.
15. Liu, Galior, Ma, Salaita 2017, *Acc Chem Res* 50:2915 — [10.1021/acs.accounts.7b00305](https://doi.org/10.1021/acs.accounts.7b00305) (PMC6066286). MTFM pN tension-probe methodology context.

**In-library cross-links (no DOI repeat):** `13_rsc-ib-c2ib20159c-1-12.md` (integrin-clutch params: k_LR ~1 pN/nm, l=30 nm, h_c=300 nm); `29_a-general-model-of-focal-adhesion-orientation-dynamics-in-re.md` (Pereverzev catch-slip FA-cluster oracle); `01_…spreading-migration-and.md` (spreading timescales ~3–10 min); `07_…tumor-microenvironment.md`, `15_…matrix-effectors.md` (MCF7 vs MDA-MB-231 EMT contrast).
