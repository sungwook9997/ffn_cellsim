# ffn_cellsim — Lab Meeting Presentation Guideline

*Ready to feed into Claude Design. Also usable as the PI's speaker guide. All figure paths are absolute and verified. Numbers are commit-grounded. The honest current state is: the fine-grained spheroid stack COMPACTS, does not yet collectively spread — presented as a rigorous, mechanistically-explained finding.*

---

## 1. Talk metadata

**Suggested titles (pick one):**
1. *"Every Filament, Explicitly: A Fine-Grained Mechanistic Simulator of MCF7 Spheroid Mechanics — and What It Tells Us About Spreading"*
2. *"From Lumped to Literal: Building a GPU Mechanistic Cell Simulator, and Falsifying Our Own First Positive Result"*
3. *"Why a Mechanistic Spheroid Compacts Instead of Spreads — A Fine-Grained Cell-Mechanics Story"*

**One-sentence thesis:** We built a fine-grained, GPU-resident mechanistic simulator in which every cytoskeletal filament, motor, adhesion clutch, and ECM cross-link is an explicit particle/bond — and using it we reproduced the *form* of the MCF7 spheroid-spreading law from mechanism alone, while rigorously localizing the *magnitude* (A/A₀ = 7–10) to an active driver the collective model cannot yet host.

**Target length:** ~20–25 min talk + 10 min Q&A. ~16 slides.

**THE single takeaway:** *Resolving cell mechanics at the filament/clutch level lets you explain — not just fit — where a behavior comes from. Our spheroid does not passively spread because it is mechanically non-wetting (Douezan S < 0) at physiological MCF7 values; the 7–10× spread requires an active protrusive driver, and the model told us exactly that rather than us assuming it.*

---

## 2. Narrative arc (the story)

We start from a deliberate architectural inversion: where most cell-mechanics models run published closed-forms (Chan-Odde, Pereverzev, Hill) *as the mechanism*, we made those closed-forms validation oracles only and resolve the actual physics fine-grained on the GPU. We built a large, lit-anchored, parity-validated mechanistic force library — Stam-Hocky myosin minifilaments, Rakshit E-cadherin catch-bonds, Pereverzev integrin clutches, SimuCell3D node-face contact, exact-volume turgor — and two engine generations (HOOMD reference → differentiable NVIDIA Warp DCM). Pointed at the PI's MCF7 spreading law A/A₀ = a + b/R + c/R², the center-based line reproduced the *form* with zero calibration (r² = 0.998) but bounded the *magnitude* as a structural limit, handing it to the fine-grained engine — which, when we looked honestly, **compacts rather than collectively spreads**, and we retracted our own "it spreads to 1.94" headline as a 7-cell peeling artifact. The payoff: this negative result is now *mechanistically explained* (non-wetting at physiological tension; spreading needs an active driver whose force scale is a deliberate, un-tuned open question) — a sharper scientific position than a fitted positive.

---

## 3. Slide-by-slide outline

---

**Slide 1 — Title**
- BULLETS: Talk title; "ffn_cellsim — fine-grained mechanistic MCF7 cell & spheroid mechanics"; Sungwook Yoon, Shin Lab, Dept. of Mechanical Engineering, KAIST; date.
- VISUAL: Full-bleed hero render — `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/ipc_sphere_clean_n1000_mesh.png` (clean 1138-cell spheroid, per-cell colored meshes + back-half cut). Title text overlaid on a darkened lower band.
- SPEAKER NOTES: "This is a ~1100-cell MCF7 spheroid where every cell is an explicit deformable mesh, not a point. Today I'll tell you how we built it, what it reproduced, and — honestly — what it does not yet do and why that's actually the interesting result."

---

**Slide 2 — The biological question**
- BULLETS:
  - MCF7 breast-cancer spheroids spread on functionalized substrates (iCVD pV4D4 + collagen-I)
  - Spreading follows an empirical law: **A/A₀ = a + b/R + c/R²** (R = initial radius)
  - 3 ligand conditions, integrin-β1 read-out: **Bare ≈ 7.2 · Pre ≈ 7.5 · Lam4 ≈ 10.0**
  - This law has **no published analog → genuinely novel**
  - Question: can we reproduce it *from mechanism*, with no fitting?
- VISUAL: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/layer2/figs/fig_layer2_pi_overlay.png` (platform curves vs PI experimental data) — or a clean recreated A/A₀-vs-R scatter with the three conditions color-coded and the median values annotated.
- SPEAKER NOTES: "This is our own lab's experiment — MCF7 spheroids, three integrin-β1 ligand presentations. Spreading scales as a + b/R + c/R²; b is a traction/curvature term, c a small-size penalty. Medians run 7 to 10×. The whole platform exists to reproduce that partition from explicit mechanism, never by fitting to these numbers — they're overlay only."

---

**Slide 3 — The core principle: mechanistic over lumped**
- BULLETS:
  - Hard rule (2026-05-19): pick the **fine-grained mechanistic** option over lumped/proxy — even at higher cost
  - v1 ran paper closed-forms **as the runtime**; v2 demotes them to **acceptance oracles only**
  - Enforced in code: oracles are *runtime-import-forbidden*
  - The ONLY sanctioned coarse-graining: ×40 mesoscale (~1,000 vs ~38,000 filaments/cell) — a hardware constraint, PI-ratified
- VISUAL: The wrong-vs-right table as a clean two-column "card stack":

  | Decision | Lumped (wrong) | Mechanistic (right) |
  |---|---|---|
  | Bond off-rate | Metropolis | **Bell-Evans** |
  | Motor F–V | Linear stall | **Hill (1938)** |
  | Myosin II | Single 2-head spring | **Stam-Hocky bipolar minifilament** |
  | Cadherin | Slip-only | **Rakshit catch-bond** |
  | Arp2/3 | Rigid 72° | **Angle-harmonic + thermal** |
  | Integrator | Euler-Maruyama | **Leimkuhler-Matthews BAOAB** |

- SPEAKER NOTES: "The design philosophy is one rule: never replace a mechanistic process with a lumped one. Where the previous project ran, say, a Bell-Evans off-rate as the mechanism, here Bell-Evans is only an oracle the emergent dynamics must match. The single exception is a ×40 filament coarse-graining, and we call that out explicitly as a hardware constraint, not a modeling convenience."

---

**Slide 4 — The physiological-baseline rule (a methods value the lab cares about)**
- BULLETS:
  - Initialize every parameter/IC at its real in-vivo value **from step 0**
  - Cytoplasm η = **65.9 Pa·s** (Dessard 2024) — *not* water
  - Baseline turgor Π₀ ≈ **40 Pa** — *not* a relaxed floppy shell
  - Anchored ECM — *not* un-anchored
  - It is **invalid** to run from a null baseline and then compare to real data — the baseline itself is wrong
  - Origin story: the "γ-floor" — Π₀=0 left the cortex un-pretensioned → myosin floored ~1000× under band
- VISUAL: Recreated "two-cell" before/after diagram — left: floppy unpressurized sphere labeled "Π₀=0, η=water (WRONG)"; right: turgor-pretensioned sphere with cortex under tension labeled "Π₀=40 Pa, η=65.9 Pa·s (physiological baseline)". Arrow between them labeled "the comparison is only meaningful from the right state."
- SPEAKER NOTES: "A rule we hold hard: you measure from the real physiological operating point. If you start a cell at water viscosity with zero turgor and then ask myosin to generate all the tension, it floors ~1000× low — and we learned that the hard way. So perturbations and measurements happen *from* the pressurized, real-viscosity baseline, not from a convenient zero."

---

**Slide 5 — The compute stack & the HOOMD → Warp pivot**
- BULLETS:
  - Python 3.13 · **HOOMD-blue 7.0.1** (fine-grained reference / parity oracle) · **NVIDIA Warp DCM** (going-forward engine)
  - DCM = Deformable Cell Model: each cell is a closed **triangulated icosphere** (nodes/faces/edges), not a point
  - Warp wins: **GPU-resident** (no per-step host sync, 30.3× single-cell), **differentiable** (autodiff vs FD = 4e-12), **parity-gated** kernel-by-kernel against committed HOOMD fixtures
  - Custom **Leimkuhler-Matthews BAOAB** integrator + M-SHAKE/Fixman constraints (~40× integrator-only, gbook A5000)
  - Hardware baseline: **RTX A5000 mandatory** for production; CPU = dev only
- VISUAL: Recreated 3-layer stack diagram (Python host → HOOMD reference / Warp DCM engine → A5000 GPU), with a side callout box "DCM cell = triangulated shell" showing a small icosphere wireframe with node dots and one face highlighted. Use the node-face mesh wireframe as the recurring visual motif.
- SPEAKER NOTES: "Two engine generations. HOOMD-blue is our fine-grained reference and standing parity oracle. The production engine is a Warp port — GPU-resident, so no per-step CPU round-trip, and crucially differentiable, so we can back-propagate to physical parameters. Every force kernel is bit-parity-checked against a committed HOOMD fixture before it counts. A cell here is a deformable mesh of ~160–640 nodes, not a particle."

---

**Slide 6 — The mechanistic component inventory (overview)**
- BULLETS:
  - A parity-validated, literature-anchored mechanistic **force library**:
  - Cortex (connected spanning actin mesh) · Myosin II (Stam-Hocky + Hill) · Cadherin catch-bond (Rakshit) · FA/integrin clutch (Pereverzev) · Turgor (exact volume) · Node-face contact (SimuCell3D) · Lamellipodium · Division/proliferation/necrosis
  - Each anchored to a named paper + KU; each kernel bit-parity-gated
- VISUAL: A 3×3 grid of "mechanism cards" (use the reusable card layout from §4) — each card = icon + mechanism name + literature anchor + a one-word STATUS chip (BUILT / VALIDATED / OPEN).
- SPEAKER NOTES: "Before results — here's the parts bin. Each of these is an explicit force law anchored to a specific paper and validated kernel-by-kernel. I'll spotlight three on the next slide, then get to what the assembled system does."

---

**Slide 7 — Three mechanism spotlights (the "right way" examples)**
- BULLETS:
  - **Myosin II — Stam-Hocky bipolar minifilament**: rigid ~700 nm backbone, ~10 heads/side, Bell-Evans head binding + Hill (1938) stepping — *not* a single 2-head spring
  - **Cadherin — Rakshit-2012 catch-bond**: explicit E-cadherin trans-dimer bonds; **de-cohesion EMERGES from bond rupture**, not a binary switch
  - **FA/integrin — Pereverzev two-pathway catch-slip clutch**: basal nodes grip substrate ligand; Chan-Odde survives only as an oracle
- VISUAL: Three side-by-side mechanism cards with small schematics: (1) bipolar minifilament = red rod with many small head springs reaching to actin; (2) two cell membranes bridged by explicit cadherin trans-dimer springs, one shown rupturing; (3) a basal node tethered by a clutch spring to a fixed ligand dot on a dish. Optionally back this slide with `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h3/figs/fig_h3_evidence_montage.png` as a small inset (cortex validation).
- SPEAKER NOTES: "Three concrete 'right-way' choices. Myosin is a multi-head bipolar minifilament with real Hill force-velocity. Cadherin is explicit trans-dimer catch-bonds — so de-cohesion is an *emergent* consequence of bonds loading past their catch peak and rupturing, never a threshold flip. And focal adhesions are Pereverzev catch-slip clutches gripping the substrate."

---

**Slide 8 — Validation: emergent behavior matches the oracles**
- BULLETS:
  - Cortex: bending equipartition 0.99 kT; 3D Boltzmann-angle KS gate PASS; all beads within the 200 nm cortex band
  - Young–Dupré doublet: measured contact angle 11°→33° with cadherin strength, matches prediction (Δ2°)
  - Engine parity: force-laws to machine-eps; kT=0 bit-for-bit on CPU; differentiable loop autodiff vs FD = 4e-12
- VISUAL: Two-panel — left `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h3/figs/fig_h3_evidence_montage.png` (cortex network validation montage), right `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/young_dupre_doublet_gate.png` (Young-Dupré contact-angle PASS).
- SPEAKER NOTES: "The validation discipline: emergent measurements from the fine-grained dynamics must hit the literature closed-forms we held back as oracles. Cortex bending recovers equipartition to 0.99 kT; a two-cell doublet recovers the Young-Dupré contact angle to within 2 degrees. And the GPU port is parity-checked to machine precision against HOOMD."

---

**Slide 9 — Line A (center-based): the FORM is reproduced**
- BULLETS:
  - Lightweight 1-particle-per-cell spheroid on the same stack — to reach the *collective* term cheaply
  - **A/A₀ = 1.000 + 47.4/R + 97.5/R², r² = 0.998** — zero calibration, full PI R₀ range
  - Signs match the law: **b > 0** (traction/curvature), **c < 0** (small-size cohesion penalty)
  - Mechanistic origin: contact-inhibited proliferating rim ∝ 1/R + Rakshit catch-bond cohesion (catch was *necessary*: static Morse fragments, r²=0.80 FAIL → catch r²=0.98 PASS)
  - Ligand ordering **Lam4 > Pre ≳ Bare** reproduced (= PI 10 > 7.5 > 7.2)
- VISUAL: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/layer2/figs/fig_layer2_prod_rlaw.png` (production R₀-law, r²=0.998) — and a small inset of `/Users/sw1/ffn_cellsim/ffn_sim/outputs/layer2/figs/fig_layer2_ligand_production.png` (Bare/Pre/Lam4 ordering).
- SPEAKER NOTES: "First line of attack: a cheap center-based model to get the collective term. With zero calibration it reproduces the law's form at r²=0.998 across the full experimental radius range, with the right signs, and the right ligand ordering. The catch-bond is load-bearing — static cohesion fragments the growing spheroid; only a force-strengthening catch-bond holds it together. So the form is genuinely mechanistic, not fitted."

---

**Slide 10 — Line A: the MAGNITUDE is a structural limit (and where it points)**
- BULLETS:
  - Platform spreads only A/A₀ ≈ 1.1–2.2; PI is 7–10 → a real ~5–9× under-spread, at *matched* R₀
  - Falsified as a tunable miss on **6 independent axes**: scale ✗ · observable ✗ · cohesion ✗ · plithotaxis ✗ · passive wetting ✗ · active wetting ✗
  - Decisive: protrusion force (9.4 nN) **exceeds** cadherin cohesion (6.5 nN) → in a point-cell model it *tears rim cells off* instead of spreading the sheet
  - Verdict: A/A₀ 7–10 is the **active complete-wetting precursor-monolayer** regime — intrinsically shape-resolved → a point-cell cannot host it → hand to fine-grained
- VISUAL: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/layer2/figs/fig_layer2_magnitude_gap.png` (the ×5.9 force-ladder, gap localized to the active side) + small inset `/Users/sw1/ffn_cellsim/ffn_sim/outputs/layer2/figs/fig_layer2_decisive_traction.png` (protrusion-anchor ejection).
- SPEAKER NOTES: "But the magnitude — the 7-to-10× — the center-based model provably cannot reach. We ruled it out as a tuning miss on six independent axes. The decisive test: the lamellipodial protrusion force exceeds the cadherin cohesion, so in a one-particle model you just rip the edge cells off rather than spread a sheet. The conclusion is clean: the magnitude belongs to an active wetting regime that's intrinsically subcellular — which is exactly what motivates the fine-grained deformable engine."

---

**Slide 11 — Line B side-finding: the γ-floor (honest negative)**
- BULLETS:
  - Built a 3-channel cortical-tension estimator (active-myosin / constrained-backbone / passive-turgor, separated)
  - Result: active-myosin γ ≈ 0 — **~1000× below the band** (REFUTE); structural tension is ~100% passive turgor
  - Diagnosis: **generation-limited** (½·n·f·ℓ at lit density is ~13–36× under band) — *not* measurement, cohesion, formula, or buckling
  - Reframe: the band [0.35–0.65 mN/m] is a non-MCF7 (rounded HeLa/L929) proxy; real MCF7 ~0.27 sits below it — partly a *wrong target*
- VISUAL: Recreated bar chart — three γ channels (passive turgor ~0.08, constrained backbone, active myosin ≈ 0) with the literature band shaded, and the active bar visibly floored against it. Annotate "active-myosin ≈ 0, ~1000× under band."
- SPEAKER NOTES: "A parallel single-cell finding worth reporting honestly: when we let cortical tension emerge, the active-myosin channel floors at essentially zero — a thousand-fold under the canonical band. We traced it to generation-limit physics and a missing cortical-myosin-density datum, and reframed: the 'band' itself is a non-MCF7 proxy. This is a negative result we surfaced rather than tuned away."

---

**Slide 12 — Line B core result: the spheroid COMPACTS, does not collectively spread**
- BULLETS:
  - Single cell **DOES** spread (fried-egg, A/A₀ → 3.8, PI-confirmed)
  - The **spheroid does NOT** — converges to A/A₀ ≈ 1.0 across every condition:
    - penalty contact 1.04 · IPC (artifact-free) 0.996 · deformable+IPC 1.02 · **cadherin-OFF control 0.90 (monotone compacts)**
  - cadherin-OFF still compacts → no-spread is **structural, not a cohesion/contact artifact**
  - Lamellipodium self-disengages in the confluent pack
- VISUAL: Recreated convergence bar/dot plot: four conditions all landing at A/A₀ ≈ 1.0 with the PI target band (7–10) shown far above as a shaded "experiment" zone — visually conveys the gap. Optionally background it with the early-phase frame from `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/two_stage_n400.png` (the 4-panel storyboard).
- SPEAKER NOTES: "Here's the headline of the fine-grained line. A single cell spreads beautifully — fried-egg, A/A₀ to 3.8. But the spheroid does not collectively spread: across contact methods, deformability, even with cadherin entirely removed, it lands at A/A₀ ≈ 1.0. The cadherin-OFF control is the key — if removing cohesion still compacts, then no-spread is structural, not an adhesion artifact."

---

**Slide 13 — We falsified our own positive: the retractions**
- BULLETS:
  - We *had* claimed "spheroid SPREADS, A/A₀ → 1.94" — **RETRACTED**
  - Re-derivation: it was **7 basal-rim cells PEELING** (move 30–54 µm) while the other **93 stayed frozen** (median 0.44 µm); maxZ held 80→81 µm (no flattening). The 1.94 was a geometric-mean artifact
  - A 35-agent adversarial review independently confirmed the retraction
  - This is **one of three** measurement/visualization artifacts we caught and corrected
  - Durable lesson: trust cross-sections + scalar geometry, never the 3D surface render or A/A₀ alone
- VISUAL: Two-panel "before/after the honest look": left = the 3D surface render that *looked* like spreading; right = `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/cohesion_fix_peeling_vs_cohesive.png` (or `A1_actual_cells_3d.png`) showing the 7 mobile vs 93 frozen cells. Add a bold "RETRACTED" stamp on the left panel.
- SPEAKER NOTES: "And this is the slide I most want you to remember methodologically. We had a positive — the spheroid spreads to 1.94 — and we killed it ourselves. It was seven rim cells peeling off while ninety-three sat frozen; the number was a geometric-mean artifact, and maxZ never dropped so nothing actually flattened. An adversarial review confirmed the retraction. Falsifying your own first positive is the result."

---

**Slide 14 — WHY it doesn't spread: a mechanistic explanation**
- BULLETS:
  - **Douezan spreading coefficient S = w_cs − 2γ < 0 → NON-WETTING** at physiological MCF7 (γ ≈ 1e-2 N/m, w_cs = 2.85e-3 J/m²; S = −0.017)
  - Basal-only wetting can't flatten a 40 µm ball; a uniform shell just minimizes area → round
  - Aggregation finding: today's cells are too **RIGID** (V/V₀ = 1.000) vs the old two-stage physics (V/V₀ = 1.12, turgor-distended, junction-flattened); fix = **replicate the old physics** (turgor + cortical γ ON at lit value + remesh + long settle + N≥400), not a new builder
  - The missing piece is a **contact-line / Young-angle model + an ACTIVE protrusive driver** — not a knob
- VISUAL: Recreated explanatory diagram: a cell-on-substrate with the Douezan force balance (cortical tension γ pulling up/round vs substrate adhesion w_cs pulling flat) and the inequality S = w_cs − 2γ < 0 annotated "non-wetting." Side inset: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/faceting_comparison.png` or `RECONSIDER_xsec_surfcore.png` (rigid vs deformed packing).
- SPEAKER NOTES: "Crucially, no-spread is now *explained*, not just observed. At physiological tension the spreading coefficient is negative — the tissue is mechanically non-wetting, so it doesn't passively spread, full stop. Separately, our aggregated cells are too rigid: volume-locked spheres can't flow into flattened junctions the way turgor-pressurized cells do. The fix is to restore the real turgor and cortical-tension physics, and then add an *active* driver — because passive wetting provably can't get you to 7–10."

---

**Slide 15 — Methodology the lab can adopt**
- BULLETS:
  - **No magic numbers** — never tune a *derived* parameter to make a result pass; legit levers only: more time, add a missing real force at its derived value, fix a real bug, or surface the limit
  - **Physiological baseline** — measure from the real operating point
  - **Visualize anomalies for the PI, don't auto-fix** — most retractions came from *looking at the picture*
  - **A/A₀ = top-down silhouette, never contact area** (contact area inflated A/A₀ ~10–20×)
  - **De-cohesion must EMERGE** from catch-bond rupture, never a binary latch
- VISUAL: A clean 5-icon "principles" strip (one icon per principle), dark cards with the accent color. No photo needed.
- SPEAKER NOTES: "These are the working rules that produced the honest verdicts instead of faked passes. The one I'd flag for any modeler: never lower a derived parameter to make a gate pass — that's the magic-number sin. When a cad-bundle sweep started chasing 'the value that spreads best,' it got stopped. The legitimate moves are more time, a real missing force at its real value, or a real bug fix — otherwise you surface the limit."

---

**Slide 16 — Where we are & what's next**
- BULLETS:
  - **Built & validated:** a parity-gated mechanistic force library + two engines (HOOMD ref ↔ differentiable Warp DCM)
  - **Reproduced:** the spreading-law FORM (r²=0.998) + ligand ordering, from mechanism
  - **Explained (honest):** the spheroid is non-wetting at physiological values → does not passively spread; magnitude needs an active driver
  - **Open / PI-gated:** (a) active protrusion as the spreading DRIVER (blocked on a force-scale anchor — deliberately *not* auto-tuned); (b) contact-line/deformable-substrate model; (c) restore turgor+γ aggregation physics, re-measure N≥400
  - **TAKEAWAY:** mechanistic resolution lets us *explain* where a behavior comes from — and tells us exactly which missing physics matters
- VISUAL: Hero closer — `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/ipc_cleave_opt_n1000_mesh.png` (the ~1128-cell engine-scale render) with the three open forks listed as a clean right-side panel.
- SPEAKER NOTES: "To close: we have a validated mechanistic force library and a differentiable GPU engine; we reproduced the spreading law's form from mechanism; and we honestly explained why the spheroid doesn't yet spread — it's non-wetting at physiological tension, and the magnitude needs an active driver. That driver is the next decision, and we're deliberately not tuning our way into it. The takeaway: fine-grained mechanism doesn't just fit the data — it tells you which physics is missing."

---

*(Optional backup slides if time/Q&A: the two-stage aggregation movie `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/two_stage_n400_agg.mp4`; the de-cohesion montage `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/a2_decohesion_n100_montage.png`; necrosis 3-zone `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_gpu_lod/figs/necrosis_3zone_spatial.png`; division mechanism `/Users/sw1/ffn_cellsim/ffn_sim/outputs/warp_decohesion/figs/CLEAVE_mechanism.png`; throughput/scaling numbers.)*

---

## 4. Design / visual direction for Claude Design

- **Overall aesthetic:** Clean scientific, **dark theme** (near-black `#0E1116` background) — the 3D mesh renders sit on dark backgrounds and full-bleed beautifully. Single accent: a **cyan/teal `#34D2C8`** (reads as "mesh wireframe / GPU") with a warm **coral `#FF6B5C`** as the *secondary* accent reserved for "active/myosin/retraction" call-outs (mirrors the red lamellipodia/myosin in the renders). Font: a clean geometric sans (Inter / IBM Plex Sans) for body; a technical mono (JetBrains Mono / IBM Plex Mono) for numbers, equations, commit hashes, parameter names.
- **Consistent visual motif:** a thin **node-face triangulated-mesh wireframe** in the accent color — corner accent on every slide, literal subject of the engine/mechanism slides. Carry the A/A₀ = a + b/R + c/R² equation as a small recurring footer glyph on the results slides so the audience always knows which question a slide answers.
- **Figure-heavy slides:** the 3D renders (1, 16) and montages (8, 12, 13) should be **full-bleed** with text in a single darkened lower/side band — do not box them in a white card. One figure per "story slide"; montages already carry a multi-frame story so don't crop them.
- **Reusable "mechanism card" layout** (slides 6, 7): a vertical card — top: small line-art schematic on a dark tile; middle: **mechanism name** (sans, bold) + **literature anchor** (mono, dim, e.g. "Rakshit 2012 · KU-4.2"); bottom: status chip (`BUILT` teal · `VALIDATED` teal-filled · `OPEN` coral-outline). Identical dimensions so the 3×3 grid reads as an inventory.
- **Movies (.mp4) — embed/live where the venue allows:**
  - Opener loop: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/two_stage_n400_agg.mp4` (cells coalescing).
  - Slide 12 (compaction story): `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_two_stage/figs/two_stage_n800_spread.mp4` (frame the EARLY phase only — see honesty note).
  - Backup/Q&A: `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_gpu_lod/figs/junction_switch_spatial.mp4` (de-cohesion rotation), `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h_dcm_gpu_lod/figs/necrosis_3zone_spatial.mp4`.
  - If movies can't embed, use the static montages (`a2_decohesion_n100_montage.png`, `two_stage_n400.png`).

---

## 5. Honesty & framing guidance

- **Frame the retraction as the scientific spine, not an apology.** Script: *"We had a positive result — the spheroid spreads, A/A₀ → 1.94 — and we falsified it ourselves. On inspection it was seven rim cells peeling while ninety-three stayed put; the number was a geometric-mean artifact and nothing actually flattened. An adversarial review confirmed the retraction."* Land it as discipline: we trust cross-sections and scalar geometry over a pretty 3D render, and we caught **three** such artifacts the same way. Then pivot immediately to slide 14 — the no-spread result is *explained* (Douezan S < 0, non-wetting at physiological tension), a stronger position than the retracted positive.
- **Scientific value of this bounded/negative result:** A fitted positive would have told the lab nothing about mechanism; instead we have a *mechanistic prediction* — the spheroid is non-wetting at physiological MCF7 tension, so passive spreading is impossible and the experimental 7–10× must come from an active protrusive driver. That converts "our model doesn't match yet" into "here is exactly which physics is missing and why," and gives the experimentalists a testable claim (perturb cortical tension / substrate adhesion to cross the wetting threshold).

---

## 6. Anticipated Q&A

- **Q: Why not just a vertex model / Cellular Potts / phase-field?**
  A: Those are lumped at exactly the level we care about — they encode tension and adhesion as energies/parameters, so they can *fit* the spreading law but can't tell you it comes from cadherin catch-bonds rupturing or integrin clutch traction. Our whole point is to make those the explicit, falsifiable mechanism. The DCM is the minimal representation that still resolves shape, contact line, and per-junction bond mechanics.

- **Q: Isn't the ×40 filament coarse-graining itself a lumped model?**
  A: It reduces *count* (~38,000 → ~1,000 filaments/cell), not *mechanism* — per-filament stiffness and every force law are unchanged (the H.3 force-constant check shows ×40 touches filament count only). It's the single sanctioned coarse-graining, called out explicitly as a hardware constraint; everything else stays fully fine-grained.

- **Q: What's your actual evidence the cortex is "right"?**
  A: Emergent measurements hit the held-back oracles: bending equipartition recovers 0.99 kT, the 3D bond-angle distribution passes a Boltzmann KS gate, all beads sit within the 200 nm cortex band, and a two-cell doublet recovers the Young-Dupré contact angle to within 2°. The honest caveat: the *emergent active cortical tension* underperforms (the γ-floor, slide 11) — which we report, not hide.

- **Q: How do you know "compacts not spreads" is real physics and not a bug?**
  A: Three ways. (1) Robust across independent contact methods (penalty, IPC, projection), deformability, and N. (2) The cadherin-OFF control still monotonically compacts to 0.90 — remove cohesion and it *still* doesn't spread, so it's not an adhesion artifact. (3) It's mechanistically explained: the Douezan spreading coefficient S = w_cs − 2γ is negative at physiological values, independently predicting non-wetting. And we held it to that standard precisely *because* we'd already been burned by the peeling artifact.

- **Q: What would falsify your model?**
  A: If, at the genuinely physiological operating point (real turgor + cortical γ ON + remesh, N≥400) with an active driver at its *derived, un-tuned* force scale, the spheroid still failed to reproduce both the form *and* magnitude of A/A₀ = a + b/R + c/R² — that falsifies the claim this mechanism set is sufficient. Conversely, if experimentally lowering cortical tension or raising substrate adhesion across the wetting threshold (S > 0) did *not* induce spreading, that falsifies the Douezan explanation.

- **Q: Can you reach experimental N and timescales?**
  A: Spheroids at N = 400–1138 cells run on one A5000 (N=400 in ~3 h); throughput *grows* with cell count because the hash-grid turns O(N²) contact into O(N·k) (10.9× at 8 cells → 36.8× at 27). For time, the constrained integrator lifts the CFL ceiling ~5000×, an implicit solver adds ~2–3× net, and `--accel-real-hours` lets one run *represent* 24–48 h of biology. Native filament count (~38k/cell) is GPU-only and not the production target — that's what the ×40 mesoscale is for.

- **Q: You said it's differentiable — what does that buy you?**
  A: The whole K-step loop runs under autodiff (validated against finite-difference to 4e-12), so we can back-propagate a loss on the final cell state to physical parameters like turgor and edge stiffness — a path to principled parameter estimation / inverse design the reference HOOMD runtime can't do. Note we keep it disciplined by the no-magic-number rule: we don't autodiff a *derived* constant to make a gate pass.

- **Q: Why is the active-driver force scale "PI-gated" instead of just measured/fit?**
  A: Because the honest literature anchor (single-stress-fiber / NMII force scale) is missing a key motor-density datum, and back-solving it from the desired spreading magnitude would be exactly the magic-number sin. Choosing the anchor is a scientific judgment call we surface explicitly rather than tune to outcome.

- **Q: The single cell spreads to 3.8 but the spheroid doesn't — contradiction?**
  A: No — that's the central finding. A single cell has a free contact line and a substrate to react traction against; in a confluent pack the lamellipodium self-disengages and basal-only wetting can't flatten a 40 µm ball. The single-cell-vs-collective divide *is* the structural limit, and it's exactly why the magnitude was handed from the center-based line to the fine-grained engine.

---

## 7. Time budget (~22 min talk)

| Slides | Section | Minutes |
|---|---|---|
| 1–2 | Title + biological question | 2.5 |
| 3–4 | Principles: mechanistic-over-lumped + physiological baseline | 3.0 |
| 5 | Compute stack + HOOMD→Warp pivot | 2.0 |
| 6–8 | Mechanism inventory + spotlights + validation | 3.5 |
| 9–10 | Line A (center-based): FORM reproduced + magnitude structural limit | 3.5 |
| 11 | γ-floor honest negative | 1.5 |
| 12–13 | Spheroid compacts + the retraction | 3.0 |
| 14 | WHY (Douezan non-wetting + rigidity) | 1.5 |
| 15 | Methodology principles | 1.0 |
| 16 | Where we are + next + takeaway | 1.0 |
| — | **Total** | **~22.5** |

*Buffer: trim slide 6→merge into 7, and slide 11 to 1 min, if running long. Slides 13–14 are the intellectual core — protect their time even if you cut elsewhere.*
