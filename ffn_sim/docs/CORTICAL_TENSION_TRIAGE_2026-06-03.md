# Cortical / membrane / tissue-tension literature triage — γ-floor diagnosis + single-cell→spheroid bridge (2026-06-03)

> **STATUS: TRIAGE scaffold — NOT a contract, NOT ratified bands.** Classifies the
> ~11 PI-supplied cell-tension papers (dropped into `ffn_sim/references/` 2026-06-03
> 15:54–15:59) against (a) the months-long KU-3.5 cortical-tension γ-floor and (b) the
> PI thesis that **single-cell cortical tension γ is the root of spheroid/aggregate
> surface tension**. Every item enters the platform **only** as an acceptance oracle,
> a literature constant/anchor, or an observable target — **never a runtime mechanism**
> (CLAUDE.md inversion rule, PI 2026-05-19). Author: Lead session on direct PI request
> ("references에 cell tension 논문들 넣었어 — 확인해봐"). Companion to
> [`v2_audit/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md`](v2_audit/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md),
> the STAGE-2 closeout (commit `7e416c1`), and the Layer-2 line
> ([`LAYER2_ROADMAP_2026-06-03.md`](LAYER2_ROADMAP_2026-06-03.md), §D2 surface-tension).

PI scope decision (2026-06-03): the two walls are **one physics** — *"둘 다 — 단일세포 γ가
spheroid 응집의 근원."* Solve the single-cell cortical tension and bridge it up to the
aggregate surface tension. This doc is organised around that unification.

---

## 0. The current wall state (latest on-disk reality — read this first)

The γ-floor diagnosis has moved **past** the `KU35_FLOOR_ROOT_CAUSE` doc. STAGE-1 + STAGE-2
(autonomous 2026-06-03, commits `5bfb146`, `9c7826e`, `02833f2`, `7e416c1`) established:

- **Two measured channels.** `g_rigid` = myosin-independent **structural/passive** cortical
  tension, native **0.57 mN/m — already IN the KU-3.5 band [0.35, 0.65]**. `g_soft` = the
  **active** channel (×340 motor response — confirmed the true active readout), floored at
  **~2×10⁻⁴ mN/m, ~1000× under band**.
- **STAGE-2 eliminated four levers** as the cause of the active floor — each FAILS to move
  `g_soft`: **percolation/connectivity** (z≈3), **turnover** (f_ss=0.591, chain-split),
  **measurement channel** (motors ON/OFF), and **coherence** (meridional alignment verified
  `|axis·meridian|=1.000`). ⇒ **The floor is UPSTREAM of transmission: force GENERATION /
  the soft-coupling CEILING (STAGE-1) / model fidelity — NOT connectivity, NOT √N dipole
  coherence, NOT transmission.**
- **The head-tension smoking gun** (`02833f2`): bound heads ARE loaded to **~44% of stall**
  (F/F_stall≈0.44), **but** the grip accumulator `s_grip≈0` (heads not actually walking) and
  **only ~3 % of heads are bound** (9–16 of 400). ⇒ Per-head force exists but **does not
  aggregate into network tension** — too few engaged motors + no sustained material transport.

**Consequence for this triage:** the cell-tension literature must be read against *this*
narrowed wall, not the older "coherence/cancellation" framing. The c3fafed force-isotropy
cancellation was a real earlier finding **but STAGE-2 empirically retired it as the residual
cause** of the active floor.

### 0.1 VERIFICATION RESULT (2026-06-03) — the binding-fraction hypothesis is REFUTED

The Track-1 diagnosis proposed the floor was a **batch_dt rate-throttle artifact**: bound
fraction sits at ~3 % because binding is slow per tick, and would reach ~0.38 at equilibrium,
lifting g_soft. PI approved the fix **pending a verification run**. The run
(`scripts/track1_gsoft_verify.py`, batch_steps 100→**7600**, CFL product 9.9e-4 < 1e-3 ✓, to
**~3 τ**; `outputs/h3/track1_gsoft_verify_VERDICT.json`, `figs/fig_track1_gsoft_verify.png`)
**REFUTES it:**

- Bound fraction climbed **0.02 → 0.165 (8.2×)** over 0.1→3.0 τ, yet **g_soft stayed flat**
  (mean 7.3e-5 mN/m, **corr(bound_frac, g_soft) = −0.10**, second-half *lower* than first).
- Binding *did* lift g_soft **~3000× over the pre-binding baseline** (2.4e-8 → ~5e-5 mN/m) **but
  SATURATED at the very onset (bf = 0.02)** and never rose further — g_soft is **decoupled** from
  how many heads are bound. Heads are bound and loaded (F/F_stall ≈ 0.21, steady) the whole time.
- g_soft remains **~4800× under the band floor.**

**The real wall (now empirically isolated): FORCE AGGREGATION, not binding count/rate/timescale.**
A small fixed set of bound, loaded heads (≈0.21 F_stall each) yields a small fixed g_soft; adding
more bound heads does not increase it. The per-head force **does not aggregate into sustained
network tension** — exactly the second floor STAGE-2's head-tension note suspected, and the
predicted failure mode of the §1/KU35_FLOOR_ROOT_CAUSE "relabel-not-transport" mechanism: heads
load momentarily but the load is not sustained/summed into shell-wide stress (grip-walk re-stretch
too slow per tick, and/or no bipolar antiparallel organization → per-head forces don't sum into
net contraction). **The batch_steps fix is therefore NOT adopted.** (Secondary finding: binding
itself plateaus at bf ≈ 0.16, below the naive two-state 0.38 — a separate, smaller issue.)

**Next single-cell γ step (redirected):** target force *aggregation/sustenance*, not binding —
(a) confirm the grip-walk actually re-stretches the head spring against load each tick (s_grip per
tick vs ℓ₀; the KU35_FLOOR_ROOT_CAUSE §3.1 transport step), and (b) enforce bipolar sidedness so
+/− head sets grip antiparallel filaments (§3.3) → a net contractile dipole. Both are PI-gated core
myosin physics.

---

## 1. Magnitude anchors — the KU-3.5 band is physically CORRECT (Chugh 2017)

`emss-72183.pdf` — **Chugh, Clark, … Salbreux, Paluch, "Actin cortex architecture regulates
cell surface tension," Nat Cell Biol 19:689–697 (2017), DOI 10.1038/ncb3525.** The single
most load-bearing paper here. **Role: magnitude oracle + architecture constants — NOT runtime.**

| Quantity | Value | Cell type / method | Provenance | Maps to |
|---|---|---|---|---|
| Cortical tension, interphase HeLa | ~few hundred pN/µm | AFM constant-height (Fischer-Friedrich Eq.1) | Fig 1c | **band sanity**: 1 mN/m = 1000 pN/µm ⇒ KU band [0.35,0.65] = **350–650 pN/µm** |
| Cortical tension, mitotic HeLa | ~1000–1500 pN/µm | AFM | Fig 1c, 3d | ≈1.0–1.5 mN/m, ~2–4× over band ceiling (rounded mitotic) |
| **Model normalization T₀** | **230 pN/µm**; coherent peak T/T₀≈1.6 ⇒ **~0.37 mN/m** | 3D actin+myosin+xlink sim, anchored to myosin stall force | Fig 4 legend, 4c/d | **the magnitude oracle** — a correctly-built mechanistic cortex peaks at our band floor (0.35–0.37). The band is right; the deficit is real. |
| Filament-length optimum for max γ | **~400–500 nm** (non-monotonic; ≈0 below ~250 nm) | sim | Fig 4d | cortex segment-length / mesh architecture target |
| Cortex thickness | interphase ~300–450 nm, mitotic ~200 nm | dual-colour linescan | Fig 1c | thickness ↑ ⇏ tension ↑ — thickness is **not** the knob |
| Two-condition rule (text, p.6) | net tension needs **(1) dense+connected enough AND (2) tensile stress asymmetry** | sim + theory | p.6 | the architectural law (see §2) |
| Myosin NOT an architecture regulator | MYH9/10, blebbistatin, TPM4 → no thickness effect; tension drops on length perturbation with **myosin unchanged** | siRNA/drug screen | Fig 2a, 3e/f | tension scale set by **stall force × connectivity**, not myosin count alone |

**Verdict:** Chugh **validates the KU-3.5 band** (a coherent intermediate-length mechanistic
cortex peaks at ~0.37 mN/m = the band floor) and tells us tension scale = myosin **stall
force × connectivity**. Cross-read with STAGE-2 (connectivity z≈3 OK, coherence OK): by
Chugh's two-condition rule, with connectivity + asymmetry already satisfied, the remaining
factor is **myosin stall-force EXPRESSION** — i.e. force delivery/aggregation. This is
exactly STAGE-2's "force generation / soft-ceiling" conclusion, independently corroborated
by the literature. The head-tension diagnostic (44 % stall, 3 % bound, s_grip≈0) is the
concrete mechanism by which stall force fails to express.

---

## 2. What the papers DO and DON'T do for the active γ-floor

- **Chugh (§1) — DOES:** validate band + magnitude (T₀=230 pN/µm), supply the architecture
  optimum (filament length ~400–500 nm), and corroborate "stall-force-expression, not
  myosin-count" — pointing at the *force-delivery* wall STAGE-2 isolated. **DOESN'T:** hand
  us a runtime fix; raising expressed force is core contractile physics, **PI-gated** (and a
  k-stiffening sweep is hard-rule forbidden without a PI nod — flagged in `7e416c1`).
- **De Belly 2023 (`S0092867423005330`) + Keren preview (`S0092867423005858`) — CONFIRMATORY,
  not the wall.** Long-range cortical-tension transmission is **near-unattenuated** (inter-point
  delay **1.2 ± 1.2 s**, no resolvable decay length) **iff force engages a CONTINUOUS cortex**;
  force on membrane alone does not propagate. STAGE-2 already ruled transmission OUT as our
  residual cause, so these **confirm** transmission is not the floor. Their lasting value is the
  **membrane↔cortex coupling oracle** (3-tier model: elastic membrane *k* — MCA/ERM friction
  *µ* — viscous contractile cortex) for the future **H.9 membrane layer** ([`cortex/erm.py`]
  already exists), not the current wall. ⚠️ These measure **membrane** tension (a different
  quantity from actomyosin cortical tension); no absolute pN/µm is given — do not borrow a band
  number from them.
- **Nat Methods 2024 (`s41592-024-02277-8`) — membrane, not cortical.** O(0.1–0.5 mN/m)
  membrane-tension *changes* on HFF; confirms cytoskeleton **confines** membrane tension
  (bleb propagates, intact cell does not — confinement length ~0.7–1.7 µm). Membrane-tension
  perturbation-scale anchor only; **do not conflate with the cortical band.**

**Net:** the new literature **narrows and corroborates** the active-γ wall (force generation /
stall-force expression) but does **not** unblock it — that remains the PI-gated core-myosin
fidelity work STAGE-2 proposed (head-binding-fraction, soft-ceiling k-sweep, buckling).

---

## 3. The single-cell γ → spheroid surface-tension BRIDGE (runnable NOW)

The PI thesis is the **literature default**, not a leap: single-cell cortical tension enters
every tissue model as the per-cell surface/perimeter tension term, and aggregates to a tissue
surface tension via Young-Laplace. **Crucially, our `g_rigid` channel is ALREADY in-band
(0.57 mN/m)** — a real, usable cortical surface tension — so the bridge to the spheroid line
is **not blocked by the active-channel diagnosis.** This is the scientifically productive
track available immediately.

The three bridge papers (all paper-models ⇒ oracle/observable/bridge-relation, never runtime):

| Paper (file) | The bridge it supplies | Classification |
|---|---|---|
| **Fastabend et al. 2026, Phys Rev E 113:024403** (`z152-x4l1`) — "Cortical tension links curvature to tissue growth in the CPM" | (i) cortical tension = the CPM perimeter term `λ_P(P−P₀)²`; (ii) **2D Young-Laplace analog `R = λ/σ`** (radius set by line/surface tension balance) → 1/R; (iii) **growth rate ∝ local interfacial curvature** (∝1/R), with cortical-tension stretch ΔP as the proliferation signal, growth only on concave surfaces. **The published mechanism for our `A/A₀ = a + b/R + c/R²`.** | oracle + bridge-relation |
| **Okuda group 2026** (`2026.03.17.712503`) — 3D vertex, interfacial-tension modulation | (i) **Γ_interface = cortical contractility − adhesion** (verbatim) → our γ minus E-cadherin catch-bond cohesion = effective interfacial tension; (ii) **monolayer→3D cap transition once free-surface (cortical) tension > ~0.2 κ₀** → underwrites the L2.6 "MCF7 = 3D cap not monolayer" finding; (iii) virial tissue-stress coarse-grain `σ_αβ = Σ σ^j V_j / V_T` (observable); (iv) differential-Γ sorting (Steinberg DAH) for heterotypic spheroids | oracle + constant-anchor |
| **Roffay, Chan, Guirao, Hiiragi, Graner 2021, Development 148:dev192773** (`DownloadCombined…`) — inferring junction tension & pressure from geometry | (i) vertex force balance `Σt_i = 0` (measurement protocol); (ii) **Young-Laplace `ΔP = t·K`, 3D `K = 1/R + 1/R'`** — the exact pressure↔tension↔curvature law for a real 3D cap; (iii) **outer (cell-medium) tension ≈ 1.6–2× interior cell-cell tension** (pipette-validated) — the empirical signature that an aggregate **surface** is a higher-tension boundary = aggregate surface tension is cortex-derived | observable + bridge-relation + constant-anchor (⚠️ mouse embryo; ratio transfers, absolute a.u. do not) |

**The explicit chain (all literature-anchored):**
```
KU-3.5 single-cell cortical tension γ  (g_rigid in-band 0.57 mN/m)
   − E-cadherin catch-bond adhesion (Rakshit, L2.5)          [Okuda: Γ = cortical − adhesion]
   = effective interfacial tension Γ
   → coarse-grained tissue surface tension σ                 [Okuda virial; Roffay outer>inner ~1.6×]
   → Young-Laplace  ΔP = σ·(1/R + 1/R')                      [Roffay 3D]
   → curvature-proportional pressure & growth (∝ 1/R)        [Fastabend growth ∝ curvature]
   → the a + b/R + c/R² dependence of A/A₀ on spheroid radius [L2.5 G3 PASS, r²=0.98]
```
Our Layer-2 line **already reproduces the endpoint** (G3 PASS). These papers supply the
**published derivation** of why it works and the **bridge equations** to make the link
quantitative + the **observables** (Roffay vertex-tension inference; virial tissue stress) to
*measure* the emergent surface tension as a validation gate (L2 D2).

---

## 4. Measurement-protocol candidate — Nishitani-Miura V = σκ (`Quantitative_Analysis…`)

Nishitani & Miura 2025 (arXiv 2504.14887): estimate an **effective surface tension from a
shape time-series** via the curvature-velocity law **V = σκ** (+ volume term α·ΔA/A₀),
regressed on **retraction (non-protrusive) regions only**. MDCK, 2D. **Role: observable
protocol (relative σ), NOT a magnitude anchor.**
- **Mirrorable now**: we already emit GSD trajectories → binarize per frame → curvature κ →
  normal velocity V → regress on retraction regions → slope = effective σ. Same shape-observable
  family as the L2 hull/core-area metrics.
- **Caveats (hard):** the extracted σ has units **m²/s** (a curvature-flow mobility×tension
  coefficient), **not mN/m** — use as a **relative/comparative** observable, not an absolute band
  anchor; it is **2D** (needs a 3D-cap generalization or equatorial cross-section); it deliberately
  removes active-protrusion regions. Take absolute magnitude anchoring from Chugh (§1), not here.

---

## 5. Open items → PI ratification

- [ ] **Confirm the two-track split** (§2 vs §3): Track-1 = active-γ force-generation diagnosis
      (PI-gated core myosin; STAGE-2's head-binding-fraction / soft-ceiling k-sweep[needs nod] /
      buckling). Track-2 = the **bridge (runnable now)** using the in-band `g_rigid` cortical
      tension to drive the L2 D2 spheroid surface-tension validation. Recommend **Track-2 first**
      (unblocked, high-science, directly serves the PI thesis) while Track-1 stays PI-gated.
- [ ] Ratify **Chugh** as the KU-3.5 **magnitude oracle** (band-validation: model peak 0.37 mN/m
      = band floor; T₀=230 pN/µm) + filament-length architecture constant (~400–500 nm optimum).
- [ ] Ratify the **bridge equations** as L2 D2 oracles: Young-Laplace `ΔP = σ(1/R+1/R')` (Roffay),
      `R = λ/σ` + growth∝curvature (Fastabend), `Γ = cortical − adhesion` + 0.2κ₀ cap-threshold
      (Okuda); Roffay vertex-tension inference + Okuda virial as **L2 surface-tension observables**.
- [ ] Ratify **Nishitani-Miura V=σκ** as a *relative* shape-derived σ observable (not a band anchor).
- [ ] De Belly/Keren membrane↔cortex coupling (µ/k/ERM) → **defer to H.9 membrane layer**, not the
      current wall (transmission already ruled out by STAGE-2).
- [ ] Incorporate all ~11 PDFs into `references/tag_corpus.json` + `analysis/` extracts + Notion
      SourceEvidence (with the membrane-vs-cortical-tension distinction flagged) — after PI confirms
      direction (closeout step 2).

## 6. The papers (bibliographic, for tag_corpus incorporation)

1. Chugh et al. 2017, *Nat Cell Biol* 19:689–697, DOI 10.1038/ncb3525 — `emss-72183.pdf` ⭐ magnitude oracle
2. De Belly et al. 2023, *Cell* 186:3049–3061, DOI 10.1016/j.cell.2023.05.014 — `1-s2.0-S0092867423005330-main.pdf`
3. Keren 2023 (preview), *Cell* 186:2956–2958, DOI 10.1016/j.cell.2023.05.033 — `1-s2.0-S0092867423005858-main.pdf`
4. Lüchtefeld et al. 2024, *Nat Methods* 21:1063–1073, DOI 10.1038/s41592-024-02277-8 — `s41592-024-02277-8.pdf` (membrane, not cortical)
5. Nishitani & Miura 2025, arXiv:2504.14887 — `Quantitative_Analysis_of_Cell_Membrane_Tension_in_.pdf` (V=σκ observable)
6. Fastabend et al. 2026, *Phys Rev E* 113:024403, DOI 10.1103/z152-x4l1 — `z152-x4l1.pdf` ⭐ bridge (curvature→growth)
7. Thiticharoentam, Fukamachi, Horiguchi, Okuda 2026, bioRxiv 2026.03.17.712503 — `2026.03.17.712503v1.full.pdf` (Γ=cortical−adhesion, 3D cap)
8. Roffay, Chan, Guirao, Hiiragi, Graner 2021, *Development* 148:dev192773, DOI 10.1242/dev.192773 — `DownloadCombinedArticleAndSupplmentPdf.pdf` ⭐ bridge (Young-Laplace inference)
9. (peripheral) "3D Simulation of Tissue Mechanics with Cell Polarization" — `2023.03.28.534574v1.full.pdf` — tissue-sim, scan at incorporation
10. (peripheral) "Tuning cell motility via cell tension…" — `847046v1.full.pdf` / `main (1).pdf` (duplicate) — migration, peripheral to γ/cohesion
11. (peripheral) Piezo1 membrane-tension MD supplement — `mmc4.pdf` — membrane MD, supplement-level
