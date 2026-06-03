# Layer-2 multicellular spheroid line — REPORT (L2.0 + L2.1)

> Status 2026-06-02. Parallel CBM (center-based, 1 particle/cell) spheroid line, isolated
> from the single-cell main line, on the shared HOOMD + frozen-BAOAB stack. Brief:
> `docs/LAYER2_MULTICELL_DESIGN.md`. Anchor provenance: `docs/LAYER2_ANCHORS_2026-06-02.md`.

## Milestones

| Phase | What | Status |
|---|---|---|
| **L2.0** | measurement + acceptance-oracle + config layer (parameter-free) | ✅ DONE, 21 unit tests green |
| **L2.1** | CBM physics builder (Morse + reused BAOAB) + **G1 stable-aggregate gate** | ✅ DONE, G1 PASS (200 cells) |
| **L2.2** | edge-directed active-wetting traction (spreading driver) | ✅ DONE |
| **L2.3** | ensemble A/A₀(R₀) sweep + fit — **minimal CBM is cohesion-locked** (Δ A/A₀≈0.02) | ✅ DONE (model-limit finding) |
| **L2.4** | **contact-inhibited proliferation** (the size-dependent driver) + **G4 gate** | ✅ DONE, G4 PASS; A/A₀ signal now Δ≈12 but fit blocked by fragmentation (below) |
| **L2.4.1** | **connected-core spread-area estimator** (fragmentation-robust A/A₀) | ✅ DONE; de-noises the signal (fit r² 0.27→**0.74**, Δ 23→**1.0**) but G3≥0.95 **still FAILs** — confirms L2.5 is needed, not optional (below) |
| **L2.4b** | **leak-free pooled growth** (`run_growth_pooled`): ONE Simulation + pre-allocated particle pool, division activates a parked particle via `set_snapshot` | ✅ DONE; fixes the HOOMD per-rebuild memory leak (peak RSS **7 GB→240 MB**, flat), same physics. All growth drivers now use it. |
| **L2.5** | **E-cadherin catch-bond cohesion** (faithful Rakshit-2012 sliding-rebinding) replaces the static Morse well | ✅ DONE; **resists proliferation fragmentation (catch 0/5 vs morse 1/5 seeds, variance halved) → G3 r² 0.74→0.98 PASS** (below) |
| **L2.6** | **substrate confinement** (z=0 adhesive Morse wall, quasi-2D wetting) | ✅ DONE; cohesive MCF7 forms a 3D **cap** (not a monolayer) — the correct low-invasion phenotype; D_sub = the Bare/Pre/Lam4 ligand axis |

## ⭐ HEADLINE (L2.5) — the PI spreading law A/A₀ = a + b/R + c/R² EMERGES (G3 PASS, r²=0.98)

With the fully mechanistic model — **contact-inhibited proliferation** (the 1/R proliferating-rim
driver, L2.4) + **faithful E-cadherin catch-bond cohesion** (Rakshit 2012 sliding-rebinding,
L2.5) + the **connected-core spread observable** (L2.4a) — the experiment's novel law extracts
cleanly:

```
A/A0 = −0.33 + (188.7 µm)/R + (−2655 µm²)/R²       r² = 0.980   (5 R₀, 3 seeds each)
```

| R₀ (µm) | 31.7 | 40.4 | 53.1 | 67.1 | 78.3 |
|---|---|---|---|---|---|
| A/A₀ (core, mean±sd) | 2.95±0.17 | 2.81±0.05 | 2.22±0.11 | 1.82±0.04 | 1.72±0.05 |

**Why catch-bond unlocked G3.** The L2.4 static Morse cohesion let a *growing* spheroid
fragment (proliferation tension > fixed cohesion), inflating/scattering A/A₀ (morse r²=0.74,
1/5 seeds fragment). The Rakshit catch bond *strengthens under tension up to f₀≈29 pN* — exactly
the proliferation regime — so it holds the spheroid together (catch **0/5** fragment, variance
**halved** ±0.10→±0.05). hull≡core (no fragments). All signs match the PI law: b>0 (traction/
curvature, the dominant term, 0.67–0.76), c<0 (the documented "Bare" small-size cohesion
penalty). Gates G3 (r²≥0.95) **PASS**, G4 (rim 0.93, sub-exponential) **PASS**.

**A/B on the IDENTICAL pooled path (decisive):** morse cohesion core fit **r²=0.804 FAIL**
(one seed fragments → an A/A₀ outlier, sd 0.81) vs catch **r²=0.980 PASS**. Same proliferation,
same observable, same seeds — only the cohesion differs. The catch bond is *necessary* for G3:
the static Morse fragments under growth tension and scatters the law; the force-strengthening
catch holds it. (Figures `_morse.png` vs `_catch.png`.)

**Mechanistic chain (all measured/derived-anchored, no tuned constants):** MCF7 doubling 30 h
(BNID 100685) → proliferating rim ∝ 1/R; cohesion = N_cad≈223 cadherins/contact (Iturri 6.5 nN
de-adhesion / Rakshit f₀) each a sliding-rebinding catch bond (Rakshit 2012 SI Table S1) →
force-strengthening to f₀; spread measured as the connected-core footprint. The PI A/A₀ values
remain overlay-only (never fit). Figure: `fig_layer2_aa0_growth_law_catch.png`.

## Anchored / derived parameters (MCF7)

| Quantity | Value | Provenance |
|---|---|---|
| cell diameter (r₀) | 15.0 µm (R=7.5 µm) | Wagner 2011, MEASURED (Coulter), PMC3147247 |
| adhesion well depth D_e | 1.20e-17 J (≈2804 kT) | DERIVED `N_cad·⟨F⟩·Δx*` (KU-4.2); MCF7 cohesion is a documented literature absence — **PI ratification pending** |
| per-cell Stokes drag γ | 9.77e-8 N·s/m | DERIVED `6πηR` (KU-1.26) |
| Morse α / r_cut | 6.67e5 1/m / 22.5 µm | 1/contact-zone (modeling choice) / r₀+5·range (numerical policy) |
| CFL timestep dt | 9.16e-4 s | DERIVED `safety·γ/k_spring`, overdamped |

All derivations computed by `ffn_sim.spheroid.params.resolve_layer2` and unit-verified.
A/A₀ is overlay-only (never a fitting target). D_e magnitude does not affect G1.

## G1 stable-aggregate gate — PASS

A loose blob (200 cells, seeded at 1.1·r₀) settles under adhesion + excluded volume + the
overdamped BAOAB (30 000 steps, no motility). Gate bands in
`validation/oracles/configs/layer2_cbm.yaml`:

| Metric | Result | Band | Verdict |
|---|---|---|---|
| nearest-neighbour median / r₀ | 0.981 | [0.90, 1.20] | PASS |
| detached ("gas-like") fraction | 0.0000 | ≤ 0.02 | PASS |
| Rg growth factor (final/settled) | 1.000 | ≤ 1.50 | PASS |

Physics read: the blob relaxes from the seeded cubic lattice to a liquid-like cohesive
packing at ≈r₀ (NN slightly < r₀ from many-body inward pull of 2nd/3rd-neighbour tails — a
correct solid-packing signature), stays fully cohesive (zero stragglers), and is stable
(no dispersal). The frozen BAOAB integrator + a single `md.pair.Morse` reproduce a stable
multicellular aggregate — the Layer-2 line now simulates.

## L2.4 — contact-inhibited proliferation (the size-dependent driver) — G4 PASS

L2.3 established that the minimal CBM is **cohesion-locked** (Δ A/A₀ ≈ 0.02, no smooth law) —
a model-limit finding, *not* a tuning miss: the law needs an **added mechanism**. The
experiment runs over **days** and MCF7 spheroid spreading is partly **proliferation-driven**,
so L2.4 adds proliferation as a fine-grained mechanism (no fitted terms):

- **Cell-cycle timer** per cell, mean = MCF7 uncrowded doubling **30 h** (BNID 100685;
  MCF7 30–40 h standard culture, density-dependent — and that density dependence is itself
  the contact-inhibition arm we model), per-cell CV 0.15 (desynchronisation; flagged).
- **Division gate = Drasdo-Höhme FREE-SPACE rule** (`spheroid/proliferation.py`): a
  timer-elapsed cell divides only if it can bud a daughter ≥ one repulsive-core radius
  (`r₀ − contact_zone_width ≈ 0.9·r₀`, **derived** from the Morse shape, not tuned) from every
  other cell — i.e. there is room. A buried bulk cell has no free face → quiescent; a rim cell
  buds outward at the rest separation. *Diagnostic that drove the design:* the settled liquid-
  like packing has bulk first-shell coordination only ≈11 (not the FCC 12), so a neighbour-
  **count** threshold mis-classifies — but bulk free-gap ≈ 0.7·r₀ vs rim ≈ 1.3·r₀ separates
  cleanly. The free-space rule (not a count) is the operative discriminator.
- Because dividing cells form a **rim of ~constant thickness**, the proliferating fraction
  ∝ surface/volume ∝ **1/R** — the EMERGENT geometric origin of the law's 1/R, 1/R² terms.

**G4 gate (oracle config `layer2_cbm.yaml` g4) — PASS:**

| Metric | Result | Band | Verdict |
|---|---|---|---|
| rim-localised division fraction (large spheroid) | 0.86 | ≥ 0.70 | PASS |
| sub-exponential growth (all R₀ below 2^(t/τ)=4) | 2.06–2.86 | < 4.0 | PASS |
| dilute-limit doubling (isolated cell, exponential) | ✓ | 2^n | PASS (test) |
| proliferation-OFF reduces to G1 (count conserved) | ✓ | factor 1.0 | PASS (test) |

**The signal is now real** — proliferation converts the cohesion-locked Δ A/A₀ ≈ 0.02 into
Δ A/A₀ ≈ **12** with the **correct sign**: a small spheroid (R₀=30 µm) grows 2.86× while a
large one (R₀=77 µm) grows 2.14× — the surface/volume 1/R effect, mechanistically.

**But the clean a+b/R+c/R² law does NOT yet extract (G3 r²=0.27 — FAIL, reported honestly,
gate NOT loosened).** Diagnosis (verified, not speculated): the convex-hull A/A₀ is dominated
by a **proliferation-driven FRAGMENTATION instability**. A connected-component diagnostic
shows the spheroid splits into compact pieces that drift apart — `core-frac` < 1 and the
largest-component area stays modest (A/A₀_core 0.8–2.7) while the whole-set convex hull spans
the inter-fragment gaps and balloons to 5–16 (with ±19 seed variance). The static Morse well
**cannot hold a *growing* spheroid together**: proliferation pressure exceeds the finite
cohesion (this is the L2.3 "size-specific fragmentation instability," now driven by growth).
Per-epoch relaxation is ample (≈2000 s ≫ τ_relax 17 s), so this is real physics, not under-
relaxation.

**Indicated next step = L2.5 (cadherin catch-bond), not tuning.** A catch bond *strengthens
under tension* — exactly the tension proliferation generates. The static Morse cohesion is the
wrong model under active growth pressure; the mechanistic KU-4.2 / Rakshit-2012 catch-bond is
the physically-correct fix for fragmentation. Secondary refinements (also non-tuning): z=0
**substrate confinement** (quasi-2D wetting keeps a monolayer, the experiment's actual
geometry) and an **outlier-robust / connected-core area** estimator (the alpha-shape refinement
already flagged in `observables.projected_area`). All preserve the literature-first discipline.

## L2.4.1 — connected-core spread-area estimator (fragmentation-robust A/A₀)

The first of the secondary refinements above is now built and measured. `observables.py` gains
`connected_components` (single-linkage KD-tree + union-find, `link_radius` a caller argument set
to 1.6·r₀ — between the 1st and 2nd coordination shell, **no baked constant**),
`largest_connected_component`, and `core_projected_area` (convex hull of the **largest connected
component** only). `proliferation.run_growth` now reports `area_core_over_a0` alongside the raw
hull; the growth sweep fits the PI law to the **core** area as the G3 headline (the hull is
fragmentation-inflated and unphysical — its fit even dips below A/A₀=1). 3 new observable tests
(two-cluster labelling, drifting-fragment rejection, single-cluster identity) → **39 green**.

**Result (5 R₀ × 3 seeds, memory-safe per-size driver — see Verification):**

| R₀ (µm) | A/A₀ core (mean±sd) | A/A₀ hull (mean±sd) |
|---|---|---|
| 29.1 | 2.53 ± 0.16 | 2.73 ± 0.44 |
| 38.0 | 2.22 ± 0.10 | 2.45 ± 0.22 |
| 52.2 | 1.78 ± 0.11 | 3.34 ± 2.11 |
| 63.7 | 2.12 ± 0.62 | 25.48 ± 21.55 |
| 75.2 | 1.49 ± 0.03 | 20.17 ± 20.24 |

- **fit CORE:** A/A₀ = 0.953 + (58.99 µm)/R + (−394.2 µm²)/R²  **r²=0.737**
- **fit HULL:** A/A₀ = 80.67 + (−5551.8 µm)/R + (95714 µm²)/R²  r²=0.750 (large `a`, sub-1 dip — unphysical)

**Honest verdict — the core estimator helps but does NOT rescue G3.** It collapses the
fragmentation inflation (signal Δ 23→**1.0**, error bars from ±21 down to ±0.03–0.16) and recovers
a near-monotone 1/R decrease (the correct sign) — fit r² rises **0.27→0.74**. But r²=0.74 is still
below the **G3 ≥0.95 band → G3 remains FAIL (gate NOT loosened)**. The residual scatter is real
physics, not estimator noise: at R₀=63.7 µm one of three seeds fragments so hard that even the
*core* inflates (core 2.12±0.62, the visible outlier), so a static Morse well cannot hold a
growing spheroid together even when measured robustly. **This confirms L2.5 (catch-bond) is
required, not optional** — the robust estimator was necessary to *see* the residual instability
cleanly, but the cohesion model itself is the remaining blocker.

## L2.6 — ligand-condition axis (Bare / Pre / Lam4) — mild passive effect (honest)

Sweeping the cell-substrate adhesion (the coarse ligand knob; Bare/Pre/Lam4 = ×0.5/×1.0/×2.0
the cohesion D_e) over R₀ gives three emergent A/A₀(R₀) curves that **nearly overlap** (e.g.
R₀≈40 µm: Bare 2.57, Pre 2.55, Lam4 2.78 — Lam4 highest, the correct direction, but small;
they converge at larger R₀). **Finding:** passive substrate adhesion alone does *not* reproduce
a strong ligand-condition separation — cohesive MCF7 forms a 3D cap regardless (L2.6). The
experiment's Bare/Pre/Lam4 differences therefore likely require the **active** ligand mechanism
(per-species integrin catch-slip kinetics in `bridge/ligand_species.py` driving edge traction),
the documented faithful upgrade — not the passive adhesion depth. (The per-condition r²=1.000 is
degenerate: 3 R₀ vs 3 coefficients; a trend comparison only.) Figure
`fig_layer2_ligand_conditions.png`. PI A/A₀ overlay-only.

**Active traction = the actual ligand driver (validated direction).** Adding edge-directed
active-wetting traction (L2.2 `spreading.edge_outward_forces`, now composed into the pooled
catch+substrate growth via `run_growth_pooled(f_traction=…)`) increases the spread monotonically
(core A/A₀ 1.53→1.56→1.81 at 0/3/6 nN) — the correct mechanism, where passive substrate adhesion
gave almost none. So the faithful Bare/Pre/Lam4 separation is **ligand-modulated active traction**
(per-species integrin catch-slip, `bridge/ligand_species.py`), the next anchoring step — not
substrate-adhesion depth.

**Traction axis (mechanism → law coefficient).** Sweeping active traction (0 vs 3 nN) over R₀
raises A/A₀ at *every* size (0 nN: 2.55/2.25/1.82 → 3 nN: 2.80/2.35/1.92), i.e. the law's
**b-coefficient (the 1/R traction term) grows with traction** — the platform doing its stated
job: mapping the a/b/c terms to mechanism. Figure `fig_layer2_traction_axis.png`. ⚠️ The 6 nN
level hit a numerical box/substrate-wall limit (a cell ejected past the box under
traction+growth+wall) — a known fix (size the pool box for the grown+spread footprint + soften
the Morse wall repulsive branch), not a physics error.

## B1 — box-sizing + ejection guard (numerics robustness, 2026-06-03)

The L2.6 traction-axis run crashed at high traction with a HOOMD C++ `RuntimeError: Particle …
is out of bounds`. Diagnosed (reproduced at f_traction=6–10 nN, Lam4 substrate, N₀=450, 2
doublings) as **two compounding causes**:

1. **Box too small for the wetting footprint (the real bug).** The pre-allocated pool box
   (L2.4b) was sized ONCE from `r_cluster_max = r0·N_max^(1/3)·1.3` — a *free-3D-ball* radius.
   But a substrate-confined spheroid (L2.6) **wets into a quasi-2D disk** whose in-plane radius
   scales as `√N_max`, far larger than `N_max^(1/3)`. The cluster spread past the box edge and
   a cell wrapped across the periodic boundary → out-of-bounds. **Fix:** `cbm.pool_cluster_radius`
   sizes the box for `max(3D-pack, 2D-wetting disk) × spread-safety` (derived close-packing
   geometry + a containment margin — numerical policy, never enters a force or a measurement),
   used by both pool builders (`build_pool_simulation`, `build_cbm_catch`).
2. **Genuine edge-cell detachment beyond regime (a model limit, now surfaced not hidden).**
   Once box-sizing was fixed, f_traction ≳ cohesion (~6.5 nN) still ejects a *boundary* cell:
   in the overdamped large-dt CBM a cell whose net outward traction exceeds its local cohesion
   detaches and, at terminal velocity F/γ over the ~2000 s epoch dt, leaves instantly (the
   adiabatic CBM cannot represent a slowly-peeling cell). The earlier 6 nN "success" was a
   **PBC self-interaction artifact** of the too-small box (the wrapped image artificially
   re-confined the cluster). **Fix:** a graceful **ejection guard** in `run_growth_pooled`
   (`try/except` the HOOMD out-of-bounds + a post-epoch finiteness/0.45·L containment check)
   stops cleanly on the last good state and returns `ejected=True` instead of crashing the
   whole sweep; the sweep/ligand scripts surface `⚠EJECTED k/n_seeds` (no silent truncation).

**Verified:** in-regime production (f_traction ≤ 3 nN, the REPORT's validated levels) runs clean
(`ejected=False`, traction direction preserved: 0→3 nN gives A/A₀ 1.88→1.95 at Lam4); 6–10 nN
no longer crashes (stops cleanly, flagged); the **G3 headline path is unchanged** (no substrate/
traction: N₀=120→A/A₀ 2.63, N₀=400→1.86, consistent with the §HEADLINE table). +4 tests
(`pool_cluster_radius` 2D-dominance/monotonicity/guards, box-contains-worst-case-spread,
`ejected=False` on a bounded run) → layer-2 suite **88 green** (1 unrelated cupy skip).

**Operating-regime conclusion (for A1).** Active edge-traction must stay **below the
detachment threshold** (≤ ~3 nN at the measured 6.5 nN cohesion); the ligand→traction anchor
(A1) should map the Bare/Pre/Lam4 conditions into that stable band. Higher traction is a genuine
detachment regime the overdamped CBM cannot resolve (a GPU sub-stepped-bond option, roadmap D3),
not a numerical bug to patch.

## A1 — ligand → ACTIVE edge-traction (the faithful Bare/Pre/Lam4 driver, 2026-06-03)

L2.6 mapped the ligand condition onto the *passive* substrate-adhesion depth and the three
A/A₀(R₀) curves barely separated (Δ≈0.03 — a cohesive MCF7 caps regardless). A1 is the faithful
upgrade: each PI condition maps onto the cell's **active edge-traction** through the per-species
integrin-clutch kinetics (`spheroid/ligand_traction.py`, anchored to `bridge/ligand_species.py`,
the literature SoT), the validated spreading knob. The passive wall is held COMMON so the *only*
thing differing is the ligand-set active traction — isolating the A1 mechanism.

**Mechanism (anchored, no fit to PI):** `f_traction = T_ref · density · clutch_strength`.
`clutch_strength(ligand) = φ·F_s` relative to col-I, with `φ = k_on/(k_on+k_off0)` (k_on =
KU-2.4 0.3 s⁻¹) and `F_s = k_BT/x_β` — both from the registry kinetics. With col-I (k_off0
1.3 s⁻¹, x_β 0.23 nm) and laminin-111 (k_off0 1.85 s⁻¹, x_β 0.28 nm, α7β1-invasin slip proxy):

| Condition | ligand | φ | F_s | clutch strength | density | **f_traction** |
|---|---|---|---|---|---|---|
| Bare | col-I (α2β1) | 0.187 | 18.6 pN | 1.000 | 0.60 | **1.50 nN** |
| Pre | col-I (α2β1) | 0.187 | 18.6 pN | 1.000 | 1.00 | **2.50 nN** |
| Lam4 | laminin-111 (α6β1) | 0.140 | 15.3 pN | **0.611** | 1.00 | **1.53 nN** (proxy) |

The col-I-vs-laminin ordering (laminin = 0.61× the col-I clutch — weaker occupancy AND smaller
per-bond force) is **the measured direction** (breast-epithelial traction lower on LN-111,
P=0.016–0.028), not guessed. All three sit in the B1 stable ≤3 nN band.

**Emergent result (5 R₀×… → here 4 R₀ × 2 seeds, catch cohesion, common substrate):**

| R₀ (µm) | 40.4 | 53.5 | 67.3 | 78.3 |
|---|---|---|---|---|
| Bare (1.50 nN) | 2.71 | 2.18 | 1.89 | 1.79 |
| Pre  (2.50 nN) | 2.90 | 2.38 | 1.95 | 1.77 |
| Lam4 (1.53 nN) | 2.68 | 2.26 | 1.92 | 1.75 |

**The active mechanism SEPARATES the conditions** — Δ(A/A₀) ≈ **0.139** at R₀≈60 µm (and ≈0.22
at R₀≈40 µm, above the ±0.05–0.09 seed sd), **~4.6× the L2.6 passive-adhesion separation
(~0.03)**. The ordering **Pre > Lam4 ≳ Bare tracks the resolved traction monotonically**
(2.50 > 1.53 > 1.50 nN) — the platform doing its job: ligand identity → clutch kinetics →
traction → spread, end-to-end mechanistic. Figure `fig_layer2_ligand_traction_conditions.png`.

**Honest caveats (reported, not smoothed over):**
- **Per-condition a/b/c are UNDER-DETERMINED.** 4 R₀ vs 3 coefficients = 1 dof; the fits are
  near-perfect r²≈1.00 by interpolation, and the *individual* a/b/c swing wildly (Bare a=1.41,
  b=5.3 µm, c=+1911 µm² with an unphysical c>0; Pre a=0.10, b=147 µm, c=−1369 µm²) — an artifact
  of the degeneracy (same caveat the REPORT flagged for L2.6's 3-point fit), NOT physics. **The
  A/A₀(R₀) CURVES and their SEPARATION are the robust A1 result; robust per-condition a/b/c with
  error bars need A3** (more R₀ / seeds / biological time, on GPU = B2).
- **Absolute traction scale is unanchored** (no MCF7 single-cell traction in the literature; the
  15–25 nN micropillar value was REFUTED). `T_ref` is set in the B1 stable band — the **relative**
  ordering is the anchored science, not the magnitude.
- **Density axis (Bare<Pre) is a flagged modeling knob** pending the collaborator's pV4D4 col-I
  adsorption-density data (PI-exp map: literature gap, route to Im Sung Gap / KAIST). The
  ligand-IDENTITY axis (col-I vs laminin) is the fully-anchored part.
- **Single-cell ↔ collective laminin split (expected).** Lam4's *single-cell* traction (1.53 nN,
  the weak laminin clutch) sits just above Bare, so via single-cell active traction Lam4 does
  **not** dramatically out-spread (it lands between Bare and Pre). The PI poster's *collective*
  Lam4 enhancement, if present, is the documented single-cell↔collective split (PI-exp map) —
  the A2 overlay question. A1 reports what the measured single-cell clutch kinetics produce and
  deliberately does NOT engineer the collective ordering (overlay-only hard rule).

8 new tests (`test_ligand_traction.py`: anchored col-I/laminin ordering, ≤3 nN band, Bare/Pre
density-only difference, sign-sense, input guards) → layer-2 suite **96 green**.

## A2 — PI poster overlay (OVERLAY-ONLY, never fit; 2026-06-03)

The PI provided the poster A/A₀ exports (`references/260313_{Bare,Pre,Lam4}.csv`, per-spheroid
time series of segmented area / effective radius; kept LOCAL/gitignored per the hard rule —
only this overlay figure + the extracted summary are committed). Each spheroid (`Series`)
gives one law point: R₀ = its effective radius at t=0, A/A₀ = its spread ratio at a common
observation time (sampled at t≈60 h to match the platform's 2-doubling biological time; 82 h
secondary). `scripts/layer2_pi_overlay.py` overlays these on the A1 emergent curves and reports
the agreements/gaps — **no model parameter is tuned to the PI data** (overlay-only).

**PI dataset (extracted):** Bare 8 spheroids R₀ 140–419 µm; Pre 25, R₀ 103–386 µm; Lam4 26,
R₀ 87–255 µm; all run ~82 h.

**⭐ The headline — the novel law's SHAPE is mechanistically reproduced.** In all three
conditions the PI A/A₀ **decreases with R₀** (corr(R₀,A/A₀) = −0.86 / −0.81 / −0.91; fitted
b > 0 dominant) — i.e. the experiment's own data carries the **1/R size dependence** that the
platform's fully-mechanistic model produces *emergently* (proliferating-rim surface/volume,
L2.4 → L2.5). The PI's genuinely novel A/A₀ = a + b/R + c/R² law (no published analog) and the
mechanistic CBM **agree on the fundamental sign/shape** — the platform reproduces *why* small
spheroids spread relatively more. This is the qualitative validation A2 set out to test.

**The three honest gaps (each points to a concrete next step):**

| Axis | PI (experiment) | Platform (A1) | Verdict |
|---|---|---|---|
| law shape / sign | A/A₀↓ with R₀ (corr ≈ −0.85, b>0) | A/A₀↓ with R₀ (b>0) | **MATCH ✓** |
| condition ordering | **Lam4 > Pre > Bare** (collective) | **Pre > Lam4 ≳ Bare** (single-cell) | **SPLIT** |
| magnitude (med A/A₀) | ≈ 7.5 (60 h) / 10 (82 h) | ≈ 2.1 | platform under-spreads **~4–5×** |
| R₀ range | 87–419 µm | 40–78 µm | **no overlap** → native-N (GPU) |

1. **Ordering = the single-cell↔collective laminin split, now CONFIRMED with data.** The
   experiment's *collective* ranking puts **Lam4 highest**; the platform's *single-cell*
   active-traction (A1) puts Lam4 ≈ Bare (laminin is the weaker single-cell clutch, 0.61×).
   This is exactly the split the PI-exp map flagged. **Implication:** the Lam4 collective
   enhancement is NOT single-cell traction — it must be a *collective* mechanism (the "uniform
   β1" → more uniform proliferation / a cohesion-modulation that lifts the small-size c-penalty,
   the documented Lam4 "c→0" phenotype). That is the next mechanistic hypothesis to test — and
   the platform predicting Lam4 ≠ single-cell-traction-driven is itself a useful, falsifiable
   result, not a failure.
2. **Magnitude ~4–5×.** The platform's connected-CORE area (deliberately conservative, L2.4a)
   + cohesion-locked catch bond vs the experiment's RAW segmented area (which includes spread
   protrusions / scattering the cohesive model resists) + the platform's 60 h vs 82 h. Expected;
   the core/raw and time axes are recoverable (a raw-area readout + longer biological time).
3. **R₀ range — no overlap.** The PI spheroids (R₀ 87–419 µm ≈ 10³–10⁴ cells) dwarf the
   platform's CPU first-pass (R₀ 40–78 µm). Matching the experiment's sizes needs native-N on
   GPU (B2) — the platform fit is shown EXTRAPOLATED (dotted) into the PI range and flagged.

**Net:** the platform reproduces the experiment's novel-law SHAPE (the science win); the
ordering, magnitude, and size-range gaps are characterized honestly and each maps to a defined
next step (collective-Lam4 mechanism; raw-area + longer time; native-N GPU). PI A/A₀ stays
overlay-only throughout. Figure `fig_layer2_pi_overlay.png`; summary
`outputs/layer2/pi_overlay_summary.json`.

## Figures

Regenerate all via `python -m ffn_sim.scripts.layer2_vis` (the one-entry-point convention);
each driver also auto-generates its own figure at run end (production-driver-auto-viz rule).

- `figs/fig_layer2_g1_stable_aggregate.png` — G1 result. **Left**: initial loose blob
  (1.1·r₀ jittered cubic lattice). **Middle**: settled aggregate (lattice → disordered
  cohesive packing, slightly compacted). **Right**: nearest-neighbour-distance histogram
  with the r₀=15 µm rest separation overlaid (median/r₀ = 0.981). No axis truncation; SI
  (µm) units; reference line shown.
- `figs/fig_layer2_l2_2_motility_mechanism.png` — L2.2 motility mechanism. **Top**: settled
  aggregate (f_active=0) vs under active traction (6 nN). **Bottom**: A/A₀ and detached
  fraction vs active traction, with the measured cohesion/detachment force (6.5 nN, Iturri
  2020) overlaid. Honest: isotropic self-propulsion does NOT spread a cohesive cluster — the
  spreading driver is edge-directed traction (active wetting), built in L2.3.
- `figs/fig_layer2_aa0_law.png` — **L2.3 emergent A/A₀(R₀), ensemble-averaged (8 seeds/R₀).**
  Edge-directed active-wetting traction (5 nN, Lp=11 µm) vs measured cohesion (6.5 nN), swept
  over R₀=29–76 µm. **Honest finding:** at the measured MCF7 scales the minimal CBM
  (cohesion + edge-traction) is **cohesion-locked** — A/A₀ ≈ 1.0 ± 0.02 at most sizes (MCF7 is
  low-invasion, barely spreads), with a size-specific fragmentation instability near R₀≈40 µm
  (high variance). The smooth PI law a + b/R + c/R² does **NOT** cleanly emerge from the
  minimal model (fit r²≈0.5). Regime sweeps (traction 5–6 nN, 3–8 seeds) do not change this —
  it is a model-limit finding, not a tuning miss: reproducing the experiment's spreading law
  needs ADDITIONAL mechanism (proliferation / longer biological timescale / the ligand
  conditions the experiment varies), the next research direction. The full L2.3 machinery
  (ensemble R₀-sweep → fit → error-bar figure) is in place for that. PI A/A₀ is overlay-only.
- `figs/fig_layer2_l2_4_proliferation.png` — **L2.4 mechanism.** **Left**: a grown spheroid
  mid-slice (N: 250→527), each cell coloured by first-shell coordination — the low-coordination
  **rim** (proliferation-competent, free space) vs the high-coordination, contact-inhibited
  **bulk**. **Right**: N(t) and A/A₀(t) over 2 doublings of biological time; rim-localised
  fraction 0.86, growth 2.06 (sub-exponential). SI units, no truncation.
- `figs/fig_layer2_aa0_growth_law.png` — **L2.4 proliferation-driven A/A₀(R₀), ensemble
  (3 seeds/R₀, 5 sizes R₀=30–77 µm).** **Left**: A/A₀ vs R₀ with the a+b/R+c/R² fit and the
  A/A₀=1 (no-spread) reference. **Honest finding:** the signal is now strong and measurable
  (Δ A/A₀≈12, correct 1/R sign) — proliferation is the right driver — but the fit is poor
  (r²=0.27) because the convex-hull A/A₀ is inflated by a **proliferation-driven fragmentation
  instability** (huge ±sd error bars). **Right**: the mechanism — growth factor falls with R₀
  (∝ surface/volume ∝ 1/R) and rim-localised fraction (≥0.70 G4) vs the 2^(t/τ) exponential
  ceiling. The clean law extraction needs L2.5 catch-bond cohesion + substrate confinement
  (REPORT §L2.4). PI A/A₀ overlay-only.
- `figs/fig_layer2_aa0_core_vs_hull.png` — **L2.4.1 connected-core vs raw-hull A/A₀(R₀)**
  (5 R₀ × 3 seeds, per-size driver). **Left (CORE, G3 headline)**: fragmentation-robust A/A₀
  with per-realisation points + ensemble mean±sd + a+b/R+c/R² fit (r²=0.74) and the A/A₀=1
  reference — tight error bars, near-monotone 1/R decrease, with the R₀=63.7 µm outlier (one
  seed fragments even the core) visible. **Right (HULL)**: the raw convex hull for contrast —
  ±20 error bars and a fit that dips below A/A₀=1 (unphysical), inflated by drifting fragments.
  Same axes, no truncation, SI units. The figure is the visual proof that the core estimator
  de-noises the signal but the residual scatter (→ L2.5) is real.
- `figs/fig_layer2_l2_5_cadherin_catch_bond.png` — **L2.5 catch-bond oracle.** **Left**: the
  faithful Rakshit-2012 sliding-rebinding lifetime τ(f) (catch peak F*≈28.5 pN ≈ f₀=29.2 pN,
  then slip) vs a pure Bell slip; **right**: the new-interaction probability Pₙ(f) ramp and the
  effective k_off(f) (dips at the catch peak, rises in the slip regime). SI units.
- `figs/fig_layer2_l2_5_fragmentation_resistance.png` — **L2.5 catch resists fragmentation.**
  **Left**: per-seed core A/A₀ (N₀=400, 5 seeds) for morse vs catch — morse 1/5 seeds fragment
  (annotated) with ±0.10 scatter, catch 0/5 with ±0.05. **Right**: why — the effective cohesion
  force law F_coh(ext) strengthens to a peak at per-cadherin f₀≈29 pN (overlaid: measured 6.5 nN
  de-adhesion) then slip-ruptures. The force-strengthening is the fragmentation fix.
- `figs/fig_layer2_pi_overlay.png` — **A2 PI poster overlay (overlay-only).** **Left**: PI
  per-spheroid points (○) + their a+b/R+c/R² fit (solid) for Bare/Pre/Lam4, with the platform's
  A1 emergent points (◇) + fit (dashed) and its extrapolation into the PI R₀ range (dotted,
  flagged) — both families DECREASE with R₀ (the shared 1/R law); the platform sits ~4–5× lower
  and at smaller R₀. **Right**: median A/A₀ per condition, PI (solid) vs platform (hatched) —
  the ordering split (PI Lam4>Pre>Bare collective vs model Pre>Lam4≳Bare single-cell) and the
  magnitude gap. A/A₀=1 reference shown, SI units, no truncation.
- `figs/fig_layer2_ligand_traction_conditions.png` — **A1 ligand→active-traction.** **Left**:
  the three emergent A/A₀(R₀) curves (Bare/Pre/Lam4) with per-realisation points + a+b/R+c/R²
  fit and the A/A₀=1 reference — Pre (highest traction) above, Lam4≈Bare, separation Δ≈0.14 at
  mid-R₀ (4.6× the L2.6 passive). **Right**: the mechanism — resolved `f_traction` per condition
  (= T_ref·density·φ·F_s) with the col-I/laminin clutch strength (0.61×) + density annotated and
  the B1 stable 3 nN ceiling overlaid. SI units, no truncation.
- `figs/fig_layer2_aa0_growth_law_catch.png` — **⭐ L2.5 G3-PASS A/A₀(R₀) law.** The catch-bond
  ensemble (5 R₀ × 3 seeds): A/A₀ vs R₀ with tight error bars + the a+b/R+c/R² fit (**r²=0.980**)
  + the A/A₀=1 reference; right panel shows growth-factor∝1/R + rim fraction (G4). The clean
  emergence of the PI's novel law from the fully mechanistic model. PI A/A₀ overlay-only.

## Verification

- `tests/test_spheroid_observables.py` (16) + `tests/test_layer2_params.py` (8) +
  `tests/test_spheroid_proliferation.py` (15) = **39 green** (synthetic clouds vs closed-form
  oracle; resolve derivations vs Magic-Number-Block; free-space gate / rim-localisation /
  contact-inhibition / dilute-doubling / G1-reduction; **+ connected-component labelling /
  drifting-fragment rejection / single-cluster identity** for the L2.4.1 core estimator).
- `scripts/layer2_g1_smoke.py` — reproducible G1 run; `scripts/layer2_growth_smoke.py` —
  reproducible single-spheroid growth + G4 verdict; `scripts/layer2_aa0_growth_sweep.py` —
  proliferation-driven A/A₀(R₀) sweep + fit + figure (`layer2_aa0_sweep.py` kept as the
  cohesion-locked baseline). Re-running the baseline confirms it is unchanged (still
  cohesion-locked); the growth sweep is the L2.4 headline.
- **L2.4.1 reproducible result:** `scripts/layer2_aa0_growth_persize.py` — memory-safe per-size
  driver (one `(N0,seed)` per process; the all-in-one sweep OOM/SIGKILLs around N0≈250–400 on a
  16 GB CPU box because the HOOMD epoch-rebuild loop accumulates resident memory across sizes).
  Canonical numbers + the core-vs-hull fit are in `outputs/layer2/growth_sweep_core_results.jsonl`
  + `growth_sweep_core.log` (run `… --fit growth_sweep_core_results.jsonl`).
- Isolation: runtime imports NO oracle (hard rule); **additive new files only**
  (`spheroid/proliferation.py`, growth scripts, `test_spheroid_proliferation.py`);
  `build_cbm_simulation` gained a backward-compatible optional `positions=` arg; single-cell
  main line (`cell/ cortex/ bridge/`) + `integrator/` freeze untouched.

## Open (PI ratification)

1. D_e: measured Iturri-2020 nN force-anchor adopted (was the retired pN seed) — FYI only.
2. Surface-tension validation target: emergent-only vs non-MCF7 proxy (MCF10DCIS ~21 mN/m).
3. Cell-size band position: 15 µm (low end) vs 17–18 µm.
4. **L2.4 / L2.4.1 result for PI direction:** proliferation is the correct size-dependent
   driver (G4 PASS, strong signal). The robust connected-core spread-area estimator (L2.4.1) is
   now built and de-noises the signal (fit r² 0.27→0.74, Δ 23→1.0) — but **G3 still FAILs
   (r²=0.74 < 0.95)** because a static Morse well cannot hold a *growing* spheroid together even
   when measured robustly (one R₀=63.7 µm seed fragments the core). So **L2.5 catch-bond is
   confirmed necessary, not optional**; substrate confinement remains the secondary refinement.
   Confirm this direction (and the L2.5 entry point).

## Next (L2.5)

Replace the static Morse well with the mechanistic KU-4.2 cadherin **catch-bond** (Rakshit
2012 PNAS, force-strengthening) so cohesion resists the proliferation-driven tension that
currently fragments the growing spheroid; add z=0 **substrate confinement** (`ecm/substrate.py`,
quasi-2D wetting = the experiment geometry). The outlier-robust connected-core spread area is
**done (L2.4.1)** — it is the area observable L2.5 will be scored on; then re-run the growth
sweep → fit A/A₀ = a + b/R + c/R² (G3) and overlay the PI poster.
