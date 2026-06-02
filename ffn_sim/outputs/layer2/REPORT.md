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
| L2.5 / L2.6 | cadherin catch-bond (KU-4.2) upgrade; 3D-Mikado invasion | ⬜ (L2.5 now the indicated next step — see L2.4.1 finding) |

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
