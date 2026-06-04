# H.6 (provisional) — Filopodium / invadopodium — FUTURE / DEFERRED unit

> **Status**: 📋 **REGISTERED FUTURE SCOPE — not implemented, not on the current
> critical path.** Registered 2026-06-04 on PI request ("등록해놔줘") so the omission is an
> explicit, evidence-grounded decision rather than a silent gap. Unit number **H.6 is
> provisional** — PI assigns the final number when the unit is activated.
>
> **Owner**: TBD (single-cell line). **Prereq**: H.5 lamellipodium (shares the barbed-end
> elongation + membrane-load machinery) + the invasion arm (L2.7 / 3D-Mikado ECM).

## Why filopodia is currently 0 (and why that is defensible)

Lamellipodium (H.5) was scoped and implemented; filopodium was never on the H.1→H.2→H.3→H.5→H.7
chain. This is a deliberate prioritization, grounded in the literature (KB query 2026-06-04):

1. **Role split (Mattila & Lappalainen 2008, Nat Rev MCB, `MattilaLappalainen2008_NRMCB`,
   the canonical filopodia review; corroborated by the new membrane-tension paper
   `quantitative-analysis-of-cell-membrane-t-7a514c`):** *"Lamellipodia are thin, wide
   structures with strong adhesion to the substrate and have been widely studied as the
   **primary protrusive structures involved in two-dimensional substrate-based motility**.
   Filopodia are long, fiber-like structures that work during cell migration, **microtubule
   guidance, and signal transduction**."* → Lamellipodia drive the **bulk area expansion**
   (the PI experiment's A/A₀ spreading, the b/R traction term); filopodia are **sensing /
   guidance / pathfinding antennae**, fascin-bundled parallel spikes.
2. **MCF7 phenotype:** cohesive luminal **epithelial** (low-invasion) — spreads as a
   sheet/monolayer-cap (L2.6), not via filopodial probing. Filopodia are prominent in
   **invasive/mesenchymal** lines (MDA-MB-231) and growth cones, not cohesive MCF7.
3. **Area contribution is negligible:** filopodia are ~100-300 nm-wide spikes; they add
   essentially **zero spread AREA** (the connected-core A/A₀ observable). On spreading
   substrates filopodia are even **converted into lamellipodia-like protrusions**
   (Mattila 2008, p1, spreading-fibroblast micropattern study).
4. **Not a cortical-tension source:** like lamellipodia, filopodia are **protrusive, not
   contractile-shell** structures — they do not generate cortical tension (KU-3.5). Their
   absence is therefore **not a confound** for the cortical-tension floor diagnosis.
5. **Shared machinery → cheap later:** filopodia reuse the H.5 **barbed-end elongation +
   membrane-load (Bieling 2016 force-velocity)** machinery; they differ only by **fascin
   parallel-bundling** (vs Arp2/3 branching). So H.6 is a *bundling extension* on top of
   the H.5 base, deferrable at low marginal cost.

## When H.6 activates (the trigger)

The **invasion arm**: roadmap **L2.7 invasion** (3D-Mikado ECM, invasive phenotype) and any
move to **mesenchymal / invasive lines** (MDA-MB-231, glioma). There, **filopodia /
invadopodia** (ECM probing, durotaxis/topo-sensing, matrix-degradation podia) become
first-class mechanisms. Until then H.6 is correctly deferred.

## Provisional scope (when built — mechanistic, per CLAUDE.md inversion rule)

- **Fascin-bundled parallel actin** finger (vs H.5 Arp2/3 branched mesh) — explicit fascin
  crosslink bonds, NOT a lumped bundle-stiffness scalar.
- **Mogilner-Rubinstein buckling-limited protrusion force** (a 10-30-filament bundle buckles
  at ~tens of pN; force is an emergent buckling output, not a prescribed value) — acceptance
  oracle only.
- **Tip complex** (formin/Ena-VASP barbed-end processivity) + convergent-elongation seeding
  from a lamellipodial network (Svitkina model).
- **Invadopodium variant**: + MMP matrix-degradation coupling against the 3D-Mikado ECM.

## Acceptance gates (provisional, to be ratified with the unit)

- Bundle force vs filament count (Mogilner-Rubinstein buckling); filopodial length/lifetime
  distribution; (invadopodium) ECM-degradation rate vs MMP density. All emergent, oracle-checked.

## Cross-references

- H.5 lamellipodium (`docs/briefs/H5_lamellipodium.md`) — the base protrusion machinery H.6 extends.
- L2.7 invasion (`docs/LAYER2_ROADMAP_2026-06-03.md` §D1) — the activation trigger.
- Quantitative-negligibility basis: KB query 2026-06-04 (`MattilaLappalainen2008_NRMCB`,
  `quantitative-analysis-of-cell-membrane-t-7a514c`, `Murrell2015_NRMCB`).
