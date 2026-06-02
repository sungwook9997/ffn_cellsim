# Layer-2 multicellular spheroid line — REPORT (L2.0 + L2.1)

> Status 2026-06-02. Parallel CBM (center-based, 1 particle/cell) spheroid line, isolated
> from the single-cell main line, on the shared HOOMD + frozen-BAOAB stack. Brief:
> `docs/LAYER2_MULTICELL_DESIGN.md`. Anchor provenance: `docs/LAYER2_ANCHORS_2026-06-02.md`.

## Milestones

| Phase | What | Status |
|---|---|---|
| **L2.0** | measurement + acceptance-oracle + config layer (parameter-free) | ✅ DONE, 21 unit tests green |
| **L2.1** | CBM physics builder (Morse + reused BAOAB) + **G1 stable-aggregate gate** | ✅ DONE, G1 PASS (200 cells) |
| L2.2 | active motility + substrate spring → spreading; A/A₀(R) sweep | ⬜ next |
| L2.5 / L2.6 | cadherin catch-bond (KU-4.2) upgrade; 3D-Mikado invasion | ⬜ |

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
- `figs/fig_layer2_aa0_law.png` — **L2.3 first emergent A/A₀(R₀) spreading law.** Edge-directed
  active-wetting traction (6 nN, Lp=11 µm) vs cohesion (6.5 nN), swept over initial radius
  R₀=29–78 µm. The CBM produces a measurable (Δ≈7.6%), non-monotonic A/A₀(R₀) — small clusters
  (high edge fraction) and intermediate sizes spread, large cohesion-dominated ones do not —
  fittable to the PI's a + b/R + c/R² form (peak ≈42 µm; r²=0.67 single-seed). This is the
  first quantitative contact with the experiment's collective term. G3-clean (r²≥0.95) needs
  ensemble averaging (multiple seeds/R₀) to suppress single-realization noise. PI A/A₀ values
  are overlay-only at comparison time (not shown — we hold only the functional form).

## Verification

- `tests/test_spheroid_observables.py` (13) + `tests/test_layer2_params.py` (8) = **21 green**
  (synthetic clouds vs closed-form oracle; resolve derivations vs Magic-Number-Block).
- `scripts/layer2_g1_smoke.py` — reproducible G1 run (`python -m ffn_sim.scripts.layer2_g1_smoke`).
- Isolation: runtime imports NO oracle (hard rule); additive new files only; single-cell
  main line + `integrator/` freeze untouched; no git commit yet.

## Open (PI ratification — not blocking L2.2)

1. D_e derivation `N_cad·⟨F⟩·Δx*` vs recovering Omidvar 2014 (gbook) for a measured value.
2. Surface-tension validation target: emergent-only vs non-MCF7 proxy (MCF10DCIS ~21 mN/m).
3. Cell-size band position: 15 µm (low end) vs 17–18 µm.

## Next (L2.2)

Add per-cell active traction (external force, magnitude from the single-cell scale-bridge)
+ the z=0 substrate spring (`ecm/substrate.py`, reused) → spreading; then the R₀ sweep →
fit A/A₀ = a + b/R + c/R² (G3) and overlay the PI poster.
