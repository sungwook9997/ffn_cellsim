# 3D Mikado fiber network — DESIGN (ECM geometry axis: 2D → 3D)

> **STATUS: DESIGN doc — 2026-06-02. Prep, NOT a contract, NO code.** Parallel session
> (`AUTONOMOUS_LOG_2026-06-02.md`). Specifies the 3D extension of the H.1 Mikado
> generator — the **geometry/dimensionality axis** of the 4-axis ECM decomposition
> (`ECM_PLATFORM_EXTENSIBILITY.md` axis 4). No edit to `ecm/mikado.py`. Default-2D ⇒
> current H.1 runs stay bit-for-bit identical.

## Motivation #1 — close the documented 2D-Mikado ⟨z⟩ gap (already in the code)

The H.1 sanity gate **already flags** a structural 2D-vs-3D defect: it reports
`⟨z⟩_measured − ⟨z⟩_KU-1.3` and warns *"2D Mikado lacks 3D …"* coordination
(`validation/oracles/common/sanity_gate.py:84,130`). The 2D thin-slab Mikado (box
`Lz=200 nm` with z-pinning; H.1 brief:43) cannot reach the biological 3D coordination
number `⟨z⟩` that KU-1.3 predicts — a known, logged gap. A genuine **3D fiber network
generator closes it** (the explicit validation win below). Motivation #2 = the ECM
condition matrix needs 3D matrices (3D col-I gel, tumor stroma, confinement) as config
presets, not rewrites.

## What exists (code-grounded)

- 2D generator `generate_2d_fiber_network` (`validation/oracles/ecm/fiber_network`),
  called by `build_mikado_state` / `resolve_derived` (`ecm/mikado.py:163,421`).
- Analytic 2D geometry oracle: Mikado segment length `ℓ_c = π/(2 ρ_L)` (Poisson
  spacing on a fiber; KU-1.27), `validation/oracles/common/derived_params.py`.
- Elastic band `G_0 ∈ [15, 200] Pa` (KU-1.30 #1; v1 32 Pa).
- Box convention: 2D thin slab (HOOMD is natively 3D; z is pinned/projected — H.1 brief).

So the platform is **3D-capable at the engine level** (HOOMD is 3D); the *generator* is
2D-only. The 3D extension is a generator + validation task, not an engine change.

## 3D generator design (additive; default-2D)

`dimensionality: 2D` (default) → exactly `generate_2d_fiber_network` (bit-for-bit).
`dimensionality: 3D` → a new 3D path:

1. **Fiber placement.** N rods of length L in a 3D box, positions uniform, orientations
   from a distribution with **order parameter S** (S=0 isotropic; S→1 aligned, for
   TACS-3 tumor stroma). (Nematic ⟨P₂⟩ control.)
2. **Crosslinks.** In 3D, exact rod–rod intersections are measure-zero — use a
   **capture radius** `r_x`: a crosslink forms where two fibers approach within `r_x`
   (segment–segment min-distance ≤ r_x). `r_x` from the physical crosslinker reach
   (e.g. the H.1 crosslink species length), NOT tuned. Crosslink areal/volumetric
   density ρ_x is a config knob (KU-1.V.2).
3. **Coordination number ⟨z⟩.** Each crosslink joins 2 fibers (local degree 4 like 2D),
   but the 3D embedding + finite `r_x` raises the *network* ⟨z⟩ toward the KU-1.3
   biological value — **the gap-closing observable** (validate that ⟨z⟩_3D > ⟨z⟩_2D and
   approaches KU-1.3).
4. **Rigidity / percolation.** Maxwell isostatic threshold: central-force `z_c = 2d`
   (2D→4, 3D→6); **bending-stabilized** fiber networks rigidify *below* `z_c`
   (bend-dominated regime). The generator must produce a **connected, above-rigidity-
   percolation** network; report the rigidity fraction. (Affine↔non-affine crossover
   sets whether G scales with stretching or bending — state which regime.)
5. **Elasticity.** 2D has the analytic `G_0(ρ,L)` oracle; **3D lacks a clean closed
   form** ("2D-oracle / 3D-runtime" — H.1 brief). Validate the 3D network modulus
   against the band `G_0 ∈ [15, 200] Pa` (KU-1.30) + 3D-collagen literature bands
   (KU-1.V.2), NOT a 2D analytic oracle. Density sweep ρ → G(ρ).

## Condition presets enabled (ECM_PLATFORM_EXTENSIBILITY condition matrix)

| Preset | dimensionality | S (alignment) | ρ, crosslink | KU |
|---|---|---|---|---|
| 2D reconstituted col-I (current) | 2D | 0 | ~1.5 mg/mL | KU-1.30 (default) |
| 3D col-I gel | 3D | 0 | density sweep | KU-1.V.2 + V.3 |
| Tumor stroma | 3D | **high (TACS-3)** | dense + LOX-crosslinked | KU-1.V.1/2 tumor |
| Confined channel | 3D + walls | — | pore 3–10 µm | KU-1.V.3 → H.9 |

## Sanity Gate (recorded now per CLAUDE.md)

1. **Dimensional.** L, r_x, ℓ_c [m]; ρ_x [1/m³] (3D) vs [1/m²] (2D — note the change);
   G [Pa].
2. **Boundary.** `dimensionality:2D` ⇒ `generate_2d_fiber_network` bit-for-bit (the
   default contract); S=0 ⇒ isotropic; S→1 ⇒ aligned (check the nematic order recovers).
3. **The ⟨z⟩ gap-closure check (the win).** 3D ⟨z⟩ > 2D ⟨z⟩ and trends toward
   `⟨z⟩_KU-1.3` — the sanity_gate.py-logged gap should *shrink* in 3D.
4. **Percolation / connectivity.** The generated network is connected and above the
   bending-stabilized rigidity threshold (report rigidity fraction; do not silently
   ship a floppy network).
5. **Numerical.** Segment–segment min-distance in 3D computed robustly (parallel /
   near-parallel rod degeneracy handled); RNG seed reproducible.
6. **Measurement consistency.** Report G at a stated density ρ + box size (finite-size
   effects); compare 3D G to 3D bands, never to the 2D analytic oracle.

## Validation acceptance (candidate — extends VG-H1 / KU-1.V.2)

| Gate | Criterion | KU |
|---|---|---|
| 2D default unchanged | `dimensionality:2D` bit-for-bit vs current | (regression) |
| ⟨z⟩ gap closes | 3D ⟨z⟩ → KU-1.3 (2D-Mikado gap shrinks) | KU-1.3 / KU-1.27 |
| 3D modulus band | G(ρ) ∈ 15–200 Pa; 3D-collagen lit bands | KU-1.30 / KU-1.V.2 |
| Alignment | nematic S recovers TACS-3 anisotropy | KU-1.V.2 tumor |

## Open items for PI / Lead

- [ ] Ratify the 3D generator + default-2D contract; pick its home (oracle `ecm/` vs
      runtime `ecm/`) consistent with where `generate_2d_fiber_network` lives.
- [ ] Confirm the crosslink **capture radius `r_x`** is physical (crosslinker reach),
      not tuned.
- [ ] 3D elasticity validation route (no clean 2D analytic oracle) — which 3D-collagen
      literature bands anchor G(ρ) in 3D.
- [ ] Sequence after H.9 (the confined-channel preset needs the nuclear pore limit) and
      the Mikado↔Cell coupling (the EXTEND "floating shell" fix) — both Lead-owned.
