# Parallel-prep INDEX — 2026-06-02 (handoff to PI / Lead)

> Single entry-point for the conflict-free prep produced by the parallel autonomous
> session of 2026-06-02 (narrative + contract: `AUTONOMOUS_LOG_2026-06-02.md`). All
> artifacts are **NEW files, untracked, import-isolated** — the Lead's working tree,
> the live runtime, and the existing test suite are **bit-for-bit unaffected**. Nothing
> was committed, no Lead-owned file edited, no `integrator/` touched, no Notion written.
> The Lead/PI integrate + commit at a sequencing point.

## Artifacts

### Runnable, tested building blocks (ready to wire)

| File | Purpose | Tests | Wires into (Lead) |
|---|---|---|---|
| `bridge/ligand_species.py` | ECM ligand identity → integrin clutch kinetics (col-I α2β1 slip, laminin-111 α6β1 slip-proxy, FN α5β1 catch-slip); `PereverzevParams` swap, `F_s=k_BT/x_β` | `tests/test_ligand_species.py` (**24**) | `bridge/fa.py` `resolve_h4` (see spec) |
| `bridge/clutch_spatial.py` | β1-distribution observable (`f_edge`/`angular_cv`/`n_engaged`) over engaged-clutch positions; layout-agnostic | `tests/test_clutch_spatial.py` (**10**) | analysis over production GSD (read-only) |

**34 tests, all GREEN in isolation** (pure-Python/numpy, no HOOMD — cannot red the
Lead's collection/regression; import-isolated — no existing file imports them).

### Design docs (read-only prep; specify a PI-sequenced build)

| Doc | What it specifies | Feeds | PI-gated build? |
|---|---|---|---|
| `H10_CYTOPLASM_DESIGN.md` | cytoplasm viscoelasticity, 3-tier (Tier-1 η_eff drag = MCF7≈5×MDA discriminator; Tier-2 GLE; Tier-3 poroelastic filler) | MC-H10 `draft→implemented`; cell-type η | Tier-1 (magnitude/CFL); Tier-2/3 (integrator) |
| `LIGAND_IDENTITY_FA_SPEC.md` | col-I/laminin/FN ligand identity; exact `resolve_h4` attach; KU-2.5 skip for slip | Layer-1 single-cell PI-exp validation | no (additive, default-FN) |
| `BETA1_DISTRIBUTION_OBSERVABLE.md` | β1-IF analog (diffuse/peripheral/uniform) metrics + acceptance | Layer-1 PI-exp validation | no (read-only analysis) |
| `CANCER_CELLTYPE_PARAM_MAP.md` | verified bands → sim knobs per cell-type (MCF10A/MCF7/MDA); **no global modulus** | cell-type instantiation | (cell-type contracts deferred) |
| `H8_HELFRICH_CURVATURE_DESIGN.md` | discrete mean-curvature operator for the deferred κ_m bending; analytic sphere test | `cell/membrane_surface.py` κ_m TODO | CFL softening (k_ERM precedent) |
| `MIKADO_3D_DESIGN.md` | 2D→3D fiber network (closes the logged 2D-Mikado ⟨z⟩ gap); alignment S; rigidity | ECM geometry axis / 3D presets | no (default-2D bit-for-bit) |
| `COMPOSITE_TENSION_REVALIDATION_DESIGN.md` | KU-3.5/3.1 composite re-attribution via A/B/C toggle of H.8/H.9; superposition check | VG-H3-composite (`blocked`) | **yes** (gate-contract relabel) |

## Suggested integration order (Lead)

1. **After the literal-v0 γ run lands** (don't perturb the active driver):
   wire `ligand_species.py` into `fa.py` `resolve_h4` (default `fibronectin` ⇒ bit-for-bit)
   → run preset (b) col-I single MCF7 → add the `clutch_spatial` β1 readout → **Layer-1
   single-cell PI-exp validation** (traction + β1 pattern vs lit bands; PI poster overlay).
2. **H.8 composite-tension run** (membrane already integrated) → the A/B toggle of
   `COMPOSITE_TENSION_REVALIDATION_DESIGN` → γ_mem isolated, KU-3.5 re-attributed.
3. **H.9 nucleus integration** (deferred Template-2) → A/B/C run completes KU-3.1
   composite + unblocks VG-H3-composite (PI gate-relabel).
4. **H.10 Tier-1** (η_eff drag, PI-gate magnitude) → cytoplasm viscosity discriminator
   → `CANCER_CELLTYPE_PARAM_MAP` per-compartment cell-type presets.
5. **ECM geometry**: 3D Mikado generator (default-2D) → 3D-gel / tumor-stroma presets.

## PI decision roll-up (aggregated from all docs)

- **H.10 Tier-1 η_eff** magnitude + a CFL/timestep gate re-run at the raised drag; confirm
  BAOAB reads γ per-type for the noise term (FDT under per-type η_eff).
- **Ligand identity**: ratify col-I default = single-cell anchor (Taubenberger 1.3 s⁻¹ /
  0.23 nm) vs locked-open; guard the KU-2.5 catch-peak gate to **skip slip** ligands.
- **Composite re-validation** (gate-contract): ratify A/B/C protocol + the superposition
  check; relabel current KU-3.5/3.1 "cortex-attributed"; unblock VG-H3-composite.
- **Cell-type**: ratify per-compartment (no global cancer modulus); fill the MCF7/MDA
  **lamin/E_nuc gap** from literature (cheap); method-tag the disputed stiffness ordering.
- **3D Mikado**: ratify default-2D contract; confirm the crosslink capture radius is
  physical; pick the 3D-elasticity validation route.
- **H.8 Helfrich**: ratify Tier-A(measure)→Tier-B(leading force) staging; CFL softening.
- **Frozen-integrator**: H.10 Tier-2/3 remain PI-gated (no autonomous build).

## Provenance / discipline

Every kinetic/modulus constant is literature-sourced (PI-exp dossiers + KU bands; cited
in each doc), never tuned. `F_s=k_BT/x_β` and the toggle superposition are physics, not
free parameters. No fitting to PI data (literature-first; PI poster is overlay-only).
Cross-references the Notion KU v2 / Contract Graph (KB-PIV, MC-H8/9/10, VG-H3-composite)
without duplicating it.
