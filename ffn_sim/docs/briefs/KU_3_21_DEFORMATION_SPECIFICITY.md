# KU-3.21 — Deformation-specificity of α-actinin (dilation) vs filamin (shear)

> **Status**: CANDIDATE GATE. **NOT** PI-ratified. Born from the 2026-05-29 literature
> review of Luo, Mohan, Iglesias & Robinson 2013 (*Nature Materials* 12:1064-1071);
> formalises a new validation oracle for the H.3 cortex crosslinker system
> using Dictyostelium micropipette-aspiration experimental data.
>
> **Owner**: Lead session.
> **Prereq**: H.3 ✅ DONE (KU-3.5 PASS via R1) + Cell.build(with_lamellipodium=False)
> wiring (already in tree, commit `3306eb1`).
> **Companion**: [R1 rigid Lagrange exposure](RIGID_LAGRANGE_TENSION_DESIGN.md) —
> the per-bond-type γ itemisation work below builds on R1's λ exposure.

## 1. Motivation

Our H.3 cortex already includes the KU-3.19 crosslinker mix (~30% α-actinin /
~70% filamin via `cortex/crosslinkers.py`), but the current validation surface
is **static density only** (`test_crosslinkers.py` STATIC tests +
`tension_method_of_planes` aggregate γ). We have no gate that validates the
**deformation-specific binding response** of the two crosslinkers, which is the
key Luo 2013 finding and the discriminating experimental signature.

Luo et al. 2013 §Fig. 3 + Supp. Tables 2-5 establish:

| Crosslinker | Deformation it responds to | Pipette location | Mechanism |
|---|---|---|---|
| **α-actinin** (rod dimer, 35-40 nm) | **dilation** (Δ_V) | tip (max ΔV) | catch-bond, monotonic |
| **Filamin** (V-dimer, 90-100°) | **shear** (Δθ) | neck (max Δθ) | catch-bond + structural cooperativity (on-rate increases with V-angle change) |
| **Myosin II** | dilation (lever-arm dependent) | tip | catch-bond, accumulation kinetics ~30-60s |

Force-sharing fit across WT + mutants: **ζ = F_myosin / F_internal = 1/7**
(myosin carries ~14% of cortical tension, crosslinkers ~86%). This is an
**external literature oracle** we can use without needing to fit anything
ourselves (CLAUDE.md "no fitting to PI experimental data" — Luo's is published
literature, fair game).

## 2. Experimental setup (in-silico)

**System**: Single cortex (H.3 `Cell.build`, no lamellipodium, no FA) with
either WT crosslinker mix OR ablation (e.g. α-actinin-only, filamin-only) for
the discrimination experiment.

**Perturbation**: Micropipette aspiration — implemented as a localised external
force well of radius `R_pipette` (suction). Pressure ramped to 0.5-2.0 nN/μm²
(matching Luo Fig. 3 protocol).

**Measurement**: For each of N seeds, T snapshots over ~60-120 s simulated
time, record per spatial bin (pipette tip vs neck):

- α-actinin bound-count density (KU-3.21A)
- Filamin bound-count density (KU-3.21B)
- Local dilation Δ_V and shear Δθ (from actin neighbourhood deformation
  tensor, computed in `cortex/cortex.py` post-processing)

Plus γ split itemised by bond type (extension of R1):
- γ_myosin (myosin head-actin attach)
- γ_xlink_actinin (α-actinin attach)
- γ_xlink_filamin (filamin attach)
- γ_erm (ERM tether)
- γ_rigid (R1 SHAKE Lagrange — already exposed)

→ Sum should equal γ_total (consistency check).

## 3. Sanity Gate (Luo 2013 oracle bounds)

| # | Test | Pass criterion |
|---|---|---|
| 1 | Static density (KU-3.19, already PASS) | ρ_actinin = 0.30 · ρ_xlink_total ± 5%, ρ_filamin = 0.70 · ρ_xlink_total ± 5% |
| 2 | α-actinin accumulation at pipette **tip** (max dilation) | (ρ_actinin_tip / ρ_actinin_baseline) > 1.5 at 60 s aspiration (Luo Fig. 3e), no monotonic accumulation at **neck** |
| 3 | Filamin accumulation at pipette **neck** (max shear) | (ρ_filamin_neck / ρ_filamin_baseline) > 1.5 at 60 s aspiration (Luo Fig. 3f), distinguishable from tip response within ensemble σ |
| 4 | Myosin II accumulation at **tip** | (ρ_myosin_tip / ρ_myosin_baseline) > 1.5 at 30-60 s peak (Luo Fig. 3g); lever-arm length insensitive at our fixed lever arm |
| 5 | Force-sharing ratio ζ (R1 itemised γ) | ⟨ζ⟩ = ⟨γ_myosin / γ_total⟩ ∈ [0.10, 0.20] (Luo WT range ~0.14, ζ=1/7) — within 50% factor; tighter band requires lever-arm-matched myosin parameters |
| 6 | Ablation: α-actinin-only run | Filamin response at neck NOT observed (gate 3 fails by construction) |
| 7 | Ablation: filamin-only run | α-actinin response at tip NOT observed (gate 2 fails by construction) |

Tests #6, #7 are **negative controls** — they confirm the gate is actually
measuring what we claim, not a generic localisation artefact.

## 4. Implementation plan (3 incremental phases)

### Phase A — γ itemisation by bond type (foundation, ~3-4h)

Extend `cortex/cortex.py` `_tension_method_of_planes` (or its R1 sibling) to
return a `dict[bond_type → γ]` instead of a single scalar. Caller (R1
`_tension_method_of_planes_rigid` in `h3_ku35_tension.py` + new KU-3.21
driver) sums or partitions as needed. CPU regression: existing KU-3.5 keeps
PASS because γ_total = Σ_type γ_type bit-for-bit (same arithmetic, different
intermediate granularity).

### Phase B — Micropipette aspiration force module (~4-6h)

New `ffn_sim/perturbation/micropipette.py` (new package directory):
- `MicropipetteAspiration(hoomd.custom.Action)` — localised external force well
  with time-varying pressure ramp
- Public API: `attach(sim, center, radius, pressure_pa, ramp_s, duration_s)`
- Sanity Gate: dimensional analysis (pressure × area = force ✓), pressure → 0
  recovers no-perturbation BAOAB bit-for-bit, ramp_s → 0 = step input

### Phase C — KU-3.21 driver + sweep + Sanity Gate tests (~4-6h)

- `ffn_sim/scripts/h3_ku321_aspiration.py` — multi-seed, multi-ablation driver
  (WT / α-actinin-only / filamin-only / no-xlinker control)
- `ffn_sim/tests/validation/test_ku3_21_deformation.py` — STATIC (Sanity
  Gate 1-7 above) + opt-in PRODUCTION (`KU321_PRODUCTION=1`)
- Auto-viz: `h3_ku321_vis.py` produces α-actinin density at tip vs neck +
  filamin density at tip vs neck + ζ stacked-bar chart vs Luo WT band

## 5. Cost / benefit

| Cost | Benefit |
|---|---|
| ~10-16h Lead work across 3 phases | **Independent literature oracle** for KU-3.19 + KU-3.5 (no fitting required) |
| 1 new module (`perturbation/`) | Discrimination test: confirms our crosslinker mix actually responds to deformation, not just sits there |
| Phase A is `cortex/`-freeze territory; needs PI sign-off | ζ comparison to Luo 2013 quantifies how close our cortex parameters are to a published WT system |
| Phase B is greenfield, no freeze | Reusable infrastructure for any future external-perturbation study (compression, shear flow, AFM indentation) |

## 6. Decision branches (for PI)

| Option | What it does | What it buys |
|---|---|---|
| **D0** Defer | KU-3.21 stays as a written candidate; no implementation | No work. Risk: H.3 ✅ accepted without literature-oracle discrimination test for the crosslinker mix |
| **D1** Phase A only | γ itemisation lands; ζ in sweep_analysis becomes a proper Luo-comparable number | ~3-4h work, immediate ζ oracle improvement |
| **D2** Phase A + B + C | Full KU-3.21 gate | ~10-16h work, new validation surface — eligible to be a **new ratified KU gate** comparable to KU-3.5 |

**Default Lead recommendation, conditional on PI**: **D1 first** (low cost, high
information value — turns the current "γ_soft/γ_total proxy" panel in
sweep_analysis into an actual Luo-comparable ζ). Defer D2 (full gate) until
after H.3 ✅ DONE so it doesn't slip the critical path. D2 is Phase 1.5 / early
Phase 2 work.

## 7. References

- **Primary**: Luo T, Mohan K, Iglesias PA, Robinson DN (2013). "Molecular
  mechanisms of cellular mechanosensing." *Nature Materials* 12(11):1064-1071.
  doi:10.1038/nmat3772. PMC3838893.
- **Catch-bond single-molecule data**: Ferrer JM et al. 2008 (cited by Luo
  for α-actinin / filamin off-rate measurements via optical trap).
- **Companion design doc**: [RIGID_LAGRANGE_TENSION_DESIGN.md](RIGID_LAGRANGE_TENSION_DESIGN.md)
  — R1 already exposes the rigid-bond γ contribution; this doc extends to
  per-bond-type γ for itemised ζ.
- **Current implementation**:
  - `cortex/crosslinkers.py` (KU-3.19 30%/70% α-actinin/filamin)
  - `cortex/myosin.py` (Stam-Hocky bipolar minifilament)
  - `scripts/h3_ku35_tension.py` (KU-3.5 method-of-planes + R1 split)
