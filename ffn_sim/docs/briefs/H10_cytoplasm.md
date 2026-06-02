# H.10 — Cytoplasm viscoelasticity (DRAFT skeleton)

> **STATUS: DRAFT skeleton — 2026-05-30. NOT a contract.** Pending PI ratification
> of the EXTEND decision (`ffn_sim/docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`). KU
> anchors: Notion KU v2 layer (KU-3.B3.x). Read-only research session; no code written.
>
> **NOTE (EXTEND §B.5):** this is the ONE compartment that is NOT "just another
> custom force." A viscoelastic/poroelastic cytoplasm is a background medium with
> memory — it likely needs integrator-adjacent design (memory-kernel drag or filler
> particles), beyond Templates 1–2. Scope separately and expect a frozen-integrator
> (`integrator/`) PI-gate.

**Budget**: TBD
**Branch**: TBD (after H.9)
**Owner**: Lead
**Prereq**: H.3 cortex; EXTEND ratified.
**Reference**: EXTEND memo §B.5; Notion KU-3.B3.

## Goal

Represent the cytoplasm as a viscoelastic / poroelastic medium (currently only a
water-viscosity Stokes drag `γ_b` exists). Cytoplasm viscosity is the most sensitive
discriminator of metastatic potential (MCF7 ≈ 5.3× MDA-MB-231).

## Deliverables (proposed)

| File | Content |
| --- | --- |
| `cell/cytoplasm.py` | Viscoelastic background: memory-kernel / per-type effective drag, or biphasic poroelastic filler. **Design TBD — integrator-adjacent.** |
| `configs/phase1_h10.yaml` | KU-3.B3 params (η, G, Dp). |
| `tests/test_cytoplasm.py` | Sanity gate + viscosity-contrast oracle. |
| `outputs/h10/REPORT.md` | cytoplasm viscoelastic evidence. |

## KU anchors (SI — KU-3.B3)

- Cytoplasm viscosity **η** (MRS, L=3±1 µm wires): MCF10A 41.6 / MCF7 56.4 / MDA 10.7 Pa·s. **MCF7 ≈ 5.3× MDA** (the discriminator). η probe-length-dependent (η ∝ L²) — state wire length.
- Cytoplasm elastic modulus **G ≈ 30–80 Pa** (MCF10A 79.3 / MCF7 32.9 / MDA 38.6) — ≈10–20× below whole-cell AFM (DISTINCT quantity).
- Poroelastic diffusion **Dp ≈ 40–60 µm²/s** (4–6×10⁻¹¹ m²/s); pore ξ≈14 nm; poroelasticity dominates for events faster than ~0.5 s.

## Attach contract

NOT Templates 1–2 alone. Options: (a) per-type / memory-kernel effective drag in a
BAOAB-compatible Action (touches the drag model), (b) explicit poroelastic filler
particles (Template 2 + cytosol-redistribution). **Requires design + likely an
`integrator/`-freeze PI-gate.**

## Validation acceptance (candidates)

| Gate | Criterion | KU |
| --- | --- | --- |
| Cytoplasm viscosity contrast | MCF7 η ≈ 5.3× MDA | KU-3.B3.1 |
| Cytoplasm modulus band | G ∈ 30–80 Pa, ≪ whole-cell AFM | KU-3.B3.1 |
| Poroelastic timescale | relaxation faster than ~0.5 s; Dp 40–60 µm²/s | KU-3.B3.2 |

## Friction

The genuinely harder module (EXTEND §B.5). Frozen-integrator changes need PI sign-off.
