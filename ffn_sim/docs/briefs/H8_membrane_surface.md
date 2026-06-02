# H.8 — Plasma-membrane surface mechanics (DRAFT skeleton)

> **STATUS: DRAFT skeleton — 2026-05-30. NOT a contract.** Pending PI ratification
> of the EXTEND decision (`ffn_sim/docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`). KU
> anchors live in the Notion KU v2 layer (KU-3.B1.x). Created in a read-only
> research session; no code written.

**Budget**: TBD
**Branch**: TBD (cut after H.4 FA integration lands)
**Owner**: Lead (post-FA)
**Prereq**: H.3 cortex (ERM tether = membrane–cortex linker seed); EXTEND ratified.
**Reference**: EXTEND memo Part B (attach contract, Template 1); Notion KU-3.B1.

## Goal

Add a plasma-membrane **surface** mechanics module (additive + default-off): in-plane
tension γ_mem, Helfrich bending κ_m, area-expansion modulus K_A, membrane–cortex
adhesion γ_MCA. Supplies the missing whole-cell surface tension so KU-3.5 becomes
**composite** (γ_total = γ_membrane + γ_cortex), correcting the current
cortex-only attribution.

## Deliverables (proposed)

| File | Content |
| --- | --- |
| `cell/membrane_surface.py` | `MembraneSurfaceTension(hoomd.md.force.Custom)` over cortex shell tags + `attach_membrane_surface(...)` with CFL gate. Fully-dynamical lipid sheet (new particle type) = later refinement (Template 2). |
| `configs/phase1_h8.yaml` | KU-3.B1 params (T, κ_m, K_A, γ_MCA). |
| `tests/test_membrane_surface.py` | Sanity gate + tether-force oracle. |
| `outputs/h8/REPORT.md` | composite-tension evidence. |

## KU anchors (SI — KU-3.B1)

- Apparent tension **T ≈ 3×10⁻⁵ N/m** (band 0.03–0.3 mN/m); decomposition **T = T_m + γ_MCA**.
- Bending **κ_m ≈ 1×10⁻¹⁹ J** (10–30 k_BT).
- Area modulus **K_A ≈ 0.24 N/m**; lysis 3–10 mN/m, max strain 2–5%.
- MCA adhesion **γ_MCA ≈ 10⁻⁵ J/m²** (10⁻⁶–10⁻⁴).
- Tether **f_t ≈ 5–40 pN**; **f_t = 2π·√(2·κ_m·(T_m+γ_MCA))**.

## Attach contract (EXTEND PoC, Template 1)

Custom force over cortex shell tags `[0, n_cortex_actin)` — mirrors
`enclosed_volume.py` / `erm.py`: read `cpu_local_snapshot`, mask shell tags, write
`cpu_local_force_arrays`; CFL gate `τ = γ_b/k_eff`. Promote ERM to the
membrane–cortex linker. NO integrator change; default-off ⇒ pre-H.8 builds
bit-for-bit identical.

## Validation acceptance (candidates — Notion Validation Matrix)

| Gate | Criterion | KU |
| --- | --- | --- |
| Apparent tension | T ∈ 0.03–0.3 mN/m | KU-3.B1.1 |
| Tether force | f_t ∈ 5–40 pN via f_t=2π√(2κ(T_m+γ_MCA)) | KU-3.B1.4 |
| Composite cortex tension | γ_total = γ_mem + γ_cortex; KU-3.5 re-derived | KU-3.5 (re-validate) |
| Cancer overlay | metastatic ~2× lower T | KU-3.B1.5 |

## Friction (EXTEND §B.5)

CFL stiffness gate may require softening (PI-gate, precedent: k_ERM 1000× softening
2026-05-26) — surface to PI, do not soften silently.
