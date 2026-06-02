# H.9 — Nucleus mechanics (DRAFT skeleton)

> **STATUS: DRAFT skeleton — 2026-05-30. NOT a contract.** Pending PI ratification
> of the EXTEND decision (`ffn_sim/docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`). KU
> anchors: Notion KU v2 layer (KU-3.B2.x). Read-only research session; no code written.

**Budget**: TBD
**Branch**: TBD (after H.8)
**Owner**: Lead
**Prereq**: H.3 cortex; EXTEND ratified. Couples to KU-1.V.3.3 (nuclear pore limit).
**Reference**: EXTEND memo Part B (Template 2 + Template 1); Notion KU-3.B2.

## Goal

Add a nucleus module (additive + default-off): a stiff nuclear element with
lamin-A/C-set shell stiffness + chromatin interior, deformable under confinement.
The nucleus is the rate-limiting mechanical element in confined migration and is
mechanically co-dominant in cancer cells.

## Deliverables (proposed)

| File | Content |
| --- | --- |
| `cell/nucleus.py` | `nucleus_bead` particle group near centroid + `NucleusConfinement(hoomd.md.force.Custom)` shell-stiffness force; optional post-BAOAB centroid pin. |
| `configs/phase1_h9.yaml` | KU-3.B2 params (E_nuc, lamin-A scaling, radius). |
| `tests/test_nucleus.py` | Sanity gate + nucleus:cytoplasm ratio + confinement limit. |
| `outputs/h9/REPORT.md` | nuclear modulus + confined-migration evidence. |

## KU anchors (SI — KU-3.B2)

- Nuclear modulus **E_nuc ≈ 1–10 kPa** (in-situ ~5, isolated ~8); cytoplasm ~0.5–1 kPa.
- Nucleus:cytoplasm ratio **1.4–5× in-situ / 3–10× isolated** (technique-dep; do not hard-code 10×).
- **Lamin-A scaling**: lamin-A ~ tissue-E^0.7; viscosity ~ [lamin-A]^(3±1); elasticity ~ [lamin-A]^(~0.5). Lamin-A → viscosity, lamin-B → elasticity.
- **Confined-migration limit**: arrest at nucleus compressed to ~10% cross-section; critical pore ~7 µm² (tumor) / 4 (T-cell) / 2 (neutrophil).
- Cancer: lamin-A/C down → softer/more-deformable → invasion↑ (context-dependent).

## Attach contract (EXTEND PoC)

Template 2 (new `nucleus_bead` type pre-snapshot + `gamma_map` entry + optional pin)
+ Template 1 (`NucleusConfinement` custom force over nucleus tags). Two-regime model:
chromatin interior (small strain) + lamin shell (strain-stiffening). NO integrator
change; default-off ⇒ pre-H.9 bit-for-bit identical.

## Validation acceptance (candidates)

| Gate | Criterion | KU |
| --- | --- | --- |
| Nuclear modulus | E_nuc ∈ 1–10 kPa | KU-3.B2.1 |
| Nucleus:cyto ratio | 1.4–5× in-situ | KU-3.B2.1 |
| Lamin-A stiffness scaling | viscosity ~ laminA^(3±1) | KU-3.B2.2 |
| Confined-migration arrest | ~10% cross-section / pore 7 µm² | KU-3.B2.3 / KU-1.V.3.3 |
| Composite rounding | KU-3.1 re-derived with nuclear term | KU-3.1 (re-validate) |

## Friction

Quasi-rigid nucleus needs post-BAOAB pin (no BAOAB integration-group exclusion).
CFL stiffness gate (PI). MMP/pore-enlargement rescue branch couples to ECM (KU-1.V.3.3).
