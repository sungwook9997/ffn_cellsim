# FF Engine Re-Validation — REPORT

## ⚠⚠ RETIRED-BUILD BANNER (added 2026-07-25 — read before citing any part result below)

> This banner is deliberately its own `##` section, not part of the title block: the numerals it quotes
> are **facts from another document and a retraction notice**, NOT headline claims of this report, and
> the `kb_coverage` results ratchet fingerprints the title block. Do not move it above the heading.
>
> Every part in this report that builds a cortex through `ff/gamma_floor.build_crosslinked_cortex`
> (P4.1 turgor, P7.F.7 native-cortex ⊗ Biot ⊗ IBM, and every γ measurement) was **measured on the
> pre-2026-07-23 cortex build**; it has **not been re-run**; see
> **`../../docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md`** (5 structural defects).
>
> **Three of the five are still live** on that builder (verified 2026-07-25): `n_xl` is a caller
> argument and `ff/architecture_spec.CORTEX` is never read (all call sites pass `n_xl = n_filaments` ⇒
> crosslink density 1.0, the under-percolating relic `bc5ff3b0` replaced with 20 in
> `architecture_spec.py` only); `resolve_overlaps` is never passed (default `False`); `length_dist`
> stays `"mono"`.
>
> Scope honesty: the parts here are **identity / conservation / shape** oracles (e.g. P4.1 asserts
> γ = ΔP·R/2 exactly; P7.F.7 asserts Σgrid·h³ = Σnode_src). Those are Π₀- and connectivity-**invariant**,
> so the PASS verdicts are not overturned by the retired build — but any **magnitude** read out of this
> report (γ values, `N_cortex_nodes`, forces) is on the retired cortex and must be re-run. Π₀ used
> here is 40 Pa, a **CONVENIENCE (HeLa-proxy)** value, now gated
> (`ffn_sim/common/params_turgor.yaml`).

Ground-up part-by-part re-validation (PI 2026-07-16). `validate` = KB-fidelity (analytic oracle) + phenomenon-completeness, NOT experimental. All native params.


| Part | Name | Result | Oracle | max rel-err |
|---|---|---|---|---|
| P7.F.7 | Native cortex ⊗ Biot ⊗ IBM composition — multiscale, momentum-conserving | ✅ PASS | `Σgrid·h³ = Σnode_src (momentum conserved across the coupling); coarse fluid grid spans cell (dx≪R)` |  |

## Details

### P7.F.7 — Native cortex ⊗ Biot ⊗ IBM composition — multiscale, momentum-conserving  ✅
- **Oracle (KB-fidelity):** Σgrid·h³ = Σnode_src (momentum conserved across the coupling); coarse fluid grid spans cell (dx≪R)
- **Phenomena:** momentum_conserved_at_native_scale ✓, spread→interp_bounded ✓, correct_multiscale(dx≪R,grid_spans_cell) ✓
- **Metrics:** `{"N_cortex_nodes": 21000, "R_um": 7.498012936753268, "fluid_grid": "40\u00b3", "dx_fluid_um": 0.6, "momentum_rel_err": 8.185656005822488e-15, "device": "cpu"}`
- the FSI infrastructure composes with the real cortex at the right multiscale (fine fibers + coarse smooth fluid grid), engine UNTOUCHED. FINAL step = wire into network_warp.py (additive, default-off→physiological-ON) + retire 6πηR — the deliberate double-count-trap-aware commit.
