# ⚠⚠ RETIRED-BUILD BANNER — `outputs/engine_reval/` (added 2026-07-25)

**Applies to every file in this directory, including `results.json` and `figs/*`.**

> Every part of this re-validation that builds a cortex through
> `ff/gamma_floor.build_crosslinked_cortex` was **measured on the pre-2026-07-23 cortex build**; it has
> **not been re-run**; see `../../docs/v2_audit/CORTEX_STRUCTURE_AUDIT_2026-07-23.md` (5 structural
> defects).

## Why `results.json` carries no inline banner

`results.json` is a JSON **array** of `PartResult` records that downstream readers iterate. Inserting an
annotation object would change its element type and could break those readers, so the file was left
**byte-untouched** and this sibling banner covers it. `REPORT.md` in this directory carries the inline
banner.

## The three still-live defects (verified 2026-07-25 against the code, not a doc claim)

`ff/gamma_floor.build_crosslinked_cortex`:

1. **`n_xl` is a caller argument and `ff/architecture_spec.CORTEX` is never read.** Every in-tree call
   site passes `n_xl = n_filaments` (`PROD_N_XL == PROD_N_FIL == 70686`) ⇒ crosslink density
   **1.0 / filament**, the ×40 relic that left the native network at mean degree 2 (giant component
   79.6 %, 11,430 fragments, below the ln(N) = 11.2 percolation threshold). Commit `bc5ff3b0` raised
   `density_per_fil` 1.0 → 20 in `architecture_spec.py` **only**; that commit's claim that "both weave
   and gamma_floor read the same spec" is **false for this builder**.
2. **`resolve_overlaps` is never passed** to `build_cortex_network`, so the cortex is built at the
   default `False`. Only `ff.weave.weave` passes it.
3. **`length_dist` is left at `"mono"`.**

`bc5ff3b0`'s own message names under-percolation as "the upstream root of … **no myosin transmission**".

## What this does and does not invalidate

* **Not invalidated:** the identity / conservation / shape oracles (P4.1 asserts γ = ΔP·R/2 exactly;
  P7.F.7 asserts Σgrid·h³ = Σnode_src). Those are connectivity- and Π₀-**invariant**, so the PASS
  verdicts stand as *identity* checks.
* **Invalidated until re-run:** every **magnitude** — γ values, `N_cortex_nodes`, forces, `R_um` —
  because they are read off the retired cortex.
* **Additionally:** the turgor used is Π₀ = 40 Pa, a **CONVENIENCE (HeLa-proxy)** value. Π₀ is now a
  gated, no-silent-default parameter — `ffn_sim/common/params_turgor.yaml` (`value: null`,
  `evidence_status: "GAP — PI"`). Three values are live in the tree (40 / 133 / ~72 Pa) and none is an
  MCF7 measurement.

## Removal condition

Delete this banner only after a native A5000 re-run of the affected parts on the **fixed** cortex
(crosslink density 20, `resolve_overlaps=True`, exponential `length_dist`). The executable check
`ffn_sim/tests/ac/test_turgor_pi0_gate_and_banners.py::test_three_2026_07_23_defects_are_still_live_in_the_gamma_floor_builder`
starts **failing** once the builder is fixed — that failure is the signal to re-run and then remove the
banner, not to weaken the test.
