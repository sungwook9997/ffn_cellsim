# FF S4 full integration — wiring the collagen ECM into the crawl loop as live DOFs (design, 2026-07-09)

**Status.** The S4 CORE is done + validated (coarse + native): `ff/fa_ecm.py` (FA↔collagen two-sided Newton
clutch) + the standalone `ff_ecm_remodel_demo.py` (FA traction recruits fibers, OFF-audit PASS). What remains is
the FULL integration — the crawling cell (S1+S2) and the collagen matrix (S4) evolving TOGETHER in the crawl
loop, so the cell crawls AND remodels its matrix simultaneously (the PI's real MCF7-on-collagen experiment). This
is the concrete blueprint; it was NOT executed overnight because it couples two networks with two integrators
inside the validated crawl loop — a change best made with PI oversight (subtle-regression risk), not unattended.

## What has to change (all additive, behind a new `--ecm` flag; default off = current crawl untouched)

1. **`build()`** — when `--ecm`: build a Mikado collagen slab UNDER the cell (`build_mikado_network`, box spanning
   the basal footprint × a depth below `z_sub`, far faces pinned Dirichlet). Store the ECM (`net.pos`, bend triples,
   crosslinks, segment pairs, pinned mask) in `S["ecm"]`. Attach the basal clutches to ECM nodes with
   `fa_ecm.attach_clutches_to_ecm(actin_pos, ecm_pos, capture)` → `S["ecm_node"]` (per-clutch ECM index).

2. **`run()`** — when `--ecm`, the basal FA clutch changes target from the FIXED substrate anchor to the LIVE ECM
   node: REPLACE `clutch_spring_kernel(... anch_d ...)` with `fa_ecm.clutch_ecm_spring_kernel(cell_pos, ac_d,
   ecm_pos, ecm_node_d, k_int, rest, cell_f, ecm_f)` — +f on the actin (unchanged traction on the cell), −f on the
   ECM node (the new reaction). The cell's implicit solve is otherwise UNCHANGED (this preserves the validated
   crawl path — the clutch force on the cell is the same magnitude; only its anchor is now a moving DOF).

3. **ECM evolution (a separate explicit sub-loop each crawl step)** — the ECM is soft + small-dt; evolve it with
   an explicit overdamped sub-step (as in the remodel demo): zero `ecm_f`; `cytosim_bending_kernel` +
   `link_spring_kernel`(crosslinks) + `link_spring_kernel`(segments) + the clutch reaction (from step 2) →
   `axpy_physical_kernel(ecm_pos, dt_ecm, ecm_gamma, ecm_f)`; pin via `ecm_gamma[pinned]=∞`. Run `n_ecm_sub`
   micro-steps per crawl step so the ECM (CFL dt ~µs) keeps pace with the crawl's large implicit dt (the same
   operator-split sub-cycling pattern as the DCM cadherin KMC). `dt_ecm` = CFL of the ECM (bending+crosslink+seg).

4. **Clutch treadmill re-attach onto ECM** — when a clutch turns over (kmc), re-run `attach_clutches_to_ecm` for
   the newly-nascent clutches (they grab whatever fiber node is now under them), so as the cell crawls it grips
   fresh collagen (the treadmill onto a live matrix). Rear clutches release their ECM node (`ecm_node=-1`).

5. **Recording + viz** — save the ECM frames alongside the cell frames; the crawl viewer gets an ECM fiber scene
   (the matrix deforming as the cell crawls over it).

## Validation gates (write before running)

- **OFF-path invariance**: `--ecm` off ⇒ byte-identical to the current crawl (the fixed-substrate clutch path).
- **Clutches-OFF audit**: with `--ecm` + clutches off ⇒ neither the cell translocates NOR the ECM remodels.
- **Traction-driven remodel**: `--ecm` + clutches on ⇒ the ECM deforms under the crawling cell (fibers recruited
  along the traction footprint), and net COM drift is still traction-driven (the S2 audit still passes).
- **Momentum/Newton**: the clutch is a Newton pair (cell +f / ECM −f), so no spurious net force enters the
  combined cell+ECM system (the `test_fa_ecm` gate, now at loop scale).
- **Coarse-first**: validate at coarse resolution (fast, mac) before any native (A5000) run — the crawl's own
  grid-dependent-drag caveat (`FF_CRAWL_DIAGNOSIS_2026-07-09`) still applies to the native COM speed, so the
  native cell-on-collagen crawl inherits that open solver item.

## Risk + why it was deferred

Two networks / two integrators (cell implicit large-dt + ECM explicit small-dt) sub-cycled and coupled inside the
one crawl loop is a real integration with regression surface on the validated crawl path. It is flag-gated so the
default path is safe, but getting the sub-cycling cadence + the treadmill re-attach + the recording right is a
multi-hour, oversight-worthy change — the natural first task for the next working session (or PI-directed), with
the S4 physics already de-risked by the primitive + the standalone remodel demo.
