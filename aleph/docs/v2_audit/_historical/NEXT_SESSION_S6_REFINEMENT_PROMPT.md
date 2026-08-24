---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Next-session boot prompt — S6 physiological matrix remodeling

*(Paste this to boot the next session. Then run the boot protocol and continue S6.)*

---

You are continuing the `ffn_cellsim` FF single-cell mechanobiology axis. The **crawl axis is CLOSED** — a real
polarized crawl is research-level and dramatic single-cell crawling is not MCF7 (epithelial, poorly migratory)
physiology (see `FF_POLARIZED_CRAWL_DESIGN_2026-07-09.md`). **Do NOT re-chase the crawl.** The live axis is **S6 —
the cell ↔ collagen-I matrix interaction (remodeling)**, which is MCF7's real mechanobiology and matches the PI's
pV4D4/collagen-I experiment.

## Boot (minute 0)
1. Read `CLAUDE.md` FULLY — especially the **HARD rule: validate/debug/visualize NATIVE + FULL compartments only,
   never coarse/stripped** (a coarse ~250-filament cell collapses on its own = a filament-count artifact that
   burned hours). Native = `--cortex-fil 38000 --from-resting --microtubules` (Nc≈266k).
2. `git log --oneline -8`; confirm branch `dcm/main`, latest FF commits `a95bf72 · 2414bcb · 1994e33`.
3. `conda activate ffn_sim`; GPU work → gbook A5000 (`~/ff_scratch`; rsync edited files; run with full path
   `~/miniconda3/envs/ffn_sim/bin/python`; **PYTHONUNBUFFERED=1** so the log isn't buffered; monitor by LOG FILE).
4. Skim `FF_TREADMILL_CRAWL_DESIGN_2026-07-09.md` §S6 + the ECM co-sim code in `ff_crawl_on_substrate.py`.

## S6 state (validated, native)
The native cell's basal clutches grip a Mikado collagen-I slab (`--ecm`); a co-sim substeps the compliant collagen
(bending + crosslink + segment-inextensibility, pinned bulk BC, 20 substeps/cell-step) under the cell's clutch
traction. **Validated native** (266k cell + 3000-fiber/33k-node collagen, 21 min):
- adheres — 10953/10953 clutches on collagen; traction deforms the matrix **gripped-node 266 nm, peak 586 nm**;
- **traction-driven** — clutch-OFF audit → 0 (`mask_ecm_by_bound_kernel` gates only BOUND clutches; coarse ON/OFF 468/0);
- STABLE ADHERED, V/V₀ = 1.000. Runs: `--ecm` prints `[ECM-REMODEL]`, npz stores `ecm_pos0/ecm_posf/ecm_node/ecm_foff`.

**Honest caveat (the S6 task):** the deformation is present but **NOT yet coherent physiological remodeling** —
the coherent inward pull is only **+15.7 nm** of the 266 nm (the collagen follows the basal clutches as the cell
settles on the compliant matrix; not inward densification/alignment).

## S6 task — make the remodeling PHYSIOLOGICAL
Goal: coherent **inward densification + fiber alignment** toward the cell, traction-driven, native, visualized —
not clutch-following. Three levers:

1. **Two-way coupling (cell FEELS the matrix).** Currently one-way (cell→collagen). Add the cell-side
   `clutch_ecm_spring` force to `f_d` at the 3 cell force-assembly sites (`clutch_spring_kernel` launches at
   ff_crawl_on_substrate.py ~364 full_force / ~429 gpu_force_fn / ~503 main-body) — gate `ecm_on and clutches`,
   inputs `[pos_d, ac_d, Ep_d, en_d, k_int, rest, f_d, Edummy_e]` (collagen side discarded there; `Edummy_e` already
   allocated). **Move the `mask_ecm_by_bound_kernel` launch to the TOP of the loop** (before the force `_zero`) so
   `en_d` is current for the force assembly. This makes the cell traction depend on the compliant collagen (→ biphasic).
2. **Coherent contractile inward pull.** The cell's myosin contraction must pull the clutches (and collagen)
   coherently INWARD, not the disorganized basal motion. Check: is myosin strong enough / is the basal cortex
   contracting inward? Consider a rear/peripheral-enriched or basal-focused contraction so the footprint pulls the
   matrix in. (Physiological values, no tuning-to-pass.)
3. **Biphasic cross-check (Chan-Odde).** Sweep collagen stiffness (κ_collagen / xl_k, or an effective substrate E)
   → traction should peak at the ~2–300 kPa optimum (Bangasser 2013). A legit *cross-check* vs lit, not a gate to pass.

## Gates (visual-first, native)
- **inward recruitment** coherent (≫ +15 nm; collagen densifies toward the cell, not just follows clutches);
- **fiber alignment** near the cell (collagen fibers reorient toward the footprint);
- **traction-driven** — clutch-OFF `[ECM-REMODEL]` ≈ 0 (add the remodel print to the `--audit` branch too);
- **biphasic** traction vs collagen stiffness (cross-check);
- **native + full** (Nc≈266k, all compartments) + **visualized** — displacement quiver + the remodel figure +
  the interactive cell/stress viewer. `browser_check.py` to verify HTML renders.

## Key files
- `ff/fa_ecm.py` — `clutch_ecm_spring_kernel` (two-sided), `mask_ecm_by_bound_kernel`, `attach_clutches_to_ecm`.
- `ff/ecm_mikado.py` — `build_mikado_network` (κ, xl, seg, pinned BC).
- `scripts/ff_crawl_on_substrate.py` — ECM co-sim setup (~308, `ecm_on`, `_Ensub=20`, `Edummy_e`), collagen substep
  (~632, before record), remodel print/save in `main()`.
- `scripts/ff_ecm_remodel_demo.py` — the validated standalone remodel (reuse helpers `_segment_pairs`, `K_SEG`).

## Native S6 run command (template)
```
cd ~/ff_scratch && nohup env PYTHONUNBUFFERED=1 PYTHONPATH=~/ff_scratch ~/miniconda3/envs/ffn_sim/bin/python \
  aleph/scripts/ff_crawl_on_substrate.py --from-resting --cortex-fil 38000 --microtubules --implicit \
  --dt-impl 0.05 --myosin-linear --ecm --ecm-fibers 3000 --steps 2000 --kmc-every 50 --record-every 500 \
  --audit --tag ff_ecm_native --device cuda:0 > ~/ff_scratch/ff_ecm_native.log 2>&1 &
```

## Closeout obligations (every session)
Commit on the feature branch (do NOT push `ffn/foundation` without PI sign-off). `make kb-check` clean before
commit. Update the Notion Dev Logs + Open items. Refresh figures at milestones (native, full, browser-verified).
Final user-facing line: **`Notion 업데이트 완료`**. Respond to the PI in Korean; code/commits/docs stay in the
existing convention. Offer real choices on PI decisions (2–4 alternatives), never one-and-approve.
