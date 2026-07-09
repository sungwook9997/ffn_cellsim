# FF S6 → migration-on-ECM — 12 h autonomous plan (2026-07-10)

**PI mandate (2026-07-10, 02:42 KST):** run autonomously ~12 h developing the FF engine in the current
direction. Goal: **the MCF7 cell moves well *on the collagen-I ECM* and visibly LEAVES its starting
position** over a (long) simulation, with the **ECM rendered like the cell** — where it is stressed, where
it deforms, where it *reorients* — in the `ffn_sim/outputs/ff/**/figs` interactive HTML the PI browses.
Check every 30 min that the work is progressing and the sim is not diverging; back up (commit) so nothing
is lost; from ~4 h on, check whether the TAG/KB refresh has completed and, if so, re-plan on the new info.

This plan is **living** — appended each loop tick with what landed and what's next.

---

## Where we are (start of the mandate)

- **S6 two-way (Option B) landed** (`606713a`): when `--ecm`, the collagen-I Mikado network **is** the
  substrate — the basal clutches pull the live fibres (two-way, cell feels the matrix), replacing the rigid
  dish pin that had absorbed the traction (the root cause of the "+15.7 nm-of-266 nm follows-not-remodels"
  caveat). Coherent-remodel metrics (recruit / densification / settling / coherence + fibre radial-alignment
  index) + biphasic knob `--ecm-lp-um` + `--audit` traction-driven verdict. Validated crawl path untouched.
- **Native two-way run in flight** on the A5000 (`ff_ecm_s6_twoway`, 266 k cell + 3000-fibre collagen,
  10953/10953 clutches, from-resting physiological baseline). CPU smoke (non-authoritative) already shows
  recruit +70–102 nm and OFF-audit → 0 (traction-driven).

## The one hard blocker between us and "the cell leaves its position at native scale"

**Native crawl SPEED is grid-drag-limited** (`FF_CRAWL_DIAGNOSIS_2026-07-09`, memory
`project-ff-crawl-mechanism`): `physical_node_gammas` gives Σγ ∝ Nc (per-node single-fibre log drag), so the
whole-cell COM drag is ~100–300× the physical Stokes 6πηR ⇒ **v ∝ 1/Nc** (native ~0.12 nm/s, ~500× too slow).
The correct Σγ = 6πηR (`--bulk-drag`) makes the per-node γ so small that γ/dt ≪ K → the cortex **rigid-body
modes go unregularised → NaN**. Two prior fixes failed: `--bulk-drag` (uniform-small-γ NaN) and `--com-drag`
(rank-3 COM correction → **runaway** 767 µm/30 s, despite the arg-help calling it "working" — the memory
diagnosis supersedes that label). **This is a PI-level numerics item: the whole-cell translation drag must be
made physical INSIDE the implicit solve with the rigid-translation mode regularised.** Solving it is the
central technical task of this mandate — without it, native migration is either artifactually frozen or blows up.

---

## Phases (each ends in a browsable HTML + a backup commit)

### Phase 1 — result + make the ECM visible (0–2 h)
1. **Native two-way remodel result.** Read `[ECM-REMODEL]`/`[ECM-ALIGN]`; compare coherent recruit +
   densification + RAI Δ vs the +15.7 nm baseline. If coherent inward ≫ +15 nm with coherence up and RAI Δ>0 →
   the caveat is resolved by two-way coupling alone (no new contraction mechanism = no tuning). If still flat →
   record it as the lever-2 fork for the PI (basal-inward contraction: body-force kernel vs build-side basal
   myosin enrichment) — do NOT silently add a tuned force.
2. **ECM stress / strain / reorientation viewer (PI's explicit ask).** Extend the collagen layer in
   `ff_crawl_viewer.py` (and/or a dedicated `ff_ecm_field_viewer.py`) so the fibres are coloured, like the
   cortex σ_vm/strain, by: (a) **fibre tension/stress** (segment + crosslink strain energy → per-node virial),
   (b) **areal/segment strain** (‖posf−pos0‖ already have; add local stretch), (c) **reorientation** (Δangle of
   each fibre tangent frame0→final, and the radial-alignment field). Save the per-frame ECM stress/strain into
   the npz from `run()` so the viewer can animate the fields, not just displacement. Toggleable scenes; full-res
   (no downsampling, PI rule). `browser_check.py` to confirm it renders.
3. **Longer remodel sim** so the matrix reorganisation is visible over time (steps ↑, record more frames).

### Phase 2 — migration on the compliant ECM (2–5 h)
4. **Does the cell migrate differently now collagen is the substrate?** The traction no longer sinks into a
   rigid pin; run the crawl (`protrude` + clutch turnover) ON the collagen and measure COM drift + whether the
   drag picture changed. Cross-seed.
5. **Attack the native drag blocker (the core task).** Implement rigid-translation-mode regularisation in the
   implicit solve so Σγ_translation = 6πηR while the deformation modes keep their well-conditioned per-node γ —
   the correct, stable version of what `--com-drag` botched. Validate: (a) no NaN/runaway, (b) v scales with the
   physical drag not 1/Nc, (c) a free (no-clutch) cell under a body force translates at 6πηR-consistent speed,
   (d) native stability (V/V0≈1) preserved. This is the unlock for native migration.
6. **Long native migration run** — the cell crawls across the collagen slab and LEAVES its start (COM trace),
   remodelling/aligning fibres in its wake. Animated HTML with the COM trajectory + the ECM fields.

### Phase 3 — TAG re-plan + keep improving (4 h onward, interleaved)
7. From ~4 h, poll whether the TAG/KB refresh finished (new SourceEvidence/claims). When complete, query the
   new info (`tag_query.py`) for MCF7-on-collagen migration/traction/remodel data and **re-plan** improvements
   on it (append here). Until then, keep executing Phase 1–2.

## Operating rules for the autonomous loop
- **Native + full only** for any conclusion; coarse/CPU = labelled non-authoritative smoke tests.
- **No tuning-to-pass / no magic numbers / no gate-loosening.** New mechanisms at physiological values; a real
  design fork → record as a PI choice with alternatives, proceed with the recommended one (PI: "내 권고
  기다리면 너 권장대로 진행").
- **Every ~30 min:** confirm the running sim is advancing and not diverging (log tail + etime + no NaN), advance
  the plan, and **commit a backup** (nothing lost). ScheduleWakeup ~1800 s cadence.
- **Visualisation is the deliverable surface** — every milestone lands as browsable HTML in `ff/**/figs`,
  browser-checked, ECM shown like the cell.
