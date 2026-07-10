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

---

## Loop log (living)

### Tick 0 (02:42–04:00 KST) — S6 caveat RESOLVED + ECM viz + drag-blocker breakthrough (coarse)
- **S6 two-way (Option B) native (`adb051e`):** collagen-as-substrate → densification **+264.9 nm** (was +15.7),
  traction-driven, STABLE. TAG cross-check OK (per-clutch 17–25 pN = KB 5–20). Caveat resolved.
- **ECM viz (`adb051e`):** render the ECM like the cell — recruit/tension/reorientation scenes. **Direct visual
  check found a real bug** (line layers ignored `color_frames` → fixed) — grep would never have caught it.
- **Longer native (400 s, `c507d2d`):** remodel is PROGRESSIVE (densification +264.9 → **+391.4 nm**, not a 100 s
  artifact); fibre reorientation genuinely slow (RAI Δ +0.001 → +0.002 — hours-scale). Frame cap ≤4–5 at native
  (8 frames = 614 MB HTML loads blank; 4 = 360 MB OK).
- **⭐ Drag-blocker insight (the migration unlock):** the `--com-drag` RUNAWAY was a *dish-context* artifact — there
  the clutch anchors follow the cell, so nothing regularises the rigid COM mode. In **S6 the collagen is PINNED**,
  so bound clutches DO resist COM translation → the rigid mode is regularised. Coarse test confirms: `--com-drag`
  in S6 is **STABLE (no runaway, V/V0=1.000)** and ~10× faster crawl (0.89 → 9.42 nm/s). `com_gamma=6πη·R` is a
  DERIVED physical drag (not tuned). **Native S6 `--com-drag` migration test running** — the decisive question is
  Nc-invariance (does native stay ~9 nm/s, so the cell visibly leaves its position over a long sim?).
- Next: read the native migration result; if stable + Nc-invariant → long migration sim + animated COM-trail viewer
  (cell crawls across + remodels the collagen). If it runs away at native → document + fall back to the rigid-mode
  regularisation (handle rotation modes too, not just translation).

### Tick 1 (04:00 KST) — native migration verdict + NEW-KB cross-check reframes the whole direction
- **Native `--com-drag` migration verdict:** STABLE at native (no runaway, V/V0=1.000 — the collagen-regularisation
  STABILITY hypothesis is CONFIRMED) BUT **still ~0 migration** (v = −0.02 nm/s). com_drag alone does NOT unlock
  native migration. Root cause: the cell is FIRMLY ANCHORED — 10953 clutches grip PINNED collagen, and each clutch
  stays TETHERED to its frame-0 collagen node (`en_base_d` never refreshed), so the cell remodels in place but can't
  translocate. **The blocker is not the drag — it's that clutches don't RE-GRIP new collagen ahead as the cell
  protrudes.**
- **A big, directly-relevant KB-1 ECM series just landed** (parallel session's ECM library). Cross-checked the S6
  remodel field against it (`ff_s6_kb_crosscheck.py`):
  - **KB-1.10** (σ(r)~r^-n; n~1 fibrous force-chains, n=3 continuum): our field decays with **n=2.24** — more
    continuum-like than the fibrous n~1 real collagen shows. Cause: our Mikado is LINEAR — it lacks the
    **strain-stiffening (KB-1.4)** that builds the long-range tensile force chains. The parallel session's
    `ecm_library`/`ecm_mechanics` already has EMERGENT strain-stiffening → **integrate it as the S6 substrate**.
  - **KB-1.9** (nematic order S=<cos2θ>; tumor stroma 0.3–0.7, healthy <0.1): our fibre alignment built only
    **ΔS=+0.005** over 400 s (S 0.033→0.038) — negligible vs tumor stroma. Alignment is a slow, hours-scale +
    reorientation-enabled process (matches the tiny RAI).
- **Re-plan (KB-grounded):**
  1. **Migration = clutch RE-GRIP (KB-1.23: FA captures fibres within R_FA≈1.5µm), not a drag fix.** As the cell
     protrudes, released (rear) clutches rebind the NEAREST CURRENT collagen node ahead → the grip rolls forward →
     the cell walks across the slab (+ com_drag keeps the rigid mode physical). Implement: at the kmc rebind, re-run
     `attach_clutches_to_ecm` against the live `Ep_d` for rebinding clutches and update `en_base_d`.
  2. **Fidelity: use the ECM library's strain-stiffening collagen (KB-1.4)** so stress propagates fibrous n~1
     (KB-1.10) — coordinate with the parallel session's `ecm_mechanics` rather than duplicate.
  3. **MMP proteolytic migration (KB-1.20)** is the deeper cancer-invasion mechanism (degrade fibres ahead, remove
     when L_f<0.5µm; k_deg~1e-3/s) — KB flags it Phase-3+, so stage it AFTER re-grip.
  4. **Longer sims** for alignment to approach the KB-1.9 tumor-stroma S band, measured with the KB nematic order.
- Next action: implement clutch re-grip (migration unlock), coarse-test that the cell now translocates, then native.

### Tick 2 (04:30 KST) — re-grip verdict: the native migration blocker is ISOTROPY, not drag/tethering
- **Native `--ecm-regrip --com-drag`:** STABLE (V/V0=1.000) but STILL ~0 migration (disp∥=+0.007µm, v=0.03 nm/s).
  Re-grip raised traction 187 → **1291 nN** (~157 pN/clutch, ABOVE physiological 5–20 — over-engagement) and
  bound 0.75, yet no net translocation. Visual check: the recruit field is symmetric/ISOTROPIC (pulled equally
  from all sides).
- **⭐ Confirmed diagnosis (two native tests):** native migration is blocked by **ISOTROPIC adhesion**, not the
  drag (com_drag stable) and not tethering (re-grip re-engages). ~10⁴ clutches grip all around the footprint →
  the traction is radially balanced → **net forward force ≈ 0**. The coarse "emergent crawl" (C2, ~13 nm/s) was a
  small-N asymmetry fluctuation that the law of large numbers washes out at native. This matches the biology:
  strongly-adherent EPITHELIAL MCF7 remodels the matrix in place but is poorly migratory. **Migration requires
  explicit POLARISATION (front-back asymmetry), which is a real Rho/Rac program, not tuning.** (The parallel DCM
  session independently reached the same place: it solved cell spreading via ACTIVE SELF-PROPULSION / jamming→
  unjamming — the multicellular analog of single-cell polarisation.)
- **PI decision — how to make the cell leave its position (recommendation first):**
  1. **(RECOMMENDED) Explicit polarisation program:** rear-enriched actomyosin CONTRACTION (pulls the body toward
     the front) + front-biased nascent adhesion (front clutches grip new collagen) + rear de-adhesion — all fixed
     to `phat` (NOT COM-relative, to avoid the rejected `--treadmill` run-down). This is the physiological
     mesenchymal-migration engine; it directly cures the isotropy.
  2. Compliant/sparse collagen (biphasic regime) — lets the matrix deform more, but the REAR stays isotropically
     anchored, so alone it won't translocate; useful as a co-lever with (1).
  3. MMP proteolysis (KB-1.20) — degrade fibres ahead to release the front / rear; cancer-invasion mechanism.
  4. Accept MCF7-as-non-migratory (remodel-in-place) as the honest epithelial result.
- Proceeding with (1): implement a minimal direction-fixed polarisation (`--polarize`), coarse-test for net
  translocation, then native. If it converges → long migration demo + COM-trail viewer.

### Tick 3 (05:00 KST) — migration verdict CONSOLIDATED: MCF7 remodels-in-place, doesn't migrate (robust)
- `--polarize` native: coherence 0.03→**0.15** (directionality emerged) but v=0.09 nm/s (no translocation).
- `--n-fa 200` (physiological FA count, the anchor fix + the 2× traction correction) native: coherence **0.29**
  (highest) but v=0.07 nm/s (still ~0). traction 0.29 nN.
- **ROBUST CONCLUSION (5 native tests):** across emergent-crawl, com_drag, re-grip, polarize, and n-fa, the
  directionality steadily improved (coherence 0.05→0.15→0.29) but **translocation never emerged** (v≈0.07–0.09
  nm/s throughout). A native MCF7 cell on stiff PINNED collagen-I **remodels the matrix in place but does NOT
  migrate.** This is physiologically CORRECT (epithelial MCF7 is poorly migratory) and was NOT tuned for — the
  model reproduces the remodel-not-migrate phenotype. The distributed grip on the stiff pinned matrix anchors the
  cell; front-back asymmetry alone is too weak to translocate. The S6 CORE deliverable (matrix-remodel caveat
  RESOLVED, +264.9 nm coherent inward) STANDS and is unaffected.
- **PI DECISION POINT (migration is a research direction, not a tweak — needs the PI's greenlight):**
  1. **(RECOMMENDED) MMP proteolytic invasion (KB-1.20):** the cell secretes MMP at the leading edge, degrades
     collagen ahead (dL_f/dt=−k_deg·ρ_MMP·L_f, remove at L_f<0.5µm), opening a channel it grips + pulls into.
     This is THE cancer-invasion mechanism and directly removes the stiff-matrix anchor. A substantial new
     mechanism (secretion + diffusion + fibre degradation + cyclic protrusion) — a proper implementation session.
  2. Compliant/unpinned collagen (biphasic regime) — a softer, yielding matrix the cell can deform through.
  3. Model a MOTILE cell type (mesenchymal / EMT-MCF7) with a single dominant leading edge, not epithelial MCF7.
  4. Accept remodel-in-place as the honest MCF7 result (it matches the biology).
- **Proceeding meanwhile (certain, KB-grounded deliverables):** the biphasic clutch cross-check (Chan-Odde /
  Bangasser, KB-1.5 stiffness, `--ecm-lp-um` sweep) and the `ecm_library` physiological-collagen integration
  (nematic-S alignment KB-1.9). These produce clean results while the migration-mechanism decision awaits the PI.
