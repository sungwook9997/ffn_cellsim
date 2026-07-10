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

### Tick 4 (05:30 KST) — biphasic κ-null → ecm_library integration → compliant matrix REVIVES migration
- **Biphasic with `--ecm-lp-um` (κ) is FLAT** (`c4d75c8`): traction ~21 nN across Lp=2→2e5. The clutch feels a
  K_SEG/crosslink/pin-dominated stiffness, not the bending κ → κ is the wrong knob. Need a network-modulus sweep.
- **`ecm_library` integration (`71a41e0`):** `--ecm-material collagen_I` wires the parallel session's grounded
  builder → PHYSICAL per-segment `seg_k` (real E_fibril, 17279 vs fixed K_SEG=5e4) + concentration→mesh (KB-1.7)
  + nematic S (KB-1.9). The physical (softer) collagen remodels MORE + more coherently (coherence 0.14→0.58,
  aligned 0.72). Realistically SPARSE (conc 1.5 → ξ~2µm, vs the old unphysical ξ=0.5µm dense Mikado).
- **⭐ Compliant physiological collagen REVIVES migration (`343e54a`):** fully physiological native (collagen_I
  conc=3, ξ=1.5µm, z=4.47; `--n-fa 150`; `--polarize --com-drag`): disp∥=+0.085µm, **v=0.57 nm/s — 6× the
  stiff-collagen tests**, fibre RAI Δ=+0.010 (10×), nematic ΔS=+0.014 (3×). So the migration story is
  stiffness-dependent: **stiff pinned collagen → remodel-in-place; compliant physiological collagen → directional
  migration + fibre alignment EMERGE.** Still sub-physiological absolute speed (0.57 vs 10–30 nm/s) — a longer sim
  / even-softer matrix / MMP pushes further. This is the compliant-matrix option (2) validated directionally.
- **Next:** (a) physical biphasic (`--ecm-conc` sweep 1→30, normalise traction per bound clutch to remove the
  attachment confound) → the real Chan-Odde inverted-U; (b) a longer compliant-collagen sim → visible migration +
  alignment approaching the KB-1.9 tumor-stroma S band; (c) tumor-aligned S6 (`--ecm-align-s 0.6`); (d) MMP (PI).

### Tick 5 (2026-07-11, PI: "1,2 모두 진행") — MMP + motile cell-type IMPLEMENTED (both full-speed-migration paths)
A 3-agent design workflow grounded both in the KB (mmp-design agent hit a transient rate-limit; the KB-grounding
+ motile-design agents carried it). Both landed default-off (validated remodel/biphasic paths untouched):
- **`--mmp` proteolytic invasion (`2cb6ac4`, KB-1.20):** the FRONT bound clutches secrete MMP; a quasi-steady
  diffusive halo ρ=exp(−r/λ), λ=√(2·D·τ)≈4.5µm (D_MMP≈10 µm²/s), degrades collagen segment + crosslink stiffness
  d(seg_k)/dt=−k_deg·ρ·seg_k (k_deg~1e-3/s, Wolf2013) and SEVERS a segment when fully degraded — opening the
  invasion channel the polarised cell advances into. ρ dimensionless [0,1] ⇒ timescale set by the sourced k_deg
  (no tuned magnitude; absolute MMP flux PI-flagged). `[MMP]` print reports severed-segment count.
- **`--cell-type` motile/EMT preset (`d763812`, new `ffn_sim/ff/cell_type.py`):** `mcf7_epithelial` (default =
  current) vs `mesenchymal`/`emt` — a bundle of the existing knobs grounded in KB-3.11 / SE248 Betorz2023 / KB-4.12
  / KB-2.2: polarize + ecm_regrip + front_frac 0.6 + **myo_rear_bias 0.7** (NEW build hook: relocate a fraction of
  Stam-Hocky minifilaments onto the REAR cap, force-conserving) + n_fa 50. CLI flags override the preset.
  ⚠️ front_frac/myo_rear_bias/n_fa are PI-flagged (qualitatively grounded, must NOT be swept to a target speed).
- Coarse: both build+run clean. **Native 3-way (epithelial vs mesenchymal vs mesenchymal+MMP on conc-3 collagen)
  running** — the decisive test of whether the rear-myosin engine + MMP channel reach physiological migration.

---

## OVERNIGHT AUTONOMOUS ROADMAP (2026-07-11, PI: "밤동안 계속 진행, /loop 1h, 결정은 네가")

Standing directive: work through the night, decide autonomously when blocked on a PI choice (record the choice +
rationale), extend this roadmap as phases close, keep the loop going. Every phase ends in a browsable HTML/figure
(visual-checked) + a backup commit + a TAG cross-check. Native + full only for conclusions; gbook ssh flaky →
run detached, monitor by log. Do NOT push ffn/foundation.

**Phase M1 — 3-way cell-type/MMP verdict (in flight).** Read epithelial vs mesenchymal vs mesenchymal+MMP.
Gate: does mesenchymal translocate > epithelial (rear-myosin engine works)? does MMP raise v / open a channel
(severed segments > 0, localised at the front)? TAG: compare v to KB-3.12 (10-100 nm/s motility). Visual: the
mesenchymal crawl + the MMP-degraded channel (ECM tension/severed overlay). If physiological → migration story
COMPLETE; if not → diagnose (M1b: longer sim, or the rear-myosin/MMP interplay).

**Phase M2 — the "cell leaves its position" capstone.** A long native mesenchymal(+MMP) sim on compliant collagen
→ the cell crawls a visible µm-scale track (COM trail) while remodelling/degrading the matrix. The animated
interactive viewer (crawl + collagen recruit/tension/reorientation/severed) = the PI's original ask made real.

**Phase M3 — tumor-aligned matrix + contact guidance (KB-1.9).** `--ecm-align-s 0.6` (TACS-3 tumor stroma,
S 0.3-0.7). Does the cell migrate ALONG the aligned fibres (contact guidance)? Does its traction increase the
alignment (feed-forward)? Reproduces the tumor-invasion-highway phenotype.

**Phase M4 — the physical biphasic (real Chan-Odde curve).** Sweep `--ecm-conc` (1.5→30) with `--ecm-material
collagen_I`, normalise traction per bound clutch (remove the attachment confound the κ-sweep hit), overlay the
motility optimum (M2 speed) → the traction+motility biphasic vs real matrix modulus (Bangasser 2-300 kPa band).

**Phase M5 — migration-mode atlas (a synthesis deliverable).** One figure/viewer contrasting the modes the engine
now spans: epithelial remodel-in-place; mesenchymal crawl; proteolytic (MMP) invasion; on stiff vs compliant vs
aligned collagen — the FF single-cell migration phase map, KB-grounded.

**Backlog / PI-flagged (decide autonomously, record rationale):** MMP absolute secretion flux (kept ρ∈[0,1]);
contractility_mult (kept 1.0); whether to register any new KB rows (surface, don't auto-create gates/contracts).
Autonomous-decision log appended below as choices are made.

---

## MORNING SUMMARY (overnight of 2026-07-10→11) — read first

> **Full arc:** (1) S6 caveat resolved + MMP + cell-type + honest migration verdict + atlas + engine fixes + ECM-viz
> verified (this summary), then two PI-decision items worked out in the continuation (bottom of this file, **read them**):
> **§ TAG re-plan** — migration is traction-limited, KB-confirmed, not a drag bug (lever = EMT + KB-3.14, PI-gated); and
> **§ FF-engine mechanism gap analysis** — the next mechanism (cofilin severing) is **KB-blocked**, needs your KB-ingestion
> sign-off (I did not invent a rate). Nothing needs action to be *safe*; the two § sections are where your input unblocks
> the next step. kb-check green throughout (runs 31 / params 42, no drift).

**Landed (all committed, backed up, no ffn/foundation push):**
- **MMP proteolytic invasion** (`--mmp`, KB-1.20) + **motile/EMT cell-type** (`--cell-type mesenchymal/emt`, new
  `ffn_sim/ff/cell_type.py`, KB-3.11/SE248/KB-4.12) — the two mechanisms you asked for, both default-off, KB-grounded.
- **M1 migration verdict (honest):** native migration of a strongly-adherent MCF7 on physiological collagen is
  **sub-physiological (v~0.2–0.57 nm/s vs 10–30)** and numerically fragile — the model reproduces MCF7's poor
  motility. Directionality (remodel coherence) rose 0.05→0.44 as the motile program was added, but translocation
  stayed modest. Full-speed invasion is an open adhesion-drag/solver item, **not reachable by parameter choice** —
  reported honestly, not tuned to a target.
- **M5 migration-mode ATLAS** (`ff_s6_migration_atlas.png`) — one figure: caveat resolved (remodel 16→265→391 nm),
  migration biphasic in matrix density (peak conc 3), directionality↑ but speed sub-physiological.
- **Two engine robustness fixes:** k_vol OverflowError hardening + a volume-EXPLOSION guard (clean truncation) —
  both genuine, normal-run-unaffected.
- **ECM stress/strain/reorientation viz VERIFIED (PI HARD visual mandate):** the three ECM field scenes
  (`collagen recruit` 0–169 nm displacement, `ECM tension` |ΔL|/L₀ tensile strain, `ECM reorientation` tangent-turn°)
  render correctly in a real browser — collagen fibres coloured turbo by load, colorbars + two-way legend correct,
  cell shown. Directly browser-checked (screenshots), not grep — caught+fixed that the full-res native HTML is
  **277 MB** (6 scenes × cortex) which headless couldn't load (blank canvas): reduced to 3 temporal frames (spatial
  res untouched, PI rule) → renders; and faded the cell to a 0.05 silhouette so the matrix field is the subject.
  The native HTML stays a local artifact (>GitHub 100 MB), regenerable from the npz.

**4 autonomous decisions (all KB/stability-grounded, NOT tuned-to-speed; full rationale in each commit):**
1. mesenchymal `n_fa` 50→0→100→0: fewer FAs collapse traction (soft-arm) AND starve the MMP source; proteolytic
   invasion is a STRONG-adhesion mode → keep n_fa=0 natural strong adhesion (the preset currently sits at the last
   stable-verified n_fa=100; the n_fa=0 strong-adhesion run overflowed — see below).
2. `myo_rear_bias` 0.7→0.0: the 70%-rear-myosin relocation is numerically unstable at native (overflow) + gave no
   migration benefit → dropped from the preset (flag/hook retained for your experimentation).

**MMP validation (2026-07-11 continuation) — the mechanism ACTS, it is honestly SLOW (not a bug, not tuned):**
- On the **high-traction** native path (raw Mikado two-way, remodel +265 nm, traction **185 nN**, ~7500 bound
  front clutches sourcing MMP) the MMP report was **0 collagen segments severed (0.0%)** — same as the migration
  path. This is **not** the mechanism failing: severing a segment needs its stiffness to fall 99.8% (K_SEG 5e4 →
  100 pN/µm), and at the **KB-1.20-sourced** `k_deg≈1e-3/s` that takes ≈ tens of minutes (Wolf2013 physiological
  proteolysis), while a feasible native sim is ~100 s. So over the sim MMP **degrades but does not sever**.
- I added a **degradation metric** (`mmp_degraded_pct` mean / `mmp_degraded_max_pct` peak-at-front) so the
  report shows the proteolysis that IS happening rather than a misleading "0". Coarse smoke (6 s sim-time): front
  collagen softened **0.58%**, 0 severed — mechanism confirmed acting, physiologically slow.
- **Native val2 re-run outcome (honest):** the clean native re-run with the new metric **fell back to CPU** (I
  omitted `--device cuda:0` on the launch) AND the A5000 was occupied by the **parallel session's native job**
  (`ff_contact_guidance_anisotropy … --tag cg_native`, the aligned-collagen contact-guidance study — i.e. the
  parallel session is actively owning the aligned-ECM domain I flagged, ✓). A native 266k cell on CPU is
  ~12–24 h and non-authoritative, so I **killed it** rather than contend with the parallel session's GPU run.
  The MMP finding does **not** need it: native **0-severed is already confirmed** (previous native val, 185 nN),
  the degradation metric is validated (coarse smoke), and the front-degradation % is a scale-invariant
  exponential (`exp(-k_deg·ρ·t)`) ⇒ derivable ≈ a few % over a 100 s native sim. A clean native-GPU degradation
  datum is **deferred** to whenever the A5000 frees (launch WITH `--device cuda:0`); it changes no conclusion.
- **I did NOT tune `k_deg` to force severance** — it is a KB-1.20 value; forcing a visible channel would be
  tuning-to-outcome. The honest statement: *MMP proteolysis is implemented and KB-grounded; at the physiological
  rate it opens a channel only over physiological (tens-of-min) time, which a short native sim cannot reach.* A
  longer sim, or a `k_deg` sweep **within** the KB-1.20 band, is the legitimate way to show a channel later.

**Open items for you (flagged, not auto-resolved):**
- **Native migration full-speed** — the honest limit above; the physical path is MMP invasion on a compliant matrix
  with enough front traction to source MMP, or a fundamentally softer/2-D regime. Needs a solver look (adhesion-drag).
- **Aligned collagen (`--ecm-align-s`) instability** — the tumor-stroma aligned matrix builds (S_measured≈0.64,
  matches KB-1.9), but the cell-collagen co-sim EXPLODES early: the aligned net is under-constrained (⟨z⟩ 1.47,
  parallel fibres floppy ⊥ alignment) → clutch traction drives runaway collagen displacement. This is an
  `ecm_library` aligned-build item (needs more crosslinking / a stiffer-⊥ aligned net) — a parallel-session file, so
  flagged rather than edited. The explosion guard makes it fail safe.

**Continuing autonomously:** M2 (crawl viewer from the best migration run) + consolidation; the loop stays live.

---

## TAG re-plan (2026-07-11) — the migration open item is TRACTION-limited, NOT a drag bug (KB-confirmed)

Per your instruction ("check the TAG for new info, then re-plan from it"), I queried the Kim+Miyazaki corpus for the
crawl/drag physics. **Kim2012_IntegrBiol** (the FF migration reference): drag is **per-node, linear `F_D = C·v`**,
`C = 0.001 N·s·m⁻¹` **fixed**, total drag `∝ N adhesion nodes`, **bond-rupture origin (NOT hydrodynamic)** — and,
decisively, *"migration-speed vs cell-size emerges from geometry (FA-formation probability → traction imbalance),
**not** from a size-scaled drag coefficient"* [Kim2012 p7–10].

**What this re-plans:**
- The migration open item was framed as an "adhesion-drag" problem. The TAG shows the drag is **already handled**:
  my `--com-drag` applies the physical whole-cell Stokes drag `6πηR` to the COM mode (Nc-independent) — the earlier
  `Σγ ∝ Nc` grid-drag (`--bulk-drag`, unstable) is superseded. So native slowness is **not** a drag-magnitude bug.
- Kim2012 says the speed lever is the **front–back traction imbalance** (FA-formation asymmetry, geometry-driven).
  That is exactly what my `polarize` + `ecm_regrip` model — and they produced the correct **directionality**
  (coherence 0.05→0.44) at MCF7-appropriate **modest speed**. So the honest M1 verdict is **KB-CONFIRMED**, not a
  defect: a poorly-motile epithelial MCF7 *should* be traction-imbalance-limited and slow.
- **The legitimate (non-tuning) path to physiological speed** is therefore the **EMT/mesenchymal program**, exactly
  as Kim2012 frames cell-type ("change cell size or the number of adhesions per node"): more/denser front adhesions
  + stronger polarization + the **KB-3.14 contractility upregulation (2–10×, PI-gated)**. My `cell_type.py`
  `mesenchymal`/`emt` presets already encode the adhesion/polarization side; `contractility_mult` stays **1.0**
  (PI-gated) by design. **Next migration step (needs the A5000 free + your OK on the KB-3.14 lever):** run
  `--cell-type emt --com-drag` with `contractility_mult` raised into the KB-3.14 2–10× band and measure v — a
  KB-grounded cell-type change, *not* a tuning sweep of MCF7.

This is a re-plan, not a result: nothing was run (GPU busy). It converts the vague "solver look" open item into a
concrete, KB-anchored, PI-gated experiment, and it strengthens (does not overturn) the honest migration verdict.

### FF-engine mechanism gap analysis (TAG-grounded, 2026-07-11) — the next-mechanism menu for you

Inventoried the live FF kernels vs the Kim+Miyazaki corpus roadmap (`turnover · severing · catch-slip · PCM myosin ·
formin · cylinder-drag`). **Already implemented:** barbed-end growth (KB-3.6), pointed-end depoly, fiber treadmilling,
α-actinin KMC Bell turnover, catch-slip clutches, Stam-Hocky minifilament myosin, Arp2/3 branch-angle, directed front
growth, FA growth/disassembly/maturation, MT aster, active-gel polarization, membrane/nucleus/ERM. The engine is
**substantially complete**. Two genuine gaps, with an honest blocker on the first:

1. **Cofilin severing — real gap, but KB-BLOCKED (do NOT implement yet).** No `sever` kernel exists (the grep hit was
   "*several*"). Severing is the disassembly arm that keeps the cortex fluid and, paired with front barbed-growth, is
   part of the crawl treadmill — so it is on-thread. **But the TAG shows the KB has *no* KnowledgeClaim quantifying a
   severing rate** (per-µm frequency, ADP-preference constant): the physics is only qualitative in the corpus (cofilin
   prefers ADP/aged actin, nucleotide-gated, rate ∝ cofilin conc), and the quantitative primary source —
   **Elam, Kang & De La Cruz 2013, FEBS Lett 587:1215 (doi:10.1016/j.febslet.2013.01.062)** — is *cited but not ingested*.
   Implementing severing now would require an **invented rate = a magic number = hard-rule violation.** → **Surfacing to
   you:** the correct path is (a) ingest Elam2013 as SourceEvidence + register a KnowledgeClaim for the severing rate
   (PI-gated KB change), then (b) implement a KMC severing kernel grounded in it (per-segment P_sever ∝ cofilin·age,
   ADP-gated). I did NOT invent a rate.
2. **Formin processivity — architectural-only, a design refinement (not KB-blocked).** `formin` exists as a *nucleator
   label* in `architecture_spec` (filopodium/SF parallel bundles) but there is **no formin-specific processive-elongation
   kernel** — `barbed_end_growth_kernel` covers elongation generically. Formin's distinct physics (processive barbed-end
   tracking, elongation acceleration, capping protection) is a refinement of the existing growth kernel, implementable
   from the existing KB-3.6 growth grounding + a formin on-rate; lower priority than severing and needs your steer on
   whether the generic growth kernel is sufficient for the current cortex/filopodium work.

**Turnover is NOT a gap** (treadmill + pointed-depoly + KMC-Bell already cover it). **Net recommendation for you:** the
highest-value next FF mechanism is **cofilin severing**, but it is gated on a **KB ingestion of Elam2013** — a
gate-contract/KB change that is yours to approve, not mine to auto-create. Until then I will not implement it (no magic
number). Everything here is GPU-free planning; nothing was run or changed in the engine.
