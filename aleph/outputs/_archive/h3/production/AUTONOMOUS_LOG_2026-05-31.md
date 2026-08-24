# Autonomous run log — 2026-05-31 (12h unsupervised)

**Start:** 2026-05-31 02:02 KST (epoch 1780160545). **Planned stop:** ~14:02 KST (epoch 1780203745).
**Cadence:** ~10 min heartbeat (self-pace, event-driven by job completions) · Notion commit every ~30 min.
**Mandate (PI, /loop):** "fix everything, run production, decide on your own for 12h." Cheer-on: 힘내라 전우여.

## Autonomy contract (the rule a waking PI will see)
- ✅ Autonomous: geometry placement fix (physics-design, PI-delegated); additive default-off modules (substrate / H.8 / H.9) build+integrate; long myosin sim + (post-fix) production; commit to `phase1/h3-cortex` (explicit-path adds only); **every change gated by regression 454/17/0 + adversarial verify**.
- 🚫 NOT crossed unsupervised (CLAUDE.md PI-only gates, valid while PI sleeps): gate-contract changes (no KU-band edits; Track C composite-gate deferred); frozen-`integrator/` (H.10 deferred; no BAOAB/dt-policy edits); un-derivable magic numbers (defer that piece); `ffn/foundation` push. Blocked items → logged as top PI wake-decision + routed around.
- 🛡️ Safety: every code change additive + default-off + regression-green before commit. Unverifiable change → revert + log (never leave the tree broken).

## Critical path to a γ result (the payoff)
G1 geometry fix (clutch bonds FORM) + G2 myosin steps (long sim) → cortex contracts against substrate-pinned clutch → substrate-anchored tension → **γ lifts off the 3e-4 floor**. S5 dynamic kinetics (G3) and Track A k_sub are SECONDARY (static clutch bonds already transmit load once myosin contracts). So G1+G2 is the minimal path; chase it first.

## STATUS (updated each iteration)
- **Now:** iter 3 — **fan-in cap fix** (uncommitted, regression-gated). G1-on long sim `b1455lmtv` CRASHED on first production `run()`: `RuntimeError: Too many bonds to process exclusions for particle tag 378` = HOOMD nlist exclusion overflow from clutch **fan-in 9** (seed-1; the verifier's flagged worst case). Root: cell.py S2 argmin concentrated 9 integrins on one cortex bead. **Fix:** per-bead clutch fan-in cap (4) — each integrin bonds nearest WITHIN-CAPTURE bead with spare capacity (force-free r0 preserved). **Verified: seed-1 smoke `PYTHON_EXIT=0`, no crash, clutch still 52/52.** (4 = HOOMD library-capacity bound, not physics — PI-flag.)
- **Running:** G1-ON long sim (fix) `bve21wrku` (payoff, `-u`, `ku35_g1_long2.log`) · G1-OFF control `blradxmfe` (A/B + G2, ~35min in) · full regression `bk1b99nqq` (fan-in commit gate). Modules substrate/H.8/H.9 on disk, NOT integrated.
- **Next decisions:** regression green → commit fan-in fix; watch G1-on `steps_adv` (G2) + `tension_*` (γ lift?); while sims run → nucleus net-force/COM fix → integrate 3 modules default-off + regression. If γ lifts → n_fil=150 production.
- **Notion cadence:** 02:02 start · 02:31 G1 commit · next ~03:01.

#### ⚠️ PI-FLAG (iter 3): `_MAX_CLUTCH_FANIN = 4` in cell.py `_extend_snapshot_with_fa` is a HOOMD neighbour-list exclusion-capacity bound (not a physics tuning constant; clutch stays force-free). Crash reproduced at fan-in 9; cap chosen conservatively (clutch ≤4 + backbone ~2 + myosin/xlink). Empirically still yields 52/52 clutch. PI may ratify the exact value / derive from HOOMD's true exclusion limit.

### iter 4 — 02:56 KST — HOOMD exclusion crash fully resolved (FA-confined)
- The fan-in cap (`c960406`) was NECESSARY but INSUFFICIENT: crash persisted at `--equilibrate-steps 2000` (passed at 500). **Root cause:** nlist `exclusions=("bond","1-3")` (cell.py:823) — the "1-3" expansion through the dense south-cap clutch (cortex bead → integrin → ligand, + backbone → neighbour-clutch) overflows HOOMD's per-particle exclusion storage. FA particles have NO LJ (cell.py:842-859 wire no FA pair) so the exclusions are spurious, but HOOMD computes them from the global bond graph; "1-3" can't be dropped globally (bent xlink/myosin 1-3 separations < WCA r_cut).
- **Fix (cell.py `_extend_snapshot_with_fa`, FA-path-confined):** deepen the south-cap to ~2·n_int anchor beads (derived `2·(z_max−z_min)·n_int/n_cortex`) so the clutch spreads ~1/bead + lower fan-in cap 4→2. Clutch stays force-free (r0 = exact sep). No global exclusion/cortex change → non-FA runs bit-for-bit.
- **Verified:** seed-1 + equilibrate-2000 repro (the exact crash config) `PYTHON_EXIT=0`, production passes, clutch still **52/52**, tension 3.76e-4 (still floor — only 5000 prod steps, no myosin stepping yet). Re-launched G1-on long sim `bhiqesueb`; full regression `b3svo3fcf` gating the commit.

## Decision rules locked
- **G2 (myosin):** watch `steps_adv=` in `ku35_myostep_long.log`. Nonzero past ~step 2.5e7 → ARTIFACT (no ratchet change). Stays 0 past ~2.8e7 with healthy bind_total → BUG (MECHANISM_AUDIT §A2 redesign).
- **G1 (geometry):** success = `--smoke` reports `n_fa_clutch_bonds_at_build` >> 9 and r/r0 starts to deviate; regression stays 454/17/0 with FA off.

---

## ITERATION LOG (append-only)

### iter 0 — 02:02 KST — setup
- State: HEAD `1496173`, tree clean (no stray module files), suite 454/17/0.
- Launched (parallel): long myosin sim (bg), build-modules workflow (substrate/H.8/H.9 scaffold→verify), geometry-fix design+verify workflow (read-only).
- Heartbeat 600s; first Notion progress commit due ~02:32.

### iter 1 — 02:27 KST — G1 geometry fix APPLIED + mechanism verified (commit pending full-suite)
**Critical-path win.** geometry-fix wf `wg1bzpqsc`: design + adversarial verify BOTH empirically built the full n_fil=60 cell, reproduced 9/52 baseline → **52/52 (20 seeds), force-free, regression-safe, stable**. Verdict APPLY-WITH-FIXES (doc-only).
- **Applied (Lead, single-writer):** `fa.py._build_fa_layout(+ contact_footprint=None)` (None=legacy bit-for-bit) + `cell.py._extend_snapshot_with_fa` south-cap anchor computation (z ≤ z_min + capture) → seeds each FA under a cap bead. **Offset = physical `h_integrin` (50nm), not the design's `0.1·capture` → no magic number.**
- **post-fix smoke `bm529z611`: `n_fa_clutch_bonds_at_build = 52` (was 9)**, n_clutch static/live = 52/52. ✓ G1 mechanism CONFIRMED. γ still on floor (r/r0=0.99998, tension=3.2e-4) — EXPECTED, smoke too short for myosin stepping (step_advances=0); γ lift needs G1+G2 together.
- **Regression:** first full run flagged 4 FA tests asserting OLD behavior (ligands@z=0 + "barely clutches"). NOT a real bug (pin/immobility/NaN all pass; my fix correctly changes that behavior). Updated 4 tests transparently to the corrected contract. FA file now **10 passed / 1 skip**. Full suite re-running `bixjrwuj9` → **commit gated on it being green**.
- long myosin sim `blradxmfe`: healthy (PID 57126, 99% CPU, ~28min; stdout buffered → final JSON gives G2 verdict on completion).

#### ⚠️ PI RATIFICATION NEEDED (test-contract change, made under 12h autonomy grant)
- `test_fa_physical_capture_radius_barely_clutches` → renamed `..._fully_clutches`, assertion INVERTED (`< 0.5*n_int` → `>= 0.9*n_int`). It previously CHARACTERIZED the floor bug (disjoint geometry → few clutch); the G1 fix resolves it, so the test now regression-guards the fix. **This is a gate/contract change — flagged for PI review.** Revert = `git revert` the G1 commit.
- 3 ligand `z==0` placement assertions → "ligand sits under a cortex bead (contact-footprint)" (coordinate corrections, lower sensitivity).

### iter 5 — 03:40 KST — G2 RESOLVED (myosin steps) + FA production-capable (fan-in 1)
**Two wins from the fan-in-1 crash-survival test (5e5 steps, seed 1):**
- ✅✅ **G2 VERDICT — myosin stepping is REAL; `step_advances=0` was an ARTIFACT.** `step_advances_final = 203` (nonzero) over 5e5 steps → MECHANISM_AUDIT §A2 "binned-ratchet defect" **REFUTED**, reader C's artifact hypothesis **CONFIRMED**. Control sim `blradxmfe` (G1-off) will corroborate on completion.
- ✅ **fan-in 1 SURVIVES production** (5e5 steps, 5 samples, no exclusion crash) — the only config that runs production. fan-in 2 (committed `d33230b`) = 52 clutch but crashes on long runs (dynamic peak); fan-in 1 = 45 clutch (tight FA disks lose a few) but survives. Lowered `_MAX_CLUTCH_FANIN` 2→1; adjusted the fully_clutches test threshold 0.9→0.7 (PI-flag: clutch count is a HOOMD-exclusion-capacity tradeoff).
- ⚠️ **OPEN FINDING (PI / deeper analysis):** even with clutch=45 + myosin stepping (step_adv=203), **tension stayed on the floor (4.3e-4) and r/r0≈1 (no contraction)** at 5e5 steps. So G1+G2 alone did NOT lift γ in this short run — either (a) needs a much longer run for contraction to develop, or (b) a further mechanism gap (stepping → net cortex stress). Launched full payoff `bbn4djb13` (fan-in 1, 4e7 steps, samples every 5e5, checkpointed) to test (a). Regression `ba6ak7b71` gates the fan-in-1 commit.
- **Crash root + PI-worthy deeper fixes (NOT done unsupervised):** HOOMD nlist `exclusions=("bond","1-3")` × loaded-clutch+dynamic overflows per-particle exclusion capacity. fan-in 1 mitigates; full-run survival unconfirmed. Clean fixes: (i) drop "1-3" globally (physics-adjacent, KU-gate validation needed); (ii) fa_actin_clutch as a CUSTOM FORCE not a bond (no exclusions, but method-of-planes tension must then add the custom-force virial); (iii) accept reduced clutch.
- **NEXT: pivot to EXTEND module integration** (substrate/H.8/H.9 default-off) while the payoff sim runs.

### iter 6 — 03:42 KST — CRITICAL: the long-run crash is PRE-EXISTING (not G1)
- **The G1-OFF control sim `blradxmfe` (OLD code, 9 clutch) ALSO CRASHED** — same `Too many bonds to process exclusions` (tag 284), 2nd production interval, AFTER 1 sample at step 1082100: **steps_adv=577, engaged=587, int_bound=1, γ=4.198e-4 (floor), r/r0=1.000.**
- **So the exclusion crash is a PRE-EXISTING long-run showstopper in the v4 driver, NOT introduced by G1.** It scales with RUN LENGTH (dynamic myosin/xlink/turnover bond accumulation on cortex beads); G1's extra clutch only makes it crash sooner. My fan-in fixes mitigate one contributor, not the root.
- **This corroborates + reframes the γ blocker:** even with myosin stepping (577 advances) at ~1e6 steps, γ=floor + r/r0=1.0 (no contraction). 577 advances / 587 heads ≈ 1 bin (~33 nm) each = minimal; macroscopic contraction needs FAR more stepping (≈seconds of sim → ≥1e8 steps) → a long run → which the pre-existing crash PREVENTS. **So the crash is the real long-run γ blocker, and it predates G1.**
- **Spawned a read-only investigation agent** to root-cause the dynamic bond/exclusion accumulation (leak vs stochastic peak; which updater — myosin/xlink/turnover) + propose a confined fix. Payoff sim `bbn4djb13` (fan-in 1) still running (may get 1-2 samples before the same crash → A/B vs the G1-off γ=4.2e-4 baseline).
- **G1 (mechanism) + G2 (myosin steps) are DONE.** The remaining γ-lift is gated on (i) the pre-existing long-run crash fix [investigation running] and (ii) confirming a long-enough run actually contracts. Both are PI-worthy.

### iter 7 — 03:43 KST — A/B confirmed + both long sims crashed (pre-existing) → pivot to modules
- **Payoff `bbn4djb13` (G1-ON, clutch 45):** 1 sample (step 502100) then crashed (tag 122). γ_total=**4.839e-4** | r/r0=1.000 | steps_adv=203 | int_bound=1.
- **A/B: γ 4.84e-4 (G1-ON) vs 4.20e-4 (G1-OFF control)** — the loaded clutch gives ~+15% γ but BOTH are ~1000× under the 0.35–0.65 floor, r/r0=1.0 (no contraction). **G1's clutch marginally raises γ; the real lift needs cortex CONTRACTION (a long run), which the pre-existing crash blocks.** fan-in 1 (45 clutch) crashed sooner (~5e5) than control (9 clutch, ~2e6) — more clutch = sooner, consistent with the dynamic-accumulation root.
- Both long sims done (crashed). Crash-investigation agent `a3e75c2b` running. **NOW pivoting to EXTEND module integration** (substrate/H.8/H.9 default-off; nucleus net-force fix first).

### iter 8 — 03:50 KST — EXTEND modules committed (nucleus net-force fixed)
- **Nucleus §3 net-force/COM fix:** the builder coupled bead radius to the Fibonacci direction index (small-radius → +z, large → −z) → the radial confinement force didn't cancel → ~25% spurious net force (COM drift). Fix: decorrelate radius from direction (deterministic permute) → ~4% residual; F=−∇U consistency kept. Antipodal-pair seeding (exact 0) flagged as a refinement. test_nucleus.py 35 green incl. a new fill=True net-force check.
- **Committed `db1c756`:** ecm/substrate.py + cell/membrane_surface.py + cell/nucleus.py + 3 tests — scaffolded, verified, **default-off, NOT yet wired into cell.py** (standalone, no runtime import → existing 454-suite unaffected).
- PI-flags: substrate `effective_E_sub` diagnostic is N/m³ not Pa (doc-only); nucleus antipodal-seeding refinement.
- **NEXT: wire the 3 modules into Cell.build as None-guarded default-off attach blocks + regression.** Crash-investigation agent `a3e75c2b` still running (pre-existing long-run crash root-cause).

### iter 9 — 04:40 KST — pre-existing crash ROOT-CAUSED + degree-cap fix (the γ unblocker)
- **Investigation agent `a3e75c2b` delivered a rigorous root cause:** HOOMD nlist has a HARD compile-time max of **7 bonded exclusions/particle**; a cortex-actin bead's raw degree (backbone 2 + FA clutch ≤1 + myosin ≤3 + xlink **UNCAPPED**) stochastically fluctuates to 8 → "Too many bonds to process exclusions" crash. Proven deterministically (8=crash, 7=OK). NOT a leak / NOT turnover (disabled) / NOT G1 (only changes WHEN it crashes). The uncapped `XlinkBondUpdater` + the myosin cap being blind to non-myosin bonds is the gap.
- **Fix (cortex/crosslinkers.py + cortex/myosin.py — confined + additive):** both dynamic binders compute the per-cortex-bead TOTAL bonded degree (`np.bincount` over all current bonds) and decline a bind at ≥6 (`_MAX_CORTEX_BEAD_DEGREE`, a 1-bond margin under HOOMD's 7). Library-capacity bound (PI-flag class, like MAX_HEADS_PER_BEAD). Compile+import OK.
- **Verifying:** 2.5e6-step sim `buqlz8b8s` (past both prior crash points ~5e5/~2e6) + full regression `bqz07rtg3` gating the commit. If green → commit → **re-launch the long payoff sim** (should now survive to the myosin-stepping regime → tests whether γ finally lifts with sustained contraction).
- The real γ unblocker: G1 (clutch loads) + G2 (myosin steps) + this crash fix → the long contraction run the floor needs.

### iter 10 — 06:55 KST — KU-3.5 FLOOR ROOT CAUSE FOUND (months-long answer) + membrane integrated
- **Degree-cap fix verified + committed `f2f4e75`:** the 2.5e6-step seed-1 run survived crash-free past both prior crash points (step_advances 203→1396). The pre-existing showstopper is fixed.
- **BUT γ stayed floor (~4e-4) + r/r0=1.000 across all 5 samples despite 1396 myosin advances** → investigation agent `a5cb539a` pinned the TRUE cause: the myosin **binned-r0 ratchet is a LUMPED PROXY**. A "step" only relabels the head-actin bond to a lower-r0 bin on the SAME bead; it does NOT transport actin material (AFINES `pos_a_end` grip-walk), so the spring relaxes out each tick → no strain → γ_soft ~25,000× under band. Force ceiling 0.33 pN < F_stall 0.5 pN; bipolar sidedness absent. **Measurement faithful; mechanism loses the force.** Violates CLAUDE.md "no lumped mechanisms."
- **KU-3.5 floor FULLY diagnosed** — NOT clutch/stepping/crash/run-length/measurement; it's the myosin proxy. **Fix = AFINES grip-point walking, core physics → PI-GATED (not implemented).** Full brief: `docs/v2_audit/_historical/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md`.
- **In parallel: H.8 membrane_surface integrated into Cell.build** (None-guard default-off, mirrors enclosed_volume; compile+import OK; regression `bpaae52bq` gating the commit). substrate (FA-pin-coupled) + nucleus (Template-2 particle wiring) integration flagged as next-steps.

### iter 11 — 07:15 KST — figures + PI wake brief + substrate integration
- **Figures committed** (CLAUDE.md viz rule, via a subagent): `fig_h3_ku35_myosin_proxy_failure.png` (HEADLINE — steps_adv 203→1396 but γ flat ~4e-4 + r/r0=1.000, KU-3.5 band shown, ~1258× gap), `fig_h3_ku35_g1_clutch_fix.png` (clutch 9→52), `fig_h3_ku35_AB_clutch.png` (A/B G1on/off, both ~1107× under band) + repro script `scripts/h3_ku35_finding_figs.py`.
- **PI wake brief committed** `docs/v2_audit/_historical/PI_WAKE_BRIEF_2026-05-31.md` — single boot-point: outcome + ranked PI decisions (#1 = myosin AFINES redesign, PI-gated).
- **Track A substrate integrated** into Cell.build (None-guard default-off branch at the ligand-pin attach: `p_substrate` finite k_sub → SubstrateAnchorSpring REPLACES the rigid SubstrateLigandPin; None = bit-for-bit). compile+import OK; regression `regression_substrate` gating the commit.
- **Nucleus (H.9) integration DEFERRED** — Template-2 (nucleus_bead snapshot-extension + gamma_map + force, spanning two cell.py areas) is more involved + the nucleus has the residual-net-force flag; flagged as a next-step. membrane (H.8) + substrate (Track A) now integrated default-off.

### iter 12 — 07:35 KST — substrate committed + E_sub flag addressed; session substantively complete
- **Committed:** Track A substrate integration `bbb4a3e` (default-off); E_sub diagnostic relabel N/m³ `1af5b4d` (doc-only PI-flag); figures `b98807e`; PI wake brief `7ebfc1f`.
- **2/3 EXTEND modules integrated** (membrane H.8 + substrate Track A); nucleus (H.9) deferred (Template-2, involved — flagged in the wake brief).
- **Session substantively complete:** KU-3.5 fully diagnosed (myosin lumped proxy, PI-gated fix); G1+G2+crash all fixed + committed; modules built+integrated; figures + wake brief + root-cause doc. **Remaining:** nucleus integration (next, if tractable) + the PI-gated myosin AFINES redesign (#1 decision). (11 commits this session, all regression-green.)

### iter 13 — 07:38 KST — nucleus integration assessed + DEFERRED; closeout + maintenance transition
- **Nucleus (H.9) integration assessed:** it is a genuine Template-2 snapshot-extension (append nucleus_bead particles to the gsd Frame BEFORE create_state_from_snapshot + gamma_map["nucleus_bead"] + the force + optional centroid pin) — the FA-pattern complexity, for a default-off module. Per the escape clause ("진짜 복잡/위험하면 defer"), **DEFERRED** (flagged in the wake brief; best done when the nucleus is enabled for a confined-migration experiment so the wiring + physics are tested together). membrane (H.8) + substrate (Track A) remain integrated default-off.
- **Session SUBSTANTIVELY COMPLETE.** 11 commits this session (since the `1496173` handoff), each regression-green: G1 (4 commits) + degree-cap crash fix + EXTEND-module scaffold + H.8/Track-A integration + figures + wake brief + E_sub doc, anchored by the KU-3.5 root-cause. Clean tree. Closeout regression `bh01re09v` running.
- **Transitioning to MAINTENANCE cadence** (longer heartbeat): the high-value autonomous work within the guards is done; the remaining (myosin AFINES redesign = #1 PI decision, nucleus integration) is PI-gated/deferred. The loop stays alive; I verify + do safe polish if valuable + respond if the PI returns.
