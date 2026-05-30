# KU-3.5 v4 floor diagnosis + EXTEND-scope plan — session handoff (2026-05-30)

> **Purpose.** Durable handoff so a fresh Lead session can continue without re-deriving.
> Produced by the Lead (ultracode) session on 2026-05-30 from a 5-reader read-only
> `understand-floor-and-extend` workflow (adversarial audit included), **before any code
> was written**. No commits this session; HEAD unchanged at `6e1c6ac`; working tree clean.
>
> **Boot state at handoff:** branch `phase1/h3-cortex` @ `6e1c6ac`; suite **454 passed /
> 17 skipped / 0 failed**; env `ffn_sim`, HOOMD 7.0.1, Python 3.13. No running procs.
>
> **Companion docs:** floor diagnosis `outputs/h3/production/ku35_v4_smoke_diagnosis.md`;
> PI queue `docs/v2_audit/PI_DECISION_QUEUE_2026-05-30.md`; mechanism audit
> `docs/v2_audit/MECHANISM_AUDIT_2026-05-30.md`; EXTEND `docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`;
> ECM platform `docs/ECM_PLATFORM_EXTENSIBILITY.md`; briefs `docs/briefs/H{8,9,10}_*.md`.

---

## 0. TL;DR (read this first)

The KU-3.5 v4 γ floor (≈3.1e-4 mN/m, ~1,600× under [0.35, 0.65]) is **two independent
mechanism gaps, both OUTSIDE the cf7ba49 B1/B2 files**. The PI's working hypothesis
("the bug is in cf7ba49's `dt_reconcile.py` / `equilibration.py` / `cell.py`") is
**refuted by an adversarial audit + arithmetic** — those three files are clean.

1. **Blocker (a) — FA clutch never loads = a geometry-PLACEMENT gap.** The cortex shell
   is built centered on the origin (bottom at z ≈ −R_cell ≈ −10 µm); substrate ligands
   sit at z = 0; **nothing translates the cell into contact**. The v4 driver *believes*
   `equilibrate=True` "settles the cell onto the substrate," but the equilibration code
   is pure **in-place** local force relaxation — no body force, no CoM translation. So
   `n_integrin_bound_final = 0`, `r/r0 = 1.0000`, no substrate→cortex tension. **This is
   the load-bearing fix and it is a physics-design decision (3 options below) → PI gate.**
   Track A's `k_sub` substrate-compliance spring is a *separate* axis (it tunes stiffness
   *once the clutch is loaded*) and does **not** fix this.

2. **Blocker (b) — myosin `step_advances = 0` = almost certainly a too-short-sim
   ARTIFACT**, not a bug. The smoke runs ~1,500 production steps (~2.2 µs simulated); the
   fractional-accumulator ratchet needs ~2.53e7 steps (≈165 ms) for the **first** bin
   advance. step_advances=0 at smoke length is arithmetically forced, not a dead gate.
   **One residual risk** keeps this honest: mean bound lifetime (~0.1 s at k_off0=10/s)
   is *comparable to* the 165 ms step-time, so a head may unbind before it ever steps —
   the long diagnostic sim (command in §1.2) is the decider. **PI explicitly asked to
   re-confirm with a longer sim first.** Not launched this session (paused per PI).

**EXTEND scope** is well-scoped and additive. H.8 (membrane surface) and H.9 (nucleus)
are CLEAN Template-1/2 modules (default-off, bit-for-bit when absent, zero integrator
change). H.10 (cytoplasm viscoelasticity) is **integrator-adjacent → frozen-`integrator/`
PI-gate**; do not build until ratified. Track C = re-derive KU-3.5/KU-3.1 as composite
gates (a gate-contract change → PI). The verified additive template is `cortex/erm.py` +
`cortex/enclosed_volume.py`.

---

## 1. Floor bug — full analysis

### 1.0 What the smoke reported (`ku35_v4_smoke.json`, diagnosis doc)
`γ_total = 3.13e-4 mN/m` (soft 2.1e-6 + rigid 3.11e-4); `n_fa_clutch_bonds_at_build = 9/52`;
`n_integrin_bound_final = 0`; `bind_total = 83` (myosin heads engage); `step_advances = 0`;
`r/r0 = 1.0000`. Reconciled dt `1.30e-8 → 6.52e-9 s` (FA `k_int_bare` is the stiffest soft bond).

### 1.1 Blocker (a): FA clutch never loads — geometry-placement gap

**Evidence chain (file:line):**
- Cortex centroid = origin, radius R_cell ~10 µm: `cortex.py:638` `_sample_sphere_surface(rng, F, p.R_cell)` (helper `cortex.py:562-564` returns points centered on 0) → cortex beads span z ∈ [−10, +10] µm. Box L = 3·R_cell = 30 µm (`cortex.py:438`).
- Substrate ligands at z=0: `fa.py:434,450-454`; integrins a few tens of nm above (`fa.py:428`); FA centres scattered across the full ±15 µm box face (`fa.py:393-395`), **not** confined under the cell.
- The disjoint-layer fact is stated by the code itself: `cell.py:690-699` ("the substrate plane (z≈0) and the cortex shell (r≈R_cell ~10 µm) are spatially DISJOINT … With the physical radius ZERO clutch bonds form at construction").
- S2 static clutch capture: `cell.py:388-403`; bond formed only if nearest cortex bead within `capture_radius` (default `capture_radius_R_FA` ≈ 1.5 µm, `cell.py:708-712`). The 9 that form are built **force-free** (`r0` = exact construction separation, `cell.py:396,440-446`) → 0 J at build, exert no load until the cell moves.
- S1 dynamic Pereverzev catch-bond count `n_integrin_bound` is **hardwired to 0** at `h3_ku35_v4_fa.py:423`; the dynamic load-and-fail kinetics is the explicit **S5 TODO** (`cell.py:700`; `h3_ku35_v4_fa.py:418`).
- **The "settle" is claimed but not implemented.** The driver `--help`/docstring asserts `equilibrate=True` is "a pre-production settle that brings the cell down onto the substrate" via "the internal substrate anchor that pulls the cortex shell onto the plane." But:
  - `equilibration.py:219-248` + `ecm/equilibrate.py::_softstart_step` (`ecm/equilibrate.py:135-173`) move each particle **along its own net_force, clipped** — local relaxation only. **No body force, no gravity, no z-translation, no CoM recentering.**
  - `SubstrateLigandPin` (`cell.py:180-239`, reset `act()` at `cell.py:230-239`) pins the **LIGANDS** to z=0 every step — it holds the *ground* immobile; it does **not** pull the *cortex* down.
  - Runtime corroboration: `r/r0 = 1.0000` (cortex radius unchanged) and `n_integrin_bound = 0`.
  - The only existing workaround is `--fa-capture-radius ~3–5e-6` to *fake* contact by inflating capture range (`h3_ku35_v4_fa.py:351-354`) — unphysical, not a fix.

**Verdict:** the cortex is never brought into contact with the substrate; the settle mechanism the driver assumes is absent. **Fix is geometry placement, a flagged physics-design decision (diagnosis doc: "not unsupervised code patches"). → PI gate.**

**Fix options (recommend PI pick before implementation):**
1. **Construction-time placement — shift the cortex centroid to z = R_cell** (or R_cell − contact_depth) so the bottom hemisphere sits at z≈0 within `capture_radius` of the ligands. Cleanest, physical, no new force. *(Lead's lean: this or #2.)*
2. **Seed FAs at the contact footprint** — restrict ligand/integrin xy to the cell's footprint and place ligands directly beneath the cortex lower pole so nearest-bead < capture_radius by construction.
3. **A real settle step** — add a downward body force / contact constraint during equilibration that translates the cell onto z=0 (matches the driver's stated intent, but is equilibration-adjacent → near the B2 surface; handle with care).
4. *(workaround, not a fix)* `--fa-capture-radius` inflation — already available, unphysical; use only to unblock a smoke, never for a production γ claim.

**Note:** even after contact, S1 dynamic clutch kinetics (S5 TODO) must be implemented for the catch bonds to actually carry/fail load. Placement fix is necessary but the S5 kinetics is the co-requisite for a *sustained* loaded clutch.

### 1.2 Blocker (b): myosin step_advances=0 — artifact (long-sim decides)

**Increment path (verbatim, `cortex/myosin.py:983-1002`, `MyosinStepUpdater.act()` step 3):**
fractional accumulator `_head_step_accum[h] += d_bin` with `d_bin = v_step·batch_dt/bin_width`;
`adv = floor(accum)`; step advances only when `adv ≥ 1`. `v_step = hill_velocity_clamped`
(`myosin.py:979`) is positive for sub-stall heads → accumulator grows monotonically. Reset to
0 on unbind (`myosin.py:837`). **No gate that structurally pins adv=0** — it just grows slowly.
Comment `myosin.py:732-737` documents the accumulator was added *because* `round(d_bin)` was
"always 0 at the model's timescale (v·batch_dt ≪ bin_width)" — step_advances=0 in a short run
is the *known, intended* behavior.

**Arithmetic (smoke is ~75,000× too short):**
- dt = 6.515e-9 s; batch_steps=100 → batch_dt = 6.515e-7 s; bin_width = 3.3e-7/10 = 3.30e-8 m; v0 = 0.2 µm/s.
- d_bin_max = 0.2e-6·6.515e-7/3.30e-8 = **3.95e-6 bins/tick** → 1 bin needs 2.53e5 ticks = **2.53e7 BAOAB steps ≈ 165 ms**. (Exactly the 단계-12 plateau figure.)
- Smoke production ≈ 1,500 steps (~14 µs) → accumulator ~6e-5 ≪ 1. step_advances=0 is forced.

**Residual risk (the real thing the long sim decides):** at k_off0=10/s, mean bound lifetime
~0.1 s ≈ 153k ticks, *comparable to* the 253k ticks needed for one advance. A head may unbind
(~100 ms) before its first step (~165 ms). The long sim tests whether **any** head survives
long enough to cross a bin boundary.

**EXACT long diagnostic command (flags validated against `--help`):**
```
cd /Users/sw1/ffn_cellsim && conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_v4_fa \
  --n-fil 60 --seed 1 --n-warmup 2000 --equilibrate-steps 2000 \
  --n-sample 40 --interval 1080000 \
  --device cpu --skip-integrity --no-resume \
  --out ffn_sim/outputs/h3/production/ku35_myostep_long.json
```
- ≈ 4.32e7 production steps (= n_sample·interval) ≈ 280 ms (×1.6 margin past the 165 ms threshold).
- **Decision rule:** watch the `steps_adv=` column in `PROGRESS sample=` lines.
  - first nonzero past sample ~24 (≈step 2.5e7) → **ARTIFACT confirmed**, MECHANISM_AUDIT §A2 "binned-ratchet defect" **refuted**; no ratchet code change warranted (165 ms-to-first-step vs 100 ms lifetime may still be a *parameter* question for PI, not a bug).
  - stays 0 past ~2.8e7 with healthy bind_total + non-trivial lifetime → **promote to bug** (A2): the lifetime-vs-step-time race kills every advance; bin_width (33 nm) too coarse or v0·batch_dt too small → re-design (finer bins, or step bond r0 continuously instead of binned).
- Wall-time: constrained build (M-SHAKE+Fixman+per-100-step Python updaters+FA+KDTree) ≈ 1–5e4 steps/s single-core → **~15–70 min on gbook (central ~30 min)**. Lighter alt: `--n-sample 28 --interval 1000000` (~2.8e7, ~10% margin) ~20–45 min.
- **Do NOT pass `--smoke`** (it clobbers n-sample/interval back to 3/500). `--n-fil` sets cortex size; equilibrate+reconcile_dt are hard-wired True in the production build (not CLI).

### 1.3 PI hypothesis refuted: cf7ba49 (B1/B2) is clean

- **`dt_reconcile.py` (B1) correct.** `compute_global_cfl_dt` (`dt_reconcile.py:230-413`) = per-element τ=γ/k list, `binding=min(active)`, `dt=safety·τ` (`:388-394`). Standard overdamped bound. Smoke's 1.30e-8→6.52e-9 is **arithmetically exact** (γ_integrin=6π·6.913e-4·5e-9=6.515e-11; τ_fa=6.515e-8; dt=0.1·τ=6.515e-9) and only a **~2× reduction** vs cortex — benign, cannot freeze contractility (would need 10³–10⁴×). dt only ever *lowers* (`if cfl_result.dt_min < dt_used`), never raises into instability; SHAKE timescale assert (`:402-411`) only fires when requested dt *exceeds* the bound (driver pre-empts it).
- **`equilibration.py` (B2) correct.** Phase-A clipped steepest-descent + Phase-B BAOAB drain; detaches/re-attaches BAOAB correctly (`:211-213`/`:246`). Cannot relax away active machinery: clutch S2 bonds are force-free by construction (nothing to relax), S1 is runtime-populated (nothing pre-formed to break), Sanity-Gate §3 guarantees energy only drains. It does **not** translate the cell (the actual gap, §1.1).
- **`cell.py` (cf7ba49 diff) correct order, default-off.** Reconcile runs *before* integrator/BAOAB attach (no "reconcile after attach" bug); equilibrate runs at the very end *after* all bonds + the ligand pin exist (no "equilibrate before placement" bug). Defaults `reconcile_dt=False`, `equilibrate=False` → bit-for-bit identical when off. The two on-path bugs the commit message cites were caught by the test + fixed pre-commit (not live).

---

## 2. EXTEND scope — analysis & build specs

### 2.0 The verified additive default-off template (copy this)
- **Template-1** (custom force over a tag range): exemplars `cortex/enclosed_volume.py` (`EnclosedVolumePressure`, `:365-470`) and `cortex/erm.py` (`ERMHarmonic`, `:184-234`). A `hoomd.md.force.Custom` subclass that masks a tag range, computes per-bead forces, writes `self.cpu_local_force_arrays` (read by BAOAB each step — no integrator registration beyond `ig.forces.append`). Attach helper carries a **CFL gate** τ=γ_b/k_eff, raises if dt>safety·τ (`erm.py:283-298`, `enclosed_volume.py:525-541`).
- **Default-OFF gate** (load-bearing, `cell.py:923-940`): the builder kwarg `p_X: ResolvedX | None = None` (`cell.py:514`). `if p_X is not None: attach...(); else skip`. When None, nothing is appended to `ig.forces` → **bit-for-bit identical run**. The None-default param *is* the gate. **Modules must have zero import side-effects / no global registration.**
- **Template-2** (new particle subsystem): append particles/bonds to the `gsd.hoomd.Frame` *before* `create_state_from_snapshot` (FA exemplar `cell.py:339-479`); add the `gamma_map[type_name]` entry *before* `make_baoab_updater` (`cell.py:945-964`); optional quasi-static pin = post-BAOAB position-reset Action (precedent `SubstrateLigandPin`, `WaveMembranePin` `cell.py:911-920`).
- **BAOAB tag contract** (`integrator/baoab.py`): attach-time coverage assert (`:289`) raises if any present type lacks a gamma entry; runtime growth hook `_extend_tag_buffers` (`:324-384`) handles appended tags (γ from `gamma_map[type_name]`, `bd_prefactor=√(kT/(2γdt))`, `prv_rnds=0`). **Append-only, dense tags; shrinkage → RuntimeError.** Forget a gamma entry → loud halt at attach (never silent).

### 2.1 Track A — `ecm/substrate.py` (tunable k_sub anchor spring)
- Class `SubstrateAnchorSpring(md.force.Custom)`: F_i = −k_sub·(r_i − anchor_i) over ligand tags; mirror `ERMHarmonic`. Helper `attach_substrate_spring(...)` with CFL gate τ=γ_ligand/k_sub.
- **Default = `pin:true` keeps the existing `SubstrateLigandPin` verbatim** (the only bit-for-bit-identical path — a finite-large k_sub is *not* identical because it adds a force). The spring attaches only for finite k_sub (compliant substrate).
- `k_sub` is a REQUIRED param (no magic default); preset layer supplies it. Dimensional bridge `E_sub ≈ k_sub·n_ligand/A_substrate` (KU-1.V.1), do not hard-code E_sub. Needs a CFL entry in `dt_reconcile.compute_global_cfl_dt` when finite.
- **This is the substrate-compliance axis, NOT the floor fix** (see §1.1). Default-rigid preserves all current FA runs.

### 2.2 Track B — H.8 / H.9 (CLEAN) + H.10 (PI-gate)
- **H.8 `cell/membrane_surface.py`** — Template-1 over cortex shell tags [0, n_cortex_actin). Class `MembraneSurfaceTension`. Minimal first version: surface tension γ_mem (Laplace dP=2γ/R) + area-elasticity K_A; Helfrich κ_m as a documented TODO. **KU-3.B1:** T∈0.03–0.3 mN/m; K_A≈0.24 N/m; κ_m≈1e-19 J (10–30 kT); γ_MCA≈1e-5 J/m²; tether f_t=2π√(2κ_m(T_m+γ_MCA))∈5–40 pN. CLEAN, no integrator change. **Name-collision warning:** `cell/membrane.py` already exists and is the H.5 leading-edge **load reservoir** (adds no particles; `membrane.py:1-94`, wired `cell.py:1075-1088`) — H.8 is a *new, different* module; do not touch `membrane.py`. Promote `cortex/erm.py` to membrane–cortex linker as part of H.8.
- **H.9 `cell/nucleus.py`** — Template-1 `NucleusConfinement` (radial harmonic F=−k_eff(r−R_nuc)r̂, two-regime: chromatin small-strain + lamin strain-stiffening) + Template-2 `nucleus_bead` group (Lead appends to snapshot + gamma_map + optional centroid pin). **KU-3.B2:** E_nuc 1–10 kPa; nucleus:cytoplasm 1.4–5× in-situ (**do not hard-code 10×**); lamin-A scaling (η~laminA³, E~laminA^0.5); critical pore ~7 µm². CLEAN, no integrator change. Couples to **KU-1.V.3.3** (pore limit) — in-vivo ECM confinement presets gate on H.9.
- **H.10 `cell/cytoplasm.py`** — viscoelastic/poroelastic background **WITH MEMORY** → NOT Template-1/2 alone; touches the drag model (memory kernel or poroelastic filler). **integrator-adjacent → frozen-`integrator/` PI-gate. DO NOT BUILD until ratified.** KU-3.B3: η contrast MCF7≈5.3×MDA; G 30–80 Pa; poroelastic τ<0.5 s, Dp≈40–60 µm²/s. EXTEND §B.5 flags this as the lone hard module.

### 2.3 Track C — composite-gate re-derivation (after Track B)
Current KU-3.5 (tension) + KU-3.1 (rounding) attribute **all** whole-cell mechanics to cortex
(EXTEND §A.3 Issue 2 — the real scientific risk). Re-derive as composite sums:
- **KU-3.5:** γ_total = γ_membrane(H.8) + γ_cortex(H.3). (Membrane further T = T_m + γ_MCA.)
- **KU-3.1:** rounding = cortex + nucleus(H.9) terms.
This is a **gate-contract change → PI gate** (no-gate-loosening rule). Until done, label current
results "cortex-attributed", not "whole-cell" (EXTEND §A.5). Code is kept; gates RE-VALIDATE
(not a teardown). Especially load-bearing for cancer cells (nucleus+cytoplasm co-dominant).

---

## 3. Governance — what the next session needs PI to ratify

| # | Decision | Type | Recommendation |
|---|---|---|---|
| G1 | **Floor (a) geometry-placement fix** — pick option 1/2/3 (§1.1) | physics-design | Lead leans #1 (centroid→z=R_cell) or #2 (footprint seeding); cleanest + physical |
| G2 | **Floor (b) verdict** — gated on the long-sim result (§1.2) | data | run long sim → artifact (no change) vs bug (A2 ratchet redesign) |
| G3 | **S1 dynamic clutch kinetics (S5 TODO)** implementation | mechanism | co-req of G1 for a sustained loaded clutch |
| G4 | **H.10 cytoplasm** integrator-adjacent design | frozen-`integrator/` | DO NOT build until PI ratifies; H.8/H.9 proceed without it |
| G5 | **Track C composite KU-3.5/KU-3.1** re-derivation | gate-contract | after Track B; re-label current as "cortex-attributed" now |
| G6 | **CFL softening** if H.8/H.9 stiffness violates the gate | magic-number/gate | surface with the τ=γ/k numbers; do not soften silently (k_ERM precedent) |

**Still-open from `PI_DECISION_QUEUE_2026-05-30.md`** (already-applied: B1/B2 `cf7ba49`; corrections `c2d0b01` = filamin catch / δ_cap / branch-angle / KU-5.1 band):
- **A2** myosin physical stepping — *gated on G2*; do not touch the ratchet until the long sim says "bug".
- **A4** KU-3.5 interphase vs metaphase regime; **A3** KU-5.1 band (done?); **C1** integrin F* 6.99→~30 pN (Kong 2009); **C3** myosin v0; **A5** Arp2/3 debranching; **D1** `ffn/foundation` push; **D2/D3** stale branches / sync-conflict shadows.
**EXTEND/ECM open items** (`docs/*EXTEND*`, `docs/ECM_PLATFORM*`): ratify EXTEND vs rebuild (PI directed proceeding in this session's boot prompt); roadmap amendment H.8/9/10 as first-class units; ratify 4-axis `ecm_condition` preset + rigid-pin→k_sub unlock + (b) col-I-coated-rigid as DEFAULT preset; ligand species set; **in-vivo ECM presets gate on H.9 pore**.

---

## 4. Execution state this session (for an accurate restart)

- **Understand workflow** `wcvcrlhts` — DONE (5 readers, ~308k tok). Full maps: `tasks/wcvcrlhts.output` (transcript dir under the session `subagents/workflows/wf_836be351-3ea`).
- **Build workflow** `wwl0i634g` (would scaffold `ecm/substrate.py`, `cell/membrane_surface.py`, `cell/nucleus.py` + tests via scaffold→adversarial-verify pipeline) — **STOPPED before any file was written** (PI said pause). Tree verified clean: none of the 6 target files exist; 0 commits; HEAD still `6e1c6ac`. Script saved at `workflows/scripts/build-extend-modules-wf_7593473b-1ec.js` (re-runnable when PI green-lights).
- **Long myosin sim** — **NOT launched** (paused). Command ready in §1.2.
- **gbook probe** — reachable, 16 cores, repo `/home/sungwook/ffn_cellsim`. ⚠️ `git HEAD = d9249e5` (OLD, pre-v4). **But Syncthing syncs working-tree files, not git** (`.stignore` excludes `.git`/`.claude`), so the driver file may already be present despite the stale git log — **verify the actual file (`ls -la ffn_sim/scripts/h3_ku35_v4_fa.py` + a quick `--help` on gbook) before launching there.** Per memory, CPU jobs default to gbook (concurrent, single-core each).
- **Guards honored:** cell.py untouched; no governed-surface edits; no git mutations by subagents; no commits.

---

## 5. Recommended next-session order

1. **Launch the long myosin diagnostic** (§1.2) in the background (gbook per default, after verifying file sync; or local). It's the slow long-pole — start it first. Grep `steps_adv=`; apply the §1.2 decision rule.
2. **Surface G1 (geometry fix) to PI** with the 3 options; on ratification, implement the placement fix (cell.py single-writer, Lead). This is the actual floor lift for blocker (a).
3. **In parallel, re-run the build workflow** (`build-extend-modules-wf_7593473b-1ec.js`) to scaffold + verify `substrate.py` / `membrane_surface.py` / `nucleus.py` (additive, default-off, cell.py untouched). Then Lead-integrate the None-guarded attach blocks serially + run regression (must stay 454/17/0 with modules off).
4. **Hold H.10** for PI (G4). **Stage Track C** (G5) re-derivation once H.8/H.9 land; re-label current KU-3.5/3.1 as "cortex-attributed" in the interim.
5. Re-measure KU-3.5 v4 only after G1 (+ S5 kinetics, G3) and the myosin verdict (G2) are resolved — not before (would reproduce the floor and waste gbook hours, per the diagnosis doc).

---

*Handoff authored 2026-05-30 by the Lead session. Analysis is read-only-derived and
adversarially audited; the placement-gap and artifact verdicts are corroborated by the
runtime smoke (`r/r0=1.0`, `n_bound=0`, `step_advances=0`) but the definitive confirmations
are (b) the long sim and (a) an empirical placement-fix re-measure. No gate was declared PASS.*
