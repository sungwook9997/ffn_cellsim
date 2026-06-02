# Autonomous PARALLEL-prep session log — 2026-06-02

**Start:** 2026-06-02 02:17 KST (epoch 1780334226).
**Mode:** PARALLEL prep session, running *alongside* an active Lead session.
**Mandate (PI, /loop, KO):** "너는 리드 세션과 병렬로 동작한다. 리드는 프로덕션으로
바쁘다. 새로 추가된 개념 + 이 프로덕션 이후 구축해야 할 것들 중 병렬 가능한 것을 모두
하라. 충돌만 안 생기게 조심해서 자율주행하라. /loop로 20분마다 진행 + 리드 상태를
같이 확인하라."

> This is **NOT** the Lead session. The Lead owns `phase1/h3-cortex`, the working
> tree's mutable runtime files, the GPU production on gbook, the PI-exp validation
> map, and the Notion record. This session does **conflict-free additive prep only**.

## What the Lead is actively doing (read-only observation, 02:17)

- **GPU-main port + NATIVE n_fil=38,000 (266k particle) transport pilot** on gbook
  (RTX A5000). cupy port of the constrained Action → GPU 11–12× at native scale,
  per-step host sync removed. TRANSPORT ENGAGED (myosin walks, adv 0→131 at native).
- **γ still ~100× below band (flat)**: `v0_accel` speeds the walk but not literal
  crosslink binding → desync → no percolation → no tension. v0_accel = transport
  diagnostic, NOT γ. True γ needs a **literal-timescale long run** (next PI-gated).
- **ACTIVE run:** `gpu_native38k_TEST_v1500.{log,grip_walk.gsd}` (grip-walk Tier-1,
  v0_accel=1500×, force-scaling ON), GSD every 40000 steps. Healthy at boot.
- Also touching: `PI_EXP_VALIDATION_MAP.md`, the Obsidian RAG mirror.
- Next (PI-gated, Lead): literal-v0 gold-standard γ run; binder GPU port;
  `ffn/foundation` push (~87 ahead).

## My autonomy contract (the rule a waking PI / the Lead will see)

✅ **Allowed (this session):**
- Write **NEW files only**, each complete-on-write (never leave a half-written file
  in the shared tree), under `ffn_sim/docs/` and — for runnable additive building
  blocks — `ffn_sim/bridge/` + `ffn_sim/tests/`.
- New runtime modules must be **import-isolated**: imported by *no existing file*, so
  the live runtime + the existing test suite are bit-for-bit unaffected.
- Any new test must be **verified GREEN in isolation** (`pytest <file>`) before it is
  left on disk — so it can never red the Lead's collection / regression gate.
- Literature-anchored constants only, with provenance + proxy flags (PI-exp map +
  KU bands). Physics conversions (e.g. F_s = k_BT/x_β) carry no free parameter.

🚫 **NOT crossed (conflict + PI-gate guards):**
- **No edits to Lead-owned tracked files**: `cell/cell.py`, `bridge/fa.py`,
  `bridge/integrin_bonds.py`, `cortex/*`, `ecm/mikado.py`, the briefs,
  `PI_EXP_VALIDATION_MAP.md`, configs in use, `STRUCTURE.md`/`README.md`/`pyproject`.
- **No `integrator/` touch** — BAOAB/dt frozen (H.10's invasive option is PI-gated).
- **No git** (add/commit/checkout/stash/branch) on the shared tree — the index/HEAD
  is shared with the Lead (memory `reference-subagent-no-shared-git`). New files stay
  **untracked**; the Lead/PI integrates + commits at a sequencing point. (Untracked
  files are NOT swept by `git commit -am`.)
- **No `outputs/h3/production/` writes**, no SSH interference with gbook production.
- **No Notion writes** (Lead's authoritative record) — read-only for context.
- **No gate-contract changes, no magic-number invention, no fitting to PI data.**

🛡️ **Safety:** every artifact additive + isolated. If a build needs to touch a
Lead-owned file or the integrator, I STOP and log it as a PI/Lead hand-off item
rather than editing.

## Why these tasks (alignment with the authoritative open items)

From the Notion Dev-Logs board + PI_EXP_VALIDATION_MAP + the EXTEND memo, the two
genuinely-open, parallelizable, non-PI-gated *prep* items are:

1. **Ligand-identity FA extension** — Notion open item, verbatim: *"sequenced after
   FA production lands → extend `bridge/fa.py` ligand identity."* This is exactly the
   PI's "이 프로덕션 이후 구축해야 할 것." The literature is fully anchored (3 dossiers).
   I prep the **ready-to-implement spec + a standalone, tested data module** the Lead
   wires in — without touching `fa.py`.
2. **H.10 cytoplasm** — the ONE EXTEND unit with no code/design artifact yet
   (membrane H.8 + substrate Track A already built+integrated; nucleus H.9 built,
   integration deferred). Notion has `MC-H10-cytoplasm` as `draft`; the brief says
   "scope separately, expect a frozen-`integrator/` PI-gate." I produce the **design +
   literature + sanity-gate doc** (read-only prep, the sanctioned working mode) that
   feeds draft→implemented — NO integrator-touching code.

Deliberately NOT done (PI-gated / deferred / Lead-owned): the literal-v0 γ run, binder
GPU port, `ffn/foundation` push, H.9 nucleus cell.py integration (Template-2, edits
cell.py), cell-type contracts ("의도적 보류 until cell-type instantiation"), KU-3.5/3.1
composite-gate re-derivation (gate-contract).

## Key technical findings (grounding the build)

- The integrin↔ligand off-rate runtime is `IntegrinBondUpdater` driving
  `pereverzev_k_off(F, PereverzevParams(k_s, F_s, k_c, F_c))`
  (`bridge/integrin_bonds.py:102,272`). `validation.pereverzev` is **runtime-imported
  already** (it is the catch-slip *model*, not a forbidden `validation/oracles/` gate).
- **Slip-only bonds are a degenerate Pereverzev** (`k_c = 0`): a pure Bell-Evans slip
  `k_off(F) = k_off0·exp(F·x_β/k_BT)` maps to `k_s = k_off0`, `F_s = k_BT/x_β`. So
  **ligand identity is a per-species `PereverzevParams` swap**, not new mechanism:
  - **col-I → α2β1** = slip-only (k_c=0), from measured x_β / k_off0.
  - **laminin-111 → α6β1** = slip-only proxy (α7β1-invasin shape), scaled weaker than FN.
  - **FN → α5β1** = full catch-slip = the platform's existing KU-2.18 default.
  This keeps the runtime updater untouched; only the param set differs per ligand.
- H.8 `membrane_surface.py` is self-contained (literature constants inline, no YAML) —
  no config gap. `test_substrate_spring.py` already covers `ecm/substrate.py`.

---

## ITERATION LOG (append-only)

### iter 0 — 02:17 KST — boot + plan
- Booted: read CLAUDE.md, the 3 EXTEND/ECM/PI-exp docs, the H8/9/10 briefs, the
  2026-05-31 autonomous log, the Notion Dev-Logs board (open items), and the
  Pereverzev/FA runtime interfaces. Confirmed Lead's active areas + a zero-conflict
  prep surface. Production (`TEST_v1500`) healthy.
- Plan: (1) ligand-identity spec + `bridge/ligand_species.py` + test;
  (2) H.10 cytoplasm design doc; then maintenance-cadence /loop (20 min) watching
  production health + Lead commits while extending prep.

### iter 1 — 02:30 KST — H.10 design + ligand-identity prep BUILT (4 new files, all isolated)
Delivered (all NEW, untracked, conflict-free — no Lead-owned file or `integrator/` touched):
- `docs/H10_CYTOPLASM_DESIGN.md` — the missing EXTEND artifact. 3-tier attach analysis
  (Tier-1 per-type η_eff drag = the MCF7≈5×MDA discriminator, gamma_map-only, FDT-safe
  by construction, PI-gate on magnitude; Tier-2 Markovian-GLE storage = frozen-
  `integrator/` PI-gate; Tier-3 poroelastic filler = Template-2). Hu 2024 bands, the
  3-distinct-moduli guard, full Sanity Gate, candidate VG-H10 gates. Feeds MC-H10
  `draft→implemented` prep.
- `docs/LIGAND_IDENTITY_FA_SPEC.md` — ready-to-implement spec for the Notion open item
  "after FA production → extend `bridge/fa.py` ligand identity". Key design fact:
  ligand identity = a per-species `PereverzevParams` swap (slip = catch disabled), NOT
  new mechanism; exact `resolve_h4` attach points; KU-2.5 catch-peak gate must skip slip.
- `bridge/ligand_species.py` — standalone, import-isolated registry: col-I α2β1 (slip,
  Taubenberger default + locked-open variant), laminin-111 α6β1 (slip PROXY, α7β1-invasin
  shape, weaker than FN), FN α5β1 (KU-2.18 catch-slip). `F_s=k_BT/x_β` conversion;
  `bell_evans_k_off`; `pereverzev_params_for`. Discovered `PereverzevParams.__post_init__`
  forbids `k_c=0` → slip catch disabled via `k_c=1e-9·k_s` (F* → NaN, correct).
- `tests/test_ligand_species.py` — **24 passed in 0.05s** (verified GREEN in isolation;
  no HOOMD, fast; import-isolated ⇒ cannot red the Lead's collection/regression).
Production check @02:30: `TEST_v1500.log` healthy, no crash signature.
**Next:** set up the 20-min /loop (production health + Lead-commit watch + extend prep:
candidate next items = β1-distribution observable design, cancer cell-type param map,
H.8 Helfrich-κ_m discrete-curvature note).

### iter 2 — 02:40 KST — loop cycle 1: cancer cell-type param map
- **Production HEALTHY:** Lead launched a new run `gpu_native38k_v300full.log` (02:38,
  v0_accel=300×, force-scaling ON, grip_walk) — just started (header only), no crash
  signature, no sync-conflict files. Lead commits unchanged (HEAD 877cf40, branch
  phase1/h3-cortex). Conflict risk nil.
- **Built `docs/CANCER_CELLTYPE_PARAM_MAP.md`** (NEW, conflict-free) — the
  how-to-parameterize layer mapping verified bands → sim knobs per cell-type
  (MCF10A/MCF7/MDA): cytoplasm η (H.10 Tier-1, the robust ~5× MCF7/MDA discriminator),
  G (Hu 2024, ordering ≠ η — independence guard), nucleus E_nuc/lamin (with the MCF7/MDA
  lamin GAP flagged, not invented), cortex/membrane (existing knobs), whole-cell E₀ as
  OVERLAY-only (emergent, not a knob). Governing principle: **no global cancer
  modulus** (per-compartment vector); the disputed whole-cell stiffness ordering kept
  as a method-tagged dual band, never "cancer = softer." Distinct from the Notion KU v2
  claims (cross-ref, not duplicate); cell-type *contracts* stay deferred.
- **Next candidates:** β1-distribution observable design (pairs with ligand spec),
  H.8 Helfrich-κ_m discrete-curvature note, 3D-Mikado design.

### iter 3 — 03:00 KST — loop cycle 2: β1-distribution observable (note + module + test)
- **Production WATCH (soft flag, NOT a crash):** `gpu_native38k_v300full.log` static at
  101 bytes (header only) since 02:38; **no `v300full` GSD yet**; nothing synced from
  gbook in ~25 min. No crash signature, no traceback. Most likely native-38k warm-up +
  buffered stdout (GPU 97% saturated; prior runs buffer stdout → flush at samples/
  completion). **Watch threshold:** if still static next cycle (≥40 min) AND no GSD AND
  no relaunch, give the PI a low-key heads-up (possible stall) — still not a crash alarm.
  Lead commits unchanged (HEAD 877cf40).
- **Built (NEW, conflict-free):** `docs/BETA1_DISTRIBUTION_OBSERVABLE.md` +
  `bridge/clutch_spatial.py` + `tests/test_clutch_spatial.py` (**10 passed in 0.07s**,
  isolated GREEN). The PI-exp β1-IF analog: a read-only geometry observable over engaged-
  clutch positions giving `f_edge` (edge-localization vs uniform-disk ref 1−ρ²),
  `angular_cv` (uniformity), `n_engaged`. Distinguishes Bare(diffuse)/Pre(peripheral)/
  Lam4(uniform). Layout-agnostic (positions + engaged mask + center) ⇒ Lead wires to the
  real clutch tags at analysis time; Clark-Evans/Ripley noted as a deferred extension.
  Completes the Layer-1 single-cell toolkit (ligand identity + traction + β1 pattern).
- **Next candidates:** H.8 Helfrich-κ_m discrete-curvature note, 3D-Mikado design.

### iter 4 — 03:20 KST — loop cycle 3: H.8 Helfrich-κ_m discrete-curvature design
- **Production HEALTHY — soft-flag RESOLVED:** `gpu_native38k_v300full.log` advanced to
  214 B at 03:20 (GSD-writer line appeared) + `v300full.grip_walk.gsd` initial frame
  (5376 B) written. Confirmed the prior static window was native-38k **warm-up**, not a
  stall (correctly did NOT cry wolf). No crash signature; Lead commits unchanged
  (HEAD 877cf40); no sync-conflicts.
- **Built `docs/H8_HELFRICH_CURVATURE_DESIGN.md`** (NEW, conflict-free) — closes the one
  documented TODO in the shipped `cell/membrane_surface.py` (κ_m bending deferred). Specs
  the mesh-free discrete-H operator options (local quadric fit for measurement +
  Laplace-Beltrami `Δr=2Hn̂` for the force; cotangent-mesh rejected for a dynamic shell),
  a staged Tier-A (measurement-only, validates κ_m vs the tether oracle) → Tier-B
  (leading-order Helfrich force, default-off) build, and the decisive **analytic sphere
  test** (per-bead H≈1/R; ∮(2H)²dA→16π ⇒ U_bend→8πκ_m, R-independent). CFL/softening
  PI-gate noted (k_ERM precedent).
- **Next candidates:** 3D-Mikado design; then shift toward a consolidated PREP INDEX +
  provenance polish (high-value safe items are nearly exhausted).

### iter 5 — 03:40 KST — loop cycle 4: 3D-Mikado design
- **Production HEALTHY + STREAMING:** `gpu_native38k_v300full.log` 352 B (03:34) —
  **sample 1/40** landed: engaged=23835 (native-scale clutch), drift 3.39e-15 (stable),
  adv=0 (early), g_tot 4.7e-3 mN/m (= the known γ-floor, expected). No crash; Lead
  commits unchanged (HEAD 877cf40); no sync-conflicts.
- **Built `docs/MIKADO_3D_DESIGN.md`** (NEW, conflict-free) — the ECM geometry/
  dimensionality axis (2D→3D). Key grounding: the H.1 sanity gate ALREADY logs a
  2D-Mikado-vs-3D coordination ⟨z⟩ gap (`sanity_gate.py:84,130`, KU-1.3) — a true 3D
  generator closes it (the explicit validation win). Specs: 3D rod placement + nematic
  order S (TACS-3 alignment), distance-threshold crosslinks (3D intersections are
  measure-zero), Maxwell rigidity percolation (z_c 4→6, bending-stabilized below),
  default-2D bit-for-bit, 3D-elasticity validated vs bands (no clean 2D analytic
  oracle). Enables the 3D-gel / tumor-stroma / confined-channel ECM presets.
- **Portfolio now spans all 6 loop candidates + the EXTEND/ECM/PI-exp prep surface.**
  Next: a consolidated PREP_INDEX for the PI/Lead + light provenance polish; high-value
  *new-build* items are essentially exhausted within the conflict-free contract.

### iter 6 — 03:48 KST — loop cycle 5: composite-tension re-validation design
- **Production HEALTHY + STREAMING:** `v300full.log` sample 2/40 (g_tot 4.18e-3 mN/m,
  engaged 23946, drift 3.39e-15); GSD grew to 18.4 MB (active trajectory writes). No
  crash; Lead commits unchanged (HEAD 877cf40); no sync-conflicts.
- **Built `docs/COMPOSITE_TENSION_REVALIDATION_DESIGN.md`** (NEW, conflict-free) — the
  EXTEND memo's flagged "real scientific cost" (Issue 2): re-derive KU-3.5/3.1 as
  composite. Key insight: it is a re-ATTRIBUTION, NOT a band change. Protocol = an
  A/B/C **factorial toggle** of the shipped default-off H.8/H.9 modules
  (A cortex-only=current; B +membrane ⇒ γ_mem = B−A; C +nucleus ⇒ KU-3.1 nuclear term),
  with a **superposition sanity check** (γ_B−γ_A ≈ 2γ_mem/R, catch coupled vs additive).
  The gate-contract relabel itself stays PI-gated (VG-H3-composite seeded `blocked`).
- Next: consolidated PREP_INDEX (capstone, includes this note) for PI/Lead handoff.

### iter 7 — 04:15 KST — loop cycle 6: PREP_INDEX capstone
- **Production HEALTHY + STREAMING:** `v300full.log` samples 3–4/40 (g_tot 3.0→2.7e-3
  mN/m, engaged ~23964, s_grip 172→227 nm walking, adv=0, drift 3.39e-15); GSD 30.9 MB.
  No crash; Lead commits unchanged (HEAD 877cf40); no sync-conflicts.
- **Built `docs/PARALLEL_PREP_INDEX_2026-06-02.md`** (NEW, conflict-free) — single
  handoff entry-point: table of all 11 artifacts (2 runnable+tested modules = 34 green;
  7 design docs), the suggested **integration order** for the Lead, and the **PI
  decision roll-up** aggregated across every doc. Cross-refs Notion (KB-PIV, MC-H8/9/10,
  VG-H3-composite) without duplicating.
- **Prep burst COMPLETE.** All 6 loop candidates + the EXTEND/ECM/PI-exp/composite
  surface covered. Within the conflict-free contract (no Lead-file edit, no git, no
  integrator, no Notion), high-value new-build items are exhausted. **Shifting to
  MAINTENANCE cadence:** the loop stays live — each cycle keeps watching production
  health + Lead commits, and only builds if a genuinely new safe item appears (else
  reports "healthy, holding"). Will surface any fresh production crash to the PI at once.

### iter 8 — 04:29 KST — loop cycle 7: ligand presentation (Bare/Pre/Lam4 availability axis)
- **Production HEALTHY:** `v300full.log` sample 5/40 (g_tot 3.55e-3 mN/m, s_grip 281.6 nm
  walking, engaged 23967, adv=0, drift 3.39e-15). ~1 sample / ~7 min. No crash; Lead
  commits unchanged (HEAD 877cf40); no sync-conflicts.
- **Built `docs/LIGAND_PRESENTATION_MECHANISM.md`** (NEW, conflict-free) — closes the
  LAST open PI-exp Layer-1 mapping TODO (Pre vs Lam4 mechanism). Key: presentation is a
  separate **availability axis** (ligand density × k_on × spatial), ORTHOGONAL to the
  identity/off-rate axis (already done). Recommends availability as primary (avidity
  optional); maps Bare(diffuse)/Pre(peripheral)/Lam4(uniform) → falsifiable against the
  `clutch_spatial` β1 metrics; honest med-confidence + the refuted soluble-FN drop + the
  single-cell↔Layer-2 split flagged. Uses only existing anchors (no new web sweep).
- **Layer-1 single-cell PI-exp mapping is now COMPLETE** (identity + presentation +
  β1 observable + traction). Firmly in MAINTENANCE: subsequent cycles report "healthy,
  holding" unless a genuinely new safe item arises.

### iter 9 — 04:57 KST — loop cycle 8: MAINTENANCE (healthy, holding)
- **Production HEALTHY:** `v300full.log` sample 7/40 (g_tot 3.79e-3 mN/m, s_grip 388 nm,
  engaged 23952, adv=0, drift 3.39e-15); GSD 43.4 MB. No crash; Lead commits unchanged
  (HEAD 877cf40); no sync-conflicts.
- **No new build (deliberate, anti-sprawl):** all candidates + the EXTEND/ECM/PI-exp/
  composite surface are covered; manufacturing filler docs would dilute the handoff.
- **Integrity re-verify:** the 2 shipped modules **34 passed in 0.18s** (GREEN after the
  Syncthing round-trip); all 12 prep files intact. Holding — will act only on a fresh
  production crash, a new Lead commit/direction, or a genuinely new safe high-value item.

### iter 10 — 05:11 KST — MAINTENANCE (healthy, holding)
- Production HEALTHY: `v300full` sample 8/40 (g_tot 4.26e-3 mN/m, s_grip 440 nm, engaged
  23945, drift 3.39e-15); GSD 55.9 MB. No crash; Lead HEAD 877cf40 unchanged; no
  sync-conflicts. No new build (anti-sprawl). Holding.

### iter 11 — 05:39 KST — MAINTENANCE (healthy; TRANSPORT ENGAGED milestone)
- Production HEALTHY + milestone: `v300full` samples 9–10/40 — **adv jumped 0→17424→
  19073** (myosin grip-walk transport engaged at native 38k); s_grip reset 440→133 nm
  (grip-walk cycle); g_tot still floor (~3–4e-3 mN/m) — consistent with the documented
  finding (v0_accel speeds transport but γ needs literal-timescale coherence). drift
  4.87e-15 (stable); GSD 68.4 MB. No crash; Lead HEAD 877cf40 unchanged; no
  sync-conflicts. No new build. Holding.
- **iter 12 — 05:53 KST — MAINTENANCE:** `v300full` sample 11/40 (g_tot 3.37e-3 mN/m,
  s_grip 190 nm, engaged 23964, adv 19390 sustained, drift 3.60e-15). Healthy; no crash;
  Lead HEAD 877cf40; no sync-conflicts. Holding (no new build).

### iter 13 — 06:01 KST — 🔴 FRESH PRODUCTION CRASH (surfaced to PI; I did NOT touch it)
- **`v300full` (v0_accel=300×, force-scaling ON, grip_walk) CRASHED** at ~06:01, right
  after sample 11/40, during `sim.run()`:
  `FloatingPointError: Constrained BAOAB _wrap_into_box_xp: |fractional coord| exceeded
  the int32-image guard (worst |round(f)|=3.255e+08 > 1e+08)`
  (`integrator/constrained_baoab.py:949 act → :192 _wrap_into_box_xp`). A particle blew
  ~3.26e8 box-lengths out in one step; the defensive int32-image guard caught it and
  failed loud (not silent corruption — the guard worked).
- **Timeline:** samples 1–8 adv=0 (no transport), drift ~3.4e-15; sample 9 transport
  ENGAGED (adv 0→17424), drift still ~3.4e-15; 10–11 adv→19390, drift ~3.6–4.9e-15;
  then a SUDDEN single-particle blow-up between 11 and 12 (not slow drift). ⇒ likely a
  momentary over-large force shortly after grip-walk transport engaged (accelerated walk
  × force-scaling ON compounding), or a stiff bond at tiny separation, or dt too large in
  the transport regime. NOTE: Notion board says **"force-scaling OFF = gold-standard"**;
  this run is force-scaling ON (a diagnostic arm) — the instability may be a known
  sensitivity of the ON arm.
- **My actions:** STOP building (already maintenance — nothing to stop). Did NOT touch
  `integrator/`, the driver, or relaunch (frozen / PI-gated / Lead-owned). Surfaced
  prominently to the PI with the traceback + timeline. Data: samples 1–11 + 68 MB GSD
  intact (transport-engaged data captured); recovery (force-scaling-OFF gold-standard /
  smaller dt / relaunch) is a Lead/PI decision. Lead HEAD 877cf40 unchanged.

### iter 14 — 06:21 KST — crash follow-up: read-only diagnostic note
- **Lead has NOT relaunched** (~25 min since crash; no new production artifact; crash log
  unchanged; HEAD 877cf40; no sync-conflicts). Lead may be away or investigating.
- **Built `docs/CRASH_v300full_BAOAB_2026-06-02.md`** (NEW, conflict-free, read-only
  analysis — did NOT touch integrator/driver/config/run). Code-grounded diagnostic:
  the SHAKE finite-guard did NOT fire (position finite-but-enormous), constraint drift
  healthy ~e-15 → a single over-forced particle, NOT slow drift; the int32 guard
  (`constrained_baoab.py:182`) worked correctly (loud-fail). Cause grounded in the driver
  (`:98-118`): v0_accel=300× accelerated walk × `mesoscale_force_scaling` ON, once
  transport engaged (sample 9), delivers an over-large single-step force — consistent
  with the documented diagnostic-accelerant caveat, NOT a new defect. Safe recovery
  options for Lead/PI (no integrator change): (1) re-run force-scaling-OFF gold-standard
  (the ratified call; this was the ON diagnostic arm), (2) lower v0_accel, (3)
  --couple-accel, (4) smaller dt. Data 1–11 + 68 MB GSD intact.
- Holding + watching for the Lead's recovery.

- **iter 15 — 06:41 KST — WATCH (production still down; Lead idle):** `v300full` crash
  unchanged (~45 min down); Lead has not relaunched (no new artifact), HEAD 877cf40, no
  sync-conflicts. Recovery is Lead/PI (I have no run/integrator/SSH access by contract).
  Crash diagnostic + recovery options already delivered (`CRASH_v300full_BAOAB_2026-06-02.md`).
  Integrity re-verify: 2 modules **34 passed**; 11 prep docs intact. No new build. Holding.
- **iter 16 — 07:01 KST — WATCH (unchanged):** production down ~1 h (crash 06:01, no
  relaunch); Lead idle (HEAD 877cf40, no new artifact, no sync-conflicts). No action
  available to me (recovery = Lead/PI). No new build. Holding at 20-min cadence.
- **iter 17+ — ROLLING IDLE WATCH (updated in place to avoid log sprawl):**
  last check **14:41 KST** · production DOWN since 06:01 (~8 h 40 m) · Lead idle, no
  relaunch · HEAD 877cf40 · no sync-conflicts · my prep (34 tests / 11 docs) intact.
  No action available (recovery = Lead/PI; crash diag delivered). Will break out of the
  rolling line + report verbosely the moment state changes (Lead relaunch / new crash /
  new commit / PI direction).
  **→ ROLLING IDLE-WATCH CLOSED at 15:15 (Lead returned).**

### iter 40 — 15:15 KST — 🟢 STATE CHANGE: Lead returned, production relaunched
- After ~9 h idle, the Lead launched **two new grip-walk runs**:
  `gpu_native38k_nxl5k_diag` (15:01, likely n_xl=5000 crosslink diagnostic — cf. commit
  80fdc33 `--n-xl overrides`) and `gpu_native38k_ALIGNED` (15:15, aligned-fiber config).
- **NOTE:** both are **v0_accel=300×, force-scaling ON** — the SAME config class as the
  crashed `v300full`, NOT the force-scaling-OFF gold-standard I flagged. The Lead appears
  to be probing the crash via n_xl / alignment variants. Both at header-only so far; no
  crash signature; HEAD 877cf40 (no new commit).
- **Active-watch resumed:** watching for a crash RECURRENCE at the transport-engage point
  (~sample 9, where v300full died). Will surface immediately if it recurs.

### iter 41 — 15:35 KST — Lead active: new commit + runs warming up
- **New Lead commit `c3fafed`:** *"native γ diagnosis: force-isotropy = the wall (budget
  closed) + meridional-align diagnostic flag + viz"* — Lead diagnosed **force-isotropy as
  the fundamental γ-floor blocker** ("the wall") and added a meridional-align diagnostic
  (the `ALIGNED` run). This is the Lead's γ territory — noted as context, not acted on.
- New runs `nxl5k_diag` / `ALIGNED` still header-only (~20–35 min in) — consistent with
  native-38k warm-up (v300full had a similar long pre-sample window); no crash sig, no
  sync-conflict. Watching for the first sample + crash recurrence at transport-engage.
- No new prep item triggered (force-isotropy is outside my EXTEND/PI-exp surface). Holding
  builds; active-watch continues.

### iter 42 — 15:57 KST — ALIGNED run warming up healthy (no crash recurrence yet)
- `ALIGNED` (active run): advanced to 213 B, GSD-writer line + initial frame (5376 B) —
  same healthy native-38k warm-up pattern as v300full; no crash. Watching for first
  sample + the transport-engage point (~sample 9, where v300full died).
- **NOTION UPLOAD (PI-authorized, overrides the no-Notion clause):** posted the prep
  portfolio (constructive output only — Lead-watch/monitoring excluded per request) as a
  new milestone child page under the Dev-Logs board:
  https://www.notion.so/373120daec5d8194a068c9f3792e9633 (created additively, no board-
  body edit ⇒ no race with the active Lead). 15 files / 34 tests + integration order +
  PI decision roll-up.
- `nxl5k_diag` dormant (header-only since 15:01) — gbook has ONE GPU ⇒ runs are
  sequential; nxl5k was likely superseded by ALIGNED (no crash sig — not a 2nd crash).
- HEAD c3fafed; no sync-conflicts. Active-watch continues; holding builds.

### iter 44 (logged early, ahead of iter 43 entry) — 16:39 KST — 🟢 Lead LANDED my prep
- **New Lead commit `7bc2ea7`: "Land parallel-prep EXTEND/PI-exp artifacts + crash
  diagnostic + gamma-breakdown viz"** — the Lead integrated + committed all 15 of my
  prep files (the EXTEND/PI-exp portfolio + the crash diagnostic). **HANDOFF COMPLETE** —
  my untracked artifacts are now in the repo history. Intended outcome achieved.
- `ALIGNED` healthy at sample 3/12 (g_tot ~3.5e-3 mN/m, engaged 23959, adv=0 pre-engage,
  drift 3.39e-15); GSD 18.5 MB; no crash; no sync-conflicts. Watching toward transport-
  engage for crash recurrence.
- **iter 45 — 17:15 KST:** ALIGNED still sample 4/12 (log+GSD frozen at 16:53, ~22 min,
  modestly over the ~14 min/sample cadence). adv=0 (pre-engage); no crash; HEAD 7bc2ea7;
  no sync-conflicts. Not a stall escalation yet (no crash sig; could be a slow sample /
  pre-engage dynamics). If still at sample 4 next cycle (~40 min+), give the PI a soft
  heads-up. Hypothesis (alignment lifts γ?) still UNVERIFIED — decisive point (transport
  engage, ~sample 9) not reached.
- **iter 46 — 17:20 KST:** plateau resolved (slow sample) — ALIGNED now sample 6/12
  (s5 g_tot 4.85e-3 soft 5.25e-4; s6 2.24e-3 soft 5.7e-5; s_grip 281→335 nm); GSD 43.4 MB;
  drift stable; no crash; HEAD 7bc2ea7; no sync-conflicts. **adv still 0 (pre-engage)** —
  s_grip still monotone (no engage-point reset yet); γ still floor, soft noisy (no lift).
  Hypothesis still UNVERIFIED. Note: 12-sample run may finish before transport engages →
  alignment test could be inconclusive in this short run. Watching.
- **iter 47 — 17:34 KST:** ALIGNED sample 7/12 (g_tot 3.08e-3, soft 5.3e-5, s_grip
  388 nm, adv=0, drift stable); no crash; HEAD 7bc2ea7; no sync-conflicts. s_grip 388 nm
  is near v300full's engage point (~440 nm, its sample 8→9) ⇒ **sample 8–9 is the
  decisive window** (= both the crash-recurrence risk AND the γ-lift test). Still
  pre-engage / γ floor / hypothesis UNVERIFIED. Watching closely next cycles.

### iter 48 — 17:55 KST — 🔴 CRASH RECURRED on ALIGNED (same point, confirms diagnosis)
- ALIGNED crashed during the sample 8→9 interval: same `FloatingPointError`
  (`constrained_baoab.py:192 _wrap_into_box_xp`), **|round(f)|=1.777e9** (worse than
  v300full's 3.255e8).
- **Reproduced at the SAME point:** both runs had sample 8 at **s_grip=440.4 nm, adv=0**,
  then blew up at/just-after that point. The crash is a **reproducible instability of
  force-scaling-ON + v0_accel=300× at the grip-walk transport-engage transition**
  (s_grip≈440 nm) — NOT a one-off. Confirms the iter-13/14 diagnosis (over-large force at
  engage) across two runs.
- **Hypothesis STILL UNVERIFIED & now BLOCKED:** ALIGNED stayed adv=0 / γ-floor through
  sample 8, then died exactly at the engage transition — so post-engage γ (the alignment
  test) was never measured. The crash is the common blocker for BOTH the γ run and the
  alignment test. Recovery (Lead/PI): force-scaling-OFF gold-standard, or lower v0_accel,
  to survive engage. Data samples 1–8 + 55.9 MB GSD intact.
- **STOP building** (already maintenance); did NOT touch integrator/run/git/Notion. HEAD
  7bc2ea7; no sync-conflicts. Surfaced to PI. (Did not edit the now-committed crash-diag
  note — recorded the recurrence here instead.)
- **iter 49 — 18:15 KST:** ALIGNED crash unchanged (17:55); Lead has not relaunched
  (~20 min, no new run/artifact); HEAD 7bc2ea7; no sync-conflicts. Crash already surfaced
  + recorded (disk + git). No new build. Holding; recovery = Lead/PI. (NOTE for PI: the
  crash findings are on disk + git but NOT yet in Notion — offered to add, awaiting OK.)
- **iter 50+ — ROLLING POST-CRASH WATCH (updated in place):** last check **18:35 KST** ·
  ALIGNED crashed 17:55, Lead not yet relaunched · HEAD 7bc2ea7 · no new run/artifact · no
  sync-conflicts. Crash surfaced + recorded (disk+git). Holding; recovery = Lead/PI. Will
  report verbosely on relaunch / new crash / new commit / PI direction.
  · update **19:15 KST** — still no relaunch (~1 h 20 m post-crash), HEAD 7bc2ea7, no
  sync-conflicts.

### iter 43 — 16:11 KST — ALIGNED sample 1/12 (healthy)
- `ALIGNED` sample 1/12: g_tot 4.46e-3 mN/m, s_grip 60.5 nm, engaged 23851, adv=0
  (pre-engage), drift 3.18e-15 (stable). Shorter run (n-sample=12 vs v300full's 40).
  No crash; HEAD c3fafed; no sync-conflicts. Watching toward the transport-engage point
  (~sample 9) where v300full crashed.
