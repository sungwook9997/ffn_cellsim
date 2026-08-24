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

# H.7 session handoff → next main session (2026-06-07)

Branch `h7/full-cell-integration`. This main session drove the H.7 γ-floor diagnosis + the
mesh (faithful cortex + surface-manifold) coding. Two OTHER sessions ran concurrently (see §5).

## 1. The scientific bottom line — SOLID, triply-confirmed
**The active cortical-tension γ floor is GENERATION-LIMITED** (per-head force / myosin engagement /
buckling), NOT a measurement/formula, connectivity, or density problem.
- `γ_soft` (method-of-planes) ALREADY captures network transmission/amplification — it reads the
  realized deformed-network bond tensions; sitting ~1000× below even the dipole ceiling proves the
  deficit is UPSTREAM (generation), not the formula. The IK whole-shell virial (`gamma_soft_ik`)
  agrees (lower variance). [`H7_SIGMA_A_NETWORK_INVESTIGATION_2026-06-07.md`]
- Correct literature anchoring consistently WIDENS the gap (NOT gate-chasing): density 3.0→0.6/µm²
  (Nie 2015; the "Salbreux 3.0" was a CONFIRMED misattribution) → M1 ceiling 0.10→**0.02 mN/m**;
  EM geometry (heads 56 Niederman, dipole 301 nm Billington). Gap to the MCF7 datum 0.27 (Hosseini,
  ⚠️PROVISIONAL — read off Fig 2g, γ-protocol unconfirmed) is ~13×.
- Two serial gates (the only levers): **Gate A** = myosin heads bind+load but don't WALK (s_grip≈0
  → no contraction → no tension); **Gate B** = M-SHAKE forbids the filament buckling that the
  literature (Ronceray-Broedersz-Lenz 2016) says amplifies active stress ~10×.
- The active-gel seam M1 (`cortex/active_gel_seam.py`) is a DIAGNOSTIC, not a fix; it separates
  timescale-gap from generation-limit and confirmed generation.

## 2. Commits this session (git log, newest first)
```
acbed26 surface-manifold prototype (geometry-only + patch-search benchmark; both gates PASS)
bac4312 faithful_connected_mesh production-reachable (Cell.build + build_baseline_cell)
84accc4 grip_walk variable-length port (faithful M0 mesh + grip_walk now works)
e443680 density re-anchor 3.0→0.6/µm² (Nie 2015, not Salbreux)
7e61ebe GPU-main: device-resident BAOAB into unconstrained build (OPT-IN, default-off)
f248fdf surface-manifold-explicit-cortex DESIGN cross-check (4-agent workflow)
786323b M0: faithful variable-length bimodal connected cortex wired into build_cortex_full_simulation
8abfd0a band re-target → MCF7 0.27 reference (band demoted to rounded-cell envelope)
9a57715 IK whole-shell virial gamma_soft_ik channel
bb7115c bimodal myosin fil_idx bug fix
```
**⚠️ UNCOMMITTED (Bash classifier outage blocked the commit — RETRY FIRST):** the manifold FA-traction
files: `aleph/bridge/manifold_traction.py` + `aleph/scripts/h7_manifold_traction.py` +
`outputs/h7/figs/h7_manifold_traction.png` + `.json`. Contact-conservation gate PASS (residual 0.0).
Commit ONLY these (surgical) — see §5 for the concurrent-session files NOT to touch.

## 3. Mesh coding — DONE (PI directions)
- **faithful M0 mesh** (bimodal variable-length cortex): DONE + production-reachable.
  `build_baseline_cell(faithful_connected_mesh=True)` builds + equilibrates the full cell (grip_walk
  + mesoscale force-scaling + compartments) on the bimodal cortex (Arp2/3 branches, Chugh lengths).
  Default False → uniform fast-hybrid, BIT-IDENTICAL (54 tests green).
- **surface-manifold** (spatial substrate, NOT mesh-as-cortex): `cortex/surface_manifold.py`
  (geometry-only icosphere) + `scripts/manifold_search_benchmark.py`. Broad-phase-NEUTRALITY +
  resolution-invariance gates PASS. HONEST finding: the cortex-search speedup is NOT realized at the
  ×40 scale in CPU-Python (scipy cKDTree already fast) — the manifold's value is the SPATIAL/contact
  infrastructure (shared frames, FA/ECM contact patches, deformable shell, Blender viz), not a cortex
  speedup. [`H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md` + the sub-session's
  `H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md` converge on this.]
- **(b) manifold as FA spatial substrate** — FIRST increment DONE (pending commit, §2): FA
  contact-patch traction FIELD. Tangential 22.6 Pa >> normal 1.94 Pa (basal shear, as expected),
  azimuthal polarity 0.884@47°, conservation gate PASS. Sparse at smoke scale (4 integrin bonds — the
  binding-throughput limit; a full-scale GPU run densifies it).

## 4. NEXT (in priority order)
1. **RETRY the manifold-traction commit** (§2) once Bash works.
2. **(b) continue** wiring the manifold spatial substrate: lamellipodium region-masks on manifold
   (basal_ring/polarized_patch as patch masks + local frames); then cell↔ECM contact manifold.
3. **Gate-A experiment** (does γ rise as s_grip develops?) — driver `scripts/h7_gate_a.py` ready, BUT
   BLOCKED on GPU throughput (<1500 steps/s, the per-step cpu_local_snapshot). Needs the GPU-main
   port (§5) to run the ~2×10⁸-step contraction.
4. **Gate B** (M-SHAKE relax → buckling) — PI PRE-AUTHORIZED a CONTROLLED experiment (keep stretch
   stability, allow bending; measure r/r0<1, Euler margin, drift, γ above the dipole ceiling).
   Needs the GPU throughput too. Integrator-freeze; the sign-off is for the experiment only.
5. faithful-mesh full-cell COMPARISON VIZ (existing uniform vs new bimodal — PI's "how it differs"
   request) — option (a), not yet done.
6. fa.py force-free-integrin fix (spec'd `H7_FA_INTEGRIN_OVERLOAD_FIX`, not implemented).
7. Hosseini-2020 0.27 γ cross-check (currently provisional, read off a figure).

## 5. Concurrent sessions — DO NOT COLLIDE (file ownership)
- **GPU-opt session** owns + is editing: `integrator/baoab_device.py`, `scripts/h7_gate_a.py`
  (added `SnapshotFreeMyosinGripClock` — a SPEED diagnostic, NOT a valid production myosin: it doesn't
  retarget bonds or do Bell-Evans, so its γ is not a production value), `tests/test_myosin.py`, and
  is ADDING `scripts/h7_force_stack_profile.py`, `scripts/h7_hotloop_scaling.py`, `outputs/cupy_cache/`,
  + GPU docs (`H7_GPU_ENGINE_OPT_FOR_CODEX`, `HOOMD_GPU_BESTPRACTICE_GAP`, `NATIVE_HOT_LOOP_MIGRATION`).
  Do NOT commit/edit these — they're the GPU session's.
- **Sub-session** is docs-only: `AGENTS.md`, `H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`,
  `H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07.md`. ∅ runtime intersection.
- **GPU-main port state**: the constrained BAOAB is GPU-ported+validated (Phase-1); the unconstrained
  device BAOAB is wired OPT-IN (`FFN_GPU_DEVICE_BAOAB=1`, default-off = zero regression) but its
  full-system validation is BLOCKED on gbook launch + cupy-compile flakiness (a Codex consult was
  prepared; the proven launch pattern is a `.sh` FILE launcher with conda INSIDE the setsid, not
  inline `bash -c` with env-var prefixes). This is THE blocker for long γ production.

## 6. Hard-won lessons / tripwires
- gbook detached launches: use a `.sh` file launcher (`setsid bash file.sh`) with `source conda;
  conda activate` INSIDE; inline `setsid bash -c "...FFN_X=1 python -u..."` with env-var + quoting
  silently fails to create the log. Foreground `... | grep` shows nothing (python block-buffers
  stdout when piped + timeout kills before flush) — use `-u` + detach-to-file + read the file.
- Surgical commits only (concurrent shared tree): `git add <explicit paths>`, never `-A`.
- `git commit -F file` (never heredoc through the persistent shell). Co-author line:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- The GPU-env test `test_ku35_patch_contractility.py::test_motors_off_baseline_has_no_myosin` FAILS on
  a CPU-only Mac (requests `--device gpu`) — that's the env, not a regression.
</content>
</invoke>
