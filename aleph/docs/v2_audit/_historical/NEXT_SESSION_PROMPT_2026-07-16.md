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

# Next-session boot prompt (2026-07-16)

Paste this to boot the next session. State lives on disk + memory + Notion — this points to it.

---

## Where we are (this session: 37 commits, branch `ff/mech-hierarchy`)

**PI directive that drove the session:** "S1-3도 믿지 말고 아예 처음부터 재검증 = FF 엔진을 부품 하나하나 처음부터
재검증하고, 신뢰가능한 기반 위에 계속 쌓아올려라" (autonomous overnight).

**Done + committed:**
1. **FF engine re-validated from scratch — 28/28 PASS** (`scripts/engine_reval.py`, ladder
   `docs/v2_audit/FF_ENGINE_REVALIDATION_LADDER_2026-07-16.md`): T0 integrator → T1 filament → T2
   connectors → T3 cortex → T4 volume → T5 membrane/MT/nucleus → **T6 pre-stressed cell** → T7 missing
   physics. Each part = analytic oracle + phenomenon-completeness + figure; native cell-scale parts on
   the gbook A5000. `python -m aleph.scripts.engine_reval --all` reproduces; `--summary` = one-page status.
2. **⭐ P6.1 trustworthy baseline** = the thing S1-3 lacked: native full cell (`resting_reval.npz`,
   Nc=494,802 cortex+MT + nucleus 3000 + MT 40, ALL compartments ON, sphericity 0.99999, pre-stressed
   ΔP=40Pa/γ=0.15mN/m). Interactive 3D HTML: `outputs/mech_hier/figs/resting_reval_cell.html`.
3. **FEM+CFD infrastructure built + validated + visible:** `ff/biot_fluid_warp.py` (device-resident Biot
   field + IBM spread/interp, parity-gated). T7.F.1 Terzaghi → F.2 3-D field → F.3 IBM (momentum/adjoint,
   1e-16) → F.4 coupled two-way FSI → F.5/F.6 Warp kernels → F.7 native composition. **⭐ CFD-ON cell is
   VISIBLE + VERIFIABLE**: `scripts/ff_cfd_cell_native.py` → `outputs/mech_hier/figs/cfd_cell_native.html`
   (pore-pressure field radial gradient, drains on poroelastic τ_p=0.11s ≈ measured 0.13s, NOT the 0-D
   efflux 30-200s; `cfd_cell_metrics.json`).
4. **Honest re-validation findings** (the value of "don't trust S1-3"): FF = mechanical-equilibrium
   solver NOT thermal MD; myosin = LINEAR not Hill; prior S1 "linear-shell r²=0.992" is REGIME-SPECIFIC
   (small-strain shell 0.986 confirmed, large-strain stiffens); membrane Helfrich bending + IF genuinely
   ABSENT (seeded); 6 bugs caught. Memory: `project-ff-engine-revalidation`.
5. **KB**: 79 SourceEvidence + 73 KnowledgeClaim drafts uploaded (CrossRef-verified, marker
   `AUDIT-2026-07-15`); gates GREEN (verify_runs/params no-drift). Roadmap
   `docs/v2_audit/DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md`.

## Immediate next task — RECOMMENDED: two-way FSI wiring commit (finish the CFD)

The CFD field is real + visible but **one-way** (compression → pressure → drainage). Make it **two-way**
= the fluid feeds back on the mechanics:
- Wire the Biot+IBM FSI into `network_warp.py`'s step as an **additive, default-OFF channel** (mirror the
  `local_load`/`_ll` hook at `network_warp.py:844`) → FSI-OFF == the current cell (the 28 parts + native
  regression PASS unchanged = the safety gate).
- FSI-ON: node divergence → grid source (`ibm_spread`), step `BiotField`, interp ∇p → node Darcy force
  (`ibm_interp`) added to `f_d`. **RETIRE 6πηR → γ_solid** (the double-count trap; RESOLVED:
  `γ_solid = units.fiber_point_drag(eta=ETA_SOLVENT≈1e-3)`, the effective 65.9 Pa·s EMERGES from the
  Biot coupling — see roadmap §"γ_solid RESOLUTION").
- **Emergence gate:** FSI-ON native cell reproduces (a) P4.2 drained↔undrained rate-dependence and (b) the
  spatial poroelastic τ_p — confirming the effective viscosity emerges, not double-imposed. Render the
  FSI-coupled cell as interactive HTML (per the viz rule).
- ⚠ This is the "biggest correctness trap" (η double-count) — do it carefully with the regression gate,
  fresh focus, not at a marathon tail.

## Alternative direction (PI's call): advance S1-S8

S1✅ S2✅ S3✅ (re-validated native) · **S4 🔶 partial** (compartments validated, not the layered
rate-dependent STAGE) · **S5-S8 ❌ not started** (spheroid / contact / active). If the PI wants S-stage
progress over the FSI wiring: build S4 (three-layer distinct rate/duration response) then S5 (spheroid of
validated S1 units) per `docs/v2_audit/HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md`.

## Hard constraints (do not violate)
- **NATIVE + FULL only** for any conclusion (`--cortex-fil 38000 --from-resting --microtubules`, gbook
  A5000; `~/miniconda3/envs/ffn_sim/bin/python`, rsync + nohup + log). Coarse = non-authoritative smoke.
- **Visualize as INTERACTIVE 3D CELL HTML** (compartments/fields visible, cut-away, colorbar) — NOT just
  matplotlib. (This session's miss — the PI had to ask; browser-verify with `scripts/browser_check.py`.)
- **No η double-count** on the FSI wiring; **physiological pre-stressed baseline**; no magic numbers /
  gate-loosening (surface to PI); Korean replies, English code/commits.
- New physics → add a part to `engine_reval.py` (oracle + phenomenon + figure), native-gated.

## Boot steps
1. Read `CLAUDE.md`, this prompt, the ladder + roadmap docs, memory `project-ff-engine-revalidation`.
2. `git log --oneline -5` on `ff/mech-hierarchy`; `conda activate ffn_sim`.
3. `python -m aleph.scripts.engine_reval --summary` (confirm 28/28) + open `cfd_cell_native.html`.
4. Confirm the direction with PI (FSI wiring vs S4-S8), then execute native-gated.
