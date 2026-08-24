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

# Overnight autonomous orchestration — live state (2026-07-23)

Lead session runs the ac/engine build autonomously to completion (PI /goal 2026-07-23). This file is the
CONTINUITY anchor — if the Lead context is summarized, resume from here. Update it every loop.

## ⭐ MORNING — PI decisions to resume (autonomous work at its productive limit)
The overnight run advanced everything advanceable WITHOUT PI input. To resume progress the PI must decide:
1. ⭐ **resting_bound_myosin_fraction + per-head force** (sourced) → activates the built+demonstrated mechanism (`9bb91e0a`, 40.74→2.17 pN, Laplace-exact) → **closes the resting gate** → unblocks the A5000 CONNECTED chain. RECOMMENDED = candidate ① (already built, opt-in). CLI: `--resting-bound-myosin-fraction/force/source`.
2. PI-GAPs (surfaced, no defaults invented): MT DI v_grow/v_shrink/catastrophe/rescue, filopodium k_fascin, IF WLC EA/x_max.
3. P2 cards review: KB-6.1.7 registered verified/High (may downgrade); run `bash outputs/tag_kb/refresh.sh` (Notion token) to sync duckdb/vault.
4. Small follow-ups: fluid_core CUDA test fixture (CPU→device grid, like the load_path oracle fix); the f_rupt diagnostic comment P2 handed to the ac/cell lane.

## Standing goal (PI /goal)
Run overnight, many parallel lanes, to **full ff-ac/engine completion**. Register to Notion + cross-check the
integration plan every loop (context is filling). Keep extending the plan continuously until the engine is
done. Keep the compartment visualization refreshed in the ff-engine HTML style (force/deformation; remove the
"point" toggle; per-compartment independent on/off as compartments grow). **/loop every 30 min**: check lane +
simulation progress; if a PI decision is pending, DECIDE per the surfaced recommendation and proceed; keep
going. When the immediate to-dos finish, pull the next work from the plan.

## Governing docs (cross-check each loop)
- `AC_ENGINE_INTEGRATION_PLAN_2026-07-22.md` (the plan + launch pack + §12 viz gate)
- `AC_DECISION_CARDS_2026-07-22.md` + `P2_RATIFICATION_PACKAGE_2026-07-22.md` (parameter/gate decisions)
- `cell_engine/ROLLING_ROADMAP.md` (R0–R15+ waves; the 8-state evidence ladder) + `COMPONENT_BUILD_MATRIX.md`
- `RESTING_BASELINE_DIAGNOSIS_2026-07-22{c,d,e}.md` (P0 critical path)

## Lanes running (update status each loop)
| lane | what | agent/task id | status |
|---|---|---|---|
| **P2 apply** | apply ratification package | agent a7123ceb | ✅ DONE — 3 cards committed (Hosseini γ KB-6.1.7 `4e648d9e`; collagen conc-gate KB-1.32 1.5mg/mL=13.14 Pa `51a4e17d`; ERM HOLD `3ef9802f`). kb-check green/no-drift. 40 Pa unchanged. HELD for PI AM: (1) f_rupt diagnostic COMMENT in ac/cell/assemble.py+preload_contract.py → P0 lead applies (out of P2 scope; correction fully in KB-3.B1.4/6); (2) run `refresh.sh` (Notion token) to sync duckdb/vault; (3) KB-6.1.7 registered verified/High + 6 new gates draft — PI may downgrade. |
| **Track-C wave1** | 8 components authored | workflow war265yg9 | ✅ DONE + integrated (`db5bc8f2`) — full engine suite GREEN. Honest states: KERNEL_BOUND(ecm/linc) · struct-bound(surface) · PARTIAL(nmii, MOTOR xbridge unbound) · lamellipodium-KB(protrusion) · MT/IF/fluid SEAMED (verify caught overstated KERNEL_BOUND — A5000 exec pending). 6 INTEGRATION_*.md patch-notes → Lead applies dispatch/world/ac-motor wiring. |
| **Track-C wave2** | shared wiring | agent a1d30fc0 | ✅ DONE (`43cee611`) — ⭐MILESTONE: whole 13-comp/32-conn composed world REGISTERED+dispatched (composition.py builder), nmii MOTOR 2-array adjoint crossbridge closed, dump_state (global unique filament IDs = no double-draw). **engine 244 passed/6 skip**. PI-GAPs surfaced (no defaults): MT DI v_grow/v_shrink, filopodium k_fascin, IF WLC EA/x_max, global boundary-ledger class (design). Handoff: INTEGRATION_LEAD_WAVE1.md. |
| **Track-C wave3** | ledger + composed-world viz | agent ace6aa06 | ✅ DONE — GlobalCellLedger + Newton-3rd balance gate (`234963f5`, 249 passed); REAL 13/32 composed world rendered in per-compartment viz (`94156518`, disjoint filament-IDs asserted, connectors=explicit joints, browser-verified). CPU breadth now largely EXHAUSTED. |
| **P0 wave1** | resting gate attempt 1 | agent a6a11ea6 | ⚠️ NOT closed (native 0.5062, gate 0.21) — honest, no loosening. backbone-aware ERM-augmented block LANDED+verified (FD 1.6e-16). ⭐KEY re-diag (23a): ERM PRELOAD IS COUNTERPRODUCTIVE (0.78 = pure membrane RADIAL turgor; 22e floors were preload artifacts). True blocker = stiff-crosslink/soft-membrane 174× ratio; ALL additive preconditioners exhausted. |
| **P0 wave2 GS** | multiplicative GS (no preload) | `64afff42` | ⚠️ NOT closed — native no-preload 0.7766→gs_tour 0.636/gs_pure 0.709 (best of 3 solver families, still 3× the 0.21 gate). ⇒ **3 solver families exhausted; residual = pure membrane radial turgor.** |
| **P0 wave3 PHYSICS** | form-finding cortex pre-tension | `e30d3257` (23c) | ⚠️ REFUTED at native — cortex pre-tension is POWERLESS on the gate (0.478 w/ vs 0.484 w/o; force-free ERM transmits nothing + sparse crosslinks can't carry balanced hoop tension). ΔP-γ: residual=0.644 matches 40 Pa proxy exactly, NOT 72 Pa. Code UNCHANGED (refuted). |
| **⭐ RESTING GATE = PI DECISION** | model completeness | — | 🔴 **The resting gate is NOT solver/passive-physics; it needs a resting cortical-tension SOURCE + force-transmitting ERM. All 3 solver families + passive pre-strain + preload REFUTED.** PI decides the mechanism: **①resting bound-myosin setpoint (fine-grained, RECOMMENDED)** / ②cortex area-tension shell (lumped-violation) / ③dense ERM. ① needs PI-GAP: resting bound-myosin FRACTION + per-head force. |
| **P0 wave4 mechanism** | resting bound-myosin | agent aa5d90e2 | ✅ DONE (`9bb91e0a`) — ⭐BREAKTHROUGH: mechanism built+DEMONSTRATED. bound NMII heads (existing crossbridge kernel, no new law) + force-transmitting ERM (no preload) → membrane residual **40.74→2.17 pN (95%)**, minimum lands EXACTLY at γ_cortex=ΔP·R/2−γ_mem (Laplace prediction, not fitting). opt-in default-OFF (parity kept). Gate needs ONLY the PI-GAP: **resting_bound_myosin_fraction + per-head force** (product = γ_cortex 140 pN/µm @40Pa). CLI ready. |
| **A5000 CUDA_UNIT verify** | run components' CUDA gates on A5000 | gbook | ✅ DONE — **254 passed / 2 failed / 1 skip**. Most components' CUDA-gated tests PASS on real hardware (→CUDA_UNIT verified). 2 fails: fluid_core (SEAMED — CUDA subslice unrun, EXPECTED/honest) + **load_path (real CUDA bug: device force/kinetics assertion fails — the one real-kernel component)**. |
| **fluid_core CUDA fix** | fix binding CUDA guard | agent acf6838c | ✅ DONE (`06220b2b`) — root cause = fluid_core binding bug (guard read grid.device STRING 'cuda:0' for .is_cuda → false-rejected a real CUDA grid; now checks grid.p array device). **⭐ FULL A5000 tests/ac/engine = 256 passed, 0 failures.** |
| **load_path CUDA fix** | fix device force/kinetics assertion | agent a6b208ef | ✅ DONE (`ec2716d3`) — root cause = TEST ORACLE bug (hardcoded [10,0,0] vs non-axis-aligned geometry; kernel CORRECT, A5000 = host oracle exactly). Replaced magic-number expectations with device-vs-host parity + Newton-3rd closure. **A5000 test_load_path 12 passed.** load_path now genuinely CUDA_UNIT. Remaining: fluid_core CUDA test (fixture provides CPU BiotSubstrate grid — SEAMED, fixture-fix like load_path, deferred to AM). |
| **viz** | per-compartment on/off | agent ae57b2125 | ✅ DONE (`3f7c4ed3`) — point slider removed; grouped layer checkboxes (COMPARTMENTS/CONNECTORS/REFERENCE, connectors per family, multi-on); scene→preset; data-driven (new component auto-toggle); browser-verified 9-comp/14-conn + native 494k, 0 JS errors. |

## Loop protocol (every ~30 min)
1. Check each lane (git log, `/workflows`, gbook sim progress). Integrate verify-PASSED component modules to
   `codex/ff-ac-codex` (Lead applies any shared dispatch/world registration from `INTEGRATION_<c>.md` patch-notes).
2. Self-decide PI-gated items per the surfaced recommendation (8h mandate) — EXCEPT never loosen a gate or
   invent a magic number; those stay documented+surfaced. P2 gate-contract apply is PI-pre-authorized.
3. Persist: update this file + the integration plan (extend it toward completion) + Notion Dev Log.
4. Launch the next wave from the plan/roadmap (deeper bindings, connector families, dispatch/world integration,
   whole-cell composition, re-verify, viz refresh). When P0 closes the gate → start R2 native binding.
5. ScheduleWakeup ~1800 s for the next check.

## Hard rules (unchanged)
Warp-CUDA only (native on gbook A5000, serial; Mac = source + CPU oracles + structural tests). No magic
numbers, no gate-loosening, physiological baseline, validate at full native. Commits end `Co-Authored-By:
Claude Opus 4.8 (1M context) <noreply@anthropic.com>`; no ffn/foundation push. gbook:
`PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim`.

## Completion definition (goal clears when)
All 13 components + 32 connectors reach NATIVE (per ROLLING_ROADMAP), the resting baseline gate is closed,
the whole-cell native baseline (R8) passes with figures, and the plan records it. Until then, keep looping.
