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

# Session handoff — 2026-07-23 (overnight run → fresh session)

Branch `codex/ff-ac-codex`. The overnight autonomous run is COMPLETE and consolidated; this hands off to a
fresh Lead session (the prior context filled). Everything is committed + verified. Boot from this file.

## Boot prompt (paste into the fresh session)

```
너는 이 세션의 Lead다. ffn_cellsim, 브랜치 codex/ff-ac-codex, HEAD 최신. 부트:
1. CLAUDE.md 정독.
2. docs/v2_audit/SESSION_HANDOFF_2026-07-23.md (이 파일) + AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md +
   OVERNIGHT_ORCHESTRATION_STATE_2026-07-23.md(§MORNING) 정독 — 상태·다음 방향·PI 결정지 전부 여기 있음.
3. git log --oneline -10; conda activate ffn_sim; CUDA Warp only 확인.
4. 다음 액션(아래 §Forward direction)에서 이어간다. gbook 실행은 항상
   PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/aleph.
현재 상태 한 줄: ac/engine 13컴포넌트/32커넥터 composed world 등록·dispatch 완료, A5000 CUDA 256/0,
resting 게이트 해법 메커니즘(bound-myosin) 실증됨 — PI-GAP(fraction/force)만 대기. 완성은 각 컴파트먼트
NATIVE + whole-cell R8.
```

## Where things stand (all committed, verified)
- ⭐ **Resting gate = SOLVED as a mechanism**, blocked on one PI-GAP. `9bb91e0a` (`ac/motor/resting_setpoint.py`)
  + `RESTING_BASELINE_DIAGNOSIS_2026-07-23d.md`. Bound-myosin cortical tension + force-transmitting ERM →
  membrane residual 40.74→2.17 pN (95%), minimum EXACTLY at the Laplace target. Opt-in, default-OFF. Activate:
  `--resting-bound-myosin-fraction/force/source`. Diagnosis chain 22c/d/e + 23a/b/c/d all committed + honest.
- **Composed world REGISTERED + dispatched** (13 comp / 32 conn; `composition.py`, `43cee611`), global
  `GlobalCellLedger` Newton-3rd gate (`234963f5`), `dump_composed_world_state` (unique filament IDs). Per-
  compartment on/off viz + composed-world render. **A5000 `tests/ac/engine` 256/0.** 8 components authored to
  honest CUDA_UNIT/SEAMED states (adversarial-verify corrected 3 overstated claims).
- **Parameters**: P2 cards applied (collagen conc-gate, Hosseini γ KB-6.1.7, ERM HOLD; kb-check green).

## Forward direction (from `AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md`)
1. ⭐ **Close the resting gate** — the single highest-leverage unblock. Either the PI supplies a sourced
   `resting_bound_myosin_fraction + per-head force`, OR a session SOURCES them (NMII resting duty ratio +
   isometric per-head tension from literature/KB — surface if unfound, do NOT invent). Their product must land
   `γ_cortex = ΔP·R/2 − γ_mem` (140 pN/µm at 40 Pa). Activate via CLI → native `<0.21` gate → R2.
2. **Phase B (advanceable NOW, no PI, parallel to CUDA_UNIT)**: cytosol/nucleus CONNECTED binding · SF disjoint
   population build · nmii per-head Bell KMC · lamellipodium Arp2/3 net · ECM GPU SoA topology (merge
   `codex/ecm-gpu-topology` +3) · MT/IF/filopodium kernels with rate/stiffness as UNSET PI-GAP slots.
3. **Phase C (post-resting, serial on one A5000)**: R2 membrane⊕cortex⊕cytosol → R3 ecm → R4 SF-FA-ECM+nmii →
   R5 MT/IF/LINC/nucleus → R6 protrusion → R8 whole-cell NATIVE = 완성.

## PI decisions pending (8 GAPs — surfaced, no defaults invented)
⭐ resting bound-myosin fraction+force · MT DI rates · IF WLC EA/x_max · filopodium k_fascin · collagen gate
(P2 ratification near done) · MCF7 ERM density · nucleus I0-B2 (5) · NMII backbone L_p. Full table + unblock
mapping in the completion roadmap §2.

## Operational notes
- gbook A5000 for native (serial, one GPU); dev Mac = source + CPU oracles + structural tests. Always
  `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim` (stale `~/ffn_cellsim` editable shadows the package).
- Shared branch (Codex bundling hazard): `git diff --stat` + EXPLICIT file staging (never `git add -A`) before
  every commit; verify md5 on gbook (Syncthing/Codex can revert). Commits end
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. No ffn/foundation push (PI-gate).
- No magic numbers, no gate-loosening, physiological baseline, validate at full native. Adversarial-verify
  every authored increment before integrating (overnight verify caught 3 overstated evidence states).
- A Codex read-only adversarial review of the overnight work was requested (prompt handed to PI); when its
  findings land, VERIFY each against code (don't rubber-stamp) before acting.
