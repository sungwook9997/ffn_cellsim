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

# Session handoff — 2026-07-23b (cortex structural fixes landed; Whole-Cell Event Runtime ratified → GATE A next)

Branch `codex/ff-ac-codex`, HEAD `2af39966` (latest). This hands off to a **fresh Lead session** for the next
phase. Everything is committed + native-verified. **Boot from this file + the ratified architecture doc.**

## Boot prompt (paste into the fresh session)

너는 이 세션의 Lead다. ffn_cellsim, 브랜치 codex/ff-ac-codex, HEAD 최신. 부트:
1. CLAUDE.md 정독.
2. `docs/v2_audit/WHOLE_CELL_EVENT_RUNTIME_2026-07-23.md` (⭐PI-RATIFIED 최종 방향) +
   `SESSION_HANDOFF_2026-07-23b.md`(이 파일) + `project-ac-cortex-structural-fix`/`project-whole-cell-event-runtime`
   메모리 정독 — 상태·다음 작업·게이트 전부 여기.
3. `git log --oneline -12`; `conda activate ffn_sim`; CUDA Warp only 확인 (Mac=CPU 오라클/구조테스트, native=gbook A5000).
4. **GATE A (cortex first vertical slice)**에서 이어간다: 아래 §다음 작업.
   gbook native는 항상 `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim`. 공유 브랜치 = 명시적 staging + md5.

현재 상태 한 줄: cortex 5-결함 중 2개(density 1→20 `bc5ff3b0`, overlap_free 기본ON `81805535`) native-verified
landing → cortex는 이제 **단일 spanning·rigid·0.2µm 3D·무관통 STATIC CONTROL**. resting residual 2634→0.78(진단
baseline)로 환원. 방향은 Whole-Cell Event Runtime으로 **PI 확정**.

## 확정 방향 (요약)
CellState(느린 맥락, force 소유 X) → conditional component event rate → explicit state change. 공유 spine
(accepted-step 트랜잭션 · population/ID ledger · conserved-pool connector[G-actin=Cytosol 소유] · scheduler 계약
propose/snapshot/rollback/commit/ledger). inter-comp event=connector 소유 = ac/engine 그래프 + event 척추.
**GATE A**(static baseline 먼저) · **GATE B**(트랜잭션=공유 인프라 먼저) 협상 불가. ensemble 검증 · sourced-only
CellState · 신중한 KMC 방법. cortex부터 컴포넌트별. ac/weave under ac/engine (feature-frozen ac/cell 신규 biology 금지).

## §다음 작업 — cortex GATE A (수렴된 static resting baseline 확보)
목표: **첫 수렴된 힘-평형 native 세포**(static connected control) = 고정 레퍼런스. "static control 수렴"까지만 주장.
- **P0.2** 강성 operator용 inner solver 재튜닝. Stage-D에서 연결+두께 cortex는 residual_start 0.7766까지 내려가나
  **강성 operator라 solver가 overshoot**(candidate 0.78→2.79 reject→rollback). step/line-search/n_inner를 강성
  operator에 맞춰 재튜닝해 0.78→0.21 게이트로. (진단 22c-23d 연장이지만 이제 well-conditioned operator 위에서.)
  native: `driver --from-resting`(기본 density=20+overlap_free) `--membrane-subdivisions 6`.
- **P0.3** resting cortical-tension source: `resting_bound_myosin`(ac/motor/resting_setpoint.py)을 **연결된 cortex**
  위에 시딩 → 전달되는지(이전 fragmented에선 inert였음) + 0.78→0.21. PI-GAP: fraction+per-head force(NM2B-pure,
  `params_i0b3.yaml`). ⚠️ head-on-actin(defect#3, capture 0.05µm서 17%만 시딩)이 아직이라 full 시딩 안 되면 그것부터
  (RESTING_SETPOINT_SOURCING_2026-07-23.md §12 + WHOLE_CELL 계획 P4/head-on-actin).
- **병행(코딩 아님):** whole-cell **공통 계약** 스펙 작성 — accepted-step 트랜잭션 추상화(GATE B) + scheduler 인터페이스
  + population/ID ledger + conserved-pool(G-actin) connector 규약. **biology phase 열기 전 PI 보고.**

## 금지/주의
density=20은 STATIC CONTROL이지 최종 cortex 아님(과대주장 금지). gate 맞춤 튜닝 금지 · CellState→force 금지 ·
population 재샘플 금지 · host-authoritative queue 금지 · coarse population 결론 금지 · gate 완화 금지. native
full-70,686 + (확률적이면) ensemble 검증. 커밋 끝: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`; no
ffn/foundation push(PI-gate).

## 이번 세션 커밋 (참고)
`233ff9c1`(감사) `bc5ff3b0`(density) `81805535`(두께) `5e6a8102`/`abefed6c`/`f1f7615c`(계획 진화)
`02111f78`(whole-cell 아키텍처) `2af39966`(RATIFY). resting 소싱: `82c263fa`/`fc84727d`/`8226f54c`.
native probe JSON: `outputs/ac/resting_native/`.
