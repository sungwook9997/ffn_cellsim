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

# 다음 세션 부팅 프롬프트 — ffn_cellsim (2026-05-31 자율세션 인계)

> PI께: 다음 리드 세션 시작 시 아래 블록을 그대로 붙여넣으세요. (이 문서 자체가 인계 프롬프트입니다.)

---

[새 리드 세션 부팅 — ffn_cellsim / KU-3.5 myosin 재설계]

CLAUDE.md 세션 부팅 프로토콜대로 시작해줘. 부팅 입력:

- **상태**: `phase1/h3-cortex` @ `11c554f`, 직전 자율세션 13커밋 전부 회귀-green, clean tree, env `ffn_sim`(HOOMD 7.0.1). cell은 이제 **2.5e6 step crash-free**.
- **부트 시 읽기** (직전 세션 산출, 전부 `aleph/docs/v2_audit/`):
  - `PI_WAKE_BRIEF_2026-05-31.md` — **단일 boot-point + PI 결정큐** (먼저 읽기)
  - `KU35_FLOOR_ROOT_CAUSE_2026-05-31.md` — KU-3.5 floor 진단
  - `KU35_GRIP_WALK_DESIGN_2026-05-31.md` — myosin fix 구현 설계안 (**#1 작업**)
  - `aleph/outputs/h3/production/AUTONOMOUS_LOG_2026-05-31.md` — 전체 trail (iter 0~13)
  - Notion Dev Logs 상단 closeout 노트 + 진행 페이지 `370120daec5d81aaa9a6c4103cb42732`

**핵심 상태 (직전 세션 결과)**:
- 🎯 **몇 달짜리 KU-3.5 floor 완전 규명**: myosin binned-r0 ratchet = **lumped proxy** (actin 물질 미운반 → 수축력 유실, γ_soft ~25,000× 미달). 측정(method-of-planes)은 정확. CLAUDE.md "no lumped mechanisms" 위반.
- **모든 path-블로커 해결·커밋**: G1 FA clutch(9→45) / G2 myosin stepping=아티팩트(실제 step함) / pre-existing HOOMD 제외목록 크래시(degree-cap fix, 2.5e6 crash-free 검증).
- **EXTEND**: 3모듈(substrate/H.8/H.9) 빌드 + H.8 membrane·Track A substrate 통합(default-off). nucleus(H.9) 통합 deferred.

[이번 세션 #1 작업 — myosin grip-walk fix 구현]
1. **먼저 PI 설계 결정 ratify 받기** (KU35_GRIP_WALK_DESIGN의 PI-decision 표):
   - (a) cortex filament **polarity 할당 방식** — cortex actin은 polarity 없음(설계 추천: label-only "minus=bead 0")
   - (b) **bipolar sidedness** binding 규칙
   - (c) **연속 sub-bead `pos_a_end`** (정수-only walk는 ~4e6 step 미만 run에서 여전히 γ≈0 — 결정적)
2. ratify 후 **opt-in `stepping_mode`(default=현행 `binned_r0`)로 구현** — 단일 핵심 edit `cortex/myosin.py:985-1022`(bin-r0 감소 → grip-point walk + r0≈0). re-target primitive는 이미 존재(set_snapshot). **additive default-off**.
3. **Tier-1 micro-diagnostic(분 단위, ~1e5 step)** 으로 sustained-tension-vs-relaxation 검증(40h sim 불필요). 회귀 green 유지.
4. A/B(binned_r0 vs grip_walk) → PI 보고 → PI가 default 전환 ratify(gate-affecting).
5. 그 후 (crash-free) v4 production 재측정 → γ가 [0.35,0.65]로 오르는지 확인. (full production ~40h/constrained dt — dt co-tune/GPU 고려.)

[병렬/후속 (PI 결재 큐 — PI_WAKE_BRIEF 참조)]
- **자율 PI-flag ratify/revert**: `barely→fully` 테스트 반전, HOOMD-용량 상수(`_MAX_CLUTCH_FANIN`/`_MAX_CORTEX_BEAD_DEGREE`/cap-depth), E_sub 라벨, nucleus antipodal-seeding refinement.
- **nucleus(H.9) Template-2 통합**(deferred) + **Track C 복합 KU-3.5/3.1 게이트**(gate-contract → PI).
- 직전 큐(`PI_DECISION_QUEUE_2026-05-30.md`): C1 integrin F*→30pN, A4 KU-3.5 regime, 등.
- **`ffn/foundation` push** (13+커밋, PI 결재).

**가드 (필수)**: cell.py = 단일 라이터(Lead 직렬 통합) / 전부 additive + default-off(켜기 전 회귀 0) / 서브에이전트 git 금지 / `git add`는 파일지정만(-A 금지, 리서치 산출물 보존) / PI-only 게이트(gate-contract·integrator-freeze·magic-number·foundation push, **그리고 myosin core 구현은 설계 ratify 후**)는 PI 결재.

부팅 끝나면 다음 작업 1-2문장으로 복창하고 시작해줘.

---

## (참고) 직전 세션 커밋 13개
```
11c554f wake brief → grip-walk 설계 포인터 + 연속-sub-bead/run-length 발견
d1fb8ac KU-3.5 myosin grip-walk fix: PI 결재 설계 brief (#1)
1af5b4d substrate effective_E_sub 라벨 N/m³ (PI-flag, doc-only)
bbb4a3e Track A substrate 통합 (default-off)
b98807e KU-3.5 finding figures (proxy-failure headline + G1 clutch + A/B)
7ebfc1f PI wake brief: 결정큐 + 핸드오프
206044b H.8 membrane_surface 통합 (default-off) + KU-3.5 floor ROOT CAUSE
f2f4e75 H.3 degree-cap (pre-existing HOOMD 제외목록 크래시 fix)
db1c756 EXTEND 3모듈 scaffold (substrate/membrane_surface/nucleus, default-off)
3d70e6b G1 fix 3: clutch fan-in 2→1 (production 생존)
d33230b G1 fix 2: 남극 cap 분산 (HOOMD 1-3 제외목록 overflow)
c960406 G1 fix: clutch fan-in cap (HOOMD 제외목록 크래시)
75f464e G1: FA contact-footprint seeding (floor blocker a)
```
