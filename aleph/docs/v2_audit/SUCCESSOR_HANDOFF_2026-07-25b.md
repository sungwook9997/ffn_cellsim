# ffn_cellsim — 인계 프롬프트 (Opus 5 세션 #1 → 다음 세션)

당신은 `ffn_cellsim`(NVIDIA Warp/CUDA 단일세포 mechanobiology 시뮬레이터)의 Lead Claude Code 세션입니다.
PI는 Sungwook, **한국어로 대화**(코드/커밋/문서는 영어 관례). 전임 세션(Opus 5 #1)이 인계합니다.

## 0. 부팅 (minute 0)

1. `CLAUDE.md` 전체 정독 — 아키텍처 원칙(fine-grained/mechanistic > lumped), I0-A 계약(Warp-CUDA가 유일 런타임,
   HOOMD 금지), 하드룰(magic-number·gate-loosening 금지 / 생리 baseline / NATIVE+FULL 검증 / co-location≠connection).
2. `git log --oneline -30` — 브랜치 `codex/ff-ac-codex`, HEAD 근처 `0ca14341`(viz+closeout) / `fd482920`
   (nmii_sf_motor 런타임). **`ffn/foundation` push 절대 금지**(PI 승인 사항).
3. auto-memory `~/.claude/projects/-Users-sw1-ffn-cellsim/memory/MEMORY.md` 인덱스 훑기. 특히
   `project-ac-nmii-sf-motor-landed.md`(직전 작업) + `[[feedback-*]]` 피드백 전부.
4. **⭐ 세 감사 결과를 먼저 읽으세요** — 전임 세션이 프로젝트 전체를 35-에이전트로 감사했고 결과를 아래에
   저장했습니다. 이게 현재 상태의 가장 정확한 스냅샷입니다(SoT 문서보다 최신):
   - `docs/v2_audit/AUDIT_AC_ENGINE_2026-07-25.md` — ac/ 엔진 6차원(컴포넌트 소유권 / 커넥터 실체 / native
     증거 / 하드룰 / frozen·드리프트 / 파라미터 GAP), 각 차원 적대적 반증 후 종합
   - `docs/v2_audit/AUDIT_WHOLE_REPO_2026-07-25.md` — ac/ 밖 전부(ff/, dcm/, KB 시스템, repo 위생,
     테스트·오라클, 문서·산출물, 엔진 간 정합성)
   - `docs/v2_audit/TRAJECTORY_HOW_WE_GOT_TO_AC_2026-07-25.md` — v1 ActiveCellSim → ffn_cellsim → ff+dcm
     → ac 결정 궤적. 각 전환이 실제로 강제된 것인지(측정된 실패 / 구조적 불가능 / PI 판단 / 근거 복원 불가)
     분류하고, **ac가 앞선 시대를 끝낸 벽들을 실제로 해결하는지** 판정
   ⚠️ 위 세 파일이 없으면 전임 세션이 워크플로 결과 회수 전에 종료된 것입니다. 그 경우 워크플로 transcript를
   직접 읽으세요(같은 세션이 아니라 resume은 불가, 완료된 결과는 파일에 남아 있음):
   `~/.claude/projects/-Users-sw1-ffn-cellsim/db982de2-47e9-4da5-a475-fb7e4bc840d9/subagents/workflows/`
   아래 `wf_ec9ae307-58c`(ac) / `wf_10f8ea01-cd5`(whole repo) / `wf_946791cc-366`(궤적) —
   각 디렉토리의 `journal.jsonl`이 에이전트별 반환값을, `agent-*.jsonl`이 전문을 담고 있습니다.
5. SoT 문서: `docs/v2_audit/cell_engine/ROLLING_ROADMAP.md` +
   `outputs/ac/gate_b_sf_motor/REPORT.md`(직전 작업 클로즈아웃) + `NATIVE_COMPOSITION_SCOPE_2026-07-25.md` +
   `INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md` + `PI_GAP_EVIDENCE_CARDS_2026-07-25.md`.
   Notion "Full-cell breadth" milestone(page 3a7120da…4deb) + 그 하위 day-log.
6. `conda activate ffn_sim`. gbook(RTX A5000)이 유일 CUDA GPU — **UP**(Tailscale, 16GB 전부 여유).
   Mac은 CPU 전용(커널 SOURCE + NumPy reference). 배포 `rsync -az --relative <path> gbook:~/ffn_ac_native/`,
   실행 `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python …`.

## 1. 전임 세션(Opus 5 #1)이 한 것

**`nmii_sf_motor` MOTOR 엣지에 REAL 런타임 착지** (`fd482920` + `0ca14341`). R0부터 그래프에 선언만 되어
있고 런타임이 없던 엣지입니다(`sf_mechanics.SF_CONNECTOR_BINDING_STATUS`가 "left honestly SEAMED"로 기록).

- `aleph/engine/sf_nmii_population.py`(신규) — SF sarcomere에 explicit head-resolved bipolar NMII 미니필라멘트를
  **straddle 배치**. `ac/motor` 원시함수(`straddle_frame` / `placed_positions`) 그대로 재사용, `ac/cell` import
  없음, `warp` import 없음(CPU-importable).
- `aleph/engine/sf_motor_slice.py`(신규) — SF bind-target port + `nmii` owner(자기 파티클 배열) + connector를
  **하나의 CellTransaction / 참여자 3개**로. `sf_arc`가 cortex 슬라이스와 달리 **port가 아니라 진짜 owner**이고,
  split ownership이 **물리적**(두 배열 절대 병합 안 됨 — cortex lane은 `cell.pos_d`를 alias).
- `aleph/engine/cortex_motor_slice.py` — `CortexMotorConnector` → `FilamentMotorConnector`로 타깃-일반화
  (`name`/`component_b` 파라미터화, cortex 기본값 byte-identical, alias 유지). 잘못된 타깃 port를 주면 raise.
- `aleph/engine/sf_population.py` — sarcomere **straddle 기하**(flag-gated, 생략 시 legacy byte-identical). 두
  길이 모두 미니필라멘트에서 **유도**(`lateral = 2·head_offset`, `overlap = backbone contour`).
- `aleph/engine/sf_mechanics.py` — per-segment identity/polarity/arclength를 링크를 emit하는 **같은 루프**에서
  내보내 port가 세그먼트 순서에서 이탈할 수 없게.

**native 게이트 (gbook A5000, `outputs/ac/gate_b_sf_motor/`)**
- **MECHANISM PASS 8/8**: heads `k_on` Poisson 이벤트로 0 → 317/320, 무결합 구간 장력 정확히 0, 장력이 단일
  head 하중 초과(섬유 따라 전달), **FA 앵커 반작용 INWARD 100%**(수축성), clock 1/step, rejected step이
  binding SoA + **양쪽 position 배열** bit-restore.
- 🔴 **QUANTITATIVE BLOCKED — 크기 인용 금지**: explicit relax 미수렴, free-node residual이 보고 장력의
  **15.5%**. 결정적 대조: relax 10× → `T_max`가 **3.4× 변함**(0.451 → 0.133 pN, bound 수 동일).
  게이트가 이제 매 스텝 residual을 출력하고 MECHANISM/QUANTITATIVE 주장을 분리하며 1% 미달 시 정량 주장 거부.

**기하 (측정, host)**: head→필라멘트선 residual **0.200 → 0.000 µm**(capture 0.210 여유 5% → 100%), straight
class bipolar dot **−1.0000**(정확), fiber 배정 1.20 → 1.00. **곡선 perinuclear cap은 구조적으로 제외**
(dot −0.600, scale-invariant — arch는 rigid bipolar 모터가 straddle 불가) → **cap 배치는 PI 결정**으로 상신.

**자체 오류 3건 정정** (전부 측정으로): ① 초기 "T 0 → 1.92 pN VERDICT PASS" → 미수렴 transient였고 **철회**.
② "legacy heads가 actin에 안 닿는다"는 첫 진단 → **세그먼트 기준 측정으로 반증**(닿음). ③ lateral 방향이
"free gauge"라는 주장 → 측정과 모순(connector endpoint를 움직이고 곡선 cap의 load path를 바꿈).

**독립 adversarial 검증**(refute 지시) 12건 지적 전부 조치. 부호/인덱스/stale 버그 없음. 특히 제 테스트
threshold가 cap 측정값(−0.600) **바로 바깥**에 있었던 것 = 결과에 맞춘 gate → 기하적 항등식 dot=−1로 교체.

1078 passed / 70 skipped (was 915) · `make kb-check` runs OK / params OK, 드리프트 없음 · viz 2 scene
browser-verified + 육안 확인(두 antiparallel 필라멘트 사이에 backbone, arms가 양쪽에 닿고, barbed 외향).

## 2. 🔴 최우선 자율 작업 — SF motor 슬라이스에 implicit/CG inner solve

**SF 정량 주장 전부의 유일한 blocker입니다. 메커니즘은 끝났고 solver만 남았습니다.**

문제는 구조적입니다: SOURCED α-actinin dorsal↔arc 가교(4.6e5 pN/µm)가 **여러 개 한 arc apex 노드에 몰려**
λ_max ≈ 3.7e6 pN/µm를 지배 → CFL-stable explicit step이 ~1.4e-7 µm/pN인데 SF 축 모드는 그보다 수십~수백배
느리게 완화. relax를 늘리는 건 답이 아닙니다(10× 늘려도 residual 55% → 17%).

**해법은 in-tree에 이미 있습니다** — GATE A에서 cortex가 같은 벽을 만났고 `ProjectedAnalyticCG` +
fiber-arclength multigrid로 뚫었습니다(`docs/v2_audit/FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md`, 커밋
`09c14eb8`). 그걸 SF 슬라이스에 붙이는 작업입니다. 시작점:
- `aleph/scripts/ac_gate_b_sf_motor_native.py`의 `inner_solve()`(현재 explicit axpy 루프)
- `aleph/engine/sf_mechanics.py`가 SF 강성/토폴로지를 이미 SoA로 갖고 있음 → 연산자 조립 가능
- pinned FA 앵커의 Dirichlet 처리가 이미 있음(`_pin_force_kernel`) — CG의 제약 처리로 옮겨야 함
- 검증: residual/T < 1%가 되면 게이트가 QUANTITATIVE PASS를 내도록 이미 배선돼 있음. 그리고 **relax A/B로
  T가 relax에 무관해지는지**(수렴 증거)를 반드시 확인 — 이번 세션의 핵심 교훈.

α-actinin 집중 자체도 별개 이슈입니다(한 apex에 4개 joint = 이산화 아티팩트일 가능성). solver를 붙이기 전에
그 기하를 먼저 볼지는 판단 사항 — 다만 **가교 강성을 낮추는 것은 금지**(SOURCED 값, Ferrer 2008).

## 3. 🔴 PI-DECISION QUEUE (전부 상신·문서화됨, 사인오프 대기)

**신규 2건 (Opus 5 #1)**
1. **SF motor 슬라이스 implicit/CG** — §2. 자율 진행 가능하지만 solver 선택은 설계 결정이라 PI가 방향을
   지정하고 싶을 수 있음. 착수 전 한 줄 확인 권장.
2. **곡선 sarcomere(perinuclear cap) 모터 배치 설계** — cap은 현재 모터 없음. 근사하지 않고 비워뒀습니다.

**기존 6건 (변동 없음)**
3. coupling-ownership transition **Card-5**: interior fluid column을 진짜 CONNECTED로. frozen
   `aleph/components/incumbent/driver.py` 계약 변경 = strangler 사인오프 필요. plan `INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md`.
   **자율로 driver 편집 금지.**
4. **α2β1-collagen catch-bond Card A1** — FA clutch KMC blocker, 소싱 필요.
5. **PI-GAP cards 33 슬롯** — NMII 클러스터(N1–N9)가 cortex γ + SF traction 동시 unblock = 최고 leverage.
6. turgor Π₀ Card-3(72Pa) · 7. `radial_span` default-on + window refinement · 8. collagen-modulus + PI-GAP
   Notion SoT 등록.

## 4. 다음 자율 작업 후보 (우선순위)

1. **§2 implicit/CG** — 정량 주장 전부를 unblock. 임계경로.
2. `nmii_sf_motor`를 `build_native_composed_cell_world`에 합치기 — 단 두 MOTOR 엣지 간 **head-exclusivity**
   (한 head는 한 타깃에만 결합)가 먼저 필요. 별도 composition 단계.
3. composed cell 확장(다른 private-array owner/connector). MT/IF/filopodium은 PI-GAP-blocked.
4. **세 감사 결과가 지시하는 것** — §0-4의 감사 문서에 임계경로와 위험 원장이 있습니다. 그것이 이 목록보다
   최신이면 그쪽을 따르세요.

## 5. 작업 방식 — PI가 중시하는 것 (반드시)

- **adversarial verify로 real vs theater 판별.** 매 substantive 작업에서 독립 검증 에이전트가
  "facade/스캐폴드/driver-deferred인가"를 **refute 시도**. 통과한 것만 커밋. 스캐폴드는 정직하게 라벨링.
  structural pass ≠ production. 이번 세션에서 이 검증이 제 gate-fitting을 잡았습니다 — 없으면 못 잡습니다.
- **production/native로 검증, 틀리면 즉시 정정.** 결론을 추론이 아니라 측정으로. 성급한 판단 금지.
  공간상관≠인과. **그리고 시간상관≠인과** — bound가 포화한 뒤에도 T가 오르면 그건 walk일 수도, 미수렴일
  수도 있습니다. residual을 보고하지 않으면 구별 불가.
- **하드룰**: magic-number 금지(derivable/sourced/grid-invariant; GAP은 required-no-default로 wire, fallback
  거부). gate-loosening 금지. PI 실험데이터 fitting 금지(literature-first). 생리 baseline(default-off도
  production선 ON). NATIVE+FULL population(70,686 cortex). Warp-CUDA만. co-location≠connection.
- **PI-gated는 자율로 안 바꾼다**: feature-frozen `ac/cell` 계약, mechanistic 모델 결정, 새 gate/contract,
  파라미터 소싱. 정밀 설계+상신하고 권고 기록. 단 **막힌 결정에 멈추지 말고 진행 가능한 건 진행**.
- **memory는 짧게**(상세는 doc/Notion/TAG). **각 milestone viz**(`browser_check.py`로 실제 렌더 + **육안 확인** —
  grep은 검증 아님). **Codex는 아이디어용, 검증 권위 아님**.
- **ETA 보고**(long job은 wall-clock + 근거). 롱-런은 SoT 신뢰 + 서브에이전트 병렬화.
- Ultracode/Max일 때: substantive는 **Workflow**(fan-out + adversarial verify). 토큰 아끼지 말고 exhaustive.
  scout(읽기/스코핑) → orchestrate 순서.
- 세션 클로즈아웃: 커밋 → Notion 3-store(status board / Open items / milestone day-log) → figures →
  `make kb-check` → 마지막 메시지 끝에 **`Notion 업데이트 완료`**.

## 6. 이번 세션에서 배운 함정 (반복 금지)

- **측정은 실제 결합 대상 기준으로.** head가 actin에 닿는지는 **노드**가 아니라 **세그먼트**까지의 수직거리로
  판단해야 합니다(attach는 세그먼트의 barycentric 점에 붙음). 노드로 재면 오진합니다.
- **verdict threshold를 측정 후에 놓지 마라.** 측정값 바로 바깥의 임계값은 gate-fitting입니다. 임계값은
  기하적 항등식이나 물리적 기준에서 유도하세요(예: straight sarcomere는 dot = −1이 정확).
- **수렴 없이 magnitude 인용 금지.** free-node residual을 매 스텝 보고하고, relax A/B로 값이 relax에
  무관한지 확인. 아니면 그 숫자는 transient입니다.
- **verdict를 MECHANISM / QUANTITATIVE로 분리**하면 정량이 막혔을 때도 메커니즘 성과를 정직하게 주장할 수
  있습니다. no-op 임계값(`> 1e-10`)과 구조적으로 보장된 상관관계는 verdict가 아닙니다.
- 도구 함정: Workflow script는 plain JS이고 **template literal 안에 백틱 금지**(파싱 실패).
  `browser_check.py --scene`은 **정확한 전체 scene 이름** 필요(부분 매칭 안 됨) — 틀리면 조용히 첫 scene을
  렌더하고 OK를 냅니다. viewer 카메라는 세포 스케일 고정이라 0.3 µm 물체는 **등방 배율**을 명시해서 확대해야
  육안 검증이 됩니다.

## 7. 한 줄 요약

cortex mesh 지적에서 셀 전체 통합까지 왔고, 직전 세션은 SF에 **진짜 모터**를 달았습니다(메커니즘 native 확정,
정량은 solver-blocked). 다음은 **그 solver** — 그게 SF 정량 전부를 엽니다. 매 결론을 production/adversarial로
확증하고, 틀리면 즉시 정정하고, PI-gated는 상신하며, 셀을 evidence-ladder 위로 **정직하게** 밀어올릴 것.
그리고 §0-4의 세 감사 문서를 먼저 읽으세요 — 이 프로젝트의 현재 상태와 여기까지 온 이유가 거기 있습니다.
— 전임 세션(Opus 5 #1) 드림.
