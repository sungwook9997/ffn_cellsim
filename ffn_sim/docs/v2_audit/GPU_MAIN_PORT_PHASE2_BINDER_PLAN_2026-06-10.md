# GPU-main port — Phase 2 PLAN: myosin/xlink binder host-sync (2026-06-10)

> Autonomous-loop authored (PI 수면 중). **계획만** — gbook이 dirty 구 브랜치라
> cupy GPU 검증 불가 → 검증 없는 GPU 포팅은 커밋 안 함(verified-only 룰). 이 문서는
> 다음 GPU 세션(gbook 정리 후)이 바로 실행할 수 있는 phased plan.

## 0. 왜 이게 next target인가
`GPU_MAIN_PORT_PHASE1_2026-06-01.md` §4: integrator(constrained BAOAB)는 cupy/native
포팅 완료(native scale 11-12×). **남은 per-step/batch host-sync = 2개:**
1. `MyosinHeadForce.set_forces` (myosin.py:1774) — `cpu_local_snapshot` **매 스텝**
   (continuous_stroke 모드의 production form). per-step sync = 가장 비싼 잔존 병목.
2. `XlinkBondUpdater.act` (crosslinkers.py:946) + `MyosinStepUpdater.act`
   (myosin.py:978) — `get_snapshot`+`scipy.cKDTree`+`set_snapshot` **매 batch(100스텝)**.

## 1. 구조적 제약 (포팅 설계를 좌우)
HOOMD에서 **bond add/remove(topology 변이)는 snapshot 연산 불가피** — 현 API상
device-resident 토폴로지 mutation 경로가 없다. 따라서 binder updater를 *완전*
GPU-resident로 만들 수는 없다. 포팅의 목표는 **host round-trip을 최소화**하는 것:
- 비싼 부분 = (a) 전체 positions host 복사, (b) scipy cKDTree neighbor query.
- 둘 다 topology 변이가 아니다 → **GPU로 옮길 수 있다.** topology delta 적용만 작은
  snapshot 연산으로 남긴다.

## 2. Phased plan (각 단계 additive, CPU bit-invariant, 단계별 검증 게이트)

### P2a — MyosinHeadForce per-step sync 제거 (최우선, 가장 큰 이득)
`set_forces`의 `cpu_local_snapshot`를 device-aware로:
- GPU device → `gpu_local_snapshot` + cupy로 positions/forces 읽고 per-head
  F=min(k·s_grip,F_stall) 벡터 계산을 cupy 커널로(데이터 device-resident 유지).
- CPU device → 기존 numpy 경로 (bit-identical, 기존 13 sanity 테스트가 reference).
- `xp=np` 백엔드 kwarg 패턴(integrator 포팅과 동일)으로 함수 순수화.
- ⚠️ s_grip 누적은 binder updater(batch)가 소유 — set_forces는 읽기만. s_grip
  버퍼를 device-resident로 두면 batch updater와의 인터페이스 정리 필요.
- **검증:** GPU-vs-CPU per-step force match(랜덤 fixture), continuous_stroke 13
  sanity 테스트 GPU에서 재통과, end-to-end g_soft GPU≈CPU(통계적, seed noise 내).

### P2b — cKDTree neighbor query → GPU (topology 변이 분리)
`XlinkBondUpdater.act` / `MyosinStepUpdater.act`의 acceptor 탐색:
- positions를 `gpu_local_snapshot`로 device에서 읽고, candidate acceptor 거리 query를
  cupy(또는 HOOMD nlist) 로 — scipy.cKDTree 대체.
- Bell-Evans break/bind **확률 판정**도 cupy로(RNG는 cupy stream, 통계적 동등).
- **결과 bond-delta(add/remove 목록)만** host로 내려 작은 `set_snapshot` topology
  연산에 적용 — 전체 positions 왕복 제거.
- ⚠️ 메모리: `reference-hoomd-rebuild-leak` — set_snapshot 토폴로지 변이가 rebuild를
  유발하면 leak. delta 적용이 in-place bond array 갱신으로 끝나는지 확인(pre-alloc pool
  패턴, L2.4b 참조). leak 나면 P2b 보류하고 P2a만.
- **검증:** 같은 seed에서 bind/break 결정 분포 GPU≈CPU, n_engaged 궤적 일치,
  degree-budget 캡(≤6) 불변.

### P2c — native-scale 벤치 + Gate 재검증 (PI-gated)
- `gpu_scaling_constrained.py` 패턴으로 binder-on 전체 cortex N-scaling 재측정.
- full-cell native ETA(현 ~11.5d/2e8) 개선폭 정량.
- H.2/H.3 contract gate(tau_min 등) GPU device 재통과 = PI sign-off 필요.

## 3. 검증 인프라 (이미 존재)
- `tests/test_constrained_baoab_gpu.py` = GPU-vs-CPU per-function match 패턴(복제).
- `scripts/gpu_smoke_constrained.py` / `gpu_scaling_constrained.py` = end-to-end +
  N-scaling. binder 버전을 sibling으로 추가.
- CPU 경로 기존 게이트가 항상 reference(additive 원칙).

## 4. 우선순위 + 과학적 맥락
- **P2a > P2b > P2c.** P2a(per-step)가 P2b(매 100스텝 batch)보다 sync 빈도 100× →
  가장 큰 이득. continuous_stroke가 production form이 되면 P2a는 필수.
- ⚠️ **과학적 우선순위 주의(2026-06-10 A/B 결과):** crosslink 1e-3 vs 1e-7 contraction
  A/B에서 γ 차이 ~5%뿐(transmission 벽은 crosslink lever 아님; AB_XLINK 결과 참조).
  즉 binder 강성은 γ를 못 움직인다 — binder GPU 포팅은 **성능 최적화 가치**(native
  long-run 가속)이지 γ-floor 해결책이 아니다. γ magnitude는 density/coherence 문제로
  남는다(PI_DECISION_BRIEF §2).

## 5. 차단자
- gbook 코드 동기화: 현재 phase1/h3-cortex(d9249e5) dirty. GPU 검증 전 PI가 gbook을
  h7/full-cell-integration로 정리(미커밋 작업 처리) 필요 = remote-overwrite, PI-gated.
- integrator/ freeze 불가침(P2는 cortex/ + cell/ 만 건드림).
