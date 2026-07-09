# Cadherin 성숙 + MCF7 σ 재보정 — 구현·검증됐으나 aggregate-σ 압축을 늦추지 못함 (2026-07-09)

**Engine:** DCM (Warp, A5000)  **Branch:** dcm/main  **Commits:** 8ec79ce (구현)
**PI 지시:** 실제 spheroid 비교 후 — "재보정하면서 cadherin 관련도 바로 넣어서 진행."

## 목표

실제 MCF7 spheroid 압축은 24–48h인데 모델은 ~100s. 문헌 비교(KB)로 두 간극을 짚었다: (1) 내 σ=5 mN/m은
embryonic Foty(KB-5.13)지 MCF7(**0.15–0.75 mN/m**, KB-PIV-9)이 아님, (2) 실제 타임스케일은 느린 cadherin
재배열(성숙, junction lifetime 5–30분 KB-4.11)이 rate-limit. → σ 재보정 + 성숙 kinetics를 넣어 진행.

## 한 것

1. **σ 재보정 노출**: MCF7 σ (0.15–0.75 mN/m) 사용 가능 (driver `SIGMA`). MCF7은 Phase B σ-sweep의 **weak
   regime**(σ=1이 −6.4%)이라 압축이 약함 — embryonic σ=5(−8.7%)는 MCF7 과대평가였다.
2. **junction 성숙 구현** (`dcm_cadherin_gpu` + host + driver, `--cad-mature`): bond가 AGE를 ping-pong으로
   운반, 성숙분율 m=1−exp(−age/τ_mature)가 off-rate를 nascent catch-slip → mature floor(k_off=1/lifetime,
   5–30분 KB-4.11)로 blend. 힘 의존성 유지. 기본 OFF(back-compat, host update() 경로·gate5 테스트 무손상).

## 3개 결정적 테스트 — 성숙은 작동하나 실제 타임스케일을 못 만든다

| 테스트 | 설정 | 결과 | 의미 |
|---|---|---|---|
| **① 성숙 작동?** | 정적 confluent, τ_mature=0.01 (즉시성숙) | cum_broken **8259→70** (junction 동결) | ✅ 메커니즘 정상(버그 아님) |
| **② lit τ에서 engage?** | 정적, τ_mature=10s (life 600) | cum_broken 8259→8332 (불변) | ❌ **성숙 전 파열** |
| **③ 압축이 cadherin-gate?** | loose 압축, junction 동결(τ=0.01) | porosity 0.236→**0.234** (불변) | ❌ **σ 압축은 cadherin 무관** |

**② 원인 — junction off-rate가 single-molecule.** node-pair junction은 bundle(cadherin ~10개)인데
**FORCE만 ×bundle_n 스케일이고 BREAK는 single-molecule catch-slip**으로 파열 → 유효 junction lifetime
~0.74s, 성숙 τ(5–30분)보다 ~800× 빠름. bond가 성숙하기 전에 파열되어 성숙이 foothold를 못 얻는다. 실제
junction(bundle)은 훨씬 안정(5–30분). 이건 기존 모델의 bundle 불일치.

**③ 원인 — aggregate-σ 압축은 drag-limited 반경방향 densification이지 cadherin-gated T1 재배열이 아니다.**
junction을 얼려도(파열 72× 감소) 압축이 그대로 0.234 → σ 힘이 cell을 안쪽으로 당기는 이동은 cytoplasm
drag(η=65.9)만이 저항하고, cadherin bond는 soft 인력 tether라 반경방향 이동을 막지 못한다. 성숙(=T1
이웃-교환을 gate)은 T1-limited 과정만 늦출 수 있는데, aggregate-σ 압축은 T1-limited가 아니다.

## 정직한 결론

**"cadherin 성숙만 추가"로는 실제 24–48h 타임스케일이 안 나온다.** 실제 spheroid의 느림은 cadherin이
gate하는 **느린 T1 재배열**이 rate-limit이기 때문인데, 모델의 aggregate-σ 압축은 그 재배열을 우회하는
**liquid-drop 반경방향 pull(drag-limited, 빠름)**이라 구조적으로 cadherin을 rate-limiter로 쓰지 못한다.
실제 타임스케일을 emergent하게 만들려면 **두 가지 재설계**가 필요:

1. **junction(bundle) collective off-rate** — node-pair를 single-molecule이 아니라 bundle_n load-sharing
   집합 rate로 파열(→ junction lifetime 5–30분 → 성숙이 engage). bundle FORCE만 스케일하던 기존 불일치 수정.
2. **rearrangement-driven 압축** — aggregate-σ 전역 pull 대신 cadherin junction 장력이 cell **T1 재배열**을
   구동(그때 junction 안정성이 rate-limit). 단 [[project-dcm-compaction-turgor-blocked]] §2e: junction lever
   단독은 전역 densify를 못 함(그래서 aggregate-σ가 도입됨) — 이 둘을 잇는 게 진짜 open 문제.

즉 PI 직관(실제 타임스케일 = cadherin)은 실제 생물학엔 옳지만, **모델의 aggregate-σ 프레임이 그걸 쓸 수
없다.** aggregate-σ는 끝상태(dense 구)와 driver(Foty σ)는 맞히지만 **rate(타임스케일)를 잘못된 이유로 빠르게**
낸다. 이건 [[reference-dcm-meters-ff-microns]]류 파라미터 슬립이 아니라 **압축 메커니즘 자체의 표현 격차**다.

## Files / 테스트
- `dcm/dcm_cadherin_gpu.py` (cad_break/form + age·maturation), `dcm/dcm_cadherin_host.py` (params, age
  ping-pong, update_gpu), `dcm/dcm_warp_decohesion.py` (--cad-mature/-tau-mature/-mature-lifetime).
- 로그: `~/ff_scratch/_prod_out/matsm_{0,1}.log`(동일 압축), `matstat_{0,1}.log`(정적 불변),
  `matbug.log`(τ=0.01 동결), `frzcomp.log`(얼린-junction 압축 0.234).

⚠️ host update() 경로는 아직 성숙 미구현(GPU 경로만) — cpu 성숙 parity는 후속. Related:
DCM_AGGREGATE_COMPACTION_LARGEDT_2026-07-08, memory [[project-dcm-compaction-turgor-blocked]].
