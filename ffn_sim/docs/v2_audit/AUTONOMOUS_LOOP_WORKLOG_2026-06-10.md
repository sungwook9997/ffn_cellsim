# Autonomous Loop Worklog — 2026-06-10

> PI가 6h 수면 중 Lead가 10분 self-paced 루프로 진행한 **PI 결정 불필요 unblocked
> 작업**의 durable 기록. core-physics/gate-contract 변경 없음(HARD 룰). 다음 세션/PI
> 부트 시 여기부터 읽으면 됨.

**범위 규칙(자기부과).** ✅ 허용: 분석/진단(샌드박스, production 미변경), GPU 포팅
(additive·CPU bit-invariant·검증 통과 시만 커밋), 문서, KB dry-run. ❌ 금지: 검증 안 된
core-physics 커밋, gate-contract 변경, magic-number, integrator/ 편집, 타 세션 lane.

---

## Iter 1 (시작) — 상태 파악 + PI 브리프 + KB manifest
- 8 워크트리 + Notion/git/RAG·TAG 상태 종합. 결론: 전부 active-myosin force-deficit
  한 벽의 단면, 모두 PI 결정 앞 HALT. 로컬·gbook 유휴(GPU 0%).
- **PI_DECISION_BRIEF_2026-06-10.md 작성** — 3개 대기 결정(crosslink k / myosin density /
  SF line) 정량증거+선택지+추천.
- KB harvest **dry-run** manifest 생성(`OPS_HARVEST_CANDIDATES_2026-06-09.md`, apply 안 함).
- 미커밋 myosin.py(loop24d M-band)·h7 script = 검증 안 됨 → 보존, 커밋 안 함.
- NEXT: task#2 crosslink γ-민감도 진단 샌드박스 + task#3 crosslinker binder GPU 포팅 스코핑.

## Iter 2 — ⭐ 브리프 결정1 정정 + force-budget A/B 플래그
- ⭐ **git 확인 결과 crosslink 재anchor(1e-7→1e-3)는 이미 loop18(236101e, PI-approved)
  완료.** 최초 브리프가 stale 설계doc §10(loop17)을 참조한 오류 → 결정1은 닫힘, "측정만
  남음"으로 정정(PI_DECISION_BRIEF §1 rewrite). cross-bridge 강성 fix(continuous_stroke
  +k=1e-3)도 loop14-17 빌드·검증 완료.
- `h7_active_force_budget.py`에 `--xlink-k` A/B override 추가(additive, sensitivity-only,
  scripts/ lane). CPU smoke 검증 PASS(continuous_stroke+1e-3 빌드·실행, exit 0).
- 단 CPU smoke = loading phase(s_grip≈0). 결정적 transmission A/B(1e-3 vs 1e-7,
  full contraction, s_grip→0.5)는 gbook GPU ~30min/config 필요 → 다음 이터에서 배포·실행.
- NEXT: gbook 코드 동기화 방식 확인 → A/B 2-config 런 launch. 병행 task#3 GPU 포팅 스코핑.

## Iter 3 — ⚠️ gbook 안전판단 + 로컬 CPU A/B launch
- ⚠️ **gbook은 오래된 다른 브랜치(phase1/h3-cortex, d9249e5)에 미커밋 변경(fa.py,
  cell.py, phase1_h3.yaml 등) dirty.** rsync/checkout 덮으면 gbook 미커밋 작업 파괴 →
  되돌리기 어려운 remote-overwrite, **PI 승인 없이 자율 안 함.** 결정적 GPU A/B는 PI-gated.
- 대안: **로컬 CPU 단일빌드 A/B 2-config 병렬 launch**(leak-safe, sweep 아님):
  continuous_stroke, n_fil=120, warmup 2000, contract 150000(≈1500 ticks ≫ Gate-A
  contraction onset 560), sample 15000. A=xlink 1e-3(prod, PID 21409), B=xlink
  1e-7(pre-re-anchor, PID 21573). outputs/h7/production/ab_xlink/. ETA 각 ~25-40min.
- tick0: A g_soft=1.379e-2, B=1.354e-2 mN/m (수축 전 거의 동일 — 예상대로 loading은
  passive prestress 지배). 관전 포인트 = 수축 진행(s_grip↑) 시 actin-network γ(WALL-A)가
  A에서 B보다 오르는가 = crosslink 재anchor의 transmission 기여.
- NEXT: heartbeat로 A/B 폴링 → 완료 시 WALL-A A/B 비교를 브리프/worklog에 기록. 병행 task#3.

## Iter 4 — A/B 예비결과(⭐ crosslink ≠ transmission lever) + GPU 포팅 계획 doc
- ⭐ **A/B contraction 예비결과(plateau, tick 120000/150000):** per-head **meanT가
  F_stall(8.42pN) 도달**(생성 fix 작동 확인), gen force ~3.4nN, 결합 406 heads >
  band-closure 예산 245 → **raw 생성력은 충분**. 그런데 s_grip/l0=0.001 정체(stiff k →
  즉시 stall → 안 걸음, 물리적으로 정확), **g_soft ~1.37e-2 mN/m 평탄**(~13× under
  band). **A(1e-3) vs B(1e-7) 차이 ~5%뿐**(g_ik 5.2e-4 vs 4.95e-4, gen 3.43 vs 3.27nN).
  ⇒ **crosslink 재anchor는 transmission 돌파구 아님.** 벽은 순수 transmission/기하 —
  생성력은 충분한데 isotropic 상쇄로 hoop tension에 안 모임(Gate-A "transmission wall"을
  교정된 작동점에서 재확인). γ magnitude = density/coherence 문제로 남음.
- task#3: `GPU_MAIN_PORT_PHASE2_BINDER_PLAN_2026-06-10.md` 작성 — P2a(MyosinHeadForce
  per-step sync 제거, 최우선) > P2b(cKDTree query GPU화, topology변이 분리) > P2c(native
  벤치+gate, PI-gated). gbook dirty라 GPU 검증 불가 → 계획만, 포팅 커밋 안 함(verified-only).
- NEXT: A/B 최종 JSON(WALL-A myosin/actin 분해+verdict) 캡처 → 브리프 §1에 수치 기록.

## Iter 5 — ⭐ A/B 확정: crosslink는 전파 2.6× 개선하나 γ는 density-bound
- 두 런 완료. **A(1e-3) vs B(1e-7):** WALL-A 전파효율 **4.56% vs 1.72% = 2.6× 개선**
  (재anchor가 transmission 고리를 정말 회복). actin-network γ 6.29e-4 vs 2.38e-4.
  **BUT 총 γ_soft 1.374e-2 vs 1.354e-2 = +1.4%만** — 지배항 myosin-dipole γ=1.378e-2
  (crosslink 무관, 양쪽 동일). per-head meanT=F_stall(8.44pN), 결합 406 heads.
- ⭐ **결론: γ magnitude는 transmission-bound 아님 = density/coherence-bound.** dipole이
  ~13× under(WALL B) = 생성/density envelope 문제. PI 결정 2(A) 데이터 뒷받침. 브리프
  §1/§2/묶음 + 결정1 "닫힘" 갱신.
- figure: outputs/h7/figs/h7_xlink_ab_2026-06-10.png (채널 분해 + WALL-A 비교).
- task#2 DONE. 남은 자율거리: SF 고밀도 테스트(결정3-A 가설검증, 자율가능) OR 워크트리 정리.

## Iter 6 — ⭐ γ-vs-density 곡선: 밴드 도달 ρ≈20/µm² (결정2 정량 확정)
- 3점 완료. **γ_soft: ρ=0.6→0.0137, 6.0→0.0742, 21.0→0.191 mN/m.** power-law
  **γ∝ρ^0.74** → band_lo(0.18) @ **ρ≈19.6/µm²**, band_hi(0.40) @ ρ≈58 = proxy 0.6의
  **~33–96×**. ⇒ γ는 명확히 density-bound, 밴드엔 33–96× density 필요(MCF7 datum 없음).
- ⚠️ mesoscale caveat: density override = force_factor(4.24→148) 스케일이지 head 수 아님
  (native density의 sanctioned ×40 표현; meanT 294pN은 meso-effective). γ∝ρ는 유효.
- PI 브리프 결정2에 곡선/fit/caveat 기록. fig h7_gamma_vs_density_2026-06-10.png. task#5 DONE.
- 누적 결론(iter5+6): **γ magnitude floor = density-bound(NOT transmission, NOT generation
  per-head: F_stall 도달). 밴드 도달엔 ρ~20/µm² 필요, lit proxy는 0.6(MCF7 datum 부재).**
  PI 결정2(A) 데이터 완비. 남은 자율거리: 워크트리 정리 또는 SF 고밀도(결정3) — 다음 이터.

## Iter 7-8 — seed-robustness: slope robust, |γ| ~1.5× scatter
- seed-2 3점 완료. **slope robust(seed1 0.739 / seed2 0.713)이나 절대 γ ~40-50% 낮음**
  (0.6: 0.0137 vs 0.0083; 21: 0.191 vs 0.109). band_lo 도달 density seed1=19.6 vs
  seed2=47.3/µm²(prefactor scatter가 멱승역산으로 ~2.4× 증폭).
- ⇒ 정직 진술: **density-bound 성격은 확실, 필요 density는 정밀값 아닌 "수십/µm²,
  proxy의 ~30-80×" order-of-magnitude.** 2-seed thin → firm 절대값은 더 많은 seed +
  native GPU(PI-gated). loop23 단일-seed 과대주장 방어 = 이번 규율의 핵심.
- 브리프 결정2에 seed caveat 기록, fig 2-seed overlay 갱신. task#6 DONE.
- ✅ γ-floor 자율 특성화 완료(transmission A/B + density 2-seed). 머신 건강 양호(mem 96%
  free, 좀비 없음). 이후: 가벼운 유지 + PI 대기. SF 고밀도(결정3)는 PI-gated 유지.

## Iter 9-10 — ⭐ γ_active 헤드라인 floor 앙상블 + isometric-stall 명확화
- seed 3/4/5 @ d=0.6 완료 → **n=5 앙상블: γ_active = 1.16e-2 ± 0.21e-2 mN/m (CV 18%)
  = 15.5× under band_lo.** per-seed 0.83–1.37e-2. 단일-seed→앙상블 격상(2-seed thin 해소).
  fig h7_gamma_active_ensemble_2026-06-10.png. task#7 DONE.
- 브리프 명확화: phase=loading 라벨은 s_grip<0.02 휴리스틱 탓, 실제는 stiff-bridge가
  ~4nm에서 즉시 F_stall stall(g_soft tick0→120000 평탄 = 진짜 steady state, 아닌
  under-equilibration). continuous_stroke = isometric stall 작동점.
- ✅✅ **γ-floor 자율 특성화 최종 완료.** 결론: active-γ floor = lit density에서 ~15×
  under, **density/coherence-bound**(transmission 풀림+per-head F_stall), 밴드엔 ρ~수십/µm²
  필요(MCF7 datum 없음). PI 결정1=닫힘, 결정2=데이터 완비, 결정3=PI-gated.

## Iter 11-12 — health-check + KB 대조 (citation-grounding)
- Iter11 health-check: 전부 안정(git clean, 0 procs, mem 96%, gbook 여전히 dirty 구 브랜치=GPU 차단).
- Iter12 KB 대조(tag_query.py, 읽기전용): 내 density 결론을 KB SoT에 대조 → ⭐**MCF7 myosin
  density datum 부재 확증** + 새 앵커 **KB-3.18 cortex-wide 분자 myosin ~100/µm²**(~30/minifil →
  minifil ~3.3/µm²). 모델 proxy 0.6의 ~5×지만 밴드 요구 ~20-58/µm²의 여전히 6-17× 아래 ⇒
  KB generic 분자 consensus조차 밴드보다 한 자릿수 낮음 → 밴드 direct-dipole 도달불가 가능
  (coherence 필요 OR 밴드=de-adhered artifact). 브리프 결정2에 KB 인용 보강.
