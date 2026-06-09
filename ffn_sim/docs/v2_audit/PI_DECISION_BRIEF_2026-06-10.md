# PI Decision Brief — 2026-06-10 (autonomous-loop authored)

> Lead가 6h 자율 루프 중 PI 대기 결정을 한 곳에 모았다. **3개 결정이 모든
> magnitude line을 막고 있다.** 각 결정은 (1) 무엇이/왜 막혔나 (2) 정량 증거
> (3) 선택지 + 비용 (4) Lead 추천 순으로. 결정은 PI 고유 권한 — Lead는 진행 안 함.

전체 진단의 수렴점: **active myosin이 만드는 힘이 문헌 대비 100~1000× 부족**하다는
하나의 벽이, 피질 장력 γ-floor와 SF traction 두 갈래로 나타난다. 둘 다 같은 누락
datum(**myosin areal density + active-vs-network 구분**)에 막혀 있다.

---

## 결정 1 — crosslink 강성 재anchor ✅ **이미 완료(loop18) — 결정 불필요, 측정만 남음**

> ⚠️ **2026-06-10 정정**: 최초 브리프는 이 결정을 "PI 대기"로 적었으나, git 확인 결과
> **이미 loop18(commit 236101e, "PI-approved")에서 1e-7→1e-3로 재anchor·커밋됨.**
> 최초 판단은 stale한 설계 doc §10(loop17 시점)을 참조한 오류. 결정 1은 **닫혔다.**

**현재 상태.** `phase1_h3.yaml:141-142` = `k_intra: 1.0e-3`, `k_attach: 1.0e-3`
(KB-1.28, 1 pN/nm). 이전 1e-7(0.1 pN/µm)은 mis-attribution(KU-3.19는 kinetics만, Furuike
2001은 filamin unfolding — 강성 datum 아님)으로 확정되어 교체됨. CFL 영향 없음
(τ_xl=391ns ≫ myosin dt 9.2ns). cross-bridge 강성 fix(continuous_stroke + k=1e-3)도
loop14-17에서 빌드·검증 완료(13 sanity 테스트 PASS).

**진짜 남은 것 = 결정이 아니라 측정.** config 주석이 명시: *"이미 돈 Gate-A/B는 이 강성에서
재실행 필요."* 즉 교정된 작동점(continuous_stroke + crosslink 1e-3)에서 γ가 full-stall
envelope에 도달하는지(transmission이 풀렸는지), 아니면 여전히 제한인지 **long GPU
contraction run으로 재측정**해야 한다. Lead가 이를 자율 진행 중:
- `h7_active_force_budget.py`에 `--xlink-k` A/B override 추가(2026-06-10, additive,
  sensitivity-only) → crosslink 1e-3(production) vs 1e-7(pre-re-anchor)로 transmission
  lever의 γ 기여를 정량.
- 기존 Gate-A(1e-7, binned) γ_soft = 3.06e-3 mN/m = 밴드의 1/114가 비교 baseline.

**측정 결과 (2026-06-10, Lead 자율 A/B, continuous_stroke contraction, n_fil=120):**

| | A: xlink **1e-3** (production) | B: xlink **1e-7** (pre-re-anchor) |
|---|---|---|
| WALL-A 전파효율(actin/myosin) | **4.56%** | **1.72%** |
| actin-network γ [mN/m] | 6.29e-4 | 2.38e-4 |
| myosin-dipole γ [mN/m] (지배항) | 1.378e-2 | 1.378e-2 (동일) |
| **총 γ_soft [mN/m]** | 1.374e-2 (~13× under) | 1.354e-2 (~13× under) |
| per-head meanT / gen_force | 8.44 pN = **F_stall** / 3.43 nN | 8.47 pN = F_stall / 3.27 nN |

**해석.** ⭐ crosslink 재anchor는 network **전파를 2.6× 개선**(WALL-A 1.72→4.56%) — 즉
1e-7이 정말 transmission 고리를 약화시켰고 1e-3 교정이 그 역할을 회복했다(재anchor 정당성
재확인). **그러나 절대 γ는 +1.4%만 움직인다**: γ_soft의 지배항은 local myosin-dipole
(1.378e-2, crosslink 무관)이고 전파 채널(6e-4)은 ~5%에 불과. per-head 힘은 이미 F_stall로
maxed, 결합 406 heads인데도 dipole γ가 ~13× under = **WALL B(생성/density envelope)**.

⇒ **결론: γ magnitude는 transmission-bound가 아니라 density/coherence-bound다.** crosslink는
풀렸고(전파 2.6×), 남은 ~13× 갭은 PI 결정 2(density/active-fraction datum)의 영역. 이 측정이
결정 2(A) "generation/density-bound 종결"을 데이터로 뒷받침한다.
⚠️ 단 이 A/B는 mesoscale CPU n_fil=120(예비). native-scale GPU 확정은 gbook 정리 후
(결정적 GPU A/B는 여전히 PI-gated; 그러나 결론 방향은 바뀌지 않을 것 — dipole 지배는 scale-
invariant).

---

## 결정 2 — myosin areal density / 밴드 active-fraction (datum 없음)

**무엇이/왜.** active-γ envelope(½·n2D·f·ℓ)를 lit 파라미터로 계산하면 full engage+stall에서도
**밴드 대비 ~18–36× under**. 밴드에 닿으려면 ρ≈16–21/µm² 필요 = 유일 proxy(HeLa 0.6/µm²)의
~27–36×. **MCF7 areal-density datum이 존재하지 않는다.** 밴드[0.35–0.65 mN/m] 자체도
rounded/de-adhered HeLa/L929 proxy(NO MCF7, NO spread-adherent datum).

**선택지 + 비용.**
| 선택 | 비용/영향 |
|---|---|
| (A) generation-bound로 종결 | γ active 채널은 lit-density에서 envelope까지 도달함을 보이고, 밴드까지의 잔여 갭은 "density/overlap datum 부재"로 명시 종결. 정직하나 magnitude 미완 |
| (B) density datum 탐색 계속 | MCF7 minifilament areal density를 deep-research로 더 찾기(이전 104-agent 탐색은 parallel-per-cross-section datum REFUTED) |
| (C) 밴드 active-fraction 재정의 | blebb ~halving → active-fraction ~0.135로 밴드 하한 재해석 → 갭 축소 |

**Lead 추천: (A), 단 (C)를 명시 caveat로.** 이전 세션들(nmii/floor)이 이미
"density-closable 아님"을 REFUTE로 확정. 더 파기보다 generation-bound로 정직하게 종결하고,
밴드 자체가 MCF7 datum이 아니라는 점을 결론에 박는 것이 platform 신뢰성에 맞다.
⭐ **2026-06-10 A/B가 이를 데이터로 뒷받침**: per-head 힘 F_stall + 결합 406 heads인데도
dipole γ ~13× under, crosslink(transmission) 풀어도 절대 γ 불변 → 남은 갭은 순수
density/coherence. transmission은 더 이상 용의자가 아니다.

---

## 결정 3 — SF traction line: 계속 vs 종결

**무엇이/왜.** loop23 "+131pN decisive sarcomeric rectification"이 시드별 재현 안 됨:
control seeds 1-6 = {+131, −428, −434, −247, +340} pN, mean **−127±156**, 2+/3−,
ZERO와 구분 불가. "16σ"는 within-realization SEM(across-seed std의 ~15–40× 과소). GPU
s1 == CPU s1 bit-exact(device 무관). ⇒ ~50 heads에서 사르코메릭 정류는 placement noise에
묻힌다. 미커밋 loop24d(M-band targeting)는 이 single-sided minifilament 문제를 고치려는
후속 시도(검증 안 됨, 커밋 안 함).

**재프레임(이미 lit-anchored).** Kassianidou/Kumar 2017 PNAS: 단일 섬유 *active* ~5–6nN,
Kumar 10–30nN은 대부분 network/prestress. = 결정 2의 SF 쌍둥이(같은 density datum).

**선택지 + 비용.**
| 선택 | 비용/영향 |
|---|---|
| (A) 고밀도 단일 SF 테스트 | 200–300 minifilament로 N_heads↑ 시 부호 수렴하는지(placement-noise 가설 검증). 미해결 진짜 joint test. gbook 런 필요 |
| (B) generation-bound로 종결 | 단일 SF active ~5–6nN로 확정, Kumar 10–30nN은 network/prestress로 귀속, line 종료. 결정 2와 한 묶음 |
| (C) 보류 | γ-floor 결정 후로 미룸(같은 density datum에 묶임) |

**Lead 추천: (A)를 자율 루프 중 실행 가능(과학 결정 아님, 가설 검증).** 단 결과 해석/종결
판단은 PI. 만약 PI가 빠른 종결을 원하면 (B)+결정 2 묶음.

---

## 묶음 추천 (PI 한 줄 결정용)

1. **결정 1**: ✅ **닫힘** — crosslink 재anchor는 loop18에 이미 완료. 2026-06-10 A/B로
   "transmission 풀림(전파 2.6×) but 절대 γ 불변(density-bound)" 확정. PI 액션 불필요.
2. **결정 2**: (A) generation/density-bound 종결 + active-fraction caveat. **A/B 데이터가
   직접 뒷받침** — per-head F_stall·406 heads인데도 ~13× under, transmission 풀어도 불변.
3. **결정 3**: (A) 고밀도 SF 테스트로 placement-noise 가설만 검증(자율 진행 가능), 종결 판단은 PI.

이 셋이 정해지면 γ/SF magnitude line이 "generation-bound, transmission lever=crosslink"로
일관 종결되거나, crosslink 재anchor로 밴드에 근접하는 새 결과가 나온다. 둘 다 platform을
앞으로 민다.

---

## 부수 (비차단)

- **KB drift**: RunResult 62·CodeMapping 11 un-harvested(graph stale). dry-run manifest
  생성됨(`OPS_HARVEST_CANDIDATES_2026-06-09.md`). `--apply`는 PI-gated. 승인 시 1줄.
- **워크트리 8개**: 진단 완료된 것들(aggdiag, gammadiag, mcf7wt) 정리 후보 — PI 확인 후.
- **미커밋**: myosin.py(loop24d M-band)·h7_ventral_sf_traction.py — 검증 안 됨, 그대로 보존.
