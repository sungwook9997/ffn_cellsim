# PI Decision Brief — 2026-06-10 (autonomous-loop authored)

> Lead가 6h 자율 루프 중 PI 대기 결정을 한 곳에 모았다. **3개 결정이 모든
> magnitude line을 막고 있다.** 각 결정은 (1) 무엇이/왜 막혔나 (2) 정량 증거
> (3) 선택지 + 비용 (4) Lead 추천 순으로. 결정은 PI 고유 권한 — Lead는 진행 안 함.

전체 진단의 수렴점: **active myosin이 만드는 힘이 문헌 대비 100~1000× 부족**하다는
하나의 벽이, 피질 장력 γ-floor와 SF traction 두 갈래로 나타난다. 둘 다 같은 누락
datum(**myosin areal density + active-vs-network 구분**)에 막혀 있다.

---

## 결정 1 — crosslink 강성 재anchor `k_intra/k_attach: 1e-7 → 1e-3 N/m`

**무엇이/왜.** cortex `dynamic_crosslinkers.k_intra = k_attach = 1.0e-7 N/m`(0.1 pN/µm)은
KB-1.28(검증됨, High: H=(k_xl/2)|r|², k_xl 1e-4–1e-2 N/m, default **1e-3 N/m = 1 pN/nm**)
보다 **~10⁴× 약하다**(range floor 1e-4보다도 10³×). config 주석의 anchor
("0.1 pN/µm, KU-3.19 Furuike 2001")는 강성 datum이 아니다 — KB-3.19는 crosslinker
**kinetics(off-rate)만** 명시(k 없음), Furuike 2001은 filamin **unfolding kinetics**.
= 2026-06-02/06-08 audit가 잡은 mis-attribution 클래스와 동일 패턴.

**왜 중요(γ lever).** Gate-A FINAL = REFUTE: γ_soft가 3.06e-3 mN/m = 밴드의 1/114에서
plateau. 진단 결론은 **transmission-limited** — heads는 bind/load/walk/step/contract를
다 하지만(γ 16.2× 상승 확인) rigid M-SHAKE backbone에서 국소 수축이 shunt되어 ~1%만
hoop tension에 도달. crosslink가 fiber↔fiber 전달 고리이고, 1e-7은 그 고리를 ~10⁴× 약하게
만든다. **이것이 가장 유력한 transmission lever.**

**증거 요약.** WALL-A 전파 smoke가 connected mesh에서 이미 66.8%; CFL은 안 조여짐
(τ_xl=391ns @ 1e-3 ≫ myosin dt 9.2ns). 정량 γ 영향은 Lead가 진단 중(task #2, 샌드박스).

**선택지 + 비용.**
| 선택 | 비용/영향 |
|---|---|
| (A) 승인 — 1e-3 재anchor | lit-first(KB-1.28). **production-wide gate-contract 변경**: 모든 cortex build + 이미 돈 Gate-A/B 무효화·재실행. CFL 영향 없음. 단, flexible crosslinker의 분자 강성 ≠ sim harmonic-bond 강성 개념(KB range 100× span) |
| (B) 보류 — 더 검증 | 1e-7의 올바른 단일분자 anchor를 KB에서 더 탐색 후 결정 |
| (C) 먼저 영향 분석 | 재anchor 시 γ가 실제 얼마나 오르는지 샌드박스 정량(task #2) → 그 수치 보고 PI 판단 |

**Lead 추천: (C) → (A).** 진단(task #2)이 "재anchor가 γ를 밴드 쪽으로 유의하게
끌어올린다"를 보이면 (A)가 lit-first로 명확. 끌어올리지 못하면 transmission이 lever가
아니라는 추가 증거 → 결정 2(density)로 무게 이동.

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

1. **결정 1**: (C) 진단 먼저 → 결과 보고 (A)/(B). *Lead가 자율 진행 중.*
2. **결정 2**: (A) generation-bound 종결 + active-fraction caveat.
3. **결정 3**: (A) 고밀도 SF 테스트로 placement-noise 가설만 검증(자율 진행), 종결 판단은 PI.

이 셋이 정해지면 γ/SF magnitude line이 "generation-bound, transmission lever=crosslink"로
일관 종결되거나, crosslink 재anchor로 밴드에 근접하는 새 결과가 나온다. 둘 다 platform을
앞으로 민다.

---

## 부수 (비차단)

- **KB drift**: RunResult 62·CodeMapping 11 un-harvested(graph stale). dry-run manifest
  생성됨(`OPS_HARVEST_CANDIDATES_2026-06-09.md`). `--apply`는 PI-gated. 승인 시 1줄.
- **워크트리 8개**: 진단 완료된 것들(aggdiag, gammadiag, mcf7wt) 정리 후보 — PI 확인 후.
- **미커밋**: myosin.py(loop24d M-band)·h7_ventral_sf_traction.py — 검증 안 됨, 그대로 보존.
