# ab_xlink transmission 테스트 재설계 — PI 브리프 (2026-06-11)

**Status:** PROPOSED (PI 결정 대기). 트리거: Slater et al. 2021 *Soft Matter* 17:10274
(DOI 10.1039/d0sm01911a; TAG key `10274-soft-matter-2021-17-10274-10285-55a265`) 검토.

## 1. 현재 ab_xlink A/B가 무엇을 테스트했나
`outputs/h7/production/ab_xlink/` (fb_cont_xl1e3 vs fb_cont_xl1e7): crosslink **강성**
`k_xl`을 1e3 vs 1e7 (N/m 스케일)로 흔들어 contraction γ 변화 측정 → **γ ~5%만 변함**
→ 결론 "transmission은 active-γ 벽이 아니다 (crosslink lever 아님)". 이 결론은
GPU_MAIN_PORT_PHASE2 §4와 loop24b HALT 서사("transmission 아님, density/coherence
문제")의 한 축.

## 2. Slater 2021이 말하는 것 — 우리가 틀린 손잡이를 돌렸을 수 있다
같은 discrete actomyosin+ECM 모델 패밀리(Taeyoon Kim 그룹)에서, 수축세포가 가교된
네트워크로 stress를 **전파**하는 거리·크기·이완을 좌우하는 변수는:
- **crosslink 밀도(fiber concentration + cross-linking density)** — 높을수록 연결성↑ →
  **장거리 전파 + 느린 이완**; 낮을수록 단거리 + 빠른 이완 (논문 Fig.2, "Lower fiber
  concentrations → shorter-range stress transmission, rapid relaxation").
- **힘의존 unbinding kinetics** — Bell's law `k_ub = k_ub,0·exp(x_ub|F|/kT)`
  (ref values `k_ub,0=1e-5 s⁻¹`, `x_ub=1e-10 m`). zero-force rate + force sensitivity
  `x_ub` 두 파라미터가 stress의 buildup/transmission/relaxation을 조절.
- contractility `κ_s,c` ↑ → 변형↑ + 전파거리↑, peak stress ∝ κ_s,c. 강한 수축에선
  force-dependent unbinding이 enhanced relaxation 유발.
- **강성은 1차 전파 손잡이가 아님.** Slater의 전파/이완 sweep은 density·k_off·x_ub
  위에서 도는 것이지 crosslink 강성 위가 아니다.

→ **우리 ab_xlink는 전파를 강성으로 프록시했는데, 전파의 실제 손잡이는 밀도+off-rate
(Bell-Evans)다.** "γ 5%, transmission 아님"은 *강성에 대해* 참일 뿐, transmission이
active-γ에 무관하다는 일반 결론의 근거로는 **불충분**할 수 있다.

## 3. 제안하는 재설계 (PI 승인 시)
crosslink **강성 대신 (density, k_off0, x_ub)**를 sweep하는 transmission 진단:
1. **축**: (a) crosslink 밀도(개수/단위면적 또는 nucleation rate), (b) Bell-Evans
   zero-force off-rate `k_off0`, (c) force sensitivity `x_ub`. 3D fine-grained cortex/SF
   에서 — 우리 α-actinin은 `k_off0=0.066 s⁻¹` (Ferrer 2008) 앵커이므로 그 ±배수 범위.
2. **관측량(2개, Slater 프로토콜 차용 — §4/(b) 참조)**: (i) **공간분해 stress(r)** =
   반경 r 면을 가로지르는 chain들의 axial spring force / 면적 (현 coherent-FA-traction
   differential과 독립적인 2차 추정기), (ii) **시간분해 relaxation** stress(t).
3. **판정**: density↑/k_off↓ 에서 stress 전파거리↑·이완 느려짐이 나오는가? active-γ가
   density 손잡이를 따라 움직이는가(움직이면 transmission이 magnitude에 기여 — loop24b
   서사 수정; 안 움직이면 "generation-limited" 확정 강화).
4. **출력**: `outputs/h7/production/ab_xlink_v2/` + figure + REPORT, AB_XLINK 기존
   결론을 "강성-불변" + "density/kinetics-감응(또는 비감응)"으로 정밀화.

## 4. 주의 (하드룰)
- Slater는 **모델링 논문**(1차 실험데이터 아님) + crosslink이 **collagen ECM 가교**
  (transglutaminase, `k_ub,0=1e-5 s⁻¹`)다. 우리 **cortex는 α-actinin**(`k_off0=0.066
  s⁻¹`). **원리(density+kinetics가 전파를 지배)는 이식되나 값은 맥락별** — cortex
  active-γ엔 α-actinin 값, ECM 전파엔 collagen 값. magic-number 금지: 모든 sweep 값은
  KU-앵커(α-actinin Ferrer 2008 등) ±배수로, 문헌 근거 명시.
- 이 재설계는 magnitude **특성화** → loop24b HALT와 같은 PI 결정 영역. 본 문서는 제안서.

## 5. 한 줄 결론
"transmission은 벽이 아니다"는 *강성*에 대해서만 검증됐다. Slater는 전파의 손잡이가
**밀도+Bell-Evans off-rate**임을 보이므로, active-γ가 transmission-limited인지 아닌지
판정하려면 그 축으로 다시 sweep해야 한다. → PI 승인 시 ab_xlink_v2 착수.
