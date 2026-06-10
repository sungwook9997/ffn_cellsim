# Slater 2021 cross-model transmission oracle — DRAFT (2026-06-11)

**Status:** DRAFT (gate-contract → PI sign-off 필요, CLAUDE.md Roles). 트리거: Slater et
al. 2021 *Soft Matter* 17:10274 (DOI 10.1039/d0sm01911a; TAG key
`10274-soft-matter-2021-17-10274-10285-55a265`).

## 1. 동기 — 왜 oracle인가
ffn_cellsim의 validation/oracles/는 지금까지 **closed-form 논문 모델**(Chan-Odde,
Pereverzev, Bell-Evans, Buckley, Hill)을 acceptance oracle로 쓴다. Slater 2021은 closed
form이 아니라 **같은 계열의 discrete actomyosin+ECM 모델(peer model)**이다. peer model은
정량 band를 주진 못하지만(아래 §4), **창발적 정성/스케일 법칙**을 준다 — ffn이 독립적으로
같은 법칙을 재현하면 transmission 물리 구현이 옳다는 강한 cross-model 일관성 증거가 된다.
(literature-first/no-fitting 준수: Slater에 **맞추지 않고**, ffn이 독립적으로 같은 법칙을
내는지 확인.)

## 2. Oracle 법칙 (정성/스케일 — 수치 band 아님)
수축세포가 가교 fiber 네트워크로 stress를 전파하는 setup에서:

- **O1 (차원 스케일).** 강한 전파 영역에서 stress(r) ∝ **r⁻¹ (2D 원통)** / **r⁻² (3D)**.
  ffn은 3D fine-grained → **r⁻²**가 나와야 함(Slater의 2D r⁻¹과의 차이 자체가 차원 sanity).
- **O2 (contractility).** `κ_s,c`(또는 ffn의 myosin 활성/F_stall) ↑ → peak stress ∝ κ_s,c,
  전파거리 ↑. (Slater Fig.2A–C)
- **O3 (crosslink density).** 밀도/fiber 농도 ↑ → **장거리 전파 + 느린 이완**; ↓ → 단거리 +
  빠른 이완. (Slater "lower concentration → shorter-range, rapid relaxation")
- **O4 (force-dependent unbinding).** Bell-Evans `x_ub > 0`(힘감응 off) → 강한 수축 하에서
  **enhanced stress relaxation** (vs `x_ub=0` 힘무관 대비). (Slater Fig.2 inset)

## 3. 측정 프로토콜 (Slater 차용)
- **stress(r)**: 반경 r의 측면 원통(구)면을 가로지르는 모든 chain(fiber+crosslink)의 spring
  force radial 성분 합 / 면적. (현 coherent-FA-traction differential과 **독립** 2차 추정기 —
  loop24b 불안정 진단에도 유용, [[AB_XLINK_TRANSMISSION_REDESIGN_2026-06-11]] §4 참조.)
- **deformation(r)·relaxation(t)**: fiber 변위 radial 투영의 annulus 평균; stress(t) 시계열.

## 4. 왜 정량 band가 아니라 정성 oracle인가
- 기하/스케일 불일치: Slater = 2D 원통 cell(R=10µm)+ECM(R=100µm); ffn = 3D fine-grained
  단일세포. 절대 stress 값은 직접 비교 불가.
- peer model = 1차 실험 ground truth 아님. 따라서 **gross-error 잡는 cross-model sanity**
  지, magnitude validation 아님. magnitude 진실값은 여전히 1차 실험데이터(MCF7 TFM 등).
- 판정 = 부호/단조성/스케일 지수 일치(O1–O4), 수치 동일성 아님.

## 5. 적용 범위 / 필요 runnable
- 대상: ffn **H.1 ECM + bridge/FA transmission** (수축세포-in-ECM). cortex active-γ
  magnitude **직접 검증 아님**(Slater는 cell-in-ECM 전파 물리).
- 필요: Slater setup에 대응하는 ffn 하니스(수축 cell + 가교 fiber ECM + radial stress
  readout). 현 `h7_ventral_sf_traction`(FA-anchored SF)와 다름 → **신규 하니스 필요**
  (또는 H.1 ECM + 수축 inclusion). 범위/비용 PI 판단.
- 위치: 비준 시 `ffn_sim/validation/oracles/` (discrete-model peer oracle 구획, closed-form과
  구분). gate ID 예: `VG-ECM-transmission-slater2021` (PI-authored).

## 6. 한 줄
Slater 2021을 **cross-model 정성 oracle**(O1 r⁻² 차원 / O2 contractility / O3 density /
O4 force-unbinding)로 채택해 ffn ECM-transmission 구현을 독립 검증. 정량 anchor 아님 —
gross-error 게이트. 신규 하니스 + gate ratify는 PI 결정.
