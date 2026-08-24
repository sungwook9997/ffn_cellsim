# ffn_cellsim — 인계 프롬프트 (Opus 4.x → Opus 5)

당신은 `ffn_cellsim`(NVIDIA Warp/CUDA 단일세포 mechanobiology 시뮬레이터)의 Lead Claude Code
세션입니다. 전임 세션(Opus 4.x)이 장시간 자율 작업을 마치고 당신에게 인계합니다. PI는 Sungwook이며
한국어로 대화합니다(코드/커밋/문서는 영어 관례). 아래를 읽고 이어가세요.

## 0. 부팅 (minute 0)
1. `CLAUDE.md` 전체 정독 — 아키텍처 원칙(fine-grained/mechanistic > lumped), I0-A 계약(Warp-CUDA가 유일
   런타임, HOOMD 금지), 하드룰(magic-number 금지·gate-loosening 금지·생리 baseline·NATIVE+FULL 검증).
2. `git log --oneline -30` + `git branch` — 현재 브랜치 `codex/ff-ac-codex`, HEAD 근처 `495ac651`
   (첫 composed native cell). **`ffn/foundation`에 절대 push 금지**(PI 승인 사항).
3. auto-memory `~/.claude/projects/-Users-sw1-ffn-cellsim/memory/MEMORY.md` 인덱스 훑기.
   특히 `project-ac-fine-mesh-multigrid.md`(fine-mesh 근본원인+fix)와 [[feedback-*]] 피드백들.
4. SoT 문서: `docs/v2_audit/cell_engine/ROLLING_ROADMAP.md`(현재 현실 반영됨, 2026-07-25 갱신) +
   `NATIVE_COMPOSITION_SCOPE_2026-07-25.md` + `INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md` +
   `FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md` + `CORTEX_KINK_AUDIT_2026-07-25.md` +
   `PI_GAP_EVIDENCE_CARDS_2026-07-25.md`. Notion "Full-cell breadth" milestone(page 3a7120da…4deb).
5. `conda activate ffn_sim`. gbook(RTX A5000)이 유일 CUDA GPU — **UP 상태**(Tailscale). Mac은 CPU 전용
   (커널 SOURCE + NumPy reference; native launch은 gbook). ssh 배포: `rsync -az <file> gbook:~/ffn_ac_native/…`,
   실행: `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python …`.

## 1. 지금까지 (이번 인계 세션의 성과, 전부 codex/ff-ac-codex, additive, adversarial-verified)
- **cortex fine-mesh 완결**: fine 75nm(2.9M) 비수렴 plateau(max|PF| 3.42)의 근본원인을 **측정+감사+production
  으로** 규명 = `overlap_free` WCA relaxation이 교차 node를 매끈 arc 밖으로 밀어 fine서만 5° kink→bending 발산.
  **solver도 steric도 아님**(둘 다 production 반증). smooth arc면 3.42→0.13 수렴. `radial_span` fix로 12×
  (3.42→0.28, steric ON). fiber-arclength multigrid 설계+구현+속도fix도 완료(부산물, plateau는 solver 문제가
  아니었음). ⚠️ 성급한 가설(grid-dependent steric)을 세웠다가 `--no-steric` production으로 스스로 반증·정정함 —
  **이 규율(공간상관≠인과, production으로 검증)이 PI가 가장 중시하는 것.**
- **native ECM 3층 CUDA-verified**: `codex/ecm-gpu-topology` clean additive merge(device SoA topology +
  remodel transaction, 22 CUDA tests) → 구성법칙 바인딩(collagen link_spring+cytosim_bending over SoA,
  relaxed force-free/shear-restoring) → owner-driven(ecm_world). R3 un-gated = SF→FA→ECM traction spine 개통.
- **sf_arc**: SEAMED→KERNEL_BOUND(자체 커널 launch, disjoint population id_base 1e6), native-closed
  (CFL-stable relax, Gershgorin tangent bound로 dt 도출). milestone viz browser-verified.
- **🏆 첫 composed native cell**(`495ac651` + gbook gate PASS): `build_native_composed_cell_world`가 REAL
  owner(sf_arc+ECM+NMII, cortex=bind-target PORT)를 **하나의 CellTransaction**에 바인딩 → full native
  population서 CUDA step: rest≈0, perturb→창발, NMII 3465/8840 heads 이벤트 결합, accepted-step clock 0→1,
  4 population disjoint. **N 슬라이스 → 하나의 composed 드라이버 = full-cell 통합 밀스톤.**
- interior-column CellTransaction scaffold(정직: coupling은 아직 driver-owned, native-runs 40 steps).
- PI-GAP evidence cards(33 슬롯/9 컴포넌트 fill-in 폼). Codex(GPT) 자문받고 판별. SoT 정정. milestone viz
  (sf_arc+ECM browser-verified).

## 2. 🔴 PI-DECISION QUEUE (전부 상신+문서화됨, PI 사인오프 대기 — PI 깨어나면 우선)
1. **coupling-ownership-transition**: interior fluid column을 진짜 CONNECTED로(connector가 coupling 소유).
   `_accumulate_all`에 flag-gated omit-mask + injected pass. **feature-frozen `aleph/components/incumbent/driver.py` 계약 변경 =
   Card-5 strangler 사인오프 필요.** 정밀 plan = `INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md`. 작은 첫
   단계: pressure만 먼저 이동. (자율로 driver 편집 금지.)
2. **α2β1-collagen catch-bond Hand card**(FA clutch KMC blocker, 소싱 필요; PI-GAP Card A1).
3. **PI-GAP cards 33 슬롯**(NMII 클러스터가 cortex+SF 동시 unblock; MT DI/IF WLC/filopodium fascin/ERM Bell/
   turgor). literature-first, PI가 소싱하거나 선택.
4. turgor Π₀ Card-3(72Pa), radial_span default+window refinement, collagen-modulus + PI-GAP Notion SoT 등록
   (PI-authored gates).

## 3. 다음 자율 작업 (non-PI-gated, code-ready)
- **composed cell 확장**: 더 많은 private-array owner/connector를 `build_native_composed_cell_world`에
  (membrane surface, lamellipodium branch-angle 등 — 단, MT/IF/filopodium은 PI-GAP-blocked). 패턴은
  `aleph/engine/composed_native.py`(495ac651). interior fluid column은 Card-5 사인오프 후.
- 각 확장은 **adversarial-verify workflow**(real-owner-not-facade + one-CellTransaction + disjoint + additive)
  로 검증하고, gbook native-gate로 확인.

## 4. 작업 방식 — PI가 중시하는 것 (반드시 지킬 것)
- **adversarial verify로 real vs theater 판별.** 매 substantive 작업(특히 workflow)에서 검증 에이전트가
  "facade/스캐폴드/driver-deferred인가"를 refute 시도. 통과한 것만 커밋. 스캐폴드는 **정직하게 라벨링**(예:
  interior-column은 "composition scaffold, coupling driver-owned"로 커밋). structural pass ≠ production.
- **production/native로 검증, 틀리면 즉시 정정.** 결론을 추론이 아니라 측정으로. 성급한 판단 금지("성급하게
  판단 내리지마"). 공간상관≠인과.
- **하드룰**: magic-number 금지(derivable/sourced/grid-invariant; GAP은 required-no-default로 wire, fallback
  거부). gate-loosening 금지(gate 틀리면 PI에 계약변경 상신). PI 실험데이터에 fitting 금지(literature-first).
  **생리 baseline**(turgor-pressurised 등 real in-vivo서 시작, default-off도 production선 ON). **NATIVE+FULL**
  population(70,686 cortex)서 검증(coarse/stripped로 결론 금지). Warp-CUDA만(HOOMD 금지). co-location≠connection.
- **PI-gated는 자율로 안 바꾼다**: feature-frozen `ac/cell` 계약, mechanistic 모델 결정, 새 gate/contract(PI-
  authored), 파라미터 소싱. 정밀 설계+상신하고 내 권고 기록. 단, **막힌 결정에 멈추지 말고 내 권고대로 진행할
  수 있는 부분은 진행**(파괴적/비가역/PI-gated 제외).
- **memory는 짧게**(상세는 doc/Notion/TAG 시스템). **각 milestone viz**(WebGL, `browser_check.py`로 실제
  렌더 검증 — grep은 검증 아님; 육안이 버그 잡음). **Codex는 아이디어용, 검증 권위 아님**(자문받되 내가 판별).
- **PI에겐 한국어로 답변**. 커밋 끝에 `Co-Authored-By: Claude <...>`(원래 Fable 5였으나 당신 버전에 맞게).
  **`ffn/foundation` push 금지·PI 승인 사항.** 롱-런은 SoT 신뢰 + 서브에이전트 병렬화로 drift 방지.
- Ultracode/Max일 때: substantive 작업은 Workflow(fan-out + adversarial verify)로. 토큰 아끼지 말고 exhaustive.
  scout(읽기/스코핑) → orchestrate(workflow) 순서.

## 5. 한 줄 요약
cortex mesh 지적에서 셀 전체 통합까지 실제로 도달했다. 다음은 composed cell 확장 + PI 사인오프 대기 항목.
매 결론을 production/adversarial-verify로 확증하고, 틀리면 즉시 정정하고, PI-gated는 상신하며, 셀을
evidence-ladder 위로 정직하게 밀어올릴 것. — 전임 세션 드림.
