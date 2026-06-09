# 새 세션 부트 프롬프트 — option 《가》 myosin cross-bridge custom-force 재설계 (2026-06-09 작성)

> 아래 블록을 새 세션 첫 프롬프트로 그대로 붙여넣으세요.

---

너는 ffn_cellsim의 단일 Lead Claude Code 세션이다. CWD=/Users/sw1/ffn_cellsim, 브랜치=h7/full-cell-integration.
다른 세션은 /Users/sw1/ffn_cellsim-platform(h7/compartment-platform)에서 구획 플랫폼을 한다 — 그쪽 파일·브랜치·configs/recipes/·8개 구획 모듈·compartment_registry.py·layer2/는 절대 건드리지 마라.

[부팅] conda activate ffn_sim → python -c "import hoomd; print(hoomd.version.version)" (7.0.1 확인).
다음을 순서대로 읽어라:
1. CLAUDE.md 전체.
2. docs/v2_audit/H7_MYOSIN_OVERLAP_MECHANISM_DESIGN_2026-06-09.md **전체, 특히 §9** (= 이번 미션의 완전한 명세).
3. docs/v2_audit/H7_ADHERENT_VENTRAL_PIVOT_2026-06-09.md (작동점/구조/lane).
4. docs/v2_audit/H7_ARCH_AUTORUN_LOG_2026-06-09.md (loop1-13 전체 궤적).
5. memory project-gamma-floor-layered-resolution (맨 위 "2026-06-09 ARCH session-2 FULL" 항목).
6. Notion Dev Logs 데이로그 37a120daec5d81afadfbc074ac6f6810.
7. git log --oneline -8 (최신=723225e) → git status (clean 확인).
8. 핵심 코드: cortex/myosin.py (grip_walk stepping + bind + cortex_myosin_attach_bin_rest_lengths),
   configs/phase1_h3.yaml (myosin block), bridge/motor.py (canonical head_spring_k=1e-3).

[배경 — 왜 이 미션인가 (이전 세션 확정)]
active-γ 바닥의 per-head force 결손(D1)의 근원은 **cortex cross-bridge stiffness k_head=1e-6이 프로젝트
자체의 canonical motor(bridge/motor.py 1e-3=1 pN/nm; 문헌 Veigel/Kaya 0.3-2 pN/nm)보다 1000× 약한
pN/µm-vs-pN/nm 단위 슬립**이다. PI가 교정을 승인했다. 그러나 단순 config 변경(k→1e-3)만 하면 grip_walk의
r0≈0 컨벤션 때문에 갓 결합한 head가 k·r≈322pN로 폭발(비물리적). 그리고 bin-양자화 attach 본드는 stiff
cross-bridge의 ~4nm stroke를 330nm 결합범위에서 분해 못 한다(n_bins≳165 필요=비현실적). ⇒ **단위슬립
교정 = myosin contractility 재설계**이다. (이전 세션이 bin 접근 시도→실패→되돌림, §9에 완전 명세.)

[미션 — option 《가》: 연속 per-head custom-force myosin contractility (§9의 atomic set 구현)]
다음을 **작은 단위로 sanity-gate 검증하며** 구현하라 (core contractility라 신중히):
1. **연속 per-head 힘** — harmonic-bin attach-FORCE를 md.force.Custom으로 교체: 결합 head마다
   F = min(k·s_grip, F_stall)를 head→bead 단위벡터로 가하고 반력은 head/backbone에. attach 본드는
   Bell-Evans off-rate 부기 용도로만 유지(이중계산 금지). s_grip은 기존 grip_walk 누적 + Hill stall 재사용.
   ⚠️md.force.Custom은 per-step host-sync(GPU-main 영역) → cupy/native 경로가 production form. opt-in
   stepping_mode로 만들어 legacy grip_walk byte-identical 보존.
2. **k_head_spring/k_head_actin 1e-6→1e-3** (phase1_h3.yaml; Magic-Number Block + canonical/Veigel/Kaya anchor).
3. **head_actin_max_bind_dist 축소** (현 330nm은 stiff bridge엔 비물리적 — stiff cross-bridge는 저변형
   ~수십 nm에서 결합; 유도근거 명시, 튜닝 금지).
4. **CFL 재유도** — k_backbone=1e-2 → τ_backbone~39ns(native)·~9ns(mesoscale). dt 안정성 실측 →
   필요시 dynamics.cfl_safety_factor(cortex config)로 축소, OR backbone을 진짜 rigid constraint로
   (브리프 의도; 더 나은 최적화). **integrator/는 PI-freeze — 절대 수정 금지(읽기만).**
5. **crosslink k_intra/k_attach 재anchor** — 1e-7(0.1 pN/µm)은 H.1 가교 1e-3(KU-1.28)·문헌보다 ~10000×
   약함 = fiber↔fiber 전달 고리. α-actinin 강성 문헌을 KB(tag_query.py)로 검증 후 교정(또는 PI surface).
[검증 게이트 — §5 sanity gates + 다음]
- per-head **DELIVERED** force가 F_stall에 cap(k·r 폭발 아님) — force-budget의 실제 g_myo로 확인(322pN 아님).
- 차원/경계(zero overlap→0; full→envelope)/보존(Newton-3)/부호(contractile)/측정일관성(g_actin가 myosin로
  RISE — crosslink fix 후 transmission)/CFL/no-gate-chasing(밴드 LOCKED; γ→full-stall envelope 18-36×
  under = density/overlap 갭은 남는 알려진 갭, 튜닝 금지).
- 끝나면 ventral SF 2b-2 재개: 교정된 myosin을 ventral_stress_fiber 위에 결선 → 2c 차분·시간평균 traction.

[소유권 — HARD] 쓰기 허용: cortex/(myosin.py, ventral_stress_fiber.py, manifold_index.py 등), cell/cell.py,
cell/manifest.py, configs/phase1_h3.yaml(myosin/crosslink 블록), scripts/h7_, outputs/h7, 관련 tests.
읽기만: bridge/(fa.py, motor.py), ecm/substrate.py(소비 인터페이스), 그 외. 절대 금지: integrator/(PI-freeze),
8개 구획 모듈+compartment_registry.py+configs/recipes/+layer2/+h7/compartment-platform 브랜치+다른 세션 파일.

[HARD RULES] no magic numbers(파라미터는 유도가능·grid-invariant·튜닝아님 증명 or PI surface). Hosseini/Gate
밴드 LOCKED — gate-loosening 금지. physiological-baseline(부착 작동점: FA-ON, 실제 점도/turgor). integrator-
freeze 편집·gate-contract 변경·magic-number 트리거·lane 밖 작업 = 중단·문서화·PI surface. k=1e-3과
head_actin_max_bind_dist 변경은 PI가 이번 미션으로 승인함(단 문헌 anchor + Magic-Number Block 필수); crosslink
k 교정은 문헌 검증 후, 불확실하면 PI surface.

[인프라] gbook(RTX A5000)=GPU 프로덕션; Mac→gbook rsync 배포는 PI 승인 remote-overwrite. Syncthing=outputs만.
재사용 도구: scripts/h7_active_force_budget.py(g_soft/g_myo/g_actin 채널 분리, mean_T=k·r은 tool 아티팩트라
오해 말 것 — 실제 force는 g_soft/g_myo로), h7_ventral_sf_traction.py(SF 스캐폴드), h7_manifold_index_validate.py.
긴 작업엔 ETA 적어라.

[커밋 규율] 파일 명시, git add -A 금지. h7/full-cell-integration에만. ffn/foundation 푸시는 PI 승인.
메시지 끝: Co-Authored-By: Claude Opus 4.8 (1M context) noreply@anthropic.com. 커밋 메시지/PR 본문은 heredoc로
persistent shell에 직접 넣지 말고 -F file 사용(stdin 손상 주의).

[루프] 가설/측정 1스텝 → 실행 → outputs/h7 + REPORT/로그 → 커밋 → docs/v2_audit/AUTORUN_LOG_<date>.md 한 줄
append → 반복. 검증 안 된 core-physics는 커밋 금지(되돌리고 surface). 클로즈아웃 시 Notion Dev Logs 데이로그 +
figs 갱신, 마지막 줄 "Notion 업데이트 완료".

[정지조건] integrator-freeze 변경 / gate-contract 변경 / 새 magic-number 트리거 / lane 밖 작업 필요 / 검증
실패한 core-mechanism → 중단·문서화·PI surface.

/goal §9의 atomic set(연속 per-head custom-force + k=1e-3 + 결합범위 축소 + CFL + crosslink)을 sanity-gate
검증하며 구현해, per-head delivered force가 F_stall에 cap되고 γ가 full-stall envelope에 도달함을 보여라.
그다음 교정된 myosin으로 ventral SF traction(2b-2→2c)을 측정하라. 막히거나 PI 결정 필요하면 surface.
