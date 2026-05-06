# ActiveCellSim 파일 구조 지도

작성: 2026-05-06 KST. **파일 → 버전/역할** 매핑 인덱스 문서. 새 세션에서 어디 무엇이 있는지 빠르게 파악할 때 사용.

이 worktree에서 옵션 B가 적용됨: `docs/`와 `tests/`가 v1/v2 sub-folder로 분리됨 (`docs/v1/`, `docs/v2/`, `tests/v1/`, `tests/v2/`). `acs/` 코드는 그대로 (`acs/v2/` 외에는 v1/공유 혼재). 코드까지 분리하는 옵션 C는 V2-2 Unit 1 lock 봉인 후 별도 schedule로 검토.

---

## 0. 한눈에 (최상위 분류)

| Bucket | 위치 | 상태 |
|---|---|---|
| **Shared infra** (version-agnostic) | `acs/{gpu,config,runner,provenance,logging_setup}.py`, `acs/io/`, `acs/analysis/`, 루트 메타파일들 | active, 양 버전이 공유 |
| **v1** (frozen, spheroid-first MPM continuum) | `acs/{physics,boundary,adhesion,visualization}/`, `docs/` 중 v1-era, `tests/` 중 v1, `configs/` 중 stage·path_c·lam4 계열, `scripts/` 중 v1 runner | frozen reference. 새 biology 추가 금지 (CLAUDE.md) |
| **v2** (active, image-constrained cell-resolved) | `acs/v2/`, `docs/v2/v2_*.md` + `00_project_vision_v2.md` + `10_dev_roadmap_v2.md`, `tests/v2/test_v2_*.py`, `configs/final_pilot_v2_*` + `imaging/`, `scripts/run_phase_*` | **active**. 모든 신규 작업은 여기 |
| **v3** | 없음 | MCP 워크룸 이름만 존재 (`v3-design-discussion`, `v3-implementation-work`) — 코드/docs 0개 |
| **데이터 (shared)** | `data/experimental/` (read-only PI 데이터), `data/literature/` (예약, 비어있음) | v1/v2 모두 read-only로 사용 |

진행 상황 % 분포는 `docs/v2/10_dev_roadmap_v2.md`와 lock 파일 상태 참조.

---

## 1. Shared infra (version-agnostic)

| 파일/폴더 | 역할 |
|---|---|
| `acs/__init__.py` | 패키지 루트 |
| `acs/gpu.py` | Taichi 백엔드 dispatcher (cuda/vulkan/opengl/metal/cpu/auto). `ACS_GPU_BACKEND` env var 라우팅 |
| `acs/gpu_profiler.py` | GPU 프로파일링 헬퍼 |
| `acs/config.py` | YAML config 로더 + 검증 |
| `acs/runner.py` | 실행 하네스 (run loop, frame scheduling) |
| `acs/provenance.py` | git hash + config dump (재현성) |
| `acs/logging_setup.py` | logging (DEBUG/INFO/…) |
| `acs/io/hdf5_writer.py` | HDF5 frame writer (양 버전 공통) |
| `acs/analysis/shape_metrics.py` | 형상 메트릭. **`top_down_projection_area`** 가 PI A/A₀ 매칭의 정답 측정 (Hard Rule 11) |
| `acs/analysis/radial_reduction.py` | 방사 ODE 환원 (v1 specific지만 위치는 shared) |
| `acs/viz/__init__.py` | placeholder (실제 v2 viz는 `acs/v2/viz/`) — **dead, 향후 정리 후보** |
| 루트 `pyproject.toml`, `requirements*.txt`, `.gitignore`, `.stignore` | 패키지/의존성 |
| 루트 `CLAUDE.md` | Claude용 프로젝트 컨텍스트 + hard rules + sanity gate protocol |
| 루트 `AGENTS.md` | Agent 역할 + 협업 구조 |
| 루트 `README.md` | 일반 overview |
| 루트 `FIRST_MESSAGE.md` | 세션 부트스트랩 안내 |
| `tools/collab_mcp/` | acs-collab MCP 서버 (PI↔Claude↔Codex 비동기 통신) |
| `logs/`, `results/` | runtime 생성 |

---

## 2. v1 (frozen, spheroid-first MPM continuum)

> 새 biology 추가 금지. 참조용 + 회귀 검증용으로만 보존. 재현 필요시 `docs/v1/v1_continuum_backup.md` 참조.

### 2.1 코드
| 파일 | 역할 |
|---|---|
| `acs/physics/mlsmpm.py` | **v1 본체**. MLS-MPM Taichi 구현 (Stage 1a~2 기반) |
| `acs/physics/__init__.py` | placeholder |
| `acs/boundary/__init__.py` | placeholder (boundary biology 미실구현) |
| `acs/adhesion/__init__.py` | placeholder (adhesion ODE 미실구현) |
| `acs/visualization/dashboard.py`, `plots.py`, `live_imaging.py`, `state_overlays.py`, `parameter_tables.py`, `figure_style.py` | v1-era 차트/대시보드 |

### 2.2 문서

**`docs/` 루트 — 아키텍처 핵심 (v1 origin, 일부는 v2도 참고)**
- `00_project_vision.md` — 프로젝트 WHY
- `01_model_architecture.md` ~ `09_visualization.md`, `11_performance_protocol.md`, `12_validation.md`, `13_data_schema.md` — 모델/검증/스키마/성능 문서 13종
- `references.bib` — 인용 BibTeX

**`docs/v1/` (frozen 기록 39개)**
- `10_dev_roadmap.md` — v1 roadmap
- `v1_continuum_backup.md` — v1 frozen 정당화 + 사용 가이드
- `v15_deferred_diagnostics_plan.md`
- `outcomes_*.md` 11종 (stage1a_plus, stage1a_plus_plus, stage1b, stage1c, stage1d, stage1dc_lam4_cellcount_sweep, stage1e, stage2, v15, path_c, layer3_audit_and_full_stack)
- `stage1*_sanity.md`, `stage1*_review.md`, `stage2_sanity.md` 14종
- `production_lam4_{finding,outcomes,sanity}.md`
- `path_c_sanity.md`
- `parameter_registry.md`, `gate_fail_taxonomy.md`
- `marangoni_review.md`, `layer3_phi_audit.md`, `option_f_gate_contract_sanity.md`
- `anchor_force_balance_investigation.md`, `horizontal_momentum_drift_investigation.md`

### 2.3 테스트 (`tests/v1/`)
- `test_physics_conservation.py` (MLS-MPM 보존)
- `test_measurement_boundary.py` (측정 경계 — v1 자료)

### 2.4 설정 (`configs/`)
- `dev.yaml`, `pilot.yaml` (v1 stage 1a target)
- `production_16gb.yaml`, `production_24gb.yaml`, `production_5k_80hr.yaml`
- `production_lam4.yaml`, `full_stack_production_5k_{bare,lam4,pre}.yaml`
- `path_c_v15_baseline.yaml`
- `stage1a_mid.yaml`, `stage1a_pilot.yaml`

### 2.5 스크립트 (`scripts/`)
- `diag_csf_penetration.py`, `diag_radial_print.py` (진단)
- `stage1e_compare.py` (v1 stage 1e 비교)
- 기타 v1-era render/movie scripts (production_comparison, cellcount_comparison)

---

## 3. v2 (active, image-constrained cell-resolved)

### 3.1 코드 (`acs/v2/`) — 7 그룹

#### (a) State schemas (root)
| 파일 | 역할 | 관련 phase |
|---|---|---|
| `data_contract.py` | imaging dataset spec, metric spec, single-cell contract | V2-1 |
| `single_cell.py` | single-cell state + protrusion event + FA observation + projected polygon area | V2-2 |
| `cell_cluster.py` | 다세포 클러스터 상태 | V2-2/3 (현재 V2-2 single-cell 모드만 사용) |
| `measurement_boundary.py` | 측정/샘플링 경계 조건 | V2-1/2 |
| `adhesion_network_state.py` | 접착 네트워크 상태 (φ 등) | V2-2 adhesion |
| `focal_adhesion.py` | FA 객체 (state level) | V2-2 |
| `protrusion.py` | protrusion 상태 | V2-2 |
| `junction.py` | cell-cell junction 상태 | V2-2/3 |
| `ecm_substrate.py` | ECM 기질 (변형 가능) | V2-2 (단일세포-기질) |

#### (b) Imaging contract (`imaging_contract/`) — V2-1
| 파일 | 역할 |
|---|---|
| `loader.py` | 실험 imaging loader, mask 읽기, pixel-size correction |
| `split_builder.py` | calibration/validation split, frame extraction |

> **현재 V2-2 Unit 1 작업이 이 모듈을 확장 중** (PI 분석 알고리즘 재현 후 `Area_px` 검증)

#### (c) Dynamics (`dynamics/`) — V2-2 핵심
| 파일 | 역할 | lock 상태 |
|---|---|---|
| `active_contour.py` | 단일세포 boundary representation 동역학 | locked + sanity gate |
| `adhesion_network_dynamics.py` | φ ODE (E-cad ↔ Int-β1) | locked |
| `focal_adhesion.py` | FA 동역학 (state machine + 힘) | locked + sanity gate |
| `fa_rate_response.py` | FA rate sensitivity | locked |
| `protrusion_coupled_focal_adhesion.py` | protrusion ↔ FA 결합 | locked + sanity gate |
| `ecm_open_loop.py` | ECM open-loop 동역학 (FA→ECM 단방향) | locked |
| `ecm_constitutive_response.py` | ECM 본구성 응답 (Neo-Hookean + viscosity) — **HB#1.2** | locked + sanity gate |
| `ecm_to_fa_bias.py` | ECM→FA 바이어스 — **HB#4** | locked + sanity gate |
| `fa_to_ecm_scattering.py` | FA→ECM 산란 — **HB#3** | locked + sanity gate |
| `ecm_lyapunov_metric.py` | 안정성 메트릭 — **HB#5** | locked + sanity gate |
| `closed_loop_phase_d.py` | 폐루프 Phase D scaffolding | locked + sanity gate |
| `closed_loop_phase_e.py`, `closed_loop_phase_e_sweep.py` | 폐루프 Phase E composition + sweep | locked |
| `phase_f_minimal_motility.py` | Phase F 최소 운동성 통합 | locked |

> "HB" = Hard Blocker. closed-loop ECM 게이트의 사전 조건들. 모두 locked 상태.

#### (d) Harnesses & runners (root, integration)
| 파일 | 역할 |
|---|---|
| `active_contour.py` | AC harness (gate 단위 테스트용) |
| `active_contour_harness.py` | AC full harness |
| `ecm_open_loop_harness.py` | ECM open-loop harness |
| `phase_e_v2_pilot_runner.py` | Phase E pilot end-to-end runner |
| `phase_f_minimal_motility_pilot_runner.py` | Phase F pilot runner |

#### (e) Output (`output/`)
- `frame_dump.py` — HDF5 frame 출력

#### (f) Visualization (`viz/`) — V2 자체 viz
- `phase_e_v2_pilot_dashboard.py`, `phase_e_v2_pilot_frame_replay.py`, `phase_e_v2_pilot_movie.py`
- `phase_f_annotated_visualization.py`
- `stub3d.py` (테스트/디버그용 3D stub)

#### (g) Metrics
- `metrics.py` — shape / ECM strain / FA force 측정

### 3.2 문서 (`docs/v2/` — 56개)

> 파일들은 모두 `docs/v2/` 안. 아래 목록은 sub-cluster별 그룹.

#### V2 핵심 framing
- `00_project_vision_v2.md` — v2 pivot rationale
- `10_dev_roadmap_v2.md` — **v2 active roadmap** (V2-1~V2-6 phase definition)
- `v2_current_build_state.md` — build state snapshot. **2026-05-03자, 현재 stale** (V2-1 lock + Phase D/E/F locked 미반영)

#### V2-1 (imaging contract)
- `v2_layer_1_imaging_input_contract_locked.md` — **방금 lock seal** (53eb6c1 → 79376e7 hotfix)
- `v2_phase1_plan_{claude,codex,consolidated}.md`, `v2_phase1_alpha_implementation_brief.md`, `v2_phase1_forward_roadmap.md`
- `v2_p0_morning_summary.md`, `v2_p1_derivation_locked.md`, `v2_p1_active_contour_sanity_gate.md`

#### V2-2 single-cell — focal adhesion
- `v2_focal_adhesion_dynamics_brief.md` (없으면 result_typed_brief 사용)
- `v2_focal_adhesion_dynamics_sanity_gate.md`
- `v2_focal_adhesion_dynamics_result_typed_brief.md`
- `v2_focal_adhesion_dynamics_result_typed_locked.md`
- `v2_focal_adhesion_dynamics_result_typed_sanity_gate.md`

#### V2-2 single-cell — protrusion coupling
- `v2_63b_protrusion_coupling_design_brief.md`
- `v2_63b_protrusion_coupling_locked.md`
- `v2_63b_protrusion_coupled_fa_sanity_gate.md`

#### V2-2 single-cell — adhesion network
- `v2_adhesion_network_layer_locked.md`

#### V2-2 single-cell — Hard Blockers (closed-loop ECM 사전 조건)
- HB#1 (constitutive direction): `v2_hard_blocker_1_constitutive_direction_brief.md`, `v2_hard_blocker_1_2_constitutive_direction_locked.md`, `v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
- HB#3 (FA→ECM 산란): `v2_hard_blocker_3_fa_to_ecm_scattering_brief.md`, `_locked.md`, `v2_fa_to_ecm_scattering_sanity_gate.md`
- HB#4 (active brief + ECM→FA bias): `v2_hard_blocker_4_active_brief.md`, `_locked.md`, `_sanity_gate.md`, `v2_hard_blocker_4_ecm_to_fa_bias_target_brief.md`, `_locked.md`, `v2_ecm_to_fa_bias_sanity_gate.md`, `v2_ecm_ol_sweep_sanity_gate.md`
- HB#5 (Lyapunov metric): `v2_hard_blocker_5_lyapunov_metric_brief.md`, `_locked.md`, `_sanity_gate.md`

#### V2-2 closed-loop ECM 통합
- `v2_closed_loop_ecm_gate_dependency_brief.md`
- `v2_closed_loop_ecm_gate_phased_plan_locked.md`
- `v2_closed_loop_ecm_gate_phase_a_evidence.md`
- `v2_closed_loop_ecm_gate_session_summary.md`

#### V2-2 Phase D / E / F (단일세포 motility 통합)
- Phase D: `v2_phase_d_no_op_scaffolding_brief.md`, `_locked.md`, `_sanity_gate.md`
- Phase E: `v2_phase_e_composition_brief.md`, `_locked.md`, `_sanity_gate.md`, `v2_phase_e_v2_composition_brief.md`, `v2_phase_e_v2_composition_step_2_locked.md`
- Phase F: `v2_phase_f_minimal_cell_motility_pilot_locked.md`, `v2_phase_f_proper_1_fa_rate_response_locked.md`

#### V2-3+ design notes (예약/스케치)
- `v2_cell_cycle_dynamics_design_note.md` — 세포 주기 (V2-2/3 확장)
- `v2_cytokinesis_design_note.md` — 세포분열

#### V2 sweep 도구
- `v2_item_5_sweep_harness_brief.md`, `_locked.md`, `_sanity_gate.md`

#### V2 process 메타
- `v2_sister_gate_mirror_audit_2026_05_05.md` — sister gate audit (어제)

#### `docs/` 루트의 cross-cutting (양 버전 + 협업 인프라)
- `codex_review_synthesis.md` — Codex review 종합
- `implementation_workflow.md` — 구현 워크플로
- `workroom_layout.md` — 워크룸 레이아웃 (collab infra)
- `SESSION_HANDOFF.md`, `SESSION_HANDOFF_archive_v10.md` — 세션 핸드오프
- (`claude_codex_log.md` — Claude/Codex shared ledger; 현재 worktree에는 없을 수 있음)

> **명명 규약**: `_brief` (디자인 스케치) → `_locked` (봉인된 합의) → `_sanity_gate` (실행 전 검증 컨트랙트). 셋이 같이 있으면 lock + execution-ready.

### 3.3 테스트 (`tests/v2/test_v2_*.py`) — 35
| 모듈 | 테스트 |
|---|---|
| state schemas | `test_v2_contract.py`, `test_v2_single_cell_slow_biology.py`, `test_v2_cell_cluster.py`, `test_v2_junction.py`, `test_v2_protrusion.py`, `test_v2_focal_adhesion.py`, `test_v2_ecm_substrate.py`, `test_v2_adhesion_network_dynamics.py` |
| imaging contract | `test_v2_imaging_input_contract.py` |
| dynamics — AC | `test_v2_active_contour.py`, `test_v2_active_contour_dynamics.py`, `test_v2_active_contour_gate.py` |
| dynamics — FA | `test_v2_focal_adhesion_dynamics.py`, `test_v2_protrusion_coupled_focal_adhesion.py` |
| dynamics — ECM | `test_v2_ecm_open_loop.py`, `test_v2_ecm_ol_harness.py`, `test_v2_ecm_ol_sweep.py`, `test_v2_ecm_constitutive_response.py`, `test_v2_ecm_to_fa_bias.py`, `test_v2_ecm_to_fa_bias_active.py`, `test_v2_fa_to_ecm_scattering.py`, `test_v2_ecm_lyapunov_metric.py` |
| closed-loop | `test_v2_closed_loop_phase_d.py`, `test_v2_closed_loop_phase_e.py`, `test_v2_closed_loop_phase_e_sweep.py`, `test_v2_closed_loop_phase_e_sweep_v2.py` |
| phase F | `test_v2_phase_f_minimal_motility.py`, `test_v2_phase_f_minimal_motility_pilot_runner.py` |
| viz/runner | `test_v2_phase_e_v2_pilot_runner.py`, `test_v2_phase_e_v2_pilot_dashboard.py`, `test_v2_phase_e_v2_pilot_frame_replay.py`, `test_v2_phase_e_v2_pilot_movie.py`, `test_v2_stub3d.py`, `test_v2_frame_dump.py` |
| metrics | `test_v2_metrics.py` |

### 3.4 설정 (`configs/`)
- `final_pilot_v2_bare.yaml`, `final_pilot_v2_lam4.yaml`, `final_pilot_v2_pre.yaml` — v2 pilot 3 조건 (Bare/Lam4/Pre)
- `v2_single_cell_data_contract_example.json` — data contract 예시
- `imaging/` — V2-1 imaging contract config (방금 lock된 contract 산출물)

### 3.5 스크립트 (`scripts/`)
- `run_phase_e_v2_pilot.py`, `run_phase_f_minimal_motility_pilot.py` — V2-2 phase runner
- `run_ecm_ol_harness.py`, `run_ecm_ol_sensitivity_sweep.py` — V2-2 ECM open-loop
- `run_p1_alpha_gate.py` — Phase 1 alpha
- `run_protrusion_coupled_fa_harness.py` — protrusion harness
- `render_phase_e_v2_pilot.py`, `movie_phase_e_v2_pilot.py`, `replay_phase_e_v2_pilot.py` — 시각화

---

## 4. v3 — empty

| 항목 | 상태 |
|---|---|
| 코드 (`acs/v3/` 등) | **없음** |
| 문서 (`docs/v3_*.md` 등) | **없음** |
| 테스트 (`tests/test_v3_*.py`) | **없음** |
| 설정 / 스크립트 / configs | **없음** |
| MCP 워크룸 | `v3-design-discussion`, `v3-implementation-work` 두 개 존재. 트래픽 거의 없음 (status에 unread 1건씩) |

> v3 시작 의도가 있으면 별도 결정 필요. 의도 없으면 워크룸 2개를 닫는 것도 정리 옵션.

---

## 5. 데이터 (shared, read-only)

- `data/experimental/260313_{Bare,Lam4,Pre}.csv` — PI 실험 데이터 (2026-03-13). **read-only overlay 전용. 절대 fitting 금지** (CLAUDE.md hard rule)
- `data/experimental/README.md`
- `data/literature/` — 비어있음, 문헌 파라미터 표 보관 예정

> 260313은 spheroid 단위 (`EffectiveRadius_um ≈ 200–265 μm`). **단일 세포 (Layer 1) 검증 데이터 아님** — V2-3/V2-5 territory. 자세한 framing 정리는 V2-2 Unit 1 lock 직전 라운드에서 진행 예정.

---

## 6. 빠른 lookup (자주 묻는 질문)

| 질문 | 답 |
|---|---|
| "어떤 phase에 어떤 코드?" | `docs/v2/10_dev_roadmap_v2.md` + 위 `acs/v2/dynamics/` 표 |
| "이 .md가 design인지 봉인된 합의인지?" | suffix 보기. `_brief`=초안, `_locked`=봉인, `_sanity_gate`=실행 전 검증 |
| "이 모듈의 테스트는?" | `tests/v2/test_v2_<module>.py` 동일 이름 매핑 (v1 모듈은 `tests/v1/`) |
| "V2 pilot 어떻게 실행?" | `scripts/run_phase_e_v2_pilot.py` + `configs/final_pilot_v2_*.yaml` |
| "새 docs 어디에 두지?" | v2 작업이면 `docs/v2/v2_<topic>_<brief|locked|sanity_gate>.md` |
| "PI A/A₀ 비교 측정?" | `acs/analysis/shape_metrics.py:top_down_projection_area` (Hard Rule 11) |
| "GPU 백엔드 변경?" | `ACS_GPU_BACKEND` env var, `acs/gpu.py` |

---

## 7. 알려진 정리 후보 (이 문서 작성 시점)

- `acs/viz/` — 빈 placeholder, 이미 모든 v2 viz는 `acs/v2/viz/`에 있음. 삭제 검토.
- `docs/v2/v2_current_build_state.md` — 2026-05-03자, 현재 stale. V2-1 lock seal + Phase D/E/F locked 반영 필요.
- v3 워크룸 2개 — 사용 의도 없으면 닫기 검토.
- 폴더 자체를 v1/v2/v3 top-level로 재배치 (옵션 C) — V2-2 Unit 1 lock 봉인 후 별도 schedule.

이 문서는 정적 매핑이라, 새 파일을 추가하면 같이 갱신할 책임은 작업자에게.
