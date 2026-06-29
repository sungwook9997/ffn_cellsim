# ffn_cellsim 파일 구조 지도

작성: 2026-05-20 KST. `~/ActiveCellSim_v2` → `~/ffn_cellsim` rename + 대규모 cleanup 직후 스냅샷.

이 문서는 **파일 → 역할** 매핑 인덱스. 새 세션에서 어디 무엇이 있는지 빠르게 파악할 때 사용.

> ⚠️ **2026-06-29 — TWO-LAYER RESTRUCTURE.** 아래 본문(2026-05-20 작성)은 Phase C 이전이라 최상위 레이아웃 기준으로 **historical**. 현재 활성 레이아웃:
>
> | dir | 역할 |
> |---|---|
> | `ffn_sim/ff/` | **Filament-FEM 엔진** (Cytosim 물리 — Nédélec & Foethke 2007, Warp). 빌드 진행 중 (γ-floor 프로토타입). |
> | `ffn_sim/dcm/` | **Deformable Cell Model 엔진** (SimuCell3D 물리, Warp). `warp_port/`에서 rename; hoomd-free 런타임. |
> | `ffn_sim/common/` | engine-agnostic 공유 (surface_manifold, filament_math, gsd_traj, …). |
> | `ffn_sim/validation/oracles/` | closed-form oracle (불변, runtime-import-forbidden). |
> | `ffn_sim/archive/hoomd_legacy/` | **퇴역 HOOMD 런타임** — cell, cortex, bridge, ecm, junction, integrator, spheroid, gpu_opt. frozen 참조 / parity oracle. |
> | `ffn_sim/{scripts,tests,docs,outputs,configs}/` | 최상위 유지 (hoomd 의존 항목은 archive 경로에서 import). |
>
> 브랜치: `dcm/main`·`ff/main`·`ffn/foundation`(보호). 과거 14개 → `archive/*` 태그. 안전 태그 `pre-restructure-2026-06-29`.

---

## 0. 최상위 (`~/ffn_cellsim/`)

| 항목 | 역할 |
|---|---|
| `ffn_sim/` | active HOOMD-blue runtime + 부속 (docs, tests, scripts, outputs, validation) |
| `README.md` | 사람 대상 프로젝트 overview |
| `CLAUDE.md` | Claude Code 세션 컨텍스트 + 아키텍처 원칙 + hard rules |
| `STRUCTURE.md` | 이 파일 |
| `pyproject.toml` | 빌드 + 의존성 + ruff 설정 |
| `requirements.txt` / `requirements-dev.txt` / `requirements.lock.txt` | HOOMD-blue 스택 의존성 |
| `.gitignore` / `.stignore` | git + Syncthing ignore 패턴 |

git remote `origin` → `~/ActiveCellSim` (local path remote; v1 + image-constrained-v2 시절 자산이 이쪽에 보존되어 있음).

---

## 1. `ffn_sim/` — active HOOMD runtime

### 1.1 runtime subpackages (Plan v2 H.1~H.7 매핑)

| 폴더 | Plan 단원 | 역할 | 현재 상태 |
|---|---|---|---|
| `ffn_sim/ecm/` | H.1 | Mikado fiber network + cross-link HOOMD topology | skeleton (`__init__.py`만) |
| `ffn_sim/cell/` | H.3 | 다중 필라멘트 cell body composition | skeleton |
| `ffn_sim/cortex/` | H.2, H.3 | actin cortex (cortex-only filament, cortical band) | skeleton |
| `ffn_sim/bridge/` | H.4 | focal adhesion + Stam-Hocky multi-head motor-clutch | skeleton |
| `ffn_sim/junction/` | H.6 | E-cadherin cell-cell adhesion (full catch-bond KU-4.2) | skeleton |
| `ffn_sim/integrator/` | H.0.4 | Leimkuhler-Matthews BAOAB-limit custom plugin | skeleton (Phase 0 산출물에 알고리즘 정의 있음) |
| `ffn_sim/common/` | shared | derived params, sanity-gate helpers, units, logging | skeleton |
| `ffn_sim/__init__.py` | — | 패키지 루트 | 비어있음 |
| `ffn_sim/warp_port/` | Phase C | NVIDIA Warp GPU-resident 미분가능 DCM 엔진 (HOOMD 대비 piece별 parity-gated; `engine.py`/`ENGINE.md`) | adopted 2026-06-21 (G1 GPU parity, G2 substrate) |

### 1.2 validation/

| 폴더 | 역할 |
|---|---|
| `ffn_sim/validation/` | runtime cross-check infra (oracle 호출, 회귀 테스트) — skeleton |
| `ffn_sim/validation/oracles/` | **v1에서 살아남은 closed-form oracles**. 아래 1.3 참조. |

### 1.3 `ffn_sim/validation/oracles/` — v1 acs_kb에서 carve해온 oracle (runtime-import-forbidden)

`acs_kb/_DEPRECATED.md`가 분류한 "v2 runtime survivor" 목록. 모두 closed-form 수학/지오메트리 / KU 상수만 — paper-model-as-mechanism wrapper는 전부 archive에서 제외됨.

| 파일 | 역할 |
|---|---|
| `oracles/__init__.py` | 패키지 docstring (origin: v1 acs_kb 2026-05-20 carve) |
| `oracles/ecm/fiber_network.py` | Mikado geometry generator (HOOMD topology seed) |
| `oracles/ecm/cross_links.py` | segment-intersection cross-link seeding; harmonic kernel oracle |
| `oracles/ecm/fiber_mechanics.py` | KU-1.24 WLC discrete H formula (numpy oracle for `md.bond.Harmonic`/`md.angle.Harmonic`) |
| `oracles/ecm/diagnostics.py` | integrator-agnostic measurement utility |
| `oracles/junction/contact_angle.py` | Maître Young-equation KU-4.4 contact angle |
| `oracles/junction/types.py` | `EcadherinJunction` dataclass (data container; Cell 타입은 `Any`로 stub) |
| `oracles/bridge/traction.py` | 1-D analysis reducer for FA force |
| `oracles/bridge/types.py` | `FocalAdhesion` dataclass |
| `oracles/common/sanity_gate.py` | KU sanity-gate helpers |
| `oracles/common/derived_params.py` | analytical derived parameters (Stokes drag, Mikado ℓ_c, τ_min) |
| `oracles/common/derived_params_cell.py` | cell-side derived parameters |
| `oracles/configs/phase1_unit*.yaml` | KU-anchored literature constants (5 YAMLs: unit1, 2_1, 2_2, 3, 4_1) |

**원칙**: 어떤 `ffn_sim/<runtime>/` 모듈도 `ffn_sim.validation.oracles.*` 를 import하면 안 됨. oracle은 오로지 `ffn_sim/tests/` + `ffn_sim/validation/` 에서만 호출됨.

### 1.4 `ffn_sim/docs/` — 활성 프로젝트 문서

| 파일 | 역할 |
|---|---|
| `PHASE_0_CLOSEOUT.md` | Phase 0 종료 보고 (2026-05-19) |
| `PHASE_0_3_DECISIONS.md` | PI 비준 7대 디자인 결정 (Bell-Evans, Hill, Stam-Hocky, KU-4.2, 각도 harmonic, Leimkuhler-Matthews, LJ-EV) |
| `AFINES_ALGORITHM_NOTES.md` | AFINES 마스터 리뷰 (897 줄) |
| `briefs/H1_ecm_mikado.md` | Worker A 디스패치 brief |
| `briefs/H2_single_filament.md` | Worker C 1단계 |
| `briefs/H3_cortex.md` | Worker C 2단계 |
| `briefs/H4_fa_motor_clutch.md` | Worker B 디스패치 brief |
| `v2_audit/CODEBASE_AUDIT.md` | v1 acs_kb → v2 reuse/port/archive 분류 (45% archive 결정 근거) |

### 1.5 `ffn_sim/tests/`, `ffn_sim/scripts/`, `ffn_sim/outputs/`

| 폴더 | 내용 |
|---|---|
| `ffn_sim/tests/` | HOOMD runtime 테스트 (현재 `__init__.py`만) |
| `ffn_sim/scripts/` | `hoomd_polymer_sanity.py` — M1 Max 벤치마크 재현 (89k steps/s) |
| `ffn_sim/outputs/h_0_2/` | Phase 0.2 HOOMD env 검증 REPORT.md |

---

## 2. 알려진 정리 후보 (2026-05-20)

- `ffn_sim/validation/oracles/junction/types.py` — Cell 타입 annotation을 `Any`로 stub함. `make_ecadherin_junction()` 헬퍼는 v1 Cell 메서드(`compute_cortex_boundary_position`)에 의존하므로 v2에서는 호출하지 말 것. 필요하면 v2 Cell 객체로 재작성.
- `oracles/configs/phase1_unit*.yaml` 안의 주석 일부에 v1 경로 (`acs_kb/outputs/...`) 흔적이 남아있음. 사용에는 지장 없음.
- `requirements.lock.txt` — v1 Taichi 스택 pin. HOOMD 스택용으로 재생성 필요 (Phase 1 시작 시 `pip freeze` 다시 캡처).

---

## 3. 빠른 lookup

| 질문 | 답 |
|---|---|
| "이 단원 무슨 코드?" | `ffn_sim/docs/briefs/H<n>_*.md` |
| "왜 이 디자인이지?" | `ffn_sim/docs/PHASE_0_3_DECISIONS.md` |
| "v1에서 뭐가 살아남았지?" | `ffn_sim/docs/v2_audit/CODEBASE_AUDIT.md` |
| "oracle 수식이 어디 있지?" | `ffn_sim/validation/oracles/<subpkg>/` |
| "벤치마크 어떻게 돌리지?" | `python ffn_sim/scripts/hoomd_polymer_sanity.py --steps 50000 --bench` |
| "v1 코드 어디?" | `~/ActiveCellSim` (이 repo의 git origin), `ffn_cellsim` 안에는 없음 |

이 문서는 정적 매핑이라, 새 파일을 추가하면 같이 갱신하는 책임은 작업자에게.
