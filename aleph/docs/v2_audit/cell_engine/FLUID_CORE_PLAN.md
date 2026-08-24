# Cell Engine — Fluid Volume + Core Body 실행 계획

> **Production-baseline override (2026-07-22):** native lamina/chromatin state remains authoritative for the
> first full-population run. The 24–64-mode Core below is an optimisation candidate and may become a runtime
> backend only after its exact-volume, LINC response, force/work, and frequency-response gates close against
> the native nucleus at the physiological operating point.

Status: **Wave-1 implementation contract — 2026-07-22**
Scope owner: `cytosol` (`FLUID_VOLUME`) + `nucleus` (`CORE_BODY`) + 두 컴포넌트의 field/boundary coupling
Common contract: `CELL_ENGINE_ARCHITECTURE.md`, `aleph/engine/contracts.py`
Goal: 저 Reynolds 수 세포질을 값싼 porous field로, 핵을 값싼 reduced deformable body로 풀되 LINC·압력·기계화학 피드백의 양방향 동역학은 보존한다.

## 1. 확정 결정

1. **Cytosol 기본식은 Biot storage + Darcy flux다.** 관성 Navier–Stokes나 전역 LBM은 사용하지 않는다.
   저 Reynolds 수만으로 Darcy가 도출되는 것은 아니며, cytosol을 cytoskeletal porous medium으로
   평균화한다는 두 번째 가정이 함께 들어간다.
2. **Brinkman은 전역 기본식이 아니다.** `ell_B=sqrt(k)`의 shear-screening layer 또는 실제
   free-fluid pocket이 격자에서 분해되고, Darcy-only 기준모델의 경계 traction/drag 오차가 등록된
   허용범위를 넘는 국소 블록에서만 활성화한다.
3. **핵은 642개 표면 노드를 전역 미지수로 풀지 않는다.** 동일한 642-node mesh는 geometry,
   collision, pressure quadrature, strain visualization용으로 유지하되 mechanics는 rigid pose +
   validated deformation basis의 generalized coordinates로 푼다.
4. **핵의 부피는 큰 penalty stiffness가 아니라 constraint로 보존한다.** 이것이 현재 `k_vol`이
   전체 pseudo-time step을 제한하는 경로를 제거한다.
5. **LINC는 핵 내부 spring 배열이 아니라 connector graph의 joint다.** 핵은 LINC socket geometry와
   force/velocity projection만 제공한다. 상대 endpoint는 SF/actin cap, MT, IF가 소유하고 결합·해리는
   `ConnectorFamily.LINC`가 소유한다.
6. **Fluid/Core는 독립 solver여도 하나의 candidate transaction 안에서 반복 수렴한다.** pressure가
   핵/표면을 움직이고, 움직인 경계가 같은 candidate의 pressure와 connector geometry를 다시 바꿔야
   한다. 한 번씩 순서대로 호출하고 끝내는 one-way split은 허용하지 않는다.

생산상 의미는 다음과 같다.

```text
Fluid Volume = pressure/flux/mass 전달 매질
Core Body    = LINC·압력 하중을 받는 reduced deformable actor
LINC graph   = SF/MT/IF에서 핵까지 이어지는 동적 joint
World solver = 이 셋과 Surface Body를 같은 accepted physical clock에서 수렴시키는 소유자
```

## 2. 현재 코드 감사 결과

2026-07-22에 다음 범위를 읽고 선택 회귀를 실행했다.

- `aleph/components/fluid/`, `aleph/components/nucleus/`
- `aleph/components/incumbent/{assemble,driver,compartments,fsi_coupling,live_mesh_domain,membrane_pressure}.py`
- 대응 `aleph/tests/ac/{fluid,nucleus,cell}`과 `outputs/ac/{fluid-spine,nucleus,fsi,foundation-hardening}`
- 결과: **122 PASS, 1 CUDA-only SKIP**

### 2.1 이미 좋은 부분

- `FieldGrid`와 `biot_pmass_update_kernel`은 masked conservative FV, double buffer, device-resident
  content reduction을 제공한다.
- `boundary.py`의 live-triangle Kedem–Katchalsky quadrature는 실제 triangle area로 flux를 보존적으로
  grid에 deposit한다.
- `live_mesh_domain.py`는 membrane/nucleus `wp.Mesh`를 device에서 refit하고 host readback 없이
  domain을 다시 분류한다.
- `Domain.remap()`은 moving-face content를 별도 ledger로 기록한다.
- `PressureCoupling`과 `SolidDilatationCoupling`은 Peskin transfer, partition of unity, work/conservation
  시험의 기반을 갖고 있다.
- `PhysicalScheduler`는 device acceptance latch, rollback, post-accept irreversible commit, zero-D2H
  transaction을 이미 구현한다.
- nucleus 쪽은 Helfrich/lamina/volume/LINC/chromatin kernel과 독립 analytic oracle, live mask provider를
  갖고 있다. 이들은 reduced basis 생성과 reference backend에 재사용할 수 있다.

### 2.2 구조적으로 잘못 묶인 부분

| 현재 상태 | 실제 영향 | target |
|---|---|---|
| 모든 위치가 `AssembledCell.pos_d` 한 배열에 있음 | 소유권과 coupling이 구분되지 않음 | Fluid/Core가 각자 state를 소유하고 connector/coupler만 교환 |
| `dx=0.5 um`, `dt_phys=0.05 s`에서 explicit Biot **67 subcycle** | 작은 43^3 grid에 불필요한 kernel launch 반복 | backward-Euler matrix-free PCG + geometric multigrid |
| 79,507 cells를 할당하지만 현재 live cytosol active cells는 약 9.8k | 메모리 문제는 아니나 전체-grid launch가 많음 | 첫 slice는 dense 유지, stencil fuse; 큰 domain에서만 active-tile compaction |
| `node_volume`은 actin에만 부여되고 nucleus/membrane은 0 | bulk `PressureCoupling`과 `div(v_s)`가 사실상 actin-only | coupling 종류별 surface/volume measure를 명시 |
| nucleus는 fluid mask의 no-flux hole일 뿐 직접 pressure traction이 없음 | pressure가 핵을 직접 변형하지 못함 | nucleus surface pressure quadrature ↔ reduced generalized force |
| cortex node velocity spread와 moving-domain remap을 별도 합산 | ALE boundary motion과 skeleton dilatation의 중복 위험 | geometric conservation과 Biot skeleton source를 별도 typed channel로 정의 |
| scheduler는 이전 accepted `div(v_s)`로 fluid를 풀고 mechanics 뒤 다음 값을 만듦 | 1-step-lagged first-order stagger | 같은 candidate 안에서 pressure/geometry interface iteration |
| nucleus 642 nodes = 1,926 global positional DOF | 수 자체는 작지만 global stiff solve와 `k_vol`에 편입 | 6 pose + 보통 24–64 deformation coordinates |
| LINC 642개를 nearest cortex node에 정적으로 생성 | SF/MT/IF 연결성이나 kinetics가 없음 | connector-owned dynamic sockets/endpoints |
| `NucleusCompartment.accumulate()`이 chromatin WLC와 nucleoplasm viscosity를 호출하지 않음 | 문서와 실제 composed physics 불일치 | ROM constitutive calibration에 포함하거나 명시적으로 격리 |
| nucleus rupture는 per-face state이나 reduced kinematics와 계약이 없음 | local tear를 ROM으로 억지 표현할 위험 | everyday regime 밖이면 local/full-shell backend로 promote |

현재 full-native 기록은 70,686 cortex fibers, 511,114 total nodes, 79,507 allocated field cells이며
0.05 s outer attempt가 A5000 laptop에서 1.1307 s였다. 이 중 큰 병목은 mechanics convergence지만,
Fluid/Core는 그 병목에 stiffness와 반복 launch를 더하지 않는 형태로 바꿔야 한다.

## 3. 코드 disposition

### 3.1 그대로 또는 얇은 adapter로 재사용

| Existing code | Disposition |
|---|---|
| `fluid/field_grid.py` | field allocation, mask codes, reduction을 `FluidVolumeState` 내부로 이동 |
| `fluid/biot_substrate.py` | conservative face flux와 analytic parity는 유지; explicit step은 oracle/backend로 격리 |
| `fluid/{fv_reference,manufactured,consolidation_analytic,greens_analytic,darcy_analytic}.py` | acceptance oracle로 그대로 유지 |
| `fluid/boundary.py`, `fluid/surface_trace.py` | membrane hydraulic BC와 pressure trace law 재사용 |
| `fluid/velocity.py` | `q`와 `v_f` derived-field kernel 재사용; 매 mechanics iteration마다 계산하지 않음 |
| `cell/live_mesh_domain.py` | live surface query와 D2D refit 재사용 |
| `fluid/scheduler.py` | snapshot/rollback/device latch를 world transaction으로 일반화 |
| `nucleus/{geometry,lamina_analytic,chromatin_analytic,linc_analytic}.py` | high-fidelity reference와 ROM basis/response oracle |
| `nucleus/envelope.py` | full-shell reference backend 및 local promotion backend |

### 3.2 migration 동안만 격리

- `cell/assemble.py`의 global node concatenation
- `cell/driver.py::_accumulate_all`와 현재 monolithic inner solve
- `cell/compartments.py::{NucleusCompartment,MeshNucleusMaskProvider}`
- `cell/fsi_coupling.py::SolidDilatationCoupling`

이들은 A/B regression adapter로 남기되 새 컴포넌트가 import하거나 private field를 공유하지 않는다.

### 3.3 새 경로에서 폐기

- nucleus↔nearest-cortex 정적 LINC pairing
- full shell의 global positional DOF 등록
- volume penalty를 global `kmax`에 넣는 방식
- 모든 solid node에 하나의 임의 `node_volume`을 붙여 pressure/divergence를 처리하는 방식
- fixed-count explicit Biot subcycle을 production 기본 solver로 쓰는 방식
- pressure field를 푼 뒤 한 번만 mechanics를 호출하는 one-way outer split

## 4. Fluid Volume state와 DOF

### 4.1 authoritative state

```text
p_excess[cell]            spatial pore-pressure deviation [Pa]
p_bar                     mean osmotic/hydrostatic channel [Pa]
fluid_fraction[cell]      cut-cell volume fraction
face_aperture[face]       open Darcy face fraction
material_id[cell]         outside / cytosol / nucleus / optional Brinkman patch
water_source[cell]        non-boundary source [1/s]
boundary_flux[face/quad]  membrane hydraulic flux [um/s]
candidate/accepted buffers and conservation ledger
```

`q`, `v_f`, transported concentrations는 필요한 컴포넌트가 있을 때 derived/cache state로 만든다.
현재 43^3 box에서는 scalar pressure unknown이 최대 79,507개이고 실제 cytosol unknown은 약 9.8k다.
따라서 첫 최적화는 AMR가 아니라 **implicit solve + launch fusion**이다.

### 4.2 governing operator

Moving finite-volume cell `c`에 대해 다음 conservative form을 사용한다.

```text
d/dt integral_Omega_c(t) S p dV
  + alpha integral_Omega_c div(v_s) dV
  + sum_faces q_f A_f
  = integral_Omega_c s_water dV

q_f = -(k/mu)_f (grad p - rho b)_f
```

- cut-cell volume/aperture 변화는 discrete geometric conservation law가 소유한다.
- porous skeleton dilatation은 `PorousSkeletonMap`이 소유한다.
- 같은 cortex motion을 moving-face와 `div(v_s)`에 두 번 넣지 않는다.
- membrane water permeation은 `FLUID_BOUNDARY` connector가 소유한다.
- nuclear envelope는 moving impermeable boundary다. 상대 normal flux는 0이며 pressure traction과
  boundary velocity work는 adjoint pair로 계산한다.

### 4.3 solver

첫 production solver는 backward Euler의 SPD pressure system이다.

```text
(S V / dt + D(K, aperture)) p_(n+1) = rhs(p_n, geometry, sources)
```

- matrix-free PCG, geometric multigrid V-cycle preconditioner
- mask/cut-cell coefficient가 변하지 않은 동안 hierarchy 재사용
- Warp CUDA graph로 fixed launch sequence capture
- residual, content balance, finite-state를 device scalar로 판정
- pure-Neumann/mean channel은 `p_bar`와 분리하여 null-space를 명시적으로 다룸

현재 67 explicit subcycle은 reference backend에 남기고, implicit 결과를 같은 Terzaghi/Green/
manufactured gates에 대조한다. 압력 변화가 매우 작으면 이전 `p`를 predictor로 쓰되 solve 자체를
frame-rate 기준으로 생략하지 않는다.

### 4.4 optional local Brinkman patch

Brinkman patch는 다음 식의 vector velocity unknown을 국소 block에서만 추가한다.

```text
-grad p + mu_eff laplacian(v_f) - (mu/k)(v_f-v_s) + f = 0
div(v_f) = prescribed storage/boundary rate
```

활성 조건은 `ell_B=sqrt(k)`가 grid에서 분해되고 Darcy-only traction/drag error estimator가 등록
허용범위를 넘는 경우다. Darcy–Brinkman interface는 normal flux와 traction을 연속시킨다. Wave 1에는
patch API와 Darcy limit test만 만들고 실제 patch 활성은 요구하지 않는다.

## 5. Core Body state와 DOF

### 5.1 reduced kinematics

```text
pose = translation(3) + rotation(3)
q_r[r]                    deformation generalized coordinates
lambda_V                  exact volume constraint multiplier
x_surface = T(pose) [X0 + Phi q_r + nonlinear_correction(q_r)]
v_surface = J(q_r) v_r
```

초기 `r`은 고정 숫자가 아니라 full-shell reference의 response convergence로 정한다. 구현 시작값은
24–64 modes이며 포함 대상은 axial flattening, 두 equatorial modes, shear, bending, low-order surface
waviness, LINC-interface modes다. Craig–Bampton/component-mode synthesis로 active LINC socket의 경계
응답을 보존하고 interior shell/chromatin DOF를 static-condense한다.

642-node mesh는 계속 매 candidate interface iteration에 GPU에서 reconstruct한다. 그러므로 fluid
mask, collision, LINC barycentric point, strain map과 rendering의 공간해상도는 잃지 않으며 global
mechanical unknown만 약 1,926에서 수십 개로 줄어든다.

### 5.2 constitutive response

- lamina bending와 areal strain-stiffening은 full-shell energy를 reduced basis에 projection한다.
- chromatin WLC/net response와 nucleoplasm viscosity는 basis 생성 및 reduced damping에 포함한다.
- exact volume constraint는 one scalar Schur complement로 푼다.
- small/moderate everyday deformation은 corotational reduced Newton으로 푼다.
- ROM residual estimator가 band를 벗어나거나 rupture/local indentation이 발생하면 해당 surface patch
  또는 full shell backend로 promote한다. rupture를 저차 mode에 억지로 smearing하지 않는다.

### 5.3 dynamic sockets와 mechanochemical outputs

Core Body가 제공하는 socket은 `(surface face, barycentric coordinate, local frame)`이다. LINC connector가
이를 SF/MT/IF endpoint와 묶고 force law, binding epoch, kinetics를 소유한다. Core는 socket load를
`J_socket^T f`로 generalized residual에 더하고 socket velocity를 `J_socket v_r`로 반환한다.

Core가 매 accepted step에 노출할 read-only response는 다음이다.

- per-socket LINC force/exposure
- nuclear surface principal strain과 lamina tension
- chromatin extension/strain proxy
- volume, aspect ratio, curvature, damage indicator
- mechanochemical feedback input vector와 epoch

실제 YAP/TAZ 또는 transcription model이 나중에 연결되더라도 geometry를 사후 측정하는 방식이 아니라
같은 accepted state의 load history를 읽게 한다.

## 6. runtime API

`contracts.py`는 build-time graph를 검증한다. 다음은 그 위에 추가할 device-runtime protocol이다.

```python
class FluidVolume:
    state: FluidVolumeState
    def snapshot_candidate(self, accepted_latch): ...
    def update_geometry(self, boundary_views): ...
    def solve_candidate(self, dt_phys, porous_sources, boundary_fluxes): ...
    def pressure_trace(self, surface_view) -> FieldTrace: ...
    def accumulate_boundary_traction(self, surface_view, out_surface_force): ...
    def derived_velocity(self, request_mask) -> VelocityView: ...
    def residual_view(self) -> DeviceResidual: ...
    def rollback_or_commit(self, accepted_latch): ...

class CoreBody:
    state: CoreBodyState
    def geometry_view(self) -> SurfaceView: ...
    def socket_view(self, family) -> SocketView: ...
    def project_surface_force(self, force_view, out_generalized_force): ...
    def accumulate_internal_residual(self, out_generalized_force): ...
    def solve_reduced_candidate(self, dt_phys): ...
    def feedback_view(self) -> MechanochemicalView: ...
    def rollback_or_commit(self, accepted_latch): ...
```

필수 algebraic contract는 다음 두 개다.

```text
surface interpolation:  v_surface = J v_generalized
force projection:       f_generalized = J^T f_surface
work identity:           f_surface dot v_surface = f_generalized dot v_generalized
```

Fluid boundary도 동일하게 pressure/velocity interpolation과 traction/source scatter가 transpose pair여야
한다. 각 컴포넌트는 다른 컴포넌트의 private kernel이나 array offset을 알지 못한다.

## 7. 양방향 coupling과 candidate 순서

한 candidate outer step의 Fluid/Core 부분은 다음 fixed-point 또는 block-Newton loop다.

1. accepted state를 device-to-device snapshot한다.
2. Surface/Core geometry와 connector socket을 현재 candidate로 reconstruct한다.
3. membrane hydraulic flux, nuclear no-flux kinematics, porous skeleton source를 만든다.
4. implicit Biot/Darcy pressure를 푼다.
5. membrane/nucleus pressure traction과 optional slender-body fluid drag를 계산한다.
6. Surface/Core/SF/FA/ECM coupled mechanics를 푼다. Core surface motion과 LINC reaction을 갱신한다.
7. geometry/cut-cell coefficients와 connector capture geometry를 갱신한다.
8. pressure residual, interface displacement, force/work, mass ledger가 모두 수렴할 때까지 3–7을 반복한다.
9. global validators가 통과하면 physical clock, connector kinetics, damage/feedback epoch를 commit한다.
10. 하나라도 실패하면 fields, reduced coordinates, geometry, connector states와 time을 bit-exact rollback한다.

첫 slice는 구현이 단순한 partitioned Picard + Aitken acceleration로 시작한다. 수렴 속도가 부족한
경우 pressure/Core reduced block은 Schur complement로 묶되 전체 cell을 거대한 monolithic matrix로
만들지 않는다.

### coupling ownership

| Coupling | Forward | Back reaction | Owner |
|---|---|---|---|
| membrane↔cytosol | surface velocity/geometry, osmotic data | pressure traction, hydraulic volume change | `FLUID_BOUNDARY` |
| cortex↔cytosol | porous skeleton kinematics | `-alpha p I` work | `IMMERSED_TRANSFER` |
| nucleus↔cytosol | moving impermeable surface | pressure traction projected to Core DOF | fluid boundary/immersed transfer adapter |
| SF/MT/IF↔cytosol | line velocity and drag support | relative drag; optional Darcy body source | later `IMMERSED_TRANSFER` |
| SF/MT/IF↔nucleus | connector endpoint motion | equal/opposite LINC force | `LINC` |

## 8. 첫 vertical slice

목표 시나리오는 **motor-driven SF load가 LINC를 통해 reduced nucleus를 변형하고, 동시에 membrane/
cytosol pressure가 그 응답을 바꾸며, 역으로 움직인 핵 경계가 같은 candidate pressure를 바꾸는 것**이다.

### Slice FC-1 — pressure-driven Core

1. 기존 43^3 live membrane-minus-nucleus domain을 `FluidVolume` adapter로 감싼다.
2. existing 642-node nucleus reference에서 reduced basis를 만들고 `CoreBody`를 생성한다.
3. nucleus pressure surface quadrature와 `J^T` projection을 구현한다.
4. prescribed membrane compression/pulse로 spatial pressure를 만들고 nucleus가 quadrupolar하게
   flatten되도록 coupled Picard loop를 닫는다.
5. nucleus motion이 mask/cut-cell content에 즉시 반영되고 relative no-flux/content ledger가 닫히는지 본다.

### Slice FC-2 — one live LINC path

1. nearest-cortex LINC 자동생성을 끈다.
2. Load-path track의 한 `sf_arc` endpoint와 Core surface socket 사이에 동적 LINC 하나를 등록한다.
3. NMII/SF load → LINC → Core deformation → socket displacement → LINC load 변화가 한 candidate에서
   수렴하게 한다.
4. LINC ablation 시 pressure-only nuclear response만 남는지 확인한다.
5. reject를 강제하여 pressure, Core `q_r`, reconstructed mesh, LINC state, feedback epoch, time이 모두
   bit-exact 복구되는지 확인한다.

이 두 slice가 통과하기 전에는 Brinkman, multi-species transport, rupture, AMR를 구현하지 않는다.

## 9. acceptance tests

### Fluid unit gates

- uniform pressure/zero source exact preservation
- impermeable content conservation과 membrane flux==content change
- Terzaghi, Green, manufactured moving-boundary convergence
- implicit vs current explicit oracle agreement; time/space refinement order
- cut-cell geometric conservation under translating and volume-preserving deforming boundaries
- Darcy slab linearity, `v_f-v_s=q/phi`
- Brinkman patch의 `k→0` Darcy limit, patch-interface flux/traction continuity
- no host readback and no per-step allocation

### Core unit gates

- rigid translation/rotation에서 internal energy와 spurious force 0
- full-shell 대비 uniform pressure, dipole, quadrupole, single/multiple LINC point-load response
- static displacement, strain energy, surface traction, aspect/volume error가 declared band 이내
- frequency-response band에서 amplitude/phase parity
- exact volume constraint residual
- `J/J^T` work identity와 connector endpoint permutation invariance
- ROM residual estimator가 의도한 out-of-band case를 promotion으로 보냄

### Coupled gates

- pressure traction↔moving-boundary work equality
- moving nucleus relative no-flux 및 content closure
- closed system net force/torque와 COM drift
- SF motor→LINC→nucleus load path; connector ablation locality
- pressure pulse→nuclear deformation→updated pressure의 same-candidate convergence
- rejected step의 bit-exact full transaction rollback
- component permutation과 solver partitioning에 대한 accepted state invariance
- A5000/4090 profiler에서 phase별 time, iteration, bytes ledger

기존 122-pass selection은 migration 내내 compatibility baseline으로 유지한다.

## 10. migration 순서와 파일 경계

제안 파일은 한 개념당 하나로 나눈다.

```text
ac/engine/fluid_volume.py             state + component facade
ac/engine/porous_pressure_solver.py   implicit SPD operator/PCG/MG
ac/engine/fluid_boundary.py           membrane/nucleus boundary connectors
ac/engine/brinkman_patch.py           optional local patch interface
ac/engine/core_body.py                reduced state/solve/geometry reconstruction
ac/engine/core_basis.py               full-shell -> ROM basis and validation metadata
ac/engine/component_runtime.py        runtime protocols/views only
tests/ac/engine/test_fluid_volume_*.py
tests/ac/engine/test_core_body_*.py
tests/ac/engine/test_fluid_core_slice.py
```

1. runtime views/protocol과 compatibility adapters
2. Core basis generator + full-shell response harness
3. reduced Core solve + exact-volume constraint
4. nucleus fluid-boundary pressure/work coupler
5. implicit pressure solver; existing explicit backend과 parity
6. coupled Picard/Aitken candidate loop
7. FC-1 gate
8. dynamic LINC connector 연결과 FC-2 gate
9. old driver의 nucleus/static-LINC/explicit-fluid path를 compatibility-only로 표기
10. profiling 후에만 active tiles, CUDA graph fusion, optional Brinkman 구현

Fluid/Core 작업은 각각 내부 파일을 병렬로 만들 수 있지만, `component_runtime.py`, connector schema,
world transaction은 Lead가 직렬 통합한다. 두 작업자가 같은 global assembler를 동시에 수정하지 않는다.

## 11. 성능 목표와 예상

이 값은 acceptance budget이며 사전 측정 없는 성능 주장이 아니다.

| 항목 | 현재 | Wave-1 목표 |
|---|---:|---:|
| pressure advance / 0.05 s | 67 explicit subcycles | 1 implicit solve per interface iteration |
| field kernel work | 67 stencil passes + source launches | 2–8 PCG/MG effective passes를 목표로 측정 |
| nucleus global mechanics DOF | 1,926 positional | 6 pose + 24–64 deformation + 1 volume multiplier |
| nucleus volume enforcement | stiff `k_vol` penalty | exact small Schur constraint |
| Fluid+Core A5000 step budget | 미분리 | **<=5 ms** at current 43^3 grid and validated ROM |
| Fluid+Core 4090 desktop budget | 미분리 | **약 2–3 ms 목표**, 반드시 profiler로 확정 |
| fluid/core device memory | 전체 cell ledger에 혼재 | component별 exact bytes, 첫 slice **<64 MB** budget |

Implicit pressure로 field launch work는 대략 8–30배, nucleus global DOF는 약 30–70배 줄 가능성이 있다.
전체 cell step의 동일 배속을 뜻하지는 않는다. 현재 1.13 s outer attempt의 주 병목은 cortex/whole-cell
mechanics이므로 Fluid/Core의 성공 기준은 (a) 자체가 millisecond budget을 지키고, (b) stiff penalty와
global node coupling을 제거하여 world solver 반복을 악화시키지 않는 것이다.

## 12. 공통 계약 대조와 수정 제안

`CELL_ENGINE_ARCHITECTURE.md`의 component/connector/transaction/LOD 정의와 이 계획은 일치한다.
`contracts.py`의 `cytosol: FLUID_VOLUME`, `nucleus: CORE_BODY`, `LINC`, `IMMERSED_TRANSFER`, 새로 추가된
`FLUID_BOUNDARY`도 그대로 사용한다. 다음은 구현 전 명확히 해야 할 계약 차이이다.

1. **Lead 통합 완료:** 핵–cytosol은 `nucleus_cytosol_boundary: FLUID_BOUNDARY` 하나로 확정됐고 porous
   skeleton transfer와 구분된다.
2. **Lead graph 통합 완료, runtime dispatch 남음:** SF/MT/IF↔cytosol `IMMERSED_TRANSFER` edge는 모두
   reference graph에 등록됐다. 이제 각 solid facade와 public cytosol residual endpoint를 실제 coupler에
   바인딩해야 한다.
3. **build-time contract만으로 rollback/work API를 검증할 수 없다.** 위 `FluidVolume`/`CoreBody`
   runtime protocol과 `StateView`, `GeometryView`, `ResidualView`를 별도 module로 추가해야 한다.
4. **`owns_geometry=True`의 cytosol 의미가 불명확하다.** cytosol은 field-grid geometry는 소유하지만
   membrane/nucleus boundary geometry는 소비한다. 이후 contract에는 `owned_geometry`와
   `boundary_dependencies`를 분리하는 편이 안전하다.
5. **비kinetic connector의 `commit_on_accept=False`는 state rollback 불필요를 뜻하면 안 된다.**
   FSI/cut-cell cache도 transaction state다. 이 flag는 irreversible kinetic commit에만 적용된다고
   docstring을 강화해야 한다.

이 다섯 항목은 물리 설계 변경이 아니라 component ownership과 double-count를 막기 위한 API 정리다.

## 13. 완료 정의

Fluid/Core의 native-first wave는 다음이 모두 참일 때 끝난다.

- cytosol이 implicit Biot/Darcy component로 독립 state와 ledger를 소유한다.
- nucleus가 native lamina/chromatin backend로 움직이며 압력/LINC/volume gate를 통과한다.
- membrane/cortex/nucleus pressure·flux transfer가 양방향 work/mass gate를 닫는다.
- 최소 한 SF LINC path가 connector graph를 통해 핵에 동적 하중을 전달한다.
- 핵 변형이 같은 candidate의 fluid domain과 pressure에 되먹임된다.
- 실패한 candidate가 fields, Core, LINC, feedback epoch와 physical time을 모두 복구한다.
- Fluid/Core의 native A5000+ wall time, memory, launch, iteration cost가 phase ledger로 측정된다. Reduced
  Core의 시간 목표는 이 native 기준과 response/force/work mapping gate가 생긴 뒤 optimisation wave에서
  별도로 정한다.

이 상태에서 핵은 단순히 찌그러지는 visual mesh가 아니고, cytosol도 배경 pressure texture가 아니다.
둘은 SF/MT/IF 연결망과 surface body가 사용하는 양방향 physics actor가 된다. Reduced Core는 이 native
closeout 뒤에만 생산 eligibility를 얻을 수 있는 최적화 후보다.
