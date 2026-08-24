# ECM World 확정 계획

> **FF reuse decision (2026-07-22):** `ff/ecm_library.py` is the canonical material/topology source and
> `ff/ecm_mikado.py` is its older collagen-only precursor. The new engine must port/reuse the former's
> collagen-I card, concentration/alignment geometry, segment WLC/EA, bending and crosslink topology rather
> than inventing a second ECM. Only its host KD-tree/build mutation and validation relaxation loop are
> quarantined from the Warp-CUDA physical-time runtime.

상태: **ECM-W1 구조 슬라이스 구현 및 구조 게이트 완료**
소유 컴포넌트: `ecm`
환경: MCF7 + collagen-I ECM + α2β1–collagen clutch
생산 런타임: Warp-CUDA 전용

## 1. 확정 결론

ECM은 고정된 바닥이나 탄성률 하나로 환산한 연속체가 아니다. 생산 모델의 ECM은 다음을 각각
명시적으로 보존하는 **live collagen fiber/crosslink world**로 확정한다.

| 항목 | 확정안 |
|---|---|
| 기본 표현 | collagen fiber별 node/segment rod network + explicit crosslink population |
| 내부 역학 | segment axial law + bending + excluded volume + crosslink mechanics |
| 세포 연결 | geometry-less FA의 두 semantic edge를 하나의 SF–FA–α2β1–collagen composite joint로 해석 |
| ECM endpoint | 고정 node가 아니라 collagen segment 위 barycentric material point |
| 유한 계산영역 폐쇄 | 첫 baseline은 far-field Dirichlet anchor; 자유 ECM은 production 금지 |
| 장기 경계 대안 | 반무한 bulk를 근사하는 continuum-impedance boundary, 검증 후 교체 가능 |
| 동적 변화 | crosslink kinetics, fiber remodelling, damage, refinement/coarsening, sleep/wake 모두 명시 상태 |
| 상태 갱신 | 모든 불가역 변화는 global accepted-step transaction에서만 commit |
| 최적화 | 생물학적 fiber 수를 줄이지 않고 connected-island sleeping + local refinement로 계산량 절감 |
| 보존 검증 | clutch의 equal/opposite force·moment·work와 far-field reaction/work를 device ledger로 검증 |

이 계획은 collagen 밀도나 fiber identity를 낮추는 coarse-graining이 아니다. 생물학적 population과
topology는 유지하고, 현재 힘 전달에 관여하지 않는 영역의 **계산 빈도와 해상도만 보수적으로
줄이는 것**이다.

## 2. 컴포넌트 경계와 상태 소유권

### 2.1 `ecm`이 소유하는 상태

`ecm` actor는 다른 컴포넌트와 분리된 CUDA 배열을 소유한다.

- `position_d`, `force_d`: live collagen node kinematics와 force accumulator
- `segments_d`: collagen segment의 node pair
- `fiber_offsets_d`: unique collagen fiber identity와 segment 구간
- fiber constitutive state: contour/rest length, stretch/bending material parameters와 필요 history
- `fiber_sleep_state_d`: fiber별 wake/sleep/pending 상태
- `segment_refinement_level_d`: segment별 현재 공간 해상도
- `segment_damage_d`: segment/fibril damage 및 회복 가능한 remodelling history
- `endpoint_refcount_d`: graph connector가 붙은 segment를 sleep/remesh에서 잠그는 보수적 cache
- `topology_epoch_d`: accepted topology의 device epoch
- neighbor/BVH와 refinement scratch: candidate-step cache이며 rejected step에서 복구

`endpoint_refcount_d`는 clutch bond의 화학 상태를 복제하지 않는다. bond identity, age, kinetic state는
connector가 소유하고, ECM은 “이 segment를 현재 이동/재분할해도 되는가”라는 topology lock만
소유한다.

### 2.2 graph connector가 소유하는 상태

- `ecm_crosslink`: collagen segment–segment association, age, kinetic state, rest geometry, generation
- `alpha2beta1_collagen_series`: SF actin endpoint부터 FA state를 거쳐 collagen endpoint까지의 하나의
  composite mechanical joint
- connector RNG epoch와 candidate/committed state
- endpoint actor/entity/topology generation과 material coordinate

따라서 ECM 내부 힘 배열에 crosslink 힘은 더해지지만, crosslink bond state를 ECM component 배열에
몰래 넣지 않는다. 같은 이유로 ECM은 FA ligand state나 integrin kinetic state를 소유하지 않는다.

### 2.3 boundary가 소유하는 상태

유한 ECM domain의 far-field boundary는 ECM constitutive state와 구분된 boundary participant다.

- pinned/far-field node mask 또는 impedance quadrature
- prescribed target 혹은 impedance history
- candidate snapshot/cache
- reaction force, moment, boundary work ledger

boundary는 ECM 좌표를 빌려 쓰되 collagen topology나 position을 소유하지 않는다.

## 3. 생리적 baseline과 경계 고정

첫 MCF7 collagen baseline은 `far_field_dirichlet`로 확정한다. 계산영역 바깥의 collagen gel이 계속
존재한다는 것을 나타내도록 세포에서 충분히 먼 외곽 band를 bulk ECM 기준 위치에 고정한다.
이는 collagen node를 FA의 고정 substrate point로 바꾸는 것이 아니다. 대부분의 collagen network는
live DOF로 변형되고, 오직 잘린 representative volume의 외곽에서만 bulk 반력이 발생한다.

필수 조건은 다음과 같다.

- production config에서 `pin_face="none"` 또는 무구속 finite network를 금지한다.
- boundary mask, target, reaction은 CUDA에 유지한다.
- 고정 node의 힘을 단순히 0으로 지워 버리기 전에 **지워질 반력**을 ledger에 먼저 누적한다.
- prescribed boundary가 움직이는 실험에서는 `reaction · boundary_velocity`를 boundary work로 기록한다.
- anchor band 두께, box 크기, cell–boundary 거리는 결과를 맞추는 tuning knob가 아니라 domain-error
  convergence gate로 정한다.
- cell traction의 순힘은 ECM 내부에서 사라지지 않고 far-field reaction과 닫혀야 한다.

장기적으로는 far-field Dirichlet과 같은 domain 크기에서 검증된 continuum-impedance boundary를
지원한다. 이는 외부 bulk만 축약하는 것으로, 세포 근방의 explicit collagen network를 대체하지
않는다. impedance constitutive card와 흡수/반사 오차 gate가 준비되기 전에는 baseline이 아니다.

## 4. collagen 역학

한 candidate mechanics pass의 순서는 다음으로 고정한다.

1. caller가 ECM force accumulator를 0으로 만든다.
2. fiber axial, bending, excluded-volume/contact force를 누적한다.
3. graph-owned collagen crosslink force를 동일 ECM 배열에 adjoint scatter한다.
4. far-field boundary residual/reaction을 누적한다.
5. 하나의 α2β1 composite clutch runtime이 cell 쪽과 ECM 쪽에 equal/opposite 힘을 한 번만 scatter한다.
6. world solver가 다른 컴포넌트와 함께 inner mechanics convergence를 반복한다.

물리 시간은 outer clock 하나뿐이다. inner Newton/relaxation iteration 횟수를 시간으로 해석하지
않는다. crosslink association, damage, remodelling, sleep/refinement 전이는 inner iteration에서
반복 commit하지 않는다.

### 4.1 첫 constitutive port 범위

- collagen axial law와 bending law는 기존 Warp force kernel을 포팅/재사용한다.
- collagen fiber와 crosslink의 서로 다른 chemistry/material parameter를 분리한다.
- collagen–collagen 및 collagen–cell excluded volume은 production baseline에서 활성화한다.
- stiffness, persistence length, crosslink off/on rate, damage threshold는 literature/KB parameter card가
  있어야 하며 이 문서에서 임의 숫자를 만들지 않는다.
- 현재 `ecm_mechanics.py`의 host relaxation loop는 constitutive runtime이 아니라 validation harness다.

## 5. α2β1–collagen endpoint와 composite clutch

### 5.1 endpoint 표현

collagen ligand endpoint는 segment `e=(i,j)`와 material coordinate `u∈[0,1]`로 표현한다.

```text
x_ligand = (1-u) x_i + u x_j
F_i     += (1-u) F_ligand
F_j     += u F_ligand
```

gather와 scatter가 정확한 transpose이므로 nodal virtual work와 endpoint work가 일치한다. endpoint는
actor generation, entity generation, segment ID, `u`, topology epoch를 함께 가진다. accepted remesh가
segment를 나눌 때 기존 material point의 물리 위치와 fiber arclength를 보존해 새 segment와 `u`로
원자적으로 remap한다.

고정 ECM node index를 영구 endpoint로 사용하는 방식은 폐기한다. node endpoint는 다음 문제를
가진다.

- refinement 때 node identity에 부당하게 묶인다.
- 힘을 segment 형상함수에 분배하지 못한다.
- nearest-node host KD-tree attachment가 production GPU-only 조건을 어긴다.
- ECM node density가 clutch mechanics를 바꾸는 grid dependence를 만든다.

### 5.2 하나의 composite series joint

`fa_actin_anchor`와 `integrin_collagen_clutch`는 생물학적 의미상 두 edge지만 mechanical group
`alpha2beta1_collagen_series` 하나로 해석한다. 두 개의 독립 spring을 만들어 stiffness나 energy를
이중 계산하지 않는다. geometry-less `focal_adhesion`은 별도 위치 배열을 갖지 않고 clutch state와
chemistry를 갖는다.

composite clutch ledger는 최소한 다음을 제공해야 한다.

- cell-side와 ECM-side resultant의 합
- 두 actor 기준 torque closure
- 양쪽 instantaneous power와 합
- spring/series energy와 energy increment
- active/bound/reinforced clutch count
- endpoint generation mismatch와 remap count

## 6. sleeping과 local refinement

### 6.1 sleeping 단위

개별 node를 독립적으로 재우지 않는다. crosslink topology를 기준으로 형성된 **connected mechanical
island**가 sleeping 판단 단위다. 구현 배열은 fiber별 상태를 갖되, 같은 island의 fiber는 동일한
wake/sleep 결정을 받는다.

sleep 가능 조건은 모두 만족해야 한다.

- island의 모든 clutch/contact endpoint refcount가 0이다.
- far-field active boundary나 현재 awake island와 전달 가능한 crosslink path가 없다.
- crosslink/damage/remodelling hazard의 다음 event horizon 전이다.
- elastic residual과 저장 에너지 변화가 등록된 discretization/convergence tolerance 이하다.
- 주변 broad-phase envelope가 다음 outer step 동안 침범될 수 없다.

sleep은 topology 삭제가 아니다. unique fiber ID, nodes, crosslinks, damage, material state는 그대로
남는다. force pass와 neighbor update cadence만 낮춘다.

wake trigger는 다음과 같다.

- α2β1 clutch capture 후보나 직접 contact broad phase 진입
- 연결된 awake island에서 crosslink force 전달
- boundary perturbation
- 예정된 kinetic/remodelling/damage event
- error estimator가 refinement를 요구하는 경우

### 6.2 refinement/coarsening

refinement는 cell/FA 근방에서 필요한 segment만 세분화한다. trigger는 curvature, axial strain/energy
gradient, contact interpolation error, live connector occupancy를 이용한다. threshold는 gate를 맞추기
위한 magic number가 아니라 coarse/fine solution 오차와 virtual-work convergence에서 유도한다.

필수 보존량은 다음과 같다.

- fiber contour/arclength와 unique biological fiber identity
- mass/material measure
- linear momentum 및 가능하면 angular momentum
- elastic energy의 허용 오차
- damage/remodelling history의 material-coordinate transfer
- crosslink와 clutch endpoint의 material-coordinate remap
- endpoint refcount와 topology generation

live endpoint가 있는 segment는 remap transaction 없이 분할/병합하지 않는다. coarsening은 연결성과
crosslink geometry를 보존할 수 있고 local fine/coarse gate가 통과할 때만 허용한다.

## 7. crosslink, remodelling, damage accepted-step transaction

한 outer physical step은 다음 transaction으로 고정한다.

1. ECM, crosslink connector, boundary cache, composite clutch가 committed state를 snapshot한다.
2. GPU KMC/hazard kernel이 crosslink association/dissociation 후보를 제안한다.
3. damage/remodelling과 sleep/refinement 후보를 scratch state에 제안한다.
4. endpoint remap과 refcount 변화도 scratch에만 기록한다.
5. 전체 cell–ECM candidate mechanics/field solve를 수렴시킨다.
6. conservation, work, topology, numerical gate가 최종 device acceptance predicate를 만든다.
7. reject면 네 participant를 모두 rollback한다.
8. accept면 동일 predicate 아래 crosslink, damage/remodelling, topology epoch, endpoint remap, clutch state,
   RNG epoch를 한 번만 commit한다.

부분 commit은 금지한다. 예를 들어 collagen segment는 분할됐는데 clutch endpoint는 옛 generation을
가리키거나, crosslink는 끊겼는데 damage energy가 복구되는 상태가 생기면 안 된다.

crosslink는 collagen chemistry를 갖는 explicit connector다. actin crosslink KMC의 상태 machine과
accepted-step pattern은 참고할 수 있지만, actin-specific rate/force law를 collagen에 그대로 쓰지
않는다.

## 8. ledger와 acceptance gate

ECM World가 global ledger에 제공할 항목은 다음과 같다.

### 8.1 역학/일

- fiber, crosslink, clutch, contact별 force resultant와 torque
- clutch `F_cell + F_ecm`, `τ_cell + τ_ecm`
- clutch 양쪽 power와 work closure
- far-field reaction resultant/torque 및 boundary work
- fiber axial/bending, crosslink, clutch elastic energy
- damage dissipation, remodelling chemical/mechanical work
- candidate residual과 inner convergence

### 8.2 topology/population

- unique active/allocated collagen fiber IDs
- nodes, segments, explicit crosslinks의 active/allocated capacity
- awake/sleeping island와 fiber 수
- refinement level별 segment 수, split/merge/remap 수
- clutch-locked segment 수와 stale generation 수
- damage state population과 irreversible event 수
- exact peak GPU bytes 및 scratch/cache bytes

첫 구조 슬라이스의 `ECMWorldDofLedger`는 host-visible allocation metadata만 반환한다. active population,
force/work, event count는 device ledger가 authoritative하다. facade는 ledger 값을 읽거나 acceptance를
결정하지 않고 모든 participant의 contribution을 빠짐없이 dispatch한다.

## 9. 구현된 ECM-W1 구조 슬라이스

파일: `aleph/engine/ecm_world.py`

구현 범위:

- `ECMWorldSettings`: 자유 ECM, kinetics/damage/remodelling/sleep/refinement 누락 build 차단
- `ECMStateOwner`: 별도 CUDA 배열, shape/device/storage ownership gate
- `CollagenMaterialPointEndpointView`: non-owning segment endpoint와 typed `PortRef`
- `ECMBoundaryAnchorFacade`: explicit far-field boundary, reaction/work ledger 강제
- `ECMWorldStepBindings`: composite clutch, internal crosslink, boundary를 외부 registry에서 주입
- `ECMWorld`: mechanics, snapshot/rollback/commit, ledger의 단일 orchestration seam
- `ECMWorldDofLedger`: nodes/segments/fibers 및 state-slot allocation
- `LegacyLocalMikadoMechanicsAdapter`: local Warp mechanics만 허용하는 좁은 재사용 경계

아직 구현하지 않은 범위:

- 새 collagen constitutive kernel
- GPU neighbor/BVH builder와 collagen-crosslink partner search
- 실제 sleeping island labelling kernel
- refinement/coarsening/remap kernel
- collagen-specific KMC와 damage/remodelling law
- continuum-impedance boundary
- global ledger field schema 및 world solver integration

구조 슬라이스는 이 미구현 물리를 host proxy로 대신하지 않는다. delegate가 준비될 때까지 production
gate는 열리지 않는다.

## 10. 기존 코드 재사용/격리/폐기 판정

| 기존 파일/기능 | 판정 | 이유와 처리 |
|---|---|---|
| `ff/ecm_library.py` material/spec/validation data | **재사용/포팅** | literature-backed material card와 acceptance band를 production config/oracle로 이전 |
| `ff/forces_warp.py` collagen에 쓰인 bending primitive | **재사용 후보** | local CUDA force contribution으로 adapter/port; collagen parameter card 분리 필요 |
| `ff/network_warp.py` link/WLC/EV primitives | **재사용 후보** | 단위·adjoint·precision gate 후 ECM local kernels로 사용 |
| `ff/ecm_mikado.py:MikadoNetwork` topology schema | **참고/포팅** | fiber offsets, links, pinned metadata는 유용; production state는 CUDA SoA로 재구성 |
| `ff/ecm_mikado.py:build_mikado_network` | **격리** | NumPy/SciPy KD-tree neighbor query와 host topology mutation은 production 금지; source/oracle용 |
| `ff/ecm_mechanics.py:_Dev/_force_pass` | **kernel만 재사용** | low-level Warp force kernel은 유용하지만 host relaxation, `.numpy()` readback, private loop는 validation 전용 |
| `ff/fa_ecm.py:clutch_ecm_spring_kernel` | **parity oracle/제한적 진단** | equal/opposite 의도는 유지; fixed node endpoint는 segment material-point joint로 교체 |
| `ff/fa_ecm.py:attach_clutches_to_ecm` | **격리** | host cKDTree attachment는 GPU-only production contract 위반 |
| `ff/fa_clutch_warp.py` fixed substrate anchor | **production 폐기** | live collagen ECM reaction/remodelling을 없애므로 collagen baseline에 사용 금지 |
| `ac/engine/load_path.py` segment joint/FA group resolution | **직접 재사용** | barycentric adjoint scatter와 두 semantic edge→한 joint guard가 ECM endpoint 계약과 일치 |
| `ac/weave/crosslink_kmc*.py` accepted-step pattern | **패턴만 재사용** | actin chemistry를 collagen에 대입하지 않고 transaction/RNG 구조만 참고 |

`LegacyLocalMikadoMechanicsAdapter`가 허용하는 범위는 `node_off=0`인 local CUDA 배열에 fiber-local
force를 더하는 것뿐이다. legacy object가 crosslink, boundary anchor, host neighbour query, private
integrator를 소유하면 adapter가 거부한다.

## 11. 공통 graph/runtime gap 제안

이번 작업에서는 공통 파일을 수정하지 않았다. 통합 단계에서 다음을 별도 승인/구현해야 한다.

1. **Lead 통합 완료:** `ecm_crosslink`는 collagen용 `FIBER_CROSSLINK` family로 승격됐다.
2. **Lead 통합 완료:** `ecm_far_field_anchor: ENVIRONMENT_BOUNDARY`와 `world_boundary`가 등록됐다.
3. **Lead 통합 완료:** FA clutch와 별도인 `membrane_ecm_contact: CONTACT`가 등록됐다.
4. **Lead 통합 완료:** connector contract에 `generation_required`, `remap_on_accept`,
   `blocks_sleep_refine`가 추가됐고 live ECM endpoint connector에 적용됐다.
5. common `LedgerContributor`는 method만 정의하고 force/work/moment/topology field schema가 없다.
   clutch pair closure와 boundary reaction을 global acceptance gate가 동일 이름으로 읽도록 typed device
   ledger schema가 필요하다.
6. transaction protocol에 propose/remap dependency ordering이 없다. crosslink proposal → topology remap
   → clutch remap → single commit 순서를 world scheduler phase로 명시해야 한다.
7. **Transaction 경로 통합 완료:** `CellActor.unique_runtime_objects()`와 world transaction이 composite
   runtime의 snapshot/rollback/commit/ledger를 identity 기준으로 한 번만 dispatch한다. Mechanics phase의
   정확히 한 번 dispatch는 rolling roadmap R1 scheduler gate로 남아 있다.

## 12. 검증 계획

### 12.0 생산 gate 충돌 — PI 결정 전 HALT

FF 자산 감사에서 collagen-I modulus gate가 두 값으로 갈라진 것이 확인됐다.

- `ff/ecm_library.py`의 canonical-looking `ECMSpec.modulus_band_Pa`는 `30–100 Pa`다.
- `scripts/ff_ecm_validate.py`의 별도 `LIT` 표는 `5–100 Pa`를 사용하고 native atlas도 이를 import한다.
- 이 때문에 약 `12–15 Pa` 결과가 후자에서는 `IN BAND`지만 전자 계약에서는 `OUT`이다.

농도, probe/measurement protocol, modulus 종류별 reference를 Notion/KB SoT의 별도 validation contract로
정규화하고 PI가 승인하기 전에는 어느 band도 production gate로 임의 선택하지 않는다. 기존 Hookean
law는 reference modulus 일부를 맞추지만 농도 지수와 `N1` 부호가 맞지 않고, WLC 후보는 stiffening은
개선되지만 절대 modulus와 `N1` 문제가 남아 있다. 따라서 WLC를 단순히 더 mechanistic하다는 이유만으로
production default로 승격하지 않는다. 이 blocker는 device schema, scheduler, topology port 작업을 막지
않지만 production collagen constitutive closeout은 막는다.

### ECM-W1 구조 게이트 — 구현 완료

`aleph/tests/ac/engine/test_ecm_world.py`가 다음을 검증한다.

- CUDA-only array metadata와 독립 storage
- node/segment/fiber/state-slot shape 일치
- actor element population과 mechanics population 일치
- collagen segment material-point endpoint와 generation 전달
- internal → crosslink → boundary → composite clutch mechanics wiring
- foreign ECM endpoint와 boundary mode mismatch 거부
- equal/opposite/adjoint/work ledger contract 강제
- remodelling/crosslink/boundary/clutch의 동일 snapshot/rollback/commit predicate
- invalid `dt_phys`/RNG seed가 mutation 전에 실패
- 모든 participant ledger와 allocation ledger
- 비생리적 자유 boundary와 누락된 damage/remodelling 거부
- legacy adapter의 local/connector-free/boundary-free 범위

### ECM-W2 단위 게이트

- collagen axial/bending single-fiber force와 finite-difference energy gradient
- crosslink pair의 equal/opposite force, moment, virtual work
- far-field reaction과 applied load closure
- dormant island가 awake full solve와 허용 오차 내 동일 결과
- wake trigger가 connector/contact event를 놓치지 않음

### ECM-W3 transaction/refinement 게이트

- rejected crosslink/damage/remesh step의 bitwise committed-state 복구
- split/merge 후 contour, endpoint 위치, force/work, history 보존
- live clutch/crosslink endpoint generation의 atomic remap
- stale endpoint가 force kernel 실행 전에 차단됨

### ECM-W4 native production gate

- physiological collagen population/geometry와 anchored baseline
- full cell load path에서 SF–FA–ECM traction 전달
- force/moment/work/energy ledger closure
- full native population, exact peak GPU-byte ledger
- physical-time loop 안 authoritative GPU→CPU roundtrip 0
- local refinement/sleeping ON/OFF parity와 wall-time/real-time factor

## 13. 구현 순서

1. **ECM-W1 — 완료:** ownership, endpoint, facade, transaction/ledger seam과 구조 테스트
2. **ECM-W2:** 기존 Warp collagen force primitives를 local owner API에 연결하고 far-field reaction ledger 구현
3. **ECM-W3:** GPU collagen crosslink KMC, damage/remodelling state와 accepted-step rollback/commit 구현
4. **ECM-W4:** GPU broad phase, connected-island sleep/wake, refinement/remap 구현
5. **ECM-W5:** composite SF–FA–ECM joint와 global world solver/typed ledger 통합
6. **ECM-W6:** native MCF7 collagen baseline, profiler/physics/visual closeout

## 14. 완료 정의

ECM World는 다음을 모두 만족할 때만 production-ready다.

- ECM이 live collagen fiber/crosslink network이며 fixed substrate proxy가 아니다.
- far-field anchor가 physiological baseline에서 켜져 있고 reaction/work가 ledger에 잡힌다.
- α2β1 clutch가 collagen segment material point에 연결되고 SF 쪽과 equal/opposite/adjoint로 닫힌다.
- 두 FA semantic edge가 정확히 하나의 mechanical joint로 평가된다.
- crosslink, damage, remodelling, sleeping, refinement, endpoint remap이 accepted-step atomic transaction이다.
- sleeping/refinement가 biological fiber identity/population을 삭제하지 않는다.
- 모든 partner search, physics, mechanics iteration, KMC/PDE update가 Warp-CUDA resident다.
- native population에서 force/moment/work/topology/bytes ledger와 시각화 gate가 통과한다.
