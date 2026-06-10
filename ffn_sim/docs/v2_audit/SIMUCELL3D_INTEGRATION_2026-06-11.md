# SimuCell3D 통합 계획 (Integration Plan) — 2026-06-11

대상: ffn_cellsim HOOMD-blue DCM 스택의 cell-cell contact 엔진 교체 및 spheroid/division 업그레이드
참조 구현: Runser, Vetter, Iber, "SimuCell3D: three-dimensional simulation of tissue mechanics with cell polarization," *Nat Comput Sci* 4:299-309 (2024), doi 10.1038/s43588-024-00620-9 (+ Publisher Correction doi 10.1038/s43588-024-00635-2). 코드: `/Users/sw1/simucell3d_ref/` (BSD-3, C++/OpenMP).

> 주의 — 본 문서는 검증 verdict가 교정한 사실을 우선 채택했고, **Lead가 C++ 원본을 직접 재유도해 한 가지 verdict 오류를 재교정**했습니다. (1) **adhesion 힘 = 논문 Eq.2의 bilinear traction-separation tent와 정확히 일치** (`F=ω·A·d` 상승 / `ω·A·(c−d)` 하강, §2.2 step 4). C++의 중간 스칼라 `force_amp`만 보면 `(c_adh/min_d−1)`로 1/d처럼 보이나, 실제 힘 `F=r_vec·force_amp`에서 **비정규화** r_vec의 크기(=min_d)가 1/d를 상쇄한다 — 본 문서 작성 워크플로의 verdict('1/d-falloff, bilinear 아님')는 이 곱셈 단계를 놓친 오류이며 **여기서 교정함**(저장 에너지 `−dE/dx`도 같은 tent를 줘 정합). (2) division 절단면 축은 covariance의 **largest** eigenvalue (논문/SI의 "smallest"는 오류). (3) edge-swap 임계값은 **0.2** (논문의 0.3은 오류). (4) `l_max = 3*l_min`은 **하드코딩**(XML 인자 아님). (5) 파라미터는 paper Table 2와 코드 example XML이 여러 곳에서 다르며, 포팅 시 **의도적으로** 골라야 합니다.

---

## 1. 요약 (TL;DR) — 권장 경로

- **물리는 포팅하고 코드는 포팅하지 않는다.** SimuCell3D를 외부 엔진으로 통째로 채택하는 것은 우리 BAOAB integrator / GPU-main 원칙 / HOOMD particle-bond 데이터 모델과 맞지 않는다. 우리가 가져올 것은 **node-vs-face penalty contact (CONTACT_MODEL_INDEX=0)** 의 힘 법칙 하나다.
- **coupling 모델(index 1,2)은 포팅하지 않는다.** 두 모델은 adhesion을 node 위치를 midpoint로 snap하고 force/momentum을 평균내는 **hard weld**로 구현한다 — force pipeline 밖에서 위치를 덮어쓰므로 HOOMD/BAOAB와 깨끗하게 매핑되지 않는다. node-face spring은 순수 가산력(additive force)이라 `md.force.Custom` Action으로 drop-in 가능하고 GPU 친화적이다.
- **HOOMD Custom force가 정답이지 SimuCell3D 엔진 직접 구동이 아니다.** node당, 다른 cell의 최근접 face까지의 signed distance를 계산 → 침투 시 repulsion(non-penetration), 외부 측이면서 cutoff 안이면 adhesion → barycentric으로 face 3노드에 분배. 이 Action을 우리 `md.Integrator`에 append하면 BAOAB가 그대로 적분한다.
- **누락된 shell 에너지 항을 추가한다.** 우리 DCM에는 explicit surface tension(γ), Helfrich bending(k_b), global area-elasticity(k_a, A0=cbrt(Q0·V0²))가 **없다**. ball→compact→spread 형태는 edge spring이 아니라 **γ vs ω 경쟁**이 지배하므로 surface tension과 (선택적) bending을 반드시 추가해야 한다. edge spring/turgor/substrate well/BAOAB는 유지한다.
- **remeshing이 숨은 필수 의존성이다.** contact 힘이 face area에 비례하고 division이 edge-length bound를 요구하므로, split/collapse/swap(`l∈[l_min, 3·l_min]`, `S_f<0.2`)이 없으면 큰 spreading에서 sliver triangle이 생겨 힘이 폭주한다. 이것이 우리 H.7 "LJ explosion"의 정체다. CPU host-side `CustomUpdater` + State rebuild로 구현(기존 `ProliferationUpdater`의 epoch 머신 재사용).
- **division은 dcm_prolif의 rim-division 트리거(convex-hull rim + necrosis gating)는 유지하고 daughter 배치 메커니즘만 mesh-split으로 교체한다.** centroid 통과 + PCA longest-axis(largest eigenvalue) 절단면 + Poisson+Delaunay septum + `target_volume/2`. 회귀 금지 기준: 7→18, A/A0 2.38.
- **하이브리드 전략을 채택한다 (reuse-the-binary로 검증, port-the-physics로 production).** 이미 빌드·검증된 `/Users/sw1/simucell3d_ref/build/simucell3d`를 oracle로 써서 우리 HOOMD contact force의 Young-Dupre / interpenetration / A/A0 결과를 대조한 뒤, production은 HOOMD GPU에 포팅한 force로 돌린다.
- **시작 파라미터는 mcf7_p0.xml을 베이스로** (단일 isotropic MCF7 cell type, γ=2.7e-4 N/m, adherence=repulsion=1e8, K=2.5e3 Pa), 단 절대값이 아니라 **turgor:tension:adhesion:repulsion 비율**을 맞추고 우리 단위계로 변환·재교정한다. adherence는 cadherin 모듈과 중복 가능하므로(아래 4e) 실험별로 결정한다.

---

## 2. SimuCell3D 모델 정수

### 2.1 Cell representation 및 내부 에너지

각 cell = **disjoint closed triangulated surface mesh** (cell 간 node 공유 없음, vertex model과 달리 polyhedral/convex 제약 없음). 우리 icosphere shell과 node/face/edge 구조가 1:1 대응한다.

총 에너지 (Eq. 1):
```
U = K·V·(ln(V/V0) − 1)              # (1) 압축성 cytoplasm 압력 (turgor)
  + (k_a/2)·(A/A0 − 1)^2            # (2) global membrane area elasticity
  + ∫_∂Ω [ γ + (k_b/2)·(2H)^2 ] dS  # (3) cortical surface tension + Helfrich bending
```

**(1) Pressure / turgor — log-bulk 법칙** (우리의 선형/harmonic turgor와 다름):
```
p = −K·ln(V/V0)   (cap p ← min{p, p_max})
f_{p,f} = p·A_f·n̂_f   →  node당  f_{p,i} = (1/3)·Σ_{f∈F_i} f_{p,f}
V = (1/6)·|Σ_f det[r_i r_j r_k]|   # divergence-theorem (우리 DcmTurgorForce와 동일)
```
우리 turgor의 enclosed-volume 계산은 이미 정확히 일치한다(유지). 단 압력 법칙은 선형(`dP0+K_vol·(V0−V)/V0`) → log-bulk 전환 검토(성장과 contact 압력의 열역학적 일관성).

**(2) Area elasticity — GLOBAL** (per-triangle 아님, verdict 교정):
```
A0 = cbrt(Q0·V0^2)   # Q0 = target isoperimetric ratio (sphere = 36π ≈ 113)
membrane_factor = −(k_a/A0)·(A_cell/A0 − 1)   # whole-cell A_cell 사용
f_{m,i} = membrane_factor · Σ_{f∈F_i} grad_i A_f
```
Q0=250(deformable)이면 비구형/spreading 허용, Q0=1(코드 filler) 또는 113(sphere)이면 round-up. **A0의 V0 결합**이 성장-형태 결합의 핵심 knob이다.

**(3a) Surface tension** — area-gradient 기반 nodal force (cheap, 가장 지배적인 round-up 항):
```
f_{s,i} = −Σ_{f∈F_i} γ_f · grad_i A_f ,   grad_i A_f = (1/2)·n̂_f × (r_k − r_j)
```
주의(verdict): 코드의 stored energy는 `0.5·γ·A_f`(논문 Eq.1 적분항은 1/2 누락)지만 **force는 양쪽 동일**(`−γ_f·grad A_f`). 동역학에는 영향 없음.

**(3b) Bending — Wardetzky discrete shell** (per interior edge / diamond i,j,k,l):
```
U_b ≈ Σ_(i,j) k̄_b·(‖e_ij‖^2/A_ij)·(2cos(θ_ij/2))^2 ,  k̄_b=(k_{b,a}+k_{b,b})/2, A_ij=A_a+A_b
θ_ij = −sgn(n̂_a·e_d)·arccos(−n̂_a·n̂_b)
```
힘 prefactor: verdict 교정 — 실제 코드는 **3·k̄_b** (논문 force eq의 `2·k_b`는 오류; energy는 일치). dihedral θ>135°이면 skip(안정성). **모든 read XML에서 bending_modulus=0** — 즉 bending은 기본 OFF이며 형태는 γ + area/volume + contact가 만든다. 따라서 **bending은 1단계에서 보류 가능**.

총 nodal force: `f_i = f_{s,i} + f_{m,i} + f_{p,i} + f_{b,i} + f_{c,i}` (tension + area + pressure + bending + contact).

### 2.2 Cell-cell CONTACT + ADHESION (node-LJ를 대체하는 핵심)

세 가지 모델이 있으나 **포팅 대상은 node-vs-face spring penalty (CONTACT_MODEL_INDEX=0)** 하나다. C++ 소스(`contact_node_face_via_spring.cpp`)가 ground truth이다. **재유도 결과 코드의 실제 힘은 논문 Eq.2의 bilinear traction-separation과 정확히 일치한다**(아래 step 4에서 도출 — 워크플로 verdict의 "1/d, bilinear 아님"은 오류).

**알고리즘 (per step, USPG로 후보 생성 후 node n(cell c1) vs face f(cell c2≠c1)):**

1. **Closest point on triangle** — Ericson 7-region barycentric test:
```
(d2, bary) = closest_point_on_triangle(n.pos; a,b,c)   # d2 = squared distance
cpa = a·bary.x + b·bary.y + c·bary.z
r_vec = n.pos − cpa      # face→node 벡터, 방향 + 크기(min_d=√d2)를 모두 보유
```

2. **Branch selector** — face normal 기준 부호:
```
adhesive  iff  r_vec · f.normal > 0   # node가 OUTWARD 측
repulsive iff  r_vec · f.normal < 0   # node가 막을 넘음 → 침투
# ECM(type1)·nucleus(type3) 쌍은 부호 반전
```

3. **REPULSION (non-penetration)** — `!adhesive` 이고 `d2 < c_rep^2`:
```
F = r_vec · repulsion_strength · A_face    # 단위: m·(Pa/m)·m^2 = N
```
침투 깊이와 face area에 선형 비례하는 penalty spring. **node-LJ가 못 막는 shell 관통을 이것이 막는다** — 힘이 face 평면에서 측정되므로 두 shell이 node 사이로 빠져나갈 수 없다.

4. **ADHESION (cohesion)** — `adhesive` 이고 `d2 < c_adh^2`, `hardening_distance = c_adh/2`. 코드 그대로:
```
min_d = √d2 ;  A = A_face
hardening (min_d ≥ c_adh/2):  force_amp = adherence_strength·(c_adh/min_d − 1)·A
softening (min_d <  c_adh/2):  force_amp = adherence_strength·A
F = r_vec · force_amp        # r_vec = n.pos − cpa : 비정규화(크기 = min_d)
```
**재유도(중요)**: `force_amp`만 보면 1/d처럼 보이지만 `F=r_vec·force_amp`에서 `|r_vec|=min_d`가 곱해져 1/d가 **상쇄**된다. 실제 힘 크기는 ω=`adherence_strength`로:
```
hardening (c_adh/2 ≤ d < c_adh):  |F| = ω·A·(c_adh − d)   # d=c_adh에서 0까지 선형 하강
softening (0 ≤ d < c_adh/2):       |F| = ω·A·d              # 접촉(d=0)에서 0까지 선형 상승
```
**즉 실제 adhesion 힘은 d=c_adh/2에서 peak(ω·A·c_adh/2)인 대칭 tent이며, 논문 Eq.2의 bilinear traction-separation `σ=ω·d / ω·(c−d)`와 정확히 일치한다(σ×A_face=F).** 워크플로 verdict가 `force_amp`를 힘으로 오독해 "bilinear 아님"이라 한 것은 오류. 저장 에너지 `E_adh`의 `−dE/d(min_d)`도 이 tent를 재현(아래 에너지 항 정합). 포팅 시 **bilinear tent를 그대로 구현**하면 된다 — repulsion도 동형의 선형 penalty `|F|=ξ·A·d`(ξ=`repulsion_strength`).

5. **Newton-3 분배** (두 branch 공통):
```
face 3노드:  +F·bary.x, +F·bary.y, +F·bary.z
node n:      −F
```

**검증용 에너지** (`FACE_STORE_CONTACT_ENERGY`):
```
E_rep        = 0.5·repulsion_strength·A·min_d^2
E_adh(soft)  = 0.5·adherence_strength·A·min_d^2
E_adh(hard)  = 0.5·A·adherence_strength·(0.5·c_adh^2 − (c_adh − min_d)^2)
```
→ 단일 정적 triangle에 node를 밀어 넣고 `−dE/dx`가 구현된 force와 일치하는지 unit test (visualize-at-closeout / sanity-gate 규칙).

**neighbor search — USPG (uniform space-partitioning grid, 매 step 재구축, Verlet 재사용 없음):**
```
aabb_padding = max(c_rep, c_adh)
voxel_size   = 3·min_edge_len + 2·aabb_padding
# face AABB = 3노드 [min,max]를 padding만큼 팽창 → face는 겹치는 모든 voxel에 삽입
# query node는 자기 voxel 1개만 조회, same-cell face 제외(c1.id≠f.owner.id), AABB pre-filter
```
비대칭 구조(face는 padding으로 smear, node는 voxel 1개 조회)가 in-range face의 co-residency를 보장. **우리 포팅에서는 이 매-step CPU 재구축 대신 HOOMD native neighbor list를 쓴다**(아래 4a).

**(미포팅) coupling 모델 (index 1=node-node, 2=face-face):** facing low-curvature epithelial node 쌍(normal anti-align within cos45°, curvature<H_max, within c_adh)을 양방향 coupling → force 평균 + midpoint snap. 위치를 force pipeline 밖에서 덮어쓰므로 HOOMD analogue 없음 → **포팅 제외**.

### 2.3 Integration / remeshing / division

**Time integration** (DYNAMIC_MODEL_INDEX): 기본 **overdamped forward-Euler** `x += F·dt/ζ` (mass-free, 결정론적, noise 없음). dynamic 모드는 semi-implicit Euler `p += (f − ζp/m)dt; r += p·dt/m`, `m = ρ·V/N_n`. dt=1e-7 s, ζ=2.5e-10(dynamic)/3e-9(overdamped) kg/s. **우리 BAOAB(stochastic Langevin)와 다름** — 우리는 고-friction/저-kT 극한으로 overdamped를 근사한다(아래 6).

**Remeshing** (`local_mesh_refiner`, 매 iteration): edge 길이를 `[l_min, l_max=3·l_min]`(하드코딩)에 유지.
- **SPLIT** `l²>l_max²`: midpoint node, 2→4 faces; A,B는 momentum 2/3 유지, E는 1/3·(A+B).
- **COLLAPSE** `l²<l_min²` 이고 `can_be_merged`(두 node가 정확히 2개 공통 이웃 — non-manifold 방지): midpoint I, momentum A+B 합산.
- **SWAP** quality `S_f = 36·A_f/(√3·P_f²) < 0.2`(verdict: 0.3 아님): 최장 edge flip.
- 수렴 실패 시 `mesh_integrity_exception`("reduce dt"). slivers가 주된 불안정 원인.

**Growth / Division:** `target_volume += dt·growth_rate` → `p=−K·ln(V/V_target)`로 팽창. division trigger `V ≥ division_volume`(매 5 iteration 확인). 절단면: origin=centroid, normal = covariance의 **largest** eigenvalue eigenvector(=PCA long axis, Hertwig). pipeline: edge-plane intersection 삽입 → 5-node face를 3 triangle로 분할 → interface polygon을 xy로 quaternion 회전 → Poisson(spacing l_min)+2D Delaunay septum → 3D 복원 → 면을 plane side로 daughter에 배정 + interface를 양 daughter에 반대 winding 복사 → `target_volume/2` + refine. 실패 시 mother 생존(try/catch). cost ~few µs.

---

## 3. 우리 DCM 현황과 격차

### 3.1 현재 보유 (production tree: `/Users/sw1/ffn_cellsim-platform/ffn_sim/cell/`)

- **cell = SurfaceManifold icosphere shell**: subdiv 1 = 42 nodes/80 tris (dcm/prolif/state 기본), subdiv 2 = 162/320 (ecm/spread). triangle은 outward로 re-wind되어 divergence-theorem 부피가 양수.
- **edge springs** (`md.bond.Harmonic`, type `dcm_edge`): in-plane cortical elasticity, k_edge=1e-3(dcm)/5e-4 N/m, r0=mean_edge (dcm_ecm는 contractility 0.7로 prestress). — **이것이 surface tension+bending의 cruder proxy**.
- **`DcmTurgorForce`** (`md.force.Custom`, 전 cell 1개): 정확한 triangulated 부피 `V_c=(1/6)Σv0·(v1×v2)` + outward face-normal 압력 `ΔP=turgor_dP0+K_vol·(V0−V_c)/V0`(선형), face당 3노드 등분. dP0=133 Pa, K_vol=1e3/5e3 Pa. **face-normal 형태라 flattened cell이 rim을 측면으로 밀어냄(wetting driver) — 유지 필수**.
- **substrate well** (`DcmSubstrateForce` / `ModulatedSubstrateForce`): capped-harmonic, z0에서 |dz|≤range(=R_cell)이면 `F_z=−k_well·dz`, durotaxis k_sub, per-cell integrin int_mult. **cell-cell과 분리 — 유지**.
- **integration**: 커스텀 Leimkuhler-Matthews **BAOAB-limit overdamped** Action. γ_node=3.9e-10 N·s/m, dt=1e-9/5e-10/3e-10 s. 유일한 time-stepper.
- **division**: contact-inhibited **rim division on pre-allocated dormant pool** (`ProliferationUpdater`/`GatedProliferationUpdater`). rim = active centroid의 convex-hull boundary; daughter = parent_centroid + outward·(2R+gap)에 **force-free 새 icosphere** 활성화(mesh split 아님). 검증: 7→18, A/A0 1.0→2.38.
- **상위 3개 calibrated layer**: necrosis 3-zone(depth d=R_cluster−r; <40µm PROLIF / 40-150µm QUIESCENT / >150µm NECROTIC), bulk-pressure crowding proxy(→[0,6] kPa), pressure-gated junction switch(onset 0.5 kPa; q*=3.81(2D)/5.4(3D) 검증 observable; switched 시 cad_mult×0.3, int_mult×3.0).

### 3.2 정확한 격차와 code seam

| 격차 | 현재 | SimuCell3D | seam |
|---|---|---|---|
| **interpenetration** | node-pair contact 3종(모두 center-line 반발) | node-vs-face signed-distance | 아래 3개 클래스 |
| **surface tension 부재** | edge spring proxy | explicit γ area-gradient force | 신규 force 추가 |
| **area elasticity 부재** | edge로 암묵 결합 | global k_a, A0=cbrt(Q0·V0²) | 신규 force 추가 |
| **bending 부재** | 없음 | Wardetzky(기본 OFF) | 신규(보류 가능) |
| **remeshing 부재** | 고정 icosphere | split/collapse/swap | 신규 CustomUpdater |
| **division = 새 shell** | force-free icosphere 활성화 | true mesh-split | daughter 배치만 교체 |
| **turgor 압력 법칙** | 선형 | log-bulk `−K·ln(V/V0)` | DcmTurgorForce 수정 |

**교체 대상 contact 3종 (정확한 seam):**
1. `cell/dcm.py:350-366` `build_dcm_simulation`의 per-cell-type `md.pair.LJ` 블록 (NxN type matrix = GPU anti-pattern도 제거됨). same-cell WCA(eps_rep_kT=5) + diff-cell well(W_cc=W_cc_Jm2·area_per_node, W_cc_Jm2=0.2e-3).
2. `cell/dcm_prolif.py:114-200` `class CellCellAdhesion(md.force.Custom)` (프롬프트의 "DcmCellCellAdhesion"). soft capped-harmonic core + linear adhesive well, scipy cKDTree pair search, mutable `cell_of_node`.
3. `cell/dcm_spheroid_state.py:121-185` `class ModulatedCellCellAdhesion` — 동일 법칙 + per-cell `mult=√(cad_mult_i·cad_mult_j)` (junction switch 의존).
4. (전역 교체 시 4번째 site) `cell/dcm_ecm.py:946-958` multicell builder의 per-cell-type LJ cross-cell well.

**보존 계약(절대 깨면 안 됨):**
- `DcmTurgorForce.faces / face_cell / n_cells` 인터페이스 — division updater가 in-place로 mutate(`_rebuild_turgor_faces`, np.add.at dense compact index). 신규 mesh가 remeshing으로 face 배열 크기를 바꾸면 이 경로와 fixed-tag-space 가정이 깨진다 — **최대 topology 리스크**.
- **fixed tag/type space** (HOOMD + BAOAB hard constraint): create_state 후 particle type/tag 추가 불가; BAOAB는 tag-space 축소 거부. runtime에 node 수가 바뀌는 remesh/mesh-split은 **node pool 사전할당** 또는 **epoch/State-rebuild** 패턴으로만 가능.
- 신규 contact force는 `cad_mult` seam(`mult=√(cad_mult_i·cad_mult_j)`)을 보존해야 junction switch가 작동. substrate의 int_mult는 분리되어 영향 없음.

---

## 4. Task 4 deliverables

### 4a. node-LJ → SimuCell3D face-based contact: 에너지 항, 파라미터, 그리고 HOOMD Custom force vs 엔진 직접 구동 결정

**결정: HOOMD `md.force.Custom` Action으로 우리 BAOAB integrator 안에서 구현. SimuCell3D 엔진 직접 구동 아님.**

근거:
- SimuCell3D 엔진은 별도 CPU C++/OpenMP 프로그램으로, 우리 BAOAB integrator·GPU-main 의무·HOOMD particle/bond 데이터 모델 어디에도 live로 붙지 않는다 (기존 baker 결합도 의도적 offline/one-shot). python binding `.so`도 빌드되어 있지 않다.
- node-face spring(index 0)은 순수 가산력이므로 기존 `DcmTurgorForce`/`CellCellAdhesion`과 똑같이 `md.Integrator`에 append → BAOAB가 적분. GPU-portable.
- HOOMD는 node-to-triangle contact pair force를 **native 제공하지 않는다**(pair force는 모두 node-node center-line). 따라서 contact는 반드시 Custom force여야 한다. 단 shell 내부 항(edge spring/volume/area/bending)은 HOOMD 7.0.1 native mesh potential(`md.mesh.bond.Harmonic`, `md.mesh.conservation.{Area,Volume}`, `md.mesh.bending.{Helfrich,BendingRigidity}`)로 GPU-resident 이전 가능.

**가져올 에너지 항 + 파라미터:**
- **repulsion (필수)**: `F = r_vec·repulsion_strength·A_face`, branch = `r_vec·n̂_f < 0`.
- **adhesion (선택, cadherin 모듈과 중복 주의)**: bilinear tent `|F|=ω·A·(c_adh−d)`(hardening) / `ω·A·d`(softening), branch = `r_vec·n̂_f > 0` 이고 `min_d<c_adh`. (코드 표현은 `F=r_vec·force_amp`, §2.2 step 4 재유도 참조.)
- (보존) surface tension γ, area k_a, turgor — 위 §2.1.

**권장 HOOMD custom-force 설계:**

```
class FaceContactForce(md.force.Custom):
    # 사전계산: cell_of_node(mutable), face→cell, face node-id triplets,
    #          face_proxy_tag(각 face centroid에 proxy particle), exclusion-by-cell
    def set_forces(self, ts):
        pos = local_snapshot.particles.position[argsort(tag)]   # tag-ordered (ParticleSorter-safe)
        # 1) face normal/area 매 step 재계산 (force ∝ A_face 이므로 deforming mesh 추적)
        #    n_f = (v1-v0)×(v2-v0); A_f = 0.5‖n_f‖; n̂_f = n_f/‖n_f‖
        # 2) neighbor: HOOMD nlist(node ↔ face-centroid proxy), r_cut = max(c_rep,c_adh)+max_edge_len
        #    (+edge term: node가 centroid가 ~1 edge 떨어진 face에 최근접일 수 있음)
        #    same-cell pair는 HOOMD exclusion(cell-id tag)로 skip
        # 3) 각 (node, face) 후보: closest_point_on_triangle (Ericson 7-region) → r_vec, min_d, bary
        # 4) branch + force_amp (§2.2 step 3,4)
        # 5) scatter: +F·bary → face 3노드, −F → node;  adhesion에 mult=√(cad_mult_i·cad_mult_j)
```

핵심 설계 결정:
- **neighbor search**: SimuCell3D의 매-step CPU USPG를 재구축하지 말고 **HOOMD native nlist**를 쓴다. 두-pass: (a) 각 face centroid에 proxy particle 배치 + (n1,n2,n3) lookup, (b) node↔proxy nlist, `r_cut = max(c_rep,c_adh)+max_edge_len`. same-cell은 cell-id 기반 exclusion list.
- **face normal/area는 매 step force 내부에서 재계산** (frozen 금지 — deforming mesh).
- **closest-point kernel**(Ericson, ~30 LOC)이 핵심 신규 geometry 코드. cheaper first cut으로 face-plane projection을 쓸 수 있으나 정확한 non-penetration엔 true closest-point 필요. interior branch의 cpa 전사 quirk는 무시(반환 barycentric은 정확).
- GPU: 처음엔 Python `md.force.Custom`(CPU local arrays)으로 정확성 검증 → 이후 closest-point kernel을 cupy/CUDA로(gpu_local) 포팅. 처음부터 gpu_local로 쓰는 것이 per-step CPU sync 병목 재유입을 막음.

**XML→우리 파라미터 시작값 매핑** (verdict: paper Table 2와 코드 example이 다르므로 **의도적 선택**; mcf7_p0.xml 베이스 권장):

| 파라미터 | 코드 example | paper Table 2 | 권장 시작 | 비고 |
|---|---|---|---|---|
| `repulsion_strength` | 1e9 | 1e9 (ξ) | **mcf7: 1e8** Pa/m | 양 소스 1e9 일치하나 mcf7은 1e8 |
| `adherence_strength` | 0 (OFF) | 1e9 (ω) | **0 또는 cadherin 대체** | cadherin과 중복(4e) |
| `c_rep = c_adh` | 7.5e-7 | 2e-7 (c) | **5e-7** m (mcf7) | 두 cutoff 독립 |
| `hardening_distance` | c_adh/2 | — | 파생 | |
| `min_edge_len` | 1e-6 | 2e-7 | **7.5e-7** m (mcf7) | |
| K (bulk) | 1e4 | 2500 | **2500** Pa | turgor log-bulk |
| Q0 | 150 | 250 | **250** (deformable) | A0=cbrt(Q0·V0²) |
| p_max | uncapped | 2500(dyn)/INF(over,mcf7) | **INF** (mcf7) | verdict: "both 2500" 오류 |
| γ | apical 1e-3 등 | 1e-3 | **2.7e-4** N/m (mcf7) | |
| dt | 1e-7 | 1e-7 | **재교정 필요** (아래 6) | |

> **단위 경고**: ξ/ω는 Pa/m = N·m⁻³(동일). raw 숫자를 우리 단위계로 복사 금지 — 변환 후, 비차원수 `γ̄=γ/(K·l)`, `ω̄=ω·l/K`(l=⟨V0⟩^(1/3))로 spreading regime을 튜닝.

### 4b. Spheroid construction (ball→compact→spread, multicellular) — 구체 레시피

SimuCell3D 자체에는 **N-cell close-pack initializer가 없다**(big_sphere 1개에서 division으로 채움). spheroid builder는 우리가 작성해야 한다. 두 단계 레시피:

**1단계 (즉시, 기존 자산 재사용):**
1. `_cluster_centers(mode='3d')`(dcm.py:222-262)로 FCC ball 중심 N개 생성. spacing은 **2.05-2.3·R**(t=0 node overlap 금지 — 강한 adhesion+t=0 overlap은 BAOAB 폭주, dcm.py:277-279).
2. 각 중심에 icosphere shell 배치(subdiv 2 권장 — basal node 충분), outward-wound.
3. 신규 `FaceContactForce`(repulsion only, adherence 약하게) + turgor + γ로 짧은 overdamped relaxation → 이웃이 Voronoi-like polyhedral contact로 정착(repulsion이 non-penetration, adhesion/γ가 gap 닫음).
4. substrate에 안착시켜 spreading phase로(`lower_to_substrate`, half-R gap, dcm_spheroid_state.py:696-717).

**2단계 (full fidelity, SimuCell3D 템플릿):**
- 432-cell **Voronoi-tessellation-of-sphere** vesicle 레시피(non-growing lumen + 선택적 ECM surface). lattice icosphere packing(`spacing_factor 2.1·R`)을 Voronoi seeding으로 교체 → space-filling cells. Poisson-disk + Ball-Pivoting watertight reconstruction은 후속.
- 문헌 대안: Thomson-dual sphere(Rozman 2020) 또는 Fibonacci-Voronoi-dilation(Laussu 2025) — 둘 다 near-uniform cell area를 줌.

**ball→compact→spread 물리**: paper의 핵심 — `ω̄=ω·l/K`(adhesion) ↑ → mutual contact area 최대화(compact/spread), `γ̄=γ/(K·l)`(cortical tension) ↑ → round-up(ball). **edge spring만으로는 안 됨** — γ vs ω 경쟁이 지배. 따라서 1단계에서도 γ는 반드시 켠다.

### 4c. Division을 dcm_prolif rim-division에 매핑 — 무엇이 바뀌나

**유지** (상위 gating):
- contact-inhibition 트리거 = convex-hull rim 검출(`ProliferationUpdater.act`, dcm_prolif.py:297-384).
- necrosis state gating(PROLIFERATING rim cell만, `GatedProliferationUpdater`).
- epoch/State-rebuild 머신(fixed-N HOOMD 대응) + `_rebuild_turgor_faces`.

**교체** (daughter 배치 메커니즘만, dcm_prolif.py:350-373 블록):
- 기존 "dormant pre-built icosphere를 parent_centroid+outward에 force-free 활성화" → **true mesh-split**:
  1. trigger를 `V ≥ avg_division_volume`으로(rim은 eligibility predicate로 위에 유지). mcf7: 3.534292e-15 m³(std 10%), 삭제 min_vol=1.767e-16.
  2. 절단면: centroid 통과, normal = node covariance의 **largest eigenvalue** eigenvector(numpy.linalg.eigh, 3×3, 신규 의존성 없음). verdict: largest(논문 "smallest" 오류).
  3. edge-plane intersection walk → 5-node face를 3 triangle 분할 → interface를 xy로 quaternion 회전 → Poisson(l_min)+2D Delaunay septum → 3D 복원 → plane side로 face 배정 + interface 양 daughter 반대 winding 복사.
  4. `target_volume = mother/2` (daughter pressure jump 방지 — 중요), refine to `[l_min,l_max]`.
- **mass/volume 보존**: daughter가 mother의 실제 sub-volume이 됨(기존엔 새 full-size shell).
- **HOOMD 장애물**: mesh split은 node/face/triangle을 **추가**하므로 fixed-N State와 충돌 → epoch/rebuild 패턴 필수(기존 budding이 쓰는 것과 동일). split/collapse 시 신규/제거 node에 일관된 velocity/force를 주어 에너지 주입 방지(BAOAB-aware; SimuCell3D는 split 시 1/3 momentum 이전).

**회귀 금지 기준**: 7→18 cells, A/A0 ∈ [2,4] (현재 2.38). 이것이 division의 dominant route(A/A0=a+b/R+c/R² 중 b/R 항). 순수 mechanical wetting은 A/A0~1.2에서 cap.

### 4d. Reuse-the-binary vs port-the-physics — 노력 추정, 장단점, 권장

**옵션 A — reuse the binary** (`/Users/sw1/simucell3d_ref/build/simucell3d`를 spheroid에 직접 실행 후 우리 분석 post-process):
- 방법: spheroid VTK(N cells, 각 type-42 polyhedron + cell_type_id) 생성 → parameter XML → `./simucell3d params.xml` → cell_data/face_data VTK frame에서 necrosis(depth>150µm 3-zone)/pressure/junction + A/A0=a+b/R+c/R² 분석.
- **노력: ~2-4일** (spheroid VTK generator + XML baker는 이미 존재, multi-cell generator만 작성; convert_to_cell_mesh.py 재사용; 분석 스크립트는 VTK reader만 신규).
- 장점: 검증된 face-contact 물리를 즉시 사용; division/remeshing이 native; 빠른 oracle(16 cells/8001 steps/~53s on Mac).
- 단점: **CPU/OpenMP only — GPU-main 의무 위반**; offline/one-shot(live force injection 불가); A5000 production path에 못 들어감; 우리 BAOAB/calibrated layer와 분리.

**옵션 B — port only contact/adhesion terms into dcm.py**:
- 방법: §4a의 `FaceContactForce` + §2.1 누락 항 + §3 remeshing.
- **노력: ~3-5주** (closest-point kernel ~3일, nlist 통합 ~3일, surface tension/area ~2일, remeshing+State-rebuild ~1-2주(최대 난이도), division mesh-split ~1주, GPU 포팅 ~1주).
- 장점: GPU-resident, BAOAB·calibrated layer 호환, production path; HOOMD native mesh potential로 shell 항 GPU 이전.
- 단점: 가장 큰 엔지니어링; remeshing이 fixed-tag-space와 충돌; sliver/CFL 안정성 디버깅.

**권장 — 하이브리드 (binary로 검증, port로 production):**
1. **즉시**: 옵션 A로 oracle 확보 — 동일 spheroid를 binary로 돌려 Young-Dupre 접촉각, interpenetration(IoU/compute_interpenetration, screening_analysis/ 재사용), A/A0=a+b/R+c/R² fit, 7→18 division을 ground truth로 기록.
2. **production**: 옵션 B로 HOOMD GPU에 포팅, 각 단계를 oracle와 대조. unit test = §2.2 에너지로 `−dE/dx` 검증, integration test = binary vs HOOMD의 A/A0 곡선 일치.
- 이 순서는 "검증된 물리 위에 GPU 구현을 쌓는다"는 무-마법넘버 원칙과 일치하며, port 디버깅 중 oracle가 항상 곁에 있다.

### 4e. Necrosis / junction / active-traction layer와의 결합 — face-based contact가 압력·응력 입력을 어떻게 바꾸나

세 layer는 **active-cell centroid + SpheroidStateArrays만 읽으며 contact force에 agnostic**하다. 단 다음을 honor하면 그대로 호환:

- **necrosis 3-zone** (depth d=R_cluster−r): centroid 기반이므로 contact 교체와 **무관**. face contact가 더 나은 compaction을 주면 R_cluster/centroid 분포가 물리적으로 정확해져 depth zoning이 오히려 개선. **변경 없음**.
- **bulk-pressure crowding proxy** (crowd→[0,6] kPa): 현재는 2.2R 내 centroid 수의 **proxy**다. face contact는 **실제 per-cell contact pressure**(turgor p=−K·ln(V/V0), contact area fraction, contact force 합)를 제공하므로, proxy를 실제 압력으로 **업그레이드 가능**(선택). SimuCell3D simulation_statistics.csv가 contact_area_fraction을 이미 출력 — oracle로 proxy 보정. junction 임계(0.5-5 kPa)는 그대로 유지(literature onset 0.5 / saturation 5 kPa). solid-stress 검증 bound: Helmlinger/Jain 6-16 kPa 성장억제 band **아래**에 있어야 하므로 face contact가 이 범위에서 stress-arrest를 보이면 안 됨(sanity bound).
- **junction switch** (P>0.5 kPa → cad_mult×0.3, int_mult×3.0): **신규 contact force가 `mult=√(cad_mult_i·cad_mult_j)` seam을 반드시 보존**해야 cell-cell adhesion 약화가 작동. substrate의 int_mult는 분리되어 그대로. q*=3.81/5.4는 외부 계산 observable이라 영향 없음.
- **active-traction / FA layer** (dcm_ecm catch-slip clutch): cell node↔fiber bead WCA + FA clutch로 **cell-cell과 직교**. FaClutchForce/EcmBondForce/DcmTurgorForce는 그대로. 단 dcm_ecm multicell builder의 per-cell-type LJ cross-cell well(946-958)은 전역 교체 시 4번째 site.

**압력 입력의 변화 요약**: turgor를 선형→log-bulk로 바꾸면 모든 layer가 소비하는 압력의 *법칙*이 바뀌므로(절대 스케일은 재교정), junction 임계와 pressure proxy를 oracle로 재보정해야 한다. contact가 node-LJ→face-based로 바뀌면 cell 간 압력 전달이 면적 적분으로 더 부드러워져 압력 분포가 매끄러워진다(layer 호환은 유지, 보정만 필요).

---

## 5. 단계별 로드맵 (phased)

> validation gate는 (G1) Young-Dupre 접촉각, (G2) interpenetration=0, (G3) A/A0=a+b/R+c/R² fit(r²≥0.95), (G4) 7→18 division A/A0∈[2,4](현재 2.38) — 회귀 금지.

**Phase 0 — Oracle 확보 (옵션 A, ~2-4일)**
- multi-cell spheroid VTK generator 작성(N type-42 polyhedra + cell_type_id) 또는 convert_to_cell_mesh.py 사용.
- binary로 7-cell→division, 432-cell vesicle 실행; A/A0 / interpenetration / Young-Dupre를 ground truth로 baker+morphology_vis로 기록.
- **Gate**: binary가 G1,G2,G4 재현 → oracle dataset 확정.

**Phase 1 — Contact force 코어 (옵션 B 시작, ~1.5주)**
- `closest_point_on_triangle`(Ericson) + `FaceContactForce`(repulsion only) Python `md.force.Custom`.
- node↔face-centroid-proxy nlist + cell-id exclusion.
- unit test: §2.2 에너지로 `−dE/dx` 일치. 단일 정적 triangle에 node 밀어넣기.
- **Gate G2**: 두 shell이 관통하지 않음(node-LJ 대비).

**Phase 2 — Shell 에너지 보강 (~1주)**
- surface tension γ(area-gradient) + global area k_a(A0=cbrt(Q0·V0²)) 추가. turgor를 log-bulk로 전환(선택).
- HOOMD native mesh potential(`md.mesh.bond.Harmonic`, conservation Area/Volume)로 가능한 항 GPU 이전.
- adhesion branch 추가, `mult=√(cad_mult_i·cad_mult_j)` seam 보존.
- **Gate G1, G3**: 2-3 cell doublet/triplet Young-Dupre; spheroid A/A0 fit.

**Phase 3 — Remeshing (~1.5-2주, 최대 난이도)**
- split/collapse/swap(`l∈[l_min,3l_min]`, `S_f<0.2`) CPU CustomUpdater + State rebuild(epoch 패턴, node pool 사전할당).
- BAOAB-aware momentum 이전(split 1/3, collapse 합산).
- turgor face-group 계약(`faces/face_cell/n_cells`) 재구축.
- **Gate**: 큰 spreading에서 sliver/blow-up 없음(H.7 LJ explosion 회귀 검사).

**Phase 4 — Mesh-split division (~1주)**
- §4c: V≥division_volume 트리거(rim+necrosis gating 유지), PCA-largest-axis cut, Poisson+Delaunay septum, target_volume/2.
- **Gate G4**: 7→18, A/A0 2.38 재현(±회귀 없음).

**Phase 5 — GPU 포팅 + spheroid builder + layer 재보정 (~1.5주)**
- closest-point + nlist를 cupy/CUDA(gpu_local), per-step CPU sync 제거.
- §4b Voronoi spheroid builder.
- pressure proxy를 실제 contact pressure로 보정(oracle 대조), junction 임계 재확인.
- visualize-at-closeout: outputs/h{X}/figs/ PNG + REPORT.md §Figures.

총 예상: oracle 2-4일 + production ~6-7주.

---

## 6. 리스크 & 미해결

- **timestep vs contact stiffness vs BAOAB**: SimuCell3D는 ξ=1e9 Pa/m을 dt=1e-7 s + ζ damping + remeshing으로 안정화한다. 우리 BAOAB overdamped CFL은 `dt < γ_node/(k_rep·A_face_max)`. **1e9를 그대로 복사 금지** — 훨씬 부드러운 k_rep로 시작해 interpenetration이 멈출 때까지 ramp. 우리 H.7 "LJ explosion"(cap/stiff-tether/box>3)이 정확히 stiff contact + 과신장의 sliver 실패이며 remeshing이 이를 막는다.
- **remeshing vs fixed-tag-space**: runtime node 수 변화가 HOOMD State와 BAOAB tag-buffer와 충돌. **가장 큰 미해결 엔지니어링**. node pool 사전할당 + epoch rebuild가 유일 경로. turgor face-group 재구축이 매번 일관되어야 함.
- **adhesion 파라미터 매핑 (verify in calibration)**: 힘 *형태*는 확정 — 논문 Eq.2 bilinear tent = 코드(재유도 일치, §2.2 step 4; 워크플로의 'bilinear 아님' verdict는 교정됨). 남은 불확실성은 *크기*: 우리 `W_cc_Jm2`(J/m²) ↔ SimuCell3D `ω`(Pa/m=N·m⁻³) 매핑을 contact triangle 적분 면적/길이 factor를 거쳐 검증해야 함(baker의 P0.1 calibration TODO와 동일). cadherin 모듈과 double-count 금지 — repulsion은 항상, adhesion은 선택.
- **파라미터 출처 충돌**: paper Table 2 vs 코드 example XML이 K(2500/1e4), Q0(250/150), c(2e-7/7.5e-7), p_max(2500/INF), ω(1e9/0)에서 다름. 절대값이 아니라 **비율**을 맞추고 우리 단위계로 재교정. mcf7_p0.xml 베이스 권장.
- **GPU implications**: SimuCell3D는 CPU/OpenMP only. 단일 ~7k-scale spheroid는 A5000을 미포화(util 0-4%, CPU 대비 ~2.8×). 432-cell·subdiv2(~50k-300k nodes)는 GPU가 paying off 시작하는 regime. 풀활용은 sweep batch 또는 대형계. closest-point/nlist를 처음부터 gpu_local로 작성해 per-step CPU sync 병목 재유입 방지.
- **closest-point 정확도**: face-centroid projection은 cheap first cut이나 true edge-edge/point-triangle closest-point가 정확한 non-penetration에 필요. Ericson interior branch의 cpa 전사 quirk는 무시(반환 barycentric 정확).
- **division 실패 처리**: degenerate cut 시 mother 생존(try/catch). subdiv 1(42 node)은 division septum의 Poisson+Delaunay에 너무 coarse — division 전 subdiv 2 또는 adaptive remesh 필요.
- **bending 보류 타당성**: 모든 read XML에서 bending=0이라 1단계 생략 가능하나, 매우 비구형 spreading에서 막 강성이 필요하면 Wardetzky(force prefactor 3, θ>135° skip) 추가. HOOMD native Helfrich로 GPU 무료 추가 가능.
- **necrosis 검증 oracle**: depth>150µm는 확산 제한 O2 침투 깊이(~100-200µm)와 일치. 연속 O2/glucose field 결합(Schaller-Meyer-Hermann)을 necrotic-core oracle로 추가 검토.

---

## 7. 참고문헌

핵심:
- Runser S, Vetter R, Iber D. *SimuCell3D: three-dimensional simulation of tissue mechanics with cell polarization.* Nat Comput Sci 4(4):299-309, 2024. https://doi.org/10.1038/s43588-024-00620-9 (full text PMC11052725; preprint bioRxiv 2023.03.28.534574)
- Runser S, Vetter R, Iber D. *Publisher Correction.* Nat Comput Sci 4(5):379, 2024. https://doi.org/10.1038/s43588-024-00635-2
- 코드: https://git.bsse.ethz.ch/iber/Publications/2024_runser_simucell3d (BSD-3); Zenodo https://zenodo.org/records/10796908. 로컬: `/Users/sw1/simucell3d_ref/` (binary `build/simucell3d`).

DCM / contact 계열:
- Vetter R, Runser S, Iber D. *PolyHoop: Soft particle and tissue dynamics with topological transitions.* Comput Phys Commun 299:109128, 2024. https://arxiv.org/abs/2307.15006 (2D sibling, contact penalty math 최경량 참조)
- Odenthal T, Smeets B, Van Liedekerke P, et al. *Analysis of Initial Cell Spreading Using Mechanistic Contact Formulations for a Deformable Cell Model.* PLoS Comput Biol 9(10):e1003267, 2013. https://doi.org/10.1371/journal.pcbi.1003267 (Maugis-Dugdale 고정밀 대안)
- Van Liedekerke P, Neitsch J, Johann T, et al. *A quantitative high-resolution computational mechanics cell model for growing and regenerating tissues.* Biomech Model Mechanobiol 18:1297-1321, 2019. https://doi.org/10.1007/s10237-019-01204-7 (triangulated-DCM spheroid 성장/division 템플릿)
- Zhao J, Manuchehrfar F, Liang J. *Cell-substrate mechanics guide collective cell migration through intercellular adhesion (DyCelFEM).* Biomech Model Mechanobiol 19:1781-1796, 2020. https://doi.org/10.1007/s10237-020-01308-5 (face-to-face cadherin spring 대안)

GPU 참조 구현:
- Madhikar P, Astrom J, Westerholm J, Karttunen M. *CellSim3D: GPU accelerated software for simulations of cellular growth and division in 3D.* Comput Phys Commun 232:206-213, 2018. https://doi.org/10.1016/j.cpc.2018.05.024 (github.com/SoftSimu/CellSim3D — GPU triangulated divided-cell 최근접 참조)
- HOOMD-blue 7.0.1 mesh module (md.mesh bending/conservation/bond; Mesh triangle types). https://hoomd-blue.readthedocs.io/en/stable/

Spheroid 구성 레시피:
- Rozman J, Krajnc M, Ziherl P. *Collective cell mechanics of epithelial shells with organoid-like morphologies.* Nat Commun 11:3805, 2020. https://doi.org/10.1038/s41467-020-17535-4 (Thomson-dual sphere seeding)
- Laussu J, et al. *Deciphering the interplay between biology and physics with a FEM-implemented vertex organoid model.* PLoS Comput Biol 21(1):e1012681, 2025. https://doi.org/10.1371/journal.pcbi.1012681 (Fibonacci-Voronoi-dilation builder)

Spheroid pressure / necrosis 검증 oracle:
- Helmlinger G, Netti PA, Lichtenbeld HC, Melder RJ, Jain RK. *Solid stress inhibits the growth of multicellular tumor spheroids.* Nat Biotechnol 15(8):778-783, 1997. https://doi.org/10.1038/nbt0897-778 (6-16 kPa 성장억제 band — junction 0.5-5 kPa는 아래)
- Schaller G, Meyer-Hermann M. *Continuum versus discrete model: a comparison for multicellular tumour spheroids.* Philos Trans A 364(1843):1443-1464, 2006. https://doi.org/10.1098/rsta.2006.1780 (necrotic-core O2/glucose oracle)
- Yan H, Ramirez-Guerrero D, Lowengrub J, Wu M. *Stress generation, relaxation and size control in confined tumor growth.* PLoS Comput Biol 17(12):e1009701, 2021. https://doi.org/10.1371/journal.pcbi.1009701

> 출처 표기 주의: 위 PubMed/문헌 재사용 시 DOI 링크 및 PubMed 귀속 표기 유지.
