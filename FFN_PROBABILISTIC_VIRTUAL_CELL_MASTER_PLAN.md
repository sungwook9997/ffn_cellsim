# FFN Probabilistic Mechanistic Virtual Cell — 통합 마스터 플랜

**상태:** **RATIFIED 2026-08-09** — PI가 개명(결정 A1)으로 IDENTITY 주장을 승인했다.
이 문서의 정체성 절은 이제 `CLAUDE.md` §What this is에 있고, 그쪽이 권위다.
아래 본문은 그 결론에 이른 설계 논증으로 보존한다.

<!-- 이전 상태: "PI 검토용 재창설 제안" (2026-07-29~2026-08-09) -->

**작성 기준일:** 2026-07-29

**이번 개정:** `CellState`를 고정된 세포종/상태/환경 축 목록이 아니라 확장 가능한
component–connector owner graph의 내부 상태로 수학화하고, state space·hybrid dynamics·
accepted transition·posterior·observation·intervention·tensor-network factorization을 하나의
계약으로 닫는다. Environment는 별도 conditioning context로 분리한다.

**범위:** 물리 런타임 · 확률동역학 · 텐서 네트워크 · 신경망 · 영상 AI · 지식베이스 ·
다세포종 아틀라스 · in-silico 실험 · 다세포/조직 확장

**비범위:** 이 문서만으로 기존 gate, parameter, runtime contract를 변경하지 않는다. 실제 wet-lab
프로토콜 설계와 실험 발주는 후행 단계로 둔다.

---

## 0. 한 화면 결정

### 프로젝트의 새 정체성

`ffn_cellsim`은 특정 세포주 몇 개의 파라미터를 맞추는 시뮬레이터가 아니다.

> **문헌·omics·현미경·역학 데이터를 물리적 세포 상태의 사후분포로 바꾸고, 그 분포를
> accepted-step mechanochemical dynamics로 전파하며, 다시 힘·형태·동역학·합성영상·개입반응으로
> 내보내는 Probabilistic Mechanistic Virtual Cell 플랫폼이다.**

최종 데이터 단위는 세포주 이름이나 단일 파라미터가 아니라 다음 객체다.

```text
CellStatePosterior(t)
  = StateSchemaManifest M
  + P(S_cell(t), θ, regime, model hypothesis
      | observations, external context, interventions, evidence, M)
  + observation/intervention operators
  + representation and approximation-error contract
  + evidence provenance
  + known model inadequacy
```

`external context`는 조건부 우변에만 있고 `S_cell` 안에는 없다. `S_cell`은 manifest가 선언한
component–connector owner graph에 분산 소유된 **세포 내부 상태**이며, posterior는 그 상태 자체가
아니라 그 상태에 대한 지식의 표현이다.

### 다섯 가지 표현의 역할

| 표현 | 역할 | 권위 |
|---|---|---|
| 명시적 mechanochemical simulator | 힘 경로와 사건을 실행하는 reference instrument | 물리 런타임 기준 |
| 고전 확률동역학 | 위치·형태·결합·파라미터 불확실성을 분포로 전파 | 기준 물리의 확률 표현 |
| 텐서 네트워크 | 고차원 분포·다축 텔레메트리·파라미터 응답의 조건부 저랭크 표현 | 오차가 계측된 reduced representation |
| 신경망 | 관측 역문제, solver 가속, surrogate, 이상 탐지, 실험 선택 | trust region 안의 학습 도구 |
| 영상 AI | 실제 현미경과 가상 세포 사이의 양방향 observation interface | 실험 관측과의 연결 |

### 양자라는 말의 정확한 위치

- 세포 노드와 막의 위치는 **복소 양자 진폭**이 아니다.
- 세포 스케일의 기준 동역학은 과감쇠 고전 SDE, KMC, PDE/FEM이다.
- Fokker–Planck, Koopman/Perron–Frobenius, Schrödinger bridge, MPS/TT, instanton/WKB 등
  양자·통계물리에서 발전한 수학은 **quantum-inspired algorithm**으로 사용할 수 있다.
- 실제 QM/MD는 ATP 반응에너지, 결합 장벽, 국소 free-energy landscape처럼 분자 수준에서
  계산 가능한 prior를 만드는 별도 레인이다.
- 외부 논문과 제품 설명에서 `quantum cell simulator`라고 부르지 않는다. 물리적 양자효과를
  직접 계산하지 않는 결과에는 `quantum-inspired` 이상을 주장하지 않는다.

### 중심 실험

당분간 실험 프로그램의 중심은 wet lab이 아니라 **재현 가능한 simulation sandbox**다.

1. 물리법칙·수치법·표현법을 서로 다른 형식으로 교차검증한다.
2. 가상 perturbation으로 인과축을 분리한다.
3. 세포 archetype 사이에서 공통 connector law와 cell-state effect를 분리한다.
4. 합성 현미경 영상을 만들어 image-to-state 역문제를 닫는다.
5. 실제 실험 제안은 sandbox에서 식별 가능성과 예상 정보량이 증명된 뒤에만 작성한다.

---

## 1. 북극성, 성공 정의, 금지선

### 1.1 북극성

하나의 세포를 정답 궤적으로 재현하는 것이 목적이 아니다. 다양한 세포가 주어진 환경과
개입 아래에서 만들 수 있는 **상태·궤적·관측량의 분포**를 예측하고, 왜 그런 분포가 나오는지
connector 수준에서 설명하는 것이 목적이다.

### 1.2 플랫폼 성공 정의

플랫폼은 다음 질문에 답할 때 성공한다.

1. 이 세포 상태를 구성하는 컴포넌트와 실제 힘 전달 경로는 무엇인가?
2. 관측과 양립하는 분자 파라미터의 사후분포는 무엇인가?
3. 이 계는 fixed point, stochastic NESS, quasicycle, limit cycle, chaos 중 어느 레짐인가?
4. 특정 약물·유전자·환경 개입은 어느 connector와 population을 통해 phenotype을 바꾸는가?
5. 실제 영상과 맞지 않는 부분은 파라미터 불확실성인가, 누락된 메커니즘인가?
6. 기존 세포에서 배운 법칙 중 어떤 것이 다른 세포 archetype으로 전이되는가?
7. 다음 simulation intervention 중 무엇이 불확실성을 가장 많이 줄이는가?

### 1.3 변하지 않는 금지선

- **Co-location is not connection.** 연결은 typed connector로만 존재한다.
- 물리 connector는 양방향이며 adjoint closure를 갖는다.
- kinetic state는 accepted physical step에서만 commit한다.
- force-accepted 궤적은 accepted physical trajectory가 아니다.
- surrogate 결과는 reference physics 결과인 것처럼 기록할 수 없다.
- 텐서 압축 오차와 neural uncertainty는 반드시 artifact에 포함한다.
- gate는 전체 field/tensor를 저장한 뒤, 사전 선언된 projection으로 만든 scalar에만 적용한다.
- 한 이미지에서 보이지 않는 분자 파라미터를 point estimate로 단정하지 않는다.
- 한 세포 archetype에서 성립한 모델을 다른 archetype의 보편법칙으로 승격하지 않는다.
- 실제 양자 coherence를 계산하지 않는 모델에 양자역학적 물리효과를 주장하지 않는다.
- simulation sandbox에서 자기 모델이 만든 정답으로 자기 모델을 채점하지 않는다.

---

## 2. 전체 시스템 아키텍처

```text
                           ┌──────────────────────────┐
 Literature / KB ─────────▶ typed priors + mechanisms│
 Omics / atlas ───────────▶ cell identity priors      │
 Microscopy ──────────────▶ geometry/state likelihood │
 Force/time-series ───────▶ dynamical likelihood      │
                           └────────────┬─────────────┘
                                        ▼
                 ExperimentContext ───▶ CellStatePosterior(t0)
                   (external, hashed)       │
                                        │
                        ┌───────────────┴───────────────┐
                        ▼                               ▼
          Explicit accepted-step runtime       Reduced representations
          FEM/PDE/KMC/connectors                TT/MPS/POD/ensembles
                        │                               │
                        └───────────────┬───────────────┘
                                        ▼
                           Distributional telemetry
                 time × entity × component × mode × seed × θ
                                        │
              ┌─────────────────────────┼─────────────────────────┐
              ▼                         ▼                         ▼
       scalar gate projections   neural inference/surrogate   image renderer
              │                         │                         │
              └─────────────────────────┼─────────────────────────┘
                                        ▼
                         intervention-response posterior
                                        │
                                        ▼
                           next in-silico experiment
```

### 2.1 Reference dynamics layer

현재 explicit Warp/CUDA 엔진은 폐기 대상이 아니라 reference instrument다. 다만 reference가 되려면
먼저 다음을 만족해야 한다.

- accepted-step predicate가 실제 production driver에 배선된다.
- alias 없이 candidate/rollback/commit을 증명할 수 있다.
- inner mechanical admissibility와 outer dynamical stationarity를 구분한다.
- connector별 양방향 force, work, event count가 lossless telemetry에 남는다.
- 수치 진동과 생물학적 진동을 `dt`, tolerance, iteration budget, seed에 대한 불변성으로 구분한다.

### 2.2 Cell State 수학 계약

이 절은 `CellState`의 계획서 수준 canonical definition이다. 현재 런타임의 고정 6축 객체나
현재 component census를 보편 정의로 승격하지 않는다.

#### 2.2.1 의미론적 경계

하나의 cell-state schema manifest를 \(M\), 그 manifest가 허용하는 typed owner graph를

\[
G_M=(V_M,E_M)
\]

로 둔다. \(V_M\)은 state-owning component instance, \(E_M\)은 두 개 이상의 endpoint 사이에
존재하는 state-owning connector instance다. `component`는 오늘의 막·피질·핵 목록에 닫힌
열거형이 아니며, 명시적인 ownership과 transition law를 가진 새로운 세포내 state owner를 계속
추가할 수 있는 문법 항이다. `connector`도 오늘의 기계적 edge에 닫히지 않으며, 기계적·화학적·
전기적·수송·정보적 coupling을 typed law로 추가할 수 있다.

다음 객체는 서로 구분한다.

| 객체 | 수학적 역할 | Cell State에 포함되는가 |
|---|---|---|
| \(M\) | subject annotation, 허용 graph와 field schema | 상태의 schema이지 동적 값은 아님 |
| \(G_t\) | 시점 \(t\)의 instantiated owner/connector topology | 포함 |
| \(S_{\mathrm{cell}}(t)\) | 모든 owner의 accepted 내부 상태 | **Cell State 본체** |
| \(U(t)\) | 외부 조건, boundary forcing, protocol | **포함하지 않음** |
| \(\theta\) | 동역학 law의 molecular parameter | state와 공동 추론할 수 있으나 state 값은 아님 |
| \(Y(t)\) | observation operator가 만든 관측 | 포함하지 않음 |
| \(\Pi_t\) | 상태·파라미터·레짐·모델에 대한 posterior | Cell State에 대한 지식 |

ECM, 유체, organelle, signalling module 같은 이름만 보고 내부/외부를 고정하지 않는다. 특정
모델 경계 안에서 dynamics와 state ownership을 갖는다면 \(V_M\)의 component이고, boundary
condition 또는 imposed input으로만 주어진다면 \(U(t)\)다. 이 판정은 이름이 아니라 manifest의
system-boundary·ownership contract가 한다.

#### 2.2.2 Owner-local state space

각 owner \(a\in V_M\cup E_M\)는 열린 field 집합 \(F_a\)를 가진다. field 하나의 schema는 최소한

\[
\sigma_{a,f}
=
(\text{owner ref},\text{field id},\text{kind},\text{support},\text{unit},
\text{basis},\text{time scale},\text{conservation class},\text{provenance})
\]

를 선언한다. `kind`는 closed enum이 아니라 versioned type URI다. 따라서 이후의 새로운
biochemical, electrical, transcriptional, epigenetic 또는 아직 이름 붙지 않은 state field가
기존 최상위 축을 깨지 않고 들어온다.

owner \(a\)의 국소 상태공간은

\[
\mathcal X_a=\prod_{f\in F_a}\mathcal X_{a,f}
\]

이며, \(\mathcal X_{a,f}\)는 필요에 따라 다음 중 하나 또는 그 조합이다.

- 연속 유한차원 공간: position, orientation, concentration, internal strain
- 이산 공간: binding state, chemical state, activation state
- 함수공간/mesh field: membrane, fluid, electrochemical or reaction–diffusion field
- 유한 집합·point process: 명시적 molecule/filament/site population
- topology state: active edge set, ownership ledger, generation/epoch
- memory state: Markov closure에 필요한 internal clock, dwell time, hereditary variable

이 목록은 schema 예시이며 보편적인 고정 compartment 목록이 아니다.

component는 자신의 population과 물리 state를 유일하게 소유한다. connector는 endpoint의 state를
복제하지 않고, 결합 자체에 필요한 occupancy·mapping·kinetic clock·committed topology만 소유한다.
완전히 stateless한 connector는 singleton state space를 가지며, 상태가 없다는 이유로 endpoint
state를 자기 필드처럼 읽어 소유하지 않는다.
whole-cell state는

\[
\mathcal X_{G_t}
=
\left(
\prod_{v\in V(G_t)}\mathcal X_v
\right)
\times
\left(
\prod_{e\in E(G_t)}\mathcal X_e
\right)
\cap \mathcal C_{G_t}
\]

의 원소다. \(\mathcal C_{G_t}\)는 ownership disjointness, endpoint compatibility, conservation,
finite-state, topology-generation, connector adjoint 조건을 모은 admissible constraint set이다.

\[
S_{\mathrm{cell}}(t)
=
\operatorname{Compose}_{G_t}
\left(
\{S_v(t)\}_{v\in V(G_t)},
\{S_e(t)\}_{e\in E(G_t)}
\right)
\in\mathcal X_{G_t}
\]

`Compose`는 새 authoritative array를 만드는 연산이 아니다. 한 accepted transaction에 속한
owner state들의 namespaced reference와 content identity를 묶는 비소유 view/snapshot이다.

#### 2.2.3 가변 population과 topology

세포의 population count, binding graph, connector instance 수는 시간에 따라 바뀔 수 있으므로
하나의 고정 Cartesian product만으로 전체 상태공간을 정의하지 않는다. manifest \(M\)이 허용하는
instantiated graph/topology sector의 집합을 \(\Gamma(M)\)라 하면

\[
\mathcal X_M
=
\bigsqcup_{g\in\Gamma(M)}\mathcal X_g
\]

이다. 여기서 \(\bigsqcup\)는 서로 다른 population/topology sector를 섞지 않는 disjoint union이다.
polymerization, severing, assembly, junction formation 같은 사건 \(r\)은

\[
R_r:\mathcal X_g\rightarrow\mathcal X_{g'}
\]

인 typed jump map이다. 새 owner나 connector가 생길 때 stable identity, unique ownership,
generation counter와 rollback state가 함께 이동해야 한다. topology sector가 다른 두 상태를
같은 고정 길이 벡터로 zero-padding해 동일한 상태라고 부르지 않는다.

#### 2.2.4 Hybrid mechanochemical dynamics

하나의 sector \(g\) 안에서 연속·확률·jump dynamics를 공통 generator로 쓴다. 충분히 매끄러운
관측함수 \(\varphi\)에 대해

\[
\begin{aligned}
\mathcal L_{\theta,U(t)}\varphi(x)
={}&
b(x;\theta,U)\cdot\nabla\varphi(x)
+\frac12\operatorname{tr}
\left[D(x;\theta,U)\nabla^2\varphi(x)\right] \\
&+\sum_r\lambda_r(x;\theta,U)
\left[\varphi(R_r x)-\varphi(x)\right]
+\mathcal L_{\mathrm{field}}\varphi(x).
\end{aligned}
\]

- \(b\): deterministic mechanics, transport, constitutive drift
- \(D=\sigma\sigma^\top\): thermal 또는 명시적으로 모델된 stochastic forcing
- \(\lambda_r\): KMC/event hazard
- \(R_r\): topology/population/chemical-state jump
- \(\mathcal L_{\mathrm{field}}\): PDE/FEM/continuum field evolution
- \(U(t)\): 외부 입력; state coordinate가 아님

graph locality를 보존하도록 generator는

\[
\mathcal L
=
\sum_{v\in V}\mathcal L_v
+\sum_{e\in E}\mathcal L_e
+\sum_{h\in H_{\mathrm{local}}}\mathcal L_h
\]

로 분해한다. \(\mathcal L_v\)는 component-local law, \(\mathcal L_e\)는 connector endpoint
coupling, \(\mathcal L_h\)는 둘 이상의 connector가 하나의 물리 joint를 구성할 때의 명시적인
hyperedge/group law다. `co-location`은 이 합에 항을 만들지 않는다.

기계적 connector \(e=(a,b)\)의 connector-quadrature traction을 \(f_e\), endpoint
interpolation을 \(J_a,J_b\)라 두면

\[
f_{e,a}+f_{e,b}=0,\qquad
F_{e\to a}=J_a^\top f_e,\qquad
F_{e\to b}=-J_b^\top f_e
\]

이고 adjoint virtual-work identity는

\[
\langle F_{e\to a},\delta x_a\rangle
+\langle F_{e\to b},\delta x_b\rangle
=
\langle f_e,J_a\delta x_a-J_b\delta x_b\rangle
\]

다. passive connector의 endpoint power는 저장에너지 변화와 dissipation으로 닫히고, active
connector는 추가 energy-injection channel과 source를 별도로 기록한다. 내부 elastic connector의
work 자체를 0으로 두지 않는다. 이 식들은 총합만 맞는 global cancellation이 아니라 connector
instance별로 평가한다.

#### 2.2.5 Accepted-step transition kernel

물리 generator와 수치 candidate를 동일시하지 않는다. 현재 accepted state \(x\)에서 candidate
\(c\)를 제안하는 kernel을 \(Q_{\Delta t}(dc\mid x,\theta,U)\), device predicate를
\(A(c,x)\in\{0,1\}\), commit map을 \(C(c)\), bit-exact rollback을 \(R(c,x)=x\)라 두면 실제
기록되는 one-step kernel은

\[
K_{\Delta t}(dy\mid x)
=
\int Q_{\Delta t}(dc\mid x)
\left[
A(c,x)\,\delta_{C(c)}(dy)
+(1-A(c,x))\,\delta_x(dy)
\right].
\]

clock, RNG, topology epoch와 connector kinetic candidate도 같은 predicate로 commit/rollback한다.
따라서 accepted timeline은 attempted-step timeline의 단순 subsampling이 아니라 이 kernel이 만든
별도 확률과정이다. tensor/Koopman/Perron–Frobenius 분석은 \(Q_{\Delta t}\)가 아니라
runtime-attested \(K_{\Delta t}\)와 accepted receipt chain에 결합한다.

#### 2.2.6 Cell State의 충분성 기준

Cell State는 “기록 가능한 값을 모두 모은 것”이 아니라 미래 예측에 충분한 내부 상태여야 한다.
같은 미래 외부 입력 \(U_{t:}\)와 parameter law 아래에서

\[
\Pr(Y_{t:}\mid \mathcal H_{<t},S_{\mathrm{cell}}(t),U_{t:},\theta,M)
=
\Pr(Y_{t:}\mid S_{\mathrm{cell}}(t),U_{t:},\theta,M)
\]

가 성립해야 한다. \(\mathcal H_{<t}\)를 추가했을 때 예측이 체계적으로 달라지면 현재 state
schema가 Markov-sufficient하지 않다는 뜻이다. 이 실패를 parameter retuning으로 숨기지 않고
새 memory field, internal clock, component 또는 connector law로 확장한다. non-Markovian model을
직접 쓰는 경우에는 필요한 history functional과 memory kernel을 state contract에 명시한다.

#### 2.2.7 실제 상태, 분포, posterior

실제 세포 상태 \(X_t\)와 그 상태에 대한 지식 \(\Pi_t\)를 구분한다. 관측 \(y_{0:t}\), 외부 조건
\(U_{0:t}\), intervention \(I_{0:t}\), evidence \(\mathcal E\)가 주어졌을 때

\[
\Pi_t
=
p\!\left(
dX_t,d\theta,dR_t,dH
\mid
y_{0:t},U_{0:t},I_{0:t},\mathcal E,M
\right)
\]

로 둔다.

- \(X_t\): whole-cell internal state
- \(\theta\): molecular-parameter posterior
- \(R_t\): fixed point/NESS/cycle/chaos 등 dynamical-regime uncertainty
- \(H\): graph 또는 governing-law hypothesis

초기 구현은 하나의 \(H\)에 condition할 수 있다. 이후 missing-law 경쟁모델은

\[
\Pi_t=\sum_{h\in\mathcal H}w_h\,\Pi_t^{(h)}
\]

인 model mixture로 확장한다. 지원하지 않는 topology를 기존 graph의 parameter 폭으로 흡수하지
않는다.

확률은 하나가 아니라 출처별로 분해한다.

```text
total uncertainty
  ├─ thermal / explicitly modelled external stochastic forcing
  ├─ stochastic molecular events
  ├─ quenched structural variability
  ├─ cell-to-cell biological variability
  ├─ experimental measurement noise
  ├─ parameter epistemic uncertainty
  ├─ model/graph uncertainty
  └─ numerical / representation error
```

각 항은 convolution, mixture, conditional branch 또는 bias bound 중 어떤 방식으로 결합되는지
선언한다. 이 항들을 출처 없는 하나의 Gaussian 폭으로 합치지 않는다.

사용 가능한 posterior 표현은 explicit ensemble/path measure, local Gaussian, Gaussian mixture,
jump-diffusion, sequential Monte Carlo, normalizing flow, tensorized density 등이다. 표현 선택은
Cell State의 의미를 바꾸지 않으며 representation descriptor와 measured error만 바꾼다.

#### 2.2.8 Observation, intervention, inference

observation operator \(O_k\)는 state의 일부가 아니라 one-way likelihood다.

\[
Y_k\sim p_{O_k}(dy\mid X_{t_k},\theta,U_{t_k},M).
\]

Bayesian update는

\[
\Pi_{t_k}^{+}(dX,d\theta)
\propto
p_{O_k}(y_k\mid X,\theta,U_{t_k},M)
\Pi_{t_k}^{-}(dX,d\theta)
\]

이며, 보이지 않는 방향의 posterior를 point로 만들지 않는다.

intervention \(I(t)\)도 state field를 임의로 덮어쓰는 값이 아니다. component/connector law,
hazard, admissible topology map 또는 외부 forcing을 변경하는 typed operator로

\[
\mathcal L_{\theta,U}^{\,I}
=
\mathcal I\!\left[\mathcal L_{\theta,U}\right]
\]

를 정의한다. intervention 이전과 이후의 state continuity, onset/offset schedule, off-target
hypothesis와 accepted-step commit behavior를 기록한다.

#### 2.2.9 전체 generative/path model

위 조각들을 하나의 joint law로 닫는다. accepted 시점 \(t_0,\ldots,t_K\)에서

\[
\begin{aligned}
&p(H,\theta,X_{0:K},R,Y_{\mathcal O}
\mid U_{0:K},I_{0:K},\mathcal E,M) \\
&\quad=
p(H\mid\mathcal E,M)\,
p(\theta\mid H,\mathcal E,M)\,
p(X_0\mid\theta,H,M) \\
&\qquad\quad\times
\prod_{k=0}^{K-1}
K_{\Delta t_k,\theta,U_k}^{H,I_k}
(dX_{k+1}\mid X_k)
\prod_{k\in\mathcal O}
p_{O_k}(Y_k\mid X_k,\theta,U_k,M)\,
p_R(R\mid X_{0:K},\mathcal C_R).
\end{aligned}
\]

\(\mathcal C_R\)는 사전 선언된 regime evidence contract다. `CellStatePosterior(t_k)`는 이 joint
law의 필요한 marginal/conditional이고, intervention-response posterior는 intervention \(I\)가
다른 두 path measure의 관측 차이다.

\[
\Delta O
=
O(X^{I}_{0:K})-O(X^{0}_{0:K}),
\qquad
(X^{I},X^{0})\sim\gamma.
\]

\(\gamma\)는 independent coupling, common-random-number coupling 또는 event-coupled path measure
중 어느 것인지 선언한다. intervention effect의 분산과 parameter uncertainty를 섞지 않는다.

현재 중심 객체의 모든 항은 다음 수학 객체에 대응한다.

| 계획서 항 | 수학 객체 | 비고 |
|---|---|---|
| component/connector graph | \(M,G_t,H\) | schema, instantiated topology, model hypothesis 분리 |
| population/topology | sector \(g\in\Gamma(M)\) | 가변 count/topology는 disjoint union |
| geometry와 내부 field | owner-local \(X_{a,f}\) | component가 유일 소유 |
| connector binding/coupling | connector-local \(X_e,\mathcal L_e\) | endpoint state를 복제하지 않음 |
| molecular parameters | \(p(\theta\mid\cdot)\) | state와 공동 추론하되 state field가 아님 |
| accepted dynamical regime | \(p_R(R\mid X_{0:K},\mathcal C_R)\) | trajectory와 contract에 조건부 |
| intervention response | paired path measure와 \(\Delta O\) | intervention operator를 명시 |
| observation operators | \(p_O(Y\mid X,\theta,U,M)\) | one-way likelihood |
| evidence provenance | \(\mathcal E\)와 source graph | prior/likelihood/applicability를 결합 |
| known model inadequacy | \(H\), missing-law/model mixture | parameter 폭으로 숨기지 않음 |
| external context | \(U\) | 조건부 우변, Cell State 밖 |
| representation error | \(\varepsilon_{\mathrm{repr}}\) | 생물학적 uncertainty와 분리 |

이 factorization은 tensor network를 적용할 수 있는 정확한 seam을 준다. time edge는 transition
factor, owner graph edge는 connector factor, observation은 likelihood factor, evidence는 prior
factor다. 서로 의미가 다른 축을 임의로 하나의 dense tensor에 넣을 이유가 없다.

#### 2.2.10 다세포·조직 확장

세포 \(i\)마다 namespaced graph \(G_i\)와 state \(S_i\)를 유지하고, 세포 사이에는 명시적
intercellular connector graph \(J\)를 둔다.

\[
G_{\mathrm{tissue}}
=
\left(\bigsqcup_i G_i\right)\cup J,
\qquad
S_{\mathrm{tissue}}
=
\operatorname{Compose}
\left(\{S_i\}_i,\{S_j\}_{j\in J},\{S_s\}_{s\in\mathrm{shared}}\right).
\]

shared state owner는 하나의 namespaced object이며 세포마다 복제하지 않는다. collective state가
single-cell state의 field를 재소유하지 않고, junction/shared-owner connector를 통해서만 결합한다.
이 정의는 세포 수가 늘어나도 single-cell Cell State의 의미를 바꾸지 않는다.

#### 2.2.11 최소 완전성 계약

어떤 `CellStateSnapshot` 또는 `CellStatePosterior`도 다음을 만족하지 않으면 불완전하다.

1. 모든 authoritative field가 정확히 한 component/connector owner를 가진다.
2. 모든 connector instance가 선언된 endpoint와 law, generation을 가진다.
3. state hash가 owner content, topology sector, accepted step과 representation을 구분한다.
4. 외부 조건은 별도 context hash로 기록되며 Cell State field로 흡수되지 않는다.
5. 누락된 field/component/connector/law는 빈 값이나 0이 아니라 first-class missing verdict다.
6. state schema는 versioned이며 새 field kind가 기존 hash 의미를 바꾸지 않는다.
7. posterior는 uncertainty source, representation kind, identifiability와 approximation error를
   보존한다.
8. observation과 intervention은 state를 소유하지 않고 typed operator로 state에 작용한다.

### 2.3 Cell State의 tensor-network 모델

텐서 네트워크는 위 상태 의미론을 대신하지 않는다. \(\mathcal X_M\), \(\Pi_t\),
\(K_{\Delta t}\), observation likelihood 또는 telemetry 중 **어떤 객체를 근사하는지 먼저
선언한 뒤** 쓰는 representation이다.

#### 2.3.1 Tensorization plan

각 tensor artifact는 계산 전에 다음 plan을 고정한다.

```text
TensorizationPlan
  subject_object: state-density / transition-operator / likelihood /
                  telemetry / parameter-response / cross-modal-map
  state_schema_hash:
  topology_sector:
  local_coordinates:
    - owner_ref
    - field_id
    - support and units
    - basis/discretization id
    - coordinate or label hash
    - truncation domain
  grouping:
  axis_order or network topology:
  symmetry/conservation sectors:
  requested_error_contract:
  fallback:
```

component 이름이나 현재 array 순서를 보고 축 의미를 추정하지 않는다. 축은
`owner_ref × field_id × basis_index`로 식별한다. entity ID, accepted time, parameter, seed,
observation channel 같은 비상태 축은 별도 role로 선언한다.

#### 2.3.2 Mixed-state basis expansion

연속·이산·field state를 tensor coefficient로 만들기 위해 각 국소 coordinate group \(q\)에
기저 \(\{\phi_{q,i_q}\}\)를 둔다.

\[
p(x)
\approx
\sum_{i_1,\ldots,i_d}
C_{i_1\cdots i_d}
\prod_{q=1}^{d}\phi_{q,i_q}(x_q).
\]

기저는 grid, finite element, spectral mode, wavelet, discrete indicator, occupation basis,
learned local basis 중 하나일 수 있다. basis truncation error와 tensor truncation error를
분리한다. full-native reference state는 그대로 보존하고, entity grouping 또는 modal reduction은
reduced representation에서만 measured-error contract 아래 수행한다.

가변 population/topology는 한 dense axis로 padding하지 않는다. sector \(g\)별 coefficient
\(C^{(g)}\)를 두어

\[
p(x)=\sum_{g\in\Gamma(M)}w_g\,p_g(x),
\qquad
p_g(x)\leftrightarrow C^{(g)}
\]

로 표현한다. topology-changing jump는 sector 간 block operator다.

#### 2.3.3 Network topology 선택

| 수학 객체 | 우선 표현 | 선택 이유 |
|---|---|---|
| 1차원 vector/series | dense | tensor overhead가 이득 없음 |
| 2차원 matrix/Jacobian | SVD 또는 sparse matrix | TT로 이름만 바꾸지 않음 |
| 자연스러운 chain 축 | TT/MPS | 순차 split과 오차 측정이 단순 |
| owner hierarchy가 tree에 가까움 | TTN/hierarchical Tucker | component 내부–connector–whole-cell 계층 보존 |
| 일반 sparse component graph | graph tensor network | endpoint locality 보존 |
| local transition/generator | MPO 또는 graph operator network | \(\sum_v\mathcal L_v+\sum_e\mathcal L_e\) 이용 |
| 고랭크 또는 global-coupled density | ensemble/sparse/latent fallback | 저랭크를 강요하지 않음 |

PEPS, MERA, DMRG 같은 이름은 대상 graph와 알고리즘이 실제로 그 구조를 사용할 때만 쓴다.
양자계 표기법을 차용해도 coefficient는 고전 확률·관측·operator의 표현이며 세포의 양자진폭이
아니다.

#### 2.3.4 Graph-aware ordering과 partition

TT를 쓸 경우 axis order는 임의 파일 순서가 아니라 owner graph에서 도출한 후보를 비교한다.

- 같은 component의 strongly coupled field를 우선 group한다.
- connector field는 가능한 한 endpoint group 사이에 둔다.
- conservation/symmetry sector를 가르는 축은 block boundary로 둔다.
- graph separator, estimated mutual information, operator coupling norm을 ordering 후보 생성에 쓴다.
- 후보를 실제 rank, byte cost, decode cost, declared-observable error로 평가한다.

rank는 ordering과 tolerance의 함수

\[
r_k=r_k(\pi,\varepsilon,\mathcal B,\mathcal P)
\]

이며 \(\pi\)는 축 순서, \(\mathcal B\)는 basis, \(\mathcal P\)는 partition이다. 따라서 bond
dimension은 우선 representation complexity이고 물리량이 아니다.

#### 2.3.5 State density와 operator의 TT/MPS 형태

chain ordering을 택한 coefficient tensor는

\[
C_{i_1\cdots i_d}
\approx
\sum_{\alpha_1,\ldots,\alpha_{d-1}}
G^{(1)}_{1,i_1,\alpha_1}
G^{(2)}_{\alpha_1,i_2,\alpha_2}
\cdots
G^{(d)}_{\alpha_{d-1},i_d,1}
\]

로 쓴다. transition/generator는

\[
K_{j_1\cdots j_d,i_1\cdots i_d}
\approx
\sum_{\beta_1,\ldots,\beta_{d-1}}
W^{(1)}_{1,j_1,i_1,\beta_1}
W^{(2)}_{\beta_1,j_2,i_2,\beta_2}
\cdots
W^{(d)}_{\beta_{d-1},j_d,i_d,1}
\]

형태의 MPO/operator network로 표현한다. 실제 예측–갱신은

\[
\tilde\Pi_{t+\Delta t}^{-}
=
\widetilde K_{\Delta t}\tilde\Pi_t^{+},
\qquad
\tilde\Pi_{t_k}^{+}
=
\operatorname{Normalize}
\left(
\tilde L_{y_k}\odot\tilde\Pi_{t_k}^{-}
\right)
\]

이다. 각 apply/compress/normalize 뒤에 mass, negativity, declared projection과 state-sector
보존을 다시 측정한다. posterior update가 rank를 폭발시키면 tolerance를 사후 완화하지 않고
mixture 분할, sector 분할, ensemble 또는 native query로 후퇴한다.

#### 2.3.6 고전 확률의 비음수성과 정규화

SVD/TT-SVD의 작은 Frobenius error는 확률 validity를 보장하지 않는다. density 표현은 최소한

\[
\int \tilde p(x)\,dx=1,\qquad
\tilde p(x)\ge 0
\]

를 검사한다. 사용할 수 있는 표현은 다음과 같다.

- 직접 density TT: reconstruction 뒤 negativity/mass를 측정
- nonnegative-core factorization: 최적화 실패와 local minimum을 보고
- square-root density \(p=\psi^2/Z\): 비음수성은 보장하지만 \(\psi\)를 양자 파동함수로 해석하지 않음
- mixture of tensor networks: multimodal/regime sector를 한 core에 억지로 합치지 않음
- sample-to-TN fit: held-out samples와 observable calibration으로 검증

조용한 clipping이나 renormalization은 원래 approximation error를 숨기므로 금지한다. repair가
필요하면 repair 전후 artifact와 추가 error budget을 별도로 기록한다.

#### 2.3.7 Tensor error contract

tensor acceptance는 단일 Frobenius norm으로 닫지 않는다. 최소 error vector는

\[
\boldsymbol\varepsilon_{\mathrm{TN}}
=
(\varepsilon_{\mathrm{basis}},
\varepsilon_{\mathrm{repr}},
\varepsilon_{\mathrm{mass}},
\varepsilon_{\mathrm{neg}},
\varepsilon_{\mathrm{obs}},
\varepsilon_{\mathrm{dyn}},
\varepsilon_{\mathrm{tail}},
\varepsilon_{\mathrm{byte/decode}})
\]

다.

- \(\varepsilon_{\mathrm{basis}}\): local discretization/modal truncation
- \(\varepsilon_{\mathrm{repr}}\): dense/ensemble reference에 대한 representation error
- \(\varepsilon_{\mathrm{mass}},\varepsilon_{\mathrm{neg}}\): 확률 정규화와 음수 질량
- \(\varepsilon_{\mathrm{obs}}\): 사전 선언된 observation/projection별 오차
- \(\varepsilon_{\mathrm{dyn}}\): 여러 step rollout 또는 invariant measure 오차
- \(\varepsilon_{\mathrm{tail}}\): rare event와 low-probability region 보존 실패
- \(\varepsilon_{\mathrm{byte/decode}}\): 실제 저장 byte와 decode 비용

각 허용치는 run 전에 정하고 artifact에 복사한다. global norm이 통과해도 하나의 필수 projection,
tail probability, conservation sector 또는 accepted-step identity가 실패하면 tensor
representation은 거부한다.

#### 2.3.8 Tensorization 우선순위

1. owner-semantic telemetry:
   `accepted time × owner ref × field/basis × entity × mode × seed`
2. graph-conditioned parameter response:
   `law parameter × intervention × observation × external-context coordinate`
3. 낮은 effective rank를 실제로 보이는 discrete CME/topology sector
4. reduced state coordinate와 parameter가 결합된 Fokker–Planck density
5. local generator/transition의 MPO 또는 graph operator network
6. observation likelihood와 mechanistic state 사이의 cross-modal tensor
7. multi-cell graph에서 repeated cell-local core와 junction core의 공유 가능성

`external-context coordinate`는 response tensor의 조건 축일 수 있지만 Cell State 축은 아니다.
telemetry tensor, state-density tensor, response tensor, image latent tensor는 서로 다른 artifact
kind이며 축 이름이 겹친다는 이유로 하나의 텐서로 합치지 않는다.

#### 2.3.9 실패와 fallback

다음 중 하나면 tensor network가 기본 표현이 아니다.

- 축 순서에 따라 결론이 바뀐다.
- requested tolerance에서 compression ratio가 1 이하이다.
- connector/global coupling으로 rank가 지속적으로 증가한다.
- rare-event mass 또는 필수 관측 projection이 global norm보다 먼저 실패한다.
- topology sector를 합칠 때 multimodality가 사라진다.
- operator apply 뒤 rank 증가를 threshold 사후 변경 없이는 제어할 수 없다.
- basis error와 tensor error를 분리할 수 없다.

fallback 순서는 문제 종류에 따라 dense/SVD, sparse operator, explicit ensemble, sector mixture,
adaptive latent coordinate, native evaluation 중에서 선택한다. fallback도 같은 state-schema hash,
observation contract와 provenance를 유지한다.

#### 2.3.10 State-to-tensor compiler 순서

향후 tensor backend는 다음 순서를 건너뛸 수 없다.

1. `CellStateManifest`에서 owner graph와 field schema를 읽는다.
2. `ExperimentContextManifest`는 별도 conditioning input으로 고정한다.
3. 모든 field를 `owner_ref × field_id`로 namespacing하고 local state space를 만든다.
4. admissible topology/population sector와 conservation block을 생성한다.
5. 각 local coordinate의 basis와 basis-error oracle을 선언한다.
6. component law, connector law, jump map을 local factor/operator로 compile한다.
7. graph separator와 measured coupling으로 TT order, TTN 또는 graph-network 후보를 만든다.
8. dense/toy/ensemble oracle에서 rank와 전체 error vector를 측정한다.
9. prediction–likelihood update–compression을 실행하고 매 단계 probability와 projection을 감사한다.
10. 실패하면 선언된 fallback을 실행하고 실패한 sector/field/observable을 verdict에 남긴다.

필수 artifact chain은 다음과 같다.

```text
state_schema_manifest
  + experiment_context_manifest
  + topology_sector_manifest
  + tensorization_plan
  + basis_error_report
  + rank/order/partition_report
  + tensor_error_report
  + projection_replay_report
  + fallback_or_acceptance_verdict
```

이 chain이 없으면 생성된 core 파일은 Cell State의 tensor representation이 아니라 의미가
검증되지 않은 숫자 배열이다.

### 2.3a Operator/spectral laboratory

양자·카오스·확률론에서 가져올 수 있는 수학을 한 실험실에 두되, 서로 같은 객체라고 부르지 않는다.

- Perron–Frobenius: 상태분포 전파
- Koopman: 관측함수 전파와 저차원 좌표 탐색
- generator spectrum: relaxation mode와 metastability
- non-Hermitian decomposition: detailed-balance 위반의 후보 진단
- Lyapunov/Floquet: 예측 지평과 주기궤도 안정성
- Schrödinger bridge: 두 관측분포 사이의 최소변형 확률경로 추론
- path action/rare-event methods: rupture, escape, detachment 확률

이들은 공통 telemetry와 state manifest를 사용하지만, `tangent operator`, `score martingale`,
`transfer operator`, `quantum Hamiltonian`을 하나로 동일시하지 않는다. 각 formalism은 analytic
positive/negative control과 applicability statement를 따로 가진다.

### 2.4 Neural layer

| 모델 | 입력 → 출력 | 우선도 | 물리 권위 |
|---|---|---:|---|
| learned preconditioner | operator/state → solver acceleration | 최우선 | 답을 바꾸지 않아야 함 |
| heterogeneous equivariant GNN | component/connector graph → latent/update proposal | 높음 | reference 검증 필요 |
| mesh GNN | surface/filament mesh → next-state proposal | 높음 | rollout trust region 필요 |
| neural operator | boundary/field/parameter → PDE field | 중간 | native PDE 비교 필요 |
| neural posterior estimator | image/force/omics → `p(θ,z₀|y)` | 높음 | posterior calibration 필요 |
| regime classifier | lossless trace → regime probabilities | 높음 | analytic/synthetic controls 필요 |
| model-inadequacy detector | observation residual → missing-law candidates | 중간 | 제안만 가능 |
| intervention policy | posterior → next simulation point | 높음 | information-gain 검증 |

네트워크가 예측하는 force는 임의 벡터가 아니라 connector contract가 요구하는 adjoint pair 또는
potential/work-compatible representation이어야 한다. active connector는 energy injection channel을
별도로 선언한다.

### 2.5 Image layer

영상 AI는 다섯 방향으로 사용한다.

1. **Segmentation/tracking:** membrane, nucleus, organelle, adhesion, filament proxy 추출
2. **Image → initial state:** geometry와 관측 가능한 population에 대한 likelihood
3. **State → image:** PSF, noise, stain, exposure, missing labels를 포함하는 synthetic microscopy
4. **Perturbation prediction:** intervention 후 image distribution 예측
5. **Residual diagnosis:** 실제 영상과 합성 영상의 차이를 누락된 physics 후보로 변환

실제 이미지와 합성 이미지는 같은 encoder와 같은 batch/confounder correction을 거친다.
이미지 생성 품질은 시각적 사실성이 아니라 held-out morphology statistics, calibration,
perturbation retrieval로 평가한다.

### 2.6 Knowledge/evidence layer

문헌은 point parameter를 공급하지 않고 typed prior와 applicability 조건을 공급한다.

```text
SourceEvidence
  → KnowledgeClaim
  → Component / Connector
  → ParameterPrior
  → CellState applicability
  → ObservationOperator
  → SimulationRun
  → Posterior
  → Retraction / model inadequacy
```

새 provenance 종류:

- `SOURCED`
- `DERIVED`
- `COMPUTED_QM_MD`
- `INFERRED_POSTERIOR`
- `CALIBRATION_ONLY`
- `CONVENIENCE`
- `PI_GAP`

---

## 3. Cell State를 다루는 방식 — 고정 분류가 아니라 확장 문법

### 3.1 CellStateManifest

`CellStateManifest`는 현재 세포종이나 component 목록을 최상위 필드로 고정하는 설문지가 아니다.
상태의 subject annotation, owner graph, state-field schema와 적용 한계를 content-addressed하게
묶는 열린 manifest다. species, lineage, disease label 같은 descriptor는 ontology-backed
annotation으로 들어오며 보편적인 state axis가 되지 않는다.

```yaml
cell_state_manifest:
  schema_id:
  subject_annotations:
    # open ontology URI -> typed value; no closed species/cell-line axis
  system_boundary:
    included_state_owners:
    external_reference_policy:
  owner_graph:
    components:
      - owner_ref:
        type_uri:
        state_schema_ref:
        population_schema_ref:
    connectors:
      - connector_ref:
        type_uri:
        endpoint_refs:
        state_schema_ref:
        law_ref:
        commit_policy:
  state_fields:
    - owner_ref:
      field_id:
      kind_uri:
      support:
      unit:
      basis_or_discretization:
      time_scale:
      conservation_class:
      provenance:
  admissible_sectors:
    population_and_topology_schema:
    invariants:
  applicability:
    missing_state_owners:
    missing_connectors:
    missing_fields:
    missing_laws:
    unsupported_claims:
  provenance:
```

환경, imposed boundary, intervention와 measurement protocol은 Cell State 안에 넣지 않고 별도
`ExperimentContextManifest`로 둔다.

```yaml
experiment_context_manifest:
  context_id:
  external_inputs:
  boundary_conditions:
  intervention_operators:
  observation_operators:
  schedules:
  protocol_provenance:
```

posterior 또는 run artifact는 `cell_state_manifest_hash`와 `experiment_context_manifest_hash`를
각각 가진다. 따라서 같은 Cell State schema를 여러 외부 조건에서 시험할 수 있고, 다른 조건에서
얻은 결과가 같은 state인 것처럼 합쳐지는 것도 막는다. context manifest가 Cell State의 일부가
아니더라도 posterior의 조건과 provenance에는 반드시 남는다.

특정 대상이 component인지 external context인지는 이름으로 정하지 않는다. system boundary 안에서
authoritative state와 transition law를 갖는 대상은 owner graph에 들어가고, 외부에서 주어지는
대상은 context manifest에 남는다. 동일한 대상도 연구 질문과 해상도에 따라 다른 model boundary를
가질 수 있으나, 한 artifact 안에서는 경계가 고정되고 hash된다.

### 3.2 계층적 파라미터 모델

파라미터의 계층을 species·cell line·ECM 같은 고정 이름의 additive table로 만들지 않는다.
parameter \(\theta_a\)는 먼저 component/connector law \(a\)에 붙고, manifest와 context에서
schema-versioned feature map을 읽는다.

\[
\theta_{a,i}
=
\theta_{\mathrm{law}(a)}
+B_a\,\phi_S(M_i,S_i)
+C_a\,\phi_U(U_i)
+b_{a,i},
\]

여기서

- \(\theta_{\mathrm{law}(a)}\): 여러 instance에서 공유할 수 있는 molecular/constitutive law
- \(\phi_S\): state schema와 현재 internal state에서 읽는 typed feature
- \(\phi_U\): 별도 external-context manifest에서 읽는 feature
- \(b_{a,i}\): instance 또는 experiment random effect

feature가 존재하지 않거나 law applicability가 성립하지 않으면 0으로 대체하지 않고
`missing-feature` 또는 `missing-law`를 낸다. pooling 구조 자체의 posterior는

\[
p(\theta,\beta,H\mid Y,U,M)
\propto
p(Y\mid \theta,H,U,M)\,
p(\theta,\beta,H\mid \mathcal E,M)
\]

로 추론한다. 어떤 descriptor가 실제로 parameter variation을 설명하는지는 사전 고정된 cell-type
delta가 아니라 hierarchical shrinkage, held-out transfer와 identifiability가 결정한다.

이 구조가 있어야 새로운 세포 archetype이 기존 표에 새 열을 추가하는 작업으로 축소되지 않는다.
공통 connector law는 graph 전반에서 부분 pooling하고, state-dependent modulation·population·
architecture·외부 forcing 효과는 서로 다른 조건부 항으로 남긴다.

### 3.2a 현재 타입의 migration 의미

이 수학 계약은 현재 타입을 이름만 바꿔 통과시키지 않는다. 승인 후 다음 delta를 명시적으로
처리해야 한다.

| 현재 타입/필드 | 현재 역할 | 계획서상 migration |
|---|---|---|
| `ac.engine.CellState` | 고정 6축 slow-context conditioner | `CellRegulatoryState` 성격의 adapter로 한정; whole-cell state로 사용 금지 |
| `virtual_cell.CellStateManifest.environment` | state manifest 안의 external context | 별도 `ExperimentContextManifest`로 이동, 두 hash를 envelope가 결합 |
| `virtual_cell.MechanisticState` | observation이 보는 host-owned field snapshot | owner-semantic `CellStateSnapshot`의 observation view로 한정 |
| `PosteriorGraph` | grammar와 connector/population wiring | `StateSchemaManifest` 및 topology-sector identity와 결합 |
| `MolecularParameterPosterior` | parameter uncertainty | 유지하되 whole-state measure \(\Pi_t\)의 한 marginal로 배치 |
| `CellStatePosterior` | graph·parameter·regime·operator를 묶는 record | owner-state measure/path marginal과 context hash를 추가하는 schema revision |

기존 artifact를 다시 해석해 hash를 바꾸지 않는다. 구 schema는 version을 유지하고,
`legacy_manifest_hash → state_schema_hash + context_hash` migration record를 별도 artifact로
만든다. 기존 `environment` 값을 새 Cell State field로 복사하지 않는다.

### 3.3 확장용 archetype ring

아래는 고정 우선순위가 아니라 물리 문법을 시험하는 확장 지도다.

#### Ring A — 현재 connector grammar를 직접 시험

- 비극성/극성 상피세포
- 간엽형 이동세포
- fibroblast / myofibroblast
- endothelial cell
- amoeboid immune cell: neutrophil, T cell, macrophage
- stem/progenitor cell
- 다양한 암세포 상태: EMT, 약물내성, 세포주기, 부착/부유

#### Ring B — 새 component가 필요하지만 같은 계약으로 확장

- neuron / growth cone: 장거리 MT transport와 국소 actin
- cardiomyocyte / smooth muscle: sarcomere와 주기적 active stress
- platelet: 빠른 activation과 극단적 shape transition
- red blood cell: spectrin–membrane mechanics
- adipocyte: 큰 inclusion과 cortical shell
- ciliated epithelial cell: axoneme와 active beating

#### Ring C — 문법 자체의 일반성을 시험

- plant cell: cell wall, turgor, plasmodesmata
- fungal cell: wall growth와 tip mechanics
- bacterial cell: wall synthesis와 훨씬 작은 scale
- synthetic/minimal cell: 제한된 component inventory
- organoid/tissue-resident cell: niche는 별도 context, cell–cell coupling은 명시적 intercellular
  connector graph로 필요한 문법

Ring B/C의 실패는 프로젝트 축소 사유가 아니라 `missing component law`를 발견한 결과다.

### 3.4 초기 benchmark panel의 선택 원리

초기 panel은 특정 두 암세포주에 고정하지 않는다. 다음 물리 축을 가능한 적은 archetype으로
직교하게 덮는다.

1. junction-dominant epithelial
2. stress-fiber/ECM-remodelling mesenchymal
3. cortex-dominant amoeboid
4. flow/junction-coupled endothelial
5. oscillatory contractile
6. membrane–spectrin-dominant
7. long-range cytoskeletal transport
8. differentiation/state-transition

기존 특정 breast-cancer pair는 이 panel 안의 regression benchmark 하나로만 유지할 수 있다.
프로젝트 목적, 데이터 스키마, 논문 제목, 진척도 분모를 그 pair로 정의하지 않는다.

### 3.5 진척도는 세포주 개수가 아니라 capability coverage다

```text
progress
  = capability
  × cell archetype
  × environment
  × observation modality
  × evidence level
```

한 세포에서 많은 파라미터를 출력하는 것보다, 서로 다른 물리 archetype에서 같은 connector law가
어디까지 전이되고 어디서 `missing-law`를 내는지가 더 높은 진척이다. coverage matrix의 빈칸은
실패가 아니라 다음 simulation campaign과 component project의 큐다.

---

## 4. 통합 작업축

### W0 — Charter, evidence, migration

**목적:** 새 정체성과 증거 경계를 저장소가 강제하게 한다.

**산출물**

- `CLAUDE.md`: probabilistic mechanistic virtual cell 헌장
- `ROADMAP.md`: capability/evidence stage로 재작성
- `README.md`: 외부용 개요
- `STRUCTURE.md`: physics/probability/tensor/learning/imaging/evidence 지도
- artifact provenance schema
- historical/output immutability 규칙

**Exit**

- fresh boot agent가 literal quantum, surrogate gate, 특정 세포주 중심화를 재제안하지 않는다.
- 모든 결과가 `representation`, `evidence_source`, `cell_state_manifest_hash`,
  `experiment_context_manifest_hash`를 구분해 가진다.

### W1 — Accepted reference dynamics

**목적:** 확률·AI가 배울 수 있는 유효한 reference trajectory를 만든다.

**산출물**

- accepted-step admissibility
- independent candidate/rollback/commit state
- connector dispatch/adjoint/work ledger
- numerical oscillation controls
- full native accepted trace

**Exit**

- deliberate bad balance가 reject되고 state/RNG/clock이 rollback된다.
- iteration budget과 `dt`를 바꿔도 물리 통계가 허용오차 내에서 유지된다.

### W2 — Distributional telemetry and regimes

**목적:** write-time loss를 제거하고 동역학적 수렴 대상을 올바르게 정의한다.

**산출물**

- lossless multi-index archive
- versioned scalar projections
- autocorrelation and effective sample size
- fixed point / NESS / limit cycle / quasicycle / quasiperiodic / chaos classifier
- irreversibility, response, active work observables

**Exit**

- 과거 scalar headline을 원시 trace 재사영으로 재계산할 수 있다.
- synthetic controls에서 각 regime을 혼동행렬과 함께 구분한다.

### W3 — Probabilistic state propagation

**목적:** 점 상태를 다양한 불확실성의 분포로 확장한다.

**산출물**

- `StateFieldSpec`, `StateOwnerSchema`, `CellStateSnapshot`의 representation-neutral contract
- topology/population sector와 typed jump map
- ensemble/path measure API
- jump-diffusion state
- structural/parametric/measurement/model/representation uncertainty 분리
- local closure와 mixture fallback
- filtering/data-assimilation prototype
- Markov-sufficiency/history-dependence test

**Exit**

- direct ensemble moments와 reduced propagation이 사전 선언된 오차 내에서 일치한다.
- bimodal event state를 단일 Gaussian으로 숨기지 않는다.
- history를 추가해야 예측되는 대조군은 missing-memory/state-variable verdict를 낸다.

### W4 — Tensor network laboratory

**목적:** 어디에 저랭크가 실제로 존재하는지 측정한다.

**산출물**

- versioned `TensorizationPlan`
- deterministic oracle 위 TT-cross/ALS harness
- tensorized telemetry compressor
- CME/Fokker–Planck toy systems
- parameter-response tensor
- state-density와 transition-operator TT/MPO positive/negative controls
- graph-derived partition/order 후보와 rank/tolerance/axis-order robustness report
- basis/representation/mass/negativity/observable/dynamics/tail error vector

**Exit**

- 알려진 low-rank/high-rank positive and negative controls를 모두 구분한다.
- 압축본에서 선언 projection을 재생했을 때 원본 오차 한계를 만족한다.
- 확률 질량·비음수성·topology sector와 connector-local observable을 보존한다.
- rank가 큰 경우 자동으로 ensemble/SVD/sparse representation으로 후퇴한다.

### W5 — Sensitivity, inference, and sweep replacement

**목적:** brute-force grid를 정보 기반 simulation campaign으로 바꾼다.

**산출물**

- continuous-mechanics tangent/adjoint
- stochastic-event score estimator
- event-time/saltation handling
- common-randomness and coupled-process controls
- Fisher/identifiability
- sequential design
- calibrated posterior

**Exit**

- tangent/score/saltation이 독립 finite-difference 또는 analytic oracle과 일치한다.
- regime가 다른 두 점의 차분을 sensitivity로 보고하지 않는다.
- posterior predictive가 held-out sandbox intervention을 재현한다.

### W6 — Neural acceleration and surrogate

**목적:** 답을 대체하지 않고 계산과 탐색을 가속한다.

**산출물**

- learned preconditioner
- connector-aware equivariant GNN
- membrane/fluid neural operator
- rollout error monitor
- trust-region/OOD detector
- surrogate-assisted acquisition

**Exit**

- preconditioner가 converged answer를 이동시키지 않는다.
- surrogate 밖으로 나가면 자동으로 native evaluation을 요청한다.
- surrogate-only artifact가 native evidence로 승격될 수 없다.

### W7 — Image AI and synthetic microscopy

**목적:** 세포 영상과 mechanistic state를 양방향으로 연결한다.

**산출물**

- segmentation/tracking adapters
- observation renderer: geometry → fluorescence/brightfield-like image
- PSF/noise/stain/batch model
- image encoder
- image-to-state posterior
- perturbation image prediction
- residual attribution

**Exit**

- renderer가 geometry truth를 알고도 segmentation pipeline을 속이는 trivial shortcut을 만들지 않는다.
- held-out synthetic scenes에서 state posterior가 calibrated된다.
- public perturbation images에서 batch가 아닌 perturbation signal을 회수한다.

### W8 — Multi-cell-state mechanome atlas

**목적:** 공통 law와 세포별 architecture/state effect를 분리한다.

**산출물**

- CellStateManifest + ExperimentContextManifest registry
- hierarchical posterior
- applicability matrix
- cross-cell transfer benchmark
- missing-mechanism ledger

**Exit**

- 한 archetype에서 배운 connector law가 held-out archetype에서 성공/실패하는 이유를 분해한다.
- unsupported archetype은 억지 prediction 대신 missing-law verdict를 낸다.

### W9 — Connector-centric knowledge compiler

**목적:** 문헌을 실행 가능한 prior와 반증 가능한 provenance로 바꾼다.

**산출물**

- component/connector typed labels
- cell-state applicability
- prior distribution extraction
- searched-not-found/proxy/wrong-quantity/paywall 구분
- posterior와 source evidence의 동일 row 연결
- automatic proposal, human/PI ratification boundary

**Exit**

- wrong quantity positive control이 등록 전에 차단된다.
- 같은 connector에 대한 상충 논문이 평균으로 사라지지 않고 mixture/context prior로 남는다.

### W10 — Molecular computed priors: QM/MD

**목적:** 문헌에 없는 일부 분자 파라미터를 계산 근거로 보완한다.

**가능 대상**

- ligand–receptor free-energy profile
- transition barrier와 Kramers prefactor에 필요한 curvature
- ATP/hydrolysis energetic constraints
- small-domain elasticity
- coarse-grained potential calibration

**불가능하거나 부적절한 대상**

- 세포 내 protein density
- cell-state-specific expression
- whole-cell stall force를 단순 DFT에서 직접 계산
- 관측 없이 세포종별 population 추정

**Exit**

- `COMPUTED_QM_MD` provenance와 uncertainty가 기록된다.
- 계산값을 cell-level posterior와 혼동하지 않는다.

### W11 — Multicellular and tissue composition

**목적:** 단일세포 connector grammar를 조직으로 확장한다.

**산출물**

- cell–cell junction connector
- shared ECM and fluid environment
- tissue-level conservation/work ledger
- collective migration/jamming/wetting regimes
- organoid/niche manifests

**Exit**

- single-cell 내부 force와 intercellular force가 ledger에서 분리된다.
- collective phenotype을 개별세포 parameter 하나로 환원하지 않는다.

### W12 — Platform and interactive tools

**목적:** 연구자가 세포 상태와 개입을 탐색하는 실행 환경을 만든다.

**산출물**

- run/config registry
- tensor/trace explorer
- cell-state comparison UI
- image ↔ simulation overlay
- posterior and provenance browser
- intervention sandbox
- reproducible paper artifact bundles

**Exit**

- UI에서 보이는 모든 숫자가 원 artifact, projection code, manifest, build로 역추적된다.

---

## 5. Simulation sandbox 실험 프로그램

### 5.1 샌드박스의 원칙

- 실제 실험을 흉내 내되, wet-lab protocol 제안을 바로 산출하지 않는다.
- intervention은 parameter 값을 임의 조절하는 것이 아니라 connector/component operation으로
  정의한다.
- 가능한 경우 같은 질문을 서로 다른 수학 표현으로 두 번 푼다.
- positive control, negative control, adversarial control을 모두 포함한다.
- train/calibration/held-out intervention을 run 전에 분리한다.
- 실패는 `model inadequacy`, `numerical inadequacy`, `non-identifiability`,
  `representation failure` 중 하나 이상으로 분류한다.

### 5.2 공통 실험 좌표 — Cell State와 외부 조건을 분리

```text
state schema
  owner graph / field schema / population-topology sector / internal timescale

external context
  geometry and boundary / imposed matrix or fluid condition / confinement / flow /
  external osmolarity / shared-object references

intervention operator
  target owner/law / onset-offset / off-target hypothesis / topology operation

observation operator
  modality / required state fields / protocol / likelihood / projection

representation
  point / ensemble / mixture / sparse / tensor network / surrogate

evidence
  source / applicability / identifiability / approximation error / missing law
```

current cell-cycle, polarity, EMT, activation 같은 값은 고정된 universal axis 목록이 아니다. 해당
모델이 실제 state field와 transition law를 선언할 때 owner-local schema로 들어간다. cortex,
ECM, junction 같은 이름도 이 표의 독립 “축”이 아니라 manifest가 instantiate한 owner/connector다.
외부 조건은 실험 좌표로 사용하되 Cell State tensor의 state coordinate로 혼입하지 않는다.

### 5.3 실험 패밀리

| ID | Sandbox experiment | 핵심 질문 | 반증 기준 | 실패가 여는 확장 |
|---|---|---|---|---|
| S0 | Numerical invariance battery | 진동과 분포가 수치법 산물인가? | `dt`/tol/budget에 레짐이 종속 | integrator/acceptance 연구 |
| S1 | Passive fluctuation canon | membrane/FEM thermal fluctuation이 analytic spectrum과 맞는가? | mode별 spectrum 불일치 | stochastic boundary law |
| S2 | Active-cortex regime map | activity와 turnover가 어떤 attractor를 만드는가? | regime classifier가 controls 혼동 | non-Markovian/state-variable 확장 |
| S3 | Adhesion–wetting map | adhesion/tension/ECM이 spreading transition을 어떻게 만드는가? | hysteresis/critical behavior 재현 실패 | active wetting law |
| S4 | Amoeboid–mesenchymal transition | 같은 cell state가 migration mode를 전환하는가? | connector intervention과 phenotype 분리 실패 | polarity/signaling component |
| S5 | Confinement–nuclear coupling | cortex–LINC–lamina가 통과/회복을 결정하는가? | LINC ablation이 무효 | nuclear/chromatin law |
| S6 | Osmotic–membrane–cortex response | volume, tension, reservoir가 닫히는가? | mass/work/fluctuation ledger 불일치 | channel/regulation component |
| S7 | Endothelial flow/junction | shear가 cell–cell/ECM coupling을 어떻게 바꾸는가? | flow reversal control 실패 | mechanosensing law |
| S8 | Contractile oscillator | spontaneous/forced oscillation을 분리할 수 있는가? | forcing frequency에만 peak 존재 | biochemical oscillator coupling |
| S9 | Virtual perturbation panel | drug/KO 효과가 connector-level intervention으로 설명되는가? | held-out perturbation 방향 실패 | missing target/off-target model |
| S10 | Cross-cell transfer | shared molecular law가 archetype 사이 전이되는가? | 모든 parameter 재학습 필요 | hierarchical/component grammar 확장 |
| S11 | Image inverse challenge | 영상만으로 무엇이 식별 가능한가? | calibrated posterior가 truth 미포함 | 추가 modality/observation 설계 |
| S12 | Model-inadequacy challenge | 의도적으로 빠진 law를 residual이 찾는가? | parameter shift로 누락법을 숨김 | inadequacy detector |
| S13 | Multicellular composition | single-cell law로 collective phenotype이 나오는가? | junction off/on causal contrast 실패 | tissue connector |
| S14 | Rare-event campaign | rupture, detachment, escape 확률을 계산할 수 있는가? | brute ensemble과 rare-event 법 불일치 | tensor/path-sampling 연구 |

### 5.4 가상 perturbation 문법

```text
Molecular intervention
  → target component/connector
  → affected rate/force/population distribution
  → off-target candidates
  → onset/offset schedule
  → accepted-step commit behavior
  → predicted observables
  → falsifying controls
```

예:

- myosin inhibitor는 단순히 `force = 0`이 아니라 motor availability, duty ratio, cycle kinetics 중
  어떤 항을 바꾸는지 분리한다.
- actin depolymerization은 filament density, length distribution, monomer pool, connector availability를
  함께 재정의한다.
- adhesion inhibition은 force law와 ligand occupancy를 구분한다.
- genetic perturbation은 즉시 parameter toggle이 아니라 발현/turnover timescale을 가진다.

### 5.5 wet-lab 제안으로 넘어가는 조건

실제 실험 제안은 다음을 모두 만족한 질문에만 작성한다.

1. sandbox에서 두 개 이상의 경쟁 메커니즘이 존재한다.
2. 현재 관측으로 두 메커니즘이 구분되지 않는다는 것이 정량적으로 확인됐다.
3. 추가 관측의 expected information gain이 계산됐다.
4. numerical/representation uncertainty보다 예상 효과가 크다.
5. 세포 archetype 하나에만 의존하지 않는다.

이전에는 simulation campaign과 public-data validation을 우선한다.

---

## 6. 영상·omics·문헌을 결합하는 추론 루프

### 6.1 관측 모델

```text
mechanistic state z, θ
  ├─ fluorescence renderer     → multi-channel microscopy
  ├─ brightfield renderer      → label-free shape/motion
  ├─ traction operator         → traction map
  ├─ AFM/operator protocol     → force–indentation curve
  ├─ sequencing likelihood     → composition/state prior
  └─ summary projection        → legacy scalar observable
```

서로 다른 modality는 같은 latent state를 보지만 같은 정보를 주지 않는다. 영상이 geometry와
organelle organization을 잘 제약해도 microscopic bond rate를 직접 식별하지 못할 수 있다.

### 6.2 학습 순서

1. analytic/toy simulator에서 observation operator 검증
2. accepted native simulation으로 synthetic paired data 생성
3. 이미지 encoder와 posterior estimator의 simulation calibration
4. public microscopy/perturbation data에 domain adaptation
5. simulation-to-real residual 분해
6. residual이 지속되면 parameter가 아니라 missing-law 후보로 승격

### 6.3 공개 데이터 연결

- Cell Painting/JUMP: 대규모 chemical/genetic perturbation morphology
- Allen Cell 계열: intracellular organization과 shape variation
- Human Cell Atlas/Tabula Sapiens: cell identity와 molecular-state prior
- HuBMAP: tissue spatial context
- live-cell imaging datasets: dynamics와 cell-to-cell variation

이 자료들은 reference mechanics를 대체하지 않는다. 각각 prior, observation, held-out phenotype,
domain-shift test로 사용한다.

---

## 7. 산출 가능한 연구 프로그램

### 7.1 기반 방법론

1. **Accepted stochastic mechanochemical runtime**
2. **Distributional telemetry and regime classification for active-cell simulators**
3. **Tensorized probability propagation for hybrid jump–mechanical systems**
4. **Event-aware sensitivity and posterior inference**
5. **Connector-constrained graph neural simulation**
6. **Learned preconditioning for full-native cellular mechanics**

### 7.2 영상·AI

7. **Mechanistic state inference from microscopy**
8. **Differentiable synthetic microscopy for virtual cells**
9. **Mechanistically constrained perturbation-image generation**
10. **Detecting missing cell physics from image residuals**

### 7.3 세포생물학

11. **Active-cortex dynamical-regime atlas across cell archetypes**
12. **Amoeboid–mesenchymal migration as connector reconfiguration**
13. **Adhesion–wetting–ECM transition atlas**
14. **Nuclear mechanotransduction under confinement**
15. **Flow–junction–traction coupling in endothelial states**
16. **Oscillatory actomyosin across contractile cell classes**
17. **Membrane–skeleton mechanics from epithelial cells to erythrocytes**

### 7.4 자원·플랫폼

18. **Executable Cell Mechanome Atlas**
19. **Connector-centric evidence and posterior graph**
20. **Mechanochemical simulation benchmark with negative results**
21. **Interactive in-silico intervention laboratory**
22. **Cross-cell transfer benchmark for virtual-cell models**

### 7.5 다세포·조직

23. **From single-cell connectors to collective migration**
24. **Mechanistic virtual tissue under shared ECM and flow**
25. **Organoid state inference from image and mechanical priors**

각 번호는 자동으로 한 편의 논문이 된다는 뜻이 아니다. 독립 반증 질문, 독립 dataset,
독립 검증이 있을 때만 논문 단위로 분리한다.

---

## 8. 단계별 실행

기간보다 exit condition을 기준으로 한다. 현재 native runtime 비용과 acceptance가 확정되지 않았기
때문에 달력 추정은 증거가 아니다.

### Wave 0 — 재창설과 보존

**한다**

- 이 계획을 PI가 승인한 뒤 `CLAUDE.md`, `ROADMAP.md`, `README.md`, `STRUCTURE.md`를 개정한다.
- 기존 outputs, historical docs, retractions를 불변 evidence로 고정한다.
- representation/evidence/state-schema/context-manifest를 서로 다른 hash로 정의한다.
- 기존 특정 세포주 기반 roadmap을 regression benchmark로 재분류한다.

**Exit**

- 새 세션이 프로젝트를 특정 세포주 library나 literal quantum simulator로 오해하지 않는다.

### Wave 1 — 유효한 reference trace

**한다**

- accepted-step/rollback/aliasing P0
- tensor-valued lossless telemetry
- numerical invariance battery
- non-authoritative regime observer

**Exit**

- 하나 이상의 full-native lane에서 accepted trace와 provenance-complete artifact가 나온다.

### Wave 2 — 확률·텐서 수학 실험실

**한다**

- analytic stochastic systems
- membrane fluctuation benchmark
- ensemble vs closure
- TT positive/negative controls
- owner-graph state density와 local generator MPO controls
- topology-sector jump와 posterior filtering control
- parameter tensor oracle

**Exit**

- 어떤 representation이 어떤 state schema와 failure에서 살아남는지 자동 선택할 수 있다.

### Wave 3 — 영상 폐루프

**한다**

- simulator state renderer
- synthetic microscopy
- segmentation/encoder
- image-to-state posterior
- held-out synthetic perturbation

**Exit**

- simulation truth를 모르는 분석 pipeline이 영상에서 calibrated state posterior를 회수한다.

### Wave 4 — 다세포 archetype sandbox

**한다**

- Ring A에서 직교적인 archetype panel
- hierarchical parameter sharing
- cross-cell transfer
- missing-law detection

**Exit**

- 공통 connector law와 cell-state-specific architecture를 분리한 결과가 나온다.

### Wave 5 — adaptive inference campaign

**한다**

- event-aware sensitivity
- Fisher/identifiability
- tensor/surrogate acquisition
- held-out virtual perturbations
- posterior predictive checks

**Exit**

- grid sweep 없이 simulation budget 안에서 posterior를 줄이고, 실패 시 이유를 분류한다.

### Wave 6 — 다세포·조직과 public-data bridge

**한다**

- junction/shared environment
- collective regimes
- public image/omics alignment
- simulation-to-real inadequacy map

**Exit**

- single-cell에서 배운 connector law의 조직 수준 적용 범위가 측정된다.

### Wave 7 — 선택적 wet-lab 제안

Wave 0–6이 구분하지 못한 메커니즘에 대해서만 최소 정보량 실험을 제안한다. 이 단계는 플랫폼
초기 개발의 중심이 아니다.

---

## 9. 파일 마이그레이션 계획

“모든 파일을 새 방향으로 고친다”는 의미를 **모든 파일을 재작성한다**로 해석하지 않는다.
그렇게 하면 historical evidence와 artifact가 현재 서사에 맞게 변조된다.

### 9.1 수정 대상

- boot/charter: `CLAUDE.md`
- roadmap/status contract: `ROADMAP.md`, `STATE.md`
- public description: `README.md`
- repository map: `STRUCTURE.md`
- canonical architecture/contract docs
- active code와 tests
- KB schema/refresh/query code

### 9.2 보존 대상

- committed run artifacts
- retraction/nonquotable ledger
- historical plans and audits
- archived runtime source
- raw literature corpus

보존 파일에는 새 본문을 덮어쓰지 않고 manifest에서 `historical`, `superseded`, `invalidated`,
`regression-only` 상태를 연결한다.

### 9.3 Strangler migration

```text
old explicit runtime
    │
    ├─ accepted reference adapter
    ├─ distributional telemetry adapter
    ├─ probabilistic state observer
    ├─ tensor reducer
    ├─ neural/image observation layer
    └─ new virtual-cell orchestration
```

기존 엔진을 한 번에 교체하지 않는다. 새 레인은 같은 artifact contract와 evidence ladder를 통해
하나씩 reference path에 연결한다.

---

## 10. 적대적 확장 감사

### 10.1 감사 목적

이 감사는 “너무 크니 무엇을 버릴까”를 묻지 않는다.

> 각 약점이 어느 추가 수학, 데이터, 세포 archetype, component law, 연구 프로젝트를 요구하는지
> 찾고, 확장 경로가 없는 과장만 제거한다.

감사 질문:

1. 이 주장이 실패하면 프로젝트 전체가 죽는가, 아니면 다른 representation으로 확장되는가?
2. 현재 문법이 지원하지 못하는 세포는 무엇이며, 어떤 새 component가 필요한가?
3. 한 modality의 비식별성을 다른 modality 또는 intervention이 해소할 수 있는가?
4. tensor/neural 실패를 수치로 탐지하고 native/ensemble로 후퇴할 수 있는가?
5. simulation-only 실험이 자기확증에 빠지는 지점은 어디인가?
6. 현재 한계가 독립 논문 질문으로 바뀔 수 있는가?

### 10.2 감사 결과와 강제 수정

#### A1. “확률적 위치”가 한 가지 분포라는 착각

**공격:** 열 요동, KMC, cell-to-cell variation, parameter uncertainty, measurement noise를 하나의
`Σ`로 합치면 원인과 시간척도가 사라진다.

**판정:** 단일 Gaussian universal state는 기각.

**확장:** source-factored uncertainty와 ensemble/mixture/flow/TT representation selector를 W3에
강제했다. 이 차이 자체가 uncertainty-decomposition 논문이 될 수 있다.

#### A2. full Fokker–Planck의 차원 폭발

**공격:** full-native continuous state 전체에 density grid를 만드는 것은 불가능하다.

**판정:** full density materialization 기각. 확률동역학 방향은 유지.

**확장:** path ensembles, local closures, latent collective coordinates, tensorized parameter/state
subspace를 병렬 표현으로 둔다. 어떤 표현이 어디서 깨지는지가 W4의 연구 질문이다.

#### A3. tensor rank를 물리량으로 과대해석

**공격:** rank는 axis order와 tolerance에 의존한다.

**판정:** rank=correlation length 같은 직접 동일시는 금지.

**확장:** representation complexity atlas, axis-order robustness, entanglement-like diagnostics를
별도 연구축으로 유지한다.

#### A4. neural surrogate가 새 lumped model이 될 위험

**공격:** 빠른 surrogate가 편리해지면 reference보다 더 많은 결론을 생산하게 된다.

**판정:** surrogate는 search/proposal/initialization에 한정하고 evidence source를 강제.

**확장:** learned preconditioner, trust-region monitor, native-query acquisition, surrogate failure
atlas를 별도 산출물로 둔다.

#### A5. image inverse problem의 비식별성

**공격:** 서로 다른 분자 파라미터가 같은 형태 이미지를 만들 수 있다.

**판정:** image→parameter point prediction 기각.

**확장:** posterior prediction, temporal images, traction/omics 결합, virtual interventions,
expected-information-gain observation design으로 확장한다.

#### A6. 합성영상으로 학습하고 합성영상으로 검증하는 순환

**공격:** renderer artifact를 encoder가 학습하면 완벽한 자기채점이 가능하다.

**판정:** 같은 renderer family 안의 test는 구조 테스트일 뿐 생물 검증이 아니다.

**확장:** renderer holdout, multiple rendering engines, style/randomization, public real-image domain
shift, adversarial artifact controls를 W7에 추가한다.

#### A7. simulation sandbox의 자기확증

**공격:** 모델 안에서 가상 knock-out을 하고 모델이 예상한 결과를 얻는 것은 생물 발견이 아니다.

**판정:** 단일 formalism intervention 결과만으로 biological claim 금지.

**확장:** analytic oracle, independent discretization, ablation, blind held-out intervention,
public-data phenotype를 계층적으로 결합한다. sandbox의 1차 산출물은 인과 구분능력과
식별 가능성이다.

#### A8. 특정 세포군의 문법을 보편 세포 문법으로 부를 위험

**공격:** actomyosin 중심 mammalian grammar는 RBC, neuron, ciliated cell, plant/fungal cell에서
불충분하다.

**판정:** 현재 component set을 universal이라고 부르지 않는다.

**확장:** Ring A/B/C와 `missing component law` verdict를 도입했다. 적용 실패가 새 component
프로젝트를 연다.

#### A9. accepted reference가 아직 완성되지 않음

**공격:** force-accepted, aliased, non-admissible trace 위에서 AI를 만들면 solver artifact를
정교하게 학습한다.

**판정:** W1은 전체 플랫폼의 P0. 다만 다른 레인을 정지시키지는 않는다.

**확장:** W4는 analytic oracle, W7은 public/synthetic toy geometry, W9는 literature graph에서
병렬 개발한다. native biological claim만 W1에 잠긴다.

#### A10. 계산 비용이 확장을 막을 위험

**공격:** archetype과 intervention을 늘리면 native budget이 선형 이상으로 커진다.

**판정:** 세포종 수를 줄이는 방식으로 해결하지 않는다.

**확장:** hierarchical sharing, adaptive design, tensor parameter response, learned preconditioner,
multi-fidelity evidence를 결합한다. 계산 예산 자체가 algorithm paper의 평가축이다.

#### A11. 문헌·omics·영상의 cell-state 불일치

**공격:** 서로 다른 세포주, passage, cycle, substrate, 시간이 섞이면 prior가 거짓 평균이 된다.

**판정:** cell type label만으로 evidence를 합치지 않는다.

**확장:** CellStateManifest의 state-schema applicability와 별도 ExperimentContextManifest의
context applicability, context-conditioned mixture prior, proxy distance를 지식그래프에 추가한다.

#### A12. QM/MD 레인의 과장

**공격:** DFT가 whole-cell density나 migration parameter를 주지 않는다.

**판정:** molecular free-energy/energetics로 범위를 명시한다.

**확장:** QM→MD/CG→connector prior의 scale-bridging 연구를 열되, cell-level posterior와 분리한다.

#### A13. 프로젝트가 서로 무관한 AI 도구 모음이 될 위험

**공격:** segmentation, GNN, TT, renderer가 각각 별도 demo로 끝날 수 있다.

**판정:** 모든 도구는 `CellStatePosterior` 또는 `ObservationOperator`를 입력/출력해야 한다.

**확장:** 공통 artifact schema와 intervention sandbox를 통합 시험대로 둔다.

#### A14. 논문 수를 늘리기 위한 salami slicing

**공격:** 하나의 pipeline을 작은 논문으로 과분할할 수 있다.

**판정:** 독립 falsifier·dataset·method contribution이 없는 항목은 논문이 아니라 module이다.

**확장:** 실패를 공유 benchmark/resource paper에 합치고, 생물 질문은 archetype 간 비교로
확장한다.

#### A15. 문서 재창설이 다시 문서 드리프트를 만들 위험

**공격:** 새 방향마다 새 master document를 만들면 이전 audit 실패가 반복된다.

**판정:** 이 문서는 승인 후 `ROADMAP.md`로 흡수되거나 그 파일을 대체해야 한다. 병렬 canonical
plan으로 영구 유지하지 않는다.

**확장:** machine-readable capability graph와 generated status view로 전환한다.

### 10.3 감사 총평

이 계획에서 제거해야 할 것은 세 가지뿐이다.

1. literal quantum position dynamics
2. universal single-Gaussian/full-density representation
3. surrogate 또는 synthetic data가 스스로 물리·생물 gate를 닫는 구조

그 외의 실패는 범위 축소가 아니라 다음 확장으로 라우팅할 수 있다.

```text
rank failure          → ensemble / sparse / latent-coordinate research
image nonidentifiability → temporal + force + intervention modalities
cell transfer failure → new component law / context-conditioned prior
native cost failure   → solver / preconditioner / adaptive design
simulation-real gap   → model-inadequacy discovery
QM scale failure      → MD/CG/experimental prior boundary
```

---

## 11. 의사결정 규칙

새 아이디어는 다음 다섯 질문으로 배치한다.

1. **Physics:** reference dynamics를 바꾸는가, 관측/축약/추론만 하는가?
2. **Representation:** point, ensemble, tensor, graph, field, image 중 무엇을 소유하는가?
3. **Evidence:** native, oracle, reduced, surrogate, public data 중 무엇이 결과를 지지하는가?
4. **Applicability:** 어떤 cell-state schema와 어떤 별도 external context에서만 유효한가?
5. **Failure expansion:** 실패 시 어느 새 component/method/data 레인을 여는가?

판정:

- reference physics 변경 → PI contract + A/B + native validation
- observer/telemetry → additive, default-on recording을 우선
- reduced/tensor → approximation error와 fallback 필수
- neural → trust region, calibration, provenance 필수
- image → observation model과 confounder model 필수
- cell expansion → missing-law verdict 허용

---

## 12. 첫 실행 묶음

이 계획이 승인되면 첫 묶음은 다음 순서다.

1. `CLAUDE.md`를 새 정체성·표현 경계·증거 계층으로 개정
2. `ROADMAP.md`를 특정 세포주 stage가 아니라 capability waves로 개정
3. 기존 모든 파일을 `active / evidence / historical / archive / generated`로 분류
4. `StateFieldSpec`·`StateOwnerSchema`·`CellStateSnapshot`·`ExperimentContextManifest` 계약 작성
5. artifact에 `representation`, `evidence_source`, `cell_state_manifest_hash`,
   `experiment_context_manifest_hash`, `projection_hash`, `trust_region` 추가
6. accepted trace 한 레인과 owner-semantic lossless tensor telemetry 연결
7. analytic oracle 위 state-density/transition-operator TT positive/negative harness 구축
8. toy/native state → synthetic microscopy → image-to-state posterior 폐루프 구축
9. Ring A archetype manifest 초안과 missing-law matrix 작성
10. sandbox S0/S1/S11/S12부터 실행
11. 감사 A1–A15를 자동/반자동 체크리스트로 변환

첫 목표는 “가상 세포 전체 완성”이 아니다.

> **하나의 accepted physical trace가 확률 artifact가 되고, 텐서로 축약되며, 합성영상으로 렌더링되고,
> 다시 calibrated state posterior로 복원되는 첫 폐루프**가 초기 통합 성공 조건이다.

이 폐루프가 닫히면 세포 archetype, intervention, modality, tissue scale을 같은 계약 아래 계속
확장할 수 있다.

---

## 13. 문헌상 위치

이 계획은 다음 이미 존재하는 연구축을 결합하지만, 어느 하나와도 동일하지 않다.

- AI virtual cell: multimodal·multiscale cell-state prediction
  (<https://doi.org/10.1016/j.cell.2024.11.015>)
- tensor network를 이용한 reaction–diffusion rare-event 계산
  (<https://doi.org/10.1103/PhysRevX.13.041006>)
- chemical master equation의 TT/QTT 직접해
  (<https://doi.org/10.1371/journal.pcbi.1003359>)
- 가변 입자수 reaction–diffusion의 classical Fock-space 확률 형식
  (<https://doi.org/10.1007/s11005-022-01539-w>)
- Fokker–Planck/parabolic problem의 TT/QTT 시간진화
  (<https://doi.org/10.1137/120864210>)
- high-dimensional Fokker–Planck tensor neural representation
  (<https://arxiv.org/abs/2404.05615>)
- Fourier neural operator
  (<https://openreview.net/forum?id=c8P9NQVtmnO>)
- MeshGraphNets
  (<https://arxiv.org/abs/2010.03409>)
- Cellpose
  (<https://doi.org/10.1038/s41592-020-01018-x>)
- JUMP Cell Painting perturbation-image resource
  (<https://doi.org/10.1038/s41592-024-02241-6>)
- perturbation-conditioned morphology generation
  (<https://doi.org/10.1038/s41467-025-63478-z>)
- stochastic actomyosin–adhesion cell migration
  (<https://arxiv.org/abs/1810.11435>)
- Tabula Sapiens multi-organ cell atlas
  (<https://doi.org/10.1126/science.abl4896>)
- HuBMAP Human Reference Atlas
  (<https://pmc.ncbi.nlm.nih.gov/articles/PMC11978508/>)
- integrated intracellular organization in human iPS cells
  (<https://pmc.ncbi.nlm.nih.gov/articles/PMC9834050/>)
- scGPT single-cell foundation model
  (<https://doi.org/10.1038/s41592-024-02201-0>)

FFN의 목표 노벨티는 개별 기술의 최초 사용이 아니다.

> **명시적 connector-resolved mechanics, source-factored probability, tensor/neural reduced models,
> microscopy observation operators, cell-state-conditioned evidence와 posterior를 하나의 falsifiable
> accepted-trajectory loop로 닫는 것**이 통합 연구 주장이다.

---

## 14. 실행 개시 기록 — 2026-07-29

이 절은 canonical engine·gate·`CLAUDE.md`를 바꾸지 않고 시작한 CPU-only 격리 tranche다.

### 구현된 공통층

- `aleph/virtual_cell/contracts.py`
  - 특정 세포주에 고정되지 않은 single-cell archetype manifest
  - uncertainty source와 representation provenance
  - approximation error를 source/reduced blob에 묶는 typed report
  - generic layer의 `NATIVE_ACCEPTED` 자가 선언 금지
  - sandbox experiment의 positive/negative/adversarial/held-out/falsifier/expansion 계약
- `aleph/virtual_cell/population.py`
  - heterotypic archetype, stable cell instance ID, namespaced interaction graph, shared environment,
    collective missing-law manifest
- `aleph/virtual_cell/step_receipt.py`
  - PI-gated predicate를 해석하지 않는 `unverified-observer` receipt
  - decision/build/run/predicate/state-coverage/state/RNG digest 결합
  - attempted/accepted clock 및 receipt-chain 연속성 검사
- `aleph/virtual_cell/tensor_train.py`, `reduction.py`, `probability.py`
  - dense/SVD/TT-SVD CPU oracle와 measured error
  - 고랭크·비압축·오차초과 시 dense fallback
  - 확률 질량의 비음수성/정규화 보존 실패 시 fallback
  - 선언 projection 실패가 global Frobenius 통과보다 우선
  - scale-safe norm으로 유한 대수 입력의 overflow/NaN 승인 차단
- `aleph/virtual_cell/sidecar.py`
  - content hash와 byte count
  - safe NPY decode 후 dtype/shape/order/axis 검증
- `aleph/virtual_cell/synthetic_microscopy.py`, `observation.py`
  - Gaussian point-PSF + photon/read-noise/saturation/quantization CPU oracle
  - world-y handedness와 pixel-edge origin 명시
  - frame crop 전에 PSF를 정규화하여 boundary ghost 차단
  - state/manifest/optics/camera/seed/channel-map/transform/sidecar를 한 composite hash에 결합

### 적대적 감사에서 실제로 잡아 고친 것

1. enum만으로 `NATIVE_ACCEPTED`를 자칭하던 평행 증거 경로
2. 임의 문자열을 rollback hash로 받아들이던 receipt
3. receipt 누락·역행·clock/state/RNG 불연속 미검출
4. 정상 확률을 SVD한 뒤 음수와 질량 손실이 생겨도 승인
5. global error가 작지만 희귀 local projection이 95% 틀린 압축 승인
6. `1e308` 유한 입력의 overflow가 NaN error를 만들고 비교를 우회
7. 임의 byte blob을 배열 sidecar로 오인
8. sidecar·manifest·optics·camera를 서로 바꿔 끼울 수 있던 cross-object swap
9. frame 밖 PSF 꼬리를 완전한 밝기로 재정규화하던 경계 ghost
10. whitespace-normalized predicate/channel/mapping key 충돌

### 기존 아카이브 rehearsal

`sf_operator_observables.npz`를 새 selector로 읽은 결과:

- eigenvalue vector 두 개: dense
- Fisher `(8, 8)`: dense
- positions `(904, 3)` 두 개: dense
- Jacobian `(2207, 8)` 두 개:
  - relative tolerance `0.001`: dense
  - `0.032`: rank 7 SVD, measured error 약 `0.0228–0.0229`, parameter-count ratio `1.138×`
  - `0.1`: rank 6 SVD, measured error 약 `0.0604–0.0605`, ratio `1.328×`

즉 첫 실제 rehearsal은 강한 저랭크 성공이 아니다. 이것이 fail-closed selector가 내야 할 결과다.
parameter-count ratio는 실제 직렬화 byte·metadata·decode 비용을 포함하지 않으므로 저장 성능 주장으로
승격하지 않는다.

### 검증 상태

- 새 virtual-cell suite: 88 tests PASS
- 기존 ensemble/operator/sensitivity/stationarity CPU 회귀와 합쳐 192 tests PASS
- Ruff lint/format PASS
- import 시 Warp 미초기화 확인
- GPU 사용 없음

### 아직 잠긴 것

- runtime-attested accepted native trace와 canonical acceptance conjunction
- tensor sidecar와 receipt chain의 production adapter
- accepted native image/state pair 기반 neural posterior
- public-data source resolver와 source-DAG cycle 검증
- pixel-integrated PSF, z projection, line/surface emitter
- 자동 axis-order search와 실제 serialized-byte benchmark
- component/connector endpoint를 canonical engine grammar에 대조하는 validator
- quantum amplitude/phase/Born map/Hamiltonian/open-system dynamics

마지막 항목은 현재 작업의 누락이 아니라 물리 경계다. 현 구현은 **고전 확률·선형대수·
quantum-inspired tensor representation**이며 양자역학적 세포 위치 동역학이 아니다.

---

## 14.5 이번 작업의 위치 — **Phase 0이지 최종 목표가 아니다** (PI 2026-07-29)

> **이 문단은 §15·§16을 읽기 전에 반드시 먼저 읽는다.** §15·§16만 보면 이 방향이 *검증 도구
> 모음*처럼 보이는데, 그것은 오독이며 실제로 한 번 발생했다(Lead가 PI에게 그렇게 보고했고
> PI가 정정했다).

**최종 목표는** 양자·통계물리에서 발전한 확률분포·연산자·텐서 방법을 **고전 accepted-step
세포역학과 결합**하여, 하나의 점 궤적이 아니라 **`CellStatePosterior`를 전파하는 Probabilistic
Mechanistic Virtual Cell**을 만드는 것이다. 이번 라운드는 그 방향을 **바로 엔진에 넣기 전에**,
거짓 수렴·거짓 압축·비식별·모델 누락·증거 세탁을 잡을 **CPU oracle과 계약을 먼저 만든 것**이다.

한 줄로: **검증 도구 프로젝트가 아니라, 명시적 세포역학을 확률적 상태 전파·Koopman 분석·텐서
네트워크·데이터 동화·영상 AI와 연결하는 장기 플랫폼이고, 지금 만든 것은 그 플랫폼이 스스로
거짓말하지 않게 하는 첫 바닥층이다.**

### 14.5.1 핵심 질문

> 세포를 하나의 결정론적 점 궤적으로만 계산하지 않고, **위치·형태·결합·파라미터를 확률분포로
> 들고 다니며**, 그 분포를 물리 엔진으로 전파하고, 텐서·연산자·신경망으로 분석·축약·관측할 수
> 있는가?

세포주를 늘리는 것이 아니라 §2의 폐루프를 닫는 것이 목적이다.

### 14.5.2 양자역학적 개념의 정확한 위치

**차용하는 것은 알고리즘이지 주장이 아니다.** 아래는 실제 알고리즘 후보이며, **세포 자체가 양자
상태라는 주장은 아니다.**

| 방향 | 차용 개념 | FFN에서의 용도 |
|---|---|---|
| 확률적 세포 상태 | Gaussian wavepacket, density operator에 가까운 **표현적 직관** | 노드를 점이 아니라 `(μ, Σ)` 또는 ensemble/mixture로 |
| 동역학 연산자 | Liouville, Koopman, Perron–Frobenius | 진동·NESS·metastability·예측 지평 분석 |
| 탄젠트 전파 | variational dynamics, Lyapunov subspace | 파라미터 감도, covariance 전파, 유효 예측 반경 |
| 텐서 네트워크 | MPS/TT/DMRG | 고차원 분포·응답 텐서·텔레메트리의 **조건부** 저랭크 |
| 경로 측도 | Girsanov/score martingale, path integral | KMC 사건의 파라미터 민감도, likelihood ratio |
| 데이터 동화 | ensemble Kalman/particle/variational | brute-force sweep → posterior 갱신 문제로 전환 |
| 희귀사건 | Kramers, WKB, instanton | 결합 탈출률·장벽 전인자의 장기 mechanistic 확장 |
| 양자화학 | DFT → MD/CG energy landscape | 문헌에 없는 분자 파라미터의 **prior** 생성 |

**감사 이후 살아남은 정확한 표현:** *FFN은 고전 mechanochemical runtime 위에 quantum-inspired
확률·연산자·텐서 방법을 얹는 플랫폼이다.*

### 14.5.3 명확히 기각된 것 (재제안 금지)

- 세포 노드의 위치를 **복소 파동함수**로 계산
- 현재 active cell 전체를 **Schrödinger 방정식**으로 전환
- `D_xΦ` · Hessian · Koopman · Fisher · score operator를 **하나의 연산자로 동일시**
- **양자컴퓨터·QAOA**로 파라미터 스윕 가속
- TT bond dimension을 **즉시** 상관길이·얽힘 같은 물리량으로 해석
- `Tr(K)`를 **엔트로피 생산률과 동일시**
- Jensen 항을 **현재 point simulation의 버그**라고 주장

### 14.5.4 살아 있는 개발 방향 8가지

1. **확률적 위치·형태**: `node = x` → empirical ensemble / local `(μ, Σ)` / mixture / latent density.
   힘을 **평균 위치에서 한 번**이 아니라 **분포의 모멘트에서** 계산해야 한다. 확률적 막과
   Helfrich covariance, cortex–membrane tether 분포 평균, FA bond ensemble, reflected SDE로 확장.
2. **한 런에서 감도와 국소 스윕 반경**: 편차 벡터를 탄젠트 연산자 위로 전파. 단 **같은 regime
   내부의 유한 반경**만 유효하고, topology나 regime이 바뀌면 pathwise derivative를 **중단**하고
   ensemble/regime-partition으로 후퇴한다. "한 런으로 전체 파라미터 공간"은 기각.
3. **스윕 → 데이터 동화**: 고정 격자 → sensitivity-guided adaptive → posterior update →
   assimilation cycle. 질문이 "몇 점 돌릴까"에서 "무엇이 식별 가능한가 / 어떤 intervention이
   posterior를 가장 줄이는가"로 바뀐다.
4. **텐서 네트워크는 압축기 이상**: 고차원 `ρ(position, orientation, state, parameter, time)`,
   response/covariance 구조, 다축 텔레메트리. 단 1D→dense, 행렬→SVD, 3축 이상→TT 후보,
   실패→dense/ensemble.
5. **Koopman/Perron–Frobenius 동역학 해석**: accepted native trajectory 위의 transfer spectrum,
   metastable state, relaxation mode, probability current, 시간역전 비대칭.
6. **합성영상–역추론–신경망 폐루프**: accepted state → observation operator → synthetic image →
   calibrated posterior encoder → `CellStatePosterior` → physics propagation.
7. **다세포종·조직 확장**: manifest + missing-law 방식. **cell–cell connector가 0개**라는 구조적
   병목이 확인됐으므로 조직 단계는 파라미터 추가가 아니라 **새 connector law**가 필요하다.
8. **분자 수준 QM/MD 경로**: DFT → atomistic/CG MD → free-energy landscape → Kramers/WKB rate →
   cellular prior. DFT가 곧바로 pN stall force나 세포 수준 밀도를 주지 않는다.

### 14.5.5 완료 / 아직 안 됨

```text
[완료 — 이번 두 라운드]
증거·표현 계약 · 수치 불변성 sandbox · regime analytic oracle
텐서/SVD 선택기 · 합성 optics oracle · inverse identifiability oracle
surrogate rejection harness · archetype/missing-law grammar · 적대적 감사 register

[아직 안 됨]
CellStatePosterior 실제 타입          ← 계획서가 선언한 중심 객체
ObservationOperator 공통 인터페이스    ← 계획서가 선언한 중심 객체
runtime tangent/score adapter · 블록 covariance 전파 · 확률 막·확률 노드
Koopman/transfer operator on native trace · data assimilation driver
neural posterior · FA ensemble dynamics · reflected cytosol boundary
Kramers/WKB molecular-prior pipeline
```

**가장 중요한 결손은 위 두 줄이다.** 감사 register(§16.5)가 A13에서 이미 기계적으로 잡았다 —
계획서가 **중심 객체로 선언한** `CellStatePosterior`와 `ObservationOperator`가 코드에 존재하지
않는다. 감사 도구가 아무리 잘 만들어졌어도 **이 두 인터페이스가 없으면 각 레인이 하나의 확률적
virtual-cell loop로 닫히지 않는다.** 지금 레인들은 좋은 부품이지만 아직 **회로가 아니다.**

---

## 15. Lane A–F 병렬 실행 — 2026-07-29 (진행 중)

### 15.0 먼저 해결한 것 — 이 레이어는 아무도 커밋할 수 없었다

`aleph/virtual_cell/**`와 그 테스트는 `ownership.yaml`에 **unclaimed**였고, 이 저장소에서
unclaimed는 free가 아니라 **모든 세션에게 refused**다. §14의 tranche는 디스크에만 있었고
자기를 선언한 어떤 세션도 커밋할 수 없는 상태였다.

`virtual-cell` 레인을 선언해 해결했다 (`a51a3c1c`). 확인된 동작:

- positive control — `aleph/virtual_cell/**`, `aleph/tests/ac/virtual_cell/**`,
  `aleph/scripts/vc_*.py`, 이 계획서 4경로 모두 `OK ... owned by 'virtual-cell'`
- negative control — `aleph/components/incumbent/driver.py`는 `REFUSED ... belongs to 'cortex'`

**이 레인은 `aleph/ac/` 아래를 하나도 claim하지 않는다.** 이 레이어는 구성상
unverified observer이므로, 엔진 안에서 필요한 변경은 소유 레인에 올리는 **제안**이지
이 레인이 직접 하는 편집이 아니다. 커밋은 전부 pathspec으로 했고, 동시 작업 중인 `cortex`
세션의 dirty file 4개는 커밋 전후로 손대지 않았음을 확인했다.

### 15.1 측정된 baseline (2026-07-29, CPU, dev machine)

계획서 §14가 적은 "192 tests"는 이 시점에 재현되지 않는다. 실제로 실행해 센 값은 다음과 같다.

| 대상 | 측정값 |
|---|---|
| `tests/ac/virtual_cell/` | **88 passed** |
| 기존 CPU 회귀 (`test_ensemble`, `test_observe_operator`, `test_observe_sensitivity_energy`, `test_observe_stationarity`, `test_ensemble_stall_oracle`, `test_gate_b_cortex_motor_stationarity_wiring`) | **122 passed** |
| 합계 | **210 passed** |

192가 아니라 210인 이유는 회귀 대상 목록이 그 사이에 늘었기 때문이며, 계획서 쪽 숫자를
그대로 인용하지 않고 다시 잰 것이 옳은 처리다. **문서의 숫자보다 실행한 숫자가 이긴다.**

### 15.2 병렬 레인 배치 — 파일 단위로 서로소

동시 세션 충돌과 레인 간 충돌을 같은 방식으로 막았다: 각 레인은 **새 모듈 파일과 새 테스트
파일만** 소유하고, 공유되는 `virtual_cell/__init__.py`는 아무 레인도 건드리지 않으며 Lead가
마지막에 export를 통합한다.

| Lane | 소유 모듈 | 주제 |
|---|---|---|
| A | `admissibility.py` | accepted-step admissibility sandbox + PI decision card |
| B | `telemetry.py` | receipt → tensor sidecar adapter, axis-order search, 실제 byte/decode 측정 |
| C | `regime.py`, `sweep.py` | regime 분류, cross-regime FD 거부, ESS, Jensen=moment-closure 비용 |
| D | `optics.py`, `source_registry.py` | pixel-integrated PSF, line/surface emitter, z projection, source-DAG |
| E | `archetypes.py`, `grammar.py` | Ring A/B/C archetype, endpoint/owner validator |
| F | `surrogate.py` | trust region / OOD / coverage / held-out intervention — **학습은 잠금** |

결과는 착지하는 대로 이 절에 **측정값과 함께** 기록한다. 아직 실행되지 않은 것은 이 표에
숫자로 적지 않는다.

### 15.3 Lane A 착지 (`2a18859a`) — 결함의 더 날카로운 형태

핸드오프가 지목한 결함은 "`balance_ok`는 adjoint force closure이지 mechanical convergence가
아니므로, inner solve가 rollback된 뒤에도 outer가 force-accept로 commit될 수 있다"였다.
반례를 쓰기 위해 코드를 읽는 과정에서 **더 정확한 진술**이 나왔고, Lead가 `7daedb3d`에서
독립 확인했다.

| 경로 | 무엇으로 commit하는가 | 나머지 절반 |
|---|---|---|
| `ac/engine/transaction.py` | `ledger.balance_ok_d` (adjoint closure) — `solve()`의 반환을 **아예 버린다** | convergence를 한 번도 참조하지 않음 |
| `ac/cell/driver.py:1207` | `outer_accepted = inner_converged` | balance gate를 **조립조차 하지 않음** |

즉 **두 경로 중 어느 쪽도 conjunction을 구현하지 않고, 서로 다른 절반씩을 구현한다.**
`transaction.py`는 그 하나의 flag로 `finalize_candidate`와 `clock.advance`를 **둘 다** 구동한다.
그리고 driver는 `:1293`에서 full conjunction을 이미 계산하지만 그것은 transaction predicate가
아니라 **사후 `stable` 보고 필드**다. `ac_gate_b_cortex_motor_native.py`는 별개로 모든 스텝을
force-accept하며 자기 소스에 그렇게 적혀 있다.

측정값 (CPU, predicate의 성질이지 물리가 아님):

- **balanced-but-stretched**: adjoint 상대잔차 `8.07e-17`로 통과하는 동안 inner 잔차는
  tolerance의 `2.63e+06`배, iteration budget `400/400` 소진. adjoint-only predicate는 **ACCEPT**,
  conjunctive는 `inner_converged`에서 reject.
- **per-connector vs global**: 두 커넥터가 `+12`/`−12` pN을 누출하면 합이 정확히 0이라
  **global 형태는 ACCEPT**한다. global closure는 per-connector보다 **엄격히 약하다**.
- **nonfinite**: `state_finite`를 **먼저** 평가하고 나머지 clause는 `NOT_EVALUATED`로 남긴다.
  `nan > tol`과 `nan <= tol`이 둘 다 False이므로 순서를 잘못 잡으면 조용히 승인된다.
- **rollback**: `tobytes()` 위의 SHA-256. clock/topology/RNG/state 4개 음성대조군이 각각 검출되고
  어느 채널이 샜는지 이름이 나온다. float가 아니라 byte인 이유는 `+0.0 == -0.0`이 True이고
  `nan == nan`이 False이기 때문.
- **tolerance는 기본값 없는 필수 키워드 인자**이며 `inspect.signature`로 강제된다. 값을 고르는
  것은 라이브러리 기본값이 아니라 **PI 결정**이다.

**BLOCKED (추정하지 않음):** production reject fraction. A5000에서 transaction 계측이 필요하다.
카드는 이를 명시하고, 사용 가능한 두 개의 100% 수치를 일반화하기를 **거부한다** — 하나는
구성상 100%(스크립트가 predicate를 평가하지 않음)이고 다른 하나는 자기 상태가
`INCOMPLETE_OR_FAIL`인 단일 아카이브 런이다.

**PI-DECISION:** conjunction 채택(P1), per-connector closure(P2), P1이 live-lock을 피하려면 필요한
dt controller(P3), reject fraction을 *얻기 위한* shadow-mode 순서(P4), 그리고 모든 tolerance 값.
카드 → `aleph/docs/v2_audit/PI_DECISION_CARD_ADMISSIBILITY_2026-07-29.md`.
**서명 전까지 아무것도 바뀌지 않으며 엔진은 한 줄도 건드리지 않았다.**

검증: 37 tests, 레인 디렉터리 전체 125 passed / 0 failed, ruff clean, import 후
`'warp' not in sys.modules` 확인.

### 15.4 Lane F 착지 (`68367d93`) — 거부 하네스를, 진리가 닫힌 형태인 곳에서 검증

학습은 잠겨 있다. 학습할 accepted native trace가 없기 때문이다. 착지한 것은 **그런 모델을
거부할 하네스**이며, 지금은 진리를 closed form으로 아는 analytic target 위에서 발화시켰다.
torch·jax·warp 없음, 세포 물리 없음. 유일하게 fit되는 객체는 numpy 최소제곱 다항식이고,
그것은 게이트가 잡을 대상이 필요해서 존재한다.

**권위 금지는 문서가 아니라 구조다.** `SurrogatePrediction`에는 `evidence_source` 필드가
없다 — `__slots__` 위의 property이므로
`dataclasses.replace(..., evidence_source=NATIVE_ACCEPTED)`가 raise한다(Lead가 정적 조회로
독립 확인). 우회 경로는 별도로 막았다: gate ledger가 `source_artifact_ids`를 **전이적으로**
따라가고, cycle을 검사하며, 미등록 조상은 fail-closed이고, 거부할 때 세탁 경로를 이름으로 찍는다.

측정값 (seed `20260729`):

| 검사 | 측정 | 대조 |
|---|---|---|
| calibration, n=20,000, nominal 95% | **0.94600** (z=−2.60) | 데이터 보기 전에 선언한 3σ 이항 밴드 통과 |
| 같은 모델이 σ를 0.3배로 보고 | **0.43460** (z=−334.4) | closed form 기대 `0.44346`와 **0.009**까지 일치 |
| trust region 밖 질의 | 2.543 box-width 밖에서 **거부** | 거부 안 했다면 값 `21198.28`, 진리 `sin(12)=−0.53657` — 유계 1인 대상에 오차 `2.12e4` |
| bounding box가 못 보는 OOD | annulus 원점: box excess **정확히 0.0** | k-NN 거리 `0.40182` vs covering radius `0.02509` → 거부 |
| held-out intervention | interpolation RMSE `1.6659e-3` vs held-out **`2.4051`** | 비율 **1443.8** → `mechanism_falsified`. 그 점은 trust-region box **안**에 있고 manifold 검사만 잡는다 |

**게이트가 이 레인 자신의 첫 구현 결함을 잡았다.** predictive sigma로 residual sigma를
보고하면 coefficient variance가 빠져 체계적으로 과신하게 된다 — 측정 coverage `0.9314`
(z=−5.12), 모든 질의를 거부했다. 수정은 **게이트를 넓히는 것이 아니라 모델 클래스를 바꾸는
것**이었다(`σ̂·√(1+leverage)`). 이 구분이 정확히 헌장의 hard rule이 말하는 것이고, 결과에
맞춰 threshold를 고르지 않았다는 실증이다.

전부 **oracle-only**. **BLOCKED**: simulation data에 fit한 surrogate, 세포 observable에 대한
예측 주장, surrogate artifact의 승격 — accepted native trace 전까지 전부. **PI-decision 없음**
(새 gate contract 없음, 선언된 3σ 밴드 외 magic number 없음).

검증: 27 tests. ⚠️ 이 시점 `ruff check`는 Lane C의 **작업 중** `regime.py`에서 F821 2건을
보고한다 — Lane F 파일이 아니며 착지 시 재검사한다.

### 15.5 Lane E 착지 — canonical census에 **cell–cell junction connector가 하나도 없다**

STATE.md에서 베낀 목록이 아니라 **코드에서 읽은** 실제 grammar에 대조했다. Lead가 독립 확인:

- **components 14 / connectors 36** — STATE.md와 정확히 일치.
- **junction 계열 커넥터 0개** — cadherin·desmosome·tight·gap 어느 것도 없다.
- `fa_actin_anchor` owner = `ECMWorld`, endpoints = `focal_adhesion` / `sf_arc`.
- census를 다 읽어도 `warp`·`aleph.ac` 둘 다 `sys.modules`에 없다.

**census를 읽는 것 자체가 함정이었다.** `import aleph.engine.contracts`는 패키지
`__init__`를 통해 **Warp를 초기화한다.** 그래서 census는 file path로 로드하고(`contracts.py`는
`dataclasses`/`enum`만 import하는 leaf), ownership은 `dispatch.py`에서 AST로 추출한다.
subprocess 테스트가 이를 강제한다. **fallback 목록은 없다** — 소스가 없거나 파싱 실패하면
`CensusUnavailableError`.

**§2.4의 구조적 결론:** 36개 edge가 전부 intra-cell이거나 cell↔shared-environment다.
따라서 **이 레이어의 모든 multicell coupling은 missing-law 프로젝트**이고, 이 갭 하나가
endothelial·epithelial·stem-colony·cardiomyocyte archetype을 **동시에** 막는다. §W11이 기다리는
것이 바로 이 항목이다.

**ownership 규칙은 추정이 아니라 측정으로 보고한다.** `ConnectorContract`에는 owner 필드가
**아예 없다.** ownership은 8개 facade에 걸친 exact-once dispatch claim이다. "owner가 endpoint
하나를 소유한다"는 **36 중 35**에서 성립하고, 유일한 예외 `fa_actin_anchor`는 결함이 아니다 —
`focal_adhesion`이 `owns_geometry=False`라 facade가 없고 edge가 mechanical group을 공유하므로
series joint가 한 번만 소유되는 하나의 런타임 객체다. 측정된 규칙을 인코딩하고 예외는
`EXPANSION`으로 보고했다(누락시키지 않았다).

**Ring B/C는 REFUSE하며, 그것이 산출물이다.** 각 verdict는 엔진이 스스로 선언한 representation
문자열을 인용하고, 두 테스트가 그 인용이 verbatim임을 강제한다.

| Archetype | 대표적 missing law | 왜 파라미터 문제가 아닌가 |
|---|---|---|
| erythrocyte | `membrane_in_plane_shear_elasticity` | 막이 *"fluid Helfrich surface FEM"* — fluid surface는 구성상 in-plane shear modulus가 **0**이다 |
| neuron-growth-cone | `processive_motor_cargo_transport` | 유일한 MT motor edge가 cortical dynein **capture** 사이트이고, cargo state를 소유하는 것이 없다 |
| ciliated-cell | `axonemal_dynein_sliding` | 한 개 빼고 모든 motor edge가 actin 위 NMII — track이 MT doublet인 motor law가 없다 |
| cardiomyocyte | `sarcomeric_registration` | sf_arc가 *"active rod/cable graph"* — register도 Z-disc도 없어 length–tension이 유도 불가 |

`MissingLawVerdict`는 빈 law 목록을 구조적으로 거부한다 — **generic failure는 verdict가 아니다.**

**정직한 부정 결과 (테스트로 고정해 썩지 않게 함):** 14-component census는 **archetype-coarse**다.
endothelial과 fibroblast가 **동일한** 전체 component set을 가져가고, Ring A 6개 중 3개가 중복이거나
component 하나 차이다. 이 archetype들을 가르는 것은 component census가 아니라 `missing_laws`에 있다.

**BLOCKED:** archetype prior 전부 `PI_GAP` 슬롯이며 구조적으로 숫자를 담을 수 없다. 읽지 않은
문헌은 하나도 인용하지 않았고, 유일한 인용은 저장소 자신의 gated HeLa turgor proxy다.
테스트가 어떤 archetype도 MCF7·MDA-MB-231·HeLa를 이름에 담지 않음을 강제한다 —
**프로그램은 다시 좁혀지지 않았다.**

**PI-decision 3건:** (a) Ring A 패널 구성과 각 archetype의 생리적 niche stiffness,
(b) **cell–cell junction connector를 canonical census에 추가할지** — 현재 가장 큰 구조적 갭,
(c) `fa_actin_anchor` owner 예외를 규칙으로 비준할지 census를 재구성할지.

검증: 55 tests, 레인 디렉터리 251 passed / 0 failed.

### 15.6 Lane B 착지 — **global 2.29%가 local 99.98%를 가리고 있었다** (실제 아카이브에서)

§14의 rehearsal은 "약한 저랭크"라고만 적었다. Lane B가 그것을 **byte로** 다시 재고, 그 과정에서
훨씬 중요한 것을 찾았다. Lead가 커밋된 아티팩트에서 직접 재현했다:

```
jacobian_fd  (2207, 8),  tol 0.032,  rank-7 SVD
  global relative Frobenius error : 0.022928   (2.29%)
  entry [521, 5]  original        : 0.972998
  entry [521, 5]  reconstructed   : 1.502616e-04
  local relative error            : 99.9846%       (max|J| = 1.0007)
```

선언된 `VALUE` projection이 이 압축을 **거부**한다. 즉 §14가 기록한 `1.138×` 축약은
**모든 global norm 검사를 통과하면서 full-magnitude Jacobian 성분 하나를 파괴했을 것이다.**
per-projection gate가 방어적 상용구가 아니라는 **실데이터 확증**이다.
→ 라우팅: `projection-failure → ensemble / regime partition / latent coordinate`.

**byte와 decode를 parameter count와 분리해 실측했다.** metadata member까지 포함한 실제
uncompressed NPZ 직렬화 기준:

| tol | rank | param-count ratio | **byte ratio** | decode dense → reduced (median) |
|---|---|---|---|---|
| 0.032 | 7 | 1.1382 | **1.1330** | 0.090 → **0.180 ms** |
| 0.100 | 6 | 1.3279 | **1.3194** | 0.090 → **0.172 ms** |

byte ratio가 parameter ratio보다 **일관되게 낮다** — 두 수는 증명 가능하게 다른 수이고 서로를
대체하지 못한다. 그리고 유일하게 압축되는 member는 byte의 13%를 사기 위해 **decode wall-clock
약 2배**를 낸다. metadata 때문에 오히려 커지는 경우도 실재하며 4×4 rank-1 행렬 단위 테스트로
고정했다. **§14의 약한-결과 판정은 유지되며, 이제 parameter가 아니라 byte로 측정됐다.**

축은 위치가 아니라 **의미**다. time 축은 accepted subsequence에서 **파생**되고 caller가 공급할 수
없다. attempted step은 개수와 인덱스를 기록하며 버린다. non-monotonic·gapped accepted clock은
재인덱싱이 아니라 **거부**한다. 한 time index에서 entity 순서가 바뀌면 raise한다 — 조용히
재바인딩되는 위치 인덱스가 축 선언이 막으려는 바로 그 결함이기 때문이다.

axis-order search는 **항상 `unsearched_fraction`을 보고한다.** 이 2-D 사례에서는 `0.0`이고 두
순서가 동률이다 — 행렬에서 axis order는 아무것도 사지 못한다는, 가정이 아니라 측정된 답이다.

**BLOCKED (즉흥 처리하지 않음):** `RepresentationKind`에 matrix-SVD 멤버가 없다. 2축 SVD를
`TENSOR_TRAIN`으로 라벨링하면 **거짓 representation 주장**이 되므로 `compression_envelope()`는
accepted SVD를 **거부한다**. `MATRIX_SVD` 추가는 다른 소유자의 `contracts.py` 변경이자
**PI decision**이다. 거부를 테스트 2개로 고정했다.

검증: 44 tests.

### 15.7 Lane D 착지 — point sampling은 광자의 절반 이상을 잃는다

`synthetic_microscopy.py`는 **수정하지 않았다.** `OpticsConfig`가 `MicroscopyConfig`를 대체하지
않고 감싸자 공유 코드 필요가 사라졌고, 기존 렌더러는 새 구성공간의 **정확히 한 모서리**로
`2.5e-16` 상대오차 재현된다 — world-y handedness, pixel-edge origin, crop 전 정규화가 재유도가
아니라 **구성상** 보존된다.

**pixel-integrated PSF**는 quadrature가 아니라 closed form(pixel edge 위 `erf` 차의 곱)이다.
σ = 0.3 px에서 point-sampled 렌더러의 총 광자수는 **emitter가 픽셀 안 어디에 있느냐만으로**
자기 값의 절반 이상 흔들리고 최솟값이 선언된 1.0의 약 `0.44`까지 내려간다. 같은 sub-pixel
스캔에서 integrated 형태는 `1e-16`으로 보존된다. Lead가 보존 결과와 최솟값을 독립 재현했다
(Lead의 재구현은 정규화 규약이 달라 최댓값은 일치시키지 않았으므로 swing 수치는 인용하지 않는다).
부수 효과 하나가 기록할 만하다: integrated kernel의 회수된 2차 모멘트는 `2(σ²+1/12)` px²이고
point sampling은 `2σ²`다 — **point sampling은 카메라가 실제로 기록하는 폭을 과소보고한다.**

line/surface emitter는 **단위 길이당·단위 면적당**이고 두 convolved field 모두 closed form이다
(선분은 lateral Gaussian × axial `erf` 쌍, 삼각형은 divergence theorem으로 3개 signed wedge —
각각 `arctan` − **Owen T**). 쓰기 전에 수치적분으로 대조했다. **per-primitive 정규화를 잡는**
subdivision 불변성이 `1e-16`에서 성립한다: 선분을 2/3/5/17로 쪼개거나 삼각형을 64개로 정제해도
이미지가 round-off만큼만 바뀐다 — wedge가 edge 반전에 정확히 반대칭이라 내부 edge가 상쇄되기 때문.

z projection은 **widefield 전용이며 그렇게 적혀 있다.** defocus에서 frame total이 `8.9e-16`로
보존되고 peak는 `15.7×` 떨어진다. **confocal과 light-sheet는 근사하지 않고 typed capability
verdict로 거부한다** — 둘 다 axial **gate**이므로 detection-side kernel을 넓혀서는 얻을 수 없다는
명시된 근거로.

**provenance DAG — A6(자기채점) 케이스가 핵심이다.** synthetic renderer 출력에서 **두 홉 아래**의
morphology statistic이 blocking source 이름과 전체 경로와 함께 배제된다. 배제는 **임의 깊이의
ancestry에 대해 닫혀 있고**, origin은 **whitelist**로만 허용되므로 새 origin 종류는 fail-closed다.
cycle은 freeze에서 거부하며 경로로 보고하고, 테스트가 연속 쌍이 전부 실제 child→parent edge임을
확인해 그것이 집합이 아니라 **경로**임을 강제한다. `CALIBRATION_ONLY`는 **따로, 전이적으로**
차단된다 — 증거로는 인용 가능하지만 held-out으로는 절대 쓸 수 없다.

**아무것도 다운로드하지 않았고 아무것도 통합되지 않았다.** JUMP·Allen hiPSC·IDR live-cell은
`blocked / not-yet-integrated`로 등록됐고, dataset contract는 14개 필수 필드 중 하나라도 없으면
`INTEGRATED` 주장을 거부한다. 유일한 worked fixture는 `SourceStatus.FIXTURE`를 달고 있어
자신과 파생물이 전부 배제된다 — **fixture는 증거가 아니다.**

검증: 72 tests, 기존 8개 파일 무손상(microscopy 16/16).

### 15.8 레인이 레인을 잡았다 — Lane F 테스트 결함 (`eb28ea9a`)

Lane B가 전체 `tests/ac` 실행에서 Lane F의
`test_importing_the_harness_does_not_initialise_warp_or_a_learning_framework` 실패를 보고했다.
원인은 이 테스트가 **자기 프로세스의 `sys.modules`**에 대해 단언한다는 것이다. 그것은 이름이
말하는 주장과 다르다 — "이 pytest 세션에서 무언가가 Warp를 import한 적이 있는가"를 묻는 셈이고,
스위트에는 정당하게 import하는 테스트가 있다. 그래서 **레인 단독 실행에서는 통과하고 전체
실행에서는 실패했다.**

Lane A와 Lane E가 이미 쓰던 subprocess probe로 교체했다. 새 인터프리터가 이름이 말하는 것을 잰다.
**이것이 병렬 레인의 값이다 — 레인을 쓴 쪽이 아니라 다른 레인이 잡았다.**

**그리고 같은 결함이 Lane C에도 두 군데 있었다** (`06f9f949`에서 수정). 통합 시점의 전체 스위트
실행이 잡았다. 즉 이 결함은 **한 레인의 실수가 아니라 이 계약을 테스트하는 자연스러운-그러나-틀린
방법**이며, 따라서 앞으로 이 레이어에 추가되는 모든 CPU-only 주장은 subprocess로 재야 한다.

### 15.9 Lane C 착지 — **Lyapunov 지수는 chaos를 NESS와 분리하지 못한다** (측정됨)

11개 analytic control에 대한 confusion matrix의 **off-diagonal이 0**이고, UNDECIDED가 서로 다른
세 음성대조군에서 각자의 이유와 함께 반환된다(white noise / 미완화 transient / 선언 최소길이 미달).

이 레인의 값은 **측정해 보고 기각한 세 통계**에 있다. 각각 모듈에 문서화했고 조용히 버리지 않았다.

| 후보 통계 | 왜 기각됐는가 (측정) |
|---|---|
| raw Lyapunov 지수 | best-R² fit-window 탐색에서 **OU가 0.68/s, logistic map이 0.64/s** — 확률 대조군이 카오스 대조군보다 **위**. 추정기를 고쳐도(fit window를 탐색이 아니라 측정된 `tau_int`에서 유도) white noise가 4.3/s로 Lorenz(0.835/s)보다 **여전히 높다**. → **λ는 보고하되 절대 gate하지 않으며**, 리포트에 그렇게 적는 필드가 있다 |
| raw recurrence determinism | Theiler 배제 없이는 **매끄러움**을 재지 결정성을 재지 않는다 — OU 0.78 vs Lorenz 0.99 |
| determinism vs phase surrogate z | logistic은 정확히 잡지만(z=14.9) **Lorenz는 z=1.08** — 매끄러운 flow의 phase surrogate는 똑같이 매끄럽기 때문 |

**작동하는 통계: surrogate-referenced nonlinear predictability.** Lorenz z=5.5, logistic z=45.7,
비카오스 대조군 전부 ≤ 0.588, 선언 임계 4.0. **누락된 통계를 이름만 부른 것이 아니라 구현했다.**

대조군이 드러낸 결함 4건 추가 수정: first-zero-crossing embedding delay가 Lorenz를 lag 177에
놓았고(1/e 규칙은 16) 그 나쁜 embedding에서는 예측이 **자기 null보다 나빴다**; 고정-bin spectral
baseline은 두 역할을 동시에 못 한다(25 bin에서 quasicycle의 hump가 자기 baseline을 삼키고,
101 bin에서 OU roll-off가 prominence 8.13의 허위 line을 낳는다); Lorenz가 수치 꼬리에서
**"peak" 191개**를 만들었고 전부 최댓값의 2% 미만인데 그중 둘로 quasiperiodic 판정을 받았다;
broad Welch hump는 비율 1.065로 갈라지는데 그건 작은 유리수가 아니므로 quasiperiodic 판정은 이제
**resolution-limited line 2개**를 요구하며 sharp/broad 임계 사이에 **의무 UNDECIDED 간격**이 있다.
실제 크래시도 하나 고쳤다: `np.argsort`는 stable하지 않아 admissible neighbour가 없는 행이 배열
밖을 인덱싱했다 — 그런 행은 이제 clip이 아니라 **버린다**. clip된 인덱스는 **날조된 예측**이기 때문.

sweep 측정값:

- **cross-regime 거부 작동**: OU(μ=−0.5, `stochastic-ness`) vs Hopf(μ=1.0, `limit-cycle`) →
  `REFUSED_CROSS_REGIME`, `value is None`. 같은 regime 쌍은 구성된 참값 0.4에 대해
  **0.37936 ± 0.079879**.
- **ESS는 일괄 팽창이 아니다**: AR(1) φ=0.98에서 N=8000, **n_eff = 107.1**, 구간 **8.644× 확대**.
  white noise에서는 확대가 **1.00**으로 붕괴한다.
- **phase decomposition**: 120주기 limit cycle에서 time average `0.001897` vs cycle amplitude
  `1.859` — **스칼라 평균은 자기가 지운 것의 약 0.1%다.**
- **Jensen 항은 point model의 버그가 아니라 분포 표현 채택의 moment-closure 비용**으로 다뤘고
  closed form에 대조했다: `f=x²`는 gap이 모집단 분산과 `1e-12`까지 일치(정확한 항등식),
  Gaussian 위 `f=exp`는 `e^μ(e^{σ²/2}−1)`에 `1.85e-3`까지 일치. 위반은
  ensemble/latent-coordinate로 라우팅되며 **절대 dynamics로 가지 않는다.**

**Provenance:** `integrated_autocorrelation_time`이 엔진의 `observe/stationarity.py`와
**3.38e-16** 상대오차로 일치한다(Lead 독립 확인). 엔진 모듈은 import할 수 없어서(패키지
`__init__`가 Warp를 초기화) file path로 로드하고 **수치적 동일성**을 단언한다 — 호출을 공유하는
것보다 강한 검사다.

**PI-decision:** `RegimeEvidenceContract`는 **어떤 필드에도 기본값이 없다**(테스트로 강제).
테스트 파일의 값은 *테스트의* 선언이지 프로젝트 기본값이 아니며 production contract는 PI 서명이
필요하다. per-observable `amplitude_scale`도 PI다 — fixed-point vs NESS가 **선언된** 스케일에
대해 결정되므로, 데이터를 본 뒤 그 스케일을 고르는 것이 요동치는 궤적을 fixed point로 만드는 방법이다.

검증: 33 tests.

### 15.10 통합 (`06f9f949`) — 재수출하면 **안 되는** 이름 하나

18개 모듈 **195개 public name**을 통합했다. 모든 모듈의 `__all__`이 재수출되고, 재수출된 모든
이름이 resolve되며, 패키지 import 후에도 Warp는 초기화되지 않는다.

**충돌은 하나뿐이고, 승자를 고르는 방식으로 풀지 않았다.** `sweep`과 `surrogate`가 각각
`ExpansionRoute`를 정의하는데 **서로 다른 enum**이고 멤버가 부분적으로만 겹친다
(`adaptive-design`·`model-inadequacy`는 양쪽, `regime-partition`·`native-evaluation`은 한쪽씩).
bare `ExpansionRoute`를 export하면 **마지막에 import된 쪽으로 조용히 resolve**되는데, 그것이
정확히 이 패키지가 잡으려는 결함 부류다. 그래서 `SweepExpansionRoute`/`SurrogateExpansionRoute`로
내보내고 bare 이름은 **의도적으로 없다.** 둘을 병합하는 것은 cross-lane 설계 변경이므로
**통합 단계가 발명하지 않는다.**

한 번도 export되지 않았던 이름 4개를 추가했다. 그중 `validate_step_receipt_chain`은
**receipt-chain 연속성 검사**로 load-bearing인데 패키지 루트에서 도달 불가였다.

**전체 `tests/ac`: 1705 passed, 87 skipped, 0 failed.**

### 15.11 boot budget이 빨갛게 되었다 — 내 변경 때문이 아니다

통합 커밋이 `check_boot.py`에 막혔다: boot total **46,152 B > 46,000 B**. 내 커밋은 boot 컨텍스트에
**한 바이트도 더하지 않는다.** 개별 예산은 셋 다 통과했고(총합만 152 B 초과) 커진 것은
`7daedb3d`가 `CLAUDE.md`에 금지선을 추가하고 `0485c6d7`가 `STATE.md`를 재작성한 결과다 —
**다른 세션의 것이고, 그 세션은 지금도 `STATE.md`를 dirty로 들고 있다.**

`TOTAL_BUDGET`을 올리는 것은 PI 결정이고 **하지 않았다.** 다른 세션의 in-flight 파일도 건드리지
않았다. 대신 내가 소유한 유일한 boot 파일인 `MEMORY.md`의 preamble을 **그 파일 자신의 규칙에
맞게** 줄였다 — 그 preamble은 "숫자·날짜·커밋·발견을 여기 적지 말라"고 하면서 스스로 철회된
magnitude 3개와 토큰 수 하나를 담고 있었다. 게이트와 무관하게 옳은 정리다. 결과 **45,670 B**.

**PI에게 올리는 조건:** 총합이 천장에 붙어 있다. 다음에 `CLAUDE.md`나 `STATE.md`가 조금이라도
커지면 **모든 세션의 커밋이 막힌다.** 예산을 올릴지, 무엇을 boot에서 빼낼지는 PI 결정이다.

> **후속 (같은 세션, 라운드 2 시작 시점):** 예측대로 **다시 빨개졌다** — 46,031 B. 이번에는
> `MEMORY.md`를 **또 줄이지 않았다.** 다른 세션의 증가를 반복해서 흡수하면 문제가 보이지 않게
> 될 뿐이다. 이후 그 세션이 자기 쪽에서 해소해 45,920 B가 되었다. **조건 자체는 그대로다.**

---

## 16. Round 2 — sandbox 실행 (진행 중)

계획서 §12의 첫 실행 묶음 중 **PI 대기에 막히지 않은** 항목만 골랐다: 항목 9(sandbox
S0/S1/S11/S12)과 항목 10(감사 A1–A15를 체크리스트로).

### 16.1 S1 착지 — **q⁻² 스펙트럼은 Helfrich fit이 절대 기각하지 못한다. 그것이 Helfrich이기 때문이다**

브리프는 "q⁻⁴와 q⁻²를 가르려면 q가 몇 decade 필요한가"를 물었다. exponent 검정은 답한다
(N=100/mode에서 **0.20 decade**, N=4000에서 **0.05 decade**, 예측 5σ span과 15% 이내 일치).
**그런데 정작 중요한 측정은 그 질문이 조금 빗나갔음을 말한다:**

| 검정 | q⁻² 데이터에 대한 기각률 | 판정 |
|---|---|---|
| 2-parameter Helfrich goodness-of-fit | 1.0 / 2.0 / 4.0 decade 전부 **0.010 = nominal α** | **구조적으로 못 잡는다** |
| exponent 구간 검정 | 0.05~0.20 decade면 200/200 | 잡는다 |
| 족(family) 밖 스펙트럼(q⁻³)에 대한 χ² | N=1000에서 **0.75 decade** 필요 | 잡는다 |

q⁻²는 **κ=0인 Helfrich 스펙트럼**이다. 따라서 fit 품질은 "tension을 bending으로 오인"을 원리적으로
검출할 수 없고, exponent 구간이나 명시적 `kappa_resolved` 플래그만이 할 수 있다.
`HELFRICH_FAMILY_ABSORPTION`으로 기록.

**그래서 native run의 진짜 구속은 exponent 검정이 아니라 두 modulus를 함께 분해하는 것**이고,
그건 crossover가 대역 안에 있어야 한다. 기하로 환산:

```
q* = 1.09897e7 m⁻¹  (λ* = 571.7 nm)   [κ=8.28e-20 J, σ=1.0e-5 N/m, 310 K]
±1 decade about q*  ⇒ 5.717 µm patch @ 28.59 nm spacing = 201×201 = 40,401 nodes
±0.5 decade         ⇒ 1.808 µm @ 90.4 nm = 441 nodes
±1.5 decade         ⇒ 18.08 µm @ 9.04 nm = 4.0M nodes
```

exponent 검정만 보고 잡았을 span의 **약 10배**다.

그 밖의 측정: κ·σ를 0.561%/0.704%로 회수(참값이 95% CI 안), 구간 coverage **114/120 = 0.950**
(nominal 0.95) · 1-D string이 mode당 `k_BT/2`를 1.002544로, 관측 산포 0.00903 vs 예측 0.01000 ·
**amplitude rescale은 (κ,σ)와 정확히 축퇴**(χ²가 상대 1.3e-14만 움직임)라 equipartition은
독립 modulus에 대해서만 검정 가능 · **correlated mode는 variance fit(p=0.521)과
equipartition(0.9992)을 둘 다 통과**하고 독립성 검사만이 잡는다(그래서 별도 검사) ·
**mesh cutoff를 관통해 fit한 κ는 참값의 33.63배**, trusted prefix에서는 1.0049배.

q의 잘못된 거듭제곱은 fit되지 않고 **raise한다** — 공식이 4-exponent 차원 대수로 조립된다.

숨기지 않고 적은 한계: κ,σ가 비음수로 제약되므로 참값이 경계에 있으면 유효 파라미터가 2보다
적은데 `dof`는 여전히 M−2를 보고한다 → mode 수가 적으면 goodness-of-fit이 **과기각**한다
(7 mode에서 0.130, 21 mode에서 nominal 0.010 회복). 오차가 안전한 방향이지만
**7-mode p-value는 calibrated p-value가 아니다.**

검증: 66 tests, 전체 `tests/ac` 1839 passed / 0 failed.

### 16.2 ⚠️ 차원 검사가 **내 소유가 아닌 파일에서** 문서 결함 2건을 찾았다 (Lead 독립 확인)

둘 다 산술로 확인했고, **소유 레인이 아니므로 고치지 않고 올린다.**

**(1) `aleph/docs/v2_audit/BLEBBING_VALIDATION_PLAN_2026-07-20.md` — 단위 환산이 4.1배 틀렸다.**

문서는 `kappa_m` 0.0828 pN·µm을 `(≈3.4e-19 J)`라고 적는다. 실제로
`0.0828 pN·µm = 0.0828 × 1e-18 J = **8.28e-20 J**`이다. 그리고 `3.4e-19 J`는 310 K에서
**≈79 k_BT**로, 같은 줄이 주장하는 `20 kBT`도 아니고 KB 자신이 말하는 밴드(0.4–1.2e-19 J)
**밖**이다.

**(2) `aleph/components/incumbent/compartments.py:122` — `20 kBT`가 어느 온도인지 적혀 있지 않고, 값은 300 K다.**

```
20·k_B·300 / 1e-18 = 0.082839 pN·µm   ← 코드의 0.0828과 일치
20·k_B·310 / 1e-18 = 0.085600 pN·µm   ← 생리 baseline이라면 이 값
```

차이는 3.3%다. 크지 않지만 **"모든 파라미터는 실제 in-vivo 값에서 시작한다"는 hard rule과
온도 규약이 어긋나 있고, 주석이 온도를 말하지 않아 검출되지 않았다.** 상수를 바꿀지 주석에
300 K를 명시할지는 소유 레인과 PI의 몫이다.

**PI-decision (S1 provenance):** κ·σ는 `MEMBRANE_DEFAULTS`에서 인용 KB-3.B1.2(Rawicz 2000) /
KB-3.B1.1(Diz-Muñoz 2013)과 함께 가져왔으나 **그 논문들을 읽지 않았다.**
`INHERITED-UNVERIFIED`로 표기했고, 모든 진입점이 κ/σ를 인자로 받으며 **어떤 것도 이 값으로
기본값을 삼지 않는다.**

### 16.3 S0 착지 — 가짜 limit cycle이 **해석적 dt\*를 품은 구간**에서 나타난다

battery를 공허하지 않게 만드는 절반은 **반드시 실패해야 하는 positive control**이고, 그것이
이론이 말하는 스텝에서 실패했다. 감쇠 진동자(f₀=0.5 Hz, ζ=0.3)에 explicit Euler:

```
dt = 0.01 / 0.02 / 0.04 / 0.15  → fixed-point
dt = 0.20                       → limit-cycle          verdict FLAGGED_REGIME_ARTEFACT
전이 구간 (0.15, 0.20) s ∋ dt* = 2ζ/ω = 0.190986 s     (Lead 확인)
per-step growth 0.969187 → 1.00886,  상대진폭 9.97e-4 → 186
```

**같은 계·같은 다섯 스텝에서 semi-implicit은 전 구간 `fixed-point`이고 `INVARIANT`를 반환한다** —
플래그가 계가 아니라 **방법**에 관한 것임을 보이는 것이 이 held-out이다.

**이 프로젝트가 실제로 밟은 함정:** per-step solver 잔차를 *physical로 선언*하면 tolerance
ladder에서 3.49e-5 → 3.49e-11(동적 범위 9.998e5, z_ess=217)로 움직이고 `FLAGGED_RESIDUAL_TRACKING`으로
라우팅된다. 정직하게 `SOLVER_DIAGNOSTIC`으로 선언하면 같은 숫자가 보고되되 **단언되지 않는다.**
ladder에 불변 anchor가 없으면 일반 statistic-artefact로 강등되고 **구분 못 함을 스스로 말한다.**

**ESS 함정(측정):** sem_ess는 8× dt ladder에서 1.05%까지 dt-불변인데 sem_naive는 반감마다 정확히
√2씩 줄어든다. held-out semi-implicit ladder에서 **실제 위양성**이 나왔다 — variance가
z_ess=1.734(정합)인데 z_naive=3.938로 **같은 선언 한계 3을 넘는다.** 둘 다 보고하고 z_ess만 행동한다.

**프로토콜 거부 2건(테스트됨):** dt ladder를 **step 수**로 맞추면(물리 duration 불일치)
`REFUSED_UNMATCHED_DURATION`으로 **어떤 rung도 재지 않으며**, 사유가 blebbistatin 철회
(`STATE.md` (c) 16)를 지목한다. 자기 τ_int의 20배 미만을 덮는 ladder는 증거 부족으로 거부.

**개정안이 필요로 하는 깨끗한 분리, 한 측정으로:** Euler–Maruyama 분산 편향이 closed form과
**0.12%** 일치(측정 1.29317 vs 해석 1.29167)한다 — 즉 **regime은 불변인데 physical statistic은
5.9σ로 dt-편향**되었다. **falsifier는 앞의 것뿐이다.**

seed: 8개 독립 run, modal fraction 0.875 ≥ 선언 0.80 · ESS calibration ratio **1.331**(선언 밴드
0.5–2.0), naive error면 **5.87**로 밴드 밖 · 평균 쌍별 궤적 차 **1.401σ ≈ √2**(독립성 증거).
**궤적 불변성은 어디에도 단언하지 않는다** — 보고 타입에 필드가 없고, seed 8개 미만이나 중복
seed는 거부한다. contract 11필드 전부 기본값 없음.

### 16.4 S11 착지 — **식별성을 사주는 modality가 곧 forward-model 오차에 노출시키는 modality다**

축퇴를 가정하지 않고 **쟀다.** 단일 widefield frame의 Fisher 고윳값
`[2.5e-12, 0.89, 1.6e3, 5.6e3, 1.8e4]`, 조건수 **7.1e15**. 최악 고유방향이 closed-form defocus
null과 **|cos| = 1.0000**으로 겹친다 — PSF 폭과 축방향 위치가 이미지에 **하나의 수로만** 도달하므로
근사 축퇴가 아니라 **정확한 구조적 축퇴**다. "0"과 "작음"은 차분 스텝을 반으로 줄여 구분했다
(null 고윳값 66% 이동, 결정된 넷은 상대 3.1e-8 이하).

`L/s0 = 0.25…8`에서 flux 방향 정보는 최선 방향 대비 **3%를 한 번도 넘지 않는다** — 단일
widefield frame에서 길이와 밀도는 **그 범위 어느 스케일에서도** 깨끗이 분리되지 않는다.

**대조(probe 방향은 augmented 결과를 보기 전에 고정):**

| 추가 modality | 조건수 | defocus null | flux 방향 |
|---|---|---|---|
| — (단일 평면) | 7.1e15 | 1.7e-12 | 1.80 |
| +2번째 초점면 | 1.7e4 | **×2.5e15 해소** | 4.10 (×2.3, **미해소**) |
| +독립 40 nm 채널 | 729 | **해소** | **×41.8 해소** |

**하나의 modality가 하나의 축퇴만 푼다.** 그 대조가 산출물이다.

**그리고 예상하지 못한 비용.** 렌더러 불일치에 대해 단일 평면 posterior는 거의 안 움직인다
(aliased point-chain truth는 coverage 0.8점 이하, 15% 축방향 오보정은 3.3점 이하) — **움직여야 할
파라미터가 어차피 축퇴이기 때문**이다. 같은 15% 오보정을 **2평면 설계**에 걸면 coverage가
`log_psf_sigma_um` **0.950 → 0.671**, `z_um` **0.929 → 0.592**로 붕괴한다.
**식별성을 사준 modality가 정확히 forward-model 오차에 노출시킨 modality다** →
observation-design 이득과 forward-model calibration 요구가 **결합되어 있다.** 이것이 이 폐루프의
정직한 천장이다.

**falsifier가 약한 형태로 발화했고 다듬지 않고 기록했다:** 1000 draw에서 nominal 95% 구간의
coverage가 `0.933…0.957` — 최대 2.8점 실제 결손, nonlinear sampler로 라우팅. 테스트의 선언
240 draw에서는 한 파라미터가 3σ 이항 밴드 **밖**이고 **테스트가 그 불리한 사실을 단언한다.**
과신 대조군(σ×0.5)은 |z| = 19.3으로 크게 잡힌다.

**A5는 구조적으로 강제된다:** `point_estimate`가 caller 지정 identifiable fraction·고윳값 floor
(둘 다 기본값 없음) 아래에서 `UnidentifiableDirectionError`를 raise한다. 2평면도 길이/밀도는
여전히 거부하고 superres 채널은 다섯 전부 허용한다. **credible interval은 항상 제공된다** —
거부되는 것은 point estimate이지 추론이 아니다.

가는 길에 고친 실제 결함 2건(실재하므로 보고): Fisher scoring이 ~1% draw에서 발산 → 선언된
prior-σ box로 제한하되 hit를 **세어서 보고**(삼키지 않음); 첫 line search가 quasi-score에 없는
`ln v` gradient까지 최소화해 자기 스텝과 싸웠다(~150/300 draw 조기 정지) → outer iteration마다
가중치 고정으로 1–9/300, coverage 결손 ~6점 → ~1점.

### 16.5 A1–A15 감사 register — **A13의 판정이 존재한 적 없는 타입을 지목한다**

§10.2가 산문에서 **기계 검사 가능한 register**로: A1–A15 + §10.3 제거 3건(R1–R3) = 18행.
**status는 선언이 아니라 파생된다** — `AuditItem.status`가 읽기 전용 property이고, 강제를
주장하는 항목은 **실제 예외를 지목하는 binding과 residual gap을 둘 다** 지녀야 하며 각각
음성대조 테스트가 있다.

**결과: 10 MECHANICALLY_ENFORCED / 2 EVIDENCED / 6 UNENFORCED.**

⭐ **A13의 판정문은 "모든 도구는 `CellStatePosterior` 또는 `ObservationOperator`를 입출력해야
한다"인데, Lead 확인 결과 두 타입 모두 `aleph/virtual_cell` 어디에도 존재하지 않는다.**
`.py`에서 그 이름이 나오는 곳은 **register 자신이 부재를 기록한 문장뿐**이다. 즉 계획서의 통합
판정이 **작성된 적 없는 인터페이스를 지목**하고 있어 그대로는 강제 불가능하다. 어느 쪽이든
나타나는 순간 실패하는 tripwire 테스트를 걸어 register가 조용히 낡지 않게 했다.
`contracts.ArtifactEnvelope`가 사실상 모든 레인이 내보내는 공유 스키마지만 **소비를 요구하는 것은
없고**, 판정문이 지목한 타입도 아니다.

강제 10건은 전부 테스트에서 **실제로 발화시켰다**(surrogate `evidence_source` 덮어쓰기에
`TypeError`, `HeldOutValidationRefused`, `CensusUnavailableError`, manifest 키 `"MCF7"`에 raise 등).
**저장소 문자열 스캔은 어디에도 없다** — 대상을 참조로 들고 `importlib`로 재해석해 **객체 동일성**을
단언하므로 rename이 조용한 강등이 아니라 **큰 소리의 실패**가 된다.

EVIDENCED 2건 중 A1은 여기서 재현했다: thermal(σ=0.1)과 biological(σ=1.0)을 moment matching으로
한 Gaussian에 합치면 **분산은 0.505로 정확히 맞고**(그래서 2차 모멘트로는 병합을 볼 수 없다)
`E[x⁴]`를 **1.9608배** 잘못 보고한다.

UNENFORCED 6건이 **다음 빌드 큐**다(부재가 아니라): A5(역방향 map이 없어 point estimate를 거부할
주체가 없다 — S11이 그 일), A10(이 레이어에서 native run이 나온 적 없어 비용 모델이 없다),
A13(위), A14·A15(PI), R1(의도적 — `RepresentationKind`를 9개 고전 멤버로 고정하면 양자 표현의
*선언*은 막지만 propagator *구현*은 못 막고, register가 그렇게 적는다).

### 16.6 라운드 2가 라운드 1의 실제 크래시를 찾았다 (`77eea85e`)

S0이 발산을 **플래그하려다** traceback을 받았다. 원인은 미묘하다:

**`_as_series`는 모든 샘플이 유한함을 보장한다. 그것은 모든 통계가 유한하다는 뜻이 아니다.**
overflow 직전의 explicit scheme은 `np.isfinite`를 전부 통과하는 1e200 근방 샘플을 만들고,
그 다음 **2차 모멘트가 overflow**한다 → variance `inf` → 정규화 ACF `nan` → `tau_int` `nan` →
60줄 아래 Theiler window의 `round(nan)`이 `cannot convert float NaN to integer`를 raise한다.
원인에서 멀고 **라이브러리 버그처럼 읽힌다.**

재현: `standard_normal(3000) * 1e200`, 전부 유한, max 3.886e200 → 수정 전 raise, 수정 후
UNDECIDED + 사유.

**내 첫 수정은 틀렸고 그 점이 기록될 가치가 있다.** 입력에 비유한 가드를 넣었는데 `_as_series`가
이미 비유한 입력을 거부하므로 **절대 발화할 수 없는 죽은 코드**였다. 가드는 통계가 실제로
비유한이 되는 지점으로 옮겼고, **회귀 테스트가 샘플이 `np.isfinite`를 통과함을 단언**하므로
훗날 누군가 입력 쪽에서 "고쳐" 실제 케이스를 공허하게 통과시킬 수 없다. 리포트가 담는 통계도
처음엔 0 배열을 넣었다가 **실제 시계열**로 바꿨다 — 0으로 대체하면 평균도 분산도 없는 샘플에
평균과 분산을 **날조**한다.

### 16.7 라운드 2 통합 (`5cf5c9f7`) — 재수출하면 안 되는 **두 번째** 이름

23모듈 **351 public name**, 완전성 검증(모든 `__all__` 재수출·전부 resolve·중복 0·Warp 미초기화).

라운드 1은 충돌 하나를 찾았다. 라운드 2는 **종류가 더 나쁜** 두 번째를 찾았다. 세 sandbox 모듈이
각각 `EVIDENCE_AUTHORITY`를 정의하는데 **값의 종류가 다르다**:

```
sandbox_invariance   'unverified-observer'          ← 권위 라벨
sandbox_fluctuation  'unverified-observer'          ← 권위 라벨
sandbox_inadequacy   'oracle-only: the truth is constructed in closed form. ...'  ← 산문 서술
```

패키지 수준의 단일 `EVIDENCE_AUTHORITY`는 **caller가 무엇을 읽고 있는지에 대한 거짓말**이 된다.
게다가 이곳은 **증거 provenance 어휘**라 모호한 이름의 피해가 가장 큰 자리다. `ExpansionRoute`와
같은 처리를 했다(모듈 한정 alias, bare 이름 부재). 세 모듈을 하나의 어휘로 정리하는 것은
cross-lane 변경이라 **통합 단계가 승자를 고르지 않는다.**

**전체 `tests/ac`: 1969 passed / 87 skipped / 0 failed.**
