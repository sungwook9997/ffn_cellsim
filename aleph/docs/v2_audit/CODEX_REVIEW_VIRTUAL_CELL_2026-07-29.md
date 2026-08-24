# Codex advisory review — `virtual_cell` (2026-07-29)

**요청:** `CODEX_REVIEW_REQUEST_VIRTUAL_CELL_2026-07-29.md`  
**대상:** `a51a3c1c..abf215e0` (`aleph/virtual_cell`, `aleph/tests/ac/virtual_cell`)  
**성격:** 결함 탐색 결과이며 검증 권위나 native evidence가 아니다.

`abf215e0..HEAD`에서 대상 두 디렉터리의 변경이 없음을 확인했으므로 현재 worktree의 해당 파일은
요청 범위와 같다. CUDA가 없는 개발 머신에서 CPU oracle과 정적 계약만 점검했으며 simulation
physics 결론은 내리지 않았다.

## 기준선

- `python -m pytest aleph/tests/ac -p no:randomly`: **2139 passed, 87 skipped**
- `ruff check aleph/virtual_cell aleph/tests/ac/virtual_cell`: **PASS**
- `make kb-check`: **EXIT=0**

## 발견

### 1

[결함]  
파일: `aleph/virtual_cell/cell_state_posterior.py:1207-1209,1446-1481`;
`aleph/tests/ac/virtual_cell/test_cell_state_posterior.py:1075-1080`  
주장: `record_hash`는 “whole record”의 SHA-256이므로 posterior의 어느 부분이 달라도 달라진다.  
반례 또는 재현: 서로 다른 `EvidenceGraph`를 넣고 `AdmissibilityVerdict.attempted_step_index`를
`7 -> 999`로 바꾼 두 posterior를 만들었다. `evidence.graph_hash`는 달랐지만
`record_hash_equal=True`였다. `to_dict()`에는 `evidence`가 전혀 없고 accepted-step verdict는
`admissible` boolean 하나만 들어간다. 이름이 “any part”인 테스트도 `posterior_id` 하나만 바꾼다.  
왜 중요한가: 다른 evidence DAG, 다른 accepted-step transaction을 같은 content-addressed posterior로
오인한다. source ID의 다른 그래프 재사용도 hash에서 보이지 않아 10-lane composition의 핵심
provenance가 명목상 필드로만 남는다.

### 2

[결함]  
파일: `aleph/virtual_cell/observation_operator.py:344-401,724-736,889-930`  
주장: 관측 identity는 operator/configuration과 “observed state”에 결합된다.  
반례 또는 재현: `state_id`, `manifest_hash`, `source_id`, field 이름과 단위를 같게 두고
`tension=[1,2,3]`과 `[1,2,30]`만 다른 두 `MechanisticState`를 같은 projection으로 관측했다.
값은 `3.0`과 `30.0`인데 `state.provenance_digest`, envelope `artifact_id`,
`observation_provenance_sha256`가 모두 같았다. state mapping은 얕게만 복사하므로 원본 ndarray를
사후 변경해도 같은 충돌을 만들 수 있다.  
왜 중요한가: 서로 다른 관측값이 같은 artifact identity를 가진다. upstream caller가 state ID를
완벽히 새로 발급한다는 문서 규율만 있고 schema가 이를 강제하지 않는다.

### 3

[결함]  
파일: `aleph/virtual_cell/cell_state_posterior.py:1301-1310,1419-1427`;
`aleph/virtual_cell/observation_operator.py:235-289,533-575,964-1023`  
주장: blocked modality는 `UnimplementedOperator`로만 표현되고 관측 시 항상 raise한다.  
반례 또는 재현: 기존 `SummaryProjectionOperator`를 subclass하여 `modality`를 `TRACTION`,
`operator_id`를 `unimplemented::traction`으로 바꿨다. runtime-checkable
`ObservationOperator`를 만족해 posterior에 등록되었고, `posterior.observe(...)`는
`reduced-model` evidence의 scalar `5.0`을 반환했다.  
왜 중요한가: capability table이 `BLOCKED`라고 하는 traction을 동일한 reserved ID의 숫자 반환
operator로 조용히 바꿀 수 있다. “unimplemented는 숫자를 만들 수 없다”는 강제가 class-local이고
registry-global이 아니다.

### 4

[결함]  
파일: `aleph/virtual_cell/cell_state_posterior.py:1321-1334`;
`aleph/tests/ac/virtual_cell/test_cell_state_posterior.py:802-811`;
`aleph/virtual_cell/surrogate.py:700-708,748-797`  
주장: `__slots__` 위 property라 `native_accepted=True`를 “spell”할 수 없다.  
반례 또는 재현: `class ForgedPosterior(CellStatePosterior)`에서 `native_accepted` property를
`True`로 override하고 base dataclass fields로 생성했다. `isinstance(..., CellStatePosterior)`가
참이고 `to_dict()["native_accepted"]`도 참이었다. 같은 방식으로 `SurrogatePrediction`의
`evidence_source` property도 override된다.  
왜 중요한가: posterior는 실제로 native self-promotion을 serialize한다. 현재 테스트는
`dataclasses.replace`만 공격한다. Surrogate는 envelope 생성 시 `SURROGATE`를 다시 hard-code하므로
최종 envelope laundering은 막지만 “그 값은 spell 불가능”이라는 더 강한 주장은 거짓이다.

### 5

[결함]  
파일: `aleph/virtual_cell/stochastic_node.py:1220-1223,1225-1298`  
주장: closure breakdown bisection은 `gap_cancellation_floor` 아래의 차이를 신뢰하지 않는다.  
반례 또는 재현: 상수 `1e6`과 `exp(-x)`를 합친 law, tolerance `0.01`, search low `1e-4`,
high `2.0`에서 `cancellation_floor_at_low=0.0887139786 > 0.01`인데도 함수는 거부하지 않고
`breakdown_sigma=0.2003352948`을 반환했다. 코드가 `floor_low/floor_high`를 계산해 report에
기록하지만 어떤 branch나 bisection 조건에도 사용하지 않는다.  
왜 중요한가: sanity floor보다 작은 상대 오차를 실제 closure 정확도로 읽고 `σ <= 0.2003 λ`
같은 안전 범위를 발급할 수 있다.

### 6

[공허한테스트]  
파일: `aleph/tests/ac/virtual_cell/test_sandbox_fluctuation.py:325,613,631`  
주장: factor-of-two energy, independently supplied moduli, wrong spectrum object의 세 negative
control이 equipartition 경로를 깨뜨리는지 검사한다.  
반례 또는 재현: 세 함수 이름이 모두 `check_...`라 pytest가 수집하지 않는다.
`pytest --collect-only -q .../test_sandbox_fluctuation.py`는 66개만 보고하며 세 함수는 목록에 없다.
대상 테스트 트리 전체에서 top-level `check_`는 이 세 개뿐이다.  
왜 중요한가: 가장 고전적인 factor-of-two 오류와 독립 moduli rescale 방어가 baseline의 2139 PASS에
전혀 참여하지 않았다.

### 7

[과일반화]  
파일: `aleph/virtual_cell/optics.py:13-26,688-718`;
`aleph/virtual_cell/synthetic_microscopy.py:32-42`;
`aleph/tests/ac/virtual_cell/test_optics_sources.py:62-76,300-326,364-392`  
주장: line/triangle extended emitter image는 subdivision에 대해 round-off까지 불변이다.  
반례 또는 재현: 테스트 helper는 실제 기본값 `psf_radius_sigma=4` 대신 `8`을 강제한다. 기본 4σ로
같은 17-piece line과 64-triangle subdivision을 재현하면 max absolute difference가 각각
`8.90e-7`, `1.20e-7`, peak-relative difference가 `7.51e-6`, `1.50e-6`였다. 8σ에서만
`3.75e-16`, `6.94e-17`로 내려간다. `_render_extended`가 primitive마다 다른 truncation window를
잡아 analytic interior cancellation 뒤에 서로 다른 tail을 버리기 때문이다.  
왜 중요한가: analytic field의 telescoping은 맞지만 default renderer의 유한 support에는 그대로
전이되지 않는다. 현재 테스트가 기본 configuration에서의 회귀를 가린다.

### 8

[결함]  
파일: `aleph/virtual_cell/sandbox_inadequacy.py:1655-1661,1681-1705`;
`aleph/tests/ac/virtual_cell/test_sandbox_inadequacy.py:64-73,344-364`  
주장: null-control replicate와 D2 permutation은 “no two replicates share a stream”이다.  
반례 또는 재현: 테스트 계약은 `seed == permutation_seed == 20260729`다. replicate 0의 training
seed와 permutation seed가 같고, replicate 0의 held-out seed `20260730`은 replicate 1의
permutation seed와 같다.  
왜 중요한가: 동일 seed의 새 generator들은 독립 stream이 아니다. 측정된 FPR 자체가 즉시
무효라는 뜻은 아니지만, 독립성에 근거한 binomial 해석의 전제가 코드가 주장한 형태로 성립하지 않는다.

### 9

[과일반화]  
파일: `aleph/tests/ac/virtual_cell/test_sandbox_fluctuation.py:737-772`  
주장: mesh cutoff를 관통한 κ 적합의 `33.63x` 손상이 대표적인 배수다.  
반례 또는 재현: 같은 설계에서 cutoff만 바꾸자 κ ratio가 `8e7 -> 193457`,
`1.0e8 -> 476.5`, `1.2e8 -> 33.63`, `1.5e8 -> 5.36`, `2e8 -> 1.81`,
`3e8 -> 1.13`으로 바뀌었다.  
왜 중요한가: `33.63`은 이 합성 cutoff 위치의 수치이지 전이 가능한 손상 배수가 아니다. 다만
커밋된 테스트 계약은 정확한 배수를 고정하지 않고 `>5`와 trusted-range 회복만 검사하므로 그
테스트 자체는 이 반례로 깨지지 않는다.

### 10

[의심스러움]  
파일: `aleph/virtual_cell/sandbox_fluctuation.py:1326-1375`;
`aleph/tests/ac/virtual_cell/test_sandbox_fluctuation.py:831-846`  
주장: `q_min=2π/L`, `q_max=π/a`로부터 201×201 = 40,401 node native patch가 나온다.  
반례 또는 재현: 산술은 `L/a = 2 q_max/q_min = 200`, 코드가 endpoint를 더해 201이다. 그러나
`q_min=2π/L`은 periodic spectrum의 fundamental이고 periodic grid는 중복 endpoint 없이
200 unique nodes를 쓴다. endpoint-inclusive 201-node grid라면 boundary condition에 따라
fundamental 관계가 달라진다.  
왜 중요한가: 5.717 µm와 28.59 nm 산술은 맞지만 40,401은 periodic/nonperiodic sampling
convention을 혼합했을 가능성이 있다. native-run 설계 입력으로 사용하기 전에 boundary condition을
선언해야 한다.

## 반박 실패 및 독립 재현

### Stochastic A/B/D

[반박실패]  
파일: `aleph/virtual_cell/stochastic_node.py:409-417,420-430`  
주장: Gaussian raw moments와 tilt identity가 맞고, quenched/thermal split의 mean 불변 및 양 극단
variance ratio `N`은 exponential artifact가 아니다.  
반례 또는 재현: SciPy `quad`로 `bσ=3.58, 6, 6`, raw order 12, tilted order 8까지 독립 적분했다.
최대 상대 오차는 raw `4.1e-15`, tilted `3.6e-15`였다. 임의 force law에서도 각 bond의 marginal
`q+t_i` 분포는 total variance만으로 같으므로 mean 불변이고, all-shared 대 all-iid 평균의 variance
ratio는 정확히 `N`이다. predicted exponent 경로와 fitted exponent 경로도 분리되어 있었다.  
왜 중요한가: 이 부분은 다음 검수에서 같은 Gaussian 항등식 반박을 반복할 이유가 없다.

### S11

[반박실패]  
파일: `aleph/virtual_cell/sandbox_inverse.py:84-113`  
주장: 가장 작은 Fisher eigenvalue는 구조적 0이고, reverse 15% miscalibration도 two-plane
coverage를 무너뜨리며, n=1000 deficit은 숨기지 않는다.  
반례 또는 재현: 차분 step `1e-2 -> 1e-6`에서 null eigenvalue는 `3.45e-5 -> -4.34e-14`로
움직였지만 analytic null overlap은 `0.9999999684 -> 1.0`이었다. 좌표 scale을 바꾸면 condition
number가 `7.1e15`에서 `~1e21`까지 변해 그 숫자 자체는 불변량이 아니나, module은 이미 step
sensitivity를 적어 두었다. 반대 방향 Rayleigh scale `0.85`에서도 two-plane coverage가
`0.629/0.688`로 붕괴했다. n=1000 coverage `0.933,0.940,0.957,0.938,0.922`와 `z=-4.1`도
module 문서가 명시적으로 falsifier firing으로 기록한다.  
왜 중요한가: 구조적 degeneracy 결론은 유지된다. 고유값과 condition number의 정확한 headline만
coordinate/step 의존 수치로 읽어야 한다.

### S1 family absorption

[반박실패]  
파일: `aleph/virtual_cell/sandbox_fluctuation.py:201-210`  
주장: pure q^-2는 two-parameter Helfrich family의 κ=0 member라 chi-square fit quality만으로는
reject할 수 없다.  
반례 또는 재현: 6/8 decade까지 확장한 200-trial 실험의 rejection rate는 각각 `0.015`, `0.0`으로
nominal sampling fluctuation 안이었다.  
왜 중요한가: “family member이므로 fit quality가 식별자가 아니다”라는 구조적 주장을 깨지 못했다.

### Telemetry

[반박실패]  
파일: `aleph/virtual_cell/telemetry.py:1068-1112,1161-1274`  
주장: byte ratio는 parameter ratio와 다른 serializer-specific 측정이고 reduced decode는 약 2배다.  
반례 또는 재현: serializer는 record에 `numpy-npz-uncompressed`로 명시된다. 10회 독립 호출에서
median slowdown은 `2.02..2.30`; 100 warm-up 후 1000회 dense-first/reduced-first/random 순서를
바꾼 결과 `2.249/2.247/2.296`이었다.  
왜 중요한가: 고정 실행 순서는 개선 여지가 있지만 인용된 “약 2배”를 뒤집을 cache/order 반례는
이번 환경에서 나오지 않았다.

### Optics 공식

[반박실패]  
파일: `aleph/virtual_cell/optics.py`  
주장: pixel-integrated erf 공식과 triangle Gaussian mass/Owen-T 식 자체가 맞다.  
반례 또는 재현: 독립 수치 적분으로 subpixel narrow-σ point PSF와 500개 무작위 triangle case를
비교했다. triangle worst absolute error는 약 `9.3e-15`였고 공식 오류는 찾지 못했다. 발견 7은
공식이 아니라 finite render-window 적용의 문제다.  
왜 중요한가: 다음 수정은 closed form을 교체하는 것이 아니라 default truncation과 subdivision
계약을 일치시키는 쪽이어야 한다.

### A13, 이름 충돌, overflow

[반박실패]  
파일: `aleph/virtual_cell/audit_checklist.py:810-868`;
`aleph/tests/ac/virtual_cell/test_audit_checklist.py:314-349,882-913`;
`aleph/virtual_cell/regime.py:1426-1448`  
주장: A13 flip과 proof-test search 확장은 상태 전이로서 정당하고, qualified aliases는 충돌을
숨기지 않으며, overflow guard는 도달 가능하다.  
반례 또는 재현: export 삭제 tripwire는 실제 회귀를 잡고 proof-test는 전체 lane에서 정확히 하나만
허용한다. A13의 universal “EVERY tool” 미강제는 `residual_gap`과 전용 테스트에 명시되어 있었다.
authority 값 불일치도 package docstring이 숨기지 않고 cross-lane reconciliation으로 남긴다.
`1e180/1e200/1e250` finite samples는 모두 correlation statistic을 nonfinite로 만들어 현재 guard에
도달했고, `tau_s`는 finite인데 후속 PSD/recurrence만 NaN이 되는 재현은 찾지 못했다.  
왜 중요한가: 이 세 항목은 현재 코드가 주장 범위를 제한해 둔 상태다. 단, 발견 3은 A13의
`UnimplementedOperator` class-local 보증보다 더 넓은 registry 우회다.

## 우선순위

1. posterior/observation content identity 결함(발견 1–2)
2. blocked modality와 native self-promotion 우회(발견 3–4)
3. cancellation floor 미적용(발견 5)
4. 수집되지 않는 negative controls와 default optics 회귀(발견 6–7)
5. RNG stream 계약과 설계 수치의 범위 한정(발견 8–10)

어떤 threshold, gate contract, census, `aleph/ac/**`, native physics도 변경하지 않았다.
