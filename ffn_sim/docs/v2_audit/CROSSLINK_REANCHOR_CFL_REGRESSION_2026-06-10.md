# 잠복 회귀: loop18 crosslink 재anchor가 crosslinker batch-CFL 테스트 2개를 깨뜨림

> Lead 자율-루프(2026-06-10 iter13) 회귀 테스트 중 발견. **수정 안 함 — CFL contract
> 테스트라 gate-loosening 금지(하드룰), PI surface.** Production은 안전(γ 결론 무관).

## 증상
`python -m pytest ffn_sim/tests/test_crosslinkers.py` → **2 FAILED**:
- `TestBoundary::test_batch_cfl_shrinks_when_violated` — `k_off_max=2.557e+39/s`라 어떤
  batch_steps로도 CFL(batch_dt·k_off_max ≤ 1e-3) 충족 불가.
- `TestKinetics::test_demo_xlink_sim_builds_and_runs_no_NaN` — `XlinkBondUpdater.__init__`가
  `RuntimeError: violates D2 batch CFL: batch_dt·k_off_max = 1.3e-8·2.557e39 = 3.3e31 > 1e-3`.

## Root cause (정밀)
batch-CFL은 envelope 힘 `_F_env = k_attach · max_bind_dist` (crosslinkers.py:277)로
Bell-Evans off-rate ceiling `k_off_max = k_off0·exp(x·F_env/kT)`을 잡는다.

| | k_attach | max_bind_dist | _F_env | k_off_max |
|---|---|---|---|---|
| **Production** (phase1_h3.yaml) | 1.0e-3 | 60 nm | **60 pN** | **18 /s** ✅ batch=100 정상 |
| **Demo/test** (test_crosslinkers.py:70, `max_bind_dist=1µm`) | 1.0e-3 | 1 µm | **1 nN** | **2.56e+39 /s** ❌ |

- 데모는 50-filament sparse 빌드용으로 bind 반경을 **1µm**로 넓힘(60nm의 16.7×).
- **OLD k=1e-7**에선 _F_env_demo = 1e-7·1µm = 1e-13 N = 0.1pN → 무해(주석 264-275가
  "k_attach is tiny, amplification≈1" 명시).
- **loop18 재anchor(k 1e-7→1e-3, ×10⁴)** 후 _F_env_demo = 1e-3·1µm = **1nN** → Bell-Evans
  exp(~6.5nm·1nN/kT)=exp(~90) → **1e39/s**. CFL 폭발.
- ⇒ `_F_env = k·max_bind_dist` envelope는 **stiff k에서 비물리적**: k=1e-3 본드는 1µm까지
  늘어나지 않는다(1nN는 불가, 수~수십 pN에서 파열). soft-k 가정으로 작성된 envelope가
  재anchor 때 갱신되지 않은 것.

## 영향 범위
- ✅ **Production 안전**: 60nm 반경에서 k_off_max=18/s, batch_steps=100 유지. **이번 세션의
  γ A/B·density·앙상블 결론은 영향 없음**(production 작동점에서 측정).
- ❌ **Test suite RED**: 2개 contract 테스트 실패. 데모-scale(넓은 bind 반경) 빌드는 현재
  재anchored k에서 빌드 불가.
- ⚠️ 잠재: 향후 어떤 production-scale 빌드가 max_bind_dist를 넓히거나 큰 x_β crosslinker
  파라미터를 쓰면 batch_steps→1 강제(100× perf 회귀) 또는 빌드 raise 가능.

## PI 결정 필요 (수정안 — Lead는 진행 안 함)
envelope 힘 모델 vs 테스트 fixture 중 무엇을 고칠지 = CFL-contract 판단:
- (A) **envelope를 물리적 rupture force로 cap** — `_F_env = min(k·max_bind_dist, F_rupture)`
  (catch-bond 특성 rupture ~수십 pN). 모델 측 수정, production-wide, CFL contract 변경.
- (B) **thermal-displacement bound 사용** — `_F_env = k·√(kT/k) = √(kT·k)` (실제 본드가
  겪는 열적 변형 규모; 주석 268이 이미 언급). k=1e-3 → √(4.28e-21·1e-3)=2.07e-12=2pN →
  k_off_max 정상. max_bind_dist(기하 인자)와 분리.
- (C) **데모 fixture의 1µm bind 반경 축소** — 테스트만 수정(production 무관). 단 contract
  테스트 편집이라 PI sign-off 필요.

Lead 추천: **(B)** — envelope를 기하 bind-반경이 아니라 실제 열적 변형 규모로 잡는 게
물리적으로 맞고(stiff/soft k 모두 견고) max_bind_dist 변화에 robust. 단 CFL contract 변경이라
PI sign-off. 그 다음 데모 fixture 재검증.

## 재현
```
conda activate ffn_sim
python -m pytest ffn_sim/tests/test_crosslinkers.py::TestBoundary::test_batch_cfl_shrinks_when_violated -q
```
