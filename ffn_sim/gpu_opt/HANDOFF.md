# HANDOFF — cupy 포트 완료, A5000 박스에서 검증만 남음

(이 문서는 개발 머신(Mac)에서 작업한 Claude 세션이 GPU 박스에서 이어받을 세션에게 남기는 인수인계입니다. 2026-06-11 작성.)

## 현재 상태: 구현·로컬 검증 완료

`PROMPT.md`의 과제 — `kernels_cpu.py`의 NumPy 커널 4개를 cupy로 포팅 — 는 **구현이 끝났고**,
Mac에서 가능한 모든 검증을 통과했습니다. 남은 것은 **A5000 박스에서 실제 cupy로
패리티 + 속도 측정 한 번**뿐입니다.

### 산출물 (이 폴더)
- `kernels_gpu.py` — cupy 포트 (최종본, 리뷰·수정 반영).
- `test_shim_parity.py` — cupy 없는 머신용 검증: `sys.modules`에 numpy 기반 가짜
  `cupy`/`cupyx`(scatter_add→`np.add.at`)를 심어서 `kernels_gpu.py`의 **수식/로직**을
  CPU 레퍼런스와 대조. (cupy API 자체는 검증 못 함 — 그건 GPU 박스 몫.)

### 구현 요약
- **K1 `mesh_pressure_forces`** — 직역. cross product는 cupy 버전 호환을 위해 수동 성분
  계산(`np.cross`와 동일한 FP 연산이라 bit-identical). per-group 부피 reduction과
  per-vertex 힘 scatter는 `cupyx.scatter_add`.
- **K2 `plane_well_forces`** — `cp.where` 체인으로 elementwise. mask 지원.
- **K3 `bond_spring_forces`** — 직역, ± `cupyx.scatter_add`. 빈 bonds 가드.
- **K4 `group_pair_forces`** — O(N²) 루프를 **균일 그리드 neighbor list**로 교체:
  빈 크기 `r_cut`, 셀 ID 정렬 + `searchsorted` 범위 조회(dense 셀 테이블 없음 —
  메모리는 입자 수에만 비례), 27-셀 스텐실, cumsum+searchsorted 세그먼트 전개로
  후보 페어 생성(`cp.repeat`에 ndarray repeats를 안 씀 — released cupy ≤14.0.x에서
  크래시하기 때문), `i<j` 중복 제거 → 그룹/컷오프 필터 → 레퍼런스와 bit-identical한
  페어 수식 → Newton-3 ± scatter. 파이썬 루프는 고정 27회뿐(N 비례 루프 없음).

### 검증 이력 (Mac)
- 멀티 에이전트 리뷰(4개 렌즈: cupy API / 수식 bit-identity / neighbor-list 완전성 /
  엣지케이스) → 확정 결함 2건 수정 완료 (위 K4 설명의 cp.repeat 회피, dense 테이블 제거).
- `test_shim_parity.py`: 11개 케이스 **전부 PASS**, max abs err ≤ 5.9e-21
  (40×42 / 200×162 / 2000pt 랜덤 stress + gid −1 / mask / 빈 bonds 포함.
  K4는 두 크기 모두 실제 O(N²) 레퍼런스와 대조).
- CPU 기준치 (Mac, 참고용): 200×162 (N=32,400)에서 `group_pair` **~7.6 s/call**,
  `mesh_pressure` 5.1 ms, `plane_well` 0.4 ms, `bond_spring` 0.07 ms.

## 해야 할 일 (GPU 박스, ~10분)

```bash
cd ~/kernel_port
# (이 폴더가 박스에 없으면 kernels_cpu.py, kernels_gpu.py, bench.py,
#  test_shim_parity.py, PROMPT.md, HANDOFF.md를 복사해 올 것)

python -c "import cupy; print(cupy.__version__); import cupy as cp; print(cp.cuda.runtime.getDeviceProperties(0)['name'])"

python bench.py --n-groups 40 --n-per 42      # 1) 패리티 (rtol 1e-9)
python bench.py --n-groups 200 --n-per 162    # 2) 대형 N 속도 측정
```

**합격 기준** (PROMPT.md): 4개 커널 모두 `parity OK` (maxabs ≲ 1e-12), 대형 N에서
cupy가 NumPy보다 빠를 것(특히 `group_pair`).

### 알아둘 것
- 대형 N 실행은 **CPU 쪽 `group_pair` 타이밍 루프 때문에 ~3분** 걸림(7.6s × 21회). 정상.
- GPU `scatter_add`는 atomic이라 누적 순서가 ULP 수준에서 비결정적 → maxabs가 0이 아니라
  1e-13~1e-12대로 나올 수 있음. rtol 1e-9 기준엔 충분히 안쪽. Newton-3 대칭은 같은 배열을
  ±로 scatter하므로 정확히 유지됨.
- 첫 호출에 CUDA 커널 JIT 컴파일 비용이 있지만 bench.py가 타이밍 전에 워밍업 호출을 함.
- **만약 GPU에서 에러가 나면** 의심 지점(numpy shim이 검증 못 한 cupy API — 리뷰에서
  이론상 OK 판정됐지만 실기기 미확인): ① `cupyx.scatter_add(F, idx, val)`의 (M,3) value +
  (M,) index 브로드캐스트, ② `cp.searchsorted`의 `side=` 인자, ③ boolean-mask 인덱싱
  (`pi[keep]` 등), ④ 0-d cupy 배열의 `int()` 변환. 수식은 건드리지 말고 해당 API 호출만
  등가 대체할 것.
- `group_pair` GPU가 기대보다 느리면: 27-오프셋 루프 안의 `int(cnt.sum())` 동기화가
  오프셋당 1회 있음(총 27회) — 그래도 CPU 7.6s 대비 수백 배 빨라야 정상.

## 마지막 deliverable (측정 후 채울 것)

PROMPT.md deliverable 3번 — 아래 문단의 빈칸을 bench 출력으로 채우면 끝:

> kernels_cpu.py의 NumPy 커널 4개를 cupy로 포팅했다(`kernels_gpu.py`). K1–K3은
> `cupyx.scatter_add` 기반 직역이고, 지배적 비용인 K4(`group_pair`)는 O(N²) 브루트포스를
> 빈 크기 `r_cut`의 균일 그리드 neighbor list(27-셀 스텐실, 정렬+searchsorted 페어 생성,
> 전 과정 벡터화)로 교체하되 페어별 수식은 bit-identical하게 유지했다. A5000에서
> `--n-groups 200 --n-per 162`(N=32,400) 기준 4개 커널 모두 parity OK(maxabs ___),
> 속도는 mesh_pressure ___×, plane_well ___×, bond_spring ___×, **group_pair ___×**
> (CPU ___ ms → GPU ___ ms).
