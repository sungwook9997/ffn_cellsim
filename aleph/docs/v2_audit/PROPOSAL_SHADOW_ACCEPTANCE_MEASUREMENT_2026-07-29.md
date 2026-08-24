# 제안 — accepted-step conjunction의 shadow 측정

**보내는 곳:** `observe-infra` (`aleph/engine/transaction.py`) · `cortex` (`aleph/components/incumbent/driver.py`)
**보내는 레인:** `virtual-cell` — 이 레인은 `aleph/ac/` 아래를 claim하지 않으므로 **구현하지 않고 제안만 한다.**
**관련:** `PI_DECISION_CARD_ADMISSIBILITY_2026-07-29.md` 의 **P4** · `STATE.md` (e) 1

## 1. 이 제안이 푸는 문제

PI 결정 (e) 1(gate 개정안)이 보류된 이유는 **정보 부족**이지 이견이 아니다. conjunction을 채택하면
지금 commit되는 스텝 중 **몇 %가 거부되는가?** 2%면 채택은 사소한 결정이고, 90%면 전혀 다른
결정이다. 그 숫자가 없어서 결정이 멈춰 있다.

카드는 이를 `BLOCKED`로 적어 두었다. 추정하지 않은 것은 옳았지만, **측정하지 않은 채로 두면
결정도 멈춘다.** 이 제안은 그 숫자를 얻는 최소 경로다.

## 2. 왜 지금 싸게 할 수 있는가

conjunction의 두 항이 **이미 둘 다 디바이스 배열로 존재한다.**

| 항 | 어디에 | 현재 쓰임 |
|---|---|---|
| adjoint closure | `ledger.balance_ok_d` | `transaction.step()`의 기본 predicate |
| inner convergence | `inner_d.converged_d` | `driver.py:1196`이 만들고 `:1203`에서 호스트로 읽음 |

둘을 결합하는 코드만 없다. 새 물리도, 새 solver도, 새 상태도 필요 없다.

## 3. 제안하는 변경 — 최소·가산적·기본 OFF

### 3.1 `transaction.step()`에 선택 인자 두 개

```python
def step(
    self,
    *,
    dt_phys: float,
    solve: Callable[[], None],
    accepted_d: wp.array | None = None,
    ...,
    shadow_converged_d: wp.array | None = None,   # (1,) int32, 호출자 소유
    shadow_counts_d: wp.array | None = None,      # (4,) int64, 호출자 소유
) -> None:
```

둘 다 `None`이면 **아무 일도 하지 않는다.** 기본값이 `None`이므로 옵트인하지 않은 런은
**byte-identical**이어야 하고, 그것이 이 제안의 수용 조건이다(D2 `_accumulate_all(omit=)`가
같은 방식으로 통과한 전례가 있다).

### 3.2 predicate 결정 직후, commit 직전에 커널 한 번

```python
predicate = accepted_d if accepted_d is not None else ledger.balance_ok_d

if shadow_counts_d is not None and shadow_converged_d is not None:
    wp.launch(_accumulate_shadow_agreement, dim=1,
              inputs=[predicate, shadow_converged_d, shadow_counts_d])

self.world.finalize_candidate(predicate, ...)   # 변경 없음
self.clock.advance(predicate, dt_phys)          # 변경 없음
```

**`predicate`는 건드리지 않는다.** commit과 clock은 지금과 완전히 동일한 값으로 구동된다.

### 3.3 커널

```python
@wp.kernel
def _accumulate_shadow_agreement(
    predicate: wp.array(dtype=wp.int32),   # 실제로 결정한 값
    converged: wp.array(dtype=wp.int32),   # inner solve가 수렴했는가
    counts: wp.array(dtype=wp.int64),      # (4,)
):
    p = predicate[0] != 0
    c = converged[0] != 0
    # 0: 둘 다 accept                        1: predicate accept / NOT converged  ← 구하는 숫자
    # 2: predicate reject / NOT converged    3: predicate reject / converged
    idx = wp.select(p, wp.select(c, 2, 3), wp.select(c, 1, 0))
    wp.atomic_add(counts, idx, wp.int64(1))
```

> ⚠️ **`wp.select`의 인자 순서는 삼항연산자와 반대다** — `wp.select(cond, a, b)`는 `cond`가
> **거짓일 때 `a`**, 참일 때 `b`를 돌려준다. 이 제안의 초안은 그 자리에서 버킷 2와 3을 뒤바꿔
> 썼고, 진리표를 손으로 돌려서 잡았다. 구하는 숫자인 bucket 1은 우연히 맞았지만 나머지 둘은
> 조용히 뒤집혀 있었을 것이다. **구현할 때 네 조합 전부에 대한 단위 테스트를 먼저 쓸 것.**

카운터는 **디바이스에 머문다.** 스텝마다 호스트로 읽지 않으므로 GPU-only 계약(프로파일러
게이트의 hot-loop D2H 0)을 깨지 않는다. 런 끝에 한 번만 읽는다.

### 3.4 드라이버 쪽

`driver.py`는 이미 `inner_d.converged_d`를 들고 있으므로 그대로 넘기면 된다. 새 상태 없음.

## 4. 이 측정이 내놓는 숫자

런 하나에 정수 4개, 그리고 그중 **하나가 P1 결정을 좌우한다.**

```
bucket 1 = predicate ACCEPT  &&  inner NOT converged
         = 지금 commit되고 있으나 conjunction이면 거부될 스텝 수

reject_fraction = bucket 1 / (bucket 0 + bucket 1)
```

`run-record@2`에 네 카운트와 이 비율을 그대로 싣는다.

## 5. 2단계 — 불일치의 *크기* (1단계가 나온 뒤에만)

비율만으로는 "아슬아슬하게 어긋난 것"과 "완전히 어긋난 것"이 구분되지 않는다. 1단계 숫자가
크면 2단계가 필요하다: 불일치 스텝에서 `residual / tolerance` 비를 디바이스에 누적(합·최댓값·
로그 히스토그램)한다. **1단계를 먼저 하고, 필요할 때만 한다** — 지금 둘 다 제안하면 결정에
필요 없는 작업까지 묶는 것이 된다.

## 6. ⚠️ 이 측정이 답하지 **못하는** 것 — 채택 전에 반드시 읽을 것

**shadow 측정은 "지금 술어가 만든 궤적 위에서" 두 판정이 얼마나 갈리는지를 잰다.**
conjunction을 실제로 채택하면 거부된 스텝이 commit되지 않으므로 **궤적 자체가 달라지고**,
그 이후의 불일치율은 여기서 잰 값이 아니다.

즉 이 숫자는 **채택 후 거부율의 추정치가 아니라, 현재 궤적이 얼마나 비수렴 스텝 위에 서 있는지의
측정**이다. 그 구분을 흐리면 이 프로젝트가 이미 두 번 밟은 함정 — 고정 스텝수에서 잰 값을 다른
조건으로 전이한 것(`STATE.md` (c) 16, (c) 17) — 과 같은 종류의 오류가 된다.

**그래도 결정에 충분한 이유:** bucket 1이 0에 가까우면 conjunction 채택의 위험이 낮다는 것은
그대로 성립한다(현재 궤적에서 거의 아무것도 안 걸리므로). 반대로 크면 **현재 결과의 상당 부분이
비수렴 스텝 위에 있다**는 뜻이고, 그것은 채택 여부와 무관하게 그 자체로 보고 가치가 있는 사실이다.
어느 쪽이든 결정이 정보 위에서 이루어진다.

## 7. 실행 비용

- 스텝당 커널 1회, 스레드 1개, atomic 1회. 측정 가능한 비용이 아니라고 **예상**하지만,
  그 예상은 이 제안의 근거가 아니다 — 옵트인 OFF에서 byte-identical A/B가 근거다.
- 새 GPU 런이 필요 없다. **다음에 어차피 돌릴 native 런에 플래그만 얹으면 된다.**
- 새 파라미터·임계값·게이트 없음. 세는 것뿐이다.

## 8. 이 제안이 요구하지 않는 것

- 술어 변경 없음. conjunction은 **채택되지 않는다** — 세어질 뿐이다.
- gate 계약 변경 없음. `STATE.md` 행 추가 없음.
- 물리 상수·census·rung 변경 없음.
- PI 서명 불필요(계측은 gate 변경이 아니다). **단, 나온 숫자로 P1을 결정하는 것은 PI 몫이다.**

## 9. 소유 레인에 묻는 것

1. `transaction.py`의 선택 인자 두 개와 커널 — `observe-infra`가 받을 만한 형태인가?
2. `driver.py`가 `inner_d.converged_d`를 넘기는 것 — `cortex` 레인에서 사소한가?
3. 다음 native 런에 얹을 수 있는가, 아니면 전용 런이 필요한가?

셋 중 하나라도 아니라면 어디가 걸리는지 알려주면 제안을 그쪽에 맞춰 다시 쓴다. 이 레인은
`ac/` 아래를 편집하지 않으므로 **구현은 소유 레인의 것**이다.
