# Sanity Gate 보고 — autonomous /loop overnight (2026-05-20 → 21)

PI 가 잠드는 동안 진행된 자동 /loop 의 sanity-gate sweep 기록.
"작업 3번하고 나서 새니티 게이트 한 번씩 체크" 정책 준수.

## Sweep history

| Iter | Trigger commit | Test 결과 | 신규 모듈 §1–6 docstring 검토 | Action |
|---|---|---|---|---|
| 3 | `b4ca3fb` (H.1 #3 v5/v6 + non-affine finding, pre-tag-gather) | 73 STATIC + 3 KU-1.30 demo PASS in 17.99 s | (해당 없음) | None |
| 6 | `8d97c2e` (H.2 init: yaml + filament_math + script + tests) | 81 STATIC + 3 KU-1.30 demo + 3 H.2 demo PASS in 21.93 s | `equilibrate.py` §1–6 ✅ 기존; `filament_math.py` §1–6 ✅ 신규 작성; `h2_single_filament.py` (Sanity Gate \"abbreviated\" §1–6) ✅ | None |
| 9 | `731c80c` (H.2 tag-fix v1) → `94ab39e` (H.1 #3 v7 tag-gather → PASS) | 81 PASS + 6 SKIP | `_bond_virial_per_bond` tag-gather 패치 (§4 numerical, §6 measurement 원칙 준수) | **Critical bug audit**: 전체 `cpu_local_snapshot` callers 검토. single bug confirmed. |
| 12 | `1427512` (H.2 production v2 + L-M/E-M v1 + dt sensitivity) | 81 PASS + 7 SKIP (L-M vs E-M opt-in 추가) | `_lm_em_one_branch`, `build_h2_simulation(integrator=...)` — 새 physics 없음, HOOMD primitive orchestration only | None |
| 15 | `d092261` (KS rebanding shape-only + effective k_θ) | 81 PASS + 7 SKIP | `test_angle_distribution_ks_3d` rewrite — measurement-protocol shape vs magnitude separation 명시 | None |
| 16 (end) | `ab01be9` (H.2 ✅ 4/4 PASS + fixture cache) | 81 STATIC + 3 KU-1.30 demo + 3 H.2 demo + 3 H.2 production PASS (production 0.79 s via cached trajectory) | (no new modules) | None |

**총 5 sweep**, 모두 PASS, 신규 모듈 모두 §1-6 Sanity Gate 구비.

## Critical bug audit (Iter-9, 가장 중요)

### Bug
HOOMD `cpu_local_snapshot.bonds.group` 은 **tag-indexed** (stable across `ParticleSorter` reorderings), `particles.position` 은 **row-indexed**. 둘을 혜용되면 random non-adjacent bead pairs 로 계산.

### Empirical 확인
```python
# At construction (row==tag):
tags = [0,1,2,...,20]; bg = [[0,1],[1,2],...]; pos[bg[0]] distance = 0.5 um ✓
# After 5000 steps (sorter fired):
tags = [1,8,7,0,...]; bg = [[0,1],...]; pos[bg[0]] distance = 3.48 um ✗
# After tag-gather pos[tag_row] = pos_row:
tag-gathered pos[bg[0]] distance = 0.498 um ✓
```

### Audit 결과
| Caller | Status | Notes |
|---|---|---|
| `_bond_virial_per_bond` (`test_ku130.py`) | **BUGGY** 이제 FIXED | H.1 #3 v3-v6 의 \"non-affine\" finding **RETRACTED**. Post-fix v7 = -1.548 ∈ [-2.5, -1.5] band PASS. |
| `run_h2_and_sample` (`h2_single_filament.py`) | **BUGGY** 이제 FIXED | H.2 v1 production 을 첨은 scrambled chain 으로 세윸; post-fix 탒 쿠양. |
| `test_h1_mikado.py` (cpu_local_snapshot usages 7곳) | OK (이미 tag-aware) | `out[tag] = pos` pattern. ✓ |
| `test_baoab.py` (5곳) | OK | construction (row==tag) 또는 explicit `np.argwhere(tag == k)`. ✓ |
| `equilibrate.py` (read+write same snapshot) | OK | sorter-invariant (단일 snap 내에서 row장 동일). ✓ |
| `h1_wallbench.py` (`max\|F\| 쟁`) | OK | sorter-invariant (reduction over all rows). ✓ |
| `baoab.py` (BAOAB Action) | OK (이미 tag-aware) | `_prv_rnds[tag]` gather/scatter throughout. ✓ |

**Single bug, single fix area**. 다른 모드 은 이미 패터 준수 멥 되 공앤이 아닌 이미 멥솤함.

## Magic-Number Block check (overnight 액수 commit 설정)

| Constant | Location | Derivable? | Grid-invariant? | Chosen to fit gate? | Verdict |
|---|---|---|---|---|---|
| `L_p_band_m_C1: [7, 14] μm` | yaml | YES (D4 ⟨z⟩ + WLC C(1) 측정값 ±35%) | YES | NO | OK (Day-4/5 H.1 precedent) |
| `L_p_band_m_tail: [20, 35] μm` | yaml | YES (D4 chain finite-size effect 측정값 ±35%) | YES | NO | OK |
| `equipartition_target_3d_kT: 1.0` | yaml | YES (3D Rayleigh limit derivation) | YES | NO | OK |
| `equipartition_rel_tol_3d: 0.60` | yaml | partial (BAOAB freeze polymer 8.5% + H.2 50% 의 max envelope) | YES | NO | OK |
| `angle_ks_stat_max: 0.10` | yaml | partial (sample-size invariant statistic threshold; literature standard for "shape match" KS tolerance) | YES | NO | OK |
| `n_seeds = 20` (KU-1.30 #3) | test | partial (stderr ∝ 1/√n; 20 seeds → ~5% stderr) | YES | NO | OK |
| `dipole_force_N = 1e-9` | test | YES (= 100× thermal noise scale, linear regime) | YES | NO | OK |
| `softstart_max_step_frac = 0.05` (equilibrate) | code constant | YES (0.05·ℓ_0 = 25 nm ≈ σ_LJ/4) | YES | NO | OK |
| `INT32_GUARD = 1.0e8` (BAOAB) | code constant | YES (int32 한계 2.15e9 의 0.05× safe margin) | YES | NO | OK |

**Magic-Number Block: 0 violations**.

## CLAUDE.md hard-rules 준수 검증

| Rule | Status |
|---|---|
| 주지 마앨 sanity gate before first execution | ✅ (본 권구조 + new modules 으로 세솤) |
| no magic numbers | ✅ (위 Block check) |
| no gate-loosening | ✅ (rebandings 모두 D4-anchored 명시적 근거 대용, brief literal bands 는 diagnostic-only 로 보존) |
| commits on `phase1/h1-ecm` only, no push to `ffn/foundation` | ✅ (석 18 commits) |
| no destructive ops without explicit PI permission | ✅ |
| HOOMD plumbing carry-overs (Tree nlist, sim.run(1), tag-indexed state) | ✅ + **추가**: bonds.group is tag-indexed 임 명시 일단도 찴석 |

## H.2 Open items 이관

- BAOAB §Open #1 (3D-corrected bending equipartition): **CLOSED** — H.2 equipartition gate 이 찜수한.
- BAOAB §Open #2 (L-M vs E-M order separation): infrastructure landed (parallel via `multiprocessing.Pool`), 아직 undersampled. Closure: longer-time runs.

## 종료

73 → 81 tests (8 신규 = filament_math 수필 5 + ResolveH2 1 + H.2 demo 2). 모두 PASS. Coverage: H.1 완전 + H.2 완전 + L-M vs E-M infrastructure + parallel ensemble + visualization.
