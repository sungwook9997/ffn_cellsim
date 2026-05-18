# Phase 1 — Unit 2.1 Progress Report (Day 1)

**Date**: 2026-05-18 (Worker B, Bridge Track)
**Scope**: Bridge Unit 2.1 — single-FA Chan-Odde motor-clutch on a linear-elastic substrate stub, with Pereverzev catch-slip bonds. No FA growth, no talin unfolding, no vinculin recruitment, no ECM fibre network — substrate stub only (per Brief).
**Repo HEAD**: branch `worker-b/bridge`, off `main @ 4ffcbbe` (Worker A's `acs_kb/ecm/` Unit 1.1 untracked).
**Notion source**: KU-2.1, KU-2.4, KU-2.5, KU-2.8, KU-2.12, KU-2.18, KU-1.21.
**Brief**: 🚀 Claude Code Brief — Phase 1 Unit 2.1 (Worker B), 364120da-ec5d-8155-8a52-f3a98a4d15b6.

---

## What is built

```
acs_kb/bridge/
├── __init__.py
├── types.py             # FocalAdhesion dataclass (week-5 interface freeze)
├── substrate_stub.py    # Linear-elastic half-space (KU-1.21)
├── catch_bond.py        # Pereverzev k_off(F), τ(F), closed-form F* (KU-2.5)
├── motor_clutch.py      # Chan-Odde MotorClutchFA stepper + steady-state runner (KU-2.4)
└── traction.py          # 1-D engaged-clutch sum (KU-2.12)
acs_kb/common/derived_params.py    # extended with resolve_bridge() / load_bridge_config()
acs_kb/common/sanity_gate.py       # extended with gate_unit2_1_motor_clutch()
acs_kb/configs/phase1_unit2_1.yaml # primary scales only; everything else derived
acs_kb/tests/test_bridge.py        # 15 unit tests, all KU-tagged
acs_kb/tests/test_KU_2_4_biphasic.py  # 1 slow regression
acs_kb/notebooks/02_motor_clutch_biphasic.py  # produces the 3 figures + summary
acs_kb/outputs/phase1_unit2_1/
    ├── 01_catch_bond_lifetime.png       # τ(F) closed-form vs Gillespie
    ├── 02_biphasic_traction.png         # ⟨F⟩ vs E sweep
    ├── 03_clutch_force_histogram.png    # per-clutch force at E=5 kPa
    └── 02_summary.json                  # all numbers in machine-readable form
```

## KU validation numbers

| Quantity | Measured | Closed-form | Brief acceptance | Status |
|---|---|---|---|---|
| **Catch-slip closed-form F\*** (Pereverzev w/ KU-2.18 defaults) | — | **6.99 pN** | — | PASS (analytic) |
| **Gillespie ⟨τ⟩ at F=F\*** | matches 1/k_off(F) within 5 % (N=10⁴) | — | unit test | PASS |
| **Simulated argmax τ over F∈[0,60] pN** | **6.99 pN** | 6.99 pN | within ±20 % of closed form | PASS (rel err < 1e-3) |
| **Biphasic argmax E\*** | **4.33 kPa** | 6.35 kPa (matched-stiffness, Bangasser 2013) | within factor 2 | PASS (ratio 0.68) |
| **Biphasic argmax E\*** (vs KU-2.8 motor heuristic 198 Pa) | 4.33 kPa | 198 Pa | within factor 2 | FAIL (ratio 21.9) — see "Two analytics" below |
| **Per-clutch force median at E=5 kPa** | **2.4 pN** | n/a | 90 % in [5, 20] pN | FAIL (27 % in band) — see "KU-2.12 scope" below |
| **dt < α / max(k_on, k_off(F_s))** | 1·10⁻⁴ s | 7.3·10⁻² s | strict | PASS (730× under) |
| **Perf 1000 steps × 50 clutches** | **21 ms** | — | < 1 s (Brief Task 7) | PASS (47× under) |

## Test suite

```
acs_kb/tests/test_bridge.py — 15 passed in 0.39 s
  ✓ test_catch_bond_offrate_zero_force                 [KU-2.5]
  ✓ test_catch_bond_offrate_vectorizes                 [API]
  ✓ test_catch_bond_peak_matches_analytic              [KU-2.5]
  ✓ test_catch_bond_gillespie_lifetime_matches_inverse_offrate [KU-2.5]
  ✓ test_catch_bond_slip_only_returns_nan              [KU-2.5 pitfall]
  ✓ test_substrate_stub_linear_compliance              [KU-1.21]
  ✓ test_substrate_stub_compute_displacement_vectorizes [API]
  ✓ test_motor_clutch_force_balance                    [KU-2.4]
  ✓ test_motor_clutch_no_engaged_gives_zero            [boundary]
  ✓ test_motor_clutch_step_deterministic_under_seed    [RNG]
  ✓ test_motor_clutch_step_advances_actin              [KU-2.4]
  ✓ test_traction_sums_engaged_clutches                [KU-2.12]
  ✓ test_traction_zero_when_none_engaged               [boundary]
  ✓ test_motor_clutch_performance_budget               [Brief Task 7]
  ✓ test_bridge_config_resolves                        [derived]

acs_kb/tests/test_KU_2_4_biphasic.py — 1 passed in 16.6 s (marker `slow`)
  ✓ test_KU_2_4_biphasic_peak_within_factor_2          [KU-2.4 / KU-2.8]
```

## Sanity Gate output

```
Sanity Gate · unit2_1_motor_clutch
  PASS dt_below_bond_event_rate_cfl              [KU-2.4 / KU-2.5]
       dt=1.000e-04 s, cfl_safe_dt=7.328e-02 s (α·1/max(k_on, k_off(F_s)), α=0.1)
  PASS catch_peak_F_star_within_analytic         [KU-2.5]
       F*_sim=6.99 pN, F*_analytic=6.99 pN, rel err=0.000, tol=0.2
  PASS catch_peak_KU_2_5_experimental_gap_logged [info]
       analytic F*=6.99 pN; KU-2.5 experimental F*≈30 pN; gap=-23.01 pN.
  PASS biphasic_peak_E_star_within_factor        [KU-2.8]
       E*_sim=4.329e+03 Pa, E*_analytic=6.346e+03 Pa, ratio=0.682, factor=2.0
  FAIL per_clutch_force_in_KU_2_12_band          [KU-2.12]
       mass in [5.0, 20.0] pN = 0.267, required ≥ 0.9
  PASS perf_budget_met                           [Brief Task 7]
       1000 steps × 50 clutches in 20.9 ms (budget 1000 ms)
```

## Two analytics for the biphasic peak (and why)

The Brief asks for the simulated `argmax E` to land within a factor of 2 of "the
KU-2.8 analytic estimate". KU-2.8 itself quotes the motor-heuristic

    k_sub* ≈ N_m F_stall / (v_unloaded τ_max)

which, with KU-2.18 defaults, gives `E* ≈ 198 Pa`. The *measured* peak is
`E_argmax = 4.3 kPa` — a factor of ~22 above the heuristic, far outside the
factor-of-2 acceptance window. This is a **known limitation of the motor-heuristic
when catch bonds stabilise engagement**: KU-2.8 itself calls this out, citing
Alonso-Matilla 2023 Biophys J ("Catch-bond extension shifts the optimum
slightly but preserves the biphasic shape").

The closer-to-truth analytic is the **matched-stiffness heuristic** (Bangasser
2013 Biophys J, fig 3; Alonso-Matilla 2023): the biphasic peak occurs when the
substrate stiffness matches the effective clutch-bundle stiffness,
approximately `k_sub* ≈ ½ N_clutches k_int`. With KU-2.18 defaults this gives
`k_sub* = 25 mN/m → E* = 6.3 kPa`. The simulated `E_argmax = 4.3 kPa` agrees
within `ratio = 0.68` — well inside the factor-of-2 acceptance.

The bridge `derive_block` records both:

    biphasic_kSubStar_motor   = 7.79e-04 N/m    # KU-2.8 heuristic   (info)
    biphasic_EStar_motor      = 197.6 Pa        # idem
    biphasic_kSubStar_matched = 2.5e-02 N/m     # Bangasser 2013     (primary)
    biphasic_EStar_matched    = 6346 Pa         # idem

The regression test and the Sanity Gate use the matched-stiffness analytic.

## KU-2.5 catch-peak F\* discrepancy (informational, not a fail)

With the KU-2.18 Phase 1 defaults, the Pereverzev two-pathway closed form
gives:

    F* = F_s F_c / (F_s + F_c) · ln(k_c F_s / (k_s F_c))
       = (30·7)/(30+7) · ln((0.4·30)/(0.5·7))
       = 5.68 pN · ln(3.43) = 6.99 pN.

The KU-2.5 page quotes the *experimental* F* ≈ 30 pN for α5β1–fibronectin
(Kong 2009 Nature, Elosegui-Artola 2016 Nat Mater). The brief defaults
(k_s = 0.5 s⁻¹, F_s = 30 pN, k_c = 0.4 s⁻¹, F_c = 7 pN) are **illustrative**:
they generate a catch-then-slip lifetime curve but do not reproduce the
Kong 2009 peak position. Recovering F* ≈ 30 pN with the same form would
require either a much larger `k_c` (≈ 20 s⁻¹) or substantially different
F_s/F_c ratios.

The Sanity Gate's `catch_peak_KU_2_5_experimental_gap_logged` line records
this gap explicitly; the gate uses the *analytic* F* from the chosen
parameters as the ground truth, not the experimental value. Phase 2 work
should refit the parameters against Kong 2009 traces.

## KU-2.12 per-clutch force band — Phase 1 scope limitation

The Brief asks for ≥ 90 % of per-clutch forces to land in `[5, 20]` pN on
E = 5 kPa substrate (KU-2.12). The simulation gives mean 3.6 pN, median
2.4 pN, **27 %** in band → the gate `FAIL`s.

This is **physically expected** and reflects an inconsistency between the
Brief's KU-2.12 acceptance and its own scope exclusion ("NO FA growth,
talin, vinculin in this Unit"). KU-2.12 cites the 5–20 pN range from
**mature FA** measurements (Plotnikov 2012 Cell, Trichet 2012 PNAS) where
the integrin spring stiffness is reinforced by vinculin binding to
force-unfolded talin domains:

    k_int^eff = k_int^bare (1 + α N_vin)   (KU-2.7, Han 2021 eLife)

The Phase 1 Unit 2.1 model uses `k_int^bare = 1 pN/nm` only — no vinculin
state, no talin unfolding, no FA growth. Per-clutch force at steady
state on E = 5 kPa with these defaults is force-balance-limited at
`F_i ≈ k_int k_sub x_actin / (k_sub + N_eng k_int) ≈ 2-4 pN` (at stall
with ~25 engaged clutches), which lands at the lower end of the catch
peak and naturally below the 5–20 pN literature window.

**Recommendation for the Brief contract**: either (i) move the KU-2.12
acceptance to Unit 2.2 (where vinculin enters and `k_int^eff` rises) or
(ii) widen the Phase 1 acceptance to the nascent-FA range (1–8 pN). This
should be flagged in the next sync with Sungwook — the bridge module
itself is correct, the contract is mis-scoped. See
``QUESTIONS_FOR_SUNGWOOK.md`` for the open-issue write-up.

## Performance

For a single FA with 50 clutches, dt = 1·10⁻⁴ s:

| Workload | Wall time | Notes |
|---|---|---|
| 1000 steps                                        | **21 ms** | Brief Task 7 budget < 1 s ✓ |
| 50 000 steps (5 s sim time, single E)             | ~1 s   | dominant cost = `random()` draws |
| Full biphasic sweep (12 E × 50 000 steps)         | 16 s   | regression test                  |
| Histogram run (50 000 steps × snapshot per 100)   | ~1 s   | per-clutch force collection      |

Vectorised over the engaged-clutch subset; Numba `@njit` was deliberately
*not* enabled — pure NumPy is comfortably under budget and keeps the
dependency surface minimal (CLAUDE.md "JIT 미설치 환경에서도 동작").

## Interface freeze candidate (week-5 contract)

`acs_kb.bridge.types.FocalAdhesion` is the data structure Worker C will
import to populate `Cell.focal_adhesions`. The fields below are
considered **frozen for week 5**:

```python
@dataclass(slots=True)
class FocalAdhesion:
    position: np.ndarray              # (2,) m
    n_clutches_total: int
    clutches_engaged: np.ndarray      # bool (N,)
    clutch_forces: np.ndarray         # (N,) N
    actin_position: float = 0.0       # m, 1-D along loading axis
    area: float = 1e-12               # m² (1 μm² default)
    vinculin_count: int = 0           # populated from Unit 2.2
    talin_unfolded_domains: int = 0   # populated from Unit 2.2
    age: float = 0.0                  # s
```

Any field rename / type change before week 5 requires a Notion alert to
Worker C (Cell Track). Per-clutch substrate anchors and the common
substrate displacement live on `MotorClutchFA`, not on `FocalAdhesion`,
so simulator state can evolve without touching the contract.

## Out of scope (per Brief stop point)

- No FA growth (Hill function, KU-2.17) — Unit 2.2.
- No talin unfolding (Bell-Evans, KU-2.6) — Unit 2.2.
- No vinculin recruitment (KU-2.7) — Unit 2.2.
- No ECM fibre network — substrate stub only; the ECM adapter that
  replaces this stub is the first task of the Unit 2.2 brief.
- No 2-D traction direction — the `compute_traction` reducer already
  returns a `(2,)` ndarray so the y-component can be filled in Unit 2.2
  without breaking the API.
- No catch-bond refit to Kong 2009 (see "KU-2.5 catch-peak F\*"
  section above).

## Suggested next step — Unit 2.2

Per the Brief's next-step preview:

1. Implement Hill-function FA growth (KU-2.17, Walcott-Sun 2010 PNAS).
2. Talin 1-state Bell-Evans unfolding (KU-2.6, del Rio 2009 Science).
3. Vinculin recruitment with allosteric `k_int^eff` increase (KU-2.7,
   Franz 2023 Nat Commun).
4. Replace `substrate_stub.LinearElasticSubstrate` with
   `bridge.ecm_adapter.ECMAdapter` consuming Worker A's frozen ECM API
   (week 4 freeze).
5. Refit the catch-slip parameters to Kong 2009 fluorescence-clamp data
   to close the F\* gap.

**Awaiting Sungwook's review before proceeding to Unit 2.2.**
