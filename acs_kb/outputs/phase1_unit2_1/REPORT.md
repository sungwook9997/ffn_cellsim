# Phase 1 — Unit 2.1 Progress Report (post-PI-review v2)

**Date**: 2026-05-18 (Worker B, Bridge Track) — revised after PI conditional-acceptance feedback the same day.
**Scope**: Bridge Unit 2.1 — single-FA Chan-Odde motor-clutch on a linear-elastic substrate stub with Pereverzev catch-slip bonds. No FA growth, no talin unfolding, no vinculin recruitment, no ECM fibre network — substrate stub only (per Brief).
**Repo HEAD**: branch `worker-b/bridge`, off `main @ 4ffcbbe`.
**Notion source**: KU-2.1, KU-2.4, KU-2.5, KU-2.8, KU-2.12, KU-2.18, KU-1.21.
**Brief**: 🚀 Claude Code Brief — Phase 1 Unit 2.1 (Worker B), 364120da-ec5d-8155-8a52-f3a98a4d15b6.

## PI conditional-acceptance fixes (this revision)

Four contract surfaces tightened per PI feedback, reflected in the current code/config/tests and in this report:

1. **KU-2.12 reclassified**. Phase 1 Unit 2.1 acceptance moves to the **nascent-FA 1–8 pN band** with a mass threshold derived from the expected exponential-tail distribution (mean ≈ 3.6 pN ⇒ `exp(-1/3.6) − exp(-8/3.6) = 0.65`, so the threshold 0.55 leaves ~10 pp stochastic margin). The KU-2.12 **mature-FA 5–20 pN band moves to the Unit 2.2 contract** where vinculin reinforcement of `k_int^eff` lets the per-clutch regime cross into the 5–20 pN window. Both bands are tracked in the Sanity Gate output; only the nascent band gates Unit 2.1.
2. **`from_config` everywhere**. `MotorClutchParams`, `CatchSlipParams`, and `LinearElasticSubstrate` each expose a `from_config(cfg)` classmethod that reads the resolved YAML. The notebook + tests no longer construct these dataclasses with class defaults — the YAML is the single source of truth. A new `test_params_from_config_round_trip` enforces this contract.
3. **Biphasic claim rephrased**. The original single-seed sweep found `argmax_E ≈ 4.3 kPa`. The PI requested replicated-seed validation. We ran 5 seeds × 12 stiffnesses × 5 s sim. Result: per-seed argmax_E spans **2.27 decades** (`[534, 4329, 12328, 12328, 100000]` Pa) and the peak prominence above the asymptote is **0.34 %**. The curve **saturates** at ≈ 97 pN ≈ 1.0 · `N_motors · F_stall`. The regression test now asserts `verdict == "saturating"`, an asymptote within `[0.7, 1.1] · N_m · F_stall`, and *no longer* asserts argmax-E vs analytic. The matched-stiffness analytic (Bangasser 2013) is logged as info but not gated.
4. **`contact_radius` clarified**. Renamed neither field; instead the docstrings and the YAML now explicitly call out that `substrate.contact_radius` (Boussinesq-form compliance length) and `motor_clutch.fa_area` are **independent model parameters**. The geometric projection `a = √(A/π) = 0.564 μm` is documented as an alternative that would shift `k_sub` by ×1.77 — Phase 1 stays with the Brief Task 2 specification `a = 1 μm`. Unit 2.2's ECM adapter obviates this ambiguity by computing `k_sub` directly from the fibre network.

## What is built

```
acs_kb/bridge/
├── __init__.py
├── types.py             # FocalAdhesion (week-5 freeze candidate) + make_focal_adhesion()
├── substrate_stub.py    # LinearElasticSubstrate (KU-1.21) with from_config
├── catch_bond.py        # Pereverzev k_off(F), τ(F), F* analytic + from_config
├── motor_clutch.py      # MotorClutchFA stepper + run_steady_state + from_config
└── traction.py          # 1-D engaged-clutch reducer (KU-2.12)
acs_kb/common/derived_params.py    # resolve_bridge() + load_bridge_config()
acs_kb/common/sanity_gate.py       # gate_unit2_1_motor_clutch() (verdict-aware)
acs_kb/configs/phase1_unit2_1.yaml # primary scales only; everything else derived
acs_kb/tests/test_bridge.py        # 16 unit tests (incl. from_config round-trip)
acs_kb/tests/test_KU_2_4_biphasic.py  # 1 multi-seed slow test (saturation verdict)
acs_kb/notebooks/02_motor_clutch_biphasic.py  # 4 figures + summary
acs_kb/outputs/phase1_unit2_1/
    ├── 01_catch_bond_lifetime.png        # τ(F) analytic vs Gillespie samples
    ├── 02_biphasic_traction.png          # per-seed + mean ⟨F⟩(E) — saturating
    ├── 03_clutch_force_histogram.png     # per-clutch force at E=5 kPa, two-band view
    ├── 04_biphasic_argmax_spread.png     # per-seed log10(argmax_E), 2.3-dec spread
    ├── 02_summary.json                   # machine-readable numbers
    └── QUESTIONS_FOR_SUNGWOOK.md         # 1 open item after PI review
```

## KU validation numbers

| Quantity | Measured | Analytic | Acceptance | Status |
|---|---|---|---|---|
| **Catch-slip F\*** (Pereverzev w/ KU-2.18 defaults) | **6.99 pN** (argmax τ(F) grid) | 6.99 pN | within ±20 % of analytic | PASS (rel err < 1e-3) |
| **Gillespie ⟨τ⟩(F=F\*)** | within 5 % of 1/k_off(F) for N=10⁴ samples | — | unit test | PASS |
| **Biphasic shape (5 seeds × 12 E × 5 s sim)** | verdict = **saturating** (prominence 0.34 %, spread 2.27 dec) | — | not a sharp peak | PASS (honest) |
| **Asymptote ⟨F⟩ vs N\_m · F\_stall** | **97 pN** / 100 pN ratio 0.970 | 100 pN | [0.7, 1.1] | PASS |
| **Per-clutch force band on E=5 kPa** | mean 3.57 pN, median 2.38 pN, **63.6 %** in [1, 8] pN | exp-tail ⇒ ~65 % | ≥ 55 % in [1, 8] pN | PASS |
| **Per-clutch force in KU-2.12 mature 5-20 pN band** | 27 % | (Unit 2.2 target) | info only | logged |
| **dt < α / max(k\_on, k\_off(F\_s))** | 1·10⁻⁴ s | 7.3·10⁻² s | strict | PASS (730× under) |
| **Perf 1000 × 50 clutches** | **21 ms** | — | < 1 s | PASS (47× under) |

## Test suite

```
acs_kb/tests/test_bridge.py — 16 passed in 0.12 s
  ✓ test_catch_bond_offrate_zero_force                          [KU-2.5]
  ✓ test_catch_bond_offrate_vectorizes                          [API]
  ✓ test_catch_bond_peak_matches_analytic                       [KU-2.5]
  ✓ test_catch_bond_gillespie_lifetime_matches_inverse_offrate  [KU-2.5]
  ✓ test_catch_bond_slip_only_returns_nan                       [KU-2.5 pitfall]
  ✓ test_substrate_stub_linear_compliance                       [KU-1.21]
  ✓ test_substrate_stub_compute_displacement_vectorizes         [API]
  ✓ test_motor_clutch_force_balance                             [KU-2.4]
  ✓ test_motor_clutch_no_engaged_gives_zero                     [boundary]
  ✓ test_motor_clutch_step_deterministic_under_seed             [RNG]
  ✓ test_motor_clutch_step_advances_actin                       [KU-2.4]
  ✓ test_traction_sums_engaged_clutches                         [KU-2.12]
  ✓ test_traction_zero_when_none_engaged                        [boundary]
  ✓ test_motor_clutch_performance_budget                        [Brief Task 7]
  ✓ test_bridge_config_resolves                                 [derived]
  ✓ test_params_from_config_round_trip                          [PI fix #2]

acs_kb/tests/test_KU_2_4_biphasic.py — 1 passed (mark: slow) in 80 s
  ✓ test_KU_2_4_biphasic_is_saturating_not_peaked               [KU-2.4 / KU-2.8]
```

## Sanity Gate output (post-fix)

```
Sanity Gate · unit2_1_motor_clutch
  PASS dt_below_bond_event_rate_cfl              [KU-2.4 / KU-2.5]
       dt=1.000e-04 s, cfl_safe_dt=7.328e-02 s (α·1/max(k_on, k_off(F_s)), α=0.1)
  PASS catch_peak_F_star_within_analytic         [KU-2.5]
       F*_sim=6.99 pN, F*_analytic=6.99 pN, rel err=0.000, tol=0.2
  PASS catch_peak_KU_2_5_experimental_gap_logged [info]
       analytic F*=6.99 pN; KU-2.5 experimental F*≈30 pN; gap=-23.01 pN.
  PASS biphasic_shape_is_saturating              [KU-2.4 / KU-2.8]
       verdict=saturating; prominence=0.0034; inter-seed log10 argmax spread=2.27 dec
  PASS biphasic_asymptote_near_motor_stall       [KU-2.4]
       asymptote ⟨F⟩=96.96 pN, N_m·F_stall=100.00 pN, ratio=0.970
  PASS biphasic_matched_stiffness_E_star_info    [info]
       matched-stiffness E*=6.35e+03 Pa (Bangasser 2013) vs argmax E=1.23e+04 Pa.
  PASS per_clutch_force_in_phase1_nascent_band   [KU-2.12 Phase 1 (nascent)]
       mass in [1.0, 8.0] pN = 0.636, required ≥ 0.55 (KU-2.12 mature 5-20 pN → Unit 2.2)
  PASS perf_budget_met                           [Brief Task 7]
       1000 steps × 50 clutches in 21.0 ms (budget 1000 ms)
```

**8 / 8 PASS** (previous revision was 7 / 8 with KU-2.12 band scope-mismatch fail).

## Biphasic shape: replicated-seed evidence

| Seed root | argmax E (Pa) |
|---|---|
| 1234        | 4 329  |
| 2234        | 12 328 |
| 3234        | 534    |
| 4234        | 100 000 (sweep boundary) |
| 5234        | 12 328 |

- **argmax_E spread**: 2.27 decades — the "peak" is noise-dominated.
- **Peak prominence above asymptote**: 0.0034 (0.34 %).
- **Asymptote** `⟨F⟩` at high E: 96.96 pN ≈ 0.970 · `N_motors · F_stall = 100 pN`.
- **Saturation knee** (E at which ⟨F⟩ first reaches 90 % of asymptote): 187 Pa.

The per-clutch force at force-balance with `k_int = 1 pN/nm`, `N_m·F_stall = 100 pN`, ~25 engaged clutches gives `F_i ≈ 4 pN`, well below the slip-pathway force `F_s = 30 pN`. With no clutch entering the slip regime, the system has **no mechanism to dump force at stiff substrates** → the curve saturates rather than declining. Unit 2.2's vinculin allostery (`k_int^eff = k_int^bare (1 + α N_vin)`, KU-2.7) is the mechanism that should recover a true peak.

## KU-2.5 catch-peak F\* discrepancy (still open — unchanged from v1)

The Pereverzev closed form with KU-2.18 defaults gives `F* = 6.99 pN`, while KU-2.5 quotes the experimental Kong 2009 α5β1 value at `F* ≈ 30 pN`. The two are **not simultaneously satisfiable** with the KU-2.18 parameter set — recovering 30 pN needs either `k_c ≈ 20 s⁻¹` or substantially different `F_s/F_c` ratios. See `QUESTIONS_FOR_SUNGWOOK.md` Q1 for the decision tree; refit deferred to Phase 2.

## Interface freeze candidate (week-5 contract)

`acs_kb.bridge.types.FocalAdhesion` is the data structure Worker C will import to populate `Cell.focal_adhesions`. The fields are considered **frozen for week 5**; any rename / type change before week 5 requires a Notion alert to Worker C.

```python
@dataclass(slots=True)
class FocalAdhesion:
    position: np.ndarray              # (2,) m
    n_clutches_total: int
    clutches_engaged: np.ndarray      # bool (N,)
    clutch_forces: np.ndarray         # (N,) N
    actin_position: float = 0.0       # m, 1-D along loading axis
    area: float = 1e-12               # m² (1 μm² default, INDEPENDENT of substrate.contact_radius)
    vinculin_count: int = 0           # populated from Unit 2.2
    talin_unfolded_domains: int = 0   # populated from Unit 2.2
    age: float = 0.0                  # s
```

Per-clutch substrate anchors and the common substrate displacement live on `MotorClutchFA` (simulator state), not on `FocalAdhesion` (contract).

## Performance

For a single FA with 50 clutches, dt = 1·10⁻⁴ s:

| Workload | Wall time |
|---|---|
| 1000 steps                                        | **21 ms** |
| 50 000 steps (5 s sim time, single E)             | ~1 s   |
| Replicated biphasic sweep (5 seeds × 12 E × 50 k) | ~80 s  |
| Histogram run (50 000 steps × snapshot/100)       | ~1 s   |

## Out of scope (per Brief stop point)

- FA growth (Hill function, KU-2.17) — Unit 2.2.
- Talin unfolding (Bell-Evans, KU-2.6) — Unit 2.2.
- Vinculin recruitment (KU-2.7) — Unit 2.2 (mechanism for entering the KU-2.12 mature force band and recovering a true biphasic peak).
- ECM fibre network — substrate stub only; ECM adapter is the first task of the Unit 2.2 brief.
- 2-D traction direction — `compute_traction` already returns `(2,)` ndarray.
- Pereverzev refit to Kong 2009 — deferred to Phase 2 (QUESTIONS Q1).

## Suggested next step — Unit 2.2

1. Implement Hill-function FA growth (KU-2.17, Walcott-Sun 2010 PNAS).
2. Talin 1-state Bell-Evans unfolding (KU-2.6, del Rio 2009 Science).
3. Vinculin recruitment with allosteric `k_int^eff` increase (KU-2.7, Franz 2023 Nat Commun).
4. Replace `substrate_stub.LinearElasticSubstrate` with `bridge.ecm_adapter.ECMAdapter` consuming Worker A's frozen ECM API.
5. Inherit the KU-2.12 **mature** 5–20 pN per-clutch force band as the Unit 2.2 acceptance; expect a true (peaked, not saturating) biphasic curve once vinculin reinforcement lets clutches enter the slip regime.
6. Refit Pereverzev parameters against Kong 2009 force-clamp data (QUESTIONS Q1).

**Status**: this revision removes the earlier overclaim around biphasic peak validation and records Unit 2.1 as a Phase 1 bridge scaffold. Q1 remains open before claiming literature-level catch-bond calibration.
