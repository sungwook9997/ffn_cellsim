# Stage 1a+ — bounded outcomes (Option α: mechanical anchor only)

This document records the **PI-approved bounded-outcome decision tree** for
the Stage 1a+ Option α pilot (γ_sub = 0, mechanical anchor only). It is
written *before* the pilot runs so the result is read against a contract
that exists prior to any data.

The companion documents are:
- `docs/stage1a_plus_substrate_sanity.md` — six-check Sanity Gate +
  Magic-Number Block + Option α resolution record
- `docs/outcomes_v15.md` — v15 baseline (R drift 24.4%, single-layer limit)
- `docs/12_validation.md` — Sanity-Gate Protocol + Magic-Number Block

## Mechanism question (one)

**Does mechanical-only substrate anchoring — a rigid reflective floor at
`z = 0` with no wetting-energy term — materially reduce the v15 single-layer
baseline R drift of 24.4%?**

The pilot is a **single simulation** producing a **single number**
(`R_drift_1aplus`), evaluated against the v15 baseline (`R_drift_v15 =
0.244`). All other gates are pass/fail correctness contracts (mass,
momentum-horizontal, energy-monotone, no-NaN, anchor force balance,
contact-band ρ_kernel) and play no role in the bounded-outcome bucketing.

## Bounded outcomes (4 buckets)

The buckets are defined on `R_drift_1aplus = max_t |R(t)/R₀ − 1|` over the
same `t ≥ t_check` window used by v15 (`radius_drift_check_after_s = 600 s`,
i.e. t* ≥ 10).

| Outcome (R drift) | Interpretation | Action |
|---|---|---|
| **< 5%** | Mechanical anchor alone restores spheroid integrity. The v15 24% drift was almost entirely a free-floating-spheroid artifact; energetic wetting is a small correction. | STOP. Surface to PI. Stage 1a+ passes. Decisions to consider: (i) move to Stage 1a++ (Layer 2 boundary biology) directly, treating the energetic-wetting question as Stage 3 future work; (ii) optionally run Option β (Ca_sub sweep) as a refinement before Stage 1a++. PI selects. **No automatic mid-run.** |
| **5 – 15%** | Mechanical anchor contributes substantially but cannot fully balance Laplace. Energetic wetting (γ_sub > 0) is the natural next missing piece. | STOP. Surface to PI. Recommended next step: **Option β** — `Ca_sub = α·Ca_cc` sweep for α ∈ {0.1, 0.3, 1.0} (3 runs anchored to Maître Science 2012's Ca_cc), reporting R drift as a response curve. PI signs off before Option β runs. **No automatic mid-run.** |
| **15 – 22%** | Mechanical anchor adds little; the v15 single-layer limit is the dominant effect. | STOP. Surface to PI. Decisions: (i) revisit the v15 deferred (i)–(iii) diagnostics from `docs/outcomes_v15.md` (Maxwell deviatoric residual / G2P kernel-density bias / small-strain effective bulk modulus) — these now have a substrate-anchored configuration to test against; (ii) consider Option β at higher α; (iii) consider architectural options from `docs/outcomes_v15.md` Path A. PI selects. **No automatic mid-run.** |
| **≥ 22%** (drift unchanged or worse) | Substrate is contributing **nothing** in the chosen scheme — or worse, destabilising. | STOP. Surface to PI as a **scheme-level failure** of mechanical-only substrate anchoring. The relative-improvement gate ("R_drift_1aplus < R_drift_v15") fails outright. This is the architectural-review territory analogous to v15 Outcome 4. PI decides among architectural alternatives or scope-level reframing. **No automatic mid-run, no within-scheme retries.** |

**Why these buckets**: Stage 1a+ Option α is a deliberate one-knob test
(mechanical anchor on/off; no continuous parameter). The bucketing reflects
the *qualitative* possible outcomes of the mechanism question, not a
parameter sweep. Each bucket corresponds to a distinct downstream path.

## Gates that must independently pass (Pillar 1 correctness)

These are correctness gates, not part of the outcome bucketing. A
violation halts the run and surfaces to PI **before** the R-drift
bucketing is read.

- mass conservation (exact, limit 1e-10)
- momentum drift — *horizontal only* (vertical absorbed by substrate per
  the Cousin-Rule contract change recorded in
  `docs/stage1a_plus_substrate_sanity.md` §3)
- energy monotone (KE + strain energy; no surface-energy term added under
  γ_sub = 0)
- no NaN / Inf
- max-speed bound
- peak VRAM ≤ ceiling
- v15-inherited gates (mass, calibration `<J>=1`, curvature, CSF
  localisation): these stay live and continue to enforce the v15 Layer-1
  contract
- **anchor force balance** (new under Option α): substrate reaction at
  contact band ≈ ∫σ_vol_zz dA over the same band; tolerance ≤ 0.20
- **contact-band ρ_kernel / ρ_ref** (new under Option α): ∈ [0.85, 1.15]
  (witnesses kernel sees substrate via reflective BC, distinct from
  free-surface ~0.5·ρ_ref)

## Diagnostics (no gate, log-only under Option α)

- substrate contact area `A_contact` (no Young analytical reference under
  γ_sub = 0)
- apparent contact angle θ (no Young analytical reference; recorded for
  comparison with future Option β runs)
- shell-averaged `<ρ_kernel>(r/R₀)` profile (carried over from v15)

## Files of record (will be created by the pilot run)

```
results/stage1a_plus_pilot/gate_report.md         — Pillar-1 gate evaluation
results/stage1a_plus_pilot/metrics.csv            — per-frame invariants + shell summary + new diagnostics
results/stage1a_plus_pilot/shell_profile.csv      — long-format <ρ_kernel>(r/R₀)
results/stage1a_plus_pilot/contact_metrics.csv    — A_contact, θ_apparent, ρ_kernel_contact, anchor force per frame
results/stage1a_plus_pilot/snapshots.h5           — particle state per frame
results/stage1a_plus_pilot/curvature_validation.json
results/stage1a_plus_pilot/reference_calibration.json
results/stage1a_plus_pilot/run_manifest.json      — git hash + config + host
```

## Stop conditions (unchanged)

- Magic number 도입 금지
- Gate semantic 변경 금지 (Cousin Rule)
- v13 anti-pattern (measurement-protocol consistency 위반) 시 즉시 STOP
- 한도 외 escalation 금지 — pilot 종료 후 PI 결정 대기, 자동 다음 cycle 진입 금지
