# Session Handoff — 2026-04-29 (overnight)

This document records where the overnight Claude Code session stopped so the next session (or PI in person) can resume cleanly. Everything below has been committed to git in the `v10 — Brackbill smoothing partial fix, reference calibration drift discovered, STOP for PI analysis` commit.

## Timeline

| Block | What happened |
|---|---|
| **Stage 0** | Conda env `activecellsim` (Python 3.11.15, Taichi 1.7.4), `acs/` editable package, `ACS_GPU_BACKEND` dispatcher (cuda/vulkan/opengl/metal/cpu/auto), structured logging, run-provenance manifest, pytest + ruff, 4 staircase configs (`stage1a_{pilot,mid,prod_burnin,production}.yaml`), `scripts/verify_env.py`. Stage 0 closed. |
| **Stage 1a — staircase + scheme decisions** | Roadmap §1a updated with 4-step staircase (pilot → mid → prod-burnin → production) and PASS/FAIL gates. Numerical scheme locked: overdamped MLS-MPM (Re ≈ 10⁻¹³), exponential Maxwell integrator, density-based boundary tag, CSF on background grid. Dimensionless units (length=R₀, time=τ_relax, stress=K). `docs/v1/stage1a_assumption_review.md` written before any code ran. |
| **Pilot v1 (Stage 0 deferred CFL fix)** | Discovered explicit-elastic CFL violation (dt=0.02s vs limit ≈1 μs, 16,000× over) before writing solver kernels. PI chose overdamped + dimensionless path. Sanity-Gate Protocol caught the issue and was added to CLAUDE.md Hard Rules. |
| **Pilot v1 (real)** | NaN at step 200 (`is_boundary` shell colour function double-edged force, plus CSF impulse applied to grid_v *before* mass division). Initial step time 137 ms/step (atomic_max global serialisation). |
| **Pilot v2** | Boundary tagger rewritten to "empty-neighbour count" (no global atomic_max), CSF moved to `_grid_op_overdamped` after mass division, ρ_local clamp. Step time 1.56 ms/step (88× faster). NaN gone. R/R₀ drifts +12.5% (expansion) — bulk pressure dominated CSF. |
| **Pilot v3 (reference calibration, Hu 2018 §4.3)** | Per-particle F = (ρ_ref/ρ_actual)^(1/3)·I after one mass-only P2G round trip. Drift improved 12.5% → 8.5% but still expanding. **Found Jensen's inequality bias**: arithmetic-mean ρ_ref gives mean(J) > 1 → outward pressure. |
| **Pilot v3b (harmonic mean)** | ρ_ref = harmonic mean of well-resolved particle densities. Drift 8.5% (expansion, slow). Diagnostic: `<J>=1.556` after calibration — clearly **NOT 1**, contradicting the harmonic-mean derivation. *Unresolved at handoff* (see "Open root cause" below). |
| **Pilot v4 (CSF sign)** | Switched colour function from boundary-only shell to all-particle "ρ/ρ_bulk inside, 0 outside" with corrected sign (impulse `+γ∇c·dt/ρ`). R/R₀ now contracts (correct physics) but to 0.824 — way past Laplace prediction R_eq ≈ 0.99. |
| **Pilot v5 (colour normalisation)** | `c = grid_m / (ρ_bulk · dx³)`. Stronger contraction (R/R₀ = 0.657, drift 34%). Ran step ordering rebuild (P2G first, then CSF). |
| **Pilot v6, v7, v8 — MAGIC NUMBER TUNING (REVERTED)** | Introduced `csf_kappa_scale` knob to compensate for missing curvature factor: 0.25, 0.05, then dropped Ca to 0.001. PI flagged this as exactly the anti-pattern Sanity-Gate Protocol exists to prevent. **All three of those changes reverted.** Rule added: **Magic-Number Block** in CLAUDE.md Hard Rules + `docs/12_validation.md` (3-test gate, plus cousin rule that gate semantics are immutable per run). |
| **Pilot v8 (clean — Brackbill curvature operator)** | Implemented full Brackbill 1992 §IV: `n̂[I] = ∇c/√(\|∇c\|²+ε²)`, `κ[I] = -∇·n̂`, impulse `dv = +γ·κ·∇c·dt/ρ_local`. Added **static-sphere curvature gate** (`κ_measured` vs analytical `2/R`, ±10%). FAIL: `κ_measured = 8.24 ± 11.1` vs `2.0` (312% over). κ encoding grid-scale (≈ 1/dx ≈ 10) instead of physical 2/R. |
| **Pilot v9 (Step 1, option C — selection criterion)** | Switched κ-band selection from "top-decile \|∇c\|" to "boundary-tagged-particle cells with \|∇c\| ≥ 25% of band max". `κ_measured = 9.72 ± 11.1` (worse). Selection isn't the issue; the κ field itself encodes grid-scale. |
| **Pilot v10 (Step 2 — Brackbill 1992 §V colour smoothing)** | 2-pass 3³ box averaging of colour before ∇c, n̂, κ. Improved: `κ_measured = 5.85 ± 6.17` (192% over vs 312% in v8) — smoothing helps but is not enough. **STOP HERE for PI analysis** (per directive). |

## Discovered root causes and fixes (cumulative)

| # | Root cause | Status | Anchor |
|---|---|---|---|
| 1 | Explicit-elastic CFL violation (dt=0.02s, K=1kPa → dt_max ≈ 1μs) | Fixed: overdamped + dimensionless reformulation | Marchetti 2013, Pérez-González 2019 |
| 2 | `is_boundary` shell colour produced double-edged ∇c force | Fixed: `c = ρ/ρ_bulk` from all particles | Brackbill 1992 |
| 3 | CSF impulse applied to grid *momentum* before mass division | Fixed: moved into `_grid_op_overdamped` after `v = grid_v / m`, with min-mass clamp | standard MPM convention |
| 4 | `atomic_max` of grid_count to single scalar (boundary tag) → 137ms/step | Fixed: empty-neighbour count tag, no global reduction | — |
| 5 | Reference state F=I gives σ_vol(t=0) ≠ 0 due to ρ_actual ≠ ρ_nominal | Fixed: harmonic-mean ρ_ref, F_p = (ρ_ref/ρ_p)^(1/3)·I | Hu et al. 2018 §4.3, Jiang et al. 2015 |
| 6 | Curvature-less CSF (`F_v = γ·∇c`) over-drives surface tension | Fixed: full Brackbill `F_v = γ·κ·∇c` with `κ = -∇·n̂` | Brackbill 1992 §IV |
| 7 | Single-cell-thin diffuse interface → κ encodes 1/dx not 2/R | Partially fixed: 2-pass colour smoothing reduced κ from 8.24 → 5.85, still 192% over | Brackbill 1992 §V |

## v10 status (current, blocked)

```
Curvature gate     FAIL  κ_measured = 5.85 ± 6.17, analytical 2/R = 2.0   (192% over, limit 10%)
Mass               PASS  |Δm/m₀| = 0  (exact)
Energy monotone    PASS  ΔE_max = 6.4e-5 ≤ tol 1.2e-3
No NaN             PASS  nan_count = 0
Sphericity ψ       PASS  min ψ = 0.983
Max speed          PASS  v_max = 1.2e-3
VRAM peak          PASS  2.23 GB
Momentum drift     FAIL  0.065 (downstream of κ — over-strong CSF)
Radius drift       FAIL  44.2% (downstream of κ)

Wall-clock 0.83 min, GPU util 22%, peak VRAM 2.23 GB.
```

## Open root cause — reference calibration drift

This is the new finding the PI flagged before going to bed. After `calibrate_reference_state()` runs, the diagnostic reports:

```
ρ_ref(harmonic) = 2.5367   ← expected 1.0 (the configured density_star)
<J> after calib = 1.5558   ← expected exactly 1.0 by construction
F_scale ∈ [0.843, 1.567]
```

Two things are wrong here, and only one (Jensen) was supposedly fixed in v3b:

1. **`ρ_ref ≈ 2.54 ≠ 1.0`** even though `density_star = 1.0`. Reason this is "okay" *within the calibration's own logic*: ρ_ref is the harmonic mean of `grid_m / dx³` interpolated back onto particles. It is NOT the configured `density_star`. So the units of "ρ_actual" used inside calibrate_reference_state() are not the same units as the configured reference. They differ by a particle-volume-vs-grid-volume scaling (each particle scatters mass over 27 grid cells, so the kernel-weighted local density is larger than the per-particle bulk).
2. **`<J> = 1.556` after harmonic-mean calibration**. Mathematically this should be exactly 1.0:
   ```
   F_scale_p = (ρ_ref / ρ_p)^(1/3)
   J_p = det(F_scale_p · I) = F_scale_p^3 = ρ_ref / ρ_p
   <J> = ρ_ref · <1/ρ_p>
   harmonic mean ρ_ref = 1 / <1/ρ_p>
   ⇒ <J> = (1 / <1/ρ_p>) · <1/ρ_p> = 1   exactly.
   ```
   The fact that we report `<J> = 1.556` means **either** the harmonic mean is being computed only over a sub-population (`well_resolved`) but `<J>` is being averaged over **all** particles (Jensen bias re-introduced), **or** the F_scale clamp (`[1/8, 8]^(1/3)`) is hiding the violation in tail particles.

This is consistent with the symptoms: even with proper Brackbill curvature and colour smoothing, R/R₀ drifts heavily because `<σ_vol>(t=0) = K · (<J> − 1) = 0.556 ≠ 0`. **The system starts pre-stressed outward by ≈ 50% K**, which dominates everything.

## Two diagnoses for next session

**(a, most likely)** The reference-calibration sub-population mismatch above. Fix candidate: compute the harmonic mean over the *same* population the F-scale will be applied to (i.e., do not use `well_resolved`-only when then applying to all particles). Or: re-derive ρ_ref by an iterative pass that drives `<J>` to 1.0 exactly. Either path is principled (Hu/Jiang approach, no magic number).

**(b)** Smoothing of 2 passes is not enough — Brackbill 1992 §V mentions 1–3 passes; the upper end may be required at this grid_n. If (a) is fixed and curvature still fails, try 3 or 4 passes and re-measure.

**(c)** A different fundamental issue still hidden in the colour/curvature pipeline (e.g., the FD divergence operator on `n̂` is biased near the interface where `\|∇c\| → 0`, even with the ε² regulariser).

## Where to start the next session

1. Read `docs/SESSION_HANDOFF.md` (this file) + `docs/v1/stage1a_assumption_review.md` + the last commit message.
2. Open `acs/physics/mlsmpm.py` and **review `calibrate_reference_state()` directly** (the section between `def calibrate_reference_state` and the end of `_set_F_isotropic_from_calib`). Specifically check:
   - Does `rho_used = well_resolved if len(well_resolved) else rho_np` introduce a population bias when the F-scale is then applied to *all* `rho_np`?
   - Does the `np.clip(rho_ref / rho_np, 1/8, 8)` clamp suppress the violation evidence?
3. Then re-derive: what does `<J>` actually equal if we compute the harmonic mean over the well-resolved subset and apply the F-scale to all particles? Is it 1.0 or biased?
4. Decide diagnosis (a) vs (b) vs (c).
5. Apply the principled fix (NOT a magic number, NOT a gate-tolerance change).
6. Re-run pilot v11.

## Bounded-protocol reminder for the next session

The PI's auto-mode rules from this session remain in force unless explicitly relaxed:

- **Stop conditions** still apply: magic number temptation → stop; gate-semantic-change temptation → stop; same root cause debugged twice → stop; ≥2 fundamental issues → stop; 60% token usage → stop.
- **Magic-Number Block**: any new empirical scaling factor blocked unless it passes the 3-test in `CLAUDE.md` Hard Rules.
- **Gate-Semantics-Immutability**: gate tolerances and check windows are part of the validation contract written before the run; do not edit them inline to make a failing run pass.
- **Two-outcomes-only for the next pilot**: either all PASS (R/R₀ ≈ 0.99 ± 0.01, κ within ±10%) → stop and surface; or FAIL → stop and surface. Do NOT enter mid-run automatically.

## Task tracker state at handoff

The active task is **#27 "Step 2 — Brackbill 1992 §V colour smoothing"** in the in-progress state. v10 result above. Step 3 (Stage 1a-mid) was queued but **NOT** entered per the Stop conditions. All earlier debug-cycle tasks (#26 "Step 1 (option C, selection criterion) — FAILED") are completed for the record.

---

## Stage 1b framing memo — Bare/Pre/Lam4 = phenotype, NOT substrate (recorded 2026-04-29)

When Stage 1b activates Layer 3 (φ-ODE for E-cad ↔ Int-β1, see
`docs/03_adhesion_dynamics.md`), the PI's three experimental conditions
enter the simulation. Use this exact framing — re-deriving it from
scratch in a new session is wasted effort and risks the substrate-
misframing trap that was caught and corrected on 2026-04-29.

**One-liner**: *Bare/Pre/Lam4 = same Col1 substrate, different spheroid
initial phenotype (formation environment dependent).*

**Mechanism**: PI seeds MCF7 single cells onto pV4D4 hydrogel under three
laminin-presentation conditions during spheroid formation:

| Condition | Formation environment |
|---|---|
| `Bare` | pV4D4 only, no protein presentation |
| `Pre` | pV4D4 + pre-adsorbed laminin |
| `Lam4` | pV4D4 + laminin in media (4 µg/mL) |

This three-condition design extends Cho et al. 2020's ULA-derived → pV4D4
framework with the PI's own laminin-presentation tuning. Different
laminin-integrin engagement during formation produces different
**E-cadherin / Integrin-β1 ratios** in the formed spheroid (Cho 2020
mechanism, well-established). After formation, all spheroids are *transferred*
onto the same Col1-coated confocal dish for the 24–96 hr spreading assay.

**Stage 1b mapping** (when implemented; not now):
- Activate Layer 3 φ-ODE.
- Three runs with three different *initial φ values*, mapped from
  formation environment:
  - `Bare` → low laminin engagement → E-cad-dominant initial state →
    low initial φ (E-cad/(E-cad+Int) ratio)
  - `Pre` → moderate laminin engagement → mid initial φ
  - `Lam4` → high laminin engagement → Int-β1-dominant initial state →
    high initial φ
- Substrate stays Col1-only across all three runs. Single γ_sub_Col1
  value (the same anchored in `docs/v1/stage1a_plus_substrate_sanity.md`).

**Anti-pattern to avoid**: do *not* parameterise Bare/Pre/Lam4 as
γ_sub_{Bare,Pre,Lam4}. The spreading-time substrate is identical in all
three cases. Inter-condition variance in `A/A₀(t)` is attributable to
the spheroid's formation-derived starting state, not to spreading-
substrate differences. This was caught and fixed during the Stage 1a+
sanity-md review on 2026-04-29; see also memory entry
`memory/experimental_design.md`.


---

