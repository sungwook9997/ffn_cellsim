# Phase C — implicit accelerator GPU-verify + decisive (i)/(ii) spread test (2026-06-22)

Continues `PHASE_C_SESSION_HANDOFF_2026-06-22.md`. Resolves its two open questions.
Branch `h7/compartment-platform`. Engine: Warp DCM, A5000 (gbook).

## Headline

**The fully-mechanistic (proxy-free) spheroid DOES spread — answer (ii) slow-transient,
not (i) no-spread-equilibrium. A/A0 climbs 1.0 → 1.94 over T≈3–18 s of physical time and
plateaus. The earlier "A/A0 = 1.001 FLAT" headline was a TIMESCALE ARTIFACT** — the explicit
production only reached T_spread ≈ 0.96 s, far below the minutes-scale onset of lamellipodial
spreading. Reaching the decisive horizon required the implicit accelerator.

## 1. Implicit driver was GPU-broken — fixed (commit 917e51d)

The wired `--integrator implicit` path crashed on the **first spread step** on BOTH GPU and CPU
(`ZeroDivisionError` on `1.0/s` in `device_cg`). The handoff's "CPU-validated" only covered the
standalone `_spring_demo`/I2–I3, never the full multicell+contact+nucleus driver. Three root-cause
fixes in `dcm_warp_implicit.py::device_cg` (all standalone I1/I2/I3 unchanged after):

1. **Crash** — near equilibrium F(x)≈0, the forward-difference JVP catastrophically cancels
   (Fp−Fx≈0); ‖p‖→inf underflows the perturbation `s=eps/vn` to exactly 0 → `1/s` crashes.
   Fix: early-out when F(x) already at tolerance (Δx≈0 correct) + degenerate-direction guard.
2. **8e13 blowup at high dt** — A = (γ/dt)·I + K. Diagonal `a·‖p‖²` is exact & positive; physical
   stiffness K is PSD at a stable equilibrium, so `pᵀAp ≥ a·‖p‖²` ALWAYS. Near equilibrium the FD
   K·p is rounding noise and `pᵀKp` goes spuriously negative → runaway α=rs/pAp. Fix: floor pAp at
   the exact diagonal (no scale-tuned constant).
3. **Intra-solve runaway** — non-symmetric, K-dominated operator at large dt (a=γ/dt ≪ ‖K‖) diverges
   to ~1e128 in one solve. Fix: divergence guard (break at 4× initial residual, return current iterate).

## 2. Real speedup — measured, honest (GPU A5000, full fine-grained N=24)

| integrator | dt | steps/s | per-step cost | **net wall-clock speedup** |
|---|---|---|---|---|
| explicit BAOAB | 8e-6 | 250.1 | 1× | 1× (ref) |
| implicit 30× | 2.4e-4 | 17.70 | 14× | **2.1×** |
| implicit 50× | 4e-4 | 10.11 | 25× | **2.0×** |
| implicit 100× | 8e-4 | 7.14 | 35× | **2.85×** |

Matrix-free CG costs ~14–35 stiff-force evals/step, so the net win is **~2–3×**, NOT the dt-ratio.
100× is the best net and is stable + physics-matched (V/V0≈1.000, aa0 agrees across dts). N=100
implicit-100× rate = 2.94 steps/s.

**Cadence validity:** lamellipodium `p_advance = v_front·(batch_steps·dt)·S/ℓ₀` and cadherin/ecm
`p_on,p_break = 1−exp(−k·(batch_steps·dt))` are dt-aware → event rates are physical-time-invariant.
Implicit-at-100× is therefore physically equivalent to explicit at matched physical time (only the
binding-batch granularity is coarser, not biased). The decisive test below is valid.

## 3. Decisive run — N=100 implicit 100×, T_spread ≈ 30 s

Config: `--n-cells 100 --subdiv 2 --warmup 200 --settle-steps 1000 --gap 2.05 --builder fcc
--cadherin --ecm-clutch --lamellipodium --nucleus --surface-tension --bending --edge-edge
--steps 37500 --frames 60 --integrator implicit --accel-dt 8e-4`. 37500 steps @ 2.94 steps/s ≈ 3.5 h.

A/A0 trajectory (top-down silhouette, the PI-mandated metric — `_topdown_area_um2`):

| T_spread (s) | A/A0 | maxZ (µm) | V/V0 | pen | drift (µm) | ecm-clutch engaged |
|---|---|---|---|---|---|---|
| 0.5  | 1.000 | 80.3 | 1.000 | 2.08 | 0.00 | 55 |
| 3.5  | 1.231 | 80.4 | 1.000 | 3.96 | 0.03 | 218 |
| 6.5  | 1.577 | 80.4 | 1.000 | 0.00 | 0.14 | 303 |
| 9.5  | 1.836 | 80.5 | 1.000 | 0.00 | 0.15 | 373 |
| 15.5 | 1.927 | 80.5 | 1.000 | 0.00 | 0.14 | 380 |
| 21.5 | 1.938 | 80.5 | 1.000 | 0.00 | 0.16 | 368 |
| 30.0 | 1.942 | 80.6 | 1.000 | 0.00 | 0.16 | 390 |

**Reading:** flat at T=0.5 s (= the explicit production's whole window), then spreads over T≈3–18 s,
plateau A/A0 ≈ 1.94. Footprint expands while maxZ holds 80→81 µm and V/V0=1.000 — genuine spreading,
not collapse. COM drift 0.16 µm (symmetric, not translation). ecm-clutch engagements 55→390 (the
Pereverzev traction driver engaging at the spreading rim); cadherin bonds hold ~2280 (1.30M formed /
1.30M broken — turnover at equilibrium). rim 7/100 cells, 93 dragged.

## 4. Resolves the handoff open questions

- **Q1 (why A/A0 flat?)** → answer **(ii)**. Not a no-spread equilibrium; a slow transient the
  explicit run was 30× too short to see. The mechanistic stack spreads to A/A0 ≈ 1.94 with the
  wetting/settling PROXY fully removed (ecm-clutch + lamellipodium do the work). No spreading driver,
  magic number, or proxy was added (PI prohibition respected) — only solver robustness + a longer run.
- **Q2 (pen=2.49 G2 FAIL)** → **transient, self-heals**. pen peaks 3.96 during the fast-spread phase
  (T≈3.5 s) then returns to **0.000** for the rest of the run (final 0.000). The interpenetration is
  a fast-spreading overshoot the contact penalty recovers from, NOT a stuck internal overlap. G2 read
  on the converged state passes.

## Figures
- `outputs/warp_decohesion/figs/n100_implicit_decisive_contact_montage.png` — 10-frame top-down +
  side montage, A/A0 0.99→1.94, footprint expands while height holds (rim cells red).
- `outputs/warp_decohesion/figs/n100_implicit_decisive_contact_surface.mp4` — 62-frame surface movie.

## Open for PI
- **Magnitude calibration** vs PI experimental A/A0=a+b/R+c/R² (overlay-only, no fitting): A/A0≈1.94
  for N=100. An R-sweep would give the a,b,c form. Is the platform now ready for the form comparison?
- Optional: fold I4 Newton into `device_cg` to push the stable dt past 100× (would raise the net
  speedup ceiling above ~3×).
- This **retracts** the memory line "spheroid COMPACTS not spreads" for the proxy-free Warp stack.

## Artifacts
`outputs/warp_decohesion/n100_implicit_decisive.{json,log,npz}` (npz on gbook). Run script:
`run_decisive_implicit.sh` (gbook). Benchmark: `bench_implicit.sh`, `bench_hi.sh` (gbook).
