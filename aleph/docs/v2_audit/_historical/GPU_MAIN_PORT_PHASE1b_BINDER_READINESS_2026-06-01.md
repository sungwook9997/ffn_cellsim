---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# GPU-main port — Phase 1b: native readiness + binder-cost data (2026-06-01)

> Autonomous follow-on to the cupy port #1 (`GPU_MAIN_PORT_PHASE1_2026-06-01.md`).
> Goal: **prepare the myosin-ON native run**. Findings: (1) the large-N build
> blow-up is fixed by soft-start equilibration; (2) the binding updaters are
> **NOT the bottleneck even at native scale** (data below) → the risky binder
> cupy port is **deferred** (data-grounded); (3) full-cortex GPU pulls ahead of
> CPU earlier than the solver-only benchmark, as predicted.

## 1. Large-N build blow-up — fixed (soft-start equilibration)

The int32-image blow-up at `n_fil≥1000` (`GPU_MAIN_PORT_PHASE1` §3 caveat) was
**not a packing bug** (Explore-confirmed: `cortex.py:590–686` places filaments
uniformly on the sphere — sparser, not denser, at large N). Cause: the profile
script's warm-up did a **raw `sim.run(n_warmup)`** with no soft-start, so an
LJ overlap drove a per-step displacement past the int32 guard.

Fix (script-only, low-risk): `scripts/h3_profile_constrained.py::_warm_and_constrain`
now uses the **build's existing equilibration** — `build_cortex_full_simulation(
equilibrate=True, equilibrate_steps=n_warmup, equilibrate_softstart_steps=…)` —
which runs `equilibrate_cell` → `ecm.equilibrate._softstart_step` (clipped-Brownian:
per-step `|dr|` capped at `frac·ℓ₀`, so overlaps drain safely regardless of N)
then a BAOAB drain. **The production build path already supports `equilibrate=True`**
(`cell.py:584`, wired at `cell.py:1323`) — so a native production driver just
passes `equilibrate=True`; no production code change was needed.

Result: `n_fil=1000` now runs clean on **both** devices, ON and OFF (§2).

## 2. Native-scale binder cost — the decisive data (gbook A5000)

`h3_profile_constrained.py --no-cprofile`, 800 timed steps, batch-binding myosin
ON vs OFF (OFF = pure constrained BAOAB; the ON−OFF gap = the binding updaters'
per-batch `get/set_snapshot` + cKDTree + grip-walk cost, amortised):

| n_fil | N_part | GPU ON | GPU OFF | GPU binder Δ | CPU ON | CPU OFF |
|---|---|---|---|---|---|---|
| 1,000 | ~7,000  | 97.3 | 103.9 | **−6%** | 84.0 | 102.4 (−18%) |
| 3,000 | ~21,000 | 90.4 |  97.3 | **−7%** | 32.7 | (M-SHAKE non-converge, §5) |

**Finding: the binding updaters are NOT the bottleneck at native scale** — they
cost ~6–7% on GPU (~18% on CPU) at n_fil=1000, and the GPU fraction does not
grow from 1k→3k. Per project culture (data-ground before porting), the
**binder cupy port is therefore DEFERRED** — porting the binding kinetics
(serial grip-walk while-loop, segment-projection acceptance, bond-topology
mutation) is large and science-critical (feeds KU-3.5), and the data says it
would buy <10% on GPU. Revisit only if a much-larger-N run (N~10⁵⁺) shows the
per-batch full-system `get/set_snapshot` growing into the bottleneck.

## 3. Full-cortex GPU crossover (with LJ + neighbour list)

Unlike the solver-isolated benchmark (Phase-1 §3, no pair force), this is the
**real cortex** (LJ + nlist + bonds + angles + binders). GPU/CPU speed-up (myosin
ON): **1.16× at n_fil=1,000 → 2.8× at n_fil=3,000.** The full-cortex GPU pulls
ahead of CPU at a *smaller* N than the solver-only crossover (~n_fil=1,000)
because the 71% C++ LJ/nlist is GPU-accelerated for free — **confirming the
Phase-1 caveat** that real large-N full-cortex GPU speed-up beats the
solver-only number. (Solver-only crossover was N~10k particles; full-cortex
crosses at ~7k and is already 2.8× by 21k.)

## 4. Binder cupy-port design (documented, NOT executed — for when warranted)

From the Explore map (`myosin.py:735–1318`, `crosslinkers.py:532–814`). Each
binder `act()` fires every `batch_steps` (~100) and does a full host
`get_snapshot()` (pos + bonds) → numpy/scipy logic → full `set_snapshot()`
(rebuilds bond topology). Classification for a future port:

| piece | class | plan |
|---|---|---|
| unbind loops (`myosin.py:1003`, `crosslinkers.py:696`) | EASY | vectorise (Python loop → boolean-mask scatter); CPU-bit-equivalent |
| Bell-Evans / Pereverzev rates | EASY | already vectorised numpy → cupy elementwise |
| cKDTree `query_ball_point` candidate search (`myosin.py:1043`, `crosslinkers.py:720`) | REPLACE | HOOMD GPU nlist as head→actin candidate source, OR cupy fixed-radius cell-list |
| segment-projection + bipolar acceptance (`myosin.py:1067`) | HARD | per-head conditional logic; CPU-gather a small candidate set + host accept |
| **grip-walk while-loop** (`myosin.py:1226`) | SERIAL | inherently sequential per head (depth 1–3 typical); **stays host** on a small gathered array |
| full bond-topology mutation (`set_snapshot`) | HARD | needs host bond table + nlist rebuild; fully device-resident requires a compiled plugin |

**Recommended (if/when revisited)**: Phase-1 = drop the full `get/set_snapshot`
for a `gpu_local_snapshot` position read + bond-array-only writes; keep cKDTree +
grip-walk host on a small gather. Additive + CPU-bit-invariant, same as the
constrained-Action port. Est. modest (<10% at current N per §2) — low ROI now.

## 5. Open items / native-run readiness

- ✅ **Native myosin-ON run is ENABLED**: `build_cortex_full_simulation(
  device=GPU(), constrained=True, equilibrate=True, with_baoab=True, p_myosin=…)`
  uses the GPU-resident constrained Action (port #1) + soft-start equilibration +
  host binders (cheap, §2). No code blocker remains for a moderate-native
  (n_fil~1,000–3,000) myosin-ON GPU run.
- ⚠️ **Large-N equilibration budget** (resolved → IC budget, not a bug): at
  n_fil=3,000 one realisation (CPU, OFF) hit an **M-SHAKE non-convergence**
  (drift 1.9e3) on the post-warmup constrained run with warmup=3,000 (n_soft=375)
  — an LJ overlap the soft-start budget didn't fully drain. **Re-test with
  warmup=10,000 → CLEAN** (CPU n_fil=3,000: ON 32.9 / OFF 35.5 steps/s, no
  failure). So it is an **equilibration-budget** issue, not a port bug:
  production native drivers should scale `equilibrate_softstart_steps` / warmup
  with N (rule of thumb: warmup ≳ 3·n_fil worked here). A constrained-run
  per-step overlap guard (cap predictor displacement, mirroring the soft-start
  clip) would make large-N starts robust without hand-tuning — a small future
  hardening, flagged for PI.
- 🔴 **PI-gated** (unchanged): full-length GPU L_p/equipartition (H.2/H.3 `tau_min`
  contract) before native production; any default-device flip in a production
  driver; ffn/foundation push.

## 6. Files touched

- `scripts/h3_profile_constrained.py` — warm-up now uses build soft-start
  equilibration (fixes large-N blow-up). Script-only; no production-path change.
- (no binder code changed — port deferred per §2.)
