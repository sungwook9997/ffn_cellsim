# H.7 — GPU-device myosin grip_walk s_grip bug (BLOCKS Gate-A); precise repro + elimination

**Date:** 2026-06-07 (night, autonomous) · **Status:** OPEN, blocks Gate-A · **Severity:** high
**Found by:** Lead, while launching the Gate-A production run on the GPU-main stack.

## TL;DR
On a **GPU** device the cortex-myosin `grip_walk` grip-stretch accumulator
`_head_grip_s` stays **exactly 0.0** indefinitely, even though heads bind (78%),
the Hill walk velocity `v_step ≈ v0 = 2e-7 m/s` is computed nonzero every tick,
and `batch_dt = 6.98e-5 s`. The **identical code + config on a CPU device grows
s_grip** (~14 pm/tick, reaching 2.65e-10 m after 1500 steps). ⇒ a GPU-device-
specific bug in the myosin walk write-back. **Gate-A ("does γ rise as s_grip →
ℓ₀?") is INVALID on GPU** until fixed — the experiment cannot develop s_grip.
The 2×10⁸-step (~5-day) Gate-A run was launched, this was caught at tick 100,
and the run was **stopped** (≈5 GPU-days saved). Do NOT relaunch Gate-A on GPU
until s_grip accumulates.

## Definitive apples-to-apples (identical code, n_fil=100, constrained, cupy integrator, NO native compartments)
| burst (500 steps each) | CPU `sgmax` | GPU `sgmax` |
|---|---|---|
| 0 | 1.26e-10 | **0.000e+00** |
| 1 | 1.95e-10 | **0.000e+00** |
| 2 | 2.65e-10 | **0.000e+00** |

Both bind heads (CPU 45/62/87, GPU 50/75/101), `n_bind_total` monotonic == n_bound
(no churn/rebind-reset), `_n_step_advances_total = 0` (s_grip never reaches ℓ₀).

## Elimination chain (what it is NOT)
1. **NOT the native constrained integrator.** GPU + *cupy* `ConstrainedLeimkuhlerMatthewsBAOAB` (no swap) also pins s_grip = 0.
2. **NOT the native compartment ForceCompute.** GPU with `FFN_GPU_DEVICE_COMPARTMENTS` unset still pins s_grip = 0.
3. **NOT bond-dropping.** On GPU `attach_snap == n_bound` (50/75/101) — the `cortex_myosin_attach*` bonds persist in the snapshot, so Step-3 (`if attach_bonds.shape[0] > 0`) runs.
4. **NOT bond column-order swap.** On BOTH CPU and GPU the attach-bond column order is identical: `col0_is_head = 1.00`, `col1_is_head = 0.00`. So `head_locals` (derived from `attach_bonds[:,0]`) is correct on both.
5. **NOT v_step = 0.** Monkeypatching `myosin.hill_velocity_clamped` on GPU shows it is called (50/75/101 rows), `F_load_max ≈ 2.4e-14 N` (sub-stall), returning `v_step_max = 2e-7 = v0`, `v_step_mean ≈ 1.97e-7`.
6. **NOT batch_dt = 0.** `cell.myosin_action.p.batch_dt = 6.98e-5 s` on GPU.
7. **NOT array reassignment.** Wrapping `myosin_action.act()` and reading `_head_grip_s` immediately after each call (same `id()`): `max = 0.0`, `nonzero = 0` after acts #1,2,5,20,60.
8. **NOT a GPU branch in the updater.** `cortex/myosin.py` has no device/`gpu_local`/cupy branch — it reads `sim.state.get_snapshot()` uniformly (line ~1029) and updates `_head_grip_s` in plain host numpy (`myosin.py:1386` `s_new = s + v_step*batch_dt`; writes at `:1396`/`:1427`).

## The remaining pinpoint (for the fix)
`grip_walk` Step-3 runs with the correct head indices and a nonzero `v_step`, yet
`_head_grip_s[head_locals] = s_new` does not persist on GPU. So `s_new` is
evaluating to 0 *for the bound heads*, or the write targets the wrong rows, in a
way that only manifests on the GPU device. The next step is **inside-`act()`
instrumentation on GPU**: print `head_locals`, `s` (= `_head_grip_s[head_locals]`),
`v_step`, `s_new`, and `_head_grip_s` immediately before/after the write loop
(`myosin.py:1386-1427`). Leading suspects to check there:
- `head_locals` dtype/values when `head_tags`/`per_motor` arithmetic meets a
  GPU-sorted snapshot (even though column order is preserved, verify `head_locals`
  are the same int indices CPU produces);
- whether `s_new[row]` aligns row-for-row with `head_locals[row]` on GPU
  (any reorder between the `hill` input rows and the write loop).
A robust hardening (mirrors `bridge/manifold_traction.py`'s by-type resolution):
resolve the head index from particle TYPE rather than column position + assert
`head_locals` are in `[0, n_heads)` and unique, so a snapshot reorder can never
silently mis-target the write.

## Reproduce (gbook)
```
source ~/miniconda3/etc/profile.d/conda.sh; conda activate ffn_sim; cd ~/ffn_cellsim
# device cpu → s_grip grows; device gpu → s_grip stays 0.0
python ffn_sim/scripts/h7_sgrip_gpu_diag.py cpu
python ffn_sim/scripts/h7_sgrip_gpu_diag.py gpu
```

## Impact / disposition
- **Gate-A (GPU) PAUSED** — not relaunched (would waste ~5 GPU-days producing a
  meaningless s_grip = 0 verdict).
- The GPU-main throughput work (native integrator 5.86×, native compartment
  ForceCompute 8.3×, γ bit-parity) is UNAFFECTED and correct — this is an
  independent myosin-walk write bug.
- A valid Gate-A could run on **CPU** today (s_grip develops there) but at the
  cupy CPU throughput it is far slower; the right fix is the GPU s_grip write-back
  so Gate-A runs at the ~4.7-day GPU-main rate.
