# Oracle-B blebbing — growth driver (readiness)

Status: **growth driver complete + hardened; native execution GATED on (b) full-native resting convergence.**

The single-shot nucleation criterion already lands on the A5000 (`bleb_70686.json`): with a 25° cap on the
70,686-filament cell the measured local ΔP = 40 Pa and γ_mem = 10 pN/µm give `a_crit = 2γ/ΔP = 0.5 µm`, and
the 3.17 µm de-adhered cap `nucleates = true`. This report covers the **multi-step growth** driver
(`ac/cell/bleb_growth.py`) that carries that nucleation into the mechanical outer loop and measures the
emergent bleb. It is finished and hardened so it runs the moment resting convergence lands; it does **not** run
now, by construction (see *Dependency*).

## What the driver does (all mechanistic, no lumped bleb)

`run_bleb_growth` reuses the frozen resting machinery (`PhysicalScheduler` + `make_inner_solve`, imported
read-only — no shared-driver edit) and adds only a perturbation and diagnostics:

1. **Equilibrate to the resting baseline.** `equilibrate_steps` outer steps relax the composed cell (cortex +
   membrane + ERM + Biot pressure at Π₀ = 40 Pa) at its physiological setpoint. This is the physiological
   baseline the bleb nucleates *from* — not a freshly-built, un-relaxed shell.
2. **Commit-safe local ERM de-adhesion.** The angular-cap ablation (`bleb_perturbation`) is irreversible
   biology, so — exactly like the ERM Bell KMC and the nuclear-envelope rupture — it fires only at an
   **accepted** outer physical-time boundary. The driver launches it from inside the scheduler's
   `commit_irreversible` hook, predicated on the authoritative device acceptance scalar, and exactly once
   (idempotent `1→0`). A rejected mechanical candidate never leaves a bleb behind.
3. **Emergent readouts, per step:** the baseline-subtracted outward displacement of the de-adhered patch vs the
   intact control membrane (the bulge); the mean cap pore ΔP vs the far field (Charras non-equilibration
   signature); and the Laplace balance `a_patch` vs `a_crit(t) = 2γ/ΔP(t)` (the cap bulges while
   `a_patch > a_crit`). The bulge, its expansion, and the pressure transient are outputs of the membrane
   bending + tension + Biot physics — never scripted.

## Dependency — why it does not run now (the (b) gate)

The driver makes the resting-convergence dependency **explicit and observable**. If any equilibration step is
rejected — i.e. full-native resting mechanics did not converge — it HALTS and returns
`resting_converged = false`, `blocked_on = "resting-convergence"`, and applies **no** de-adhesion. Today the
70,686-filament resting state does not converge (implicit outer step rejected; residual ~O(10³) pN from
long-wavelength fiber modes + WCA ill-conditioning — the WCA session's critical path), so this run would stop
at step 1 rather than measure a bleb off an unphysical baseline. This is intentional: per the physiological
baseline hard rule, comparing a bleb grown from a non-converged resting state to reality is meaningless.

The growth driver takes the **same solver-selection flags** as `run_from_resting`
(`--inner-solver`, `--max-inner-retries`, `--implicit-*`, `--rkc-stages`), so once the WCA session's
conditioning fix makes resting converge under some config, the growth run uses that identical config.

## Run procedure once (b) lands (gbook A5000, Warp-CUDA only)

Preconditions: `run_from_resting` reports `outer_accepted = true` / `inner_converged = true` at the full native
70,686-filament population under some `--inner-solver …` config. Then, on the gbook (log to a file; monitor the
log, not the terminal):

```bash
# on gbook (RTX A5000); <CFG> = the flags that made resting converge
~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.ac.cell.bleb_growth \
    --n-filaments 70686 --patch-deg 30 \
    --equilibrate-steps 8 --growth-steps 40 \
    <CFG e.g. --inner-solver contact_tournament --max-inner-retries 4> \
    --out ffn_sim/outputs/ac/bleb/growth.npz \
    > ffn_sim/outputs/ac/bleb/growth.log 2>&1 &
# monitor:
tail -f ffn_sim/outputs/ac/bleb/growth.log
```

Then render the real figure on any machine (CPU, no CUDA):

```bash
python -m ffn_sim.scripts.ac_bleb_growth_vis --npz ffn_sim/outputs/ac/bleb/growth.npz
```

### What validates the bleb

Read `growth.json`. The run **passes** the growth signature when, with `resting_converged = true`:

- `bleb_deadhered = true` and `n_tethers_deadhered` ≈ the cap-membership count (`n_tethers_patch`);
- `all_growth_steps_accepted = true` (the post-de-adhesion cell stays a converged transaction);
- `bleb_grows = true` — the patch ends outward beyond the control AND beyond its own de-adhesion-step value;
- `cap_expands_last = true` — `a_patch` (µm-scale) stays above `a_crit(t) = 2γ/ΔP` (sub-µm), the Laplace
  criterion for continued expansion;
- `dp_cap_last_Pa` shows the local pore ΔP evolving (the non-equilibration transient) rather than instantly
  re-equilibrating to the far field.

If `bleb_grows = false` while resting converged and de-adhesion committed, that is a **finding** to surface to
PI (e.g. quasi-static per-step growth too small at the sourced k_erm CFL — extend `--growth-steps` / `--dt-phys`
and re-examine, do not re-tune sourced parameters to force an outcome).

## Analysis-only vs native-execution

- **Validated now (CPU, dev Mac, no CUDA):** the host analysis pipeline — cap membership, patch/control outward
  displacement, the ΔP-set Laplace critical radius, and the bulge/grows verdict — via
  `tests/ac/cell/test_bleb_growth_reference.py` (8 tests, pure NumPy on analytic geometry). The Laplace anchor
  is checked against the committed nucleation numbers (γ=10, ΔP=40 → `a_crit = 0.5 µm`). These functions are
  the *same* ones the device path calls on its post-step host reads.
- **NOT claimed now:** any physical bleb outcome (bulge magnitude, growth rate, ΔP transient). That requires the
  native Warp-CUDA run above, which is gated on (b).

## Figures

- `figs/bleb_growth_analysis_validation.png` — ANALYSIS validation (synthetic area-uniform membrane sphere,
  CPU): the +z 30° cap is 7.0% of the sphere (matches the (1−cos θ½)/2 solid-angle fraction), and the analysis
  pipeline recovers an imposed cap bulge exactly (patch on y=x) while the control stays at ≈0. Not a physical
  bleb claim.
- `figs/bleb_nucleation.png` — the single-shot nucleation criterion from `bleb_experiment` (pre-existing).
- `figs/bleb_growth.png` — the REAL growth figure (patch/control outward displacement + cap ΔP with `a_crit(t)`
  overlay); generated by `ac_bleb_growth_vis.py --npz growth.npz` after the gbook run lands.
