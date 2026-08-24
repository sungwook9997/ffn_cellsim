# Resting-baseline (2026-07-22e) — fork-2 ERM-Schwarz preconditioner: PARTIAL improvement, gate NOT closed

Branch `codex/ff-ac-codex`. Executes PI-approved fork-2 (the coupled membrane↔ERM↔cortex Schwarz-block
preconditioner recommended by `_historical/RESTING_BASELINE_DIAGNOSIS_2026-07-22d.md`). New module
`aleph/components/incumbent/erm_schwarz.py` (`ERMPairSchwarz`), wired into the driver as inner-solvers `erm_schwarz` (standalone
candidate) and `erm_tournament` (block_descent + anderson + **erm_schwarz** + rkc). Probes:
`scripts/ac_erm_schwarz_probe.py`, `ac_erm_solver_sweep.py`. Full native (70,686) on the gbook A5000.

## What was built (faithful to the fork-2 spec)
An overlapping additive Schwarz over the **static resting load-path spring graph** — the union of actin
**crosslinks** (`k~8e5`), membrane↔cortex **ERM** tethers (`k_erm=4600`), and membrane **area-tension** edges
(`γ_mem`). Each edge owns a 2-node 6×6 block `[[D_i,-H],[-H,D_j]]` (`D` = regularizer + every spring tangent at
the node; `H` = that edge's tangent), degree^-1/2-weighted → SPD; the one-sweep direction `P M⁻¹ PF` is a
projected descent for the driver's exact residual-monotone line search. Mirrors `ContactPairSchwarz`, extended
to solve crosslinks as their **own** blocks (a first ERM-only-pairs attempt that folded crosslinks into the
diagonal detonated the 8e5 spring — two crosslinked cortex nodes reached by different ERM blocks stepped
inconsistently; recorded so it is not re-chased).

## Result — a real but PARTIAL improvement; the 0.21 gate is NOT closed

| solver (subdiv 6, ERM-preload, `n_inner=200`) | 8k `residual_candidate` | native 70,686 `residual_candidate` |
|---|---|---|
| `analytic_implicit` (Newton) | 2.166 | 0.783 |
| `erm_schwarz` (standalone) | 2.166 | 0.783 |
| **`erm_tournament`** (fiber-block + erm_schwarz + anderson + rkc) | **1.374** | **0.620** |

`erm_tournament` improves the native plateau **0.776 → 0.620** — the coupled membrane-ERM-cortex blocks capture
conditioning the Newton misses. But **it does not reach the 0.21 gate.** The *isolated* `erm_schwarz` direction
is actually an **ascent** direction for the nonlinear residual (probe: one sweep sends max|PF| 0.80 → 172 at
`t=1` and no line-search fraction descends) — the gain comes only from the tournament combining it with the
per-fiber block. **Gate stays OPEN honestly (no loosening).**

## Refined diagnosis — why a spring-pair Schwarz is insufficient (two residual gaps)

**(a) The resting cortex is BACKBONE-inextensibility-dominated + sparsely crosslinked.** Measured at the
preloaded state: the ERM-attached cortex nodes have a diagonal trace of **7608 = 3a(3000) + one ERM(4600)** —
i.e. the *median* ERM-cortex node carries **NO crosslink** (`n_xl=8000` for ~72k actin nodes; ~0.45
crosslinks per ERM-cortex node). The stiffness resisting a cortex node's motion is therefore mostly the
**projected backbone inextensibility (P) + bending**, not a local crosslink neighbourhood. A spring-pair
Schwarz cannot represent the backbone constraint, so its block direction over-steps. **Correct preconditioner:
BACKBONE-AWARE — extend the existing per-fiber Cholesky block (which already captures the backbone) to include
the ERM-tethered membrane nodes**, solving membrane→ERM→cortex→backbone consistently per fiber. This is a
different, larger design than the `ContactPairSchwarz` mirror.

**(b) A separate membrane floor.** The membrane carries an intrinsic ~**0.38 pN** (native) tangential
mesh/area-tension residual that neither the ERM preload nor the cortex-focused step clears; it sits above the
0.21 gate on its own. Needs a membrane tangential relaxation and ties to **Card 3** (`ΔP=2γ/R` consistency,
Hosseini γ).

## Recommendation (surfaced to PI)
1. **Cortex:** backbone-aware fiber-block extended with ERM-tethered membrane DOFs (not a spring-pair Schwarz).
2. **Membrane:** tangential area-tension relaxation to clear the ~0.38 floor.
`erm_schwarz`/`erm_tournament` are kept as **additive, documented options** (defaults unchanged) — they help
(0.776→0.62) and are the harness the next attempt builds on, but they are **not** a gate closure.

## Artifacts
- Module: `aleph/components/incumbent/erm_schwarz.py`; probes `scripts/ac_erm_schwarz_probe.py`, `ac_erm_solver_sweep.py`.
- Data: `outputs/ac/implicit/erm_schwarz_{native,8k}_subdiv6.json`, `erm_sweep_{native,8k}.log`.
- Reproduce: `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim python scripts/ac_erm_solver_sweep.py --native`.
