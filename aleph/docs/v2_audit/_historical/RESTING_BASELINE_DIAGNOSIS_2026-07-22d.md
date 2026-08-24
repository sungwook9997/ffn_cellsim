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

# Resting-baseline diagnosis (2026-07-22d) — FD-vs-operator probe: the operator is EXACT; the blocker is conditioning, not assembly

Branch `codex/ff-ac-codex`. Executes the PI P0 mandate (AC_DECISION_CARDS_2026-07-22.md Card 4): the
finite-difference-vs-analytic-operator comparison, run before any new solver/preconditioner code. Probe:
`aleph/scripts/ac_fd_operator_probe.py`. Raw artifacts: `outputs/ac/implicit/fdprobe_{native,8k}_subdiv6.json`,
`fdprobe_native_subdiv6.log`. Run on the gbook A5000 at the fixed state **subdiv 6 + overlap_free + radial ERM
+ ERM preload**, frozen state/active-set. **Confirmed at FULL NATIVE (70,686 filaments, n_total = 551,434)**
per the PI native-debug rule — and native is materially WORSE than 8k (below), so 8k under-represented it.

## The decisive result — every prior "operator gap" root cause is REFUTED

| probe question | measurement (native 70,686) | verdict |
|---|---|---|
| (a) per-family analytic tangent `K·v` vs true-force FD `-[F(x+εv)−F(x−εv)]/2ε` | bending rel **9.9e-6**, crosslink rel **9.9e-6**, ERM rel **7.0e-6** (ε=1e-8) | tangents are the **EXACT** Hessians of the true force kernels |
| operator stiffness in the residual direction `s = v·(Kv)` | `s_operator=57537`, `s_fd_true=57537`, **ratio 1.0000** | operator delivers **exactly** the true measured stiffness |
| (b) `‖PF‖` vs `‖F‖` (residual in projection nullspace?) | ratio **0.9946** (PF≈F) | residual is **NOT** in the projector nullspace; P is not eating it |
| crosslink compression / clamped negative curvature | frac_compressed 0.000, `|ext|≈0` | no indefiniteness from the SPD transverse clamp |

**This refutes `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md` §Update 2's confirmed lead** ("the analytic
Jacobian `PKP` does not deliver the measured stiffness into the step for the turgor/ERM-reaction mode"). It
does deliver it — to 5 significant figures (`ratio_op/fd = 1.0000`), at full native, per force family. The
`add_pair_stiffness_kernel`/`add_bending_stiffness_kernel`/`add_erm_stiffness_kernel` tangents ARE the exact
linearizations of `link_spring_kernel`/`cytosim_bending_kernel`/`erm_tether_force_kernel` (verified by code:
same `kxl_d/r0xl_d/alpha_d` arrays, quadratic bending ⇒ FD = analytic exactly). **There is no first
mismatching force family — there is no mismatching family at all.** The ERM-tangent-zero-at-rest lead
(46aee6de) is likewise not the blocker: with the preload the ERM is stretched (`ext>0`) so its tangent is
present, and it matches FD to 7e-6.

## What IS the blocker — the exact operator produces an OVERSHOOTING, ill-conditioned Newton step

The step `dx = (aI + PKP)⁻¹ PF` (RHS `PF`, exactly as the driver forms it, `a=1000`):

| quantity | 8k | native 70,686 |
|---|---|---|
| residual after ERM preload (max\|PF\|) | 2.21 | **0.800** |
| ONE full Newton step (t=1): max\|PF\| after | 2.21 → **4.45** (2.0×) | 0.80 → **6.59** (8.2×) |
| best line-search fraction `t*` and its residual | t=0.5 → 1.39 (37% ↓) | **t=0.0625 → 0.778 (3% ↓)** |
| cortex localized stiffness `s` (matches FD) | 20,458 | 57,537 |

The full Newton step **overshoots** (native by 8×), and the damped-Newton REFERENCE — the exact operator +
a proper backtracking line search over `t∈{1,½,…}` (a diagnostic reference, not production solver) —
**STALLS**:

```
native:  step0 0.800→t=0.0625→0.778   step1 0.778→t=0.031→0.776   step2 STALLED (best_t=0)   gate=0.21
8k:      step0 2.21 →t=0.5   →1.39    …  step8 0.919 STALLED
```

Native stalls at **0.776 pN — 3.7× the 0.21 gate** — after moving essentially nothing (max drift 0.046 nm).
The stalled max is on the **cortex** (0.776; membrane 0.377, other-actin 0.163). No crosslink compression of
consequence (`|ext|≈0` even at the stall).

### Mechanism — density-scaling conditioning of the coupled soft-membrane / stiff-cortex system

The severity scales with the physical density: 8k overshoots 2× and descends to 0.92; native overshoots 8×
and is stuck at 0.78 after one effective step. The difference is the cortex stiffness (20,458 → 57,537) — the
denser the true cortex, the larger the stiffness ratio between the stiff crosslinked cortex and the soft
membrane/ERM manifold (`k_erm=4600`, area-tension = 0 transverse at rest, membrane bending parked in the
scalar `a`), and the worse `A⁻¹` conditions. `A = PKP + aI` is SPD and locally tangent-exact, but `PKP` is
near-singular on the soft coupled manifold; the `a=1000` Tikhonov floor makes the soft-mode step
`≈PF/a ≈ sub-nm` which, coupled through the ERM/crosslinks into the stiff cortex, **overshoots** the stiff
residual. A single scalar line-search cannot rescale per-mode ⇒ stall. **This is exactly the "stiff soft/stiff
coupling" 22c named — but it is a conditioning problem of an EXACT operator, not a missing/wrong stiffness.**
It is invisible at coarse density, which is why the PI native-debug rule was decisive here.

## Fork decision (the Card-4 P0 question: conditioned-implicit vs targeted per-node preload)

Both naïve shortcuts are now eliminated by measurement:
- **Targeted rest-length preload (fork 1)** reduces to `K δr = −F` with the **same** ill-conditioning (22c
  §Update: stuck at 2.22 / diverges under Jacobi). Not fixable by the preload itself.
- **Conditioned-implicit as scalar-damped Newton (fork 2, cheap form)** STALLS at 0.78 native (this doc).

⇒ **The operator is proven correct, so the justified path is fork 2 in its REAL form: a preconditioner that
captures the coupled membrane↔ERM↔cortex soft/stiff modes.** The existing pieces do NOT: the per-fiber block
Cholesky is per-fiber only (no ERM coupling), node-Jacobi is diagonal, and the `l≤2` rigid/strain coarse space
was already shown useless (22c §Update 2 — the residual is not a global `l≤2` mode). Leading concrete,
bounded candidate: an **ERM-pair (and cortex-crosslink-neighbourhood) Schwarz block preconditioner**, mirroring
the existing `ContactPairSchwarz` that already does exact small-block solves on overlapping WCA pairs. This is
a **solver-architecture change, not a minimal family fix**, so it is surfaced to PI for sign-off rather than
landed unilaterally under the "minimal fix only" clause.

**Gate status: OPEN, honestly.** No minimal fix closes it because the FD probe proves there is no family to
fix; the residual is real and no gate was loosened. The next investment is the coupled-block preconditioner
(PI-gated). Secondary open item surfaced by the probe: the membrane carries an intrinsic ~0.38 pN residual
(native) that the radial ERM preload only reduces ~20% and that the cortex-focused step leaves inert — a
membrane initialization/Young-Laplace-consistency question (ties to Card 3 `ΔP=2γ/R`).

## Reproduce
```
# gbook A5000
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
  ~/miniconda3/envs/ffn_sim/bin/python aleph/scripts/ac_fd_operator_probe.py --native --subdiv 6 \
  --json /tmp/fdprobe_native.json
```
