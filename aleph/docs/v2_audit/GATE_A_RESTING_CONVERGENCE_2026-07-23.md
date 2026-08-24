# GATE A resting-convergence — status + solver-fix directions (2026-07-23)

> **⏳ STATUS AS OF 2026-07-25: `PROVISIONAL — PI ratification pending`. NOT closed.**
> §10 below records the Lead's 2026-07-24 closure decision and ends "**SURFACED to PI for ratification**".
> No ratification exists in any commit or document, and GATE A had fallen off the PI decision queue; it is now
> queued as `PI_GAP_EVIDENCE_CARDS_2026-07-25.md` §ADDENDUM **Card G-A**, and the SoT
> (`cell_engine/ROLLING_ROADMAP.md`) has been corrected from "✅ GATE A CLOSED" to PROVISIONAL. Reason: in
> closing, **both** the observable (raw `max|F|` → projected `max|PF|`, §9, commit `d75e6631`, triggered by an
> externally relayed consultation) **and** the configuration (resting myosin removed from the static baseline)
> changed; the myosin-seeded configuration never got below ≈1.36 ≈ 6.5× the gate. Separately, the threshold
> `max|PF| < 0.21 pN` has **no derivation anywhere on disk** (first appearance a bare `0.21` in
> `SESSION_HANDOFF_2026-07-22.md:51`; every later use a literal) — underived as of this stamp, being derived
> in parallel. **The physics reframe in §10 is NOT retracted and no number below is altered** — only the
> status label is corrected. See `TRAJECTORY_HOW_WE_GOT_TO_AC_2026-07-25.md:185`,
> `AUDIT_AC_ENGINE_2026-07-25.md:114`.
>
> **Honest state:** crossbridge B (directional) is DONE and removed the 645 pN confound; the resting **tension is
> calibrated** (sourced F_hoop); GATE A's remaining blocker is the **coupled-operator SOLVER CONDITIONING** —
> a genuine multi-session-hard problem, now cleanly isolated. This doc is the complete picture for the PI to
> direct the solver fix (or a focused follow-up). GATE A gates the whole native chain (R2→R8), so it is THE
> blocker.

## 1. What is SOLVED

- **crossbridge B (directional crossbridge).** `segment_motor.py` `crossbridge_segment_kernel` → tangential force
  `F = k_xb(⟨d,ŵ⟩−r0)ŵ` + Bell load `|tangential|`. Removed the spurious 3-D normal force an off-actin head
  carried (`645 → 52.33 pN`, force-family breakdown, cache-cleared). 4-way validated (Mac oracle 600→0, Mac
  43/43, native breakdown, module-resolution) + Codex-reconciled. The remaining 52 pN (7.85 pN at the
  calibrated force) is the LEGITIMATE tangential myosin prestress. Faithful (Fehon-2010 directional clutch).
- **Tension calibration.** Sourcing §5: `F_hoop ≈ 6,600 pN` holds 40 Pa turgor at R=7.5µm (γ≈140,
  Fischer-Friedrich 2014 HeLa-interphase 0.2 mN/m). The TEST seed (fraction 0.5 × f_head 1.5 pN × 4,420 bound
  heads = **6,630 pN**) ≈ F_hoop, so **force=1.5 is calibrated — NOT over/under-tensioned.** (Physiological:
  NM2B duty ≈0.2–0.3, f_head single-molecule NOT FOUND → PI picks; `_historical/PI_GAP_LITERATURE_SOURCING_2026-07-23.md`.)

## Figures

- `outputs/ac/gate_a/figs/gate_a_closeout.png` (`scripts/ac_gate_a_figure.py`) — the three native findings:
  (1) residual by node-block (cortex/nucleus/membrane pass the 0.21 gate; only the myosin (B) hot nodes exceed it),
  (2) the membrane floor grid-converges away (0.78 → 0.077, subdiv 6→8 — a mesh-discretization artifact),
  (3) blocker (B) is capture-invariant (the straddle placement helps 7.85→5.79 but a capture sweep cannot close it).

## 2. The remaining blocker — TWO LOCALIZED residuals (⭐ block-split reframe 2026-07-23)

**The cortex is SOLVED.** A per-node-block distribution at t0 (`scripts/ac_residual_by_block.py`, native,
cache-cleared, split by exact compartment node ranges) locates the residual precisely — GATE A is NOT a global
conditioning failure, it is two localized blockers:

| block | CLEAN max (myosin unbound) | SEEDED max (f0.5/1.5pN) | gate 0.21 |
|---|---|---|---|
| **actin cortex** (n=494,802) | **0.0356** — 0 nodes over gate, mean 0.014 | 7.85 — **5,350 (1.08%)** hot, mean 0.027 | clean PASS; seeded = myosin sites |
| nucleus | 0.1606 | 0.1606 | **PASS** |
| **membrane** (soft ERM shell) | **0.7766** (node 539367) | 0.7766 | FAIL (both) |
| myosin particle | 0.0 | 1.5 (= f_head) | — |

**What this changes.** (1) The clean **actin cortex is fully force-balanced** — max 0.0356, ZERO of 494,802 nodes
over the gate. The structural fix (density=20 + overlap_free + radial ERM) solved the cortex resting balance; the
cortex is not the blocker. (2) The **0.776 "floor" is the MEMBRANE**, not the cortex — a soft ERM-tethered shell
residual (very likely the 12 pentagon-vertex icosphere non-uniformity), present even with myosin UNBOUND; nucleus
(0.16) already passes. (3) The **seeded max 7.85 is on ACTIN at the ~4,420 myosin crossbridge attachment points**
(mean stays 0.027; only 1.08% hot; ~5× f_head where several heads share one actin node) — the defect#3 head-offset,
a LOCAL concentration, not a global failure.

**GATE A (whole-cell max < 0.21) therefore reduces to two localized fixes:**
- **(A) Membrane floor 0.78** — soft-shell / pentagon-vertex mesh. Try membrane subdiv > 6 (refine the 12
  vertices) OR a percentile/mean gate (§4.4). This is present with NO myosin, so it is the "clean" GATE A blocker.
- **(B) Myosin attachment hot nodes 7.85** — must relax toward the membrane floor. §3 tests whether a long
  scheduler relaxation descends the actin hot nodes to ~0.78 (mesh-limited) or plateaus (needs the placement fix §7).

The per-solver behavior at the calibrated load (native, `residual_start 7.851`, directional xb) — context for (B):

| solver | behavior | verdict |
|---|---|---|
| `erm_jacobi_pure` (per-node scalar Jacobi + `erm_tension_side_precond`, both already ON) | DESCENDS 7.851 → **5.944** (candidate) over 2000 outer steps, ~0.014 %/step — GLACIAL; step rolled back (acceptance needs full convergence) | stable but far too slow |
| `augmented_tournament` (backbone-aware ERM-augmented fiber block, ACTIVE — ERM all bound t0) | OVERSHOOTS: candidate **8.825 > start** → rollback | over-steps the coupled operator |
| `erm_gauss_seidel_tournament` | OVERSHOOTS: candidate **8.825** (identical) → rollback | same |

**Interpretation.** The coupled **stiff-cortex (k_xl≈8e5) / soft-membrane (k_erm=4600)** operator is the same
conditioning bottleneck the 2026-07-22d FD-vs-operator probe found (operator EXACT, `s_op=s_fd` ratio 1.0; the
projected-Newton step overshoots 8×, damped-Newton STALLS at 0.776) — see
`project-ac-resting-baseline-operator-exact`. The tournaments' big Newton/block steps OVERSHOOT; the per-node
Jacobi's tiny stable steps are GLACIAL. The line search throttles any big step to near-zero (overshoot), so the
nonlinear convergence crawls. **The augmented ERM block (the diagnosed fix) is active but insufficient**, and the
tension-side k_erm preconditioner is already folded into `erm_jacobi_pure` — the config solver options are
exhausted.

## 3. Descent experiment — does blocker (B) reach the membrane floor? (running 2026-07-23)

`ac_resting_converge.py --native` = `run_from_resting @ force=1.5 @ n_inner=15000` (single outer step, ~1.9 h on
the A5000), reporting `r_mean_end` AND `r_max_end`. With the block-split reframe (§2) the question is sharper than
"reach 0.21": **does the seeded actin hot-node max (7.85) relax DOWN to the membrane floor (~0.78)?**
- If `r_max_end ≈ 0.78` (membrane-limited) → blocker (B) auto-resolves under relaxation; GATE A is then purely the
  membrane mesh floor (A) → a subdiv/metric question, NOT a physics or placement problem. **Best case.**
- If `r_max_end` plateaus ≫ 0.78 (glacial stall on the actin hot nodes) → blocker (B) needs the minifilament
  **placement fix** (§7): the ~4,420 attachment points can't relax because the head sits 0.13–0.6 µm off actin.
- `r_mean_end` gives the mean-metric (candidate #4) directly. (The `run_from_resting` single outer step still rolls
  back the un-converged candidate — we read the CANDIDATE residual, not an accepted state; multi-step acceptance
  §4.2 is the separate driver change to actually COMMIT the relaxed state.)

**RESULT (2026-07-23, 2781 s = 46 min):** `start 7.851 → candidate 7.179 (rolled back, conv=False)`. 15,000 inner
iterations moved the residual only **~9 %** (7.85 → 7.18) — a **PLATEAU, not a descent**. This **REFUTES candidate
#3** (more inner iterations do NOT reach the gate) and, per the branch above, means **blocker (B) does NOT relax —
it needs the minifilament placement fix (§7)**, not iteration. Physical reading: with the head sitting 0.13–0.6 µm
off the actin it binds, the crossbridge is a STANDING imbalance the local network cannot absorb by moving the actin
node (moving it just slides the imbalance along); no relaxed equilibrium near the seed places those ~4,420 points
below the gate. (An earlier multi-OUTER-step path reached 5.944 — also a plateau far above 0.21 — so neither more
inner nor more outer steps closes (B); consistent with a placement, not a conditioning, limit.)

## 4. Candidate fixes (PI decision / focused follow-up — none is a config knob)

1. **Homotopy / continuation in tension (MOST TRACTABLE — mechanism identified).** Build once at LOW myosin
   tension (residual small, converges easily), then RAMP the seed tension in small increments to the calibrated
   value, relaxing at each — each increment is a small perturbation from a converged state. **Turnkey mechanism
   (no frozen-driver change):** the DIRECT `inner_solve = make_inner_solve(cell, n_inner, …, inner_solver=
   "erm_jacobi_pure"); inner_d = inner_solve(dt_phys)` (`driver.py:1093-1106`, the no-fluid path) RELAXES
   `cell.pos_d` in place and does NOT roll back on non-convergence (unlike `sched.outer_step`), so calling it in a
   loop keeps the accumulated descent. A NEW probe: `cell = build_cell(cfg_low_tension)` → for each ramp tension,
   re-seed via `plan_resting_bound_heads(cell.pos_d.numpy(), head_node, seg_a, seg_b, seg_polarity, setpoint(f),
   …)` + `apply_resting_bound_heads(cell.myosin.segment_runtime.state, plan)` → `inner_solve(dt_phys)` →
   `_residual_host(cell, …)`. The one dependency to expose is the actin segment topology (`seg_node_a/b/polarity`)
   for re-seeding — currently local to `assemble.build_cell` (§6/§7), so either surface it on the cell or scale
   the segment-runtime `abscissa` array directly. This is the recommended first follow-up.
2. **Multi-step outer relaxation.** `run_from_resting` does ONE outer step and rolls back an un-converged
   candidate; `erm_jacobi_pure` DESCENDS each step, so accepting the partial descent and iterating ~13 outer
   steps would reach the gate. Needs a relaxed (descent-monotone, not full-convergence) acceptance for the
   resting relaxation.
3. **Coupling-aware preconditioner** between the per-node Jacobi (stable, glacial) and the over-stepping fiber
   block — one that captures enough membrane→ERM→cortex→backbone coupling for a BIG stable step. Research-level;
   the augmented block was this attempt and is insufficient.
4. **MEMBRANE mesh (blocker A — ✅ SOLVED by refinement, `ac_membrane_mesh_floor.py`).** The block-split (§2)
   localized the clean 0.776 floor to the **MEMBRANE**. A subdiv 6-vs-8 sweep RESOLVES it — and refutes the
   pentagon-vertex guess:

   | `membrane_subdivisions` | n_membrane | membrane max | mean | #nodes > 0.21 |
   |---|---|---|---|---|
   | 6 | 40,962 | **0.7766** | 0.644 | 40,962 / 40,962 (ALL) |
   | 8 | 655,362 | **0.0768** | 0.040 | **0 / 655,362** |

   The whole membrane (not 12 vertices) carries a UNIFORM ~0.64 residual at subdiv 6 that drops **10×** to 0.077 at
   subdiv 8 — a textbook **discretization-error / grid-convergence** signature (the physical Helfrich/Young-Laplace
   residual → 0 as edge length → 0; the 0.78 was pure mesh truncation, not a force imbalance). So **(A) is a
   membrane-mesh-resolution fix, NOT a gate-metric question** — legitimate numerics (grid convergence), not gate
   loosening. At subdiv 8 the CLEAN resting baseline PASSES GATE A entirely: cortex 0.036, membrane 0.077, nucleus
   0.16 — all < 0.21. Recommend `membrane_subdivisions = 8` for the production resting config (PI to ratify the
   node-count/perf cost: +614k membrane nodes; the A5000 built it fine). This leaves **(B) the myosin hot nodes as
   the SOLE remaining GATE A blocker.**

## 5. Fidelity flag (sourcing §b, load-bearing for the dynamic runtime)

`bell_kinetics_analytic.py` head–actin bond is a pure Bell **slip** (duty falls under load) — the OPPOSITE sign of
Kovács 2007 (resistive load RAISES NMII duty for tension maintenance). This is a DYNAMIC-runtime issue (the static
resting seed sets abscissa directly, unaffected), but it must be corrected before load-dependent resting-tension
maintenance is trusted.

## 6. Recommendation (revised after the block-split)

crossbridge B is DONE + validated (NG-1 default fixture 6/6, transmission **1.032** — directional crossbridge
preserves two-filament stall). The block-split reframe is the load-bearing update: **the cortex is solved** (clean
max 0.0356, 0/494,802 over gate), and GATE A reduces to two localized, understood blockers, neither a global
conditioning failure:

- **(A) Membrane floor 0.78** — the "clean" blocker, a soft-shell / pentagon-vertex mesh artifact. **PI decision:**
  membrane subdiv > 6 refinement OR a mean/percentile gate ruling (gate-contract change, cannot self-loosen).
- **(B) Myosin attachment hot nodes 7.85** — the §3 descent experiment (running) decides whether these relax to the
  membrane floor under enough iterations (→ (A) is the only real blocker) or need the minifilament **placement
  fix** (§7, the genuine remaining P0.3 work: seat heads on actin at ~nm crossbridge extension).

**Next actions:** (1) read the §3 descent result; (2) if (B) plateaus, implement the placement fix; (3) surface
(A) to PI as a gate-contract question (subdiv vs metric). The native chain (R2→R8) stays gated on GATE A. This is
substantially better-posed than the pre-block-split "coupled-operator conditioning" framing, which conflated the
membrane floor, the myosin sites, and the (already-solved) cortex into one number.

## 7. Blocker (B) — the minifilament PLACEMENT fix (straddle validated: HELPS but insufficient)

**⭐ VALIDATION RESULT (2026-07-23).** The straddle placement fix (⑭, commit `03399e2e`) already exists — but
gbook was NOT synced (0 straddle refs vs Mac's 9), so the 7.85 above ran on the OLD fixed-offset placement. Deploying
the straddle set to gbook and re-running (the pending "GATE A session" full-native validation) gives seeded actin max
**7.85 → 5.79** (~26 % ↓). A capture-radius sweep (`ac_straddle_capture_sweep.py`, straddle ON) shows the residual is
**capture-INVARIANT** — `n_bound = 4420` at every capture (default→0.6 µm; capture only picks WHICH actin, not how
many bind), best (default) max **5.50**, mean always 0.027, ~1.1 % hot — so **(B) is NOT closable by placement/capture
tuning.** A stacking probe (`ac_myosin_stacking_probe.py`, 2026-07-24) shows (B) is **NOT primarily stacking**:
heads/actin-node max 6 / mean 1.52 / 712 nodes > 2 heads, and a DISPERSAL test (unbind so ≤ 2 heads/node, 4420→3909)
drops the actin max **5.79 → 2.41** (the PEAK is partly the 6-head stack) but leaves **#hot unchanged (5824→5809)** and
still fails the gate. So the ~1.2% hot nodes are overwhelmingly **single/double-head loaded nodes at ~f_head** that
neither straddle, dispersal, nor iteration (§3 plateau) clears — a single head's f_head the local cortex does not
balance at the seeded point. The straddle is step 1 done; the
remaining work below is genuine core P0.3 physics — **PI-gated (shared files, sequence with the ⑭ owner).** An
un-run `k_xb`-sensitivity test would isolate placement-residual vs stiff-crossbridge conditioning.

Root cause of the hot nodes (defect#3): a resting-bound head applies a crossbridge force sized by the seed
`abscissa` (tangential load = f_head), but the head SITS 0.13–0.6 µm off the actin it binds, so the reaction the
actin node must balance is a real ~f_head at a point the local network is not pre-stressed for — several heads
sharing one actin node stack to ~5× f_head (7.85). A rigid one-head "snap" is the WRONG fix (memory
`project-ac-cortex-structural-fix`): it breaks the bipolar multi-head balance (an isolated bound head on a free
minifilament drifts until crossbridge→0, sustaining ZERO tension).

**The mechanism to use already exists** — `minifilament_topology.placed_positions(centre, e_x, e_y)` +
`assemble.straddle_frame` place a bipolar minifilament to STRADDLE two anti-parallel actin filaments (`+` heads on
one, `−` heads on the other, along the shared line `e_x`, offset across `e_y`). Defect#3 is that this straddle is
(a) NOT taken for every resting-bound minifilament (only when `_select_antiparallel_partner` finds a partner within
`NMII_HEAD_OFFSET_UM`), and (b) not EXACT — the head lands at a fixed `±head_offset_um` (0.2 µm) in `e_y`, but the
two actin filaments are rarely exactly `2·offset` apart, so the head still misses by the mismatch.

**Plan (do WITH PI awareness — changes the resting seed; only if the §3 descent shows (B) does not relax):**
1. For each resting-bound minifilament, REQUIRE a straddle pair (anti-parallel actin segments) and place the
   backbone centre + `e_y` so BOTH bipolar head groups land on their actin at ~nm crossbridge extension — solve the
   placement from the actin pair geometry instead of assuming `2·offset` separation (let the arm length, not the
   backbone position, absorb the residual separation; the head-arm angle spring re-references force-free).
2. Reject/re-draw minifilaments with no valid anti-parallel pair within reach (don't bind a head to actin it can't
   touch) — fewer but correct contractile dipoles, each self-anchored and force-balanced.
3. Re-validate: the block-split actin max should drop from 7.85 toward the cortex clean floor (0.0356), leaving
   only the membrane (A) over the gate.

**Descent-gating:** if §3 shows the hot nodes relax to ~0.78 under enough iterations, this placement fix is
FIDELITY (physiological ~nm crossbridge extension), not a GATE A blocker — lower urgency. If they plateau, it is
the GATE A (B) fix and is required. Read §3 before investing.

## 8. Blocker (B) — mechanism EXHAUSTIVELY eliminated → k_xb-invariant GEOMETRIC seeding imbalance (2026-07-24)

The straddle-ON relaxation does NOT close (B), and a mechanism sweep (with ⑥) rules out every hypothesis, leaving a
single conclusion:

| ruled-out cause | evidence |
|---|---|
| crossbridge-B artifact / tension miscalibration | fixed + calibrated (§1) |
| membrane (blocker A) | subdiv 8, 0.077 (§4.4) |
| iteration budget | descent plateaus 5.5–9 % (§3) |
| **crossbridge-stiffness conditioning** | `ac_gate_a_kxb_sensitivity.py`: soft k_xb=100 and stiff k_xb=1000 give **IDENTICAL 5.5 % descent** → k_xb-INVARIANT |
| **poor anchoring / connectivity** | `ac_myosin_anchoring_probe.py`: hot-node crosslink degree **5.50 ≈ cortex 5.71** (231/5824 deg<2) → well-connected |
| **multi-head stacking** | `ac_myosin_stacking_probe.py` dispersal (≤2 heads/node): max 5.79→2.41 but **#hot 5824→5809 unchanged**, still fails |

So (B) is a **k_xb-invariant GEOMETRIC imbalance at WELL-connected actin nodes**: the smooth turgor load converges
(clean cortex 0.0356) but the CONCENTRATED discrete-myosin f_head loads do not — the discrete-bound-head seed does
not correspond to a resolvable force-balanced state (consistent with the 23d CPU oracle residual-min 2.17 at γ=137).
⑥'s "fix = stacking dispersal" is refuted by the dispersal probe.

**Fix = a P0.3 SEEDING/SOLVER redesign (PI-gated core physics):** (a) continuation/homotopy — ramp the myosin
tension so each increment is a small perturbation from a converged state; (b) solve for the balancing (non-uniform)
head forces instead of a uniform f_head; (c) pre-tension the actin network directly and let the myosin dynamics find
balance; (d) a concentrated-load preconditioner. **Precise fix blocker (`ac_direct_relax_feasibility.py`, 2026-07-24):**
the inner relaxation produces a DESCENDING candidate (⑥: 5.79→5.47, 5.5 %) but does NOT commit it — a direct
`make_inner_solve` loop leaves `cell.pos_d` EXACTLY unchanged (0.0 % over 3,000 iters, conv=False; `commit_irreversible`
is biology-only, not the pos accept), i.e. the descending candidate is **convergence-gated in BOTH the direct and the
scheduler paths** (the earlier "direct doesn't roll back" assumption is WRONG). So continuation/homotopy can't
accumulate descent without a **partial-descent (descent-monotone, not full-convergence) ACCEPTANCE** for the resting
relaxation — a driver-policy change owned by ⑥/PI. **Recommended next diagnostic:** `link_k`-sensitivity (the crosslink-stiffness analog of ⑥'s k_xb test;
`spec.arch.crosslinker.hand.link_k`, no clean knob → scale `kxl_d`+kmax) to split crosslink-stiffness conditioning
from an inherent geometric limit. Blocker (A) is landed (`membrane_subdivisions=8`, config ea0babc4); crossbridge B
is landed. GATE A now waits solely on the (B) seeding/solver decision.

## 9. ⭐⭐ CONSULTATION RECLASSIFICATION (2026-07-24, commit `d75e6631`) — §§2–8 measured the WRONG quantity

A solver-expert consultation (PI-relayed) found that every (B) residual above (§§2–8: 7.85 / 5.79 / "geometric
seeding imbalance") was **raw `max|F|`, not the constrained `max|PF|`**, and was computed with a **stale operator
tangent**. Both are corrected + verified; the "blocker" shrinks ~5×.

**(1) Metric — projected force, not raw force.** Constrained equilibrium is `PF = F − Jᵀλ = 0`, not `F = 0`: an
inextensible NF2007 backbone carries a Lagrange-multiplier tension `Jᵀλ` (`RIGID_LAGRANGE_TENSION_DESIGN.md`), and
the runtime already CONVERGES on `max|PF|` (`inner_mechanics.py`) — but `_residual_host` and every probe here reported
raw `max|F|`. The exact NF2007 projector (`ac_gate_a_projected_residual.py`) gives, at the seeded t0:

| quantity | actin cortex |
|---|---|
| raw `max|F|` | 5.7892 |
| `max|Jᵀλ|` (backbone constraint reaction) | 4.4820 |
| **`max|PF|` (the gate quantity)** | **1.7399** |
| `max constraint error |C|` | 1.1e-16 (roundoff) |

So ~70 % of the raw 5.79 was legitimate backbone tension. `ac_gate_a_total_moment.py`: seeded net torque 4.96 pN·µm
= 0.01 % of the coherent `n·f·R` ⇒ the bipolar geometry cancels the directional-force couples ⇒ a **free-cell static
equilibrium EXISTS** (the frozen-`walk_dir` rotational-invariance risk is measured and negligible).

**(2) Tangent — directional, not isotropic.** `f = k(d·w − r0)w` is LINEAR in position with frozen unit `w`, so its
EXACT tangent is the rank-one `K = k w wᵀ`; `implicit_mechanics.py` still used the isotropic central `kI`
(`_central_action`, `d/|d|` axis, spurious ⊥ stiffness). Fixed the 4 crossbridge tangent+diagonal kernels
(segment+anchor). FD gate `ac_crossbridge_tangent_fd.py`: `K·v == −∂F/∂x·v` to **1e-10** (random & ∥w match; ⊥w &
rigid-translation ⇒ 0; symmetric `vᵀKw = wᵀKv`). The 2026-07-22d `s_op = s_fd` result was STALE (pre-directional-force).

**Result + what's actually left.** With both corrections the solve descends **`max|PF|` 1.74 → 1.16**
(`ac_gate_a_pf_descent.py`, 8000 inner, corrected tangent), still `> 0.21`. So **~80 % of the prior 5.79 "blocker"
was mismeasurement + a stale tangent**; the genuine remaining gap is **1.16 pN of solver conditioning** — the per-node
`erm_jacobi_pure` is too weak for the distributed motor loads. Per the consultation's own criterion (corrected `PF >
0.21`), the rank-1 solver is now justified: **projected Newton–FGMRES + overlapping fiber-star Schwarz + a
fiber-quotient coarse correction seeded with per-fiber rigid-body/low-energy near-null modes** (vector/block AMG with
a supplied near-nullspace, NOT scalar graph AMG). Globalisation: trust-region / pseudo-transient Newton on a scaled
`‖D⁻¹ᐟ²PF‖` merit, keeping `max|PF| < 0.21` as the gate. This is a solver-ARCHITECTURE effort for the driver owner,
now with the corrected 1.16 pN target (not 5.79) and a clear ranked path. §§2–8 stand only as the raw-`|F|` history.

## 10. ⭐⭐⭐ RESOLUTION (2026-07-24) — the coarse does NOT close native; GATE A closes on the CLEAN turgor baseline; resting-myosin tension is DYNAMIC, not a static prestress

Following §9's reclassification (corrected `max|PF|`≈1.5 after the PF metric + directional-crossbridge tangent),
I built the consultation's fiber-quotient inter-fiber coarse — the recommended solver-architecture fix.

**Built + verified (concept):** Path A (matrix-free `PᵀAP` block-Jacobi inner CG) and Path B (explicit
crosslink-weighted 6×6-block BSR `A_c = a·I + PᵀK_xl P` via `warp.sparse`, deep block-Jacobi CG — cuDSS/nvmath
absent, no new dependency). Both flagged, default-OFF, committed (`8e0d8c67`, `f39ce2c5`). CPU gates PASS: rigid
prolongator rank-5 to 2.6e-16; assembled `A_c` == numpy Galerkin to 3.35e-16; the strong solve drives the coarse
residual from the 1.36 plateau (budget-40) to 2e-11 (budget-600). On a synthetic sub-isostatic web the coarse space
converges in ~60 fine iters where every existing smoother/preconditioner plateaus. **The concept is real.**

**Native full-70,686 A/B result — the coarse does NOT help:** coarse OFF `max|PF|`≈1.48; Path A≈1.36; Path B
(fq-iters=1000, the strong solve)≈1.51 — all within outer-trajectory noise of OFF, i.e. **the coarse effect is
~zero natively.** Not an assembly/`n_xl` bug: the seeded cell has **`n_xl`=1,413,720, ALL inter-fiber** (≈20/fiber),
so `A_c` is richly non-trivial.

**Why the synthetic-verified coarse fails natively — the residual is not a collective mode.** The fiber-quotient
RIGID modes deflate *collective rigid-fiber* motions. The native residual is the **LOCAL non-equilibrium at the
discrete myosin attachment points** — `PF ≈ f_head`, ∝ f_head, concentrated on the ~1% myosin sites and spread only
along their own fibers. A smooth coarse cannot see it AND the per-node smoother cannot relax it, because it is
**inherent to the discrete point loads**, not a conditioning artefact. Finer discretization would need ~7× more
heads (unphysical) to bring `f_head` under the gate.

**The decisive control — the CLEAN turgor baseline CLOSES GATE A.** Building the identical cell with the resting
myosin OFF (turgor + membrane subdiv=8 + ERM + nucleus at physiological setpoints, no discrete heads) and running
the projected quasi-static descent: **`max|PF|` = 0.1606 → 0.1398 over 30 steps — CONVERGED, < 0.21.** So the static
turgor-pressurised cortex+membrane+nucleus is a genuine converged equilibrium; seeding 4,420 discrete isometric
heads is the *only* thing that injects the ~1.5 residual, and no solver removes it. **Ensemble-validated (the
non-negotiable gate):** three independent cortex-construction seeds all close — seed 0 `0.1398`, seed 1 `0.1416`,
seed 2 `0.1586` (mean ≈0.147, max 0.1586, ≥24% margin under 0.21) — so the closure is seed-robust, not one lucky
realization.

**Reframe (physically grounded, not a workaround).** A physiological resting cortex's tension is **maintained
dynamically** by myosin turnover — it is a *dynamic steady-state*, not a static equilibrium. Forcing that tension
into a STATIC equilibrium (a fixed isometric seed) is the category error the 1.5 residual reports. This matches the
project's own framing: CLAUDE.md's physiological-baseline rule says *measure on the turgor-pressurised cell, then
**add myosin as a modulator***; and the PI-ratified (2026-07-23) state-conditioned Event Runtime has cortical
tension **emerge from myosin binding/unbinding EVENTS**, not a static prestress.

**Decision (Lead, 2026-07-24 — SURFACED to PI for ratification).**
1. **GATE A (static turgor baseline) is CLOSED** on the clean turgor-pressurised cortex+membrane(subdiv8)+nucleus,
   `max|PF|`=0.1398<0.21. This is the correct static baseline: a real equilibrium the runtime starts from.
   — **[status stamp 2026-07-25: this item is a PROPOSAL, status `PROVISIONAL — PI ratification pending`; see
   §11. Do not quote this line as a closure.]**
2. The resting-myosin cortical **tension is deferred to GATE B (the dynamic runtime)**, where it emerges from
   explicit NMII binding EVENTS (Bell kinetics + Hill FV on Stam-Hocky heads) — the fine-grained, physically correct
   home for a dynamically-maintained tension.
3. The fiber-quotient coarse **stays in-tree, flagged, default-OFF** — a correct, CPU-verified preconditioner for
   genuinely inter-fiber-mode-limited solves; it is simply the wrong tool for a discrete-point-load residual.
4. §§2–9 stand as the diagnostic history; §10 is the resolution.

## 11. Status stamp (2026-07-25) — §10 item 1 is `PROVISIONAL`, not `CLOSED`

§10 item 1's "GATE A … is CLOSED" is the **Lead's proposal**, surfaced above for PI ratification and **not yet
ratified**. Between 2026-07-24 and 2026-07-25 it was nevertheless recorded as settled fact in the SoT
(`cell_engine/ROLLING_ROADMAP.md:5`, "✅ GATE A CLOSED") while being absent from the PI decision queue. Corrected
2026-07-25 (PI-approved status correction):

- **Status: `PROVISIONAL — PI ratification pending`.** Queued as `PI_GAP_EVIDENCE_CARDS_2026-07-25.md`
  §ADDENDUM **Card G-A**. SoT headline corrected to PROVISIONAL.
- **What the PI is being asked to sign is a gate-contract change, not a result.** Between the FAIL and the PASS
  the **observable** changed (raw `max|F|` → projected `max|PF|`; §9, `d75e6631`, prompted by an externally
  relayed consultation) **and** the **configuration** changed (resting myosin removed from the static baseline).
  The myosin-seeded configuration never got below ≈**1.36** ≈ **6.5×** the gate.
- **The threshold itself is underived.** `max|PF| < 0.21 pN` has no derivation on disk (first appearance a bare
  `0.21` in `SESSION_HANDOFF_2026-07-22.md:51`; ~12 later hand-copied literals; absent from runtime code, where
  the executing criterion is the `sqrt(eps)·ℓ` triple test). Underived as of this stamp; a derivation is being
  attempted in parallel and is **not** performed here.
- **Physics unchanged.** §§2–10 stand verbatim, including the discrete-myosin-local-non-equilibrium reframe and
  the ~7×-heads unphysicality argument. No number in this document was altered by this stamp.
