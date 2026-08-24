# Resting-baseline diagnosis (2026-07-22c) — geometry fix landed + f_rupt/subdiv root cause

Branch `codex/ff-ac-codex`, geometry fix committed **`23022295`**. This continues handoff 2026-07-22b
(§DO THIS FIRST tether-geometry check). The check is done and it reframes the blocker: it is **not one
thing**. Three stacked constraints, two now resolved, one remaining.

## 1. Geometry: ERM tethers were tangential (FIXED — `23022295`)

The overlap-free cortex dispersed filament radii symmetrically `±t/2` about `R_CORTEX_UM=7.40`, so outer
nodes reached `7.50 = R_mem`. Nearest/radial ERM pairing then picked cortex partners at ~the membrane
radius → tethers ~zero-length / **tangential** → cannot carry the radial turgor load, and `u=d/|d|` in the
preload is degenerate.

Fix: disperse **inward** `[R-t, R]` — cortex outer surface stays at `R_CORTEX_UM`, thickness extends toward
the interior. Steric decorrelation preserved; γ-parity (`t=0`) path untouched.

Full-native (70,686) ERM tether **radial fraction** (1=radial, 0=tangential), overlap_free=True:

| pairing | before fix | after fix |
|---|---|---|
| nearest | mean 0.51, 189/642 tangential, 6 degenerate, 24637 cortex nodes touching membrane | mean **0.91**, 0 tangential, 0 degenerate, 0 touching |
| radial  | — | mean **0.994**, 0 tangential, 0 degenerate |

Uniform 0.10 µm radial gap restored. Necessary but **not sufficient** (below).

## 2. f_rupt: the coarse resting membrane has NO single-ERM equilibrium (build at subdiv ≥6)

> **⚠️ CORRECTION (PI 2026-07-22).** `f_rupt` here is a **continuum membrane tube/tether-extraction force
> scale** (bilayer bending + membrane–cortex adhesion), NOT the single-ERM rupture force. Labeling the
> `if tension > f_rupt: return` cap as "the single ERM ruptures" below is a **mis-attribution**: the actual
> single ezrin–F-actin unbinding force is ~**50 pN** (Braunger 2014), and `k_erm=4600 pN/µm` IS Braunger's
> single-molecule stiffness — so the code mixes a single-molecule stiffness with a continuum rupture threshold.
> This means the "subdiv-3 has NO equilibrium because 44 > 11.4" conclusion is **contingent on the wrong
> threshold** and must be re-examined in P0/P2: a single ERM can hold ~50 pN, and with the (production-HOLD)
> MCF7 ERM *density* each linker carries far less. The subdiv-≥6 recommendation is still a safe build choice,
> but the *reason* is the per-node/per-linker load vs the CORRECT single-molecule capacity + real density, not
> this continuum `f_rupt`. See `AC_DECISION_CARDS_2026-07-22.md` Card 2. The measurements below stand; only the
> `f_rupt` interpretation is corrected.

`f_rupt = 2π√(2κ_m(γ_mem+γ_MCA)) = ` **11.43 pN** (KB-3.B1.4 band 5–40). The ERM force kernel
(`erm_tether_force_kernel`) drops any tether whose tension exceeds `f_rupt` (super-threshold ⇒ bleb, no
candidate force). Membrane node turgor = `ΔP·A_node`:

| subdiv | verts | turgor/node | vs f_rupt (11.4) |
|---|---|---|---|
| **3 (resting default)** | 642 | **44.0 pN** | **EXCEEDS** → single ERM ruptures → load dropped |
| 4 | 2562 | 11.0 pN | ~at threshold |
| **6 (ratified dynamic)** | 40962 | **0.69 pN** | far under → holdable |
| 7 (validation) | 163842 | 0.17 pN | under |

**This is why no solver ever converged:** at subdiv 3 there is literally no force-balanced resting state —
one ERM per node cannot hold 44 pN, so the membrane is unheld regardless of solver. The recurring "≈11 pN
transmitted" in prior handoffs is exactly the `f_rupt` cap. **The resting baseline must be built at membrane
subdiv ≥6** (the already-ratified dynamic resolution), not subdiv 3.

Preload transmission probe (8k, overlap_free, radial pairing), membrane max|F| / whole residual:

| subdiv | before preload | after ERM preload | note |
|---|---|---|---|
| 3 | 47.4 / 47.4 | 47.4 / 47.4 (0% ) | ERM > f_rupt → dropped, nothing transmits |
| 4 | 12.3 / 12.3 | mean 10.3→3.6, cortex reaction 11.5 | transmits except few nodes at cap |
| 6 | 0.78 / 0.78 | membrane mean 0.64→**0.18**, cortex reaction 2.22 | **load path works**; residual moves membrane→cortex |

## 3. Remaining blocker: FIRE is CFL-throttled at the reachable subdiv-6 equilibrium

At subdiv 6 an equilibrium EXISTS (residual_start 0.78 pN, well-conditioned) but FIRE barely moves:
`inner_max_displacement ≈ 1.7e-11 µm/step`. `dt_mu` is set by the stiffest bond (crosslink k≈8e5 pN/µm), so
the **soft membrane** (k_erm 4600) needs to drift `turgor/k_erm ≈ 0.15 nm` but crawls at 1e-11 µm/step —
~1e7 iters. The residual lives on the slow mode; the step is limited by the fast mode. Classic stiff
soft/stiff coupling.

- **ERM preload** holds the membrane (rest pre-stretched) but transfers the reaction to the force-free
  cortex → residual 0.78→2.22 (now a cortex-node reaction needing ~2.8e-6 µm cortex motion).
- **Global cortex prestrain** (scale all crosslink rests by 1−ε) does **not** produce clean hoop tension in
  the irregular network — it injects unbalanced local forces: residual rises monotonically with ε
  (ε=1e-4 → +55 pN; ε=1e-3 → 554 pN). This mechanism is wrong for cortex pretension.

**FIRE-iteration escalation (subdiv 6, 8k, no preload) — the decisive result.** Best candidate residual vs
n_inner (start = 0.77656 pN, gate ≈ 0.21):

| n_inner | candidate | drop | maxdisp (µm) |
|---|---|---|---|
| 200 | 0.776554 | +4.5e-6 | 1.7e-11 |
| 2000 | 0.775972 | +5.9e-4 | 1.9e-10 |
| 50000 | **1.181** | **−0.40 (WORSE)** | 4.5e-9 |

FIRE does not slowly converge — it **diverges at long horizons**: the adaptive `dt` accelerates
(`maxdisp` grows 1.7e-11 → 4.5e-9) until it overshoots and the residual climbs from 0.78 → **1.18**.
`residual_end` stays pinned at `residual_start` only because the non-converged step is rolled back. The
residual force lives on a mode (soft membrane / ERM-reacted cortex vs the stiff-crosslink-throttled `dt_mu`)
that explicit/FIRE cannot descend. **Explicit/FIRE structurally cannot close the resting baseline**, even
with the geometry fixed and the equilibrium existing at subdiv 6. Any solution must therefore either start
essentially AT the balanced state (near-exact coupled preload — FIRE can't be trusted to descend a residual)
or use a conditioned implicit stiff step.

## The fork (for PI)

The resting baseline needs the soft membrane + its cortex reaction to reach balance without a 1e7-iter
explicit crawl. Options:

1. **Consistent coupled preload** — pre-stretch ERM (membrane held) AND give the paired cortex nodes a
   *targeted, per-node* reaction preload (NOT global prestrain) so the whole system starts at ~0 residual.
   Needs the cortex to carry Laplace hoop tension `T = ΔP·R/2` distributed correctly; the ε→T mapping is a
   derivable numerical calibration of the built network (not a magic number). MCF7 `γ_cortex` is a PI-gap
   only for the *consistency gate* `ΔP=2γ/R`, not for computing T from the sourced `dP_hyd`.
2. **Preconditioned/implicit stiff step for statics** — the soft/stiff split is exactly what an implicit
   solve (or a diagonal-mobility preconditioner rescaling soft DOFs) handles. Handoff 22b retired implicit
   for statics; this finding argues for a *conditioned* stiff step, reconciling with 22a's "coupled
   multilevel preconditioner is the dominant work."

Recommendation: pursue (1) first — a targeted coupled preload is physical (it IS the resting pre-stress) and
avoids reopening the implicit-vs-FIRE decision. If the targeted preload can't be made self-consistent
cheaply, (2) is the fallback.

## Update — option (1) tried and RULED OUT; (2) is unavoidable

Implemented the coupled **rest-length preload** (`scripts/ac_coupled_restlength_preload_probe.py`): hold the
built positions fixed and adjust rest lengths (ERM + crosslinks) so the as-built geometry is force-balanced —
the physiological pre-stress, done out-of-hot-loop with the exact GPU force kernels (accumulate → read →
adjust r0 → assign), a Jacobi relaxation `r0 += ½(f_i−f_j)·û / k` per link. Result at subdiv 6, 8k:

- **Undamped** (both signs): diverges — sign=+1 flat at the 2.22 ERM-floor for ~10 sweeps then blows up
  (2.22 → 2468 pN over 40 sweeps); sign=−1 explodes immediately (→1e17).
- **Damped** (0.1, 0.01): stable but **stuck at 2.22 pN** — the ERM preload holds the membrane, but the
  resulting cortex reaction cannot be absorbed by local rest-length adjustment. No stable regime *reduces* the
  residual.

**Conclusion: rest-length space has the SAME stiff-connected-network conditioning problem as position space.**
Adjusting `r0` to cancel the cortex reaction is solving `K δr = −f` by unpreconditioned Jacobi, which
stuck-or-diverges exactly like FIRE. Option (1) reduces to option (2). Both cheaper shortcuts are now
eliminated:

| approach | result |
|---|---|
| geometry fix alone (22b hypothesis) | necessary, not sufficient |
| explicit/FIRE position solve | diverges (0.78 → 1.18) |
| coupled rest-length preload (option 1) | stuck at 2.22 / diverges |

## Update 2 (PI "ㄱㄱ") — implicit-CG solver chain → the TRUE root cause: relaxed cortex has no radial stiffness

Discovered the coarse-mode multilevel preconditioner (22a) is **already built** in this branch
(`implicit_mechanics.py`: `rigid_strain_coarse_basis` l≤2 = 3 translations + 3 rotations + 6 strains,
Galerkin `A_c=BᵀAB` + Cholesky + restrict/prolong; `add_erm_stiffness_kernel`, crosslink/bending/LINC
tangents, membrane area-tension edge springs). gbook had been reverted to a pre-coarse-mode version
mid-session (Codex/Syncthing) — re-aligned by rsync of `implicit_mechanics{,_analytic}.py`.

Solver sweep at subdiv 6, 8k, overlap_free + radial + ERM-preload (residual after preload = **2.22 pN**,
the cortex ERM-reaction):

| solver knob | candidate | reading |
|---|---|---|
| `analytic_implicit` cm=0 / 3 / 12 | **2.12** (identical) | coarse deflation useless — residual is NOT a global l≤2 mode |
| cg_max_iterations 96 / 300 / 800 | **2.134** (identical) | the CG **linear solve already converges** — not a CG-conditioning problem |
| n_inner 80 / 300 / 1000 | 2.12 → 2.15 → **1.74** | each Newton step reduces the residual only ~4–7%; descends but ~impractically slowly |
| regularizer a = 1000 → 0.042 | **2.15** (identical) | the scalar `aI` (nucleus k_vol=1000) is NOT the over-damper |

Every lever fails identically ⇒ **the analytic Jacobian `PKP` itself does not stiffly resist the residual's
mode**, and it is not the regularizer.

### A soft-cortex hypothesis — RAISED then REFUTED by measurement (recorded so it is not re-chased)

I first hypothesised the cortex was radially *soft* at rest (force-free crosslinks have zero transverse
tangent `k·(L−r0)/L·(I−ûû)`) and that hoop pre-tension was the fix. **Measurement refutes this:**
- Cortex uniform breathing stiffness ≈ **71.5 pN/µm/node**, and it is **tension-independent** (unchanged from
  0 to 736 pN mean crosslink tension) — because breathing is *axial* (`k·ûû`), not transverse.
- The **localized** radial stiffness of the ERM-paired cortex nodes is ≈ **17,653 pN/µm/node** — very STIFF,
  also tension-independent. A 2.22 pN ERM reaction needs only **0.13 nm** of cortex drift to balance.

So the cortex is NOT soft, and hoop pre-tension does NOT stiffen the loaded mode. The 2.22 pN residual has a
stiff equilibrium 0.13 nm away that the solvers simply are not reaching.

### Confirmed operator gap (the concrete lead): the ERM tangent vanishes at the force-free rest

`add_erm_stiffness_kernel` (implicit_mechanics.py:245) returns **zero tangent when `extension = length − rest
≤ 0`**. At the resting build every ERM tether is force-free (`extension = 0`), so **the membrane nodes get NO
ERM stiffness in the analytic operator** — the solver has nothing to step the membrane against, and the
turgor residual (0.78 pN, no-preload) stays put. Evidence: analytic_implicit, subdiv 6, 500 iters —
no-preload 0.78 → 0.75 (3.3%, ERM tangent absent, extension=0); with-preload 2.22 → 2.03 (8.4%, ERM tangent
present because the preload stretches it). The unilateral spring's tangent at `extension=0` should take the
tension-side stiffness `k` when the residual force would stretch it, not 0.

This is a concrete gap — **but fixing it did NOT move the residual** (recorded so it is not re-chased):
changing the guard to `extension < 0` (so the axial ERM tangent `k·ûû` applies at the extension==0 rest, SPD
via `_central_action`) and clearing the Warp kernel cache left the no-preload result **bit-identical**
(0.7766 → 0.7507) and the with-preload result identical (2.22 → 2.03). So the missing ERM tangent is real but
is **not** what caps the descent — the change was reverted (unproven, shared validated operator).

**Net: still an OPEN diagnosis after eliminating five levers and three mechanisms.** What is airtight: the
elimination table (coarse-modes / cg-iters / regularizer / FIRE / rest-length all fail identically), the two
cortex stiffness measurements (breathing 71.5, localized 17,653 pN/µm/node — cortex is NOT soft), and that the
implicit-CG **linear solve converges** yet each **Newton step** only cuts the residual ~4–8%. That combination
means the analytic Jacobian `PKP` is not delivering the (measured, available) stiffness into the step for the
turgor/ERM-reaction mode — but the specific missing coupling is not yet identified. Next session should
instrument the CG at the Newton-step level: dump `‖P F‖` vs `‖F‖` (is the residual in the projection
nullspace?), dump the proposed `dx` on the membrane/loaded-cortex nodes vs the 0.13 nm expected, and check the
`_central_action` crosslink/ERM tangent assembly against a finite-difference of the true force. That
finite-difference-vs-operator check is the decisive next probe; do it before writing any new solver code.

## Operational note (bit PythonPath)
gbook `~/ffn_cellsim` (stale June clone) carries the `-e` editable install and shadows the `ffn_sim` package;
absolute `from aleph...` imports resolve stale. Run every gbook job with
`PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim` (parent first) so path-finding beats the editable
meta-finder. Verified fresh resolution. Likely the source of prior "partial-sync" confusion.
