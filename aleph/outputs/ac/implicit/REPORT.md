# Active Cell analytic implicit mechanics audit

Status: **CUDA implementation and safety guards PASS; analytic PCG, full-native convergence, and production
wall-time acceleration FAIL.** The nonlinear tournament sharply lowers the early residual but reaches a
force-gated stationary point rather than convergence. A rigid-fiber contact Schwarz extension lowers that
stationary force from `38.8567` to `18.2513 pN`, but remains `86.63x` above the invariant force gate and is
slower, so it is diagnostic rather than production. The default inner solver remains `explicit`. The optional `analytic_implicit` path may install a trial
displacement only when its PCG solve reaches the predeclared relative `sqrt(eps64)` gate and one exact nonlinear
projected-residual trial beats the matched explicit candidate. No full-native trial satisfied both conditions.

> **ERM interpretation correction (2026-07-21):** these sweep artifacts predate the Bell on/off integration.
> Their 422,094-tether cross-cell proxy is diagnostic only, and their `preload_*capacity*` fields multiply a
> continuum membrane tube-extraction force by linker count. Those fields are not a single-ERM molecular
> capacity contract and cannot satisfy NG-3.

## Scope

The Warp-CUDA matrix-free operator solves

`(a I + P K_analytic P) dx = P F`,

where `P` is the exact NF2007 inextensibility projector. `K_analytic` contains the exact actin bending Hessian
and current symmetric-positive tangents for crosslinks, WCA radial contact, explicit NMII backbone/arms/bound
crossbridges, NMII backbone and head-arm angle harmonics, LINC, and unilateral ERM. It is a numerical
accelerator only and does not replace any mechanistic force.

The topology-aware preconditioner adds two actin-fiber levels to the Cartesian node diagonal:

1. an exact per-fiber Cholesky solve for `D_external + D2^T alpha D2`; and
2. a per-fiber translation coarse space, with its initial correction advanced by a fixed `1/2`-damped Jacobi
   schedule on the represented contact graph.

The final audit uses an isotropic trace majorizer for the vector-valued coarse graph. That bound is deliberately
conservative: it prevents an unproved componentwise diagonal from being treated as a safe graph spectral bound.

The positive regularizer `a` is derived from already-declared stiffnesses omitted from the analytic operator
(membrane/nuclear surface mechanics, nuclear volume, and live pressure), with a `sqrt(eps64) * k_max` rigid-mode
floor. A fixed device-resident line search evaluates scales `1, 1/2, 1/4, ...` against the exact nonlinear
projected residual. The matched explicit-CFL state is always the fallback. No fitted trust radius, relaxation
factor, iteration-as-time interpretation, or loosened gate was introduced.

## Native gate configuration

- GPU: NVIDIA RTX A5000 (`gbook`), Warp CUDA, float64 mechanics.
- Cortex: 70,686 unique active/allocated fibers, 494,802 actin nodes.
- Whole mechanical state: 511,114 nodes, including 15,028 NMII particles, 642 nuclear nodes, and 642 membrane
  nodes.
- Explicit molecular state: 422,094 ERM tethers and 8,840 NMII heads.
- WCA production stiffness ledger: `k_max = 2,828,099.2886817832 pN/um`.
- ERM density: 600 /um2 cross-cell diagnostic proxy, **not an MCF7 production value**.
- Final majorized audit Warp mempool high-water after solver allocation: `671,676,007 bytes`.
- Exact whole-device process peak: **OPEN**. NVML lifetime accounting returned `ACCOUNTING_DISABLED` on the
  A5000, while the Warp mempool excludes non-mempool CUDA allocations; neither is an exact peak-byte ledger.

## Original safety sweep

All rows use the full native population and 20 nonlinear inner iterations. The convergence tolerance is
`sqrt(eps64) * segment_length = 7.4491633985e-9 um`.

| Configuration | Final max projected force [pN] | Last displacement / tolerance | Verdict |
|---|---:|---:|---|
| WCA off, explicit | 21.2935671 | 362.755 | reference |
| WCA off, unpreconditioned CG8 | 146.100563 | 60,254.7 | rejected |
| WCA off, PCG8 + CFL trust | 43.1940099 | 709.589 | rejected |
| WCA off, PCG32 + CFL trust | 35.1030367 | 583.126 | rejected |
| WCA on, explicit | 2,204.600875 | 9,455.427 | production reference |
| WCA on, scalar PCG8 + local residual guard | 2,215.302370 | 9,457.101 | rejected |
| WCA on, vector PCG8 + local residual guard | 2,216.277567 | 9,457.251 | rejected |
| WCA on, PCG8 + linear-convergence gate | **2,204.600875** | **9,455.427** | safe fallback; 0 implicit accepts |

The final convergence-gated run is bit-identical to its explicit reference in projected residual, last maximum
displacement, and NF2007 constraint residual (`1.33226762955e-15 um`). PCG did not reach the derived relative
`sqrt(eps64)` residual target even with a separate 128-iteration native probe, so no implicit candidate was
authorized. This is the intended failure mode: extra compute, but no mutation of the production trajectory.

## Full-native topology-aware audit

Each row is a one-inner-iteration diagnostic from the same 70,686-fiber WCA-on state. It is deliberately not a
production trajectory. The linear authorization target is `sqrt(eps64) = 1.490116119e-8` in relative residual
norm.

| Preconditioner / budget | Relative PCG residual | Last implicit / explicit exact residual | Implicit accepts | Verdict |
|---|---:|---:|---:|---|
| vector Jacobi, CG128 | 0.952879 | 39.85 | 0 | linear FAIL |
| exact fiber block, CG128 | 0.952518 | 39.82 | 0 | linear FAIL |
| two-level, coarse 8, CG128 | 0.0902172 | 39.72 | 0 | linear FAIL |
| two-level, coarse 32, CG256 | 0.00136620 | 9,981,033 | 0 | linear + nonlinear FAIL |
| additive two-level, coarse 32, CG256 | 0.000766206 | 2,054.85 | 0 | linear + nonlinear FAIL |
| two-level, coarse 32, CG768 | 0.000107395 | 196.29 | 0 | linear + nonlinear FAIL |
| **trace-majorized additive two-level, coarse 32, CG256** | **0.458052** | **0.998915** | **0** | **linear FAIL; fallback** |

The non-majorized topology variants show that the missing long-wavelength fiber modes are a real conditioning
problem, but they do not establish an admissible solver: even the best linear ratio remains about 7,200 times
above the fixed gate, and its nonlinear candidate is worse than explicit. In the conservative final variant,
the smallest tested implicit scale reduced the exact residual from `2515.8572 pN` to `2513.1275 pN`, but the
linear solve remained unconverged (`0.458052`), so the driver correctly rejected it. **There were zero implicit
accepts in every full-native row.**

This is not yet an *active* full-native convergence result. The physiological preload itself remains
unconverged, so the transactional scheduler cannot commit irreversible NMII KMC head updates. Claiming active
acceleration from these one-iteration probes would invert the engine's physical-time contract.

## Residual-monotone nonlinear accelerator tournament

Three source-independent numerical candidates now share one exact nonlinear acceptance mechanism:

- `block_descent`: the projected live fiber-block preconditioner direction `P M^-1 P F`;
- `anderson`: depth-one Type-II Anderson acceleration of that block map, with an exact singular-history
  fallback; and
- `rkc1`: a first-order Runge--Kutta--Chebyshev recurrence composed from the validated explicit step.

`tournament` generates all three from the same state and tests every declared geometric line-search scale. The
unchanged current state and the matched explicit-CFL update are also candidates. Consequently, an accelerated
run can neither increase the exact projected nonlinear residual nor lose the explicit fallback through a
non-finite sibling candidate. Candidate validity latches remain independent and GPU resident.

The RKC stage count is a numerical compute budget, not a biological parameter. A power-of-two A5000 sweep
(`s=2,4,8,16`) selected `s=8` for this diagnostic because it produced the lowest 20-iteration projected force;
it was not made the production default.

| Full-native WCA-on path, 20 iterations | Max projected force [pN] | Last displacement / tolerance | A5000 s/iteration |
|---|---:|---:|---:|
| explicit | 1,991.6517 | 9,455.43 | 0.02271 |
| block descent | 1,075.3322 | 36,627.29 | 0.27968 |
| Anderson depth 1 | 1,072.5532 | 44,661.11 | 0.28221 |
| RKC `s=2` | 2,423.3120 | 19.74 | 0.12798 |
| RKC `s=4` | 1,150.1761 | 507.97 | 0.16312 |
| **RKC `s=8`** | **421.8935** | 2,184.07 | 0.23237 |
| RKC `s=16` | 720.2443 | 8,222.94 | 0.37115 |
| **tournament, RKC `s=8`** | **45.3542** | 19,206.89 | **0.82235** |

These are the re-run residual-monotone numbers. Earlier standalone RKC diagnostics that compared only against
the next explicit candidate understated the final residual (for example, `s=8` reported `168.1090 pN`). They
are superseded by this table. The small last displacement of `s=2` is a stationary hold, not convergence; its
projected force remains above the 20-step explicit reference and far above the invariant force gate.

The 200-iteration residual-monotone run accepted 26 accelerated candidates: three block, one Anderson, and 22
RKC. It accepted zero explicit candidates, then retained the unchanged state for the remaining 174 iterations.
Projected force fell from `45.3542 pN` at iteration 20 to `38.8567 pN` at iteration 40 and remained exactly on
that plateau through iteration 200. NF2007 constraint residual remained `1.39e-15 um`. Wall time was
`163.586 s` (`0.81793 s/iteration`), about 36 times the timed explicit per-iteration cost. This is strong early
nonlinear progress, but neither a converged solve nor a production wall-time win.

### Accelerator-safe convergence invariant

Variable backtracking invalidates displacement-only convergence: a solver can make its update arbitrarily
small while retaining a large force. The convergence latch therefore now requires all three conditions:

`max|dx| <= tol`, `max constraint error <= tol`, and `dt_mu * max|P F| <= tol`.

The third condition is exactly the old explicit-step displacement criterion expressed as force, so it adds no
empirical threshold and does not loosen a gate. At the stationary tournament state, displacement was only
`5.62e-7` times tolerance, but projected force was `38.8567 pN` against the derived `0.210670 pN` gate
(`184.44x` too high). The run correctly remained `converged=false`; without this invariant, backtracking would
have produced a false success.

## Contact-cluster Schwarz audit

The plateau was measured before adding a solver. At the exact 200-iteration candidate, the live WCA graph had
31,360 undirected contacts over 58,904 nodes. Only 74 contacts were force-capped. The maximum-residual node had
one uncapped WCA neighbor (`r/sigma = 0.759996`, radial curvature `4.82726e5 pN/um`) and one crosslink to a
*different* node. Its fiber carried a `45.2713 pN` translation residual. This rejects the hypothesis that the
plateau is primarily a force-cap artifact and explains why scalar/node-local corrections miss the coupled
fiber mode.

Three nested, parameter-free additive Schwarz spaces were tested. Every pair block is FP64, device resident,
degree-weighted symmetrically, projected through the NF2007 tangent operator, and accepted only when it wins the
same exact nonlinear residual competition:

1. exact 6-DoF two-node WCA blocks;
2. exact 6-DoF two-fiber translation blocks; and
3. exact 12-DoF two-fiber rigid blocks (translation plus rotation about the live fiber center).

The local Cartesian/Galerkin matrices contain the live WCA and crosslink diagonal tangents plus the already
derived omitted-family regularizer. There is no damping factor or cluster radius. A `sqrt(eps64)` rotational
pivot handles the unobservable axial spin of a perfectly straight fiber and prolongs to zero nodal motion in
that limit.

| Full-native 200-iteration path | Final max projected force [pN] | Force/gate | A5000 s/iteration | Verdict |
|---|---:|---:|---:|---|
| original tournament | 38.8567 | 184.44x | 0.81793 | plateau |
| + node-pair contact block | 38.8567 | 184.44x | 1.03945 | 0 contact wins; reject |
| + fiber-translation contact block | 38.5330 | 182.91x | 1.03938 | small improvement; plateau |
| **+ rigid-fiber contact block** | **18.2513** | **86.63x** | **1.32612** | **best diagnostic; plateau** |
| rigid block + separate crosslink-edge overlap | 48.2600 | 229.08x | 1.32729 | trajectory worsened; reverted |

The retained rigid candidate won ten iterations in the 200-step tournament. It reduced the residual through
iteration 60, then the unchanged state won the remaining stationary tail. Constraint residual stayed at
roundoff and the final convergence latch remained false. Thus the contact-cluster hypothesis is supported—the
best residual is 2.13x lower than the prior plateau—but full-native convergence and wall-time optimization are
still failures. The next block must couple the WCA neighbor and the distinct crosslink neighbor in one
multi-fiber star, or add the fiber's internal deformation modes; merely adding crosslink pairs as independent
overlaps was experimentally falsified.

## Validation

- Pure NumPy/static wiring: exact bending stencil, symmetry/PSD, central-spring axial/transverse tangent, dense
  projected solve, public option forwarding, and removal of the explicit-displacement cap from the runtime path.
- CUDA: matrix-free bending solve vs dense oracle, WCA radial curvature vs closed form, full `P K P` solve vs
  independent dense NF2007 projector, rejection of a lower-residual candidate from an unconverged PCG solve,
  NMII angle tangent vs dense Gauss--Newton/straight-rod oracles, exact short-fiber Cholesky block, and
  fiber-translation coarse correction.
- A5000 focused nonlinear/contact/wiring suite: 23/23 PASS; full executable-scope `tests/ac`: 444/444 PASS.
- Full-native implicit profiler control (`70,686` fibers, `494,802` actin nodes, one diagnostic inner
  iteration): `0` hot-loop DtoH copies and `1` intentional post-loop DtoH positive-control copy. The run did
  not converge and rolled back, as expected from the one-iteration budget; it is residency evidence, not a
  mechanics pass.
- The default `make_inner_solve` selection remains `explicit`; no production default changed.

## Figures

- `figs/fig_ac_implicit_sweep.png` — full-native residual and displacement-to-tolerance comparison, with the
  matched explicit reference overlaid from zero-based axes.
- `figs/fig_ac_implicit_twolevel_audit.png` — relative PCG convergence and exact nonlinear residual ratios for
  every topology-aware variant; all failed authorization states and zero accept counts are visible.
- `figs/fig_ac_nonlinear_accelerator_tournament.png` — 20-iteration progress/cost comparison, 200-iteration
  force-gated plateau, and candidate/line-scale winner distribution; logarithmic force axis is marked.
- `figs/fig_ac_contact_cluster_schwarz.png` — node, fiber-translation, and rigid-fiber Schwarz trajectories,
  cost/final-residual tradeoff, measured contact-graph facts, and rigid-candidate winner distribution.

## Open numerical work

1. Extend the proved pair Schwarz space into a multi-fiber star that simultaneously contains each WCA edge and
   its distinct crosslink neighbors, or add the exact constrained internal deformation modes of both fibers.
2. Add analytic membrane/nuclear surface and pressure tangents so the regularizer does not stand in for those
   omitted numerical couplings.
3. Converge the physiological preload first, then execute a transactional multi-step run with committed NMII
   KMC state and repeat the native comparison.
4. Use the retained `18.2513 pN` rigid-cluster plateau as the next target; the independently-overlapped
   crosslink variant is falsified, and backtracking must not be extended to manufacture a smaller displacement.
5. Do not raise the PCG budget, weaken `sqrt(eps64)`, or lower the biological population merely to authorize a
   step. Re-run this exact full-native audit before changing the production default.
