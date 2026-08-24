# T10 · `extracellular_medium` — exterior Stokes medium

**Lane** `ecm-adhesion` · **Build** `a9861bed` (branch `codex/ff-ac-codex`) · **Date** 2026-07-28

## Status, stated before anything else

| Axis | Value | Why |
|---|---|---|
| `EvidenceRung` | **`CUDA_UNIT`** | The operator executed on the A5000 over the full native membrane surface (642 quadrature points / 1,926 DOF), build `576ebf42`. Not `CONNECTED`: `membrane_medium_traction` is not dispatched. |
| Gate verdict | **FAIL — 3 of 10 mechanism predicates** | And the FAIL stands. Diagnosis below: the three are mis-specified, the fix is a gate-contract change, and that is **PI's to sign** (`GC-0002`, `pi_signature: PENDING`), not this lane's to apply. |
| `QuantitativeClaim` | **`BLOCKED`** | `mu_medium` is a PI-GAP — no KnowledgeClaim, no SourceEvidence. No drag magnitude from this work is quotable. |
| `BLOCKS_CRAWL` | **still standing** | The medium loads the rigid modes, but the world boundary is still symmetric (no wall) and the connector is not dispatched. |

**No tier-(a) `STATE.md` row was added.** The run's gate verdict is FAIL, and a repository that adds a
headline row off a failing run has learned nothing from its own `(c)` list. The row comes after `GC-0002`
is signed and the gate is re-run — in that order, which is what the contract ledger checks.

## The native run

`ffn_gpu.py run` → gbook A5000, lease acquired after a 15 s wait, 6.62 s wall including Warp compilation,
11 device resistance solves. The runner's own guard reported *"remote tree matches local for 12 first-party
files"* — the hash check proposed in the previous revision of this report now exists and passed.

What the operator did, all of it viscosity-free or a ratio:

- **All six rigid modes carry strictly positive resistance.** Eigenvalues `0.147244 ×3` (translation) and
  `11.8502 ×3` (rotation) — three-fold degenerate in each block, which is the isotropy the geometry
  demands. This is the BLOCKS_CRAWL statement in its falsifiable form: an `aI` regulariser cannot produce
  two distinct three-fold-degenerate blocks whose ratio is set by the cell's radius.
- **Stokes' law recovered** at the surface quadrature's own discretisation error: translation
  `measured/6πμa = 1.041539` and rotation `measured/8πμa³ = 1.116898`, matching the host convergence
  study for this refinement level to the digit.
- **Rest is exactly force-free**; dissipated work strictly negative; force exactly linear in viscosity;
  net force on a *spinning* sphere `|Σf|/Σ|f| = 2.15e-14` (it must vanish, and it does); torque sign −z,
  opposing. CG converged in 55 / 79 iterations, condition number 2157.

## The FAIL, and why the threshold is not being touched

Three predicates missed: translation isotropy `2.77e-13`, translation–rotation decoupling `1.25e-12`,
time-reversal asymmetry `4.74e-13` — each against a declared floor of `64·eps64 = 1.42e-14`.

All three assert properties that are **exact for the operator**, and all three are measured **through a
CG solve**. The fourth member of that family, `operator_symmetry`, is the only one with no solve in its
path, and it came out **exactly `0.0`**. That is circumstantial, so it was tested rather than asserted.

**The discriminating measurement.** A CG-tolerance sweep at `1e-6 / 1e-8 / 1e-10 / 1e-12` (four native
runs). An operator defect is tolerance-**independent**; solver error is not:

| cg_tol | iters | isotropy | decoupling | time-reversal | transverse | drag ratio |
|---|---|---|---|---|---|---|
| 1e-6 | 22 | 1.06e-10 | 1.14e-12 | 2.18e-11 | 9.32e-14 | 1.041539447 |
| 1e-8 | 38 | 3.93e-10 | 1.53e-12 | 2.47e-10 | 1.54e-10 | 1.041539447 |
| 1e-10 | 55 | 2.77e-13 | 1.25e-12 | 4.74e-13 | 6.84e-14 | 1.041539447 |
| 1e-12 | 76 | 4.21e-13 | 1.02e-14 | 8.64e-16 | 2.09e-16 | 1.041539447 |

Spans across the sweep: ×932, ×150, ×2.9e5, ×7.3e5. **None is horizontal**, every value stays under the
`cond·tol` bound at every tolerance, and the physical answer — the drag ratio — is **bit-identical across
all six decades**. The operator-defect hypothesis is refuted; the residuals are the solver's.

**Stated against my own case, because it is the honest limit:** the dependence is *not* a clean power law.
`decoupling` sits near-flat over three decades before dropping two. That is expected — CG error depends on
which Krylov iterate the stopping test lands on, not smoothly on the tolerance — but it means the sweep
refutes tolerance-independence rather than demonstrating a clean `O(tol)` law, and it should not be
reported as the latter. A second symptom of the same mis-specification: at `1e-12` two of the three
failures clear the *unchanged* floor while `linear_in_viscosity`, which passed at `1e-10`, newly fails.
Which predicates pass is being decided by solver noise around a direct-solve-sized threshold — the
clearest possible sign that the floor, not the operator, is wrong.

**What is NOT being done.** The floors are not edited. Under the no-gate-loosening rule a wrong gate is
surfaced, not fixed inline, so `GC-0002` is filed in `docs/v2_audit/gate_contracts/contract_changes.yaml`
with `pi_signature: PENDING`, which authorises nothing by that file's own rules. It proposes replacing the
direct-solve floor with the CG error bound `cg_relative_tolerance · cond(M)` — **derived** from quantities
the run already measures, tightening automatically whenever the solve is tightened, and leaving
`operator_symmetry`'s exact-equality predicate alone. The row is written **before** any superseding run,
and no re-run to a PASS has been attempted. T10 also has no `GATE_*.yaml`, the same hole `GC-0001` exposed,
so even a signed row would authorise nothing machine-checkable until that file exists.

## What the component is

PI decision D5-A (2026-07-28) declared the 14th component and its `membrane_medium_traction` connector while
stating that the declaration bought none of the physics. The gap it named: the cell's six whole-body rigid
modes carry zero internal stiffness by construction, so every implicit path in the tree holds them with a
numerical term — `aI` with `a = 1/(dt/γ)` in `ac/engine/sf_implicit.py`, `compute_regularization_kernel` in
`ac/cell/`. That term is a mobility the CFL bound already computed. It has no relation to the cell's size or
shape, so it can neither be checked against a measurement nor **load** a mode.

`ac/engine/medium_exterior.py` puts a physical object in that slot: a regularised-Stokeslet (Cortez) single
layer over the cell's outer surface quadrature, with the resistance `f_medium = -M⁻¹ v` applied matrix-free
by device-scalar CG. Its tangent `M⁻¹/dt` is symmetric positive definite on the whole space including the
rigid block, and the sphere limits `6πμa` and `8πμa³` are exact references for it.

**Why a boundary integral and not a grid.** Outside the cell `k → ∞`, so the Brinkman screening length
`ℓ_B = √k → ∞` and Darcy is *wrong* there rather than approximate. What remains is Stokes, and the
single-layer integral representation **is** the exterior Stokes solution — the only approximations are the
surface quadrature and the blob width, both refinable. A grid solve would add a truncated far field and an
immersed-boundary blob on top, and at any affordable `dx` could not resolve the cell–substrate lubrication
film that crawl traction lives in.

## What is measured, and where

All host-side, on the dev Mac, at radius 7.5 µm with `ε = h` held fixed across refinements. The convergence
table and its provenance live in **`ffn_sim/ac/engine/medium_params_t10.yaml` §derived.blob_epsilon**, which
is the parameter card for this component; the gates that assert it are
`ffn_sim/tests/ac/engine/test_medium_exterior.py` (25 host gates green, 4 CUDA-gated and skipping).

Exact structural properties, verified and population-independent because they are properties of the
operator's *form*:

- the assembled mobility is **bit-exactly symmetric** (the kernel is even in the separation vector, so the
  two blocks are the same IEEE products) — Lorentz reciprocity as a discrete identity, not a tolerance;
- it is positive definite, so dissipation is strictly positive for every non-zero surface velocity;
- on an icosphere the translation resistance block is **exactly isotropic** and decouples from rotation, both
  at the round-off floor, because icosahedral symmetry forces any invariant rank-2 tensor to be isotropic —
  a check that needs no reference value at all;
- the response is exactly linear in the viscosity and exactly reversed under time reversal.

`ε` is **derived**, not tuned: it is the quadrature's own mean nearest-neighbour spacing, so it falls under
refinement, and the gate asserts the error **sequence** plus the fitted order rather than an absolute
tolerance — which is precisely how a refinable approximation avoids being frozen into a magic constant.

## On "native" for this component, because the word is doing real work here

The exterior medium addresses
the cell's outer *surface* only, and the native membrane surface **is** the icosphere at subdivision 3
(`ac/cell/compartments.py:104,641`). So the gate's 642 quadrature points are this component's *complete*
native population, not a slice of it — unlike the `nmii_sf_motor` slice at 0.18% of the cell. What is not
native is the **composition**: no composed cell, no connector dispatch.

## Figures

Rendered from the committed record by `scripts/ac_gate_t10_medium_vis.py` (no GPU, no re-run), committed
beside the data. This is the ISOLATION render the ladder gate requires of a component at `CUDA_UNIT`.

- `figs/traction_field.png` — per-node traction under rigid translation and rotation, each drawn in ITS OWN
  plane, with the net force and the net torque judged separately. Shows `|Σf|/Σ|f| = 2.15e-14` on the
  rotation panel, i.e. a spinning sphere is visibly not pushed.
- `figs/rigid_modes.png` — the six resistance eigenvalues against `6πμa` and `8πμa³`.
- `figs/grid_invariance.png` — Stokes-law error vs surface spacing, log-log, with the fitted order.
- `figs/dissipation.png` — per-node traction·velocity: strictly negative at every node, not merely on average.
- `figs/solver_tolerance_sweep.png` — **the figure that carries the diagnosis**: the four residuals against
  CG tolerance, with the `cond·tol` bound and the declared direct-solve floor on the same axes.

Three defects were caught by the figures *before* any GPU time was spent, on a host-reference fixture: the
rotation panel drawn in the wrong plane for a z-rotation; a torque ratio labelled onto a force arrow whose
value vanishes by construction; and an `aI` reference line drawn at a value no run had measured. The third
recurred while building the sweep figure — a "what an operator defect would look like" line drawn at a
measured value — and was removed the same way. Drawing an un-measured reference on a measurement is the
specific failure this project's visualization rule names, and it takes deliberate effort not to repeat.

## Trap #4 — the obligation changed shape, and that is the finding

`AC_EXECUTION_PLAN_2026-07-25.md` §T10 item 2 states it as a deletion: remove `γ_node = 6πηR/Nc` in the same
change. The audit (`docs/v2_audit/T10_GAMMA_NODE_AUDIT_2026-07-28.md`) measured where that term is and found
it has **never existed in `ac/`**. Both live sites are in `ff/` — `network_warp.py`'s whole-cell compression
driver and `implicit_ff.py`'s modal COM drag — plus `motility_warp.physical_node_gammas`, and every caller is
a legacy `scripts/ff_*` driver. Nothing in `ac/**` imports any of them.

So `ac/` goes from **zero** velocity-proportional drag owners to one, and the obligation is an *invariant*
rather than a deletion: `medium_exterior.assert_single_dissipation_owner` refuses those four symbols by AST
walk. AST rather than a source regex, because this repository has twice been bitten by regex static guards
firing on docstrings, and the medium module's own docstring discusses `γ_node` at length. A guard is strictly
stronger than the deletion that was asked for — a deletion is discharged once and cannot prevent a re-import,
and it would also have broken four `ff_*` drivers whose results are already invalidated-pending-rerun for
unrelated reasons.

Not discharged, and named so nobody reads this as closing them: `ff/` still carries all three terms (correct
*for* `ff/`, which has no exterior bath — `engine-library`'s call); the shear-dissipation hole
(`COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md` §2c hazard 9) is untouched, and the exterior medium does not
fill it because a surface operator loads *rigid* modes, which is a different functional support from internal
shear; and T6b's repo-wide dissipation census is still to come, of which this guard is the first instance.

## Open items, in the order they block things

1. **`GC-0002` → PI signature.** The gate stays FAIL until then, and no re-run to a PASS will be attempted
   before it — that ordering is the only thing separating a fix from a fit, and this ledger checks it by
   commit order. Paired with authoring `GATE_T10_exterior_medium.yaml`, without which even a signed row
   authorises nothing machine-checkable.
2. **`mu_medium` → PI.** Needs a SourceEvidence row for a measured culture-medium viscosity at 37 °C with an
   OK verdict. Blocks only magnitudes; the mechanism is already verified without it.
3. **`HALF_SPACE_BLAKE`.** T10 item 3, the asymmetric world boundary (basal 2D collagen + free media face).
   Exact and closed-form, with the Brenner wall correction as its oracle. Spreading and crawl stay gated on
   it — a symmetric boundary can do neither.
4. **Dispatch `membrane_medium_traction`.** Cross-lane: the claim belongs to `SurfaceBody` (`cortex`). Until
   then the rung cannot exceed `CUDA_UNIT`.
5. **Retire the `aI` regularisers** in `sf_implicit.py` and `ac/cell/` once (4) lands — `sf-motor` and
   `cortex` lanes. Doing it before the medium is dispatched would leave the rigid modes held by nothing.

## Remote-tree provenance — a real hazard, but not the one this session first wrote down

**Corrected.** This report initially blamed Syncthing for `~/ffn_ac_native` being stale. That was wrong, and
the correction is the useful part: **Syncthing on this repo shares only `ffn_sim/outputs/` — source code has
never been synced, by design.** So `ffn_sim/scripts/` not carrying a fresh commit is the *expected*
behaviour, the five minutes spent polling for it were spent waiting on something that was never going to
happen, and the files deleted in the 2026-07-28 cleanup are still on gbook for the same ordinary reason: the
tree has not been rsync'd since. The two concurrently running Syncthing daemons are a known half-finished
upgrade and are irrelevant to code, because code was never their job.

**What genuinely does not have a guard, and is worth PI attention.** `ffn_gpu.py run` stamps the *local*
`HEAD` onto whatever code the remote tree happens to hold. Since that tree is updated by a manual rsync, a
forgotten rsync silently produces a run record whose build commit is wrong — the exact defect `STATE.md` (b)
names as the audit's top finding, and it cannot be caught by reading the record afterwards. Two mitigations,
both cheap:

1. **Standing habit (no code change):** rsync the source tree to gbook as the first step of any native
   launch, then verify — `md5sum` the modules the driver imports on both sides — rather than trusting that
   the remote copy is current. That is what was done here before launching.
2. **Proposed guard (`observe-infra`'s file, not this lane's):** have `ffn_gpu.py run` hash the remote copies
   of the driver and its imported modules against the local `HEAD` blobs and refuse the launch on a
   mismatch. The stamp is only meaningful if the code it names is the code that ran.
