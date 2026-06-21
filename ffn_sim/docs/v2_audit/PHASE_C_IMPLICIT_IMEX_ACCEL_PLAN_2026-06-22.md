# Phase C — Implicit / IMEX integrator: acceleration strategy (PLAN)

**Status: DESIGN (not implemented). Authored 2026-06-22.** Drives the "reach biological
timescales" question. NOT a parameter change — a new integrator path alongside the explicit
BAOAB. No physics is altered; only HOW the same force law is integrated in time.

## 0. The problem this solves (why explicit dt can't be "swept" up)

The Warp engine integrates overdamped dynamics `γ ẋ = F(x)` with an explicit BAOAB-limit step.
The stable timestep is capped by the **stiffest** force (a CFL limit, `dt < ~γ/k_stiffest`), NOT
by the dynamics we care about:

| force term | stiffness role | timescale it imposes |
|---|---|---|
| node-face / node-node contact repulsion | stiffest (excluded volume) | the CFL bottleneck |
| cadherin trans-dimer spring (k_trans) | stiff | small dt |
| ECM clutch spring (k_fa) | stiff | small dt |
| turgor (k_vol), bending (k_bend), edges (k_edge) | stiff–moderate | small dt |
| **spreading / de-cohesion (the OBSERVABLE)** | **soft, slow** | **minutes (what we want)** |

Measured headroom (production cfl≈0.01–0.09 at dt=8e-6 → tunneling at ~0.3): **explicit dt can
rise only ~3×** (≈2.7e-5) before contact tunnels. Biological spreading (~minutes ⇒ ~10⁸ steps,
~35 days wall) needs **~100–1000×**, which explicit integration **cannot** give. A dt *sweep* only
locates this wall; it does not move it. The wall is moved by treating the stiff terms IMPLICITLY.

## 1. Strategy — linearly-implicit (semi-implicit) overdamped step + IMEX split

Deterministic spread runs use `kT=0`, so the integrator is overdamped first-order
`γ ẋ = F(x)` — no thermostat to preserve, which makes the implicit step clean (no
implicit-Langevin subtlety; that is a later concern for `kT>0`).

**IMEX split** — implicit for the stiff, explicit for the slow/host-managed:
- **Implicit (stiff, every step):** contact repulsion + cadherin/clutch/turgor/edge/bending
  springs + nucleus. These set the CFL.
- **Explicit (soft / low-cadence, unchanged):** substrate wetting/clutch *traction*,
  lamellipodium tether, and ALL host updaters (bond form/break, remesh, division, necrosis) —
  already at ≥50-step cadence, far below any stability concern.

**Linearly-implicit step** (treat the stiff force's linearization implicitly):
```
(γ/dt · I + K) · Δx = F_total(xₙ)        with   Δx = xₙ₊₁ − xₙ ,  K = −∂F_stiff/∂x
xₙ₊₁ = xₙ + Δx
```
- `K` (the stiff Hessian) is large (3N×3N), sparse, **changes every step** (contact neighbours,
  bonds). We never form it — **matrix-free**: a CG/BiCGStab solve needs only `K·v` products,
  and `K·v ≈ −[F_stiff(x+εv) − F_stiff(x)]/ε` = **run the existing stiff force kernels on a
  perturbed position** (or analytic JVPs where cheap). The kernels we already have BECOME the
  linear operator.
- `(γ/dt I + K)` is symmetric-positive-definite-dominant for repulsion+springs (γ/dt diagonal
  regularises it), so **conjugate gradient** converges fast; matrix-free CG is GPU-ideal (Warp
  reductions + the force kernels). Target ≤ ~50 CG iters/step.

Net: cost/step ↑ (CG iters) but dt ↑ 100–1000× ⇒ large net win **if CG stays cheap** (the
γ/dt regulariser keeps the system well-conditioned at moderate dt).

## 2. Integrator order / accuracy

- **Implicit Euler** (above) is unconditionally stable but numerically DAMPED → perfect for
  reaching the **quasi-static equilibrium** (the A/A0 steady footprint — likely the validation
  observable, per the time-scale discussion), where over-damping is harmless.
- For **trajectory** accuracy (rate-resolved), upgrade to **trapezoidal / BDF2** (2nd-order,
  less damping). Decision deferred to whether the observable is equilibrium (Euler ok) or
  transient (need BDF2).

## 3. GPU implementation sketch (Warp, alongside explicit BAOAB)

1. `stiff_force(x)` = a composed launch of the stiff kernels into a force buffer (contact,
   cohesion-repulsion, cadherin/clutch springs, turgor, edges, bending, nucleus) — reuse the
   existing kernels, gated to the stiff subset.
2. `apply_operator(v) → (γ/dt) v + K v`, with `K v` = matrix-free JVP via a perturbed
   `stiff_force`. All on device.
3. Matrix-free **CG** in Warp (axpy + dot reductions + `apply_operator`); preconditioner =
   Jacobi (diagonal γ/dt + diagonal stiffness) for speed.
4. One linear solve → `Δx` → update `pos_d`. Soft/explicit terms (wetting, lamellipodium) added
   to the RHS `F_total` explicitly. Host updaters unchanged.
5. `--integrator {baoab, implicit}` flag; `--accel-dt` for the (now larger) implicit dt.

## 4. Validation plan (MUST — before trusting any speedup)

1. **Convergence to explicit:** at the SAME small dt, implicit must reproduce the explicit
   BAOAB trajectory to tolerance (correctness of the operator + CG).
2. **Conservation invariants:** turgor volume V/V0, manifold (no inversion), interpenetration
   `pen` must hold as dt is pushed up — same gates as the explicit path (D11).
3. **dt ramp:** push dt 10×, 100×, 1000×; record where the EQUILIBRIUM A/A0 (or trajectory)
   stops matching the converged explicit reference — that is the implicit accuracy limit (NOT
   a stability limit). Report it; never push past it silently.
4. **CG cost:** iters/step vs dt — confirm the net wall-clock win is real (cost/step × steps).
5. GPU parity (gbook A5000) on every step, per the project rule.

## 5. Expected payoff + risks

- **Payoff:** if CG ~30–50 iters/step and dt ↑ ~500×, a biological-minute run drops from ~35
  days (explicit) to ~hours — i.e., the A/A0 law could be evaluated at the assay's real
  timepoint, OR (if A/A0 is quasi-static) the equilibrium reached far faster.
- **Risks:** (a) K's per-step topology change (contact/bonds) can spike CG iters at events —
  mitigate with a warm-started CG + the γ/dt regulariser; (b) implicit-Euler damping can erase
  fast transients — use BDF2 if the observable is rate-sensitive; (c) `kT>0` (thermal) needs an
  implicit-Langevin scheme (BBK-implicit / geometric Langevin) — out of scope until the
  deterministic path is validated.

## 6. Sequencing (when work resumes; PI-gated)

I1. Composed `stiff_force` + matrix-free `apply_operator` + Jacobi-preconditioned CG (Warp).
I2. Implicit-Euler step + `--integrator implicit`; validation §4.1–§4.2 at matched dt.
I3. dt ramp §4.3 → find the accuracy ceiling; report the achievable speedup.
I4. (cond.) BDF2 upgrade if the observable is transient-sensitive.
I5. (cond.) implicit-Langevin for kT>0.

> This is the modern replacement for a dt *sweep*: the sweep finds the explicit CFL wall; the
> implicit/IMEX path moves the wall (stability → accuracy-bound dt). No physics changes — only
> the time-integration of the same force law.
