# Rigid-constraint Lagrange shell tension — design + Sanity Gate

> **Status**: EXPERIMENTAL. **NOT** yet PI-ratified. Touches
> `ffn_sim/integrator/constrained_baoab.py` (integrator-freeze) for a
> *diagnostic exposure*, not a behavioural change. Requires PI sign-off.
>
> **Author**: Lead session, 2026-05-28 autonomous /loop.
> **Reference commit (baseline)**: `0f65e75`.

## 1. Problem

The current KU-3.5 cortical tension measurement
(`_tension_method_of_planes` in `h3_ku35_tension.py`) sums **soft-bond
contributions only** — xlinks, myosin head-actin attach, ERM, myosin
internal bonds. The dominant tension source in the cortex network is the
**rigid actin backbone** under M-SHAKE constraint. That contribution is
currently **not** in the measurement; it is documented but missing.

Quantitatively: with 333 engaged motor heads each pulling up to ~0.5 pN
(F_stall), total soft-bond contractile force ≤ 166 pN, giving cortical
tension ≤ 2.7 × 10⁻⁶ N/m = 2.7 × 10⁻³ mN/m. Far below the KU-3.5 target
of 0.5 mN/m. The factor-200 gap is the **rigid-bond stress** that
propagates through the actin backbone but is invisible to a soft-bond
plane-method sum. Without including this, KU-3.5 will systematically
under-report tension and FAIL the gate even when the physics is correct.

## 2. Theory

For a rigid bond enforced by M-SHAKE, the Lagrange multiplier λᵢⱼ is the
scalar force-along-bond that the constraint applies each step to
maintain |r_i − r_j| = r₀. For a bond between particles i,j:

```
F_constraint_on_i = +λᵢⱼ · (r_j − r_i) / r₀
F_constraint_on_j = −λᵢⱼ · (r_j − r_i) / r₀
```

The total scalar bond tension is **T_ij = λᵢⱼ** (signed: +T for the
bond pulling its endpoints together, −T for compression).

M-SHAKE in our plugin solves the per-chain tridiagonal Newton system
inside `shake_project_chains` (constrained_baoab.py:266). Each Newton
iteration computes a Lagrange-multiplier vector λ for the chain's
M bonds. The final λ at convergence is the per-step Lagrange multiplier
needed to satisfy the constraint.

**The Lagrange multiplier λ is ALREADY computed every step — it is
just discarded after projection.** We need only **expose** it.

## 3. Algorithm — diagnostic-only modification

Modify `shake_project_chains` to also return the converged λ vector
(currently it returns only the projected positions):

```python
def shake_project_chains(pred_pos, ref_pos, chains, rest_length,
                          inv_mass, box, *, tol=1e-10, max_iter=100):
    # ... existing M-SHAKE solver ...
    return projected_pos, lambda_vec  # NEW: also return λ
```

The custom-Action's `act()` method captures λ each step into a buffer:

```python
class ConstrainedLeimkuhlerMatthewsBAOAB(hoomd.custom.Action):
    def __init__(...):
        ...
        self._lambda_buf = None   # most-recent step's λ vector
        self._lambda_history = [] # if self.record_lambda is True

    def act(self, timestep):
        ...
        new_pos, self._lambda_buf = shake_project_chains(...)
        if self.record_lambda:
            self._lambda_history.append((timestep, self._lambda_buf.copy()))
```

This **does NOT change** any positions, velocities, forces, or integrator
behaviour. It only exposes a quantity that was already computed and
discarded. The behavioural assertions in `tests/integrator/` should pass
bit-for-bit.

## 4. Tension measurement (consumer side)

New `_tension_method_of_planes_with_rigid` in `h3_ku35_tension.py`:

```python
def _tension_method_of_planes_with_rigid(sim, R_cell, n_planes=12):
    # 1. Soft-bond contribution (existing).
    gamma_soft = _tension_method_of_planes(sim, R_cell, n_planes)

    # 2. Rigid-bond contribution.
    act = <find ConstrainedLeimkuhlerMatthewsBAOAB action>
    lam = act._lambda_buf      # shape (n_filaments, n_bonds_per_filament)
    chains = act._chains       # particle-tag pairs per chain bond
    # ... convert chain-local bond index → particle tag pair → bond vector
    #     × λ → contribution to plane-cut tension γ.
    gamma_rigid = sum_over_crossing(...)

    return gamma_soft + gamma_rigid
```

## 5. Sanity Gate (must pass before PI sign-off)

| # | Test                                  | Pass criterion |
|---|---------------------------------------|----------------|
| 1 | M1 dimer/trimer empty actuation       | All previous BAOAB-bit-for-bit tests still pass. |
| 2 | Lambda dimensional analysis           | λ has units of force (N) — verify via `pN / dt` etc. |
| 3 | Lambda sign convention                | At equilibrium with random thermal noise, ⟨λ⟩ should be of order 0; pulled (stretched) bond → +λ. |
| 4 | Reproducibility                        | Same seed → same λ trajectory (record_lambda=True), within machine ε. |
| 5 | Symmetry across plane choice           | Method-of-planes γ_total = γ_soft + γ_rigid is independent of plane orientation (variance over plane orientations ≤ existing soft-bond variance). |
| 6 | Plateau magnitude vs literature         | KU-3.5 plateau γ_total ∈ [0.35, 0.65] mN/m for full assembly with k_ERM = 5.6e-5 and 100 motors. (This is the gate itself.) |
| 7 | k_record overhead                      | record_lambda=False → 0% wall overhead; record_lambda=True every step → ≤5% overhead (just a copy + append). |

## 6. Open question for PI

Even though this is a diagnostic exposure (no behaviour change), it
adds a return value to the public `shake_project_chains` signature.
That could break downstream callers if any exist outside the
constrained_baoab module.

**Search outcome**: the only caller of `shake_project_chains` is
`ConstrainedLeimkuhlerMatthewsBAOAB.act()` inside the same file. No
external imports. Safe to extend the return tuple.

**Test impact**: `tests/integrator/test_constrained_baoab.py` mocks /
calls the projector? Need to audit before merge. (Lead has not yet
audited; defer to PI sign-off cycle.)

## 7. Decision branches (for PI)

| Option | What it does | What it buys |
|--------|--------------|--------------|
| **R0** Accept soft-only γ | KU-3.5 measured as soft-bond only; documented in result file. | No code change. KU-3.5 BANDS will need adjustment OR a separate "rigid contribution" reported. PI sign-off on band adjustment required either way. |
| **R1** Diagnostic Lagrange exposure (this doc) | Modify `shake_project_chains` return + Action buffer + new measurement function. ~50 lines, 1 file in integrator/ + 1 in scripts/. | KU-3.5 fully-implemented per literature theory; γ_total reportable. Sanity Gate 1-7 above. |
| **R2** Re-derive bands at soft-bond scale | Theoretical exercise: what γ_soft would be expected if literature's reported γ is the full quantity. Probably ≪0.5 mN/m. | No code change but PI band adjustment needed. Less rigorous. |

Default Lead recommendation, **conditional on PI**: **R1**. The
information is already computed; we just need to expose it. Without R1,
KU-3.5 cannot be honestly compared to literature.
