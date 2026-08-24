# `ac/engine/observe` — the first three measurements of the 2026-07-25 execution plan

Build `db841185` · A5000 `cuda:0` · 2026-07-28 · record `sf_operator_observables.json` + `.npz`
Evidence: **`CUDA_UNIT` × `QuantitativeClaim.BLOCKED`** (two axes, D1-B) · gate **PASS** (D7-classified)

> ## ⚠️ CORRECTION, 2026-07-28 (PI pushback, same day this landed)
>
> **904 nodes is 0.18% of the 511,114-node native cell, and four of this report's headline numbers do not
> survive that.** Retired to `STATE.md` (c) 15 and marked `retract` in `results_manifest.yaml`:
> the **1.871× Gershgorin factor** (a function of the connectivity graph), the **condition number 2.03e11 /
> 11.31 decades** (`λ_min` is a GLOBAL mode, so a small system *understates* it — the implicit-solver
> argument survives, its number does not), the **floppy fraction 439/2,712** (its denominator is this
> slice's own composition), and the **entire Fisher block** (identifiability is a function of what is in
> the system; the α-actinin ranking may be an artifact of this slice carrying only 4 dorsal-arc joints).
> This repo already has the precedent: §(b)'s fiber-quotient coarse operator passed synthetically and
> failed at native.
>
> **What survives is law-level**, and is marked as such below: the crossbridge is not a gradient field
> (a property of the force law's differential form — one bound head settles it), and the tangent is
> multilinear + symmetric (per-contribution algebra). The instrument is sound; only the denominator was
> wrong. Sections 1 and 2 are kept verbatim as the record of what was measured, **not** as claims.
>
> Rule this produced, now standing: **develop on a slice, conclude at native.** Column probing refuses
> above 8,192 DOF by design, so going there needs Lanczos first.
>
> **Follow-up landed the same day** (`49d2215a`): the floppiness is no longer merely retired but
> **restated in the form that transfers** — see **§1b**, a per-minifilament count on a 102-DOF block that
> was verified decoupled first. And matrix-free Lanczos now exists
> (`observe/spectrum.py::lanczos_extremal_spectrum`), so the remaining three retired numbers are
> re-measurable at native rather than merely withdrawn.

The plan's §2 asks for three measurements — the λ spectrum, the Fisher information, and an energy
ledger — on one small slice, with no new physics and no new parameter. All three ran. Nothing here
computes a force: every number comes from applying the tangent and the force that the GATE-B
`nmii_sf_motor` run already launches.

**Slice.** 360 `sf_arc` nodes + 544 `nmii` particles = 904 nodes / **2,712 DOF**; 320 heads on 16
minifilaments over 16 straight motor stations; 20 pinned FA anchors. Two configurations are probed:
**t0** (no head bound, so the field is purely passive) and **bound** (318/320 heads bound after 10
accepted steps). The bound configuration is *not* an equilibrium — a tangent is perfectly well defined
there, but nothing about the equilibrium may be inferred from its spectrum.

---

## 1. The λ spectrum — the operator had never been looked at

The tangent is applied matrix-free by every inner solve in this repo, so its spectral properties were
unmeasured. Column probing assembles it exactly (a tangent is linear in its argument), and `eigh` does
the rest.

| Quantity | t0 (passive) | bound (motor engaged) |
|---|---|---|
| true `λ_max` [pN/µm] | 1.857e6 | **1.863e6** |
| Gershgorin row-sum bound of the SAME operator [pN/µm] | 3.469e6 | 3.485e6 |
| **bound / truth** | 1.868× | **1.871×** |
| stable explicit mobility step `2/λ_max` [µm/pN] | 1.077e-6 | 1.074e-6 |
| smallest resolvable `λ` [pN/µm] | 1.000e-2 | 9.16e-6 |
| condition number | 1.86e8 | **2.03e11** |
| time-scale separation [decades] | 8.27 | **11.31** |
| modes below the resolvability floor | 596 | 505 |
| rigid-body lower bound (6 × unanchored components) | 192 | 66 |
| **excess (internal floppy) modes** | **404** | **439** |
| negative (saddle) modes | 0 | 0 |

Four things follow, none of which needed a threshold:

1. **Every explicit step in this repo is sized ~1.87× smaller than it needs to be.** The Gershgorin
   bound and the truth are computed from the same assembled operator, so this measures the bound's
   conservatism and not two different formulas.
2. **11.3 decades of time-scale separation** at the engaged configuration. That is the case for an
   implicit solve stated as a number rather than as experience. Caveat, recorded: the smallest
   resolvable eigenvalue sits within a small multiple of the eigensolve's own float64 floor, so the
   condition number is a **lower** bound on the true one.
3. **~16% of the slice's degrees of freedom are floppy beyond its rigid-body content** (439 of 2,712).
   The measured null space is compared against `6 ×` the number of connected components carrying no
   Dirichlet anchor, computed independently from the bond graph; it never falls below that bound
   (which would be a hard defect), and the excess is the built model's own internal softness. Most of
   it is minifilament-internal: a head held to its backbone by one spring is free to swing about it.
   **The bipolar minifilament as built is not a rigid body.**
4. The softest passive mode is `λ = 1.000e-2 pN/µm`, occupying **5.3 of 904 nodes** (0.59%) with an
   energy-weighted gyration radius of **2.60 µm** — a mode that is sparse in node count yet extended in
   space. The full per-mode table is in the record; the full spectrum is in the `.npz`, not summarised.

The tangent is also **exactly multilinear in its stiffness parameters**: assembling each family
separately and summing gives `‖K − Σ_j K_j‖ / ‖K‖ = 9.35e-17`. That identity is checked, not assumed,
and it is what licenses the analytic gradient below.

### 1b. The floppiness, restated as a per-object property (added 2026-07-28, build `49d2215a`)

Row 3 above was retired because its denominator is the slice's own composition. The transferable form of
the same finding is a count **per minifilament**, and a minifilament is 102 DOF — so it is settled by an
eigenproblem, not by scale. Measured on all 16 blocks of the same t0 tangent:

| Quantity | Measured |
|---|---|
| DOF per minifilament (34 particles: 14 backbone + 20 heads) | 102 |
| coupling to everything outside the block, at t0 | **exactly 0.000e+00 pN/µm** |
| zero modes per minifilament | 25 |
| rigid-body content (3 translations + 3 rotations, built explicitly) | 6 |
| **internal zero-stiffness directions per minifilament** | **19** |
| — same number by projecting out the rigid-body basis, not by subtracting | **19** ✓ |
| — same on every one of the 16 blocks | ✓ |
| fraction of that subspace lying on the HEAD particles | **95.6 %** |

> **A built bipolar minifilament does not resist 19 internal directions, and 95.6 % of that subspace
> lives on its heads.** The object is the same object at any population, and the block was verified
> decoupled before the count was taken — which is the check the retired fraction never had.

Three honesty notes, all of which the record carries:

* The count is taken **two independent ways** — counting eigenvalues under the derived resolvability
  floor, and projecting the measured null space onto the complement of the six explicitly constructed
  rigid-body motions. They agree, so "25 − 6" is a measurement rather than an assumption.
* The null space is **degenerate**, so it has no canonical basis: the per-mode participation figure
  (8.46 particles) describes the eigensolver's arbitrary basis, **not the object**. What is
  basis-independent — the dimension, and the subspace's weight on the heads — is what is quoted.
* 19 internal directions against 20 heads is *consistent with* one free rotation per head about its own
  arm, and the head weight supports it — but **this run does not test that**, because a degenerate
  subspace cannot be decomposed into per-head modes without a further measurement.

## 2. Fisher information — which parameters this observable could ever infer

Observable: the **2,207 resolvable eigenvalues** (a field, per the plan's §3 rule), differentiated in
log-log form so every entry is dimensionless. Step: the derived `ε^(1/3)` central-difference optimum.

- **Effective dimension 5.42 of 8** stiffness parameters; all 8 directions are formally resolvable;
  **sloppiness 2.72 decades** — comparatively *stiff* by the usual standard, where sloppy models show
  five decades and up. This is the project's first quantitative identifiability statement.
- Response ranking `‖∂ log λ / ∂ log θ‖`: `sf_bending` 22.1 · `nmii_head_arm` 18.6 · `nmii_backbone`
  17.0 · `sf_axial` 16.5 · `nmii_angle_arm` 15.1 · `crossbridge` 14.5 · **`alpha_actinin_arc` 2.40** ·
  `nmii_angle_backbone` 0.99.
- **The α-actinin dorsal↔arc crosslink is the least identifiable of the mechanically significant
  parameters** — 6–9× below the rest. It is also the *only* SOURCED stiffness here (Ferrer 2008,
  4.6e5 pN/µm), the one that dominates the Gershgorin bound and the reason this slice needed an
  implicit solve at all. It sets the numerics and barely touches the observable.
- **The analytic gradient is validated.** Multilinearity makes first-order perturbation theory exact,
  `∂ log λ_i / ∂ log θ_j = θ_j u_iᵀ K_j u_i / λ_i`, and the finite difference agrees with it: median
  |difference| 2.0e-9, **88.3% of entries within 1e-6**, p99 1.2e-3. The worst entry (6.5e-2) sits on
  near-degenerate eigenvalues, which are not individually differentiable — expected, and reported as a
  distribution rather than as one worst case. This is the plan's "FD reference that validates any
  adjoint gradient", available now.

## 3. Energy — a conservativity test that needs no energy function

No module in this engine returns a potential, so the balance the plan wants cannot yet be formed. But
`∮F·dx = 0` for **any** conservative field whatever its potential, so the loop integral tests directly
whether an energy exists for the forces the slice launches. A nonzero residual has exactly two possible
origins and they scale differently, so the residual identifies itself — **no tolerance is applied**:

| origin | amplitude `a` | segments `N` |
|---|---|---|
| midpoint-rule truncation | `a³` | `N⁻²` |
| genuine circulation (curl ≠ 0) | `a²` | independent of `N` |

| configuration | `∮F·dx` [pN·µm] | amplitude exponent | segment exponent | verdict |
|---|---|---|---|---|
| t0 — passive | 4.187e-12 | 3.976 | 1.998 | **quadrature** |
| **bound — motor engaged** | **3.425e-08** | **1.955** | **0.00013** | **circulation** |
| control — same configuration, heads detached | 4.192e-12 | 3.974 | 2.001 | **quadrature** |

**The control is what makes this an attribution.** Between t0 and bound *two* things changed — the
configuration moved and heads bound — so neither series alone attributes anything. Detaching every head
at the **same** configuration changes exactly one of the two (the binding SoA is snapshotted, zeroed,
walked, and restored bit-for-bit, verified not assumed). The circulation vanishes, dropping by a factor
of **8,170** back onto the passive curve. So:

> The bound `nmii_sf_motor` force field is **non-conservative, and the non-conservative content is the
> crossbridge**. Measured: 3.425e-08 pN·µm per circuit of 1e-3 µm edge, scaling as `a²` and exactly
> independent of the quadrature.

Two consequences that were previously assumptions:

- **An energy ledger for this lane must carry an active-work term.** No potential exists for the
  composed field, so "monotone energy decrease" can never be this lane's convergence criterion — the
  plan already says the gate is the balance *closing*, and this is the direct measurement of the term
  that makes it not close otherwise.
- **The implicit solver's tangent is a quasi-Newton operator by construction.** It is symmetric to
  round-off at both configurations (`‖K − Kᵀ‖` max entry 1.02e-10 / 5.82e-11 against a derived floor of
  3.04e-7 / 3.05e-7), i.e. it *cannot* represent the non-conservative content it is differentiating —
  it omits the crossbridge's dependence on the segment direction `ŵ`. That is legitimate for locating a
  fixed point (`F = 0 ⇒ dx = 0` for any non-singular `A`) but it **bounds the convergence rate**, and it
  is a candidate explanation for the recorded Newton overshoot on this family of lanes.

## 4. The gate, and one honest correction

Gate: **does the PASSIVE assembled field admit an energy?** — the t0 configuration contains no active
element, so a circulation there would be an unambiguous defect. Two independent tests must agree: the
tangent's antisymmetry (exact) and the loop's own scaling verdict. Verdict **PASS**, residual/signal
`7.31e-17` against a void ceiling of 1e-2 declared before the run.

The predicate **first** declared for this observer applied that same expectation to the *bound*
configuration and evaluated **false** there. It was mis-specified: a myosin crossbridge is an active
element and its force is not a gradient, so demanding a vanishing circulation of it asks the physics to
be something it is not. The predicate is therefore scoped to the passive configuration, the bound
circulation is reported as a **measurement with its control**, and the original verdict is recorded in
the artifact under `gate_as_first_declared`. **Narrowing a gate's scope after a run is a gate-contract
change. PI-SIGNED 2026-07-28** and recorded as **GC-0001** in
`ffn_sim/docs/v2_audit/gate_contracts/contract_changes.yaml` — the ledger's first authorising row. It was not
applied retroactively. The row also names the defect it ran into: `observe` has no `GATE_*.yaml`, so there is no
contract hash for the change to be pinned to and the row authorises nothing a machine can check. Authoring that
contract file is the follow-up.

## 5. What may and may not be quoted

*Rewritten by the 2026-07-28 correction above. The original text listed the Gershgorin factor and the
effective dimension as quotable "STRUCTURE"; they are not — being structural is not the same as being
population-independent, and that conflation is exactly what the correction caught.*

**Quotable** — only what does not depend on how big the system is:

- the **crossbridge is not a gradient field**, with its attribution control and the scaling exponents that
  identified it (`a²` / N-independent vs `a³` / `N⁻²`);
- its two consequences: an energy ledger for this lane needs an active-work term, and the symmetric tangent
  is a quasi-Newton operator;
- the tangent is **multilinear** and **symmetric**, for the families this slice launches;
- the FD-vs-analytic gradient agreement;
- **19 internal zero-stiffness directions per minifilament, 95.6 % of that subspace on the heads** (§1b) —
  a single-object property, taken only after verifying the block is exactly decoupled, and cross-checked
  by rigid-body projection against the eigenvalue count.

**Not quotable** — the 1.871× Gershgorin factor, the condition number / decades, the floppy-mode fraction,
and the whole Fisher block (all four population-dependent, see the banner); any eigenvalue as a
stress-fiber material property (every mechanical parameter here is an unresolved KB/PI GAP, cards N1–N9 and
S1, except the α-actinin crosslink); and any statement about the equilibrium, which this run does not
measure.

**The transferable form of the floppiness result** was measured on 2026-07-28 (build `49d2215a`) and is
§1b above: **19 internal zero-stiffness directions per minifilament**, not a fraction of this slice's DOF.
The per-mode participation figure inside it is **not** quotable — the subspace is degenerate, so that
number describes the eigensolver's arbitrary basis rather than the object.

## Figures

Regenerate everything with `python ffn_sim/scripts/ac_observe_vis.py`.

| Figure | What it shows |
|---|---|
| `figs/spectrum.png` | the full eigenvalue field at both configurations on a log axis, with the resolvability floor, the true `λ_max` and the Gershgorin bound drawn on the same axes; the near-null block and its rigid-body lower bound are marked, nothing is truncated |
| `figs/localization.png` | participation ratio and energy-weighted gyration radius per softest mode — the stress decay length, with "fully delocalised" and "single node" overlaid as references |
| `figs/sensitivity.png` | the Fisher eigenvalue spectrum against its resolvability floor (the sloppiness plot), and the per-parameter response of the observable |
| `figs/loop_work.png` | the loop residual vs amplitude and vs segment count for all three configurations, with **both** predicted power laws drawn; the passive and heads-detached series coincide on the `a³`/`N⁻²` line while the motor-engaged series sits four decades above it, flat in `N` |

## Reproduce

```
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
  ffn_sim/scripts/ac_observe_sf_operator_native.py \
  --k-axial 1000 --k-xb 1000 --k-backbone 1000 --k-head-arm 100 --backbone-lp 1.0 \
  --k-on 50 --f-stall 0.5 --v0 0.12 --kappa 0.5 --capture 0.21 \
  --n-bb 14 --n-side 10 --backbone-len 0.301 --head-offset 0.2 \
  --bind-steps 10 --build-commit db841185 \
  --out ffn_sim/outputs/ac/observe/sf_operator_observables.json
```

`--build-commit` is required on the native machine: its tree is synced rather than checked out, so a
git-only stamp degrades silently to `"unknown"` — the missing-build defect the record exists to close.
The record marks it `source: "declared"` so it is never mistaken for one git verified.
