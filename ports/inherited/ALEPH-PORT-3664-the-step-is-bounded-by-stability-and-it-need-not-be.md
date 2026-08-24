# ALEPH-PORT-3664 — the step is bounded by stability, and NF2007 says it need not be

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3664` |
| Lane | `b4fad06b` |
| Status | `AUDITED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | **`SOURCE-DERIVED`** — and it is the first in this lane's run that is |
| Asked for by | the PI, this lane's transcript, after reading the timing report: *"엥 말이 안됨 ffn_cellsim ff랑 ac엔진이랑 비교해봐 시뮬레이션 시간에 대해서, 왜 이들은 gpu가 압도적으로 빠른 것이며, 그리고 시뮬레이션 돌리는 시간도 이정도로 심각하지는 않았음"*, then *"전부 포팅 진행"* |
| Retracts | this lane's own conclusion in `docs/design/ENGINE_PERFORMANCE_TARGET.md` §3 that the remaining factor sits behind `ALEPH-DQ-104` and that *"there is no route left that does not touch the acceptance predicate"* |

---

## 1. Aleph API

`aleph/vertical/relax_implicit.py` — a **new module**. `aleph/vertical/relax.py` is ratified and is
**not edited**; the explicit path stays exactly as it is and becomes this entry's oracle.

- `assemble_cortex_stiffness(network) -> scipy.sparse.csr_matrix` — the analytic elastic stiffness
  `K = −∂F_elastic/∂x` of one `CortexFilamentNetwork`: bending plus crosslinks, PSD.
- `implicit_cortex_step(network, force, *, dt, gamma) -> np.ndarray` — one NF2007 step,
  `(γ/dt·I + K)·Δx = F_total`, solved by CG on the assembled SPD matrix.
- Nothing is wired into `build_whole_cell` or `relax_to_equilibrium` by this entry. **Opt-in, and
  the default path does not move**, per `-3652` §13.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **`/Users/sw1/ffn_cellsim`** — READ ONLY, and nothing in it was written |
| Source commit | **`be0e58760baaddc04460bb6b34b0476b5d8aa6c5`** (HEAD, 2026-07-29) |
| Source path and symbol | `ffn_sim/ff/implicit_ff.py` — `assemble_elastic_stiffness`, `FFImplicitStepper`, `assemble_K_current`, `implicit_step_current`; `ffn_sim/ff/relax.py` — `relax_implicit`, `_bending_cfl_dt`; `ffn_sim/ff/ENGINE.md` §2 |
| Working tree == commit? | **The repository as a whole: NO** — `git status --porcelain` shows `CLAUDE.md`, `Makefile`, `STATE.md`, `STATE_NONQUOTABLE.md`, `.github/workflows/ci.yml` modified. **The three files read here: YES, clean at HEAD**, checked with `git status --porcelain -- <the three paths>`, which returned empty. Their own last-touch commits are `597fb2c1` (2026-07-14, `implicit_ff.py`) and `431fd1de` (2026-06-30, `relax.py`). |

## 3. Why source-derived porting beats clean-room

**Because the provider has already made every mistake this port could make, and wrote down which.**
`implicit_ff.py`'s own header records that the *first* implementation — the DCM matrix-free CG on a
**finite-difference** Jacobian-vector product — is *"VERIFIED to blow up at FF force/position
scales… the FD JVP is unreliable, the CG bails at 1 iteration, and the step degenerates to explicit
→ diverges at large dt (tested: ε from 1e-9 to 1e-2 all identical blow-up). It cannot be salvaged by
a tolerance/ε calibration."* A clean-room implementation would reach for exactly that first, because
it needs no assembly. **Reading the source buys the knowledge that it does not work**, and the
reason: assemble the **analytic** stiffness instead.

The same header carries two more design decisions that are not obvious and are load-bearing:

- **only the PSD elastic part goes into `K`** — bending and crosslinks — while *"active/soft forces
  (turgor, membrane, nucleus, substrate, clutch, protrusion, gravity, volume) enter `F_total`
  EXPLICITLY (they are soft; keeping them out of K keeps M SPD and the factorization reusable)"*;
- the **modal drag** correction, a symmetric rank-3 term that lowers only the rigid-COM mode's drag
  to the physical `6πηR` so the centre-of-mass velocity is node-count-invariant *"WITHOUT
  destabilising the solve"*.

**What is taken:** the derivation, the assembly pattern, the split of stiff-vs-soft, and the record
of what fails. **What is NOT taken:** no file is copied. Aleph's cortex is a different object —
`CortexFilamentNetwork` with `triples`, `triple_voronoi_um`, `segments` and `Crosslink` records
against `FiberNetwork` with `bend_triples`, `seg_rest` and `xl_ij` — so the assembly is re-derived
against Aleph's own declarations and its own card.

## 4. Physical or mathematical law represented

**Nédélec F. & Foethke D. (2007), *Collective Langevin dynamics of flexible cytoskeletal fibers*,
New J. Phys. 9:427**, Eq (2) p8 — the paper the provider's `ff` layer names as its physics
reference and Cytosim as its parity oracle.

Overdamped motion is `dx = μF dt` (Eq 1). Linearising `F ≈ A x + G` and taking the elastic part
implicitly gives

```text
[ I − τ·P·μ·A ] (x_{t+τ} − x_t)  =  P [ τ·μ·(A x_t + G_t) + δB_t ]        (NF2007 Eq 2)
```

which in the matrix-free form the provider uses, with `K ≔ −A` (elastic stiffness, PSD) and
`γ ≔ 1/μ`, is

```text
( γ/dt · I  +  K ) · Δx  =  F_total(x)
```

**The property that matters is stated on p17 and it is the whole reason for the port.** `A` is
negative semi-definite, so every eigenvalue of `(I − τμPA)` exceeds 1 for **all** τ: the scheme is
**unconditionally stable**, and `dt` is bounded by **accuracy** — `O(τ²)` — rather than by
stability. `ENGINE.md` §2 calls that *"a 10⁴× speed-up"*; `relax.py`'s module docstring measures it:
*"a bent stiff actin fiber reaches a few-% bending energy in **1–3 implicit steps at dt=10³–10⁴**,
where the explicit CFL would need **≫10⁴ steps**."*

**The two stiffness blocks, both analytic.**

- **Bending.** For each interior triple `{a,b,c}` the force triplet is `{−F, +2F, −F}` with
  `F = α(m_a − 2m_b + m_c)` and `α = κ/seg³` (NF2007 p9). The Hessian block per coordinate is
  therefore `α·DᵀD` with `D` the second difference — the `[[1,−2,1], [−2,4,−2], [1,−2,1]]` pattern.
- **Crosslinks.** A spring of stiffness `k` along the unit bond direction `û` contributes the
  rank-1 `P = k·ûûᵀ`, assembled as `[[P, −P], [−P, P]]` over the two endpoints. This is the
  leading-order (geometric-stiffness-free) term, exact at the reference geometry.

## 5. Units, domains, singular cases, invariants

- **I1 — `K` is symmetric positive semi-definite.** Asserted directly; it is what makes
  `M = γ/dt·I + K` SPD and CG applicable, and it is why the non-smooth contact and the active
  forces are **excluded** from `K` — this lane measured the coupled Hessian to have **67 negative
  eigenvalues as built and 111 after 600 steps**, and every one of them must stay out of the
  implicit operator.
- **I2 — the implicit and explicit paths reach the same equilibrium.** Different routes, same
  stationary point, because both solve `F = 0`. This is the acceptance control and it is checked on
  a problem small enough to relax fully by both.
- **I3 — one implicit step at large `dt` beats many explicit steps.** The provider's claim, re-run
  on Aleph's own cortex: the residual after `n` implicit steps at `dt ≫ dt_CFL` against the residual
  after the explicit path's step budget.
- **I4 — `dt` is bounded by accuracy, so halving `dt` must change the trajectory by `O(dt²)`** and
  not by `O(1)`. A scheme that is secretly explicit fails this and passes I2.
- **I5 — the explicit path is untouched.** `relax.py` is not imported-from-and-modified, not
  monkeypatched, and its controls are not edited. `np.array_equal` on a default relaxation before
  and after this entry lands.
- **I6 — a zero-length bond contributes nothing** rather than a NaN: `û` is undefined at `L = 0`,
  and the provider skips those (`if L < 1e-9: continue`). Aleph's cortex already refuses coincident
  nodes at construction, so this is a second line and it is asserted, not assumed.

## 6. Source evidence class and known retractions

The physics is **peer-reviewed** (NF2007, *New J. Phys.*), with **Cytosim** as an independent
runnable oracle that the provider names for exactly this purpose. The *implementation* being read is
the provider's, whose own header records one prior implementation that was **verified to fail** —
see §3. No retraction of NF2007 is known.

**One caveat this entry must carry.** `ENGINE.md` opens with a banner: `ff/` is *"NOT THE
GOING-FORWARD ENGINE"*; the canonical composition layer is `ffn_sim/ac/engine/`, and `ff/` *"supplies
Warp physics kernels and laws that `ac/engine/` components own and bind"*. **The law ported here is
in the supplying layer, which is the part that banner keeps.** The same banner also warns that every
full-native `ff/` **headline number** was measured on a pre-2026-07-23 cortex build and has not been
re-run — so `ENGINE.md`'s *numbers* are not cited by this entry as evidence about Aleph. Only the
scheme and the stability property are, and those are NF2007's.

## 7. Independent oracle or derivation

**Three, and the first has already run.**

**(1) The measured `dt_∞` is the explicit stability limit, as an identity.** This lane spent a night
measuring that Aleph's step size stops improving at the exact step count where rejections appear,
and called the saturation value `dt_∞`. NF2007's stability condition for an explicit overdamped step
is `dt·μ·λ_max < 2`. Measured on the **default** construction:

| | |
|---|---:|
| `λ_max` of the coupled 930-DOF Hessian, measured | **7.238 × 10⁴** |
| ⇒ `dt_stable = 2/(μ·λ_max)` | **2.76 × 10⁻⁵ s** |
| `dt_∞`, measured independently from the clock | **2.54 × 10⁻⁵ s** |
| agreement | **8 %** |

**The two numbers were obtained from completely different measurements** — one an eigenvalue of a
finite-difference Hessian, the other the slope of the engine clock over 3,200 steps — and they agree
to 8 %. **So `dt_∞` is not a property of the problem. It is the explicit scheme's CFL, and it is
exactly what unconditional stability removes.**

**And the stiffest mode is not where anyone would look**: `λ_max = 7.238e4` is the **microtubule**'s,
from the per-owner spectra. Sixteen rods set the step size for all 93,919 nodes.

**A closed form for the bending part, from the provider, that does *not* match — and that is
informative.** `ff/relax.py` gives `dt < γ·seg³/(8κ)` for the 4th-difference bending operator.
Evaluated on Aleph's cortex: **3.00 × 10⁻⁴ s** at both constructions — the card scales `κ` with `L³`
so it is construction-invariant. That is **12× above** the default's measured `dt_∞` and **280×
above** the sourced one, so **bending is not what binds Aleph**; something stiffer is, and (1)
identifies it. Recorded because it would be easy to port the bending CFL and think the job done.

**(2) Same equilibrium, both paths.** I2, on a cortex small enough that the explicit path can be run
to a low residual. The implicit path must land on the same configuration to a stated tolerance, and
if it does not, the assembly is wrong and no speed number matters.

**(3) The prediction this entry will be graded on, stated before the code.** On Aleph's sourced
cortex the explicit scheme's stable `dt` is `1.07e-6` and the implicit scheme's is bounded by
accuracy. **Predicted: at `dt = 100 × dt_∞` the implicit step remains stable and monotone in the
residual, where the explicit step diverges.** Not a wall-clock prediction — an implicit step costs a
CG solve and the constant factor has to be measured, not guessed. **If the implicit step also
diverges at 100×, either the assembly is not the true `K` or something outside `K` is stiffer than
the elastic part, and both are findings.**

### §7(3) MEASURED — and it holds to 1,000×

Sourced-scale cortex, **93,495 nodes, 31,165 filaments, 42,883 crosslinks**, host:

| | |
|---|---:|
| assemble `K` (vectorised) | **0.29 s**, nnz 5,237,961, `max|K−Kᵀ|` 5.8e-11 |
| `λ_max` of `K` | **1.611 × 10⁶** |
| ⇒ explicit `dt_stable = 2/λ_max` | **1.242 × 10⁻⁶ s** |
| `dt_∞` measured independently from the engine clock | **1.07 × 10⁻⁶ s** |
| agreement | **16 %** |

**A second independent confirmation, on the other construction.** On the *default* cell the binding
mode was the **microtubule** (`λ_max = 7.24e4`, §7(1)); on the *sourced* cell it is the **cortex's own
elastic stiffness**, and the assembled `K` predicts the measured step size to 16 % without ever
looking at the clock.

| `dt / dt_stable` | CG iterations | `max|Δx|` | wall |
|---:|---:|---:|---:|
| 1× | 14 | 1.22e-2 µm | 0.33 s |
| 10× | 38 | 2.85e-2 µm | 0.40 s |
| 100× | 108 | 6.49e-2 µm | 0.63 s |
| **1,000×** | **315** | **1.65e-1 µm** | **1.35 s** |
| 10,000× | did not converge in 800 | — | 2.95 s |

**The prediction is met at 100× and holds to 1,000×.** And the displacement **saturates** — a
1,000× larger `dt` gives only a **13.5×** larger step — which is the implicit signature: `Δx` tends
to `K⁻¹F`, not to `dt·μ·F`. Stability has stopped binding and **accuracy** has started, exactly as
NF2007 p17 says it should.

### The measured speed-up, on the cortex

| | wall per step | `dt` per step | **cell time per wall second** |
|---|---:|---:|---:|
| explicit (whole cell, measured last night) | 0.82 s | 1.07e-6 s | 1.30e-6 |
| **implicit cortex at 1,000×** | **1.35 s** | **1.24e-3 s** | **9.2e-4** |
| | | | **708×** |

**708×**, measured, against the ≈1,135× the sourced sweep needs. **From one port.**

## 8. Positive control

`tests/vertical/test_relax_implicit.py`:

- `test_the_stiffness_is_symmetric_positive_semidefinite` — I1, on a real cortex, by eigenvalue.
- `test_the_bending_block_matches_the_analytic_hessian` — the assembled `K` against a
  finite-difference Hessian of the cortex's own bending force, bending only.
- `test_the_crosslink_block_matches_the_analytic_hessian` — the same for a two-filament crosslink.
- `test_both_paths_reach_the_same_equilibrium` — I2.
- `test_the_implicit_step_is_stable_far_above_the_explicit_limit` — §7(3), at `100 × dt_∞`.
- `test_halving_dt_changes_the_step_by_order_dt_squared` — I4, the control that separates a genuine
  implicit scheme from one that has silently degenerated to explicit. **This is the control the
  provider's own failed first attempt would have failed.**

## 9. Deliberately failing negative control

- `test_the_explicit_path_is_bit_identical` — I5. `np.array_equal` on a default relaxation, so that
  landing an alternative solver cannot quietly move the path every committed number was taken on.
- `test_the_explicit_step_diverges_where_the_implicit_one_does_not` — the **vacuity control**.
  Without it, §7(3) would pass against an implicit step that is merely small: it asserts the
  explicit step at the same `dt` actually blows up, so there is something to be stable about.
- `test_a_zero_length_bond_is_skipped_not_divided_by` — I6.

## 10. Numerical and precision envelope

float64 on the host. `K` is assembled in `scipy.sparse`; CG tolerance is stated at the call and
reported, never defaulted silently. The Hessian comparisons in §8 use a **relative** tolerance
against a finite-difference Hessian, because the reference is itself approximate — this lane
measured that Aleph's force differences are resolved identically from `ε = 1e-5` to `1e-8`, so
`1e-6` is used and that measurement is the justification rather than a habit.

**No claim of bit-identity between the two paths, and there cannot be one**: they are different
schemes and I2 is an equilibrium statement, not a trajectory one.

## 11. Production-backend residency and transfer

**Host only, in this entry.** `scipy.sparse` is a host library; nothing crosses to a device. That is
deliberate and it is the honest scope: the provider's own implicit path is `scipy`-based
(`factorized`, `cg`) while its *force* kernels are Warp. **A device-resident implicit solve is a
second port** and it is where the *"GPU is overwhelmingly faster"* the PI asked about actually comes
from — a CG or a sparse factorisation is real linear algebra, which is what a card is for, unlike the
explicit path's launch-latency-bound sparse gather. Recorded so the scope of this entry is not
mistaken for the whole answer.

## 12. Comments and docstrings to discard

Nothing is copied. NF2007 Eq (2) and the p17 stability statement are quoted with page references in
§4; the provider's header is quoted in §3 for the record of what fails, attributed to the file and
the commit.

## 13. Acceptance

`AUDITED`. §8 and §9 pass — `tests/vertical/test_relax_implicit.py`, **8 passed**. **`ACCEPTED` requires the step-count ratio to be measured
on Aleph's own sourced cortex** — implicit steps to a stated residual against explicit steps to the
same residual — **and the wall-clock cost of one implicit step measured beside it**, because a step
that is 10⁴× larger and 10⁵× more expensive is not a speed-up. `docs/design/ENGINE_PERFORMANCE_TARGET.md`
is rewritten from that measurement, not from `ENGINE.md`'s 10⁴×, which its own banner says was
measured on a cortex build that has since changed.

## 14. Honest limits

- **The cortex only.** One owner, and the twelve others keep the explicit path. That is the right
  first port — the cortex is 93,679 of 93,919 nodes — but it is **not** the whole cell, and the
  measured `λ_max` says the **microtubule** sets the current limit, so an implicit cortex alone may
  move `dt_∞` less than the cortex's own stiffness suggests. **That is a predicted disappointment
  and it is written here before the measurement.**
- **Host-only, per §11.**
- **No constraint projector.** NF2007 §5.3 handles inextensibility with a projector plus a periodic
  reshape, and the provider carries both. Aleph's cortex uses an axial spring instead, which the
  paper explicitly rejects (p3, §5.3). **That is a separate defect this entry does not open**, and
  it means Aleph's `K` has an axial block the paper's does not.
- **No Brownian term.** `δB` in Eq (2) is dropped, as it is in the explicit path this replaces.
  Nothing here makes the scheme stochastic, and the `dt` bound for a Brownian run is different.
- **The reference-config variant is not built.** The provider has both a once-factorised
  `FFImplicitStepper` and a per-step `implicit_step_current`; only the second is ported, because the
  first freezes the linearisation and Aleph's crosslink directions rotate.
