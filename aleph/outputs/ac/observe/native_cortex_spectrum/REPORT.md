# Native cortex tangent — family certification, then a matrix-free spectrum

Record `record.json` · driver `ffn_sim/scripts/ac_observe_native_spectrum.py` · figures `figs/`
Evidence axes: **`EvidenceRung` × `QuantitativeClaim`** (D1-B) — the record self-stamps both, and the
values below are read from it rather than restated here.

> **STATUS: RUN LANDED AT FULL NATIVE.** `551,434` nodes / `1,654,302` DOF on `cuda:0`, build
> `a44a1de3` (declared — the run host is not a git checkout), gate **PASS**. The build stamp is
> trustworthy for a reason this session had to build: `ffn_gpu.py` now verifies the driver's whole
> first-party import closure against the remote tree before launch, and reported 79/79 matching.
>
> **The headline: the production explicit inner step is stable at native — and the constant it is
> sized by was never a bound.** `λ_max = 2.946× kmax`, so the largest SINGLE stiffness in the build
> under-states the operator by exactly the accumulation over incident bonds. The `0.1` safety factor in
> `dt_mu = 0.1/kmax` absorbs that and leaves **6.79×** headroom. See §5.2.
>
> **The small end did NOT converge, and that is reported rather than papered over.** One of the three
> numbers `STATE.md` (c) 14 retired is now re-measured at native; the condition number and its decades
> are **not**, because the folded pass did not resolve `λ_min` at this Krylov dimension. §5.3.
>
> **The two smoke runs during development are not committed and nothing from them may be cited.** They
> did their job — one of them refuted this driver's own first method (§2.3) — but one ran the superseded
> differencing construction and the other ran with `overlap_free_cortex` still defaulted ON.

---

## 1. What this exists to do

On 2026-07-28 four spectral numbers were retired to `STATE.md` (c) 14 one session after they landed.
The reason was not that the instrument was wrong — it was that the instrument had been pointed at 904
nodes, which is **0.18% of the 511,114-node native cell**, and all four numbers are functions of the
connectivity graph or of system size:

| Retired number | Why a slice cannot carry it |
|---|---|
| Gershgorin / `λ_max` = 1.871× | a row sum is a property of node degree, i.e. of the graph |
| condition number 2.03e11 / 11.31 decades | `λ_min` is a GLOBAL mode; a small system **understates** the separation |
| floppy fraction 439 / 2,712 DOF | its denominator is the slice's own composition |
| the Fisher block (effective dimension 5.42 / 8, …) | identifiability depends on what is in the system |

The floppiness was restated the same day as a **per-minifilament** count, which is a single-object
property and does transfer — that row is in `STATE.md` (b). The other three needed an instrument that
works where column probing refuses (above 8,192 DOF, by `assemble_dense_operator`'s own design), and
`observe/spectrum.py::lanczos_extremal_spectrum` landed to be that instrument. **This run is its first
caller.**

## 2. Method

### 2.1 Two operators, and why neither substitutes for the other

The engine's inner iteration is `pos += dt_mu · P F` (`ac/cell/inner_mechanics.py`). Linearised, and
using that `range(P)` is invariant under it — for `v ∈ range(P)`, `(I − dt·P K)v = v − dt·P K P v` —
the iteration matrix on the constraint tangent space is `I − dt_mu · P K P`. So:

* **`P K P` restricted to `range(P)`** is what the solve's stability and conditioning are properties
  of. Its raw small end is identically zero by construction (the constraint directions), so the
  measurement is taken **on the restriction**: both Lanczos passes start from `P(random)` and the
  folded operator is `σ·P − P K P`, not `σ·I − P K P`. Without that second substitution the identity
  term carries the constraint directions at eigenvalue `σ` — the folded operator's *largest* — and the
  pass converges to the projector's null space, reporting `λ_min = 0` for any operator whatsoever.
  That failure mode is covered by a test with its own negative control
  (`test_without_the_projector_the_same_operator_reports_its_constraint_zeros`).
* **`K` unprojected** is what conservativity is a statement about: symmetry of `∂F/∂x` is exactly the
  integrability condition for `F = −∇U`. `P K P` can be symmetric while `K` is not, so measuring only
  the projected operator would launder an asymmetric assembly.

### 2.2 The folding shift is measured, not chosen

`λ_min = σ − λ_max(σP − K)` must not depend on `σ`. Rather than pick one shift and defend it, the run
measures the small end at **`σ = 2θ` and `σ = 4θ`** (`θ` from a short probing pass) and reports the
disagreement. A chosen constant becomes a measured invariance. The two shifts also use independent
start-vector seeds, so the run reports `λ_max` twice and the seed-to-seed agreement says whether the
large end converged or whether it is reporting its start vector.

Every Ritz value carries its own **rigorous** Parlett residual bound `|β_m · s_{m,i}|`: some true
eigenvalue lies within that distance. So convergence announces itself and nothing needs a chosen
tolerance to be believed. What Lanczos does **not** give is the interior of the spectrum, so a mode
*count* — the retired floppy fraction — stays a dense-probe question and is not attempted here.

### 2.3 Symmetry is Lanczos's unchecked premise, so it is certified separately

A Krylov method handed an asymmetric operator returns confident numbers about the symmetric part and
issues no warning. Symmetry therefore has to be established where the operator can be assembled. The
certification stage builds a **probeable** cell, assembles `K` by column probing, and then assembles
each stiffness family **alone** by disabling every other family's live gating attribute.

Two things make this more than a slice result:

1. The statements taken from it are **per-contribution algebraic** — is family *j*'s tangent symmetric;
   is the assembly the sum of its families. Both are properties of the differential forms being
   summed, and one instance settles them. This is the same argument `STATE.md` (b) already accepts for
   the `nmii_sf_motor` multilinearity row. **No count, ratio, spectrum or magnitude from the
   certification build is quoted**, and the record says so in its own `scope` field.
2. Certification is **by isolation, not by launch dimension**. The WCA steric kernel launches over
   every node and contributes exactly nothing when no pair is inside its cutoff — which is the
   *designed* state of an `overlap_free_cortex` build. A launch-dimension census would count that
   family as covered. A family whose isolated tangent is exactly zero is reported `exercised: false`
   and named in `families_unexercised`, which is a coverage hole stated rather than hidden.

**The obvious construction was tried first and is wrong, which the run measured.** Isolating a family
as `K_on − K_off` made *every* family report the identical asymmetry — one half-ulp of the **full**
operator's largest entry, not of the family's. A difference of two matrices each accumulated at the
scale of the whole assembly carries the whole assembly's round-off, so the smallest families sat under
the crosslink family's noise while being scored against floors derived from their own much smaller
magnitudes, and most of them "failed" for that reason alone. Assembling each family alone removes the
subtraction: its floor is genuinely its own, and `exercised` becomes an exact test rather than a
thresholded one. It is also the **stronger** multilinearity test — if one family's tangent depended on
another's presence, the leave-one-in sum would not reproduce the assembly, whereas under leave-one-out
the sum reproduces it by algebra no matter what the physics does. The failure is kept as a regression
test at the measured magnitudes
(`test_a_small_family_is_not_buried_under_a_large_one_s_round_off`).

`STATE.md` (b) flags the cortex-only paths (multigrid, fiber-block, ERM, membrane) as untested for
exactly the multilinearity + symmetry properties this stage measures, so the stage also closes a
recorded gap.

Two families are expected to come back **uncertified**, and each carries its reason in the record
rather than being quietly absent. `steric_wca` needs interpenetrations, so `--certify-overlap-free`
defaults OFF. `myosin` cannot be certified on *any* probeable cortex build — its sites are generated
from cortex node proximity at native areal density, four decades above the column-probing ceiling — so
the record instead carries a bounded cross-reference: `sf_implicit._stiffness` launches the same
uniform-pair, angle-Gauss-Newton and segment-crossbridge kernels and was measured symmetric at build
`db841185`. The one thing that cross-reference does **not** cover is named too: the non-segment
`add_crossbridge_stiffness_kernel` that `ac/cell` takes when `myosin.segment_runtime is None`.

### 2.4 The instrument is scored against ground truth before it is trusted

On the certification build both a dense eigensolve and Lanczos are possible, so the run compares
`λ_max` from the two **on the same Warp operator**. That check has no value at native, which is
precisely why it is done where it has value.

## 3. The gate, declared before the run

Three conjuncts, each closing a way this measurement could be confidently wrong:

1. the assembled `K` admits a potential to within the ledger's own float64 accumulation floor
   (`assemble_balance_tolerance_ratio`, PI D8 — the same function the accepted-step balance gate uses,
   never a second independently-written tolerance);
2. **every family that actually contributed is individually symmetric.** An assembly can be symmetric
   overall while one family's asymmetry cancels another's, and that cancellation is a property of the
   configuration, not of the forms — it would not survive a population change;
3. the large end announced its own convergence through Parlett's bound rather than through an
   iteration count.

`VoidCeiling(0.01)`: an antisymmetric part above 1% of the tangent's own largest entry leaves no
separable conservative structure, so a PASS there would be a statement about round-off and a FAIL a
statement about nothing. Residual/signal = `max|K − Kᵀ|` / `max|K_ij|`.

A run given `--stage native` alone **refuses to score itself** (`form_verdict` raises) rather than
emitting an unqualified record.

## 4. Population and configuration

From `record.json:census`. The spectrum stage ran at **full native**; the certification stage is a
separate, deliberately tiny build whose only outputs are the algebraic statements of §5.1.

| | spectrum stage | certification stage |
|---|---|---|
| nodes / DOF | **551,434 / 1,654,302** | 604 / 1,812 |
| cortical F-actin fibers | **70,686** | 40 |
| actin nodes | 494,802 | — |
| inter-fiber crosslinks | **1,413,720** | — |
| ERM tethers / membrane faces | 40,962 / 81,920 | — |
| NMII head springs / crossbridges | 8,840 / 8,840 | 0 |
| device | `cuda:0` (A5000) | `cuda:0` |

Krylov dimension 120 per pass (1.59 GB host basis, full reorthogonalization), probing pass 40.
Regularization exactly `0.0`. Configuration probed is the **as-built** state, not a relaxed one — a
tangent is well defined at any configuration, but nothing about an equilibrium follows from this one.

## 5. Results

### 5.1 The cortex tangent is multilinear, and every exercised family admits a potential

This closes the gap `STATE.md` (b) flags on its multilinearity row — *"only for the families this slice
launches: the cortex-only paths (multigrid, fiber-block, ERM, membrane) are untested."*

| family | `‖K_j‖_F / ‖K‖_F` | `max\|K_j − K_jᵀ\|` | its own floor | verdict |
|---|---|---|---|---|
| `crosslink` | 1.000 | 4.657e-10 | 2.774e-06 | symmetric |
| `membrane_erm` | 1.913e-04 | 4.547e-13 | 7.971e-10 | symmetric |
| `membrane_area_edges` | 8.934e-06 | 3.553e-15 | 9.838e-12 | symmetric |
| `nucleus_linc` | 5.094e-06 | 7.105e-15 | 2.822e-11 | symmetric |
| `bending` | 1.995e-06 | **0 exactly** | 1.338e-12 | symmetric to the last bit |
| `steric_wca` | — | — | — | **NOT exercised** (launch dim 604, tangent exactly 0) |
| `myosin` | — | — | — | **NOT exercised** (launch dim 0) |

**Multilinearity residual `‖K − Σ_j K_j‖/‖K‖ = 1.0160e-16`.** Because each family is assembled ALONE
rather than differenced out, this also establishes that no family's tangent depends on another being
present — which the leave-one-out construction cannot test at all.

Assembled `K`: `max|K − Kᵀ| = 4.657e-10` against a derived floor of `2.774e-06`, i.e. **four decades
under**. The projected operator `P K P` is separately symmetric (`6.985e-10`), which also exercises the
runtime per-fiber Thomas projector.

**Instrument check**: on the same Warp operator, Lanczos and the dense eigensolve agree on `λ_max` to
**5.29e-16 relative**. The iterative path is therefore scored against ground truth before being used
where no ground truth exists.

**The two uncertified families are a real hole, stated as one.** `myosin` cannot be certified on *any*
probeable cortex build — its sites come from cortex node proximity at native areal density, four decades
above the column-probing ceiling — and the record carries a kernel-level cross-reference to `db841185`
plus the one variant that cross-reference does **not** cover. `steric_wca` was expected to be fixed by
running the certification with `--certify-overlap-free` OFF; **it was, and the family still assembled to
exactly zero**, because 40 great-circle arcs on a 7.5 µm shell put no pair inside the WCA cutoff. So the
cross-reference note's suggested remedy is, on this evidence, insufficient — closing it needs a build
dense enough to have contacts yet small enough to probe, and whether that window exists is not yet
established.

### 5.2 The production explicit step is stable, and `kmax` was never a bound

Measured on `P K P` restricted to `range(P)` — the operator the inner iteration `pos += dt_mu·P F`
actually linearises to:

| quantity | value |
|---|---|
| `λ_max(P K P \|range(P))` | **8.331183e+06 pN/µm** (Parlett bound 6.92e-18; two independent seeds agree to 2.24e-16) |
| `λ_max(K)` unprojected | 9.599162e+06 pN/µm (bound 1.07e-22) |
| build's own `kmax` | 2.828099e+06 pN/µm |
| **`λ_max / kmax`** | **2.9459** |
| `dt_mu` | 3.535944e-08 µm/pN |
| **`dt_mu · λ_max`** | **0.2946** — stable (forward-Euler limit 2) |
| headroom at the production step | **6.79×** |
| true stability limit `2/λ_max` | 2.4006e-07 µm/pN |

**What this says.** `kmax` is the largest *single* stiffness constant in the build, not a row sum, so it
bounds nothing: a node accumulates one contribution per incident bond, and `2.9459` **is** that
accumulation, measured. The `0.1` in `dt_mu = 0.1/kmax` therefore does exactly the work of covering it,
and leaves 6.79× beyond. The production step is stable — this is the first time that was checked rather
than assumed.

**FIRE, reported and deliberately not adjudicated.** FIRE raises the step to `10·dt_mu`, putting the
product at `2.9459`, above the plain-descent limit. That number is *the same number* as `λ_max/kmax` —
an identity (`10 · (0.1/kmax) · λ_max`), not a second independent fact, and it should not be read as
one. Whether it is unstable is **not** established here: FIRE is damped MD with a fictitious mass and a
velocity, and its bound scales differently. What this run supplies is the measured `λ_max` that FIRE's
own analysis would have to use.

### 5.3 The small end did not converge, and no condition number is reported

| | `λ_min` (folded) | its bound | resolved? | shift disagreement (`σ=2θ` vs `4θ`) |
|---|---|---|---|---|
| `P K P \|range(P)` | 52.11 pN/µm | 1.315e+03 | **no** | 4.84% |
| `K` | 39.59 pN/µm | 1.184e+03 | **no** | 0.39% |

The Parlett bound is **25× larger than the value** in both cases, so the reported `λ_min` is a property
of the start vector and the iteration count, not of the operator. `condition_number` is therefore
`null` in the record, with its reason, rather than a number built from an unresolved end.

This is the honest limit the method declared in advance: the solve-free folded small end converges as
`(λ₂ − λ₁)/(σ − λ₁)`, so **the more ill-conditioned `K` is, the slower its own small end converges** —
the measurement gets hardest exactly where it matters most. Reaching it needs either far more Krylov
iterations or shift-and-invert, and shift-and-invert requires a linear solve, which would make the
measurement depend on the very solver it is meant to diagnose.

### 5.4 Cost

`145.87 s` wall, 560 operator applications, `0.2605 s` per Krylov iteration; `comparable: true` because
this run's gate PASSED. Not a benchmark: the driver takes no physical-time step, and Warp kernel
compilation is included on a cold cache.

## 6. What may be quoted, and what may not

The record self-stamps **`CUDA_UNIT` × `QuantitativeClaim.BLOCKED`**.

**Quotable**, all measured at the native population stated in §4:

* the **stability verdict and its structure** — `λ_max/kmax = 2.9459`, `dt_mu·λ_max = 0.2946 < 2`,
  headroom 6.79×, and the statement that `kmax` is not an upper bound on the operator;
* **multilinearity `1.0160e-16`** and the per-family symmetry verdicts of §5.1 — per-contribution
  algebraic properties, so these transfer with population;
* the **instrument agreement** `5.29e-16` between Lanczos and the dense eigensolve;
* the **negative result** of §5.3: at Krylov dimension 120 the folded small end does not resolve.

**Not quotable:**

* **`λ_min`, the condition number, the decades** — unresolved (§5.3). `STATE.md` (c) 14's condition
  number is re-measur*able* but still not re-measur*ed*;
* **any eigenvalue as a cortical material property.** `k_xb`, `k_backbone` and `k_linc` are unresolved
  PI-GAPs and the myosin family is present in the native build, so `8.331e6 pN/µm` is a property of this
  parameter set, not of a cortex;
* **anything about equilibrium.** The tangent was taken at the as-built configuration;
* **the certification build's counts, ratios or spectrum** — it is not native, and §2.3's scope says so.

## Figures

Regenerate with `python ffn_sim/scripts/ac_observe_native_spectrum_vis.py` (no GPU, reads the record).

| Figure | What it shows, and the reference drawn on the same axes |
|---|---|
| `figs/family_certification.png` | per-family share of `‖K‖`, and each family's asymmetry against **its own** derived round-off floor |
| `figs/native_ritz.png` | Ritz values with Parlett bounds as error bars, for `P K P` and `K` |
| `figs/explicit_step.png` | `dt_mu·λ_max` at the production step and at FIRE's maximum, with the forward-Euler limit 2 drawn — and only the plain-descent bar judged against it, because FIRE is damped MD with a differently-scaled bound |
