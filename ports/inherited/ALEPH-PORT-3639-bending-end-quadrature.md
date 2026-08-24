# ALEPH-PORT-3639 — the bending sum leaves out the ends, and the cortex is at the worst possible n

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3639` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **INTERNAL FIX with an external anchor.** No code is read from the provider. The *finding* is anchored to the provider's own Cytosim parity record, which is an independent C++ implementation, and reproduced here from scratch. |
| Provider | `/Users/sw1/ffn_cellsim/ffn_sim/ff/cytosim_parity.py` — read for its recorded finding only, READ ONLY, untouched |
| Aleph target | `aleph/vertical/cortex_filaments.py`, `aleph/vertical/cortex_surface_coupling.py` |
| Exists because | The PI, 2026-08-06 15:44 KST: *"일단 cytosim 정도는 되야지 우리 코르텍스가"*. This is the first item in that gap and the largest. |

---

## 0. The finding

**The cortex's bending energy is 49.9% of the continuum**, and the deficit is a systematic
quadrature error, not a discretisation error that shrinks with resolution in the usual way.

Measured on a circular arc of known radius, `κ = 0.07 pN·µm²`, `L = 0.3 µm`, `R = 1 µm`, against
`E = κL/2R²`:

| n nodes | interior triples | `E / E_continuum` | `(n−2)/(n−1)` |
|---:|---:|---:|---:|
| **3 — the cortex** | **1** | **0.49953138** | 0.50000 |
| 4 | 2 | 0.66638935 | 0.66667 |
| 5 | 3 | 0.74982423 | 0.75000 |
| 7 | 5 | 0.83326388 | 0.83333 |
| 9 | 7 | 0.87494873 | 0.87500 |
| 17 | 15 | 0.93748627 | 0.93750 |
| 65 | 63 | 0.98437410 | 0.98438 |

The ratio is `(n−2)/(n−1)` to within the same `sin(x)/x` term at every n. **A cortical filament has three nodes, which
is the worst value the factor can take**: the card declares `κ = 0.07 pN·µm²` and the filament bends
as though it were 0.035.

### Why, geometrically

A filament of `n` nodes has `n−1` segments of length `s = L/(n−1)`. A turning angle exists only at an
**interior** node, so there are `n−2` curvature samples, and each represents the arc from the
midpoint of one segment to the midpoint of the next — exactly `s`. The sum therefore integrates over

```
(n−2)·s   of a filament that is   (n−1)·s   long
```

and the two **end half-segments** — `s/2` each, `s` in total — carry no curvature sample at all. The
missing arc is one whole segment, independent of `n`, which is why the deficit is a ratio of counts
and not a power of `h`.

### The docstring says otherwise, and has since the law was written

`cortex_filaments.py:910-912`:

> *"Continuum limit: … `E_i -> (kappa/2) C^2 h` and the sum tends to `integral (kappa/2) C^2 ds`
> **with no coefficient left over**."*

True as `n → ∞`. **False at every finite `n`, and 2× false at the `n` this project uses.** The
sentence is about the *integrand*; the deficit is in the *domain*.

## 1. Independently anchored

`ffn_sim/ff/cytosim_parity.py:16-22` records the same factor, found by driving **Cytosim** — an
independent C++ implementation by Nédélec and Foethke — on identical inputs:

> *"on the bending-energy of a fixed circular arc, **Cytosim matches the continuum κL/2R² at every
> resolution, while FF's NF2007 interior-triple discrete energy under-counts by the end-factor
> (n−2)/(n−1)** (exact: FF·(n−1)/(n−2) recovers the continuum; converges as n→∞; ~17 % low at the
> cortex's n=7)."*

Their cortex is at `n = 7` and 17% low. **Aleph's is at `n = 3` and 50% low.** The two projects
reached the factor independently — theirs against Cytosim, this one against the closed-form
continuum — which is the strongest form this finding could take.

## 2. What is fixed, and what it is NOT

**It is a quadrature weight, not a material constant.** That distinction is the whole design.

With `E_i = (κ / 2h_i)·|t₂ − t₁|²` and `|t₂ − t₁| ≈ θ_i = C·s`, node `i` contributes
`(κ/2h_i)·C²s²`, and the continuum wants `(κ/2)·C²·w_i` where `w_i` is the arc that node represents.
So

```
h_i  =  s² / w_i
```

For an interior node `w_i = s`, giving `h_i = s` — **which is exactly what the code already does**,
and is why nothing looked wrong. The error is only at the two ends, and only because their arc is
unclaimed.

Give each end half-segment to its nearest interior triple — a zeroth-order extrapolation of the
curvature into the one place where none was sampled:

| triple | arc `w` | `h = s²/w` |
|---|---|---|
| interior | `s` | `s` |
| first, last | `3s/2` | `2s/3` |
| the only one, at `n = 3` | `2s = L` | `s/2` |

`Σ w_i = (n−2)s + s = (n−1)s = L` **exactly, at every n ≥ 3.** The weights become a partition of the
filament, which they were not before.

Measured after the change:

| n | before | **after** |
|---:|---:|---:|
| **3** | 0.49953138 | **0.99906276** |
| 5 | 0.74982423 | 0.99976564 |
| 9 | 0.87494873 | 0.99994141 |
| 129 | 0.99218727 | 0.99999977 |

The residual 0.094% at `n = 3` is the genuine `2 sin(θ/2)` vs `θ` term and converges as `O(θ²)`. It is
a discretisation error. What is removed is the systematic one.

**Nothing is done to `κ`.** No material number changes; the card still declares 0.07 pN·µm² and now
the filament bends as though it were 0.07.

## 3. This is a choice, and it is recorded as one

The provider surfaced the same finding to its PI rather than changing it silently, and gave the
reason: *"Both are valid discretisations of the same operator (Cytosim end-corrects; FF is faithful
to the paper's published interior sum)."*

That is right, and it applies here. The two positions are:

* **faithful to the published interior sum** — what the tree does today, correct as written, and 2×
  low against the continuum at `n = 3`;
* **continuum-exact at every resolution** — what Cytosim does, and what this entry implements.

**The PI directed the second** (2026-08-06 15:44 KST, *"cytosim 정도는 되야지"*). The direction is
cited, not transcribed as an authorisation, and **the first behaviour is kept as the default** so
that every committed number stays reproducible and the choice stays visible in the tree rather than
only in this file. No `decided_by` field appears anywhere.

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **B-1** | **Continuum exactness.** `E / (κL/2R²)` on a circular arc, corrected, at n = 3, 5, 9, 17, 65. | `≥ 0.998` at every n |
| **B-2** | **The weights partition the filament.** `Σ w_i` against `L`. | equal to 1e-12 relative |
| **B-3** | **The uncorrected path is unchanged.** Default arguments reproduce `K_A = 1081.819088` at ρ = 0.5, level 2. | bitwise |
| **B-4** | **The gradient still is one.** `F = −∂E/∂x` under the corrected weights, against a central difference. | `≤ 1e-6` relative |
| **B-5** | **Bending stays blind to uniform stretch.** A uniformly scaled straight filament has exactly zero bending energy and force. | exactly `0.0` |
| **B-6** | **The cortex's areal modulus moves, and by how much is reported not predicted.** `K_A` at ρ = 8, corrected vs not. | measured, no threshold |

**B-4 and B-6 are the ones that can fail honestly.** The weights enter the force through the same
`κ/h` coefficient, so if the gradient stops matching a finite difference the change is not a
reweighting but a different law, and it would be withdrawn. B-6 has **no threshold on purpose**:
this lane does not know which way `K_A` moves — bending is not the only term, and predicting the
answer before measuring it is the error `docs/results/2026-08-06-*` already records four times.

## 5. What is NOT claimed

- **Not that Cytosim was consulted here.** No Cytosim binary was run by this lane. The anchor is the
  provider's recorded parity finding plus this lane's own closed-form comparison; both agree, and
  neither is a run of the C++ code.
- **Not that `n = 3` becomes adequate.** It becomes *unbiased*. Three nodes still give one hinge per
  filament and cannot represent a shape, which is `ALEPH-PORT-3638` §C's problem and is untouched.
- **Not that the zeroth-order extrapolation is unique.** A higher-order end treatment would also
  work; this one is chosen because it is exact for uniform curvature at every `n` and is one line.
- **Not applied outside the cortex.** `intermediate_filament`, `microtubule` and `sf_arc` carry their
  own bending stencils and were **not** measured. If they share the form, they share the defect —
  that is a question for whichever lane owns them, and it is raised here rather than fixed silently.

---

## Amendment, 2026-08-06, after the gates ran: B-6, and a retraction

### The gates

| # | outcome |
|---|---|
| B-1 continuum exactness ≥ 0.998 | **PASS** — 0.99906 at n = 3, 0.99999908 at n = 65 |
| B-2 the weights partition the filament | **PASS** — `Σw = 0.299982422184` against `L = 0.299982422184`; uncorrected `Σℓ = 0.262485` |
| B-3 the retained control untouched | **PASS** — `K_A = 1081.819088`, bitwise |
| B-4 the force is still `−∂E/∂x` | **PASS** — 5.734e-10 relative against a central difference |
| B-5 blind to uniform stretch | **PASS** — energy exactly `0.0`, max force exactly `0.0` |
| B-6 what `K_A` does | **measured, and it is not what this lane said** |

### B-6, and the retraction it forces

| ρ = 8 | `K_A` affine | `K_A` relaxed |
|---|---:|---:|
| `INTERIOR_ONLY` | 14,544.8825 | **465.3530** |
| `END_CORRECTED` | 14,544.8825 | **468.7343** |
| ratio | **1.000000** | **1.00727** |

**The bending energy is 2× low and the areal modulus barely notices: 0% affine, +0.73% relaxed.**

This lane told the PI (2026-08-06 15:47 KST) that the correction *"오늘까지의 모든 피질 계수를 2배
바꿉니다"* — that it would double every cortex number. **That is retracted; it is wrong by two orders
of magnitude**, and B-6 was written with no threshold precisely because the answer was not known.

Two reasons, and both are about `n = 3`:

1. **Affine: exactly zero, by construction.** A cortical filament is built straight — three collinear
   points at `midpoint ± half·tangent` — and an affine dilation maps collinear points to collinear
   points. The bending energy of the strained configuration is therefore exactly `0.0`, and twice
   zero is zero. The affine estimator **cannot see this defect at all.**
2. **Relaxed: 0.73%.** Under tangential relaxation the filaments do bend, but a three-node filament
   has **one** hinge, so there is very little shape for the correction to act on, and the axial and
   crosslink terms dominate the virial.

**So the defect is real and currently almost invisible**, and it will stop being invisible exactly
when `ALEPH-PORT-3638` §C gives a filament enough nodes to have a shape. Fixing it first is still
right — it is a bias, and a bias that grows into relevance is worse than one that is already
visible — but it is not the lever this lane claimed it was.

### What this changes about the order of work

The crosslink term dominates the relaxed virial, so **the crosslinker's identity (`§A2`) and the
node count (`§A3`) are the levers; the bending quadrature was a correctness fix, not a lever.**
Stated here rather than left implicit, because this lane has now twice guessed which change would
move a number and been wrong both times.


---

## Amendment 2, 2026-08-06 — every number in this file was from a script the tree does not run

A verifier checking this entry found the tables **twice as far from the continuum as the shipped
code actually is**, and it was right.

**The cause.** The tables were produced by a scratch script that built a circular arc and passed
`h = L/(n−1)` — the **arc** length — as the Voronoi length. `CortexFilamentNetwork` takes its
`triple_voronoi_um` from `step`, which is `np.linalg.norm` of the difference between adjacent node
positions: the **chord**. Chord and arc differ by exactly `sin(x)/x` with `x = θ/2`, and the energy
carries that factor once, so the two conventions differ by exactly one power of it:

| n | code (chord) | this file, before | `sin(x)/x` | `[sin(x)/x]²` |
|---:|---:|---:|---:|---:|
| 3 | **0.99906276** | 0.99812641 | 0.99906276 | 0.99812641 |
| 129 | **0.99999977** | 0.99999954 | 0.99999977 | 0.99999954 |

The tables above are corrected to the values `scripts/_ws_bending_gate.py` prints, which are the
values the tree produces. `aleph/vertical/cortex_filaments.py`'s `BendingQuadrature` docstring
carried the same wrong pair and is corrected with them.

**What this says about the practice, which is the part worth keeping.** The gate script and the
ledger were written from two different programs, and only one of them is the engine. The gate
passed — B-1 asked for `≥ 0.998` and the code gives 0.99906 — so **nothing failed**, and the error
would have propagated into every artefact citing 3639 without anything going red. The rule that
would have caught it is the one this project already has and this entry did not follow: **a number
in a ledger must come from the script the ledger names.** From here the tables are pasted from the
gate's own stdout rather than retyped from a scratch run.

Found by an independent verifier, not by this lane.
