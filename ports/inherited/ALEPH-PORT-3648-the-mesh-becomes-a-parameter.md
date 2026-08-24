# ALEPH-PORT-3648 — the filament stops being three nodes, and two defaults would hide what that costs

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3648` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-07` — **AFTER the code, by roughly twenty minutes. See §0.** |
| Port class | **INTERNAL.** One builder argument and the arc that generates a filament's nodes. Nothing is read from a provider. |
| Aleph target | `aleph/vertical/cortex_surface_coupling.py` |
| Depends on | `ALEPH-PORT-3643` (which blocked this), `-3638` §C (which asked for it) |
| Exists because | `Filament(..., 3, ...)` was hard-coded at two sites and the cortex could not be refined at all. |

---

## 0. A process failure, recorded rather than backdated

`CLAUDE.md` §3: *"Every port needs a ledger entry in `ports/ledger/` written **before** the code."*
**This entry was written after its code**, by about twenty minutes, during an unattended run.

It is recorded here instead of being quietly dated earlier, because the value of the rule is exactly
that a gate written afterwards is a gate written knowing the answer — and every other entry in this
lane's arc (`-3643` through `-3647`) was written before its code and two of them failed their own
gates as a result. That is the rule working. This one did not get that protection, and a reader
should discount its gates accordingly.

The one mitigation that is real: **§2's measurements were taken by an independent survey agent before
this lane wrote a line**, and they are what the design follows.

## 1. What changed

`build_radial_cortex_network` gains `filament_node_count: int = 3`. Both node-block generators —
the mesh-parameterised retained control and the density-parameterised population — now place nodes on
a shared arc `linspace(-1, 1, n)` scaled by the half-length, so a refinement means the same thing in
each. `3` is the **retained control** and the default does not move: `K_A` = 1081.819088 through
`build_vertical`, bitwise.

It also gains `steric_law: StericLaw = StericLaw.NODE_PAIR`, and the reason that is in this entry
rather than its own is that **it was not settable from anywhere in `aleph/` at all**.
`ALEPH-PORT-3643` landed `StericLaw.SEGMENT_PAIR` and the only code that ever reached it was three
`scripts/_ws_*` files assigning the attribute after construction. A law nothing in the package can
select is a law the package does not have.

## 2. The two defaults that would make a refinement study measure something else

Measured by an independent survey before the change:

**`BendingQuadrature.INTERIOR_ONLY`, the shipped default**, is `(n−2)/(n−1)` of the continuum. On a
circular arc, `E / E_continuum`:

| `n` | INTERIOR_ONLY | END_CORRECTED |
|---:|---:|---:|
| 3 | **0.49869893** | 0.99739787 |
| 5 | 0.74951181 | 0.99934909 |
| 9 | 0.87485759 | 0.99983725 |
| 17 | 0.93746185 | 0.99995931 |

**Refining under it silently doubles the effective bending rigidity**, from `κ/2` toward `κ`. A
refinement study that leaves it on moves the mesh and `κ_eff` at once and can attribute neither.

**`StericLaw.NODE_PAIR`, also the shipped default**, is not refinement invariant. On this builder's
own shell geometry:

| ρ | `n` = 3 → 9, node-pair energy | node pairs | segment-pair energy |
|---:|---|---|---|
| 24 | 2.620487e+08 → 2.985348e+08 (**+13.9%**) | 19 → 1,302 | identical to every printed digit |
| 100 | 5.973331e+18 → 7.876258e+18 (**+31.9%**) | 2,880 → 25,087 | identical to every printed digit |

That is `-3643`'s F-2 confirmed on real shell geometry rather than on a test crossing, and it makes
`-3643` §3's ordering claim a measurement instead of an argument. It also shows the block was only
half lifted: **§3 was satisfied at module level and not at the builder**, because until now the
builder could not select the law it says must be used first.

## 3. The gates

Written after the code — see §0 — so they are reported with that discount.

| # | Gate | Outcome |
|---|---|---|
| **K-1** | The retained control is untouched: `K_A` = 1081.819088 at ρ = 0.5, level 2. | **PASS**, bitwise |
| **K-2** | `n` = 3, 5, 9, 17 all build and carry the right node count. | **PASS** — 7,479 / 12,465 / 22,437 / 42,381 nodes at ρ = 8 |
| **K-3** | **Does the affine `K_A` converge under refinement?** ρ = 8, `n` = 3 → 17. | **measured: 0.03% drift** |
| **K-4** | How much of any drift is the mesh and how much is a law that was never mesh-independent? Four variants: both defaults, each alone, both off. | **measured, and the answer is a null with a cause** |
| **K-5** | The relaxed modulus under refinement. | **not run here** — needs the solver, and `-3647` J-5 is still open |

### K-3 and K-4, and the null result has a cause worth stating

| variant | `n` = 3 | 5 | 9 | 17 | drift |
|---|---:|---:|---:|---:|---:|
| both defaults | 14274.7 | 14277 | 14278.2 | 14278.9 | **0.03%** |
| end-corrected | 14274.7 | 14277 | 14278.2 | 14278.9 | **0.03%** |
| segment steric | 14274.7 | 14277 | 14278.2 | 14278.9 | **0.03%** |
| both fixed | 14274.7 | 14277 | 14278.2 | 14278.9 | **0.03%** |

**All four are identical, and that is not a coincidence.** The *affine* modulus is a virial at fixed
affine displacement, and on this construction:

* the filaments are built **straight**, so an affine dilation keeps them collinear and the bending
  energy is **exactly 0** — the quadrature has nothing to be wrong about (`-3639`'s amendment records
  the same fact retracting a 2× claim);
* the shell is built **overlap-free**, so the segment-pair count is **0 at every `n`** and the steric
  term is exactly 0 under either law.

So the affine estimator cannot see either of the two things §2 says would corrupt a refinement study
— which is the **fourth** time this instrument has turned out to be blind to something a reader
would assume it measured (`-3642` §0 connectivity, `-3644`'s retraction interpenetration, `-3645`'s
H-2 the shell span, and now the mesh).

**The affine `K_A` converges to 0.03% across a 5.7× refinement.** That is a real W-3 pass for that
observable and it is a weak one, because the observable is nearly insensitive to the discretisation
by construction.

## 4. What is NOT claimed

- **Not that the cortex has a mesh-independent modulus.** K-3 is the *affine* bound. The relaxed
  modulus is where bending and contact actually act, and it is not measured here.
- **Not that `n` = 3 was wrong.** It was unmeasurable. It now has a control at 5, 9 and 17 beside it.
- **Not that the two defaults are fixed.** Both stay where they are; this entry makes them
  *selectable* and measures what selecting them costs.
- **Not that every hard-coded 3 is gone.** `aleph/scenarios/spread.py` takes
  `cortex_map.filament_nodes[:, 1]`, which means *the middle node of a three-node filament* and is
  wrong at any other count. That file is this lane's and the fix is a separate change with its own
  control; it is named here so it is not discovered again.

---

## Amendment, 2026-08-07 — **K-5 ran, and it cannot answer its own question yet**

`4090-1`, citation `5f969e5e`, campaign job 60 study B. END_CORRECTED + device-side filament-pair
steric, so neither of §2's two defaults is in the way. CG at 1e-3 pN, budget 200,000.

| ρ | `n` | nodes | `K_A` affine | **`K_A` relaxed** | iterations | **rms** | seg pairs | closest/σ |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 3 | 7,479 | 15011.08 | **175.2095** | 11,465 | 0.00085 | 1 | 1.0996 |
| 8 | 5 | 12,465 | 15013.74 | **176.2494** | 3,495 | 0.00267 | 1 | 1.1034 |
| 8 | 9 | 22,437 | 15015.13 | **176.5053** | 11,955 | 0.00242 | 3 | 1.1007 |
| 24 | 3 | 22,440 | 49118.84 | **1436.165** | 34,583 | 0.00086 | 13 | 1.0111 |
| 24 | 5 | 37,400 | 49129.53 | **1448.994** | 11,359 | 0.00108 | 10 | 1.0341 |
| 24 | 9 | 67,320 | 49138.71 | **1496.489** | **4,357** | **0.00847** | **40** | **0.7814** |

Drift `n` = 3 → 9: **+0.740%** at ρ = 8 and **+4.200%** at ρ = 24, both above the ±0.5% run-to-run
spread `ALEPH-PORT-3645` measured — so both are *resolved* in the sense that the solver's own noise
does not explain them.

### The confound, and it is in the column nobody would look at

**The iteration count FALLS as the mesh is refined.** 34,583 → 11,359 → 4,357 at ρ = 24, on problems
of 22,440, 37,400 and 67,320 nodes. A three-times-larger problem converging in an eighth of the
iterations is not a problem getting easier; **it is the line search giving up earlier**, and the RMS
residual confirms it — 0.00086 → 0.00108 → **0.00847**, ten times worse at the finest mesh.

And the ρ = 24, `n` = 9 row has **40 segment pairs at 0.7814 σ**: the shell it produced is
interpenetrating, which the other five rows are not. **A `K_A` measured on that configuration is not
comparable to one measured on the others**, and it is the row carrying most of the +4.2%.

### So K-5 is recorded as **inconclusive**, not as a convergence result

The affine modulus converges (K-3, 0.03% across a 5.7× refinement). Whether the **relaxed** modulus
does is still unknown, because this measurement varied the mesh *and* the solver's stopping point
together and cannot attribute the difference.

**The honest statement:** the relaxed areal modulus rises by 0.7–4.2% between `n` = 3 and `n` = 9,
and at least the larger figure is contaminated by a finer-mesh run that stopped 8× earlier at 10×
the residual with a shell that had begun to interpenetrate.

### What would answer it, and it is the method `-3645` already used

Drive each node count over the same three initial steps and three tolerances that `-3645`'s W-3 run
used, and compare the **best-converged run at each `n`** rather than one run each. That is what made
the `-3644` G-6 correction possible — the published `+2.785%` became `−0.838%` once the two sides
were compared at matched RMS.

**This is not a re-gate.** K-5 asked "does the relaxed `K_A` converge under refinement" and carries no
threshold; it has not been answered, and saying so is the answer this pass produces.
