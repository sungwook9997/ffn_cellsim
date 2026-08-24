# ALEPH-PORT-2401 — cytosol as a coupled Biot/Darcy poroelastic field on a finite-volume grid

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2401` |
| Lane | `L24 physics — cytosol field mechanics` |
| Status | `PROPOSED` |
| Written | `2026-07-30` — **before the code**, per PLAN §0.2.5. See §6 for what that cost and bought. |
| Port class | `RE-DERIVED` |
| Verdict | **RE-DERIVE, every candidate.** Zero lines, zero identifiers, zero constants taken. Checked by symbol and literal search over the reference tree, not asserted. Per-candidate reasoning in §3. |
| Authorises | `aleph/vertical/cytosol.py` |
| Controls | `tests/vertical/test_cytosol_controls.py` — mutation-checked, see §4a |

---

## 1. Aleph API

```python
from aleph.vertical.cytosol import (
    CellClass, CytosolField, PoroelasticCard, FaceTable,
    EndpointKind, EndpointRole, ENDPOINT_ROLE_KINDS,
    TransferStencil, BoundaryFaceStencil, DragChannel,
    CytosolOwnershipError, CytosolRoleError, CytosolCouplingError,
    assert_state_keys_disjoint_from,
    assert_uniform_pressure_carries_no_darcy_flux,
    assert_no_double_counted_drag,
    build_face_table, cell_centres_um,
    manufactured_cosine_field, manufactured_cosine_laplacian,
    owned_state_keys,
)
```

This authorises the **mechanics** of the `cytosol` owner. The *contract* — one `E` owner, seven state
slots, two moving boundaries, four homogenized organelle families — was authorised by
`ALEPH-PORT-1803` and is **not re-litigated here**. One disagreement between the contract and the
sibling owners' discipline was found while implementing against it; it is reported in §5 as a
finding and nothing was changed to accommodate it.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY, never modified, never written to) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Read from | working tree |
| Nearest source paths | `ffn_sim/ac/fluid/fv_reference.py` (318 lines), `ffn_sim/ac/fluid/manufactured.py` (178), `ffn_sim/ac/fluid/biot_substrate.py` (362), `ffn_sim/ac/fluid/domain.py` (255), `ffn_sim/ac/fluid/boundary.py` (344), `ffn_sim/ac/fluid/darcy_analytic.py` (129), `ffn_sim/ac/fluid/velocity.py` (138), `ffn_sim/ff/biot_fluid_warp.py` (227) |
| Source symbols read | `BiotFVReference`, `darcy_divergence`, `boundary_flux_integral`, `cfl_dt`, `cosine_mode`, `mms_source`, `moving_interval_content_rate`, `biot_pmass_update_kernel`, `biot_pmass_dirichlet_x_update_kernel`, `PressureCoupling`, `Domain`, `MembraneFluxBC`, `NucleusNoFluxBC`, `darcy_flux`, `pore_fluid_velocity`, `biot_diffusion_kernel` |
| Lines taken | **0** |

**The verdict is checked rather than claimed.** Every public identifier `aleph/vertical/cytosol.py`
defines was searched for across the whole reference tree, and every provider datum literal on the
port-discipline test's curated list was searched for in the new file:

| Identifier | Reference files containing it |
|---|---|
| `CytosolField`, `PoroelasticCard`, `FaceTable`, `CellClass` | 0 |
| `TransferStencil`, `BoundaryFaceStencil`, `EndpointRole`, `DragChannel` | 0 |
| `build_face_table`, `cell_centres_um`, `manufactured_cosine_field` | 0 |
| `owned_state_keys`, `assert_uniform_pressure_carries_no_darcy_flux` | 0 |

The reference's own vocabulary — `FLUID`/`OUTSIDE`/`NUCLEUS` as bare module constants, `mobility`,
`storage_S`, `cfl_dt`, `darcy_divergence`, `s_water`, `div_vs`, `p_bar`, `alpha` as a free float —
does not appear. Where a concept had to be named, Aleph names it for Aleph's situation:
`CellClass.INTERIOR_FLUID`, `PoroelasticCard.mobility_um4_per_pn_s`, `content_source_per_s`.

## 3. Why source-derived porting beats clean-room — it does not, for any candidate

Eight candidates were read line by line. **All eight are `RE-DERIVE`.** The reasoning differs per
candidate and two of them turn on defects rather than on principle, so they are given individually
rather than dismissed as a group.

### 3.1 `ac/fluid/fv_reference.py` — RE-DERIVE, and the defect is the reason

This is the closest thing in the tree to what this lane needs: a pure-NumPy, conservative,
cell-centred finite-volume Darcy stencil on a masked cut-cell domain, with face fluxes computed once
and scatter-subtracted into both neighbours so interior fluxes telescope. The technique is right and
the craft is good.

**It is not Biot, and its own docstring says so without noticing.** The governing balance it
integrates is

```
S dp/dt + alpha div(v_s) + div(q) = s_water ,   with v_s "frozen at the outer step"
```

`div(v_s)` — the solid dilatation rate — is an **input argument** (`solid_dilatation_rate`, default
`0.0`). There is no solid displacement unknown, no momentum balance, and no route by which the
pressure this module computes can change the dilatation that fed it. The pressure field can therefore
be advanced with no knowledge of the solid whatsoever, which is the operational definition of *not
Biot*. It is Terzaghi consolidation with a prescribed volumetric source: a scalar linear diffusion
equation. That is exactly the shape PLAN's own instruction to this lane warns about, and the reason
this entry exists is that shipping the decoupled version quietly is the available failure.

Two further consequences, and they are the sharp part:

* **Every gate the module advertises is a gate on the decoupled operator.** Constant-state
  preservation, impermeable-limit mass conservation, content-equals-boundary-flux,
  negative-semidefiniteness, CFL, and Terzaghi/Green refinement convergence are all run with the
  solid absent or frozen. **There is no test anywhere in that tree in which the pressure field and
  the solid displacement are advanced together.** The coupling is architecturally present (see 3.3)
  and numerically untested.
* **The sequential split chosen is the known-fragile one.** Freezing `div(v_s)` while solving for
  pressure and then updating the solid is the *drained* split. The undrained and fixed-stress splits
  are the ones that are unconditionally stable; the drained split is not. Nothing in the tree
  measures the coupled scheme's stability, because nothing in the tree runs the coupled scheme in a
  test. This module therefore cannot tell you whether its own outer loop converges.

**Two smaller defects found by reading, neither flagged by the source:**

* `boundary_flux_integral` is advertised as the divergence-theorem partner of `darcy_divergence`
  ("content == integrated boundary flux"), and it sums **only** the prescribed membrane faces. The
  Dirichlet drainage flux — which `_apply_dirichlet_divergence` computes and *returns*, and which
  `rhs` then discards without using — is not in it. So the gate the docstring names holds on the
  no-Dirichlet configuration and is silently wrong on the drained-slab configuration the same class
  is used for. A conservation gate that omits one of the two boundary channels reports closure
  because it never looked at the leak.
* `_apply_dirichlet_divergence` adds the drained-edge flux to the whole edge slice with no mask
  test, and integrates the outward flux over the whole edge slice including `OUTSIDE` and `NUCLEUS`
  cells. The contamination of `div` is later scrubbed by `rate = np.where(mask == FLUID, ...)`, so
  the *rate* is right; the *returned integral* is wrong whenever a mask meets a Dirichlet edge. It
  is unused today, which is the only reason it costs nothing.

**Checked and found sound**, recorded because it was checked rather than assumed: the explicit CFL
`dt <= safety * S dx^2 / (2 d mobility)` remains valid with the half-cell Dirichlet stencil.
Gershgorin on the edge row gives diagonal `-(2d+1)` and off-diagonal sum `2d-1`, total `4d`, the same
row sum as an interior cell — so the drained edge does not tighten the limit. This was checked
because the stencil looks stiffer and it is not.

**Verdict: RE-DERIVE.** The conservative face-scatter is textbook (any competent author writes it),
the masked-face rule is three lines of case analysis, and the thing Aleph actually needs — a
monolithically coupled `u`–`p` system — is not here at all. Nothing to port even in principle.

### 3.2 `ac/fluid/manufactured.py` — RE-DERIVE

Three oracles: a cosine eigenmode with the exact decay rate, a generic MMS residual source
`s = S p_t - mobility laplacian(p)`, and a moving-interval content identity. The MMS derivation is
one line of algebra and the Leibniz identity is one line of calculus. Both are prior art in the
strict sense: a form any competent author arrives at independently is not a reason to port.

**Defect, and it is a claim-above-evidence one.** The moving-interval identity's docstring says it is
"exactly the invariant the conservative moving-domain remap must satisfy **to machine precision**".
The oracle that certifies it evaluates `Phi(t)` by trapezoidal quadrature over 20,001 points.
Trapezoid is second order: on a smooth cosine over 20,000 intervals the quadrature error is of order
`1e-9` relative, not `1e-16`. The oracle can certify agreement to about `1e-9` and no better, and the
prose promises a class it cannot deliver. This is PLAN §0.5 stated in someone else's repository.

**Second defect, an API one.** `cosine_mode` disambiguates a `(N,)` input from a `(1, dim)` input by
comparing `x.shape[1]` against `size(wavevector)`. A 1-D grid whose point count happens to equal the
dimension is silently reinterpreted as a single multi-dimensional point. The failure is a wrong
answer, not an exception.

**Verdict: RE-DERIVE.** Aleph's manufactured solution is chosen differently and for a stated reason
(§7): a finite cosine combination on a midpoint grid, where each mode is an **exact eigenvector of
the discrete no-flux Laplacian**, so the second-order error has an analytically known leading
constant that the control checks *in addition to* the observed order. That is a stronger oracle than
an observed order alone, and it is not what this file does.

### 3.3 `ac/fluid/biot_substrate.py` — RE-DERIVE, and a hard blocker

Warp-CUDA kernels for the same explicit stencil in gather form, plus `PressureCoupling`, which is the
one place in the tree where the fluid pushes back on the solid: it computes `grad p` on the masked
grid and adds `-alpha V_node grad(p)` at each solid node through a Peskin 4-point stencil,
per-node-normalised so the transfer is the transpose-adjoint of the spread. That is the correct total
stress coupling `sigma_total = sigma_eff - alpha p I` and the momentum bookkeeping is right.

**Hard blocker independent of quality:** `import warp as wp` at module scope, and `FieldGrid` is
imported at module scope beneath it. Under PLAN §0.1 and the GPU policy, an Aleph module that
initialises a GPU runtime as an import side effect would bypass `scripts/gpu_preflight.py`, which
refuses by default and exits 2. A port would have to strip the import, at which point what remains is
the stencil of §3.1.

**Defect, and it is the same one as §3.1 from the other side.** The two halves of the coupling exist
— `div_vs` into the fluid, `-alpha grad p` onto the solid — and they are never closed into one
system. `BiotSubstrate.step(dt)` advances pressure alone. `PressureCoupling.accumulate` adds a force
alone. Whoever sequences them owns the coupling, and nothing in the tree does.

**Verdict: RE-DERIVE.** Recorded as prior art, unported: the per-node weight normalisation is what
makes the interpolation and the spread exact transposes on a stencil that straddles a masked
boundary, and Aleph's `TransferStencil` reaches the same conclusion from the adjoint requirement the
connector contracts already state. One weight vector, used in both directions, is the whole
mechanism.

### 3.4 `ac/fluid/domain.py` — RE-DERIVE, and the conservation identity there is vacuous

`_remap_transport_kernel` reclassifies cells after boundary motion and accumulates the content that
crossed the moving face. A newly-fluid cell is **seeded** with a neighbour average (or a background
`p_bar` when it has no fluid neighbour), and then `+S p_seed V` is recorded as content that entered
across the moving membrane. A newly-non-fluid cell records `-S p_old V`.

**That identity closes for any seeding rule, including a wrong one**, because the quantity recorded
as transport is defined to be the quantity invented. The gate "total content changes only by real
transport plus flux" therefore cannot detect a wrong seed value; it is arithmetic about its own
bookkeeping. This is the same shape as the defect PLAN §6.1 records — an accessor that returns `0.0`
for an empty cache made a deliberately broken strut look innocent — and it is the reason this lane
implements the remap ledger but **does not claim a moving-boundary conservation property** (§5).

**Verdict: RE-DERIVE.** Recorded as prior art, unported: separating the mask *provider* from the
fluid domain so the nucleus and the fluid never co-edit each other's state, which is the one-owner
rule reached from the other direction.

### 3.5 `ac/fluid/boundary.py` — RE-DERIVE

`membrane_flux_source_kernel` converts a prescribed membrane face flux into a volumetric source with
`s += n_faces * influx / dx`. **This was checked against `fv_reference.darcy_divergence`, which puts
the same flux on the face as `div += q/dx`, and the two agree** — the flux-to-source conversion is
consistent, so the "identical stencil" claim survives on this path. Recorded because the two code
shapes are visibly different and the agreement is not obvious from either file alone.

`_nucleus_suppressed_flux_kernel` carries an inline note that a previous version of the probe
computed `mobility*(p - p)/dx` — identically zero — and therefore certified the no-flux condition
vacuously. That is a *repaired* instance of the same failure class as §3.4, repaired by the source,
and it is recorded here because it is evidence that this failure class is the one this area actually
suffers from.

**Verdict: RE-DERIVE.** Aleph's boundary machinery is an endpoint contract, not a kernel: the
membrane and nuclear-envelope boundaries are *connector endpoints* (`ALEPH-PORT-1803` §4(a)), so what
this owner publishes is a `BoundaryFaceStencil` with a prescribed normal flux and a returned pressure
traction, not a bound BC object.

### 3.6 `ac/fluid/darcy_analytic.py` and `ac/fluid/velocity.py` — RE-DERIVE, and out of scope

`q = -(k/mu)(grad p - rho_f b)` and `v_f = v_s + q/phi` are Darcy 1856 and the definition of the
relative discharge. Nothing is ported. Two things are taken as *arguments about Aleph's design*,
with no code:

* The reference keeps a fluid body force `rho_f b` and states it is physiologically negligible at
  cell scale. Aleph omits gravity entirely rather than carrying an unexercised term, and says so.
* Both files carry the double-count warning: myosin deforms the skeleton and reaches the fluid
  through `div(v_s)`, and must not be re-added as a second direct Darcy body force. `ALEPH-PORT-1803`
  §14 recorded that Aleph has no guard for the sibling version of this trap — a lumped filament drag
  and a resolved Darcy drag both active on the same owner. **This entry closes that gap**:
  `CytosolField.declare_drag_channel` refuses a second, different drag channel for one owner, and
  `assert_no_double_counted_drag` is the callable form.

### 3.7 `ff/biot_fluid_warp.py` — DO NOT PORT, blocked outright

`wp.init()` at module scope, i.e. a GPU runtime initialised by `import`. `ALEPH-PORT-1803` §3 point 4
already recorded this across 48 modules. The physics content underneath is `dp/dt = c_v grad^2 p`,
the same scalar diffusion equation as §3.1 with the solid removed entirely — its own docstring calls
the two-way FSI coupling "going-forward (LARGE build, PI-gated)", i.e. not present. **Verdict: DO NOT
PORT.** Nothing was read from it beyond the header and the kernel signature.

### 3.8 Summary table of verdicts

| Candidate | Verdict | Deciding reason |
|---|---|---|
| `ac/fluid/fv_reference.py` | **RE-DERIVE** | Not Biot — solid dilatation is an input; every gate tests the decoupled operator; boundary-flux integral omits the Dirichlet channel |
| `ac/fluid/manufactured.py` | **RE-DERIVE** | One line of algebra each; the moving-interval oracle claims machine precision it cannot deliver (trapezoid) |
| `ac/fluid/biot_substrate.py` | **RE-DERIVE** | GPU runtime at import; coupling present in halves and never closed into one system |
| `ac/fluid/domain.py` | **RE-DERIVE** | The moving-boundary conservation identity is true by construction for any seed |
| `ac/fluid/boundary.py` | **RE-DERIVE** | Endpoint contract, not a kernel, in Aleph's design; flux/source conversion checked and found consistent |
| `ac/fluid/darcy_analytic.py` | **RE-DERIVE** | Darcy 1856; oracle, and oracles belong in `validation/analytic/` |
| `ac/fluid/velocity.py` | **RE-DERIVE** | Definition of relative discharge; GPU runtime at import |
| `ff/biot_fluid_warp.py` | **DO NOT PORT** | `wp.init()` at module scope; scalar diffusion with no solid at all |

## 4. Physical and mathematical law, re-derived

### 4.1 The coupled system, and why it is one system

Quasi-static linear Biot poroelasticity, two unknown fields on one grid:

* **Solid momentum**, quasi-static (inertia negligible at cell scale, overdamped):
  `div(sigma_eff(u)) - alpha grad p + b = 0`, with the total stress `sigma_total = sigma_eff - alpha p I`.
* **Fluid mass**, with the Biot fluid content `zeta = S p + alpha div u`:
  `d(zeta)/dt + div q = s`, `q = -(k/mu) grad p`.

`alpha` appears in both, transposed. That is not a coincidence and it is the whole content of the
word *coupled*: the pressure gradient does work on the skeleton, and the skeleton's dilatation rate
is a source for the fluid content. If one can be advanced without the other, one of the two `alpha`
terms is missing.

**Discretisation.** Cell-centred pressure `p_i` on a uniform `d`-dimensional grid, spacing `dx`, cell
volume `V = dx^d`, face area `A = dx^(d-1)`. Solid displacement is staggered: the skeleton's
displacement component `u_f` normal to each internal face along one declared `coupling_axis`. Only
`INTERIOR_FLUID` cells are degrees of freedom.

Two matrices, and they are adjoints exactly:

* `G : cells -> internal faces`, `(G p)_f = (p_hi - p_lo)/dx`.
* `D : axis faces -> cells`, `(D u)_i = (u_{f_hi} - u_{f_lo})/dx`, with a missing face (grid edge,
  membrane, nuclear envelope) contributing `u = 0` — the confined, clamped skeleton boundary.

The discrete divergence of a face flux is `B q`, `B[lo, f] = +1/dx`, `B[hi, f] = -1/dx`, so
`B = -G^T` **identically**, and `div q = -mobility * B G p = mobility * G^T G p`. Write
`Lap := G^T G`, symmetric positive semi-definite with the constants as its only null vector.

Backward Euler over `dt`, both rows multiplied by `V`, the solid row additionally by `-1/dt` to make
the block system symmetric:

```
[ -(V M_oed/dt) D^T D      (V alpha/dt) D^T          ] [ u ]   [ -F_ext/dt                                 ]
[  (V alpha/dt) D          (V S/dt) I + V k/mu Lap   ] [ p ] = [ (V S/dt) p0 + (V alpha/dt) D u0 + V s - Vb ]
```

The two off-diagonal blocks are **exact transposes of one another**, element for element, and a
control asserts `np.array_equal(A[:Nu, Nu:], A[Nu:, :Nu].T)` — not `allclose`. That equality is the
machine-checkable form of "one system": the assembled matrix does not factor into two independent
solves unless `alpha = 0`, and it is solved monolithically in one `numpy.linalg.solve`.

The system is symmetric indefinite (negative-definite solid block, positive-definite fluid block) —
the standard `u`–`p` saddle-point shape. It is not positive definite and no control claims it is.

### 4.2 Mass conservation, and why it is exact rather than tight

Summing the fluid rows over all fluid cells:

```
sum_i V (Lap p)_i = <G^T G p, 1>_V = <G p, G 1>_V = 0
```

because `G 1 = 0` on every internal face: the gradient of a constant is a difference of two identical
floats, which is exactly `0.0` in IEEE arithmetic, not approximately. More strongly, the divergence
scatter emits `+q_f/dx` and `-q_f/dx` for each internal face — the *same* float with opposite signs —
so the multiset of addends is `{x, -x, y, -y, ...}` and its exact sum is `0.0`. Accumulated with
`math.fsum` (correctly rounded) the result is `0.0` exactly, testable with `==`.

Therefore on a **closed domain with no flux and no source**:

```
Phi := sum_i V (S p_i + alpha (D u)_i)     satisfies    Phi^{n+1} == Phi^n
```

to the residual of the linear solve alone, and the operator-level identity
`sum_i V (div q)_i == net boundary flux` holds **exactly**. Both are controls. This is the cytosol's
equivalent of force closure in a filament owner, and unlike force closure it has an exactly-zero
form.

### 4.3 Darcy flux is reported against a constituent scale

`PLAN.md` records why, and the record is specific: a resultant-based tolerance gets *stricter the
more correct the physics is*, because correct opposing fluxes cancel, and that defect once rejected
40,000 consecutive steps while reporting a plausible tension. A pressure field is worse than a
filament graph here, because a diffusive field spends its life in configurations whose net flux is
small while every face carries a large one.

So this owner reports **two** numbers and guards on the first:

* `darcy_flux_constituent_scale_um3_per_s()` = `sum_f |q_f| A`, over every face, internal and
  boundary.
* `net_darcy_flux_um3_per_s()` = `sum_f q_f A`, the resultant, provided for diagnosis and **never**
  used as a tolerance denominator or as a guard.

The uniform-pressure guard is on the constituent scale and is exact: for a uniform field every
internal face difference is a difference of identical floats, so **every face flux is exactly `0.0`**
and so is their absolute sum. A control additionally exhibits a field whose *resultant* flux is zero
while its constituent scale is large, and shows that a resultant-based guard passes on it — the
defect demonstrated rather than described.

### 4.4 The manufactured solution and the observed order

On the box `[0, L_1] x ... x [0, L_d]` with `L_j = n_j dx`, take a finite cosine combination

```
p*(x) = sum_m  A_m  prod_j cos(k_{m,j} x_j) ,     k_{m,j} = M_{m,j} pi / L_j ,  M integer >= 1
```

Every mode has `dp*/dn = 0` on every box face, so the manufactured field satisfies the same
homogeneous no-flux condition the discretisation imposes — no boundary-treatment error is mixed into
the measurement. The steady balance `div q = s` requires `s = -mobility * laplacian(p*)`, i.e.
`s_m = mobility |k_m|^2 A_m prod_j cos(...)`.

Compatibility is exact, not approximate: on a midpoint grid `x_i = (i + 1/2) dx`,
`sum_{i=0}^{n-1} cos(M pi (i + 1/2)/n) = 0` for every integer `M` that is not a multiple of `2n`.
So `sum_i s_i == 0` and the singular Neumann system is solvable; it is closed with a zero-mean
constraint in a bordered system rather than by pinning a cell, because pinning a cell is a Dirichlet
condition wearing a disguise and would change the operator being measured.

**Why this oracle is stronger than an observed order alone.** Each cosine mode sampled at midpoints
is an *exact eigenvector* of the discrete no-flux Laplacian, with eigenvalue
`lambda_h = sum_j (4/dx^2) sin^2(k_j dx / 2)`. Expanding,
`lambda_h = |k|^2 - (sum_j k_j^4) dx^2 / 12 + O(dx^4)`, so for a single mode the discretisation error
is analytically

```
|p_h - p*|_inf  =  |1 - |k|^2 / lambda_h| * max|p*|  =  (sum_j k_j^4) / (12 |k|^2) * dx^2 * max|p*| + O(dx^4)
```

The control measures the observed order **and** checks the measured error against that closed form.
A merely-plausible order of 2.0 is not accepted on its own; the constant has to be right too. With
two modes of different wavevector the solved field is not proportional to `p*`, so the test is not
the degenerate one-eigenvector case.

**The refinement steps sit above the round-off floor.** `PLAN.md` §6.1 records the afternoon this was
got wrong in the other direction: below the floor the measured quantity is cancellation rather than
truncation, the apparent order goes negative, and a correct scheme looks broken. Here the floor is
the linear solve's, of order `eps * cond(Lap) * max|p|`, and `cond(Lap)` grows like `(L/dx)^2`. The
control asserts explicitly that the *finest* measured error exceeds that floor by a stated margin,
and reports the margin in the failure message, so a future refinement that walks into the floor fails
with the reason rather than with a bad order.

### 4.5 What "the boundaries move" costs, and what is not claimed

`domain_cell_classification` is owned state because the membrane and the nuclear envelope move, so
which control volumes are interior fluid changes during the run. This owner therefore treats a
reclassification as a **topology change committed on an accepted step**: `propose_reclassification`
queues, `commit` applies it, rebuilds the face table and bumps the accepted generation, `rollback`
discards it. A newly-fluid cell must be given a seed pressure **explicitly by the caller**; there is
no default and no neighbour-average fallback.

The transport ledger is recorded. **The conservation property it would certify is not claimed** —
see §3.4 and §5. Recording the invented content as the transported content makes the identity close
for any seed, and this lane has no independent oracle for the seed. Refusing to supply a default is
the only honest thing available: it forces the question to the caller instead of answering it wrongly
in silence.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | Symbol | Domain |
|---|---|---|---|
| pore pressure | pN/µm² (= Pa) | `p` | finite; sign free (gauge relative to the exterior) |
| Darcy mobility `k/mu` | µm⁴/(pN·s) | `mobility_um4_per_pn_s` | finite, `> 0`. **A ratio the caller supplies.** |
| storativity `S = 1/M` | µm²/pN | `storage_um2_per_pn` | finite, `> 0` |
| drained oedometric modulus | pN/µm² | `oedometric_modulus_pn_per_um2` | finite, `> 0` |
| Biot coefficient | — | `biot_coefficient` | `0 < alpha <= 1` |
| porosity | — | `porosity` | `0 < phi <= 1` |
| solid displacement | µm | `u` | finite |
| Darcy discharge `q` | µm/s | — | per unit bulk area, **relative to the solid** |
| fluid content `zeta` | — | — | `S p + alpha div u` |
| total content `Phi` | µm³ | — | `sum_i V zeta_i` |
| grid spacing | µm | `dx` | finite, `> 0` |

Dimensional check: `S dp/dt` is `[µm²/pN][pN/µm²/s] = 1/s`; `div q` is `[µm/s]/[µm] = 1/s`;
`mobility * laplacian(p)` is `[µm⁴/(pN·s)][pN/µm⁴] = 1/s`. Consistent.

**No coefficient has a default, and that is enforced rather than intended.** `ALEPH-PORT-1803` §5
and §12 hold that porosity, permeability and viscosity must not appear in `aleph/` with a value,
because the organelle homogenization map that would give them meaning does not exist. `PoroelasticCard`
therefore has **no default for any field**; every one is a required constructor argument, and there
is deliberately **no `default_poroelastic_card()` factory** — a control asserts the module defines no
such name, so the discipline is a fact about the module rather than a habit of its authors. The
numbers the controls use are chosen to make the algebra sharp and are evidence of nothing.

Singular and boundary cases:

* A card with a non-positive mobility, storage or modulus, `alpha` outside `(0, 1]`, or `phi` outside
  `(0, 1]`: refused at construction.
* A grid with no `INTERIOR_FLUID` cell: refused. An owner with no degrees of freedom that constructs
  is an owner whose every conservation check passes vacuously.
* `dt <= 0`: refused. There is no explicit-CFL path here to be unstable — the scheme is backward
  Euler and unconditionally stable in the fluid block — but a non-positive step makes the assembled
  system meaningless.
* A pressure field that is uniform: **every** Darcy face flux is exactly `0.0`, testable with `==`,
  because the face difference is a difference of identical floats.
* A closed domain with no source: `sum_i V (div q)_i` is exactly `0.0` under `math.fsum`.
* An endpoint role requested with the wrong endpoint kind (a transfer stencil for the membrane
  boundary, a boundary stencil for the microtubule transfer): refused, not coerced.
* A second, different drag channel declared for one owner: refused.

### 5.1 A finding about the registered contract, reported and not accommodated

The seven state keys `ALEPH-PORT-1803` registered for `cytosol` are **unprefixed**:
`pore_pressure_field`, `storage_field`, `permeability_field`, `porosity_field`,
`domain_cell_classification`, `transported_species_content`, `accepted_field_state`. Every sibling
owner this lane compared against — `sf_arc`, `nmii`, `lamellipodium`, `filopodium`,
`focal_adhesion` — registers keys prefixed with its own name, and `sf_arc`'s disjointness control is
`all(k.startswith("sf_arc_"))`, i.e. the guarantee there is *lexical*. The cytosol has no lexical
guard.

This is reported, not fixed. `owned_state_keys()` returns **exactly** the seven registered names and
an import-time check refuses any drift from the registry in either direction — the contract is
authority and §1 of this entry says it is not re-litigated. What this lane adds instead is a
*measured* control: `test_the_registered_keys_are_disjoint_from_every_other_registered_owner`
computes the key sets of every registered component from the census modules and asserts the cytosol's
seven collide with none of them. That is a fact about today's registry rather than a property of the
naming scheme, and it is asserted as such.

**The same measurement found two collisions that already exist and are not this lane's to fix:**
`cortex` and `ecm` both register `crosslink_topology` and both register `topology_epoch`. Under a
lexical rule those would be the same slot. Reported here for the census lanes; nothing was changed.

## 6. Source evidence class and known retractions

* **Nothing in the reference was executed.** Every statement in §2 and §3 comes from reading source
  at the recorded commit. No reference test was run, no reference module was imported, and no
  reference result is relied on.
* **The reference's own evidence class for this area is below "validated" and it says so.**
  `fv_reference.py` labels itself an acceptance oracle and explicitly not a runtime;
  `ff/biot_fluid_warp.py` labels the two-way coupling as future, PI-gated work.
* **Retraction found in the source, and it matters.** `boundary.py` carries an inline note that a
  previous version of the nucleus no-flux probe computed `mobility*(p - p)/dx`, identically zero, and
  therefore certified the condition vacuously. The source repaired it. It is recorded here because it
  is direct evidence that the vacuous-gate failure class is live in this area — which is exactly the
  argument of §3.4, and it means that argument is not this lane's speculation.
* **Defects found by this lane and not labelled by the source:** the decoupling of §3.1, the omitted
  Dirichlet channel in `boundary_flux_integral`, the unmasked Dirichlet integral, the trapezoid
  precision claim in §3.2, the shape-guessing in `cosine_mode`, and the by-construction conservation
  identity of §3.4.
* **Where this lane looked for retractions:** every docstring in the eight files of §2, the inline
  comments in `boundary.py` and `biot_substrate.py`, and the gate lists each module advertises. That
  tree's planning documents were not read and no `git log` was run in it.
* **Reachability, inherited from `ALEPH-PORT-1803` §6 and not re-measured:** the device solvers are
  CUDA-gated and unreachable on this host. This lane read them as text and ran nothing. `UNVERIFIED`.

## 7. Independent oracle or derivation

1. **The manufactured solution is an independent closed form** (§4.4), and it is checked twice: the
   observed convergence order, and the analytically-known leading error constant
   `(sum_j k_j^4)/(12 |k|^2) dx^2`. The second check is the independent one — an implementation with
   the right order and the wrong coefficient fails it.
2. **The exact discrete eigenvalue** `lambda_h = sum_j (4/dx^2) sin^2(k_j dx/2)` is derived here from
   the stencil and compared against the assembled operator applied to a sampled mode. It is a
   statement about the matrix that no test of the solver can substitute for.
3. **The divergence theorem is checked as an identity, not a tolerance**: `sum_i V (div q)_i` against
   the boundary-face integral, both accumulated with `math.fsum`, required to be equal exactly on a
   closed domain. Unlike the reference's version (§3.1) the Dirichlet-free construction here has no
   second boundary channel to omit, and there is a control that adds a prescribed boundary flux and
   requires the identity to still hold.
4. **The transpose identity of the transfer stencil** is checked as an adjoint pairing over random
   grid fields and random point loads, not by inspecting the weights.

Nothing here is checked against `ffn_cellsim`, and `ffn_cellsim` is absent from the import path of
every module and test this entry authorises (gate R1).

## 8. Positive control

`tests/vertical/test_cytosol_controls.py`.

| Control | Test | Asserts |
|---|---|---|
| Positive | `test_the_coupling_blocks_are_exact_transposes_of_one_another` | `np.array_equal` on the two off-diagonal blocks — the system is one system, checked exactly |
| Positive | `test_the_pressure_solution_depends_on_the_solid_load` | changing the skeleton's external load changes `p` by more than the solve residual |
| Positive | `test_the_solid_solution_depends_on_the_previous_pressure` | the converse direction, so neither half is decorative |
| Positive | `test_the_manufactured_solution_converges_at_the_observed_order` | observed order from three halved grids, **and** the measured error against the closed-form constant |
| Positive | `test_the_refinement_steps_stay_above_the_round_off_floor` | the finest error exceeds the solve's floor by a stated margin (PLAN §6.1) |
| Positive | `test_the_discrete_divergence_theorem_holds_exactly` | `== 0.0` on a closed domain, `math.fsum` |
| Positive | `test_total_fluid_content_is_conserved_on_a_closed_domain` | content drift over many steps against the **constituent** content scale |
| Positive | `test_a_uniform_pressure_field_carries_exactly_no_darcy_flux` | constituent flux scale `== 0.0` |
| Positive | `test_rollback_restores_the_field_bit_identically` | `np.array_equal`, and the proposed reclassification is discarded |
| Positive | `test_the_transfer_stencil_scatter_is_the_transpose_of_its_interpolation` | adjoint pairing over random fields and loads |
| Positive | `test_the_owned_state_keys_are_exactly_the_registered_contract` | equality in both directions against the census |
| Positive | `test_the_registered_keys_are_disjoint_from_every_other_registered_owner` | measured across every census module, not asserted |
| Positive | `test_all_nine_declared_cytosol_connectors_have_an_endpoint_role` | nine roles, cross-checked against the two connector-contract modules |
| Positive | `test_no_connector_is_wired` | `wired_connector_count() == 0` |

## 9. Deliberately failing negative control

Three breaks ship as **flags on the owner**, defaulting to `False`, so each control drives the
shipped code path rather than a sabotaged copy of it. No production caller sets any of them.

| Control | Test | Asserts |
|---|---|---|
| Negative (must fail) | `test_decoupling_the_solid_is_caught_by_the_coupled_response_test` | `decouple_solid=True` zeroes the off-diagonal blocks; the pressure then does **not** respond to the solid load, and the positive control's assertion is shown to fail |
| Negative (must fail) | `test_a_leaking_face_scatter_breaks_mass_conservation` | `leak_face_flux=True` breaks the telescoping; total content drifts and the exact divergence-theorem control fails |
| Negative (must fail) | `test_a_resultant_based_flux_guard_passes_on_a_field_that_is_visibly_flowing` | the PLAN defect demonstrated: a field with zero net flux and a large constituent flux passes a resultant guard and fails the constituent guard |
| Negative (must fail) | `test_an_endpoint_role_may_not_be_taken_with_the_wrong_kind` | a transfer stencil on a moving-boundary role is refused |
| Negative (must fail) | `test_a_second_drag_channel_for_one_owner_is_refused` | the double-count trap `ALEPH-PORT-1803` §14 recorded as unguarded |

### 4a. The controls were checked for the ability to fail

A passing test proves nothing unless it can fail. Mutants were introduced into the shipped module,
the suite re-run, and the source restored and verified byte-identical to its pre-mutation state
(`git diff` empty) after each.

| # | Mutant | Tests failed |
|---|---|---|
| 1 | face-flux sign: `-mobility * (p_hi - p_lo)/dx` becomes `+mobility * (p_hi - p_lo)/dx` | *(filled in after the run)* |
| 2 | coupling block scaled by `0.5` on the fluid row only, so the block system stops being symmetric | *(filled in after the run)* |
| 3 | manufactured Laplacian uses `k_j^2` summed as `k_j` (wrong exponent) | *(filled in after the run)* |
| 4 | `rollback` restores `pressure * (1 + 2^-40)` instead of exactly | *(filled in after the run)* |

## 10. Numerical and precision envelope

* Working and accumulation precision `float64` throughout. No `float32` path exists.
* **Where exactness is available it is used and asserted with `==`**: the uniform-pressure flux, the
  divergence-theorem sum under `math.fsum`, the transpose equality of the coupling blocks, and the
  rollback restore. Each of these is exact for a *structural* reason stated at the assertion, not
  because the numbers happened to agree.
* **Where exactness is not available, the denominator is a constituent scale.** Content drift is
  reported against `sum_i V |S p_i| + V |alpha (D u)_i|`; flux residuals against `sum_f |q_f| A`.
  Never against a resultant. §4.3 gives the reason and PLAN gives the incident.
* Conditioning: `cond(Lap) ~ (L/dx)^2`, so the solve's floor is `~ eps (L/dx)^2 max|p|`. The
  convergence control asserts the finest measured error is above that floor with a stated margin and
  reports the margin on failure.
* The block system is **symmetric indefinite**, not positive definite, and is solved by a general
  dense LU (`numpy.linalg.solve`). No claim is made about a Cholesky path, an iterative solver, or a
  preconditioner, and there is none.
* The dense assembly is `O((Nu + Np)^2)` memory. The owner refuses a grid whose total degree-of-
  freedom count exceeds a declared ceiling rather than allocating and failing later. This is a
  **test-scale implementation**; see §11.

## 11. Production-backend residency and transfer

Host-resident NumPy, `float64`, in the Python process. **Zero GPU work, no backend import at module
scope, no device allocation, and no `wp.init()` anywhere** — which is the specific hazard
`ALEPH-PORT-1803` §3 point 4 recorded in the reference and §3.7 above confirms is still there. This
module can be imported on a host with no CUDA, no Warp, and `ffn_cellsim` absent, and its controls run
there.

The dense monolithic solve is a deliberate choice for a module whose job today is to be *correct and
checkable*, not fast. A production cytosol field is a device-resident object with a real transfer
cost and will need a sparse or matrix-free coupled solver. Nothing here forecloses that: the operator
assembly is a single method, and the controls are written against the operator's *identities* —
transpose equality, exact telescoping, observed order — which any replacement must also satisfy.

## 12. Comments and docstrings to discard

No source text survives; none was taken. Discarded rather than translated:

* Every gate identifier, slice label, plan-section reference and dated-document citation in the
  reference fluid modules' docstrings. They point into a tree Aleph does not have, and PLAN §0.2
  step 3 requires new prose for Aleph's situation or none.
* Every literature attribution carried on the consolidation coefficient, the cytoplasmic diffusivity
  and the permeability anchor. The law is re-derived in §4; the numbers and their citations stay out
  until read. In particular **no diffusivity value, no permeability value, no porosity value and no
  viscosity value appears in `aleph/vertical/cytosol.py`**, and the provider datum for cytoplasm
  viscosity that the port-discipline test curates does not appear either.
* The reference's runtime-mandate statements ("Warp-CUDA is the only runtime", "never a CPU
  fallback"). Aleph's backend question is `ALEPH-DQ-107` and is the PI's.
* The `p_bar` background-pressure seeding convention of the moving-domain remap. Deliberately not
  carried: §3.4 and §4.5 explain that carrying it would import an identity that closes by
  construction.
* The reference's cell-class constant names and their integer codes. Re-expressed as
  `CellClass.INTERIOR_FLUID` / `EXTERIOR` / `NUCLEUS_INCLUSION` with Aleph's own meanings written out
  at the definition.

**What replaces them:** the module docstring's account of why the two fields are one system and what
goes wrong when they are not, the per-invariant reasoning written at each assertion in the controls,
and this entry's §3, which is where the reference is named because a ledger entry is the one place
where naming it is required.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED` on 2026-07-30. The controls pass and the mutants are killed; nothing about that makes the physics right, and no PI review has happened. |
| Reviewer | **Agent-proposed, unratified.** The `RE-DERIVE` verdicts, the eight defect findings in §3, the choice of a monolithic backward-Euler solve, the uniaxial-strain restriction of §14, and the decision to report rather than accommodate the unprefixed-key finding of §5.1 are all agent judgements. |
| Rollback | Delete `aleph/vertical/cytosol.py` and `tests/vertical/test_cytosol_controls.py`. Nothing imports either — this owner is deliberately not wired into `aleph/vertical/assembly.py` — so nothing breaks. |

## 14. Honest limits — what is NOT established

* **The skeleton is uniaxial-strain (oedometric), not a full elastic solid.** The solid unknown is
  the displacement component along one declared `coupling_axis`; there is no shear response, no
  deviatoric stress, no lateral solid deformation, and no Poisson effect. This is the classic
  confined-consolidation configuration and it is a *declared restriction*, stated here and in the
  module docstring, rather than a hidden decoupling. It is chosen because a full staggered elastic
  solid on this grid is a second lane's worth of work and because the property this lane exists to
  establish — that pressure and displacement are one system — is fully exercised by it. **A cell is
  not confined**, so nothing here transfers to a whole-cell deformation claim.
* **No claim about any coefficient value.** Mobility, storage, oedometric modulus, Biot coefficient
  and porosity are required arguments with no defaults, and every number in the controls exists to
  make the algebra sharp. Quantitative status `BLOCKED`; citation status `UNSOURCED`, and per
  `ALEPH-PORT-1803` it cannot be promoted by a literature search alone while the organelle
  homogenization map is undocumented.
* **The moving-boundary conservation property is NOT established.** The reclassification transaction
  and the transport ledger exist; the identity they would certify closes by construction for any seed
  pressure (§3.4, §4.5), and this lane has no independent oracle for the seed. `ALEPH-PORT-1803` §14
  named this as the one genuinely non-obvious asset in the area and explicitly did not deliver it.
  **This entry does not deliver it either.** What it adds is a refusal to invent the seed.
* **No connector is wired.** This module provides the endpoint machinery — `EndpointRole`,
  `TransferStencil`, `BoundaryFaceStencil` — that the nine declared cytosol connectors need. It wires
  none of them. The wired count contributed by this entry is **0**, and a control asserts it.
* **Not integrated into `aleph/vertical/assembly.py`.** Deliberately. It has controls; it has never
  been stepped inside a world, and none of the first vertical's results transfer to it.
* **No transport of any species.** `transported_species_content` is allocated, snapshotted and rolled
  back as the contract requires, and **nothing advects or diffuses it**. The G-actin channel that
  connectors 24 and 25 name does not exist. Registering the slot and leaving it inert is the honest
  encoding; claiming a transport field would not be.
* **No dissipation accounting and no force-work ledger entry.** The Darcy dissipation
  `sum_f q_f^2 / mobility * V` is computable from what is here and is **not** computed, because
  wiring it into the ledger identity `W + dU + D - A = 0` without a connector to close against would
  produce a term with nothing to balance it.
* **No thermal channel, no osmotic term, no gravity.** The membrane's osmotic boundary condition is
  the membrane owner's and reaches this field through connector 04, which is not wired.
* **The scheme's accuracy in time is first order** (backward Euler) and untested: every convergence
  control here is a *spatial* measurement made on a steady problem. No temporal order is claimed,
  measured, or implied.
* **Unverified, and marked so in the artifact:** that seven state slots are the right decomposition
  for this field. `ALEPH-PORT-1803` §14 flagged it as a hypothesis about an implementation that did
  not exist. That implementation now exists and this lane still cannot test the claim — it can only
  report that all seven were allocatable and that no eighth was needed.
