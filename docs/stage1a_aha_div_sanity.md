# AHA 2010 §4 reproducing divergence — verification before code

Before modifying `_build_curvature` to use the Adami-Hu-Adams 2010 §4
reproducing divergence operator, this document records the three verification
items the PI requested, plus the resulting design choice.

The paper:
> Adami, S., Hu, X. Y., Adams, N. A. (2010), "A new surface-tension
> formulation for multi-phase SPH using a reproducing divergence
> approximation", *J. Comp. Phys.* **229**, 5011–5021.

---

## 1. What does AHA §4 actually do?

§4 of the paper introduces a "reproducing kernel" SPH divergence that is
exact for first-order polynomial vector fields (Bonet-Lok 1999 lineage).
For a vector field A at particle i:

```
(∇·A)_i = L_i^{−1} : Σ_j (A_j − A_i) ⊗ ∇W_ij · V_j
```

where the renormalisation tensor is

```
L_i = Σ_j (x_j − x_i) ⊗ ∇W_ij · V_j
```

The point of L_i^{−1}: standard SPH divergence has truncation errors that
scale with the (finite) kernel-density variance. Multiplying by L_i^{−1}
absorbs that variance, restoring linear-reproducing accuracy.

## 2. Does this translate literally to our MPM-grid setup?

Our solver does NOT use SPH; it uses an MPM-Eulerian background grid with
quadratic B-spline kernels. The "divergence" we currently use is

```
∇·n̂[I] = (n̂[i+1,j,k][0] − n̂[i−1,j,k][0]) / (2·dx) + … (y, z)
```

which is **already** first-order linear-reproducing on a uniform Cartesian
grid (central FD has zero truncation error for linear fields).

If we apply the AHA §4 formula treating each grid cell as a "particle" with
volume dx³ and the unit-cell connectivity as the SPH kernel:

```
L_I = Σ_J∈neigh (x_J − x_I) ⊗ ∇W_IJ · dx³
    = (dx · I) · (constants from B-spline)
    = constant · I
```

so L_I^{−1} is a constant matrix, and the reproducing-divergence formula
collapses back to a constant-scaled central difference — **no improvement
over what we already do**. This is consistent with the well-known fact that
Bonet-Lok renormalisation matters for irregular SPH point distributions but
is identity-equivalent on a uniform grid.

**Conclusion of verification 1+2:** Literal AHA §4 on our grid gives no
gain. The 40% κ error in v12 is NOT from divergence-operator order; it is
from the n̂ field being ill-defined where |∇c| → 0 inside the smoothed band.

## 3. AHA-inspired alternative consistent with the paper's spirit

The deeper purpose of AHA §4 is a *robust* curvature near the free surface.
The grid-native equivalent is a long-standing CSF refinement
(Brackbill 1992 §V, Sussman-Smereka-Osher 1994, Kang-Fedkiw-Liu 2000) that
exploits the identity

```
κ = −∇·n̂ = −∇·(∇c / |∇c|) = − Δc / |∇c| + (∇c · ∇|∇c|) / |∇c|²
```

For a sphere of radius R with a smooth, radially-symmetric diffuse colour
c(r), the second term vanishes at the gradient peak and the leading term
recovers κ = 2/R *exactly* (analytical: c″(r)=0 at the peak, Δc = (2/r)c′(r),
|∇c| = |c′(r)|, κ = −Δc/|∇c| = +2/R).

This **Laplacian-form curvature** is the grid-native AHA §4 analogue:

```
κ_I = − Δc[I] / max(|∇c[I]|, ε_div)
```

where Δc is the standard 7-point Laplacian and ε_div is a numerical-safety
floor (preventing 0/0 in pure-vacuum cells; not a tuning parameter — passes
all three Magic-Number Block tests).

## Sanity-Gate checks for the new operator

### 1. Dimensional analysis
- Δc has units of [c]/[L]² (= 1/[L]² since c is dimensionless after AHA §3).
- |∇c| has units of 1/[L].
- κ = −Δc/|∇c| has units of 1/[L]. ✓ (matches geometric curvature 2/R)

### 2. Boundary cases
- **Bulk (after AHA §3 normalisation):** c ≡ 1, ∇c ≡ 0, Δc ≡ 0. κ = 0/ε = 0. ✓
- **Vacuum:** c = 0, ∇c = 0, Δc = 0. κ = 0. ✓
- **Surface band:** c varies smoothly from 1 to 0, |∇c| > 0, Δc captures
  the curvature signature. ✓
- **Single-grid-cell c-jump (n_smoothing_passes=0):** Δc has a singular
  discrete value, |∇c| is large, κ ≈ Δc/|∇c| ≈ ±1/dx (grid-scale). The
  smoothing passes (kept) widen the interface to ≈ 4·dx, recovering κ ~ 2/R.

### 3. Conservation invariants
The operator only computes a diagnostic / impulse coefficient field κ;
it does not change conservation properties of the solver. CSF impulse
itself (`dv = γ·κ·∇c·dt/ρ`) is structurally identical to v12 — only the
*value* of κ at each grid cell changes.

### 4. Numerical sanity
- 7-point Laplacian: standard 2nd-order central FD on a uniform grid,
  truncation error O(dx² · c″″) for smooth c.
- For a sphere: at the gradient peak we have c″ = 0 and the error is
  O(dx² · c″″) ≈ O((dx/R)² · 1/R) ≈ 1% at our resolution (dx=0.094, R=1).
- f32 precision: each Δc term is one subtraction-of-near-equals (c[i+1] +
  c[i−1] − 2·c[i]); cancellation is at the |∇c|·dx ≈ 0.3 scale, well
  inside f32 range when c is O(1). Justified.

### 5. Sign / sense check
- For a sphere c(r) decreasing from 1 to 0 across the surface:
  ∇c points inward, |∇c| > 0.
  Δc = (2/r) · c′(r) at the peak; c′ < 0 ⇒ Δc < 0.
  κ = −Δc/|∇c| = −(negative)/positive = +positive = +2/R. ✓
  The CSF impulse `dv = γ·κ·∇c·dt/ρ` then has γ > 0, κ > 0, ∇c inward
  ⇒ inward impulse, physically correct restoring surface tension. ✓

## Cost estimate

- Drop the n̂ field computation (one 4-vec field × n_g³ writes, eliminated).
- Add Δc field (one f32 field × n_g³ writes, computed via the existing
  smoothed colour). Net field memory: change of one Vec3 → one f32, **net
  saving** 8 bytes/cell ≈ 2 MB at n_g=64.
- One additional grid pass for Δc; but we save the n̂ pass. **Net step time
  ≈ unchanged** (predicted ±5%).

## Step-ordering consistency

Replace `_build_curvature` only; everything else is unchanged:

```
_clear_grid → _tag_boundary → _p2g
            → _build_csf_field   (unchanged: AHA-§3 colour seed + smoothing + ∇c)
            → _build_curvature   (AHA-§4 form: κ = −Δc / max(|∇c|, ε_div))
            → _grid_op_overdamped (unchanged: dv = γ·κ·∇c·dt/ρ)
            → _g2p_and_constitutive
```

The new kernel reads `grid_color` (smoothed) and `grid_color_grad` and
writes only `grid_kappa`. The `grid_normal` field is no longer used; it
will be retained as an optional diagnostic (cleared but not consumed by the
CSF impulse).

## What this verifies, in one sentence

The literal AHA §4 SPH formula collapses to central FD on a uniform grid
and would not change the v12 result; the grid-native reproducing-divergence
analogue (`κ = −Δc/|∇c|`, citing Brackbill 1992 §V / Sussman-Smereka-Osher
1994 / Kang-Fedkiw-Liu 2000 alongside AHA 2010 §4) is the principled
alternative — passes the Sanity-Gate checks above, costs ≈ 0 extra, and is
expected to drive κ from ≈ 2.80 toward 2.0 for our sphere.
