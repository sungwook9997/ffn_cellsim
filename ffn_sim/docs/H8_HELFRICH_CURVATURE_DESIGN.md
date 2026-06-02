# H.8 Helfrich bending κ_m on a particle shell — discrete-curvature DESIGN

> **STATUS: DESIGN doc — 2026-06-02. Prep, NOT a contract, NO code.** Parallel session
> (`AUTONOMOUS_LOG_2026-06-02.md`). Closes the one documented TODO in the *shipped*
> H.8 module `cell/membrane_surface.py` (lines 61–73): the Helfrich bending term κ_m is
> deliberately deferred there because mean curvature on an unstructured particle shell
> is non-trivial. This specifies the discrete-curvature operator + the staged build, so
> a later (PI-sequenced) refinement starts from a vetted plan. No `integrator/` touch.

## What exists vs the gap (code-grounded)

`cell/membrane_surface.py` implements the in-plane membrane mechanics as a Young-Laplace
shell force: bare tension `γ_mem` + area-elasticity `K_A·(A−A0)/A0`, both inward
`ΔP = 2γ_tot/R`. It **defers** the Helfrich curvature-elastic term

    U_bend = ½ κ_m ∮ (2H − c0)² dA                                          [J]

because the per-bead **mean curvature H** must be estimated on a point-cloud shell with
**no fixed mesh connectivity** at the ×40 coarse scale. Today κ_m enters NO runtime
force; it is only (a) a literature constant and (b) the **tether-force oracle**
`f_t = 2π√(2κ_m(T_m + γ_MCA))` used TEST-ONLY to check that κ_m / T land f_t in the
KU-3.B1 band. This doc is the design to turn κ_m into a real (default-off) force.

**KU anchors (KU-3.B1):** κ_m ≈ 1×10⁻¹⁹ J (10–30 k_BT); tether f_t ∈ 5–40 pN.

## The hard part — discrete H on an unstructured shell

Options for a mesh-free per-bead mean-curvature estimate (cortex shell tags
`[0, n_cortex_actin)`, centroid `c`, outward normal estimate `n̂_i`):

| Operator | Idea | Fit for a *dynamic* particle shell | Cost |
|---|---|---|---|
| **Local quadric fit** (recommended for measurement) | per bead: PCA normal from k-NN → fit `w = a u² + b uv + c v² + d u + e v` in the tangent frame → `H = (a(1+e²) − bde + c(1+d²))/(1+d²+e²)^{3/2}` (≈ a+c for small slope) | **good** — mesh-free, robust, k-NN only | O(N·k) |
| **Laplace-Beltrami of the embedding** | `Δ_LB r = 2H n̂` (mean-curvature normal); discrete via a Gaussian/heat-kernel graph Laplacian (Belkin-Niyogi) | **good** — mesh-free; gives the curvature *normal* directly (force-ready) | O(N·k) |
| Cotangent Laplacian (Meyer 2003) | exact on a triangle mesh | **poor** — needs a Delaunay re-triangulation each step on a moving shell | O(N) + mesh maint. |
| Osculating-sphere fit | local sphere radius → H = 1/R_local | cheap but **noisy** | O(N·k) |

The **Laplace-Beltrami** route is the most mechanistically natural for the *force*: the
mean-curvature normal `Δ_LB r = 2H n̂` is exactly what a bending energy differentiates
into. The **local quadric fit** is the most robust for the *measurement* (per-bead H,
principal curvatures) used to validate κ_m against the tether oracle.

## Staged build (default-off, mirroring the existing H.8 shell-force shape)

**Tier-A — measurement only (no runtime force; pure validation).** Add a per-bead H
estimator (local quadric fit) as an analysis function. Validate against the analytic
sphere (below) and feed the tether-force oracle. This alone lets H.8 *report* the
bending state without changing any force — zero runtime risk, immediate value.

**Tier-B — leading-order bending force (default-off).** Add the Helfrich force as a
Template-1 custom force over the shell tags (same shape as the existing
`MembraneSurfaceTension`): `F_i = −δU_bend/δr_i`. For a first mechanistic version use
the Laplace-Beltrami curvature normal to build a **bending-stiffness force that
penalizes local curvature deviation from `c0`** (a discrete Willmore/Canham-Helfrich
leading term). Flag clearly that the FULL Helfrich force is high-order (it involves the
surface Laplacian of curvature, `Δ_LB(2H) + 2H(2H²−2K)`) and is a research-grade
discrete-differential-geometry task — Tier-B ships the leading bending resistance, not
the complete Willmore flow. Default-off ⇒ pre-Tier-B builds bit-for-bit identical.

## Sanity Gate (recorded now per CLAUDE.md)

1. **Dimensional.** `H` [1/m]; `κ_m` [J]; `U_bend = ½κ_m∮(2H−c0)²dA` [J·(1/m)²·m² = J] ✓;
   `F = −∇U` [N] ✓.
2. **Boundary — the analytic sphere test (the key one).** A perfect sphere of radius R
   has uniform `H = 1/R`, so (c0=0) `U_bend = ½κ_m·(2/R)²·4πR² = 8π κ_m`, **independent
   of R**. Any discrete H operator MUST reproduce (i) per-bead `H ≈ 1/R` on a sphere
   shell within tolerance, and (ii) `∮(2H)²dA → 16π` (⇒ U_bend → 8πκ_m). This is the
   single decisive validation of the discrete operator. (Add Gaussian κ_G later — its
   integral is the topological 4πκ_G, constant for a sphere.)
3. **Conservation.** A uniform sphere ⇒ zero *net* bending force (uniform curvature, no
   gradient) — only non-spherical perturbations feel a restoring force. Check ΣF ≈ 0 on
   a sphere (like the nucleus net-force check, AUTONOMOUS_LOG 2026-05-31 iter 8).
4. **Numerical / CFL.** κ_m bending is a stiff term → the Template-1 `τ = γ_b/k_eff`
   stiffness gate applies; precedent k_ERM 1000× softening (PI-gate, 2026-05-26). The
   effective bending stiffness `k_bend ~ κ_m/⟨ℓ⟩²` over the bead spacing ⟨ℓ⟩ — surface
   to PI if it forces `dt` down; do not soften silently.
5. **Sign-sense.** Positive κ_m resists curvature change: a dimple (locally higher H)
   feels a force flattening it toward c0. Wrong sign = a runaway crumple.
6. **Measurement consistency.** Report H at a stated k-NN count; the operator must be
   grid-/density-invariant (no magic neighbor count chosen to pass the sphere test).

## Validation acceptance (candidate — extends VG-H8)

| Gate | Criterion | KU |
|---|---|---|
| Discrete-H sphere test | per-bead H ≈ 1/R; ∮(2H)²dA → 16π (U_bend → 8πκ_m) | (operator correctness) |
| Tether force | f_t = 2π√(2κ_m(T_m+γ_MCA)) ∈ 5–40 pN with κ_m ≈ 1e-19 J | KU-3.B1.4 |
| Bending modulus band | κ_m ∈ 10–30 k_BT | KU-3.B1 |

## Open items for PI / Lead

- [ ] Ratify the staging (Tier-A measurement-only first; Tier-B leading-order force).
- [ ] Pick the operator: local quadric fit (measurement) + Laplace-Beltrami (force).
- [ ] Confirm Tier-B ships the *leading* bending term (not full Willmore) — scope the
      complete Helfrich force as a separate research item if needed.
- [ ] CFL stiffness gate decision for `k_bend` (softening precedent k_ERM).
- [ ] Sequence after the membrane composite-tension re-validation (KU-3.5) lands.
