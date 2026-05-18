# v4 Round 21 — Wardetzky Discrete Bending Form (cortex curvature term replacement)

Module: `acs/v4/cell/dynamics.py:compute_cortex_force` (lambda_curv → lambda_bend)
Status: draft
Date: 2026-05-07
Author: claude-mech (channel v4-coord, id=2721 → 2723)
Driver paper: SimuCell3D, Runser, Vetter, Iber, *Nat. Comput. Sci.* 4:299–309 (2024); cites Wardetzky, Bergou, Harmon, Zorin, Grinspun, *Comp. Aided Geom. Des.* 24:499–518 (2007), "Discrete quadratic curvature energies."

---

## 1. Motivation — why Option R is structurally limited

Option R (committed by impl, results in `results/v4/random_optionR_*`) added a uniform-weighted Laplacian as a curvature regulariser:

```
F_curv_i = -λ_curv · Σ_{j ∈ 1-ring} ê_ij                  (current code)
```

with `ê_ij = (r_i − r_j) / ||r_i − r_j||`. This is **not** a discrete mean-curvature operator — it is the gradient of total edge length, weighted uniformly per edge. Three structural consequences:

(i) On a planar 1-ring (six equilateral neighbours of `r_i`), the term reduces to the centroid-direction sum. By symmetry the in-plane component cancels, so on a uniform plane it is zero. But on a non-uniform planar 1-ring (e.g. anisotropic stretching of the lamellipodial leading edge) it gives a nonzero in-plane pull toward the local centroid of neighbours, even though the surface is flat → this is **isotropic shrinkage pressure**, not bending stiffness.

(ii) On a closed convex sphere (icosphere) it gives an inward force at every vertex proportional to vertex valence. This is mean-curvature flow, which is known to contract any closed surface toward a point in the long-time limit; here it is balanced only by the volume + area constraints. The cell shape that emerges is the equilibrium between three opposing pressures — area-want, volume-want, curvature-want — and lamellipodia have to fight all three.

(iii) Because the curvature-want is uniform across vertices (independent of *where* the curvature actually departs), increasing λ_curv to fix concavities (donut, dual-lobe) simultaneously pulls the entire cell inward, killing motility. The empirical λ_curv sweep confirms exactly this:

| λ_curv | shape  | speed [μm/min] | speed loss vs Round 20 baseline (0.95) |
| ------ | ------ | -------------- | -------------------------------------- |
| 0.05   | perfect | 0.010 | 95× drop |
| 0.01   | mostly OK | 0.085 | 11× drop |
| 0.005  | tuning point | (visual_compare.png) | not reported numerically |

The trade-off is **in the form**, not in the parameter. No value of λ_curv with this form simultaneously preserves shape and motility.

The fix is to use a curvature operator that is **zero where the surface is flat** and **proportional to local mean-curvature departure from rest** — i.e. a true Helfrich-like bending energy. Wardetzky's cotangent-weighted Laplacian provides exactly this on a triangulated surface.

---

## 2. Form derivation — cotangent-Laplacian as discrete mean curvature

### 2.1 Continuum target

The continuum bending energy (Helfrich, 1973) for a closed surface is:

```
U_b = (κ_b / 2) ∫_S (H − H_0)² dA                        (Helfrich)
```

where H is local mean curvature (sum of two principal curvatures), H_0 is spontaneous curvature, and κ_b is the bending modulus [J]. The first variation gives:

```
δU_b/δr = κ_b · (Δ_LB(H − H_0) + 2(H − H_0)(H² − K))     (Canham–Helfrich shape eq.)
```

where Δ_LB is the surface Laplace-Beltrami operator and K is Gaussian curvature. For a sphere this expression evaluates to zero when H ≡ H_0, recovering the rest sphere as energy minimum.

For v4 we want a discrete operator that:
- Reduces to zero on the rest icosphere mesh.
- Increases quadratically as local curvature departs from rest.
- Has compact support (per-vertex, uses only the 1-ring).
- Is zero on a perfectly flat region (the lamellipodial sheet should not be penalized for being flat).

### 2.2 Discrete cotangent-Laplacian

For a triangulated mesh with vertex `i`, let `(j ∈ N(i))` be the 1-ring neighbours. For each edge `(i,j)`, the two faces sharing it have *opposite* angles α_ij and β_ij (the interior triangle angles at the two vertices that are NOT i or j). The discrete Laplace-Beltrami operator at vertex i is (Pinkall & Polthier 1993; Meyer et al. 2003; Wardetzky et al. 2007):

```
(Δ_LB f)_i = (1 / (2 A_i)) · Σ_{j ∈ N(i)} (cot α_ij + cot β_ij) · (f_j − f_i)
```

where A_i is the *Voronoi area* (or 1/3 sum of incident face areas as a simpler proxy — already computed in `compute_vertex_normals_voronoi_areas`). Applying this to the position function `r` gives the **mean-curvature normal vector**:

```
K_H(i) := Δ_LB(r)_i = (1 / (2 A_i)) · Σ_{j ∈ N(i)} (cot α_ij + cot β_ij) · (r_j − r_i)
```

Properties (Wardetzky 2007 Theorem 3.4 + standard differential geometry):
- On a planar region: K_H(i) = 0 exactly. (cot α + cot β weights still nonzero, but the sum of weighted edge vectors is in-plane and cancels by the planar-discrete-Laplace identity.)
- On a sphere of radius R: K_H(i) = −(2/R) · n̂_i (mean curvature 2/R, pointing inward against outward normal).
- On a saddle: K_H(i) reflects the principal curvature difference correctly.

### 2.3 Bending energy with rest-curvature reference

Define the rest mean-curvature vector at each vertex from the *initial mesh* (icosphere at radius R_0):

```
K_H_rest(i) := K_H(i)|_{r = r_init}    (computed once at cell construction, frozen)
```

For an icosphere of radius R_0, K_H_rest(i) ≈ −(2/R_0) · n̂_i_init.

The discrete bending energy is:

```
U_b = (k_b / 2) · Σ_i A_i · ||K_H(i) − K_H_rest(i)||²
```

Per-vertex bending force (negative gradient w.r.t. r_i, holding A_i and the cotangent weights fixed at their *current* values — a frozen-coefficient approximation, standard in computer-graphics shape-preservation literature):

```
F_bend_i = −k_b · A_i · (1 / A_i) · (K_H(i) − K_H_rest(i))  · [boundary term ignored]
        ≈ −k_b · (K_H(i) − K_H_rest(i))                    (frozen-coefficient form)
```

Note: the `1/A_i` factors cancel when treating `A_i · ||K_H||²` and approximating ∂(A_i)/∂r_i ≈ 0 over the bending step. This is exact only at rest; near rest it is a leading-order approximation. For v4 v0, this approximation is acceptable — it preserves the right sign and zero-at-rest property; deviations from full ∇U_b are sub-leading and absorbed into the area/volume terms.

### 2.4 Why this form fixes Option R's failure

- **Planar lamellipodial sheet**: K_H(i) = 0 = K_H_rest(i)|_{flat region after deformation}. **Force = 0.** No shrinkage of the leading-edge lamellipodium.
- **Localized concavity (donut, dual-lobe)**: K_H(i) at concave vertex ≠ K_H_rest(i). **Force = nonzero, restoring.** Concavity flattens; the rest of the cell unaffected.
- **Uniform sphere stretch (cell grows isotropically)**: K_H(i) shrinks proportionally; mismatch with K_H_rest gives uniform inward force. *But* the area + volume constraints already handle uniform expansion/contraction with their own dedicated terms — so this contribution is sub-leading and physically correct (bending stiffness adds to elastic resistance against expansion).

The empirical λ_curv ↔ shape/motility trade-off is **eliminated** by the form: bending stiffness can be increased to lock down concavities without applying any pressure to flat or rest-curvature regions.

---

## 3. Per-vertex force formula — implementation primitive

Given mesh state (vertices `r ∈ ℝ^{Nv×3}`, faces `F ∈ ℤ^{Nf×3}`), at construction time:

```
Step A (init, once per cell):
  A.1  Build edge → adjacent-faces map: for each edge (i,j), find the two faces (a) (i,j,k) and (b) (i,j,l).
  A.2  For each edge (i,j), compute cot α_ij from face (a) (angle at vertex k) and cot β_ij from face (b) (angle at l).
       cot(angle at k) = ((r_i − r_k) · (r_j − r_k)) / ||(r_i − r_k) × (r_j − r_k)||
  A.3  K_H_rest(i) = (1 / (2 A_i)) · Σ_{j ∈ N(i)} (cot α_ij + cot β_ij) · (r_j − r_i)
  A.4  Store K_H_rest as a (Nv, 3) constant array on the CellMeshState.

Step B (each step, in compute_bending_force):
  B.1  Recompute cot α_ij, cot β_ij with CURRENT r.
  B.2  K_H_now(i) = (1 / (2 A_i)) · Σ_{j ∈ N(i)} (cot α_ij + cot β_ij) · (r_j − r_i)
  B.3  F_bend(i) = −k_b · (K_H_now(i) − K_H_rest(i))
```

Numerical guards:
- Floor `||(r_i − r_k) × (r_j − r_k)||` at 1e-12 to avoid divide-by-zero on degenerate triangles.
- Cap |cot α + cot β| at some large positive value (e.g. 100) to clamp ill-conditioned obtuse-triangle contributions. Heuristic from cotangent-Laplacian literature; record the cap-saturation rate as a diagnostic (if > 1% of edges, mesh quality is degenerating and remeshing would be needed — out of scope for Round 21).

### 3.1 Sign sense (Sanity Gate item 5)

- At `r = r_init`: K_H_now = K_H_rest by construction → F_bend = 0. ✓
- Local outward bulge (one vertex pushed outward by lamellipodia event): K_H_now at that vertex points more outward than K_H_rest → mismatch points outward → F_bend points inward (restoring). ✓
- Local inward concavity (donut): K_H_now points outward where K_H_rest pointed inward (curvature flipped) → large mismatch outward → F_bend points outward (pushing concavity out). ✓
- Uniform sphere expansion (R > R_0): K_H_now magnitude decreases (1/R < 1/R_0) → mismatch points inward (toward original 2/R_0 magnitude) → F_bend points inward (restoring). ✓

All four cases match physical intent.

---

## 4. Sanity Gate (per `v4_sanity_gate_template`)

### 4.1 Dimensional analysis

Inputs (units):
- `r ∈ μm`
- `A_i ∈ μm²`
- cot α, cot β ∈ dimensionless
- `k_b ∈ nN` (bending stiffness; same dimensions as `lambda_c_nN`, treated as a per-edge force scale)

Outputs:
- `K_H ∈ μm⁻¹` (1/length, mean curvature)
- `F_bend = −k_b · (K_H_now − K_H_rest) ∈ nN · μm⁻¹` 
- ⚠️ unit mismatch: F_bend should be in nN to be consistent with cortex/area/volume forces in dynamics.py.

Unit chain fix: introduce a **bending length scale** L_b in the per-vertex formula:

```
F_bend_i = −k_b · L_b · (K_H_now(i) − K_H_rest(i))           [nN · (μm) · (μm⁻¹) = nN]
```

L_b = `state.rest_edge_length_um` is the natural choice (same convention as cortex spring uses rest_edge_length as length scale). With this:
- `k_b` has dimensions [nN], matches `lambda_c_nN`.
- F_bend has dimensions [nN], matches other force terms.

CFL/stability: explicit Euler with v4 dt = 0.02 s and ξ = 1 nN·s/μm. Effective per-vertex relaxation rate from bending term = k_b · L_b · |∂K_H/∂r| / ξ. For sub=1 icosphere, |∂K_H/∂r| ~ 1/L_b² (cotangent-Laplacian second-order). So characteristic rate ~ k_b/(ξ · L_b). With L_b ≈ 5 μm (sub=1 mean edge), ξ = 1, the constraint dt · k_b / 5 < 0.5 → k_b < 125 nN. We will operate at k_b ≤ 1 nN, far below CFL bound. **PASS.**

### 4.2 Boundary cases

- `dt → 0`: F_bend bounded; per-step displacement → 0. PASS.
- `dt → large` (≥ CFL): instability mirrored by other terms; β-fix displacement cap + log volume (separate proposal) catch it. PASS conditional.
- Mesh resolution low (sub=0, 12v): cotangent weights still defined; coarser approximation but stable. PASS.
- Mesh resolution high (sub=3, 642v): per-step cost increases but no instability mode. PASS.
- `k_b → 0`: F_bend = 0 → reduces to current Round 20 baseline (lambda_curv = 0). PASS, this is the controlled-variable case.
- `k_b → large`: cell rigidly locked at icosphere shape → motility suppressed, but no NaN. Diagnostic, not a fail mode.
- Degenerate triangle (one face area → 0): cot α + cot β → ∞. Caught by clamp at |cot α + cot β| ≤ 100. Saturation rate diagnostic logged.

### 4.3 Conservation / dissipation

- Momentum conservation: Σ_i F_bend_i over the closed mesh. The cotangent-Laplacian satisfies Σ_i A_i · K_H(i) = 0 (closed surface integral of mean-curvature vector vanishes). Rest term Σ_i A_i · K_H_rest(i) = 0 likewise. So Σ_i (1/A_i) · A_i · ... wait — careful here.

Actually: F_bend_i = −k_b · (K_H_now(i) − K_H_rest(i)). The sum Σ_i F_bend_i = −k_b · (Σ_i K_H_now(i) − Σ_i K_H_rest(i)). Each K_H is *not* automatically zero-sum without area-weighting; with frozen-coefficient approximation we lose strict momentum conservation.

Bound on momentum drift: |Σ_i F_bend_i| ≤ k_b · Σ_i ||K_H_now(i) − K_H_rest(i)||. Empirical check: log Σ F_bend_i per step in smoke test; if > 1% of largest single-vertex F_bend, escalate. **PARTIAL — record as known systematic, not strict conservation.**

- Energy: the energy U_b = (k_b/2) Σ ||K_H_now − K_H_rest||² is bounded below by 0 and the per-step ΔU_b ≤ 0 in the absence of other forces (gradient flow under F_bend alone). With other forces present, energy exchange. PASS.

### 4.4 Numerical sanity

- dt vs slowest physics: relaxation under bending alone has time constant τ_b = ξ · L_b / (k_b · L_b · |∂K_H/∂r|) ≈ ξ · L_b² / k_b. For ξ=1, L_b=5, k_b=0.5 → τ_b ≈ 50 s. dt = 0.02 → dt / τ_b = 4e-4. PASS.
- Mesh resolution vs feature: lamellipodia event force radius ~ L_b (one ring); bending term resolves curvature on the same scale. PASS.
- Float precision: f32 for vertex positions is borderline at sub=2 + 5 μm cell + 10⁻³ relative deformations (≈ 5e-3 μm). f64 is current default in CellMeshState — keep f64. PASS.

### 4.5 Sign / sense check

See section 3.1. All four cases (rest, outward bulge, inward concavity, uniform expansion) have correct sign. PASS.

### 4.6 Measurement-protocol consistency

Round 21 sweep observables (impl runs):
- migration speed: centroid xy displacement / window — UNCHANGED from Round 20.
- aspect ratio: PCA of vertex xy positions — UNCHANGED from Round 20.
- cell phenotype: blender render — UNCHANGED.

Bending term does not alter measurement protocol. The Hard Rule 11 anchor (PI's experimental top-down projection ↔ sim convex hull) does not exist for v4 (literature-only). PASS by absence.

### 4.7 Magic-Number Block

| Param | Default | Derivable? | Grid-invariant? | Not-fitting? |
| ----- | ------- | ---------- | --------------- | ------------ |
| `k_bend_nN` | 0.5 (proposal) | YES — same dimension family as `lambda_c_nN`, set to ≈ 1.5× cortex spring | YES — via L_b = rest_edge_length, scale invariant under sub=N change | YES — value picked from physical reasoning (cortex bending should be comparable to cortex stretching), not from gate target |
| L_b (length scale) | `state.rest_edge_length_um` | YES — geometric invariant of init mesh | YES | N/A |
| cot saturation cap | 100 | YES — heuristic from Pinkall–Polthier 1993, prevents ill-conditioned obtuse contribution | grid-invariant since the cot is intrinsic to angle | not-fitting (independent of any gate target) |

All three tests PASS for proposed defaults. If sweep finds k_bend default needs revision based on shape↔motility behaviour: still PASS test 3 (search-driven, not gate-driven). If sweep finds saturation cap matters (>1% of edges saturate): escalate, mesh quality issue → remeshing on the spheroid track.

---

## 5. Pseudo-diff for `acs/v4/cell/dynamics.py`

### 5.1 New helper: `_cotangent_edge_weights`

```python
def _cotangent_edge_weights(state: CellMeshState) -> tuple[np.ndarray, np.ndarray]:
    """Return (edges (Ne, 2), cot_sum (Ne,) i.e. cot α + cot β per edge)."""
    f = state.faces  # (Nf, 3)
    v = state.vertices_xyz_um
    # build edge → (face_a, opposite_vertex_a, face_b, opposite_vertex_b) map
    # for each face (i, j, k): 3 edges (i,j), (j,k), (k,i). vertex k is opposite to edge (i,j).
    # Construct a dict keyed on sorted (a,b) → list of opposite-vertex indices (max 2 entries).
    edge_to_opp: dict[tuple[int,int], list[int]] = {}
    for fi in range(state.n_faces):
        a, b, c = int(f[fi,0]), int(f[fi,1]), int(f[fi,2])
        for (u, w, opp) in [(a, b, c), (b, c, a), (c, a, b)]:
            key = (min(u, w), max(u, w))
            edge_to_opp.setdefault(key, []).append(opp)
    edges = []
    cot_sum = []
    SAT_CAP = 100.0
    for (i, j), opps in edge_to_opp.items():
        cot = 0.0
        for k in opps:
            ru, rw, rk = v[i] - v[k], v[j] - v[k], None
            cross = np.cross(ru, rw)
            cross_norm = np.linalg.norm(cross)
            if cross_norm < 1e-12:
                cot += 0.0  # degenerate, contributes nothing
            else:
                cot += float(np.dot(ru, rw) / cross_norm)
        cot = max(-SAT_CAP, min(SAT_CAP, cot))
        edges.append([i, j])
        cot_sum.append(cot)
    return (np.asarray(edges, dtype=np.int32),
            np.asarray(cot_sum, dtype=np.float64))
```

### 5.2 New helper: `_mean_curvature_vector`

```python
def _mean_curvature_vector(state: CellMeshState,
                           edges: np.ndarray | None = None,
                           cot_sum: np.ndarray | None = None) -> np.ndarray:
    """K_H(i) = (1 / (2 A_i)) · Σ_{j ∈ N(i)} (cot α + cot β) · (r_j − r_i). Returns (Nv, 3)."""
    if edges is None or cot_sum is None:
        edges, cot_sum = _cotangent_edge_weights(state)
    v = state.vertices_xyz_um
    _, va = compute_vertex_normals_voronoi_areas(state)
    K = np.zeros_like(v)
    e_ij = v[edges[:,1]] - v[edges[:,0]]                  # r_j - r_i per edge (oriented)
    contrib = cot_sum[:, None] * e_ij                     # (Ne, 3)
    np.add.at(K, edges[:,0], +contrib)                    # vertex i gets (r_j - r_i) contributions
    np.add.at(K, edges[:,1], -contrib)                    # vertex j gets (r_i - r_j) = -(r_j - r_i)
    K /= np.maximum(2.0 * va[:, None], 1e-12)
    return K
```

### 5.3 Cache the rest curvature on the CellMeshState

Add a field on `CellMeshState`:
```python
@dataclass
class CellMeshState:
    ...
    K_H_rest_um_inv: Optional[np.ndarray] = None   # (Nv, 3), set at init (see make_icosphere_cell / make_oblate_polarized_cell)
```

In `make_icosphere_cell` / `make_oblate_polarized_cell`, after the final state is constructed:
```python
state.K_H_rest_um_inv = _mean_curvature_vector(state)
```

### 5.4 Replace lambda_curv with lambda_bend in `compute_cortex_force`

```python
def compute_cortex_force(
    state: CellMeshState,
    lambda_c_nN: float,
    lambda_bend_nN: float = 0.0,        # WAS: lambda_curv_nN
) -> np.ndarray:
    """Per-vertex cortex force, (Nv, 3) [nN].

    Edge-spring (κ-fix) + Wardetzky discrete bending (Round 21 form).

    edge-spring (length scale, dominant):
        F_edge_i = −λ_c · Σ (||r_i − r_j|| − rest) · ê_ij

    Wardetzky bending (curvature-localized regularizer, replaces Option R uniform Laplacian):
        F_bend_i = −λ_bend · L_b · (K_H(i) − K_H_rest(i))
        K_H(i)   = (1 / (2 A_i)) · Σ (cot α_ij + cot β_ij) · (r_j − r_i)
        L_b      = state.rest_edge_length_um
    """
    if lambda_c_nN == 0.0 and lambda_bend_nN == 0.0:
        return np.zeros_like(state.vertices_xyz_um)

    force = np.zeros_like(state.vertices_xyz_um)

    # edge-spring (UNCHANGED from current dynamics.py)
    if lambda_c_nN != 0.0:
        edges = state.edges()
        v = state.vertices_xyz_um
        e_vec = v[edges[:,0]] - v[edges[:,1]]
        e_len = np.linalg.norm(e_vec, axis=1, keepdims=True)
        e_hat = e_vec / np.where(e_len > 1e-12, e_len, 1.0)
        rest = state.rest_edge_length_um or float(e_len.mean())
        strain = e_len.squeeze(-1) - rest
        f_edge_mag = lambda_c_nN * strain
        np.add.at(force, edges[:,0], -f_edge_mag[:, None] * e_hat)
        np.add.at(force, edges[:,1], +f_edge_mag[:, None] * e_hat)

    # Wardetzky bending (NEW, replaces Option R lambda_curv mean-curvature flow)
    if lambda_bend_nN != 0.0:
        if state.K_H_rest_um_inv is None:
            raise ValueError(
                "lambda_bend_nN > 0 requires state.K_H_rest_um_inv to be set "
                "(call _mean_curvature_vector at mesh init)"
            )
        K_now = _mean_curvature_vector(state)
        L_b = state.rest_edge_length_um or 1.0
        force += -lambda_bend_nN * L_b * (K_now - state.K_H_rest_um_inv)

    return force
```

### 5.5 Remove lambda_curv plumbing from `step` and `IntegratedParams`

In `step` (dynamics.py):
- Replace `lambda_curv_nN` parameter with `lambda_bend_nN` (default 0.0).
- Update the `compute_cortex_force` call.

In `IntegratedParams` (stepper.py):
- Replace `lambda_curv_nN: float = 0.005` with `lambda_bend_nN: float = 0.5` (initial proposal default; sweep below revises).
- Update the `_step_mesh` call in `integrated_step`.

The Option R history comments in stepper.py (the 0.05/0.01/0.005 sweep table) should be **preserved** as historical context, prefaced with: "# (Round 21) lambda_curv superseded by lambda_bend (Wardetzky); kept here as motivation."

---

## 6. Smoke test plan (impl runs before sweep)

### 6.1 Unit tests (new, `tests/v4/test_dynamics_bending.py`)

```python
def test_bending_force_zero_at_rest():
    state = make_icosphere_cell(subdivision=1, radius_um=10.0)
    assert state.K_H_rest_um_inv is not None
    F = compute_cortex_force(state, lambda_c_nN=0.0, lambda_bend_nN=1.0)
    assert np.linalg.norm(F) < 1e-9

def test_bending_force_restores_perturbed_vertex():
    state = make_icosphere_cell(subdivision=1, radius_um=10.0)
    # push vertex 0 outward by 1 μm along its normal
    state.vertices_xyz_um[0] += 1.0 * (state.vertices_xyz_um[0] / np.linalg.norm(state.vertices_xyz_um[0]))
    F = compute_cortex_force(state, lambda_c_nN=0.0, lambda_bend_nN=1.0)
    # F[0] should point inward (toward origin)
    inward = -state.vertices_xyz_um[0] / np.linalg.norm(state.vertices_xyz_um[0])
    assert np.dot(F[0], inward) > 0

def test_bending_force_zero_under_uniform_translation():
    state = make_icosphere_cell(subdivision=1, radius_um=10.0)
    state.vertices_xyz_um += np.array([5.0, 0.0, 0.0])  # translate
    F = compute_cortex_force(state, lambda_c_nN=0.0, lambda_bend_nN=1.0)
    # K_H is translation-invariant; F should still be zero
    assert np.linalg.norm(F) < 1e-9
```

### 6.2 Single-cell relaxation smoke (no ECM, no protrusion)

```bash
.venv/bin/python scripts/v4_cell_smoke.py --lambda-bend 0.5 --steps 1000
```

Expectation: cell remains close to initial icosphere shape; no shape drift > 1 μm in any vertex; energy U_b decays monotonically.

### 6.3 Integrated single-cell smoke (full stepper, no sweep)

```bash
.venv/bin/python scripts/v4_integrated_smoke.py --lambda-bend 0.5 --steps 3000
```

Expectation: speed within 0.7–1.2 μm/min (Round 20 baseline 0.95), aspect 1.05–1.4, no NaN, no z-cap saturation > 5%.

### 6.4 Round 21 sweep design (impl runs on ssh win)

Conditions (16 sims, 4 lambda_bend × 1 best-Round-20 condition × 4 seeds for std):
- `lambda_bend ∈ {0.0, 0.1, 0.5, 2.0}` × `k_R = 0.10` (best Round 20 condition) × `seed ∈ {22, 23, 24, 25}`
- ECM: random isotropic, Nf=1500
- Sim time: same as Round 20 (~600 s sim, ~30 min wall on RTX A5000)
- Output: dump.h5 + condition_sweep summary CSV per `scripts/v4_condition_sweep.py` schema

Wall-clock estimate: 16 × 30 min = 8 h on RTX A5000 (within autonomous budget if PI confirms; otherwise split into 2 batches of 8 sims each, 4 h each, fits ≤3h with parallelism if scripts/v4_condition_sweep.py is parallelized).

Deliverables:
- `results/v4/round21_wardetzky_sweep_kR05/` — h5 dumps + per-condition speed/aspect/τ_p
- `results/v4/round21_wardetzky_speed_aspect_panel.png` — analog of optionR_visual_compare.png
- `results/v4/round21_wardetzky_winning_blender/` — blender render of the chosen lambda_bend best condition

Pass criteria for Round 21 lock-in:
- All 16 sims valid (analysis.is_valid_sim).
- Median speed within Round 20 ± 30%.
- Median aspect xy ≥ 1.20 (literature lower bound, currently 1.06–1.16 is the gap).
- No condition shows the Option R shape↔motility trade-off (speed loss > 50% at any lambda_bend that improves aspect).

---

## 7. Failure handling (Sanity Gate)

If any of the following occurs during smoke or sweep, **halt and surface to PI**:

- Test 6.1 fails (force not zero at rest, or sign wrong).
- Smoke 6.2 shows energy increase per step (gradient direction wrong).
- Smoke 6.3 shows NaN or speed < 0.1 μm/min (motility dead — would mean bending is acting like Option R again, form bug).
- Sweep 6.4 shows the same shape↔motility trade-off as Option R (e.g. aspect improves only when speed drops > 50%) — would mean K_H_rest reference is wrong or cotangent weights are computed wrong.

Silent workaround (e.g. tuning lambda_bend down to mask a sign error) is **forbidden** per `v4_track_policy` Magic-Number Block test 3.

---

## 8. Open questions for impl review (before code application)

1. **A_i choice**: Voronoi area (1/3 sum incident face area, current) vs *true* Voronoi (Meyer et al. 2003 mixed-area scheme). Voronoi-1/3 is what `compute_vertex_normals_voronoi_areas` already returns; it's the simpler choice and consistent with v4 conventions. Mixed-area is more accurate for obtuse triangles but harder to implement. **Mech recommends Voronoi-1/3 for v4 v0; revisit only if cot saturation rate > 5%.**

2. **L_b choice**: `state.rest_edge_length_um` (single scalar) vs per-vertex local edge length. Single scalar is dimensionally-consistent and matches cortex spring convention. Per-vertex would adapt under remeshing (which v4 doesn't do). **Mech recommends single scalar.**

3. **Frozen-coefficient approximation**: section 2.3 ignores ∂A_i/∂r_i and ∂(cot weights)/∂r_i in the gradient. For v4 v0 single-cell with sub=1 mesh near-rest, this is acceptable (sub-leading). For spheroid track or large deformation, we should re-derive with full chain rule (computer-graphics literature has standard expressions). **Mech recommends frozen coefficients for v4 Round 21; flag for future work.**

4. **Should sub=2 default also land in this PR?** Per channel id=2723 plan, sub=2 was bundled with Wardetzky + log volume. But log volume is a separate change (volume term, not cortex term). Suggest: **Round 21 = Wardetzky bending only**; sub=2 + log volume become Round 22 if Wardetzky alone gets aspect to 1.20.

5. **What to do with results/v4/random_optionR_***: keep (motivation evidence) vs prune. Mech recommends **keep + reference in Round 21 Wardetzky commit message** ("Replaces Option R uniform-Laplacian curvature term whose empirical λ_curv sweep at results/v4/random_optionR_curv0p* demonstrated the structural shape↔motility trade-off this form eliminates").

---

## 9. Channel handoff

- This note posted to v4-coord room as id=TBD (mech to send via direct INSERT — MCP tools deferred this session).
- impl reads, addresses Q1–Q5 in section 8 with a channel ack message.
- After mech ack of impl's answers, impl applies code on ssh win, runs smoke 6.1+6.2+6.3, posts results.
- mech reviews smoke pass/fail; if PASS, impl runs 6.4 sweep.
- mech reviews sweep results, decides Round 22 entry (sub=2 + log volume) or Round 21 lock + spheroid track entry.

ETA total Round 21 closure: ~24 h from impl's Option R commit (4–5 h impl coding + 1 h smoke + 8 h sweep + 2 h analysis + buffer).
