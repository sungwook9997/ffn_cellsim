# V2 Phase 1 P1 Active Contour — Locked Design Derivation

**Date**: 2026-05-04 KST
**Authors**: Claude + Codex design-discussion (5-round adversarial lock,
PI id=809 aggressive debate posture applied)
**Source unit**: design-discussion `topic=v2-layer-2-p1-derivation`,
MCP id 943-964
**PI ratify status**: full delegation per PI id=939 (overnight automation).
This lock is the design-team source of truth for P1 active contour
implementation. impl-work uses this for the Sanity Gate doc rewrite.

---

## 0. Scope

This document locks the six design questions that gated P1 active contour
executable code (per consolidated plan §13 P1 conditional). Lock is
sufficient for impl-work to rewrite the P1 Sanity Gate doc from scratch
and, on PASS, enter executable code (`acs/v2/active_contour.py` + `acs/v2/dynamics/active_contour.py`).

P1 alpha scope (from forward roadmap): single-cell free-boundary 2D
active contour with cortex + area force only. NO protrusion, NO FA,
NO ECM coupling, NO junction, NO cell cycle in P1.

---

## 1. Q1 Cortex line tension — energy-gradient form

### Locked

**Energy**: `E_c = λ_c * P` where `P = Σ_i |r_{i+1} - r_i|` is perimeter.

**Per-vertex force**: `F_i^c = -∂E_c/∂r_i = λ_c * (t_i - t_{i-1})` where
`t_i = (r_{i+1} - r_i) / |r_{i+1} - r_i|` is the unit tangent on edge `i→i+1`.

**Units**:
- `r_i` in μm
- `P` in μm
- `λ_c` in **nN** (note: NOT nN/μm; this is the perimeter coefficient,
  not the literature line tension density)
- `E_c` in nN·μm
- `F_i^c` in nN

**Why this form**: variational, energy-monotone, avoids the April 2026
curvature anti-pattern (proof-at-peak vs measurement-band-averaged
gap). In the smooth limit equivalent to `λ_c · κ · n̂ · ℓ`, but the code
must NOT estimate curvature first.

### Literature mapping

Literature reports cortex/surface tension as `σ_c [nN/μm]`. The reduction
to the implemented coefficient is:

`λ_c = σ_c * H_eff`

where `H_eff [μm]` is the effective vertical integration height of the
2D contour (top-down projection of a 3D cell with thickness).

**Anchors**:
- Biological `σ_c` value: Salbreux/Charras/Paluch 2012 *Trends Cell Biol*,
  Chugh et al. 2017 *Nat Cell Biol*, Bambardekar et al. 2015 *PNAS*.
  Range: `σ_c ~ 0.1-1 nN/μm`.
- Variational vertex perimeter energy form: Farhadifar et al. 2007
  *Current Biology*, Alt-Ganguly-Salbreux 2017 *Phil Trans B*.

### Reject for P1

- Edge-spring with rest length: introduces extra knob (k + L_rest); not
  derived as exact perimeter gradient.
- 2.5D cortex sheet: requires height as state variable; out of P1 scope.
- Curvature-first: April anti-pattern risk.

---

## 2. Q2 Area-restoring force — normalized energy

### Locked

**Energy**: `E_A = 0.5 * K_A * (A - A_0)² / A_0`

**Pressure-like scalar**: `p_A = ∂E_A/∂A = K_A * (A - A_0) / A_0`

**Per-vertex force**: `F_i^A = -p_A * ∇_i A`
where for CCW polygon `∇_i A = 0.5 * (y_{i+1} - y_{i-1}, x_{i-1} - x_{i+1})`.

**Units**:
- `A`, `A_0` in μm²
- `(A - A_0)² / A_0` in μm²
- `K_A` in **nN/μm** (force per length, not nN/μm³)
- `E_A` in nN·μm
- `p_A` in nN/μm
- `∇_i A` in μm
- `F_i^A` in nN

**Sense**:
- `A > A_0`: `p_A > 0`, `F_i^A` inward (compression)
- `A < A_0`: `p_A < 0`, `F_i^A` outward (expansion)

### Runtime config

Single entry: `k_a_nN_per_um`. No H_eff pattern, no derived form.

### Literature mapping (docs only, not code)

Farhadifar 2007 vertex model uses unnormalized energy `E = 0.5 * K * (A - A_0)²`
with `K [nN/μm³]`. Conversion:

`K_A = K * A_0`

Example (literature K=100 pN/μm³, A_0=800 μm²):
`K_A = 100 pN/μm³ * 800 μm² = 80,000 pN/μm = 80 nN/μm`

Document `A_0` used in any literature conversion.

---

## 3. Q3 Drag — line drag density

### Locked

**Per-vertex drag**: `ζ_i = ξ_line * ℓ_i` where `ℓ_i = (|r_{i+1} - r_i| + |r_i - r_{i-1}|) / 2`
is the vertex control length (half-edges sum).

**Equation of motion (overdamped)**: `dr_i/dt = F_i / ζ_i`

**Units**:
- `ξ_line` in **nN·s/μm²**
- `ℓ_i` in μm
- `ζ_i` in nN·s/μm
- `F_i / ζ_i` in μm/s

### Why line drag density

Uniform per-vertex mobility `μ` is grid-dependent: doubling vertex
count keeps `μ` per vertex constant but doubles total drag for the
same physical contour, slowing the contour artifically. Line drag
density is grid-invariant: `ζ_i` halves when `ℓ_i` halves, so
`F_i / ζ_i` velocity is invariant under refinement.

### Literature mapping

Substrate friction is typically reported as areal: `ξ_areal [nN·s/μm³]`.
Reduction to line drag:

`ξ_line = ξ_areal * H_eff`

**Anchor**: Trichet et al. 2012 *PNAS* (substrate stiffness sensing,
effective viscous regime). Range: `ξ_areal ~ 10-100 nN·s/μm³` (broad,
cell-type and substrate dependent).

### Runtime config — two-path with provenance

Either:
- (a) Direct: `xi_line_nN_s_per_um2` supplied
- (b) Derived: `xi_areal_nN_s_per_um3 + effective_height_um` supplied,
  code computes `xi_line = xi_areal * H_eff`

If both supplied: verify agreement within float64 roundoff, accept,
record provenance "both supplied, verified".

If neither: reject with explicit error.

If both `cortex (lambda_c derived)` and `drag (xi_line derived)` use
derived path, they share the SAME `effective_height_um` — single
canonical H per `ActiveContourParameters` instance.

---

## 4. Q4 Timestep derivation — analytic Gershgorin row-sum

### Locked rule

**Stability bound**: For overdamped explicit Euler, `dt < 2 / rate_max`
where `rate_max = max_i (1/ζ_i) * Σ_j |∂F_i/∂r_j|`.

**Runtime requirement**: `dt_cell * rate_max ≤ 0.5` (named
accuracy/stability margin, **not** calibrated to pass any gate, **not**
adjustable to make tests pass).

If `dt * rate_max` violates 0.5 at runtime, fail with clear error message
suggesting smaller `dt_cell` or different integrator. Do NOT auto-adjust.

### Component bounds (force-Jacobian rows vs rate rows)

**Strict separation** (per impl-Codex review id=983 → 985):

- The **force-Jacobian row sum** `J_i = Σ_j |∂F_i/∂r_j|` has units
  `[nN/μm]`. It is what the analytic Gershgorin bound aggregates,
  and it is computed *before* dividing by drag.
- The **rate row** `r_i = J_i / ζ_i` has units `[1/s]`. It is the
  per-vertex rate that enters `rate_max = max_i r_i`. The division
  by `ζ_i` happens **once**, at the rate-row step.

**Cortex Jacobian row**:
- `∂F_i^c / ∂r_i` block norm ≤ `1/ℓ_i + 1/ℓ_{i-1}` per axis component
- Row sum: `J_i^c = Σ_j |∂F_i^c/∂r_j| ≤ 2 λ_c (1/ℓ_i + 1/ℓ_{i-1})`
  in `[nN/μm]`
- Globally: `J^c_global ≤ 4 λ_c / ℓ_min` in `[nN/μm]`

**Area Jacobian row** (exact since area is bilinear in `i±1`
neighbors):
- `dF_i/dr_j = -(K_A/A_0) g_i ⊗ g_j - p_A * H^A_{ij}` where
  `g_i = ∇_i A` and `H^A_{ij} = ∂²A/∂r_i ∂r_j`
- `H^A` only touches `i-1, i, i+1` (bilinear)
- Row sum: `J_i^A = (K_A/A_0) |g_i| Σ_j |g_j| + |p_A| · H^A_row_i`
  in `[nN/μm]`
- `H^A_row_i` is exact, computable in closed form per polygon.

**Combined rate**: `rate_max = max_i (J_i^c + J_i^A) / ζ_i`. Note
`(J_i^c + J_i^A)` is summed in the Jacobian-row units `[nN/μm]`
*then* divided by `ζ_i [nN·s/μm]` once to give `[1/s]`. Do not
divide either component by `ζ_i` before the sum.

### Metadata contract

Record per simulation step (or per chunk):
- `dt_cell`
- `rate_max` (analytic bound used)
- `dt_cell * rate_max`
- ratio `(dt_cell * rate_max) / 0.5` (must be ≤ 1)

Production tuning (Lanczos eigenvalue estimate for tighter bound) is a
separate post-P1-alpha unit.

---

## 5. Q5 Vertex spacing — fixed uniform, no remeshing

### Locked

P1 alpha:
- Fixed uniform `N` vertices.
- No insertion / deletion / adaptive refinement.
- Constructor from regular polygon or ellipse, fixed `N`.

### Validation per timestep (or every K timesteps)

- `ℓ_min` above floor (per `MeasurementBoundary` `min_edge_um`)
- Polygon simple (no self-intersection)
- Vertex coordinates finite
- No vertex coalescence below `min_edge_um`

If validation fails, reject step (do not auto-repair) — surface as
runtime error for design-discussion to triage.

### Adaptive vertex density

Out of P1 alpha. Belongs to a separate unit after `ProtrusionEventScheduler`
exists (so adaptive refinement is filopodia-driven, not arbitrary).

---

## 6. Q6 Tests — derived, no arbitrary percent gates

### Test 1: zero-force baseline

`λ_c = 0`, `K_A = 0` (or equivalently no force computation enabled).

Pass criteria:
- Vertices do not move except float64 roundoff
- Export through `MeasurementBoundary` round-trip
- Shoelace area / perimeter match the raw state round-trip

### Test 2: area-only

`λ_c = 0`, `K_A > 0`.

Pass criteria:
- If `A_init = A_0`: no motion (within float64)
- If `A_init ≠ A_0`: area force pushes correct direction (Sanity Gate §5)
- Energy `E_A` is **non-increasing each step, strictly decreasing
  only while `||F||_2 > 0`** (equilibrium steps must satisfy
  `E_A(t+1) == E_A(t)` to within float64 round-off, not strict
  inequality)
- Equilibrium: `A → A_0` exponentially

### Test 3: cortex-only

`λ_c > 0`, `K_A = 0`.

Pass criteria:
- Regular polygon shrinks self-similarly (no shape distortion)
- Energy `E_c` is **non-increasing each step, strictly decreasing
  only while `||F||_2 > 0`**
- No area preservation (`A → 0` over long time)

### Test 4: coupled ellipse → discrete equilibrium residual — DIAGNOSTIC

**Status amendment** (impl-Codex review id=985, before code lands):
This test was originally drafted with `ε_float64_envelope` and a
`(ℓ/R)²` anisotropy tolerance as gating pass criteria. On
implementation review, those tolerances were found to combine three
independent error orders (float64 round-off accumulation across
`T_steps · N` vertex evaluations, forward-Euler local truncation
across `T_steps`, and `O((ℓ/R)²)` polygon-vs-circle finite-N
discretisation) into a single named constant whose O(1) coefficients
were not derivable in the alpha window without introducing a
hidden tolerance knob — itself another Sanity Gate per Plan §11
"requires Sanity Gate derivation before execution". Test 4 was
therefore **downgraded to DIAGNOSTIC (non-gating)** in
`docs/v2/v2_p1_active_contour_sanity_gate.md` §4. Promotion to a
PASS/FAIL gate is a separate post-alpha unit (analytic envelope
or Lanczos eigenvalue estimate).

`λ_c > 0`, `K_A > 0`. Initial: ellipse polygon (e.g., 1.5×1 aspect, N=64).

**Reference equilibrium**: continuum solution `R*` from
`K_A * (πR² - A_0)/A_0 + λ_c/R = 0` (quartic, solve via numpy.roots,
take the positive real root closest to `R_init`).

**Discrete baseline**: build a regular N-gon at radius `R*`, evaluate
the IMPLEMENTED force residual `||F||` on that discrete polygon.
This is the `finite_N_residual` baseline (recorded as a diagnostic).

**Diagnostic checks (no PASS/FAIL gate on Test 4 alone)**:
- Energy `E_c + E_A` is non-increasing each step (strictly decreasing
  only while `||F||_2 > 0`).
- Final `||F||_final` and `finite_N_residual` are recorded; their
  ratio and difference are reported as diagnostics (no tolerance
  comparison gates the run).
- Final shape anisotropy `(R_max - R_min) / R̄` recorded as diagnostic.
- Final area is recorded alongside `A* = πR*²`.

**No handpicked tolerance** survives in the gate as a result of the
downgrade.

---

## 7. Implementation contract

`ActiveContourParameters` (in `acs/v2/active_contour.py`):

```python
@dataclass(frozen=True, slots=True)
class ActiveContourParameters:
    # Cortex (one of):
    lambda_c_nN: Optional[float] = None          # direct
    sigma_c_nN_per_um: Optional[float] = None    # derived path component
    # Drag (one of):
    xi_line_nN_s_per_um2: Optional[float] = None  # direct
    xi_areal_nN_s_per_um3: Optional[float] = None # derived path component
    # Shared (required if any derived path is used):
    effective_height_um: Optional[float] = None
    # Area:
    k_a_nN_per_um: float = ...     # required
    target_area_um2: float = ...   # required
    # Timestep:
    dt_cell_s: float = ...         # required, validated against rate_max at runtime
    # Vertex:
    n_vertices: int = ...          # required
    min_edge_um: float = 1e-3      # default carried from MeasurementBoundary

    def validate(self) -> ActiveContourParameters:
        # cortex path:
        cortex_direct = self.lambda_c_nN is not None
        cortex_derived = (self.sigma_c_nN_per_um is not None and
                          self.effective_height_um is not None)
        if not (cortex_direct or cortex_derived):
            raise V2ContractError("supply lambda_c_nN OR (sigma_c_nN_per_um + effective_height_um)",
                                  failure_kind="cortex_underspecified")
        if cortex_direct and cortex_derived:
            derived_lambda = self.sigma_c_nN_per_um * self.effective_height_um
            if abs(self.lambda_c_nN - derived_lambda) > 1e-9 * abs(self.lambda_c_nN + derived_lambda):
                raise V2ContractError("lambda_c_nN and sigma_c*H_eff disagree",
                                      failure_kind="cortex_overspecified_disagreement")
        # ... (same pattern for drag)
        # k_a, target_area, dt_cell, n_vertices basic checks
        return self

    @property
    def lambda_c_resolved_nN(self) -> float:
        if self.lambda_c_nN is not None:
            return self.lambda_c_nN
        return self.sigma_c_nN_per_um * self.effective_height_um

    # ... similar for xi_line_resolved_nN_s_per_um2
```

`ActiveContourState` (mutable simulation state):
- `vertices_xy_um: np.ndarray (N, 2)`, float64
- `to_measurement_boundary(source_modality)` → frozen `MeasurementBoundary`

`acs/v2/dynamics/active_contour.py`:
- `compute_cortex_forces(state) -> np.ndarray (N, 2)`
- `compute_area_forces(state, params) -> np.ndarray (N, 2)`
- `compute_rate_max(state, params) -> float` (analytic Gershgorin bound)
- `step(state, params) -> ActiveContourState` (overdamped Euler with
  rate_max validation)

### P1 alpha test isolation

Active contour unit tests MUST use:
- `cell.protrusions = ()`
- `cell.adhesions = ()`
- `cluster.junctions = ()` (if cluster-level test)

Active contour code MUST assert these are empty in P1 alpha or skip
without consuming them.

---

## 8. Sanity Gate § cross-reference (impl-work rewrites Sanity Gate doc against this)

| Sanity Gate Section | Item from this lock |
|---|---|
| §1 dimensional | unit chains for λ_c, K_A, ξ_line, dt — all derived above |
| §2 boundary cases | N≥3, finite, simple polygon, ℓ_min above floor (Q5) |
| §3 conservation | energy monotonicity (tests 2/3/4), area conservation N/A in P1 |
| §4 numerical | dt·rate_max ≤ 0.5 with metadata, float64 (Q4) |
| §5 sign | cortex contractile, area inward when A>A_0, drag dissipative |
| §6 measurement protocol | `to_measurement_boundary` round-trip (test 1), exported polygon CCW canonical |

---

## 9. References (locked)

- Salbreux G, Charras G, Paluch E. *Trends Cell Biol* 2012; 22(10): 536-545.
- Chugh P et al. *Nat Cell Biol* 2017.
- Bambardekar K et al. *PNAS* 2015; 112(5): 1416-1421.
- Farhadifar R et al. *Current Biology* 2007; 17(24): 2095-2104.
- Alt S, Ganguly P, Salbreux G. "Vertex models" *Phil Trans B* 2017.
- Trichet L et al. *PNAS* 2012; 109(18): 6933-6938.
- CLAUDE.md Rule 10 (dimensional comparison verification, anti-pattern
  Path C `g_star=0.1` per-volume vs per-area mismatch)
- CLAUDE.md Rule 11 (sim-experiment measurement matching, MeasurementBoundary
  top-down convention)

---

## 10. Cross-room dispatch

This file is the design-team input to implementation-work for:
1. P1 Sanity Gate doc rewrite (impl Claude, fresh from this lock, NOT
   patching the prior stale draft)
2. impl Codex review (verify unit reductions, no hidden defaults, no
   arbitrary tolerances, dt metadata contract)
3. On Sanity Gate PASS: executable code in `acs/v2/active_contour.py` +
   `acs/v2/dynamics/active_contour.py` + tests
4. On Sanity Gate FAIL: blocker to PI approval queue + design-discussion
   re-open for fallback options

Rounds 1-5 of the design lock are MCP id 943-964.
