# WCA steric ill-conditioning — diagnosis and remediation

**Date:** 2026-07-21 · **Branch:** `ac/wca-conditioning` (worktree `ffn_ac-wca`, isolated from the
shared-tree WCA session; no `ac/cell/` edits, no commits) · **Scope:** why full-native `analytic_implicit`
never converges the resting preload with WCA on (blocker (b) on the magnitude / cortical-tension / blebbing
critical path).

> **Execution provenance (dev Mac = CPU, no CUDA).** Every quantitative claim below is either (a) a
> **closed-form** result from the repo's own `ffn_sim/ac/solid/wca_analytic.py` oracle, or (b) **CPU-executed
> NumPy** on a tractable 2–16-fiber model of the exact `(a I + P K P)` operator
> (`scratchpad/wca_conditioning_v2.py`). No Warp/CUDA kernel was run here. The native-scale confirmations are
> listed as **gbook A5000 runs to execute** in §5 — they are proposed, not performed.

---

## 0. Verdict

The WCA-on resting blocker is **not one problem, it is three coupled defects, all rooted in a single design
choice**: *node-node* WCA at the *physical σ = 7 nm* on a *0.5 µm node mesh*, closed by a *hard force cap*.

| # | Defect | Symptom in REPORT.md | Mechanism (this analysis) |
|---|---|---|---|
| **D1** | **Coverage / fidelity failure** | (silent) | node-node at σ=7 nm detects **0.1 %** of real fiber crossings; the contacts that fire are accidental node coincidences → grid-variant, not physical (§4). |
| **D2** | **Force-cap residual trap** | `2204 pN` resting residual (vs `21 pN` WCA-off); explicit stuck | the hard cap zeroes the tangent **exactly where the force is largest**; a node with ~2–3 deep capped overlaps carries ~2× f_cap of standing force that has **no curvature to relax** (§3.2). |
| **D3** | **Stiffness-contrast / κ floor** | best PCG `~7,200×` above gate; `0` implicit accepts; `k_max=2.83e6` | the 7 nm well forces `k_max = 2.83e6 pN/µm`; the regularizer `a = √eps64·k_max` pins `κ ≈ 1/√eps64 ≈ 10⁸`; and there are `O(#fibers) ≈ 10⁵` near-null modes (§3.1, §3.3). |

**Recommendation (§5): reformulate the contact — (ii) segment-segment (edge-edge) EV at a resolution-faithful
σ, with the hard force-cap replaced by a C¹ finite-stiffness linear core** — retaining the existing two-level
coarse space for the near-null modes and adding stiffness continuation to land the resting solve. The
**single highest-value, lowest-risk first step** is the linear-core replacement of the cap (D2), which is
CFL-neutral and unblocks the resting residual on its own. σ_EV is surfaced to PI as the discretization
decision `params_i0b2b.yaml` already flags — **no constant is tuned to pass a gate.**

---

## 1. The operator, and where k_max = 2.83e6 comes from

The inner solve is `(a I + P K P) dx = P F`, `P` the NF2007 inextensibility projector,
`K = K_bend + K_wca + K_xlink + …`. From `assemble.py`: production runs `σ_EV = 0.007 µm` (7 nm physical
actin), `k_ev = 1000 pN/µm`, `force_cap = 1000 pN`, node spacing `ℓ = 0.5 µm`.

Closed form (`wca_analytic`, reproduced exactly):

```
σ = 7 nm     r_c = 2^(1/6)σ = 7.857 nm     ε = 8.574e-4 pN·µm
force-cap boundary  r_eq(f_cap=1000) = 4.698 nm   (overlap = 40.2 % of r_c)
k_max = U''(r_eq)  = 2.8281e6 pN/µm      ← matches REPORT k_max = 2,828,099
bending  α = EI/ℓ³ = 0.583 pN/µm  (EI = k_BT·L_p = 4.114e-3·17.7)
WCA / bending stiffness contrast = 4.85e6
```

`k_max` is **the WCA radial tangent `U''` at the cap boundary** — the stiffest *live* (uncapped) contact the
operator can assemble. It is `4.85 × 10⁶` times the bending stiffness. This single number drives both the
explicit CFL (`dt_µ = 0.1/k_max = 3.5e-8 s`) and, through the regularizer, the linear conditioning.

**Why the well is so stiff:** at σ = 7 nm the WCA well is ~1 nm wide, so resisting even the f_cap = 1000 pN
operating load demands curvature `~10⁶ pN/µm`. This is purely a consequence of the tiny σ — a wider effective
tube needs far less curvature to hold the same load (§4, C3): σ=50 nm → k_max=3.5e5 (8× softer, 8× larger dt).

---

## 2. Two hypotheses, decided

The REPORT poses: is conditioning dominated by **(A) a few near-contact node pairs (large local Hessian)** or
**(B) long-wavelength fiber modes**? The answer from the spectrum decomposition (`range_spec` on `P K P`):
**both are real, they live at opposite ends of the spectrum, and they are DIFFERENT bottlenecks needing
different fixes.**

- **λ_max is set by the live WCA contact (Hyp A).** In the uncapped band the top eigenvalue tracks the WCA
  tangent one-for-one (§3.1): as overlap goes 2 %→39 %, `k_ax` rises `1.5e3 → 2.1e6` and `λ_max` rises
  `3.0e3 → 4.3e6`. So the stiff, localized contacts own the top of the spectrum, exactly to `k_max`.
- **λ_min is the near-null fiber mode at the regularizer floor (Hyp B), and it is PRE-EXISTING.** With WCA
  entirely off, `κ` is still `1.58e8` (§3.3) — identical to WCA-on — because `λ_min` is a rigid/translation
  fiber mode sitting at `a = √eps64·k_max`. WCA raises `λ_max` *and* `a` proportionally, so it does not change
  κ; it changes the CFL and the residual (D2).

**Consequence for the fix:** a preconditioner that only attacks the soft end (Hyp B — the coarse space) cannot
finish the job, because (i) the residual trap D2 has *zero* tangent and no preconditioner inverts a zero, and
(ii) `λ_max = k_max` remains. The critical-path fix must target the **WCA contact formulation** (D1+D2+the
`k_max` magnitude), *then* lean on the coarse space for the near-null modes.

---

## 3. Evidence (CPU-executed NumPy on the exact operator)

Model: perpendicular crossing fibers, 7 nodes each, exact NF2007 projector, discrete bending stencil
`α(v_a−2v_b+v_c)` (identical to `add_bending_stiffness_kernel`), and the axial WCA tangent identical to
`add_wca_stiffness_kernel` (capped pair → zero tangent). `a = √eps64·max|diag K|`, mirroring
`omitted_regularization_base`. Spectra are eigenvalues **restricted to range(P)** — the space CG lives in.

### 3.1 λ_max tracks the WCA tangent up to k_max, then the cap drops it to zero

```
overlap%  r[nm]  Fpair[pN] capped?    k_ax     lam_min     lam_max     kappa
   2%     7.700  1.94e-01   False   1.50e+03  2.25e-05  3.011e+03  1.34e+08
  10%     7.072  2.41e+00   False   6.76e+03  1.01e-04  1.353e+04  1.34e+08
  20%     6.286  1.76e+01   False   4.23e+04  6.31e-04  8.462e+04  1.34e+08
  30%     5.500  1.19e+02   False   2.99e+05  4.46e-03  5.985e+05  1.34e+08
  39%     4.793  7.67e+02   False   2.13e+06  3.18e-02  4.266e+06  1.34e+08
  60%     3.143  1.94e+05    True   8.05e+08  5.21e-08  8.219e+00  1.58e+08   ← capped: tangent = 0
  95%     0.393  1.07e+17    True   3.55e+21  5.21e-08  8.219e+00  1.58e+08   ← capped: tangent = 0
```

`λ_max` follows the live tangent exactly (Hyp A). **κ is pinned near `1.3–1.6e8 ≈ 1/√eps64` regardless of
overlap depth** — because `a = √eps64·k_max` scales `λ_min` up in lockstep with `λ_max`. This is a
*structural* condition-number floor, not a magnitude accident. Note the projector does **not** rescue us: a
transverse single-node displacement is first-order inextensible, so the stiff WCA mode survives projection and
`P K P` keeps `λ_max ≈ k_max` (the projector's own spectrum is a benign orthogonal projection).

### 3.2 The force-cap residual trap — the 2204 pN mechanism

The cap sets `F = f_cap`, `tangent = 0` for any overlap deeper than 40.2 %. Deep overlaps therefore carry a
standing force with **exactly zero curvature**:

```
1 capped contact              |net force| =  1000.0 pN   assembled tangent = 0   → Newton step = 0
2 capped contacts (0.3 rad)   |net force| =  1977.5 pN   assembled tangent = 0   → Newton step = 0
3 capped contacts (0.5 rad)   |net force| =  2755.2 pN   assembled tangent = 0   → Newton step = 0
```

A worst-case node with ~2–3 deep capped cross-fiber overlaps has a standing residual of ~2000–2800 pN and the
solver has **no gradient to push it out** — the explicit path is additionally CFL-throttled to `dt_µ=3.5e-8`.
This is precisely the REPORT's `2204 pN` WCA-on resting residual (≈ 2.2 × f_cap) vs `21 pN` WCA-off. **D2 is
the direct cause of the unconverged resting preload**, and it is a defect no linear preconditioner can fix
(you cannot precondition a zero tangent).

### 3.3 λ_min is a pre-existing rigid-mode floor; decoupling the regularizer collapses κ

```
bending only                              a=5.21e-08  lam_min=5.21e-08  lam_max=8.22   kappa=1.58e8
bending + WCA (uncapped)                  a=5.21e-08  lam_min=5.21e-08  lam_max=8.22   kappa=1.58e8
bending + WCA, a fixed = 1.0 (decoupled)  a=1.00e+00  lam_min=1.00e+00  lam_max=9.22   kappa=9.22
```

κ is the same with and without WCA (Hyp B is pre-existing), and setting `a` to a *physically derived* value
instead of the `√eps64·k_max` floor collapses κ from `1.6e8` to `9.2`. **The floor coupling `a ∝ k_max` is
itself a conditioning lever** — when `k_max` is a 7 nm-well artifact, the floor drags the whole spectrum into a
`1/√eps64` condition number. Lowering `k_max` (via a faithful, wider contact) or decoupling `a` from it both
help; the former is the mechanistic fix, the latter needs analytic omitted-stiffness tangents (REPORT open
item 2) to be legitimate rather than a tuned floor.

### 3.4 Near-null modes scale with fiber count → native CG stalls

```
 2 fibers ( 14 nodes): near-floor modes= 10  CG(jac)= 262  CG(+coarse)= 197
 4 fibers ( 28 nodes): near-floor modes= 20  CG(jac)= 723  CG(+coarse)= 504
 8 fibers ( 56 nodes): near-floor modes= 40  CG(jac)=1269  CG(+coarse)= 961
16 fibers (112 nodes): near-floor modes= 80  CG(jac)=1947  CG(+coarse)=1533
```

Near-floor modes grow `~5 per fiber`; Jacobi-CG iteration count grows with them; the per-fiber translation
coarse space deflates ~25 %. Extrapolated to **70,686 fibers → ~10⁵ near-null modes**, CG768 is orders short —
consistent with REPORT (two-level coarse-32 CG768 reaching only `1.07e-4`, still `~7,200×` above the
`√eps64 = 1.49e-8` gate). This is why the coarse space *helps* (0.95 → 0.0011) but cannot *finish*: it deflates
the translation part of an O(10⁵)-dimensional near-null space while `λ_max = k_max` holds the top fixed.

---

## 4. node-node vs segment-segment

### 4.1 Coverage — node-node barely enforces excluded volume at all (randomized node phase)

```
sigma=physical 7nm  r_c= 7.86nm : real crossings=6000  node-node HITS=  6   coverage= 0.1%
sigma=MT 25nm       r_c=28.06nm : real crossings=6000  node-node HITS= 72   coverage= 1.2%
sigma=coarse 50nm   r_c=56.12nm : real crossings=6000  node-node HITS=304   coverage= 5.1%
```

With nodes 0.5 µm apart and `r_c ≈ 8 nm`, two fiber centerlines can cross within the excluded shell yet the
nearest *nodes* are ~250 nm apart, so node-node fires on **0.1 %** of real crossings. The 99.9 % it misses are
genuine steric contacts silently ignored; the 0.1 % it fires are **accidental node coincidences** — which
contacts exist depends on where nodes happen to land, i.e. it is **grid-variant in the physical sense** and
violates the mechanistic-fidelity hard rule. This is D1, and it is independent of any conditioning concern:
node-node EV at physical σ is not enforcing excluded volume.

Segment-segment (edge-edge) tests the **centerline minimum distance** (the physically meaningful quantity), so
it detects 100 % of crossings by construction and is grid-invariant.

### 4.2 Stiffness concentration

```
node-node sigma= 7nm   k_max=2.828e6   peak_nodal=2.828e6   dt_mu=3.54e-8 s
node-node sigma=50nm   k_max=3.538e5   peak_nodal=3.538e5   dt_mu=2.83e-7 s
seg-seg   sigma= 7nm   k_max=2.828e6   peak_nodal=7.070e5   dt_mu=3.54e-8 s   (4× lower peak nodal)
seg-seg   sigma=50nm   k_max=3.538e5   peak_nodal=8.846e4   dt_mu=2.83e-7 s
```

node-node deposits the entire rank-1 `k_max` on one node pair. edge-edge distributes the contact to the four
segment endpoints via barycentric weights `w ≤ 1`, so peak nodal stiffness drops by up to `4×` (`w=0.5` mid-
segment → `w²=0.25`). Combined with a resolution-faithful (wider) effective σ, this lowers both peak nodal
stiffness *and* `k_max` itself (σ=50 nm → 8× softer, 8× larger dt). It does **not** by itself remove the
`κ ≈ 1/√eps64` floor (that needs the regularizer decoupling / coarse space of §3.3–3.4), but it cuts `λ_max`,
the CFL cost, and the residual-trap depth.

---

## 5. Recommendation

**Primary: (ii) segment-segment (edge-edge) WCA at a resolution-faithful σ, with the hard force-cap replaced
by a C¹ finite-stiffness linear core**, keeping the two-level coarse space (Hyp B) and adding a stiffness
continuation schedule to land the resting solve.

Ranked by impact/effort — each is a **structural** change, none tunes a constant to an outcome:

### Fix 1 (immediate, minimal, CFL-neutral) — replace the hard cap with a C¹ linear core

The hard cap `F = min(F_WCA, f_cap)` gives `tangent = 0` in the deep region (D2). Replace it, below the cap
radius `r_eq(f_cap)`, with the C¹ linear continuation
`F(r) = f_cap + k_core·(r_eq − r)`, `k_core = U''(r_eq)` — matching the WCA force **and** slope at `r_eq`. The
tangent is then `k_core` (constant, **positive**, equal to the current `k_max`) instead of zero. Effects:
`k_max` is unchanged → **dt_µ and CFL are unchanged**; the deep-overlap tangent becomes non-zero → **the 2204
pN residual becomes relaxable**. This alone should unblock the resting preload (D2) and is a ~15-line kernel
edit. SPD-preserving: the deep-region tangent is `+k_core·(û⊗û)`, the same rank-1 PSD form as the uncapped
branch.

Warp-CUDA-guarded tangent (edit to `add_wca_stiffness_kernel` / `add_wca_preconditioner_kernel`, and the
matching force in `wca_steric_kernel`):

```python
# device (Warp), f64; replaces the `capped -> contribute zero` branch
r_eq   = _r_eq_cap          # precomputed host-side: compressed_pair_separation(f_cap, sigma, eps)
k_core = _k_core_cap        # precomputed host-side: wca_local_stiffness(r_eq, sigma, eps)  == k_max
if force_magnitude >= force_cap:            # deep (was: contribute nothing)
    axial_k = k_core                        # C1 linear core: constant POSITIVE tangent, not 0
else:
    axial_k = (24.0*epsilon/(length*length))*(26.0*sr6*sr6 - 7.0*sr6)   # unchanged WCA branch
unit   = delta / length
action = action + axial_k * wp.dot(unit, vi - vector[j]) * unit
# force kernel: fmag = f_cap + k_core*(r_eq - length)  for length < r_eq   (C1 continuation, was clamp to f_cap)
```

NumPy reference gate (extends `steric_reference.py`, dev-Mac runnable): `F` and `F'` are **C¹ at `r_eq`**
(finite-difference of the new force = the new tangent to `O(h²)`), `tangent > 0` for all `0 < r < r_c` (no zero
region), force is monotone decreasing, and `k_max` (= `max tangent`) equals the pre-change value so the CFL
ledger is unchanged. Native gate (§5.4): resting max projected residual must fall from `2204 pN` toward the
`21 pN` WCA-off reference.

### Fix 2 (structural, PI-gated σ) — edge-edge geometry at resolution-faithful σ

Detect contacts on the **segment–segment minimum distance**, not node coincidence (D1). The production motor
already computes point-to-segment nearest points in `ffn_sim/ac/motor/segment_query.py` — reuse that geometry
(extend to segment-segment / clamped-parameter closest points) so the machinery is device-proven. Distribute
force and the rank-1 tangent to the four endpoints via the barycentric weights `[(1−s), s, −(1−t), −t] ⊗ û`
(§4.2 — 4× lower peak nodal stiffness). This is the mechanistically-faithful EV the CLAUDE.md worked-example
table and `params_i0b2b.yaml` (`reason_gap`, option (a)) both point to.

**σ_EV is a PI decision, surfaced not chosen.** `params_i0b2b.yaml` already flags σ_EV as a discretization GAP.
This analysis supplies the evidence PI needs: node-node at physical σ enforces ~0 % of contacts (§4.1), so the
faithful resolutions are (a) edge-edge at physical per-type σ (7/25 nm) — recommended — or (b) node-node at a
declared coarse σ ≈ node spacing. Do **not** pick silently.

### Fix 3 (solver) — keep the coarse space; add continuation; SPD contact-graph block later

Retain and extend the two-level coarse space for the O(10⁵) near-null modes (Hyp B, §3.4). To land the resting
solve from the current overlapped checkpoint, add a **continuation ramp** on the core stiffness `k_core`
(option iii): solve at a soft `k_core`, then step it up to the CFL value across a few outer iterations
(barrier-style), each stage warm-started — this globalizes the descent without a fitted trust radius. The SPD
contact-graph preconditioner (REPORT open item 1, option i) becomes worthwhile **after** Fixes 1–2 lower the
stiffness contrast; leading with it treats the symptom (κ) while leaving the physics unfaithful (D1) and the
residual trap intact (D2).

### What must NOT be done (hard rules)

- Do **not** lower `k_ev`, raise `tol`/loosen `√eps64`, cap the PCG budget as a pass, or reduce fiber density
  to converge — all forbidden (CLAUDE.md; REPORT open item 4). Fixes 1–3 are contact-model/solver
  reformulations, not constant tuning.
- Do **not** treat `a = √eps64·k_max` as physical; §3.3 shows it is a conditioning lever. Its principled
  removal is the analytic omitted-stiffness tangents (REPORT open item 2), not a hand-set floor.

### 5.4 Native gbook A5000 runs to execute (validation — NOT run here)

1. **Fix 1 alone** — full-native resting, WCA on, linear-core, explicit + `analytic_implicit`: expect resting
   max projected residual `2204 → O(21)` pN (D2 resolved); `k_max`, `dt_µ` unchanged. This is the go/no-go for
   the resting-preload unblock.
2. **Fix 1 + coarse-32 CG768** — same state: expect the relative PCG residual to fall well below the current
   `1.07e-4` now that deep overlaps carry curvature; report iterations to `√eps64` and implicit-accept count.
3. **Fix 2 (edge-edge, σ per PI)** — full-native: certify Warp edge-edge vs a NumPy segment-segment brute-force
   reference (grid-invariance re-cert, analogue of the existing hash-grid parity gate); report contact count
   (expect ≫ node-node), peak nodal stiffness, `k_max`, `dt_µ`.
4. **Fix 3 continuation** — full-native resting with the `k_core` ramp: report per-stage residual and the first
   authorized implicit accept, then a transactional multi-step run (REPORT open item 3).

Re-run the exact §"Full-native topology-aware audit" table after each, per REPORT open item 4.

---

## 6. Files

- Analysis script (CPU NumPy, reproducible):
  `ffn_sim/outputs/ac/implicit/wca_conditioning_scripts/wca_conditioning_v2.py`
  (run: `PYTHONPATH=<repo> python .../wca_conditioning_v2.py`).
- Read (not edited): `ffn_sim/ac/solid/{wca_analytic,steric_reference,steric_warp}.py`,
  `ffn_sim/ac/solid/params_i0b2b.yaml`, `ffn_sim/ac/cell/{implicit_mechanics,assemble}.py`,
  `ffn_sim/outputs/ac/implicit/REPORT.md`.
