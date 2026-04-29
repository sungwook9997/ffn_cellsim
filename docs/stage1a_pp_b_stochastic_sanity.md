# Stage 1a++.b sanity gate — Stochastic boundary events

PI directive 2026-04-29 (autonomous sequential implementation).
Implements `docs/marangoni_review.md` Option β (Codex review item 1):
discrete lamellipodia / filopodia / leader-cell events at the
boundary band, complementing the continuum Layer 2 active stress.

## Scope of code changes

1. Add SolverConfig fields:
   - `layer2_b_enabled: bool = False`
   - `lambda_lam_star: float = 0.0` (event rate, dimensionless 1/τ_relax)
   - `impulse_lam_star: float = 0.0` (per-event tangential impulse)
2. New kernel `_apply_stochastic_events` called once per step after
   `_g2p_and_constitutive`:
   ```
   for each particle p:
     if p is in contact band (z_p < h_band):
       u ~ uniform(0, 1)
       if u < lambda_lam · dt:
         # Outward tangential direction: project (x_p − COM) onto
         # the local tangent plane (perpendicular to grid_normal).
         radial_dir = (x_p − COM_xy)  # in xy plane, outward
         apply Δv_p = impulse_lam · radial_dir / ||radial_dir||
   ```
3. New invariants in `invariants()`:
   - `n_lam_events_step`: # events fired this step (diagnostic)
   - `lam_event_count_total`: cumulative count
4. New gate:
   - `lambda_lam · dt ≤ 0.5` (Poisson stability)
   - `n_lam_events_step` finite
   - `<lam_event_rate>` matches `lambda_lam` × n_band over time

## Sanity Gate checks

### 1. Dimensional analysis
- `lambda_lam_star` has dimension 1/time (in star units). Multiplied
  by `dt` gives dimensionless probability.
- `impulse_lam_star` has dimension velocity (in star units; same as v).
- COM = (mean x, mean y) over all particles, computed host-side per
  step (small overhead).
- **CFL check**: Poisson stability requires `lambda_lam · dt ≤ 0.5`
  (otherwise multiple events per step possible; Bernoulli underestimate).
  Default lambda_lam=0 trivially satisfies. For non-zero lambda_lam,
  runtime enforce.

### 2. Boundary cases
- `lambda_lam = 0`: no events ever fire → backwards-compat with
  legacy continuum-only Layer 2.
- `n_band = 0`: nothing in contact band → no events fire (loop skips).
- `radial_dir = 0` (particle exactly at COM xy): degenerate; skip
  with `||radial_dir|| < ε` guard. Numerically vanishing measure.
- Particle migrating out of contact band: stops getting events,
  c_act preserved per Stage 1b.b semantics.

### 3. Conservation invariants
- Stochastic events inject *momentum* (impulse) into the particle.
  Total system momentum is NOT conserved by these events — they are
  external active-matter forcing, analogous to the continuum Layer 2
  active stress (which also injects momentum).
- However, the events should be *isotropic* in the xy plane on
  average (each particle has a different radial direction, summing
  to ~0 over a symmetric initial pack). The horizontal momentum
  drift gate F2 (Option F Week 2 ACCEPTED-LIMITATION) covers this.
- Vertical (z) component of impulse: by construction, radial_dir is
  in the xy plane (z=0), so no vertical contribution. Substrate
  reaction is unaffected.
- **PASS** — no conservation regression beyond what's already accepted.

### 4. Numerical sanity
- `Δv = impulse_lam · dir`: bounded by `impulse_lam_star` per event.
  For `impulse_lam = 0.01` and at most ~lambda_lam·dt fraction of
  band particles per step, the per-step velocity injection is small
  compared to the continuum spreading dynamics.
- Random number sampling: `ti.random()` in Taichi is per-thread
  PRNG; deterministic given seed.
- **PASS**.

### 5. Sign / sense check
- Outward radial impulse on contact-band particles: pushes the
  spheroid edge outward, increasing the spreading area. Matches the
  experimental observation that lamellipodia extend the leading edge.
  **PASS**.

### 6. Measurement-protocol consistency
- Gate `lambda_lam · dt ≤ 0.5`: numerical stability, well-defined.
- Gate `<lam_event_rate>` ≈ `lambda_lam` × <n_band>: a self-
  consistency check on the Poisson process (mean event rate equals
  the configured rate × eligible-particle count). Measured over
  the entire run.
- **PASS**.

## Magic-Number Block

| Constant | Default | Test 1 | Test 2 | Test 3 |
|---|---|---|---|---|
| `lambda_lam_star` | 0.0 | YES — non-zero requires PI sweep + literature anchor (Mattila & Lappalainen 2008 *NRMCB* IF 113: lamellipodia events 1-10/min/cell → in star units: 1-10 / 60s / τ_relax ≈ 1.7e-2 to 1.7e-1) | YES — rate per particle, scale-invariant | NO — default off |
| `impulse_lam_star` | 0.0 | YES — non-zero requires PI sweep + literature anchor (lamellipodia force ~50 pN, cell mass ~10⁻¹² kg, Δv ~ 50pN·dt/m ~ 5 mm/s for dt~1ms; in star units, dimensionless O(1e-3 to 1e-2)) | YES — per-event impulse, scale-invariant | NO — default off |

Both default 0; non-zero values must be Magic-Number-Block PARTIAL
justified per ζ_star Option α' precedent.

## Configuration

```yaml
layer2_b:
  enabled: true               # NEW: enable stochastic events
  lambda_lam_star: 0.05       # NEW: rate per particle per τ_relax
  impulse_lam_star: 0.005     # NEW: per-event Δv
```

## Cross-references

- `docs/marangoni_review.md` Option β
- `docs/codex_review_synthesis.md` item 1
- `docs/layer3_phi_audit.md` (c_act_p semantics under stochastic events)
- Mattila & Lappalainen 2008 *Nat Rev Mol Cell Biol* **9** 446
  (lamellipodia / filopodia mechanism review)
- Khalilgharibi et al 2019 *Nat Phys* **15** 839 (collective tissue
  mechanics with discrete events)
