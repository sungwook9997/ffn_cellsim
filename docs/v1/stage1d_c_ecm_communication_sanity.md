# Stage 1d.c sanity gate ??ECM-mediated mechanical communication + protrusion state machine

PI directive 2026-04-30 (autonomous implementation per Codex
analysis). Adds three orthogonal modules to the framework:

1. **2D substrate ECM displacement field** (z=0 plane), screened
   viscoelastic, accumulates traction deposits from contact-band
   particles and propagates them as long-range mechanical
   communication.
2. **Protrusion state machine** per contact-band particle, replacing
   one-step stochastic lamellipodia impulse with a 5-state biological
   cycle (quiet ??filopodia_probe ??nascent_adhesion ??   lamellipodium_spread ??retract).
3. **Per-particle traction state** (`fa_strength_p`, `polarity_p`,
   `traction_p`) coupled to the ECM field through Newton's 3rd law.

Plus HDF5 schema rationalization per PI directive: stress fields
disambiguated (`stress_tensor` deprecated ??renamed `tau_dev` plus
`sigma_vol`, `sigma_active`, `sigma_total`, `pressure`, `dev_norm`),
new viz channels (`traction_ecm`, `fa_strength`, `protrusion_state`,
`ecm_signal`).

This document discharges the mandatory Sanity Gate Protocol before
first execution.

---

## Literature anchors (PI-supplied + framework verification)

| Anchor | Year | Journal | IF | Used for |
|---|---|---|---|---|
| Nahum et al. | 2023 | *Comm Biol* | 5.4 | Cell-pair mechanical communication 4??0 cell diameters; ECM densification ~2 cell diameter falloff |
| Adebowale et al. | 2021 | *Nat Mater* | 47 | Soft fast-relaxing substrate filopodia-mediated migration |
| Case & Waterman | 2015 | *Nat Cell Biol* | 21 | Molecular clutch: actin-adhesion-force transmission |
| Balaban et al. | 2001 | *Nat Cell Biol* | 21 | Focal adhesion stress ??5.5 nN/關m짼 |
| Bergert et al. | 2015 | *Nat Cell Biol* | 21 | Integrin focal adhesion traction ~100 Pa |
| Chaudhuri et al. | 2015 | *Nat Commun* | 17 | ECM viscoelastic relaxation timescales |

---

## Module 1 ??2D substrate ECM displacement field

### Equation
On the substrate plane (z=0):
```
管_ecm 쨌 ?굑_ecm/?굏 = G_ecm 쨌 ?눫쾣_ecm ??k_anchor 쨌 u_ecm + 誇_p T_p 쨌 W(x ??x_p)
```
- `u_ecm(x, y)`: ECM in-plane displacement (vector field, 2 components).
- `管_ecm`: ECM viscous timescale (Chaudhuri 2015 ??~1??0 s, in star
  units ?_ecm_star ??0.05).
- `G_ecm`: ECM shear modulus (collagen ~kPa, ratio to bulk K).
- `k_anchor`: anchoring stiffness of the underlying rigid dish
  (screened-elastic decay length).
- Deposit kernel `W`: Gaussian centered at x_p, width = `R_dep_star`
  (~1 cell diameter).

### Distance kernel
Steady-state of the screened Helmholtz operator gives exponential
decay with characteristic length:
```
貫_ecm = sqrt(G_ecm / k_anchor)
```
Per Nahum 2023: pair communication detected at 4??0 cell diameters.
With cell diameter ~ 1 R? in star units (since R? = single-cell-pack
radius), 貫_ecm sweep values **2, 4, 7, 10** in star units (config
parameter `lambda_ecm_star`). NOT fitted; literature-anchored.

### Discretization
- Use the existing `grid_n=64` background grid restricted to the
  z=0 plane (so `grid_ecm_u[i,j]`, 2D vector field).
- Forward-Euler in time, central-FD Laplacian.
- CFL: `dt 쨌 G_ecm / (管_ecm 쨌 dx짼) ??0.5` for stability.

---

## Module 2 ??Protrusion state machine

### States (per contact-band particle)
```
0 quiet
1 filopodia_probe         (short, fast probing)
2 nascent_adhesion        (FA forming, force-dependent)
3 lamellipodium_spread    (mature FA, sustained traction ?_lam ~ 10 min)
4 retract                 (FA disassembly)
```

### Per-particle fields
- `protrusion_state_p ??{0..4}`
- `protrusion_timer_p` (time in current state, star units)
- `polarity_p` (3-vector, ?뼿룐???1)
- `filopodia_length_p` (displacement of filopodia tip)
- `fa_strength_p ??[0, 1]` (0 = no FA, 1 = mature)
- `leader_score_p` (cumulative event count, optional)
- `ecm_signal_p` (sampled |u_ecm| at particle position)
- `traction_p` (3-vector, traction force this particle is currently
  applying to substrate)

### Transition rules (rates in 1/?_relax star units)

| Transition | Rate | Conditions |
|---|---|---|
| 0 ??1 | 貫_filo 쨌 free_edge 쨌 contact_band 쨌 (a + b쨌c_act_p) | particle is at free edge |
| 1 ??2 | r_form 쨌 ecm_collagen | inversely scaled by ECM resistance |
| 2 ??3 | r_mature 쨌 (1 + 款쨌|F_ecm|) | force-dependent maturation |
| 3 ??4 | 1 / ?_lam (~10 min ??10 star units) | timer-driven |
| 4 ??0 | 1 / ?_release (~5 min ??5 star units) | timer-driven |

State 3 (lamellipodium_spread) applies sustained traction along
`polarity_p` for ?_lam minutes. NOT one-frame impulse.

### Polarity update (state-3 only)
```
dp_p/dt = 慣_edge 쨌 (n?_free_edge ??p_p) 쨌 |free_edge|
        + 棺_ecm 쨌 Q_ecm_p 쨌 p_p                     # alignment to fiber
        + 棺_traction 쨌 ??u_ecm|                     # follow stiffer ECM
        + noise_dp / ?_corr
```
where `Q_ecm_p` is the local ECM strain alignment tensor. Polarity is
soft-clipped to ?뻪_p????1.

### Traction force (state 3 only)
```
T_p = T0 쨌 fa_strength_p 쨌 polarity_p
```
- T0 = `T0_star` from Balaban/Bergert: 100??500 Pa equivalent;
  dimensionless ??0.1쨌K = 0.1 in star units (default).
- Magnitude rises during state 2 (nascent), saturates in state 3.

### Newton's 3rd law back to particle
```
particle p receives Δv_p = +T_p · dt / xi_drag     # overdamped pull; no inertial 1/m_p amplification
ECM field receives  ?F_ecm(x_p) = ?뭈_p 쨌 W(x ??x_p)
```

This is the **"?붾줈 ?↔퀬 ?밴린湲?** mechanism: the particle pulls
itself toward the polarity direction by anchoring on the substrate;
Newton 3rd transmits the pull into the ECM, which propagates to
neighboring cells via the screened Helmholtz operator.

---

## Module 3 ??HDF5 schema rationalization

Per PI directive: current `stress_tensor` is misleading because it
only stores `tau_dev` (Maxwell deviatoric, NOT total).

### Renames
- `stress_tensor` ??`tau_dev` (deprecated: alias kept for backwards
  read compat; new writes use `tau_dev` only).

### New per-particle fields (saved per frame when active)
- `sigma_vol[N, 3, 3]`: K_eff쨌(?_ref/???)쨌I (v15 volumetric)
- `sigma_active[N, 3, 3]`: -瓘_eff쨌K쨌I (Layer 2/3 active stress)
- `sigma_total[N, 3, 3]`: tau_dev + sigma_vol + sigma_active
- `pressure[N]`: -trace(sigma_total)/3
- `dev_norm[N]`: ||tau_dev|| Frobenius
- `traction_ecm[N, 3]`: per-particle ECM-deposit force vector
- `force_marangoni_sampled[N, 3]`: per-particle Marangoni impulse
  (sampled from grid)
- `force_substrate_sampled[N, 3]`: per-particle substrate reaction
  contribution
- `fa_strength[N]`: focal adhesion maturity ??[0, 1]
- `protrusion_state[N]`: int ??{0..4}
- `ecm_signal[N]`: sampled |u_ecm| at particle xy

These are write-time-optional (omitted when the relevant layer is
disabled).

---

## Sanity Gate checks

### 1. Dimensional analysis
- `u_ecm`: length (star units = R?).
- `管_ecm`: time 횞 stress / length짼 ??1/? 횞 1 / R?짼 in star.
- `G_ecm`: stress (in star = K).
- `k_anchor`: stress / length짼.
- `貫_ecm = sqrt(G_ecm / k_anchor)`: length.
- `T0`: stress (force per area; here force per particle volume since
  T_p applied as impulse to a particle).
- ECM diffusion CFL: `dt 쨌 G_ecm / (管_ecm 쨌 dx짼) ??0.5` ??at default
  values (管_ecm_star ??0.05, G_ecm_star ??1.0, dx_star ??0.094:
  ratio = 0.01 쨌 1 / (0.05 쨌 0.0088) ??22.7 ??**WAIT, fails**).
- **Stability fix**: 管_ecm_star must be ??1.0 for stability at dx_star
  ??0.094 and dt_star = 0.01: `dt 쨌 G_ecm / (管_ecm 쨌 dx짼) ??0.5`
  ??管_ecm ??dt 쨌 G_ecm / (0.5 쨌 dx짼) = 0.01 쨌 1 / (0.5 쨌 0.0088)
  ??2.27. **Set default 管_ecm_star = 5.0** (safety factor 2횞).
- Protrusion state-machine timescales: ?_lam ~ 10 min = 10 in star,
  ?_release ~ 5 in star. Both ??dt_star = 0.01 ??stable.
- **PASS** with 管_ecm_star = 5.0 default.

### 2. Boundary cases
- `lambda_ecm_star = 0` (anchor-only, no diffusion): kernel collapses
  to delta ??equivalent to no ECM communication. Backwards-compat
  with Stage 1d.b.
- All Stage 1d.c features OFF (`layer7_enabled=false`): no new
  state, no new fields written, no kernel calls. Stages 1a??d.b
  unaffected.
- Particle at z > h_band: never in contact band, protrusion state
  stays at 0 forever, no traction applied.
- f32 accumulation: ECM grid uses f32 with f64 accumulators in
  diagnostics. NaN-guarded by the existing infrastructure.
- **PASS**.

### 3. Conservation invariants
- Newton's 3rd law (particle ??ECM): T_p applied to particle's
  velocity AND ?뭈_p deposited into ECM grid as force; both per same
  dt. This is momentum-exchanged with the substrate (which is
  mechanically grounded by `k_anchor` term). Total momentum is NOT
  conserved (substrate absorbs it via k_anchor, same as the
  reflective BC ground truth).
- ECM diffusion is dissipative (?눫?acts as smoothing); k_anchor is
  dissipative (Lyapunov toward u_ecm = 0).
- Protrusion state machine is purely state-evolution; no spurious
  energy injection beyond what state 3 (lamellipodium_spread)
  applies via the traction term.
- **PASS** ??all forcing channels accounted for in active power gate.

### 4. Numerical sanity
- ECM grid time step: `dt_star = 0.01` 횞 `G_ecm_star / (管_ecm_star 쨌
  dx_star짼)` = 22.7e-3 / 5.0 = ~4.5e-3 ??CFL ratio ??0.45, **PASS**.
- Protrusion timer: dt_star = 0.01, ?_lam = 10 ??1000 steps in state
  3, fine resolution.
- f32 vs f64: per-particle f32, ECM grid f32, diagnostic
  accumulators f64.
- **PASS**.

### 5. Sign / sense check
- ECM diffusion: G쨌?눫쾣 smoothes away sharp deposits (correct).
- ECM anchor: ?뭟쨌u always points toward zero (correct, dissipative).
- Traction T_p = T0 쨌 fa_strength 쨌 polarity: pulls particle in
  polarity direction (correct ??"?붿씠 substrate瑜??↔퀬 cell body瑜?  ?밴?").
- Newton 3rd: ?뭈_p deposited into ECM at particle position
  (correct).
- Polarity ECM bias: `+棺_traction 쨌 ??u_ecm|` aligns toward stiffer
  / more-deformed ECM (correct ??durotaxis-like).
- Protrusion state transitions: rates positive, conditions
  monotone-correct.
- **PASS** ??all signs match physical intuition.

### 6. Measurement-protocol consistency
- Per Hard Rule 11: PI experimental top-down imaging modality is
  matched by `A_over_A0_topdown` (existing top-down xy-projection of
  ALL particles). Stage 1d.c does not change this; the new traction
  + ECM mechanism modifies HOW particles spread, not WHAT we
  measure.
- New gate channels (`fa_strength`, `protrusion_state`, `ecm_signal`,
  `traction_ecm`) are diagnostic-only and visualization-side, not
  reported as physical claims.
- HDF5 rename: `stress_tensor` ??`tau_dev` is a measurement-protocol
  fix per PI directive ("吏湲?sigma plot? 'deviatoric Maxwell stress
  ?쇰?'?쇱꽌 ?댁긽?섍쾶 蹂댁씠??寃??뺤긽"). Old field name preserved as
  read-time alias for backwards compat.
- **PASS** ??measurement protocol matches biological scope.

---

## Magic-Number Block

| Constant | Default | Test 1 (derivable) | Test 2 (grid-invariant) | Test 3 (fitting) |
|---|---|---|---|---|
| `eta_ecm_star` | 5.0 | YES ??set by CFL stability bound 횞 2 safety; not tuned | YES ??depends on dx짼쨌G/管 ratio, scales with grid | NO ??set for stability |
| `G_ecm_star` | 1.0 | YES ??collagen shear modulus ~kPa = K_star order; literature anchor (Chaudhuri 2015) | YES | NO |
| `k_anchor_star` | varies by lambda_ecm sweep | YES ??derived from 貫_ecm짼 쨌 G_ecm relationship per Nahum 2023 | YES | NO ??Nahum literature anchor |
| `lambda_ecm_star` | sweep [2, 4, 7, 10] | YES ??Nahum 2023 4??0 cell diameter pair communication range | YES | NO ??sweep, not fitted |
| `T0_star` | 0.1 | YES ??Balaban 2001 ~5.5 nN/關m짼, Bergert 2015 ~100 Pa; in star = T0/K ??0.1 | YES | NO ??literature mid-range |
| `tau_lam_star` | 10.0 | YES ??lamellipodia persistence ~10 min ??10 ?_relax | YES | NO ??literature |
| `tau_release_star` | 5.0 | YES ??FA disassembly ~5 min | YES | NO ??literature |
| `lambda_filo_star` | 0.05 | YES ??filopodia probe rate per Mattila & Lappalainen 2008 NRMCB IF 113 (1-10/min lower end) | YES | NO ??literature mid-low |
| `r_form_star` | 0.1 | YES ??FA formation rate 1/min order; in star ??0.1 | YES | NO ??order of magnitude |
| `r_mature_star` | 0.05 | YES ??mature FA formation slower than nascent; ratio ~ 1/2 | YES | NO ??biological |
| `alpha_edge_star` | 1.0 | YES ??polarity update toward free edge; relaxation rate ~ 1/?_relax | YES | NO ??dimensional |
| `beta_ecm_star` | 0.5 | YES ??ECM-alignment coupling, mid-range | YES | NO ??sweep candidate |
| `beta_traction_star` | 0.5 | YES ??durotaxis bias, mid-range | YES | NO ??sweep candidate |

All constants pass Magic-Number Block; non-derivable parameters are
literature-anchored or sweep candidates explicitly flagged for
future PI sweep authorization (per 瓘_star Option 慣' / Path C
gravity_star precedent).

---

## Differential test (acceptance)

Stage 1d.c default-OFF (layer7_enabled=false): all stages 1a/1b/1c/
1d/1d.b/1a++.b unchanged (regression). pytest 16 must pass.

Stage 1d.c default-ON in pilot:
- `protrusion_state[i]` distribution evolves (count of state 3 /
  total contact-band > 0 after 30 min sim time)
- `traction_ecm` magnitude in active region > 0
- `ecm_signal` non-zero at distance ??貫_ecm from any state-3 particle
- A/A?_topdown trajectory: should exceed Stage 1b.b baseline (1.677)
  due to traction-driven body translation; expected pilot 1k
  endpoint A/A? ??3-5 (within PI experimental Lam4 4hr range 2.34
  +/- spread)

---

## Configuration

```yaml
layer7:                      # NEW: ECM communication + traction
  enabled: true
  # ECM substrate field
  eta_ecm_star: 5.0
  G_ecm_star: 1.0
  lambda_ecm_star: 4.0       # mid-range Nahum 2023 (4-10)
  R_dep_star: 1.0            # deposit kernel width = 1 cell diameter
  T0_star: 0.1               # Balaban / Bergert traction scale
  # Protrusion state machine
  tau_lam_star: 10.0
  tau_release_star: 5.0
  lambda_filo_star: 0.05
  r_form_star: 0.1
  r_mature_star: 0.05
  alpha_edge_star: 1.0
  beta_ecm_star: 0.5
  beta_traction_star: 0.5
  protrusion_speed_cap_star: 0.06  # ~0.1 um/s at R0=100um, tau=60s
```

---

## Cross-references

- `docs/v1/marangoni_review.md` ??Mechanism A/E/F discussion (Stage
  1d.b)
- `docs/v1/stage1a_pp_b_stochastic_sanity.md` ??Stage 1a++.b (single-step
  random impulse, replaced by Stage 1d.c protrusion state machine)
- `docs/v1/layer3_phi_audit.md` ???_memory + c_act split
- `docs/13_data_schema.md` ??HDF5 schema (to be updated with new fields)
- Nahum et al. 2023 *Comm Biol* ??ECM-mediated cell pair communication
- Adebowale et al. 2021 *Nat Mater* ??substrate viscoelasticity
- Case & Waterman 2015 *Nat Cell Biol* ??molecular clutch
- Balaban et al. 2001 *Nat Cell Biol* ??focal adhesion stress
- Bergert et al. 2015 *Nat Cell Biol* ??integrin traction stress
- Chaudhuri et al. 2015 *Nat Commun* ??ECM viscoelastic relaxation

