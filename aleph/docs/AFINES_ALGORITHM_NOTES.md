# AFINES Algorithm Notes — Phase 0.3 deliverable

*ActiveCellSim v2 HOOMD-blue migration. Compiled 2026-05-19 from primary
source review of [github.com/Simfreed/AFINES](https://github.com/Simfreed/AFINES)
(master @ 2025-03-31) and the Freedman 2017 BPJ / Bieling 2016 Cell /
Funk 2021 + Li/Bieling 2022 eLife papers. Source files inspected listed at the end.*

> **Audience**: PI + future Worker A/B/C who will implement Phase 1 Units
> H.1–H.5 (ECM Mikado / single filament / cortex multi-filament / FA +
> motor-clutch / lamellipodium dendritic). This document is the design
> reference each Worker must read before opening a HOOMD file.

---

## 0. Executive summary — six things to know before any HOOMD code is written

These contradict pieces of the v2 Build Plan and need explicit PI sign-off.
Each is expanded in the relevant section below.

| # | Finding | Plan v2 assumption | Implication |
| --- | --- | --- | --- |
| 1 | Canonical AFINES repo is `github.com/Simfreed/AFINES` | `github.com/Shibalab-Lehigh/AFINES` (flagged as needing reconfirmation) | Update Plan §2 H.0.3 URL. `Shibalab-Lehigh/AFINES` does not exist. |
| 2 | AFINES motors/xlinks use **Metropolis** (Glauber detailed-balance) off-acceptance, with constant `k_off` × ΔU-dependent acceptance | "Bell-Evans break events" (force-dependent k_off(F) = k_0·exp(F·xβ/kT)) | Different physical models. Bell-Evans is the v2 *upgrade*, not a port. Need explicit decision: stay-with-AFINES (Metropolis) or upgrade to Bell-Evans. |
| 3 | **AFINES has no Arp2/3 branching code** anywhere — not in master, not in any of the 9 forks. The closest feature is the `dinner-group/AFINES@growing` branch which does **linear barbed-end elongation only** | Plan H.5 = "AFINES branched mode port" | Plan H.5 is **greenfield design** from Bieling 2016 / Funk 2021 + Li/Bieling 2022 phenomenology, not a port. Adds ~2 weeks scope vs the original "port AFINES branched mode" estimate. |
| 4 | AFINES integrator is **Leimkuhler-Matthews BAOAB-limit** (O(Δt²) on harmonic systems) | not specified in Plan (just "HOOMD Brownian") | HOOMD's `md.methods.Brownian` is Euler-Maruyama (O(Δt)). Two options: (a) accept E-M and halve Δt; (b) write a custom L-M HOOMD updater (~1–2 days). Recommend (a) for Phase 1 baseline. |
| 5 | AFINES default filament: N=11 beads, l_link=1.0 μm, R_bead=0.5 μm — these are **decoupled** (beads don't touch) | Plan v2: "ℓ₀ = 0.5 μm, 21 beads per fiber" | Probably a length-scale clarification: with L_f=10 μm and ℓ₀=0.5 μm bond rest length, N=21 beads is consistent (twice the AFINES bead density). But "ℓ₀ = 0.5 μm" might be conflated with AFINES's `R_bead = 0.5 μm`. Disambiguate. |
| 6 | AFINES is **2D, CPU-only, single-threaded** (`#pragma omp parallel for` lines are commented out in all hot loops) | Plan v2 = HOOMD GPU + 3D | We are gaining GPU + 3D + parallel for free; we are losing every AFINES-side performance optimisation we might have implicitly counted on. Wall-time projections in Plan §11 are HOOMD-from-scratch, not AFINES-port. |

The remaining sections give the per-subsystem AFINES details and the
proposed HOOMD port.

---

## 1. Repository, language, scope

| Repo | URL | State | Use it for |
| --- | --- | --- | --- |
| **Canonical** | <https://github.com/Simfreed/AFINES> | master, last commit 2025-03-31, 9 forks, GPLv3, ≈80 KB of C++ | All reference work — the only repo whose algorithm matches the published paper. |
| Dinner-lab fork | <https://github.com/dinner-group/AFINES> | 28 feature branches; default `sorting_guv`; last pushed 2021-09-20 | Source for features not in canonical: `exv.cpp` (excluded volume), `box.cpp`, refactored `quadrants.cpp`, `growing` branch (treadmilling), `growing_and_spacer` and `polymerize` branches. |
| BanerjeeLab fork | <https://github.com/BanerjeeLab/AFINES> | identical to canonical, stale (2018) | Skip. |

AFINES is **C++14, single-threaded, 2D, CPU only**, configured via plaintext
`.cfg` files and a custom main `prog/network.cpp` that hard-codes the loop
order. Typical workstation runs simulate 10²–10³ filaments × 10⁵ steps
over hours. `prog/network.cpp:407-496` is the master loop.

The per-step main loop is:

1. Trajectory write if `count % n_bw_print == 0`.
2. Apply shear delta if oscillatory/linear strain is configured.
3. `net->update()` — per filament: stretching forces → bending forces →
   Brownian-Langevin position update.
4. `net->quad_update_serial()` — rebuild neighbor-list cells
   (every `quad_update_period` steps; default 1).
5. `crosslks->motor_walk(t)` — for each xlink: Brownian-relax unbound
   heads, attach/detach (Metropolis), step bound heads (skipped because
   `p_motor_v = 0`).
6. `myosins->motor_walk(t)` — same machinery as xlinks, with stepping.
7. `net->clear_broken()` — discard fractured-filament index list.
8. `t += dt; count++;`

Note 4: AFINES quadrant rebuild is the dominant cost at large N. A HOOMD
GPU port automatically wins here.

---

## 2. Filament (bead–spring + bending)

### 2.1 AFINES implementation

Headers `include/filament.h`, `include/bead.h`, `include/spring.h`;
implementations `src/filament.cpp` (651 lines), `src/bead.cpp`,
`src/spring.cpp`. Ensemble: `include/filament_ensemble.h`,
`src/filament_ensemble.cpp`.

State per filament (`filament.h:25–137`):

```
filament                          # one filament
├── vector<bead*>   beads         # N = nmonomer; default 11
├── vector<spring*> springs       # N-1 inter-bead springs
├── vector<array<double,2>> prv_rnds   # previous Gaussian draws (Leimkuhler memory)
└── doubles: kb, temperature, dt, damp, bd_prefactor, ubend, fracture_force, ...

filament_ensemble                 # whole network
├── vector<filament*> network
├── springs_per_quad[nq[0]][nq[1]]
└── n_springs_per_quad
```

Bead (`bead.cpp:19–28`): `(x, y, rad, visc, friction = 6πη·rad, force[2])`.
`friction` is the Stokes drag coefficient stored per-bead.

Spring (`spring.cpp:21–42`): `(l0, kl, max_ext, aindex[2], hx[2], hy[2],
llen, llensq, disp[2], direc[2], force[2])`. `aindex` are the two bead
indices; `disp` is the periodic-image displacement; `force = kl*(llen-l0)*direc`.

#### Force pipeline (`filament_ensemble::update()`, `filament_ensemble.cpp:468–485`)

```
for each filament f:
    update_filament_stretching(f)   # spring forces -> beads, fracture check
    f.update_bending(t)             # 3-body angle forces (LAMMPS-style)
    f.update_positions()            # Brownian-Langevin integration
update_energies()
t += dt
```

#### Stretch (`spring::update_force`, `spring.cpp:74–78`)

```cpp
double kf = kl * (llen - l0);
force = {{kf*direc[0], kf*direc[1]}};
```

Pure Hookean: `U_stretch = ½ kl (|r_{i+1} − r_i| − l0)²`. FENE and
Marko-Siggia variants exist (`update_force_fraenkel_fene`,
`update_force_marko_siggia`) but are **not** called by `update()` —
dead-ish helpers.

#### Bend (`filament::lammps_bending_update`, `filament.cpp:507–564`)

`U_bend = (κ_B / 2 l_0) Σᵢ θᵢ²` where θᵢ = π − (angle between consecutive
springs). Forces are derived as in LAMMPS `angle_harmonic.cpp`. Comment at
line 540 states "in this implementation, Lp = kb/kT" — i.e., the class
member `kb` equals **κ_B**, not the LAMMPS coefficient. Paper form:
`κ_B = L_p · k_B T`. With default `polymer_bending_modulus = 0.068 pN·μm²`
and kT = 0.004 pN·μm: **L_p = 17 μm** (matches actin literature).

Numerical guard: divides by `sin(θ)` with `s = max(s, maxSmallAngle)`
(`maxSmallAngle = 0.001`, line 534) — stability floor for nearly-straight
filaments.

#### Brownian-Langevin position update (`filament::update_positions`, `filament.cpp:200–230`)

```cpp
for each bead i:
    new_rnds = {rng_n(), rng_n()};                // two unit Gaussians
    vx = bead.force[0]/damp + bd_prefactor*(new_rnds[0] + prv_rnds[i][0]);
    vy = bead.force[1]/damp + bd_prefactor*(new_rnds[1] + prv_rnds[i][1]);
    prv_rnds[i] = new_rnds;                       // store for next step
    newpos = pos_bc(BC, delrx, dt, fov, {vx,vy}, {bead.x + vx*dt, bead.y + vy*dt});
    bead.set_xcm(newpos[0]); bead.set_ycm(newpos[1]);
    bead.reset_force();
```

with `bd_prefactor = sqrt(temperature / (2*dt*damp))`. This is the
**Leimkuhler-Matthews BAOAB-limit** for overdamped Langevin
(`filament.cpp:78` and `filament_ensemble.cpp:466` comments both reference
"Leimkuhler, 2013"). Concretely:

> r(t+Δt) = r(t) + (F/γ) · Δt + √(kT / (2 γ Δt)) · (W_n + W_{n−1}) · Δt
>         = r(t) + (F/damp) · Δt + sqrt(T/(2 Δt damp)) · (W_n + W_{n-1}) · Δt

Averaging two Gaussians is what makes this Leimkuhler-Matthews rather than
Euler-Maruyama. For harmonic systems it gives "minimal canonical
deviations" — important because actin springs are stiff.

#### Boundary conditions

`pos_bc()` in `globals.cpp` dispatches `PERIODIC` (default), `REFLECTIVE`,
`NONE`, `LEES-EDWARDS` (when `delrx ≠ 0` for shear). Minimum-image
distance via `rij_bc()`.

#### Fracture (`filament.cpp:267–288, 408–436`)

If any spring's stretch² > `fracture_force_sq`, the filament splits at
that node; new filaments are spawned from the bead vector. Default
`fracture_force = 10⁸ pN`, effectively disabled.

#### Subtleties

- `y_thresh` (`filament.cpp:206–211`) skips bead updates outside ±y_thresh
  of mid-y. Default 2 (never triggers). If set < 1, pinned-bead band near
  the box wall.
- `prv_rnds` is initialised to `{0,0}` for new beads (e.g. after
  fracture) — first integration step after fracture is effectively
  Euler-Maruyama, then L-M from step 2.
- `bead::update_force` (`bead.cpp:54–63`) **aborts the simulation** on
  NaN/inf force — useful invariant to replicate.

### 2.2 Paper parameters (Freedman 2017 Methods)

| Quantity | Symbol | Default | Notes |
| --- | --- | --- | --- |
| Bead radius | R | 0.5 μm | `actin_length` flag (misnomer; radius) |
| Spring rest length | l_link | 1.0 μm | `link_length` flag |
| Beads per filament | N_B | 11 | so filament contour ≈ 10·l_link = 10 μm |
| Bending modulus | κ_B | 0.068 pN·μm² | from L_p = 17 μm |
| Spring stiffness | k_a | 1 pN/μm (network sims); paper notes "70 pN/μm corresponds to actin" but used 1 to allow larger Δt | |
| Temperature | k_B T | 0.004 pN·μm | room temp |
| Viscosity | η | 0.001 mg/(μm·s) | water |
| Box | L_x × L_y | 50 × 50 μm | typical |
| Δt | | 2×10⁻⁵ s (network); 5×10⁻⁵ s (motility assay) | |
| BC | | PERIODIC + Lees-Edwards for shear | |

### 2.3 HOOMD port design — filament

| AFINES concept | HOOMD 7.0.1 mapping |
| --- | --- |
| `bead` (2D point with friction) | one HOOMD particle of type `actin`. HOOMD is 3D; constrain z by using a slab box with `box.Lz = 0` (HOOMD does not support strict 2D; 3D-with-thin-z is the convention) or pin z via a stiff `md.external.field`. |
| `vector<bead*> beads` per filament | a HOOMD particle group, with `body` tag = filament id (useful for analysis grouping). |
| `vector<spring*> springs` | `md.bond.Harmonic` with bond type `actin-actin`. HOOMD uses `U = ½ k (r − r0)²`; AFINES uses `U = ½ k_l (r − l_0)²`. Conventions match: `k_HOOMD = k_AFINES`, `r0 = l_link`. |
| 3-body bending (LAMMPS-style) | `md.angle.Harmonic` with `U = ½ k (θ − t0)²`. AFINES form is `(κ_B / 2 l_0) Σ θ²` (equilibrium at θ=0, i.e., straight; angle measured as deviation from π between consecutive bond vectors). HOOMD definition uses interior angle, so set `k_HOOMD = κ_B / l_0` and `t0 = π`. |
| Brownian-Langevin (L-M) integrator | `md.methods.Brownian(filter=actin_filter, kT=kT, default_gamma=6π η R)`. **HOOMD `Brownian` is Euler-Maruyama (O(Δt)), not BAOAB-limit (O(Δt²))**. Two options: (a) accept E-M and halve dt to compensate (~2× wall-clock), (b) write a custom HOOMD `Updater` (Python or C++) that stores `prv_rnds` per particle. Recommend (a) for Phase 1, defer (b) to performance tuning. |
| `damp = 6πη·R` | `default_gamma` in `Brownian` takes γ directly (NOT viscosity — convert). |
| `pos_bc()` periodic wrap | HOOMD wraps automatically when `Box` has periodic=(True, True, False). |
| Lees-Edwards shear | `hoomd.update.BoxResize` + custom shear updater. Exact API name in 7.0.1 needs verification when implementing H.1. |
| Fracture | not native; implement as a Python `Updater` that monitors bond forces and rebuilds the bond table. Expensive. AFINES default disables it; recommend the same. |

#### Concrete topology recipe (AFINES default config: npolymer=3, nmonomer=11)

- 33 particles, all type `actin`
- 30 bonds, type `actin-actin`, `k=1.0 pN/μm`, `r0=1.0 μm`
- 27 angles, type `actin-bend`, `k = κ_B/l_0 = 0.068 pN·μm`, `t0 = π`

#### Bead radius vs spring rest length

AFINES has them decoupled (`R = 0.5 μm`, `l_link = 1.0 μm`) — adjacent beads
do **not** touch. This is a coarse-graining choice with **no excluded
volume by default**. For dense networks (post-branching), enable
`md.pair.LJ` truncated to repulsive only at ε ~ 0.5 kT, σ ~ 2R.

#### Open questions

- **Q1.1** L-M custom plugin vs vanilla `md.methods.Brownian` (E-M).
  Recommendation: start with `Brownian` @ Δt = 1×10⁻⁵ s (half of AFINES);
  upgrade only if validation requires.
- **Q1.2** Excluded-volume? AFINES master ignores it; `dinner-group/AFINES`
  `exv` branch adds it. Required post-branching; defer for Phase 1 ECM.
- **Q1.3** N=21 beads/filament (v2 plan) or N=11 (AFINES default)? With
  L_f=10 μm and ℓ₀=0.5 μm bond rest length, N=21 is consistent (double
  AFINES bead density). Confirm intent.
- **Q1.4** Plan v2 `ℓ₀ = 0.5 μm` vs AFINES `R_bead = 0.5 μm` and
  `link_length = 1.0 μm` — disambiguate which is intended.

---

## 3. Crosslinkers (passive motors)

### 3.1 AFINES implementation

AFINES has **no separate xlink class**: both active motors and passive
crosslinkers are `motor` instances in two `motor_ensemble` containers
("amotors" / "pmotors"). A motor becomes an "xlink" by config:
`p_motor_v = 0` (no stepping), typically shorter rest length
(`p_motor_length = 0.150 μm` for filamin), lower stiffness.

Files: `include/motor.h` (lines 24–136), `src/motor.cpp` (581 lines,
the most complex AFINES file), `src/motor_ensemble.cpp` (269 lines).

#### Per-motor state (`motor.h:118–134`)

```
hx[2], hy[2]               # two head positions
disp[2], direc[2], len     # head-head displacement, normalised direction, length
force[2]                   # spring force on head 0 (head 1 gets -force)
tension                    # scalar = mk*(len − mld)
state[2]                   # head state: 0=unbound, 1=bound, -1=dead
f_index[2]                 # filament index per head (-1 if unbound)
l_index[2]                 # spring index within that filament
pos_a_end[2]               # arc-length from head to the j+1-th bead of its spring
mld, mk                    # rest length, stiffness
kon, koff, kend            # base attach / detach / barbed-end-detach rates (× dt)
vs, stall_force            # walking velocity at zero load, stall force
max_bind_dist, max_bind_dist_sq    # rcut; default a_m_cut = 0.063 μm
damp = 6π η mld            # head friction (uses motor length as the "size")
bd_prefactor = sqrt(T/(2·damp·dt))
prv_rnd_x[2], prv_rnd_y[2] # Leimkuhler memory
ldir_bind[2][2], bind_disp[2][2]   # filament orientation + offset at binding time
at_barbed_end[2]           # latched flag
```

#### Per-step pipeline (`motor_ensemble::motor_walk`, `motor_ensemble.cpp:159–200`)

```cpp
check_broken_filaments();
for each motor m:
    s = m.get_states();
    if s[0] == 0: m.brownian_relax(0);
    if s[1] == 0: m.brownian_relax(1);
    m.update_angle();
    m.update_force();
    m.filament_update();
    if s[0] == 0: m.attach(0); else m.step_onehead(0);
    if s[1] == 0: m.attach(1); else m.step_onehead(1);
update_energies();
```

#### Binding (`motor::attach`, `motor.cpp:249–302`)

For unbound head `hd`:

1. Query the neighbor-list quadrant containing `(hx[hd], hy[hd])` for all
   springs, sorted by distance² (`filament_ensemble::get_dist`).
2. For each spring (in order):
   - Reject if past `max_bind_dist_sq` (break, since sorted).
   - Reject if the other head is already bound to this same spring
     (`allowed_bind`).
   - Compute closest point on the spring (`spring::calc_intpoint` —
     projects + clamps to endpoints).
   - Compute Metropolis acceptance: `not_off_prob += metropolis_prob(hd,
     fl_idx, intPoint, kon)`.
   - If `mf_rand < not_off_prob`: snap head to `intPoint`, record indices,
     store `ldir_bind` and `bind_disp`, set `pos_a_end`. Return true.
3. If no spring accepted, return false.

#### Metropolis acceptance (`motor::metropolis_prob`, `motor.cpp:231–242`)

```cpp
stretch = dist(newpos, hx[other_head]) − mld;
delE    = 0.5 * mk * stretch * stretch − this->get_stretching_energy();
prob    = maxprob;                          // = kon*dt or koff*dt
if (delE > 0) prob *= exp(−delE/temperature);
return prob;
```

Glauber detailed balance: base rate `kon` (or `koff`), accept with
probability 1 if ΔU ≤ 0, else `exp(−ΔU/kT)`. **No Bell-Evans.** No
force-dependent k_off.

#### Off-position rotation (`motor::generate_off_pos`, `motor.cpp:379–393`)

Rotates `bind_disp` by the rotation that maps `ldir_bind` (filament
direction at binding time) onto the current filament direction; applies
the rotated offset to the current head position. This preserves detailed
balance even when filaments rotate while motors are bound — a subtle
correctness point.

#### Force back-propagation (`motor::filament_update_hd`, `motor.cpp:472–484`)

```cpp
pos_ratio = pos_a_end[hd] / spring_length;
filament.update_forces(f_index, l_index,   force * pos_ratio);
filament.update_forces(f_index, l_index+1, force * (1 − pos_ratio));
```

Standard Nedelec 2002 lever rule — partition motor force linearly between
the two endpoint beads of its current spring. Head 0 applies `+force`,
head 1 applies `−force` (Newton's third law for the motor's own spring).

#### Brownian relax for unbound heads (`motor::brownian_relax`, `motor.cpp:331–346`)

```cpp
vx = (−1)^hd * force[0] / damp + bd_prefactor*(new_rnd_x + prv_rnd_x[hd]);
vy = (−1)^hd * force[1] / damp + bd_prefactor*(new_rnd_y + prv_rnd_y[hd]);
pos = boundary_check(hd, hx[hd] + vx*dt, hy[hd] + vy*dt);
hx[hd] = pos[0]; hy[hd] = pos[1];
```

Sign `(−1)^hd`: head 0 gets +force, head 1 gets −force.

#### Subtleties

- `kon, koff, kend` are **pre-multiplied by dt** in the constructor
  (lines 49–51) — they are per-step probabilities by the time `event()`
  sees them.
- Quadrant rebuild dominates cost at scale.
- `allowed_bind` only checks the **spring**, not the filament — a single
  head may bind to the same filament its partner is on (just not the same
  spring). Self-binding artefact loophole.
- Cutoff `r_c = a_m_cut = p_m_cut = 0.063 μm` ≈ `√(kT/k_motor)` at
  k = 1 pN/μm — binding is restricted to a thin Boltzmann-shell around
  the rest length.

### 3.2 Paper parameters (Freedman 2017)

| Quantity | Crosslinker | Active motor | Notes |
| --- | --- | --- | --- |
| Length | 0.150 μm | 0.4 μm | filamin / myosin minifilament |
| Stiffness k | 1 pN/μm | 1 pN/μm | low for stability |
| k_on | 1 s⁻¹ | 1 s⁻¹ | base attach |
| k_off | 0.1 s⁻¹ | 0.1 s⁻¹ | base detach internal |
| k_end | 0.1 s⁻¹ | 1 s⁻¹ | barbed-end detach (10× for motors) |
| Stall F_s | 0 (no walking) | 0.5 pN | for `p_motor_v = 0` it's moot |
| v_0 | 0 μm/s | 1 μm/s | |
| Binding cutoff r_c | 0.063 μm | 0.063 μm | ≈ √(kT/k) |

### 3.3 HOOMD port design — crosslinkers

This is the most architecturally different subsystem from native HOOMD.
HOOMD has no concept of "dynamic bonds with stochastic create/destroy
events" out of the box — bonds are static topology, modified only by
custom plugins or `Snapshot` manipulation.

#### Design A (recommended — Python periodic updater)

```python
class XlinkUpdater(hoomd.custom.Action):
    def __init__(self, kon, koff, mld, k, kT, max_bind_dist):
        ...
    def act(self, timestep):
        with self._simulation.state.cpu_local_snapshot as snap:
            # 1. Unbinding pass: for each xlink particle bonded as harmonic-spring
            #    to actin, compute current stretch ΔU; Metropolis acceptance; if
            #    accepted, remove the bond from bond.group / bond.typeid arrays.
            # 2. Binding pass: for each unbound xlink head, query neighbor list
            #    (md.nlist on actin); for each candidate spring point compute
            #    Metropolis prob; on acceptance, add bond.
```

- **Pros**: clean Python expression matching AFINES one-to-one; debuggable.
- **Cons**: Python GIL; bond mutation requires snapshot upload (~µs per call).
- **Frequency**: every step at AFINES Δt is too slow. Batch every N steps
  where N·Δt ≪ 1/k_on. With AFINES Δt=2e-5s and k_on=1/s, per-step bind
  prob is 2e-5 — batching every 100 steps gives 0.002 events/xlink/check,
  fine.

#### Design B (C++ plugin)

Re-implement the Metropolis kernel as a HOOMD C++/CUDA plugin. Big effort,
only worth it once Design A is validated as the bottleneck. The
[glotzerlab/dybond](https://github.com/glotzerlab/dybond) plugin is a
reference for dynamic bonding in HOOMD; note it implements Bell-Evans, not
Metropolis.

#### Design C (continuous-bond, finite-lifetime approximation)

Always-bonded WCA pairs that effectively detach when stretch > r_c.
Fundamentally different physics; do not use unless the PI explicitly wants
a soft-bond model.

#### Particle topology for xlinks

- One particle per xlink head, type `xlink_head`, `body` tag = xlink id.
- One bond between paired heads, type `xlink_spring`,
  `k = p_motor_stiffness`, `r0 = p_motor_length`.
- When a head is "bound" to an actin spring, add an additional bond
  `xlink-actin` between the head and the chosen actin bead (or to a
  virtual point — but HOOMD bonds are particle-to-particle, so either
  bind to the nearer of the two endpoint beads (approximate) or insert a
  ghost particle on the spring constrained via `md.constrain.Rigid`
  (faithful, more complex).

#### Open questions

- **Q3.1** Bell-Evans vs Metropolis — **physically different**. Bell-Evans
  fits experimental single-molecule data better; Metropolis preserves
  detailed balance with the spring U. Pick one explicitly.
  Recommendation: Bell-Evans, document AFINES deviation in validation.
- **Q3.2** Faithful continuous-spring binding (ghost particles) vs
  snap-to-bead approximation. Worst-case force error of snap = ~0.5 pN
  with default k=1 pN/μm, l_link=1 μm.
- **Q3.3** Run xlink updater every step or batched? Batching by factor B
  introduces O(B·Δt·k_on) error in binding statistics — acceptable up
  to ~10⁻³.

---

## 4. Active motors (myosin minifilaments)

### 4.1 AFINES implementation

**Same `motor` class as xlinks**, differing only by config:

- `a_motor_v = 1 μm/s` (vs `p_motor_v = 0` for xlinks)
- `a_m_stall = 0.5 pN`
- `a_m_kend = 1 s⁻¹` (10× internal off-rate)
- `a_motor_length = 0.4 μm`

The real-myosin **bipolar minifilament** structure (~30 heads per side,
shared backbone, titin-like elasticity) is **not modeled**. AFINES motors
are two-headed "minifilament-as-spring": one Hookean spring between two
heads, each head independently binding/walking. No backbone, no
cross-bridge spring vs rigid-rod distinction.

#### Stepping kinetics (recap of §3 plus the stepping kernel)

Per step, for each bound head:

1. Compute proposed unbinding position (rotate `bind_disp` to current
   filament frame).
2. Compute off-probability via Metropolis with stretch-energy ΔU and base
   rate `kend` (if at barbed end) or `koff`, sample.
3. If not detached: advance `pos_a_end` by `vm · dt` where vm comes from
   the linear stall force-velocity, using the projection of the spring
   force onto the filament direction with sign `(−1)^hd`.
4. `pos_a_end` overflowing the current spring → transition to next spring
   toward barbed end (`update_pos_a_end`, `motor.cpp:424–455`).
5. Reaching the last spring latches `at_barbed_end = true` —
   subsequent stepping blocked; only detachment via `kend`.

#### Force-velocity (`globals::my_velocity`)

Standard piecewise-linear `v(F) = v0 · max{1 + F·r̂/F_s, 0}`, truncated
above at `2·v0`. At zero load `vm = vs`; at stall load `vm = 0`. **No
catch bond, no Hill function — just linear stall.**

#### Force back-propagation

Same lever rule as §3 — linear partition between the two endpoint beads
of the spring the head is currently on.

#### Subtleties

- Force-velocity uses **only** `(−1)^hd · dot(force, filament_direction)`
  — the force projected along the filament axis. Negative projection
  (assisting force) gives `vm > vs`, capped at `2·vs`.
- A head can step "off the end" of its current spring onto the next
  spring of the same filament. Reaching `l_index = 0` (barbed end)
  latches `at_barbed_end`; no mechanism to reset — escape only by
  detachment.
- Heads **cannot** step toward the pointed end. The check `if pos < 0`
  in `update_pos_a_end` is marked "shouldn't be possible if vm > 0"
  (line 440). If `vm < 0` were allowed (assisting force pushing the
  motor backwards), the code handles it, but it's not the typical regime.
- "Dead" motors (`dead_head_flag`): head `dead_head` (default 0) is
  killed via `kill_head` (`state = −1`); that head never participates.
  Models treadmilling-arrested or surface-tethered motors.

### 4.2 Paper parameters (Freedman 2017)

| Quantity | Value | Comment |
| --- | --- | --- |
| v_0 | 1 μm/s | non-muscle myosin II in vitro |
| F_s | 0.5 pN | per-head stall, small relative to 3–4 pN for cardiac myosin |
| k_on | 1 s⁻¹ | base |
| k_off | 0.1 s⁻¹ | base |
| k_end | 1 s⁻¹ | 10× boost at barbed end |
| Motor length | 0.4 μm | |
| Motor stiffness | 1 pN/μm | |

### 4.3 HOOMD port design — motors

Architecturally **identical** to xlinks (§3.3), with one additional kernel:
the **stepping update** (advance bound head along the filament).

```python
def act(self, timestep):
    super().act(timestep)   # do the Metropolis bind/unbind from xlink logic
    for motor in self.motors:
        for hd in (0, 1):
            if motor.state[hd] == 1 and not motor.at_barbed_end[hd]:
                f_dot_dir = dot(motor.force * (-1)**hd, motor.fil_direction[hd])
                vm = max(self.vs * (1 + f_dot_dir / self.Fs), 0.0)
                vm = min(vm, 2 * self.vs)
                motor.pos_a_end[hd] += vm * dt
                if motor.pos_a_end[hd] >= motor.spring_length[hd]:
                    if motor.l_index[hd] == 0:
                        motor.at_barbed_end[hd] = True
                        motor.pos_a_end[hd] = motor.spring_length[hd]
                    else:
                        motor.l_index[hd] -= 1
                        motor.pos_a_end[hd] -= motor.spring_length[hd]
                # update bonded-actin-bead index in HOOMD topology
```

Expensive part is recomputing `f_dot_dir` (needs current filament
tangent). With N_motors ~10⁴ and one call per step, ~10⁴ × 3 ops/step in
Python = ~ms/step, tolerable for AFINES-scale runs. For 10⁵ motors at
large scale, write a C++/CUDA plugin.

#### Minifilament representation choice

If v2 wants a more realistic minifilament than AFINES single-spring:

- **Option A (AFINES-equivalent)**: one 2-head bipolar bridge per motor.
  Simple, matches AFINES paper-for-paper.
- **Option B (Stam-Hocky 2017 PRL style)**: a chain of beads making up
  the rod (~14 beads), each with two "heads" reaching out via auxiliary
  springs. ~30× more particles per motor; matches published minifilament
  models from the Hocky group.
- **Option C (rigid-body backbone)**: rigid body for the bipolar rod via
  `md.constrain.Rigid`, with multiple cross-bridge springs hanging off.

#### Open questions

- **Q4.1** Single-spring motor (AFINES) or multi-head minifilament (Stam-Hocky)?
  Former matches AFINES paper-for-paper; latter is what cellular
  biophysicists usually want for contractility.
- **Q4.2** Catch-bond detail? AFINES has no catch bond. Cardiac myosin
  exhibits catch behaviour (k_off ↓ with load). v2 says "Bell-Evans"
  which is slip. Confirm slip vs catch.
- **Q4.3** Hill force-velocity (`v = v0·(F_s − F)/(F_s + F/a)`) vs linear
  stall (AFINES). Hill is more biologically accurate but adds a free
  parameter `a`.
- **Q4.4** What happens when the motor reaches the barbed end? AFINES
  latches and only detaches via `kend`. v2 may want explicit
  polymerisation (treadmilling) or filament-end interactions.

---

## 5. Arp2/3 dendritic branching — **greenfield (no AFINES reference)**

### 5.1 AFINES status

**There is no Arp2/3 branching code in AFINES.** Searching
`Simfreed/AFINES`, all 9 forks, and all 28 dinner-group branches for
`arp`, `branch`, `nucleat` returns **0 hits in source**.

The closest existing feature is the `dinner-group/AFINES@growing` branch's
`filament::grow(double dL)` (`filament.cpp:705–761`), which performs
**linear barbed-end elongation only**: extends `links[0].l0` and inserts
a new bead between the two leading beads when `l0` exceeds `l0_max`. No
lateral branch nucleation.

Plan v2 Unit H.5 ("AFINES branched mode port") therefore needs
**greenfield design** from experimental phenomenology, not a port. Adds
~2 weeks scope vs the original estimate.

### 5.2 Bieling 2016 + Funk 2021 + Li/Bieling 2022 phenomenology (what HOOMD must reproduce)

From [Bieling et al. 2016 Cell](https://doi.org/10.1016/j.cell.2015.11.057)
[PMC5033619] + [Li, Bieling et al. 2022 eLife](https://doi.org/10.7554/eLife.73145)
[PMC9328761]:

#### Single-molecule nucleation

- WAVE1ΔN per-molecule rate at zero force: **k_b⁰ = 0.037 s⁻¹**.
- Decreases by ~20 % from zero to stall load (~1200 Pa).
- Functional form: **NOT exponential**; per-network nucleation falls
  **linearly** with load to ~50 % at stall. The discrepancy (single
  ~20 % vs bulk ~50 %) is explained by "abortive nucleation" — Arp2/3
  complexes that bind and rapidly fall off under high load.
- Qualitative fit: `R_nucleate(F) ≈ R_0 · (1 − α F)` with
  `α ≈ 0.4 / F_stall`, plus an abortive-failure population above ~500 Pa.

#### Per-filament rates (Bieling Fig. 5C)

Both elongation and capping are well-fit by **single exponential decay in
force per filament**:

- `k_elong(f) = k_elong⁰ · exp(−f · δ_elong / kT)` — Brownian ratchet.
- `k_cap(f)   = k_cap⁰   · exp(−f · δ_cap · sin θ / kT)` — Li/Bieling 2022.

Funk fits: `δ_cap ≈ 0.3 pN` (tethering force), `sin θ` accounts for
barbed-end geometry. Zero-force capping: `k_cap⁰ ≈ 3 s⁻¹ · (100 nM CP)⁻¹`.

#### Density / architecture

- Filament density `n(F)` increases ~8× across functional load range; not
  a power law (saturates near stall).
- Free barbed-end density: ~160/μm² (low F) → ~550/μm² (high F).
- Mean filament length **invariant**: ~300 nm (~110 monomers) — because
  elongation and capping respond identically to per-filament force.
- Mean attack angle: ~54° (low F, Y-branch geometry) → ~16° (high F,
  compression).
- Elastic modulus ∝ n^0.6 (weaker than n² for isotropic gels) — branched
  anisotropic networks have different physics.

#### Mechanism (Funk 2022 model)

- Branching occurs from the side of an existing filament near its
  barbed end.
- Arp2/3 binds the mother filament; NPF (WAVE/WASp) at the membrane
  recruits actin-profilin; daughter initiates at ~72° from mother.
- Force on the network slows mother elongation → daughter Arp2/3 binding
  more probable per unit time (autocatalytic). Force also reduces NPF
  availability via "barbed-end interference" (free barbed ends near
  membrane sequester WH2 domains needed for monomer recruitment) → the
  ~20 % per-Arp2/3 nucleation drop.
- Net: density UP because capping slows MORE than nucleation slows.

### 5.3 HOOMD port design — branching

Greenfield. Recommended sketch:

#### Particle topology additions

- New type `arp23` — small particle representing the branch point.
- New bond type `arp-mother` (arp23 ↔ actin bead at branch point).
- New bond type `arp-daughter` (arp23 ↔ first actin bead of daughter).
- New angle constraint `arp-daughter-actin0-actin1` with t0 = 72° ≈ 1.257 rad.

#### Branching kernel (Python Updater, runs every N steps)

```python
class BranchUpdater(hoomd.custom.Action):
    def __init__(self,
                 k_b_zero=0.037, F_per_filament_scale=2.0,  # pN
                 branch_angle=72*pi/180, daughter_length_init=0.5,
                 npf_zone_thickness=0.5):  # μm — only branch within this distance of "membrane"
        ...
    def act(self, timestep):
        # 1. Identify candidate mother filaments — barbed end within
        #    npf_zone_thickness of the membrane surface (e.g. +y wall).
        # 2. Compute force-per-filament: total reaction force at the WAVE
        #    patch normal / # free barbed ends. Simpler proxy: load
        #    force / N_growing.
        # 3. Per-mother branching rate:
        #    k_b(F) = k_b_zero * (1 - 0.2 * F/F_stall)        # Bieling
        #    p_branch_per_step = k_b * N * dt
        # 4. For each candidate firing the Bernoulli:
        #    - Pick a random position near the mother barbed end (within 0.5 μm).
        #    - Spawn an arp23 particle there.
        #    - Spawn a daughter filament (5-7 actin beads); first bead bonded
        #      to arp23, second bead constrained at 72°.
        #    - Mutate snapshot to add particles + bonds + angle.
        # 5. Capping: per-step, for each free barbed end, sample with
        #    k_cap(F) = k_cap_zero * exp(-F * delta / kT), delta = 0.3 pN.
        #    Set is_capped flag; treadmilling skips capped ends.
```

#### Elongation kernel

Independent of branching — extend the barbed-end spring of each
*non-capped* filament at rate `k_elong(F)` with Bell-Evans form. The
`dinner-group/AFINES@growing` branch has the bead-insertion logic to
port.

#### Subtleties to anticipate

- **Bookkeeping cost**: branching + capping grows N_filaments linearly
  until steady state. Snapshot mutation is expensive (~ms). Batch
  topology updates every ~10³ steps.
- **"Membrane" geometry**: AFINES has no membrane. Decide: (a) flat
  surface at `y = Y_max` with WAVE density σ_NPF, (b) Gaussian zone,
  (c) explicit polymer-coated bead surface. (a) simplest.
- **Per-filament force partition**: load divided by free barbed ends in
  the WAVE zone. Count + partition each step.
- **Volume exclusion** becomes critical as density grows 8×. Recommend
  enabling `md.pair.LJ` repulsive-only between all actin beads at
  moderate ε ~ 0.5 kT.

#### Open questions

- **Q5.1** Branch geometry: rigid 72° Y-branch (Arp2/3 crystal) or 70°±5°
  angular fluctuation? Different papers use different prescriptions.
- **Q5.2** k_b functional form: Bieling network-level `(1 − 0.4 F/F_stall)`
  OR per-Arp2/3 `(1 − 0.2 F/F_stall)` + abortive failure? Latter is
  mechanistically richer.
- **Q5.3** Capping: Bell-Evans `k_cap = k_0 · exp(−F·0.3/kT)` per Funk
  2022, or constant `k_cap` per AFINES no-capping behaviour? v2 doesn't
  mention capping, but Bieling phenomenology requires it.
- **Q5.4** Treadmilling at pointed ends? Bieling assumes no pointed-end
  dynamics (Arp2/3 protects); fine for second-scale sims.
- **Q5.5** Daughter initial length: 1 bead (and elongate) or seed at
  ~3 beads (1.5 μm) to avoid instant capping? Recommend 3.
- **Q5.6** Debranching? AFINES doesn't model it; real `k_debranch ~ 10⁻⁵
  s⁻¹` negligible for second-scale sims. Defer.

---

## 6. Overdamped Brownian integrator

### 6.1 AFINES implementation

- **Scheme**: Leimkuhler-Matthews BAOAB-limit. Labeled at
  `filament_ensemble.cpp:466` ("Overdamped Langevin Dynamics Integrator
  (Leimkuhler, 2013)") and stated in Freedman 2017 Methods.
- **Update**:
  - `v = F/γ + sqrt(kT / (2γΔt)) · (W_n + W_{n−1})`
  - `r(t+Δt) = r(t) + v Δt`
  - `W_n, W_{n−1}` independent unit Gaussians; previous draw stored as
    `prv_rnds`.
- **Hardware**: CPU only, serial. `#pragma omp parallel for` commented
  out in all hot loops.
- **Δt**: 2×10⁻⁵ s for network sims, 5×10⁻⁵ s for motility assay.
  Springs at k=1 pN/μm with γ=6π·η·R = 6π·0.001·0.5 ≈ 0.0094 mg/s give
  τ = γ/k ≈ 0.01 s; Δt < τ/100 needed for accuracy.
- **Reference**: Leimkuhler & Matthews, *J. Chem. Phys.* 138, 174102
  (2013), [doi:10.1063/1.4802990](https://doi.org/10.1063/1.4802990).

### 6.2 HOOMD port design — integrator

- `hoomd.md.methods.Brownian(filter=All(), kT=variant_kT, default_gamma=γ)`
  — Euler-Maruyama, O(Δt).
- `hoomd.md.methods.Langevin` is **full Langevin (underdamped)** — not the
  AFINES overdamped scheme.
- **Δt**: with E-M O(Δt) vs L-M O(Δt²) for harmonic systems, halve Δt as
  starting point — use **Δt = 1×10⁻⁵ s**.
- **Custom L-M plugin**: if benchmarking shows Δt=1×10⁻⁵ is the
  bottleneck, write a HOOMD `Updater` (Python or C++) that does:
  1. Compute forces.
  2. Per particle: draw `W_n`, read `prv_rnds`, compute
     `v = F/γ + bd_pre·(W_n + W_{n−1})`, integrate `r += v·Δt`, store
     `W_n` as `prv_rnds`.
- **GPU benefit**: HOOMD `Brownian` is a CUDA kernel — expect ~10–100×
  speedup vs AFINES single-thread CPU for >10⁴ particles. On M1 Max
  CPU, the polymer sanity benchmark (§7) shows we are at ~90k steps/s
  for 100 beads; HOOMD scales sub-linearly with N up to GPU saturation.

#### Open questions

- **Q6.1** Accept HOOMD `Brownian` Euler-Maruyama with Δt = 1×10⁻⁵ s, or
  implement L-M to match AFINES at Δt = 2×10⁻⁵ s? Plugin cost ~1–2 days;
  Δt halving costs ~2× wall-clock per simulated second.
- **Q6.2** GPU vs CPU: HOOMD supports both. For typical ActiveCellSim
  scale (10⁴ beads), GPU recommended once available.
- **Q6.3** Thermostat consistency between actin beads and motor heads:
  AFINES uses the same `bd_prefactor` form for both. In HOOMD, both
  particle types must use `md.methods.Brownian` with the same kT but
  different `default_gamma` per type. Motor head γ_motor = 6π·η·motor_length
  (NOT bead_radius). Configure per-type gamma table.

---

## 7. Per-subsystem effort + integration matrix

| Subsystem | AFINES code | HOOMD primitives available | Port effort |
| --- | --- | --- | --- |
| Filament bead-spring | `filament.cpp`, `bead.cpp`, `spring.cpp` | `md.bond.Harmonic` + `md.angle.Harmonic` | Low — 1 day |
| Overdamped integrator | `filament::update_positions` (L-M) | `md.methods.Brownian` (E-M) | Low if accept E-M; medium (1–2 d) for custom L-M plugin |
| Crosslinker dynamics | `motor.cpp` (Metropolis bind/unbind, no walking) | none — Python custom Action | High — 1–2 weeks |
| Motor stepping | `motor.cpp` (linear stall F-V + Metropolis off) | none — Python custom Action | High — 1–2 weeks (extends xlink) |
| Arp2/3 branching | **none in AFINES** | none — full design from Bieling/Funk | Highest — 2–3 weeks |
| Filament fracture | `filament::fracture` | none — Python | Medium; usually disabled |
| Excluded volume | `dinner-group@exv` (not master) | `md.pair.LJ` repulsive-only | Low |
| Shear / Lees-Edwards | `filament_ensemble::update_d_strain`, `pos_bc` | `BoxResize` / custom box updater | Medium |
| Filament growth (treadmilling) | `dinner-group@growing` | Python Action | Medium |

Total Phase 1 effort estimate (Plan v2 said 3 months for Phase 1 full;
remains consistent with the above table assuming Workers A/B/C in
parallel).

---

## 8. Consolidated open questions for PI

Listed in priority order (those affecting Plan-scope first):

1. **Branching mechanism not in AFINES** (§5, all of Q5). Plan H.5 must
   be greenfield from Bieling/Funk phenomenology — needs sign-off on
   `k_b` form, capping inclusion, branch geometry.
2. **Bell-Evans vs Metropolis** for xlink/motor off-rate (Q3.1, Q4.2).
   v2 says Bell-Evans; AFINES uses Metropolis. **Physically different
   models**; pick one. Recommendation: Bell-Evans.
3. **Integrator: L-M custom plugin vs vanilla Brownian** (Q1.1, Q6.1).
   2× Δt cost penalty for vanilla; ~1–2 days dev for plugin.
4. **Bead resolution N=11 (AFINES) vs N=21 (v2)** + `ℓ₀=0.5 μm` vs
   `link_length=1.0 μm` clarification (Q1.3, Q1.4). v2 numbers may
   confuse AFINES `link_length` and `actin_length` (= bead radius).
5. **Motor minifilament representation** (Q4.1). Single 2-head spring
   (AFINES) vs multi-head bipolar (Stam-Hocky). Important for
   contractility realism.
6. **Force-velocity form** (Q4.3). Linear stall (AFINES) vs Hill.
   Probably stick with AFINES form for parity.
7. **Excluded volume between filaments** (Q1.2). AFINES master ignores;
   matters for dense post-branching networks.

---

## 9. References

### Primary code

- Simon Freedman et al., AFINES repository (canonical),
  <https://github.com/Simfreed/AFINES>. master @ 2025-03-31. GPLv3.
- Dinner-lab AFINES fork (extended features),
  <https://github.com/dinner-group/AFINES>. Branches: `master`, `growing`,
  `polymerize`, `exv`, `growing_and_spacer`, `sorting_guv`, ...
  Last pushed 2021-09-20.

### Primary papers

- Freedman, Banerjee, Hocky, Dinner. "A versatile framework for simulating
  the dynamic mechanical structure of cytoskeletal networks." *Biophysical
  Journal* 113, 448–460 (2017). doi:[10.1016/j.bpj.2017.06.003](https://doi.org/10.1016/j.bpj.2017.06.003).
  [PMC5529201](https://pmc.ncbi.nlm.nih.gov/articles/PMC5529201/).
- Bieling, Li, Risca, Sun, Fletcher, Mullins. "Force feedback controls
  motor activity and mechanical properties of self-assembling branched
  actin networks." *Cell* 164, 115–127 (2016).
  doi:[10.1016/j.cell.2015.11.057](https://doi.org/10.1016/j.cell.2015.11.057).
  [PMC5033619](https://pmc.ncbi.nlm.nih.gov/articles/PMC5033619/).
- Funk, Bieling, Mueller, Mullins, Fletcher. "The molecular mechanism of
  load adaptation by branched actin networks." *eLife* 11, e73145 (2022).
  doi:[10.7554/eLife.73145](https://doi.org/10.7554/eLife.73145).
  [PMC9328761](https://pmc.ncbi.nlm.nih.gov/articles/PMC9328761/).

### Algorithmic references

- Leimkuhler, Matthews. "Robust and efficient configurational molecular
  sampling via Langevin dynamics." *J. Chem. Phys.* 138, 174102 (2013).
  doi:[10.1063/1.4802990](https://doi.org/10.1063/1.4802990).
- Nedelec, F. "Computer simulations reveal motor properties generating
  stable antiparallel microtubule interactions." *J. Cell Biol.* 158,
  1005–1015 (2002).
- Bell, G. "Models for the specific adhesion of cells to cells." *Science*
  200, 618–627 (1978). (Bell-Evans formalism foundation.)

### Adjacent / context

- Chandrasekaran, Edelmaier, Voth, Hocky. "Toward the cellular-scale
  simulation of motor-driven cytoskeletal assemblies." *eLife* 11, e74160
  (2022). doi:[10.7554/eLife.74160](https://doi.org/10.7554/eLife.74160).
  aLENS framework: independent post-AFINES re-implementation in modular
  C++ with OpenMP+MPI, no GPU, no HOOMD. Useful comparison point if we
  ever consider an aLENS+HOOMD hybrid.
- HOOMD-blue 7.0.1 documentation, <https://hoomd-blue.readthedocs.io/>.
  `hoomd.md.methods.Brownian` API.

---

## Appendix A — source files inspected (verbatim list)

`include/`: `filament.h`, `filament_ensemble.h`, `motor.h`,
`motor_ensemble.h`, `spring.h`, `bead.h`, `globals.h`.

`src/`: `filament.cpp`, `filament_ensemble.cpp`, `motor.cpp`,
`motor_ensemble.cpp`, `spring.cpp`, `bead.cpp`.

`prog/`: `network.cpp`, `motors_on_struct.cpp`.

`misc/`: `makefile`, `README.md`, example `.cfg`.

Dinner-group fork (excerpts): `motor.cpp`, `exv.cpp`, `quadrants.cpp`,
`bundles.cpp`, and `growing` branch `filament.cpp` for treadmilling logic.
