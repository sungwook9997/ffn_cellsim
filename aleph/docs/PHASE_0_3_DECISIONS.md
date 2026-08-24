# Phase 0.3 — PI decisions ratified (2026-05-19)

**Principle (PI 2026-05-19)**: "full resolution simulation 방향에 맞게
너의 옵션을 믿음. 처음 프레임워크에 추상화나 다른 과정이 들어가면 안됨."

Translation: pick the **most physically faithful, fine-grained,
mechanistic** option for every Phase 1 design decision, even at higher
implementation or wall-time cost. The whole point of v2 is to build
behaviours emergently from particle/bond dynamics rather than wrap
literature closed-forms as runtime mechanism (the v1 paradigm).

**Exception**: the ×40 mesoscopic coarse-graining (1000 effective filaments
per cell instead of ~38,000 native) is an explicit Plan v2 §3 H.3
hardware constraint already ratified in the 2026-05-19 v3.1 revision.
Within that fixed scale, every other decision goes to full-fidelity.

---

## D1 — Arp2/3 dendritic branching (H.5)

**Decision**: greenfield implementation (no AFINES code to port). Per-Arp2/3
mechanistic form with thermal angle fluctuation and Bell-Evans capping.

**Concrete spec**:

| Component | Implementation |
| --- | --- |
| Per-Arp2/3 nucleation rate | `k_b(F) = k_b⁰ · (1 − 0.2 · F / F_stall)` with abortive-failure population above ~500 Pa (Funk 2021 mechanism); k_b⁰ = 0.037 s⁻¹ per WAVE molecule (Bieling 2016). **Not** the lumped network-level `(1 − 0.4 F/F_stall)` form. |
| Branch angle | `md.angle.Harmonic` with `t0 = 72°` (Arp2/3 crystal structure) and `k_angle` from Arp2/3 binding stiffness (≈ 100 pN·μm/rad²; literature TBD). Thermal fluctuation emerges naturally — no manual `±5°` jitter. |
| Capping | Bell-Evans: `k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sinθ / kT)`, δ_cap = 0.3 pN (Li/Bieling 2022). k_cap⁰ ≈ 3 s⁻¹ at 100 nM CP. NOT constant rate. |
| Daughter seed | 1 actin bead at nucleation, elongates by Bell-Evans `k_elong(F) = k_elong⁰ · exp(−F · δ_elong / kT)`. NOT pre-seeded at 3 beads. Accept the early-capping rate; that is the mechanism. |
| Force partition | Network load divided by N_free_barbed_ends in the WAVE zone, per timestep. |
| Membrane geometry | Flat surface at `y = Y_max` with WAVE density σ_NPF (simplest). |

**Replaces** `acs_kb/ [v1, deleted] /cell/lamellipodia.py` (archived). HOOMD updaters live
in `aleph/cell/lamellipodium.py`.

---

## D2 — Bond off-rates (xlinks / motors / junctions / FA clutch)

**Decision**: Bell-Evans force-dependent off-rates per bond type, with
catch-bond formulations where the experimental literature supports them.
**Not** AFINES Metropolis (Glauber detailed-balance is a detailed-balance
convenience, not the experimentally-measured mechanism).

**Concrete spec**:

| Bond type | Form | Parameters | Source |
| --- | --- | --- | --- |
| Filamin / α-actinin xlink | Bell-Evans slip: `k_off(F) = k_0 · exp(F·x_β/kT)` | `k_0 = 0.1 s⁻¹`, `x_β ≈ 0.3 nm` (filamin literature) | Furuike 2001, Ferrer 2008 |
| Myosin motor head | Bell-Evans slip (NMII) | `k_0 = 0.1 s⁻¹` internal, `k_end = 1 s⁻¹` barbed-end | AFINES defaults preserved |
| E-cadherin trans-bond | **Full catch-bond (KU-4.2 Buckley 2014)**, NOT slip-only (KU-4.17). Two-pathway: slip exponent + catch exponent. | Buckley 2014 Fig 4 fit | Buckley 2014 Science |
| Integrin (FA clutch) | Catch-slip Pereverzev (KU-2.5, KU-2.18) | `k_s = 0.5 s⁻¹, F_s = 30 pN, k_c = 0.4 s⁻¹, F_c = 7 pN` | Pereverzev 2005, Kong 2009 |

**Cadherin Phase 1/2 split override**: Plan v2 originally placed KU-4.2 full
catch-bond as a Phase 2 follow-up, with KU-4.17 slip-only as the Phase 1
baseline. **Override per "no abstractions"**: full catch-bond from Phase 1.
The slip-only formulation is exactly the kind of v1 simplification PI is
warning against.

**Implementation**: Python `hoomd.custom.Action` per bond family,
batched every ~100 steps (`100·Δt·k_max ≈ 10⁻³` events per bond per
batch — acceptable error vs single-step evaluation). Bond
addition/removal via `snapshot.bonds.group` mutation.

---

## D3 — Overdamped integrator

**Decision**: Custom HOOMD `Updater` implementing Leimkuhler-Matthews
BAOAB-limit overdamped Langevin (matches AFINES, O(Δt²) error on harmonic
systems). **Not** vanilla `md.methods.Brownian` (Euler-Maruyama, O(Δt)).

**Concrete spec**:

```python
class LeimkuhlerMatthewsBAOAB(hoomd.custom.Action):
    """r(t+Δt) = r(t) + (F/γ)·Δt + sqrt(kT/(2γΔt)) · (W_n + W_{n-1}) · Δt
    
    Stores per-particle prv_rnds (W_{n-1}) between calls.
    Matches AFINES filament::update_positions (filament.cpp:200-230).
    """
```

**Δt**: 2×10⁻⁵ s (AFINES default), not halved.

**Effort**: 1-2 days dev, then verify against the §7 polymer sanity bench
to confirm equipartition convergence rate matches AFINES theory.

**Lives in**: `aleph/integrator/baoab.py`.

---

## D4 — Filament discretisation (bead resolution)

**Decision**: Plan v2 spec (N=21 beads/filament, ℓ₀=0.5 μm) confirmed as
**bond rest length**. AFINES `link_length` (1.0 μm) and `actin_length`
(bead radius, 0.5 μm) are decoupled; v2 doubles the bond density.

**Concrete spec**:

| Parameter | Value | Note |
| --- | --- | --- |
| Beads per filament | **N = 21** | doubles AFINES density |
| Bond rest length `ℓ₀` | **0.5 μm** | bond `r0` in `md.bond.Harmonic` |
| Contour length `L_f` | **10 μm** | = (N−1) · ℓ₀ |
| Bead radius (actin bundle) | **R = 30 nm** | ×40 mesoscopic bundle cross-section |
| Bead radius (collagen) | **R = 50 nm** | KU-1.2 |
| Bead radius (decoupled from ℓ₀) | ✓ | AFINES convention |

Per-bead Stokes drag `γ_b = 6π·η·R` uses the bead radius, not ℓ₀.

---

## D5 — Motor minifilament representation

**Decision**: Stam-Hocky bipolar multi-head minifilament. **Not** AFINES
single 2-head Hookean spring (that is the abstraction PI is warning
against).

**Concrete spec**:

| Component | Value |
| --- | --- |
| Backbone | rigid-rod approximation via `md.constrain.Rigid`, ~14 beads, length ≈ 700 nm |
| Heads per side | ≈ 10 (Hocky-group convention; real NMII ~30 but ~10 is computational standard) |
| Cross-bridge spring | Hookean, `k = 1 pN/μm` (AFINES motor stiffness), `r0 = 200 nm` (head-to-backbone) |
| Per-head dynamics | independent Bell-Evans bond + Hill stepping (D2 + D6) |
| Architecture | bipolar — 10 heads each side, opposite polarity |

**Plan v2 §3 H.3 says "100 myosin minifilaments per cell"**. With ~14
backbone beads + 20 heads per minifilament = 34 particles per minifilament
× 100 = **3,400 extra particles per cell** for motors. Within Plan §11
budget (~10K beads/cell Phase 1 total).

**Lives in**: `aleph/cortex/myosin.py`.

---

## D6 — Motor force-velocity

**Decision**: Hill form `v(F) = v0 · (F_s − F) / (F_s + F/a)` with
`a/F_s = 0.5` (Kovács 2003 for non-muscle myosin II). **Not** AFINES
piecewise-linear stall.

**Rationale**: Hill 1938 is the canonical biophysical model fit to muscle
force-velocity data; AFINES piecewise-linear is a computational
convenience. For non-muscle myosin specifically, Kovács 2003 J Biol Chem
fits give `a ≈ 0.5 F_s`.

**Concrete spec**:

```python
def hill_velocity(F: float, v0: float, F_s: float, a: float) -> float:
    """Hill (1938) force-velocity for myosin head stepping.
    
    F: load force (pN), F_s: stall force (pN), v0: unloaded velocity (μm/s),
    a: Hill shape parameter (pN), default 0.5·F_s for NMII (Kovács 2003).
    Returns step velocity (μm/s), clamped to [0, 2·v0] for assisting loads.
    """
    if F >= F_s:
        return 0.0
    if F <= -F_s:
        return 2.0 * v0  # assisting-load cap
    return v0 * (F_s - F) / (F_s + F / a)
```

Parameters per minifilament head: `v0 = 1 μm/s, F_s = 0.5 pN, a = 0.25 pN`.

---

## D7 — Excluded volume

**Decision**: `md.pair.LJ` repulsive-only (WCA-style) **ON from Phase 1**.
Not deferred to post-branching. Excluded volume is a physical fact, not
an optional feature.

**Concrete spec**:

```python
lj = md.pair.LJ(nlist=nlist, default_r_cut=0)  # disabled by default per-pair
# Enable per pair-type, truncate to repulsive only via WCA convention
lj.params[('actin', 'actin')] = dict(epsilon=kT * 0.5, sigma=2 * R_actin)
lj.r_cut[('actin', 'actin')] = 2**(1/6) * 2 * R_actin   # WCA cutoff
# Repeat for arp23, xlink_head, motor_head pair types
lj.mode = 'shift'
```

WCA cutoff `2^(1/6)·σ` truncates LJ at its minimum, keeping only the
repulsive branch — purely steric, no attractive tail. ε ~ 0.5 kT is the
standard "soft" repulsion that avoids overlap without producing
artificial crystallisation.

Per-pair-type table needed for: actin-actin, actin-arp23, actin-xlink,
actin-motor_head, motor_head-motor_head. Phase 1 default is to enable all
intra-cell pairs.

**Cost**: HOOMD `md.pair.LJ` is GPU/CPU-optimised; expect ≤ 20% wall-time
overhead at Phase 1 density.

---

## Summary table

| # | Decision | v1 / AFINES default | **v2 full-fidelity choice** |
| --- | --- | --- | --- |
| 1 | Branching | none in AFINES | per-Arp2/3 `(1 − 0.2 F/F_s)` + abortive + Li/Bieling 2022 capping + 1-bead daughter seed |
| 2 | Bond off-rates | Metropolis (AFINES); slip-only cadherin (Plan Phase 1) | **Bell-Evans per bond type; full catch-bond cadherin from Phase 1** |
| 3 | Integrator | E-M (HOOMD vanilla) | **Leimkuhler-Matthews BAOAB-limit (custom Updater)** |
| 4 | Bead resolution | N=11, l_link=1.0 μm (AFINES) | **N=21, ℓ₀=0.5 μm** (Plan v2 confirmed) |
| 5 | Motor architecture | single 2-head spring (AFINES) | **Stam-Hocky bipolar multi-head minifilament** |
| 6 | Force-velocity | piecewise-linear stall (AFINES) | **Hill 1938** with `a/F_s = 0.5` (Kovács 2003 NMII) |
| 7 | Excluded volume | off (AFINES master) | **`md.pair.LJ` WCA repulsive-only ON from Phase 1** |

## Consequences for the dispatch plan

The 7 decisions feed back into Phase 1 Worker dispatch (Plan v2 §3):

- **Worker A (H.1 ECM Mikado)**: D3 (integrator), D4 (bead resolution), D7 (EV). +2-3 days for L-M plugin + EV table — within H.1's 2-week budget.
- **Worker B (H.4 FA + motor-clutch)**: D2 (Bell-Evans + Pereverzev for integrin), D5 (multi-head minifilament for motor side), D6 (Hill). +1 week vs AFINES single-spring port; H.4 budget extends from 2 → 3 weeks.
- **Worker C (H.2 single filament)**: D3 + D4 only. Same 1-week budget.
- **Worker C (H.3 cortex)**: D2 (filamin Bell-Evans), D5 (motor architecture), D7. Same 3-week budget.
- **Worker C (H.5 lamellipodium)**: D1 (greenfield), D7. Plan said 2 weeks for "AFINES branched mode port"; greenfield is **4 weeks** because the design itself is new.
- **Worker D (H.8 junction Phase 2)**: D2 cadherin full catch-bond from Phase 1 means the Phase 2 work shifts from "add catch detail to slip-only baseline" to "validate the always-catch implementation against Maître two-cell rounding".

**Net Phase 1 schedule impact**: +2 to +3 weeks (mostly H.5 greenfield).
Phase 1 total: ~3 months → ~3.5-4 months.

## What was NOT decided (still open)

- Integrator GPU vs CPU choice (Plan v2 §0: M1 Max CPU through Phase 1, GPU at Phase 2 entry; no change).
- Notion `aleph/configs/` YAML schemas — will be drafted as part of H.1/H.3/H.5 Worker dispatch briefs.
- `aleph/validation/` oracle library structure — will be drafted as v2 Bell-Evans/Hill/Buckley validation tests come online.
