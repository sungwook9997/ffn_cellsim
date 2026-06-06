# H.7 multiscale active-gel γ-seam — design + Sanity Gate

> **Status**: DESIGN (proposal, NOT ratified). PI-chosen path 2026-06-07 for the
> active cortical-tension floor. Awaiting PI sign-off before implementation.
> **Author**: H.7 DESIGN subagent, 2026-06-07.
> **Supersedes nothing**; it is the structural answer to the timescale gap
> diagnosed in `H7_GATEB_FINDINGS_2026-06-07.md` and
> `H7_OVERNIGHT_SUMMARY_2026-06-07.md`.
> **Scope**: design ONE seam module. No code is written here; no git is run.

---

## 0. The problem this seam solves (one paragraph)

The full physiological MCF7 cell assembles and runs at full ×40 scale, but its
**active** cortical tension floors at `γ_active ≈ 0.0001 mN/m` — ~6–10× below the
overlay band, robustly, across 8 seeds, constrained and unconstrained, turnover
ON and OFF, and ×40-converged (`H7_GATEB_FINDINGS_2026-06-07.md`). The root cause
is **not** a missing mechanism: it is a **timescale gap**. The tension-generating
and -remodelling processes are *seconds*-scale (myosin contraction builds over
~seconds; actin turnover τ½ = 10 s, Chugh 2017), while the fine-grained
constrained-MD step is ~µs, so a feasible run (10³–10⁴ steps) reaches only ~ms.
The cortex never reaches its active-tension steady state. This is exactly the
regime the active-gel / active-shell continuum (Salbreux–Prost–Joanny 2009;
Jülicher–Grill–Salbreux 2018; Borja da Rocha 2022) works in *directly* — it does
not integrate µs MD; it writes the *coarse-grained* constitutive law of the
actomyosin gel and relaxes it over seconds. The seam keeps the fine-grained cell
as the **mechanism source** (it measures the *instantaneous* local active stress
honestly, at F/F_stall ≤ 1, Hill-valid) and hands that coarse stress to a
seconds-scale active-gel continuum that performs the **slow relaxation** the MD
cannot reach. This is the *inverse* of the planned Layer-2 cortical→tissue seam
(there the fine-grained cell exports γ *up* to a tissue CBM; here it exports a
*stress* up to a continuum that returns the relaxed γ).

---

## 1. The two scales + what crosses the seam

### 1.1 The two scales

| | **Fine-grained (FG)** | **Active-gel continuum (AG)** |
|---|---|---|
| object | HOOMD particle/bond cortex (`build_cortex_full_simulation`) | active-gel constitutive PDE/ODE for the cortical shell |
| state | bead positions, bonds, M-SHAKE λ, myosin head binding | stress field σ(x,t), tension γ(t), strain-rate field v |
| timestep | dt ≈ 1.5e-9 s (cortex CFL) | dt_AG ≈ 0.01–0.1 s (seconds-scale relaxation) |
| reach | ~ms (10³–10⁴ steps) | seconds → steady state (cheap ODE/PDE) |
| role | **mechanism**: measures the *instantaneous* local active stress | **slow relaxation**: integrates that stress to steady γ |
| Hill validity | read at F/F_stall ≤ 1 (the active stress is a snapshot, not accelerated) | n/a — no per-head dynamics; ζΔμ is an input |

The FG cell is the only place where the mechanism (myosin force dipoles, catch
bonds, M-SHAKE backbone stress, turnover) lives at full fidelity (CLAUDE.md
architectural principle). The AG continuum is a *thin, parameter-poor* slow-time
relaxer whose parameters are **all set from the FG run** (§2), not free-fit.

### 1.2 What the FG cell EXPORTS (FG → AG), with units

The export is a **coarse-grained active stress** for the cortical shell, built
from the *same* three channels `cortical_tension.py` already measures, plus the
myosin force-dipole density that sets the active-stress drive. Concretely, the FG
cell exports a small bundle of seam variables, per region r (see §3 on
one-region vs binned):

1. **Active stress drive `σ_a` [Pa]** — the instantaneous actomyosin active
   stress, i.e. the myosin **force-dipole (stresslet) density**. Two equivalent
   readouts, cross-checked:
   - *Direct dipole sum*: `σ_a = (1/(A·h)) · Σ_heads (f_head · ℓ_dipole)` where
     `f_head` is the per-head contractile force [N] (= `k_head·(L−r₀)` of the
     head-actin attach bond, the *soft-channel* tension `cortical_tension._gamma_soft`
     already computes per bond), `ℓ_dipole` is the head-to-head minifilament span
     [m], `A` is the region cortical area [m²] and `h` the cortex thickness [m]
     (≈ 200 nm, KU-3.17). Units: `N·m / m³ = N/m² = Pa`. ✓
   - *Tension-equivalent*: `σ_a ≈ γ_soft / h` — the soft-MOP tension γ_soft [N/m]
     divided by the shell thickness h [m] gives the volumetric stress [Pa]. This
     reuses `measure_cortical_tension(...)['gamma_soft']` verbatim. ✓
2. **Backbone (rigid) stress `σ_r` [Pa]** = `γ_rigid / h` — the M-SHAKE Lagrange
   bond stress (`cortical_tension._gamma_rigid`, `γ_rigid·r₀/Δt`) divided by h.
   This is the elastic/structural load path; it sets the *prestress* the AG
   relaxes from (prestress matters for the viscoelastic response — Mokbel/
   Fischer-Friedrich; see §2).
3. **Active-stress coefficient inputs** (so ζΔμ is *derived*, not fit, §2.2):
   the engaged-head count `N_eng`, per-head stall force `F_stall_per_head` [N],
   minifilament dipole length ℓ_dipole [m], region area A [m²], thickness h [m].
4. **Turnover relaxation time `τ` [s]** — from the actin turnover τ½ (Chugh 2017,
   τ½ = 10 s ⇒ τ = τ½/ln2 ≈ 14.4 s) and/or the myosin head off-rate; this sets
   the Maxwell relaxation of the gel (§2.3). Read from `p_turnover` if present,
   else the literature anchor.
5. **Effective viscosity `η` [Pa·s]** and **shear modulus `G` [Pa]** — η from
   the cortex/cytoplasm (the FG run already runs at η_eff = 65.9 Pa·s MCF7
   cytoplasm, `cytoplasm.py`; the *cortical* η is η = G·τ, §2.4); G from the FG
   cortex's measured elastic response (small-strain stress/strain, or the
   literature cortical G).

All five are scalars (one-region first milestone) or per-bin fields (later). The
**load-bearing variable that crosses the seam is `σ_a` [Pa]** (the myosin
force-dipole / stresslet density); σ_r, τ, η, G are the constitutive parameters
that the *same* FG run supplies so the AG model is closed without free fitting.

### 1.3 What the AG continuum RETURNS (AG → FG / to the gate), with units

- **Primary return — the relaxed steady-state cortical tension `γ_ss` [N/m]**
  (report in mN/m for the band overlay). This is the quantity the GATE-B band
  `[0.35, 0.65] mN/m` is compared to. `γ_ss = σ_zz_relaxed · h`, the relaxed
  in-plane active stress integrated through the shell thickness.
- **Secondary (two-way only, §3.2) — a strain-rate field `v` [s⁻¹]** (or a
  velocity field) that the continuum predicts the cortex *would* flow at as it
  relaxes; fed back to bias the FG dynamics (e.g. an imposed slow compaction
  drift). The recommended minimal first version does **not** use this.

---

## 2. The active-gel continuum model

The minimal active-gel constitutive law is the **Maxwell (viscoelastic) active
fluid** of Jülicher–Grill–Salbreux 2018 (Reports on Progress in Physics) /
Joanny–Prost 2009 (HFSP J), specialised to a thin cortical shell à la
Salbreux–Prost–Joanny 2009 (PRL) and Borja da Rocha 2022 (J Mech Phys Solids).
It is deliberately the *minimal* closure — one relaxation time, one active term —
so every parameter is set by the FG run.

### 2.1 Constitutive law (the equation the seam integrates)

Deviatoric stress σ of an active Maxwell gel (Jülicher 2018, Eq. for the
Maxwell deviatoric stress with an active term):

```
τ · (Dσ/Dt) + σ = 2 η ε̇ + ζΔμ                      (1)
```

with, for the cortical shell, the in-plane (tension-relevant) projection. Here:

- `σ` = deviatoric stress tensor [Pa] (the in-plane component sets γ);
- `τ` = viscoelastic relaxation time [s] (set from turnover τ½, §2.3);
- `η` = cortical shear viscosity [Pa·s] (§2.4);
- `ε̇` = strain-rate tensor [s⁻¹] (zero in the 0-D relaxation-to-rest first
  milestone; the cortex relaxes at fixed area, ε̇ = 0);
- `ζΔμ` = **active stress** [Pa] — `ζ` the active-stress coefficient, `Δμ` the
  ATP-hydrolysis chemical-potential difference; this term is the myosin
  force-dipole density and is **the FG export σ_a** (§2.2). Sign: contractile ⇒
  ζΔμ < 0 in the SPJ sign convention (here we carry magnitude and the contractile
  sign explicitly so γ_ss > 0).

The relation between the moduli (Jülicher 2018; Maxwell):

```
η = G · τ        ⇔        τ = η / G                  (2)
```

where `G` is the elastic shear modulus [Pa]. Only two of {η, G, τ} are
independent; the seam fixes τ (turnover) and G (FG elastic response), and derives
η = G·τ — so η is **not** a free knob.

### 2.2 ζΔμ is SET from the FG myosin force-dipole density (not free-fit)

The active stress entering Eq. (1) is *literally* the FG export σ_a:

```
ζΔμ ≡ σ_a = (N_eng · ⟨f_head⟩ · ℓ_dipole) / (A · h)          (3)
```

- `N_eng` = engaged myosin heads in the region (FG: count of `myosin_action.
  _head_bound_to_actin ≥ 0`);
- `⟨f_head⟩` = mean per-head contractile force [N] (FG: head-actin attach bond
  tension `k_head·(L−r₀)`, the soft-channel per-bond tension; bounded by
  `F_stall_per_head` so the export is Hill-valid by construction, F ≤ F_stall);
- `ℓ_dipole` = minifilament head-to-head span [m] (FG layout geometry);
- `A·h` = region cortical volume [m³].

Equivalently and as a cross-check, `ζΔμ = γ_soft / h` using the already-validated
`gamma_soft` channel. **No tuning constant is introduced** — this is the
Magic-Number Block: ζΔμ is derived from measured FG quantities, grid-invariant
(it is intensive, per unit volume), and is NOT chosen to make the band pass.

### 2.3 τ is SET from turnover / off-rate (not free-fit)

```
τ = τ½ / ln 2                                                 (4)
```

with `τ½ = 10 s` (actin turnover, Chugh 2017 — already the project anchor) or the
myosin head off-rate `1/k_off` when that is the faster relaxation channel; the
seam takes `τ = min(τ_turnover, τ_myosin)` (the fastest relaxation dominates the
Maxwell time). Read from `p_turnover` if the FG run has turnover ON, else the
literature anchor with provenance recorded.

### 2.4 η, G

- `G` (cortical shear modulus) — from the FG cortex's *measured* small-strain
  elastic response (apply a small area strain, read Δσ/Δε from the soft+rigid
  channels), OR the literature cortical modulus anchor (G ~ 10²–10³ Pa, cortex)
  with provenance. Surfaced to PI if the FG-measured value is unavailable.
- `η = G·τ` from Eq. (2) — derived, not independent.

The cytoplasm viscosity η_eff = 65.9 Pa·s (MCF7, Hu 2024, `cytoplasm.py`) is the
*bulk drag* the FG run already uses; it is **not** the cortical shear viscosity η
of Eq. (1) — the seam must not conflate them (the cortical η is the actomyosin
network's own viscosity, η = G·τ). This distinction is a Sanity-Gate item (§5).

### 2.5 The 0-D steady-state solution (first milestone closed form)

At fixed area (ε̇ = 0) and steady state (Dσ/Dt = 0), Eq. (1) collapses to

```
σ_ss = ζΔμ = σ_a        ⇒        γ_ss = σ_a · h = γ_soft      (5)
```

— which would just return γ_soft and gain nothing. **The seam's value is the
transient, prestressed relaxation, not the trivial fixed point.** The minimal
*non-trivial* 0-D model the first milestone integrates is the active Maxwell
element relaxing from the FG prestress σ_r toward the active drive over τ, at
fixed area, capturing the *buildup the MD cannot reach*:

```
τ · dσ/dt + σ = ζΔμ ,     σ(0) = σ_r (FG backbone prestress) (6)
⇒  σ(t) = ζΔμ + (σ_r − ζΔμ) · e^{−t/τ}
⇒  γ(t) = σ(t) · h ,   γ_ss = ζΔμ · h     (t ≫ τ)
```

Critically, the *informative* first-milestone observable is **γ at the
seconds-scale steady state with the full active drive resolved**, where the FG
export σ_a is taken at the **Hill-valid quasi-steady local state** (heads engaged,
F/F_stall ≤ 1) rather than the MD-truncated ~ms value. The non-triviality comes
from (a) the FG σ_a being read at a *converged local* engaged-head population
(not the transient first-ms one), and (b) the prestress σ_r entering the
transient. If 0-D returns γ_ss = γ_soft and γ_soft is still ≈ 0, that is itself
the decisive result: it isolates whether the floor is a *timescale* problem (the
AG buildup lifts it) or a *generation* problem (σ_a itself is ≈ 0 even at the
converged engaged state) — see §5 honesty. The richer 1-D shell version
(SPJ/Borja da Rocha, advective contractile instability, flow-built tension) is
the *second* milestone and is where genuine flow-amplified tension can emerge.

---

## 3. The coupling protocol

### 3.1 Recommended minimal first version — ONE-WAY (FG → AG → γ)

```
   FG cell (HOOMD)                         AG continuum
   ───────────────                         ────────────
   build_cortex_full_simulation
        │  run to Hill-valid quasi-steady local state (~ms, engaged heads)
        ▼
   measure_cortical_tension(...)  ──►  extract seam bundle (§1.2):
        gamma_soft, gamma_rigid             σ_a = γ_soft/h   (= ζΔμ)
        myosin_action (N_eng, f_head)       σ_r = γ_rigid/h  (prestress)
        p_turnover / Chugh anchor           τ   = τ½/ln2
        FG elastic response / lit           G, η=G·τ
        cortex h, region A                  geometry
                                              │
                                              ▼
                                   integrate Eq. (6) over t ≫ τ
                                              │
                                              ▼
                                   γ_ss [N/m]  ──►  compare to band
```

One-way is the **right minimal first version**: it is the smallest thing that
tests the hypothesis ("the seconds-scale relaxation of the FG-measured active
stress reaches the band"), it cannot create feedback artifacts, and it keeps the
FG run untouched (the seam is a pure post-processor of a built cell, like the
GATE-B estimator). Data flow is a single hand-off bundle; no FG re-run inside the
AG loop.

### 3.2 Two-way (deferred) — continuum strain-rate → FG bias

The full multiscale loop feeds the AG-predicted strain-rate field `v` back into
the FG dynamics (e.g. impose a slow area-compaction drift, or bias myosin binding
toward the contracting region), re-measure σ_a, re-relax, iterate to a fixed
point. This couples flow-built tension (SPJ contractile instability) back to the
mechanism. It is **deferred**: it needs an FG re-run per AG iteration (expensive),
a careful timescale-separation argument (the FG must equilibrate faster than the
AG step), and a feedback-stability gate. Do NOT build it first.

**Recommendation: ship one-way 0-D first (milestone M1), then one-way 1-D shell
(M2), and only consider two-way (M3) if M2 shows flow-amplification is the
missing physics.**

---

## 4. Implementation plan (this codebase)

### 4.1 Where it slots in

New module **`ffn_sim/cortex/active_gel_seam.py`** — a *pure post-processor* of a
built cell, exactly like `cortex/cortical_tension.py` (one file = one concept).
It does **not** modify the FG runtime, the integrator, or any frozen core. It
reads from:

- `cortex/cortical_tension.py::measure_cortical_tension(...)` → `gamma_soft`,
  `gamma_rigid` (the σ_a and σ_r drives — reused verbatim, no duplication);
- the cell's `extras["handles"]`:
  - `handles["myosin_action"]` → engaged-head count, per-head force, dipole
    geometry (`_head_bound_to_actin`, `_head_bound_filament`, F_stall_per_head);
  - `handles["baoab_action"]` (+ `dt_used`) → the rigid Lagrange prestress (via
    `cortical_tension`'s rigid channel — passed through, not re-derived);
  - `handles["turnover_action"]` / `p_turnover` → τ½ for τ;
- `cell.p_cortex` → `R_cell`, `cortex_thickness` (h), region area A = 4πR²;
- `cell.p_myosin` → `F_stall_per_head`, `n_motors_per_cell`, dipole length;
- `cell.p_cytoplasm` → η_eff (bulk; for the documented distinction from cortical η).

### 4.2 Concrete functions

```python
# ffn_sim/cortex/active_gel_seam.py  (proposed signatures)

@dataclass(slots=True)
class SeamExport:
    """Coarse-grained FG → AG bundle (all SI). One region (M1)."""
    sigma_a: float        # active stress drive ζΔμ [Pa]  (= γ_soft / h)
    sigma_r: float        # backbone/rigid prestress  [Pa]  (= γ_rigid / h)
    tau: float            # Maxwell relaxation time   [s]   (= τ½ / ln2)
    G: float              # cortical shear modulus    [Pa]
    eta: float            # cortical shear viscosity  [Pa·s] (= G·τ, derived)
    h: float              # cortex thickness          [m]
    A: float              # region cortical area       [m²]
    n_engaged: int        # engaged myosin heads (provenance)
    provenance: dict      # which anchor each param came from (turnover, lit, FG)

@dataclass(slots=True)
class ActiveGelParams:
    """Minimal active Maxwell gel constitutive parameters (from SeamExport)."""
    sigma_a: float; sigma_r: float; tau: float; eta: float; G: float; h: float

def export_seam_state(
    cell, *, gamma: dict, region: str = "whole_cortex",
    tau_half_s: float | None = None, G_Pa: float | None = None,
) -> SeamExport:
    """Build the FG → AG seam bundle from a built cell + a measured `gamma`
    dict (the output of measure_cortical_tension). σ_a, σ_r from γ_soft/h,
    γ_rigid/h; τ from turnover τ½; η = G·τ. No tuning constant introduced."""

def relax_active_gel_0d(p: ActiveGelParams, *, t_max_s: float,
                        dt_s: float) -> dict:
    """Integrate the 0-D active Maxwell element Eq. (6) from prestress σ_r
    toward the active drive σ_a over τ, at fixed area (ε̇ = 0). Returns the
    σ(t) trace and the steady γ_ss = σ_ss · h [N/m]."""

def predicted_gamma(cell, *, gamma: dict, **kw) -> dict:
    """One-call: export_seam_state → relax_active_gel_0d → {gamma_ss [N/m],
    band overlay, channel provenance}. The seam's GATE-B-comparable output."""

# M2 (second milestone, same file or a sibling):
def relax_active_gel_1d_shell(...):  # SPJ/Borja-da-Rocha thin-shell, flow-built
```

A **coarse-graining step** lives inside `export_seam_state`: for the M1 single
region it is the whole-cortex average (A = 4πR², all heads); for M2 it bins the
FG stresslets into a small set of latitudinal shell patches (Fibonacci or
polar bins) so σ_a(x), σ_r(x) become fields the 1-D shell relaxes — the binning
reuses the same per-bond positions `cortical_tension` already reads.

### 4.3 First milestone (M1) — 0-D single-region relaxation → predicted γ

**Goal**: a 0-D active-gel relaxation, *driven by the FG-measured myosin
force-dipole density*, predicting a steady γ_ss, compared to the band — the
smallest end-to-end seam.

Steps:
1. Build a full physiological MCF7 cell (the GATE-B build) and run it to a
   Hill-valid quasi-steady *local* state (engaged heads, F/F_stall ≤ 1).
2. `gamma = measure_cortical_tension(...)` (reuse the GATE-B call).
3. `export = export_seam_state(cell, gamma=gamma, tau_half_s=10.0, ...)` →
   σ_a = γ_soft/h, σ_r = γ_rigid/h, τ = 10/ln2 s, η = G·τ.
4. `out = relax_active_gel_0d(ActiveGelParams(**export), t_max=5·τ, dt=0.01)`.
5. Report γ_ss [mN/m] vs the band; record provenance of every parameter.

**Driver**: `ffn_sim/scripts/h7_active_gel_seam.py` (auto-viz at run end per the
production-driver rule → `outputs/h7/figs/h7_active_gel_relaxation.png`:
σ(t)/γ(t) trace from σ_r toward γ_ss, with the band overlaid — visualization-
integrity: no axis truncation, band overlaid, SI units annotated).

**Tests**: `ffn_sim/tests/test_active_gel_seam.py` — Sanity Gate §5 below
(dimensional, t=0 = σ_r, t→∞ = σ_a·h, τ-monotone, ε̇=0 fixed point, η=G·τ).

### 4.4 What it explicitly does NOT touch

The frozen integrator/, the FG builder, the GATE-B estimator, the brief specs.
The seam is additive and read-only over the FG cell, exactly like the γ estimator.

---

## 5. Honesty / validity

### 5.1 How this AVOIDS the accelerated-dynamics Hill-invalidity trap

The discarded alternative (PI path option (a)) is to **raise the myosin v0 /
turnover rates** so the FG MD reaches steady state in feasible steps. That is
invalid because Hill force–velocity is only meaningful at F/F_stall ≤ 1; cranking
v0 drives heads super-stall and fabricates an in-band number that is **not** a
real tension (the KU-3.5 v0_accel memory flags exactly this). The seam **never
accelerates the per-head dynamics**: the FG cell runs at its *physiological*
rates and exports the active stress σ_a measured at a **Hill-valid** local state
(F ≤ F_stall by construction — σ_a is built from head forces *bounded by
F_stall*, §2.2). The seconds-scale integration happens in the *continuum*, where
there are no per-head dynamics to invalidate — ζΔμ is a steady input, not an
accelerated motor. The slow time is carried by the *constitutive relaxation time
τ* (a real, literature-anchored turnover time), not by a sped-up motor. So the
seam buys the seconds-scale steady state **without** the Hill-invalidity that
accelerated dynamics incurs.

### 5.2 What it CAN claim

- A *physically-grounded* prediction of the seconds-scale steady cortical tension
  γ_ss from the FG-measured instantaneous active stress, with **every continuum
  parameter set from the FG run or a cited literature anchor** (no free fit, no
  band-chasing) — i.e. it tests whether the floor is a *timescale* artifact.
- A clean separation of the two competing explanations: if γ_ss reaches the band,
  the floor was a timescale gap (the mechanism is fine; the MD just couldn't
  relax it). If γ_ss is *still* below band even after the seconds-scale
  relaxation, the floor is a **generation** limit (σ_a itself is too small at the
  physiological engaged-head population) — a sharper, more honest result than the
  current "active ≈ 0 at ms".
- A standard, defensible method (active-gel theory is the canonical cortex
  continuum) — it directly answers the novelty-analysis threat ("emergent γ is a
  promise, not a result") by making the seconds-scale claim explicit and oracle-
  checked.

### 5.3 What it CANNOT claim

- It is **not** a first-principles fine-grained derivation of γ at seconds scale
  — it is a *bridge*: the FG cell provides the active stress, the continuum
  provides the relaxation. The seconds-scale tension is a *continuum* result fed
  by a fine-grained drive, not a single fully-fine-grained number. State this
  plainly in any deliverable.
- The 0-D M1 cannot produce *flow-amplified* tension (the SPJ contractile
  instability that builds extra tension via cortical flow) — that needs the 1-D
  shell (M2). If the band needs flow amplification, M1 will under-predict and
  that is the honest M1 read.
- It cannot rescue a *wrong baseline*: the FG export must come from a cell at its
  physiological operating point (turgor ON, η_eff = 65.9 Pa·s, FA at physiological
  force — the physiological-baseline HARD rule). Exporting σ_a from an unphysical
  null baseline poisons γ_ss.
- γ_passive (turgor) stays a SEPARATE channel (GATE-B B3 rule) — the seam relaxes
  the **active/structural** stress only; it must NOT fold turgor into γ_ss.

### 5.4 Acceptance oracle

The seam's correctness oracle is **NOT** the band (band-chasing is forbidden).
The oracle is the **analytic active-Maxwell solution** Eq. (6) and the
active-gel limits:

| # | Oracle check | Pass criterion |
|---|---|---|
| O1 | Dimensional | σ_a, σ_r in Pa; τ in s; η in Pa·s; γ_ss = σ·h in N/m. ✓ each. |
| O2 | t = 0 | σ(0) = σ_r (the FG prestress) exactly. |
| O3 | t ≫ τ | σ(t) → ζΔμ = σ_a; γ_ss = σ_a·h = γ_soft (the §2.5 fixed point). |
| O4 | Maxwell relaxation | σ(t) = σ_a + (σ_r−σ_a)e^{−t/τ}; matches the closed form to machine ε. |
| O5 | η = G·τ | derived η equals G·τ; τ derived from τ½/ln2; no independent η knob. |
| O6 | ε̇ = 0 fixed point | at fixed area the steady stress is exactly ζΔμ (active fluid limit). |
| O7 | Hill validity | every per-head force entering σ_a satisfies F ≤ F_stall_per_head. |
| O8 | Provenance | every continuum parameter carries a source tag (FG-measured / lit anchor); none is a tuned magic number (Magic-Number Block). |
| O9 | Turgor separation | γ_ss excludes γ_passive (B3). |

The *scientific* comparison to the `[0.35,0.65] mN/m` band is an **overlay**, not
a gate the seam is tuned to (it is itself a non-MCF7 rounded-cell proxy under PI
review). The seam PASSES on O1–O9; the band tells us *which explanation*
(timescale vs generation) the data supports.

---

## 6. References (verified 2026-06-07)

All four anchor papers verified via web search / publisher this session:

- **Salbreux, Prost & Joanny 2009** — "Hydrodynamics of Cellular Cortical Flows
  and the Formation of Contractile Rings", *Phys. Rev. Lett.* **103**, 058102
  (2009). Thin-shell active-gel cortex; flow couples to filament orientation,
  driven by myosin density. (verified — APS DOI 10.1103/PhysRevLett.103.058102)
- **Joanny & Prost 2009** — "Active gels as a description of the actin–myosin
  cytoskeleton", *HFSP Journal* **3**(2), 94–104 (2009). Hydrodynamic theory of
  active polar gels; active stress ζΔμ; viscoelastic (Maxwell) response.
  (verified — PMID 19794818; PMC2707794; DOI 10.2976/1.3054712)
- **Jülicher, Grill & Salbreux 2018** — "Hydrodynamic theory of active matter",
  *Reports on Progress in Physics* **81**(7), 076601 (2018). Irreversible-
  thermodynamics derivation of active-gel constitutive equations; the active
  Maxwell form τσ̇ + σ = 2ηε̇ + ζΔμ and η = Gτ used here. (verified — PMID
  29542442; DOI 10.1088/1361-6633/aab6bb; full text on pks.mpg.de, fetched.)
- **Borja da Rocha, Bleyer & Turlier 2022** — "A viscous active shell theory of
  the cell cortex", *J. Mech. Phys. Solids* **164**, 104876 (2022). Reduces 3-D
  incompressible viscous active gel to a Koiter-like thin active shell (stretch +
  bend), the M2 1-D shell target. (verified — arXiv:2110.12089; HAL hal-03403141;
  DOI 10.1016/j.jmps.2022.104876)

Supporting (already project anchors, not re-verified here): Chugh et al. 2017
(actin turnover τ½ ≈ 10 s, cortical tension); Salbreux, Charras & Paluch 2012
(cortical-tension band, `cortical_tension.py`); Hu et al. 2024 (MCF7 cytoplasm
η_eff = 65.9 Pa·s, `cytoplasm.py`).
