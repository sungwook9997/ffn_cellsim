# ff — Filament-FEM engine (build spec)

*Stage-6 of the 2026-06-29 two-layer restructure. The going-forward fine-grained engine:
a discrete cross-linked filament network solved as a **mechanical (FEM-like) system** —
energy assembly + implicit nonlinear solve + autodiff — instead of explicit thermostatted MD.*

Physics reference: **Cytosim** (Nédélec & Foethke 2007, *New J. Phys.* **9**:427) — exactly as
the `dcm/` layer references SimuCell3D. Cytosim itself (runnable external C++) is the independent
**parity oracle**, so the active layer is never validated against only itself.

---

## 1. Why MD-free is faithful here (the design rationale, PI-discussed 2026-06-29)

A crosslinked filament network is structurally a **beam/truss finite-element network** (crosslinks =
nodes, filament segments = elastic elements, bending = rotational stiffness). The question is whether
to integrate it as Langevin MD (current `archive/hoomd_legacy`) or solve it as mechanical equilibrium.
Decision rule — **does kT enter the observable through the mean, or through the variance / rare tail?**

| regime | mechanism | treatment | cost |
|---|---|---|---|
| **A — integrable-out** | entropic (WLC) elasticity | closed-form constitutive force (Marko–Siggia); the thermal effect becomes a deterministic term | **free** |
| **B — kinetic layer** | catch/slip bonds, myosin stepping (force-dependent, ATP-driven) | event-driven KMC on the mechanical solve (not per-step thermostat) | cheap; carries its own kinetic timescale |
| **C — must sample** | *pure-thermal* rare-event configurational remodeling | cannot be integrated out; keep stochastic | the only thing that needs sampling |

Key point for an **active** cell: most remodeling is ATP-rate-driven (regime B, clean kinetic layer),
not pure-thermal exploration (regime C). So the MD-free mechanical solve is *more* applicable here than
for a passive thermal gel. Quenched/structural disorder (initial network architecture, bound states)
moves from the time-axis to a **parallel realization ensemble** (Monte-Carlo at init) — exact for slow
disorder, and embarrassingly parallel. See the γ-floor protocol (§4).

> ⚠️ NOT lumping: regimes A and B keep every mechanistic force law (WLC, Bell-Evans, Hill) — only the
> *integration method* changes. This satisfies the project's mechanistic-over-lumped hard rule.

---

## 2. Cytosim physics — GROUNDED (Nédélec & Foethke 2007; PDF in `references/`)

Verbatim from the paper with page refs (implemented values are inputs, never invented).

**Fiber + bending (p9, fig.6) — implemented + validated in [`forces_warp.py`](forces_warp.py) ✓.**
A fiber = `p+1` equidistant model-points, segment length `seg = L/p`. For each interior point the
consecutive triplet {m_{i-1}, m_i, m_{i+1}} carries the force triplet {−F, +2F, −F} with

    F = α (m_{i-1} − 2 m_i + m_{i+1}),   α = κ (p/L)³ = κ / seg³            [NF2007 p9]

κ = bending modulus = `k_B·T·L_p = E·I` (standard identity the paper assumes; it supplies
κ = 20 pN·µm² for a microtubule). Discrete energy `E = (α/2) Σ |m_{i-1}−2m_i+m_{i+1}|²`; assembled
stiffness = symmetric banded `α·E` (interior 4th-diff stencil `[−1,4,−6,4,−1]`, reduced ends,
row/col sums 0 → internal torque only). Validation anchor: buckling threshold = Euler's `π²κ/L²`.
Tests (`tests/ff/`): straight→0, restoring, **force = −∇E (FD)**, explicit relaxation straightens.

**No axial spring (p3, §5.3).** Non-extensibility is a hard length constraint
`C_k = (m_{k+1}−m_k)² − seg² = 0` + a "reshape" projection — NOT a penalty spring (the paper
explicitly rejects the spring). → the constraint projector `P = I − Jᵀ(JJᵀ)⁻¹J` is Stage-6b's
remaining piece; bending-only relaxation currently holds segment length only approximately.

**Integration — Eq (2), p8 (the key contribution).** Overdamped `dx = μF dt + dB` (Eq 1), linearise
`F = A x + G`, solve semi-implicitly:

    [I − τ Pₜ μ Aₜ] (x_{t+τ} − xₜ) = Pₜ [ τ μ (Aₜ xₜ + Gₜ) + δBₜ ]          …Eq (2)

A = elastic stiffness (bending + attractive links), implicit; G + Brownian explicit; P = constraint
projector. **Unconditionally stable** (p17): A negative-semidefinite (compression→constraint) ⇒
eigenvalues of (I − τμPA) > 1 ∀τ ⇒ dt bounded by *accuracy* (O(τ²)) not stability — a 10⁴× speed-up.
Non-symmetric after PμA ⇒ paper uses **BiCGStab** (tol = 0.1·min Brownian). Our
`dcm.dcm_warp_implicit.implicit_overdamped_step` is the matrix-free `(γ/dt·I + K)Δx = F` CG form —
**one implicit step relaxes a bent fiber to ~1.4 % bending energy ✓** (the large-step payoff).
  > ✅ **Scale finding RESOLVED (Stage 6c-a):** that solver's absolute thresholds (`newton_tol=1e-10`,
  > `cg_tol=1e-8`) are DCM-scale-calibrated; FF single-filament forces are ~1e-12 N (SI) → instant
  > false convergence (Δx=0). Fix = **nondimensionalize to pN·µm·s** ([`units.py`](units.py), forces
  > O(0.1–100)). Validated: at FF scale the *unchanged* solver relaxes a bent actin fiber with its
  > DEFAULT tolerances, while the identical SI-scale fiber false-converges (`tests/ff/test_units.py`).

**Mobility (p10, §5.2).** Isotropic scalar per point (deliberately *not* anisotropic):
`μ = log(L_h/δ) / (3π η L)`, per-point `μ_p = (p+1)μ` (L_h = min(L, hydro cutoff), δ = diameter).
Brownian: `δB = β θ`, `θ∼N(0,1)`, `β = √(2Dτ)`, `D = μ_p k_B T` (p8, p11).

**Hand binding (p20–21) — Stage 6c.** Attach prob/step = `τ k_on` to the closest site within capture
radius ε; active step `δa = τ v_max (1 − f/f_stall)`; Bell off-rate `p_off = p₀ exp(|f|/f₀)`. Kinesin
(p19): v_max=0.4 µm/s, f_stall=5 pN, k_on=10 s⁻¹, p₀=0.5 s⁻¹, f₀=2.5 pN, k=200 pN/µm, ε=10 nm.

---

## 3. Reuse map — ~80% already built in `dcm/` (verified 2026-06-29)

| FF piece | reuse | adaptation |
|---|---|---|
| element force assembly | [`dcm/network_warp.py`](network_warp.py) — harmonic bond + harmonic angle + LJ, one-thread-per-element, atomic accumulate | swap harmonic-bond → **WLC (Marko–Siggia)** constitutive; harmonic-angle → **Cytosim bending** |
| **implicit "FEM solver"** | [`dcm/dcm_warp_implicit.py`](dcm_warp_implicit.py) — `implicit_overdamped_step(pos, force_fn, γ, dt, …)`, matrix-free CG + Newton, **force-agnostic**, precond opt-in | quasi-static γ: `dt→∞` ⇒ Newton on `F(x)=0`; pass the FF force as `force_fn` |
| autodiff (inverse) | [`dcm/dcm_warp_diff.py`](dcm_warp_diff.py) — `wp.Tape` over the whole loop, FD-validated 4e-12 | back-prop γ-loss → fiber κ, prestress |
| geometry | `dcm/geometry.py`, `common/surface_manifold.py` | cortex shell seeding |
| kinetic layer template | `archive/hoomd_legacy` cadherin / junction-switch hosts (force-dependent bond update) | Hand KMC on the mechanical solve |

New code needed: WLC + Cytosim-bending kernels, the Hand KMC layer, the quenched-ensemble harness,
the fiber-network data model ([`fiber_network.py`](fiber_network.py), scaffolded).

---

## 4. First milestone — γ-floor prototype (the decisive experiment)

γ (cortical tension) is **quasi-static** → the load barely changes during measurement → even the
load-coupled myosin has a well-defined steady distribution that can be sampled once (quenched). Protocol:

1. MC-sample N realizations of the cortex fiber network + myosin bound-state at init (regime A/B disorder).
2. Each realization: WLC-beam network + myosin active prestress → **implicit mechanical equilibrium**
   (`implicit_overdamped_step`, `dt→∞`), no explicit BAOAB.
3. Run the N realizations **in parallel** → distribution of γ.
4. **Decisive validation:** compare MD-free γ vs the BAOAB-MD γ from the archived
   `cortical_tension.py` (the γ-floor finding). "Does an MD-free mechanical solve reproduce the same
   γ?" → confirms FF is MD-equivalent for γ, *or* localizes exactly where in-time load-coupling is
   load-bearing.

Ties to the open **γ anchoring** question (Moazzeni MCF7 1e-2 N/m vs SimuCell3D 1e-3 vs emergent;
no clean adherent-MCF7 datum) — γ should be **swept as a controlled variable**, faceting/tension
reported as a function of it (NOT tuned to a target). See `docs/v2_audit/PARAM_AUDIT_SIMUCELL3D_2026-06-25.md`.

---

## 5. Build stages

- **6a** ✅ this spec + fiber-network data model scaffold.
- **6b** ✅ Cytosim-bending Warp kernel (grounded) + wired to `implicit_overdamped_step`; relaxation
  validated (`forces_warp.py`, `tests/ff/test_cytosim_bending.py`).
- **6b-iv** ✅ inextensibility = constraint projector `P = I − Jᵀ(JJᵀ)⁻¹J` + exact reshape (NF2007
  §5.3, NOT a spring) + per-segment axial tension from the multipliers (`constraints.py`).
- **6c-a** ✅ nondimensionalize to **pN·µm·s** (`units.py`) — resolves the implicit-solver scale
  finding (default tolerances now work at FF scale).
- **6c-b** ✅ cortex fiber-network on the sphere, lit-anchored (`cortex_assembly.py`, `viz_cortex.py`).
- **6c-c** ✅ Hand KMC kinetic layer (NF2007 §10.1: attach / δa=τv(1−f/f_stall) / Bell + Pereverzev
  off-rate) with kinesin/NMIIA/α-actinin/filamin presets (`hand_kmc.py`).
- **6d** ✅ γ-floor prototype: method-of-planes γ estimator (`gamma_estimator.py`, faithful HOOMD
  port) + crosslinked+myosin+turgor cortex equilibrium + quenched ensemble + prestress sweep
  (`gamma_floor.py`, `gamma_floor_sweep.py`). **Result:** the MD-free solve REPRODUCES the BAOAB-MD
  γ-floor (actomyosin ~370× under band at the lit NMIIA prestress; turgor at-band) → the floor is
  STRUCTURAL (force-magnitude/transmission), not a dynamic MD artifact. See
  `docs/v2_audit/FF_STAGE6D_GAMMA_FLOOR_2026-06-29.md`. **PI decision open** (missing motor-density
  datum — same as the SF/NMII line).
- **next** Cytosim runnable parity oracle (the independent external check); WLC (Marko–Siggia)
  constitutive option; implicit-accelerated equilibrium for the loaded (not just resting) shell.

Full `tests/ff` suite: **43/43**.

## 6. Open / PI-gated before deep build

- **Cytosim paper grounding:** Nédélec & Foethke 2007 is now in `references/`
  (`Nedelec_2007_New_J._Phys._9_427.pdf`, the published IOP version; arXiv 0903.5178 mirror also
  present) — ground the §2 formulae (discrete bending operator + implicit scheme + Hand model) from
  it before any kernel coding (hard rule). SE candidate written
  (`references/SE_REGISTRATION_CANDIDATES_2026-06-29.md`).
- **γ target band:** unresolved (Moazzeni vs SimuCell3D vs emergent) — sweep, don't tune.
- **Active-driver force anchor:** still REFUTED/HALTED to PI (the SF/NMII missing motor-density datum) —
  the spreading magnitude driver remains PI-gated, not auto-tunable.
