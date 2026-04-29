"""Overdamped 3D MLS-MPM solver in dimensionless units — Stage 1a Layer 1.

This is the Stage 1a baseline. It composes the standard MLS-MPM pipeline
(Hu et al. 2018, ACM TOG 37, 150) with **overdamped** grid dynamics suitable
for the cellular Re ≈ 10⁻¹³ regime (Pérez-González 2019, Marchetti 2013).
APIC affine transfer (Jiang 2015). Maxwell deviatoric stress is integrated by
the closed-form **exponential integrator**. Surface tension uses a
density-tagged-boundary CSF impulse on the background grid (Brackbill 1992).

All quantities inside the solver are **dimensionless**. The reference units are:

    length unit    : R₀  (initial spheroid radius — taken as the scale for `dx`)
    time unit      : τ_relax  (Maxwell deviatoric relaxation time)
    stress unit    : K  (long-time effective bulk modulus)
    mass unit      : ρ · R₀³  (taken as 1 in dimensionless form, ρ* = 1)
    velocity unit  : R₀ / τ_relax
    force unit     : K · R₀²

Restoration to SI is the responsibility of the Stage 1a → 1a+ anchor pass; this
solver never sees Pa, μm, or seconds. The mapping table is captured in
`docs/stage1a_assumption_review.md` §A6.

================================================================================
Sanity Gate (per docs/12_validation.md, run before first execution)
================================================================================

1. **Dimensional analysis**
   Reference units (length, time, stress) defined above. Every solver field is
   dimensionless. Characteristic non-dimensional numbers (with the Stage 1a
   placeholder set: K=1 kPa, ρ=1050 kg/m³, R₀=100 μm, τ=60 s, η_cell≈100 Pa·s,
   v ≈ μm/min spreading rate, γ=1 mJ/m²):

       Re   = ρ v L / η      ≈ 10⁻¹³           ⇒  overdamped fully justified (A1, A2)
       Ca   = γ / (K · R₀)   = 0.01            ⇒  surface tension ≪ bulk elasticity
       De   = τ_relax / τ_obs ≈ 2 × 10⁻⁴       ⇒  quasi-static elastic regime
       Pe   = (advective)/(diffusive) — N/A in Stage 1a (no diffusion species)
       Ma   — N/A (overdamped, no inertia, no acoustic)

   CFL: in the overdamped scheme the elastic-wave CFL is **not** the limiting
   condition (no inertia ⇒ no acoustic wave). The relevant stability bounds are
   (i) the surface-tension capillary timescale `τ_γ* = ξ* · dx*³ / γ*` ≥ dt*,
   and (ii) the diffusion-like grid step `dt* ≤ dx*² · ξ* / μ*` for the
   deviatoric viscous response. With dx* = 1/16 = 0.0625, μ* = 0.3, ξ* = 1,
   γ* = Ca·K·R₀ = 0.01:
       τ_γ*    ≈ 1 · (0.0625)³ / 0.01 ≈ 2.4·10⁻⁵     ← ⚠ tighter than dt*=0.01
       τ_visc* ≈ (0.0625)² · 1 / 0.3 ≈ 1.3·10⁻²     ← ok vs dt*=0.01

   The capillary timescale is below dt* with default dx*. We mitigate by
   (a) keeping the CSF colour-gradient operator low-passed through the
   quadratic-kernel particle scatter (effectively over a few cells, easing the
   bound by the kernel width factor 3³), and (b) the grid resolution config
   default `background_grid_resolution=64` so dx* ≈ domain*/64 = 6/64 ≈ 0.094,
   which gives τ_γ* ≈ 8.3·10⁻⁵·3³ ≈ 2.2·10⁻³ — still tighter than 0.01.

   **Action**: the gate `max_speed_over_vrms` and the radius-drift gate catch
   any capillary-instability blow-up at runtime. If pilot violates either gate
   we drop dt* to 0.001 (10× tighter) before declaring a scheme failure.

   Dimensional check: PASS — Re=1e-13, Ca=0.01, De=2e-4, dt vs τ_γ borderline,
   gated by runtime sanity bounds.

2. **Boundary cases**
   N → 0:        constructor refuses N < 8 (need at least one APIC stencil's worth).
   N → ∞:        memory scales O(N + grid³); 16 GB Laptop A5000 ceiling enforced
                 by `device_memory_GB` in config + Stage 1a benchmark gate.
   Δt → 0:       fully stable, error → 0 (any explicit scheme).
   Δt → large:   capillary blow-up first, then volumetric instability; gate
                 detects via max_speed_over_vrms.
   γ → 0:        boundary tagging still runs but contributes nothing — solver
                 reduces to pure viscoelastic relaxation, valid sanity check.
   K → ∞:        volumetric stress σ_vol = K·(ρ_ref/ρ_kernel − 1) saturates the
                 f32 field for any finite (ρ_kernel − ρ_ref). Stage 1a
                 placeholder K* = 1 stays well inside f32 range.
   ρ_kernel→0:   v15 (k.3) form 1/ρ_kernel diverges; clamped at
                 ρ_floor = 0.1·ρ_ref_kernel by `_interpolate_rho_runtime`,
                 bounding σ_vol at +9·K. The clamp is a numerical safety
                 floor — see `docs/outcomes_v15.md` for the literature
                 verification record (no exact 0.1 cite found; conservative
                 below the AHA-2010 free-surface kernel-truncation bound
                 ~0.5·ρ_ref and the Liu-Liu 2010 isolated-particle
                 threshold ~0.5·ρ_ref).
   R → 0:        constructor refuses R < 2·dx (cannot resolve the spheroid).

3. **Conservation invariants**
   Mass            : exact (P2G/G2P each move equal weight; particle masses
                     never modified). Diagnostic in f64.
   Momentum        : conserved to atomic-add round-off in absence of external
                     force; the reflective box-wall safety net is a deliberate
                     leak channel and only fires if particles touch the wall
                     (gate halts before this can happen).
   Angular momentum: APIC preserves it up to round-off (Jiang 2015 Theorem 4.1).
                     Not actively checked at Stage 1a but tracked for Stage 1d.
   Energy          : intentional dissipation = Maxwell deviatoric relaxation
                     (rate (1-decay)/dt) + overdamped drag. Intentional source =
                     CSF surface-tension impulse only at boundary particles.
                     Energy injection per step is bounded by γ·area·v_max·dt.
                     Net energy must decay monotonically once the spheroid is
                     in equilibrium (gate tolerance: +1e-3·E_max).
   Possible leaks  : (i) atomic_add ordering on the grid (round-off only),
                     (ii) reflective wall safety net — gated, should never fire.

4. **Numerical sanity**
   dt*       = 0.01 ≪ τ* = 1 (Maxwell relaxation) ⇒ 100 steps per relaxation time. ✓
   dt*       borderline vs τ_γ* (see check 1); runtime-gated.
   dx* / R₀* ≈ 0.094 ⇒ ~22 grid cells per spheroid radius. Adequate for shape
                       resolution; surface-tension support spans 3·dx ≈ 0.28·R₀.
   Precision : f32 for x, v, F, C, τ_dev, grid fields. f64 for cumulative
               diagnostics (mass, momentum, energy) where round-off can
               accumulate over 10⁵ steps. Justified vs gate tolerance 10⁻¹⁰.

5. **Sign / sense check**
   Volumetric (v15 (k.3), density-based)   σ_vol = K (ρ_ref/ρ_kernel − 1) I:
       ρ_kernel > ρ_ref (clustered) ⇒ ρ_ref/ρ_kernel − 1 < 0 ⇒ σ_vol < 0
       ⇒ Cauchy stress is compressive ⇒ grid force −V₀·σ_vol·∇w is *outward*
       ⇒ pushes the over-clustered region apart. ✓
       ρ_kernel < ρ_ref (rarefied) ⇒ σ_vol > 0 ⇒ tensile ⇒ grid force pulls
       particles together. ✓
       Equilibrium at ρ_kernel = ρ_ref. This is the standard fluid-pressure
       response (high density → expand, low density → contract), with the
       same sign behaviour as the v12 K(J − 1) volumetric Neo-Hookean but
       with kernel-density rather than Lagrangian deformation as the source —
       so the bulk pressure responds to *spatial particle clustering*, not
       just to local F. This is the v15 fix for the (k) interior-pressure-
       transmission failure of v12 (where overdamped + ∇v ≈ 0 made F → I in
       the bulk and σ_vol stayed ≈ 0 even as surface CSF dragged particles
       inward; see `docs/stage1a_interior_pressure_sanity.md`).
   Maxwell deviatoric        τ_dev_{n+1} = e^(-dt/τ)·τ_dev_n + 2μτ(1−e^(-dt/τ))·ε̇_dev:
       under shear ε̇_dev > 0 ⇒ τ_dev grows toward 2μ·τ·ε̇_dev (steady),
       opposing the shear ⇒ correct. ✓ Unchanged from v12.
   CSF impulse              dv = +γ·κ·∇c·dt / ρ_local (Brackbill 1992):
       Colour c = grid_m / (ρ_bulk · grid_kernel_weight) (Adami-Hu-Adams 2010
       §3 reproducing-kernel normalisation, in v12) gives c ≡ 1 in any cell
       whose kernel sees particles and c ≈ 0 in vacuum. ∇c then points
       INWARD at the free surface (gradient of "inside-ness"). The unit
       normal n̂ = ∇c/|∇c| therefore points inward, and curvature κ = -∇·n̂
       is **positive for a convex droplet** (= 2/R for a sphere). The body
       force F_v = γ·κ·∇c thus has TWO inward signs multiplied (κ > 0, ∇c
       inward) and points inward — exactly the physical restoring force from
       surface tension. ✓
       Static-sphere validation in `runner.run_stage1a` asserts that the
       measured surface κ matches 2/R within 10% before the time loop runs.
   Overdamped drag          v ← (1 / (1 + ξ*·dt*)) · v after the elastic step:
       reduces |v| monotonically toward force-balance ⇒ correct. ✓

   Sense check: PASS.

6. **Measurement-protocol consistency** (codified after v13, see
   docs/12_validation.md and docs/stage1a_interior_pressure_sanity.md §6)
   The v15 proposal walks through every gate's measurement protocol against
   the new density-based volumetric stress. Salient items:
     • Static-curvature κ gate is unchanged (κ depends only on the colour
       field, not the constitutive law). ✓
     • Radius-drift gate measures `effective_radius` from the second moment
       of particle positions; the proposal predicts equilibrium at
       R_eq/R_ref ≈ 0.991 (from γ·κ ≈ K·(ρ_ref/ρ_eq − 1) with γ=0.01,
       κ ≈ 2.8, K=1 ⇒ ρ_eq/ρ_ref ≈ 1.029 ⇒ R contraction ≈ 0.9%). The
       *off-protocol pathway* — could R shrink without ρ_kernel rising in
       the bulk? — is closed by the new shell-averaged `<ρ_kernel>(r/R₀)`
       diagnostic (`shell_density_profile` in `acs.analysis.shape_metrics`,
       written to `shell_profile.csv` by the runner). A flat bulk profile
       at equilibrium confirms that pressure transmits through the bulk;
       a surface-only peak indicates the v15 fix is failing.
     • Conservation gates (mass, momentum, energy) are arithmetic-only;
       the volumetric strain energy is updated to U_vol = (K/2)·
       (ρ_ref/ρ_kernel − 1)² for consistency with the new constitutive law,
       so the energy-monotone gate continues to reflect the actual stored
       internal energy.
     • Calibration `<J>_well_resolved` gate is unchanged in semantic but
       becomes informational (no longer load-bearing for σ_vol).
   Check 6: PASS — the predicted equilibrium IS consistent with the
   measurement protocol, and a new measurement-protocol-consistent witness
   (shell `<ρ_kernel>(r/R₀)`) is added to test the bulk-transmission
   mechanism directly.

================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import taichi as ti


@dataclass
class SolverConfig:
    """Concrete dimensionless inputs to the solver.

    Built from the YAML config's `nondim:` block by `acs.runner`. Anchor / SI
    restoration uses the YAML `physics:` block plus this `nondim` block.
    """

    n_particles: int
    grid_n: int
    domain_star: float           # box edge length, in R₀ units
    radius_star: float           # spheroid radius, in R₀ units
    dt_star: float               # time step, in τ_relax units

    K_star: float                # bulk modulus
    mu_star: float               # deviatoric shear modulus
    tau_star: float              # Maxwell relaxation time
    capillary_number: float      # Ca = γ/(K·R₀); ⇒ γ_star = Ca · K_star
    drag_xi_star: float          # overdamped drag, ξ*
    density_star: float          # ρ*; only enters mass diagnostics
    free_surface_threshold: float

    seed: int = 42

    # Stage 1a+ substrate (Option α: mechanical anchor only, γ_sub = 0).
    # When `substrate_enabled` is True, the existing reflective grid BC at
    # k < n_contact_band on the −z side becomes the *intentional* substrate
    # (rather than a never-fired box-wall safety net), and per-step impulse
    # delivered into z = 0 is accumulated in `diag_substrate_impulse_z`
    # for the anchor-force-balance gate. See
    # `docs/stage1a_plus_substrate_sanity.md` and `docs/outcomes_stage1a_plus.md`.
    # When False, the −z wall is the same 3-cell safety net as the other five
    # box faces (Stage 1a baseline behaviour, no substrate).
    substrate_enabled: bool = False
    n_contact_band: int = 3
    # Stage 1a+ Option β: substrate adhesion energy density `γ_sub_star`
    # (dimensionless, in K·R₀ units = stress·length). When > 0 and
    # `substrate_enabled` is True, the substrate CSF impulse
    #   dv_z = −γ_sub_star · κ_sub_proxy · dt / ρ_local
    # (with κ_sub_proxy = 1/dx*, n̂_sub = −ẑ) is applied to grid cells in
    # the contact band, pulling particles toward z = 0. Anchored to
    # `Ca_cc` (Maître Science 2012, IF 47) via the runner's
    # `gamma_sub_alpha` sweep variable: γ_sub_star = α · γ_cc_star =
    # α · Ca_cc · K_star · radius_star. See
    # `docs/stage1a_plus_substrate_sanity.md` §"Option β addendum" and
    # `docs/outcomes_stage1a_plus.md` §"Option β addendum". Default 0
    # recovers Option α (mechanical anchor only).
    gamma_sub_star: float = 0.0

    # Stage 1a++ Layer 2 boundary-cell active stress (continuum entry; no
    # lamellipodia / filopodia / leader / FA stochastic events here, those
    # are Stage 1a++.b deferred). When `zeta_star` > 0, every particle with
    # `is_boundary[p] == 1` receives an additional Cauchy stress
    #   σ_act_p = -ζ_star · K_star · I
    # (sign: contractile cortex; compressive stress opposes the surface-
    # tension contraction pulling the spheroid inward). Reported as a
    # *dimensionless ratio* per PI Option α' resolution
    # 2026-04-29 — no Pa claim is made; the result is a response curve in
    # ζ/K. K is anchored to Fischer-Friedrich Nat Cell Biol 2014 IF 30
    # (already cited in `docs/02_force_models.md` §1.1); Marchetti
    # Rev Mod Phys 2013 IF 50 retained as the framework reference. See
    # `docs/stage1a_plus_plus_layer2_sanity.md` and
    # `docs/outcomes_stage1a_plus_plus.md`. Default 0 reproduces Stage 1a+
    # Option β (carrier baseline).
    zeta_star: float = 0.0

    # Stage 1b Layer 3 φ-ODE (Cho et al. 2020 mechanism, see
    # `docs/03_adhesion_dynamics.md` and `docs/stage1b_layer3_sanity.md`).
    # When `layer3_enabled` is True, every particle carries a φ ∈ [0, 1]
    # state evolved by the simplified ODE (S ≡ 1 since substrate is single-
    # condition Col1):
    #   dφ/dt = k_+_star · (1 − φ) − k_-_star · φ
    # and the Layer 2 coupling becomes φ-dependent:
    #   ζ_star_per_p = ζ_min · (1 − φ_p) + ζ_max · φ_p
    # replacing the constant `zeta_star` (which is ignored when Layer 3
    # is on). PI Option α'-style framing — bounds anchored to Stage 1a++
    # Track 1 stable regime; φ_initial mapped per Bare/Pre/Lam4 phenotype
    # by PI domain expertise (PI full authorisation 2026-04-29). All
    # values default to disable Layer 3 (k_+ = k_- = 0; ζ_min = ζ_max =
    # zeta_star) so existing Stage 1a / 1a+ / 1a++ configs are unaffected.
    layer3_enabled: bool = False
    phi_initial: float = 0.05
    k_plus_star: float = 0.0
    k_minus_star: float = 0.0
    zeta_min: float = 0.0
    zeta_max: float = 0.0

    # Path C effective gravity (buoyancy-corrected). Per
    # `docs/path_c_sanity.md` and `docs/outcomes_path_c.md`. When > 0,
    # every grid cell with mass receives a per-step downward velocity
    # impulse `Δv_z = −gravity_star · dt`, applied in
    # `_grid_op_overdamped` after CSF / before drag / before reflective
    # wall clamp. Default 0 disables Path C (Stage 1a+ Option α etc.
    # behaviour preserved). Anchored to Stewart Nature 2011 cell-density
    # framework (IF 65, already cited in `docs/02_force_models.md` §1.7);
    # specific value framed as a body-force coefficient (PARTIAL
    # Magic-Number Block per ζ_star Option α' precedent — see sanity-md).
    # Provisional default for the v15-Path-C baseline pilot: 0.1.
    gravity_star: float = 0.0

    # Stage 1c Layer 5 mechano-osmotic Tier 2 (per
    # `docs/08_mechano_osmotic.md` framework citing Guo PNAS 2017 IF 12 +
    # Venkova eLife 2022). When `layer5_enabled` is True, every particle
    # carries a `rho_osm_p ∈ [rho_osm_min, rho_osm_max]` osmotic state
    # evolved by forward Euler:
    #   dρ_osm/dt = α_osm · max(0, -tr(C_p)) − β_osm · (ρ_osm − 1)
    # ε̇^spreading = max(0, -tr(C)): positive when cell volume contracts,
    # which is when water effluxes out of the cell. Coupling to the
    # bulk modulus: K_eff = K_star · ρ_osm_p (replaces constant K in the
    # v15 σ_vol formula). All values default to disable Layer 5
    # (α_osm = β_osm = 0; ρ_osm fixed at rho_osm_initial = 1.0). PI full
    # authorisation 2026-04-29 ("원래 framework standard 풀 적용").
    # Magic-Number Block PARTIAL on the ODE constants — analogous to ζ_star
    # Option α' precedent; honest disclosure in `docs/stage1c_sanity.md`.
    layer5_enabled: bool = False
    rho_osm_initial: float = 1.0
    alpha_osm_star: float = 0.0
    beta_osm_star: float = 0.0
    rho_osm_min: float = 0.5
    rho_osm_max: float = 1.6

    # Stage 1d Layer 4 cellular Marangoni (per
    # `docs/07_internal_flow_dynamics.md` framework citing Pajic-Lijakovic
    # & Milivojevic Eur Biophys J 2022 + Fütterer Phys Rev Fluids 2022).
    # When `layer4_enabled` is True:
    #   γ(φ_p) = γ_max·(1−φ_p) + γ_min·φ_p   per particle (Maître IF 47
    #     anchor for γ_max/γ_min range; γ_1 ≡ γ_min − γ_max < 0 per
    #     Cho 2020 mechanism)
    #   γ scattered to grid via mass-weighted P2G → grid_gamma
    #   ∇γ via central differences → grid_gamma_grad
    #   Tangent projection ∇_s γ = (I − n̂⊗n̂)·∇γ using existing grid_normal
    #   Marangoni impulse Δv = ∇_s γ · dt / ρ_local at boundary cells
    # γ_max_star, γ_min_star are dimensionless (γ / (K · R₀)). Default 0
    # disables Layer 4. PI full authorisation 2026-04-29 covers PARTIAL
    # Magic-Number Block per ζ_star Option α' precedent (Pajic-Lijakovic
    # 2022 framework anchored, Maître IF 47 γ-range anchored, specific γ_max
    # / γ_min values inherit project framework defaults from
    # `03_adhesion_dynamics.md`).
    layer4_enabled: bool = False
    gamma_max_star: float = 0.0
    gamma_min_star: float = 0.0
    # Stage 1d also enables a Layer 3 spatial extension: per-particle
    # substrate-contact indicator S_p (binary; 1 if z_p < n_contact_band·dx,
    # 0 otherwise). When True, Layer 3 ODE differentiates contact-band
    # particles (φ→φ_eq) from interior (φ→0), creating ∇φ → ∇γ → Marangoni
    # driving force. Without this, single-phenotype Layer 3 has uniform φ
    # → ∇γ ≈ 0 → no Marangoni effect.
    #
    # Stage 1b.b (PI directive 2026-04-29) DEPRECATES this flag: the
    # φ_memory + c_act split intrinsically encodes the contact-band
    # spatial gating in the c_act ODE while preserving formation
    # phenotype memory in φ_memory. The flag is read but ignored;
    # `layer3_split` (default True) is the new control.
    layer3_spatial_S: bool = False

    # Stage 1b.b Layer 3 φ split (per `docs/layer3_phi_audit.md` §4 +
    # `docs/stage1b_b_phi_split_sanity.md`). When True (default if
    # layer3_enabled), φ is split into φ_memory_p (formation phenotype
    # memory, ε-decay toward phi_initial) and c_act_p (contact
    # activation, Cho 2020 ODE on contact-band only). Effective
    # variable for downstream (Layer 4 γ): φ_eff = φ_memory + κ·c_act·
    # (1 − φ_memory). When False, falls back to legacy single-φ ODE.
    layer3_split: bool = True
    layer3_kappa_act: float = 1.0      # κ in φ_eff formula (PI default 1.0)
    layer3_memory_eps_star: float = 0.0  # ε memory decay rate (PI default 0.0 = exactly fixed)

    # Stage 2 Layer 6 chemistry / ECM remodeling — minimal scope per
    # `docs/stage2_sanity.md`: MMP secretion + ECM degradation only,
    # de novo ECM secretion deferred to Stage 2.b.
    # Per-step ODE (forward Euler):
    #   d(mmp_total)/dt = α_MMP · n_contact_band(t)
    #   d(ecm_strength)/dt = -β_deg · mmp_total · ecm_strength
    # γ_sub_eff(t) = γ_sub_star · ecm_strength(t) — substrate CSF
    # impulse magnitude scaled by current ecm_strength. ecm_strength
    # clamped to [ecm_strength_min, 1.0] (defensive floor against
    # runaway full degradation). Default 0 disables Layer 6.
    # Magic-Number Block PARTIAL anchored to Egeblad-Werb 2002
    # Nat Rev Cancer IF 70 + Lu 2011 Nat Rev Mol Cell Biol IF 113
    # framework references; specific dimensionless rates are order-of-
    # magnitude derivations matching the project's τ_relax = 60 s
    # calibration. PI full authorisation 2026-04-29.
    layer6_enabled: bool = False
    alpha_mmp_star: float = 0.0
    beta_deg_star: float = 0.0
    ecm_strength_initial: float = 1.0
    ecm_strength_min: float = 0.1

    @property
    def dx_star(self) -> float:
        return self.domain_star / self.grid_n

    @property
    def gamma_star(self) -> float:
        return self.capillary_number * self.K_star * self.radius_star

    @property
    def particle_volume_star(self) -> float:
        sphere_vol = (4.0 / 3.0) * np.pi * self.radius_star**3
        return sphere_vol / self.n_particles

    @property
    def particle_mass_star(self) -> float:
        return self.density_star * self.particle_volume_star


@ti.data_oriented
class MLSMPMSolver:
    """Overdamped MLS-MPM in dimensionless units, vacuum boundary, Stage 1a Layer 1."""

    def __init__(self, cfg: SolverConfig):
        if cfg.n_particles < 8:
            raise ValueError("N must be ≥ 8 for a single APIC stencil to be defined.")
        if cfg.radius_star < 2.0 * cfg.dx_star:
            raise ValueError(
                f"radius* ({cfg.radius_star}) must be ≥ 2·dx* ({2.0 * cfg.dx_star}); "
                "increase grid_n or shrink the domain."
            )
        # Stage 1b Layer 3 φ-ODE forward-Euler stiffness invariant
        # (per docs/stage1b_layer3_sanity.md check 2): factor-2 safety
        # margin against `dt · (k_+ + k_-) ≤ 1`.
        if cfg.layer3_enabled:
            stiffness = cfg.dt_star * (cfg.k_plus_star + cfg.k_minus_star)
            if stiffness > 0.5:
                raise ValueError(
                    f"Layer 3 φ-ODE stiffness violation: dt·(k_+ + k_-) = "
                    f"{stiffness:.4f} > 0.5. Reduce dt_star or k rates "
                    f"(forward Euler stability bound)."
                )
            if not (0.0 <= cfg.phi_initial <= 1.0):
                raise ValueError(
                    f"phi_initial ({cfg.phi_initial}) must be ∈ [0, 1]."
                )
        # Stage 1c Layer 5 forward-Euler stability invariant
        # (per docs/stage1c_sanity.md check 2): dt · β_osm ≤ 0.5
        # (factor-2 safety; α-driver bounded by typical |tr(C)| ≤ 1).
        if cfg.layer5_enabled:
            stiffness_osm = cfg.dt_star * cfg.beta_osm_star
            if stiffness_osm > 0.5:
                raise ValueError(
                    f"Layer 5 ρ_osm-ODE stiffness violation: dt·β_osm = "
                    f"{stiffness_osm:.4f} > 0.5. Reduce dt_star or β_osm "
                    f"(forward Euler stability bound)."
                )
            if not (cfg.rho_osm_min <= cfg.rho_osm_initial <= cfg.rho_osm_max):
                raise ValueError(
                    f"rho_osm_initial ({cfg.rho_osm_initial}) must be ∈ "
                    f"[{cfg.rho_osm_min}, {cfg.rho_osm_max}]."
                )

        self.cfg = cfg
        n_p, n_g = cfg.n_particles, cfg.grid_n

        self.x = ti.Vector.field(3, dtype=ti.f32, shape=n_p)
        self.v = ti.Vector.field(3, dtype=ti.f32, shape=n_p)
        self.C = ti.Matrix.field(3, 3, dtype=ti.f32, shape=n_p)
        self.F = ti.Matrix.field(3, 3, dtype=ti.f32, shape=n_p)
        self.tau_dev = ti.Matrix.field(3, 3, dtype=ti.f32, shape=n_p)
        self.is_boundary = ti.field(dtype=ti.i32, shape=n_p)

        self.grid_v = ti.Vector.field(3, dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_m = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_count = ti.field(dtype=ti.i32, shape=(n_g, n_g, n_g))
        # Reproducing-kernel volume sum (Adami, Hu & Adams 2010 §3 / Shepard
        # normalisation): grid_kernel_weight[I] = Σ_p w_pI · V₀, the
        # kernel-volume coverage at grid cell I. Used to normalise the
        # colour function so c_norm = grid_m / (ρ_bulk · grid_kernel_weight)
        # is identically 1 in bulk wherever grid_kernel_weight > 0
        # (partition-of-unity), eliminating the random-pack radial trend
        # that produced CSF interior penetration in v11.
        self.grid_kernel_weight = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_color = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        # Smoothed colour for Brackbill (1992) §V: a thin diffuse interface
        # (~1 grid cell) makes ∇c picks up grid-scale variation rather than
        # physical surface curvature. Smoothing the colour with a few passes
        # of a 3³ averaging kernel widens the interface to ~5–7 dx and lets
        # the FD curvature operator resolve 2/R rather than 1/dx.
        self.grid_color_smooth = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_color_grad = ti.Vector.field(3, dtype=ti.f32, shape=(n_g, n_g, n_g))
        # Brackbill (1992) curvature operator fields:
        #   n̂[I] = ∇c[I] / sqrt(|∇c|² + ε²)   (smoothed unit normal, ε regularises f32)
        #   κ[I] = -∇·n̂[I]                   (mean curvature; positive for convex outward surface)
        self.grid_normal = ti.Vector.field(3, dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_kappa = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        # Stage 1d Layer 4 Marangoni grid scratch fields. grid_gamma is the
        # mass-weighted γ_eff per cell; grid_gamma_grad is its central-FD
        # gradient. Same allocation pattern as grid_color / grid_color_grad.
        self.grid_gamma = ti.field(dtype=ti.f32, shape=(n_g, n_g, n_g))
        self.grid_gamma_grad = ti.Vector.field(3, dtype=ti.f32, shape=(n_g, n_g, n_g))

        # Reference-calibration scratch fields (allocated once, reused on every
        # call to `calibrate_reference_state`).
        self._calib_rho = ti.field(dtype=ti.f32, shape=n_p)
        self._calib_W = ti.field(dtype=ti.f32, shape=n_p)
        self._calib_scale = ti.field(dtype=ti.f32, shape=n_p)

        # v15 (k.3) density-based volumetric stress fields.
        # `_rho_kernel_p` is the kernel-interpolated density at particle p,
        # refreshed at the start of every `step()` from the post-`_p2g_mass`
        # `grid_m`. `_rho_ref_kernel` and `_rho_floor` are scalar fields set
        # once by `calibrate_reference_state` (rho_ref = harmonic mean over
        # well-resolved particles; rho_floor = 0.1·rho_ref). See
        # `docs/stage1a_interior_pressure_sanity.md` for the proposal and
        # `docs/outcomes_v15.md` for the literature-verification record on
        # the 0.1 floor value.
        self._rho_kernel_p = ti.field(dtype=ti.f32, shape=n_p)
        self._rho_ref_kernel = ti.field(dtype=ti.f32, shape=())
        self._rho_floor = ti.field(dtype=ti.f32, shape=())
        # Stage 1b Layer 3 φ field declared early so the constructor
        # initialiser block below can populate it together with the
        # density fields. Diagnostic accumulators are declared later.
        self.phi_p = ti.field(dtype=ti.f32, shape=n_p)
        # Stage 1b.b Layer 3 φ split (per docs/layer3_phi_audit.md §4 +
        # docs/stage1b_b_phi_split_sanity.md). φ_memory_p stores the
        # formation phenotype memory (slow/fixed); c_act_p stores the
        # contact-activation state (Cho 2020 ODE on contact-band only).
        # Effective variable phi_eff_p = φ_memory + κ · c_act ·
        # (1 − φ_memory) is computed each step and assigned into
        # `phi_p` for downstream consumers (Layer 4 γ scatter, Layer 2
        # ζ_eff, diagnostics). When `layer3_split=False`, these are
        # zero-init'd and unused.
        self.phi_memory_p = ti.field(dtype=ti.f32, shape=n_p)
        self.c_act_p = ti.field(dtype=ti.f32, shape=n_p)
        # Stage 1c Layer 5 osmotic state ρ_osm field (per particle).
        # Initialised to `cfg.rho_osm_initial` in the constructor. Evolved
        # by `_integrate_osmotic_ode` once per step when Layer 5 enabled.
        self.rho_osm_p = ti.field(dtype=ti.f32, shape=n_p)
        # Stage 2 Layer 6 scalar state fields (global; not per-particle).
        # mmp_total accumulates MMP secretion over time; ecm_strength
        # decays exponentially as MMP·ecm_strength accumulates.
        # diag_n_contact_band_layer6 is a per-step accumulator for the
        # contact-band particle count needed by the Layer 6 ODE
        # (separate from the diagnostic-frame `diag_n_boundary` which
        # only updates per `_compute_invariants` call).
        self.mmp_total_field = ti.field(dtype=ti.f64, shape=())
        self.ecm_strength_field = ti.field(dtype=ti.f32, shape=())
        self.diag_n_contact_band_layer6 = ti.field(dtype=ti.i32, shape=())

        # Initialise v15 density fields to a self-consistent default
        # (σ_vol = K·(ρ_ref/ρ_kernel − 1) = 0 for ρ_kernel = ρ_ref = density_star)
        # so that any code path that calls `step()` or `invariants()` without
        # first calling `calibrate_reference_state` has well-defined behaviour.
        # `calibrate_reference_state` later overwrites all three with the
        # harmonic-mean-derived values from the actual particle pack.
        rho_default = float(cfg.density_star)
        self._rho_ref_kernel[None] = rho_default
        self._rho_floor[None] = 0.1 * rho_default
        self._rho_kernel_p.from_numpy(
            np.full(cfg.n_particles, rho_default, dtype=np.float32)
        )
        # Stage 1b Layer 3 φ field initialised to `phi_initial` (uniform
        # across all particles; per-phenotype mapping is set by the runner
        # via the YAML `layer3.phi_initial` field).
        self.phi_p.from_numpy(
            np.full(cfg.n_particles, float(cfg.phi_initial), dtype=np.float32)
        )
        # Stage 1b.b φ split initial state: φ_memory := phi_initial,
        # c_act := 0 (no contact activation at t=0). Per
        # `docs/stage1b_b_phi_split_sanity.md` §"Initialization".
        self.phi_memory_p.from_numpy(
            np.full(cfg.n_particles, float(cfg.phi_initial), dtype=np.float32)
        )
        self.c_act_p.from_numpy(
            np.zeros(cfg.n_particles, dtype=np.float32)
        )
        # Stage 1c Layer 5 osmotic state ρ_osm initialised to
        # `rho_osm_initial` (= 1.0 = full hydration by default).
        self.rho_osm_p.from_numpy(
            np.full(cfg.n_particles, float(cfg.rho_osm_initial), dtype=np.float32)
        )
        # Stage 2 Layer 6 scalar state initialisation.
        self.mmp_total_field[None] = 0.0
        self.ecm_strength_field[None] = float(cfg.ecm_strength_initial)

        self.diag_mass = ti.field(dtype=ti.f64, shape=())
        self.diag_momentum = ti.Vector.field(3, dtype=ti.f64, shape=())
        self.diag_kinetic_energy = ti.field(dtype=ti.f64, shape=())
        self.diag_strain_energy = ti.field(dtype=ti.f64, shape=())
        self.diag_surface_energy = ti.field(dtype=ti.f64, shape=())
        self.diag_max_speed = ti.field(dtype=ti.f32, shape=())
        self.diag_max_count = ti.field(dtype=ti.i32, shape=())
        self.diag_nan_count = ti.field(dtype=ti.i32, shape=())

        # Stage 1a+ Option α: per-step accumulator of substrate reaction
        # impulse (in z direction). Reset in `_clear_grid`, accumulated in
        # `_grid_op_overdamped` whenever the −z wall clamps a downward grid
        # velocity, and read by `substrate_diagnostics` to compute
        # F_substrate = diag_substrate_impulse_z / dt for the anchor-force-
        # balance gate.
        self.diag_substrate_impulse_z = ti.field(dtype=ti.f64, shape=())

        # Stage 1b Layer 3 φ field is declared earlier (next to
        # `_rho_kernel_p`) so the constructor's `from_numpy` initialiser
        # block can populate it before `_zero_state` runs.
        # Path C gravitational potential energy (sum over all particles
        # of `ρ_p · g_star · z_p · V₀`). Reset per `_compute_invariants`
        # call. Reported in `invariants()` for the energy-monotone gate
        # extension.
        self.diag_grav_pe = ti.field(dtype=ti.f64, shape=())
        # Path C centre-of-mass z (informational sedimentation depth).
        self.diag_com_z_sum = ti.field(dtype=ti.f64, shape=())

        # Stage 1c Layer 5 ρ_osm aggregation diagnostics: bulk-shell
        # <ρ_osm> via sum / N (host-side); min/max for the per-particle
        # invariant gate. Sum-of-squares for variance.
        self.diag_rho_osm_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_rho_osm_sum_sq = ti.field(dtype=ti.f64, shape=())
        self.diag_rho_osm_min = ti.field(dtype=ti.f32, shape=())
        self.diag_rho_osm_max = ti.field(dtype=ti.f32, shape=())

        # Stage 1b Layer 3 diagnostics: bulk-shell <φ>, boundary-shell <φ>,
        # contact-band <φ>, plus min/max for the φ ∈ [0,1] invariant gate.
        self.diag_phi_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_phi_sum_sq = ti.field(dtype=ti.f64, shape=())
        self.diag_phi_min = ti.field(dtype=ti.f32, shape=())
        self.diag_phi_max = ti.field(dtype=ti.f32, shape=())
        self.diag_phi_boundary_sum = ti.field(dtype=ti.f64, shape=())
        # Stage 1b.b Layer 3 split diagnostics: per-frame aggregate of
        # phi_memory and c_act for new gates 5a + 5b.
        # φ_memory: global mean for gate 5a preservation check.
        # c_act: contact-band mean for gate 5b trajectory check; plus
        # global min/max for the [0, 1] invariant.
        self.diag_phi_memory_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_phi_memory_min = ti.field(dtype=ti.f32, shape=())
        self.diag_phi_memory_max = ti.field(dtype=ti.f32, shape=())
        self.diag_c_act_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_c_act_min = ti.field(dtype=ti.f32, shape=())
        self.diag_c_act_max = ti.field(dtype=ti.f32, shape=())
        self.diag_c_act_band_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_c_act_band_count = ti.field(dtype=ti.i32, shape=())

        # Stage 1a++ Layer 2 diagnostics:
        #  • diag_active_power: instantaneous power delivered by the active
        #    boundary stress, P_act = Σ_{p ∈ boundary} -ζ·K · tr(C_p) · V₀.
        #    Positive when the cortex does positive work on the bulk (cortex
        #    expanding) — under contractile cortex (ζ > 0) this is normally
        #    negative (cortex compresses, work flows out via overdamped drag).
        #  • diag_n_boundary: count of boundary-tagged particles for the
        #    boundary-tag-flicker stability gate (host-side comparison
        #    across frames).
        self.diag_active_power = ti.field(dtype=ti.f64, shape=())
        self.diag_n_boundary = ti.field(dtype=ti.i32, shape=())

    # ------------------------------------------------------------- init ----
    def initialize_sphere(self, center_star, radius_star: float | None = None) -> None:
        """Random pack inside a sphere using rejection sampling."""
        radius_star = self.cfg.radius_star if radius_star is None else radius_star
        rng = np.random.default_rng(self.cfg.seed)
        pts = np.empty((self.cfg.n_particles, 3), dtype=np.float32)
        n = 0
        cx = np.asarray(center_star, dtype=np.float32)
        while n < self.cfg.n_particles:
            batch = rng.uniform(-1.0, 1.0, size=(self.cfg.n_particles * 2, 3)).astype(np.float32)
            inside = (batch[:, 0] ** 2 + batch[:, 1] ** 2 + batch[:, 2] ** 2) <= 1.0
            keep = batch[inside][: self.cfg.n_particles - n]
            pts[n : n + len(keep)] = keep * radius_star + cx
            n += len(keep)
        self.x.from_numpy(pts)
        self._zero_state()

    @ti.kernel
    def _zero_state(self):
        for p in self.v:
            self.v[p] = ti.Vector.zero(ti.f32, 3)
            self.C[p] = ti.Matrix.zero(ti.f32, 3, 3)
            self.F[p] = ti.Matrix.identity(ti.f32, 3)
            self.tau_dev[p] = ti.Matrix.zero(ti.f32, 3, 3)
            self.is_boundary[p] = 0

    # ----------------------------------------------- reference calibration ----
    @ti.kernel
    def _scatter_mass_only(self):
        """One P2G pass that scatters mass and the unit-volume kernel weight
        (no momentum, no stress) to the grid.

        Used by `calibrate_reference_state` to measure local density and by
        `measure_surface_curvature` (read-only static measurement). The
        kernel-volume sum `Σ_p w_pI · V₀` enables the Adami-Hu-Adams 2010 §3
        reproducing-kernel normalisation of the colour function: in any bulk
        cell where particles fully populate the kernel support, the ratio
        `grid_m / (ρ_bulk · grid_kernel_weight) = m_p / (ρ_bulk · V₀) = 1`
        identically, with NO random-pack noise.
        """
        m_p = self.cfg.particle_mass_star
        V0 = self.cfg.particle_volume_star
        dx = self.cfg.dx_star
        for I in ti.grouped(self.grid_m):
            self.grid_m[I] = 0.0
            self.grid_kernel_weight[I] = 0.0
        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    ti.atomic_add(self.grid_m[idx], weight * m_p)
                    ti.atomic_add(self.grid_kernel_weight[idx], weight * V0)

    @ti.kernel
    def _interpolate_density_to_particles(self):
        """G2P-style interpolation of `grid_m / dx³` (= local kernel density)
        into `_calib_rho`, and of `grid_kernel_weight / dx³` (= dimensionless
        kernel coverage, ≈ 1 in fully-resolved bulk) into `_calib_W`.

        `_calib_rho` retains the v11 semantics for backwards-compatible
        well-resolved selection. `_calib_W` is a diagnostic for the
        v12 task-7 consistency check (does W-based well-resolved coincide
        with ρ-based?).
        """
        dx = self.cfg.dx_star
        inv_vol = 1.0 / (dx ** 3)
        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]
            rho = 0.0
            wsum = 0.0
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    rho += weight * self.grid_m[idx] * inv_vol
                    wsum += weight * self.grid_kernel_weight[idx] * inv_vol
            self._calib_rho[p] = rho
            self._calib_W[p] = wsum

    @ti.kernel
    def _set_F_isotropic_from_calib(self):
        """Set F[p] = _calib_scale[p] · I per particle."""
        for p in self.x:
            self.F[p] = self._calib_scale[p] * ti.Matrix.identity(ti.f32, 3)

    def calibrate_reference_state(self) -> dict:
        """Population-aware reference-state calibration.

        Per Hu et al. 2018 §4.3 / Jiang et al. 2015 for the F-rescaling itself,
        and Adami, Hu & Adams 2010 (J. Comp. Phys. 229, 5011) for the
        population partition rationale.

        Random rejection sampling produces a particle pack whose
        kernel-interpolated density ρ_kernel differs from the nominal
        `density_star`. **The deviation has two distinct physical origins:**

          (i) Bulk packing noise: interior particles see a complete 3³ kernel
              and ρ_kernel fluctuates around the true bulk density due to the
              discrete random pack. This IS a real elastic mismatch and we
              want to absorb it into a non-trivial reference state F[p] ≠ I.

         (ii) Kernel-truncation artifact at the free surface: boundary
              particles see a kernel that is partially in vacuum, so
              ρ_kernel is systematically biased low (by up to ~50% in 3D)
              relative to the true cellular density. This is NOT a real
              elastic mismatch — it is a numerical artifact of the SPH-style
              kernel-density estimator (Adami, Hu & Adams 2010).

        A naive uniform calibration absorbs the kernel artifact (ii) into
        ρ_ref and biases the volumetric stress for ALL particles in a
        grid-resolution-dependent way. This was the v10 failure mode:
        ρ_ref was the harmonic mean over a "well-resolved" subset but the
        F-scale was applied to ALL particles, giving <J>_all = 1.556 ≠ 1
        (population-mismatch Jensen residual).

        Population-aware fix:
          - identify well-resolved particles by ρ_kernel > 0.5·ρ_max (interior
            kernel sees a roughly complete neighbourhood),
          - compute ρ_ref = harmonic mean over THIS subset only,
          - apply F[p] = (ρ_ref / ρ_p)^(1/3) · I to THIS subset only,
          - leave F[p] = I for boundary particles (set by _zero_state). Their
            σ_vol(t=0) = K·(J−1) = 0 by construction, contributing zero
            spurious elastic strain energy. The kernel-truncation bias stays
            where it physically belongs — in the kernel — instead of being
            absorbed into a bulk reference parameter.

        By construction <J>_well_resolved ≡ 1.0 (the harmonic-mean identity
        holds because ρ_ref and the F-scale-application population coincide),
        and <J>_boundary ≡ 1.0 (F = I trivially). The runner asserts
        |<J>_well_resolved − 1| < 1e-6 as a hard scheme-correctness gate.

        Returns a diagnostics dict for the runner to log and gate.
        """
        self._scatter_mass_only()
        self._interpolate_density_to_particles()
        rho_np = self._calib_rho.to_numpy()

        rho_max = float(rho_np.max())
        rho_arith_mean = float(rho_np.mean())
        well_resolved_mask = rho_np > 0.5 * rho_max
        n_resolved = int(well_resolved_mask.sum())

        # Default F-scale = 1 everywhere (i.e. F = I); only well-resolved
        # particles are modified.
        scale_np = np.ones_like(rho_np, dtype=np.float32)

        if n_resolved > 0:
            rho_used = rho_np[well_resolved_mask]
            # Harmonic mean over the SAME population to which the F-scale will
            # be applied — this is the only way <J>_population = 1 exactly
            # (closes the Jensen identity: ρ_ref · <1/ρ_p>_resolved ≡ 1).
            inv_rho = 1.0 / np.clip(rho_used, 1e-6, None)
            rho_ref = float(1.0 / inv_rho.mean())
            ratio = rho_ref / np.clip(rho_used, 1e-6, None)
            # Clamp guards pathological outliers; for any sane pack the
            # well-resolved ratio sits comfortably inside [1/8, 8].
            scale_np[well_resolved_mask] = np.cbrt(
                np.clip(ratio, 1.0 / 8.0, 8.0)
            ).astype(np.float32)
        else:
            # Degenerate (no well-resolved particles): leave F = I everywhere
            # and report ρ_ref = arithmetic mean for the log only.
            rho_ref = rho_arith_mean

        self._calib_scale.from_numpy(scale_np)
        self._set_F_isotropic_from_calib()

        # v15 (k.3): record ρ_ref_kernel and the numerical-safety floor
        # ρ_floor = 0.1·ρ_ref for the density-based volumetric stress, and
        # initialise `_rho_kernel_p` to ρ_ref so the very first step's σ_vol
        # is zero by construction (it is overwritten anyway by
        # `_interpolate_rho_runtime` at the top of step 1, before
        # `_p2g_momentum_and_stress` reads it; this initialisation keeps the
        # field in a defined state for any out-of-loop diagnostics that read
        # it before `step` is called).
        self._rho_ref_kernel[None] = float(rho_ref)
        self._rho_floor[None] = 0.1 * float(rho_ref)
        self._rho_kernel_p.from_numpy(
            np.full(self.cfg.n_particles, float(rho_ref), dtype=np.float32)
        )

        # Diagnostic: <J> over the three populations.
        J_per_p = scale_np.astype(np.float64) ** 3
        boundary_mask = ~well_resolved_mask
        n_boundary = int(boundary_mask.sum())
        J_well_resolved = float(J_per_p[well_resolved_mask].mean()) if n_resolved else 1.0
        J_boundary = float(J_per_p[boundary_mask].mean()) if n_boundary else 1.0
        J_all = float(J_per_p.mean())

        # Task-7 consistency check: does the W-based (kernel-coverage)
        # well-resolved subset coincide with the ρ-based one we use for
        # calibration? In Adami-Hu-Adams 2010 §3 framing the kernel-coverage
        # W_p ≈ Σ_p w_pI · V₀ / dx³ is the cleaner indicator of "kernel sees
        # full neighbourhood" — bulk-noise-immune, surface-truncation-direct.
        # We *report* the W-based mask alongside the ρ-based one. If they
        # agree (high Jaccard), the existing calibration logic is consistent;
        # if not, the v12 commit should record the discrepancy and motivate
        # a unified definition in v13.
        W_np = self._calib_W.to_numpy()
        W_max = float(W_np.max())
        W_well_resolved_mask = W_np > 0.5 * W_max if W_max > 0 else np.zeros_like(W_np, dtype=bool)
        n_W_resolved = int(W_well_resolved_mask.sum())
        intersection = int(np.logical_and(well_resolved_mask, W_well_resolved_mask).sum())
        union = int(np.logical_or(well_resolved_mask, W_well_resolved_mask).sum())
        jaccard = float(intersection / union) if union > 0 else 1.0

        # Reset grid mass so the next real step starts from a clean slate.
        self._clear_grid()
        return {
            "rho_ref_harmonic": rho_ref,
            "rho_arith_mean": rho_arith_mean,
            "rho_actual_min": float(rho_np.min()),
            "rho_actual_max": rho_max,
            "W_kernel_min": float(W_np.min()),
            "W_kernel_max": W_max,
            "W_kernel_mean": float(W_np.mean()),
            "n_W_well_resolved": n_W_resolved,
            "rho_W_jaccard": jaccard,
            "F_scale_min": float(scale_np.min()),
            "F_scale_max": float(scale_np.max()),
            "F_scale_mean": float(scale_np.mean()),
            "J_mean_well_resolved": J_well_resolved,
            "J_mean_boundary_subset": J_boundary,
            "J_mean_all": J_all,
            # Backwards-compatible alias for the runner log line.
            "J_mean_after_calib": J_all,
            "n_well_resolved": n_resolved,
            "n_boundary_subset": n_boundary,
            "calibration_population": (
                "well_resolved (ρ_kernel > 0.5·ρ_max); boundary subset retains F = I"
            ),
        }

    # ----------------------------------------------------- Layer 3 φ-ODE ----
    @ti.kernel
    def _integrate_phi_ode(self):
        """Stage 1b Layer 3 forward-Euler φ-ODE integration with Stage 1d
        spatial S_p extension.

        Per `docs/03_adhesion_dynamics.md` and `docs/stage1b_layer3_sanity.md`:

            dφ/dt = k_+ · S_p · (1 − φ) − k_- · φ

        Stage 1b default: S_p ≡ 1 (substrate-induced signal uniform; single-
        condition Col1 always present from t=0). Stage 1d extension
        (`layer3_spatial_S`): S_p = 1 if z_p < n_contact_band·dx, else 0
        — substrate-engaged particles see S=1 (φ → φ_eq = 0.75); interior
        particles see S=0 (φ → 0 via the k_- relaxation alone). This
        creates ∇φ → ∇γ → Marangoni driving force. Required for Stage 1d
        (Layer 4) since otherwise ∇γ ≈ 0 with single-phenotype.

        Forward-Euler update with stiffness checked at constructor time
        (`dt · (k_+ + k_-) ≤ 0.5`). The `clip(0, 1)` defends against any
        per-step overshoot at the boundaries of the [0, 1] domain.

        Stage 1b.b (PI directive 2026-04-29 per docs/layer3_phi_audit.md):
        when `layer3_split=True` (default for layer3_enabled), φ is
        split into φ_memory_p (formation phenotype memory, ε-decay
        toward phi_initial) and c_act_p (contact activation, Cho 2020
        ODE on contact-band only with c_act preserved in interior on
        migration). φ_eff_p = φ_memory + κ · c_act · (1 − φ_memory) is
        assigned into `phi_p` for downstream consumers (Layer 4 γ
        scatter, Layer 2 ζ_eff, diagnostics).

        When `layer3_split=False`, falls back to the legacy single-φ
        ODE (kept for backwards-compat tests).
        """
        dt = self.cfg.dt_star
        k_plus = self.cfg.k_plus_star
        k_minus = self.cfg.k_minus_star
        h_band = self.cfg.n_contact_band * self.cfg.dx_star
        if ti.static(self.cfg.layer3_split):
            # Stage 1b.b new path: φ_memory + c_act split.
            kappa = self.cfg.layer3_kappa_act
            eps_mem = self.cfg.layer3_memory_eps_star
            phi_init_const = self.cfg.phi_initial
            for p in self.phi_p:
                phi_mem = self.phi_memory_p[p]
                c_act = self.c_act_p[p]
                # φ_memory: slow decay toward phi_initial (ε=0 default
                # makes this no-op; the Lyapunov form is robust to
                # any small ε > 0 that might be set by future audits).
                dphi_mem = -dt * eps_mem * (phi_mem - phi_init_const)
                new_phi_mem = phi_mem + dphi_mem
                if new_phi_mem < 0.0:
                    new_phi_mem = 0.0
                if new_phi_mem > 1.0:
                    new_phi_mem = 1.0
                # c_act: substrate-band-only Cho 2020 ODE. In interior
                # (S=0), dc/dt = 0 — c_act *preserved* (not decayed).
                # This contrasts with the v11 spatial S=0 → k_- decay
                # which erased formation memory; here memory is in
                # phi_memory_p, and c_act represents the integrated
                # contact-engagement state.
                S_p = 1.0
                if self.x[p][2] >= h_band:
                    S_p = 0.0
                dc_act = dt * S_p * (k_plus * (1.0 - c_act) - k_minus * c_act)
                new_c_act = c_act + dc_act
                if new_c_act < 0.0:
                    new_c_act = 0.0
                if new_c_act > 1.0:
                    new_c_act = 1.0
                # φ_eff: phenotype-modulated contact-activation.
                phi_eff = new_phi_mem + kappa * new_c_act * (1.0 - new_phi_mem)
                if phi_eff < 0.0:
                    phi_eff = 0.0
                if phi_eff > 1.0:
                    phi_eff = 1.0
                self.phi_memory_p[p] = new_phi_mem
                self.c_act_p[p] = new_c_act
                self.phi_p[p] = phi_eff
        else:
            # Legacy single-φ path (Stage 1b backwards compat).
            for p in self.phi_p:
                phi = self.phi_p[p]
                S_p = 1.0
                if ti.static(self.cfg.layer3_spatial_S):
                    if self.x[p][2] >= h_band:
                        S_p = 0.0
                dphi = dt * (k_plus * S_p * (1.0 - phi) - k_minus * phi)
                new_phi = phi + dphi
                if new_phi < 0.0:
                    new_phi = 0.0
                if new_phi > 1.0:
                    new_phi = 1.0
                self.phi_p[p] = new_phi

    # ------------------------------------------------ Layer 4 Marangoni ----
    @ti.kernel
    def _scatter_gamma_to_grid(self):
        """Stage 1d Layer 4: scatter per-particle γ(φ_p) to grid_gamma.

        Per `docs/07_internal_flow_dynamics.md` §1: γ_eff per cell =
        Σ_p w_pI · m_p · γ(φ_p) / Σ_p w_pI · m_p (mass-weighted average).
        For simplicity and matching the existing CSF colour pattern, we
        scatter the mass-weighted γ contribution; the runner does NOT
        renormalise (the gradient operator is invariant under uniform
        scaling, and absolute γ value is informational only — the
        Marangoni FORCE depends only on ∇γ).

        Assumes _clear_grid has already zeroed grid_gamma. Reads
        per-particle φ_p; computes γ(φ_p) = γ_max·(1−φ) + γ_min·φ
        inline. When Layer 4 disabled this kernel is not called.
        """
        m_p = self.cfg.particle_mass_star
        dx = self.cfg.dx_star
        gamma_max = self.cfg.gamma_max_star
        gamma_min = self.cfg.gamma_min_star
        for p in self.x:
            phi_p = self.phi_p[p]
            gamma_p = gamma_max * (1.0 - phi_p) + gamma_min * phi_p
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    ti.atomic_add(self.grid_gamma[idx], weight * m_p * gamma_p)

    @ti.kernel
    def _compute_gamma_grad(self):
        """Central-difference ∇γ on grid → grid_gamma_grad.

        Same pattern as `_color_gradient_from_smoothed`. The gradient
        is then used in `_grid_op_overdamped` to compute the tangential
        Marangoni impulse via projection with `grid_normal` (already
        populated by `_build_curvature` from CSF).
        """
        dx = self.cfg.dx_star
        inv_2dx = 1.0 / (2.0 * dx)
        n_g = self.cfg.grid_n
        for I in ti.grouped(self.grid_gamma_grad):
            i, j, k = I[0], I[1], I[2]
            gx = 0.0
            gy = 0.0
            gz = 0.0
            if 0 < i < n_g - 1:
                gx = (self.grid_gamma[i + 1, j, k] - self.grid_gamma[i - 1, j, k]) * inv_2dx
            if 0 < j < n_g - 1:
                gy = (self.grid_gamma[i, j + 1, k] - self.grid_gamma[i, j - 1, k]) * inv_2dx
            if 0 < k < n_g - 1:
                gz = (self.grid_gamma[i, j, k + 1] - self.grid_gamma[i, j, k - 1]) * inv_2dx
            self.grid_gamma_grad[I] = ti.Vector([gx, gy, gz])

    def _build_marangoni_field(self) -> None:
        """Build the γ_eff field on grid + its gradient. Called per step
        when Layer 4 is enabled, after `_p2g_mass` (so particles exist on
        grid) and after `_build_curvature` (so grid_normal is populated)."""
        self._scatter_gamma_to_grid()
        self._compute_gamma_grad()

    # ---------------------------------------------- Layer 6 chemistry ----
    @ti.kernel
    def _count_contact_band_layer6(self):
        """Count particles in the substrate contact band (z < n_contact_band·dx).

        Stage 2 Layer 6 needs this per step to update mmp_total ODE.
        Writes the count into `diag_n_contact_band_layer6` field (avoids
        Taichi kernel-return-type compatibility issues across versions).
        """
        h_band = self.cfg.n_contact_band * self.cfg.dx_star
        self.diag_n_contact_band_layer6[None] = 0
        for p in self.x:
            if self.x[p][2] < h_band:
                ti.atomic_add(self.diag_n_contact_band_layer6[None], 1)

    def _integrate_layer6_ode(self) -> None:
        """Stage 2 Layer 6 forward-Euler ODE update (host-side, per step).

        Per `docs/stage2_sanity.md` Tier 2 spec:

            d(mmp_total)/dt   = α_MMP · n_contact_band(t)
            d(ecm_strength)/dt = -β_deg · mmp_total · ecm_strength

        with ecm_strength clamped to [ecm_strength_min, 1.0]. Forward
        Euler with dt = cfg.dt_star. Rates are slow (~ 1e-3 per τ_relax)
        so dt · rate ≪ 1 (massively stable per check 4).
        """
        dt = self.cfg.dt_star
        self._count_contact_band_layer6()
        n_contact = int(self.diag_n_contact_band_layer6[None])

        mmp = float(self.mmp_total_field[None])
        ecm = float(self.ecm_strength_field[None])

        # MMP secretion: monotone-increasing as long as n_contact > 0.
        mmp_new = mmp + dt * self.cfg.alpha_mmp_star * float(n_contact)

        # ECM degradation: exponential decay with rate β·mmp.
        ecm_new = ecm + dt * (-self.cfg.beta_deg_star * mmp * ecm)

        # Clamp ecm_strength to defensive bounds.
        if ecm_new < self.cfg.ecm_strength_min:
            ecm_new = self.cfg.ecm_strength_min
        if ecm_new > 1.0:
            ecm_new = 1.0

        self.mmp_total_field[None] = mmp_new
        self.ecm_strength_field[None] = ecm_new

    # ----------------------------------------------------- Layer 5 ρ_osm ----
    @ti.kernel
    def _integrate_osmotic_ode(self):
        """Stage 1c Layer 5 forward-Euler ρ_osm ODE integration.

        Per `docs/08_mechano_osmotic.md` Tier 2 framework:

            dρ_osm/dt = α_osm · ε̇^spreading − β_osm · (ρ_osm − 1)

        with ε̇^spreading ≈ max(0, -tr(C_p)) (positive when cell volume
        contracts → water leaves → ρ_osm rises). Forward-Euler update with
        stiffness checked at constructor time (`dt · β_osm ≤ 0.5`). The
        clamp `[ρ_osm_min, ρ_osm_max]` defends against transient
        overshoots at the bounds.
        """
        dt = self.cfg.dt_star
        alpha_osm = self.cfg.alpha_osm_star
        beta_osm = self.cfg.beta_osm_star
        rho_osm_min = self.cfg.rho_osm_min
        rho_osm_max = self.cfg.rho_osm_max
        for p in self.rho_osm_p:
            C_p = self.C[p]
            tr_C = C_p[0, 0] + C_p[1, 1] + C_p[2, 2]
            eps_spread = ti.max(0.0, -tr_C)
            rho = self.rho_osm_p[p]
            drho = dt * (alpha_osm * eps_spread - beta_osm * (rho - 1.0))
            new_rho = rho + drho
            if new_rho < rho_osm_min:
                new_rho = rho_osm_min
            if new_rho > rho_osm_max:
                new_rho = rho_osm_max
            self.rho_osm_p[p] = new_rho

    # --------------------------------------------------------- one step ----
    def step(self) -> None:
        """Single MLS-MPM step under v15 (k.3) density-based volumetric stress.

        Pipeline (post-v15):
          1. _clear_grid              — zero grid fields
          2. _tag_boundary            — populate grid_count, classify particles
          3. _p2g_mass                — scatter particle mass + kernel weight
                                        to grid (no momentum, no stress)
          4. _interpolate_rho_runtime — G2P → `_rho_kernel_p` per particle,
                                        with floor clamp ρ_floor = 0.1·ρ_ref
          5. _build_csf_field         — colour, smoothing, ∇c (reads grid_m)
          6. _build_curvature         — n̂ and κ from ∇c
          7. _p2g_momentum_and_stress — scatter momentum + v15 stress
                                        σ = K·(ρ_ref/ρ_kernel_p − 1)·I + τ_dev
          8. _grid_op_overdamped      — CSF impulse, drag, walls, m → v
          9. _g2p_and_constitutive    — gather to particles, advect, update F, τ_dev

        The split P2G (steps 3, 7) is the price of the v15 fix: σ_vol now
        depends on the kernel-interpolated density at p, which itself depends
        on grid_m built from all particles. So mass must scatter first, ρ_kernel
        must be interpolated next, and only then can stress be scattered.
        Estimated step-time cost: ~30–50% over v12.
        """
        self._clear_grid()
        self._tag_boundary()
        self._p2g_mass()                    # populates grid_m + grid_kernel_weight
        self._interpolate_rho_runtime()     # populates _rho_kernel_p (with floor)
        self._build_csf_field()             # reads grid_m to compute colour ρ/ρ_bulk and ∇c
        self._build_curvature()             # n̂ = ∇c/|∇c|;  κ = -∇·n̂
        if self.cfg.layer4_enabled:
            # Stage 1d Layer 4: scatter γ(φ) and compute its gradient on
            # the grid; the Marangoni impulse is applied in
            # `_grid_op_overdamped` using both grid_gamma_grad and the
            # existing grid_normal (from `_build_curvature`).
            self._build_marangoni_field()
        self._p2g_momentum_and_stress()     # scatters momentum + density-based stress
        self._grid_op_overdamped()          # converts to velocity, applies CSF impulse γ·κ·∇c
        self._g2p_and_constitutive()
        if self.cfg.layer3_enabled:
            # Stage 1b Layer 3: integrate the per-particle φ-ODE after the
            # mechanical step (φ is decoupled from positions; the ζ(φ)
            # coupling enters next step's `_p2g_momentum_and_stress`).
            self._integrate_phi_ode()
        if self.cfg.layer5_enabled:
            # Stage 1c Layer 5: integrate the per-particle ρ_osm ODE after
            # the mechanical step (ρ_osm is driven by the post-step
            # tr(C_p) deformation rate; the K(ρ_osm) coupling enters
            # next step's `_p2g_momentum_and_stress`).
            self._integrate_osmotic_ode()
        if self.cfg.layer6_enabled:
            # Stage 2 Layer 6: integrate the Layer 6 ODE (mmp_total + ecm_strength).
            # Updates global scalars; substrate CSF impulse magnitude in
            # next step's `_grid_op_overdamped` reads ecm_strength_field.
            self._integrate_layer6_ode()

    @ti.kernel
    def _clear_grid(self):
        for I in ti.grouped(self.grid_m):
            self.grid_v[I] = ti.Vector.zero(ti.f32, 3)
            self.grid_m[I] = 0.0
            self.grid_kernel_weight[I] = 0.0
            self.grid_count[I] = 0
            self.grid_color[I] = 0.0
            self.grid_color_smooth[I] = 0.0
            self.grid_color_grad[I] = ti.Vector.zero(ti.f32, 3)
            self.grid_normal[I] = ti.Vector.zero(ti.f32, 3)
            self.grid_kappa[I] = 0.0
            self.grid_gamma[I] = 0.0
            self.grid_gamma_grad[I] = ti.Vector.zero(ti.f32, 3)
        self.diag_max_count[None] = 0
        # Stage 1a+ substrate reaction-impulse accumulator: reset per step.
        self.diag_substrate_impulse_z[None] = 0.0

    @ti.kernel
    def _tag_boundary(self):
        """Tag boundary particles by counting *empty* neighbour cells in the 3³ stencil.

        A particle is boundary iff ≥ `n_empty_threshold` of its 27 neighbour
        cells are empty (count == 0). Equivalent in spirit to the density-ratio
        threshold but avoids the global atomic_max bottleneck of the previous
        formulation. Threshold of 5/27 ≈ 0.18 is a robust heuristic in 3D.
        """
        n_empty_threshold = 5  # 5 of 27 neighbour cells empty ⇒ boundary

        # Pass 1: count particles per cell (small, local atomic_add — fine).
        for p in self.x:
            base = ti.cast(self.x[p] / self.cfg.dx_star, ti.i32)
            if (
                0 <= base[0] < self.cfg.grid_n
                and 0 <= base[1] < self.cfg.grid_n
                and 0 <= base[2] < self.cfg.grid_n
            ):
                ti.atomic_add(self.grid_count[base], 1)

        # Pass 2: tag particles by empty-neighbour count.
        for p in self.x:
            base = ti.cast(self.x[p] / self.cfg.dx_star, ti.i32)
            empty = 0
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                idx = base + ti.Vector([i - 1, j - 1, k - 1])
                if (
                    idx[0] < 0 or idx[0] >= self.cfg.grid_n
                    or idx[1] < 0 or idx[1] >= self.cfg.grid_n
                    or idx[2] < 0 or idx[2] >= self.cfg.grid_n
                ) or self.grid_count[idx] == 0:
                    empty += 1
            self.is_boundary[p] = 1 if empty >= n_empty_threshold else 0

    @ti.kernel
    def _seed_color_from_mass(self):
        """Reproducing-kernel-normalised colour (Adami, Hu & Adams 2010 §3).

        c[I] = grid_m[I] / (ρ_bulk · max(grid_kernel_weight[I], ε_div))

        The kernel-volume sum `grid_kernel_weight[I] = Σ_p w_pI · V₀` and the
        kernel-mass sum `grid_m[I] = Σ_p w_pI · m_p` share the SAME particle
        weights, so any random-pack packing-density variation cancels in the
        ratio: `m_p / V₀ = ρ_bulk` is a constant. By construction:
            c = 1   in any cell whose kernel support contains particles,
            c = 0   in pure-vacuum cells (grid_m = 0; ε_div only prevents 0/0).

        This eliminates the smooth radial trend in c that v11's diagnostic
        (`scripts/diag_csf_penetration.py`) showed as the source of CSF
        interior penetration: bulk |∇c| ≈ 60% of surface peak |∇c| dropped
        to a thin diffuse interface localised at the actual free surface.

        ε_div is a numerical-safety floor preventing 0/0 in vacuum cells —
        not a physical scale or a fitting parameter (passes Magic-Number
        Block tests 1, 2, and 3).
        """
        rho_bulk = self.cfg.density_star
        eps_div = 1.0e-6
        for I in ti.grouped(self.grid_color):
            denom = rho_bulk * ti.max(self.grid_kernel_weight[I], eps_div)
            self.grid_color[I] = self.grid_m[I] / denom

    @ti.kernel
    def _smooth_color_pass(self):
        """One pass of 3³ box averaging from `grid_color` into `grid_color_smooth`,
        then copy back. Brackbill (1992) §V recommends a few passes to widen
        the diffuse interface to ~5–7 dx, so the FD curvature operator
        `κ = -∇·n̂` resolves the physical surface curvature 2/R rather than
        the grid-scale 1/dx that a single-cell interface produces.

        Box averaging is the simplest such filter; it is rotationally symmetric
        only at low order, but for the curvature-magnitude check at the
        spheroid surface it is adequate and matches the prescription of the
        original Brackbill paper. Each pass corresponds to one application of
        a discrete diffusion step (forward Euler) of the colour field.
        """
        n_g = self.cfg.grid_n
        for I in ti.grouped(self.grid_color_smooth):
            i, j, k = I[0], I[1], I[2]
            s = 0.0
            count = 0
            for di, dj, dk in ti.static(ti.ndrange((-1, 2), (-1, 2), (-1, 2))):
                ii = i + di
                jj = j + dj
                kk = k + dk
                if 0 <= ii < n_g and 0 <= jj < n_g and 0 <= kk < n_g:
                    s += self.grid_color[ii, jj, kk]
                    count += 1
            self.grid_color_smooth[I] = s / ti.cast(count, ti.f32)
        # Copy smoothed back into grid_color so a subsequent pass diffuses further.
        for I in ti.grouped(self.grid_color):
            self.grid_color[I] = self.grid_color_smooth[I]

    @ti.kernel
    def _color_gradient_from_smoothed(self):
        """Central-difference ∇c from the smoothed colour field."""
        dx = self.cfg.dx_star
        inv_2dx = 1.0 / (2.0 * dx)
        for I in ti.grouped(self.grid_color_grad):
            i, j, k = I[0], I[1], I[2]
            gx = 0.0
            gy = 0.0
            gz = 0.0
            if 0 < i < self.cfg.grid_n - 1:
                gx = (self.grid_color[i + 1, j, k] - self.grid_color[i - 1, j, k]) * inv_2dx
            if 0 < j < self.cfg.grid_n - 1:
                gy = (self.grid_color[i, j + 1, k] - self.grid_color[i, j - 1, k]) * inv_2dx
            if 0 < k < self.cfg.grid_n - 1:
                gz = (self.grid_color[i, j, k + 1] - self.grid_color[i, j, k - 1]) * inv_2dx
            self.grid_color_grad[I] = ti.Vector([gx, gy, gz])

    def _build_csf_field(self) -> None:
        """Build the smoothed CSF colour field and its gradient.

        Pipeline:
          1. seed `grid_color = ρ/ρ_bulk` from `grid_m`
          2. apply N smoothing passes (Brackbill 1992 §V; default 2)
          3. compute ∇c from the smoothed field

        Number of smoothing passes is fixed at 2 in Stage 1a baseline. Higher
        N would lose surface localisation; lower N (= 0) reproduces the v8
        failure where κ encodes grid-scale variation.
        """
        n_smoothing_passes = 2
        self._seed_color_from_mass()
        for _ in range(n_smoothing_passes):
            self._smooth_color_pass()
        self._color_gradient_from_smoothed()

    @ti.kernel
    def _build_curvature(self):
        """Brackbill (1992) curvature operator — restored after v13 revert.

        Two-pass finite-difference scheme:
            (1) smoothed unit normal:  n̂[I] = ∇c[I] / sqrt(|∇c[I]|² + ε²)
            (2) mean curvature:        κ[I] = -∇·n̂[I]   (central differences)

        Sign convention: with c≈1 inside, c≈0 outside, ∇c points INWARD, so n̂
        points inward. Then -∇·n̂ is **positive for a convex droplet**. For a
        sphere of radius R: κ = 2/R. The static-sphere validation in the
        runner asserts this. ε² = 1e-6 keeps the f32 division stable in the
        bulk where |∇c| ≈ 0.

        The ε² regulariser inside the sqrt has a non-obvious robustness
        consequence: it bounds the magnitude of n̂ at every grid cell, so
        the subsequent central-FD divergence is shielded from the off-peak
        f″/f′ asymmetry of the smoothed colour profile that broke the v13
        Laplacian-form attempt (commit 3ba63c5). See v13 commit message
        and `docs/stage1a_aha_div_sanity.md` for the failure analysis.

        References:
        - Brackbill, Kothe, Zemach (1992) "A continuum method for modeling
          surface tension", J. Comp. Phys. 100, 335.
        - Adami, Hu, Adams (2010) "A new surface-tension formulation for
          multi-phase SPH using a reproducing divergence approximation",
          J. Comp. Phys. 229, 5011 (informs §3 colour normalisation, used;
          §4 reproducing divergence shown not to translate to uniform grid).
        """
        eps2 = 1.0e-6
        n_g = self.cfg.grid_n
        # Pass 1: unit normals from ∇c.
        for I in ti.grouped(self.grid_normal):
            g = self.grid_color_grad[I]
            mag = ti.sqrt(g.dot(g) + eps2)
            self.grid_normal[I] = g / mag

        # Pass 2: κ = -∇·n̂ via central difference.
        inv_2dx = 1.0 / (2.0 * self.cfg.dx_star)
        for I in ti.grouped(self.grid_kappa):
            i, j, k = I[0], I[1], I[2]
            div = 0.0
            if 0 < i < n_g - 1:
                div += (self.grid_normal[i + 1, j, k][0] - self.grid_normal[i - 1, j, k][0]) * inv_2dx
            if 0 < j < n_g - 1:
                div += (self.grid_normal[i, j + 1, k][1] - self.grid_normal[i, j - 1, k][1]) * inv_2dx
            if 0 < k < n_g - 1:
                div += (self.grid_normal[i, j, k + 1][2] - self.grid_normal[i, j, k - 1][2]) * inv_2dx
            self.grid_kappa[I] = -div

    @ti.kernel
    def _p2g_mass(self):
        """First P2G pass: scatter particle mass + kernel volume weight only.

        v15 (k.3) splits the legacy `_p2g` into two passes so that the
        density-based volumetric stress σ_vol = K·(ρ_ref/ρ_kernel − 1)·I can
        be evaluated using the post-scatter `grid_m` *before* momentum/stress
        is scattered. Pass 1 (this kernel) writes `grid_m` and the
        Adami-Hu-Adams 2010 §3 reproducing-kernel volume sum
        `grid_kernel_weight`. `_clear_grid` runs first, so we accumulate into
        already-zeroed fields.
        """
        m_p = self.cfg.particle_mass_star
        V0 = self.cfg.particle_volume_star
        dx = self.cfg.dx_star

        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    ti.atomic_add(self.grid_m[idx], weight * m_p)
                    # Reproducing-kernel volume sum (Adami-Hu-Adams 2010 §3).
                    ti.atomic_add(self.grid_kernel_weight[idx], weight * V0)

    @ti.kernel
    def _interpolate_rho_runtime(self):
        """G2P interpolation of grid kernel density into `_rho_kernel_p`.

            ρ_kernel_p[p] = max( Σ_I w_pI · grid_m[I] / dx³ ,  ρ_floor )

        Reads the post-`_p2g_mass` `grid_m` and writes the per-particle
        kernel-interpolated density used as the source of the v15 (k.3)
        density-based volumetric stress

            σ_vol_p = K · (ρ_ref_kernel / ρ_kernel_p − 1) · I.

        ρ_floor (= 0.1 · ρ_ref_kernel by `calibrate_reference_state`) is a
        numerical-safety clamp on the 1/ρ_kernel asymptote at ρ_kernel → 0.
        The exact value 0.1 has **no specific literature reference**: a
        search of Becker-Teschner (2007) (WCSPH origin, Tait EOS, no floor
        used), Adami-Hu-Adams (2010) §3 (Shepard normalisation, free-surface
        truncation bound ~0.5·ρ_ref), and Liu-Liu (2010) (review,
        isolated-particle threshold ~0.5·ρ_ref) found no exact 0.1·ρ_ref
        value. The choice is conservative below those bounds (5× below the
        free-surface kernel-truncation limit) so that the floor is *only*
        active for particles in the rarefied vacuum tail and never at the
        equilibrium ρ_eq ≈ 1.03·ρ_ref that the radius-drift gate measures.
        Magic-Number Block (`docs/12_validation.md`): derivable from the
        AHA-2010 truncation bound, grid-invariant (a fraction of a
        per-pack-calibrated quantity), not chosen to fit any gate value.
        See `docs/outcomes_v15.md` for the full literature-verification
        record.

        Structurally identical to `_interpolate_density_to_particles`
        (calibration scratch, validated through v11/v12); separated into a
        dedicated step-time kernel because writing into `_rho_kernel_p`
        (per-particle, runtime-hot) is conceptually distinct from writing
        into `_calib_rho` (per-particle, calibration-only).
        """
        dx = self.cfg.dx_star
        inv_vol = 1.0 / (dx ** 3)
        rho_floor = self._rho_floor[None]
        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]
            rho = 0.0
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    rho += weight * self.grid_m[idx] * inv_vol
            self._rho_kernel_p[p] = ti.max(rho, rho_floor)

    @ti.kernel
    def _p2g_momentum_and_stress(self):
        """Second P2G pass: scatter momentum + v15 density-based stress.

        Volumetric stress (v15 (k.3), density-based, hybrid Eulerian-Lagrangian):

            σ_vol_p = K · (ρ_ref_kernel / ρ_kernel_p − 1) · I,

        where ρ_kernel_p is the kernel-interpolated density at p (computed in
        `_interpolate_rho_runtime`, with floor clamp). This replaces the
        previous (v12) deformation-based form K·(det F − 1)·I that produced
        no bulk pressure response in the overdamped/no-flow interior (where
        ∇v ≈ 0 ⇒ F → I), the root cause of the v12 R/R₀ ≈ 0.67 contraction.

        Deviatoric stress (unchanged from v12):  τ_dev_p, integrated by the
        closed-form Maxwell exponential integrator in `_g2p_and_constitutive`.

        Momentum scatter follows the standard MLS-MPM affine form (Hu et al.
        2018). Mass and kernel-weight scatter happen in `_p2g_mass`, not here.

        Reference: Becker & Teschner (2007) for the WCSPH density-driven
        volumetric stress idea; Stomakhin et al. (2014) §3 for the
        elastic/plastic-split MPM precedent for splitting volumetric and
        deviatoric responses; Adami-Hu-Adams (2010) §3 (already used in v12)
        for the reproducing-kernel normalisation underlying the kernel-density
        estimator.
        """
        K = self.cfg.K_star
        V0 = self.cfg.particle_volume_star
        m_p = self.cfg.particle_mass_star
        dx = self.cfg.dx_star
        dt = self.cfg.dt_star
        rho_ref = self._rho_ref_kernel[None]
        # Stage 1a++ Layer 2 / Stage 1b Layer 3 active-stress coupling.
        # When Layer 3 is OFF: ζ(φ) = ζ_star (constant across all particles;
        # Stage 1a++ behaviour). When Layer 3 is ON: ζ(φ) is computed
        # per-particle from φ_p inside the loop.
        zeta_K = self.cfg.zeta_star * K
        zeta_min_K = self.cfg.zeta_min * K
        zeta_max_K = self.cfg.zeta_max * K

        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]

            rho_p = self._rho_kernel_p[p]
            # Stage 1c Layer 5 K(ρ_osm) coupling: K_eff per particle =
            # K_star · ρ_osm_p. When Layer 5 disabled, ρ_osm_p == 1.0 by
            # constructor + ODE inactivity; K_eff == K (Stage 1b
            # behaviour). Compile-time gating via Python attribute access:
            # if layer5 disabled, multiplication by 1.0 is a no-op.
            K_eff = K
            if ti.static(self.cfg.layer5_enabled):
                K_eff = K * self.rho_osm_p[p]
            stress_vol = K_eff * (rho_ref / rho_p - 1.0) * ti.Matrix.identity(ti.f32, 3)
            # Stage 1a++ Layer 2 / Stage 1b Layer 3 boundary-cell active
            # stress. Per-particle effective ζ_K depends on whether Layer 3
            # is enabled:
            #   Layer 3 OFF: ζ_K_eff = zeta_K (constant, Stage 1a++)
            #   Layer 3 ON : ζ_K_eff = ζ_min·K·(1−φ_p) + ζ_max·K·φ_p
            # Both branches collapse to the same expression when zeta_min ==
            # zeta_max == zeta_star, so the runtime conditional is at the
            # config level (`layer3_enabled`), evaluated at trace time.
            zeta_K_eff = zeta_K
            if ti.static(self.cfg.layer3_enabled):
                phi_p = self.phi_p[p]
                zeta_K_eff = zeta_min_K * (1.0 - phi_p) + zeta_max_K * phi_p
            stress_act = (
                -zeta_K_eff * ti.cast(self.is_boundary[p], ti.f32)
                * ti.Matrix.identity(ti.f32, 3)
            )
            stress = stress_vol + self.tau_dev[p] + stress_act

            stress_term = -dt * V0 * stress * (4.0 / (dx * dx))
            affine = stress_term + m_p * self.C[p]

            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    dpos = (ti.Vector([i, j, k]).cast(ti.f32) - fx) * dx
                    ti.atomic_add(self.grid_v[idx], weight * (m_p * self.v[p] + affine @ dpos))

    @ti.kernel
    def _grid_op_overdamped(self):
        """Convert momentum to velocity, apply CSF impulse, drag, and box safety net.

        CSF is applied here (not in P2G) so it sees the mass-divided velocity and
        scales by the local density. This was the source of the first-pilot NaN:
        applying `dt·γ·∇c` directly to the *momentum* grid before mass division
        produced runaway velocities at low-mass boundary cells.
        """
        dx = self.cfg.dx_star
        dt = self.cfg.dt_star
        gamma = self.cfg.gamma_star
        rho_ref = self.cfg.density_star
        # Only inject CSF where there is a reasonable amount of mass — otherwise
        # the inverse density is unphysical. Threshold = 10% of the reference
        # cell mass ρ·dx³.
        min_cell_mass = 0.1 * rho_ref * (dx ** 3)

        damp = 1.0 / (1.0 + self.cfg.drag_xi_star * dt)
        # Compile-time constants for the substrate / box-wall reflective BC.
        # When `substrate_enabled` is True the −z wall uses `n_contact_band`
        # (intentional substrate); otherwise it uses 3 (box-wall safety net).
        sub_band = self.cfg.n_contact_band if self.cfg.substrate_enabled else 3
        # Stage 1a+ Option β substrate CSF impulse coefficient. When > 0,
        # cells in the contact band receive an attractive impulse toward
        # z = 0 of magnitude γ_sub_eff · κ_sub_proxy · dt / ρ_local with
        # κ_sub_proxy = 1/dx* and n̂_sub = −ẑ.
        # Stage 2 Layer 6: γ_sub_eff = γ_sub_star · ecm_strength(t),
        # where ecm_strength(t) decays from 1.0 as MMP accumulates.
        # When Layer 6 disabled, ecm_strength stays at its initial value
        # (1.0 by default) → γ_sub_eff = γ_sub_star (Stage 1d behaviour).
        gamma_sub_star = self.cfg.gamma_sub_star * self.ecm_strength_field[None]
        inv_dx = 1.0 / dx
        # Path C effective gravity impulse coefficient. When > 0, every
        # grid cell with mass receives a per-step downward velocity
        # impulse `Δv_z = −gravity_star · dt` (no division by ρ_local
        # because gravity is force-per-unit-mass — the impulse is just
        # an acceleration · dt; in the f = m·a sense the per-cell force
        # is m · g_star, but in MPM we apply the velocity impulse
        # directly via overdamped balance ξ·v = m·g/m = g).
        gravity_star = self.cfg.gravity_star

        for I in ti.grouped(self.grid_m):
            m = self.grid_m[I]
            if m > 1e-30:
                v = self.grid_v[I] / m
                # CSF (Brackbill 1992): body force per unit volume is
                # F_v = γ · κ · ∇c, where κ = -∇·n̂ and n̂ = ∇c/|∇c|. The
                # impulse on grid velocity is dv = F_v · dt / ρ_local.
                # ∇c points INWARD at the free surface (c≈1 inside, c≈0
                # outside) and κ > 0 for a convex droplet, so dv pulls
                # the surface inward as physically required.
                # ∇c vanishes in the bulk (uniform colour) so this only acts
                # in a thin layer at the interface.
                if m > min_cell_mass:
                    rho_local = m / (dx ** 3)
                    v += (
                        dt * gamma * self.grid_kappa[I] * self.grid_color_grad[I]
                        / rho_local
                    )
                # Stage 1a+ Option β substrate CSF impulse: cells in the
                # contact band are pulled toward the substrate (n̂_sub = −ẑ,
                # κ_sub_proxy = 1/dx*). When γ_sub_star = 0 (Option α) the
                # impulse is identically zero by arithmetic, so the branch is
                # kept guard-less for Taichi-trace simplicity (no ti.static
                # on a kernel-local Expr). The min_cell_mass guard mirrors
                # the free-surface CSF block above so a low-mass cell does
                # not blow up via 1/ρ_local. Sign verified in
                # `docs/stage1a_plus_substrate_sanity.md` §6/§"Option β
                # addendum" check 5.
                i, j, k = I[0], I[1], I[2]
                if k < sub_band and m > min_cell_mass:
                    rho_local_sub = m / (dx ** 3)
                    v[2] -= dt * gamma_sub_star * inv_dx / rho_local_sub

                # Path C effective gravity impulse: every grid cell with
                # mass receives a per-step downward velocity impulse
                # `Δv_z = −gravity_star · dt`. When gravity_star = 0 (Path
                # C disabled) the impulse is identically zero by arithmetic.
                # Sign convention: gravity_star > 0 ⇒ dv_z < 0 ⇒ pulls
                # toward substrate at z = 0. See `docs/path_c_sanity.md`
                # check 5 for sign verification.
                v[2] -= dt * gravity_star

                # Stage 1d Layer 4 Marangoni impulse: at boundary cells
                # (m > min_cell_mass; same band as the existing CSF
                # impulse), apply tangential surface-tension-gradient
                # force. Tangent projection: ∇_s γ = (I − n̂⊗n̂) · ∇γ.
                # When Layer 4 disabled, grid_gamma_grad is zero
                # (uncleared in cleared-state) so the impulse is zero by
                # arithmetic. Sign per docs/stage1d_sanity.md check 5
                # (Pajic-Lijakovic 2022): velocity flows from low γ
                # toward high γ. ∇γ points toward higher γ; impulse
                # adds in +∇γ direction → flow toward high γ. ✓
                if ti.static(self.cfg.layer4_enabled):
                    if m > min_cell_mass:
                        rho_local_M = m / (dx ** 3)
                        n_hat = self.grid_normal[I]
                        grad_g = self.grid_gamma_grad[I]
                        # Tangent projection: g_t = g − (g·n̂)·n̂
                        g_dot_n = grad_g.dot(n_hat)
                        grad_g_tan = grad_g - g_dot_n * n_hat
                        v += dt * grad_g_tan / rho_local_M

                # Overdamped per-step damping (mild; physics enters via slow-time interpretation).
                v *= damp

                # Reflective box-wall safety net (3-cell margin).
                if i < 3 and v[0] < 0:
                    v[0] = 0.0
                if i > self.cfg.grid_n - 3 and v[0] > 0:
                    v[0] = 0.0
                if j < 3 and v[1] < 0:
                    v[1] = 0.0
                if j > self.cfg.grid_n - 3 and v[1] > 0:
                    v[1] = 0.0
                # −z wall: under Stage 1a+ Option α this is the *intentional*
                # rigid substrate (γ_sub = 0, mechanical anchor only); when
                # `substrate_enabled` is False it is the same 3-cell box-wall
                # safety net as the +x, −x, +y, −y, +z faces. Either way the
                # clamp logic is identical; the only Stage-1a+-specific
                # behaviour is recording the per-step reaction impulse for
                # the anchor-force-balance gate.
                if k < sub_band and v[2] < 0:
                    if ti.static(self.cfg.substrate_enabled):
                        ti.atomic_add(
                            self.diag_substrate_impulse_z[None],
                            ti.cast(m * (-v[2]), ti.f64),
                        )
                    v[2] = 0.0
                if k > self.cfg.grid_n - 3 and v[2] > 0:
                    v[2] = 0.0

                self.grid_v[I] = v

    @ti.kernel
    def _g2p_and_constitutive(self):
        dt = self.cfg.dt_star
        dx = self.cfg.dx_star
        mu = self.cfg.mu_star
        tau = self.cfg.tau_star
        decay = ti.exp(-dt / tau)
        drive = 2.0 * mu * tau * (1.0 - decay)

        for p in self.x:
            base = ti.cast(self.x[p] / dx - 0.5, ti.i32)
            fx = self.x[p] / dx - ti.cast(base, ti.f32)
            w = [
                0.5 * (1.5 - fx) ** 2,
                0.75 - (fx - 1.0) ** 2,
                0.5 * (fx - 0.5) ** 2,
            ]

            new_v = ti.Vector.zero(ti.f32, 3)
            new_C = ti.Matrix.zero(ti.f32, 3, 3)
            for i, j, k in ti.static(ti.ndrange(3, 3, 3)):
                weight = w[i][0] * w[j][1] * w[k][2]
                idx = base + ti.Vector([i, j, k])
                if (
                    0 <= idx[0] < self.cfg.grid_n
                    and 0 <= idx[1] < self.cfg.grid_n
                    and 0 <= idx[2] < self.cfg.grid_n
                ):
                    dpos = ti.cast(ti.Vector([i, j, k]), ti.f32) - fx
                    g_v = self.grid_v[idx]
                    new_v += weight * g_v
                    new_C += 4.0 / dx * weight * g_v.outer_product(dpos)

            self.v[p] = new_v
            self.C[p] = new_C
            self.x[p] += dt * new_v
            self.F[p] = (ti.Matrix.identity(ti.f32, 3) + dt * new_C) @ self.F[p]

            sym_C = 0.5 * (new_C + new_C.transpose())
            tr = sym_C[0, 0] + sym_C[1, 1] + sym_C[2, 2]
            eps_dot_dev = sym_C - (tr / 3.0) * ti.Matrix.identity(ti.f32, 3)
            self.tau_dev[p] = decay * self.tau_dev[p] + drive * eps_dot_dev

    # ------------------------------------------------------ diagnostics ----
    @ti.kernel
    def _compute_invariants(self):
        self.diag_mass[None] = 0.0
        self.diag_momentum[None] = ti.Vector.zero(ti.f64, 3)
        self.diag_kinetic_energy[None] = 0.0
        self.diag_strain_energy[None] = 0.0
        self.diag_surface_energy[None] = 0.0
        self.diag_max_speed[None] = 0.0
        self.diag_nan_count[None] = 0
        # Stage 1a++ Layer 2 diagnostics, reset per call.
        self.diag_active_power[None] = 0.0
        self.diag_n_boundary[None] = 0
        # Stage 1b Layer 3 diagnostics, reset per call.
        self.diag_phi_sum[None] = 0.0
        self.diag_phi_sum_sq[None] = 0.0
        self.diag_phi_min[None] = 1.0e10
        self.diag_phi_max[None] = -1.0e10
        self.diag_phi_boundary_sum[None] = 0.0
        # Stage 1b.b Layer 3 split diagnostics, reset per call.
        self.diag_phi_memory_sum[None] = 0.0
        self.diag_phi_memory_min[None] = 1.0e10
        self.diag_phi_memory_max[None] = -1.0e10
        self.diag_c_act_sum[None] = 0.0
        self.diag_c_act_min[None] = 1.0e10
        self.diag_c_act_max[None] = -1.0e10
        self.diag_c_act_band_sum[None] = 0.0
        self.diag_c_act_band_count[None] = 0
        # Path C diagnostics, reset per call.
        self.diag_grav_pe[None] = 0.0
        self.diag_com_z_sum[None] = 0.0
        # Stage 1c Layer 5 ρ_osm diagnostics, reset per call.
        self.diag_rho_osm_sum[None] = 0.0
        self.diag_rho_osm_sum_sq[None] = 0.0
        self.diag_rho_osm_min[None] = 1.0e10
        self.diag_rho_osm_max[None] = -1.0e10

        m_p = ti.cast(self.cfg.particle_mass_star, ti.f64)
        K = ti.cast(self.cfg.K_star, ti.f64)
        mu = ti.cast(self.cfg.mu_star, ti.f64)
        V0 = ti.cast(self.cfg.particle_volume_star, ti.f64)
        gamma = ti.cast(self.cfg.gamma_star, ti.f64)
        gamma_sub = ti.cast(self.cfg.gamma_sub_star, ti.f64)
        h_band = self.cfg.n_contact_band * self.cfg.dx_star

        for p in self.x:
            v = self.v[p]
            speed = v.norm()
            ti.atomic_max(self.diag_max_speed[None], speed)

            bad = 0
            for d in ti.static(range(3)):
                # NaN-check via the IEEE-754 self-inequality identity.
                if v[d] != v[d] or ti.abs(v[d]) > 1e30:
                    bad = 1
                if self.x[p][d] != self.x[p][d] or ti.abs(self.x[p][d]) > 1e30:
                    bad = 1
            if bad == 1:
                ti.atomic_add(self.diag_nan_count[None], 1)

            ti.atomic_add(self.diag_mass[None], m_p)
            for d in ti.static(range(3)):
                ti.atomic_add(self.diag_momentum[None][d], m_p * ti.cast(v[d], ti.f64))

            ti.atomic_add(self.diag_kinetic_energy[None], 0.5 * m_p * ti.cast(speed * speed, ti.f64))

            # Path C gravitational potential energy: U_grav = ρ_p · g_star
            # · z_p · V₀ (positive above z=0; decreases as spheroid sinks).
            # Sign chosen so dU_grav/dz > 0 (need to do work to lift).
            # When gravity_star = 0 (Path C off) the contribution is 0 by
            # arithmetic.
            z_p = self.x[p][2]
            ti.atomic_add(
                self.diag_grav_pe[None],
                ti.cast(self.cfg.density_star * self.cfg.gravity_star * z_p, ti.f64) * V0,
            )
            ti.atomic_add(self.diag_com_z_sum[None], ti.cast(z_p, ti.f64))

            # v15 (k.3): volumetric strain energy uses the density-based form
            #   U_vol = (1/2) K_eff (ρ_ref/ρ_kernel − 1)²
            # to match the v15 constitutive law σ_vol = K_eff (ρ_ref/ρ_kernel − 1) I.
            # K_eff under Stage 1c Layer 5 is K · ρ_osm_p; under Layer 5 OFF
            # ρ_osm_p == 1.0 so K_eff = K (Stage 1b behaviour).
            # F is still updated for diagnostics + Maxwell deviatoric coupling,
            # but no longer drives the volumetric channel. The deviatoric
            # contribution τ²/(4μ) is unchanged from v12.
            rho_p_d = ti.cast(self._rho_kernel_p[p], ti.f64)
            rho_ref_d = ti.cast(self._rho_ref_kernel[None], ti.f64)
            vol_strain = rho_ref_d / rho_p_d - 1.0
            K_eff_d = K
            if ti.static(self.cfg.layer5_enabled):
                K_eff_d = K * ti.cast(self.rho_osm_p[p], ti.f64)
            tau = self.tau_dev[p]
            tau_norm_sq = 0.0
            for ii, jj in ti.static(ti.ndrange(3, 3)):
                tau_norm_sq += tau[ii, jj] * tau[ii, jj]
            U = 0.5 * K_eff_d * vol_strain * vol_strain + ti.cast(tau_norm_sq, ti.f64) / (4.0 * (mu + 1e-12))
            ti.atomic_add(self.diag_strain_energy[None], U * V0)

            if self.is_boundary[p] == 1:
                # Surface energy proxy: γ × per-particle surface element ≈ γ × V₀^(2/3).
                ti.atomic_add(self.diag_surface_energy[None], gamma * ti.cast(V0 ** (2.0 / 3.0), ti.f64))
                ti.atomic_add(self.diag_n_boundary[None], 1)
                # Stage 1a++ Layer 2 active power per boundary particle:
                #   P_act,p = -ζ_eff·K · tr(C_p) · V₀
                # (σ_act_p = -ζ_eff·K·I, ε̇_p ≈ sym(C_p),
                #  σ:ε̇ = -ζ_eff·K · tr(C_p)).
                # ζ_eff is constant under Layer 2 only; under Layer 3 it
                # depends on φ_p (per-particle φ-modulated coupling).
                C_p = self.C[p]
                tr_C = C_p[0, 0] + C_p[1, 1] + C_p[2, 2]
                zeta_eff = self.cfg.zeta_star
                if ti.static(self.cfg.layer3_enabled):
                    phi_p = self.phi_p[p]
                    zeta_eff = (
                        self.cfg.zeta_min * (1.0 - phi_p)
                        + self.cfg.zeta_max * phi_p
                    )
                ti.atomic_add(
                    self.diag_active_power[None],
                    ti.cast(-zeta_eff * self.cfg.K_star * tr_C * V0, ti.f64),
                )
                # Stage 1b Layer 3 boundary-shell <φ> accumulator.
                if ti.static(self.cfg.layer3_enabled):
                    ti.atomic_add(
                        self.diag_phi_boundary_sum[None],
                        ti.cast(self.phi_p[p], ti.f64),
                    )

            # Stage 1b Layer 3 per-particle φ aggregation (over all
            # particles, not just boundary; bulk avg comes from <total>−
            # <boundary> via the runner).
            if ti.static(self.cfg.layer3_enabled):
                phi_p_all = self.phi_p[p]
                ti.atomic_add(
                    self.diag_phi_sum[None],
                    ti.cast(phi_p_all, ti.f64),
                )
                ti.atomic_add(
                    self.diag_phi_sum_sq[None],
                    ti.cast(phi_p_all * phi_p_all, ti.f64),
                )
                ti.atomic_min(self.diag_phi_min[None], phi_p_all)
                ti.atomic_max(self.diag_phi_max[None], phi_p_all)
                # Stage 1b.b: aggregate phi_memory + c_act per gates 5a/5b.
                if ti.static(self.cfg.layer3_split):
                    phi_mem_p = self.phi_memory_p[p]
                    c_act_p_all = self.c_act_p[p]
                    ti.atomic_add(
                        self.diag_phi_memory_sum[None],
                        ti.cast(phi_mem_p, ti.f64),
                    )
                    ti.atomic_min(self.diag_phi_memory_min[None], phi_mem_p)
                    ti.atomic_max(self.diag_phi_memory_max[None], phi_mem_p)
                    ti.atomic_add(
                        self.diag_c_act_sum[None],
                        ti.cast(c_act_p_all, ti.f64),
                    )
                    ti.atomic_min(self.diag_c_act_min[None], c_act_p_all)
                    ti.atomic_max(self.diag_c_act_max[None], c_act_p_all)
                    # Boundary band only (z_p < h_band) for gate 5b.
                    if self.x[p][2] < h_band:
                        ti.atomic_add(
                            self.diag_c_act_band_sum[None],
                            ti.cast(c_act_p_all, ti.f64),
                        )
                        ti.atomic_add(self.diag_c_act_band_count[None], 1)

            # Stage 1c Layer 5 per-particle ρ_osm aggregation.
            if ti.static(self.cfg.layer5_enabled):
                rho_osm_p_all = self.rho_osm_p[p]
                ti.atomic_add(
                    self.diag_rho_osm_sum[None],
                    ti.cast(rho_osm_p_all, ti.f64),
                )
                ti.atomic_add(
                    self.diag_rho_osm_sum_sq[None],
                    ti.cast(rho_osm_p_all * rho_osm_p_all, ti.f64),
                )
                ti.atomic_min(self.diag_rho_osm_min[None], rho_osm_p_all)
                ti.atomic_max(self.diag_rho_osm_max[None], rho_osm_p_all)

            # Stage 1a+ Option β substrate-adhesion energy: per particle in
            # the contact band, subtract γ_sub · V₀^(2/3) (adhesion *reduces*
            # total energy). Gated on `substrate_enabled` (compile-time
            # Python bool, ti.static-safe); under Option α (γ_sub = 0) the
            # contribution is arithmetically zero, so no inner γ_sub-guard
            # is needed. The diagnostic is informational; the existing
            # energy-monotone gate in the runner uses KE + U_strain only and
            # is unaffected by this addition.
            if ti.static(self.cfg.substrate_enabled):
                if self.x[p][2] < h_band:
                    ti.atomic_add(
                        self.diag_surface_energy[None],
                        -gamma_sub * ti.cast(V0 ** (2.0 / 3.0), ti.f64),
                    )

    def substrate_diagnostics(self) -> dict:
        """Stage 1a+ Option α (γ_sub = 0) substrate gates + diagnostics.

        Computes:
        - n_contact_band_particles: # particles in the contact band z* < n·dx*
        - rho_kernel_contact_over_ref: <ρ_kernel>_band / ρ_ref
            Gate window broadened from [0.85, 1.15] to [0.65, 1.15] in
            Option F Week 2 (per docs/anchor_force_balance_investigation.md):
            near the −z reflective substrate boundary the SPH/MPM kernel
            truncation (Adami-Hu-Adams 2010 §3) biases the kernel-density
            estimate low by ≈ 30%, so [0.65, 1.15] is the principled
            truncation-aware envelope.
        - F_substrate_per_step: substrate reaction force from accumulated
            impulse / dt (over the most recent step).
        - F_pressure_compressive: compressive part of the contact-band
            volumetric stress integral, ≈ Σ_{P>0} P · V₀/h_band. This
            excludes the tensile contribution caused by kernel truncation
            (P_per_p < 0 when ρ_kernel < ρ_ref), reported separately as
            F_pressure_tensile_artefact for diagnostic transparency.
        - F_pressure_down: legacy (signed) integrand, kept for backwards
            comparison only; this is the value reported in pre-Week-2
            gate reports.
        - F_gravity_band: gravity contribution acting on contact-band
            particles, ρ · V₀ · g_star · n_contact (Path C term).
        - F_total_required = F_pressure_compressive + F_gravity_band.
        - anchor_force_balance_rel_err: |F_substrate − F_total_required|
            / max(|F_total_required|, ε); gate ≤ 0.20 (Option F Week 2
            new contract: substrate provides the compressive bulk +
            gravity load, ignoring the tensile truncation artefact).
        - contact_area_xy_hull: convex-hull area of contact-band particle
            (x, y) projections, in dimensionless area units.
        - apparent_contact_angle_deg: geometric angle of the spheroid
            surface near the substrate; log-only diagnostic under γ_sub = 0
            (no Young analytical reference).

        Returns {"valid": False, ...} if substrate is disabled or there are
        no contact-band particles.
        """
        if not self.cfg.substrate_enabled:
            return {"valid": False, "reason": "substrate not enabled"}

        x_np = self.x.to_numpy()
        rho_p = self._rho_kernel_p.to_numpy()
        z = x_np[:, 2]
        h_band = self.cfg.n_contact_band * self.cfg.dx_star
        in_band = z < h_band
        n_contact = int(in_band.sum())
        if n_contact == 0:
            return {"valid": False, "reason": "no contact-band particles"}

        rho_ref = float(self._rho_ref_kernel[None])
        rho_band = rho_p[in_band].astype(np.float64)
        rho_band_mean = float(rho_band.mean())

        K = float(self.cfg.K_star)
        V0 = float(self.cfg.particle_volume_star)
        # Phase 1.1 anchor-force-balance integrand update (Stage 1c sanity-md
        # contract change): K → K_eff per particle. Under Layer 5 active,
        # K_eff_p = K · rho_osm_p (Stage 1c K(ρ_osm) coupling). When Layer 5
        # is off, rho_osm_p == 1.0 and K_eff_p == K (Stage 1a+ behaviour).
        rho_osm_p = self.rho_osm_p.to_numpy()
        rho_osm_band = rho_osm_p[in_band].astype(np.float64)
        K_eff_band = K * rho_osm_band  # per-particle effective bulk modulus
        # Hydrostatic pressure per particle. σ_vol = K_eff(ρ_ref/ρ − 1)·I;
        # P_per_p = -tr(σ_vol)/3 = K_eff(1 − ρ_ref/ρ).
        # Note: in the kernel-truncation regime near the substrate boundary
        # (Adami-Hu-Adams 2010 §3), ρ_kernel < ρ_ref so P_per_p < 0
        # (tensile). The down-pushing force on the substrate is the
        # *compressive* part only — see anchor-force-balance gate below.
        P_per_p = K_eff_band * (1.0 - rho_ref / np.clip(rho_band, 1e-30, None))
        # Option F Week 2 contract change (per
        # docs/anchor_force_balance_investigation.md, PI-authorised
        # 2026-04-29). The legacy comparison failed systematically because
        # (1) the F_pressure_down formula used `(1 − ρ_ref/ρ_kernel)` which
        # goes negative in the kernel-truncation regime (ρ_kernel < ρ_ref
        # near the −z reflective BC, Adami-Hu-Adams 2010 §3) → the bulk-
        # side estimator predicted the wrong sign of normal force; and
        # (2) F_gravity used only contact-band particles whereas the
        # substrate must hold up the *entire* spheroid (gravity acts on
        # all particles). The F_substrate_up impulse is unaffected.
        #
        # New balance: F_substrate_up ≈ M_total · g_star + F_compressive,
        # where F_compressive is the compressive (positive-P) contribution
        # of the volumetric stress integrand. Tensile contributions
        # (negative-P due to kernel truncation) are diagnostic-only.
        rho_per_p = float(self.cfg.density_star)  # uniform density_star
        n_total = self.cfg.n_particles
        F_gravity_total = rho_per_p * V0 * float(self.cfg.gravity_star) * float(n_total)
        # Compressive part only (mask off tensile contributions from kernel
        # truncation). Tensile contributions are reported separately for
        # diagnostic transparency.
        compressive_mask = P_per_p > 0.0
        F_compressive = float((P_per_p[compressive_mask] * V0 / h_band).sum())
        F_tensile_artefact = float(((-P_per_p)[~compressive_mask] * V0 / h_band).sum())
        # Legacy band-only F_gravity (kept for diagnostic continuity).
        F_gravity_band = rho_per_p * V0 * float(self.cfg.gravity_star) * float(in_band.sum())
        # Legacy F_pressure_down (kept for backwards comparison only):
        # full signed integrand + band-only gravity.
        F_pressure_down_legacy = float((P_per_p * V0 / h_band).sum()) + F_gravity_band
        F_total_required = F_compressive + F_gravity_total
        F_substrate_up = float(self.diag_substrate_impulse_z[None]) / self.cfg.dt_star

        denom = max(abs(F_total_required), 1e-12)
        balance_err = abs(F_substrate_up - F_total_required) / denom

        # Contact area in (x, y) plane via convex hull (host-side, scipy).
        from scipy.spatial import ConvexHull
        xy_band = x_np[in_band, :2].astype(np.float64)
        A_contact_xy = float("nan")
        if n_contact >= 3:
            try:
                hull = ConvexHull(xy_band)
                A_contact_xy = float(hull.volume)  # 2D ConvexHull.volume = area
            except Exception:
                pass

        # Apparent contact angle (geometric, log-only): linear fit of the
        # spheroid radial extent r(z) over the lowest 0.2·R₀ slab gives
        # dr/dz, and the apparent angle is atan(-dr/dz) (radians) measured
        # from the substrate plane (90° = straight wall, 0° = parallel film).
        # This is a *geometric* quantity (not a curvature), so the v13
        # f″/f′ pathology does not apply.
        from numpy.polynomial import polynomial as Pnumpy
        cx = x_np[:, 0].mean()
        cy = x_np[:, 1].mean()
        r_all = np.sqrt((x_np[:, 0] - cx) ** 2 + (x_np[:, 1] - cy) ** 2)
        z_lo, z_hi = 0.0, 0.2 * self.cfg.radius_star
        slab = (z >= z_lo) & (z <= z_hi)
        theta_deg = float("nan")
        if int(slab.sum()) >= 8:
            # Per-z radial outer extent: take the 90th percentile of r within
            # narrow z bins to get the effective surface r(z).
            n_bins = 5
            z_edges = np.linspace(z_lo, z_hi, n_bins + 1)
            zc, rc = [], []
            for b in range(n_bins):
                m = slab & (z >= z_edges[b]) & (z < z_edges[b + 1])
                if int(m.sum()) >= 3:
                    zc.append(0.5 * (z_edges[b] + z_edges[b + 1]))
                    rc.append(float(np.percentile(r_all[m], 90)))
            if len(zc) >= 3:
                coefs = np.polyfit(zc, rc, 1)  # rc = coefs[0]*zc + coefs[1]
                drdz = float(coefs[0])
                theta_deg = float(np.degrees(np.arctan2(1.0, max(-drdz, -1e6))))
                # arctan2(1, -dr/dz): when -dr/dz → +∞ (vertical wall) → 90°,
                # when -dr/dz → 0 (flat film) → 90° too — degenerate. Use
                # a more direct formula:
                theta_deg = float(np.degrees(np.arctan(-drdz)) + 90.0)
                # = 90° at dr/dz = 0 (flat near-substrate surface, hemispherical),
                # > 90° at dr/dz < 0 (wider at z=0 than above, dewetting),
                # < 90° at dr/dz > 0 (narrowing toward z=0, beading).

        return {
            "valid": True,
            "n_contact_band_particles": n_contact,
            "h_band_star": h_band,
            "rho_kernel_contact_mean": rho_band_mean,
            "rho_kernel_contact_over_ref": rho_band_mean / rho_ref,
            "F_substrate_per_step": F_substrate_up,
            "F_pressure_down": F_pressure_down_legacy,
            "F_pressure_compressive": F_compressive,
            "F_pressure_tensile_artefact": F_tensile_artefact,
            "F_gravity_band": F_gravity_band,
            "F_gravity_total": F_gravity_total,
            "F_total_required": F_total_required,
            "anchor_force_balance_rel_err": balance_err,
            "contact_area_xy_hull": A_contact_xy,
            "apparent_contact_angle_deg": theta_deg,
        }

    def measure_surface_curvature(self) -> dict:
        """Static curvature validation: build a CSF + curvature pass on the
        current particle configuration and measure mean κ at the free surface.

        For a sphere of radius R*, the analytical surface curvature is
        κ = 2/R*. We localise the surface band by gathering cells that
        contain at least one **boundary-tagged particle** (the density-based
        free-surface tag), as opposed to a top-decile |∇c| selection. This
        avoids picking up "ghost" cells slightly outside the spheroid where
        |∇c| is small but κ from the FD operator can be noisy due to f32
        cancellation in `∇·n̂`.

        Returns a dict with the measured κ, the analytical 2/R, the relative
        error, and counts. The runner asserts |κ_measured / κ_analytical − 1|
        ≤ 0.10 (10% — generous for f32 finite differences on a 64³ grid).
        """
        # We don't move particles; just build the colour + curvature fields
        # using the current state.
        self._clear_grid()
        self._tag_boundary()
        # Mass-only scatter so colour is well-defined; momentum kernels are
        # not needed for κ measurement.
        self._scatter_mass_only()
        self._build_csf_field()
        self._build_curvature()
        ti.sync()

        kappa_np = self.grid_kappa.to_numpy()                 # (n_g, n_g, n_g)
        grad_np = self.grid_color_grad.to_numpy()             # (n_g, n_g, n_g, 3)
        grad_mag_3d = np.sqrt((grad_np ** 2).sum(axis=-1))    # (n_g, n_g, n_g)

        # Localise the surface via boundary-tagged particle cells.
        x_np = self.x.to_numpy()
        is_b = self.is_boundary.to_numpy().astype(bool)
        bdy_pos = x_np[is_b]
        if bdy_pos.shape[0] == 0:
            return {"valid": False, "reason": "no boundary-tagged particles"}

        cells = np.floor(bdy_pos / self.cfg.dx_star).astype(np.int32)
        n_g = self.cfg.grid_n
        in_range = (
            (cells[:, 0] >= 0) & (cells[:, 0] < n_g)
            & (cells[:, 1] >= 0) & (cells[:, 1] < n_g)
            & (cells[:, 2] >= 0) & (cells[:, 2] < n_g)
        )
        cells = cells[in_range]
        if cells.shape[0] == 0:
            return {"valid": False, "reason": "no in-range boundary cells"}

        # Deduplicate cells.
        flat = np.unique(cells[:, 0] * n_g * n_g + cells[:, 1] * n_g + cells[:, 2])
        ix = flat // (n_g * n_g)
        iy = (flat // n_g) % n_g
        iz = flat % n_g
        kappa_surface_all = kappa_np[ix, iy, iz]
        grad_at_band = grad_mag_3d[ix, iy, iz]

        # Within the boundary-cell set, keep only cells where |∇c| is at least
        # 25% of the maximum: this filters out boundary-tagged particles that
        # happen to sit in cells with near-bulk colour (a few isolated
        # interior particles tagged as boundary by the empty-neighbour
        # heuristic). 25% is a band selector, not a fitted parameter.
        if grad_at_band.max() <= 0:
            return {"valid": False, "reason": "no |∇c| > 0 in boundary band"}
        cutoff = 0.25 * float(grad_at_band.max())
        m = grad_at_band >= cutoff
        kappa_surface = kappa_surface_all[m]
        if kappa_surface.size == 0:
            return {"valid": False, "reason": "empty surface band after |∇c| filter"}

        # ---------- bulk / surface |∇c| ratio (CSF localisation gate) -----
        # Measure the spatial localisation of the colour-gradient field. With
        # Adami-Hu-Adams 2010 reproducing-kernel normalisation, c ≡ 1 in any
        # cell whose kernel sees particles, so |∇c| should be confined to the
        # surface band. The bulk shell (r/R₀ ∈ [0.2, 0.7]) provides a
        # representative "deep interior" sample where |∇c| should be
        # numerically negligible compared to the surface peak.
        n_g_local = self.cfg.grid_n
        ix_arr = np.arange(n_g_local)
        Ix_arr, Iy_arr, Iz_arr = np.meshgrid(ix_arr, ix_arr, ix_arr, indexing="ij")
        Xc = (Ix_arr + 0.5) * self.cfg.dx_star
        Yc = (Iy_arr + 0.5) * self.cfg.dx_star
        Zc = (Iz_arr + 0.5) * self.cfg.dx_star
        centre_d = self.cfg.domain_star * 0.5
        rr = np.sqrt((Xc - centre_d) ** 2 + (Yc - centre_d) ** 2 + (Zc - centre_d) ** 2)
        r_over_R0 = rr / self.cfg.radius_star
        bulk_mask = (r_over_R0 >= 0.2) & (r_over_R0 <= 0.7)
        n_bulk = int(bulk_mask.sum())
        bulk_grad_mean = float(grad_mag_3d[bulk_mask].mean()) if n_bulk > 0 else 0.0
        surface_grad_peak = float(grad_mag_3d.max())
        bulk_to_surface_grad_ratio = (
            bulk_grad_mean / surface_grad_peak
            if surface_grad_peak > 0
            else 0.0
        )

        return {
            "valid": True,
            "kappa_measured_mean": float(kappa_surface.mean()),
            "kappa_measured_median": float(np.median(kappa_surface)),
            "kappa_measured_std": float(kappa_surface.std()),
            "kappa_analytical": 2.0 / float(self.cfg.radius_star),
            "n_surface_cells": int(kappa_surface.size),
            "n_boundary_cells_total": int(flat.size),
            "selection": "boundary_particle_cells_with_grad_filter",
            "bulk_grad_c_mean": bulk_grad_mean,
            "surface_grad_c_peak": surface_grad_peak,
            "bulk_to_surface_grad_ratio": bulk_to_surface_grad_ratio,
            "n_bulk_cells_sampled": n_bulk,
        }

    def invariants(self) -> dict:
        self._compute_invariants()
        ti.sync()
        p = self.diag_momentum.to_numpy()
        return {
            "mass_star": float(self.diag_mass[None]),
            "momentum_star": np.array(p, dtype=float),
            "kinetic_energy_star": float(self.diag_kinetic_energy[None]),
            "strain_energy_star": float(self.diag_strain_energy[None]),
            "surface_energy_star": float(self.diag_surface_energy[None]),
            "max_speed_star": float(self.diag_max_speed[None]),
            "nan_count": int(self.diag_nan_count[None]),
            # Stage 1a++ Layer 2 diagnostics.
            "active_power_star": float(self.diag_active_power[None]),
            "n_boundary": int(self.diag_n_boundary[None]),
            # Stage 1b Layer 3 diagnostics. Under Layer 3 OFF these are
            # populated only with the boundary-tag count of n_boundary
            # particles having phi_initial; bulk/std are trivial. Under
            # Layer 3 ON they reflect the live φ field state.
            "phi_sum": float(self.diag_phi_sum[None]),
            "phi_sum_sq": float(self.diag_phi_sum_sq[None]),
            "phi_min": float(self.diag_phi_min[None]) if self.cfg.layer3_enabled else float("nan"),
            "phi_max": float(self.diag_phi_max[None]) if self.cfg.layer3_enabled else float("nan"),
            "phi_boundary_sum": float(self.diag_phi_boundary_sum[None]),
            # Stage 1b.b Layer 3 split diagnostics. Under layer3_split OFF
            # (legacy single-φ) these are sentinel NaN.
            "phi_memory_sum": (
                float(self.diag_phi_memory_sum[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "phi_memory_min": (
                float(self.diag_phi_memory_min[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "phi_memory_max": (
                float(self.diag_phi_memory_max[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "c_act_sum": (
                float(self.diag_c_act_sum[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "c_act_min": (
                float(self.diag_c_act_min[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "c_act_max": (
                float(self.diag_c_act_max[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "c_act_band_sum": (
                float(self.diag_c_act_band_sum[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else float("nan")
            ),
            "c_act_band_count": (
                int(self.diag_c_act_band_count[None])
                if (self.cfg.layer3_enabled and self.cfg.layer3_split) else 0
            ),
            # Path C diagnostics.
            "grav_pe_star": float(self.diag_grav_pe[None]),
            "com_z_star": (
                float(self.diag_com_z_sum[None]) / float(self.cfg.n_particles)
                if self.cfg.n_particles > 0 else float("nan")
            ),
            # Stage 1c Layer 5 ρ_osm diagnostics.
            "rho_osm_sum": float(self.diag_rho_osm_sum[None]),
            "rho_osm_sum_sq": float(self.diag_rho_osm_sum_sq[None]),
            "rho_osm_min": float(self.diag_rho_osm_min[None]) if self.cfg.layer5_enabled else float("nan"),
            "rho_osm_max": float(self.diag_rho_osm_max[None]) if self.cfg.layer5_enabled else float("nan"),
            # Stage 2 Layer 6 diagnostics.
            "mmp_total": float(self.mmp_total_field[None]),
            "ecm_strength": float(self.ecm_strength_field[None]),
        }
