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
        # Stage 1b Layer 3 diagnostics: bulk-shell <φ>, boundary-shell <φ>,
        # contact-band <φ>, plus min/max for the φ ∈ [0,1] invariant gate.
        self.diag_phi_sum = ti.field(dtype=ti.f64, shape=())
        self.diag_phi_sum_sq = ti.field(dtype=ti.f64, shape=())
        self.diag_phi_min = ti.field(dtype=ti.f32, shape=())
        self.diag_phi_max = ti.field(dtype=ti.f32, shape=())
        self.diag_phi_boundary_sum = ti.field(dtype=ti.f64, shape=())

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
        """Stage 1b Layer 3 forward-Euler φ-ODE integration.

        Per `docs/03_adhesion_dynamics.md` and `docs/stage1b_layer3_sanity.md`:
        with substrate-induced signalling intensity simplified to S ≡ 1
        (Stage 1b uses single-condition Col1 substrate, always present from
        t=0 per Stage 1a+ initial-position scheme), the φ-ODE reduces to

            dφ/dt = k_+ · (1 − φ) − k_- · φ

        Forward-Euler update with stiffness checked at constructor time
        (`dt · (k_+ + k_-) ≤ 0.5`). The `clip(0, 1)` defends against any
        per-step overshoot at the boundaries of the [0, 1] domain.
        """
        dt = self.cfg.dt_star
        k_plus = self.cfg.k_plus_star
        k_minus = self.cfg.k_minus_star
        for p in self.phi_p:
            phi = self.phi_p[p]
            dphi = dt * (k_plus * (1.0 - phi) - k_minus * phi)
            new_phi = phi + dphi
            # Defense-in-depth clamp on [0, 1].
            if new_phi < 0.0:
                new_phi = 0.0
            if new_phi > 1.0:
                new_phi = 1.0
            self.phi_p[p] = new_phi

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
        self._p2g_momentum_and_stress()     # scatters momentum + density-based stress
        self._grid_op_overdamped()          # converts to velocity, applies CSF impulse γ·κ·∇c
        self._g2p_and_constitutive()
        if self.cfg.layer3_enabled:
            # Stage 1b Layer 3: integrate the per-particle φ-ODE after the
            # mechanical step (φ is decoupled from positions; the ζ(φ)
            # coupling enters next step's `_p2g_momentum_and_stress`).
            self._integrate_phi_ode()

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
            stress_vol = K * (rho_ref / rho_p - 1.0) * ti.Matrix.identity(ti.f32, 3)
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
        # z = 0 of magnitude γ_sub_star · κ_sub_proxy · dt / ρ_local with
        # κ_sub_proxy = 1/dx* and n̂_sub = −ẑ.
        gamma_sub_star = self.cfg.gamma_sub_star
        inv_dx = 1.0 / dx

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

            # v15 (k.3): volumetric strain energy uses the density-based form
            #   U_vol = (1/2) K (ρ_ref/ρ_kernel − 1)²
            # to match the v15 constitutive law σ_vol = K (ρ_ref/ρ_kernel − 1) I.
            # F is still updated for diagnostics + Maxwell deviatoric coupling,
            # but no longer drives the volumetric channel. The deviatoric
            # contribution τ²/(4μ) is unchanged from v12.
            rho_p_d = ti.cast(self._rho_kernel_p[p], ti.f64)
            rho_ref_d = ti.cast(self._rho_ref_kernel[None], ti.f64)
            vol_strain = rho_ref_d / rho_p_d - 1.0
            tau = self.tau_dev[p]
            tau_norm_sq = 0.0
            for ii, jj in ti.static(ti.ndrange(3, 3)):
                tau_norm_sq += tau[ii, jj] * tau[ii, jj]
            U = 0.5 * K * vol_strain * vol_strain + ti.cast(tau_norm_sq, ti.f64) / (4.0 * (mu + 1e-12))
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
            Gate ∈ [0.85, 1.15]: witnesses the kernel sees the substrate via
            the reflective BC (substrate-side cells contribute particle mass
            normally), distinct from the free surface where ρ_kernel ≈ 0.5·ρ_ref
            (Adami-Hu-Adams 2010 §3 truncation bound).
        - F_substrate_per_step: substrate reaction force from accumulated
            impulse / dt (over the most recent step).
        - F_pressure_down: bulk pressure pushing down on the substrate,
            ≈ Σ_p P_p · V₀/h_band over contact-band particles, where
            P_p = K · (1 − ρ_ref/ρ_kernel_p) is the v15 (k.3) hydrostatic
            pressure and h_band = n_contact_band · dx*.
        - anchor_force_balance_rel_err: |F_substrate − F_pressure| /
            max(|F_substrate|, ε); gate ≤ 0.20 (finite-relaxation tolerance).
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
        # Hydrostatic pressure per particle (positive when ρ > ρ_ref, i.e.
        # compressed). σ_vol = K(ρ_ref/ρ − 1)·I; pressure = -tr(σ_vol)/3 =
        # K(1 − ρ_ref/ρ).
        P_per_p = K * (1.0 - rho_ref / np.clip(rho_band, 1e-30, None))
        # Force pushing down on the substrate from each contact-band
        # particle: P · (V0 / h_band). Sum gives total downward bulk force.
        F_pressure_down = float((P_per_p * V0 / h_band).sum())
        F_substrate_up = float(self.diag_substrate_impulse_z[None]) / self.cfg.dt_star

        denom = max(abs(F_substrate_up), 1e-12)
        balance_err = abs(F_substrate_up - F_pressure_down) / denom

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
            "F_pressure_down": F_pressure_down,
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
        }
