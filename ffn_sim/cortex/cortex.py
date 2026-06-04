"""H.3 Cortex multi-filament topology + HOOMD wiring.

Phase 1 Unit H.3 (PHASE_0_3_DECISIONS D3 + D4 + D7; Plan v2 §3 H.3 v3.1
ratified ×40 mesoscopic convention).

Builds a 3D HOOMD-blue cortex of ``n_filaments`` effective F-actin
filaments placed on a spherical shell of radius ``R_cell``. Each
effective filament represents ``×40`` native actin filaments bundled
(the only sanctioned coarse-graining per CLAUDE.md, ratified in Plan v2
§3 H.3 v3.1, 2026-05-19): the per-filament Hamiltonian constants are
the H.2 single-actin-bundle values (KU-1.1 κ, AFINES μ) — what is
coarse-grained is the *filament count* (1000 effective instead of
~38,000 native), not the per-filament stiffness.

Topology layout:

* ``n_filaments`` filaments, each with ``beads_per_filament`` beads at
  rest length ``ℓ₀`` (Plan v2 §3: 1000 × 7 ≈ 7,000 cortex beads).
* Filament center-of-mass placed uniformly on the sphere ``r = R_cell``
  (Marsaglia 1972 surface-uniform sampling).
* Tangent direction sampled uniformly within the local tangent plane.
* Beads laid along the tangent at ℓ₀ spacing. Bonds are at exact ℓ₀
  (force-free at construction); the small radial drift at filament
  endpoints is ``ℓ_end²/(2 R_cell)`` ≈ 110 nm for a 3 μm filament on a
  10 μm shell — within the biological 200 nm cortex thickness band
  (KU-3.17) so no radial projection is needed.

Per-bead Hamiltonian wired into HOOMD (identical to H.1 / H.2 except
optional crosslinker bonds, see §Crosslinkers below)::

    H = ½ k_bond Σ (|b_i| − ℓ₀)²              (md.bond.Harmonic, k = μ/ℓ₀)
      + ½ k_θ   Σ (θ_j − π)²                  (md.angle.Harmonic, t0 = π)
      + Σ_{i<j} LJ_WCA(r_ij; ε, σ, r_cut)    (D7, repulsive only)
      + with_crosslinkers ? ½ k_xl Σ (|r_ab| − r0_bin)²  (KU-3.19 xl bonds)

This unit's deliverable is the TOPOLOGY GENERATOR + AFINES-coarsened
per-filament force constants + optional static crosslinker bonds.
Dynamic D2 Bell-Evans slip kinetics, D5 Stam-Hocky myosin bipolar
minifilaments, and D7 ERM tether are LATER H.3 deliverables in the
brief (`ffn_sim/docs/briefs/H3_cortex.md` §Implementation spec) that
will live in sibling modules (`crosslinkers.py`, `myosin.py`, `erm.py`,
`cell/cell.py`). This module is the load-bearing topology backbone
that those modules will mount onto.

Box geometry — σ_z vs L_z pre-flight (H.2 slab lesson)
-------------------------------------------------------
H.2 (2026-05-21) re-discovered that artificial slab confinement (small
``L_z`` relative to the natural thermal extent ``σ_perp``) biases the
measured persistence length by tens of percent. For H.3 cortex:

* Per-filament natural thermal transverse extent (free 3D, WLC short
  limit since L < L_p):

      σ_perp ≈ √(L³ / (3 L_p))

  For ``L = 3 μm`` (Plan v2 mean), ``L_p = 17 μm`` (KU-1.1):
  ``σ_perp ≈ 730 nm``. For ``L = 5 μm`` (Plan v2 max):
  ``σ_perp ≈ 1.56 μm``.

* Cortex thickness 200 nm (KU-3.17) is a BIOLOGICAL confinement
  enforced physically by the ERM tether (later H.3 deliverable), NOT a
  box-geometry constraint. The TOPOLOGY GENERATOR places filaments on
  a shell but does not pin them to the shell at runtime.

* Box: ``L_box = 3 · R_cell = 30 μm`` cube. Half-box = 15 μm; the cell
  + thermal margin (R_cell + σ_perp_max ≈ 10 + 1.56 = 11.56 μm) sits
  comfortably inside. ``L_z / σ_perp_max ≈ 19×`` — well above the H.2
  failure ratio (Lz=0.2 μm → L_z/σ_perp ≈ 0.14 was the slab failure;
  L_z/σ_perp ≥ 10 is the H.2 strict-PASS safe regime).

* Conclusion: H.2 slab artifact CANNOT occur here. The 200 nm cortex
  thickness emerges from physics (ERM) once attached, not from box
  geometry.

Sanity Gate
-----------
*Per CLAUDE.md hard rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module." STATIC checks live in
``ffn_sim/tests/test_cortex.py``; RUNTIME checks are enforced in
``resolve_h3_derived`` and ``build_cortex_state``.*

1. **Dimensional analysis**

   - ``k_bond = μ / ℓ₀`` → [N/m]. HOOMD ``md.bond.Harmonic.params['k']``
     in SI is J/m² = N/m. ✓
   - ``k_θ = κ_B / ℓ₀`` → [N·m] = [J/rad²] (angle dimensionless). ✓
   - ``γ_b = 6π η R_bead`` → [Pa·s·m] = [N·s/m]. ✓
   - ``ε_LJ = 0.5 kT`` → [J]; ``σ_LJ = 2 R_bead`` → [m]; ``r_cut =
     2^(1/6) σ_LJ`` → [m]. ✓
   - ``k_xl = 0.1 pN/μm = 1e-7 N/m`` (KU-3.19, Furuike 2001/Ferrer
     2008) → [N/m]. ✓
   - ``dt_CFL = α · min(τ_xl, τ_stretch, τ_bend)`` → [s]. ✓
   - STATIC: ``test_cortex.py::TestDimensional`` recomputes each
     derived value from primaries and asserts SI consistency.

2. **Boundary cases**

   - ``beads_per_filament < 2``: forbidden. RUNTIME raise.
   - ``n_filaments == 0``: degenerate (no cortex). RUNTIME raise.
   - ``R_cell ≤ 0``: invalid geometry. RUNTIME raise.
   - ``cortex_thickness ≤ 0``: invalid. RUNTIME raise.
   - ``L_filament = (beads_per_filament − 1) · ℓ₀ > 2 R_cell`` (a
     filament longer than the cell diameter cannot lie tangent on the
     shell). RUNTIME raise.
   - ``L_box < 2 · (R_cell + σ_perp_max + ℓ₀)`` (cell + thermal margin
     escapes the periodic box). RUNTIME raise unless ``demo_mode``.
   - ``n_filaments > acceptance.n_filaments_max`` (cost ceiling, Plan
     v2 §11). RUNTIME raise unless ``demo_mode``.

3. **Conservation invariants**

   - **At construction (all bonds at rest length, all chains straight
     along tangent vectors)**: HOOMD ``potential_energy ≈ 0`` (each
     bond ½·k·0² = 0; each angle ½·k_θ·0² = 0 at θ = π straight
     chain). STATIC test asserts.
   - **Newton's 3rd law**: HOOMD bond/angle/LJ forces are conservative
     by construction; ``net_force.sum() ≈ 0`` to float64 at the
     initial state. STATIC.
   - **Particle count** = ``n_filaments · beads_per_filament``
     (uniform per-filament count). STATIC.
   - **Bond count** = ``n_filaments · (beads_per_filament − 1)``
     exactly (cortex backbone bonds only; xl bonds appended on top
     when ``with_crosslinkers=True``). STATIC.
   - **Angle count** = ``n_filaments · (beads_per_filament − 2)``
     exactly. STATIC.
   - **Crosslinker count** (when enabled) ≤ ``cortex.crosslinkers.n_xl``
     target (some xl requests may fail to find an acceptor within
     ``max_bind_dist`` and are dropped). The actual count is reported
     in ``ResolvedH3.extras['n_xl_realised']``.

4. **Numerical sanity**

   - ``dt_CFL > 0`` and ≤ ``α · τ_min``. RUNTIME assertion.
   - All derived values finite. RUNTIME ``np.isfinite``.
   - Positions ``float64`` (HOOMD CPU default). STATIC.
   - All bead positions within ``[-L_box/2, L_box/2)`` per axis. The
     cell ``r ≤ R_cell + σ_perp_max ≈ 11.6 μm`` < half-box 15 μm so
     there is no wrap event for any bead at construction. RUNTIME
     assertion.
   - Per-filament bond length at construction within
     ``[ℓ₀ − 1e-12, ℓ₀ + 1e-12]`` m (exact-ℓ₀ tangent placement).
     STATIC.

5. **Sign / sense**

   - Bond stretched by +δ along ``b̂``: HOOMD harmonic force on the
     two beads points along ±``b̂`` with magnitude ``k_bond · δ``
     (attractive). STATIC.
   - Bent triplet (perturb interior bead by +δ ŷ, central angle
     θ_HOOMD drops below π): angle force restores toward straight.
     STATIC.
   - LJ WCA repulsive only: two cortex beads at r < σ feel force
     pushing them apart. STATIC.
   - Crosslinker (when enabled): bond stretched by +δ from r0_bin:
     force ``k_xl · δ`` attractive between the two heads, identical
     sense to backbone bonds. STATIC.

6. **Measurement protocol**

   - **Topology smoke** (STATIC): particle / bond / angle counts match
     §3; all bead positions within ±200 nm cortex thickness band at
     construction (with small radial drift ``ℓ_end²/(2 R_cell)`` from
     tangent-plane placement, < 200 nm for ``L ≤ 4 μm``); sphere
     coverage area density ``n_filaments / (4π R_cell²) ≈ 0.8 /μm²``
     consistent with KU-3.17 cortex density.
   - **Per-filament L_p measurement**: aggregate over all filaments
     (same beads_per_filament) via
     ``ffn_sim.common.filament_math.fit_persistence_length`` on the
     ``(n_filaments, beads_per_filament, 3)`` snapshot array.
     First-principles band ``L_p ∈ [15.3, 18.7] μm`` (KU-1.1 actin
     anchor, identical to H.2). NOT a fitted band — the literature
     reference is the gate.
   - **Per-bond equipartition (3D)**: ⟨½ k_θ (θ − π)²⟩ averaged over
     all interior beads × snapshots should equal ``kT/2`` to within
     ±5 % (same first-principles tolerance as H.2 ✅ strict-PASS at
     equipartition_target_3d_kT). Aggregate sample size for 1000
     filaments × 5 interior beads × 100 snapshots = 500,000 ≫ 1/0.05²
     = 400 so the band is statistically tight by 1000×.
   - **3D Boltzmann angle KS**: per-filament angle histograms
     aggregated, KS test against
     ``boltzmann_angle_density_3d``. Same first-principles gate as
     H.2 (KS_stat ≤ angle_ks_stat_max, χ²/df ≤
     angle_chi2_per_df_max, KL ≤ angle_KL_nats_max). NOT a fitted
     gate — the same H.2 thresholds carry over because the underlying
     physics (single-filament bond + angle Brownian dynamics) is
     identical at the per-filament level.

References
----------
- Brief: ``ffn_sim/docs/briefs/H3_cortex.md`` §Cortex topology,
  §Crosslinkers (D2), §Validation acceptance.
- PHASE_0_3_DECISIONS.md §§D3, D4, D7.
- Plan v2 §3 H.3 v3.1 (2026-05-19): ×40 mesoscopic ratification.
- KU-1.1 (κ_B, ℓ_p actin), KU-1.26 (CFL τ_min, η_water), KU-3.5
  (γ_cortex), KU-3.17 (cortex thickness, density, R_cell),
  KU-3.19 (α-actinin / filamin xl mix, lengths, k_xl).
- ``ffn_sim/integrator/baoab.py`` (D3 BAOAB frozen 2026-05-20).
- ``ffn_sim/common/filament_math.py`` (per-filament L_p /
  equipartition / 3D Boltzmann gates — reused from H.2).
- AFINES μ_actin_bundle = 1.5 nN (Freedman 2017 Soft Matter) — same
  as H.2 since H.3 uses the H.2 single-filament force constants per
  effective ×40 bundle.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedH3:
    """Fully-resolved H.3 cortex parameter set (all SI).

    Populated by :func:`resolve_h3_derived` from the raw YAML ``cortex:``
    block. Every primary scale is verbatim from config; every
    ``derived`` value is computed from primaries.
    """

    # Primary biological scales (KU-1.x / KU-3.x anchors)
    R_cell: float                  # m   cell radius, KU-3.17
    cortex_thickness: float        # m   200 nm shell, KU-3.17
    n_filaments: int               # ×40 mesoscopic count, Plan v2 §3 H.3 v3.1
    beads_per_filament: int        # uniform per-filament (Plan v2 mean 7)
    L_filament: float              # m   per-filament contour length (mean 3 μm)
    bead_radius: float             # m   30 nm ×40 actin bundle cross-section
    persistence_length: float      # m   ℓ_p, KU-1.1 actin
    bending_modulus: float         # N·m²  κ_B = ℓ_p k_B T at 310 K, KU-1.1
    stretching_modulus: float      # N    μ_actin_bundle (AFINES Freedman 2017)

    temperature: float             # K   37 °C
    kT: float                      # J   k_B · T
    water_viscosity: float         # Pa·s  η_water at 310 K, KU-1.26

    seed: int

    # Box geometry (cube; σ_z vs L_z pre-flight in module docstring)
    box_factor: float              # L_box = box_factor · R_cell (default 3)

    # D7 LJ (WCA repulsive only, same as H.1)
    lj_enabled: bool
    lj_epsilon_kT: float           # ε / kT (= 0.5 per D7 Phase 1)
    lj_sigma_factor: float         # σ / R_bead (= 2 per D7 Phase 1)

    # Crosslinkers (KU-3.19 — α-actinin / filamin mix, static bonds)
    xl_enabled: bool
    xl_n_target: int               # KU-3.19 target xl count per cortex
    xl_stiffness: float            # N/m   k_xl, KU-3.19 (0.1 pN/μm = 1e-7 N/m)
    xl_alpha_fraction: float       # KU-3.19 (0.3)
    xl_alpha_length: float         # m   α-actinin xl length, KU-3.19 (35 nm)
    xl_filamin_length: float       # m   filamin xl length, KU-3.19 (150 nm)
    xl_max_bind_dist: float        # m   acceptor search radius (60 nm, brief)
    xl_n_bins: int                 # per-r0 binning count (matches H.1 cross_links)

    # D3 dynamics
    integrator_name: str
    cfl_safety_factor: float

    # Acceptance (first-principles bands ONLY — no measurement-anchored tuning)
    L_p_band_m: tuple[float, float]               # KU-1.1 literal [15.3, 18.7] μm
    equipartition_rel_tol: float                  # H.2 first-principles (0.05)
    equipartition_target_3d_kT: float             # 3D solid-angle exact (H.2 ✅)
    angle_chi2_per_df_max: float                  # H.2 strict gate, first-principles
    angle_KL_nats_max: float                      # H.2 strict gate, first-principles
    angle_ks_stat_max: float                      # H.2 diagnostic, first-principles
    radial_drift_band_m: tuple[float, float]      # construction radial position band
    n_filaments_max: int                          # cost ceiling, Plan v2 §11

    demo_mode: bool

    # Derived (filled by resolve_h3_derived)
    L_box: float = 0.0
    rest_length: float = 0.0       # ℓ₀  [m]
    bond_k: float = 0.0            # μ / ℓ₀  [N/m]
    angle_k: float = 0.0           # κ_B / ℓ₀  [N·m / rad²]
    angle_t0: float = math.pi
    lj_epsilon: float = 0.0        # [J]
    lj_sigma: float = 0.0          # [m]
    lj_r_cut: float = 0.0          # [m]
    gamma_b: float = 0.0           # 6π η R_bead  [N·s/m]
    tau_xl: float = 0.0            # γ_b / k_xl  [s]
    tau_stretch: float = 0.0       # γ_b ℓ₀ / μ  [s]
    tau_bend: float = 0.0          # γ_b ℓ₀³ / κ_B  [s]
    tau_min: float = 0.0
    dt_cfl: float = 0.0

    sigma_perp_per_filament: float = 0.0   # √(L³/(3 L_p)) thermal transverse extent
    cortex_surface_area: float = 0.0       # 4π R_cell²  [m²]
    filament_areal_density: float = 0.0    # n_filaments / area  [1/m²]
    n_xl_alpha_target: int = 0             # round(xl_n_target · xl_alpha_fraction)
    n_xl_filamin_target: int = 0

    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Crosslinker helpers (per-r0 binning, same pattern as H.1 cross_links)
# ---------------------------------------------------------------------------
def xl_bin_type_names(n_bins: int) -> list[str]:
    """HOOMD bond type names for per-r0 crosslinker binning."""
    return [f"xl_b{i}" for i in range(n_bins)]


def xl_bin_rest_lengths(n_bins: int, max_bind_dist: float) -> np.ndarray:
    """Bin-center rest lengths uniformly across [0, max_bind_dist].

    Per-bin r0 keeps the construction state approximately force-free
    for crosslinker bonds (H.1 cross_links.py PI-ratified pattern).
    """
    edges = np.linspace(0.0, max_bind_dist, n_bins + 1, dtype=np.float64)
    return 0.5 * (edges[:-1] + edges[1:])


# ---------------------------------------------------------------------------
# Resolver
# ---------------------------------------------------------------------------
def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0, got {x!r}")


def resolve_h3_derived(cfg: dict) -> ResolvedH3:
    """Resolve all derived parameters from a raw ``cortex:`` config block.

    Input ``cfg`` is the YAML root dict (as returned by ``yaml.safe_load``)
    or its ``cortex:`` sub-dict. RUNTIME §1, §2, §4 gate checks fire here.
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]

    accept = cfg["acceptance"]
    xl_cfg = cfg.get("crosslinkers", {})

    p = ResolvedH3(
        # Primary
        R_cell=float(cfg["R_cell"]),
        cortex_thickness=float(cfg["cortex_thickness"]),
        n_filaments=int(cfg["n_filaments"]),
        beads_per_filament=int(cfg["beads_per_filament"]),
        L_filament=float(cfg["L_filament"]),
        bead_radius=float(cfg["bead_radius"]),
        persistence_length=float(cfg["persistence_length"]),
        bending_modulus=float(cfg["bending_modulus"]),
        stretching_modulus=float(cfg["stretching_modulus"]),
        temperature=float(cfg["temperature"]),
        kT=float(cfg["kT"]),
        water_viscosity=float(cfg["water_viscosity"]),
        seed=int(cfg["seed"]),
        box_factor=float(cfg["box"]["L_box_over_R_cell"]),
        # D7 LJ
        lj_enabled=bool(cfg["excluded_volume"]["enabled"]),
        lj_epsilon_kT=float(cfg["excluded_volume"]["epsilon_kT"]),
        lj_sigma_factor=float(cfg["excluded_volume"]["sigma_factor"]),
        # Crosslinkers
        xl_enabled=bool(xl_cfg.get("enabled", False)),
        xl_n_target=int(xl_cfg.get("n_xl", 1000)),
        xl_stiffness=float(xl_cfg.get("k_xl", 1.0e-7)),
        xl_alpha_fraction=float(xl_cfg.get("alpha_fraction", 0.3)),
        xl_alpha_length=float(xl_cfg.get("alpha_length", 35.0e-9)),
        xl_filamin_length=float(xl_cfg.get("filamin_length", 150.0e-9)),
        xl_max_bind_dist=float(xl_cfg.get("max_bind_dist", 60.0e-9)),
        xl_n_bins=int(xl_cfg.get("n_bins", 10)),
        # Dynamics
        integrator_name=str(cfg["dynamics"]["integrator"]),
        cfl_safety_factor=float(cfg["dynamics"]["cfl_safety_factor"]),
        # Acceptance (first-principles bands)
        L_p_band_m=tuple(accept["L_p_band_m"]),
        equipartition_rel_tol=float(accept["equipartition_rel_tol"]),
        equipartition_target_3d_kT=float(accept["equipartition_target_3d_kT"]),
        angle_chi2_per_df_max=float(accept["angle_chi2_per_df_max"]),
        angle_KL_nats_max=float(accept["angle_KL_nats_max"]),
        angle_ks_stat_max=float(accept["angle_ks_stat_max"]),
        radial_drift_band_m=tuple(accept["radial_drift_band_m"]),
        n_filaments_max=int(accept["n_filaments_max"]),
        demo_mode=bool(cfg.get("demo_mode", False)),
    )

    # ---- §2 boundary checks on primaries ----
    if p.beads_per_filament < 2:
        raise ValueError(
            f"beads_per_filament must be ≥ 2 (got {p.beads_per_filament})"
        )
    if p.n_filaments <= 0:
        raise ValueError(f"n_filaments must be > 0 (got {p.n_filaments})")
    if p.integrator_name not in ("leimkuhler_matthews_baoab", "lm_baoab"):
        raise ValueError(
            "dynamics.integrator must be the D3 canonical "
            f"'leimkuhler_matthews_baoab' or 'lm_baoab' "
            f"(got {p.integrator_name})."
        )
    for name, x in [
        ("R_cell", p.R_cell), ("cortex_thickness", p.cortex_thickness),
        ("L_filament", p.L_filament), ("bead_radius", p.bead_radius),
        ("persistence_length", p.persistence_length),
        ("bending_modulus", p.bending_modulus),
        ("stretching_modulus", p.stretching_modulus),
        ("temperature", p.temperature), ("kT", p.kT),
        ("water_viscosity", p.water_viscosity),
        ("box_factor", p.box_factor),
        ("cfl_safety_factor", p.cfl_safety_factor),
    ]:
        _require_finite_positive(name, x)
    if p.L_filament > 2.0 * p.R_cell:
        raise ValueError(
            f"L_filament = {p.L_filament:.3e} m exceeds cell diameter "
            f"2·R_cell = {2 * p.R_cell:.3e} m — a single tangent filament "
            "cannot wrap the cell."
        )

    if p.xl_enabled:
        for name, x in [
            ("xl_stiffness", p.xl_stiffness),
            ("xl_alpha_length", p.xl_alpha_length),
            ("xl_filamin_length", p.xl_filamin_length),
            ("xl_max_bind_dist", p.xl_max_bind_dist),
        ]:
            _require_finite_positive(name, x)
        if not (0.0 <= p.xl_alpha_fraction <= 1.0):
            raise ValueError(
                f"xl_alpha_fraction must be in [0, 1] "
                f"(got {p.xl_alpha_fraction})"
            )
        if p.xl_n_bins < 1:
            raise ValueError(f"xl_n_bins must be ≥ 1 (got {p.xl_n_bins})")

    # ---- §1 derived: geometry + force constants ----
    p.L_box = p.box_factor * p.R_cell
    p.rest_length = p.L_filament / (p.beads_per_filament - 1)
    p.bond_k = p.stretching_modulus / p.rest_length
    p.angle_k = p.bending_modulus / p.rest_length
    p.angle_t0 = math.pi
    p.lj_sigma = p.lj_sigma_factor * p.bead_radius
    p.lj_epsilon = p.lj_epsilon_kT * p.kT
    p.lj_r_cut = (2.0 ** (1.0 / 6.0)) * p.lj_sigma
    p.gamma_b = 6.0 * math.pi * p.water_viscosity * p.bead_radius

    # ---- §1 derived: relaxation times + CFL ----
    p.tau_stretch = p.gamma_b * p.rest_length / p.stretching_modulus
    p.tau_bend = p.gamma_b * p.rest_length**3 / p.bending_modulus
    if p.xl_enabled:
        p.tau_xl = p.gamma_b / p.xl_stiffness
        p.tau_min = min(p.tau_xl, p.tau_stretch, p.tau_bend)
    else:
        p.tau_xl = float("inf")
        p.tau_min = min(p.tau_stretch, p.tau_bend)
    p.dt_cfl = p.cfl_safety_factor * p.tau_min

    # ---- §1 derived: physical predictions ----
    p.sigma_perp_per_filament = math.sqrt(
        p.L_filament**3 / (3.0 * p.persistence_length)
    )
    p.cortex_surface_area = 4.0 * math.pi * p.R_cell**2
    p.filament_areal_density = p.n_filaments / p.cortex_surface_area
    if p.xl_enabled:
        p.n_xl_alpha_target = int(round(p.xl_n_target * p.xl_alpha_fraction))
        p.n_xl_filamin_target = p.xl_n_target - p.n_xl_alpha_target

    # ---- §2 box / cost-ceiling acceptance ----
    half_box = 0.5 * p.L_box
    cell_extent = p.R_cell + p.sigma_perp_per_filament + p.rest_length
    if half_box < cell_extent and not p.demo_mode:
        raise ValueError(
            f"L_box/2 = {half_box:.3e} m is smaller than the cell + "
            f"thermal margin = R_cell + σ_perp + ℓ₀ = "
            f"{cell_extent:.3e} m. Increase box.L_box_over_R_cell or "
            "set demo_mode (CLAUDE.md no-gate-loosening)."
        )
    if p.n_filaments > p.n_filaments_max and not p.demo_mode:
        raise ValueError(
            f"n_filaments = {p.n_filaments} exceeds cost ceiling "
            f"{p.n_filaments_max}; set demo_mode to override."
        )

    # ---- §4 finite checks on derived ----
    for name in [
        "L_box", "rest_length", "bond_k", "angle_k", "lj_epsilon",
        "lj_sigma", "lj_r_cut", "gamma_b", "tau_stretch", "tau_bend",
        "tau_min", "dt_cfl", "sigma_perp_per_filament",
        "cortex_surface_area", "filament_areal_density",
    ]:
        v = getattr(p, name)
        if not (math.isfinite(v) and v > 0.0):
            raise RuntimeError(
                f"Derived {name}={v!r} is not finite-positive; "
                "H.3 resolve produced a degenerate value."
            )

    p.extras["sigma_perp_to_Lz_ratio"] = (
        p.L_box / max(p.sigma_perp_per_filament, 1.0e-30)
    )
    p.extras["radial_drift_construction_m"] = (
        ((p.beads_per_filament - 1) * p.rest_length / 2.0) ** 2
        / (2.0 * p.R_cell)
    )
    p.extras["filament_areal_density_per_um2"] = (
        p.filament_areal_density * 1.0e-12
    )

    return p


# ---------------------------------------------------------------------------
# Topology generation
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class CortexTopology:
    """Cortex shell topology — per-filament bead coordinates + bonds + angles.

    Attributes
    ----------
    positions : ndarray, shape (n_filaments, beads_per_filament, 3)
        SI [m]. Each filament i has bead j at
        ``com_i + (j − (N−1)/2) · ℓ₀ · tangent_i``.
    centers_of_mass : ndarray, shape (n_filaments, 3)
        Filament CoM on the sphere ``r = R_cell``.
    tangents : ndarray, shape (n_filaments, 3)
        Unit tangent vectors in the local tangent plane at the CoM.
    bond_groups : ndarray, shape (n_bonds, 2)
        Per-bond (i, i+1) flat indices in the (F·N,) bead array.
    angle_groups : ndarray, shape (n_angles, 3)
        Per-angle (i, i+1, i+2) flat indices.
    """

    positions: np.ndarray
    centers_of_mass: np.ndarray
    tangents: np.ndarray
    bond_groups: np.ndarray
    angle_groups: np.ndarray


def _sample_sphere_surface(rng: np.random.Generator, n: int, R: float) -> np.ndarray:
    """Uniform points on the sphere of radius ``R`` (Marsaglia 1972).

    Returns shape (n, 3), each |p|=R to float64.
    """
    # Marsaglia: sample u ∈ U(-1, 1), v ∈ U(-1, 1) inside unit disk; map to sphere.
    pts = np.empty((n, 3), dtype=np.float64)
    drawn = 0
    while drawn < n:
        batch = max(2 * (n - drawn), 64)
        u = rng.uniform(-1.0, 1.0, batch)
        v = rng.uniform(-1.0, 1.0, batch)
        s = u * u + v * v
        mask = s < 1.0
        u, v, s = u[mask], v[mask], s[mask]
        k = min(len(u), n - drawn)
        if k == 0:
            continue
        u, v, s = u[:k], v[:k], s[:k]
        factor = 2.0 * np.sqrt(1.0 - s)
        pts[drawn:drawn + k, 0] = R * u * factor
        pts[drawn:drawn + k, 1] = R * v * factor
        pts[drawn:drawn + k, 2] = R * (1.0 - 2.0 * s)
        drawn += k
    return pts


def _tangent_plane_basis(normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return two orthonormal tangent vectors per unit normal.

    Uses the deterministic "minimum component" axis trick to avoid the
    degenerate case where the chosen reference axis is collinear with
    the normal. Returns (e1, e2) with shape (n, 3) each.
    """
    n = normals.shape[0]
    e1 = np.empty_like(normals)
    # Pick the world-axis with smallest |component| of normal, cross with normal.
    abs_n = np.abs(normals)
    min_axis = np.argmin(abs_n, axis=1)
    ref = np.zeros_like(normals)
    ref[np.arange(n), min_axis] = 1.0
    e1 = np.cross(ref, normals)
    e1_norm = np.linalg.norm(e1, axis=1, keepdims=True)
    e1 = e1 / np.maximum(e1_norm, 1.0e-30)
    e2 = np.cross(normals, e1)
    return e1, e2


def generate_cortex_topology(
    p: ResolvedH3,
    rng: np.random.Generator | None = None,
    *,
    tangent_bias_axis: np.ndarray | tuple[float, float, float] | None = None,
    tangent_bias_kappa: float = 0.0,
) -> CortexTopology:
    """Place ``n_filaments`` tangent-plane filaments on the R=R_cell shell.

    Topology generation algorithm:

    1. Sample ``n_filaments`` centers-of-mass uniformly on the sphere
       (Marsaglia 1972).
    2. For each CoM, compute the local outward normal ``n̂ = CoM/R_cell``.
    3. Build a deterministic tangent-plane basis (e1, e2) at the CoM
       (deterministic given seed — so the topology is reproducible).
    4. Sample one azimuth ``φ`` → tangent direction
       ``t̂ = cos(φ) e1 + sin(φ) e2``.  Without a bias, ``φ ∈ U(0, 2π)``
       → isotropic tangent field (KU-3.20 isotropic case, S → 0).
    5. Lay ``beads_per_filament`` beads along ``t̂`` at spacing ℓ₀,
       centered on the CoM. Bonds at exact ℓ₀ → bond Harmonic 0 at
       construction; collinear beads → angle Harmonic 0.

    The construction-time radial drift at filament endpoints is
    ``ℓ_end²/(2 R_cell)`` ≈ 110 nm for L = 3 μm on R = 10 μm shell —
    within the KU-3.17 200 nm cortex thickness, no projection needed.

    Nematic alignment (KU-3.20 aligned case)
    ----------------------------------------
    When ``tangent_bias_axis`` is given with ``tangent_bias_kappa > 0``,
    each filament's azimuth is drawn from a von-Mises distribution
    concentrated about the projection of the global bias axis ``g`` onto
    the local tangent plane: ``φ ~ vonMises(φ₀, κ)`` with
    ``φ₀ = atan2(g·e2, g·e1)``. This realises a meridionally-combed
    tangent field (cf. cortical actin alignment along a stress axis). The
    spherical geometry caps the achievable scalar nematic order at
    ``S = 0.5`` in the κ → ∞ limit (a perfectly combed sphere has
    ⟨t_∥²⟩ = 2/3 → Q_max eigenvalue 0.5), comfortably above the KU-3.20
    aligned threshold S > 0.3. ``κ = 0`` (default) reproduces the
    isotropic ``U(0, 2π)`` behaviour exactly.
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)

    F = p.n_filaments
    N = p.beads_per_filament
    L0 = p.rest_length

    centers = _sample_sphere_surface(rng, F, p.R_cell)
    normals = centers / p.R_cell
    e1, e2 = _tangent_plane_basis(normals)
    if tangent_bias_axis is not None and tangent_bias_kappa > 0.0:
        g = np.asarray(tangent_bias_axis, dtype=np.float64)
        g = g / np.maximum(np.linalg.norm(g), 1.0e-30)
        # Preferred azimuth = projection of the global axis into each plane.
        phi0 = np.arctan2(e2 @ g, e1 @ g)
        phi = rng.vonmises(phi0, tangent_bias_kappa)
    else:
        phi = rng.uniform(0.0, 2.0 * math.pi, F)
    tangents = (
        np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
    )
    # Sanity: ensure tangent unit norm to float64.
    t_norm = np.linalg.norm(tangents, axis=1, keepdims=True)
    tangents = tangents / np.maximum(t_norm, 1.0e-30)

    # Bead offsets along tangent: j − (N−1)/2 for j ∈ 0..N-1.
    offsets = (np.arange(N, dtype=np.float64) - 0.5 * (N - 1)) * L0
    # positions[f, j, :] = centers[f] + offsets[j] * tangents[f]
    positions = (
        centers[:, None, :]
        + offsets[None, :, None] * tangents[:, None, :]
    )

    # Bond / angle groups (flat indices, identical pattern to H.1 mikado).
    f_idx = np.arange(F, dtype=np.int64).reshape(F, 1)
    i_b = np.arange(N - 1, dtype=np.int64).reshape(1, N - 1)
    base_b = f_idx * N + i_b
    bond_groups = np.stack([base_b, base_b + 1], axis=-1).reshape(
        F * (N - 1), 2
    )
    if N >= 3:
        i_a = np.arange(N - 2, dtype=np.int64).reshape(1, N - 2)
        base_a = f_idx * N + i_a
        angle_groups = np.stack(
            [base_a, base_a + 1, base_a + 2], axis=-1
        ).reshape(F * (N - 2), 3)
    else:
        angle_groups = np.empty((0, 3), dtype=np.int64)

    return CortexTopology(
        positions=positions,
        centers_of_mass=centers,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
    )


# ---------------------------------------------------------------------------
# Static crosslinker bond generation (optional, KU-3.19)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class CrosslinkerBonds:
    """Static initial crosslinker bonds at construction.

    Attributes
    ----------
    bond_pairs : ndarray, shape (n_realised, 2)
        Flat-index bead pairs for each xl bond.
    bond_types : ndarray, shape (n_realised,)
        Per-bond r0-bin type index in [0, n_bins).
    species : ndarray of "alpha" or "filamin" strings, shape (n_realised,)
        Diagnostic species label (KU-3.19 30 / 70 split); not used by
        the HOOMD bond force (which is r0-binned, single k for all xl).
    n_alpha : int
    n_filamin : int
    """

    bond_pairs: np.ndarray
    bond_types: np.ndarray
    species: np.ndarray
    n_alpha: int
    n_filamin: int


def generate_static_xl_bonds(
    topology: CortexTopology, p: ResolvedH3,
    rng: np.random.Generator | None = None,
) -> CrosslinkerBonds:
    """Place ``xl_n_target`` static crosslinker bonds, per-bin r0 binning.

    Per-xl algorithm:

    1. Sample a head A: random bead index in [0, F·N).
    2. Find all beads B from DIFFERENT filaments with
       ``|r_A − r_B| ≤ max_bind_dist``.
    3. If at least one acceptor exists, sample one uniformly; else drop
       this xl request (counted in extras['n_xl_dropped']).
    4. Bin the observed |r_AB| into one of ``xl_n_bins`` bins across
       [0, max_bind_dist]. The xl bond type is that bin index.
    5. Label the species as "alpha" or "filamin" by random with
       fraction ``xl_alpha_fraction`` (KU-3.19 30 / 70).

    The per-bin r0 binning ensures construction-time xl forces are
    bounded by ``½ k_xl (bin_width/2)²`` ≪ kT (for k_xl = 1e-7 N/m and
    bin_width = max_bind_dist / n_bins = 6 nm, max xl energy at
    construction = ½ · 1e-7 · (3e-9)² = 4.5e-25 J ≪ kT = 4.28e-21 J).
    """
    if not p.xl_enabled:
        return CrosslinkerBonds(
            bond_pairs=np.empty((0, 2), dtype=np.int64),
            bond_types=np.empty((0,), dtype=np.uint32),
            species=np.empty((0,), dtype="<U7"),
            n_alpha=0, n_filamin=0,
        )

    if rng is None:
        rng = np.random.default_rng(p.seed + 1)

    F = p.n_filaments
    N = p.beads_per_filament
    flat_pos = topology.positions.reshape(F * N, 3)
    # Filament index per flat bead (for "different filament" filter).
    filament_idx = np.repeat(np.arange(F, dtype=np.int64), N)

    bin_width = p.xl_max_bind_dist / p.xl_n_bins
    pairs: list[tuple[int, int]] = []
    types: list[int] = []
    species: list[str] = []
    n_alpha = 0
    n_filamin = 0
    n_dropped = 0

    # For efficiency at larger F, batch random head A draws and check
    # neighbours via direct distance compute (cortex is sparse enough at
    # 60 nm bind radius that brute force is fine for n_xl ≤ 5000).
    max_attempts = p.xl_n_target * 3  # allow up to 3× over-sampling for drops
    attempts = 0
    while len(pairs) < p.xl_n_target and attempts < max_attempts:
        attempts += 1
        a = int(rng.integers(0, F * N))
        a_fil = filament_idx[a]
        r_a = flat_pos[a]
        # Different-filament beads within max_bind_dist.
        diff = flat_pos - r_a
        dist = np.linalg.norm(diff, axis=1)
        mask = (filament_idx != a_fil) & (dist <= p.xl_max_bind_dist)
        candidates = np.nonzero(mask)[0]
        if candidates.size == 0:
            n_dropped += 1
            continue
        b = int(rng.choice(candidates))
        d_ab = dist[b]
        # Bin index: floor(d_ab / bin_width), clamped to [0, n_bins-1].
        bin_idx = int(min(p.xl_n_bins - 1, max(0, int(d_ab / bin_width))))
        species_label = (
            "alpha" if rng.random() < p.xl_alpha_fraction else "filamin"
        )
        if species_label == "alpha":
            n_alpha += 1
        else:
            n_filamin += 1
        pairs.append((a, b))
        types.append(bin_idx)
        species.append(species_label)

    bond_pairs = np.array(pairs, dtype=np.int64).reshape(-1, 2)
    bond_types = np.array(types, dtype=np.uint32)
    species_arr = np.array(species, dtype="<U7")

    p.extras["n_xl_realised"] = int(bond_pairs.shape[0])
    p.extras["n_xl_dropped"] = int(n_dropped)
    p.extras["n_xl_attempts"] = int(attempts)

    return CrosslinkerBonds(
        bond_pairs=bond_pairs,
        bond_types=bond_types,
        species=species_arr,
        n_alpha=n_alpha,
        n_filamin=n_filamin,
    )


# ---------------------------------------------------------------------------
# GSD frame + Simulation builders
# ---------------------------------------------------------------------------
def build_cortex_state(
    p: ResolvedH3, *,
    with_crosslinkers: bool | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[gsd.hoomd.Frame, CortexTopology, CrosslinkerBonds]:
    """Build a HOOMD GSD frame from a resolved H.3 cortex config.

    Returns
    -------
    snap : gsd.hoomd.Frame
        HOOMD initial frame with cortex beads + cortex-bond + cortex-angle
        + optional xl-bin bond types.
    topology : CortexTopology
        Per-filament geometry (for downstream measurement, visualization).
    xl_bonds : CrosslinkerBonds
        Per-xl bond pairs / types / species. Empty when disabled.
    """
    enable_xl = p.xl_enabled if with_crosslinkers is None else bool(with_crosslinkers)

    topology = generate_cortex_topology(p, rng=rng)
    # If xl override forces enabled but config says disabled (or vice versa),
    # honour the override on the realised bonds only — but the xl bond
    # type registration in the frame depends on the resolved p.
    if enable_xl and not p.xl_enabled:
        # Override case: generate with a local temp-enabled copy.
        p_xl = ResolvedH3(**{
            f.name: getattr(p, f.name)
            for f in p.__dataclass_fields__.values()
            if f.name != "extras"
        })
        p_xl.xl_enabled = True
        p_xl.extras = {}
        xl_bonds = generate_static_xl_bonds(topology, p_xl, rng=rng)
    elif (not enable_xl) and p.xl_enabled:
        xl_bonds = CrosslinkerBonds(
            bond_pairs=np.empty((0, 2), dtype=np.int64),
            bond_types=np.empty((0,), dtype=np.uint32),
            species=np.empty((0,), dtype="<U7"),
            n_alpha=0, n_filamin=0,
        )
    else:
        xl_bonds = generate_static_xl_bonds(topology, p, rng=rng)

    F = p.n_filaments
    N = p.beads_per_filament
    n_part = F * N
    flat_pos = topology.positions.reshape(n_part, 3)

    # §4 box-wrap sanity: all coords within [-L_box/2, L_box/2).
    half = 0.5 * p.L_box
    if (np.abs(flat_pos) > half + 1.0e-9).any():
        raise RuntimeError(
            f"Bead positions outside HOOMD box [-{half:.3e}, {half:.3e}); "
            "increase box.L_box_over_R_cell."
        )

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_part
    snap.particles.types = ["actin_cortex"]
    snap.particles.typeid = np.zeros(n_part, dtype=np.uint32)
    snap.particles.position = flat_pos
    snap.particles.mass = np.ones(n_part, dtype=np.float64)

    n_backbone_bonds = topology.bond_groups.shape[0]
    n_xl_bonds = xl_bonds.bond_pairs.shape[0]
    total_bonds = n_backbone_bonds + n_xl_bonds

    # Bond type registration: cortex-bond + per-r0-bin xl types.
    bond_types: list[str] = ["cortex-bond"]
    if enable_xl:
        bond_types.extend(xl_bin_type_names(p.xl_n_bins))

    bond_groups = np.empty((total_bonds, 2), dtype=np.uint32)
    bond_typeids = np.zeros(total_bonds, dtype=np.uint32)
    bond_groups[:n_backbone_bonds] = topology.bond_groups.astype(np.uint32)
    if n_xl_bonds > 0:
        bond_groups[n_backbone_bonds:] = xl_bonds.bond_pairs.astype(np.uint32)
        # xl typeids offset by 1 (cortex-bond occupies typeid 0).
        bond_typeids[n_backbone_bonds:] = xl_bonds.bond_types.astype(np.uint32) + 1

    snap.bonds.N = total_bonds
    snap.bonds.types = bond_types
    snap.bonds.typeid = bond_typeids
    snap.bonds.group = bond_groups

    n_angles = topology.angle_groups.shape[0]
    if n_angles > 0:
        snap.angles.N = n_angles
        snap.angles.types = ["cortex-angle"]
        snap.angles.typeid = np.zeros(n_angles, dtype=np.uint32)
        snap.angles.group = topology.angle_groups.astype(np.uint32)

    snap.configuration.box = [p.L_box, p.L_box, p.L_box, 0.0, 0.0, 0.0]
    return snap, topology, xl_bonds


def build_cortex_simulation(
    p: ResolvedH3, *,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
    with_crosslinkers: bool | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[hoomd.Simulation, Any, Any, CortexTopology, CrosslinkerBonds]:
    """Construct and wire a HOOMD Simulation for the H.3 cortex.

    Returns
    -------
    sim, updater, action, topology, xl_bonds

    ``updater`` and ``action`` are None when ``with_baoab=False``.
    """
    snap, topology, xl_bonds = build_cortex_state(
        p, with_crosslinkers=with_crosslinkers, rng=rng
    )
    enable_xl = (
        p.xl_enabled if with_crosslinkers is None else bool(with_crosslinkers)
    )

    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p.seed
    )
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    bond.params["cortex-bond"] = dict(k=p.bond_k, r0=p.rest_length)
    if enable_xl:
        bin_names = xl_bin_type_names(p.xl_n_bins)
        bin_r0 = xl_bin_rest_lengths(p.xl_n_bins, p.xl_max_bind_dist)
        for i in range(p.xl_n_bins):
            bond.params[bin_names[i]] = dict(
                k=p.xl_stiffness, r0=float(bin_r0[i])
            )

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p.angle_k, t0=p.angle_t0)

    # D7 LJ WCA repulsive only on actin_cortex × actin_cortex.
    # Cortex is denser than ECM (filament areal density ~0.8 /μm² on the
    # shell), so md.nlist.Cell is fine here (box 30 μm, r_cut 60 nm →
    # ~500 cells/axis = 1.25e8 cells — modest). Keep Tree for parity
    # with H.1 anyway: identical neighbour-list semantics, no measurable
    # cost difference for n_part ~ 7e3.
    nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    lj.params[("actin_cortex", "actin_cortex")] = dict(
        epsilon=p.lj_epsilon, sigma=p.lj_sigma
    )
    lj.r_cut[("actin_cortex", "actin_cortex")] = (
        p.lj_r_cut if p.lj_enabled else 0.0
    )
    lj.mode = "shift"

    ig = md.Integrator(dt=p.dt_cfl)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    # methods=[] required for L-M Action contract (D3 BAOAB freeze).
    sim.operations.integrator = ig

    if with_baoab:
        action, updater = make_baoab_updater(
            kT=p.kT,
            gamma={"actin_cortex": p.gamma_b},
            dt=p.dt_cfl,
            seed=p.seed,
        )
        sim.operations.updaters.append(updater)
        return sim, updater, action, topology, xl_bonds

    return sim, None, None, topology, xl_bonds


# ===========================================================================
# Variable-length filament distribution (additive — fixed-N functions above
# are unchanged; new functions below are opt-in via config block
# cortex.variable_length.enabled or direct function call).
# ===========================================================================
#
# Brief §Cortex topology specifies: "Filament length distribution: uniform
# 1–5 μm (mean 3 μm) → average 7 beads per filament at ℓ₀=0.5 μm.  Total
# cortex beads ≈ 7,000 per cell (Plan v2 §3 H.3 v3.1 estimate)."
#
# The fixed-N implementation (above) uses L = 3 μm (the brief's MEAN) for
# all filaments and gives the same total bead count and per-filament force
# constants.  The variable-length implementation below realises the full
# distribution: per-filament L_i ~ Uniform(L_min, L_max), quantized to ℓ_0
# multiples → N_beads_i = round(L_i/ℓ_0) + 1.
#
# Sanity Gate (additive — inherits §1-5 from the fixed-N module-level
# docstring; the only new sense to verify is §3 conservation + §6
# measurement-protocol consistency at the variable-N level.):
#
# §3 Conservation (variable-length):
#   - Particle count = Σ N_beads_i (per-filament sum, not F · N).
#   - Bond count    = Σ (N_beads_i − 1).
#   - Angle count   = Σ (N_beads_i − 2)  with each term ≥ 0 (filaments
#     with N=2 contribute zero angles, as expected by the harmonic
#     angle compute which needs three consecutive beads).
#
# §6 Measurement (variable-length):
#   - Per-filament L_i empirical mean approaches (L_min+L_max)/2 with
#     sample-size N error √(var(L)/F).
#   - Total bead count near n_filaments × ((L_min+L_max)/2 / ℓ_0 + 1).
#   - L_i range fills [L_min, L_max] without clipping (each ℓ_0
#     quantum bin sampled approximately uniformly).


@dataclass(slots=True)
class VariableLengthCortexLayout:
    """Variable-N cortex shell topology — flat-indexed layout.

    Attributes
    ----------
    positions_flat : ndarray, shape (n_total_beads, 3)
        Bead positions in flat order: filament 0 beads first, then
        filament 1, etc.  ``positions_flat[filament_starts[i]:
        filament_starts[i] + n_beads_per_filament[i], :]`` is filament i.
    n_beads_per_filament : ndarray, shape (n_filaments,), dtype int64
        Per-filament bead counts (varies per filament).
    filament_starts : ndarray, shape (n_filaments,), dtype int64
        Flat start indices for each filament.
        ``filament_starts[i] = sum(n_beads_per_filament[0:i])``.
    L_per_filament : ndarray, shape (n_filaments,), dtype float64
        Per-filament contour length actually realised (after ℓ_0
        quantization): ``L_i = (n_beads_per_filament[i] - 1) · ℓ_0``.
    centers_of_mass : ndarray, shape (n_filaments, 3)
    tangents : ndarray, shape (n_filaments, 3)
    bond_groups : ndarray, shape (n_total_bonds, 2)
    angle_groups : ndarray, shape (n_total_angles, 3)
    """

    positions_flat: np.ndarray
    n_beads_per_filament: np.ndarray
    filament_starts: np.ndarray
    L_per_filament: np.ndarray
    centers_of_mass: np.ndarray
    tangents: np.ndarray
    bond_groups: np.ndarray
    angle_groups: np.ndarray
    # Bimodal-architecture diagnostic: per-filament nucleator class (formin =
    # long backbone vs Arp2/3 = short infill).  Empty for the uniform layout.
    is_formin: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=bool)
    )
    # Arp2/3 70° BRANCH topology (mother-daughter junctions).  branch_bonds are
    # (mother_bead, daughter_base_bead) flat-index pairs — the Arp2/3 link that
    # makes a branched daughter BORN attached to its mother (real connectivity,
    # independent of crosslinkers).  branch_angles are (mother_neighbour,
    # mother_bead, daughter_base) triplets constrained to the 70° branch angle.
    # is_branched[f] = True for daughters nucleated off a mother.  Empty when
    # arp_branch_fraction = 0.
    branch_bonds: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=np.int64)
    )
    branch_angles: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=np.int64)
    )
    is_branched: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=bool)
    )

    @property
    def filament_idx(self) -> np.ndarray:
        """Per-bead filament index (shape (n_total_beads,)), derived from the
        per-filament bead counts.  Used by the crosslinker bridge-different-
        filament rule and connectivity measurement."""
        return np.repeat(
            np.arange(self.n_beads_per_filament.shape[0], dtype=np.int64),
            self.n_beads_per_filament,
        )


def generate_variable_length_cortex_layout(
    p: ResolvedH3,
    *,
    L_min: float,
    L_max: float,
    n_filaments: int | None = None,
    rng: np.random.Generator | None = None,
    n_beads_min: int = 2,
    n_beads_max: int | None = None,
) -> VariableLengthCortexLayout:
    """Place ``n_filaments`` variable-length filaments on the R=R_cell shell.

    Each filament:
    1. Center sampled uniformly on the sphere (Marsaglia).
    2. Random tangent-plane axis (uniform azimuth).
    3. ``L_i`` drawn from ``Uniform(L_min, L_max)``, quantized to ℓ_0
       multiples → ``N_beads_i = round(L_i/ℓ_0) + 1`` clipped to
       ``[n_beads_min, n_beads_max]``.
    4. ``N_beads_i`` beads laid along tangent at ℓ_0 spacing, centered
       on the CoM.

    Parameters
    ----------
    L_min, L_max : float
        Filament length range [m].  Brief §Cortex topology says 1-5 μm.
    n_filaments : int, optional
        Defaults to ``p.n_filaments``.
    n_beads_min : int, default 2
        Lower clamp on per-filament bead count (need ≥ 2 for a bond).
    n_beads_max : int, optional
        Upper clamp.  Default ``round(L_max/ℓ_0) + 1``.
    rng : np.random.Generator, optional
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)

    if not (math.isfinite(L_min) and L_min > 0.0):
        raise ValueError(f"L_min must be finite and > 0; got {L_min}")
    if not (math.isfinite(L_max) and L_max >= L_min):
        raise ValueError(f"L_max must be finite and ≥ L_min; got {L_max}")
    if L_max > 2.0 * p.R_cell:
        raise ValueError(
            f"L_max = {L_max:.3e} m exceeds cell diameter 2 R_cell = "
            f"{2 * p.R_cell:.3e} m — a single tangent filament cannot "
            "wrap the cell."
        )

    F = int(n_filaments) if n_filaments is not None else int(p.n_filaments)
    if F <= 0:
        raise ValueError(f"n_filaments must be > 0; got {F}")

    L0 = p.rest_length
    if n_beads_max is None:
        n_beads_max = int(round(L_max / L0)) + 1
    if n_beads_min < 2:
        raise ValueError(f"n_beads_min must be ≥ 2; got {n_beads_min}")
    if n_beads_max < n_beads_min:
        raise ValueError(
            f"n_beads_max ({n_beads_max}) must be ≥ n_beads_min ({n_beads_min})"
        )

    # Per-filament N_beads via Uniform L sampling + ℓ_0 quantization.
    L_samples = rng.uniform(L_min, L_max, F)
    n_beads_per = np.clip(
        np.round(L_samples / L0).astype(np.int64) + 1,
        n_beads_min, n_beads_max,
    )
    # Realised contour length per filament (after quantization).
    L_realised = (n_beads_per - 1).astype(np.float64) * L0

    # CoM on sphere (Marsaglia) + tangent basis.
    centers = _sample_sphere_surface(rng, F, p.R_cell)
    normals = centers / p.R_cell
    e1, e2 = _tangent_plane_basis(normals)
    phi = rng.uniform(0.0, 2.0 * math.pi, F)
    tangents = (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    tangents = tangents / np.linalg.norm(
        tangents, axis=1, keepdims=True
    ).clip(min=1e-30)

    # Flat layout.
    n_total_beads = int(n_beads_per.sum())
    filament_starts = np.zeros(F, dtype=np.int64)
    filament_starts[1:] = np.cumsum(n_beads_per[:-1])
    positions_flat = np.empty((n_total_beads, 3), dtype=np.float64)
    bond_pairs: list[tuple[int, int]] = []
    angle_triplets: list[tuple[int, int, int]] = []

    for f in range(F):
        N_f = int(n_beads_per[f])
        start = int(filament_starts[f])
        offsets = (np.arange(N_f, dtype=np.float64) - 0.5 * (N_f - 1)) * L0
        positions_flat[start:start + N_f] = (
            centers[f] + offsets[:, None] * tangents[f]
        )
        for j in range(N_f - 1):
            bond_pairs.append((start + j, start + j + 1))
        for j in range(N_f - 2):
            angle_triplets.append((start + j, start + j + 1, start + j + 2))

    bond_groups = np.array(bond_pairs, dtype=np.int64).reshape(-1, 2)
    angle_groups = (
        np.array(angle_triplets, dtype=np.int64).reshape(-1, 3)
        if angle_triplets else np.empty((0, 3), dtype=np.int64)
    )

    return VariableLengthCortexLayout(
        positions_flat=positions_flat,
        n_beads_per_filament=n_beads_per,
        filament_starts=filament_starts,
        L_per_filament=L_realised,
        centers_of_mass=centers,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
    )


# ===========================================================================
# Bimodal-exponential cortex topology (CORTEX CONSTRUCTION REBUILD 2026-06-04)
# ===========================================================================
#
# The prior cortex topology (uniform L = 3 μm, fixed 7 beads) is the root cause
# of the fragmented mesh (z = 1.3, giant 7 %, outputs/h3/production/
# cortex_network.json) that underlies the γ-floor: a uniform single-length field
# has NO long connecting backbone, so the crosslinked network does not percolate
# and force is not transmitted (Kadzik-Munro 2026: connectivity is the
# prerequisite for force transmission).
#
# The architecture lit-study (docs/ACTIN_ARCHITECTURE_NOTES.md) converged on the
# faithful cortical-filament length law:
#   * BIMODAL-EXPONENTIAL length (Fritzsche 2016 Sci Adv; Fritzsche 2017 Nat
#     Commun): a SHORT Arp2/3-branched subpopulation (native ~120 nm, ~80-90 %
#     by COUNT) + a LONG formin-linear subpopulation (native ~1200 nm, ~10-20 %
#     by count but ~10× longer).  The long formin filaments are the
#     MECHANICAL + PERCOLATION BACKBONE — a single long filament spans many mesh
#     holes and ties together many distinct crosslink partners.
#   * DISORDERED ISOTROPIC orientation (Li-Gao-Xu 2022: power-law cortex
#     rheology ORIGINATES from exponential disorder; a regular lattice gives the
#     wrong rheology).  Uniform-azimuth tangent placement (S → 0).
#   * 200 nm thick shell BAND (Clark 2013 / Flormann 2024; KU-3.17): every bead
#     is projected onto the band so even the long formin backbone curves ALONG
#     the membrane and the radial extent stays < 200 nm for ANY length.
#
# At the ×40 mesoscale the absolute native nm are inflated with the mesh (the
# only sanctioned coarse-graining is the filament COUNT, CLAUDE.md); what is
# preserved faithfully is the bimodal SHAPE and the long-backbone role.  The
# default means (L_short = ℓ_0, L_long = 10·L_short) keep the native ~10× ratio.
#
# Sanity Gate (additive — inherits §1-5 from the module docstring; new senses):
#   §3 Conservation: particle/bond/angle counts = Σ over variable N_i (as the
#      uniform variable-length layout); is_formin mask sums to ~formin_fraction.
#   §6 Measurement: per-subpopulation mean length recovers L_short / L_long;
#      every bead radius within [R_cell − 200 nm, R_cell] (shell-band project);
#      tangent field isotropic (nematic S → 0).


def _project_to_shell_band(
    positions: np.ndarray, R_cell: float, thickness: float,
    rng: np.random.Generator,
    *,
    bead_depth: np.ndarray | None = None,
) -> np.ndarray:
    """Project bead positions onto the cortex shell BAND [R_cell − thickness,
    R_cell] by radial renormalisation to a per-bead target radius.

    The tangent-plane placement gives long filaments a chord that bulges
    radially inward by ``(L/2)²/(2 R)`` — for the long formin backbone this can
    exceed the 200 nm band.  Renormalising each bead's radius curves the
    filament ALONG the shell (geodesic-like) so it stays within the KU-3.17
    cortex thickness for ANY length, while the disordered isotropic tangent
    field is preserved (only the radial coordinate is touched).

    ``bead_depth`` (per-bead radial depth in [0, thickness]) is supplied by the
    caller so that a whole filament shares ONE depth — projecting it to a single
    sphere of radius ``R_cell − depth`` as a SMOOTH arc (no per-bead radial
    zigzag, which would distort bonds and let crossing filaments coincide).
    Different filaments draw different depths → they are radially separated
    within the band, avoiding construction WCA overlaps.  Falls back to an
    independent per-bead uniform draw when ``bead_depth`` is None.
    """
    r = np.linalg.norm(positions, axis=1, keepdims=True)
    if bead_depth is None:
        depth = rng.uniform(0.0, thickness, (positions.shape[0], 1))
    else:
        depth = np.asarray(bead_depth, dtype=np.float64).reshape(-1, 1)
    r_target = R_cell - depth
    return positions / np.maximum(r, 1.0e-30) * r_target


def _nucleate_arp_branches(
    positions_flat: np.ndarray,
    starts: np.ndarray,
    nbeads: np.ndarray,
    tangents: np.ndarray,
    is_formin: np.ndarray,
    *,
    ell0: float,
    R_cell: float,
    thickness: float,
    arp_branch_fraction: float,
    branch_angle_deg: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Nucleate Arp2/3 daughters as BRANCHES off mother filaments at ~70°.

    Faithful Arp2/3 topology (Fritzsche 2016/2017; Garlick 70° ground-truth):
    a fraction ``arp_branch_fraction`` of the SHORT (Arp2/3) filaments are
    daughters BORN attached to a mother — the daughter is rigidly re-placed so
    its base bead sits one segment from a mother bead, growing at the branch
    angle (≈70°) to the mother's local tangent in the shell tangent plane.

    A branch contributes:
      * a BRANCH BOND ``(mother_bead, daughter_base)`` — the Arp2/3 link (the
        daughter is connected to the network independent of crosslinkers, so no
        branched daughter is ever isolated), and
      * a BRANCH ANGLE ``(mother_neighbour, mother_bead, daughter_base)`` held
        at ``branch_angle_deg`` (the dendritic 70° geometry).

    Mothers are restricted to NON-daughter filaments (one branch level), chosen
    as the nearest such filament bead to the daughter's original location.

    Returns (positions_flat, branch_bonds (B,2), branch_angles (B,3),
    is_branched (F,)).
    """
    from scipy.spatial import cKDTree

    F = nbeads.shape[0]
    is_branched = np.zeros(F, dtype=bool)
    short = np.flatnonzero(~is_formin)
    if arp_branch_fraction <= 0.0 or short.size == 0:
        return (positions_flat, np.empty((0, 2), dtype=np.int64),
                np.empty((0, 3), dtype=np.int64), is_branched)
    n_branch = int(round(arp_branch_fraction * short.size))
    daughters = rng.choice(short, size=n_branch, replace=False)
    is_branched[daughters] = True

    # Mother pool = beads of NON-daughter filaments (formin + unbranched Arp2/3).
    mother_fil = np.flatnonzero(~is_branched)
    mother_beads = np.concatenate([
        np.arange(int(starts[f]), int(starts[f]) + int(nbeads[f]))
        for f in mother_fil
    ])
    tree = cKDTree(positions_flat[mother_beads])

    cos_t = math.cos(math.radians(branch_angle_deg))
    sin_t = math.sin(math.radians(branch_angle_deg))
    branch_bonds: list[tuple[int, int]] = []
    branch_angles: list[tuple[int, int, int]] = []
    fil_of_bead = np.repeat(np.arange(F), nbeads)

    for d in daughters:
        s0, nd = int(starts[d]), int(nbeads[d])
        com = positions_flat[s0:s0 + nd].mean(axis=0)
        # nearest mother bead to the daughter's location (query a few, skip any
        # accidentally on a daughter via the fil-of-bead check).
        _, idxs = tree.query(com, k=min(8, mother_beads.size))
        m_bead = None
        for j in np.atleast_1d(idxs):
            cand = int(mother_beads[int(j)])
            if not is_branched[fil_of_bead[cand]]:
                m_bead = cand
                break
        if m_bead is None:
            is_branched[d] = False
            continue
        fm = int(fil_of_bead[m_bead])
        P_m = positions_flat[m_bead]
        n_hat = P_m / max(np.linalg.norm(P_m), 1e-30)        # shell normal
        t_m = tangents[fm]
        t_m = t_m - np.dot(t_m, n_hat) * n_hat               # in tangent plane
        t_m /= max(np.linalg.norm(t_m), 1e-30)
        w = np.cross(n_hat, t_m)                             # lateral in-plane
        # Branch at the 70° angle from the mother axis, at a RANDOM azimuth ψ
        # around that axis (the Arp2/3 dendritic cone) — so multiple daughters
        # off the same mother point in DIFFERENT directions (no coincident
        # bases).  The out-of-plane component is pulled back by the shell
        # projection below, keeping the quasi-2D cortex.
        psi = rng.uniform(0.0, 2.0 * math.pi)
        branch_dir = cos_t * t_m + sin_t * (
            math.cos(psi) * w + math.sin(psi) * n_hat)
        branch_dir /= max(np.linalg.norm(branch_dir), 1e-30)
        # rigidly re-place the daughter: base one segment from the mother bead,
        # beads laid straight along branch_dir, then projected back onto the band.
        base = P_m + ell0 * branch_dir
        for k in range(nd):
            positions_flat[s0 + k] = base + k * ell0 * branch_dir
        # keep daughter in the band at the mother's depth (smooth, in-shell)
        depth_m = R_cell - np.linalg.norm(P_m)
        sub = positions_flat[s0:s0 + nd]
        r = np.linalg.norm(sub, axis=1, keepdims=True)
        positions_flat[s0:s0 + nd] = sub / np.maximum(r, 1e-30) * (R_cell - depth_m)
        # branch bond + 70° branch angle
        branch_bonds.append((m_bead, s0))
        m_nbr = m_bead + 1 if (m_bead + 1) < (int(starts[fm]) + int(nbeads[fm])) \
            else m_bead - 1
        if m_nbr != m_bead:
            branch_angles.append((m_nbr, m_bead, s0))

    return (
        positions_flat,
        np.array(branch_bonds, dtype=np.int64).reshape(-1, 2),
        np.array(branch_angles, dtype=np.int64).reshape(-1, 3),
        is_branched,
    )


def generate_bimodal_cortex_layout(
    p: ResolvedH3,
    *,
    formin_fraction: float = 0.12,
    L_short_mean: float | None = None,
    L_long_mean: float | None = None,
    n_filaments: int | None = None,
    project_to_shell: bool = True,
    arp_branch_fraction: float = 0.0,
    branch_angle_deg: float = 70.0,
    rng: np.random.Generator | None = None,
    n_beads_min: int = 2,
    n_beads_max: int | None = None,
) -> VariableLengthCortexLayout:
    """Place ``n_filaments`` cortical filaments with a BIMODAL-EXPONENTIAL
    length distribution (the cortex-construction-rebuild topology).

    Per filament:
      1. Nucleator class: formin (long backbone) with prob ``formin_fraction``,
         else Arp2/3 (short).  Fritzsche 2016/2017: ~10-20 % formin by count.
      2. Length ``L_i ~ Exponential(mean)`` with mean ``L_long_mean`` (formin)
         or ``L_short_mean`` (Arp2/3), quantized to ℓ_0 multiples →
         ``N_beads_i = round(L_i/ℓ_0) + 1`` clipped to [n_beads_min, n_beads_max].
      3. Center uniform on the sphere (Marsaglia); tangent in the local plane at
         a UNIFORM azimuth (disordered isotropic, S → 0).
      4. Beads laid at ℓ_0 spacing about the CoM, then (if ``project_to_shell``)
         projected onto the 200 nm shell band so the long backbone curves along
         the membrane.

    Defaults: ``L_short_mean = ℓ_0`` (short Arp2/3 ≈ one segment at the
    mesoscale), ``L_long_mean = 10·L_short_mean`` (formin ~10× longer, Fritzsche).
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)
    if not (0.0 <= formin_fraction <= 1.0):
        raise ValueError(
            f"formin_fraction must be in [0, 1]; got {formin_fraction}"
        )

    L0 = p.rest_length
    if L_short_mean is None:
        L_short_mean = L0
    if L_long_mean is None:
        L_long_mean = 10.0 * L_short_mean
    if not (math.isfinite(L_short_mean) and L_short_mean > 0.0):
        raise ValueError(f"L_short_mean must be finite > 0; got {L_short_mean}")
    if not (math.isfinite(L_long_mean) and L_long_mean >= L_short_mean):
        raise ValueError(
            f"L_long_mean must be finite ≥ L_short_mean; got {L_long_mean}"
        )

    F = int(n_filaments) if n_filaments is not None else int(p.n_filaments)
    if F <= 0:
        raise ValueError(f"n_filaments must be > 0; got {F}")
    if n_beads_min < 2:
        raise ValueError(f"n_beads_min must be ≥ 2; got {n_beads_min}")
    # Upper clamp: a tangent chord can be at most the cell diameter; with shell
    # projection the geodesic arc is bounded by π·R_cell.  Cap conservatively so
    # no filament wraps past the cell.
    if n_beads_max is None:
        L_cap = min(2.0 * p.R_cell, max(L_long_mean * 6.0, 12.0 * L0))
        n_beads_max = max(n_beads_min, int(round(L_cap / L0)) + 1)

    # ---- per-filament nucleator class + exponential length ----
    is_formin = rng.random(F) < formin_fraction
    means = np.where(is_formin, L_long_mean, L_short_mean)
    L_samples = rng.exponential(means)
    n_beads_per = np.clip(
        np.round(L_samples / L0).astype(np.int64) + 1, n_beads_min, n_beads_max
    )
    L_realised = (n_beads_per - 1).astype(np.float64) * L0

    # ---- disordered isotropic placement (uniform azimuth) ----
    centers = _sample_sphere_surface(rng, F, p.R_cell)
    normals = centers / p.R_cell
    e1, e2 = _tangent_plane_basis(normals)
    phi = rng.uniform(0.0, 2.0 * math.pi, F)
    tangents = np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
    tangents = tangents / np.linalg.norm(
        tangents, axis=1, keepdims=True
    ).clip(min=1e-30)

    # ---- flat layout ----
    n_total_beads = int(n_beads_per.sum())
    filament_starts = np.zeros(F, dtype=np.int64)
    filament_starts[1:] = np.cumsum(n_beads_per[:-1])
    positions_flat = np.empty((n_total_beads, 3), dtype=np.float64)
    bond_pairs: list[tuple[int, int]] = []
    angle_triplets: list[tuple[int, int, int]] = []
    for f in range(F):
        N_f = int(n_beads_per[f])
        start = int(filament_starts[f])
        offsets = (np.arange(N_f, dtype=np.float64) - 0.5 * (N_f - 1)) * L0
        positions_flat[start:start + N_f] = (
            centers[f] + offsets[:, None] * tangents[f]
        )
        for j in range(N_f - 1):
            bond_pairs.append((start + j, start + j + 1))
        for j in range(N_f - 2):
            angle_triplets.append((start + j, start + j + 1, start + j + 2))

    if project_to_shell:
        # One depth per FILAMENT (smooth arc on a sphere of radius R−depth);
        # different filaments sit at different depths within the 200 nm band so
        # crossing filaments are radially separated (no construction overlap).
        filament_depth = rng.uniform(0.0, p.cortex_thickness, F)
        bead_depth = np.repeat(filament_depth, n_beads_per)
        positions_flat = _project_to_shell_band(
            positions_flat, p.R_cell, p.cortex_thickness, rng,
            bead_depth=bead_depth,
        )

    # Arp2/3 70° branching: re-place a fraction of short filaments as daughters
    # born attached to a mother (real connectivity → no isolated short filament).
    positions_flat, branch_bonds, branch_angles, is_branched = (
        _nucleate_arp_branches(
            positions_flat, filament_starts, n_beads_per, tangents, is_formin,
            ell0=L0, R_cell=p.R_cell, thickness=p.cortex_thickness,
            arp_branch_fraction=arp_branch_fraction,
            branch_angle_deg=branch_angle_deg, rng=rng,
        )
    )

    bond_groups = np.array(bond_pairs, dtype=np.int64).reshape(-1, 2)
    angle_groups = (
        np.array(angle_triplets, dtype=np.int64).reshape(-1, 3)
        if angle_triplets else np.empty((0, 3), dtype=np.int64)
    )

    return VariableLengthCortexLayout(
        positions_flat=positions_flat,
        n_beads_per_filament=n_beads_per,
        filament_starts=filament_starts,
        L_per_filament=L_realised,
        centers_of_mass=centers,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
        is_formin=is_formin,
        branch_bonds=branch_bonds,
        branch_angles=branch_angles,
        is_branched=is_branched,
    )


def build_variable_length_cortex_state(
    p: ResolvedH3, layout: VariableLengthCortexLayout,
):
    """HOOMD GSD frame builder for variable-length cortex.

    Mirrors :func:`build_cortex_state` but consumes the variable-N
    ``VariableLengthCortexLayout`` instead of the fixed-N CortexTopology.
    No xlinks at this layer (xlinks need fixed per-filament topology to
    look up filament indices — extend xlinks support in a follow-up).
    """
    import gsd.hoomd

    n_total = layout.positions_flat.shape[0]
    half = 0.5 * p.L_box
    if (np.abs(layout.positions_flat) > half + 1.0e-9).any():
        raise RuntimeError(
            f"Variable-length bead positions outside HOOMD box "
            f"[-{half:.3e}, {half:.3e}); increase box.L_box_over_R_cell."
        )

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_total
    snap.particles.types = ["actin_cortex"]
    snap.particles.typeid = np.zeros(n_total, dtype=np.uint32)
    snap.particles.position = layout.positions_flat
    snap.particles.mass = np.ones(n_total, dtype=np.float64)

    n_bonds = layout.bond_groups.shape[0]
    n_angles = layout.angle_groups.shape[0]
    # Arp2/3 70° branch bonds (mother↔daughter) + branch angles, appended on top
    # of the backbone bonds/angles as their own types.
    n_branch_b = int(layout.branch_bonds.shape[0])
    n_branch_a = int(layout.branch_angles.shape[0])

    snap.bonds.N = n_bonds + n_branch_b
    snap.bonds.types = (["cortex-bond"] + (["arp_branch"] if n_branch_b else []))
    bt = np.zeros(n_bonds + n_branch_b, dtype=np.uint32)
    bg = layout.bond_groups.astype(np.uint32)
    if n_branch_b:
        bg = np.concatenate([bg, layout.branch_bonds.astype(np.uint32)], axis=0)
        bt[n_bonds:] = 1                          # arp_branch typeid
    snap.bonds.typeid = bt
    snap.bonds.group = bg

    if n_angles + n_branch_a > 0:
        snap.angles.N = n_angles + n_branch_a
        snap.angles.types = (["cortex-angle"]
                             + (["arp_branch_angle"] if n_branch_a else []))
        at = np.zeros(n_angles + n_branch_a, dtype=np.uint32)
        ag = layout.angle_groups.astype(np.uint32)
        if n_branch_a:
            ag = (layout.branch_angles.astype(np.uint32) if n_angles == 0
                  else np.concatenate([ag, layout.branch_angles.astype(np.uint32)], axis=0))
            at[n_angles:] = 1 if n_angles else 0  # arp_branch_angle typeid
            if n_angles == 0:
                snap.angles.types = ["arp_branch_angle"]
        snap.angles.typeid = at
        snap.angles.group = ag

    snap.configuration.box = [p.L_box, p.L_box, p.L_box, 0.0, 0.0, 0.0]
    return snap


def build_variable_length_cortex_simulation(
    p: ResolvedH3, layout: VariableLengthCortexLayout,
    *,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
):
    """Variable-length cortex HOOMD Simulation builder.

    Mirrors :func:`build_cortex_simulation` (no xlinks; flat layout).
    """
    snap = build_variable_length_cortex_state(p, layout)
    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p.seed
    )
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    bond.params["cortex-bond"] = dict(k=p.bond_k, r0=p.rest_length)

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p.angle_k, t0=p.angle_t0)

    # Use nlist with bond exclusions (matches build_cortex_full_simulation
    # 단계 6 lesson: bonded WCA neighbours blow up when bond length
    # < r_cut).
    nlist = md.nlist.Tree(
        buffer=0.5 * p.lj_sigma, exclusions=("bond", "1-3"),
    )
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    lj.params[("actin_cortex", "actin_cortex")] = dict(
        epsilon=p.lj_epsilon, sigma=p.lj_sigma,
    )
    lj.r_cut[("actin_cortex", "actin_cortex")] = (
        p.lj_r_cut if p.lj_enabled else 0.0
    )
    lj.mode = "shift"

    ig = md.Integrator(dt=p.dt_cfl)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    sim.operations.integrator = ig

    if with_baoab:
        action, updater = make_baoab_updater(
            kT=p.kT,
            gamma={"actin_cortex": p.gamma_b},
            dt=p.dt_cfl,
            seed=p.seed,
        )
        sim.operations.updaters.append(updater)
        return sim, updater, action, layout

    return sim, None, None, layout
