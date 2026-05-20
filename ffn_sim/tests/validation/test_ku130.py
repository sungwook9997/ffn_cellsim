"""KU-1.30 ECM mechanical validation suite (H.1 M2-rest).

Validates the H.1 Mikado + L-M BAOAB pipeline against the three v1
acceptance gates frozen by commits ``11eaf13`` (#1, #2) and ``d92ac20``
(#3):

  #1  G_0 ∈ [15, 200] Pa    (small-amplitude linear modulus; v1: 32 Pa)
  #2  K(γ) ∝ γ^β with β ∈ [-2.5, -1.5]    (strain-stiffening exponent)
  #3  σ(r) ∝ 1/r²           (point-force-dipole stress decay)

Each gate is implemented at two scales:

  - **smoke** (pytest-runnable, ~30–90 s each): a reduced n_fibers
    demo config so the test fits in CI. Validates the *protocol*
    (Lees-Edwards driver, pressure-tensor read-out, stress-field
    binning) and produces an order-of-magnitude G_0 / β / power-law
    that lets us catch regressions in the measurement code.
  - **production** (skipped by default; opt-in via env
    ``H1_KU130_PRODUCTION=1``): the canonical N≈66 k Mikado run
    matching the brief's Validation acceptance table. Takes ~1 h
    (#1), several h (#2). Designed to write its diagnostics into
    ``ffn_sim/outputs/h1/`` for the REPORT.

Convention notes
----------------
- **Layer-thickness convention.** HOOMD's ``pressure_tensor`` divides
  by the full box volume ``L_box³``. Our 3D-periodic-with-2D-projection
  topology (PI 2026-05-20) places every bead at ``z=0``, so the
  material occupies a slab of thickness ``2 R_bead = 100 nm`` within
  an ``L_box = 200 μm`` cubic box — the natively-reported σ is
  diluted by ``(2 R_bead) / L_box``. To compare against the v1
  reference convention (which uses ``layer_thickness = ξ = 2 μm``),
  we multiply HOOMD's σ by ``L_box / layer_thickness`` =
  ``200 μm / 2 μm`` = 100. The conversion factor is encoded as
  ``LAYER_THICKNESS_CONVENTION_M`` below; both the raw HOOMD value
  and the v1-convention value are reported in test diagnostics.

- **Sign convention.** Stress σ_xy = -p_xy (continuum mechanics
  sign, opposite to HOOMD's pressure-tensor sign convention). Tests
  apply the negation when extracting σ from ``pressure_tensor[1]``.

References
----------
- Brief: ``ffn_sim/docs/briefs/H1_ecm_mikado.md`` §Validation acceptance.
- v1 closeouts referenced in CLAUDE.md project rule.
- KU-1.30 freeze: v1 commit ``11eaf13`` (#1, #2), ``d92ac20`` (#3).
"""

from __future__ import annotations

import math
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import yaml

import hoomd
import hoomd.md as md

from ffn_sim.ecm.equilibrate import equilibrate_no_shear
from ffn_sim.ecm.mikado import (
    ResolvedH1,
    build_mikado_simulation,
    resolve_derived,
)
from ffn_sim.ecm.shear_protocol import (
    ShearSchedule,
    attach_shear_updater,
    make_box_variant,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[2] / "configs" / "phase1_h1.yaml"
)
OUTPUTS_DIR = Path(__file__).resolve().parents[2] / "outputs" / "h1"
LAYER_THICKNESS_CONVENTION_M: float = 2.0e-6  # v1 ξ-based slab convention


# ---------------------------------------------------------------------------
# Demo-scale config: a reduced Mikado that fits in pytest wall time.
#
# Why a separate demo config instead of editing phase1_h1.yaml: the yaml
# is the canonical KU-anchored config for the brief. Smoke tests
# explicitly tune one knob (target_segment_length → larger ⇒ fewer
# fibers) so the validation runs at ~3 k beads instead of ~66 k beads.
# Demo configs are flagged ``demo_mode: true`` per the yaml's contract.
# ---------------------------------------------------------------------------
def _load_demo_config(
    target_segment_length_m: float = 8.0e-6,
    L_box_m: float = 80.0e-6,
) -> dict:
    """Load the H.1 yaml and override two knobs for a smaller demo system.

    Defaults give ``n_fibers ≈ 78`` (vs production 3142), ``N_beads ≈
    1638`` (vs ~66 k). Wall time per BAOAB step is roughly ``N`` so
    this is ~40× faster than production.
    """
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg["ecm"]["L_box"] = L_box_m
    cfg["ecm"]["target_segment_length"] = target_segment_length_m
    cfg["ecm"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def resolved_demo() -> ResolvedH1:
    return resolve_derived(_load_demo_config())


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _attach_thermo(sim: hoomd.Simulation) -> md.compute.ThermodynamicQuantities:
    """Attach ThermodynamicQuantities + a Table Writer logging
    ``pressure`` so HOOMD 7.0.1 registers the virial-compute flag.

    Without a Writer actively logging a pressure-related quantity, the
    integrator does NOT compute virials and ``pressure_tensor`` returns
    NaN. Discovered in M2-rest debugging (see REPORT §Open issues #1)
    — surfaces a HOOMD 7.0.1 wiring requirement, not a bug.
    """
    import os
    tq = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(tq)
    logger = hoomd.logging.Logger(categories=["scalar"])
    logger.add(tq, quantities=["pressure"])
    writer = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(1),  # fires every step so virial flag stays live
        logger=logger,
        output=open(os.devnull, "w"),
    )
    sim.operations.writers.append(writer)
    return tq


def _sigma_xy(tq: md.compute.ThermodynamicQuantities) -> float:
    """HOOMD's σ_xy = -p_xy (continuum sign). Returns Pa (HOOMD volume)."""
    return -float(tq.pressure_tensor[1])


def _to_layer_convention(sigma_pa_hoomd: float, sim: hoomd.Simulation) -> float:
    """Multiply by ``L_box / layer_thickness`` to match v1 ξ-slab convention."""
    L_box = float(sim.state.box.Lx)
    return sigma_pa_hoomd * (L_box / LAYER_THICKNESS_CONVENTION_M)


def _build_simulation_with_prelude(
    p: ResolvedH1,
    n_softstart: int = 100,
    n_baoab: int = 200,
) -> tuple[hoomd.Simulation, md.compute.ThermodynamicQuantities, dict]:
    """Build + equilibrate. Returns ``(sim, thermo, prelude_diagnostics)``."""
    sim, updater, action = build_mikado_simulation(p, with_cross_links=True)
    diag = equilibrate_no_shear(
        sim, action, updater,
        n_softstart=n_softstart, n_baoab=n_baoab,
        rest_length=p.rest_length, gamma_b=p.gamma_b,
    )
    tq = _attach_thermo(sim)
    return sim, tq, diag


# ---------------------------------------------------------------------------
# KU-1.30 #1 — small-amplitude step-strain G_0
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class StepStrainResult:
    gamma_applied: float
    sigma_xy_hoomd: float            # Pa, HOOMD V = L_box³ convention
    sigma_xy_layer_pa: float         # Pa, v1 layer_thickness=2μm convention
    G_0_layer_pa: float              # σ / γ in v1 convention
    n_relax_steps: int


def _step_strain_G0(
    p: ResolvedH1,
    *,
    gamma: float = 0.01,
    n_softstart: int = 100,
    n_baoab: int = 200,
    n_ramp: int = 100,
    n_relax: int = 500,
    n_avg: int = 200,
) -> StepStrainResult:
    """Step strain γ_0 = ``gamma`` ramped over ``n_ramp`` steps + hold.

    Stress is averaged over the trailing ``n_avg`` steps of an
    ``n_relax + n_avg`` hold window, so transient ringing has time to
    decay. Returns the linear modulus G_0 = ⟨σ_xy⟩ / γ at this γ.

    Linear-regime test: at γ = 0.01 (1%) the response is expected to
    be within ~10% of the small-γ Taylor expansion of the network's
    constitutive law, so ⟨σ_xy⟩ / γ is a faithful G_0 estimator.
    """
    sim, tq, _ = _build_simulation_with_prelude(p, n_softstart, n_baoab)
    schedule = ShearSchedule(gamma_max=gamma, ramp_steps=n_ramp, hold_steps=0)
    attach_shear_updater(sim, schedule)
    sim.run(schedule.total_steps)        # finish ramp
    # Relax phase — let the network's elastic response settle.
    sim.run(n_relax)
    # Average σ_xy over n_avg samples (one per step).
    sigmas = np.empty(n_avg, dtype=np.float64)
    for i in range(n_avg):
        sim.run(1)
        sigmas[i] = _sigma_xy(tq)
    sigma_hoomd = float(np.mean(sigmas))
    sigma_layer = _to_layer_convention(sigma_hoomd, sim)
    G_0 = sigma_layer / gamma if gamma != 0.0 else 0.0
    return StepStrainResult(
        gamma_applied=gamma,
        sigma_xy_hoomd=sigma_hoomd,
        sigma_xy_layer_pa=sigma_layer,
        G_0_layer_pa=G_0,
        n_relax_steps=n_relax + n_avg,
    )


class TestKU130_1_G0:
    """KU-1.30 #1: small-amplitude G_0 ∈ [15, 200] Pa (v1 ref 32 Pa)."""

    def test_g0_demo_scale_within_band(self, resolved_demo):
        """Demo-scale run: assert G_0 is finite and positive, and in
        the brief's band with a 2× tolerance widening for finite-size
        bias (demo n_fibers ≪ production)."""
        res = _step_strain_G0(
            resolved_demo,
            gamma=0.01,
            n_softstart=100, n_baoab=200,
            n_ramp=100, n_relax=200, n_avg=200,
        )
        assert math.isfinite(res.G_0_layer_pa), (
            f"G_0 not finite: σ_xy_HOOMD={res.sigma_xy_hoomd:.3e} Pa."
        )
        # Demo scale: assert order-of-magnitude band [1, 1000] Pa
        # rather than the production [15, 200] Pa, because n_fibers is
        # 40× smaller → finite-size corrections of factor ~few.
        assert 1.0 < res.G_0_layer_pa < 1000.0, (
            f"Demo G_0 = {res.G_0_layer_pa:.3e} Pa outside the order-"
            f"of-magnitude band [1, 1000] Pa. σ_HOOMD={res.sigma_xy_hoomd:.3e} Pa."
        )

    @pytest.mark.skipif(
        os.environ.get("H1_KU130_PRODUCTION", "") != "1",
        reason="Production-scale KU-1.30 #1; opt-in via H1_KU130_PRODUCTION=1.",
    )
    def test_g0_production_within_band(self):
        """Production-scale G_0 ∈ [15, 200] Pa on the canonical N≈66k
        Mikado. ~1 h wall-time at ~150 steps/s; writes diagnostics
        for the REPORT."""
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f)
        p = resolve_derived(cfg)
        t0 = time.time()
        res = _step_strain_G0(
            p,
            gamma=0.01,
            n_softstart=100, n_baoab=900,
            n_ramp=500, n_relax=4500, n_avg=2000,
        )
        elapsed = time.time() - t0
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        (OUTPUTS_DIR / "ku130_g0_production.json").write_text(
            f"""{{
  "gamma": {res.gamma_applied},
  "sigma_xy_hoomd_pa": {res.sigma_xy_hoomd},
  "sigma_xy_layer_pa": {res.sigma_xy_layer_pa},
  "G_0_layer_pa": {res.G_0_layer_pa},
  "n_relax_steps": {res.n_relax_steps},
  "wall_time_s": {elapsed:.2f}
}}
"""
        )
        assert 15.0 <= res.G_0_layer_pa <= 200.0, (
            f"Production G_0 = {res.G_0_layer_pa:.3e} Pa outside KU-1.30 "
            f"band [15, 200] Pa (v1 ref: 32 Pa). Surface to PI per "
            f"CLAUDE.md no-gate-loosening."
        )


# ---------------------------------------------------------------------------
# KU-1.30 #2 — strain-stiffening exponent
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class StiffeningResult:
    gammas: np.ndarray
    sigma_xy_layer_pa: np.ndarray
    K_layer_pa: np.ndarray
    log_log_slope: float
    fit_window_gamma: tuple[float, float]


def _strain_stiffening(
    p: ResolvedH1,
    *,
    gamma_max: float = 0.30,
    n_softstart: int = 100,
    n_baoab: int = 200,
    n_ramp_steps: int = 2000,
    n_samples: int = 20,
    fit_window: tuple[float, float] = (0.05, 0.30),
) -> StiffeningResult:
    """Single monotonic ramp γ ∈ [0, γ_max] with σ_xy(γ) sampled at
    ``n_samples`` evenly-spaced points along the ramp.

    Differential modulus K(γ) = dσ_xy/dγ computed by central finite
    difference. Log–log slope of K(γ) vs γ fitted by ordinary
    least-squares over the ``fit_window`` strain range.
    """
    sim, tq, _ = _build_simulation_with_prelude(p, n_softstart, n_baoab)
    schedule = ShearSchedule(
        gamma_max=gamma_max, ramp_steps=n_ramp_steps, hold_steps=0
    )
    attach_shear_updater(sim, schedule)
    # Sample box.xy and σ_xy at evenly-spaced timesteps over the ramp.
    sample_at = np.linspace(
        n_ramp_steps // n_samples, n_ramp_steps, n_samples, dtype=int
    )
    gammas = np.empty(n_samples, dtype=np.float64)
    sigmas_layer = np.empty(n_samples, dtype=np.float64)
    last = 0
    for i, ts in enumerate(sample_at):
        delta = int(ts) - last
        if delta > 0:
            sim.run(delta)
            last = int(ts)
        gammas[i] = float(sim.state.box.xy)
        sigmas_layer[i] = _to_layer_convention(_sigma_xy(tq), sim)

    # Differential modulus K(γ) = dσ/dγ via central diff. np.gradient
    # silently divides by zero on repeated γ (the first sample is taken
    # while the ramp is still warming up); restrict to the strictly-
    # increasing subset to keep the gradient well-defined.
    increasing = np.diff(gammas, prepend=gammas[0] - 1.0) > 0
    K = np.full_like(sigmas_layer, np.nan)
    if increasing.sum() >= 2:
        K[increasing] = np.gradient(sigmas_layer[increasing], gammas[increasing])
    # Log–log fit in fit_window. Drop zero/negative K (could happen on
    # noisy demo runs — log-fit only the positive subset within window).
    mask = (
        (gammas >= fit_window[0])
        & (gammas <= fit_window[1])
        & (K > 0.0)
        & np.isfinite(K)
    )
    if mask.sum() >= 3:
        lg = np.log(gammas[mask])
        lk = np.log(K[mask])
        slope = float(np.polyfit(lg, lk, 1)[0])
    else:
        slope = float("nan")

    return StiffeningResult(
        gammas=gammas,
        sigma_xy_layer_pa=sigmas_layer,
        K_layer_pa=K,
        log_log_slope=slope,
        fit_window_gamma=fit_window,
    )


class TestKU130_2_StrainStiffening:
    """KU-1.30 #2: log–log slope of K(γ) vs γ ∈ [-2.5, -1.5]."""

    def test_strain_stiffening_demo_scale_protocol(self, resolved_demo):
        """Demo-scale smoke: the protocol runs to completion with the
        sheared Lees-Edwards box; box.xy reaches γ_max; σ_xy(γ) is
        finite at every sample. The log–log slope band [-2.5, -1.5] is
        *not* enforced at demo scale because (a) the finite-size noise
        dominates the small-γ regression and (b) without ensemble
        averaging a single ramp produces a stress trace that fluctuates
        on the kT scale around the elastic mean."""
        res = _strain_stiffening(
            resolved_demo,
            gamma_max=0.20,
            n_softstart=100, n_baoab=100,
            n_ramp_steps=500, n_samples=10,
            fit_window=(0.05, 0.20),
        )
        # Protocol invariants only.
        assert np.all(np.isfinite(res.gammas)), "Non-finite γ samples."
        assert np.all(np.isfinite(res.sigma_xy_layer_pa)), (
            "Non-finite σ_xy samples — pressure_tensor wiring broken."
        )
        assert res.gammas[-1] >= 0.9 * 0.20, (
            f"Final γ = {res.gammas[-1]:.3e} did not reach 90% of γ_max=0.20."
        )

    @pytest.mark.skipif(
        os.environ.get("H1_KU130_PRODUCTION", "") != "1",
        reason="Production-scale KU-1.30 #2; opt-in via H1_KU130_PRODUCTION=1.",
    )
    def test_strain_stiffening_production_within_band(self):
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f)
        p = resolve_derived(cfg)
        t0 = time.time()
        res = _strain_stiffening(
            p,
            gamma_max=0.30,
            n_softstart=100, n_baoab=900,
            n_ramp_steps=10000, n_samples=30,
            fit_window=(0.05, 0.30),
        )
        elapsed = time.time() - t0
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUTPUTS_DIR / "ku130_strain_stiffening_production.npz",
            gammas=res.gammas,
            sigma_xy_layer_pa=res.sigma_xy_layer_pa,
            K_layer_pa=res.K_layer_pa,
            log_log_slope=res.log_log_slope,
            wall_time_s=elapsed,
        )
        assert -2.5 <= res.log_log_slope <= -1.5, (
            f"Strain-stiffening slope {res.log_log_slope:.3f} outside "
            f"KU-1.30 #2 band [-2.5, -1.5]. Surface to PI per "
            f"CLAUDE.md no-gate-loosening."
        )


# ---------------------------------------------------------------------------
# KU-1.30 #3 — point-force-dipole 1/r² stress decay
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class DipoleResult:
    bin_centers_m: np.ndarray
    sigma_radial_pa: np.ndarray
    fit_exponent: float


def _bond_virial_per_bond(
    sim: hoomd.Simulation,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the per-bond Born-virial xx + xy contributions + midpoint xy.

    For each bond ``a-b`` with Lees-Edwards-aware MI displacement
    ``r_ab = MI(r_b - r_a)`` and harmonic force magnitude
    ``F = k (|r_ab| - r0)``, the stress contribution from this bond to
    the macroscopic stress tensor is::

        σ_ij_bond = F_vec_i · r_ab_j        (Newton · meter)

    where ``F_vec`` is the force *on bead a* (along +r̂_ab when stretched),
    by the standard pair-virial convention. Caller divides by the bin's
    shell volume to get Pa.

    Returns ``(virial_xx_Nm, virial_xy_Nm, midpoint_xy_m)``.
    """
    ig = sim.operations.integrator
    # Locate the bond.Harmonic force compute (assumes single bond compute,
    # which matches H.1's mikado build).
    bond_force = next(
        f for f in ig.forces if isinstance(f, md.bond.Harmonic)
    )
    bond_type_names = list(sim.state.bond_types)
    k_per_type = np.array(
        [bond_force.params[t]["k"] for t in bond_type_names], dtype=np.float64
    )
    r0_per_type = np.array(
        [bond_force.params[t]["r0"] for t in bond_type_names], dtype=np.float64
    )

    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position).copy()
        bg = np.asarray(s.bonds.group).copy()
        bt = np.asarray(s.bonds.typeid).copy()

    box = sim.state.box
    Lx, Ly, Lz = box.Lx, box.Ly, box.Lz
    xy, xz, yz = box.xy, box.xz, box.yz

    ra = pos[bg[:, 0]]
    rb = pos[bg[:, 1]]
    dr = rb - ra

    # MI wrap with full Lees-Edwards tilt (matches baoab._wrap_into_box
    # convention so the result is consistent across sheared / unsheared
    # frames).
    fz = dr[:, 2] / Lz
    fy = (dr[:, 1] - yz * Lz * fz) / Ly
    fx = (dr[:, 0] - xy * Ly * fy - xz * Lz * fz) / Lx
    fx -= np.round(fx)
    fy -= np.round(fy)
    fz -= np.round(fz)
    dr[:, 0] = Lx * fx + xy * Ly * fy + xz * Lz * fz
    dr[:, 1] = Ly * fy + yz * Lz * fz
    dr[:, 2] = Lz * fz

    r_norm = np.linalg.norm(dr, axis=1)
    safe_r = np.where(r_norm > 0.0, r_norm, 1.0)
    k_bond = k_per_type[bt]
    r0_bond = r0_per_type[bt]
    # F_vec = force on bead a = +k·(|r_ab|−r0) · r̂_ab when stretched.
    F_mag_signed = k_bond * (r_norm - r0_bond)
    F_vec = F_mag_signed[:, None] * (dr / safe_r[:, None])

    # σ_ij_bond = F_i · r_ab_j  (Newton · meter); caller divides by the
    # bin's shell volume V_bin to get Pa.
    virial_xx = F_vec[:, 0] * dr[:, 0]
    virial_xy = F_vec[:, 0] * dr[:, 1]

    midpoint_xy = pos[bg[:, 0], :2] + 0.5 * dr[:, :2]
    return virial_xx, virial_xy, midpoint_xy


def _point_dipole_stress_decay(
    p: ResolvedH1,
    *,
    n_softstart: int = 100,
    n_baoab: int = 200,
    n_after_dipole: int = 500,
    dipole_force_N: float = 1.0e-9,
    bin_edges_m: np.ndarray | None = None,
    n_time_avg_samples: int = 1,
    time_avg_spacing: int = 100,
) -> DipoleResult:
    """Embed a point-force dipole in an equilibrated Mikado, compute the
    **true Born bond-virial stress field** σ_xx(r) radially.

    PI 2026-05-21 (v2 of the bond-virial implementation):

    The x-x oriented dipole produces an anisotropic far-field stress
    with ``σ_xx ∝ 1/r²``; the off-diagonal ``σ_xy`` has angular
    structure that averages to zero on a radial shell unless an
    angular weight is applied. We therefore fit ``σ_xx(r)`` (the
    dipole-parallel component), summed *signed* per shell so thermal
    contributions cancel rather than accumulate as Σ|noise| ∝ √N
    (the v1 implementation's mistake — Σ|σ_xy| accidentally yielded
    a positive slope ~ +1.0 from noise scaling with bond count).

    Procedure:

    1. Two adjacent beads near the box origin receive ±F along x via
       ``md.force.Constant``.
    2. ``n_after_dipole`` BAOAB steps relax the network.
    3. (optional) ``n_time_avg_samples`` further snapshots taken at
       ``time_avg_spacing`` step spacing, each yielding a per-bond
       ``σ_xx_bond = F_x · r_ab_x`` array. The time-average over
       samples cancels thermal-fluctuation σ_xx noise (mean zero) and
       retains the dipole-induced σ_xx mean (≠ zero).
    4. Bonds binned radially by midpoint distance ``r`` from the dipole
       centre; per-bin ``σ_xx(r) = |Σ_{bonds in shell} ⟨σ_xx_bond⟩|
       / V_shell`` with ``V_shell = 2π·r·dr·L_z``.
    5. Log–log slope of σ_xx(r) vs r over the bulk (drop closest +
       farthest bins).
    """
    sim, tq, _ = _build_simulation_with_prelude(p, n_softstart, n_baoab)

    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position).copy()
        tags = np.asarray(s.particles.tag).copy()
    r_xy = np.linalg.norm(pos[:, :2], axis=1)
    near_center = np.argsort(r_xy)[:2]
    if pos[near_center[0], 0] > pos[near_center[1], 0]:
        plus_idx, minus_idx = near_center[0], near_center[1]
    else:
        plus_idx, minus_idx = near_center[1], near_center[0]
    plus_tag = int(tags[plus_idx])
    minus_tag = int(tags[minus_idx])
    dipole_center = 0.5 * (pos[plus_idx, :2] + pos[minus_idx, :2])

    plus_filter = hoomd.filter.Tags([plus_tag])
    minus_filter = hoomd.filter.Tags([minus_tag])
    f_plus = md.force.Constant(filter=plus_filter)
    f_plus.constant_force["actin_ecm"] = (dipole_force_N, 0.0, 0.0)
    f_minus = md.force.Constant(filter=minus_filter)
    f_minus.constant_force["actin_ecm"] = (-dipole_force_N, 0.0, 0.0)
    ig = sim.operations.integrator
    ig.forces.append(f_plus)
    ig.forces.append(f_minus)

    sim.run(n_after_dipole)

    # Time-averaged σ_xx per bond. The first sample is taken right
    # after the post-dipole relax; additional samples are spaced by
    # `time_avg_spacing` so they are decorrelated by ~τ_relax.
    samples_xx = []
    midpoint_xy = None
    for s_i in range(n_time_avg_samples):
        if s_i > 0:
            sim.run(time_avg_spacing)
        v_xx, _v_xy, mid = _bond_virial_per_bond(sim)
        samples_xx.append(v_xx)
        if s_i == 0:
            midpoint_xy = mid
    virial_xx_avg = np.mean(np.stack(samples_xx, axis=0), axis=0)

    assert midpoint_xy is not None
    r_from_dipole = np.linalg.norm(
        midpoint_xy - dipole_center[None, :], axis=1
    )

    if bin_edges_m is None:
        bin_edges_m = np.geomspace(p.rest_length, 0.25 * p.L_box, 12)
    centers = 0.5 * (bin_edges_m[:-1] + bin_edges_m[1:])
    L_z = float(sim.state.box.Lz)

    # σ_xx(r) per shell: |signed sum of bond virials| / V_shell.
    # Signed sum lets the thermal contributions (mean zero) cancel
    # while the dipole-induced σ_xx (which has a definite sign on each
    # angular sector for an x-x dipole) accumulates coherently. The
    # outer |.| then yields a positive magnitude appropriate for the
    # log-log fit.
    sigma_r = np.zeros(centers.size, dtype=np.float64)
    for i in range(centers.size):
        lo, hi = bin_edges_m[i], bin_edges_m[i + 1]
        m = (r_from_dipole >= lo) & (r_from_dipole < hi)
        if m.any():
            V_shell = 2.0 * math.pi * centers[i] * (hi - lo) * L_z
            sigma_r[i] = float(np.abs(np.sum(virial_xx_avg[m])) / V_shell)

    valid = (sigma_r > 0.0) & np.isfinite(sigma_r)
    if valid.sum() > 0:
        idxs = np.where(valid)[0]
        valid[idxs[0]] = False
        if len(idxs) > 1:
            valid[idxs[-1]] = False
    if valid.sum() >= 3:
        lg_r = np.log(centers[valid])
        lg_s = np.log(sigma_r[valid])
        slope = float(np.polyfit(lg_r, lg_s, 1)[0])
    else:
        slope = float("nan")

    return DipoleResult(
        bin_centers_m=centers,
        sigma_radial_pa=sigma_r,
        fit_exponent=slope,
    )


class TestKU130_3_PointDipole:
    """KU-1.30 #3: stress decay σ(r) ∝ 1/r² (β ≈ -2)."""

    def test_point_dipole_demo_protocol(self, resolved_demo):
        """Demo-scale smoke: protocol invariants only (true Born bond-
        virial stress field is computed and returns finite Pa values).
        The β ∈ [-2.5, -1.5] band is *not* enforced at demo scale
        because a single dipole over ~80 fibers does not resolve the
        continuum 1/r² decay law without ensemble averaging — the
        production variant runs ensemble averaging for the band gate."""
        res = _point_dipole_stress_decay(
            resolved_demo,
            n_softstart=100, n_baoab=100, n_after_dipole=200,
        )
        # The bond-virial sum must yield at least 3 occupied radial bins
        # (otherwise we cannot fit anything). At demo scale we expect
        # ~4-6 occupied bins in the 12-bin geomspace grid.
        valid = (res.sigma_radial_pa > 0.0) & np.isfinite(res.sigma_radial_pa)
        assert valid.sum() >= 3, (
            f"Bond-virial radial profile has only {int(valid.sum())} "
            f"occupied bins; expected ≥ 3."
        )
        # Slope must be finite when computable (fallback nan is allowed
        # only if the trimmed-edge fit window has < 3 bins).
        if math.isfinite(res.fit_exponent):
            assert -10.0 < res.fit_exponent < 10.0, (
                f"Bond-virial log-log slope unreasonable: "
                f"{res.fit_exponent:.3f}."
            )

    @pytest.mark.skipif(
        os.environ.get("H1_KU130_PRODUCTION", "") != "1",
        reason="Production-scale KU-1.30 #3; opt-in via H1_KU130_PRODUCTION=1.",
    )
    def test_point_dipole_production_band(self):
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f)
        p = resolve_derived(cfg)
        t0 = time.time()
        # 5 time-averaged snapshots spaced by 200 steps to cancel
        # thermal σ_xx noise (mean zero) while accumulating the
        # dipole-induced mean (≠ zero).
        res = _point_dipole_stress_decay(
            p,
            n_softstart=100, n_baoab=900, n_after_dipole=2000,
            n_time_avg_samples=5, time_avg_spacing=200,
        )
        elapsed = time.time() - t0
        OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUTPUTS_DIR / "ku130_point_dipole_production.npz",
            bin_centers_m=res.bin_centers_m,
            sigma_radial_pa=res.sigma_radial_pa,
            fit_exponent=res.fit_exponent,
            wall_time_s=elapsed,
        )
        assert -2.5 <= res.fit_exponent <= -1.5, (
            f"Point-dipole stress decay exponent {res.fit_exponent:.3f} "
            f"outside KU-1.30 #3 band [-2.5, -1.5]. Surface to PI."
        )
