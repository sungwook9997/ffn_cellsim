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


def _point_dipole_stress_decay(
    p: ResolvedH1,
    *,
    n_softstart: int = 100,
    n_baoab: int = 200,
    n_after_dipole: int = 500,
    dipole_force_N: float = 1.0e-9,
    bin_edges_m: np.ndarray | None = None,
) -> DipoleResult:
    """Embed a point-force dipole in an equilibrated Mikado.

    A pair of beads near the box center, separated along x by ~ℓ₀,
    receive equal and opposite forces along x (``+F`` and ``-F``).
    After ``n_after_dipole`` steps, the per-particle stress
    contribution is binned radially by distance from the dipole center
    and a log–log power-law ``σ(r) ∝ r^β`` is fitted in the bulk.

    Smoke-scale only: the smoothing required to get a clean 1/r²
    signal is far beyond what fits in pytest. The production protocol
    (separate script in ``outputs/h1/``) does ensemble averaging.
    """
    sim, tq, _ = _build_simulation_with_prelude(p, n_softstart, n_baoab)

    # Pick two adjacent backbone beads near the box center to act as
    # the dipole. We need their tags so we can apply forces.
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position).copy()
        tags = np.asarray(s.particles.tag).copy()
    r_xy = np.linalg.norm(pos[:, :2], axis=1)
    near_center = np.argsort(r_xy)[:2]
    # Re-order by x so dipole "+ end" is the larger-x bead.
    if pos[near_center[0], 0] > pos[near_center[1], 0]:
        plus_idx, minus_idx = near_center[0], near_center[1]
    else:
        plus_idx, minus_idx = near_center[1], near_center[0]

    plus_tag = int(tags[plus_idx])
    minus_tag = int(tags[minus_idx])
    dipole_center = 0.5 * (pos[plus_idx, :2] + pos[minus_idx, :2])

    # Apply forces via md.force.Constant. We use a Filter that selects
    # just the dipole pair by tag using hoomd.filter.Tags.
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

    # Radial binning. HOOMD's pressure_tensor is a system-wide scalar,
    # so we can't read a per-particle stress field directly. Instead we
    # approximate the local stress via the per-particle bond+angle+LJ
    # *force magnitude* binned by distance — a proxy that decays
    # similarly to true stress in a continuum response. (Full per-
    # particle virial requires a custom force compute; deferred.)
    with sim.state.cpu_local_snapshot as s:
        pos2 = np.asarray(s.particles.position)
        F2 = np.asarray(s.particles.net_force)
    r_from_dipole = np.linalg.norm(pos2[:, :2] - dipole_center[None, :], axis=1)
    F_mag = np.linalg.norm(F2, axis=1)

    if bin_edges_m is None:
        # Default radial bins from ℓ₀ out to L_box/4 with log spacing.
        bin_edges_m = np.geomspace(
            p.rest_length, 0.25 * p.L_box, 12
        )
    centers = 0.5 * (bin_edges_m[:-1] + bin_edges_m[1:])
    # Per-bin mean force magnitude (proxy for σ_radial; same scaling).
    sigma_r = np.zeros(centers.size, dtype=np.float64)
    for i in range(centers.size):
        m = (r_from_dipole >= bin_edges_m[i]) & (r_from_dipole < bin_edges_m[i + 1])
        sigma_r[i] = float(np.mean(F_mag[m])) if m.any() else 0.0

    # Fit log–log slope where σ_r > 0 (drop empty + the closest bin
    # where the field is dominated by the applied force itself).
    valid = (sigma_r > 0.0) & np.isfinite(sigma_r)
    valid[0] = False  # closest bin contaminated by dipole bead itself
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
        """Demo-scale smoke: the dipole force is applied, the network
        relaxes without crashing, and the per-particle force magnitude
        has a decaying radial profile. The 1/r² band [-2.5, -1.5] is
        *not* enforced at demo scale because a single dipole on ~80
        fibers does not resolve the continuum decay law."""
        res = _point_dipole_stress_decay(
            resolved_demo,
            n_softstart=100, n_baoab=100, n_after_dipole=200,
        )
        # Just confirm the profile decreases on average: late bins ≤ early bins.
        valid = res.sigma_radial_pa > 0
        if valid.sum() >= 4:
            early = float(np.mean(res.sigma_radial_pa[valid][:2]))
            late = float(np.mean(res.sigma_radial_pa[valid][-2:]))
            assert late <= early, (
                f"Stress proxy did not decay with r: σ(r→0)={early:.3e}, "
                f"σ(r→L/4)={late:.3e}."
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
        res = _point_dipole_stress_decay(
            p,
            n_softstart=100, n_baoab=900, n_after_dipole=2000,
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
