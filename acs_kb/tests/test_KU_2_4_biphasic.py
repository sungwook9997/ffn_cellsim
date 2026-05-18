"""KU-2.4 / KU-2.8 biphasic regression — replicated-seed verdict.

The Chan-Odde / Bangasser model with the KU-2.18 illustrative parameters
in this Unit is **not** strongly biphasic in the simulated stiffness
window. Multi-seed evidence (n_seeds = 5):

* per-seed argmax_E wanders over ≈ 2 decades across the sweep
  (533 Pa to 100 kPa with seed-of-the-day);
* peak prominence above the asymptote is < 1 % of the asymptote;
* ⟨F_total⟩ saturates near N_motors · F_stall ≈ 100 pN for any
  E above ~ 200 Pa.

This is physically correct: the per-clutch force at force-balance is
limited at ≈ N_m · F_stall / N_eng ≈ 4 pN, which never reaches the
slip-pathway characteristic force F_s = 30 pN. With no clutch entering
the slip regime, the system has no mechanism to dump force at very stiff
substrates → the curve plateaus, not peaks. Phase 2 work (vinculin
reinforcement of k_int^eff, FA growth) is the mechanism that raises
F_per_clutch into the slip regime and recovers a true peak.

This test consequently asserts only the **honest verdict** ("saturating"
with prominence < 10 %) and the asymptote-sanity check (⟨F⟩_plateau ∈
[0.7, 1.1] · N_motors · F_stall). The matched-stiffness analytic
``E*_matched`` is recorded in the verdict dict as info but is **not**
asserted — it predicts the peak of slip-bond systems and is not the
right reference for the catch-stabilised saturating curve we have here.
The PI is asked to re-scope the KU-2.4 / KU-2.8 peak-position acceptance
into the Unit 2.2 contract, where reinforcement enables the slip regime.

The test is marked ``slow`` so the default ``pytest`` invocation does
not run it. Use ``pytest -m slow``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from acs_kb.bridge.motor_clutch import (
    MotorClutchFA,
    MotorClutchParams,
    run_steady_state,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.types import make_focal_adhesion
from acs_kb.common.derived_params import load_bridge_config

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_unit2_1.yaml"


def _mean_traction_at_E(
    E: float,
    *,
    cfg: dict,
    n_steps: int,
    rng: np.random.Generator,
) -> float:
    sub_template = LinearElasticSubstrate.from_config(cfg)
    sub = LinearElasticSubstrate(
        young_modulus=float(E),
        poisson_ratio=sub_template.poisson_ratio,
        thickness=sub_template.thickness,
        contact_radius=sub_template.contact_radius,
    )
    mc_params = MotorClutchParams.from_config(cfg)
    fa = make_focal_adhesion(np.array([0.0, 0.0]),
                             n_clutches_total=mc_params.n_clutches)
    mc = MotorClutchFA(fa, sub, mc_params)
    dt = float(cfg["bridge"]["dynamics"]["dt"])
    result = run_steady_state(mc, dt=dt, n_steps=n_steps, rng=rng,
                              burn_in_fraction=0.5)
    return result["mean_force_total"]


def _sweep_one_seed(
    cfg: dict, E_grid: np.ndarray, seed: int, n_steps_per_E: int
) -> np.ndarray:
    out = np.empty_like(E_grid)
    for i, E in enumerate(E_grid):
        rng = np.random.default_rng(seed + i)
        out[i] = _mean_traction_at_E(float(E), cfg=cfg,
                                     n_steps=n_steps_per_E, rng=rng)
    return out


def biphasic_shape_verdict(
    E_grid: np.ndarray,
    mean_F_per_seed: np.ndarray,
    *,
    prominence_threshold: float = 0.10,
) -> dict:
    """Classify a multi-seed biphasic sweep as 'peaked' or 'saturating'.

    Verdict logic:
      - peak_prominence = (max - asymptote) / asymptote, where the
        asymptote is the mean of ⟨F⟩(E) over the top third of the E
        sweep (the would-be plateau).
      - inter_seed_argmax_spread = log10(max(argmax_E_per_seed) /
        min(argmax_E_per_seed)). A real peak gives a small spread; a
        saturating curve gives a large one (argmax is noise-dominated).
      - knee_E: smallest E at which ⟨F⟩_mean ≥ 0.9 · asymptote.
    """
    seed_argmax_E = E_grid[np.argmax(mean_F_per_seed, axis=1)]
    log10_spread = float(np.log10(seed_argmax_E.max() / seed_argmax_E.min()))

    mean_F = mean_F_per_seed.mean(axis=0)
    n = len(E_grid)
    asymptote = float(mean_F[-(n // 3):].mean())
    F_peak = float(mean_F.max())
    prominence = (F_peak - asymptote) / max(asymptote, 1e-30)

    above_knee = np.flatnonzero(mean_F >= 0.9 * asymptote)
    knee_E = float(E_grid[above_knee[0]]) if above_knee.size else float("nan")

    verdict = "peaked" if (
        prominence > prominence_threshold and log10_spread <= 1.0
    ) else "saturating"

    return dict(
        verdict=verdict,
        prominence=prominence,
        prominence_threshold=prominence_threshold,
        log10_seed_argmax_spread=log10_spread,
        seed_argmax_E=seed_argmax_E.tolist(),
        asymptote=asymptote,
        F_peak=F_peak,
        knee_E=knee_E,
        E_grid=E_grid.tolist(),
        mean_F_per_seed=mean_F_per_seed.tolist(),
        mean_F=mean_F.tolist(),
    )


@pytest.mark.slow
def test_KU_2_4_biphasic_is_saturating_not_peaked():
    """Multi-seed verdict: 'saturating' with asymptote at motor stall.

    The Phase 1 / KU-2.18 parameter set produces a curve that *rises and
    saturates* in the simulated stiffness window, not a true biphasic
    peak. Per-seed argmax_E spread of > 1 decade and peak prominence
    < 10 % of the asymptote are the operational signatures of
    saturation; the asymptote should sit near N_motors · F_stall.
    """
    cfg = load_bridge_config(CONFIG_PATH)
    b = cfg["bridge"]
    sw = b["biphasic_sweep"]
    mc_params = MotorClutchParams.from_config(cfg)

    E_grid = np.logspace(sw["log10_E_min"], sw["log10_E_max"],
                         int(sw["n_points"]))
    n_seeds = int(sw["n_seeds"])
    seed_root = int(sw["seed"])
    n_steps = int(sw["n_steps_per_E"])

    per_seed = np.stack([
        _sweep_one_seed(cfg, E_grid, seed_root + 1000 * s, n_steps)
        for s in range(n_seeds)
    ])
    v = biphasic_shape_verdict(E_grid, per_seed)

    # 1. Honest verdict: with Phase 1 / KU-2.18 parameters we expect
    #    'saturating'. If the simulation ever flips to 'peaked', surface
    #    that as a finding — do not silently re-validate.
    assert v["verdict"] == "saturating", (
        f"Unexpected verdict={v['verdict']} (prominence={v['prominence']:.3f}, "
        f"spread={v['log10_seed_argmax_spread']:.2f} dec). Investigate "
        f"before claiming a real peak. Verdict={v}"
    )

    # 2. Saturation level is near motor stall. Acceptance ±30 % covers
    #    the fact that ⟨n_engaged⟩ < N_clutches at finite k_on (some
    #    clutches are disengaged at any instant so total < N_m · F_stall).
    F_stall_total = mc_params.n_motors * mc_params.F_stall_per_motor
    ratio_to_stall = v["asymptote"] / F_stall_total
    assert 0.7 <= ratio_to_stall <= 1.1, (
        f"Asymptote ⟨F⟩={v['asymptote']*1e12:.2f} pN vs motor stall "
        f"N_m·F_stall={F_stall_total*1e12:.2f} pN (ratio {ratio_to_stall:.3f}). "
        f"Expected 0.7-1.1; verdict={v}"
    )

    # 3. Sanity: peak should not be at the boundary (would indicate the
    #    sweep is too narrow to expose saturation).
    mean_F = np.array(v["mean_F"])
    assert int(np.argmax(mean_F)) > 0, "argmax at left boundary"
