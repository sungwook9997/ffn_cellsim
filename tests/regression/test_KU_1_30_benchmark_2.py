"""KU-1.30 benchmark #2 — Storm-MacKintosh strain-stiffening exponent.

Protocol
--------
Quasi-static affine shear ramp under Lees-Edwards PBC + Langevin
relaxation at every γ, with time-averaging of σ_xy over many
post-equilibration snapshots to drive the per-link thermal
fluctuation below the shear signal.

A sub-isostatic 2-D Mikado at the resolved Phase 1 defaults shows
two regimes (Broedersz & MacKintosh RMP 2014, §IV):

  - small γ (≲ 0.2): bend-dominated, soft. Non-affine modes absorb the
    shear and G(γ) is approximately flat or weakly decreasing.
  - large γ (≳ 0.2): stretch-dominated. The network enters the
    Storm-MacKintosh divergent branch
        G(γ)/G_0 = (1 − γ/γ_c)^(−p),     p ≈ 2.

The KU-1.30 #2 acceptance is the exponent ``p ∈ [1.5, 2.5]`` on the
stiffening branch (γ ∈ [γ_soft, γ_max]). G_0 here is the modulus at
the threshold of the stiffening regime, NOT the small-γ linear
modulus.

Signal/noise
------------
A single end-of-relaxation snapshot of σ_xy is dominated by per-link
thermal fluctuation; the shear-driven mean is recovered by averaging
~30 snapshots after equilibration at every γ.

Dimensional check
-----------------
[σ_xy] = N/m for the 2-D virial; [γ] dimensionless ⇒ [G_2D] = N/m;
exponent ``p`` dimensionless ⇒ band check unit-safe.

KU tags: KU-1.4 / KU-1.27 / KU-1.29 / KU-1.30 #2.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import curve_fit

from acs_kb.common.derived_params import load_config
from acs_kb.ecm.cross_links import generate_cross_links
from acs_kb.ecm.fiber_network import generate_2d_fiber_network
from acs_kb.ecm.integrator import EulerMaruyama
from acs_kb.ecm.shear_lees_edwards import (
    apply_affine_kick,
    compute_le_forces,
    compute_virial_shear_stress_2d_le,
    le_run,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[2] / "acs_kb" / "configs" / "phase1_unit1.yaml"
)

N_FIBERS = 2000
N_EQUILIBRATE_PER_GAMMA = 800
N_SAMPLE_PER_GAMMA = 30
SAMPLE_INTERVAL = 30
# The first GAMMA_SWEEP entry may sit on the soft/stiff transition bump
# (Broedersz & MacKintosh RMP 2014); the fit therefore starts from the
# argmin(G) onward so the Storm-MacKintosh divergent branch is fit on
# its own monotonically-rising data.
GAMMA_SWEEP = np.array([0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65])


def _build_network(ecm, n_fibers: int, seed: int):
    net = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=n_fibers, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=seed,
    )
    links = generate_cross_links(net, stiffness=ecm["xl_stiffness"])
    return net, links


def _sweep_sigma_time_averaged(net, links, ecm, gammas, seed=11):
    mu = ecm["stretching_modulus"]
    kappa = ecm["bending_modulus"]
    d = ecm["derived"]
    dyn = ecm["dynamics"]
    gamma_b, kT, dt = d["gamma_b"], ecm["kT"], dyn["dt"]
    rng = np.random.default_rng(seed)
    integ = EulerMaruyama()

    p = net.bead_positions.copy()
    g_curr = 0.0
    sigma_arr = np.zeros_like(gammas, dtype=np.float64)
    for i, g in enumerate(gammas):
        dg = g - g_curr
        p = apply_affine_kick(p, dg)
        g_curr = g

        def forces_fn(positions, _g=g_curr):
            return compute_le_forces(
                positions, net.rest_length, mu, kappa, net.box_size, _g,
                cross_links=links,
            )

        p = le_run(p, forces_fn, integ, gamma_b, kT, dt, net.box_size,
                   g_curr, N_EQUILIBRATE_PER_GAMMA, rng)
        samples = []
        for _ in range(N_SAMPLE_PER_GAMMA):
            p = le_run(p, forces_fn, integ, gamma_b, kT, dt, net.box_size,
                       g_curr, SAMPLE_INTERVAL, rng)
            samples.append(compute_virial_shear_stress_2d_le(
                p, net.rest_length, mu, net.box_size, g_curr, cross_links=links,
            ))
        sigma_arr[i] = float(np.mean(samples))
    return sigma_arr


def _fit_storm_mackintosh(gammas, G):
    """Fit G(γ) = G_0 · (1 − γ/γ_c)^(−p); return (p, γ_c, G_0)."""
    G0 = G[0]

    def model(gamma, p, gamma_c):
        return G0 * (1.0 - gamma / gamma_c) ** (-p)

    popt, _ = curve_fit(
        model, gammas, G,
        p0=[2.0, 1.3 * gammas[-1]],
        bounds=([0.1, gammas[-1] + 1e-6], [6.0, 5.0]),
        maxfev=10000,
    )
    return float(popt[0]), float(popt[1]), float(G0)


def test_KU_1_30_2_storm_mackintosh_exponent():
    """KU-1.30 #2: Storm-MacKintosh fit exponent p ∈ [1.5, 2.5]."""
    cfg = load_config(CONFIG_PATH)
    ecm = cfg["ecm"]
    net, links = _build_network(ecm, N_FIBERS, seed=7)
    assert len(links) > 100, f"Too few XLs ({len(links)}) for a meaningful fit."

    sigma = _sweep_sigma_time_averaged(net, links, ecm, GAMMA_SWEEP)
    G = np.abs(sigma) / GAMMA_SWEEP
    assert np.all(np.isfinite(G)) and (G > 0.0).all(), (
        f"Non-finite or zero G. σ_xy: {sigma.tolist()}, G: {G.tolist()}"
    )

    # Trim to the rising branch: start from argmin(G) so the Storm-MacKintosh
    # divergent form is fit on its own data, never on the soft/stiff
    # transition bump.
    i_min = int(np.argmin(G))
    if len(G) - i_min < 5:
        pytest.fail(
            f"Too few stiffening points (n={len(G)-i_min}) after argmin "
            f"at γ={GAMMA_SWEEP[i_min]:.3f}; extend GAMMA_SWEEP upward. "
            f"Full sweep — σ: {sigma.tolist()}, G: {G.tolist()}"
        )
    g_fit = GAMMA_SWEEP[i_min:]
    G_fit = G[i_min:]

    assert G_fit[-1] / G_fit[0] >= 1.5, (
        f"G not stiffening on the trimmed window "
        f"γ ∈ [{g_fit[0]:.3f}, {g_fit[-1]:.3f}]: "
        f"G[-1]/G[0] = {G_fit[-1]/G_fit[0]:.2f}. "
        f"Full sweep σ: {sigma.tolist()}, G: {G.tolist()}"
    )

    p, gamma_c, G0 = _fit_storm_mackintosh(g_fit, G_fit)

    assert 1.5 <= p <= 2.5, (
        f"Storm-MacKintosh exponent p={p:.3f} outside KU-1.30 #2 band "
        f"[1.5, 2.5]. Fit window γ ∈ [{g_fit[0]:.3f}, {g_fit[-1]:.3f}], "
        f"γ_c={gamma_c:.3f}, G_0={G0:.3e} N/m. "
        f"Full γ: {GAMMA_SWEEP.tolist()}, σ_xy: {sigma.tolist()}, G: {G.tolist()}"
    )
    assert gamma_c <= 2.0, (
        f"Fitted γ_c={gamma_c:.3f} unphysically large; the network is "
        f"not actually diverging in this window."
    )
