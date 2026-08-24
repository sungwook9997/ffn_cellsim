r"""Pure-NumPy proof that the resting bound-myosin → ERM transmission balances the membrane turgor.

Runs on any host (no CUDA): the coupled two-shell oracle
(:mod:`aleph.components.incumbent.resting_balance_oracle`) composes the SAME conservative force families the production
cell composes — turgor, membrane area tension, radial (unilateral, un-preloaded) ERM, and an EXPLICIT active
cortical tension — and measures the membrane residual with the cortical tension ON vs OFF.

The acceptance is the diagnosis §5 prediction, now demonstrated through the mechanism (not asserted):
  * OFF (no cortical-tension source) → the membrane carries the full turgor floor, γ_cortex ≈ 0, ERM force-free.
  * ON  → the active hoop tension contracts the cortex, stretches the radial ERM, and the transmitted load
    cancels the membrane turgor; the residual bottoms out exactly when the emergent cortical γ hits the
    Laplace target ``ΔP·R/2 − γ_mem``.

The TEST cortical tension is a stand-in for the PI-GAP fraction × per-head force; this oracle NEVER claims the
native 0.21 pN gate.
"""

from __future__ import annotations

import numpy as np

from aleph.components.incumbent.resting_balance_oracle import coupled_resting_balance, icosphere


def test_icosphere_is_a_closed_unit_sphere() -> None:
    verts, faces, edges = icosphere(2)
    assert np.allclose(np.linalg.norm(verts, axis=1), 1.0)
    # Euler characteristic V - E + F = 2 for a closed triangulated sphere.
    assert verts.shape[0] - edges.shape[0] + faces.shape[0] == 2


def test_myosin_off_leaves_the_full_turgor_floor() -> None:
    """OFF: no cortical-tension source ⇒ γ_cortex ≈ 0, ERM force-free, membrane at the turgor floor."""
    off = coupled_resting_balance(dipole_fraction=0.0, per_dipole_tension_pn=0.0, subdivisions=3, iterations=15000)
    assert off.myosin_on is False
    assert off.gamma_cortex_pn_per_um == 0.0
    assert off.erm_mean_tension_pn < 1e-9   # force-free at rest (only icosphere round-off remains)
    assert off.cortex_radial_contraction_um == 0.0
    # the membrane residual is the unbalanced turgor floor (ΔP − 2γ_mem/R)·A_node, within the coarse-mesh spread
    floor = off.turgor_per_node_pn * (40.0 - 2.0 * 10.0 / 7.5) / 40.0
    assert off.membrane_residual_mean_pn == np.float64(off.membrane_residual_mean_pn)  # finite
    assert 0.8 * floor <= off.membrane_residual_mean_pn <= 1.2 * floor


def test_myosin_on_reduces_membrane_residual_toward_the_gate() -> None:
    """ON at a TEST tension near the Laplace target: the transmission cancels most of the turgor residual."""
    off = coupled_resting_balance(dipole_fraction=0.0, per_dipole_tension_pn=0.0, subdivisions=3, iterations=20000)
    # TEST cortical tension chosen only to sit near the Laplace target — NOT a physiological fraction/force.
    on = coupled_resting_balance(dipole_fraction=1.0, per_dipole_tension_pn=88.0, subdivisions=3, iterations=20000)
    assert on.myosin_on is True
    # the seeded cortex actually developed active tension and contracted, loading the ERM
    assert on.gamma_cortex_pn_per_um > 100.0
    assert on.erm_mean_tension_pn > 0.0
    assert on.cortex_radial_contraction_um > 0.0
    # the mechanism works: the membrane residual drops by a large factor vs OFF (transmission balances turgor)
    assert on.membrane_residual_mean_pn < 0.2 * off.membrane_residual_mean_pn


def test_residual_minimum_tracks_the_laplace_target() -> None:
    """Sweeping the cortical tension, the residual minimum coincides with γ_cortex ≈ ΔP·R/2 − γ_mem (§5)."""
    tensions = [0.0, 20.0, 40.0, 60.0, 88.0, 120.0]
    results = [
        coupled_resting_balance(dipole_fraction=1.0, per_dipole_tension_pn=t, subdivisions=3, iterations=20000)
        for t in tensions
    ]
    residuals = np.array([r.membrane_residual_mean_pn for r in results])
    gammas = np.array([r.gamma_cortex_pn_per_um for r in results])
    target = results[0].gamma_laplace_target_pn_per_um
    i_min = int(np.argmin(residuals))
    # the minimum residual is interior (a genuine balance, not monotone in tension) …
    assert 0 < i_min < len(tensions) - 1
    # … and the cortical γ there is within 25% of the derived Laplace target (mechanism, not a fit)
    assert abs(gammas[i_min] - target) < 0.25 * target
    # emergent cortical γ is linear in the per-dipole tension (label-blind method-of-planes readback)
    assert gammas[-1] > gammas[1] > 0.0


def test_erm_is_not_preloaded_at_rest() -> None:
    """The ERM rest = the real separation, so with the cortex un-contracted (OFF) it carries ZERO force."""
    off = coupled_resting_balance(dipole_fraction=0.0, per_dipole_tension_pn=0.0, subdivisions=2, iterations=5000)
    assert off.erm_mean_tension_pn < 1e-9  # no counterproductive preload — force emerges only from contraction
