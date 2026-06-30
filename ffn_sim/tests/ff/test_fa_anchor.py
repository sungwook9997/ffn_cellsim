"""FA / traction integrin-clutch anchor layer (ff/fa_anchor) — increment-3 validation.

The FA clutch is an integrin α5β1–fibronectin CATCH-SLIP bond (Kong 2009 / KB-2.5), NOT cadherin;
de-adhesion is emergent from load-dependent rupture, not a binary latch. Validates the catch-slip
lifetime signature, the traction readout, and the turnover steady state on a stress fiber.
"""

import numpy as np

from ffn_sim.ff.architecture_spec import STRESS_FIBER
from ffn_sim.ff.fa_anchor import (
    attach_fa_clutches,
    clutch_forces,
    fa_end_nodes,
    step_clutch_turnover,
    traction_force_nN,
)
from ffn_sim.ff.hand_kmc import INTEGRIN_A5B1
from ffn_sim.ff.weave import weave


def test_integrin_catch_slip_lifetime_has_interior_peak():
    """The integrin α5β1 bond is a CATCH-slip bond: the lifetime τ(F)=1/k_off RISES then FALLS (an
    interior peak F*>0), unlike a pure slip bond. (Recorded KB-2.5 params give F*≈7 pN; the catch-bond
    behaviour is the assertion, not the exact peak — the KB's '30 pN' claim is flagged inconsistent.)"""
    F = np.linspace(0.0, 80.0, 400)
    tau = 1.0 / INTEGRIN_A5B1.off_rate(F)
    fstar = F[int(np.argmax(tau))]
    assert fstar > 2.0                                       # interior peak = catch present (not slip-only)
    assert tau[0] < tau.max() and tau[-1] < tau.max()       # rises from F=0 then falls (catch→slip)


def test_fa_clutch_traction_and_turnover():
    """Attach integrin clutches to a stress-fiber's FA end nodes: force-free at rest (0 traction),
    positive traction under load, and the bound fraction reaches the catch-slip+Poisson steady state."""
    sf = weave(STRESS_FIBER, rng=np.random.default_rng(0))
    pop = attach_fa_clutches(sf, rng=np.random.default_rng(0))
    assert len(pop.actin_nodes) == len(fa_end_nodes(sf)) and pop.bound.all()
    assert traction_force_nN(pop, sf.net.pos) < 1e-6        # force-free at the binding geometry

    loaded = sf.net.pos.copy(); loaded[pop.actin_nodes, 2] *= 1.001
    assert traction_force_nN(pop, loaded) > 0.0             # load develops clutch tension
    assert np.all(clutch_forces(pop, loaded) >= 0.0)

    rng = np.random.default_rng(1)
    for _ in range(400):
        step_clutch_turnover(pop, sf.net.pos, tau=0.01, rng=rng)
    bf = pop.bound.mean()
    koff0 = INTEGRIN_A5B1.off_rate(0.0)                      # force-free off-rate
    steady = INTEGRIN_A5B1.k_on / (INTEGRIN_A5B1.k_on + koff0)
    assert abs(bf - steady) < 0.12                          # matches the analytic steady bound fraction
