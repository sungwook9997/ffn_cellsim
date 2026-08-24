"""FF fiber buckling capability — the enabler for the decisive γ-amplification test (Stage 6h).

The cortical-tension active-fraction review (FF_STAGE6D §3.0c) showed the γ-floor was measured in the
DCM with buckling structurally FORBIDDEN (M-SHAKE, r/r0=1.000) — so network stress amplification
(Ronceray 2016) could never operate. This test confirms the **FF engine is buckling-CAPABLE**
(explicit bending + reshape inextensibility, NOT M-SHAKE): a fiber with excess contour length between
fixed ends bows out (Euler buckling), with the transverse amplitude growing as the end-span is
compressed, while a straight uncompressed fiber stays straight. This is what makes the decisive
buckling-enabled γ_active test runnable in FF without the DCM integrator-freeze sign-off.
"""

import numpy as np
import pytest

from aleph.laws import units as U
from aleph.laws.constraints import project_constraint_forces, reshape
from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import make_bending_force_fn


def _buckle(compress: float, n: int = 15, seg: float = 0.1, n_steps: int = 8000) -> float:
    """Pin a fiber's ends at span D=L·(1−compress) with rest contour L; relax; return transverse bow."""
    kappa = U.KAPPA_ACTIN
    L = (n - 1) * seg
    D = L * (1.0 - compress)
    x0 = np.linspace(0.0, D, n)
    pts = np.stack([x0, 1e-3 * np.sin(np.pi * x0 / D), np.zeros(n)], axis=1)   # tiny transverse seed
    net = build_fiber_network([pts], kappa=kappa)
    net.seg_rest[:] = seg                                # rest contour = uncompressed L (excess length)
    bfn = make_bending_force_fn(net)
    x = net.pos.reshape(-1, 3).copy()
    iend = [0, n - 1]
    ends = x[iend].copy()
    dt_mu = 0.05 / (kappa / seg**3)
    for it in range(n_steps):
        net.pos = x
        F = project_constraint_forces(net, bfn(x.reshape(-1)).reshape(n, 3))
        F[iend] = 0.0                                    # pin both ends
        x = x + dt_mu * F
        if (it + 1) % 50 == 0:
            net.pos = x
            xr = reshape(net, n_iter=2)
            xr[iend] = ends                              # hold ends through reshape
            x = xr
    return float(np.max(np.abs(x[:, 1])))


def test_uncompressed_fiber_stays_straight():
    assert _buckle(0.0, n_steps=4000) < 5e-3            # no excess length → no buckle


def test_buckling_amplitude_grows_with_compression():
    """Transverse bow increases with end-compression — Euler buckling in the FF engine."""
    bows = [_buckle(c) for c in (0.0, 0.08, 0.15)]
    assert bows[0] < bows[1] < bows[2]                  # monotone
    assert bows[2] > 0.1                                # clearly buckled (~0.3 µm at 15% compression)
    assert bows[0] < 0.01                               # straight baseline
