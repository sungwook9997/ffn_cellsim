"""FF relaxation to mechanical equilibrium — explicit (CFL-bound) + implicit (NF2007 Eq 2).

Two overdamped relaxers for a :class:`~ff.fiber_network.FiberNetwork`, both honouring inextensibility
via the constraint projector + periodic reshape (NF2007 §5.3):

* **explicit** — ``x ← x + (dt/γ)·P·F``. Simple, but CFL-bound: the 4th-difference bending operator
  has λ_max = 16·κ/seg³, so stability needs ``dt·μ·λ_max < 2`` ⇒ ``dt < γ·seg³/(8κ)``. For a finely
  segmented (stiff) fiber this dt is tiny and equilibration takes very many steps.

* **implicit** — the **linearly-implicit overdamped step of NF2007 Eq 2** (``[γ/dt·I + K]Δx = P·F``,
  matrix-free CG; reused from ``dcm.dcm_warp_implicit``). Unconditionally stable (NF2007 p17): dt is
  bounded by ACCURACY, not stability, so a STIFF fiber relaxes in O(1) implicit steps at huge dt —
  the paper's 10⁴×-class speed-up. Validated here: a bent stiff actin fiber reaches a few-% bending
  energy in 1–3 implicit steps at dt=10³–10⁴, where the explicit CFL would need ≫10⁴ steps.

Use explicit for the soft resting-cortex settle (not stiff — explicit is already cheap; the γ-floor
prototype uses it). Use implicit when the configuration is stiff or loaded (fine segmentation, taut
crosslinkers, the loaded shell) — same equilibrium, far fewer force evaluations.
"""

from __future__ import annotations

import numpy as np

from aleph.laws import units as U
from aleph.laws.constraints import make_projected_force_fn, project_constraint_forces, reshape
from aleph.laws.fiber_network import FiberNetwork
from aleph.laws.forces_warp import make_bending_force_fn


def _default_gamma(net: FiberNetwork) -> float:
    """Per-point overdamped drag γ [pN·s/µm] from the network's mean fiber (NF2007 §5.2)."""
    off = net.fiber_offsets
    n_per = np.diff(off)
    L_mean = float((net.seg_rest.sum() / max(net.n_fibers, 1)) if net.seg_rest.size else 1.0)
    p_mean = max(int(round(n_per.mean())) - 1, 1)
    return U.fiber_point_drag(max(L_mean, 1e-3), p_mean)


def _bending_cfl_dt(net: FiberNetwork, gamma: float, *, safety: float = 0.5) -> float:
    """Explicit-stable dt for the bending relaxation: dt < γ·seg³/(8κ) (λ_max = 16κ/seg³)."""
    seg = float(net.seg_rest.mean()) if net.seg_rest.size else 1.0
    kmax = float(net.kappa.max())
    return safety * gamma * seg**3 / (8.0 * kmax)


def relax_explicit(net: FiberNetwork, force_fn=None, *, gamma: float | None = None,
                   n_steps: int = 2000, reshape_every: int = 25, dt: float | None = None,
                   in_place: bool = True) -> np.ndarray:
    """Explicit projected overdamped relaxation + periodic reshape (CFL-bound). Returns positions."""
    if force_fn is None:
        force_fn = make_bending_force_fn(net)
    gamma = gamma if gamma is not None else _default_gamma(net)
    dt = dt if dt is not None else _bending_cfl_dt(net, gamma)
    x = net.pos.reshape(-1, 3).copy()
    for step in range(n_steps):
        net.pos = x
        F = np.asarray(force_fn(x.reshape(-1)), dtype=np.float64).reshape(net.n_nodes, 3)
        x = x + (dt / gamma) * project_constraint_forces(net, F)
        if (step + 1) % reshape_every == 0:
            net.pos = x
            x = reshape(net, n_iter=2)
    net.pos = reshape(net, n_iter=4)
    if not in_place:
        return net.pos.copy()
    return net.pos


def relax_implicit(net: FiberNetwork, force_fn=None, *, gamma: float | None = None,
                   dt: float = 1.0e3, n_steps: int = 20, n_newton: int = 2,
                   reshape_every: int = 5, cg_maxiter: int = 120, project: bool = False,
                   in_place: bool = True) -> np.ndarray:
    """Implicit overdamped relaxation (NF2007 Eq 2) + periodic reshape — unconditionally stable.

    Wraps ``dcm.dcm_warp_implicit.implicit_overdamped_step`` on the FF force. Because the step is
    unconditionally stable, ``dt`` is large (default 10³) and only a handful of steps are needed even
    for stiff fibers; reshape restores exact segment lengths periodically (inextensibility).

    ``project`` (default **False**): with it on, ``force_fn`` is wrapped with the inextensibility
    projector — but the projector's per-fiber ``(JJᵀ)⁻¹`` can go SINGULAR under the matrix-free CG's
    position perturbations (a momentarily degenerate segment), so the robust default uses the raw
    force + reshape (the reshape enforces inextensibility geometrically between steps; validated to
    reach the same equilibrium as the projected explicit relaxation to ~1e-6 µm).

    Returns the equilibrium positions.
    """
    from aleph.dcm.dcm_warp_implicit import implicit_overdamped_step

    if force_fn is None:
        force_fn = make_bending_force_fn(net)
    eval_fn = make_projected_force_fn(net, force_fn) if project else force_fn
    gamma = gamma if gamma is not None else _default_gamma(net)
    x = net.pos.reshape(-1).copy()
    for step in range(n_steps):
        x, _ = implicit_overdamped_step(x, eval_fn, gamma=gamma, dt=dt, n_newton=n_newton,
                                        cg_maxiter=cg_maxiter)
        net.pos = x.reshape(-1, 3)
        if (step + 1) % reshape_every == 0:
            x = reshape(net, n_iter=2).reshape(-1)
            net.pos = x.reshape(-1, 3)
    net.pos = reshape(net, n_iter=4)
    if not in_place:
        return net.pos.copy()
    return net.pos
