r"""Power-stroke -> force coupling — the analytic ground truth for the ACTIVE directed contraction (P2 root).

Host-side acceptance oracle (pure NumPy, NO Warp) for the ONE piece the earlier oracle suite could not gate:
that the myosin power stroke (the walked ``abscissa``) actually ENTERS the crossbridge force, so a stepping
head generates a DIRECTED contractile force and force-velocity self-limits at the stall force. It mirrors, bit
-for-formula, the device kernels :func:`~aleph.components.motor.minifilament_warp.crossbridge_kernel` /
:func:`~aleph.components.motor.minifilament_warp.compute_head_loads_kernel` (which the dev Mac cannot launch), so the
coupling is a LOCAL gate rather than something only the gbook native gate can find.

THE MECHANISM (Huxley/AFINES/Cytosim cross-bridge with a walking attachment point). A bound head grabbed an
actin node at ``anchor`` and has since walked a displacement ``s = abscissa`` (>= 0) toward that filament's
barbed end ``walk_dir`` (unit). Its zero-strain attachment reference therefore ADVANCES along actin::

    x_att = anchor + s * walk_dir                                            (the power stroke)

The crossbridge is a compliance spring (stiffness ``k_xb``, rest ``r0_xb`` ~ 0) between the head particle and
``x_att``::

    F_on_head  = k_xb (|x_att - head| - r0_xb) * (x_att - head)/|x_att - head|      (pulls head toward x_att)
    F_on_actin = -F_on_head                                                          (Newton's 3rd law)

As the head steps (``s`` grows) ``x_att`` moves away from the head, the spring stretches, and F_on_actin points
along ``-walk_dir`` — the actin is dragged the way OPPOSITE the walk: the directed contraction. With ``s = 0``
(or ``walk_dir = 0``) the active term vanishes and only the passive spring toward the grabbed node remains.

FORCE-VELOCITY EMERGES (not imposed). The load the Hill law sees is the TANGENTIAL (along-walk) tension::

    load = k_xb ( dot(x_att - head, walk_dir) - r0_xb )   ->  d(load)/d(abscissa) = k_xb

so the walked strain lifts the load. The abscissa advances by ``ds/dt = hill_velocity(load)``; in an isometric
condition (head + anchor held) it grows until ``load = F_stall`` where ``hill_velocity = 0``. Hence the settled
crossbridge force -> ``F_stall`` and the settled strain -> ``F_stall/k_xb + r0_xb`` (the physical working stroke,
5-20 nm — the same quantity ``minifilament_topology.working_stroke_strain`` certifies). Nothing is tuned: the
stall force is where the Hill velocity vanishes.

Sanity Gate (of the oracle itself — tests/ac/motor/test_powerstroke_coupling_oracle.py):
  * directed sign: a walking head (s > 0) drags actin along ``-walk_dir`` (contractile) with |F| = k_xb(s-r0_xb);
  * Newton's 3rd law: F_on_head = -F_on_actin exactly;
  * abscissa enters the force: d(load)/d(abscissa) = k_xb; s = 0 => zero active force (decoupling arbiter);
  * force-velocity / stall: isometric stepping settles at load -> F_stall, v -> 0, strain -> F_stall/k_xb + r0_xb;
  * monotone approach: the load rises monotonically toward F_stall and never overshoots it (clamped Hill).

⚠ ``walk_dir`` is the bound actin's barbed-end polarity — a GEOMETRY input, not an I0-B3 magnitude. The SOURCE
builder (``build_minifilament_nodes``) derives a bipolar-geometry default (outward along the backbone axis); the
production engine overwrites it per head from the true actin polarity (I4-weave hand-off). ``r0_xb ~ 0`` is the
physical "a bound head sits on its actin site" zero-strain length (grounded, not tuned); ``k_xb`` / ``F_stall``
/ ``v0`` / ``kappa`` remain the I0-B3 GAPs — this oracle SWEEPS them, it never chooses a value.

References:
  A.F. Huxley 1957 (cross-bridge sliding filament); Hill 1938 (force-velocity); Freedman 2017 (AFINES walking
  motor + spring); Nedelec & Foethke 2007 §10 (Cytosim Hand abscissa/step). Stam-Hocky 2015 (minifilament).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.components.motor.hill_fv_analytic import hill_velocity

__all__ = [
    "attachment_point",
    "crossbridge_force",
    "stall_abscissa",
    "step_engaged_head",
]


def attachment_point(
    anchor_pos: npt.ArrayLike,
    abscissa: float,
    walk_dir: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Power-stroke-advanced crossbridge attachment point ``x_att = anchor + abscissa * walk_dir`` [µm].

    Mirrors the device ``@wp.func`` :func:`aleph.components.motor.minifilament_warp.attachment_point`.
    """
    return np.asarray(anchor_pos, dtype=np.float64) + float(abscissa) * np.asarray(walk_dir, dtype=np.float64)


def crossbridge_force(
    head_pos: npt.ArrayLike,
    anchor_pos: npt.ArrayLike,
    abscissa: float,
    walk_dir: npt.ArrayLike,
    k_xb: float,
    r0_xb: float = 0.0,
) -> dict[str, object]:
    """One bound head's power-stroke crossbridge force + tangential load (host mirror of the Warp kernels).

    Args:
        head_pos: head particle position ``(3,)`` [µm].
        anchor_pos: grabbed actin node position ``(3,)`` [µm].
        abscissa: walked displacement ``s`` along actin [µm] (>= 0).
        walk_dir: unit barbed-end/walk direction ``(3,)``.
        k_xb: crossbridge stiffness [pN/µm].
        r0_xb: crossbridge rest length [µm] (~0, the head-sits-on-actin zero-strain length).

    Returns:
        dict with ``x_att`` (advanced attachment point), ``f_on_head`` / ``f_on_actin`` (``(3,)`` [pN],
        equal-and-opposite), ``load`` (tangential resisting tension [pN] fed to Hill), ``newton3rd_residual``
        (|f_head + f_actin|, must be 0), and ``f_actin_along_walk`` (signed projection: < 0 = contractile).
    """
    if k_xb <= 0.0:
        raise ValueError("k_xb must be positive")
    head = np.asarray(head_pos, dtype=np.float64)
    anchor = np.asarray(anchor_pos, dtype=np.float64)
    walk = np.asarray(walk_dir, dtype=np.float64)
    x_att = attachment_point(anchor, abscissa, walk)
    d = x_att - head
    length = float(np.linalg.norm(d))
    if length < 1e-12:
        f_on_head = np.zeros(3)
    else:
        f_on_head = k_xb * (length - r0_xb) * (d / length)     # head pulled toward the advanced attachment
    f_on_actin = -f_on_head                                     # Newton's 3rd law
    load = k_xb * (float(np.dot(d, walk)) - r0_xb)              # tangential (along-walk) resisting tension
    return {
        "x_att": x_att,
        "f_on_head": f_on_head,
        "f_on_actin": f_on_actin,
        "load": load,
        "newton3rd_residual": float(np.linalg.norm(f_on_head + f_on_actin)),
        "f_actin_along_walk": float(np.dot(f_on_actin, walk)),
    }


def stall_abscissa(f_stall: float, k_xb: float, r0_xb: float = 0.0) -> float:
    """Isometric settled abscissa ``s* = F_stall/k_xb + r0_xb`` [µm] — the working-stroke strain at stall.

    The fixed point of the coupled walk where the tangential load equals ``F_stall`` (so ``hill_velocity`` = 0).
    Equals ``minifilament_topology.working_stroke_strain(f_stall, k_xb)`` when ``r0_xb = 0``.
    """
    if k_xb <= 0.0 or f_stall <= 0.0:
        raise ValueError("f_stall and k_xb must be positive")
    return f_stall / k_xb + r0_xb


def step_engaged_head(
    v0: float,
    f_stall: float,
    kappa: float,
    k_xb: float,
    *,
    r0_xb: float = 0.0,
    walk_dir: npt.ArrayLike = (1.0, 0.0, 0.0),
    tau: float = 1.0e-4,
    n_steps: int = 200_000,
    anchor_pos: npt.ArrayLike = (0.0, 0.0, 0.0),
) -> dict[str, npt.NDArray[np.float64]]:
    """Isometric single-engaged-head stepping (head + anchor held) — the force-velocity / stall demonstration.

    A freshly-bound head (``abscissa = 0``, head on its anchor) walks by ``ds = tau * hill_velocity(load)`` while
    the crossbridge tension it builds feeds back as the Hill load. With head + anchor held fixed the load rises
    monotonically toward ``F_stall``; the walk (hence the force) self-limits at stall. This is the exact 1-D
    reduction of the coupled device loop (compute_head_loads_kernel -> step_detach_kernel).

    Args:
        v0: unloaded stepping velocity [µm/s]. f_stall: per-head stall [pN]. kappa: Hill curvature.
        k_xb: crossbridge stiffness [pN/µm]. r0_xb: crossbridge rest [µm] (~0).
        walk_dir: unit walk direction. tau: KMC tick [s]. n_steps: ticks. anchor_pos: held actin node.

    Returns:
        dict of ``(n_steps + 1,)`` trajectories: ``t``, ``abscissa`` [µm], ``load`` [pN], ``velocity`` [µm/s],
        ``f_actin_along_walk`` [pN] (signed; contractile is negative).
    """
    walk = np.asarray(walk_dir, dtype=np.float64)
    walk = walk / (np.linalg.norm(walk) + 1e-30)
    anchor = np.asarray(anchor_pos, dtype=np.float64)
    head = anchor.copy()                                       # head held on its anchor (isometric)

    s = np.zeros(n_steps + 1)
    load = np.zeros(n_steps + 1)
    vel = np.zeros(n_steps + 1)
    f_par = np.zeros(n_steps + 1)
    for i in range(n_steps + 1):
        xb = crossbridge_force(head, anchor, s[i], walk, k_xb, r0_xb)
        load[i] = xb["load"]
        f_par[i] = xb["f_actin_along_walk"]
        vel[i] = float(hill_velocity(load[i], v0, f_stall, kappa, clamp=True))
        if i < n_steps:
            s[i + 1] = s[i] + tau * vel[i]
    return {
        "t": np.arange(n_steps + 1) * tau,
        "abscissa": s,
        "load": load,
        "velocity": vel,
        "f_actin_along_walk": f_par,
    }
