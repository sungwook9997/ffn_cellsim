r"""Arp2/3 branch-angle harmonic force — Warp CUDA kernel SOURCE (I4).

Device realization of the 3-body angle-harmonic Arp2/3 junction (``ac.weave.branch_angle`` is the host
oracle). One thread per branch triple ``(mother, branch, daughter)`` adds the harmonic-angle force
``U = 1/2 k_theta (theta - theta0)^2`` to the global node-force accumulator; the angle FLUCTUATES thermally
(it is a soft potential, NOT a rigid 72 deg constraint — CLAUDE.md worked-example).

This is a NET-NEW port of ``ff.network_warp.branch_angle_kernel`` into the ``ac/`` engine (same standard
harmonic-angle expression, so the native gate reproduces the CPU oracle). It is authored on the dev Mac but
NEVER launched here (I0-A: Warp runs on CUDA GPU only; the Mac has no CUDA). The lead wires it into the
inner force-assembly loop on the gbook A5000; ``ff/network_warp.py`` is untouched from this worktree — the
retire/replace note lives in ``ac/weave/INTEGRATION.md``.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). HOOMD is never imported.

Sanity Gate (native, lead): device branch-angle force == host ``ac.weave.branch_angle.angle_forces`` to
tolerance; the three nodal forces sum to zero (internal potential); a bent junction restores toward theta0.
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

__all__ = ["branch_angle_kernel"]


@wp.kernel
def branch_angle_kernel(
    pos: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),   # (B, 3) [mother i, branch/center j, daughter k]
    active: wp.array(dtype=wp.int32),            # (B,) 1 = live branch (dormant daughters => 0, no force)
    theta0: wp.float64,                          # rest branch angle [rad] (Arp2/3 ~ 70 deg = 1.2217 rad)
    k_theta: wp.float64,                         # angular stiffness [pN*um/rad^2] = kT/sigma_theta^2
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Harmonic angle force ``U = 1/2 k_theta (theta - theta0)^2`` at each ACTIVE branch triple.

    theta is the angle at branch node j between (i-j) and (k-j). Standard angle force
    ``F_i = pref (r2/(d1 d2) - c r1/d1^2)``, ``F_k = pref (r1/(d1 d2) - c r2/d2^2)``, ``F_j = -(F_i+F_k)``
    with ``pref = k_theta (theta - theta0)/sin(theta)`` — identical to the host oracle and to the ff kernel.
    Dormant (inactive) branches contribute no force so the N-fixed dormant-daughter pool is force-neutral
    until a nucleation event activates it (the dendritic topology emerges without allocating new nodes).
    """
    t = wp.tid()
    if active[t] == 0:
        return
    i = triples[t, 0]
    j = triples[t, 1]
    k = triples[t, 2]
    r1 = pos[i] - pos[j]
    r2 = pos[k] - pos[j]
    d1 = wp.length(r1)
    d2 = wp.length(r2)
    if d1 < wp.float64(1e-12) or d2 < wp.float64(1e-12):
        return
    c = wp.dot(r1, r2) / (d1 * d2)
    c = wp.clamp(c, wp.float64(-1.0), wp.float64(1.0))
    s = wp.sqrt(wp.max(wp.float64(1.0) - c * c, wp.float64(1e-12)))
    theta = wp.acos(c)
    pref = k_theta * (theta - theta0) / s
    fi = pref * (r2 / (d1 * d2) - c * r1 / (d1 * d1))
    fk = pref * (r1 / (d1 * d2) - c * r2 / (d2 * d2))
    wp.atomic_add(force, i, fi)
    wp.atomic_add(force, k, fk)
    wp.atomic_add(force, j, -(fi + fk))
