"""End-to-end DIFFERENTIABLE DCM simulation loop (optimization #8 — Warp's real prize).

The B4 spot-check differentiated a single force evaluation. This differentiates the
WHOLE K-step simulation loop: run K deterministic overdamped steps (turgor + cortex
edge springs + integrator) under a ``wp.Tape``, define a scalar loss on the FINAL
state, and back-propagate to the physical PARAMETERS (turgor dP0, edge stiffness
k_edge). That is what the native HOOMD runtime cannot do at all, and the basis for
parameter estimation / inverse design / learnable constants on the DCM.

Design for clean reverse-mode through a time loop:
  * parameters are length-1 ``requires_grad`` arrays read as ``p[0]`` in the kernels
    (autodiff differentiates w.r.t. array inputs).
  * each step writes a FRESH ``pos`` array (no in-place overwrite) so every
    intermediate state the backward pass needs is preserved — the textbook
    Warp differentiable-sim pattern. (Memory = K x per-step arrays; gradient
    checkpointing is the scaling path for large K, noted below.)
  * kT = 0 (deterministic): the gradient is exact; thermal noise would add a
    parameter-independent term with zero gradient anyway.

Validation: autodiff d(loss)/d{dP0,k_edge} vs central finite differences.
"""

from __future__ import annotations

import numpy as np

import warp as wp

from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM

wp.init()


@wp.kernel
def _reduce_centroid(pos: wp.array(dtype=wp.vec3d), acc: wp.array(dtype=wp.float64)):
    i = wp.tid()
    p = pos[i]
    wp.atomic_add(acc, 0, p[0])
    wp.atomic_add(acc, 1, p[1])
    wp.atomic_add(acc, 2, p[2])
    wp.atomic_add(acc, 3, wp.float64(1.0))


@wp.kernel
def _reduce_radius(pos: wp.array(dtype=wp.vec3d), acc: wp.array(dtype=wp.float64)):
    i = wp.tid()
    cnt = acc[3]
    dx = pos[i][0] - acc[0] / cnt
    dy = pos[i][1] - acc[1] / cnt
    dz = pos[i][2] - acc[2] / cnt
    wp.atomic_add(acc, 4, wp.sqrt(dx * dx + dy * dy + dz * dz))


@wp.kernel
def _turgor_force(pos: wp.array(dtype=wp.vec3d), acc: wp.array(dtype=wp.float64),
                  dP0: wp.array(dtype=wp.float64), K_vol: wp.float64, V0: wp.float64,
                  force: wp.array(dtype=wp.vec3d)):
    i = wp.tid()
    cnt = acc[3]
    cx = acc[0] / cnt
    cy = acc[1] / cnt
    cz = acc[2] / cnt
    R_mean = acc[4] / cnt
    PI = wp.float64(3.14159265358979323846)
    V = (wp.float64(4.0) / wp.float64(3.0)) * PI * R_mean * R_mean * R_mean
    S = wp.float64(4.0) * PI * R_mean * R_mean
    dP = dP0[0] - K_vol * (V - V0) / V0
    A_i = S / cnt
    Fmag = dP * A_i
    dx = pos[i][0] - cx
    dy = pos[i][1] - cy
    dz = pos[i][2] - cz
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    force[i] = wp.vec3d(Fmag * dx / r, Fmag * dy / r, Fmag * dz / r)


@wp.kernel
def _edge_force(pos: wp.array(dtype=wp.vec3d), edges: wp.array(dtype=wp.int32, ndim=2),
                kedge: wp.array(dtype=wp.float64), r0: wp.array(dtype=wp.float64),
                force: wp.array(dtype=wp.vec3d)):
    e = wp.tid()
    i = edges[e, 0]
    j = edges[e, 1]
    dx = pos[i][0] - pos[j][0]
    dy = pos[i][1] - pos[j][1]
    dz = pos[i][2] - pos[j][2]
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    fmag = -kedge[0] * (r - r0[e]) / r
    wp.atomic_add(force, i, wp.vec3d(fmag * dx, fmag * dy, fmag * dz))
    wp.atomic_add(force, j, wp.vec3d(-fmag * dx, -fmag * dy, -fmag * dz))


@wp.kernel
def _integrate(pos_in: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
               inv_gamma: wp.float64, dt: wp.float64, pos_out: wp.array(dtype=wp.vec3d)):
    i = wp.tid()
    f = force[i]
    p = pos_in[i]
    pos_out[i] = wp.vec3d(p[0] + f[0] * inv_gamma * dt,
                          p[1] + f[1] * inv_gamma * dt,
                          p[2] + f[2] * inv_gamma * dt)


@wp.kernel
def _radius_loss(pos: wp.array(dtype=wp.vec3d), acc: wp.array(dtype=wp.float64),
                 inv_n: wp.float64, loss: wp.array(dtype=wp.float64)):
    """loss = mean bead radius about the centroid (a scalar shape observable)."""
    i = wp.tid()
    cnt = acc[3]
    dx = pos[i][0] - acc[0] / cnt
    dy = pos[i][1] - acc[1] / cnt
    dz = pos[i][2] - acc[2] / cnt
    wp.atomic_add(loss, 0, wp.sqrt(dx * dx + dy * dy + dz * dz) * inv_n)


def _forward(pos0, edges, r0, dP0v, kedgev, *, K, dt, K_vol, V0, inv_gamma,
             device, tape=None):
    N = pos0.shape[0]
    ne = edges.shape[0]
    dP0 = wp.array([dP0v], dtype=wp.float64, device=device, requires_grad=True)
    kedge = wp.array([kedgev], dtype=wp.float64, device=device, requires_grad=True)
    pos = wp.array(pos0, dtype=wp.vec3d, device=device, requires_grad=True)
    edges_d = wp.array(edges, dtype=wp.int32, device=device)
    r0_d = wp.array(r0, dtype=wp.float64, device=device)
    loss = wp.zeros(1, dtype=wp.float64, device=device, requires_grad=True)

    def body():
        nonlocal pos
        for _ in range(K):
            acc = wp.zeros(5, dtype=wp.float64, device=device, requires_grad=True)
            force = wp.zeros(N, dtype=wp.vec3d, device=device, requires_grad=True)
            wp.launch(_reduce_centroid, dim=N, inputs=[pos, acc], device=device)
            wp.launch(_reduce_radius, dim=N, inputs=[pos, acc], device=device)
            wp.launch(_turgor_force, dim=N,
                      inputs=[pos, acc, dP0, wp.float64(K_vol), wp.float64(V0), force],
                      device=device)
            wp.launch(_edge_force, dim=ne, inputs=[pos, edges_d, kedge, r0_d, force],
                      device=device)
            pos_next = wp.zeros(N, dtype=wp.vec3d, device=device, requires_grad=True)
            wp.launch(_integrate, dim=N,
                      inputs=[pos, force, wp.float64(inv_gamma), wp.float64(dt), pos_next],
                      device=device)
            pos = pos_next
        acc_f = wp.zeros(5, dtype=wp.float64, device=device, requires_grad=True)
        wp.launch(_reduce_centroid, dim=N, inputs=[pos, acc_f], device=device)
        wp.launch(_radius_loss, dim=N,
                  inputs=[pos, acc_f, wp.float64(1.0 / N), loss], device=device)

    if tape is not None:
        with tape:
            body()
    else:
        body()
    return loss, dP0, kedge


def run_diff(*, subdiv: int = 1, K: int = 40, dt: float = 2.0e-8,
             device: str = "cpu") -> dict:
    """Differentiate the K-step loss w.r.t. (dP0, k_edge) via autodiff + FD check."""
    p = ResolvedDCM(subdivisions=subdiv)
    verts, edges, tris = icosphere_mesh(p.R_cell, subdiv)
    edges = edges.astype(np.int32)
    r0 = np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1).astype(np.float64)
    R0 = float(np.linalg.norm(verts - verts.mean(0), axis=1).mean())
    V0 = (4.0 / 3.0) * np.pi * R0 ** 3
    inv_gamma = 1.0 / p.gamma_node
    common = dict(K=K, dt=dt, K_vol=p.K_vol, V0=V0, inv_gamma=inv_gamma, device=device)

    tape = wp.Tape()
    loss, dP0, kedge = _forward(verts, edges, r0, p.turgor_dP0, p.k_edge, tape=tape, **common)
    tape.backward(loss=loss)
    g_dP0 = float(dP0.grad.numpy()[0])
    g_kedge = float(kedge.grad.numpy()[0])
    L0 = float(loss.numpy()[0])

    def fd(which, rel=1e-4):
        h = rel * abs(getattr(p, which))
        kw = dict(dP0v=p.turgor_dP0, kedgev=p.k_edge)
        key = "dP0v" if which == "turgor_dP0" else "kedgev"
        kw[key] = getattr(p, which) + h
        lp, _, _ = _forward(verts, edges, r0, kw["dP0v"], kw["kedgev"], **common)
        kw[key] = getattr(p, which) - h
        lm, _, _ = _forward(verts, edges, r0, kw["dP0v"], kw["kedgev"], **common)
        return (float(lp.numpy()[0]) - float(lm.numpy()[0])) / (2.0 * h)

    fd_dP0 = fd("turgor_dP0")
    fd_kedge = fd("k_edge")
    return {
        "device": device, "subdiv": subdiv, "K": K, "N": int(verts.shape[0]),
        "loss_mean_radius": L0,
        "grad_dP0_autodiff": g_dP0, "grad_dP0_fd": fd_dP0,
        "grad_kedge_autodiff": g_kedge, "grad_kedge_fd": fd_kedge,
        "rel_err_dP0": abs(g_dP0 - fd_dP0) / (abs(fd_dP0) + 1e-300),
        "rel_err_kedge": abs(g_kedge - fd_kedge) / (abs(fd_kedge) + 1e-300),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(run_diff(), indent=2))
