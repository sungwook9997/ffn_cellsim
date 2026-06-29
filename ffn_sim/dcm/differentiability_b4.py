"""B4 differentiability spot-check — Warp reverse-mode autodiff through a force law.

This is the Warp spike's actual PRIZE (plan §B4): the native CUDA compartment
plugin (``radial_shell_force.cu``) is NOT differentiable, so a clean reverse-mode
gradient w.r.t. a physical parameter is the thing Warp buys that the native path
cannot. Here we differentiate a scalar loss through the REAL membrane radial-shell
per-bead force kernel + an atomic-add reduction, w.r.t. the bare membrane tension
``γ_mem`` (a "shell stiffness"-class parameter), and cross-check the autodiff
gradient against BOTH a closed-form analytic gradient AND a central
finite-difference — so the check is self-validating, not graded on a fresh oracle.

Setup (centroid + R_mean precomputed and held fixed, so the gradient is purely
w.r.t. the parameter, isolating the autodiff of the force law):

    F_i  = -(ΔP · A_i) · n̂_i ,  ΔP = 2·(γ_mem + K_A·(S-A0)/A0)/R_mean ,  A_i = S/n
    loss = Σ_i |F_i|²  = (ΔP · A_i)² · n   (|n̂_i| = 1 for active beads)

⇒ ΔP is linear in γ_mem, loss is quadratic ⇒ closed form:

    d(loss)/dγ_mem = 2·(ΔP·A_i)·A_i·(2/R_mean)·n
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def _membrane_force_sq_loss(
    pos: wp.array(dtype=wp.vec3d),
    cx: wp.float64,
    cy: wp.float64,
    cz: wp.float64,
    R_mean: wp.float64,
    cnt: wp.float64,
    gamma_mem: wp.array(dtype=wp.float64),  # length-1, requires_grad
    K_A: wp.float64,
    A0: wp.float64,
    loss: wp.array(dtype=wp.float64),       # length-1, accumulator
):
    """Per-bead |F|² for the membrane Laplace law, summed into loss[0]."""
    i = wp.tid()
    PI = wp.float64(3.14159265358979323846)
    S = wp.float64(4.0) * PI * R_mean * R_mean
    gamma_tot = gamma_mem[0] + K_A * (S - A0) / A0
    dP = wp.float64(2.0) * gamma_tot / R_mean
    A_i = S / cnt
    Fmag = -(dP * A_i)

    p = pos[i]
    dx = p[0] - cx
    dy = p[1] - cy
    dz = p[2] - cz
    r = wp.sqrt(dx * dx + dy * dy + dz * dz)
    rs = r
    if r <= wp.float64(0.0):
        rs = wp.float64(1.0)
    nx = dx / rs
    ny = dy / rs
    nz = dz / rs
    fx = Fmag * nx
    fy = Fmag * ny
    fz = Fmag * nz
    wp.atomic_add(loss, 0, fx * fx + fy * fy + fz * fz)


def membrane_loss_and_grad(
    *, pos: np.ndarray, gamma_mem: float, K_A: float, A0: float, device: str = "cpu",
) -> tuple[float, float]:
    """Return (loss, d loss / d gamma_mem) via Warp reverse-mode autodiff."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    centroid = pos.mean(axis=0)
    radii = np.linalg.norm(pos - centroid, axis=1)
    R_mean = float(radii.mean())
    cnt = float(N)

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    gm = wp.array([gamma_mem], dtype=wp.float64, device=device, requires_grad=True)
    loss = wp.zeros(1, dtype=wp.float64, device=device, requires_grad=True)

    tape = wp.Tape()
    with tape:
        wp.launch(
            _membrane_force_sq_loss, dim=N,
            inputs=[pos_d,
                    wp.float64(centroid[0]), wp.float64(centroid[1]),
                    wp.float64(centroid[2]),
                    wp.float64(R_mean), wp.float64(cnt),
                    gm, wp.float64(K_A), wp.float64(A0)],
            outputs=[loss],
            device=device,
        )
    tape.backward(loss=loss)
    return float(loss.numpy()[0]), float(gm.grad.numpy()[0])


def analytic_grad(*, pos: np.ndarray, gamma_mem: float, K_A: float, A0: float) -> float:
    """Closed-form d(loss)/dγ_mem = 2·(ΔP·A_i)·A_i·(2/R_mean)·N."""
    N = pos.shape[0]
    centroid = pos.mean(axis=0)
    R_mean = float(np.linalg.norm(pos - centroid, axis=1).mean())
    S = 4.0 * np.pi * R_mean * R_mean
    gamma_tot = gamma_mem + K_A * (S - A0) / A0
    dP = 2.0 * gamma_tot / R_mean
    A_i = S / N
    return 2.0 * (dP * A_i) * A_i * (2.0 / R_mean) * N


def fd_grad(*, pos, gamma_mem, K_A, A0, h_rel=1e-6, device="cpu") -> float:
    """Central finite-difference d(loss)/dγ_mem (independent numeric check)."""
    h = h_rel * max(abs(gamma_mem), 1e-30)
    lp, _ = membrane_loss_and_grad(pos=pos, gamma_mem=gamma_mem + h, K_A=K_A, A0=A0, device=device)
    lm, _ = membrane_loss_and_grad(pos=pos, gamma_mem=gamma_mem - h, K_A=K_A, A0=A0, device=device)
    return (lp - lm) / (2.0 * h)


def run_check(device: str = "cpu") -> dict:
    """Differentiate the membrane force-law loss w.r.t. γ_mem three ways."""
    import os
    fix = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
    fx = dict(np.load(os.path.join(fix, "radial_ref_membrane.npz")))
    pos = fx["pos"]
    gamma_mem = float(fx["pa"])  # pa == gamma_mem for the membrane law
    K_A = float(fx["pb"])
    A0 = float(fx["pc"])

    loss, g_ad = membrane_loss_and_grad(pos=pos, gamma_mem=gamma_mem, K_A=K_A, A0=A0, device=device)
    g_an = analytic_grad(pos=pos, gamma_mem=gamma_mem, K_A=K_A, A0=A0)
    g_fd = fd_grad(pos=pos, gamma_mem=gamma_mem, K_A=K_A, A0=A0, device=device)

    rel_an = abs(g_ad - g_an) / (abs(g_an) + 1e-300)
    rel_fd = abs(g_ad - g_fd) / (abs(g_fd) + 1e-300)
    return {
        "loss": loss,
        "grad_autodiff": g_ad,
        "grad_analytic": g_an,
        "grad_finite_diff": g_fd,
        "rel_err_vs_analytic": rel_an,
        "rel_err_vs_finite_diff": rel_fd,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(run_check(), indent=2))
