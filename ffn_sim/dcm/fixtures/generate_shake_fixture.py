"""Generate + commit the HOOMD reference fixture for the M-SHAKE Warp parity test.

Mirrors the native SHAKE parity test config (``test_shake_parity.py``): F randomly
oriented straight chains at rest length, perturbed off the constraint manifold,
projected by the COMMITTED Python reference
``integrator.constrained_baoab.shake_project_chains`` (rigid, return_lambdas=True).
The projected positions + accumulated Lagrange multipliers are the bit-parity
TARGET for ``tests/warp_port/test_shake_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_shake_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd

from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import shake_project_chains

HERE = os.path.dirname(os.path.abspath(__file__))

F, M = 1000, 6          # 1000 chains, 6 bonds (7 beads) each
N = F * (M + 1)
L = 1.0e6               # huge box → min-image inactive (isolates the solver)
R0 = 100.0              # rest length
TOL, MAX_ITER = 1.0e-10, 100
SEED = 3


def main() -> None:
    rng = np.random.default_rng(SEED)
    ref = np.zeros((N, 3))
    chains = np.zeros((F, M + 1), dtype=np.int32)
    for f in range(F):
        base = f * (M + 1)
        center = (rng.random(3) - 0.5) * (L * 1e-3)
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        for b in range(M + 1):
            ref[base + b] = center + (b - M / 2.0) * R0 * d
            chains[f, b] = base + b
    pred = ref + rng.normal(0.0, R0 * 0.05, size=(N, 3))
    inv_mass = np.ones(N)
    box = hoomd.Box(Lx=L, Ly=L, Lz=L)

    proj, lam = shake_project_chains(
        pred.copy(), ref.copy(), chains.copy(), R0, inv_mass.copy(), box,
        tol=TOL, max_iter=MAX_ITER, return_lambdas=True,
    )
    proj = np.asarray(proj)
    lam = np.asarray(lam)

    # constraint satisfaction of the reference projection (sanity)
    s = proj[chains[:, :-1]] - proj[chains[:, 1:]]
    bond = np.sqrt(np.sum(s * s, axis=2))
    max_bond_err = float(np.max(np.abs(bond - R0)) / R0)

    out = os.path.join(HERE, "shake_ref.npz")
    np.savez(
        out,
        F=F, m=M, N=N, L=L, r0=R0, tol=TOL, max_iter=MAX_ITER,
        pred=pred, ref=ref, chains=chains, inv_mass=inv_mass,
        ref_proj=proj, ref_lam=lam,
    )
    print(f"wrote {out}  F={F} m={M}  max bond rel-err={max_bond_err:.2e}  "
          f"max|proj|={np.abs(proj).max():.3e}  max|lam|={np.abs(lam).max():.3e}")


if __name__ == "__main__":
    main()
