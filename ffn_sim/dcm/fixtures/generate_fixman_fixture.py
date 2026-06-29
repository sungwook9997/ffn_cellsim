"""Generate + commit the HOOMD reference fixture for the Fixman Warp parity test.

Mirrors the native Fixman parity test config (``test_fixman_parity.py``): F random
chains (varying bond lengths + directions, varying per-bead mobility), reference
force + U_F from the COMMITTED Python ``fixman_logdet_and_force``. That is the
parity TARGET for ``tests/dcm/test_fixman_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_fixman_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd

from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import fixman_logdet_and_force

HERE = os.path.dirname(os.path.abspath(__file__))

F, M = 1000, 6
N = F * (M + 1)
L = 1.0e6
R0 = 100.0
KT = 1.0
SEED = 5


def main() -> None:
    rng = np.random.default_rng(SEED)
    pos = np.zeros((N, 3))
    chains = np.zeros((F, M + 1), dtype=np.int32)
    for f in range(F):
        base = f * (M + 1)
        p = (rng.random(3) - 0.5) * (L * 1e-3)
        for b in range(M + 1):
            chains[f, b] = base + b
            pos[base + b] = p
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            p = p + R0 * (1.0 + 0.1 * rng.normal()) * d
    inv_gamma = 1.0 + rng.random(N)  # varying per-bead mobility
    box = hoomd.Box(Lx=L, Ly=L, Lz=L)

    U_py, F_py = fixman_logdet_and_force(
        pos.copy(), chains.copy(), KT, inv_gamma.copy(), box)
    F_py = np.asarray(F_py)

    out = os.path.join(HERE, "fixman_ref.npz")
    np.savez(
        out,
        F=F, m=M, N=N, L=L, r0=R0, kT=KT,
        pos=pos, chains=chains, inv_gamma=inv_gamma,
        ref_force=F_py, ref_U=float(U_py),
    )
    print(f"wrote {out}  F={F} m={M}  U_F={float(U_py):.6e}  "
          f"max|F|={np.abs(F_py).max():.3e}")


if __name__ == "__main__":
    main()
