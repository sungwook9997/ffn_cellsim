"""Native Fixman pseudo-force (Stage 1b) parity vs fixman_logdet_and_force.

Random-orientation, varying-length chains with varying per-bead mobility →
exercises the full tridiagonal Gram metric. Compares the native CUDA Fixman
force + U_F against the Python (numpy linalg) reference. Parity is tolerance-
based: the native uses tridiagonal continuants (θ/φ) while numpy uses LAPACK
LU, so they agree to a tight relative tolerance, not bitwise.

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python test_fixman_parity.py
"""

from __future__ import annotations

import numpy as np
import hoomd

from ffn_sim.integrator.constrained_baoab import fixman_logdet_and_force
import _ffn_native


def main() -> int:
    F, m = 3000, 6
    N = F * (m + 1)
    L = 1.0e6
    r0 = 100.0
    kT = 1.0
    rng = np.random.default_rng(5)

    pos = np.zeros((N, 3))
    chains = np.zeros((F, m + 1), dtype=np.int32)
    for f in range(F):
        base = f * (m + 1)
        p = (rng.random(3) - 0.5) * (L * 1e-3)
        for b in range(m + 1):
            chains[f, b] = base + b
            pos[base + b] = p
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            p = p + r0 * (1.0 + 0.1 * rng.normal()) * d
    inv_gamma = 1.0 + rng.random(N)  # varying per-bead mobility
    box = hoomd.Box(Lx=L, Ly=L, Lz=L)

    U_py, F_py = fixman_logdet_and_force(pos.copy(), chains.copy(), kT, inv_gamma.copy(), box)
    F_nat, U_nat, bad = _ffn_native.fixman_force(
        pos.copy(), inv_gamma.copy(), chains.copy(), kT, L, L, L
    )
    F_nat = np.asarray(F_nat)
    bad = np.asarray(bad)
    F_py = np.asarray(F_py)

    dF = float(np.max(np.abs(F_nat - F_py)))
    Fscale = float(np.max(np.abs(F_py))) + 1e-30
    dU = abs(float(U_nat) - float(U_py))
    Uscale = abs(float(U_py)) + 1e-30

    f_ok = dF < 1e-7 * Fscale
    u_ok = dU < 1e-9 * Uscale
    sign_ok = int(bad.sum()) == 0

    print(f"[fixman] F={F} m={m}  bad_sign={int(bad.sum())}")
    print(f"[fixman] max|F_nat-F_py| = {dF:.3e}  (scale {Fscale:.3e}, rel {dF/Fscale:.2e})  "
          f"{'PASS' if f_ok else 'FAIL'}")
    print(f"[fixman] |U_nat-U_py|    = {dU:.3e}  (scale {Uscale:.3e}, rel {dU/Uscale:.2e})  "
          f"{'PASS' if u_ok else 'FAIL'}")

    if f_ok and u_ok and sign_ok:
        print("native-fixman parity PASS")
        return 0
    print("native-fixman parity FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
