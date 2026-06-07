"""Native M-SHAKE (Stage 1b) parity vs constrained_baoab.shake_project_chains.

Builds F randomly-oriented straight chains at rest length, perturbs off the
constraint manifold, and projects with BOTH the Python reference (numpy fast
path) and the native CUDA kernel. Checks: projected positions + accumulated
Lagrange multipliers agree, and the projected bonds satisfy the constraint.

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python test_shake_parity.py
"""

from __future__ import annotations

import numpy as np
import hoomd

from ffn_sim.integrator.constrained_baoab import shake_project_chains
import _ffn_native


def main() -> int:
    F, m = 3000, 6          # 3000 chains, 6 bonds (7 beads) each
    N = F * (m + 1)
    L = 1.0e6               # huge box → min-image inactive (isolates the solver)
    r0 = 100.0              # rest length
    tol, max_iter = 1.0e-10, 100
    rng = np.random.default_rng(3)

    ref = np.zeros((N, 3))
    chains = np.zeros((F, m + 1), dtype=np.int32)
    for f in range(F):
        base = f * (m + 1)
        center = (rng.random(3) - 0.5) * (L * 1e-3)
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        for b in range(m + 1):
            ref[base + b] = center + (b - m / 2.0) * r0 * d
            chains[f, b] = base + b
    pred = ref + rng.normal(0.0, r0 * 0.05, size=(N, 3))
    inv_mass = np.ones(N)
    box = hoomd.Box(Lx=L, Ly=L, Lz=L)

    proj_py, lam_py = shake_project_chains(
        pred.copy(), ref.copy(), chains.copy(), r0, inv_mass.copy(), box,
        tol=tol, max_iter=max_iter, return_lambdas=True,
    )
    proj_nat, lam_nat, nc = _ffn_native.shake_project(
        pred.copy(), ref.copy(), inv_mass.copy(), chains.copy(),
        r0, L, L, L, tol, max_iter,
    )
    proj_nat = np.asarray(proj_nat)
    lam_nat = np.asarray(lam_nat)
    nc = np.asarray(nc)

    dpos = float(np.max(np.abs(proj_nat - np.asarray(proj_py))))
    dlam = float(np.max(np.abs(lam_nat - np.asarray(lam_py))))
    pscale = float(np.max(np.abs(np.asarray(proj_py))))
    lscale = float(np.max(np.abs(np.asarray(lam_py)))) + 1e-30

    # native constraint satisfaction (bonds == r0)
    s = proj_nat[chains[:, :-1]] - proj_nat[chains[:, 1:]]
    bond = np.sqrt(np.sum(s * s, axis=2))
    max_bond_err = float(np.max(np.abs(bond - r0)) / r0)

    pos_ok = dpos < 1e-8 * max(pscale, 1.0)
    lam_ok = dlam < 1e-8 * lscale
    sat_ok = max_bond_err < 1e-8
    conv_ok = int(nc.sum()) == 0

    print(f"[shake] F={F} m={m}  nonconverged={int(nc.sum())}")
    print(f"[shake] max|pos_nat-pos_py| = {dpos:.3e}  (scale {pscale:.3e})  "
          f"{'PASS' if pos_ok else 'FAIL'}")
    print(f"[shake] max|lam_nat-lam_py| = {dlam:.3e}  (scale {lscale:.3e})  "
          f"{'PASS' if lam_ok else 'FAIL'}")
    print(f"[shake] max bond rel-err    = {max_bond_err:.3e}  "
          f"{'PASS' if sat_ok else 'FAIL'}")

    if pos_ok and lam_ok and sat_ok and conv_ok:
        print("native-shake parity PASS")
        return 0
    print("native-shake parity FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
