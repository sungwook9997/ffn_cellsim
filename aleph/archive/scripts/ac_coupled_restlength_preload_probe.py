"""Coupled rest-length preload prototype (option a) — make the AS-BUILT geometry force-balanced.

FIRE diverges trying to MOVE nodes to equilibrium (soft membrane vs stiff cortex CFL). Instead, hold the
positions fixed and adjust REST LENGTHS so the built state is already balanced — this IS the physiological
resting pre-stress (ERM tension holds membrane turgor; cortex crosslinks carry the Laplace hoop reaction).
The ERM half already exists (_preload_erm_resting_balance); this extends it to the crosslink network by a
Jacobi relaxation of r0 that zeroes the projected nodal force along each link.

Out-of-hot-loop host relaxation using the exact GPU force kernels (accumulate -> read -> adjust rest ->
assign). Prototype only; if it drives the projected residual under gate at subdiv 6, integrate into driver.
CUDA-gated -> gbook only.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.incumbent.assemble import build_cell, CellConfig
from aleph.components.incumbent.driver import _accumulate_all, _preload_erm_resting_balance, _residual_host


def _force(cell):
    cell.f_d.zero_()
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    return cell.f_d.numpy()


def coupled_preload(cell, sweeps: int = 40, sign: float = +1.0, damp: float = 1.0, verbose=True):
    """Alternating ERM + crosslink rest-length relaxation to zero the residual at fixed positions."""
    xl = cell.xl_d.numpy().astype(np.int64)
    k = cell.kxl_d.numpy()
    i, j = xl[:, 0], xl[:, 1]
    res0, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    if verbose:
        print(f"    sweep  0: residual {res0:.5f}")
    for s in range(1, sweeps + 1):
        _preload_erm_resting_balance(cell, sweeps=1)          # hold membrane (existing)
        f = _force(cell)
        pos = cell.pos_d.numpy()
        d = pos[j] - pos[i]
        L = np.linalg.norm(d, axis=1)
        ok = L > 1e-12
        u = np.zeros_like(d); u[ok] = d[ok] / L[ok, None]
        # projected antisymmetric residual along each link: 0.5*(f_i - f_j)·û
        g = 0.5 * (np.einsum("nc,nc->n", f[i] - f[j], u))
        r0 = cell.r0xl_d.numpy()
        r0_new = np.clip(r0 + sign * damp * g / np.maximum(k, 1e-30), 1e-6, None)
        cell.r0xl_d.assign(np.ascontiguousarray(r0_new, np.float64))
        if verbose and (s <= 5 or s % 5 == 0):
            res, _ = _residual_host(cell, cell.pos_d, cell.f_d)
            print(f"    sweep {s:2d}: residual {res:.5f}")
    res, _ = _residual_host(cell, cell.pos_d, cell.f_d)
    return res


def probe(n_filaments=8000, subdiv=6):
    for sign in (+1.0, -1.0):
        cfg = CellConfig(n_filaments=n_filaments, overlap_free_cortex=True,
                         erm_radial_pairing=True, membrane_subdivisions=subdiv)
        cell = build_cell(cfg)
        print("=" * 80)
        print(f"n={n_filaments} subdiv={subdiv} sign={sign:+.0f}")
        res = coupled_preload(cell, sweeps=40, sign=sign)
        print(f"  >> final residual {res:.5f} pN")


if __name__ == "__main__":
    import sys
    nf = 70686 if "--native" in sys.argv else 8000
    probe(nf, 6)
