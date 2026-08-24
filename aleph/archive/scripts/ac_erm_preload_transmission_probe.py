"""ERM preload force-transmission probe (does the membrane<->cortex load path now carry turgor?).

Handoff 2026-07-22b: after the inward-dispersion geometry fix (cortex outer surface at R_CORTEX, ERM
tethers radial/non-degenerate), verify that `_preload_erm_resting_balance` now transmits the full ~40 pN
membrane turgor to the cortex (before the fix it moved only ~11 of 40 pN because the tethers were tangential).

Reports, per config, the membrane-node max/mean |F| BEFORE and AFTER the ERM preload, the cortex reaction
load, and the whole-cell residual. CUDA-gated (build_cell) -> gbook only.
"""
from __future__ import annotations

import numpy as np

from aleph.components.incumbent.assemble import build_cell, CellConfig
from aleph.components.incumbent.driver import _accumulate_all, _preload_erm_resting_balance, _preload_cortex_pretension
import warp as wp


def _force_now(cell):
    cell.f_d.zero_()
    _accumulate_all(cell, cell.pos_d, cell.f_d)
    wp.synchronize_device(cell.device)
    return cell.f_d.numpy()


def _stats(fmag):
    return f"max={fmag.max():8.3f}  mean={fmag.mean():8.3f}  p95={np.percentile(fmag,95):8.3f}"


def probe(n_filaments: int, radial: bool, prestrain: float, membrane_subdivisions: int = 3) -> None:
    cfg = CellConfig(n_filaments=n_filaments, overlap_free_cortex=True, erm_radial_pairing=radial,
                     membrane_subdivisions=membrane_subdivisions)
    cell = build_cell(cfg)
    mem = cell.membrane
    erm_m = mem.erm_m_d.numpy().astype(np.int64)
    mem_nodes = np.unique(erm_m)
    n_actin = cell.n_actin

    tag = f"n={n_filaments:>6d} subdiv={membrane_subdivisions} radial_pair={radial!s:5} prestrain={prestrain}"
    print("=" * 96)
    print(tag, f" | n_verts={mem.n_verts} n_erm={mem.n_erm} f_rupt={mem.f_rupt:.2f}pN k_erm={mem.k_erm:g}")

    f0 = _force_now(cell)
    fm0 = np.linalg.norm(f0, axis=1)
    print(f"  BEFORE preload | membrane |F|: {_stats(fm0[mem_nodes])}")
    print(f"                 | cortex   |F|: {_stats(fm0[:n_actin])}")
    print(f"                 | whole    |F|: max={fm0.max():.3f}")

    if prestrain:
        info = _preload_cortex_pretension(cell, prestrain)
        print(f"  cortex pretension applied: {info}")
    info = _preload_erm_resting_balance(cell)
    print(f"  ERM preload: n_erm={info['n_erm']} res_before={info['res_before']:.3f} res_after={info['res_after']:.3f}")

    f1 = _force_now(cell)
    fm1 = np.linalg.norm(f1, axis=1)
    print(f"  AFTER  preload | membrane |F|: {_stats(fm1[mem_nodes])}")
    print(f"                 | cortex   |F|: {_stats(fm1[:n_actin])}")
    print(f"                 | whole    |F|: max={fm1.max():.3f}")
    drop = 100.0 * (1.0 - fm1[mem_nodes].max() / max(fm0[mem_nodes].max(), 1e-12))
    print(f"  >> membrane max|F| reduction: {drop:5.1f} %  ({fm0[mem_nodes].max():.2f} -> {fm1[mem_nodes].max():.2f} pN)")


if __name__ == "__main__":
    import sys
    # (n_filaments, radial_pair, prestrain, membrane_subdivisions)
    configs = [
        (8000, True, 0.0, 3),
        (8000, True, 0.0, 4),
        (8000, True, 0.0, 6),
    ]
    if "--native" in sys.argv:
        configs = [
            (70686, True, 0.0, 3),
            (70686, True, 0.0, 6),
        ]
    for nf, rad, ps, sub in configs:
        probe(nf, rad, ps, membrane_subdivisions=sub)
