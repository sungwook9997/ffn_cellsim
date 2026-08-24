"""ERM tether geometry probe (membrane<->cortex load-path diagnosis).

Handoff 2026-07-22b §DO THIS FIRST. Hypothesis: the overlap-free cortex radial thickness (0.2 um)
spreads cortex nodes across R~7.30-7.50 up to R_mem=7.5, so membrane<->cortex ERM tethers become
near-degenerate (length ~0) and/or purely TANGENTIAL, so they cannot transmit the radial turgor load.

Measures, per (overlap_free, n_filaments) config:
  - tether length distribution |pos[erm_m]-pos[erm_c]|  (== erm_rest at build; force-free rest).
  - radial fraction of each tether = |d . r_hat| / |d|  (1 = purely radial/load-bearing, 0 = tangential).
  - cortex node radius min/mean/max ; membrane node radius min/mean/max.
  - how many cortex nodes reach out to R_mem (radial gap to membrane ~ 0).

CUDA-gated (build_cell, I0-A) -> gbook only.
"""
from __future__ import annotations

import numpy as np

from aleph.components.incumbent.assemble import build_cell, CellConfig, R_CELL_UM, R_CORTEX_UM


def _pct(a, ps=(0, 1, 5, 50, 95, 99, 100)):
    return {p: float(np.percentile(a, p)) for p in ps}


def probe(n_filaments: int, overlap_free: bool, radial: bool) -> None:
    cfg = CellConfig(
        n_filaments=n_filaments,
        overlap_free_cortex=overlap_free,
        erm_radial_pairing=radial,
    )
    cell = build_cell(cfg)
    pos = cell.pos_d.numpy()[: cell.n_actin]        # actin/cortex nodes
    allpos = cell.pos_d.numpy()
    mem = cell.membrane
    erm_m = mem.erm_m_d.numpy().astype(np.int64)
    erm_c = mem.erm_c_d.numpy().astype(np.int64)
    erm_rest = mem.erm_rest_d.numpy()

    c = np.zeros(3)
    r_cortex = np.linalg.norm(pos - c, axis=1)
    pm = allpos[erm_m]
    pc = allpos[erm_c]
    r_mem = np.linalg.norm(pm - c, axis=1)
    r_c_paired = np.linalg.norm(pc - c, axis=1)

    d = pm - pc
    L = np.linalg.norm(d, axis=1)
    mid = 0.5 * (pm + pc)
    rhat = mid / np.linalg.norm(mid, axis=1, keepdims=True)
    # radial fraction: |component of tether along radial| / length
    with np.errstate(invalid="ignore", divide="ignore"):
        radial_comp = np.abs(np.einsum("ij,ij->i", d, rhat))
        radial_frac = np.where(L > 1e-12, radial_comp / L, np.nan)
    # radial gap between the two shells for each pair (signed: membrane outside cortex => +)
    radial_gap = r_mem - r_c_paired

    tag = f"n={n_filaments:>6d} overlap_free={overlap_free!s:5} radial_pair={radial!s:5}"
    print("=" * 100)
    print(tag)
    print(f"  pairing_mode      : {mem.erm_pairing_mode}")
    print(f"  n_erm tethers     : {L.shape[0]}   (membrane nodes = {mem.n_verts})")
    print(f"  R_CELL/mem={R_CELL_UM}  R_CORTEX_centre={R_CORTEX_UM}")
    print(f"  cortex node R     : min={r_cortex.min():.4f} mean={r_cortex.mean():.4f} max={r_cortex.max():.4f}  (span={r_cortex.max()-r_cortex.min():.4f})")
    print(f"  membrane node R   : min={r_mem.min():.4f} mean={r_mem.mean():.4f} max={r_mem.max():.4f}")
    print(f"  paired-cortex R   : min={r_c_paired.min():.4f} mean={r_c_paired.mean():.4f} max={r_c_paired.max():.4f}")
    print(f"  tether length [um]: min={L.min():.5f} mean={L.mean():.5f} max={L.max():.5f}")
    print(f"    length pcts     : {_pct(L)}")
    print(f"  radial GAP  [um]  : min={radial_gap.min():.5f} mean={radial_gap.mean():.5f} max={radial_gap.max():.5f}")
    print(f"    radial_gap pcts : {_pct(radial_gap)}")
    print(f"  radial FRAC (1=radial,0=tangential): min={np.nanmin(radial_frac):.4f} mean={np.nanmean(radial_frac):.4f} max={np.nanmax(radial_frac):.4f}")
    print(f"    radial_frac pcts: {_pct(radial_frac[np.isfinite(radial_frac)])}")
    # degeneracy counts
    print(f"  tethers < 0.010 um (near-degenerate len): {int((L < 0.010).sum())} / {L.shape[0]}")
    print(f"  tethers radial_frac < 0.30 (mostly tangential): {int((radial_frac < 0.30).sum())} / {L.shape[0]}")
    print(f"  tethers radial_frac > 0.70 (mostly radial)    : {int((radial_frac > 0.70).sum())} / {L.shape[0]}")
    # cortex nodes reaching membrane
    near_mem = (r_cortex > (R_CELL_UM - 0.010)).sum()
    print(f"  cortex nodes with R > R_mem-0.01 (touch membrane): {int(near_mem)} / {r_cortex.shape[0]}")


if __name__ == "__main__":
    import sys
    configs = [
        (8000, False, False),
        (8000, True, False),
        (8000, True, True),
    ]
    if "--native" in sys.argv:
        configs += [
            (70686, True, False),
            (70686, True, True),
        ]
    for nf, of, rad in configs:
        probe(nf, of, rad)
