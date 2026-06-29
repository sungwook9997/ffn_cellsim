"""Quantify whether a DCM aggregate spheroid is GENUINELY cohesive (PI verify gate).

The in-run contact metric uses c_adh=5um (=0.67 R) and over-reports — node-NODE
aggregates score ~0.88 there yet have ZERO true surface contact (cells float ~4um
apart). This checker computes contact at PHYSICAL thresholds (down to 0.5um), per-cell
turgor V/V0 (collapse/inflation), isolated-cell count, and packing fraction, from a
saved spheroid .npy (+ the cell topology). Genuine cohesion = contact survives <=2um.

Usage:  python -m ffn_sim.scripts.dcm_agg_check SPHEROID.npy [--nv 42]
"""
from __future__ import annotations

import argparse

import numpy as np
from scipy.spatial import cKDTree, ConvexHull

UM = 1e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npy")
    ap.add_argument("--nv", type=int, default=42)
    ap.add_argument("--R", type=float, default=7.5e-6)
    a = ap.parse_args()
    pos = np.load(a.npy)
    nv, R = a.nv, a.R
    nc = pos.shape[0] // nv
    cof = np.repeat(np.arange(nc), nv)

    tree = cKDTree(pos)
    d, idx = tree.query(pos, k=16)
    nf = []
    for i in range(pos.shape[0]):
        for j, dist in zip(idx[i, 1:], d[i, 1:]):
            if cof[j] != cof[i]:
                nf.append(dist); break
    nf = np.array(nf) * UM

    # per-cell turgor volume (divergence theorem) needs topology; use the canonical
    # icosphere tris for nv (rebuild from dcm.icosphere_mesh to match the run).
    from ffn_sim.archive.hoomd_legacy.cell.dcm import icosphere_mesh
    sub = 1 if nv == 42 else 2
    _v, _e, tris0 = icosphere_mesh(R, sub)
    V0 = (4.0 / 3.0) * np.pi * R ** 3

    def cellV(v):
        aa, bb, cc = v[tris0[:, 0]], v[tris0[:, 1]], v[tris0[:, 2]]
        return abs(np.einsum("ij,ij->i", aa, np.cross(bb, cc)).sum()) / 6.0
    vols = np.array([cellV(pos[c * nv:(c + 1) * nv]) / V0 for c in range(nc)])
    cents = np.array([pos[c * nv:(c + 1) * nv].mean(0) for c in range(nc)])
    dd, _ = cKDTree(cents).query(cents, k=2)
    nn = dd[:, 1] / (2 * R)
    pack = (vols * V0).sum() / ConvexHull(pos).volume

    print(f"AGGREGATE CHECK  {a.npy}  (N={nc} cells, nv={nv})")
    print(f"  GENUINE contact vs threshold (node-NODE collapses to 0 below 3um):")
    for t in (5, 4, 3, 2, 1, 0.5):
        print(f"     {t:>4}um -> {(nf < t).mean():.3f}")
    genuine = (nf < 2).mean()
    print(f"  nearest foreign-node gap [um]: p0={nf.min():.2f} median={np.median(nf):.2f} mean={nf.mean():.2f}")
    print(f"  per-cell V/V0: mean={vols.mean():.3f} min={vols.min():.3f} max={vols.max():.3f}")
    print(f"  collapsed cells (V/V0<0.5): {(vols < 0.5).sum()}   over-inflated (>1.3): {(vols > 1.3).sum()}")
    print(f"  packing fraction = {pack:.3f}   isolated cells (NN>1.5*2R): {int((nn > 1.5).sum())}/{nc}")
    verdict = ("GENUINE cohesion (surfaces touch)" if genuine > 0.3
               else "FAKE cohesion (cells float apart)")
    print(f"  => contact@2um={genuine:.3f} -> {verdict}")


if __name__ == "__main__":
    main()
