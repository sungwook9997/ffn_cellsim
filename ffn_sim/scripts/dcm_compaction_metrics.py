"""Compaction metrics from a DCM frames npz: Rg, asphericity, porosity, V/V0, per frame.

Rg = radius of gyration of all nodes (compaction: lower = tighter).
asph = asphericity of the whole aggregate (gyration-tensor eigenvalue anisotropy).
porosity = 1 - Σ(cell volumes)/V_hull(all nodes)  (void fraction; compaction eliminates voids).
V/V0 = saved per-frame mean cell volume ratio.
"""
import sys
import numpy as np
from scipy.spatial import ConvexHull


def cell_volumes(pos, faces, cof):
    """Per-cell volume via divergence theorem V_c = (1/6)Σ v0·(v1×v2) over the cell's faces."""
    fc = cof[faces[:, 0]]                       # cell-of-face (from first node)
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    tv = np.einsum("ij,ij->i", v0, np.cross(v1, v2)) / 6.0
    ncell = int(cof.max()) + 1
    V = np.zeros(ncell)
    np.add.at(V, fc, tv)
    return np.abs(V)


def analyze(npz_path):
    z = np.load(npz_path)
    frames = z["frames"]; faces = z["faces"].astype(int); cof = z["cof"].astype(int)
    vv0 = z["vv0"] if "vv0" in z.files else np.full(len(frames), np.nan)
    steps = z["step"] if "step" in z.files else np.arange(len(frames))
    out = []
    for t in range(len(frames)):
        p = frames[t].astype(float)
        c = p.mean(0)
        d = p - c
        Rg = float(np.sqrt((d**2).sum(1).mean()))
        # asphericity from gyration tensor
        G = (d.T @ d) / len(p)
        ev = np.sort(np.linalg.eigvalsh(G))[::-1]
        asph = float((ev[0] - 0.5 * (ev[1] + ev[2])) / ev.sum()) if ev.sum() > 0 else 0.0
        try:
            Vhull = float(ConvexHull(p).volume)
            Vcells = float(cell_volumes(p, faces, cof).sum())
            poro = 1.0 - Vcells / Vhull
        except Exception:
            poro = np.nan
        out.append((int(steps[t]), Rg, asph, poro, float(vv0[t])))
    return out


if __name__ == "__main__":
    print(f"{'label':>14} {'step':>7} {'Rg':>7} {'asph':>7} {'porosity':>9} {'V/V0':>6}")
    for path in sys.argv[1:]:
        label = path.split("/")[-1].replace(".npz", "")
        rows = analyze(path)
        for (s, rg, a, po, vv) in [rows[0], rows[-1]]:
            print(f"{label:>14} {s:>7} {rg:7.2f} {a:7.4f} {po:9.4f} {vv:6.3f}")
        r0, rN = rows[0], rows[-1]
        print(f"{'  Δ '+label:>14} {'':>7} {100*(rN[1]-r0[1])/r0[1]:+6.2f}% "
              f"{'asph '+format(rN[2],'.4f'):>13} poro {r0[3]:.3f}->{rN[3]:.3f}")
