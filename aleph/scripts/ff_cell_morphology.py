"""Cell MORPHOLOGY meshes for the FF viewer — membrane surface + filopodia FINGERS + nucleus.

PI 2026-07-02 ("이건 필로포디움이 아니라 그냥 작대기"): a filopodium is a membrane-wrapped actin FINGER, and
cells extend SEVERAL — not a single bare actin centerline. This builds triangle-mesh surfaces (the cell's
morphology, what the DCM mesh viewer shows) so filopodia render as fingers: a semi-transparent plasma-membrane
sphere, capped membrane TUBES for each filopodium (radius ~0.1–0.2 µm, the real diameter; Mattila-Lappalainen
2008), and the nucleus. The actin bundle cores can be overlaid as lines inside the fingers.

Meshes plug into ff_viewer_html.build_viewer as ``{"kind":"mesh","verts":(N,3),"faces":(M,3),"opacity":…}``.
"""
from __future__ import annotations

import numpy as np


def _frame(axis):
    """Orthonormal frame (axis, e1, e2)."""
    a = axis / (np.linalg.norm(axis) + 1e-12)
    t = np.array([0.0, 1.0, 0.0]) if abs(a[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, t); e1 /= np.linalg.norm(e1) + 1e-12
    e2 = np.cross(a, e1)
    return a, e1, e2


def uv_sphere(center, R, nu=28, nv=18):
    """UV sphere mesh (verts, faces)."""
    verts, faces = [], []
    for j in range(nv + 1):
        v = np.pi * j / nv
        for i in range(nu):
            u = 2 * np.pi * i / nu
            verts.append(center + R * np.array([np.sin(v) * np.cos(u), np.sin(v) * np.sin(u), np.cos(v)]))
    def idx(j, i): return j * nu + (i % nu)
    for j in range(nv):
        for i in range(nu):
            faces.append([idx(j, i), idx(j, i + 1), idx(j + 1, i + 1)])
            faces.append([idx(j, i), idx(j + 1, i + 1), idx(j + 1, i)])
    return np.array(verts, np.float32), np.array(faces, np.uint32)


def finger_tube(base, axis, length, radius, n_ring=16, n_cyl=8, n_cap=6):
    """A filopodium finger: a capped membrane cylinder from ``base`` along ``axis`` (verts, faces)."""
    a, e1, e2 = _frame(np.asarray(axis, float))
    base = np.asarray(base, float)
    verts, rings = [], []

    def ring(center, r):
        idxs = []
        for k in range(n_ring):
            th = 2 * np.pi * k / n_ring
            idxs.append(len(verts)); verts.append(center + r * (np.cos(th) * e1 + np.sin(th) * e2))
        return idxs

    for i in range(n_cyl + 1):                                    # cylinder shaft
        rings.append(ring(base + (length * i / n_cyl) * a, radius))
    for j in range(1, n_cap + 1):                                 # hemispherical tip
        phi = (np.pi / 2) * j / n_cap
        rings.append(ring(base + (length + radius * np.sin(phi)) * a, max(radius * np.cos(phi), 1e-4)))
    tip = len(verts); verts.append(base + (length + radius) * a)

    faces = []
    for i in range(len(rings) - 1):
        r0, r1 = rings[i], rings[i + 1]
        for k in range(n_ring):
            k2 = (k + 1) % n_ring
            faces.append([r0[k], r0[k2], r1[k2]]); faces.append([r0[k], r1[k2], r1[k]])
    last = rings[-1]
    for k in range(n_ring):
        faces.append([last[k], last[(k + 1) % n_ring], tip])
    return np.array(verts, np.float32), np.array(faces, np.uint32)


def merge(meshes):
    """Concatenate (verts, faces) meshes with index offsets → one (verts, faces)."""
    V, F, off = [], [], 0
    for v, f in meshes:
        V.append(v); F.append(np.asarray(f) + off); off += v.shape[0]
    return (np.concatenate(V, 0).astype(np.float32),
            np.concatenate(F, 0).astype(np.uint32) if F else np.zeros((0, 3), np.uint32))


def build_cell_mesh(center, R_cell, filopodia, *, R_nuc=None):
    """Assemble the ff_viewer scene layers for a cell morphology.

    Args:
        center: (3,) cell centroid [µm].
        R_cell: plasma-membrane radius [µm].
        filopodia: list of dicts ``{"base","axis","length","radius"}`` (µm) — the finger protrusions.
        R_nuc: nucleus radius [µm] (default 0.4·R_cell).

    Returns:
        list of ff_viewer layer dicts (membrane, filopodia, nucleus meshes).
    """
    center = np.asarray(center, float)
    R_nuc = R_nuc if R_nuc is not None else 0.4 * R_cell
    mem_v, mem_f = uv_sphere(center, R_cell)
    fingers = merge([finger_tube(f["base"], f["axis"], f["length"], f["radius"]) for f in filopodia]) \
        if filopodia else (np.zeros((0, 3), np.float32), np.zeros((0, 3), np.uint32))
    nuc_v, nuc_f = uv_sphere(center, R_nuc, nu=22, nv=14)
    return [
        {"name": "plasma membrane", "kind": "mesh", "color": "#7db8e8", "opacity": 0.30,
         "verts": mem_v, "faces": mem_f},
        {"name": f"filopodia ({len(filopodia)} fingers)", "kind": "mesh", "color": "#ffb000", "opacity": 1.0,
         "verts": fingers[0], "faces": fingers[1]},
        {"name": "nucleus", "kind": "mesh", "color": "#c65b7c", "opacity": 0.95,
         "verts": nuc_v, "faces": nuc_f},
    ]


def sample_surface_filopodia(center, R_cell, n=10, *, lengths=(1.5, 3.5), radius=0.15,
                             lead_axis=(1, 0, 0), spread_deg=55.0, rng=None):
    """Place ``n`` filopodia emerging from the membrane, clustered toward ``lead_axis`` (a leading edge)."""
    rng = rng or np.random.default_rng(0)
    lead = np.asarray(lead_axis, float); lead /= np.linalg.norm(lead)
    out = []
    for _ in range(n):
        # direction within spread_deg of the lead axis (a filopodial leading-edge fan)
        d = lead + np.tan(np.deg2rad(spread_deg)) * rng.standard_normal(3)
        d /= np.linalg.norm(d)
        base = center + (R_cell - 0.15) * d
        out.append({"base": base, "axis": d, "length": float(rng.uniform(*lengths)), "radius": radius})
    return out


if __name__ == "__main__":                                       # synthetic sanity render
    from aleph.scripts.ff_viewer_html import build_viewer
    c = np.zeros(3); R = 7.5
    filo = sample_surface_filopodia(c, R, n=12, lead_axis=(1, 0, 0), rng=np.random.default_rng(3))
    layers = build_cell_mesh(c, R, filo)
    build_viewer({"MCF7 cell with filopodia (synthetic)": layers},
                 out="aleph/outputs/ff/figs/_cell_morphology_test.html", title="FF cell morphology test")
    # a matplotlib 3D check too
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(7, 7)); ax = fig.add_subplot(111, projection="3d")
    for L in layers:
        v, f = np.asarray(L["verts"]), np.asarray(L["faces"])
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        tri = v[f]; pc = Poly3DCollection(tri, alpha=L["opacity"], facecolor=L["color"], edgecolor="none")
        ax.add_collection3d(pc)
    ax.set_xlim(-11, 11); ax.set_ylim(-11, 11); ax.set_zlim(-11, 11); ax.set_box_aspect((1, 1, 1))
    ax.set_title("MCF7 cell + 12 filopodia (synthetic morphology)")
    fig.savefig("aleph/outputs/ff/figs/_cell_morphology_test.png", dpi=120)
    print("wrote _cell_morphology_test.{html,png}")
