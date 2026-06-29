"""Generate + commit the reference fixture for the DCM node-face contact parity test.

Two interpenetrating icosahedral "cells" (closed 20-face shells with consistent
outward normals) whose overlap region produces BOTH penetrating (repulsion, sign<0)
and nearby-outside (adhesion, sign>0) node-vs-face pairs spanning both the soft
(min_d<c_adh/2) and hard (min_d≥c_adh/2) adhesion branches — asserted below. A
per-cell cadherin multiplier exercises the √(cad·cad) adhesion weighting.

The per-node reference force is the COMMITTED numpy ground truth
``dcm_face_contact.node_face_contact_forces`` (brute-force O(N·M)) — the
tolerance-parity TARGET for ``tests/dcm/test_dcm_contact_warp_parity.py``
(anti-drift guard-rail 2: grade against the committed reference, not one authored
here). The contact force is on a FIXED mesh — no dynamic topology (Warp's
no-dynamic-topology risk lives in the remesh, a later piece). SI units.

    python ffn_sim/warp_port/fixtures/generate_dcm_contact_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.dcm_face_contact import (
    closest_point_on_triangle, node_face_contact_forces)

HERE = os.path.dirname(os.path.abspath(__file__))

R = 1.0e-6                   # cell radius (µm scale)
SEP = 0.85e-6               # ± centre offset along x → shells interpenetrate (1.7<2R)
OFF_Y = 0.3e-6              # cell-B lateral offset (breaks the axial symmetry)
SEED = 1                     # cell-B relative rotation (so facets interleave at
                             #   varying depth → pairs span rep + both adh branches)
REP, ADH = 1.0e8, 1.0e9     # ξ, ω  (from tests/test_dcm_face_contact.py)
C_REP = C_ADH = 5.0e-7
CAD = np.array([0.8, 1.2])   # per-cell cadherin multiplier


def _icosahedron():
    """Unit icosahedron: 12 vertices (on the unit sphere), 20 outward-wound faces."""
    t = (1.0 + 5.0 ** 0.5) / 2.0
    v = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ], dtype=np.float64)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    f = np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ], dtype=np.int64)
    return v, f


def main() -> None:
    v, f = _icosahedron()
    nv = v.shape[0]
    # cell B is rotated + laterally offset so its facets meet cell A at a RANGE of
    # depths → contact pairs span repulsion + both adhesion branches (see audit).
    rng = np.random.default_rng(SEED)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    pos = np.vstack([v * R + [-SEP, 0.0, 0.0],
                     (v @ q.T) * R + [SEP, OFF_Y, 0.0]])
    cof = np.array([0] * nv + [1] * nv, dtype=np.int64)
    faces = np.vstack([f, f + nv]).astype(np.int64)
    fcell = np.array([0] * f.shape[0] + [1] * f.shape[0], dtype=np.int64)

    F = node_face_contact_forces(
        pos, cof, faces, fcell, rep_strength=REP, adh_strength=ADH,
        c_rep=C_REP, c_adh=C_ADH, cad_mult=CAD)

    # --- branch-coverage audit (so the fixture provably exercises the whole law) ---
    n_rep = n_adh_hard = n_adh_soft = 0
    half = 0.5 * C_ADH
    for ni in range(pos.shape[0]):
        c1 = int(cof[ni])
        for fi in range(faces.shape[0]):
            if int(fcell[fi]) == c1:
                continue
            a, b, c = pos[faces[fi]]
            cpa, _ = closest_point_on_triangle(pos[ni], a, b, c)
            r = pos[ni] - cpa
            d = float(np.linalg.norm(r))
            fn = np.cross(b - a, c - a)
            nn = float(np.linalg.norm(fn))
            s = float(r @ fn) / nn if nn > 0 else 0.0
            if s < 0 and d < C_REP:
                n_rep += 1
            elif s > 0 and d < C_ADH:
                n_adh_hard += 1 if d >= half else 0
                n_adh_soft += 1 if d < half else 0
    assert n_rep > 0 and n_adh_hard > 0 and n_adh_soft > 0, (
        f"fixture does not cover all branches: rep={n_rep} "
        f"adh_hard={n_adh_hard} adh_soft={n_adh_soft}")

    out = os.path.join(HERE, "dcm_contact_ref.npz")
    np.savez(
        out,
        N=pos.shape[0], M=faces.shape[0],
        rep_strength=REP, adh_strength=ADH, c_rep=C_REP, c_adh=C_ADH,
        pos=pos, cell_of_node=cof, faces=faces, face_cell=fcell, cad_mult=CAD,
        ref_force=F,
    )
    print(f"wrote {out}  N={pos.shape[0]} M={faces.shape[0]}  "
          f"max|F|={np.abs(F).max():.4e}  beads_with_force={int((np.abs(F).sum(1)>0).sum())}  "
          f"branches rep={n_rep} adh_hard={n_adh_hard} adh_soft={n_adh_soft}")


if __name__ == "__main__":
    main()
