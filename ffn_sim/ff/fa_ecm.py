r"""FF focal-adhesion ↔ ECM-fiber coupling (Warp, µm·pN·s) — the clutch binds a LIVE collagen-I (Mikado) fiber
node instead of a FIXED substrate point, so traction REMODELS the matrix (Newton reaction on the fiber).

The crawl/adhesion driver anchors each FA clutch to a fixed dish point (`fa_clutch_warp.clutch_spring_kernel`,
one-sided — the reaction goes to the immovable substrate). The PI's actual matrix is collagen-I: a disordered
Mikado fiber network (`ecm_mikado`) whose nodes are live DOFs. This module is the S4 primitive: a TWO-sided
clutch spring between an actin end node (in the cell's position array) and a collagen fiber node (in the ECM's
position array). The same Hookean law `f = k_int·(L−rest)/L·(ECM−actin)` now applies **+f to the actin (traction
on the cell) and −f to the fiber node (the Newton reaction)** — so the cell pulls its collagen substrate toward
it, recruiting/aligning fibers (the emergent 1/r remodel validated for DCM⊗Kim). No lumped traction law: the
force is the fine-grained clutch spring, and matrix remodelling is the reaction, molecule by molecule.

The cell and the ECM are SEPARATE position/force arrays (distinct networks, each with its own bending +
crosslink + drag dynamics); the clutch is the only mechanical bridge. `attach_clutches_to_ecm` forms the
nascent bonds (each clutch grabs the nearest fiber node within a capture radius); `clutch_ecm_spring_kernel`
applies the coupled force each step. Both are additive — a cell with no ECM (`ecm_node<0` everywhere) is
untouched, and the fixed-substrate path is unchanged.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from scipy.spatial import cKDTree


@wp.kernel
def clutch_ecm_spring_kernel(
    cell_pos: wp.array(dtype=wp.vec3d),          # cell node positions [µm]
    actin: wp.array(dtype=wp.int32),             # (M,) actin end-node index per clutch (into cell_pos)
    ecm_pos: wp.array(dtype=wp.vec3d),           # ECM (collagen) node positions [µm]
    ecm_node: wp.array(dtype=wp.int32),          # (M,) bound ECM node index per clutch (into ecm_pos); <0 = unbound
    k_int: wp.float64,                           # clutch stiffness [pN/µm]
    rest: wp.float64,                            # clutch rest length [µm]
    cell_force: wp.array(dtype=wp.vec3d),        # += traction on the cell
    ecm_force: wp.array(dtype=wp.vec3d),         # += Newton reaction on the collagen fiber (remodels the matrix)
):
    """Two-sided Hookean clutch spring actin↔collagen-node: f = k_int·(L−rest)/L·(ECM−actin). +f on the actin
    (traction), −f on the fiber node (reaction). Unbound clutch (ecm_node<0) exerts nothing on either side."""
    t = wp.tid()
    j = ecm_node[t]
    if j < 0:
        return
    a = actin[t]
    d = ecm_pos[j] - cell_pos[a]
    L = wp.length(d)
    if L > wp.float64(1e-12):
        f = (k_int * (L - rest) / L) * d
        wp.atomic_add(cell_force, a, f)          # traction on the cell (toward the fiber)
        wp.atomic_add(ecm_force, j, -f)          # Newton pair: pull the collagen fiber toward the cell


def attach_clutches_to_ecm(actin_pos: np.ndarray, ecm_pos: np.ndarray, capture_um: float) -> np.ndarray:
    """Form nascent FA↔collagen bonds: each clutch (at ``actin_pos[i]``) binds the NEAREST ECM node within
    ``capture_um``; returns the per-clutch ECM node index (``-1`` = no fiber in reach). Host, at attach time
    (occasional), via a KD-tree — cheap relative to the per-step force kernel."""
    actin_pos = np.ascontiguousarray(actin_pos, np.float64)
    ecm_pos = np.ascontiguousarray(ecm_pos, np.float64)
    if ecm_pos.shape[0] == 0 or actin_pos.shape[0] == 0:
        return np.full(actin_pos.shape[0], -1, np.int64)
    dist, idx = cKDTree(ecm_pos).query(actin_pos, k=1)
    return np.where(dist <= float(capture_um), idx, -1).astype(np.int64)


__all__ = ["clutch_ecm_spring_kernel", "attach_clutches_to_ecm"]
