"""Per-node virial (Cauchy) von-Mises STRESS + areal STRAIN for the FF cortex — the FEM-style field.

Stress: the virial/Cauchy tensor assembled from the ACTUAL FF cortex force laws (same construction as
``dcm/dcm_virial_stress.py``, red-flag-#1 fix), so it is a REAL internal stress (nonzero at equilibrium),
NOT the net residual force (which →0 at convergence):

    σ_node = (1/V_node)·Σ_{bonds on node} ½·(r_ij ⊗ f_ij)  +  turgor hydrostatic (−ΔP·I)

Deviatoric carriers (→ von-Mises σ_vm): cortex crosslinks (α-actinin/filamin, per-bond k·(L−r0)) + myosin
minifilament links (contractile f_myo). Hydrostatic carrier (→ pressure p): osmotic turgor ΔP. Bending is a
3-body term and (like the DCM dominant-term virial) is excluded — label the field accordingly. FF units are
µm·pN·s, and **pN/µm² = Pa exactly** (1e-12 N / 1e-12 m²), so the virial comes out in Pa with no conversion.

Strain: per-node areal strain from the fixed cortex surface triangulation — each face's current area vs its
frame-0 rest area, averaged to its 3 nodes. Pure geometry (no force model): shows where the cortex surface is
locally stretched (>0) or compressed (<0) under indentation / crawl.
"""
from __future__ import annotations

import numpy as np


def _accum_bond_virial(sig, i, j, pos, fvec):
    """Add ½·(r_ij ⊗ f_ij) to both endpoints (symmetric split)."""
    rij = pos[i] - pos[j]
    outer = 0.5 * (rij[:, :, None] * fvec[:, None, :])       # (M,3,3)
    for a in range(3):
        for b in range(3):
            np.add.at(sig[:, a, b], i, outer[:, a, b])
            np.add.at(sig[:, a, b], j, outer[:, a, b])


def cortex_node_stress(pos, node_vol, *, xl_ij=None, k_xl=None, r0_xl=None,
                       myo_ij=None, f_myo=0.0, dP=0.0):
    """Per-node (σ_vm, pressure) in Pa for the cortex nodes ``pos`` (µm). ``node_vol`` (µm³) per node.

    xl_ij (E,2)/k_xl (E,)/r0_xl (E,): crosslink bonds, Hookean F=k·(L−r0). myo_ij (My,2)/f_myo: myosin links,
    constant contractile f_myo. dP [Pa]: osmotic turgor → isotropic −ΔP·I (hydrostatic → pressure, not σ_vm)."""
    N = pos.shape[0]
    sig = np.zeros((N, 3, 3), np.float64)
    if xl_ij is not None and len(xl_ij):
        i, j = xl_ij[:, 0], xl_ij[:, 1]
        rij = pos[i] - pos[j]; L = np.linalg.norm(rij, axis=1)
        u = rij / np.maximum(L, 1e-30)[:, None]
        Fmag = np.asarray(k_xl) * (L - np.asarray(r0_xl))     # >0 stretched → pulls together
        _accum_bond_virial(sig, i, j, pos, -Fmag[:, None] * u)
    if myo_ij is not None and len(myo_ij) and f_myo != 0.0:
        i, j = myo_ij[:, 0], myo_ij[:, 1]
        rij = pos[i] - pos[j]; L = np.linalg.norm(rij, axis=1)
        u = rij / np.maximum(L, 1e-30)[:, None]
        _accum_bond_virial(sig, i, j, pos, -float(f_myo) * u)  # contractile (pulls together)
    sig /= np.maximum(node_vol, 1e-30)[:, None, None]
    if dP != 0.0:
        for a in range(3):
            sig[:, a, a] += -float(dP)                        # interior turgor pressure (hydrostatic)
    sxx, syy, szz = sig[:, 0, 0], sig[:, 1, 1], sig[:, 2, 2]
    sxy, syz, szx = sig[:, 0, 1], sig[:, 1, 2], sig[:, 2, 0]
    vm = np.sqrt(np.maximum(0.0, 0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                            + 3.0 * (sxy ** 2 + syz ** 2 + szx ** 2)))
    pressure = -(sxx + syy + szz) / 3.0
    return vm, pressure


def cortex_node_areal_strain(pos, faces, area0):
    """Per-node areal strain (current face area / rest area0 − 1), averaged to the 3 nodes of each face.
    ``faces`` (F,3) fixed triangulation; ``area0`` (F,) frame-0 rest areas. Returns (N,)."""
    N = pos.shape[0]
    a = pos[faces[:, 0]]; b = pos[faces[:, 1]]; c = pos[faces[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    fstrain = area / np.maximum(area0, 1e-30) - 1.0
    acc = np.zeros(N); cnt = np.zeros(N)
    for k in range(3):
        np.add.at(acc, faces[:, k], fstrain); np.add.at(cnt, faces[:, k], 1.0)
    return acc / np.maximum(cnt, 1.0)


def face_areas(pos, faces):
    """Triangle areas (F,) for the fixed triangulation — call on frame 0 to get rest areas."""
    a = pos[faces[:, 0]]; b = pos[faces[:, 1]]; c = pos[faces[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


__all__ = ["cortex_node_stress", "cortex_node_areal_strain", "face_areas"]
