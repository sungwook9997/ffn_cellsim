"""Per-node virial (Cauchy) stress for the DCM — the REAL stress field.

Fixes fidelity red-flag #1 (2026-07-07 audit): the viewer's ``fmag`` channel is the |net RESIDUAL
force| per node (∝ node velocity, →0 at convergence), NOT a stress — and it does not concentrate at
junctions. This module assembles the virial/Cauchy tensor from the ACTUAL pairwise force laws:

    σ_node = (1/V_node) · Σ_{pairs on node} ½·(r_ij ⊗ f_ij)  +  (turgor isotropic −ΔP·I)

for the dominant stress carriers — cadherin trans-dimer bonds (deviatoric, concentrates at
junctions), cortex edge springs (deviatoric, membrane tension), and turgor (hydrostatic). Returns
the von-Mises invariant (σ_vm, the deviatoric magnitude a real stress map shows) and the hydrostatic
pressure (p = −tr σ / 3) per node, in Pa.

Force laws (must match the sim kernels exactly):
  cadherin (``cadherin_bond_force_kernel``): F = k_trans_eff·max(0, L−r0), k_trans_eff = k_trans·bundle_n.
  cortex edge (``_bond_accumulate``):        F = k_edge·(L−r0_edge), Hookean along the edge.
  turgor (``dcm_turgor_force_kernel``):       ΔP_c = dP0 + K_vol·(V0−V_c)/V0, isotropic per cell.

NOT included (smaller / localised): bending, nucleus, node-face contact. So this is the dominant-term
virial, honest as such — surface-tension γ enters via the edge springs' rest-length tension, not as a
separate area term. Callers should label the field "virial stress (cadherin+cortex+turgor)".
"""
from __future__ import annotations

import numpy as np


def _accumulate_pair_virial(sig, fnode, i, j, pos, fvec):
    """Add ½·(r_ij ⊗ f_ij) to both endpoints' node tensors, and f to the net-force check array."""
    rij = pos[i] - pos[j]                                  # (P,3)
    # outer product r_ij ⊗ f_ij per pair → (P,3,3)
    outer = 0.5 * (rij[:, :, None] * fvec[:, None, :])
    for a in range(3):
        for b in range(3):
            np.add.at(sig[:, a, b], i, outer[:, a, b])
            np.add.at(sig[:, a, b], j, outer[:, a, b])     # symmetric split, same contribution
    np.add.at(fnode, i, fvec)
    np.add.at(fnode, j, -fvec)


def node_virial_stress(pos, cof, node_vol,
                       bonds=None, k_trans_eff=0.0, r0_cad=0.0,
                       edges=None, k_edge=0.0, r0_edge=None,
                       dP_cell=None):
    """Per-node virial stress → (sigma_vm, pressure) in Pa. `pos` in metres.

    Args:
      pos:      (N,3) node positions [m].
      cof:      (N,) compact cell id per node (<0 = dormant, excluded).
      node_vol: (N,) representative volume per node [m³] (= V_cell/npc).
      bonds:    (Mc,2) cadherin bonded node-index pairs, or None.
      k_trans_eff: cadherin effective stiffness = k_trans·bundle_n [N/m].
      r0_cad:   cadherin rest length [m].
      edges:    (Me,2) cortex edge node-index pairs, or None.
      k_edge:   cortex edge stiffness [N/m].
      r0_edge:  (Me,) cortex edge rest lengths [m].
      dP_cell:  (n_cells,) per-cell turgor pressure ΔP [Pa], or None.

    Returns:
      (sigma_vm, pressure): each (N,) float64 [Pa]. Dormant nodes → 0.
    """
    N = pos.shape[0]
    sig = np.zeros((N, 3, 3), np.float64)
    fnode = np.zeros((N, 3), np.float64)                   # reconstructed net force (validation)

    # --- cadherin bonds: F = k_trans_eff·max(0, L−r0) attractive along the bond ---
    if bonds is not None and len(bonds):
        i, j = bonds[:, 0], bonds[:, 1]
        rij = pos[i] - pos[j]
        L = np.linalg.norm(rij, axis=1)
        u = rij / np.maximum(L, 1e-30)[:, None]
        Fmag = k_trans_eff * np.maximum(0.0, L - r0_cad)   # >0 → pulls i toward j
        fvec = -Fmag[:, None] * u                          # force on i
        _accumulate_pair_virial(sig, fnode, i, j, pos, fvec)

    # --- cortex edges: F = k_edge·(L−r0_edge) Hookean along the edge ---
    if edges is not None and len(edges) and k_edge > 0.0 and r0_edge is not None:
        i, j = edges[:, 0], edges[:, 1]
        rij = pos[i] - pos[j]
        L = np.linalg.norm(rij, axis=1)
        u = rij / np.maximum(L, 1e-30)[:, None]
        Fmag = k_edge * (L - r0_edge)                      # >0 stretched → pulls together
        fvec = -Fmag[:, None] * u
        _accumulate_pair_virial(sig, fnode, i, j, pos, fvec)

    # normalise the pairwise virial by node volume → stress [Pa]
    live = cof >= 0
    vol = np.maximum(node_vol, 1e-30)
    sig /= vol[:, None, None]

    # --- turgor: isotropic σ += −ΔP·I per cell (hydrostatic, adds to pressure, not von Mises) ---
    if dP_cell is not None:
        dPn = np.zeros(N)
        dPn[live] = dP_cell[cof[live]]
        for a in range(3):
            sig[:, a, a] += -dPn                            # compressive interior pressure

    # invariants
    sxx, syy, szz = sig[:, 0, 0], sig[:, 1, 1], sig[:, 2, 2]
    sxy, syz, szx = sig[:, 0, 1], sig[:, 1, 2], sig[:, 2, 0]
    vm = np.sqrt(np.maximum(0.0, 0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                            + 3.0 * (sxy ** 2 + syz ** 2 + szx ** 2)))
    pressure = -(sxx + syy + szz) / 3.0
    vm[~live] = 0.0
    pressure[~live] = 0.0
    return vm, pressure, fnode
