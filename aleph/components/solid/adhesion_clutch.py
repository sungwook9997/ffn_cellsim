r"""ECM adhesion-clutch compartment for the Active Cell — the substrate traction load path (Warp-CUDA).

The molecular clutch (Chan & Odde 2008, `10.1126/science.1163595`) is how the cell GRIPS the ECM and pulls:
a Hookean linkage (integrin–talin–actin) from a basal actin node to a bound collagen node, with a
force-dependent bond that de-adheres EMERGENTLY under load — never a latch. This wraps the ff-audited
`ff.fa_ecm.clutch_ecm_spring_kernel` (two-sided actin↔collagen traction + Newton reaction that remodels the
matrix) and `ff.fa_ecm.attach_clutches_to_ecm` (nascent FA↔collagen capture) into an ac/ §1.4 primitive with
the ``accumulate(cell_pos, cell_force)`` contract, and adds the α2β1–collagen off-rate as the sourced
de-adhesion kinetics.

Sourced parameters (2026-07-22 DOI dossier; PI-registration pending, not tuned):
    * **Clutch linkage stiffness k_int = 0.8 pN/nm = 800 pN/µm** (Chan-Odde 2008 / Bangasser 2013; band 0.5–5).
    * **α2β1–collagen is a SLIP bond** (no catch data): k_off(F) = k_off0·exp(F·x_β/k_BT) with
      **k_off0 = 0.44 s⁻¹, x_β = 0.7 nm** (Attwood 2013, `10.3390/ijms14022832`) ⇒ f_β = k_BT/x_β ≈ 5.9 pN.
      (The α5β1–fibronectin catch bond, Kong 2009, is a shape PROXY only — not used here.)
    * k_on (2D α2β1 on-rate) is a hard **PI GAP** — the capture is geometric (nearest fiber within a radius)
      pending a sourced on-rate; a weak grip at the sourced k_int is a FINDING, not a knob.

Sanity Gate:
    * Dimensions: k_int [pN/µm], L,rest [µm], force [pN]; off-rate [s⁻¹].
    * Boundary: a clutch at rest length exerts zero force; an unbound clutch (ecm_node<0) exerts nothing;
      an engaged, stretched clutch pulls the actin node TOWARD its collagen node (traction) with the exact
      Newton reaction on the fiber (matrix remodelling).
    * Kinetics: k_off increases monotonically with load (slip) — de-adhesion is emergent from force.
    * Residency: clutch/ECM arrays are CUDA; ``accumulate`` launches one kernel, no host readback (attach is
      an occasional host KD-tree step, out of the hot loop).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.fa_ecm import attach_clutches_to_ecm, clutch_ecm_spring_kernel
from aleph.laws.units import KBT

__all__ = [
    "AdhesionClutchCompartment",
    "build_adhesion_clutch_compartment",
    "attwood_slip_off_rate",
    "CLUTCH_K_INT_PN_UM",
    "A2B1_KOFF0_PER_S",
    "A2B1_XBETA_UM",
]

CLUTCH_K_INT_PN_UM = 800.0     # integrin-talin-actin linkage stiffness = 0.8 pN/nm (Chan-Odde 2008; band 0.5-5)
A2B1_KOFF0_PER_S = 0.44        # α2β1-collagen unloaded off-rate [1/s] (Attwood 2013, I-domain-GFOGER SMFS)
A2B1_XBETA_UM = 0.7e-3         # Bell slip distance x_β = 0.7 nm [µm] ⇒ f_β = k_BT/x_β ≈ 5.9 pN


def attwood_slip_off_rate(force_pn, k_off0: float = A2B1_KOFF0_PER_S, x_beta_um: float = A2B1_XBETA_UM):
    """α2β1–collagen Bell slip off-rate k_off(F) = k_off0·exp(F·x_β/k_BT) [1/s] (host oracle)."""
    f_beta = KBT / x_beta_um                                    # pN·µm / µm = pN ≈ 5.9 pN
    return k_off0 * np.exp(np.asarray(force_pn, np.float64) / f_beta)


@dataclass
class AdhesionClutchCompartment:
    """ECM-clutch force primitive — ``accumulate(cell_pos, cell_force)`` adds substrate traction (§1.4).

    ``ecm_pos`` are the collagen node positions (the substrate); ``actin_d`` maps each clutch to a basal actin
    node (global index into the cell array); ``ecm_node_d`` is its bound collagen node (<0 = unbound). The
    two-sided kernel also accumulates the Newton reaction into ``ecm_force`` (matrix remodelling); with a rigid
    substrate that reaction is discarded.
    """

    device: str
    n_clutch: int
    ecm_pos_d: wp.array            # (E,3) collagen node positions [µm]
    actin_d: wp.array             # (M,) int32 basal actin node index per clutch (GLOBAL)
    ecm_node_d: wp.array          # (M,) int32 bound ECM node per clutch (<0 unbound)
    k_int: float
    rest_um: float
    capture_um: float
    _ecm_force_d: wp.array = field(default=None)   # scratch reaction accumulator (rigid substrate → discarded)

    def __post_init__(self) -> None:
        if self._ecm_force_d is None:
            with wp.ScopedDevice(self.device):
                self._ecm_force_d = wp.zeros(int(self.ecm_pos_d.shape[0]), dtype=wp.vec3d, device=self.device)

    def accumulate(self, cell_pos: wp.array, cell_force: wp.array) -> None:
        """Add the two-sided Hookean clutch traction (actin ← collagen) to the cell node force."""
        if self.n_clutch:
            self._ecm_force_d.zero_()
            wp.launch(clutch_ecm_spring_kernel, dim=self.n_clutch,
                      inputs=[cell_pos, self.actin_d, self.ecm_pos_d, self.ecm_node_d,
                              wp.float64(self.k_int), wp.float64(self.rest_um)],
                      outputs=[cell_force, self._ecm_force_d], device=self.device)


def build_adhesion_clutch_compartment(
    *, actin_pos: npt.NDArray[np.float64], actin_global_idx: npt.NDArray[np.integer],
    ecm_pos: npt.NDArray[np.float64], capture_um: float = 0.1, rest_um: float = 0.0,
    k_int: float = CLUTCH_K_INT_PN_UM, device: str | None = None,
) -> AdhesionClutchCompartment:
    """Form nascent clutches from the given basal actin nodes to the nearest collagen node within ``capture_um``.

    ``actin_pos`` / ``actin_global_idx`` are the candidate basal actin nodes' positions and their GLOBAL cell
    indices; ``ecm_pos`` the collagen substrate nodes. One clutch per candidate actin node; unbound if no fiber
    is in reach. CUDA required (I0-A).
    """
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError("AdhesionClutchCompartment is Warp-CUDA-only (I0-A)")
    dev_s = str(dev)
    ecm_node = attach_clutches_to_ecm(np.asarray(actin_pos, np.float64), np.asarray(ecm_pos, np.float64),
                                      float(capture_um))
    with wp.ScopedDevice(dev_s):
        ecm_pos_d = wp.array(np.ascontiguousarray(ecm_pos, np.float64), dtype=wp.vec3d, device=dev_s)
        actin_d = wp.array(np.ascontiguousarray(actin_global_idx, np.int32), dtype=wp.int32, device=dev_s)
        ecm_node_d = wp.array(np.ascontiguousarray(ecm_node, np.int32), dtype=wp.int32, device=dev_s)
    return AdhesionClutchCompartment(
        device=dev_s, n_clutch=int(actin_global_idx.shape[0]), ecm_pos_d=ecm_pos_d, actin_d=actin_d,
        ecm_node_d=ecm_node_d, k_int=float(k_int), rest_um=float(rest_um), capture_um=float(capture_um))
