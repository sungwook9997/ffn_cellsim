r"""Intermediate-filament (IF) perinuclear cage compartment for the Active Cell (Warp-CUDA).

Reuses the ff-audited :mod:`ff.intermediate_filaments` radial-spoke cage — the nucleus↔cortex load path that
fixes the mechanically-DECOUPLED nucleus (INTERNAL_DISPLACEMENT_DIAGNOSIS 2026-07-16) — and wraps it into an
ac/ §1.4 force primitive with the ``accumulate(pos, f)`` contract. Each spoke runs radially from just outside
the nucleus to just inside the cortex; coupling bonds are LINC (spoke inner bead ↔ nucleus) + cortex anchor
(spoke outer bead ↔ cortex) + tangential plectin crosslinks. In the LINEAR regime every IF bond is Hookean,
so the compartment rides the reused ``ff.network_warp.link_spring_kernel`` — no new kernel.

The physiological contrast is the EMT switch: **keratin (K8/K18, epithelial/MCF7)** vs **vimentin
(mesenchymal/MDA-MB-231)**. 2026-07-22 DOI sourcing: keratin single-filament l_p ≈ 0.30 µm (Lichtenstern 2012,
`10.1016/j.jsb.2011.11.003`), vimentin ≈ 0.50 µm (Lin 2010, `10.1016/j.jmb.2010.04.054`); network stretch
modulus ~6–9 MPa; keratin loss softens epithelial cells ~30–60% + drives invasion (Seltmann 2013,
`10.1073/pnas.1310493110`), vimentin doubles cytoplasmic G′ and forms the large-strain nuclear-protective
safety net (Guo 2013 `10.1016/j.bpj.2013.08.037`; Patteson 2019 `10.1083/jcb.201902046`). Only the effective
Young's modulus ``E_if_Pa`` differs between the two presets; the geometry/kernel are identical.

Sanity Gate:
    * Dimensions: E [Pa=pN/µm²], A [µm²], k = E·A/l_seg [pN/µm]; force [pN].
    * Boundary: every bond rest length = its construction distance ⇒ ZERO force at build (physiological
      baseline; the cage is a passenger until the cell deforms).
    * Coupling (the whole point): displacing the cortex anchors transmits a nonzero force to the nucleus beads
      through the spokes — the isolated-cage test asserts the nucleus feels the cortex, unfreezing it.
    * Residency: bond arrays are CUDA; ``accumulate`` launches one Hookean kernel, no host readback.

DEFERRED (behind PI-authored KB registration, never faked with a magic-number fit): the nonlinear 3-regime
strain-stiffening (soft coil → α→β unfolding plateau → re-stiffening; Kreplak 2005 `10.1016/j.jmb.2005.09.092`,
Block 2017 `10.1103/PhysRevLett.118.048101`) — the ``ff.wlc`` extensible-WLC law is the ready primitive for
it. GAP magnitudes (native cage count ``n_fil``, k_linc/k_anchor/crosslink ratio, IF cage mesh size) are
provisional PI GAPs — a weak coupling at the KB-anchored k_bb is a FINDING to surface, not a knob to tune.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.intermediate_filaments import build_if_cage, resolve_intermediate_filaments
from aleph.laws.network_warp import link_spring_kernel

__all__ = ["IntermediateFilamentCage", "build_if_compartment", "IF_E_KERATIN_PA", "IF_E_VIMENTIN_PA"]

# Effective network Young's modulus presets (2026-07-22 sourcing, ~6–9 MPa network band; the EMT split is
# carried by E only). Keratin (MCF7) is the more flexible/softer single filament; vimentin (MDA-231) the
# large-strain-tough one. Values are provisional network-modulus anchors → PI ratification pending.
IF_E_KERATIN_PA = 6.0e6        # Pa — keratin K8/K18 network modulus (epithelial / MCF7)
IF_E_VIMENTIN_PA = 9.0e6       # Pa — vimentin network modulus (mesenchymal / MDA-MB-231)


@dataclass
class IntermediateFilamentCage:
    """IF-cage force primitive — ``accumulate(pos, f)`` adds the reused Hookean bond force of the cage."""

    device: str
    cell_type: str            # "keratin" (MCF7) | "vimentin" (MDA-231)
    n_if: int
    pos0: npt.NDArray[np.float64]   # (n_if,3) spoke bead reference positions [µm] (LOCAL, for isolated runs)
    bonds_d: wp.array          # (Nb,2) int32 GLOBAL bond endpoints
    k_d: wp.array              # (Nb,) float64 stiffness [pN/µm]
    rest_d: wp.array           # (Nb,) float64 rest length [µm]
    n_backbone: int
    n_linc: int
    n_anchor: int
    n_xl: int
    E_if_Pa: float
    n_bonds: int = field(init=False)

    def __post_init__(self) -> None:
        self.n_bonds = int(self.bonds_d.shape[0])

    def accumulate(self, pos: wp.array, f: wp.array) -> None:
        """Add the reused Hookean IF-bond force (backbone + LINC + anchor + crosslink) to the node force."""
        if self.n_bonds:
            wp.launch(link_spring_kernel, dim=self.n_bonds,
                      inputs=[pos, self.bonds_d, self.k_d, self.rest_d], outputs=[f], device=self.device)


def build_if_compartment(
    *, centre, R_nuc_um: float, R_cortex_um: float, nuc_pos: npt.NDArray[np.float64],
    cortex_pos: npt.NDArray[np.float64], n_fil: int = 60, cell_type: str = "keratin",
    nuc_offset: int, cortex_offset: int, if_offset: int, crosslink: bool = True,
    device: str | None = None,
) -> IntermediateFilamentCage:
    """Build the radial-spoke IF cage between the nucleus and cortex, wired to global indices.

    ``cell_type`` selects the keratin(MCF7)/vimentin(MDA-231) modulus preset. ``*_offset`` place the IF beads,
    nucleus beads, and cortex nodes in the combined node array; CUDA required (I0-A).
    """
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError("IntermediateFilamentCage is Warp-CUDA-only (I0-A)")
    dev_s = str(dev)
    e_if = IF_E_VIMENTIN_PA if cell_type == "vimentin" else IF_E_KERATIN_PA
    res = resolve_intermediate_filaments(n_fil=n_fil, E_if_Pa=e_if)
    cage = build_if_cage(centre, R_nuc_um, R_cortex_um, np.asarray(nuc_pos, np.float64),
                         np.asarray(cortex_pos, np.float64), res, nuc_offset=nuc_offset,
                         cortex_offset=cortex_offset, if_offset=if_offset, crosslink=crosslink)
    bonds = np.ascontiguousarray(np.stack([cage.bond_i, cage.bond_j], axis=1), np.int32)
    with wp.ScopedDevice(dev_s):
        bonds_d = wp.array(bonds, dtype=wp.int32, device=dev_s)
        k_d = wp.array(np.ascontiguousarray(cage.bond_k, np.float64), dtype=wp.float64, device=dev_s)
        rest_d = wp.array(np.ascontiguousarray(cage.bond_rest, np.float64), dtype=wp.float64, device=dev_s)
    return IntermediateFilamentCage(
        device=dev_s, cell_type=cell_type, n_if=int(cage.n_if),
        pos0=np.ascontiguousarray(cage.pos, np.float64), bonds_d=bonds_d, k_d=k_d, rest_d=rest_d,
        n_backbone=int(cage.n_backbone), n_linc=int(cage.n_linc), n_anchor=int(cage.n_anchor),
        n_xl=int(cage.n_xl), E_if_Pa=float(e_if),
    )
