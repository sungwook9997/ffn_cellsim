r"""Intermediate-filament (IF) perinuclear cage — the nucleus↔cortex load path (FF/Warp, µm·pN·s).

The dominant mechanical coupling that transmits cortex/AFM deformation to the nucleus. The
internal-displacement diagnosis (docs/v2_audit/INTERNAL_DISPLACEMENT_DIAGNOSIS_2026-07-16.md) proved the
nucleus is mechanically DECOUPLED (RAW displacement 0.0 nm under AFM), and the spatial-Biot FSI is too weak
to move it (secondary). The IF cage is the fix: real cells transmit deformation to the nucleus through the
vimentin/keratin cytoskeletal cage + LINC, not the cytoplasm pressure.

GEOMETRY (design-verified 2026-07-16): a RADIAL-SPOKE cage — the tangential perinuclear shell of the older
plumbing does NOT couple nucleus↔cortex (it is a decoration). Each spoke is a bead chain running RADIALLY
from just outside the nucleus (R_nuc) to just inside the cortex (R_cortex≈R0) along a Fibonacci direction.
Coupling bonds: **LINC** (spoke inner bead ↔ nearest nucleus surface bead — nesprin/SUN) + **cortex anchor**
(spoke outer bead ↔ nearest cortex node — plectin/desmosome). Optional tangential **crosslinks** (plectin
cytolinker) knit spokes so load spreads. This is the hard-fact reframe: the coupling is the radial spokes +
anchors, NOT the cage shell.

MECHANICS — LINEAR-first (this module). IFs are soft-then-strain-stiffening (Kreplak 2005; Block 2018), but
the LINEAR regime alone (backbone slope k_bb = E_if·A_if/l_seg, KB-anchored G1 gate) already flips the
falsifier (nucleus unfreezes). Because every IF bond (backbone + LINC + anchor + crosslink) is Hookean in
this regime, they all ride the existing ``link_spring_kernel`` — NO new kernel needed. The nonlinear
3-regime strain-stiffening (unfolding plateau + re-stiffening) is deferred behind KB registration of
F_yield/λ1/λ2/k_stiff (PI-authored) — do NOT fake it with a magic-number fit.

⚠ MAGIC-NUMBERS → PI (never tuned to make coupling appear; no-param-tuning HARD rule):
  - ``n_fil`` (native cage count) — no clean MCF7 datum; a density anchor is PI-authored. Provisional here.
  - ``ratio_xl`` (crosslink/backbone stiffness), ``k_linc``, ``k_anchor`` — provisional, PI registration pending.
If the coupling is weak at the KB-anchored k_bb, that is a FINDING to surface, not a knob to turn.

Units: FF is µm, pN, s. E in Pa (= pN/µm²), A in µm², k in pN/µm.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# ── Magic-Number Block (literature anchors + range-check bands) ───────────────────────────────────
E_IF_PA = 6.0e6            # IF axial small-strain Young's modulus [Pa=pN/µm²] (Kreplak 2005 / Guo 2013; band [1,10] MPa)
E_IF_BAND = (1.0e6, 1.0e7)
D_IF_UM = 0.010           # assembled IF diameter 10 nm (Mücke 2004, canonical TEM)
L_SEG_UM = 0.5            # backbone bead spacing ~ Lp (vimentin Lp 0.4–1 µm Mücke 2004; keratin 0.3–0.5 Lichtenstern 2012)
LP_BAND_UM = (0.3, 1.0)
RATIO_XL = 0.1            # ⚠ crosslink/backbone stiffness ratio (plectin) — MAGIC-NUMBER → PI
# LINC + cortex-anchor stiffness: provisional (PI). Kept O(k_bb) so the load path is stiff enough to transmit
# but not a rigid pin (over-stiff guard: the native gate checks the nucleus deforms LESS than the cortex).
RATIO_LINC = 1.0         # ⚠ LINC (nesprin/SUN) stiffness / k_bb — MAGIC-NUMBER → PI (LINC tension oracle ~8 pN)
RATIO_ANCHOR = 1.0       # ⚠ cortex-anchor (plectin/desmosome) stiffness / k_bb — MAGIC-NUMBER → PI


@dataclass
class ResolvedIF:
    """Resolved IF cage parameters in FF units (µm, pN)."""
    k_bb: float           # backbone harmonic stiffness E_if·A_if/l_seg [pN/µm]
    k_xl: float           # crosslink stiffness ratio_xl·k_bb [pN/µm]
    k_linc: float         # LINC (nucleus-envelope) tether stiffness [pN/µm]
    k_anchor: float       # cortex-anchor tether stiffness [pN/µm]
    l_seg_um: float       # backbone rest spacing [µm]
    n_fil: int            # number of radial spokes
    E_if_Pa: float


def resolve_intermediate_filaments(*, n_fil: int, E_if_Pa: float = E_IF_PA, d_if_um: float = D_IF_UM,
                                   l_seg_um: float = L_SEG_UM, ratio_xl: float = RATIO_XL) -> ResolvedIF:
    """Resolve the IF cage: k_bb = E_if·A_if/l_seg (grid-invariant rod-segment axial spring), derived — no
    magic number in the stiffness. Bands range-checked. n_fil / the LINC & anchor ratios stay PI-flagged."""
    if not (E_IF_BAND[0] <= E_if_Pa <= E_IF_BAND[1]):
        raise ValueError(f"E_if {E_if_Pa} Pa outside band {E_IF_BAND} (Kreplak/Guo)")
    if n_fil < 1 or d_if_um <= 0 or l_seg_um <= 0:
        raise ValueError(f"n_fil={n_fil}, d_if={d_if_um}, l_seg={l_seg_um} must be positive")
    A_if = np.pi * (d_if_um / 2.0) ** 2                    # µm²
    k_bb = E_if_Pa * A_if / l_seg_um                       # pN/µm (Pa·µm² / µm = pN/µm)
    return ResolvedIF(k_bb=float(k_bb), k_xl=float(ratio_xl * k_bb), k_linc=float(RATIO_LINC * k_bb),
                      k_anchor=float(RATIO_ANCHOR * k_bb), l_seg_um=float(l_seg_um), n_fil=int(n_fil),
                      E_if_Pa=float(E_if_Pa))


def _fibonacci_directions(n: int) -> np.ndarray:
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)


@dataclass
class IFCage:
    """Built IF cage: tail-block positions + Hookean bond arrays (GLOBAL indices), ready to append to xl_*_all."""
    pos: np.ndarray        # (n_if,3) spoke bead positions [µm]
    n_if: int
    bond_i: np.ndarray     # (Nb,) global i
    bond_j: np.ndarray     # (Nb,) global j
    bond_k: np.ndarray     # (Nb,) stiffness [pN/µm]
    bond_rest: np.ndarray  # (Nb,) rest length [µm] (force-free at construction)
    n_backbone: int
    n_linc: int
    n_anchor: int
    n_xl: int


def build_if_cage(centre, R_nuc_um: float, R_cortex_um: float, nuc_pos: np.ndarray, cortex_pos: np.ndarray,
                  res: ResolvedIF, *, nuc_offset: int, cortex_offset: int, if_offset: int,
                  crosslink: bool = True) -> IFCage:
    """Build the radial-spoke IF cage. Spokes run from R_nuc·(just outside) to R_cortex·(just inside) along
    ``res.n_fil`` Fibonacci directions. Returns the tail-block positions + Hookean bonds in GLOBAL indices:
    backbone (consecutive spoke beads, k_bb) + LINC (inner bead ↔ nearest nucleus bead, k_linc) + cortex
    anchor (outer bead ↔ nearest cortex node, k_anchor) + optional tangential crosslinks (k_xl). All rest
    lengths = the construction distance (force-free at build, physiological-baseline rule)."""
    from scipy.spatial import cKDTree
    c = np.asarray(centre, np.float64)
    dirs = _fibonacci_directions(res.n_fil)
    r0, r1 = R_nuc_um + 0.15, R_cortex_um - 0.15          # spoke spans the cytoplasm gap (just off both surfaces)
    nb = max(2, int(round((r1 - r0) / res.l_seg_um)) + 1)  # beads per spoke
    radii = np.linspace(r0, r1, nb)                        # (nb,)
    pos = (c[None, None, :] + radii[None, :, None] * dirs[:, None, :]).reshape(-1, 3)  # (n_fil*nb, 3)
    n_if = pos.shape[0]
    gidx = if_offset + np.arange(n_if)                     # global indices of IF beads

    bi, bj, bk, br = [], [], [], []
    # backbone: consecutive beads within each spoke
    for a in range(res.n_fil):
        base = a * nb
        for b in range(nb - 1):
            i, j = base + b, base + b + 1
            bi.append(gidx[i]); bj.append(gidx[j]); bk.append(res.k_bb)
            br.append(float(np.linalg.norm(pos[i] - pos[j])))
    n_backbone = len(bi)

    inner = (np.arange(res.n_fil) * nb)                    # inner bead of each spoke (near nucleus)
    outer = (np.arange(res.n_fil) * nb + (nb - 1))         # outer bead of each spoke (near cortex)
    # LINC: inner bead ↔ nearest NUCLEUS surface bead
    nuc_tree = cKDTree(nuc_pos)
    _, nn = nuc_tree.query(pos[inner])
    for s in range(res.n_fil):
        bi.append(gidx[inner[s]]); bj.append(nuc_offset + int(nn[s])); bk.append(res.k_linc)
        br.append(float(np.linalg.norm(pos[inner[s]] - nuc_pos[int(nn[s])])))
    n_linc = res.n_fil
    # cortex anchor: outer bead ↔ nearest CORTEX node
    cortex_tree = cKDTree(cortex_pos)
    _, cn = cortex_tree.query(pos[outer])
    for s in range(res.n_fil):
        bi.append(gidx[outer[s]]); bj.append(cortex_offset + int(cn[s])); bk.append(res.k_anchor)
        br.append(float(np.linalg.norm(pos[outer[s]] - cortex_pos[int(cn[s])])))
    n_anchor = res.n_fil
    # tangential crosslinks: knit neighbouring spokes at matching radius (plectin), one per bead to the nearest
    # bead on a DIFFERENT spoke within reach — spreads load so the cage fails globally, not locally.
    n_xl = 0
    if crosslink and res.n_fil > 1:
        reach = 1.5 * res.l_seg_um
        spoke_of = np.repeat(np.arange(res.n_fil), nb)
        tree = cKDTree(pos)
        pairs = tree.query_pairs(reach, output_type="ndarray")
        for (u, v) in pairs:
            if spoke_of[u] != spoke_of[v]:                # cross-spoke only
                bi.append(gidx[u]); bj.append(gidx[v]); bk.append(res.k_xl)
                br.append(float(np.linalg.norm(pos[u] - pos[v]))); n_xl += 1

    return IFCage(pos=np.ascontiguousarray(pos, np.float64), n_if=n_if,
                  bond_i=np.asarray(bi, np.int64), bond_j=np.asarray(bj, np.int64),
                  bond_k=np.asarray(bk, np.float64), bond_rest=np.asarray(br, np.float64),
                  n_backbone=n_backbone, n_linc=n_linc, n_anchor=n_anchor, n_xl=n_xl)


__all__ = ["ResolvedIF", "resolve_intermediate_filaments", "IFCage", "build_if_cage",
           "E_IF_PA", "D_IF_UM", "L_SEG_UM"]
