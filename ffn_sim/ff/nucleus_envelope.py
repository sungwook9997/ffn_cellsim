r"""Nuclear envelope as a real surface (FF, Warp, µm·pN·s) — Phase 3, nucleus-stack piece 1.

The live FF nucleus is a bead CLOUD held by an independent radial bilinear spring (``nucleus_shell_kernel``) —
so the nucleus cannot behave as a coherent shell (no lateral area coupling, no bending, no rupture). This makes
the nuclear envelope a real closed triangulated surface (reusing the plasma-membrane mesh + area-tension infra)
with the strain-stiffening lamina and EMERGENT rupture — the first fine-graining of the lumped nucleus.

Mechanics (KB-3.B2, PI-ratified E_nuc = 399 Pa MCF7 in-situ):
  - Bilinear STRAIN-STIFFENING areal tension: at low areal strain the soft CHROMATIN governs (2D modulus
    K_chrom = E_nuc·h_lamina); past the lamin-engagement knee (~10% strain, KB-3.B2.2) the stiff LAMIN-A/C
    meshwork engages (K_lamin = ratio·K_chrom). Chromatin-low-strain → lamin-high-strain is the measured nuclear
    response (KB-3.B2.1/.2).
  - Helfrich bending κ_NE resists high curvature (confined-migration / blebbing regime).
  - EMERGENT nuclear-envelope RUPTURE (KB-3.B2.4: cancers down-regulate lamin-A/C → NE rupture + DNA damage):
    a face whose areal strain exceeds the rupture strain ``eps_rupture`` tears (its tension drops to 0).

No magic numbers on the mechanics: the areal moduli are E_nuc·(lamina thickness)·ratio, all KB-3.B2. ⚠ SURFACED:
the lamina thickness h_lamina (~15 nm) and the rupture strain eps_rupture are not MCF7-specific single values
(the lit-ingest flagged NE rupture as "curvature+lamin-context, may be irreducible to one number") → both are
flagged, PI-proxy defaults, NOT tuned to a gate.

Gate (analytic): (1) Young-Laplace — a pressurised nuclear shell at tension σ balances internal ΔP = 2σ/R_nuc
(reuses the plasma-membrane Laplace machinery, already verified); (2) strain-stiffening KNEE — the tangent areal
modulus jumps ×ratio at 10% strain; (3) rupture EMERGES above eps_rupture and is absent below.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from ffn_sim.ff.membrane_surface import build_membrane_mesh, MembraneMesh

# KB-3.B2 nuclear-envelope constants
E_NUC_PA = 399.0               # MCF7 in-situ (Fischer 2020; PI-ratified 2026-07-07) [Pa = pN/µm²]
H_LAMINA_UM = 0.015            # ⚠ lamina thickness ~15 nm (surfaced proxy, not MCF7-specific)
RATIO_LAMIN = 3.0             # lamin-engaged / chromatin tangent stiffness (KB-3.B2.1 in-situ 1.4–5×)
KNEE_STRAIN = 0.10            # areal strain at lamin engagement (KB-3.B2.2 ~10%)
KAPPA_NE_PN_UM = 0.0828       # envelope bending rigidity ~20 kBT (lamina+bilayer; same order as plasma membrane)
EPS_RUPTURE = 0.50            # ⚠ areal strain for NE rupture (surfaced proxy; KB-3.B2.3/.4 lamin-context)


@dataclass
class NuclearEnvelope:
    mesh: MembraneMesh
    K_chrom: float               # chromatin-regime areal modulus [pN/µm]
    K_lamin: float               # lamin-engaged areal modulus [pN/µm]
    knee_strain: float
    eps_rupture: float
    R_nuc: float


def resolve_nuclear_envelope(R_nuc_um: float, *, subdivisions: int = 3, centre=(0, 0, 0),
                             E_nuc_pa: float = E_NUC_PA) -> NuclearEnvelope:
    """Nuclear-envelope surface + KB-3.B2 areal moduli. ``K_chrom = E_nuc·h_lamina`` (soft chromatin regime),
    ``K_lamin = ratio·K_chrom`` (stiff lamin regime past the knee)."""
    mesh = build_membrane_mesh(R_nuc_um, subdivisions=subdivisions, centre=centre)
    K_chrom = E_nuc_pa * H_LAMINA_UM                       # [pN/µm² · µm] = pN/µm
    return NuclearEnvelope(mesh=mesh, K_chrom=float(K_chrom), K_lamin=float(RATIO_LAMIN * K_chrom),
                           knee_strain=KNEE_STRAIN, eps_rupture=EPS_RUPTURE, R_nuc=float(R_nuc_um))


def lamina_tension(area: float, A0: float, ne: NuclearEnvelope) -> float:
    """Bilinear strain-stiffening areal tension σ(ε) [pN/µm]: soft chromatin below the knee, stiff lamin above.
    Ruptured (ε > eps_rupture) → 0. ``ε = (A−A0)/A0`` areal strain."""
    eps = (area - A0) / A0
    if eps > ne.eps_rupture:
        return 0.0                                         # torn envelope
    if eps <= ne.knee_strain:
        return ne.K_chrom * eps                            # chromatin regime
    return ne.K_chrom * ne.knee_strain + ne.K_lamin * (eps - ne.knee_strain)   # lamin engaged (stiffer tangent)


def tangent_modulus(eps: float, ne: NuclearEnvelope) -> float:
    """dσ/dε — the tangent areal modulus (K_chrom below the knee, K_lamin above; 0 if ruptured)."""
    if eps > ne.eps_rupture:
        return 0.0
    return ne.K_chrom if eps <= ne.knee_strain else ne.K_lamin


@wp.kernel
def nuclear_envelope_kernel(npos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                            sigma: wp.float64, ruptured: wp.array(dtype=wp.int32),
                            nforce: wp.array(dtype=wp.vec3d)):
    """Areal-tension force of the nuclear-envelope surface (same area-gradient form as the plasma membrane),
    skipping faces flagged ``ruptured`` (their tension is 0 = the tear). ``sigma`` = the bilinear tension for the
    current global areal strain (computed once per step)."""
    t = wp.tid()
    if ruptured[t] == 1:
        return
    i0 = faces[t, 0]; i1 = faces[t, 1]; i2 = faces[t, 2]
    p0 = npos[i0]; p1 = npos[i1]; p2 = npos[i2]
    nrm = wp.cross(p1 - p0, p2 - p0)
    a2 = wp.length(nrm)
    if a2 < wp.float64(1e-18):
        return
    nhat = nrm / a2
    wp.atomic_add(nforce, i0, -sigma * wp.float64(0.5) * wp.cross(nhat, p2 - p1))
    wp.atomic_add(nforce, i1, -sigma * wp.float64(0.5) * wp.cross(nhat, p0 - p2))
    wp.atomic_add(nforce, i2, -sigma * wp.float64(0.5) * wp.cross(nhat, p1 - p0))


__all__ = ["NuclearEnvelope", "resolve_nuclear_envelope", "lamina_tension", "tangent_modulus",
           "nuclear_envelope_kernel", "E_NUC_PA", "H_LAMINA_UM", "RATIO_LAMIN", "KNEE_STRAIN", "EPS_RUPTURE"]
