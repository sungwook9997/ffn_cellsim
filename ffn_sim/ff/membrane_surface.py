"""Independent plasma-membrane surface — H.8 Template 2 (fully-dynamical lipid sheet), FF/Warp, µm·pN·s.

The plasma membrane as its OWN entity, not the convex-hull boundary of the cortex filaments (PI 2026-07-07:
"플라즈마 멤브레인은 필라멘트의 경계가 아닌 독립적인 개체로서 존재해야"). A separate closed triangulated
sheet (its own node set) that wraps the cortex, carries lipid-bilayer mechanics, is tethered to the cortex by
ERM linkers, and is the true osmotic envelope that CONTAINS the cortex (so escaping filament tips cannot leak
through the bilayer). This is the H.8 brief's Template-2 refinement (Template 1 = a lumped tension on cortex
nodes, which is what the model had).

Mechanics (KU-3.B1, all lit-anchored — see ``ffn_sim/common/compartments.py::resolve_membrane``):
  - in-plane tension γ_mem  = 10  pN/µm   (bilayer-only tension T_m; KB-3.B1.1)
  - area-expansion K_A       = 2.35e5 pN/µm (0.235 N/m; Rawicz 2000, KB-3.B1.3), reservoir-buffered:
    constant γ_mem until areal strain exceeds the reservoir capacity f_excess (folds/microvilli unfold —
    supplies area for spreading), then the steep K_A upturn (lysis at τ_lysis).
  - Helfrich bending κ_m     = 0.0828 pN·µm (20 kBT; KB-3.B1.2)
  - ERM tether (membrane↔cortex): breakable spring, rupture ≈ f_t = 2π√(2κ_m(T_m+γ_MCA)) ≈ 5–40 pN
    (KB-3.B1.4) — de-tachment (blebbing) is EMERGENT from rupture, not scripted.

Units: FF is µm, pN, s. γ,K_A,γ_MCA in pN/µm; κ in pN·µm; f_t in pN.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from ffn_sim.common.surface_manifold import icosphere
from ffn_sim.common.compartments import resolve_membrane, ResolvedMembrane

# Membrane–cortex adhesion (KU-3.B1.4): γ_MCA ≈ 1e-5 J/m² = 10 pN/µm (band 1e-6–1e-4 J/m² = 1–1000 pN/µm).
GAMMA_MCA_PN_UM = 10.0
# Reservoir capacity: apparent-area strain the fold/microvillus reservoir absorbs at constant tension before the
# K_A elastic upturn. The bilayer itself lyses at ~2–5% areal strain (KB-3.B1.3); the whole-cell membrane
# RESERVOIR buffers far more apparent-area change (Raucher-Sheetz plateau) — cells spread to 2–3× footprint by
# unfolding. Set to the spreading-relevant reservoir, NOT the bilayer lysis strain.
RESERVOIR_STRAIN = 0.60          # apparent-area reservoir before K_A engages (fold/microvillus unfolding)


@dataclass
class MembraneMesh:
    """A closed triangulated lipid sheet (its own nodes), FF units (µm)."""
    verts: np.ndarray            # (Nm,3) node positions [µm]
    faces: np.ndarray            # (Ntri,3) int outward-wound triangles
    edges: np.ndarray            # (Ne,2) unique undirected edges (bending/diagnostics)
    A0: float                    # reference total area 4πR_mem² [µm²]
    R_mem: float                 # membrane radius [µm]


def build_membrane_mesh(R_mem: float, subdivisions: int = 4, centre=(0.0, 0.0, 0.0)) -> MembraneMesh:
    """Icosphere lipid sheet of radius ``R_mem`` (subdivisions=4 → 5120 tri / 2562 nodes: smooth, light —
    the membrane is a continuum sheet, it does NOT need the cortex's filament-node resolution)."""
    v, f = icosphere(subdivisions, R_mem)
    v = v + np.asarray(centre, np.float64)[None, :]
    e = np.unique(np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1), axis=0)
    A0 = float(4.0 * np.pi * R_mem ** 2)
    return MembraneMesh(verts=np.ascontiguousarray(v, np.float64), faces=np.ascontiguousarray(f, np.int64),
                        edges=np.ascontiguousarray(e, np.int64), A0=A0, R_mem=float(R_mem))


def resolve_membrane_full() -> ResolvedMembrane:
    """The KU-3.B1 membrane parameter set (delegates to the validated ``compartments.resolve_membrane``)."""
    return resolve_membrane()


# ----------------------------------------------------------------------------------------------------------
# Warp kernels — the membrane is nodes [0, Nm); ``mpos`` is its own position array.
# ----------------------------------------------------------------------------------------------------------
@wp.kernel
def membrane_area_kernel(mpos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                         sigma: wp.float64, mforce: wp.array(dtype=wp.vec3d)):
    """In-plane surface TENSION σ on each triangle (lipid bilayer = fluid → area-elastic, ~no shear). The
    area-gradient force on vertex i is ``f_i = −σ·∂A/∂x_i = −σ·½·(n̂ × e_opp)`` (e_opp = the opposite edge),
    which minimises area ⇒ a taut, pressure-balanced sheet. ``sigma`` is the reservoir-buffered tension
    (γ_mem, upturning to K_A past the reservoir) computed once per step from the global areal strain."""
    t = wp.tid()
    i0 = faces[t, 0]; i1 = faces[t, 1]; i2 = faces[t, 2]
    p0 = mpos[i0]; p1 = mpos[i1]; p2 = mpos[i2]
    nrm = wp.cross(p1 - p0, p2 - p0)                       # 2A·n̂
    a2 = wp.length(nrm)
    if a2 < wp.float64(1.0e-18):
        return
    nhat = nrm / a2
    # ∂A/∂p_i = ½ (n̂ × edge_opposite_to_i); tension force = −σ·∂A/∂p_i
    f0 = -sigma * wp.float64(0.5) * wp.cross(nhat, p2 - p1)
    f1 = -sigma * wp.float64(0.5) * wp.cross(nhat, p0 - p2)
    f2 = -sigma * wp.float64(0.5) * wp.cross(nhat, p1 - p0)
    wp.atomic_add(mforce, i0, f0)
    wp.atomic_add(mforce, i1, f1)
    wp.atomic_add(mforce, i2, f2)


@wp.kernel
def membrane_area_reduce_kernel(mpos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                                area_out: wp.array(dtype=wp.float64)):
    """Σ triangle areas → area_out[0] (on-device total membrane area for the reservoir/K_A tension)."""
    t = wp.tid()
    p0 = mpos[faces[t, 0]]; p1 = mpos[faces[t, 1]]; p2 = mpos[faces[t, 2]]
    wp.atomic_add(area_out, 0, wp.float64(0.5) * wp.length(wp.cross(p1 - p0, p2 - p0)))


@wp.kernel
def erm_tether_kernel(mpos: wp.array(dtype=wp.vec3d), cpos: wp.array(dtype=wp.vec3d),
                      m_idx: wp.array(dtype=wp.int32), c_idx: wp.array(dtype=wp.int32),
                      bound: wp.array(dtype=wp.int32), k_erm: wp.float64, rest: wp.array(dtype=wp.float64),
                      f_rupt: wp.float64,
                      mforce: wp.array(dtype=wp.vec3d), cforce: wp.array(dtype=wp.vec3d)):
    """ERM linker (ezrin/radixin/moesin) — a breakable spring tethering membrane node ``m_idx[k]`` to cortex
    node ``c_idx[k]``. Equal-and-opposite Hookean force k_erm·(L−rest[k]) along the link (rest[k] = the tether's
    own formation length, so it is force-free at construction); if the tension exceeds the rupture force f_rupt
    (≈ f_t, KU-3.B1.4) the tether UNBINDS (bound→0) ⇒ the membrane locally detaches from the cortex = a BLEB,
    emergent not scripted. Bound tethers keep the membrane riding on the cortex."""
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = m_idx[k]; ci = c_idx[k]
    d = mpos[mi] - cpos[ci]
    L = wp.length(d)
    if L < wp.float64(1.0e-12):
        return
    u = d / L
    tension = k_erm * (L - rest[k])                       # >0 when stretched (membrane pulled off cortex)
    if tension > f_rupt:                                  # ERM rupture → local bleb
        bound[k] = 0
        return
    f = tension * u                                       # on membrane node toward cortex when stretched
    wp.atomic_add(mforce, mi, -f)
    wp.atomic_add(cforce, ci, f)


@wp.kernel
def membrane_containment_kernel(cpos: wp.array(dtype=wp.vec3d), centre: wp.vec3d,
                                r_env: wp.array(dtype=wp.float64), k_wall: wp.float64,
                                cforce: wp.array(dtype=wp.vec3d)):
    """The membrane is the OUTER envelope: any cortex node whose radius exceeds the local membrane radius
    ``r_env`` (an angular lookup of the sheet) is pushed back INWARD with a stiff one-sided penalty
    ``k_wall·(r−r_env)``. This is what stops filament tips leaking OUT through the bilayer (the straggler
    defect) — the lipid sheet is impermeable to actin. One thread per cortex node."""
    i = wp.tid()
    d = cpos[i] - centre
    r = wp.length(d)
    re = r_env[0]                                          # conservative single-radius envelope (min membrane r)
    if r > re:
        wp.atomic_add(cforce, i, -k_wall * (r - re) * (d / r))


def reservoir_tension(area: float, A0: float, mem: ResolvedMembrane) -> float:
    """Reservoir-buffered membrane tension σ(A) [pN/µm]: constant γ_mem while the apparent-area strain is
    within the fold/microvillus reservoir; past it the steep K_A elastic upturn (Rawicz), capped at τ_lysis.
    Below A0 the sheet is slack → σ = γ_mem baseline (Raucher-Sheetz plateau)."""
    strain = (area - A0) / A0
    if strain <= RESERVOIR_STRAIN:
        return mem.gamma_mem
    excess = strain - RESERVOIR_STRAIN
    return min(mem.gamma_mem + mem.K_A * excess, mem.tau_lysis)


__all__ = ["MembraneMesh", "build_membrane_mesh", "resolve_membrane_full", "reservoir_tension",
           "membrane_area_kernel", "membrane_area_reduce_kernel", "erm_tether_kernel",
           "membrane_containment_kernel", "GAMMA_MCA_PN_UM", "RESERVOIR_STRAIN"]
