"""Deformable-mesh nuclear envelope — Warp kernel SOURCE (authored on the dev Mac, run on gbook CUDA).

The GPU-resident force primitives of the I2 nucleus. Every kernel is one-thread-per-element,
``atomic_add``, float64 — mirroring ``ff/forces_warp.py`` and ``ff/membrane_surface.py`` idioms.
Per §1.2 these are NET-NEW under ``ac/nucleus/`` and are NOT executed on the Mac (no CUDA device);
they are gated by the pure-NumPy oracles in ``*_analytic.py`` (analytic-first, §1.5), and the lead runs
the native CUDA gates on gbook.

Force families (each an ``accumulate``-style contribution the lead-owned integrator sums, §1.4):
  * BENDING — REUSES ``ff.membrane_surface.helfrich_bending_kernel`` read-only (it only touches the 4
    hinge vertices it names, so it composes additively; κ̃ = 8π·κ_NE/Σ_ref from ``calibrate_kappa_tilde``).
  * LAMINA areal tension — NEW: PER-FACE local areal strain → framework-#6 3-regime tension
    (chromatin+lamin-B soft → lamin-A/C stiff past the knee), skipping ruptured faces. Fine-grained
    per-face (each face tears independently), NOT the ff single global-σ scaffold.
  * NUCLEOPLASM volume — NEW: incompressible (ν→½) volume penalty p_vol=−K_vol(V−V0)/V0 applied along
    the divergence-theorem volume gradient ∂V/∂x_i.
  * NUCLEOPLASM viscosity — NEW: overdamped viscous drag f=−γ_nuc·v on envelope/chromatin nodes.
  * LINC tether — NEW: nonlinear tension-only nesprin link (nucleus node ↔ cytoskeleton anchor).
  * CHROMATIN net — NEW: internal WLC (Marko–Siggia) polymer links (finite-extensible, small-strain soft).

Host helpers (pure NumPy, testable) compute the per-face reference areas + the rupture update; the
kernels consume them. This module imports ``warp`` (allowed — Warp is the mandated runtime) but is
import-only on the Mac; kernel launches happen on the gbook A5000.

Units: FF µm·pN·s. κ̃ [pN·µm]; areal moduli [pN/µm]; K_vol [pN/µm²]; γ_nuc [pN·s/µm]; k_linc [pN/µm].
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.nucleus.geometry import build_hinges, build_oblate_mesh, face_normals_areas
from aleph.components.nucleus.lamina_analytic import LaminaParams, calibrate_kappa_tilde

# Reuse the validated membrane Helfrich bending kernel read-only (do NOT edit ff/; §1.2/§2).
from aleph.laws.membrane_surface import (
    helfrich_bending_kernel,  # noqa: F401  (re-exported for the integrator)
)

__all__ = [
    "NucleusMesh",
    "build_nucleus",
    "helfrich_bending_kernel",
    "lamina_areal_tension_kernel",
    "rupture_update_kernel",
    "conditional_rupture_update_kernel",
    "nucleoplasm_volume_reduce_kernel",
    "nucleoplasm_volume_force_kernel",
    "nucleoplasm_viscosity_kernel",
    "linc_tether_kernel",
    "chromatin_wlc_kernel",
    "face_reference_areas",
]


@dataclass
class NucleusMesh:
    """The deformable nuclear envelope + its precomputed reference state (host build → device upload)."""
    verts: npt.NDArray[np.float64]        # (Nv,3) envelope nodes [µm]
    faces: npt.NDArray[np.int64]          # (Nf,3) outward-wound triangles
    hinges: npt.NDArray[np.int64]         # (Nh,4) dihedral hinges [i,j,k,l]
    A0_face: npt.NDArray[np.float64]      # (Nf,) per-face reference areas [µm²]
    V0: float                             # reference nucleoplasm volume [µm³]
    kappa_tilde: float                    # calibrated dihedral coupling [pN·µm]
    params: LaminaParams
    R_eq: float
    aspect: float


def face_reference_areas(
    verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    """Per-face reference areas A0[t] [µm²] of the resting envelope (the per-face strain denominator)."""
    _, areas = face_normals_areas(verts, faces)
    return np.ascontiguousarray(areas, np.float64)


def build_nucleus(
    r_eq: float, params: LaminaParams, *, aspect: float = 1.0, subdivisions: int = 3, centre=(0, 0, 0)
) -> NucleusMesh:
    """Assemble the resting deformable nucleus: oblate mesh + hinges + per-face reference areas +
    reference volume + κ̃ calibrated so the resting bending energy equals the continuum 8πκ_NE."""
    verts, faces = build_oblate_mesh(r_eq, aspect, subdivisions, centre)
    hinges = build_hinges(faces)
    A0 = face_reference_areas(verts, faces)
    from aleph.components.nucleus.geometry import mesh_volume
    V0 = mesh_volume(verts, faces)
    kt = calibrate_kappa_tilde(params.kappa_ne, verts, hinges)
    return NucleusMesh(verts=verts, faces=faces, hinges=hinges, A0_face=A0, V0=float(V0),
                       kappa_tilde=float(kt), params=params, R_eq=float(r_eq), aspect=float(aspect))


# ---------------------------------------------------------------------------------------------------
# LAMINA areal tension — per-face framework-#6 3-regime split, skipping ruptured faces
# ---------------------------------------------------------------------------------------------------
@wp.func
def _lamina_sigma(eps: wp.float64, k_soft: wp.float64, k_ac: wp.float64,
                  knee: wp.float64, eps_rupt: wp.float64) -> wp.float64:
    """σ(ε): soft (chromatin+lamin-B) below the knee, stiff lamin-A/C above, 0 past rupture."""
    if eps > eps_rupt:
        return wp.float64(0.0)
    if eps <= knee:
        return k_soft * eps
    return k_soft * knee + k_ac * (eps - knee)


@wp.kernel
def lamina_areal_tension_kernel(npos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                                a0: wp.array(dtype=wp.float64), ruptured: wp.array(dtype=wp.int32),
                                k_soft: wp.float64, k_ac: wp.float64, knee: wp.float64,
                                eps_rupt: wp.float64, nforce: wp.array(dtype=wp.vec3d)):
    """Per-face areal-tension force (area-gradient form ``f_i=−σ·½·(n̂×e_opp)``) with a LOCAL per-face
    strain ε=(A−A0)/A0 → framework-#6 3-regime σ(ε). Ruptured faces (flag=1) contribute 0 = the tear.
    Fine-grained: each face carries its own strain/tension/rupture (not one global σ)."""
    t = wp.tid()
    if ruptured[t] == 1:
        return
    i0 = faces[t, 0]; i1 = faces[t, 1]; i2 = faces[t, 2]
    p0 = npos[i0]; p1 = npos[i1]; p2 = npos[i2]
    nrm = wp.cross(p1 - p0, p2 - p0)
    a2 = wp.length(nrm)                                    # = 2·area
    if a2 < wp.float64(1.0e-18):
        return
    area = wp.float64(0.5) * a2
    eps = (area - a0[t]) / a0[t]
    sigma = _lamina_sigma(eps, k_soft, k_ac, knee, eps_rupt)
    nhat = nrm / a2
    wp.atomic_add(nforce, i0, -sigma * wp.float64(0.5) * wp.cross(nhat, p2 - p1))
    wp.atomic_add(nforce, i1, -sigma * wp.float64(0.5) * wp.cross(nhat, p0 - p2))
    wp.atomic_add(nforce, i2, -sigma * wp.float64(0.5) * wp.cross(nhat, p1 - p0))


@wp.kernel
def rupture_update_kernel(npos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
                          a0: wp.array(dtype=wp.float64), eps_rupt: wp.float64,
                          ruptured: wp.array(dtype=wp.int32)):
    """EMERGENT per-face rupture: once a face's areal strain exceeds ``eps_rupt`` it tears
    IRREVERSIBLY (ruptured→1). Rupture appears past the sourced threshold and is absent below —
    the falsifiable gate; the threshold is GAP (I0-B2), never tuned to a gate."""
    t = wp.tid()
    if ruptured[t] == 1:
        return
    p0 = npos[faces[t, 0]]; p1 = npos[faces[t, 1]]; p2 = npos[faces[t, 2]]
    area = wp.float64(0.5) * wp.length(wp.cross(p1 - p0, p2 - p0))
    if (area - a0[t]) / a0[t] > eps_rupt:
        ruptured[t] = 1


@wp.kernel
def conditional_rupture_update_kernel(
    npos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    a0: wp.array(dtype=wp.float64),
    eps_rupt: wp.float64,
    accepted: wp.array(dtype=wp.int32),
    ruptured: wp.array(dtype=wp.int32),
) -> None:
    """Commit nuclear-envelope rupture only when the device outer-acceptance latch is set."""
    t = wp.tid()
    if accepted[0] == 0 or ruptured[t] == 1:
        return
    p0 = npos[faces[t, 0]]
    p1 = npos[faces[t, 1]]
    p2 = npos[faces[t, 2]]
    area = wp.float64(0.5) * wp.length(wp.cross(p1 - p0, p2 - p0))
    if (area - a0[t]) / a0[t] > eps_rupt:
        ruptured[t] = 1


# ---------------------------------------------------------------------------------------------------
# NUCLEOPLASM volume — incompressible (ν→½) penalty along the divergence-theorem volume gradient
# ---------------------------------------------------------------------------------------------------
@wp.kernel
def nucleoplasm_volume_reduce_kernel(npos: wp.array(dtype=wp.vec3d),
                                     faces: wp.array(dtype=wp.int32, ndim=2),
                                     vol_out: wp.array(dtype=wp.float64)):
    """Σ signed tetra volumes → vol_out[0] = enclosed nucleoplasm volume (divergence theorem)."""
    t = wp.tid()
    p0 = npos[faces[t, 0]]; p1 = npos[faces[t, 1]]; p2 = npos[faces[t, 2]]
    wp.atomic_add(vol_out, 0, wp.dot(p0, wp.cross(p1, p2)) / wp.float64(6.0))


@wp.kernel
def nucleoplasm_volume_force_kernel(npos: wp.array(dtype=wp.vec3d),
                                    faces: wp.array(dtype=wp.int32, ndim=2),
                                    p_vol: wp.float64, nforce: wp.array(dtype=wp.vec3d)):
    """Incompressible nucleoplasm force: uniform pressure ``p_vol=−K_vol(V−V0)/V0`` [pN/µm²] pushing
    each vertex along its volume gradient ∂V/∂x_i=(1/6)Σ(x_j×x_k). Positive p_vol (V<V0) inflates."""
    t = wp.tid()
    i0 = faces[t, 0]; i1 = faces[t, 1]; i2 = faces[t, 2]
    p0 = npos[i0]; p1 = npos[i1]; p2 = npos[i2]
    wp.atomic_add(nforce, i0, p_vol * wp.cross(p1, p2) / wp.float64(6.0))
    wp.atomic_add(nforce, i1, p_vol * wp.cross(p2, p0) / wp.float64(6.0))
    wp.atomic_add(nforce, i2, p_vol * wp.cross(p0, p1) / wp.float64(6.0))


@wp.kernel
def nucleoplasm_viscosity_kernel(nvel: wp.array(dtype=wp.vec3d), gamma_nuc: wp.float64,
                                 nforce: wp.array(dtype=wp.vec3d)):
    """Overdamped nucleoplasm viscous drag ``f=−γ_nuc·v`` on each envelope/chromatin node (γ_nuc from
    the nucleoplasm viscosity, I0-B2 GAP). One thread per node."""
    i = wp.tid()
    wp.atomic_add(nforce, i, -gamma_nuc * nvel[i])


# ---------------------------------------------------------------------------------------------------
# LINC tether — nonlinear tension-only nesprin link (nucleus node ↔ cytoskeleton anchor)
# ---------------------------------------------------------------------------------------------------
@wp.kernel
def linc_tether_kernel(npos: wp.array(dtype=wp.vec3d), apos: wp.array(dtype=wp.vec3d),
                       n_idx: wp.array(dtype=wp.int32), a_idx: wp.array(dtype=wp.int32),
                       rest: wp.array(dtype=wp.float64), k_linc: wp.float64, stiffening: wp.float64,
                       nforce: wp.array(dtype=wp.vec3d), aforce: wp.array(dtype=wp.vec3d)):
    """Nonlinear LINC tether: tension ``f=k_linc·Δ·(1+stiffening·Δ²)`` (Δ=L−rest, Δ≥0 tension-only)
    along the link between nucleus node ``n_idx[k]`` and cytoskeleton anchor ``a_idx[k]``. Force-free
    at the formation length. ⚠ k_linc is GAP (I0-B7: "8 pN is a TENSION not a stiffness")."""
    k = wp.tid()
    ni = n_idx[k]; ai = a_idx[k]
    d = apos[ai] - npos[ni]
    L = wp.length(d)
    if L < wp.float64(1.0e-12):
        return
    delta = L - rest[k]
    if delta <= wp.float64(0.0):                          # tether does not push
        return
    u = d / L
    f = k_linc * delta * (wp.float64(1.0) + stiffening * delta * delta)
    wp.atomic_add(nforce, ni, f * u)                      # nucleus pulled toward the anchor
    wp.atomic_add(aforce, ai, -f * u)


# ---------------------------------------------------------------------------------------------------
# CHROMATIN internal polymer net — Marko–Siggia WLC links (finite-extensible, small-strain soft)
# ---------------------------------------------------------------------------------------------------
@wp.kernel
def chromatin_wlc_kernel(pos: wp.array(dtype=wp.vec3d), link_a: wp.array(dtype=wp.int32),
                         link_b: wp.array(dtype=wp.int32), contour_L: wp.float64,
                         persistence: wp.float64, kbt: wp.float64,
                         force: wp.array(dtype=wp.vec3d)):
    """Chromatin internal-net WLC link: Marko–Siggia tension ``f=(kT/ℓp)[r+1/(4(1−r)²)−¼]`` (r=x/L)
    along each polymer link, pulling its two nodes together (entropic, finite-extensible → never
    interpenetrates). Governs the SOFT small-strain nuclear response (framework #6)."""
    k = wp.tid()
    ia = link_a[k]; ib = link_b[k]
    d = pos[ib] - pos[ia]
    x = wp.length(d)
    if x < wp.float64(1.0e-12):
        return
    r = x / contour_L
    if r > wp.float64(0.999999):
        r = wp.float64(0.999999)                          # guard the x→L pole
    u = d / x
    f = (kbt / persistence) * (r + wp.float64(1.0) / (wp.float64(4.0) * (wp.float64(1.0) - r) * (wp.float64(1.0) - r)) - wp.float64(0.25))
    wp.atomic_add(force, ia, f * u)                       # a pulled toward b (contractile entropic)
    wp.atomic_add(force, ib, -f * u)
