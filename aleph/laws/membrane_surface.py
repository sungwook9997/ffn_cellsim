"""Independent plasma-membrane surface — H.8 Template 2 (fully-dynamical lipid sheet), FF/Warp, µm·pN·s.

The plasma membrane as its OWN entity, not the convex-hull boundary of the cortex filaments (PI 2026-07-07:
"플라즈마 멤브레인은 필라멘트의 경계가 아닌 독립적인 개체로서 존재해야"). A separate closed triangulated
sheet (its own node set) that wraps the cortex, carries lipid-bilayer mechanics, is tethered to the cortex by
ERM linkers, and is the true osmotic envelope that CONTAINS the cortex (so escaping filament tips cannot leak
through the bilayer). This is the H.8 brief's Template-2 refinement (Template 1 = a lumped tension on cortex
nodes, which is what the model had).

Mechanics (KU-3.B1, all lit-anchored — see ``aleph/laws/compartments.py::resolve_membrane``):
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

from aleph.laws.surface_manifold import icosphere
from aleph.laws.compartments import resolve_membrane, ResolvedMembrane

# Membrane–cortex adhesion (KU-3.B1.4): γ_MCA ≈ 1e-5 J/m² = 10 pN/µm (band 1e-6–1e-4 J/m² = 1–1000 pN/µm).
GAMMA_MCA_PN_UM = 10.0
# Reservoir capacity: apparent-area strain the fold/microvillus reservoir absorbs at constant tension before the
# K_A elastic upturn (Raucher-Sheetz plateau). ⚠ MAGIC-NUMBER — SURFACED (PI 2026-07-07): MCF7 is caveolae-
# DEFICIENT (Cav-1 negative; Lavie1998/Fiucci2002), so the reservoir here is fold/microvillus/ruffle unfolding,
# NOT caveolae — and the 0.60 value is NOT yet sourced (no MCF7 fold-reservoir inventory registered). The
# caveolae compartment is retained but DEFAULT-OFF for MCF7 (enableable for caveolae-competent cell types). This
# constant must be re-derived from a registered MCF7 membrane-fold reservoir; treat runs as provisional until then.
# ⚠ SUPERSEDED for `ac/` (PI `AC_EXECUTION_PLAN` §9 D2 option B, 2026-07-28). An unsourced placeholder for the
# fold/microvillus reservoir before K_A engages. `ac/engine/membrane_area.py` replaces it with a
# REQUIRED, provenance-carrying `capacity` that is a sweep axis and an inference OUTPUT, precisely so
# this kind of silent constant cannot govern the physics. Note the `ac/` path never reached this number
# anyway — it runs the pure gamma_mem plateau with the K_A upturn default-OFF ("hard-truth #8"), so the
# real defect was zero area stiffness, not this value. Kept only for the `ff/` callers that predate the
# decision; do not use it in new `ac/` code.
RESERVOIR_STRAIN = 0.60


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


# ----------------------------------------------------------------------------------------------------------
# Helfrich bending — the DEFINING lipid-membrane mechanic (absent until now; the lumped 2γ/R tension has no
# bending term). Discrete dihedral-hinge form: per interior edge shared by the two outward-wound triangles
# T1=(i,j,k), T2=(j,i,l) with unit normals n̂₁,n̂₂,
#     E_e = κ̃·(1 − n̂₁·n̂₂),   (flat reference C0=0),
# summed over edges = a discrete Helfrich/Willmore energy (∮(κ/2)(2H)²dA → 8πκ on a sphere). Force = −∂E/∂x
# (Σf=0, flat→0; the FD-gradient test in test_membrane_bending.py is the sign arbiter). One thread per hinge,
# atomic_add, float64 — mirrors forces_warp.cytosim_bending_kernel.
#
# ⚠ CALIBRATION CORRECTION (2026-07-16, measured — see test_membrane_bending.py). The design doc claimed the
# Seung–Nelson κ̃=2κ_m/√3 gives sphere energy 8πκ. That holds only for the IDEAL regular (all-degree-6,
# equilateral) triangular lattice, whose sphere dihedral sum Σ(1−cosθ)=4√3π≈21.77. Our sheet is an ICOSPHERE
# (12 degree-5 defects; not the ideal lattice) whose Σ converges to **7.50** (resolution- AND R-independent;
# icosahedron 7.64 → subdivided 7.497…7.519). So κ̃=2κ/√3 under-represents κ_m by 21.77/7.50≈2.9×. The honest,
# grid-invariant fix is to CALIBRATE κ̃ so the discrete reference sphere reproduces the continuum 8πκ_m:
#     κ̃ = 8π·κ_m / Σ_ref      (Σ_ref = the dihedral sum of the RESTING sheet ≈ SIGMA_ICOSPHERE)
# Σ_ref is a geometric constant of the icosphere topology (not a fit to pass a gate) → the Magic-Number Block
# is satisfied (derived, resolution-invariant). ``calibrate_kappa_tilde`` measures Σ from the actual mesh.
# ----------------------------------------------------------------------------------------------------------
@wp.kernel
def helfrich_bending_kernel(pos: wp.array(dtype=wp.vec3d), hinges: wp.array(dtype=wp.int32, ndim=2),
                            kappa_tilde: wp.float64, force: wp.array(dtype=wp.vec3d)):
    """Dihedral Helfrich bending force on the four hinge vertices [i,j,k,l] (i,j = shared edge; k = flap of
    T1=(i,j,k); l = flap of T2=(j,i,l)). ``pos`` may be the GLOBAL node array and ``hinges`` GLOBAL indices —
    the kernel only touches the four vertices it names, so it composes additively with the cortex/nucleus
    forces in ``f_d`` (the membrane is a tail block). ``kappa_tilde`` = 2κ_m/√3."""
    h = wp.tid()
    i = hinges[h, 0]; j = hinges[h, 1]; k = hinges[h, 2]; l = hinges[h, 3]
    xi = pos[i]; xj = pos[j]; xk = pos[k]; xl = pos[l]
    # T1=(i,j,k): N1 = (xj−xi)×(xk−xi);  T2=(j,i,l): N2 = (xi−xj)×(xl−xj). |N| = 2·area; both point outward
    # for a consistently outward-wound sheet ⇒ n̂₁·n̂₂ ≈ 1 (E ≈ 0) when flat/convex.
    N1 = wp.cross(xj - xi, xk - xi)
    N2 = wp.cross(xi - xj, xl - xj)
    a1 = wp.length(N1); a2 = wp.length(N2)                 # a = 2·area
    if a1 < wp.float64(1.0e-18) or a2 < wp.float64(1.0e-18):
        return
    n1 = N1 / a1; n2 = N2 / a2
    c = wp.dot(n1, n2)
    q1 = n2 - c * n1                                        # ⊥-component of n̂₂ in the n̂₁ tangent frame
    q2 = n1 - c * n2
    inv1 = kappa_tilde / a1                                 # κ̃/(2A₁)
    inv2 = kappa_tilde / a2                                 # κ̃/(2A₂)
    # f = −∂E/∂x. i,j on the shared edge feel both triangles; k,l only their own (design formula, FD-verified).
    fi = -inv1 * wp.cross(xk - xj, q1) - inv2 * wp.cross(xj - xl, q2)
    fj = -inv1 * wp.cross(xi - xk, q1) - inv2 * wp.cross(xl - xi, q2)
    fk = -inv1 * wp.cross(xj - xi, q1)
    fl = -inv2 * wp.cross(xi - xj, q2)
    wp.atomic_add(force, i, fi)
    wp.atomic_add(force, j, fj)
    wp.atomic_add(force, k, fk)
    wp.atomic_add(force, l, fl)


SIGMA_ICOSPHERE = 7.50   # converged dihedral sum Σ(1−cosθ) of a closed icosphere (measured; resolution- and
#                          R-independent — icosahedron 7.64 → subdivided 7.497…7.519). Geometric constant of the
#                          icosphere topology, NOT a tuned value. See test_membrane_bending.py::sphere.


def kappa_tilde(kappa_m: float, sigma_ref: float = SIGMA_ICOSPHERE) -> float:
    """Calibrated dihedral coupling κ̃ = 8π·κ_m/Σ_ref [pN·µm] so the discrete reference sphere reproduces the
    continuum Willmore energy 8πκ_m (≠ the ideal-lattice Seung–Nelson 2κ/√3 — the icosphere Σ_ref≈7.50, not
    4√3π; see the module header). ``sigma_ref`` defaults to the icosphere constant; pass the mesh-measured Σ
    (``calibrate_kappa_tilde``) for an exact per-mesh calibration."""
    return 8.0 * np.pi * float(kappa_m) / float(sigma_ref)


def calibrate_kappa_tilde(kappa_m: float, verts: np.ndarray, hinges: np.ndarray) -> float:
    """κ̃ calibrated to the ACTUAL resting mesh: measure Σ=Σ(1−n̂₁·n̂₂) on ``verts``/``hinges`` (should be a
    near-sphere) and return κ̃=8πκ_m/Σ so ``membrane_bending_energy`` at rest equals exactly 8πκ_m. This makes
    the calibration exact for whatever subdivision is used, with no hard-coded constant in the physics path."""
    sigma = membrane_bending_energy(verts, hinges, 1.0)     # κ̃=1 → the raw dihedral sum
    if sigma <= 1e-9:
        return kappa_tilde(kappa_m)                          # degenerate/flat → fall back to the constant
    return 8.0 * np.pi * float(kappa_m) / float(sigma)


def build_membrane_hinges(faces: np.ndarray) -> np.ndarray:
    """Hinge adjacency (Nh,4)=[i,j,k,l] from DIRECTED half-edges — the winding is the correctness-critical
    detail. Each face (v0,v1,v2) contributes directed edges v0→v1 (opp v2), v1→v2 (opp v0), v2→v0 (opp v1).
    For a manifold edge {i,j} the two half-edges i→j and j→i live in T1=(i,j,k) and T2=(j,i,l) respectively,
    so the two flap vertices come out with the outward winding each triangle already carries (⇒ consistent
    normals). Closed icosphere: Nh = 3·Ntri/2 (every interior edge shared by exactly two triangles)."""
    opp = {}                                               # directed (a,b) → opposite vertex of its triangle
    for f in faces:
        v0, v1, v2 = int(f[0]), int(f[1]), int(f[2])
        opp[(v0, v1)] = v2
        opp[(v1, v2)] = v0
        opp[(v2, v0)] = v1
    hinges = []
    seen = set()
    for (a, b), k in opp.items():
        key = (a, b) if a < b else (b, a)
        if key in seen:
            continue
        if (b, a) not in opp:
            continue                                       # boundary edge (open sheet) — no hinge
        seen.add(key)
        i, j = key                                         # canonical i<j
        # T1 must contain the directed edge i→j; T2 the directed edge j→i (so normals stay outward)
        k1 = opp[(i, j)]                                   # flap of the triangle carrying i→j
        l2 = opp[(j, i)]                                   # flap of the triangle carrying j→i
        hinges.append((i, j, k1, l2))
    return np.ascontiguousarray(hinges, np.int64)


def membrane_bending_energy(verts: np.ndarray, hinges: np.ndarray, kappa_tilde_val: float) -> float:
    """Discrete Helfrich energy E = Σ_hinges κ̃·(1 − n̂₁·n̂₂) [pN·µm] (host reference for the sphere→8πκ_m and
    FD-gradient acceptance tests; the kernel force must equal −∂E/∂x)."""
    if hinges.shape[0] == 0:
        return 0.0
    xi = verts[hinges[:, 0]]; xj = verts[hinges[:, 1]]; xk = verts[hinges[:, 2]]; xl = verts[hinges[:, 3]]
    N1 = np.cross(xj - xi, xk - xi)
    N2 = np.cross(xi - xj, xl - xj)
    n1 = N1 / (np.linalg.norm(N1, axis=1, keepdims=True) + 1e-300)
    n2 = N2 / (np.linalg.norm(N2, axis=1, keepdims=True) + 1e-300)
    c = np.sum(n1 * n2, axis=1)
    return float(kappa_tilde_val * np.sum(1.0 - c))


def erm_rupture_force(kappa_m: float, gamma_mem: float, gamma_mca: float = GAMMA_MCA_PN_UM) -> float:
    """ERM tether rupture force f_t = 2π·√(2κ_m·(T_m+γ_MCA)) [pN] (KB-3.B1.4; band 5–40 pN). DERIVED from the
    bilayer bending rigidity + membrane–cortex adhesion — not a magic number. Blebbing (local detachment) is
    then EMERGENT when the tether tension exceeds this in ``erm_tether_kernel``."""
    return float(2.0 * np.pi * np.sqrt(2.0 * float(kappa_m) * (float(gamma_mem) + float(gamma_mca))))


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
           "membrane_containment_kernel", "helfrich_bending_kernel", "kappa_tilde",
           "calibrate_kappa_tilde", "build_membrane_hinges", "membrane_bending_energy",
           "erm_rupture_force", "GAMMA_MCA_PN_UM", "RESERVOIR_STRAIN", "SIGMA_ICOSPHERE"]
