r"""Angle-harmonic bending for the NMII minifilament — the F6 transmission fix (Warp CUDA source).

Round-2 motor-fix (``ac/mag-f-motor-fix``). NG-1 (two-filament isometric stall) FAILED on the A5000 with
transmission |B|/|A| ≈ 0.61 and per-head tangential load ≈ 0.0004 pN (should → f_stall): the crossbridge force
was borne TRANSVERSE to ``walk_dir`` and never reached the opposite filament. Root cause (F6):

  * the head↔backbone ARM is a single central *distance* spring — it has ZERO transverse stiffness, so a head
    pulled along ``walk_dir`` (the actin/x axis) simply SWINGS the arm about the backbone bead instead of
    loading it; and
  * the backbone is *distance* springs only — no bending term — so the rod BOWS under the collected arm loads
    instead of holding straight and transmitting the axial tension collinearly.

Both are the classic "harmonic bond has no angular stiffness" failure. The mechanistic fix (CLAUDE.md PI rule
— the SAME angle-harmonic-with-thermal choice held up for the Arp2/3 branch, over a rigid/lumped proxy) adds a
harmonic-angle bending potential ``E = ½ k_θ (θ − θ₀)²`` on two triple families:

  1. BACKBONE bending — consecutive backbone triples ``(i−1, i, i+1)``, rest ``θ₀ = π`` (straight rigid rod);
  2. HEAD-ARM orientation — triples ``(backbone_neighbor, backbone_bead, head)``, rest ``θ₀ = π/2`` (the lever
     arm stands perpendicular to the backbone). This is the "angle harmonic" head-arm option of the two the
     Round-2 boot allows (the other being a 2nd arm anchor); it is the more mechanistic — it represents the
     real lever-arm bending rigidity, not a geometric truss.

With both terms the head can no longer swing freely: pulling it along ``walk_dir`` now rotates the arm against
``k_θ,arm`` (first-order transverse stiffness ``k_θ,arm / r0_head²``), the load rises with the walked abscissa
toward ``f_stall`` (the Hill self-limit), and the straight rod carries that tension to the antiparallel side —
transmission → 1 (the collinear-static ideal ``two_filament_reference.ideal_transmission_ratio`` proves).

⭐ NOT tuned to a force (Magic-Number-Block, §"No empirical magic numbers"). Both stiffnesses are DERIVED from a
grid-invariant material property in the RIGID-ROD limit, never chosen to hit a transmission band:

  * ``k_θ,backbone = κ_bb / a``  (the standard bead-angle → continuum map: ``E_bend = ∫ ½κ C² ds`` with joint
    angle ``θ ≈ a·C`` gives ``½ k_θ θ² = ½ (κ/a) θ²``), with the flexural rigidity ``κ_bb = L_p · k_BT`` from a
    backbone persistence length ``L_p`` and ``a`` the backbone-bead spacing. ``L_p`` is an unresolved
    mature-minifilament material GAP (PI/KB). Single-coiled-coil values are diagnostic proxies only; production
    bending remains blocked until an applicable source and paralog/topology contract are ratified.
  * ``k_θ,arm = k_xb · r0_head²``  — the lever arm's transverse stiffness ``k_θ,arm / r0_head²`` is set to the
    crossbridge stiffness ``k_xb`` it must transmit (rigid-lever limit, structural-consistency), so the head
    loads instead of swings. Also GAP-flagged and convergence-checked.

CFL: bending adds an effective LINEAR stiffness (``4 k_θ / a²`` at a backbone bead, ``k_θ / r0_head²`` at a
head); the composed engine must fold these into ``kmax`` when picking ``dt_mu`` (see
:func:`backbone_bending_cfl_stiffness` / :func:`arm_bending_cfl_stiffness`, wired through
``MyosinForce.cfl_stiffness``).

⚠ Runtime: NVIDIA Warp on CUDA GPU only (I0-A). Authored on the dev Mac (no CUDA) but NOT launched here — the
transmission fix is proven by the native NG-1 gate on the lead's gbook A5000. The :func:`angle_harmonic_forces`
NumPy mirror is the host acceptance reference (``tests/ac/motor/test_backbone_angle_mechanics.py``) — it is the
BIT-FOR-FORMULA twin of :func:`angle_harmonic_kernel`, so the host gate certifies the exact mechanics the device
kernel evaluates (zero force/torque at rest, restoring moment under bend, first-order transverse stiffness).

References: Landau-Lifshitz elastic rod (κ = EI); bead-angle discretisation (κ = k_θ·a); LAMMPS ``angle_harmonic``
force expansion; Billington 2013 (stiff bipolar backbone); Stam-Hocky 2015. kBT = 4.28e-3 pN·µm (= NMII_F0 kBT).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import warp as wp

# engine units: length µm, force pN, energy pN·µm, angle rad, k_theta pN·µm/rad²
KBT_PN_UM: float = 4.28e-3        # thermal energy [pN·µm] (= NMII_F0 numerator kBT; Bell f0 = kBT/x_beta)
_ANGLE_SIN_FLOOR = 1.0e-8         # sinθ floor near the straight/degenerate configuration (LAMMPS SMALL)

__all__ = [
    "KBT_PN_UM",
    "angle_harmonic_kernel",
    "angle_harmonic_forces",
    "backbone_bending_kappa",
    "backbone_bending_k_theta",
    "arm_orientation_k_theta",
    "backbone_bending_cfl_stiffness",
    "arm_bending_cfl_stiffness",
]


# ── Warp device kernel (float64) ─────────────────────────────────────────────────────────────────
@wp.kernel
def angle_harmonic_kernel(
    pos: wp.array(dtype=wp.vec3d),
    triples: wp.array(dtype=wp.int32, ndim=2),
    k_theta: wp.float64,
    theta0: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    r"""Harmonic-angle bending ``E = ½ k_θ (θ − θ₀)²`` on each triple ``triples[t] = (i, j, k)``, vertex ``j``.

    ``θ`` is the angle at ``j`` between ``r1 = pos[i] − pos[j]`` and ``r2 = pos[k] − pos[j]``. Adds the
    RESTORING force ``F = −∂E/∂x`` to the three particles (Newton's 3rd law: ``f_i + f_j + f_k = 0`` exactly).
    The LAMMPS ``angle_harmonic`` force expansion (with the ½-prefactor convention here, ``a = −k_θ Δθ / sinθ``):

        a11 = a·cosθ / |r1|²,  a12 = −a / (|r1||r2|),  a22 = a·cosθ / |r2|²
        f_i = a11·r1 + a12·r2,  f_k = a22·r2 + a12·r1,  f_j = −(f_i + f_k)

    Zero force at the rest angle (Δθ = 0 ⇒ a = 0). The ``sinθ`` floor guards the exactly-straight (θ = π) and
    exactly-folded (θ = 0) configurations where the ratio ``Δθ / sinθ`` is finite but 0/0 in floating point.
    """
    t = wp.tid()
    i = triples[t, 0]
    j = triples[t, 1]
    k = triples[t, 2]
    r1 = pos[i] - pos[j]
    r2 = pos[k] - pos[j]
    n1 = wp.length(r1)
    n2 = wp.length(r2)
    if n1 < wp.float64(1.0e-12) or n2 < wp.float64(1.0e-12):
        return
    c = wp.dot(r1, r2) / (n1 * n2)
    c = wp.clamp(c, wp.float64(-1.0), wp.float64(1.0))
    s = wp.sqrt(wp.float64(1.0) - c * c)
    if s < wp.float64(_ANGLE_SIN_FLOOR):
        s = wp.float64(_ANGLE_SIN_FLOOR)
    theta = wp.acos(c)
    dtheta = theta - theta0
    a = -k_theta * dtheta / s
    a11 = a * c / (n1 * n1)
    a12 = -a / (n1 * n2)
    a22 = a * c / (n2 * n2)
    f_i = a11 * r1 + a12 * r2
    f_k = a22 * r2 + a12 * r1
    wp.atomic_add(force, i, f_i)
    wp.atomic_add(force, k, f_k)
    wp.atomic_add(force, j, -(f_i + f_k))


# ── NumPy mirror — the host acceptance reference (BIT-FOR-FORMULA twin of the kernel) ─────────────
def angle_harmonic_forces(
    pos: npt.NDArray[np.float64],
    triples: npt.NDArray[np.int64],
    k_theta: float,
    theta0: float,
) -> npt.NDArray[np.float64]:
    r"""NumPy mirror of :func:`angle_harmonic_kernel`: return the per-particle bending force ``(n, 3)`` [pN].

    Same energy ``E = ½ k_θ (θ − θ₀)²`` and the same LAMMPS force expansion the device kernel evaluates — the
    host gate certifies the mechanics without launching a CUDA kernel (I0-A dev-Mac discipline).

    Args:
        pos: particle positions ``(n, 3)`` [µm].
        triples: angle triples ``(m, 3)`` as ``(i, j, k)`` with ``j`` the vertex.
        k_theta: angular stiffness [pN·µm/rad²].
        theta0: rest angle [rad].

    Returns:
        Per-particle force ``(n, 3)`` [pN] (the restoring ``−∂E/∂x``, summed over all triples).
    """
    pos = np.asarray(pos, dtype=np.float64)
    triples = np.asarray(triples, dtype=np.int64)
    force = np.zeros_like(pos)
    for i, j, k in triples:
        r1 = pos[i] - pos[j]
        r2 = pos[k] - pos[j]
        n1 = float(np.linalg.norm(r1))
        n2 = float(np.linalg.norm(r2))
        if n1 < 1.0e-12 or n2 < 1.0e-12:
            continue
        c = float(np.dot(r1, r2)) / (n1 * n2)
        c = min(1.0, max(-1.0, c))
        s = np.sqrt(1.0 - c * c)
        if s < _ANGLE_SIN_FLOOR:
            s = _ANGLE_SIN_FLOOR
        theta = np.arccos(c)
        dtheta = theta - theta0
        a = -k_theta * dtheta / s
        a11 = a * c / (n1 * n1)
        a12 = -a / (n1 * n2)
        a22 = a * c / (n2 * n2)
        f_i = a11 * r1 + a12 * r2
        f_k = a22 * r2 + a12 * r1
        force[i] += f_i
        force[k] += f_k
        force[j] += -(f_i + f_k)
    return force


# ── k_θ derivation (Magic-Number-Block: DERIVED, grid-invariant, convergence-checked; not fit) ────
def backbone_bending_kappa(persistence_length_um: float, kbt_pn_um: float = KBT_PN_UM) -> float:
    r"""Backbone flexural rigidity ``κ = L_p · k_BT`` [pN·µm²] from a persistence length (grid-invariant).

    Args:
        persistence_length_um: backbone persistence length ``L_p`` [µm] (GAP → PI/KB). A value from an
            isolated coiled coil is not automatically applicable to a mature multi-tail bipolar minifilament;
            the caller must supply a source-ratified material contract.
        kbt_pn_um: thermal energy [pN·µm].

    Returns:
        Flexural rigidity ``κ`` [pN·µm²].
    """
    if persistence_length_um <= 0.0:
        raise ValueError("persistence_length_um must be positive")
    return persistence_length_um * kbt_pn_um


def backbone_bending_k_theta(kappa_pn_um2: float, segment_len_um: float) -> float:
    r"""Discrete backbone angle stiffness ``k_θ = κ / a`` [pN·µm/rad²] (bead-angle → continuum map).

    ``E_bend = ∫ ½ κ C² ds`` with a joint angle ``θ ≈ a·C`` (``a`` = bead spacing) discretises to
    ``Σ ½ (κ/a) θ²`` — so ``k_θ = κ / a``. This is grid-invariant: refine the backbone (smaller ``a``) and
    ``k_θ`` rises to keep the physical rigidity ``κ`` fixed.

    Args:
        kappa_pn_um2: flexural rigidity ``κ`` [pN·µm²] (:func:`backbone_bending_kappa`).
        segment_len_um: backbone-bead spacing ``a = L_bb/(n_bb−1)`` [µm].

    Returns:
        Angular stiffness ``k_θ`` [pN·µm/rad²].
    """
    if kappa_pn_um2 <= 0.0 or segment_len_um <= 0.0:
        raise ValueError("kappa and segment_len must be positive")
    return kappa_pn_um2 / segment_len_um


def arm_orientation_k_theta(k_xb_pn_um: float, r0_head_um: float) -> float:
    r"""Head-arm angular stiffness ``k_θ,arm = k_xb · r0_head²`` [pN·µm/rad²] (rigid-lever limit).

    Sets the arm's effective TRANSVERSE stiffness at the head ``k_θ,arm / r0_head² = k_xb`` — i.e. the lever
    arm is as rigid as the crossbridge it transmits, so a head pulled along ``walk_dir`` LOADS the backbone
    instead of swinging (structural-consistency, not tuned to a transmission band). GAP-flagged (rides on the
    ``k_xb`` MASTER-knob GAP + the ``r0_head`` geometry GAP); convergence-checked at the native gate.

    Args:
        k_xb_pn_um: crossbridge stiffness [pN/µm].
        r0_head_um: head↔backbone arm rest length [µm].

    Returns:
        Angular stiffness ``k_θ,arm`` [pN·µm/rad²].
    """
    if k_xb_pn_um <= 0.0 or r0_head_um <= 0.0:
        raise ValueError("k_xb and r0_head must be positive")
    return k_xb_pn_um * r0_head_um * r0_head_um


def backbone_bending_cfl_stiffness(k_theta_bb: float, segment_len_um: float) -> float:
    r"""Effective LINEAR stiffness a backbone-bending joint adds at the middle bead ``= 4 k_θ / a²`` [pN/µm].

    Displacing the vertex bead transversely by ``δ`` gives ``Δθ ≈ 2δ/a`` ⇒ ``E ≈ (2 k_θ/a²) δ²`` ⇒ stiffness
    ``4 k_θ/a²``. Folded into ``kmax`` so ``dt_mu = 0.1/kmax`` stays stable (CFL note per the Round-2 boot).

    Args:
        k_theta_bb: backbone angular stiffness [pN·µm/rad²].
        segment_len_um: backbone-bead spacing [µm].

    Returns:
        Effective linear stiffness [pN/µm].
    """
    return 4.0 * k_theta_bb / (segment_len_um * segment_len_um)


def arm_bending_cfl_stiffness(k_theta_arm: float, r0_head_um: float) -> float:
    r"""Effective LINEAR stiffness a head-arm angle adds at the head ``= k_θ,arm / r0_head²`` [pN/µm].

    Displacing the head transverse to the arm by ``δ`` gives ``Δθ ≈ δ/r0_head`` ⇒ stiffness ``k_θ,arm/r0_head²``.

    Args:
        k_theta_arm: arm angular stiffness [pN·µm/rad²].
        r0_head_um: arm rest length [µm].

    Returns:
        Effective linear stiffness [pN/µm].
    """
    return k_theta_arm / (r0_head_um * r0_head_um)
