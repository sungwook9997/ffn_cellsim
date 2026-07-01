"""FF cortex fiber-network assembly on a sphere (Stage 6c-b).

Builds the actin cortical shell as a :class:`~ff.fiber_network.FiberNetwork` of fibers lying on a
sphere of radius ``R`` — the FF-engine counterpart of the archived HOOMD cortex
(``archive/hoomd_legacy/cortex/cortex.py``), so the eventual MD-free γ can be compared
apples-to-apples with the BAOAB-MD γ (the γ-floor experiment, Stage 6d / ENGINE.md §4).

Every geometric/material constant is a LITERATURE-ANCHORED input carried over verbatim from the
H.3 cortex config (``configs/phase1_h3.yaml``); none is invented here:

    R_cell                = 10.0 µm     cell radius (KU-3.17; MCF7, Wagner 2011)
    n_filaments           = 1000        ×40 mesoscopic count (Plan v2 §3 H.3 v3.1 — the ONLY
                                        sanctioned coarse-graining; native ≈38 000)
    beads_per_filament    = 7           mean (= L/ℓ₀ + 1 = 6 + 1; H.3 v3.1)
    seg = ℓ₀              = 0.5 µm      segment rest length (H.3 box derivation)
    L_filament            = 3.0 µm      = (beads−1)·ℓ₀
    persistence_length    = 17.0 µm     actin ℓ_p (KU-1.1; Gittes et al. 1993)
    κ = k_B·T·ℓ_p        ≈ 0.073 pN·µm² bending modulus (ff.units.KAPPA_ACTIN)
    areal density         = n/(4πR²)    ≈ 0.80 µm⁻² (cortex.py coverage check)

Fibers are placed with COM uniformly on the sphere and a random tangent orientation, then laid as
a **great-circle arc** so every model-point sits exactly on radius ``R`` (the cortex is a thin
shell; radial drift is below the mesh scale). seg lengths follow from the arc spacing (chord
≈ ℓ₀). This is the fiber GEOMETRY only — crosslinkers / myosin (the Hand kinetic layer) are added
separately in Stage 6c-c; without them the network is a set of unconnected fibers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ffn_sim.ff import units as U
from ffn_sim.ff.fiber_network import FiberNetwork, build_fiber_network


@dataclass(frozen=True, slots=True)
class CortexParams:
    """Literature-anchored cortex assembly parameters, in FF units (µm, pN·µm²).

    Defaults reproduce the H.3 production cortex (``configs/phase1_h3.yaml``). All values are
    grounded inputs (see module docstring) — not tuned.
    """

    R_um: float = 7.5                       # MCF7 cell radius [µm] (Wagner 2011; PI 2026-07-01 — was 10µm
                                            # generic, mislabeled "MCF7"; corrected to the real MCF7 value)
    n_filaments: int = 70686                # actin areal density ~100/µm² (KB-3.18; = 100·4πR²) — supersedes
                                            # the CLAUDE.md ~38000 (=30/µm², 3× sparse vs KB; ⚠️CLAUDE.md↔KB
                                            # conflict surfaced to PI). ×40 mesoscale RETIRED, FF GPU-native
    beads_per_filament: int = 7             # mean (H.3 v3.1)
    seg_um: float = 0.5                     # ℓ₀ segment rest length [µm]
    kappa: float = U.KAPPA_ACTIN            # bending modulus [pN·µm²] = k_B·T·ℓ_p
    persistence_length_um: float = U.LP_ACTIN_UM   # actin ℓ_p [µm] (KU-1.1)

    @property
    def L_filament_um(self) -> float:
        """Contour length L = (beads−1)·ℓ₀ [µm]."""
        return (self.beads_per_filament - 1) * self.seg_um

    @property
    def surface_area_um2(self) -> float:
        """Cortex shell area 4πR² [µm²]."""
        return 4.0 * np.pi * self.R_um**2

    @property
    def areal_density_um2(self) -> float:
        """Filament areal density n/(4πR²) [µm⁻²]."""
        return self.n_filaments / self.surface_area_um2


def _random_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    """``n`` uniformly-distributed unit vectors on S² (Marsaglia: normalize Gaussians)."""
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _local_tangent_frame(com_dirs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-COM orthonormal tangent frame on S²: (ê_φ azimuthal/circumferential, ê_θ meridional/polar).

    ê_φ = ẑ × ĉ (points along lines of latitude); ê_θ = ĉ × ê_φ (points pole-to-pole). At the poles
    (ĉ ∥ ẑ) ê_φ is degenerate → filled with an arbitrary but consistent tangent so the frame stays finite.
    """
    n = com_dirs.shape[0]
    z = np.array([0.0, 0.0, 1.0])
    ephi = np.cross(np.broadcast_to(z, (n, 3)), com_dirs)
    nrm = np.linalg.norm(ephi, axis=1, keepdims=True)
    fallback = np.cross(np.broadcast_to(np.array([1.0, 0.0, 0.0]), (n, 3)), com_dirs)
    ephi = np.where(nrm > 1e-9, ephi / np.maximum(nrm, 1e-12),
                    fallback / np.maximum(np.linalg.norm(fallback, axis=1, keepdims=True), 1e-12))
    etheta = np.cross(com_dirs, ephi)
    etheta = etheta / np.maximum(np.linalg.norm(etheta, axis=1, keepdims=True), 1e-12)
    return ephi, etheta


def _tangent_field(com_dirs: np.ndarray, rng: np.random.Generator, orientation: str,
                   nematic_S: float) -> np.ndarray:
    """Fiber tangent directions for a chosen cortex ARRANGEMENT (the 'filament alignment' axis).

    orientation:
      'isotropic'       — random tangent (nematic order S≈0; the default / historical cortex).
      'circumferential' — director = ê_φ (azimuthal; filaments wrap like lines of latitude).
      'meridional'      — director = ê_θ (polar; filaments run pole-to-pole).
    ``nematic_S`` ∈ [0,1] sets the alignment strength for the non-isotropic cases (S=1 fully aligned to the
    director, S=0 recovers isotropic): t = normalize(S·director + (1−S)·random_tangent).
    """
    rand = _random_unit_vectors(com_dirs.shape[0], rng)
    if orientation == "isotropic":
        return rand
    ephi, etheta = _local_tangent_frame(com_dirs)
    if orientation == "circumferential":
        director = ephi
    elif orientation == "meridional":
        director = etheta
    else:
        raise ValueError(f"orientation must be isotropic/circumferential/meridional, got {orientation!r}")
    S = float(np.clip(nematic_S, 0.0, 1.0))
    t = S * director + (1.0 - S) * rand
    return t / np.maximum(np.linalg.norm(t, axis=1, keepdims=True), 1e-12)


def nematic_order(net: FiberNetwork, orientation: str = "circumferential") -> dict:
    """Measure the cortex filament ALIGNMENT: mean |t̂·director̂| per fiber (1=aligned, isotropic≈0.5-0.64)
    plus the global Q-tensor scalar order S (largest eigenvalue of ⟨(3/2)t⊗t − ½I⟩). Diagnostic only."""
    off = net.fiber_offsets
    tans = []
    coms = []
    for f in range(len(off) - 1):
        seg = net.pos[off[f]:off[f + 1]]
        if seg.shape[0] < 2:
            continue
        t = seg[-1] - seg[0]
        t = t / max(np.linalg.norm(t), 1e-12)
        tans.append(t)
        coms.append(seg.mean(0) / max(np.linalg.norm(seg.mean(0)), 1e-12))
    tans = np.asarray(tans)
    coms = np.asarray(coms)
    Q = (1.5 * np.einsum("fi,fj->ij", tans, tans) / len(tans)) - 0.5 * np.eye(3)
    S_global = float(np.linalg.eigvalsh(Q).max())
    ephi, etheta = _local_tangent_frame(coms)
    director = ephi if orientation == "circumferential" else etheta
    align = float(np.mean(np.abs(np.einsum("fi,fi->f", tans, director))))
    return {"align_to_director": align, "S_global_Q": S_global, "n_fibers": len(tans)}


def _great_circle_arc(com_dir: np.ndarray, tangent: np.ndarray, R: float, n_beads: int,
                      seg: float) -> np.ndarray:
    """Lay ``n_beads`` model-points as a centred great-circle arc on the sphere of radius ``R``.

    Point at arc length ``s``: ``p(s) = R[cos(s/R)·ĉ + sin(s/R)·t̂]`` with ``ĉ`` the (unit) COM
    direction and ``t̂`` a unit tangent ⊥ ĉ. Beads span ``s ∈ [−(n−1)seg/2, +(n−1)seg/2]`` →
    chord spacing ≈ ``seg``; all points lie exactly on radius ``R``.
    """
    c = com_dir / np.linalg.norm(com_dir)
    t = tangent - np.dot(tangent, c) * c            # project tangent into the plane ⊥ c
    t = t / np.linalg.norm(t)
    s = (np.arange(n_beads) - (n_beads - 1) / 2.0) * seg
    ang = s / R
    return R * (np.cos(ang)[:, None] * c[None, :] + np.sin(ang)[:, None] * t[None, :])


def build_cortex_network(params: CortexParams | None = None, *,
                         rng: np.random.Generator | None = None,
                         n_filaments: int | None = None,
                         orientation: str = "isotropic",
                         nematic_S: float = 1.0) -> tuple[FiberNetwork, dict]:
    """Assemble the cortex fiber network on the sphere (FF units, µm).

    Args:
        params: cortex parameters (defaults = H.3 production cortex).
        rng: random generator (default seed 0 for reproducibility).
        n_filaments: optional override of ``params.n_filaments`` (e.g. a small prototype).

    Returns:
        ``(net, meta)`` — the assembled :class:`FiberNetwork` (κ per fiber set from ``params``) and
        a metadata dict (areal density, segment-length stats, on-shell residual, contour length).
    """
    if params is None:
        params = CortexParams()
    if rng is None:
        rng = np.random.default_rng(0)
    F = int(n_filaments) if n_filaments is not None else params.n_filaments
    R, nb, seg = params.R_um, params.beads_per_filament, params.seg_um

    com_dirs = _random_unit_vectors(F, rng)
    tangents = _tangent_field(com_dirs, rng, orientation, nematic_S)   # ARRANGEMENT (alignment) axis
    fibers = [_great_circle_arc(com_dirs[f], tangents[f], R, nb, seg) for f in range(F)]

    net = build_fiber_network(fibers, kappa=params.kappa)

    radii = np.linalg.norm(net.pos, axis=1)
    meta = {
        "n_filaments": F,
        "n_nodes": net.n_nodes,
        "beads_per_filament": nb,
        "R_um": R,
        "L_filament_um": params.L_filament_um,
        "areal_density_um2": F / params.surface_area_um2,
        "seg_len_mean_um": float(net.seg_rest.mean()),
        "seg_len_std_um": float(net.seg_rest.std()),
        "on_shell_residual_um": float(np.max(np.abs(radii - R))),
        "kappa_pN_um2": params.kappa,
        "orientation": orientation,
        "nematic_S": float(nematic_S) if orientation != "isotropic" else 0.0,
    }
    return net, meta


def equatorial_circumference_um(params: CortexParams) -> float:
    """Equatorial circumference 2πR [µm] — the denominator for the method-of-planes γ (Stage 6d)."""
    return 2.0 * np.pi * params.R_um
